from __future__ import print_function
import sys
import optparse

import ase.io
from ase.db import connect
from ase.db.summary import Summary
from ase.db.table import Table, all_columns
from ase.calculators.calculator import get_calculator

import numpy as np


def plural(n, word):
    if n == 1:
        return '1 ' + word
    return '%d %ss' % (n, word)

    
description = """Selecton is a comma-separated list of
selections where each selection is of the type "ID", "keyword" or
"key=value".  Instead of "=", one can also use "<", "<=", ">=", ">"
and  "!=" (these must be protected from the shell by using quotes).
Special keys: id, user, calculator, age, natoms, energy, magmom,
and charge.  Chemical symbols can also be used to select number of
specific atomic species (H, He, Li, ...)."""

examples = ['calculator=nwchem',
            'age<1d',
            'natoms=1',
            'user=alice',
            '2.2<bandgap<4.1',
            'Cu>=10']


def main(args=sys.argv[1:]):
    if isinstance(args, str):
        args = args.split(' ')
    parser = optparse.OptionParser(
        usage='Usage: %prog db-name [selection] [options]',
        description=description,
        epilog='Selection examples: ' + ', '.join(examples) + '.')
    
    add = parser.add_option
    add('-v', '--verbose', action='store_true', default=False)
    add('-q', '--quiet', action='store_true', default=False)
    add('-n', '--count', action='store_true',
        help='Count number of selected rows.')
    add('-l', '--long', action='store_true',
        help='Long description of selected row')
    add('-i', '--insert-into', metavar='db-name',
        help='Insert selected rows into another database.')
    add('-a', '--add-from-file', metavar='[type:]filename',
        help='Add results from file.')
    add('-k', '--add-keywords', metavar='word1,word2,...',
        help='Add keywords to selected rows.  Keywords can only contain ' +
        'letters, numbers and the underscore character and the first ' +
        'character can not be a number.')
    add('-K', '--add-key-value-pairs', metavar='key1=val1,key2=val2,...',
        help='Add key-value pairs to selected rows.  Values must be numbers ' +
        'or strings and keys must follow the same rules as keywords.')
    add('--limit', type=int, default=500, metavar='N',
        help='Show only first N rows (default is 500 rows).  Use --limit=0 ' +
        'to show all.')
    add('--delete', action='store_true',
        help='Delete selected rows.')
    add('--delete-keywords', metavar='keyword1,keyword2,...',
        help='Delete keywords for selected rows.')
    add('--delete-key-value-pairs', metavar='key1,key2,...',
        help='Delete key-value pairs for selected rows.')
    add('-y', '--yes', action='store_true',
        help='Say yes.')
    add('--explain', action='store_true',
        help='Explain query plan.')
    add('-c', '--columns', metavar='col1,col2,...',
        help='Specify columns to show.  Precede the column specification ' +
        'with a "+" in order to add columns to the default set of columns.  ' +
        'Precede by a "-" to remove columns.')
    add('-s', '--sort', metavar='column', default='id',
        help='Sort rows using column.  Default is to sort after ID.')
    add('--cut', type=int, default=35, help='Cut keywords and key-value ' +
        'columns after CUT characters.  Use --cut=0 to disable cutting. ' +
        'Default is 35 characters')
    add('-p', '--python-expression', metavar='expression',
        help='Examples: "id,energy", "id,mykey".')
    add('--csv', action='store_true',
        help='Write comma-separated-values file.')
    add('-w', '--open-web-browser', action='store_true',
        help='Open results in web-browser.')

    opts, args = parser.parse_args(args)

    if not args:
        parser.error('No database given')

    verbosity = 1 - opts.quiet + opts.verbose

    try:
        run(opts, args, verbosity)
    except Exception as x:
        if verbosity < 2:
            print('{0}: {1}'.format(x.__class__.__name__, x))
            sys.exit(1)
        else:
            raise

        
def run(opts, args, verbosity):
    filename = args.pop(0)
    query = ','.join(args)

    if query.isdigit():
        query = int(query)
    
    if opts.add_keywords:
        add_keywords = opts.add_keywords.split(',')
    else:
        add_keywords = []

    if opts.delete_keywords:
        delete_keywords = opts.delete_keywords.split(',')
    else:
        delete_keywords = []

    add_key_value_pairs = {}
    if opts.add_key_value_pairs:
        for pair in opts.add_key_value_pairs.split(','):
            key, value = pair.split('=')
            for type in [int, float]:
                try:
                    value = type(value)
                except ValueError:
                    pass
                else:
                    break
            add_key_value_pairs[key] = value

    if opts.delete_key_value_pairs:
        delete_key_value_pairs = opts.delete_key_value_pairs.split(',')
    else:
        delete_key_value_pairs = []

    con = connect(filename)
    
    def out(*args):
        if verbosity > 0:
            print(*args)
            
    if opts.add_from_file:
        filename = opts.add_from_file
        if ':' in filename:
            calculator_name, filename = filename.split(':')
            atoms = get_calculator(calculator_name)(filename).get_atoms()
        else:
            atoms = ase.io.read(filename)
        con.write(atoms, add_keywords, key_value_pairs=add_key_value_pairs)
        out('Added {0} from {1}'.format(atoms.get_chemical_formula(),
                                        filename))
        return
        
    if opts.count:
        n = 0
        for dct in con.select(query):
            n += 1
        print('%s' % plural(n, 'row'))
        return

    if opts.explain:
        for dct in con.select(query, explain=True,
                              verbosity=verbosity, limit=opts.limit):
            print(dct['explain'])
        return

    if opts.insert_into:
        con2 = connect(opts.insert_into)
        nkw = 0
        nkvp = 0
        nrows = 0
        for dct in con.select(query):
            keywords = dct.get('keywords', [])
            for keyword in add_keywords:
                if keyword not in keywords:
                    keywords.append(keyword)
                    nkw += 1

            kvp = dct.get('key_value_pairs', {})
            nkvp = -len(kvp)
            kvp.update(add_key_value_pairs)
            nkvp += len(kvp)
            con2.write(dct, keywords, data=dct.get('data'), **kvp)
            nrows += 1
            
        out('Added %s and %s (%s updated)' %
            (plural(nkw, 'keyword'),
             plural(nkvp, 'key-value pair'),
             plural(len(add_key_value_pairs) * nrows - nkvp, 'pair')))
        out('Inserted %s' % plural(nrows, 'row'))
        return

    if add_keywords or add_key_value_pairs:
        ids = [dct['id'] for dct in con.select(query)]
        nkw, nkv = con.update(ids, add_keywords, **add_key_value_pairs)
        out('Added %s and %s (%s updated)' %
            (plural(nkw, 'keyword'),
             plural(nkv, 'key-value pair'),
             plural(len(add_key_value_pairs) * len(ids) - nkv, 'pair')))
        return

    if opts.delete:
        ids = [dct['id'] for dct in con.select(query)]
        if ids and not opts.yes:
            msg = 'Delete %s? (yes/No): ' % plural(len(ids), 'row')
            if raw_input(msg).lower() != 'yes':
                return
        con.delete(ids)
        out('Deleted %s' % plural(len(ids), 'row'))
        return

    if delete_keywords or delete_key_value_pairs:
        ids = [dct['id'] for dct in con.select(query)]
        nkw, nkv = con.delete_keywords_and_key_value_pairs(
            ids, delete_keywords, delete_key_value_pairs)
        print('Removed %s and %s' %
              (plural(nkw, 'keyword'),
               plural(nkv, 'key-value pair')))
        return

    if opts.python_expression:
        for dct in con.select(query):
            row = eval(opts.python_expression, dct)
            if not isinstance(row, (list, tuple, np.ndarray)):
                row = [row]
            print(', '.join(str(x) for x in row))
        return

    if opts.long:
        dct = con.get(query)
        summary = Summary(dct)
        summary.write()
    else:
        if opts.open_web_browser:
            import ase.db.app as app
            app.connection = con
            app.app.run(host='0.0.0.0', debug=True)
        else:
            columns = list(all_columns)
            c = opts.columns
            if c:
                if c[0] == '+':
                    c = c[1:]
                elif c[0] != '-':
                    columns = []
                for col in c.split(','):
                    if col[0] == '-':
                        columns.remove(col[1:])
                    else:
                        columns.append(col.lstrip('+'))
        
            table = Table(con, verbosity, opts.cut)
            table.select(query, columns, opts.sort, opts.limit)
            if opts.csv:
                table.write_csv()
            else:
                table.write()
