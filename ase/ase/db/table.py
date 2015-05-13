from __future__ import print_function

from ase.data import atomic_masses, chemical_symbols
from ase.db.core import float_to_time_string, now, dict2constraint

import numpy as np


all_columns = ['id', 'age', 'user', 'formula', 'calculator',
               'energy', 'fmax', 'pbc', 'volume',
               #'keywords', 'keys',
               'charge', 'mass', 'smax', 'magmom']


def plural(n, word):
    if n == 1:
        return '1 ' + word
    return '%d %ss' % (n, word)

    
def cut(txt, length):
    if len(txt) <= length or length == 0:
        return txt
    return txt[:length - 3] + '...'


def cutlist(lst, length):
    if len(lst) <= length or length == 0:
        return lst
    return lst[:9] + ['... ({0} more)'.format(len(lst) - 9)]

    
def hill(numbers):
    d = {}
    for Z in numbers:
        s = chemical_symbols[Z]
        d[s] = d.get(s, 0) + 1
    result = [(s, d.pop(s)) for s in 'CH' if s in d]
    result += [(s, d[s]) for s in sorted(d)]
    return ''.join('{0}{1}'.format(symbol, n) if n > 1 else symbol
                   for symbol, n in result)
    
    
def dict2forces(d):
    forces = d.get('forces')
    if forces is None:
        return None
        
    constraints = [dict2constraint(c) for c in d.get('constraints', [])]
    if constraints:
        forces = forces.copy()
        for constraint in constraints:
            constraint.adjust_forces(d.positions, forces)
            
    return forces

    
class Table:
    def __init__(self, connection, verbosity=1, cut=35):
        self.connection = connection
        self.verbosity = verbosity
        self.cut = cut
        self.rows = []
        self.columns = None
        self.id = None
        self.right = None
        self.keys = None
        self.keywords = None
        
    def select(self, query, columns, sort, limit):
        self.limit = limit
        
        if sort != 'id':
            limit = 0
           
        self.rows = [Row(d, columns, self.cut)
                     for d in self.connection.select(
                         query, verbosity=self.verbosity, limit=limit)]

        delete = set(range(len(columns)))
        for row in self.rows:
            for n in delete.copy():
                if row.values[n] is not None:
                    delete.remove(n)
        delete = sorted(delete, reverse=True)
        for row in self.rows:
            for n in delete:
                del row.values[n]
        for n in delete:
            del columns[n]
            
        self.columns = columns
        
        if sort != 'id':
            reverse = sort[0] == '-'
            n = columns.index(sort.lstrip('-'))
            
            def key(row):
                return row.values[n]
                
            self.rows = sorted(self.rows, key=key, reverse=reverse)
            
            if self.limit:
                self.rows = self.rows[:self.limit]
                
    def format(self, subscript=None):
        right = set()
        allkeys = set()
        allkeywords = set()
        for row in self.rows:
            numbers = row.format(self.columns, subscript)
            right.update(numbers)
            allkeys.update(row.dct.key_value_pairs)
            allkeywords.update(row.dct.keywords)
            
        right.add('age')
        self.right = [column in right for column in self.columns]
        
        self.keys = sorted(allkeys)
        self.keywords = sorted(allkeywords)

    def write(self):
        self.format()
        L = [[len(s) for s in row.strings]
             for row in self.rows]
        L.append([len(c) for c in self.columns])
        N = np.max(L, axis=0)

        fmt = '{0:{align}{width}}'
        print('|'.join(fmt.format(c, align='<>'[a], width=w)
                       for c, a, w in zip(self.columns, self.right, N)))
        for row in self.rows:
            print('|'.join(fmt.format(c, align='<>'[a], width=w)
                           for c, a, w in
                           zip(row.strings, self.right, N)))

        if self.verbosity == 0:
            return
            
        print('Rows:', len(self.rows), end='')
        if self.limit:
            print(' (limited to first {0})'.format(self.limit))
        else:
            print()

        if self.keys:
            print('Keys:', ', '.join(cutlist(self.keys, self.cut)))
        if self.keywords:
            print('Keywords:', ', '.join(cutlist(self.keywords, self.cut)))
            
    def write_csv(self):
        print(', '.join(self.columns))
        for row in self.rows:
            print(', '.join(str(val) for val in row.values))        

            
class Row:
    def __init__(self, dct, columns, cut=40):
        self.dct = dct
        self.cut = cut
        self.values = None
        self.strings = None
        self.more = False
        self.set_columns(columns)
        if 'key_value_pairs' not in dct:
            dct['key_value_pairs'] = {}
        if 'keywords' not in dct:
            dct['keywords'] = []
        
    def set_columns(self, columns):
        self.values = []
        for c in columns:
            f = getattr(self, c, None)
            if f is None:
                value = getattr(self.dct, c, None)
            else:
                try:
                    value = f(self.dct)
                except (AttributeError, TypeError):
                    value = None
            self.values.append(value)
            
    def toggle(self):
        self.more = not self.more
        
    def format(self, columns, subscript=None):
        self.strings = []
        numbers = set()
        for value, column in zip(self.values, columns):
            if column == 'formula' and subscript:
                value = subscript.sub(r'<sub>\1</sub>', value)
            elif isinstance(value, int):
                value = str(value)
                numbers.add(column)
            elif isinstance(value, float):
                numbers.add(column)
                value = '{0:.3f}'.format(value)
            elif value is None:
                value = ''
            self.strings.append(value)
        
        return numbers
        
    def age(self, d):
        return float_to_time_string(now() - d.ctime)

    def formula(self, d):
        return hill(d.numbers)

    def volume(self, d):
        return abs(np.linalg.det(d.cell))

    def pbc(self, d):
        return ''.join('-P'[p] for p in d.pbc)

    def fmax(self, d):
        forces = dict2forces(d)
        return (forces**2).sum(1).max()**0.5

    def keywords(self, d):
        return cut(','.join(d.keywords), self.cut)

    def keys(self, d):
        return cut(','.join(['{0}={1}'.format(*item)
                             for item in d.key_value_pairs.items()]), self.cut)

    def mass(self, d):
        if 'masses' in d:
            return d.masses.sum()
        return atomic_masses[d.numbers].sum()

    def smax(self, d):
        return (d.stress**2).max()**0.5
