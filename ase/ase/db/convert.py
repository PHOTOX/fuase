import os
import sys

from ase.db import connect
from ase.db.sqlite import index_statements


def convert(name):
    con1 = connect(name, use_lock_file=False)
    con1._allow_reading_old_format = True
    newname = name[:-2] + 'new.db'
    con2 = connect(newname, create_indices=False, use_lock_file=False)
    for dct in con1.select():
        keywords = dct.get('keywords', [])
        kvp = dct.get('key_value_pairs', {})
        con2.write(dct, keywords, data=dct.get('data'), **kvp)
        
    c = con2._connect()
    for statement in index_statements.split(';'):
        c.execute(statement)
    c.commit()

    os.rename(name, name[:-2] + 'old.db')
    os.rename(newname, name)
    
        
if __name__ == '__main__':
    for name in sys.argv[1:]:
        convert(name)
