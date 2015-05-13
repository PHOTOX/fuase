import os
import tempfile

from ase.io import write
import ase.parallel as parallel


def view(atoms, data=None, viewer='ase-gui', repeat=None, block=False):
    # Ignore for parallel calculations:
    if parallel.size != 1:
        return

    vwr = viewer.lower()
    
    if vwr == 'ase-gui':
        format = 'traj'
        if repeat is None:
            command = 'ase-gui'
        else:
            command = 'ase-gui --repeat=%d,%d,%d' % tuple(repeat)
            repeat = None
    elif vwr == 'vmd':
        format = 'cube'
        command = 'vmd'
    elif vwr == 'rasmol':
        format = 'pdb'
        command = 'rasmol -pdb'
    elif vwr == 'xmakemol':
        format = 'xyz'
        command = 'xmakemol -f'
    elif vwr == 'gopenmol':
        format = 'xyz'
        command = 'rungOpenMol'
    elif vwr == 'avogadro':
        format = 'cube'
        command = 'avogadro'
    elif vwr == 'sage':
        from ase.visualize.sage import view_sage_jmol
        view_sage_jmol(atoms)
        return
    else:
        raise RuntimeError('Unknown viewer: ' + viewer)

    fd, filename = tempfile.mkstemp('.' + format, 'ase-')
    fd = os.fdopen(fd, 'w')
    if repeat is not None:
        atoms = atoms.repeat()
    if data is None:
        write(fd, atoms, format=format)
    else:
        write(fd, atoms, format=format, data=data)
    fd.close()
    if block:
        os.system('%s %s' % (command, filename))
        os.remove(filename)
    else:
        if os.name in ['ce', 'nt']:  # Win
            # XXX: how to make it non-blocking?
            os.system('%s %s' % (command, filename))
            os.remove(filename)
        else:
            os.system('%s %s & ' % (command, filename))
            os.system('(sleep 60; rm %s) &' % filename)
