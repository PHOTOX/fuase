from cStringIO import StringIO
from ase.io.xyz import read_xyz
from ase.units import Bohr


def read_orca(filename):
    """Method to read geometry from a ORCA output
    """

    f = filename
    if isinstance(filename, str):
        f = open(filename)

    lines = f.readlines()

    done = False
    i = 0
    xyzstring = ''

    for l, line in enumerate(lines):
        if line.find('CARTESIAN COORDINATES (ANGSTROEM)') >= 0:
            i += 1
            while not done:
                i += 1

                if not lines[l+i] == '\n':
                    xyzstring += lines[l+i]
                    sym = lines[l+i].strip().split()[0]
                else:
                    done = True
        if done:
            break

    xyzstring = str(i-2) + '\n\n' + xyzstring
    atoms = read_xyz(StringIO(xyzstring))

    if type(filename) == str:
        f.close()

    return atoms


def read_orca_input(filename):
    """Method to read geometry from an ORCA input file."""
    f = filename
    if isinstance(filename, str):
        f = open(filename)
    lines = f.readlines()

    # Find geometry region of input file.

    done = False
    i = 0
    xyzstring = ''

    for l, line in enumerate(lines):
        if line.find('xyz') > -1 and line.find('*') > -1:
            while not done:
                i += 1

                if not (lines[l+i].find('*') > -1):
                    xyzstring += lines[l+i]
                    sym = lines[l+i].strip().split()[0]
                else:
                    done = True
        if done:
            break

    xyzstring = str(i-1) + '\n\n' + xyzstring
    atoms = read_xyz(StringIO(xyzstring))

    if type(filename) == str:
        f.close()

    return atoms


def write_orca(filename, atoms, charge=0, multiplicity=1):
    """Method to write ORCA coord file
    """

    if isinstance(filename, str):
        f = open(filename, 'w')
    else:  # Assume it's a 'file-like object'
        f = filename

    f.write('*xyz %i %i\n' %(charge, multiplicity))

#FUDO| geometry string should contain charge and multiplicity; should we also try to get it from the atoms object if it doesn't exist?
#    else:
#        f.write('C1\n')

#FU| we'll only take integer nuclear charges here and put external charges in $FFDATA
    for a, atom in enumerate(atoms):
        if atom.tag == -71:  # 71 is ascii G (Ghost)
            symbol = 'bq' + '0'
        else:
            symbol = atom.symbol 
        f.write('  ' + symbol + ' ' +
                str(atom.position[0]) + ' ' +
                str(atom.position[1]) + ' ' +
                str(atom.position[2]) + '\n')
    f.write('*\n')
