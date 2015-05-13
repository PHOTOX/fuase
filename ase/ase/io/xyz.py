from ase.atoms import Atoms
from ase.parallel import paropen


def read_xyz(fileobj, index=-1):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    if not isinstance(index, int) and not isinstance(index, slice):
        raise TypeError('Index argument is neither slice nor integer!')

    fileobj.seek(0)
    natoms = int(fileobj.readline())

    for i, l in enumerate(fileobj):
        pass

    ln = i + 2
    lnsnp = natoms + 2
    lastsnap = ln // lnsnp

    rvrs = False

    if isinstance(index, int):
        if index < 0:
            tmpsnp = lastsnap + index
            trbl = range(tmpsnp, tmpsnp + 1, 1)
        else:
            trbl = range(index, index + 1, 1)
        rtnndx = -1
    elif isinstance(index, slice):
        start = index.start
        stop = index.stop
        step = index.step

        if start is None:
            start = 0
        elif start < 0:
            start = lastsnap + start

        if step is None:
            step = 1

        if stop is None:
            stop = lastsnap
        elif stop < 0:
            stop = lastsnap + stop

        trbl = range(start, stop, step)

        if step < 0:
            rvrs = True
            trbl.reverse()

        rtnndx = slice(len(trbl))

    images = []
    current = 0

    fileobj.seek(0)

    for index in trbl:
        for lnndx in range(current, index * lnsnp):
            line = fileobj.readline()
        
        line = fileobj.readline()
        line = fileobj.readline()

        symbols = []
        positions = []
        for ln in range(natoms):
            line = fileobj.readline()

            symbol, x, y, z = line.split()[:4]
            symbol = symbol.lower().capitalize()
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])

        current = (index + 1) * lnsnp

        images.append(Atoms(symbols=symbols, positions=positions))

    if rvrs:
        images.reverse()

    return images[rtnndx]


def write_xyz(fileobj, images, comment=''):
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)
    for atoms in images:
        fileobj.write('%d\n%s\n' % (natoms, comment))
        for s, (x, y, z) in zip(symbols, atoms.get_positions()):
            fileobj.write('%-2s %22.15f %22.15f %22.15f\n' % (s, x, y, z))

def read_xyz_quicker(fileobj, index=-1):

    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    if not isinstance(index, int) and not isinstance(index, slice):
        print 'Index argument is neither slice nor integer! Cannot go on!'
        exit (1)

    fileobj.seek(0)

    natoms = fileobj.readline()
    bytelen = len(natoms)
    natoms = int(natoms)

    bytelen += len(fileobj.readline())
    bytelen += len(fileobj.readline()) * natoms

    fileobj.seek(0,2)
    fsize = fileobj.tell()

    if fsize % bytelen != 0:
        print "Your file is either broken or not suited for qxyz"
        exit ( 1 )
    else:
        lastsnap = fsize / bytelen

    if isinstance(index, int):
        if index < 0:
            tmpsnp = lastsnap + index
            trbl = range(tmpsnp, tmpsnp + 1, 1)
        else:
            trbl = range(index, index + 1, 1)

        rtnndx = -1

    elif isinstance(index, slice):

        start = index.start
        stop = index.stop
        step = index.step

        if start is None:
            start = 0
        elif start < 0:
            start = lastsnap + start

        if step is None:
            step = 1

        if stop < 0:
            stop = lastsnap + stop

        trbl = range(start, stop, step)

        rtnndx = slice(len(trbl))

    images = []

    for index in trbl:

        positions = []
        symbols = []

        fileobj.seek(index*bytelen)

        fileobj.readline()
        fileobj.readline()

        for i in range(natoms):
            line = fileobj.readline()

            symbol, x, y, z = line.split()[:4]
            symbol = symbol.lower().capitalize()
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])

        images.append(Atoms(symbols=symbols, positions=positions))

    return images[rtnndx]

