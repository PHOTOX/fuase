from ase.atoms import Atoms
from ase.parallel import paropen

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

