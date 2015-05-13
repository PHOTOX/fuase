from ase import Atoms
from ase.db import connect
from ase.structure import molecule
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms, FixBondLength


for name in ['y2.json', 'y2.db']:
    c = connect(name)
    print(name, c)

    id = c.reserve(abc=7)
    c.delete([d.id for d in c.select(abc=7)])
    id = c.reserve(abc=7)
    assert c[id].abc == 7
    
    a = c.get_atoms(id)
    c.write(Atoms())
    ch4 = molecule('CH4', calculator=EMT())
    ch4.constraints = [FixAtoms(indices=[1]),
                       FixBondLength(0, 2)]
    f1 = ch4.get_forces()
    print(f1)
    
    c.delete([d.id for d in c.select(C=1)])
    id = c.write(ch4)
    f2 = c.get(C=1).forces
    assert abs(f2.sum(0)).max() < 1e-14
    f3 = c.get_atoms(C=1).get_forces()
    assert abs(f1 - f3).max() < 1e-14

    try:
        c.update(id, abc={'a': 42})
    except ValueError:
        pass
    else:
        2 / 0

    c.update(id, grr='hmm')
    assert c.get(C=1).id == id
