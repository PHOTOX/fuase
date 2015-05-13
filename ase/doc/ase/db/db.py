import ase.db
c = ase.db.connect('abc.db')

from ase import Atoms
from ase.calculators.emt import EMT
h2 = Atoms('H2', [(0, 0, 0), (0, 0, 0.7)])
h2.calc = EMT()
h2.get_forces()

c.write(h2, ['molecule'], relaxed=False)

from ase.optimize import BFGS
BFGS(h2).run(fmax=0.01)
c.write(h2, ['molecule'], relaxed=True, data={'abc': [1, 2, 3]})

for d in c.select('molecule'):
    print d.forces[0, 2], d.relaxed

h = Atoms('H')
h.calc = EMT()
h.get_potential_energy()
c.write(h)

d = c.get(relaxed=1, calculator='emt')
for k, v in d.items():
    print '%-25s: %s' % (k, v)

print d.data.abc

e2 = d.energy
e1 = c.get(H=1).energy
ae = 2 * e1 - e2
print ae

id = c.get(relaxed=1).id
c.update(id, atomization_energy=ae)

del c[c.get(relaxed=0).id]
