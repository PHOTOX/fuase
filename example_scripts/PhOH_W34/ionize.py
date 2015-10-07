#!/usr/bin/env python

from ase.io import read
from ase.calculators.orca import ORCA
from ase.calculators.gebf import GEBF
import numpy as np

atoms = read('cut-4.xyz')

natphe = 13

fragments = []
fragments.append(range(natphe))

for i in range(34):
    fragments.append(range(natphe+i*3,natphe+i*3+3))

print fragments

around_solute = 0

neu_raw = ''
xc = 'HF'
basis = 'HF-3c'
charges = 'Mulliken'

neu_label = 'large-cut_neutral_'
cat_label = 'large-cut_cation_'

frgchg = np.zeros(len(fragments), dtype='int')
frgmlt = np.ones(len(fragments), dtype='int')

frgchg[0] = 1
frgmlt[0] = 2
    
cutoff = 0.0

neu_label = 'large-cut_neutral_'

for basis in ['STO-3G', '6-31G*']:
    for xc in ['PBE', 'BHandHLYP'] :
        #### neutral #####
        
        calc = ORCA(xc=xc, basis=basis, charges=charges, raw=raw)
        
        gebf = GEBF(atoms, molecules=True, fragments=fragments, subsyscalc=calc, cutoff=cutoff, label=neu_label, printlevel=1, econv=1.e-4, qconv=1.e-2, around_solute=around_solute)
        
        atoms.set_calculator(gebf)
        
        neutral = atoms.get_potential_energy()
        print 'energy: ', neutral
        
        #### cation ####
        
        
        calc = ORCA(xc=xc, basis=basis, charges='Mulliken', raw=raw)
        
        gebf = GEBF(atoms, molecules=True, fragments=fragments, subsyscalc=calc, cutoff=cutoff, label=cat_label, printlevel=1, frgchg=frgchg, frgmlt=frgmlt, econv=1.e-4, qconv=1.e-2, around_solute=around_solute)
        
        atoms.set_calculator(gebf)
        
        cation = atoms.get_potential_energy()
        print 'energy: ', cation

        diff = cation - neutral
        print 'ionization energy: ', diff
