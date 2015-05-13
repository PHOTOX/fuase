"""This module defines an ASE interface to NWchem

http://www.nwchem-sw.org/
"""
import os, sys

import numpy as np

from warnings import warn
from ase.atoms import Atoms
from ase.units import Hartree, Bohr
from ase.io.orca import write_orca
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError

all_changes = ['positions', 'numbers', 'cell', 'pbc',
               'initial_charges', 'initial_magmoms']


class KPoint:
    def __init__(self, s):
        self.s = s
        self.eps_n = []
        self.f_n = []


class ORCA(FileIOCalculator):
#FU| added charges to possible properties calculation
#FU| need more variety on the charge calculation
    implemented_properties = ['energy', 'forces', 'dipole', 'magmom', 'charges']
    command = 'orca PREFIX.inp > PREFIX.out'
    #for testing purposes
    #command = 'mpirun -np 4 nwchem PREFIX.nw > PREFIX.out'

    default_parameters = dict(
        xc='DFT', #referes to VWN5 and Slater exchange, for some cocked-up reason keyword LDA turns on the RI approximation (v3.0.2)
        wfn='RHF',
        exet='run',
        smearing=None,
        charge=None,
        task='energy',
        geometry=None,
        convergence={'energy': None,
                     'density': None,
                     'gradient': None,
                     'lshift': None,
                     # set lshift to 0.0 for nolevelshifting
                     },
        basis='3-21G',
        basispar=None,
        ecp=None,
        bq=None,
        charges='Mulliken',
        raw='')  # additional outside of dft block control string

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='orca', atoms=None, **kwargs):
        """Construct ORCA-calculator object."""
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def check_state(self, atoms):
        system_changes = FileIOCalculator.check_state(self, atoms)
        # Ignore unit cell and boundary conditions:
        if 'cell' in system_changes:
            system_changes.remove('cell')
        if 'pbc' in system_changes:
            system_changes.remove('pbc')
        return system_changes

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters
        p.initial_magmoms = atoms.get_initial_magnetic_moments().tolist()
        p.write(self.label + '.ase')
        del p['initial_magmoms']
        f = open(self.label + '.inp', 'w')

        #FUDO| we'll construct a simple input line by default and everything special should be supplied under 'raw' input
        contrl_temp = '! %s %s'

        contrl = contrl_temp %(p.wfn, p.basis)

        cnadd = ''

        if p.xc is not 'HF': 
            cnadd += ' %s' %p.xc

        contrl += cnadd
        contrl += '\n'

        if p.task == 'gradient':
            contrl += '!EnGrad\n'

        if 'mult' not in p:
            # Obtain multiplicity from magnetic momenta:
            tot_magmom = atoms.get_initial_magnetic_moments().sum()
            if tot_magmom < 0:
                mult = tot_magmom - 1  # fill minority bands
            else:
                mult = tot_magmom + 1
        else:
            mult = p.mult

        f.write(contrl)

        if 'raw' in p:
            f.write(p.raw)

#FUDO| assuming a neutral system
        if p.charge is not None:
            charge = p.charge
        else:
            charge = 0

        write_orca(f, atoms, charge, mult)

#-------orca also can use ECPs in principle, but I don't know how...
#-------old NWChem input generation---------#

        #if p.ecp is not None:
        #    basis += format_basis_set(p.ecp, 'ecp')
        #if p.so is not None:
        #    basis += format_basis_set(p.so, 'so')
        #f.write(basis)

###        for key in p.convergence:
###            if p.convergence[key] is not None:
###                if key == 'lshift':
###                    if p.convergence[key] <= 0.0:
###                        f.write('  convergence nolevelshifting\n')
###                    else:
###                        f.write('  convergence %s %s\n' %
###                                (key, p.convergence[key] / Hartree))
###                else:
###                    f.write('  convergence %s %s\n' %
###                            (key, p.convergence[key]))
###        if p.smearing is not None:
###            assert p.smearing[0].lower() == 'gaussian', p.smearing
###            f.write('  smear %s\n' % (p.smearing[1] / Hartree))
###
###        if p.odft:
###            f.write('  odft\n')  # open shell aka spin polarized dft
###        for key in sorted(p.keys()):
###            if key in ['charge', 'geometry', 'basis', 'basispar', 'ecp',
###                       'so', 'xc', 'spinorbit', 'convergence', 'smearing',
####FU| needed to add a key exception for bq, and charges...
###                           'raw', 'mult', 'task', 'odft', 'bq', 'charges']:
####FU|
###                    continue
###                f.write(u"  {0} {1}\n".format(key, p[key]))
###            f.write('end\n')

#FU| for now we will implement external charges within the FFDATA/QUANPOL module, just without any parameters but charges defined

        if p.bq is not None:

            f.write('%%pointcharges "%s"' %(self.label + '.bq'))

            bqf = open(self.label + '.bq', 'w')

            bqf.write('%i\n' %len(p.bq))

            for q in p.bq:
                bqf.write('%12.8f %12.8f %12.8f %12.8f\n' %(q[3], q[0], q[1], q[2]))

#FU| careful in earlier NWChems this might be the wrong order:

#FU| maybe it's better to put this into the raw input
        if p.charges is not None:

            if p.charges.find('ESP') > -1:
                print 'ESP or RESP charges not implemented for orca calculator'
                sys.exit ( 1 )

#FU| Mulliken and Lowdin charges are printed out by default
#FU| and RESP is not implemented directly into orca, one would need an external program for that
#FUDO| check, (R)ESP fits might actually be implemented into orca

        f.close()

    def read(self, label):
        FileIOCalculator.read(self, label)
        if not os.path.isfile(self.label + '.out'):
            raise ReadError

        self.parameters = Parameters.read(self.label + '.ase')
        self.atoms = read_orca_input(self.label + '.inp')

        self.read_results()

    def read_results(self):
        self.read_energy()
        if self.parameters.task.find('gradient') > -1:
            self.read_forces()
        self.niter = self.read_number_of_iterations()
        self.nelect = self.read_number_of_electrons()
        self.nvector = self.read_number_of_bands()
        self.results['magmom'] = self.read_magnetic_moment()
        self.results['dipole'] = self.read_dipole_moment()
        if self.parameters.charges is not None:
            #print 'hello'
            self.results['charges'] = self.read_charges()
            #print self.results['charges']

    def get_ibz_k_points(self):
        return np.array([0., 0., 0.])

    def get_number_of_bands(self):
        return self.nvector

    def read_number_of_bands(self):
        nvector = 0
#FUDO| still needs to be done
        ###for line in open(self.label + '.out'):
        ###    if line.find('Vector ') != -1:  # count all printed vectors
        ###        nvector += 1
        ###if not nvector:
        ###    nvector = None
        return nvector

    def get_number_of_electrons(self):
        return self.nelect

    def read_number_of_electrons(self):
        nelect = None
        for line in open(self.label + '.out'):  # find last one
            if line.find(' Number of Electrons') != -1:
                nelect = float(line.split()[-1].strip())
        return nelect

    def get_number_of_iterations(self):
        return self.niter

    def read_number_of_iterations(self):

        niter = 0
        for line in open(self.label + '.out'):
            if line.find('SCF CONVERGED AFTER') != -1:  # count all iterations
                niter = int(line.strip().split()[-3])

        if not niter:
            niter = None

        return niter

    def read_magnetic_moment(self):
        magmom = None
        for line in open(self.label + '.out'):
            if line.find('Multiplicity           Mult') != -1:  # last one
                magmom = float(line.split()[-1].strip())
                if magmom < 0:
                    magmom += 1
                else:
                    magmom -= 1
        return magmom

    def read_dipole_moment(self):
#FUDO| still needs to be done
        dipolemoment = []
        ###for line in open(self.label + '.out'):
        ###    for component in [
        ###        '1   1 0 0',
        ###        '1   0 1 0',
        ###        '1   0 0 1'
        ###        ]:
        ###        if line.find(component) != -1:
        ###            value = float(line.split(component)[1].split()[0])
        ###            value = value * Bohr
        ###            dipolemoment.append(value)
        ###if len(dipolemoment) == 0:
        ###    assert len(self.atoms) == 1
        ###    dipolemoment = [0.0, 0.0, 0.0]
        return np.array(dipolemoment)

    def read_energy(self):
        """Read Energy from nwchem output file."""
        text = open(self.label + '.out', 'r').read()
        lines = iter(text.split('\n'))

        #hartree fock column

        estring = 'Total Energy       :'
        col = -4

        # Energy:
        if self.parameters.xc == 'MP2':
            estring = 'MP2 TOTAL ENERGY:'
            col = -2

        for line in lines:
            if line.find(estring) >= 0:
                energy = float(line.split()[col])
                break

        self.results['energy'] = energy * Hartree

#FUDO| this is not of immediate concern, but we can figure it out later
        self.kpts = []

        # Eigenstates
###        spin = -1
###        kpts = []
###        for line in lines:
###            if line.find('Molecular Orbital Analysis') >= 0:
###                last_eps = -99999.0
###                spin += 1
###                kpts.append(KPoint(spin))
###            if spin >= 0:
###                if line.find('Vector') >= 0:
###                    line = line.lower().replace('d', 'e')
###                    line = line.replace('=', ' ')
###                    word = line.split()
###                    this_occ = float(word[3])
###                    this_eps = float(word[5])
###                    kpts[spin].f_n.append(this_occ)
###                    kpts[spin].eps_n.append(this_eps)
###                    if this_occ < 0.1 and this_eps < last_eps:
###                        warn('HOMO above LUMO - if this is not an exicted ' +
###                            'state - this might be introduced by levelshift.',
###                                    RuntimeWarning)
###                    last_eps = this_eps
###        self.kpts = kpts

    def read_forces(self):
        """Read Forces from nwchem output file."""
        file = open(self.label + '.out', 'r')
        lines = file.readlines()
        file.close()

        if self.parameters.xc == 'MP2':
            fstring = 'The final MP2 gradient'
            skip = 1
            col = 1
        else:
            fstring = 'CARTESIAN GRADIENT'
            skip = 3
            col = 3

        for i, line in enumerate(lines):
            if line.find(fstring) >= 0:
                gradients = []
                for j in range(i + skip, i + skip + len(self.atoms)):
                    word = lines[j].split()
                    gradients.append([float(word[k]) for k in range(col, col+3)])

        self.results['forces'] = -np.array(gradients) * Hartree / Bohr

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalue array."""
        return np.array(self.kpts[spin].eps_n) * Hartree

    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array."""
        return self.kpts[spin].f_n

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2."""
        return len(self.kpts)

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""
        return len(self.kpts) == 2

    def read_charges(self, atoms=None):
        """read output from charge calculation"""

        #need to check status of calculator, if it was run, otherwise re-run it
        #dirty work-around is just to run it, because that should check by itself everything

        #self.get_potential_energy()

        if self.parameters.charges.find('Mulliken') > -1:
            charges = self.read_mulliken_charges()
        elif self.parameters.charges.find('Lowdin') > -1:
            charges = self.read_lowdin_charges()

        self.results['charges'] = np.array(charges)[:,-1]
        return np.array(charges)[:,-1]

    def read_mulliken_charges(self):

        file = open(self.label + '.out', 'r')
        lines = file.readlines()
        file.close()

        core = self.atoms.get_atomic_numbers()

        for i, line in enumerate(lines):
            if line.find('MULLIKEN ATOMIC CHARGES') >= 0:

                if self.parameters.wfn.upper() == 'RHF':
                    col = -1 # this needs confirmation #3 only works for mono-letter symbols #this is the column of the charge magnitude, the positions are to be taken from the atoms positions
                else:
                    col = -2

                charges = []

                for anum, j in enumerate(range(i + 2, i + 2 + len(self.atoms))):
                    word = lines[j].split()
                    arr = self.atoms.positions[anum,:].tolist() #.append(float(word[col]))
                    arr.append(float(word[col]))
                    charges.append(arr)

        return charges

    def read_lowdin_charges(self):

        file = open(self.label + '.out', 'r')
        lines = file.readlines()
        file.close()

        core = self.atoms.get_atomic_numbers()

        for i, line in enumerate(lines):
            if line.find('LOEWDIN ATOMIC CHARGES') >= 0:

                col = 3 #this is the column of the charge magnitude, the positions are to be taken from the atoms positions

                charges = []

                for anum, j in enumerate(range(i + 2, i + 2 + len(self.atoms))):
                    word = lines[j].split()
                    arr = self.atoms.positions[anum,:].tolist() #.append(float(word[col]))
                    arr.append(float(word[col]))
                    charges.append(arr)

        return charges

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):

        if 'forces' in properties:
            self.parameters.task = 'gradient'

        FileIOCalculator.calculate(self, atoms, properties, system_changes)

