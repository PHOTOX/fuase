"""This module defines an ASE interface to NWchem

http://www.nwchem-sw.org/
"""
import os

import numpy as np

from warnings import warn
from ase.atoms import Atoms
from ase.units import Hartree, Bohr
from ase.io.nwchem import write_nwchem
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError

all_changes = ['positions', 'numbers', 'cell', 'pbc',
               'initial_charges', 'initial_magmoms']


class KPoint:
    def __init__(self, s):
        self.s = s
        self.eps_n = []
        self.f_n = []


class NWChem(FileIOCalculator):
#FU| added charges to possible properties calculation
#FU| need more variety on the charge calculation
    implemented_properties = ['energy', 'forces', 'dipole', 'magmom', 'charges']
    command = 'nwchem PREFIX.nw > PREFIX.out'
    #for testing purposes
    #command = 'mpirun -np 4 nwchem PREFIX.nw > PREFIX.out'

    default_parameters = dict(
        xc='LDA',
        smearing=None,
        charge=None,
        task='energy',
        #task='gradient',
        # Warning: nwchem centers atoms by default
        # see ase-developers/2012-March/001356.html
        geometry='nocenter noautosym',
        convergence={'energy': None,
                     'density': None,
                     'gradient': None,
                     'lshift': None,
                     # set lshift to 0.0 for nolevelshifting
                     },
        basis='3-21G',
        basispar=None,
        ecp=None,
        so=None,
        spinorbit=False,
        odft=False,
        bq=None,
        #FU| no default mulliken charges for DFT, so we can't do anything here by default
        charges=None,
        wfn='RHF',
        raw='')  # additional outside of dft block control string

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='nwchem', atoms=None, **kwargs):
        """Construct NWchem-calculator object."""
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
        f = open(self.label + '.nw', 'w')
        if p.charge is not None:
            f.write('charge %s\n' % p.charge)

        #check here, we could give here only the list of atoms that are not bq's 
        #and write an additional file containing the bq's
        write_nwchem(f, atoms, p.geometry)

        f.write('start\n')

        if p.basispar is not None:
            basispar = 'basis ' + p.basispar
        else:
            if p.basis.find('cc-p') > -1:
                basispar = 'basis "ao basis" spherical '
            else:
                basispar = 'basis'

        def format_basis_set(string, tag=basispar):
            formatted = tag + '\n'
            lines = string.split('\n')
            if len(lines) > 1:
                formatted += string
            else:
                formatted += '  * library ' + string
            return formatted + '\nend\n'

        basis = format_basis_set(p.basis)
        if p.ecp is not None:
            basis += format_basis_set(p.ecp, 'ecp')
        if p.so is not None:
            basis += format_basis_set(p.so, 'so')
        f.write(basis)

        if p.xc == 'HF':
            task = 'scf'
            f.write('\nscf\n')
            if 'mult' in p:
                if p.mult > 1:
                    if p.wfn != 'ROHF' and p.wfn != 'UHF':
                        print 'now setting your wavefunction to unrestricted as you did not specify anything else'
                        self.parameters.wfn = 'UHF'
                        p.wfn = 'UHF'

                    f.write('\nnopen %i\n' %(p.mult-1))
                
                f.write('\n%s\nend\n' %p.wfn)
            else:
                f.write('\nend\n')

        elif p.xc == 'MP2':
            task = 'mp2'
            #FUDO| the whole thing depends on the basis set, i think (Dunning depends, Pople I don't know)
            #FUDO| so we need to check that and also that it does not get set if parameters.corecol is True
            f.write('\nmp2\nfreeze atomic\nend\n')
        else:
            if p.spinorbit:
                task = 'sodft'
            else:
                task = 'dft'
            xc = {'LDA': 'slater pw91lda',
                  'PBE': 'xpbe96 cpbe96',
                  'revPBE': 'revpbe cpbe96',
                  'RPBE': 'rpbe cpbe96'}.get(p.xc, p.xc)
            f.write('\n' + task + '\n')
            f.write('  xc ' + xc + '\n')
            for key in p.convergence:
                if p.convergence[key] is not None:
                    if key == 'lshift':
                        if p.convergence[key] <= 0.0:
                            f.write('  convergence nolevelshifting\n')
                        else:
                            f.write('  convergence %s %s\n' %
                                    (key, p.convergence[key] / Hartree))
                    else:
                        f.write('  convergence %s %s\n' %
                                (key, p.convergence[key]))
            if p.smearing is not None:
                assert p.smearing[0].lower() == 'gaussian', p.smearing
                f.write('  smear %s\n' % (p.smearing[1] / Hartree))
            if 'mult' not in p:
                # Obtain multiplicity from magnetic momenta:
                #FUDO| how to handle if we have an open-shell system? Does not seem to work here...
                tot_magmom = atoms.get_initial_magnetic_moments().sum()
                if tot_magmom < 0:
                    mult = tot_magmom - 1  # fill minority bands
                else:
                    mult = tot_magmom + 1
            else:
                mult = p.mult
            if mult != int(mult):
                raise RuntimeError('Noninteger multiplicity not possible. ' +
                                   'Check initial magnetic moments.')
            if 'mult' in p:
                if p.mult > 1:
                    if p.wfn != 'ROHF' and p.wfn != 'UHF':
                        print 'now setting your wavefunction to unrestricted as you did not specify anything else'
                        self.parameters.wfn = 'UHF'
                        p.wfn = 'UHF'

                    f.write('  mult %i\n' %(mult))
                
                    if p.wfn == 'ROHF':
                        f.write('  rodft\n  cgmin\n')

            if p.odft:
                f.write('  odft\n')  # open shell aka spin polarized dft
            for key in sorted(p.keys()):
                if key in ['charge', 'geometry', 'basis', 'basispar', 'ecp',
                           'so', 'xc', 'spinorbit', 'convergence', 'smearing',
#FU| needed to add a key exception for bq, and charges...
                           'raw', 'mult', 'task', 'odft', 'bq', 'charges', 'wfn']:
#FU|
                    continue
                f.write(u"  {0} {1}\n".format(key, p[key]))
            f.write('end\n')

#FU| with NWChem 6.3 both variants are equivalent, something was weird with a previous version
        if p.bq is not None:
            #f.write('\nbq\n load %s format 1 2 3 4\nend\n' %('.bq'))
#FUDO| this label splitting kinda sucks, but it seems to work right now if the whole thing
#FUDO| is run from (a) directory(ies) above
            f.write('\nbq\n load %s format 1 2 3 4\nend\n' %(self.label.split('/')[-1] + '.bq'))
            bqf = open(self.label + '.bq', 'w')
            for q in p.bq:
                bqf.write('%21.10f %21.10f %21.10f %21.10f\n' %(q[0], q[1], q[2], q[3]))

            bqf.close()

#FU| careful in earlier NWChems this might be the wrong order:
            #f.write('\nbq\n')
            #for q in p.bq:
            #    f.write('%21.10f %21.10f %21.10f %21.10f\n' %(q[0], q[1], q[2], q[3]))
            #f.write('\nend\n')

#FU| maybe it's better to put this into the raw input
        if p.charges is not None:

            if p.charges.find('ESP') > -1:
                f.write('\nesp\n')

                if p.charges.find('RESP') > -1:
                    f.write('restrain\n')

                f.write('end\n')
            if p.charges.find('Mulliken') > -1:
                f.write('\nProperty\nMulliken\nend\n')
#FU|

        if p.raw:
            f.write(p.raw + '\n')
        f.write('\ntask ' + task + ' ' + p.task + '\n')

        if p.charges is not None:
            if p.charges.find('ESP') > -1:
                f.write('\ntask esp\n')

            if p.charges.find('Mulliken') > -1:
                f.write('\ntask ' + task + ' property\n')

        f.close()

    def read(self, label):
        FileIOCalculator.read(self, label)
        if not os.path.isfile(self.label + '.out'):
            raise ReadError

        f = open(self.label + '.nw')
        for line in f:
            if line.startswith('geometry'):
                break
        symbols = []
        positions = []
        for line in f:
            if line.startswith('end'):
                break
            words = line.split()
            symbols.append(words[0])
            positions.append([float(word) for word in words[1:]])

        self.parameters = Parameters.read(self.label + '.ase')
        self.atoms = Atoms(symbols, positions,
                           magmoms=self.parameters.pop('initial_magmoms'))
        self.read_results()

    def read_results(self):
        self.read_energy()
        if self.parameters.task.find('gradient') > -1:
            self.read_forces()
        self.niter = self.read_number_of_iterations()
        self.nelect = self.read_number_of_electrons()
        self.nvector = self.read_number_of_bands()
        self.results['magmom'] = self.read_magnetic_moment()
        #the following will not work for regular, _unresctricted_ calculations
        #self.results['dipole'] = self.read_dipole_moment()
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
        for line in open(self.label + '.out'):
            if line.find('Vector ') != -1:  # count all printed vectors
                nvector += 1
        if not nvector:
            nvector = None
        return nvector

    def get_number_of_electrons(self):
        return self.nelect

    def read_number_of_electrons(self):
        nelect = None
        for line in open(self.label + '.out'):  # find last one
            if line.find('of electrons') != -1:
                nelect = float(line.split(':')[1].strip())
        return nelect

    def get_number_of_iterations(self):
        return self.niter

    def read_number_of_iterations(self):
        niter = 0
        for line in open(self.label + '.out'):
            if line.find('d= ') != -1:  # count all iterations
                niter += 1
        if not niter:
            niter = None
        return niter

    def read_magnetic_moment(self):
        magmom = None
        for line in open(self.label + '.out'):
            if line.find('Spin multiplicity') != -1:  # last one
                magmom = float(line.split(':')[-1].strip())
                if magmom < 0:
                    magmom += 1
                else:
                    magmom -= 1
        return magmom

    def read_dipole_moment(self):
        dipolemoment = []
        for line in open(self.label + '.out'):
            for component in [
                '1   1 0 0',
                '1   0 1 0',
                '1   0 0 1'
                ]:
                if line.find(component) != -1:
                    value = float(line.split(component)[1].split()[0])
                    value = value * Bohr
                    dipolemoment.append(value)
        if len(dipolemoment) == 0:
            assert len(self.atoms) == 1
            dipolemoment = [0.0, 0.0, 0.0]
        return np.array(dipolemoment)

    def read_energy(self):
        """Read Energy from nwchem output file."""
        text = open(self.label + '.out', 'r').read()
        lines = iter(text.split('\n'))

        # Energy:
        estring = 'Total '
        if self.parameters.xc == 'HF':
            estring += 'SCF'
        elif self.parameters.xc == 'MP2':
            estring += 'MP2'
        else:
            estring += 'DFT'
        estring += ' energy'
        for line in lines:
            if line.find(estring) >= 0:
                energy = float(line.split()[-1])
                break
        self.results['energy'] = energy * Hartree

        # Eigenstates
        spin = -1
        kpts = []
        for line in lines:
            if line.find('Molecular Orbital Analysis') >= 0:
                last_eps = -99999.0
                spin += 1
                kpts.append(KPoint(spin))
            if spin >= 0:
                if line.find('Vector') >= 0:
                    line = line.lower().replace('d', 'e')
                    line = line.replace('=', ' ')
                    word = line.split()
                    this_occ = float(word[3])
                    this_eps = float(word[5])
                    kpts[spin].f_n.append(this_occ)
                    kpts[spin].eps_n.append(this_eps)
                    if this_occ < 0.1 and this_eps < last_eps:
                        warn('HOMO above LUMO - if this is not an exicted ' +
                            'state - this might be introduced by levelshift.',
                                    RuntimeWarning)
                    last_eps = this_eps
        self.kpts = kpts

    def read_forces(self):
        """Read Forces from nwchem output file."""
        file = open(self.label + '.out', 'r')
        lines = file.readlines()
        file.close()

        for i, line in enumerate(lines):
            if line.find('ENERGY GRADIENTS') >= 0:
                gradients = []
                for j in range(i + 4, i + 4 + len(self.atoms)):
                    word = lines[j].split()
                    gradients.append([float(word[k]) for k in range(5, 8)])

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

        if self.parameters.charges.find('ESP') > -1:
            if self.parameters.charges.find('RESP') > -1:
                charges = self.read_esp_charges(kind='RESP')
            else:
                charges = self.read_esp_charges()

        if self.parameters.charges.find('Mulliken') > -1:
            charges = self.read_mulliken_charges()

        if self.parameters.charges.find('Mulldef') > -1:
            charges = self.read_mulldef_charges()

        self.results['charges'] = np.array(charges)[:,-1]
        return np.array(charges)[:,-1]

    def read_esp_charges(self, kind='ESP'):
        #there is an easier way, because NWChem produces the output file self.label + '.q.' which contains all the necessary information

        file = open(self.label + '.q', 'r')

        lines = file.readlines()
        file.close()

        #FUDO| here there is too much silly wisdom, RESP charges will of course only be generated if we ask for them, nevertheless we shouldn't do it like this here:

        cols = range(1,4) #these are the coordinates of the charges
        if kind == 'RESP':
            cols.append(5)
        else:
            cols.append(4)

        charges = []

        for j in range(1, 1 + len(self.atoms)):
            word = lines[j].split()
            charges.append([float(word[k]) for k in cols])

        #file = open(self.label + '.out', 'r')
        #lines = file.readlines()
        #file.close()

        ##FUDO| here there is too much silly wisdom, RESP charges will of course only be generated if we ask for them, nevertheless we shouldn't do it like this here:
        #for i, line in enumerate(lines):
        #    if line.find('ESP') >= 0:

        #        cols = range(2,5) #these are the coordinates of the charges
        #        if kind == 'RESP':
        #            cols.append(6)
        #        else:
        #            cols.append(5)

        #        charges = []

        #        for j in range(i + 3, i + 3 + len(self.atoms)):
        #            word = lines[j].split()
        #            charges.append([float(word[k]) for k in cols])

        return charges

    def read_mulldef_charges(self):
        #there is an easier way, because NWChem produces the output file self.label + '.q.' which contains all the necessary information

        file = open(self.label + '.out', 'r')
        lines = file.readlines()
        file.close()

        core = self.atoms.get_atomic_numbers()

        for i, line in enumerate(lines):
            if line.find('  Mulliken analysis of the total density') >= 0:
                #print 'hello'

                col = 3 #this is the column of the charge magnitude, the positions are to be taken from the atoms positions

                charges = []

                for anum, j in enumerate(range(i + 5, i + 5 + len(self.atoms))):
                    word = lines[j].split()
                    arr = self.atoms.positions[anum,:].tolist() #.append(float(word[col]))
                    #print core[anum], float(word[col])
                    arr.append(float(core[anum]) - float(word[col]))
                    #charges.append(self.atoms.positions[anum,:].tolist().append(float(word[col])))
                    charges.append(arr)
                    #print charges

        return charges

    def read_mulliken_charges(self):
        #there is an easier way, because NWChem produces the output file self.label + '.q.' which contains all the necessary information

        file = open(self.label + '.out', 'r')
        lines = file.readlines()
        file.close()

        core = self.atoms.get_atomic_numbers()

        for i, line in enumerate(lines):
            if line.find('Total      gross population on atoms') >= 0:
                #print 'hello'

                col = 3 #this is the column of the charge magnitude, the positions are to be taken from the atoms positions

                charges = []

                for anum, j in enumerate(range(i + 2, i + 2 + len(self.atoms))):
                    word = lines[j].split()
                    arr = self.atoms.positions[anum,:].tolist() #.append(float(word[col]))
                    #print core[anum], float(word[col])
                    arr.append(float(core[anum]) - float(word[col]))
                    #charges.append(self.atoms.positions[anum,:].tolist().append(float(word[col])))
                    charges.append(arr)
                    #print charges

        return charges

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):

        if 'forces' in properties:
            self.parameters.task = 'gradient'

        if 'charges' in properties and self.parameters.charges is None:
            self.parameters.charges = 'Mulliken'

        FileIOCalculator.calculate(self, atoms, properties, system_changes)

