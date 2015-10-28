"""This module defines an ASE interface to NWchem

http://www.nwchem-sw.org/
"""
import os, sys

import numpy as np

from warnings import warn
from ase.atoms import Atoms
from ase.units import Hartree, Bohr
from ase.io.gamess_us import write_gamess_us
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError

all_changes = ['positions', 'numbers', 'cell', 'pbc',
               'initial_charges', 'initial_magmoms']


class KPoint:
    def __init__(self, s):
        self.s = s
        self.eps_n = []
        self.f_n = []


class GAMESS_US(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'dipole', 'magmom', 'charges', 'polarizability']
    command = 'rungms PREFIX.inp 2015 > PREFIX.out'

    default_parameters = dict(
        xc='SVWN', #Slater exchange plus VWN5 correlation
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
        pol=None,
        #charges=None,
        #FUDO| do we want to set default for charges to Mulliken? It works and is always given as default output
        charges='Mulliken',
        vec=None,
        mcscf=None,
        norbs=None,
        raw='')  # additional outside of dft block control string

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='gamess_us', atoms=None, delold=False, **kwargs):
        """Construct GAMESS-US-calculator object."""
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        self.delold = delold

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

        #CONTRL SECTION, COLLECT INPUT PARAMETER
        #scftyp: uhf, rohf, rhf
        #runtyp energy, gradient
        #maxit up to 199
        #mult
        #icharg
        #dfttyp none, or functional, svwn (lda)
        #exetyp run, check, debug

        contrl_temp = ' $contrl scftyp=%s runtyp=%s exetyp=%s $end\n'
        # $contrl maxit=%i icharg=%i mult=%i $end
        #"""

        contrl = contrl_temp %(p.wfn, p.task, p.exet)
        cnadd = ' $contrl'

        if p.charge is not None:
            cnadd += ' icharg=%i' %p.charge

        if 'mult' not in p:
            # Obtain multiplicity from magnetic momenta:
            tot_magmom = atoms.get_initial_magnetic_moments().sum()
            if tot_magmom < 0:
                mult = tot_magmom - 1  # fill minority bands
            else:
                mult = tot_magmom + 1
        else:
            mult = p.mult

        cnadd += ' mult=%i' %mult

        if p.xc is not 'HF':
            #FUDO| include different levels of Moller-Plesset (if possible)
            if p.xc == 'MP2':
                cnadd += ' mplevl=2'
            else:
                cnadd += ' dfttyp=%s' %p.xc

        cnadd += ' $end\n'

        contrl += cnadd
        cnadd = ''

#-------here we still need to insert more options for GAMESS-US calculator

        #basis set, problem ispher options goes into general contrl options

        #now we can add additional diffuse and polarization functions
        #heavy atom polarization functions: ndfunc
        #number of heavy atom f-type polarization: nffunc
        #number of light atom p-type polarization: npfunc
        #diffuse heavy atom: DIFFSP
        #diffuse light atom: DIFFS

#-------Figure out if we need ispher?!
#-------Pople-style basis sets start with either 3-21, 6-31, can in principle start with any number in a certain range
        #they should all start with an integer if they are pople's style

        #or just create a dictionary with all valid identifiers

        pople = False
        dunning = False
        semiemp = False

        ndfunc = 0
        nffunc = 0
        npfunc = 0
        diffsp = False
        diffs = False
        polar='common'
        tc = p.basis.strip()[0]

        try:
            ngauss = int(tc)
            pople = True
            basispar = ' $basis gbasis=%s ngauss=%i $end\n'

            if p.basis.find('-311') > -1:
                popn='N311'
                polar='popn311'
            elif p.basis.find('-31') > -1:
                popn='N31'
                polar='popn31'
            elif p.basis.find('-21') > -1:
                popn='N21'
                polar='common'

        except ValueError:
            pople = False
            if p.basis.find('cc-p') > -1:
                dunning = True
            if p.basis.find('AM1') > -1 or p.basis.find('PM3') > -1:
                semiemp = True

        if p.basis.find('+') > -1:
            diffsp = True

        if p.basis.find('++') > -1:
            diffs = True

        if p.basis.find('*') > -1:
            ndfunc = 1

        if p.basis.find('**') > -1:
            npfunc = 1

        if p.basis.find('(2d,2p)') > -1:
            ndfunc = 2
            npfunc = 2

        if p.basis.find('(2df,2pd)') > -1:
            ndfunc = 2
            npfunc = 2
            nffunc = 1
            print "Please note, that the (2df, 2pd) polarization functions don't seem to work like they should and I don't know why"

        if p.basis.find('(3df,3pd)') > -1:
            ndfunc = 3
            npfunc = 3
            nffunc = 1
            print 'Please note, that if you have Hydrogen in the system, currently one D function is not accounted for...'
            #FUDO| check, for some reason now a D function is missing for the hydrogen atoms

        #still need 6-311++G(3df,3dp), or however they are called

        cnapp = ''
        bsapp = ''
        basis = ''

        if pople:
            basis = ' $basis gbasis=%s ngauss=%i polar=%s $end\n' %(popn, ngauss, polar)

            bsapp = ''
            if diffs:
                bsapp += ' diffs=.true.'
            
            if diffsp:
                bsapp += ' diffsp=.true.'

            if ndfunc > 0:
                bsapp += ' ndfunc=%i' %ndfunc
                
            if npfunc > 0:
                bsapp += ' npfunc=%i' %npfunc

            if nffunc > 0:
                bsapp += ' nffunc=%i' %nffunc

        aug = False
        core = False
        wcr = False

        if dunning:
            ngauss = p.basis[-2]
            polar = 'dunning'
            
            bsapp = ''
            bstmp = ''

            if p.basis.startswith('aug'):
                aug = True
                bstmp += 'ACC'
            else:
                bstmp += 'CC'

            bstmp += ngauss

            if p.basis.find('pCV') > -1:
                core = True
                bstmp += 'C'

            if p.basis.find('pwCV') > -1:
                #FUDO| check the validity of the below statements and compare to GAMESS-US manual
                #if ngauss == 'T':
                wcr = True
                bstmp += 'WC'
                #else:
                #    print 'omega form of core augmentation only available for triple-zeta basis set'
                #    sys.exit ( 1 )

            basis = ' $basis gbasis=%s polar=%s $end\n' %(bstmp,polar)

        if semiemp:
            basis = ''
            bsapp = 'gbasis=%s' %p.basis

        if bsapp != '':
            basis += ' $basis %s $end\n' %bsapp
        
        if dunning and p.basispar is None:
            #print 'this is basispar: ', p.basispar
            cnapp = ' $contrl ispher=1 $end\n'
            contrl += cnapp

#-------Not exactly sure how to handle this here, because cannot just append anything to the basis set information

        cnapp = ''

        if p.basispar is not None:
            cnapp = p.basispar

        contrl += cnapp

        f.write(contrl)
        f.write(basis)

        if 'raw' in p:
            f.write(p.raw)

        write_gamess_us(f, atoms, p.geometry)

#-------GAMESS-US also can use ECPs in principle, but I don't know how...
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

            f.write(' $ffdata\n\nCOORDINATES\n')

            i = 0
            for q in p.bq:
                f.write('%i 0 %12.8f %12.8f %12.8f\n' %(i, q[0], q[1], q[2]))
                i += 1

            f.write('STOP\nPARAMETERS\n')

#FUDO| maybe at some point we should get the actual QuanPol parameters here, that'd be awefully nice :)
            if p.pol is not None:
                #FUDO| could do a sanity check here and see if p.pol has the same length as p.bq

                if len(p.pol) != len(p.bq):
                    print 'array with polarizabilities does not have the same length as the one with point charges'
                    sys.exit ( 1 )

                i = 0
                for l, q in enumerate(p.bq):
                    f.write('%i 0 %12.8f %12.8f 0 0 0 0\n' %(i, q[3], p.pol[l]))
                    i += 1
            else:
                i = 0
                for q in p.bq:
                    f.write('%i 0 %12.8f 0 0 0 0 0\n' %(i, q[3]))
                    i += 1

            f.write('STOP\n $end\n')

#FU| maybe it's better to put this into the raw input
        if p.charges is not None:

            if p.charges.find('ESP') > -1:
                if p.charges.find('RESP') > -1:
                    print 'RESP charges not implemented for GAMESS-US calculator'
                    sys.exit ( 1 )

                f.write(' $elpot iepot=1 where=pdc $end\n $pdc ptsel=chelpg constr=charge $end\n')

#FU| Mulliken and Lowdin charges are printed out by default
#FU| and RESP is not implemented directly into GAMESS-US, one would need an external program for that

        if p.wfn == 'MCSCF':
            det = ' $DET NCORE=%i NACT=%i NELS=%i STSYM=%s SZ=%3.1f $END\n' %(p.mcscf['ncore'], p.mcscf['nact'], p.mcscf['nels'], p.mcscf['stsym'], p.mcscf['sz'])
            det += ' $DET NSTATE=%i' %p.mcscf['nstate']
            try:
                det += ' IROOT=%i'  %p.mcscf['iroot']
            except KeyError:
                det += ' IROOT=1'

            det += ' $END\n'

#FU| check and maybe re-normalize wstate array at this point
            try:
                wstate = np.array(p.mcscf['wstate'])
            except KeyError:
                wstate = None

            if wstate is not None:

                #10 is kind of the maximum number of state weights that we can put on one line
                ncnt = len(wstate) / 10

                det += ' $DET WSTATE(1)='
                for w in range(0, ncnt):
                    det += '%3.1f,%3.1f,%3.1f,%3.1f,%3.1f,%3.1f,%3.1f,%3.1f,%3.1f,%3.1f,\n' %(wstate[w*10+0],wstate[w*10+1],wstate[w*10+2],wstate[w*10+3],wstate[w*10+4],wstate[w*10+5],wstate[w*10+6],wstate[w*10+7],wstate[w*10+8],wstate[w*10+9])

                resid = len(wstate) % 10

                for w in range(resid):
                    det += '%3.1f' %wstate[ncnt*10+w]
                    if w != resid-1:
                        det += ','
                
                det += ' $END\n'

            f.write(det)

            #try:
            #    det += ' $DET WSTATE(1)=
            #wstate

        if p.mcscf is not None and p.vec is not None:
            if p.norbs is None:
                norbs=int(p.vec[-2].split()[0])
            else:
                norbs=p.norbs

            f.write(' $GUESS GUESS=MOREAD NORB=%i $END\n' %norbs)
            f.writelines(p.vec)

        f.close()

        #FUDO| this way it won't work, because we only want to sort the lines not containing any atoms data
        #f = open(self.label + '.inp', 'r')
        #lines = f.readlines()
        #f.close

        #f = open(self.label + '.inp', 'w')
        #lines.sort()
        #for line in lines:
        #    f.write(line)
        #f.close()

    def read(self, label):
        FileIOCalculator.read(self, label)
        if not os.path.isfile(self.label + '.out'):
            raise ReadError

        self.parameters = Parameters.read(self.label + '.ase')
        self.atoms = read_gamess_us_input(self.label + '.inp')

        self.read_results()

    def read_results(self):
        if self.parameters.exet == 'check':
            self.results['energy'] = 0.00

            if self.parameters.task.find('gradient') > -1:
                self.results['force'] = np.zeros((3,3))

            return

        self.read_energy()
        if self.parameters.task.find('gradient') > -1:
            self.read_forces()
        self.niter = self.read_number_of_iterations()
        self.nelect = self.read_number_of_electrons()
        self.nvector = self.read_number_of_bands()
        self.results['magmom'] = self.read_magnetic_moment()
        self.results['dipole'] = self.read_dipole_moment()
        self.results['vectors'] = self.read_vectors()

        self.parameters.vec = self.results['vectors']

        #can we set here the updated vectors into something like parameters.output_vectors, because we want to use them as input for the next MCSCF calculation?

        if self.parameters.charges is not None:
            #print 'hello'
            self.results['charges'] = self.read_charges()
            #print self.results['charges']

        if self.parameters.task.find('ffield') > -1:
            self.results['polarizability'] = self.read_polarizability()

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
            if line.find('NUMBER OF ELECTRONS') != -1:
                nelect = float(line.split()[-1].strip())
        return nelect

    def get_number_of_iterations(self):
        return self.niter

    def read_number_of_iterations(self):
        niter = 0
        for line in open(self.label + '.out'):
            if line.find('ITERATIONS') != -1:  # count all iterations
                niter = int(line.strip().split()[-2])

                if self.parameters.wfn.upper() == 'MCSCF' or self.parameters.xc.upper() == 'MP2':
                    break

        if not niter:
            niter = None

        return niter

    def read_magnetic_moment(self):
        magmom = None
        for line in open(self.label + '.out'):
            if line.find('SPIN MULTIPLICITY') != -1:  # last one
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
        """Read Energy from GAMES-US output file."""
        text = open(self.label + '.out', 'r').read()
        lines = iter(text.split('\n'))

        # Energy:

        if self.parameters.xc.upper() == 'MP2':
            estring = '             E(MP2)='
        #elif self.parameters.wfn.upper() == 'MCSCF':
        #    estring = '             TOTAL POTENTIAL ENERGY = '
        elif self.parameters.wfn.upper() == 'MCSCF':
            estring = 'FINAL MCSCF ENERGY'
        else:
            estring = '                       TOTAL ENERGY = '

        if self.parameters.basis == 'AM1':
            estring = 'FINAL R-AM1 ENERGY IS'
        elif self.parameters.basis == 'PM3':
            estring = 'FINAL R-PM3 ENERGY IS'
        elif self.parameters.basis == 'MNDO':
            estring = 'FINAL R-MNDO ENERGY IS'
        elif self.parameters.basis == 'RM1':
            estring = 'FINAL R-RM1 ENERGY IS'

#FU| check here, we should also have a parameter just for the averaged energy?!
        if self.parameters.wfn.upper() == 'MCSCF':
            try:
                avgrad = self.parameters.mcscf['avgrad']
                energy = []
                if avgrad:
                    estring = ' ENERGY='
            except KeyError:
                avgrad = False
        else:
            avgrad = False

        for line in lines:
            if self.parameters.basis == 'PM3' or self.parameters.basis == 'AM1' or self.parameters.basis == 'MNDO' or self.parameters.basis == 'RM1' or self.parameters.wfn.upper() == 'MCSCF':

                if avgrad:
                    if line.startswith(estring):
                        energy.append(float(line.split()[-5]))
                else:
                    if line.find(estring) >= 0:
                        energy = float(line.split()[-4])
                        break

            elif line.find(estring) >= 0:
                energy = float(line.split()[-1])
                break

        if avgrad:
            #wstate = self.parameters.mcscf['wstate']
            wstate = self.parameters.mcscf['wstate'].copy()
            nnzr = np.where(wstate > 1.e-6)[0]

            wstate[nnzr] /= wstate[nnzr].sum()

            ener = np.array(energy)
            for i, n in enumerate(nnzr):
                ener[i] *= wstate[n]

            energy = np.sum(ener)

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
        """Read Forces from GAMES-US output file."""
        file = open(self.label + '.out', 'r')
        lines = file.readlines()
        file.close()

        if self.parameters.wfn.upper() == 'MCSCF':
            try:
                avgrad = self.parameters.mcscf['avgrad']
            except KeyError:
                avgrad = False

            if avgrad:
                gradients = []
                for i, line in enumerate(lines):
                    if line.find(' STATE-SPECIFIC GRADIENT OF STATE') >= 0:
                        for j in range(i + 4, i + 4 + len(self.atoms)):
                            word = lines[j].split()
                            gradients.append([float(word[k]) for k in range(2, 5)])
                            
                grad = np.array(gradients)
                #np.savetxt('grad.txt', grad)

                natom = len(self.atoms)
                nstate = len(grad) / natom

                grad = np.reshape(grad, (nstate, natom, 3))

                #wstate = self.parameters.mcscf['wstate']
                wstate = self.parameters.mcscf['wstate'].copy()
                nnzr = np.where(wstate > 1.e-6)[0]

                wstate[nnzr] /= wstate[nnzr].sum()

#                pref = self.label
#                f = open('%s.wstate' %pref, 'a')
#                for i in range(len(wstate)):
#                    f.write(' %20.8f' %wstate[i])
#
#                f.write('\n')
#                f.close()
#
#                print 'WSTATE: ', wstate
#
                for i, n in enumerate(nnzr):
                    grad[i,:,:] *= wstate[n]

                gradients = grad.sum(axis=0)
                #np.savetxt('avgrad.txt', gradients)

                self.results['forces'] = -np.array(gradients) * Hartree / Bohr
                return

        #it is the same for HF and MP2
        for i, line in enumerate(lines):
            if line.find('                         GRADIENT OF THE ENERGY') >= 0:
                gradients = []
                for j in range(i + 4, i + 4 + len(self.atoms)):
                    word = lines[j].split()
                    gradients.append([float(word[k]) for k in range(2, 5)])

        self.results['forces'] = -np.array(gradients) * Hartree / Bohr

    def read_polarizability(self):
        """Read polarizability from GAMES-US output file."""

        file = open(self.label + '.out', 'r')
        lines = file.readlines()
        file.close()

        polstr = '  ALPHA #        X                  Y                  Z   (A.U.)'
        skip = 2

        #FUDO| should we take either dipole or eneryg-based results?
        #it is the same for HF and MP2
        for i, line in enumerate(lines):
            if line.find(polstr) >= 0:
                polarizability = []

                for j in range(i + skip, i + skip + 3):
                    word = lines[j].split()
                    #print 'WORD: ', word
                    polarizability.append([float(word[k]) for k in range(2, 5)])

                break

        return np.array(polarizability) * Bohr**3

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

        if self.parameters.wfn.upper() == 'MCSCF':
            charges = np.zeros((len(self.atoms), 4))
        elif self.parameters.charges.find('Mulliken') > -1:
            charges = self.read_mulliken_charges()
        elif self.parameters.charges.find('MOPAC') > -1:
            charges = self.read_mopac_charges()
        elif self.parameters.charges.find('Lowdin') > -1:
            charges = self.read_lowdin_charges()
        elif self.parameters.charges.find('ESP') > -1:
            charges = self.read_esp_charges()

        self.results['charges'] = np.array(charges)[:,-1]
        return np.array(charges)[:,-1]

    def read_mulliken_charges(self):

        file = open(self.label + '.out', 'r')
        lines = file.readlines()
        file.close()

        core = self.atoms.get_atomic_numbers()

        for i, line in enumerate(lines):
            if line.find('          TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS') >= 0:

                col = 3 #this is the column of the charge magnitude, the positions are to be taken from the atoms positions

                charges = []

                for anum, j in enumerate(range(i + 2, i + 2 + len(self.atoms))):
                    word = lines[j].split()
                    arr = self.atoms.positions[anum,:].tolist() #.append(float(word[col]))
                    arr.append(float(word[col]))
                    charges.append(arr)

        return charges

    def read_mopac_charges(self):

        file = open(self.label + '.out', 'r')
        lines = file.readlines()
        file.close()

        core = self.atoms.get_atomic_numbers()

        for i, line in enumerate(lines):
            if line.find('MOPAC CHARGES') >= 0:

                col = 2 #this is the column of the charge magnitude, the positions are to be taken from the atoms positions

                charges = []

                for anum, j in enumerate(range(i + 3, i + 3 + len(self.atoms))):
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
            if line.find('          TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS') >= 0:

                col = 5 #this is the column of the charge magnitude, the positions are to be taken from the atoms positions

                charges = []

                for anum, j in enumerate(range(i + 2, i + 2 + len(self.atoms))):
                    word = lines[j].split()
                    arr = self.atoms.positions[anum,:].tolist() #.append(float(word[col]))
                    arr.append(float(word[col]))
                    charges.append(arr)

        return charges

    def read_esp_charges(self):

        file = open(self.label + '.out', 'r')
        lines = file.readlines()
        file.close()

        core = self.atoms.get_atomic_numbers()

        for i, line in enumerate(lines):
            if line.find(' NET CHARGES:') >= 0:

                col = 1 #this is the column of the charge magnitude, the positions are to be taken from the atoms positions

                charges = []

                for anum, j in enumerate(range(i + 4, i + 4 + len(self.atoms))):
                    word = lines[j].split()
                    arr = self.atoms.positions[anum,:].tolist() #.append(float(word[col]))
                    arr.append(float(word[col]))
                    charges.append(arr)

        return charges

    def read_vectors(self):

        file = open(self.label + '.dat', 'r')

        lines = file.readlines()

        start = -1
        end = -1

        for i, line in enumerate(lines):
            if line.find('$VEC') >= 0:
                start = i
                for j, nline in enumerate(lines[i+1:]):
                    if nline.find('$END') >= 0:
                        end = j+i+2
                        break

        return lines[start:end]

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):

        if os.path.isfile(self.label + '.dat') and self.delold:
            os.remove('%s.dat' %self.label)
        if os.path.isfile(self.label + '.trj') and self.delold:
            os.remove('%s.trj' %self.label)
        if os.path.isfile(self.label + '.rst') and self.delold:
            os.remove('%s.rst' %self.label)
        if os.path.isfile(self.label + '.efp') and self.delold:
            os.remove('%s.efp' %self.label)

        if 'forces' in properties:
            self.parameters.task = 'gradient'

        if 'polarizability' in properties:
            self.parameters.task = 'ffield'
            #FUDO| Do we want to set some default parameters like:, I'll just do it for my sake right now
            self.parameters.raw += ' $scf conv=1.0d-7 fdiff=.false. $end\n $contrl icut=20 itol=30 $end\n'

        FileIOCalculator.calculate(self, atoms, properties, system_changes)

