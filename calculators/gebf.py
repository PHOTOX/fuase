"""Calculator for General Energy-based Fragmentation"""

import os
import copy
import subprocess
import numpy as np
from ase.data import vdw
from ase.calculators.calculator import Calculator, all_changes
from sys import stdout, exit
from os import system
from ase.units import Ha, Bohr
import ase.units as u
from ase.io import write
#from future import printf
from copy import deepcopy
from ase.calculators.nwchem import NWChem
from ase.calculators.gamess_us import GAMESS_US
from ase.calculators.orca import ORCA

class GEBF(Calculator):
    """Initializes the General Energy-based Fragmentation object

    Parameters:

    atoms: atoms (object)

    fragments: list of lists containing indices defining moleculesor generally fragments, can be optional and
               molecular fragments will be generated based on the covalent radii
               of the atoms in the atoms object (latter option not implemented, yet)

               if a single number is provided it is used as a fudge factor for (option not implemented, yet)
               the molecule generation
    
    gradients: boolean determining whether gradients are to be calculated (these are more or less meaningless)

    fullgrad: boolean determining whether fully correct gradients are to be
              calculated or only the ones based on molecule by molecule embedding
    """

    #incorporate this later
    #these are currently just the same as for NWChem except for the included bq
    #    default_parameters = dict(
    #        xc='LDA',
    #        smearing=None,
    #        charge=None,
    #        task='gradient',
    #        geometry='nocenter noautosym',
    #        convergence={'energy': None,
    #                     'density': None,
    #                     'gradient': None,
    #                     'lshift': None,
    #                     # set lshift to 0.0 for nolevelshifting
    #                     },
    #        basis='3-21G',
    #        basispar=None,
    #        ecp=None,
    #        so=None,
    #        spinorbit=False,
    #        odft=False,
    #        bq=None,
    #        raw='')  # additional outside of dft block control string
    #

    implemented_properties = ['energy', 'forces', 'charges']

    def __init__(self, atoms=None, molecules=False, fragments=None, gradients=False, fullgrad=False, subsyscalc=None, frgchg=None, frgmlt=None, around_solute=None, printlevel=0, cutoff=None, prisub=None, restart=None, bq_initial_guess=None, max_iter=100, ignore_bad_restart_file=False, label=os.curdir, econv=1.e-6, qconv=1.e-4, **kwargs):
        """Initializes calculator.
        The idea is to either use pure molecular fragments or not.
        If molecules is set to True fragments can still be read in, but are assumed to be molecules.
        If molecules is True and fragments is none, molecular fragments will be generated. (not implemented currently)
        If molecules is True and fragments is a float, molecular fragments will be generated with the value of float as fudge factor for molecule generation. (not implemented currently)

        subsys calc is the calculator for each specific subsystem
        prisub are optionally provided primary subsystems. If They are provided, they will not be changed/updated. (can be used for validation purposes against original GEBF literature, w/o having capped fragments, but calculation will not work)
        """

    #FUDO| what about charged systems???? do we need to provide charges of each fragment?!?!?!?!?!
    #FUDO| otherwise just set if self.charges is None\n self.charges=np.zeros(len(self.fragments), dtype='int')
    #FUDO| and further down we need to incorporate that into the subsystem information

        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            label, atoms, **kwargs)

        self.molecules = molecules
        self.plvl = printlevel
        self.distmat = None

        self.econv = econv
        self.qconv = qconv

        self.subsyscalc = subsyscalc

        self.bq_initial_guess = bq_initial_guess
        self.max_iter = max_iter

        self.set_atoms(atoms)

        if self.subsyscalc is None:
            print  "No subsys Calculator found, cannot continue"
            exit ( 1 )

        #FUDO| make it work for Gaussian calculator as well
        if not any(isinstance(self.subsyscalc, e) for e in [NWChem, GAMESS_US, ORCA]):
            print  "GEBF calculator is currently only supported with NWChem, GAMESS-US and ORCA"
            exit ( 1 )

        self.selfenercorr = True

        if isinstance(self.subsyscalc, NWChem):
            self.selfenercorr = False
        elif isinstance(self.subsyscalc, GAMESS_US):
            self.selfenercorr = True
        elif isinstance(self.subsyscalc, ORCA):
            self.selfenercorr = False

        if self.molecules and not fragments:
            print "The option to generate arbitrary molecular fragments is currently not implemented"
            exit ( 1 )

            #if fragments is None:
            #    self.fudge = 1.0

            #elif isinstance(fragments, float):
            #    self.fudge = fragments
            #    
            #self.fragments = self.generate_molecules(fudge=self.fudge)

        elif not self.molecules and not fragments:
            print "The option to generate arbitrary non-molecular fragments is currently not implemented"
            exit ( 1 )

            #self.fragments = self.generate_fragments()

        else:
            self.fragments = fragments

        self.prisub = prisub

        if self.prisub is not None:
            print 'Will use user-provided primary subsystems. These will not be changed under any circumstances (like MD)!!!'

        self.cutoff = cutoff

        if frgchg is not None:
            self.frgchg = frgchg
        else:
            self.frgchg = np.zeros(len(self.fragments), dtype='int')

        if frgmlt is not None:
            self.frgmlt = frgmlt
        else:
            self.frgmlt = np.ones(len(self.fragments), dtype='int')

        self.around_solute = around_solute
        self.coeffs, self.subsys, self.subsyschg, self.subsysmlt = self.update_subsystems(cutoff=self.cutoff, prisub=self.prisub, frgchg=self.frgchg, frgmlt=self.frgmlt, around_solute=self.around_solute)

        #self.subsyschg = self.get_subsystem_charges(frgchg=self.frgchg, subsys=self.subsys)
        #self.subsysmlt = self.get_subsystem_multiplicities(frgmlt=self.frgmlt, subsys=self.subsys)

    def get_subsystem_charges(self, frgchg=None, subsys=None):
        if frgchg is None:
            print 'cannot calculate subsystem charges without fragment charges'
            exit ( 1 )

        if subsys is None:
            print 'cannot calculate subsystem charges without subsystems'
            exit ( 1 )

        subsyschg = []

        for sub in subsys:
            tmp = 0
            for frg in sub:
                tmp += frgchg[frg]

            subsyschg.append(tmp)

        return subsyschg

    def get_subsystem_multiplicities(self, frgmlt=None, subsys=None, highspin=False):
        if frgmlt is None:
            print 'cannot calculate subsystem multiplicities without fragment multiplicities'
            exit ( 1 )

        if subsys is None:
            print 'cannot calculate subsystem charges without subsystems'
            exit ( 1 )

        subsysmlt = []

        for sub in subsys:
            nprd = 0
            for frg in sub:
                tmp = ( frgmlt[frg] - 1 ) / 2.
                nprd += tmp

            if highspin:
                subsysmlt.append(int(2*nprd+1))
            else:
                subsysmlt.append(int(2*(nprd % 2) + 1))

        return subsysmlt

    def update_subsystems(self, cutoff=3.0, prisub=None, frgchg=None, frgmlt=None, around_solute=None):

        if self.cutoff is None and self.prisub is None:
            return self.generate_subsystems(frgchg=frgchg, frgmlt=frgmlt, around_solute=around_solute)
        elif self.cutoff is None and self.prisub is not None:
            return self.generate_subsystems(prisub=self.prisub, frgchg=frgchg, frgmlt=frgmlt, around_solute=around_solute)
        elif self.cutoff is not None and self.prisub is None:
            return self.generate_subsystems(cutoff=self.cutoff, frgchg=frgchg, frgmlt=frgmlt, around_solute=around_solute)
        else:
            return self.generate_subsystems(cutoff=self.cutoff, prisub=self.prisub, frgchg=frgchg, frgmlt=frgmlt, around_solute=around_solute)

    def generate_molecules(self, fudge=1.0):
        """Generate molecules based on distance comparison to covalent radii
           Todo: add consistency check with database molecules
        """

        #we can steal this routine from CCUBE, it should be fairly easy to do

        molecules = []

        return molecules

    def generate_primary_subsystems(self, cutoff=3.0, around_solute=None):
        """Generate primary subsystems according to distance cutoff 
           Currently only molecules are considered as basic building blocks, but extension to arbitrary
           functional groups is straight forward. We then need to take care of capping atoms.
           Further on, better distance search could be done, currently only checking for the center of mass of molecules.
        """

        frgs = self.fragments
        coms = []

        nfrg = len(frgs)

        for frag in frgs:
            tmp = self.atoms[frag].get_center_of_mass()
            coms.append(tmp)

        coms = np.array(coms)

        frgdst = np.zeros((nfrg, nfrg))

        for i in range(nfrg):
            for j in range(i+1, nfrg):
                d = coms[i] - coms[j]
                #FUDO| keep to minimum image convention here?
                #FUDO| include via periodicity of self.atoms object
                frgdst[i,j] = np.linalg.norm(d)
                frgdst[j,i] = frgdst[i,j]

        #check which molecules are within range of cutoff from each fragment (i.e., construct fragment domain to later construct primitive subsystems)
        prisub = []

        eta = 0

        if around_solute is None:
            for i in range(nfrg):
                tmp = [i]
                for j in range(nfrg):
                    if i == j:
                        continue
                    if frgdst[i,j] <= cutoff:
                        tmp.append(j)
    
                tmp.sort()
                if len(tmp) > eta:
                    eta = len(tmp)
    
                prisub.append(tmp)
                #print 'test: ' ,tmp
    
            prisub.sort()
    
            doubles, prisub = self.remove_duplicates_lol(prisub)
    
            tor = []
            print prisub
            for i, s1 in enumerate(prisub):
                for j, s2 in enumerate(prisub[i+1:]):
    
                    lc = list(set(s2).intersection(s1))
    
                    if (lc == s1):
                        tor.append(s1)
                    elif (lc == s2):
                        tor.append(s2)
    
            if len(tor) > 1:
                doubles, tor = self.remove_duplicates_lol(tor)
            for rem in tor:
                print 'to remove: ', rem
                prisub.remove(rem)
        
        else:

            slute = []
            slvnt = []

            prisub = [[around_solute]]

            for i in range(nfrg):
                if frgdst[around_solute, i] <= cutoff:
                    if i != around_solute:
                        prisub[0].append(i)
                else:
                    prisub.append([i])

        if self.plvl > 0:
            print 'largest primitive subsys (in terms of \# fragments) has %i blocks in it' %eta

        return prisub

    def generate_subsystems(self, cutoff=3.0, prisub=None, frgchg=None, frgmlt=None, around_solute=None):
        """Generate subsystems according to original GEBF literature. JPCA 114, 8126 (2010)
           Currently only molecules are considered as basic building blocks, but extension to arbitrary
           functional groups is straight forward. We then need to take care of capping atoms.
           Further on, better distance search could be done, currently only checking for the center of mass of molecules.
        """

        frgs = self.fragments
        nfrg = len(frgs)

        if prisub is None:
            prisub = self.generate_primary_subsystems(cutoff=cutoff, around_solute=around_solute)

        eta = 0
        for sub in prisub:
            if len(sub) > eta:
                eta = len(sub)

        #FU| remove small subsystems that are completely contained in bigger subsystems
        tor = []
        for i, s1 in enumerate(prisub):
            for j, s2 in enumerate(prisub[i+1:]):

                lc = list(set(s2).intersection(s1))

                if (lc == s1):
                    tor.append(s1)
                elif (lc == s2):
                    tor.append(s2)

        for rem in tor:
            prisub.remove(rem)

        if self.plvl > 0:
            print 'largest primitive subsys (in terms of \# fragments) has %i blocks in it' %eta

        #FU| now compare fragments and remove obsolete ones
        nsub = len(prisub)

        if len(prisub) > 1:
            doubles, prisub = self.remove_duplicates_lol(prisub)
            if self.plvl > 0:
                print 'doubles: ', doubles

        if self.plvl > 0:
            print 'primitive subsystems: '

        if self.plvl > 0:
            for pri in prisub:
                print pri
        
        sbs = prisub
        self.curprisub = prisub
        cff = np.ones(len(sbs))

        #FU| now all primitive subsystems should be unique
        #FU| we can compare between all primitive subsystems to get the derivative subsystems
        #FU| again, compare each and every one subsystem with all others... (see JPCA 114, 8126 (2010))

        dersub = []
        dercoe = []

        wrketa = eta - 1

        while wrketa > 0:
            fndsth = False
            sth = []

            for i, s1 in enumerate(sbs):
                for s2 in sbs[i+1:]:
    
                    tmp = list(set(s1).intersection(s2))
                    tmp.sort()

                    lnlst = len(tmp)

                    if lnlst == wrketa:
                        fndsth = True
                        sth.append(tmp)
                    else:
                        continue

            newsbs = []
            newcff = []

            if not fndsth:
                wrketa -= 1
            else:
                if len(sth) > 1:
                    bla, sth = self.remove_duplicates_lol(sth)

                #FUDO| possible to mod the thing into one function (the below appears similarly in other places)?
                tmpc = []
                for j, new in enumerate(sth):
                    tmpc.append(1)
                    for k, sub in enumerate(sbs):
                        intr = list(set(new).intersection(sub))

                        #FUDO| int - float
                        if len(intr) == wrketa:
                            tmpc[j] -= cff[k]

                    if tmpc[j] != 0:
                        newsbs.append(sth[j])
                        newcff.append(tmpc[j])
                        print 'found new subsystem: ', sth[j], 'with coeff: ', tmpc[j]
                    else:
                        print 'found subsystem with 0 coeff: ', sth[j], tmpc[j]

                if len(newsbs) > 0:
                    for k, sub in enumerate(newsbs):
                        sub.sort()
                        sbs.append(sub)
                        cff = np.append(cff, newcff[k])

                wrketa -= 1

        subsys= []

        if self.plvl > 0:
            print 'all subsystems expressed in atoms:'

        #FUDO| make highspin optional
        subsysmlt = self.get_subsystem_multiplicities(frgmlt, sbs, highspin=False)
        subsyschg = self.get_subsystem_charges(frgchg, sbs)

        #FU| go from list of fragments to list of list of atoms
        for i, mols in enumerate(sbs):
            subsys.append([])
            for mol in mols:
                for j in frgs[mol]:
                    subsys[i].append(j)

            syms = self.atoms.get_chemical_symbols()
            atcp = self.atoms.copy()

            xat = [elem for k, elem in enumerate(range(len(self.atoms))) if k not in subsys[i]]

            for x in xat:
                syms[x] = 'X'

            atcp.set_chemical_symbols(syms)

            write(self.label + 'frag-%i.xyz' %i, atcp)

            if self.plvl > 0:
                print subsys[i]

        coeffs = cff
        np.savetxt(self.label + 'coefficients.dat', coeffs)

        return coeffs, subsys, subsyschg, subsysmlt

    def run(self, atoms=None, gradients=False):

        if (self.subsys is None) or (self.subsys == []):
            print "Something is wrong, no subsystems defined"
            #printf("Something is wrong, no subsystems defined\n")
            exit ( 1 )

        itr = 0

        if self.bq_initial_guess is None:
            if gradients:
                charges, energies, tmpfrc = self.one_iteration(atoms=atoms, dobqs=False, gradients=gradients, chgs=self.subsyschg, mlts=self.subsysmlt)
            else:
                charges, energies = self.one_iteration(atoms=atoms, dobqs=False, gradients=gradients, chgs=self.subsyschg, mlts=self.subsysmlt)
        else:
            charges = self.bq_initial_guess.copy()
            energies = np.zeros(len(self.subsys))

            if gradients:
                forces = np.zeros((len(self.subsys), len(atoms), 3))

        energies /= Ha

        if gradients:
            tmpfrc /= (Ha / Bohr)
            forces = np.zeros((len(atoms),3))

            for k in range(len(self.subsys)):
                forces += self.coeffs[k] * tmpfrc[k,:,:]

        if self.plvl > 0:
            stdout.write('\n')

        atoms.set_array('charges', charges)

        if self.plvl > 1:
            print charges
        
        np.savetxt(self.label + 'initial-charges.dat', charges, fmt='%21.10f')
        np.savetxt(self.label + 'initial-energies.dat', energies, fmt='%21.10f')

        if gradients:
            np.savetxt(self.label + 'initial-forces.dat', forces, fmt='%21.10f')

            for i, subsys in enumerate(self.subsys):
                np.savetxt(self.label + 'initial-forces-subsys-%i.dat' %i, tmpfrc[i,:,:], fmt='%21.10f')

        converged = False

        tmpen = (energies * self.coeffs).sum()
        totnr = [tmpen]

        if self.plvl == -1 or self.plvl > 0:
            stdout.write("Energy without surrounding point charges: %21.10f\n" %totnr[itr])

    #FU| careful, if we only have one subsystem (for whatever reason), we can just put converged to true

        if len(self.subsys) == 1:
            print 'only found one subsystem, so we are done!'
            converged = True
        
        if self.max_iter == 0:
            print 'no self-consistent iterations requested'
            converged = True

        #os.system('rm inpsubsys-*')
        while not converged:

            if self.plvl == -1 or self.plvl > 0:
                stdout.write('\n')

            if gradients:
                nchg, nenr, tmpfrc = self.one_iteration(atoms=atoms, gradients=gradients, chgs=self.subsyschg, mlts=self.subsysmlt)
            else:
                nchg, nenr = self.one_iteration(atoms=atoms, gradients=gradients, chgs=self.subsyschg, mlts=self.subsysmlt)

    #FU| convert to a.u., because ASE returns eV
    #FUDO| change back, because will be internally inconsistent with ASE, but need it now, for double-checking
            nenr /= Ha

            if gradients:
                tmpfrc /= (Ha / Bohr)
                nfrc = np.zeros((len(atoms),3))
                for k, subsys in enumerate(self.subsys):
                    np.savetxt(self.label + 'forces-%i-subsys-%i.dat' %(itr,k), tmpfrc[k,:,:], fmt='%21.10f')
                    nfrc += self.coeffs[k] * tmpfrc[k,:,:]

    #FUDO| probably need to convert forces as well, but need to check appropriate units first

    #FU| just the sum of the energies of primitive and derivative fragments
            tmpen = (self.coeffs * nenr).sum()

    #FU| self-energy correction (depends on calculator at hand) it is actually not tested right now, because NWChem does not need it

            if self.selfenercorr:
                slfnr = self.get_mm_self_energy(atoms, charges)
                if self.plvl == -1 or self.plvl > 0:
                    stdout.write('self-energy of point charges: %21.10f\nWhich appears %i times\n' %(slfnr, self.coeffs.sum()-1))
            else:
                slfnr = 0.
                if self.plvl == -1 or self.plvl > 0:
                    stdout.write('no self-energy correction for point charges applied\n')

    #FU| for some reason we will do this here, maybe we can puch something else into slfnr later
    #FUDO|
            tmpen -= (self.coeffs.sum()-1) * slfnr

            totnr.append(tmpen)

            totdel = totnr[-1] - totnr[itr]

            if self.plvl == -1 or self.plvl > 0:
                stdout.write('total energy in iteration %i: %21.10f\n' %(itr, totnr[-1]))


            dchg = nchg - charges
            denr = nenr - energies
            if gradients:
                dfrc = nfrc - forces

            if self.plvl > 1:
                print charges
                print dchg

    #perform additional checks here, whether change in all energies and in all charges is smaller than convergence criterion

            energies = nenr
            charges = nchg

            if gradients:
                forces = nfrc

            atoms.set_array('charges', charges)

            np.savetxt(self.label + 'charges-%i.dat' %itr, charges, fmt='%21.10f')
            np.savetxt(self.label + 'energies-%i.dat' %itr, energies, fmt='%21.10f')
            if gradients:
                np.savetxt(self.label + 'forces-%i.dat' %itr, forces, fmt='%21.10f')

    #FU| convergence checks, we could also add individual check for subsystem energies
            if (np.abs(totdel) < self.econv) and all(np.abs(dchg) < self.qconv):
                if gradients:
                    if np.all(np.abs(dfrc) < 1.e-6):
                        converged = True
                else:
                    converged = True

            itr += 1

            if (not converged) and (itr == self.max_iter):
                print "Maximum number of iterations reached. Will stop now, although not converged"
                converged = True

            continue

        #if gradients:
        #    self.results['forces'] = -np.array(forces) * Ha / Bohr

        #FUDO| correctly set the energies here as well
        if self.plvl == -1 or self.plvl > 0:
            stdout.write('\nFinal total energy: %21.10f\n' %(totnr[-1]))

    #FUDO| make that available in a full calculator object
        self.results['energy'] = totnr[-1] * Ha
        self.results['charges'] = charges

    #FUDO| somehow make it work, that we get the charges later on also without re-calculating everything
        #atoms.set_initial_charges(charges)

    #FUDO| check units below
        if gradients:
            self.results['forces'] = -np.array(forces) * Ha / Bohr

    #FUMPORTANT| if this is included the whole thing works like a charm, if atoms and not self.atoms is passed in GEBF.calculate
        #self.atoms.set_initial_charges(self.results['charges'])

    def one_iteration(self, atoms=None, dobqs=True, gradients=False, chgs=None, mlts=None):

        #FUDO| figure out if we need to really copy the calculator or if we can just pass it?!
        #calc = self.subsyscalc.copy()
        #calc = atoms.get_calculator()

        nsub = len(self.subsys)

        if chgs is None:
            chgs = np.zeros(nsub)

        if mlts is None:
            mlts = np.ones(nsub)

        calc = deepcopy(self.subsyscalc)

#FUDO| this shit should not be necessary anymore
#FUDO| change later to passable parameter
        if gradients:
            calc.parameters['task'] = 'gradient'
        else:
            calc.parameters['task'] = 'energy'

        if dobqs:
            ochg = atoms.get_array('charges')

        charges = np.zeros(len(atoms))
        energies = np.zeros(len(self.subsys))
        forces = np.zeros((nsub, len(atoms), 3))
    
        for i, subsys in enumerate(self.subsys):
    #FUDO| option for this needed?
        #    os.system('rm inpsubsys-*')

            if self.plvl == -1 or self.plvl > 0:
                stdout.write(' %i / %i\r' %(i+1, nsub))
                stdout.flush()

    #FU| get indices of all that is not in subsys (works, checked it)
    #FU| and set the charges correspondingly
    #FUDO| check how general these parameters are to other calculators (if we'll ever use them)
            if dobqs:
                bqs = [elem for k, elem in enumerate(range(len(atoms))) if k not in subsys]
                calc.parameters.bq = np.append(atoms.positions[bqs], np.transpose(np.array(ochg[bqs], ndmin=2)), axis=1)

            calc.set_label(self.label + 'inpsubsys-%i' %i)
    #FUDO| check if the clean-up is really needed or not
    #FUDO| actually, right now everything is kept for the sake of simplicity unless the calculator is updated
            #os.system('rm %s ' %self.label + 'inpsubsys-%i' %i)
            os.system('rm %sinpsubsys-%i.*' %(self.label, i))
            calc.parameters['charge'] =  chgs[i]
            calc.parameters['mult'] = mlts[i]

            if mlts[i] != 1:
                calc.parameters['wfn'] = 'UHF'

    #FU| create temporary atoms object which contains only the subsystems
    #FU| and attach calculator

            tmpat = atoms[subsys].copy()
            tmpat.set_calculator(calc)

    #FU| this is not really needed, because it is done in order to get charges anyway
            #tmpat.get_potential_energy()
    #FU| assign new charges to final array with correct indices
            charges[subsys] = tmpat.get_charges()

            energies[i] = tmpat.get_potential_energy()
            
            if gradients:
                forces[i,subsys,:] = tmpat.get_forces()

    #FU| check how we set labels, it seems to converge quicker, but also to slightly different energies
            #system('rm -rf nwchem.*')

        if self.plvl > 0:
            stdout.write('\nDone!\n')

        if self.plvl > 0:
            stdout.write('\n')
        
        if gradients:
            return charges, energies, forces
        else:
            return charges, energies

    def get_mm_self_energy(self, atoms, charges):
        """calculate self-energy of the mm charges"""

        #we could do this here also with the atoms object

        numchg = len(charges)
        numatm = len(atoms)

        if numatm != numchg:
            print  "Number of charges and atoms not the same, cannot continue"
            #printf("Number of charges and atoms not the same, cannot continue\n")
            exit ( 1 )

        slfnr = 0.

    #check here, we need to reset distmat if we are doing more than only one GEBF run, e.g., for an MD run
    #FUDO| this is not perfect, because it demands, that self.atoms exists
        if self.distmat is None:
            self.distmat = atoms.get_distances_all()

            #self.distmat = np.zeros((numchg, numchg))
            #for i in range(numchg):
            #    for j in range(i+1, numchg):
            #        self.distmat[i,j] = atoms.get_distance(i, j)
            #        #distmat[j,i] = distmat[i,j]

        for i in range(numchg):
            for j in range(i+1, numchg):
                rij = self.distmat[i, j]
                #tmp = (charges[i] * charges[j]) # * u._e**2 / ( 4 * np.pi * u._eps0 )

                slfnr += ( charges[i] * charges[j] ) / rij * Bohr

                #this should be correct for atomic units

        return slfnr

    def remove_duplicates_lol(self, list):

        #although we do sort the list in the code above, we wanna be sure, that nothing goes wrong
        list.sort()

        duplicates = True
        cnt = 0
        cntprv = cnt-1
        uniq = []

        dblcnt = []

    #the loop below will break if we only have one atom fragments

        #print 'removal routine'
        #print list

        while duplicates:

            #print 'loop start: ', cnt, cnt+2, len(list)
            #here, it should be cnt+2, because we need to be able to access the element cnt+1
            if cnt+2 == len(list):
                #print 'We reached the end'
                duplicates = False

            if list[cnt] == list[cnt+1]:
                list.remove(list[cnt])
                #print 'counts: ', cnt, cntprv
                if cntprv == cnt:
                    #print cnt, len(dblcnt)
                    dblcnt[cnt] += 1
                    #print 'doubles: ', dblcnt
                else:
                    dblcnt.append(2)
                    cntprv = cnt
                    #print 'double count: ', dblcnt

            else:
                #print 'hello'
                dblcnt.append(1)
                cntprv = cnt
                cnt += 1

        #print dblcnt
        return dblcnt, list

#FU| we need the following methods:

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        """GEBF Calculator
        """

        Calculator.calculate(self, atoms, properties, system_changes)

        if len(system_changes) > 0:  # something wrong with this way
            self.update(self.atoms)
            if 'energy' in properties:
                self.calculate_energy(self.atoms)

            if 'forces' in properties:
                self.calculate_forces(self.atoms)

            if 'charges' in properties:
                self.calculate_charges(self.atoms)

        # check we have all the properties requested
        for property in properties:
            if property not in self.results:
                if property is 'energy':
                    self.calculate_energy(self.atoms)

                if property is 'forces':
                    self.calculate_forces(self.atoms)

                if property is 'charges':
                    self.calculate_charges(self.atoms)

#FUMPORTANT|
        atoms.set_initial_charges(self.results['charges'])
#FUNOMATTER|
        #self.atoms.set_initial_charges(self.results['charges'])

    def calculate_energy(self, atoms):
        self.run(atoms, gradients=False)

    def calculate_forces(self, atoms):
        self.run(atoms, gradients=True)

    def calculate_charges(self, atoms):
        self.calculate_energy(atoms)

    def update(self, atoms):
        self.pbc = atoms.get_pbc()
        os.system('rm %sinpsubsys-*' %self.label)

    #FU| the primary subsystems are determined by the cutoff, the derivative subsystems are defined by the primary subsystems
    #FU| so only need to the check the primary subsystems
        if self.prisub is None:
            prisub = self.generate_primary_subsystems(cutoff=3.0, around_solute=self.around_solute)

            if prisub != self.curprisub:
                self.coeffs, self.subsys, self.subsyschg, self.subsysmlt = self.update_subsystems(cutoff=self.cutoff, prisub=self.prisub, frgchg=self.frgchg, frgmlt=self.frgmlt, around_solute=self.around_solute)

    def set_atoms(self, atoms):
        if (atoms != self.atoms):
            self.converged = None

        self.atoms = atoms.copy()

        #atoms.set_array('charges', atoms.get_initial_charges())
        #self.atoms.set_array('charges', atoms.get_initial_charges())

    #FUDO| we might also need to update the actual subsystem definitions (completely)
    #FUDO| to this end split "generate_subsystems" into two functions. one for generation of primitive subsystems
    #FUDO| and one for the generation of derivative subsystems
#FUDO| we could include this function here from the original Calculator class, and make it remove 'initial_charges' from system_changes and substitute it by 
#FUDO| the actual charges among things to check
#    def check_state(self, atoms, tol=1e-15):
#        """Check for system changes since last calculation."""
#        if self.atoms is None:
#            system_changes = all_changes
#        else:
#            system_changes = []
#            if not equal(self.atoms.positions, atoms.positions, tol):
#                system_changes.append('positions')
#            if not equal(self.atoms.numbers, atoms.numbers):
#                system_changes.append('numbers')
#            if not equal(self.atoms.cell, atoms.cell, tol):
#                system_changes.append('cell')
#            if not equal(self.atoms.pbc, atoms.pbc):
#                system_changes.append('pbc')
#            if not equal(self.atoms.get_initial_magnetic_moments(),
#                         atoms.get_initial_magnetic_moments(), tol):
#                system_changes.append('initial_magmoms')
#            if not equal(self.atoms.get_initial_charges(),
#                         atoms.get_initial_charges(), tol):
#                system_changes.append('initial_charges')
#
#        return system_changes

