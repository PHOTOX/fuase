"""Calculator for General Energy-based Fragmentation"""

from sys import stdout
from os import system
from ase.units import Ha, Bohr
import ase.units as u
#from future import printf
from copy import deepcopy
from ase.calculators.nwchem import NWChem
from numpy import zeros, append, array, transpose, savetxt, abs, pi

#maybe later we can make let it inherit from the Calculator base class
class GEBF(object):
    """Initializes the General Energy-based Fragmentation object

    Parameters:

    atoms: atoms (object)

    fragments: list of lists containing indices defining moleculesor generally fragments, can be optional and
               molecular fragments will be generated based on the covalent radii
               of the atoms in the atoms object

               if a single number is provided it is used as a fudge factor for
               the molecule generation
    
    gradients: boolean determining whether gradients are to be calculated

    fullgrad: boolean determining whether fully correct gradients are to be
              calculated or only the ones based on molecule by molecule embedding
    """

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
#    implemented_properties = ['energy'] #, forces]

    def __init__(self, atoms, fragments=None, gradients=False, fullgrad=False, pricalc=None):
        self.fragments = fragments

        self.pricalc = pricalc

        #check if GEBF is implemented for the calculator attached to the atoms object
        #also check if there is a caclulator object attached

#FU| that is probably the wrong way of doing it
        self.atoms = atoms
#just for now, to make it work below
#is there an ASE function to get distance matrix?
#does that yield the same distance matrix than our method?
        #self.distmat = None
        self.distmat = self.atoms.get_distances_all()
        #self.distmat /= Bohr

        if self.pricalc is None:
        #if atoms.get_calculator() is None:
            print  "Atoms object has no Calculator attached, cannot continue"
            #printf("Atoms object has no Calculator attached, cannot continue\n")
            exit ( 1 )

        if not any(isinstance(self.pricalc, e) for e in [NWChem]):
            print  "GEBF calculator is currently only supported with NWChem"
            #printf("GEBF calculator is currently only supported with NWChem\n")
            exit ( 1 )

        if fragments is None:
            self.fudge = 1.0
            self.generate_molecules(fudge=self.fudge)

        elif isinstance(fragments, float):
            self.fudge = fragments
            self.generate_molecules(fudge=self.fudge)

        else:
            self.fragments = fragments

    
    def generate_molecules(self, fudge=1.0):
        """Generate molecules based on distance comparison to covalent radii
           Todo: add consistency check with database molecules
        """

        self.fragments = []

    def generate_fragments(self):
        """Generate fragments according to original GEBF literature.
           Not yet clear how exactly but can do later
        """

#and maybe, but only maybe, we should here also use self.molecules, nah always call it fragments
        self.fragments = []

#or add atoms to GEBF, i.e., GEBF.atoms
#do we really need the fragments here a second time, after the GEBF has already been initialized with fragments argument?
    def run(self):

        if (self.fragments is None) or (self.fragments == []):
            print "Something is wrong, no fragments defined"
            #printf("Something is wrong, no fragments defined\n")
            exit ( 1 )

        #get initial guess on charges
        #there also exists a function atoms.set_initial_charges, but I don't know what it actually does... (need to check that, I mean ;-)

        charges = zeros(len(self.atoms))
        energies = zeros(len(self.fragments))

        #calc = self.atoms.get_calculator()
        #calc = self.pricalc.copy()
        calc = deepcopy(self.pricalc)

        nfrg = len(self.fragments)

#unit trouble, for now we'll do everythin in atomic units, i.e., Ha
        for i, frag in enumerate(self.fragments):
            stdout.flush()
            stdout.write('working on %i / %i\r' %(i+1, nfrg))
        
            tmpat = self.atoms[frag].copy()
            tmpat.set_calculator(calc)

            #tmpat.get_potential_energy()
            charges[frag] = tmpat.get_charges()
            energies[i] = tmpat.get_potential_energy() / Ha

        stdout.write('\n')

        self.atoms.set_array('charges', charges)
        
        savetxt('initial-charges.dat', charges)
        savetxt('initial-energies.dat', energies)

        #basically something like this:
        converged = False

        totnr = [energies.sum()]
        itr = 0

        print "Energy without surrounding point charges: ", totnr[itr]

        while not converged:
            #keep running

            stdout.write('\n')
            nchg, nenr = self.one_iteration()

#convert to a.u., because ASE returns eV
            nenr /= Ha

#here, I guess, we should use the newly acquired charges. They should be the correct thing to remove...
#var 2
#            slfnr = self.get_mm_self_energy(self.atoms, nchg)
#var 1
            slfnr = self.get_mm_self_energy(self.atoms, charges)

            stdout.write('self-energy of point charges: %21.10f\n' %slfnr)
            #stdout.write('self-energy of point charges: %e\n' %slfnr)

#fix issue with self-energy of charges
#var 2
#            totnr.append(nenr.sum() + slfnr)
#var 1
            totnr.append(nenr.sum() - (len(self.fragments)-1) * slfnr)
            totdel = totnr[-1] - totnr[itr]

            stdout.write('total energy in iteration %i: %21.10f\n' %(itr, totnr[-1]))

            dchg = nchg - charges
            denr = nenr - energies

#perform additional checks here, whether change in all energies and in all charges is smaller than convergence criterion

            energies = nenr
            charges = nchg
            self.atoms.set_array('charges', charges)

            savetxt('charges-%i.dat' %itr, charges)
            savetxt('energies-%i.dat' %itr, energies)

            if (abs(totdel) < 1.e-6) and all(abs(dchg) < 1.e-4):
                converged = True

            itr += 1

            continue

    def one_iteration(self):

        #calc = self.pricalc.copy()
        #is this really needed?
        calc = deepcopy(self.pricalc)
        #calc = self.atoms.get_calculator()
        ochg = self.atoms.get_array('charges')

        charges = zeros(len(self.atoms))
        energies = zeros(len(self.fragments))

        nfrg = len(self.fragments)

        for i, frag in enumerate(self.fragments):
#get indices of all that is not in frag (works, checked it)
            bqs = [elem for k, elem in enumerate(range(len(self.atoms))) if k not in frag]

#and set the charges correspondingly
            calc.parameters.bq = append(self.atoms.positions[bqs], transpose(array(ochg[bqs], ndmin=2)), axis=1)

#create temporary atoms object which contains only the fragments
            tmpat = self.atoms[frag].copy()
#attach calculator
            tmpat.set_calculator(calc)

#this is not really needed, because it is done in order to get charges anyway
            #tmpat.get_potential_energy()
#assign new charges to final array
            charges[frag] = tmpat.get_charges()

#use temporary array of charges for calculation of self-energy
            chgtmp = ochg.copy()
#now is the question which charges we should use for the fragment at hand
#if any, at all
#var 2
#            chgtmp[frag] = charges[frag]
#var 1
            chgtmp[frag] = 0.

            slfnr = self.get_mm_self_energy(self.atoms, chgtmp)

#remove only the self-energy of the interacting bq, but not the QM density with the Bq's
#this leads to double-counting
            #energies[i] = tmpat.get_potential_energy() - slfnr
            energies[i] = tmpat.get_potential_energy()

#check here, and do this somehow differently and also possibly while still keeping separate files for each fragment to make initial guesses better and quicker
            system('rm -rf nwchem.*')
            stdout.flush()
            stdout.write(' %i / %i done\r' %(i+1, nfrg))
        
        stdout.write('\n')
        
        #self.atoms.set_array('charges', charges)

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
        if self.distmat is None:
            self.distmat = zeros((numchg, numchg))
            for i in range(numchg):
                for j in range(i+1, numchg):
                    self.distmat[i,j] = atoms.get_distance(i, j)
                    #distmat[j,i] = distmat[i,j]

        for i in range(numchg):
            for j in range(i+1, numchg):
                rij = self.distmat[i, j]
                #tmp = (charges[i] * charges[j]) # * u._e**2 / ( 4 * pi * u._eps0 )

                slfnr += ( charges[i] * charges[j] ) / rij * Bohr

                #this should be correct for atomic units

        return slfnr
