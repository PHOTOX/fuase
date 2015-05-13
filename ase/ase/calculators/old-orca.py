import os, sys
import numpy as np
from ase.units import Hartree, Bohr
import inspect

class ORCA:
    """class for performing ORCA calculations, without inclusion of the
    ORCA python wrapper"""

    def __init__(self,label='orca',**kwargs):
        """input file hack will be put in here right now"""

        self.label = label
        self.converged = False

        self.stress = np.empty((3, 3))

        self.base_dir = os.getcwd()

        try:
            pathtoorca = kwargs.pop('pathtoorca')
            self.pathtoorca = pathtoorca
        except:
            self.pathtoorca = '/home/frank/prog/orca/orca_x86_64_exe_r2360/orca'

        try:
            calc_dir = kwargs.pop('calc_dir')
            self.calc_dir = self.base_dir + '/' + calc_dir
        except:
            self.calc_dir = self.base_dir

        try:
            prefix = kwargs.pop('prefix')
            self.prefix = prefix
        except:
            self.prefix = 'tmp'

        try:
            basis_set = kwargs.pop('basis_set')
            self.basis_set = basis_set
        except:
            self.basis_set = 'TZVP'
            #nwchem has no stress ;-)
        
        try:
            method = kwargs.pop('method')
            self.method = method
        except:
            self.method = 'DFT'

        try:
            riapprox = kwargs.pop('riapprox')
            self.riapprox = riapprox
        except:
            self.riapprox = False

        try:
            printoptions = kwargs.pop('printoptions')
            self.printoptions = printoptions
        except:
            self.printoptions = ''

        try:
            ncpu = kwargs.pop('ncpu')
            self.ncpu = int(ncpu)
        except:
            self.ncpu = 1

        try:
            charge = kwargs.pop('charge')
            self.charge = int(charge)
        except:
            self.charge = 0

        try:
            multiplicity = kwargs.pop('multiplicity')
            self.multiplicity = int(multiplicity)
        except:
            self.multiplicity = 1

        try:
            dummies = kwargs.pop('dummies')
            self.dummies = dummies
        except:
            self.dummies = None

        try:
            stuff = kwargs.pop('stuff')
            self.stuff = stuff
        except:
            self.stuff = ''

        try:
            spindens = kwargs.pop('spindens')
            self.spindens = spindens
        except:
            self.spindens = False

    def update(self, atoms):
        if (not self.converged or
            len(self.numbers) != len(atoms) or
            (self.numbers != atoms.get_atomic_numbers()).any()):
            self.initialize(atoms)
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.calculate(atoms)

    def initialize(self, atoms):
        self.numbers = atoms.get_atomic_numbers().copy()
        self.species = []
        for a, Z in enumerate(self.numbers):
            self.species.append(Z)
        self.converged = False
        
    def get_potential_energy(self, atoms):
        self.energy = False
        self.gradient = True
        #self.gradient = False
        self.update(atoms)
        return self.etotal

    def get_total_energy(self, atoms):
        self.energy = True
        self.gradient = False
        self.update(atoms)
        return self.etotal

    def get_forces(self, atoms):
        self.energy = False
        self.gradient = True
        self.update(atoms)
        return self.forces.copy()
    
    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress.copy()

    def calculate(self, atoms):

        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()

        if self.calc_dir != self.base_dir:
            try:
                os.mkdir(self.calc_dir)
            except:
                print "dir '%s' already exists - continuing" %self.calc_dir

        if not self.spindens:
            os.chdir(self.calc_dir)
            self.write_orca_input(atoms)
            os.system('%s %s.input > %s.output' %(self.pathtoorca, self.prefix, self.prefix))
            os.system('cat %s.output >> all.output' %self.prefix)
            os.chdir(self.base_dir)
        else:
            self.calculate_spin_density(atoms, "sdens.cube")

        os.chdir(self.calc_dir)
        self.read_energy()
        if self.gradient:
            self.read_forces(atoms)

        os.chdir(self.base_dir)
        self.converged = True

    def calculate_spin_density(self, atoms, filename):
        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()

        if self.calc_dir != self.base_dir:
            try:
                os.mkdir(self.calc_dir)
            except:
                print "dir '%s' already exists - continuing" %self.calc_dir

        os.chdir(self.calc_dir)
        interm = self.stuff
        self.stuff += '\n%%plots\nFormat Gaussian_Cube\nSpinDens("%s");\nend\n' %filename
        self.write_orca_input(atoms)
        os.system('%s %s.input > %s.output' %(self.pathtoorca, self.prefix, self.prefix))
        self.stuff = interm

        os.chdir(self.base_dir)

    def read_energy(self):
        f = open('%s.output' %self.prefix,'r')
        lines = f.readlines()
        f.close()

        for line in range(len(lines)):
            if lines[line].startswith('FINAL SINGLE POINT ENERGY'):
                self.etotal = float(lines[line].split()[-1])*Hartree

    def read_forces(self,atoms): #todo
        f = open('%s.output' %self.prefix,'r')
        lines = f.readlines()
        f.close()

        self.forces = np.array([],dtype='float64')
        if (self.method == 'HF') or (self.method == 'DFT'):
            grad = lines.index('CARTESIAN GRADIENT\n')
        elif self.method.endswith('MP2'):
            grad = lines.index('The final MP2 gradient\n')-2

        forces = lines[grad+3:(grad+3+np.shape(atoms.positions)[0])]
        #print forces[0].split()

        self.forces = np.empty((len(forces),3))

        i = 0
        for line in forces:
            if (self.method == 'HF') or (self.method == 'DFT'):
                grads = line.split()[3:6]
            elif self.method.endswith('MP2'):
                grads = line.split()[1:4]
            #print grads
            gradx = np.float(grads[0])
            grady = np.float(grads[1])
            gradz = np.float(grads[2])

            self.forces[i,:] = -(np.array([gradx, grady, gradz]))*Hartree/Bohr
            #self.forces[i,:] = np.array([gradx, grady, gradz]) #pure forces from ORCA
            i += 1

        #print self.forces

        #check for units of forces and sign (prefactor TM=Hartree/Bohr)

    def read(self):
        """Dummy stress for NWChem"""
        self.stress = np.empty((3, 3))

    def write_orca_input(self, atoms):
        """writes orca input inputfile consisting of method, basis set and print options"""

        inputfile = open('%s.input' %self.prefix, 'w')

        if self.riapprox:
            riprefix = 'RI'
        else:
            riprefix = ''

        inputfile.write('! %s %s %s\n' %(riprefix, self.method, self.basis_set))
        
        #if self.gradient:
        print 'running gradient'
        inputfile.write('! EnGrad TightSCF\n')

        if self.ncpu > 1:
            inputfile.write('%s \n nprocs %s \n end \n' %(str('%pal'), self.ncpu))

        inputfile.write('\n'+self.stuff+'\n')

        symbols = atoms.get_chemical_symbols()

        inputfile.write('*xyz %s %s \n' %(str(self.charge), str(self.multiplicity)))
        for i in range(np.shape(atoms.positions)[0]):
            inputfile.write('%s %f %f %f\n' %(symbols[i], atoms.positions[i,0], atoms.positions[i,1], atoms.positions[i,2]))

        if self.dummies != None:
            inputfile.write(self.dummies)

        inputfile.write('*')

        inputfile.close()
