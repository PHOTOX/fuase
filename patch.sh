#!/bin/bash

OWD=`pwd`
ASEDIR=/scratch/fuhlig/opt/ase/trunk/ase

cp calculators/{gebf.py,orca.py,nwchem.py,gamess_us.py} ${ASEDIR}/calculators/
cp io/{orca.py,gamess_us.py,nwchem.py,qxyz.py} ${ASEDIR}/io/

cd ${ASEDIR}
patch -p0 <${OWD}/atoms.patch

cd ${OWD}

#to clean up the just messed up repository go to $ASEDIR and run:
#   svn revert --recursive ./
#in the top level ase directory. 
#Be aware that this will revert all changes that have not been commited, also those that have not been done by this file.
