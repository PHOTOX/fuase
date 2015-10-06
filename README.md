Some files to patch and mod the Atomic Simulation Environment

Modifications include ORCA, NWChem, GAMESS-US, and GEBF calculator.
Also the quick-xyz routine, which should be used with utter care. XYZ files for this reader need to have a fixed amount of characters/bytes per snapshot included in the XYZ file.

To install this patch, download first the ASE repository from their web page:
https://wiki.fysik.dtu.dk/ase/

set the variable ASEDIR to the location of your ASE directory (i.e., the one that contains the atoms.py file) and then run the patch.sh file from within the directory of this FUASE repository (important, otherwise file locations will be wrong)

For licensing information see file COPYING (GNU GPL v2)
