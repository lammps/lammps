The simulation consists of a mixture of isobutane and water.
Over time (less than 1 ns), the two molecules phase-separate.

The GAFF parameters are applied only to the isobutane molecule.
(The water molecule paramters are defined explicitly in
 src/moltemplate_force_fields/tip3p_2004.lt)

WARNING: THIS IS A PRELIMINARY EXAMPLE WHICH USES AMBER'S GAFF FORCE FIELD.
         AS OF 2014-4-19, these features have not been tested.
         THE ABILITY TO DETECT AND ASSIGN GAFF FORCE FIELD PARAMETERS
         MOLECULES ACCORDING TO ATOM TYPE IS AN EXPERIMENTAL FEATURE
         AND CURRENTLY PROBABLY HAS BUGS (IN THE DIHEDRALS AND IMPROPERS).
         PLEASE REPORT BUGS AND/OR SEND CORRECTIONS.  -A 2014-4-19

----------------- CHARGE ----------------------

NOTE: The GAFF force-field DOES NOT ASSIGN ATOM CHARGE.
      In this example, atom charges were taken from the OPLSAA force field file:
      http://dasher.wustl.edu/tinker/distribution/params/oplsaa.prm
      This is not the charge in AMBER simunlations is typically assigned.
      (As of 2014, it is assigned using the "HF/6-31G* RESP2" or "AM1-BCC3"
       methods using AmberTools (which are not available in moltemplate).
       http://ambermd.org/doc6/html/AMBER-sh-19.4.html
       http://ambermd.org/tutorials/basic/tutorial4b/)


-------- REQUIREMENTS: ---------

  This example requires building LAMMPS with the "USER-MISC" package.
  (because it makes use of "gaff.lt" which uses dihedral_style fourier)
   To do this, type "make yes-user-misc" before compiling LAMMPS.
  http://lammps.sandia.gov/doc/Section_start.html#start_3

More detailed instructions on how to build LAMMPS input files and 
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step 2)
README_run.sh
