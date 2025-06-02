EAM database tool

Fortran version (create.f) by Xiaowang Zhou (Sandia), xzhou at sandia.gov
with revisions by Lucas Hale lucas.hale at nist.gov from https://www.ctcms.nist.gov/potentials/entry/2004--Zhou-X-W-Johnson-R-A-Wadley-H-N-G--Al/

Python version (create_eam.py) by Germain Clavier germain.clavier at gmail.com

Most parameters based on this paper:

X. W. Zhou, R. A. Johnson, and H. N. G. Wadley, Phys. Rev. B, 69,
144113 (2004).

Parameters for Cr were taken from:
Lin Z B, Johnson R A and Zhigilei L V, Phys. Rev. B 77 214108 (2008)

This tool can be used to create an DYNAMO-formatted EAM
setfl file for alloy systems, using any combination
of the elements discussed in the paper and listed in
the EAM_code file, namely:

Cu, Ag, Au, Ni, Pd, Pt, Al, Pb, Fe, Mo, Ta, W, Mg, Co, Ti, Zr, Cr

WARNING: Please note that the parameter sets used here are all optimized
for the pure metals of the individual elements and that mixing rules will
be applied for creating the inter-element interactions.  Those are inferior
to models where the mixed terms were specifically optimized for particular
alloys.  Thus any potential files created with this tool should be used
with care and test calculations (e.g. on multiple binary mixtures) performed
to gauge the error.

Steps (create.f):

1) compile create.f -> a.out  (e.g. gfortran create.f)
2) edit the input file EAM.input to list 2 or more desired elements to include
3) a.out < EAM.input will create an *.eam.alloy potential file

Steps (create_eam.py):

Usage: create_eam.py [-h] [-n NAME [NAME ...]] [-nr NR] [-nrho NRHO]

options:
  -n NAME [NAME ...], --names NAME [NAME ...]
                        Element names
  -nr NR                Number of point in r space [default 2000]
  -nrho NRHO            Number of point in rho space [default 2000]

1) you must have numpy installed
2) run "python create_eam.py -n Ta Cu" with the list of desired elements
3) this will create an *.eam.alloy potential file

in DYNAMO or LAMMPS context the created file is referred to as a setfl file
  that can be used with the LAMMPS pair_style eam/alloy command
