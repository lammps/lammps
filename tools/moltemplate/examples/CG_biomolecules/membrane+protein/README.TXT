 This example shows how to put a protein (inclusion) in a 
 lipid bilayer mixture composed of two different lipids (DPPC and DLPC).
 The DPPC lipid model is described here:
      G. Brannigan, P.F. Philips, and F.L.H. Brown,
      Physical Review E, Vol 72, 011915 (2005)
 (The DLPC model is a truncated version of DPPC. Modifications discussed below.)
 The protein model is described here:
      G. Bellesia, AI Jewett, and J-E Shea, 
      Protein Science, Vol19 141-154 (2010)

--- PREREQUISITES: ---

1) This example requires the "dihedral_style fourier", which is currently 
in the USER-MISC package.  Build LAMMPS with this package enabled using
   make yes-user-misc
before compiling LAMMPS. 
(See http://lammps.sandia.gov/doc/Section_start.html#start_3 for details.)

2) This example may require additional features to be added to LAMMPS.
If LAMMPS complains about an "Invalid pair_style", then
 a) download the "additional_lammps_code" from 
    http://moltemplate.org     (upper-left corner menu)
 b) unpack it
 c) copy the .cpp and .h files to the src folding of your lammps installation.
 d) (re)compile LAMMPS.


----- Details --------

This example contains a coarse-grained model of a 4-helix bundle protein
inserted into a lipid bilayer (made from a mixture of DPPC and DLPC).

    -- Protein Model: --

The coarse-grained protein is described in:
   G. Bellesia, AI Jewett, and J-E Shea, Protein Science, Vol19 141-154 (2010)
Here we use the "AUF2" model described in that paper.
(The hydrophobic beads face outwards.)

    -- Memebrane Model: --

The DPPC lipid bilayer described in:
     G. Brannigan, P.F. Philips, and F.L.H. Brown,
     Physical Review E, Vol 72, 011915 (2005)
and:
     M.C. Watson, E.S. Penev, P.M. Welch, and F.L.H. Brown
     J. Chem. Phys. 135, 244701 (2011)

As in Watson(JCP 2011), rigid bond-length constraints 
have been replaced by harmonic bonds. 

A truncated version of this lipid (named "DLPC") has also been added. 
The bending stiffness of each lipid has been increased to compensate 
for the additional disorder resulting from mixing two different types 
of lipids together.  (Otherwise pores appear.)
Unlike the original "DPPC" molecule model, the new "DPPC" and "DLPC" models 
have not been carefully parameterized to reproduce the correct behavior in 
a lipid bilayer mixture.


-------------

Instructions on how to build LAMMPS input files and 
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step2)
README_run.sh
