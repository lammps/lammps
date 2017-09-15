 This example shows how to build a multicomponent spherical vesicle.
 The lipid bilayer is composed of two different lipids (DPPC and DLPC).
 The vesicle also contains trans-membrane protein inclusions.

 The coordinates for the vesicle are constructed by PACKMOL (see below).

 The DPPC lipid model is described here:
      G. Brannigan, P.F. Philips, and F.L.H. Brown,
      Physical Review E, Vol 72, 011915 (2005)
 (The DLPC model is a truncated version of DPPC. Modifications discussed below.)
 The protein model is described here:
      G. Bellesia, AI Jewett, and J-E Shea,
      Protein Science, Vol19 141-154 (2010)

--- PREREQUISITES: ---

1) This example requires PACKMOL.  You can download PACKMOL here:

  http://www.ime.unicamp.br/~martinez/packmol/

 (Moltemplate does not come with an easy way to generate spherically-symmetric
 structures, so I used the PACKMOL program to move the molecules into position.)

2) This example requires the "dihedral_style fourier", which is currently
in the USER-MISC package.  Build LAMMPS with this package enabled using
   make yes-user-misc
before compiling LAMMPS.
(See http://lammps.sandia.gov/doc/Section_start.html#start_3 for details.)

3) This example may require additional features to be added to LAMMPS.
If LAMMPS complains about an "Invalid pair_style", then
 a) download the "additional_lammps_code" from
    http://moltemplate.org     (upper-left corner menu)
 b) unpack it
 c) copy the .cpp and .h files to the src folding of your lammps installation.
 d) (re)compile LAMMPS.

------ Details -------

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

    -- Interactions between the proteins and lipids --

This is discussed in the "system.lt" file.

--- Building the files necessary to run a simulation in LAMMPS ---

step 1) Run PACKMOL

        Type these commands into the shell.
        (Each command could take several hours.)

cd packmol_files
  packmol < step1_proteins.inp   # This step determines the protein's location
  packmol < step2_innerlayer.inp # this step builds the inner monolayer
  packmol < step3_outerlayer.inp # this step builds the outer monolayer
cd ..

step 2) Run MOLTEMPLATE
        Type these commands into the shell.
        (This could take up to 10 minutes.)

cd moltemplate_files
  moltemplate.sh system.lt -xyz ../system.xyz
  mv -f system.in* system.data ../
  cp -f table_int.dat ../
cd ..

--- Running LAMMPS ---

step3) Run LAMMPS:
        Type these commands into the shell.
        (This could take days.)

lmp_linux -i run.in.min  # Minimize the system (important, and very slow)

lmp_linux -i run.in.nvt  # Run a simulation at constant volume

If you have compiled the MPI version of lammps, you can run lammps in parallel:

mpirun -np 4 lmp_linux -i run.in.min
  or
mpirun -np 4 lmp_linux -i run.in.nvt

(Assuming you have 4 cores, for example.)

