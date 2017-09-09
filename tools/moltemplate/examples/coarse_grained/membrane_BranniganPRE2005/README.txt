Note:
   This example may require additional features to be added to LAMMPS.
If LAMMPS complains about an "Invalid pair_style", then copy the code
in the "additional_lammps_code" directory into your LAMMPS "src" directory
and recompile LAMMPS.

----- Description --------

This example contains an implementation of the DPPC lipid bilayer described in:
     G. Brannigan, P.F. Philips, and F.L.H. Brown,
     Physical Review E, Vol 72, 011915 (2005)
and:
     M.C. Watson, E.S. Penev, P.M. Welch, and F.L.H. Brown
     J. Chem. Phys. 135, 244701 (2011)

As in Watson(JCP 2011), rigid bond-length constraints
have been replaced by harmonic bonds.

A truncated version of this lipid (named "DLPC") has also been added.
Unlike the original "DPPC" molecule model, "DLPC" has not been carefully
parameterized to reproduce the correct behavior in a lipid bilayer.


-------------

Instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step2)
README_run.sh
