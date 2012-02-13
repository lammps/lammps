/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ERROR/WARNING messages:

E: Bond atoms %d %d missing on proc %d at step %ld

One or both of 2 atoms needed to compute a particular bond are
missing on this processor.  Typically this is because the pairwise
cutoff is set too short or the bond has blown apart and an atom is
too far away.

E: Angle atoms %d %d %d missing on proc %d at step %ld

One or more of 3 atoms needed to compute a particular angle are
missing on this processor.  Typically this is because the pairwise
cutoff is set too short or the angle has blown apart and an atom is
too far away.

E: Dihedral atoms %d %d %d %d missing on proc %d at step %ld

One or more of 4 atoms needed to compute a particular dihedral are
missing on this processor.  Typically this is because the pairwise
cutoff is set too short or the dihedral has blown apart and an atom is
too far away.

E: Improper atoms %d %d %d %d missing on proc %d at step %ld

One or more of 4 atoms needed to compute a particular improper are
missing on this processor.  Typically this is because the pairwise
cutoff is set too short or the improper has blown apart and an atom is
too far away.

*/
