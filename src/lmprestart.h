/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
// clang-format off
#ifndef LMP_RESTART_H
#define LMP_RESTART_H

#define MAGIC_STRING "LammpS RestartT"
#define ENDIAN 0x0001
#define ENDIANSWAP 0x1000
#define FORMAT_REVISION 3

enum{VERSION,SMALLINT,TAGINT,BIGINT,
     UNITS,NTIMESTEP,DIMENSION,NPROCS,PROCGRID,
     NEWTON_PAIR,NEWTON_BOND,
     XPERIODIC,YPERIODIC,ZPERIODIC,BOUNDARY,
     ATOM_STYLE,NATOMS,NTYPES,
     NBONDS,NBONDTYPES,BOND_PER_ATOM,
     NANGLES,NANGLETYPES,ANGLE_PER_ATOM,
     NDIHEDRALS,NDIHEDRALTYPES,DIHEDRAL_PER_ATOM,
     NIMPROPERS,NIMPROPERTYPES,IMPROPER_PER_ATOM,
     TRICLINIC,BOXLO,BOXHI,XY,XZ,YZ,
     SPECIAL_LJ,SPECIAL_COUL,
     MASS,PAIR,BOND,ANGLE,DIHEDRAL,IMPROPER,
     MULTIPROC,MPIIO,PROCSPERFILE,PERPROC,
     IMAGEINT,BOUNDMIN,TIMESTEP,
     ATOM_ID,ATOM_MAP_STYLE,ATOM_MAP_USER,ATOM_SORTFREQ,ATOM_SORTBIN,
     COMM_MODE,COMM_CUTOFF,COMM_VEL,NO_PAIR,
     EXTRA_BOND_PER_ATOM,EXTRA_ANGLE_PER_ATOM,EXTRA_DIHEDRAL_PER_ATOM,
     EXTRA_IMPROPER_PER_ATOM,EXTRA_SPECIAL_PER_ATOM,ATOM_MAXSPECIAL,
     NELLIPSOIDS,NLINES,NTRIS,NBODIES,ATIME,ATIMESTEP,LABELMAP};

#define LB_FACTOR 1.1

#endif
