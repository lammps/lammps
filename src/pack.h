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

// loop counters for doing a pack/unpack 

struct pack_plan_3d {
  int nfast;                 // # of elements in fast index 
  int nmid;                  // # of elements in mid index 
  int nslow;                 // # of elements in slow index 
  int nstride_line;          // stride between successive mid indices 
  int nstride_plane;         // stride between successive slow indices
  int nqty;                  // # of values/element 
};

// function prototypes 

void pack_3d(double *, double *, struct pack_plan_3d *);
void unpack_3d(double *, double *, struct pack_plan_3d *);
void unpack_3d_permute1_1(double *, double *, struct pack_plan_3d *);
void unpack_3d_permute1_2(double *, double *, struct pack_plan_3d *);
void unpack_3d_permute1_n(double *, double *, struct pack_plan_3d *);
void unpack_3d_permute2_1(double *, double *, struct pack_plan_3d *);
void unpack_3d_permute2_2(double *, double *, struct pack_plan_3d *);
void unpack_3d_permute2_n(double *, double *, struct pack_plan_3d *);
