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

// NeighBin child classes
// customize when adding a new class
// order does not matter, but numbers must be unique
// binname[] numbers must be consistent with enum

enum{
  NEIGH_BIN_STANDARD=1,

  NEIGH_BIN_INTEL=201,

  NEIGH_BIN_SSA=301,
};

struct BinName {
  const char *word;
  int flag;
} binname[] = {
  "standard", 1,

  "intel", 201,

  "ssa", 301,

  "", -1
};

// NeighStencil child classes
// customize when adding a new class
// order does not matter, but numbers must be unique
// stencilname[] numbers must be consistent with enum

enum{
  NEIGH_STENCIL_FULL_BIN_2D=1,
  NEIGH_STENCIL_FULL_BIN_3D=2,
  NEIGH_STENCIL_FULL_GHOST_BIN_2D=3,
  NEIGH_STENCIL_FULL_GHOST_BIN_3D=4,
  NEIGH_STENCIL_FULL_MULTI_2D=5,
  NEIGH_STENCIL_FULL_MULTI_3D=6,
  NEIGH_STENCIL_HALF_BIN_2D_NEWTON=7,
  NEIGH_STENCIL_HALF_BIN_2D_NEWTON_TRI=8,
  NEIGH_STENCIL_HALF_BIN_2D_NO_NEWTON=9,
  NEIGH_STENCIL_HALF_BIN_3D_NEWTON=10,
  NEIGH_STENCIL_HALF_BIN_3D_NEWTON_TRI=11,
  NEIGH_STENCIL_HALF_BIN_3D_NO_NEWTON=12,
  NEIGH_STENCIL_HALF_GHOST_BIN_2D_NO_NEWTON=13,
  NEIGH_STENCIL_HALF_GHOST_BIN_3D_NO_NEWTON=14,
  NEIGH_STENCIL_HALF_MULTI_2D_NEWTON=15,
  NEIGH_STENCIL_HALF_MULTI_2D_NEWTON_TRI=16,
  NEIGH_STENCIL_HALF_MULTI_2D_NO_NEWTON=17,
  NEIGH_STENCIL_HALF_MULTI_3D_NEWTON=18,
  NEIGH_STENCIL_HALF_MULTI_3D_NEWTON_TRI=19,
  NEIGH_STENCIL_HALF_MULTI_3D_NO_NEWTON=20,

  NEIGH_STENCIL_HALF_BIN_2D_SSA=301,
  NEIGH_STENCIL_HALF_BIN_3D_SSA=302,
};

struct StencilName {
  const char *word;
  int flag;
} stencilname[] = {
  "full_bin_2d", 1,
  "full_bin_3d", 2,
  "full_ghost_bin_2d", 3,
  "full_ghost_bin_3d", 4,
  "full_multi_2d", 5,
  "full_multi_3d", 6,
  "half_bin_2d_newton", 7,
  "half_bin_2d_newton_tri", 8,
  "half_bin_2d_no_newton", 9,
  "half_bin_3d_newton", 10,
  "half_bin_3d_newton_tri", 11,
  "half_bin_3d_no_newton", 12,
  "half_ghost_bin_2d_no_newton", 13,
  "half_ghost_bin_3d_no_newton", 14,
  "half_multi_2d_newton", 15,
  "half_multi_2d_newton_tri", 16,
  "half_multi_2d_no_newton", 17,
  "half_multi_3d_newton", 18,
  "half_multi_3d_newton_tri", 19,
  "half_multi_3d_no_newton", 20,

  "half_bin_2d_ssa", 301,
  "half_bin_3d_ssa", 302,

  "", -1
};

// NeighPair child classes
// customize when adding a new class
// order does not matter, but numbers must be unique
// pairname[] numbers must be consistent with enum

enum{
  NEIGH_PAIR_COPY_FROM=1,
  NEIGH_PAIR_FULL_BIN=2,
  NEIGH_PAIR_FULL_BIN_GHOST=3,
  NEIGH_PAIR_FULL_MULTI=4,
  NEIGH_PAIR_FULL_NSQ=5,
  NEIGH_PAIR_FULL_NSQ_GHOST=6,
  NEIGH_PAIR_GRANULAR_BIN_NEWTON=7,
  NEIGH_PAIR_GRANULAR_BIN_NEWTON_ONESIDED=8,
  NEIGH_PAIR_GRANULAR_BIN_NEWTON_TRI=9,
  NEIGH_PAIR_GRANULAR_BIN_NO_NEWTON=10,
  NEIGH_PAIR_GRANULAR_NSQ_NEWTON=11,
  NEIGH_PAIR_GRANULAR_NSQ_NEWTON_ONESIDED=12,
  NEIGH_PAIR_GRANULAR_NSQ_NO_NEWTON=13,
  NEIGH_PAIR_HALF_BIN_NEWTON=14,
  NEIGH_PAIR_HALF_BIN_NEWTON_TRI=15,
  NEIGH_PAIR_HALF_BIN_NO_NEWTON=16,
  NEIGH_PAIR_HALF_BIN_NO_NEWTON_GHOST=17,
  NEIGH_PAIR_HALF_FROM_FULL_NEWTON=18,
  NEIGH_PAIR_HALF_FROM_FULL_NO_NEWTON=19,
  NEIGH_PAIR_HALF_MULTI_NEWTON=20,
  NEIGH_PAIR_HALF_MULTI_NEWTON_TRI=21,
  NEIGH_PAIR_HALF_MULTI_NO_NEWTON=22,
  NEIGH_PAIR_HALF_NSQ_NEWTON=23,
  NEIGH_PAIR_HALF_NSQ_NO_NEWTON=24,
  NEIGH_PAIR_HALF_NSQ_NO_NEWTON_GHOST=25,
  NEIGH_PAIR_RESPA_BIN_NEWTON=26,
  NEIGH_PAIR_RESPA_BIN_NEWTON_TRI=27,
  NEIGH_PAIR_RESPA_BIN_NO_NEWTON=28,
  NEIGH_PAIR_RESPA_NSQ_NEWTON=29,
  NEIGH_PAIR_RESPA_NSQ_NO_NEWTON=30,
  NEIGH_PAIR_SKIP_FROM=31,
  NEIGH_PAIR_SKIP_FROM_GRANULAR=32,
  NEIGH_PAIR_SKIP_FROM_GRANULAR_OFF2ON=33,
  NEIGH_PAIR_SKIP_FROM_GRANULAR_OFF2ON_ONESIDED=34,
  NEIGH_PAIR_SKIP_FROM_RESPA=35,

  NEIGH_PAIR_FULL_BIN_GHOST_OMP=101,
  NEIGH_PAIR_FULL_BIN_OMP=102,
  NEIGH_PAIR_FULL_MULTI_OMP=103,
  NEIGH_PAIR_FULL_NSQ_GHOST_OMP=104,
  NEIGH_PAIR_FULL_NSQ_OMP=105,
  NEIGH_PAIR_GRANULAR_BIN_NEWTON_OMP=106,
  NEIGH_PAIR_GRANULAR_BIN_NEWTON_TRI_OMP=107,
  NEIGH_PAIR_GRANULAR_BIN_NO_NEWTON_OMP=108,
  NEIGH_PAIR_GRANULAR_NSQ_NEWTON_OMP=109,
  NEIGH_PAIR_GRANULAR_NSQ_NO_NEWTON_OMP=110,
  NEIGH_PAIR_HALF_BIN_NEWTON_OMP=111,
  NEIGH_PAIR_HALF_BIN_NEWTON_TRI_OMP=112,
  NEIGH_PAIR_HALF_BIN_NO_NEWTON_GHOST_OMP=113,
  NEIGH_PAIR_HALF_BIN_NO_NEWTON_OMP=114,
  NEIGH_PAIR_HALF_FROM_FULL_NEWTON_OMP=115,
  NEIGH_PAIR_HALF_FROM_FULL_NO_NEWTON_OMP=116,
  NEIGH_PAIR_HALF_MULTI_NEWTON_OMP=117,
  NEIGH_PAIR_HALF_MULTI_NEWTON_TRI_OMP=118,
  NEIGH_PAIR_HALF_MULTI_NO_NEWTON_OMP=119,
  NEIGH_PAIR_HALF_NSQ_NEWTON_OMP=120,
  NEIGH_PAIR_HALF_NSQ_NO_NEWTON_GHOST_OMP=121,
  NEIGH_PAIR_HALF_NSQ_NO_NEWTON_OMP=122,
  NEIGH_PAIR_RESPA_BIN_NEWTON_OMP=123,
  NEIGH_PAIR_RESPA_BIN_NEWTON_TRI_OMP=124,
  NEIGH_PAIR_RESPA_BIN_NO_NEWTON_OMP=125,
  NEIGH_PAIR_RESPA_NSQ_NEWTON_OMP=126,
  NEIGH_PAIR_RESPA_NSQ_NO_NEWTON_OMP=127,

  NEIGH_PAIR_FULL_BIN_INTEL=201,
  NEIGH_PAIR_HALF_BIN_NEWTON_INTEL=202,
  NEIGH_PAIR_HALF_BIN_NEWTON_TRI_INTEL=203,
  NEIGH_PAIR_HALF_BIN_NO_NEWTON_INTEL=204,

  NEIGH_PAIR_HALF_BIN_NEWTON_SSA=301,
  NEIGH_PAIR_HALF_FROM_FULL_NEWTON_SSA=302,
};

struct PairName {
  const char *word;
  int flag;
} pairname[] = {
  "copy_from", 1,
  "full_bin", 2,
  "full_bin_ghost", 3,
  "full_multi", 4,
  "full_nsq", 5,
  "full_nsq_ghost", 6,
  "granular_bin_newton", 7,
  "granular_bin_newton_onesided", 8,
  "granular_bin_newton_tri", 9,
  "granular_bin_no_newton", 10,
  "granular_nsq_newton", 11,
  "granular_nsq_newton_onesided", 12,
  "granular_nsq_no_newton", 13,
  "half_bin_newton", 14,
  "half_bin_newton_tri", 15,
  "half_bin_no_newton", 16,
  "half_bin_no_newton_ghost", 17,
  "half_from_full_newton", 18,
  "half_from_full_no_newton", 19,
  "half_multi_newton", 20,
  "half_multi_newton_tri", 21,
  "half_multi_no_newton", 22,
  "half_nsq_newton", 23,
  "half_nsq_no_newton", 24,
  "half_nsq_no_newton_ghost", 25,
  "respa_bin_newton", 26,
  "respa_bin_newton_tri", 27,
  "respa_bin_no_newton", 28,
  "respa_nsq_newton", 29,
  "respa_nsq_no_newton", 30,
  "skip_from", 31,
  "skip_from_granular", 32,
  "skip_from_granular_off2on", 33,
  "skip_from_granular_off2on_onesided", 34,
  "skip_from_respa", 35,

  "full_bin_ghost_omp", 101,
  "full_bin_omp", 102,
  "full_multi_omp", 103,
  "full_nsq_ghost_omp", 104,
  "full_nsq_omp", 105,
  "granular_bin_newton_omp", 106,
  "granular_bin_newton_tri_omp", 107,
  "granular_bin_no_newton_omp", 108,
  "granular_nsq_newton_omp", 109,
  "granular_nsq_no_newton_omp", 110,
  "half_bin_newton_omp", 111,
  "half_bin_newton_tri_omp", 112,
  "half_bin_no_newton_ghost_omp", 113,
  "half_bin_no_newton_omp", 114,
  "half_from_full_newton_omp", 115,
  "half_from_full_no_newton_omp", 116,
  "half_multi_newton_omp", 117,
  "half_multi_newton_tri_omp", 118,
  "half_multi_no_newton_omp", 119,
  "half_nsq_newton_omp", 120,
  "half_nsq_no_newton_ghost_omp", 121,
  "half_nsq_no_newton_omp", 122,
  "respa_bin_newton_omp", 123,
  "respa_bin_newton_tri_omp", 124,
  "respa_bin_no_newton_omp", 125,
  "respa_nsq_newton_omp", 126,
  "respa_nsq_no_newton_omp", 127,

  "full_bin_intel", 201,
  "half_bin_newton_intel", 202,
  "half_bin_newton_tri_intel", 203,
  "half_bin_no_newton_intel", 204,

  "half_bin_newton_ssa", 301,
  "half_from_full_newton_ssa", 302,

  "", -1
};
