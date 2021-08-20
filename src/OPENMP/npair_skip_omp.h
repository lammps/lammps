/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// There is no benefit from multi-threading for skip lists, so we
// just forward the requests to the corresponding non-omp versions.

#ifdef NPAIR_CLASS
// clang-format off
NPairStyle(skip/omp,
           NPairSkip,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_OMP);

NPairStyle(skip/half/respa/omp,
           NPairSkipRespa,
           NP_SKIP | NP_RESPA | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_OMP);

NPairStyle(skip/half/size/omp,
           NPairSkipSize,
           NP_SKIP | NP_SIZE | NP_HALF | NP_FULL | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_OMP);

NPairStyle(skip/size/off2on/omp,
           NPairSkipSizeOff2on,
           NP_SKIP | NP_SIZE | NP_OFF2ON | NP_HALF |
           NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_OMP);

NPairStyle(skip/size/off2on/oneside/omp,
           NPairSkipSizeOff2onOneside,
           NP_SKIP | NP_SIZE | NP_OFF2ON | NP_ONESIDE | NP_HALF |
           NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_NEWTON | NP_NEWTOFF |
           NP_ORTHO | NP_TRI | NP_OMP);

NPairStyle(skip/ghost/omp,
           NPairSkip,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_OMP | NP_GHOST);
// clang-format off
#endif

/* ERROR/WARNING messages:

*/
