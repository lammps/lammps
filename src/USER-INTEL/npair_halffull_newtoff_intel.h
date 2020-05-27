/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

// Only used for hybrid to generate list for non-intel style. Use
// standard routines.

#ifdef NPAIR_CLASS

NPairStyle(halffull/newtoff/intel,
           NPairHalffullNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_HALF |
           NP_ORTHO | NP_TRI | NP_INTEL)

NPairStyle(halffull/newtoff/skip/intel,
           NPairHalffullNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_INTEL)

NPairStyle(halffull/newtoff/ghost/intel,
           NPairHalffullNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_HALF |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_INTEL)

NPairStyle(halffull/newtoff/skip/ghost/intel,
           NPairHalffullNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_GHOST | NP_INTEL)

#endif

/* ERROR/WARNING messages:

*/
