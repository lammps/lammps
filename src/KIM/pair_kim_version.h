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
   Contributing authors: Ryan S. Elliott,
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_KIM_VERSION_H
#define LMP_PAIR_KIM_VERSION_H

//
// Release: This file is part of the pair-kim-v1.7.2+1 package.
//

//
// This file defines the version information for the pair-kim package.
// The values specified here must conform to the Semantic Versioning
// 2.0.0 specification.
//
// Generally the version numbering for the pair-kim package will
// parallel the numbering for the kim-api package.  However, if
// additional versioning increments are required for the pair-kim
// package, the build-metatdata field will be used to provide a
// "sub-patch" version number.
//
// The PATCH value should be incremented IMMEDIATELY after an official
// release.
//
// The MINOR value should be incremented AND the PATCH value reset to
// zero as soon as it becomes clear that the next official release
// MUST increment the MINOR version value.
//
// The MAJOR value should be incremented AND the MINOR and PATCH
// vaules reset to zero as soon as it becomes clear that the next
// official release MUST increment the MAJOR version value.
//
// The PRERELEASE value can be set to any value allowed by the
// Semantic Versioning specification.  However, it will generally be
// empty.  This value should be quoted as a string constant.
//
// The BUILD_METADATA value can be set to any value allowed by the
// Semantic Versioning specification.  However, it will generally be
// emtpy; Except for when "sub-patch" versioning of the pair-kim
// package is necessary.  This value should be quoted as a string
// constant.
//

#define PAIR_KIM_VERSION_MAJOR 1
#define PAIR_KIM_VERSION_MINOR 7
#define PAIR_KIM_VERSION_PATCH 2
//#define PAIR_KIM_VERSION_PRERELEASE
#define PAIR_KIM_VERSION_BUILD_METADATA "1"

#endif  /* PAIR_KIM_VERSION_H */
