/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ------------------------------------------------------
    This file is part of the BOCS package for LAMMPS.
   ------------------------------------------------------ */

/* ----------------------------------------------------------------------
   Registration table for LDD indicator sub-styles of pair_style ldd.

   LDD indicator sub-styles cannot be loaded as plugins, so this is a
   plain checked-in table rather than an auto-generated style header.
   It replaces the "ldd_indicator_styles.h" include list.

   To add a new LDD indicator sub-style, only this file needs to be edited:

     (1) write the sub-style class in an ldd_indicator_<name>.{h,cpp} file
         pair, deriving it from LddIndicator,
     (2) add an "#include" for its header in the 1st customization section
         below, and
     (3) add one row to the ldd_indicator_table[] array in the 2nd
         customization section below: the input-script keyword and the
         class name.
------------------------------------------------------------------------- */

#include "ldd_indicator.h"
#include "pair_ldd.h"

// 1st customization section: add the #include for a new indicator header here

#include "ldd_indicator_dpd.h"
#include "ldd_indicator_lucy.h"
#include "ldd_indicator_shell.h"
#include "ldd_indicator_smooth.h"
#include "ldd_indicator_sphere.h"

namespace LAMMPS_NS {

// factory function: one instance is created for each entry in the table below

template <typename T> static LddIndicator *creator(LAMMPS *lmp)
{
  return new T(lmp);
}

// clang-format off

// 2nd customization section: add a row for a new indicator here.
// The columns are the input-script keyword and the indicator class name.

LMP_REGISTRY_CONST LddIndicatorInfo ldd_indicator_table[] = {
  { "dpd",    &creator<LddIndicatorDpd>    },
  { "lucy",   &creator<LddIndicatorLucy>   },
  { "shell",  &creator<LddIndicatorShell>  },
  { "smooth", &creator<LddIndicatorSmooth> },
  { "sphere", &creator<LddIndicatorSphere> },
};
// clang-format on

LMP_REGISTRY_CONST int num_ldd_indicator = sizeof(ldd_indicator_table) / sizeof(LddIndicatorInfo);

}    // namespace LAMMPS_NS
