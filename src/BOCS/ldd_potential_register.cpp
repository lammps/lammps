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
   Registration table for LDD potential sub-styles of pair_style ldd.

   LDD potential sub-styles cannot be loaded as plugins, so this is a
   plain checked-in table rather than an auto-generated style header.
   It replaces the "ldd_potential_styles.h" include list.

   To add a new LDD potential sub-style, only this file needs to be edited:

     (1) write the sub-style class in an ldd_potential_<name>.{h,cpp} file
         pair, deriving it from LddPotential,
     (2) add an "#include" for its header in the 1st customization section
         below, and
     (3) add one row to the ldd_potential_table[] array in the 2nd
         customization section below: the input-script keyword and the
         class name.
------------------------------------------------------------------------- */

#include "ldd_potential.h"
#include "pair_ldd.h"

// 1st customization section: add the #include for a new potential header here

#include "ldd_potential_constant.h"
#include "ldd_potential_linear.h"
#include "ldd_potential_mdpd.h"
#include "ldd_potential_noforce.h"
#include "ldd_potential_quadratic.h"
#include "ldd_potential_tablegradlin.h"
#include "ldd_potential_tablegradspline.h"
#include "ldd_potential_tablelin.h"
#include "ldd_potential_tablespline.h"

namespace LAMMPS_NS {

// factory function: one instance is created for each entry in the table below

template <typename T> static LddPotential *creator(LAMMPS *lmp)
{
  return new T(lmp);
}

// clang-format off

// 2nd customization section: add a row for a new potential here.
// The columns are the input-script keyword and the potential class name.

const LddPotentialInfo ldd_potential_table[] = {
  { "constant",         &creator<LddPotentialConstant>       },
  { "linear",           &creator<LddPotentialLinear>         },
  { "mdpd",             &creator<LddPotentialMdpd>           },
  { "noforce",          &creator<LddPotentialNoForce>        },
  { "quadratic",        &creator<LddPotentialQuadratic>      },
  { "table/gradlin",    &creator<LddPotentialTableGradLin>   },
  { "table/gradspline", &creator<LddPotentialTableGradSpline> },
  { "table/lin",        &creator<LddPotentialTableLin>       },
  { "table/spline",     &creator<LddPotentialTableSpline>    },
};
// clang-format on

const int num_ldd_potential = sizeof(ldd_potential_table) / sizeof(LddPotentialInfo);

}    // namespace LAMMPS_NS
