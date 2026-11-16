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

/* ------------------------------------------------------------------------
   Contributing author: Zhengtao Huang (The University of Hong Kong)
                        hzt990224@gmail.com
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nph/tspin,FixNPHTSpin);
// clang-format on
#else

#ifndef LMP_FIX_NPH_TSPIN_H
#define LMP_FIX_NPH_TSPIN_H

#include "fix_nh_tspin.h"

namespace LAMMPS_NS {

class FixNPHTSpin : public FixNHTSpin {
 public:
  FixNPHTSpin(class LAMMPS *, int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif
