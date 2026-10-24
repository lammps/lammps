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
FixStyle(langevin/tspin,FixLangevinTSpin);
// clang-format on
#else

#ifndef LMP_FIX_LANGEVIN_TSPIN_H
#define LMP_FIX_LANGEVIN_TSPIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLangevinTSpin : public Fix {
 public:
  FixLangevinTSpin(class LAMMPS *, int, char **);
  ~FixLangevinTSpin() override;

  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void reset_dt() override;

 protected:
  double t_start, t_stop, t_period, t_target;
  int zeroflag;
  int index_vs, index_sm;    // custom per-atom properties, see tspin.h
  int ilevel_respa;
  double gamma_drag, gamma_random;

  class RanMars *random;

  void compute_target();
};

}    // namespace LAMMPS_NS

#endif
#endif
