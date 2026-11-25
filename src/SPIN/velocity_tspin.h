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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(velocity/tspin,VelocityTSpin);
// clang-format on
#else

#ifndef LMP_VELOCITY_TSPIN_H
#define LMP_VELOCITY_TSPIN_H

#include "command.h"

namespace LAMMPS_NS {

class VelocityTSpin : public Command {
 public:
  VelocityTSpin(class LAMMPS *);

  void command(int, char **) override;

 private:
  int igroup, groupbit;
  int momentum_flag;
  int index_vs, index_sm;    // custom per-atom properties, see tspin.h
  int spinmass_flag;
  double spinmass;

  void create(double, int);
  void zero_mean(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
