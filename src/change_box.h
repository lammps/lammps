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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(change_box,ChangeBox);
// clang-format on
#else

#ifndef LMP_CHANGE_BOX_H
#define LMP_CHANGE_BOX_H

#include "command.h"

namespace LAMMPS_NS {

class ChangeBox : public Command {
 public:
  ChangeBox(class LAMMPS *);
  void command(int, char **) override;

 private:
  int scaleflag;
  double scale[3];

  struct Operation {
    int style, flavor;
    int dim, boundindex;
    int vdim1, vdim2;
    double flo, fhi, ftilt;
    double dlo, dhi, dtilt;
    double scale;
  };

  Operation *ops;
  int nops;

  double boxlo[3], h_inv[6];

  void options(int, char **);
  void save_box_state();
  void volume_preserve(int, int, double);
};

}    // namespace LAMMPS_NS

#endif
#endif
