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

#ifndef LMP_SERVER_MD_H
#define LMP_SERVER_MD_H

#include "pointers.h"

namespace LAMMPS_NS {

class ServerMD : protected Pointers {
 public:
  ServerMD(class LAMMPS *);
  ~ServerMD();
  void loop();

 private:
  int units;
  double fconvert, econvert, pconvert;
  double **fcopy;

  void box_change(double *, double *);
  void send_fev(int);
};

}    // namespace LAMMPS_NS

#endif
