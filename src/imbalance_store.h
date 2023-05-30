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

#ifndef LMP_IMBALANCE_STORE_H
#define LMP_IMBALANCE_STORE_H

#include "imbalance.h"

namespace LAMMPS_NS {

class ImbalanceStore : public Imbalance {
 public:
  ImbalanceStore(class LAMMPS *);
  ~ImbalanceStore() override;

 public:
  // parse options, return number of arguments consumed
  int options(int, char **) override;
  // compute per-atom imbalance and apply to weight array
  void compute(double *) override;
  // print information about the state of this imbalance compute (required)
  std::string info() override;

 private:
  char *name;    // property name
};

}    // namespace LAMMPS_NS

#endif
