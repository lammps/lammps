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

#ifndef LMP_INTEGRATE_H
#define LMP_INTEGRATE_H

#include "pointers.h"
#include "compute.h"

namespace LAMMPS_NS {

class Integrate : protected Pointers {
 public:
  Integrate(class LAMMPS *, int, char **);
  virtual void init();
  virtual void setup(int flag) = 0;
  virtual void setup_minimal(int) = 0;
  virtual void run(int) = 0;
  virtual void force_clear() = 0;
  virtual void cleanup() {}
  virtual void reset_dt() {}
  virtual double memory_usage() { return 0; }

 protected:
  int eflag, vflag;            // flags for energy/virial computation
  int virial_style;            // compute virial explicitly or implicitly
  int external_force_clear;    // clear forces locally or externally

  // lists of PE,virial Computes
  std::vector<Compute *> elist_global, elist_atom, vlist_global, vlist_atom, cvlist_atom;

  int pair_compute_flag;      // 0 if pair->compute is skipped
  int kspace_compute_flag;    // 0 if kspace->compute is skipped

  void ev_setup();
  void ev_set(bigint);
};

}    // namespace LAMMPS_NS

#endif
