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

#ifdef FIX_CLASS
// clang-format off
FixStyle(plumed,FixPlumed);
// clang-format on
#else

#ifndef LMP_FIX_PLUMED_H
#define LMP_FIX_PLUMED_H

#include "fix.h"

// forward declaration
namespace PLMD {
class Plumed;
}

namespace LAMMPS_NS {

class FixPlumed : public Fix {
 public:
  FixPlumed(class LAMMPS *, int, char **);
  ~FixPlumed() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  void reset_dt() override;
  int modify_param(int narg, char **arg) override;
  double memory_usage() override;

 private:
  PLMD::Plumed *p;           // pointer to plumed object
  int nlocal;                // number of atoms local to this process
  int natoms;                // total number of atoms
  int *gatindex;             // array of atom indexes local to this process
  double *masses;            // array of masses for local atoms
  double *charges;           // array of charges for local atoms
  int nlevels_respa;         // this is something to enable respa
  double bias;               // output bias potential
  class Compute *c_pe;       // Compute for the energy
  class Compute *c_press;    // Compute for the pressure
  int plumedNeedsEnergy;     // Flag to trigger calculation of the
                             // energy and virial
  char *id_pe, *id_press;    // ID for potential energy and pressure compute
};

};    // namespace LAMMPS_NS

#endif
#endif
