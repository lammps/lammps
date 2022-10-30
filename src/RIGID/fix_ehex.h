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
FixStyle(ehex,FixEHEX);
// clang-format on
#else

#ifndef LMP_FIX_EHEX_H
#define LMP_FIX_EHEX_H

#include "fix.h"
#define EHEX_DEBUG 0

namespace LAMMPS_NS {

class FixEHEX : public Fix {

 public:
  FixEHEX(class LAMMPS *, int, char **);
  ~FixEHEX() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;
  void rescale();
  double compute_scalar() override;
  double memory_usage() override;
  void update_scalingmask();
  void com_properties(double *, double *, double *, double *, double *, double *);
  bool rescale_atom(int i, class Region *region);
  void grow_arrays(int nmax) override;
  bool check_cluster(tagint *shake_atom, int n, class Region *region);

 private:
  double heat_input;
  double masstotal;
  double scale;
  class Region *region;
  char *idregion;
  int me;

  double **x;                // coordinates
  double **f;                // forces
  double **v;                // velocities
  double *mass;              // masses
  double *rmass;             // reduced masses
  int *type;                 // atom types
  int nlocal;                // number of local atoms
  class FixShake *fshake;    // pointer to fix_shake/fix_rattle
  int constraints;           // constraints (0/1)
  int cluster;               // rescaling entire clusters (0/1)
  int hex;                   // HEX mode (0/1)
  bool *scalingmask;         // scalingmask[i] determines whether
                             // the velocity of atom i is to be rescaled
};

}    // namespace LAMMPS_NS

#endif
#endif
