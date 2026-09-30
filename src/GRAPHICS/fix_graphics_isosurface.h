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
FixStyle(graphics/isosurface,FixGraphicsIsosurface);
// clang-format on
#else

#ifndef LMP_FIX_GRAPHICS_ISOSURFACE_H
#define LMP_FIX_GRAPHICS_ISOSURFACE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGraphicsIsosurface : public Fix {
 public:
  FixGraphicsIsosurface(class LAMMPS *, int, char **);
  ~FixGraphicsIsosurface() override;

  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  void end_of_step() override;

  int image(int *&, double **&) override;

 private:
  double iso;
  double rad;
  double *pdata;    // holding space for per-atom property data
  int quality;
  int nlevels_respa;
  int nmax;

  int pflag;               // type of property data
  int pindex;              // 1-based index if data is vector, else 0
  char *pstr;              // compute/fix/variable ID
  class Compute *pcomp;    // pointer to per-atom compute
  class Fix *pfix;         // pointer to per-atom fix
  int pvar;                // property variable index

  int binary;
  int pad;
  std::string filename;

  int numobjs;
  int *imgobjs;
  double **imgparms;
};
}    // namespace LAMMPS_NS
#endif
#endif
