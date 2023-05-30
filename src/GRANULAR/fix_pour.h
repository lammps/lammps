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
FixStyle(pour,FixPour);
// clang-format on
#else

#ifndef LMP_FIX_POUR_H
#define LMP_FIX_POUR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPour : public Fix {
 public:
  FixPour(class LAMMPS *, int, char **);
  ~FixPour() override;
  int setmask() override;
  void init() override;
  void setup_pre_exchange() override;
  void pre_exchange() override;
  void reset_dt() override;
  void *extract(const char *, int &) override;

 private:
  int ninsert, ntype, seed;
  int mode, idnext, dstyle, npoly, rigidflag, shakeflag;
  int ignoreflag, ignoreline, ignoretri;
  double radius_one, radius_max;
  double radius_lo, radius_hi;
  double *radius_poly, *frac_poly;
  double density_lo, density_hi;
  double volfrac;
  int maxattempt;
  int region_style;
  double rate;
  double vxlo, vxhi, vylo, vyhi, vy, vz;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double xc, yc, rc;
  double grav;
  char *idrigid, *idshake;
  char *idregion;
  class Region *region;

  class Molecule **onemols;
  int nmol, natom_max;
  double molradius_max;
  double *molfrac;
  double **coords;
  imageint *imageflags;
  class Fix *fixrigid, *fixshake;
  double oneradius;

  int me, nprocs;
  int *recvcounts, *displs;
  int nfreq, ninserted, nper;
  bigint nfirst;
  double lo_current, hi_current;
  tagint maxtag_all, maxmol_all;
  class RanPark *random, *random2;

  void find_maxid();
  int overlap(int);
  bool outside(int, double, double, double);
  void xyz_random(double, double *);
  double radius_sample();
  void options(int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif
