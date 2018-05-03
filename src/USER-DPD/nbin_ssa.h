/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NBIN_CLASS

NBinStyle(ssa,
          NBinSSA,
          NB_SSA)

#else

#ifndef LMP_NBIN_SSA_H
#define LMP_NBIN_SSA_H

#include <cmath>
#include "nbin_standard.h"
#include "error.h"

namespace LAMMPS_NS {

class NBinSSA : public NBinStandard {
 public:

  int gairhead_ssa[8];       // index of 1st ghost atom in each AIR

  // Bounds of the local atoms in the binhead array
  int lbinxlo;               // lowest local bin x-dim coordinate
  int lbinylo;               // lowest local bin y-dim coordinate
  int lbinzlo;               // lowest local bin z-dim coordinate
  int lbinxhi;               // highest local bin x-dim coordinate
  int lbinyhi;               // highest local bin y-dim coordinate
  int lbinzhi;               // highest local bin z-dim coordinate

  NBinSSA(class LAMMPS *);
  ~NBinSSA();

  void bin_atoms_setup(int);
  void bin_atoms();

  bigint memory_usage();

  inline
  int coord2bin(const double & x,const double & y,const double & z, int &ixo, int &iyo, int &izo) const
  {
    int ix,iy,iz;

    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
      error->one(FLERR,"Non-numeric positions - simulation unstable");

    if (x >= bboxhi_[0])
      ix = static_cast<int> ((x-bboxhi_[0])*bininvx) + nbinx;
    else if (x >= bboxlo_[0]) {
      ix = static_cast<int> ((x-bboxlo_[0])*bininvx);
      ix = MIN(ix,nbinx-1);
    } else
      ix = static_cast<int> ((x-bboxlo_[0])*bininvx) - 1;

    if (y >= bboxhi_[1])
      iy = static_cast<int> ((y-bboxhi_[1])*bininvy) + nbiny;
    else if (y >= bboxlo_[1]) {
      iy = static_cast<int> ((y-bboxlo_[1])*bininvy);
      iy = MIN(iy,nbiny-1);
    } else
      iy = static_cast<int> ((y-bboxlo_[1])*bininvy) - 1;

    if (z >= bboxhi_[2])
      iz = static_cast<int> ((z-bboxhi_[2])*bininvz) + nbinz;
    else if (z >= bboxlo_[2]) {
      iz = static_cast<int> ((z-bboxlo_[2])*bininvz);
      iz = MIN(iz,nbinz-1);
    } else
      iz = static_cast<int> ((z-bboxlo_[2])*bininvz) - 1;

    ixo = ix - mbinxlo;
    iyo = iy - mbinylo;
    izo = iz - mbinzlo;

    return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
  }

 private:
  int coord2ssaAIR(const double *);  // map atom coord to an AIR number
  double bboxlo_[3],bboxhi_[3];

};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
