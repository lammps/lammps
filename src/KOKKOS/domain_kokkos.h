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

#ifndef LMP_DOMAIN_KOKKOS_H
#define LMP_DOMAIN_KOKKOS_H

#include "domain.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagDomain_remap_all{};
struct TagDomain_image_flip{};
struct TagDomain_lamda2x{};
struct TagDomain_x2lamda{};

class DomainKokkos : public Domain {
 public:
  DomainKokkos(class LAMMPS *);
  ~DomainKokkos() {}
  void init();
  void reset_box();
  void pbc();
  void remap_all();
  void image_flip(int, int, int);
  void x2lamda(int);
  void lamda2x(int);

  int closest_image(const int, int) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_remap_all, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_image_flip, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_lamda2x, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_x2lamda, const int&) const;

 private:
  double lo[3],hi[3],period[3];
  int n_flip, m_flip, p_flip;
  ArrayTypes<LMPDeviceType>::t_x_array x;
  ArrayTypes<LMPDeviceType>::t_imageint_1d image;
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal simulation box

The lower bound of the simulation box is greater than the upper bound.

*/
