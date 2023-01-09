// clang-format off
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

#ifndef LMP_DOMAIN_KOKKOS_H
#define LMP_DOMAIN_KOKKOS_H

#include "domain.h"             // IWYU pragma: export
#include "kokkos_type.h"
#include "kokkos_few.h"

namespace LAMMPS_NS {

struct TagDomain_remap_all{};
struct TagDomain_image_flip{};
struct TagDomain_lamda2x{};
struct TagDomain_x2lamda{};

class DomainKokkos : public Domain {
 public:
  DomainKokkos(class LAMMPS *);
  ~DomainKokkos() override = default;
  void reset_box() override;
  void pbc() override;
  void remap_all();
  void image_flip(int, int, int);
  void x2lamda(int) override;
  void lamda2x(int) override;
  // forward remaining x2lamda() and lambda2x() variants to parent class
  void x2lamda(double *a, double *b) override { Domain::x2lamda(a,b); }
  void lamda2x(double *a, double *b) override { Domain::lamda2x(a,b); }
  void x2lamda(double *a, double *b, double *c, double *d) {
    Domain::x2lamda(a,b,c,d);
  }

  int closest_image(const int, int) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_remap_all, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_image_flip, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_lamda2x, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_x2lamda, const int&) const;

  static KOKKOS_INLINE_FUNCTION
  Few<double,3> unmap(Few<double,3> prd, Few<double,6> h, int triclinic,
      Few<double,3> x, imageint image);

 private:
  double lo[3],hi[3],period[3];
  int n_flip, m_flip, p_flip;
  ArrayTypes<LMPDeviceType>::t_x_array x;
  ArrayTypes<LMPDeviceType>::t_imageint_1d image;
};

KOKKOS_INLINE_FUNCTION
Few<double,3> DomainKokkos::unmap(Few<double,3> prd, Few<double,6> h,
    int triclinic, Few<double,3> x, imageint image)
{
  int xbox = (image & IMGMASK) - IMGMAX;
  int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
  int zbox = (image >> IMG2BITS) - IMGMAX;
  Few<double,3> y;
  if (triclinic == 0) {
    y[0] = x[0] + xbox*prd[0];
    y[1] = x[1] + ybox*prd[1];
    y[2] = x[2] + zbox*prd[2];
  } else {
    y[0] = x[0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
    y[1] = x[1] + h[1]*ybox + h[3]*zbox;
    y[2] = x[2] + h[2]*zbox;
  }
  return y;
}

}

#endif

