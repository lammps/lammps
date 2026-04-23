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

using Kokkos::fma;

struct TagDomain_remap_all{};
struct TagDomain_image_flip{};
struct TagDomain_lamda2x{};
struct TagDomain_lamda2x_group{};
struct TagDomain_x2lamda{};
struct TagDomain_x2lamda_group{};

class DomainKokkos : public Domain {
 public:
  DomainKokkos(class LAMMPS *);
  ~DomainKokkos() override = default;
  void reset_box() override;
  void pbc() override;
  void remap_all();
  void image_flip(int, int, int);
  void x2lamda(int) override;
  void x2lamda(int,int) override;
  void lamda2x(int) override;
  void lamda2x(int,int) override;
  // forward remaining x2lamda() and lambda2x() variants to parent class
  void x2lamda(double *a, double *b) override { Domain::x2lamda(a,b); }
  void lamda2x(double *a, double *b) override { Domain::lamda2x(a,b); }
  void x2lamda(double *a, double *b, double *c, double *d) {
    Domain::x2lamda(a,b,c,d);
  }

  int closest_image(const int, int) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_remap_all, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_image_flip, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_lamda2x, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_lamda2x_group, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_x2lamda, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDomain_x2lamda_group, const int&) const;

// NOLINTNEXTLINE
  template <typename T>
  static KOKKOS_INLINE_FUNCTION
  Few<T,3> unmap(Few<T,3> prd, Few<T,6> h, int triclinic, Few<T,3> x, imageint image);

 private:
  int groupbit;
  KK_FLOAT lo[3],hi[3],period[3];
  int n_flip, m_flip, p_flip;
  ArrayTypes<LMPDeviceType>::t_kkfloat_1d_3_lr x;
  ArrayTypes<LMPDeviceType>::t_imageint_1d image;
  ArrayTypes<LMPDeviceType>::t_int_1d mask;
};

// T2 can be Few<T,3> or anything with [] operator for convenience, eg. float*
// BEFORE
//  Few<KK_FLOAT,3> x_i;
//  x_i[0] = l_x(i,0); x_i[1] = l_x(i,1); x_i[2] = l_x(i,2);
//  Few<KK_FLOAT,3> unwrap = DomainKokkos::unmap(l_prd, l_h, l_triclinic, x_i, l_image(i));
// AFTER
//  auto unwrap = DomainKokkos::unmap(l_prd, l_h, l_triclinic, &l_x(i,0), l_image(i));

// NOLINTNEXTLINE
template <typename T, class T2>
KOKKOS_INLINE_FUNCTION
Few<T,3> DomainKokkos::unmap(Few<T,3> prd, Few<T,6> h, int triclinic,
                             const T2 &x, imageint image)
{
  const T xbox = static_cast<T>((image & IMGMASK) - IMGMAX);
  const T ybox = static_cast<T>((image >> IMGBITS & IMGMASK) - IMGMAX);
  const T zbox = static_cast<T>((image >> IMG2BITS) - IMGMAX);
  Few<T,3> y;
  if (triclinic == 0) {
    y[0] = fma(xbox, prd[0], x[0]);
    y[1] = fma(ybox, prd[1], x[1]);
    y[2] = fma(zbox, prd[2], x[2]);
  } else {
    // x[0] + h[0]*xbox + h[5]*ybox + h[4]*zbox
    y[0] = fma(h[5], ybox, fma(h[4], zbox, fma(h[0], xbox, x[0])));
    // x[1] + h[1]*ybox + h[3]*zbox
    y[1] = fma(h[3], zbox, fma(h[1], ybox, x[1]));
    // h[2]*zbox
    y[2] = fma(h[2], zbox, x[2]);
  }
  return y;
}

}

#endif // !LMP_DOMAIN_KOKKOS_H

