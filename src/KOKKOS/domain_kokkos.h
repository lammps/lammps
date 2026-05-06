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

  template <typename T, class T2>
// NOLINTNEXTLINE
  static KOKKOS_INLINE_FUNCTION
  Few<T,3> unmap(Few<T,3>, Few<T,6>, int, const T2&, imageint);

  // NOLINTNEXTLINE
  static KOKKOS_INLINE_FUNCTION
  Few<int,3> image_flags(imageint image);

  template <bool TRICLINIC, typename T>
// NOLINTNEXTLINE
  static KOKKOS_INLINE_FUNCTION
  void minimum_image_big(const Few<bool,3>&, const Few<T,3>&, const Few<T,3>&, const Few<T,6>&, Few<T,3>&, T&);

  using x_t = typename DAT::t_kkfloat_1d_3_lr_randomread;
  using sametag_t = typename DAT::t_int_1d;

// NOLINTNEXTLINE
  static KOKKOS_INLINE_FUNCTION
  int closest_image( x_t x, sametag_t sametag, const int, int);


 private:
  int groupbit;
  KK_FLOAT lo[3],hi[3],period[3];
  int n_flip, m_flip, p_flip;
  ArrayTypes<LMPDeviceType>::t_kkfloat_1d_3_lr x;
  ArrayTypes<LMPDeviceType>::t_imageint_1d image;
  ArrayTypes<LMPDeviceType>::t_int_1d mask;
};


/* ----------------------------------------------------------------------
   static methods
------------------------------------------------------------------------- */

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

/* ----------------------------------------------------------------------
   minimum image convention in periodic dimensions
   use 1/2 of box size as test
   for triclinic, also add/subtract tilt factors in other dims as needed
   allow multiple box lengths to enable distance to
     far-away ghost atom returned by atom->map() to be wrapped back into box
     could be problem for looking up atom IDs when cutoff > boxsize
   this should be used when there is a large image count difference possible
     this applies for example to fix rigid/small

  my friend claude said about cpu version...
  "This ancient code smells like late 90s / early 2000s HPC. Very Sandia.
  Very “don’t trust libm.” 😄
------------------------------------------------------------------------- */

template <bool TRICLINIC, typename T>
KOKKOS_INLINE_FUNCTION
void DomainKokkos::minimum_image_big(
  const Few<bool,3>& periodic,
  const Few<T,3>& prd,
  const Few<T,3>& invprd,
  const Few<T,6>& h,
  Few<T,3>& delta,
  T& dflag)
{
  auto periodic_shift = [](const T d, const T invp, T &l_dflag) -> T {
    const T dfactor = Kokkos::round(d * invp);
    if (Kokkos::abs(dfactor) > static_cast<T>(MAXSMALLINT)) {
      l_dflag = d;
      return static_cast<T>(MAXSMALLINT);
    }
    return dfactor;
  };
  if constexpr (TRICLINIC) {
    if (periodic[2]) {
      const T fd = periodic_shift(delta[2], invprd[2], dflag);
      delta[2] = fma(-prd[2], fd, delta[2]);
      delta[1] = fma(-h[5], fd, delta[1]);
      delta[0] = fma(-h[4], fd, delta[0]);
    }
    if (periodic[1]) {
      const T fd = periodic_shift(delta[1], invprd[1], dflag);
      delta[1] = fma(-prd[1], fd, delta[1]);
      delta[0] = fma(-h[3], fd, delta[0]);
    }
  } else {
    if (periodic[2])
      delta[2] = fma(-prd[2], periodic_shift(delta[2], invprd[2], dflag), delta[2]);
    if (periodic[1])
      delta[1] = fma(-prd[1], periodic_shift(delta[1], invprd[1], dflag), delta[1]);
  }
  if (periodic[0])
    delta[0] = fma(-prd[0], periodic_shift(delta[0], invprd[0], dflag), delta[0]);
}


/* ----------------------------------------------------------------------
   return ix iy iz image flags
   convenience function for printf debugging
------------------------------------------------------------------------- */

// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
Few<int,3> DomainKokkos::image_flags(imageint image) {
  return Few<int,3>(
    (image & IMGMASK) - IMGMAX,
    (image >> IMGBITS & IMGMASK) - IMGMAX,
    (image >> IMG2BITS) - IMGMAX
  );
}

/* ----------------------------------------------------------------------
   return local index of atom J or any of its images that is closest to atom I
   if J is not a valid index like -1, just return it
------------------------------------------------------------------------- */

// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
int DomainKokkos::closest_image(x_t x, sametag_t sametag, const int i, int j)
{
  if (j < 0) return j;
  int closest = j;
  const KK_FLOAT xi0 = x(i,0);
  const KK_FLOAT xi1 = x(i,1);
  const KK_FLOAT xi2 = x(i,2);
  KK_FLOAT delx = xi0 - x(j,0);
  KK_FLOAT dely = xi1 - x(j,1);
  KK_FLOAT delz = xi2 - x(j,2);
  KK_FLOAT rsqmin = delx*delx + dely*dely + delz*delz;
  while (sametag(j) >= 0) {
    j = sametag(j);
    delx = xi0 - x(j,0);
    dely = xi1 - x(j,1);
    delz = xi2 - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < rsqmin) {
      rsqmin = rsq;
      closest = j;
    }
  }
  return closest;
}

}

#endif // !LMP_DOMAIN_KOKKOS_H

