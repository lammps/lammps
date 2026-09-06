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
struct TagDomain_lamda2x_group{};
struct TagDomain_x2lamda{};
struct TagDomain_x2lamda_group{};

class DomainKokkos : public Domain {
 public:
  DomainKokkos(class LAMMPS *);
  ~DomainKokkos() override = default;
  void reset_box() override;
  void pbc() override;
  void remap_all() override;
  void image_flip(int, int, int) override;
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

  int detached_atom_x() const;

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
  static KOKKOS_INLINE_FUNCTION
  Few<double,3> unmap(Few<double,3> prd, Few<double,6> h, int triclinic,
      Few<double,3> x, imageint image);
  static KOKKOS_INLINE_FUNCTION
  Few<double,3> lamda2x(Few<double,3> boxlo, Few<double,6> h,
      Few<double,3> lamda);
  static KOKKOS_INLINE_FUNCTION
  Few<double,3> x2lamda(Few<double,3> boxlo, Few<double,6> h_inv,
      Few<double,3> x);
  static KOKKOS_INLINE_FUNCTION
  Few<double,3> remap(Few<double,3> lo, Few<double,3> hi, Few<double,3> period, int xperiodic, int yperiodic, int zperiodic,
      Few<double,3> coord, imageint &image);

 private:
  int groupbit;
  KK_FLOAT lo[3],hi[3],period[3];
  int n_flip, m_flip, p_flip;
  ArrayTypes<LMPDeviceType>::t_kkfloat_1d_3_lr x;
  ArrayTypes<LMPDeviceType>::t_imageint_1d image;
  ArrayTypes<LMPDeviceType>::t_int_1d mask;
};

// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
Few<double,3> DomainKokkos::lamda2x(Few<double,3> boxlo, Few<double,6> h,
    Few<double,3> lamda)
{
  Few<double,3> x;
  x[0] = h[0]*lamda[0] + h[5]*lamda[1] + h[4]*lamda[2] + boxlo[0];
  x[1] = h[1]*lamda[1] + h[3]*lamda[2] + boxlo[1];
  x[2] = h[2]*lamda[2] + boxlo[2];
  return x;
}

KOKKOS_INLINE_FUNCTION
Few<double,3> DomainKokkos::x2lamda(Few<double,3> boxlo, Few<double,6> h_inv,
    Few<double,3> x)
{
  double delta[3];
  delta[0] = x[0] - boxlo[0];
  delta[1] = x[1] - boxlo[1];
  delta[2] = x[2] - boxlo[2];

  Few<double,3> lamda;
  lamda[0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
  lamda[1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
  lamda[2] = h_inv[2]*delta[2];
  return lamda;
}

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

/*
 * Should be called with x for non-triclinic and lamda for triclinix
 */
KOKKOS_INLINE_FUNCTION
Few<double,3> DomainKokkos::remap(Few<double,3> lo, Few<double,3> hi, Few<double,3> period, int xperiodic, int yperiodic, int zperiodic,
    Few<double,3> xorlamda, imageint &image)
{
  Few<double,3> coord(xorlamda);

  imageint idim, otherdims;
  if (xperiodic) {
    while (coord[0] < lo[0]) {
      coord[0] += period[0];
      idim = image & IMGMASK;
      otherdims = image ^ idim;
      idim--;
      idim &= IMGMASK;
      image = otherdims | idim;
    }
    while (coord[0] >= hi[0]) {
      coord[0] -= period[0];
      idim = image & IMGMASK;
      otherdims = image ^ idim;
      idim++;
      idim &= IMGMASK;
      image = otherdims | idim;
    }
    coord[0] = MAX(coord[0],lo[0]);
  }

  if (yperiodic) {
    while (coord[1] < lo[1]) {
      coord[1] += period[1];
      idim = (image >> IMGBITS) & IMGMASK;
      otherdims = image ^ (idim << IMGBITS);
      idim--;
      idim &= IMGMASK;
      image = otherdims | (idim << IMGBITS);
    }
    while (coord[1] >= hi[1]) {
      coord[1] -= period[1];
      idim = (image >> IMGBITS) & IMGMASK;
      otherdims = image ^ (idim << IMGBITS);
      idim++;
      idim &= IMGMASK;
      image = otherdims | (idim << IMGBITS);
    }
    coord[1] = MAX(coord[1],lo[1]);
  }

  if (zperiodic) {
    while (coord[2] < lo[2]) {
      coord[2] += period[2];
      idim = image >> IMG2BITS;
      otherdims = image ^ (idim << IMG2BITS);
      idim--;
      idim &= IMGMASK;
      image = otherdims | (idim << IMG2BITS);
    }
    while (coord[2] >= hi[2]) {
      coord[2] -= period[2];
      idim = image >> IMG2BITS;
      otherdims = image ^ (idim << IMG2BITS);
      idim++;
      idim &= IMGMASK;
      image = otherdims | (idim << IMG2BITS);
    }
    coord[2] = MAX(coord[2],lo[2]);
  }

  return coord;
}

}

#endif

