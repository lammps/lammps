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
#include "kokkos_few.h"

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
  // forward remaining x2lamda() and lambda2x() variants to parent class
  void x2lamda(double *a, double *b) { Domain::x2lamda(a,b); }
  void lamda2x(double *a, double *b) { Domain::lamda2x(a,b); }
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

	// Needed by fix_rigid_kokkos and maybe others
  static KOKKOS_INLINE_FUNCTION
  Few<double, 3> remap(Few<double,3> prd, Few<double,6> h, Few<double, 6> h_inv,
                       int triclinic,
                       Few<double,3> boxlo, Few<double,3> boxhi,
                       Few<double,3> boxlo_lambda, Few<double,3> boxhi_lambda,
                       Few<double,3> prd_lambda,
                       Few<double,3> x, imageint &image, int periodic_bits);

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


KOKKOS_INLINE_FUNCTION
Few<double, 3> DomainKokkos::remap(Few<double,3> prd,
                                   Few<double,6> h, Few<double, 6> h_inv,
                                   int triclinic,
                                   Few<double,3> boxlo, Few<double,3> boxhi,
                                   Few<double,3> boxlo_lambda,
                                   Few<double,3> boxhi_lambda,
                                   Few<double,3> prd_lambda,
                                   Few <double,3> x, imageint &image,
                                   int periodic_bits)
{
  Few<double, 3> lo, hi, period, lambda, coord;

  imageint idim, otherdims;
  // Sorry for making this look like Java. :(
  auto add_few_3 = [](Few<double, 3> a, Few<double, 3> b) {
                     return Few<double, 3>{a[0]+b[0], a[1]+b[1], a[2]+b[2]}; };

  int xperiodic = (periodic_bits & 1);
  int yperiodic = (periodic_bits & 2);
  int zperiodic = (periodic_bits & 4);


  if (triclinic == 0) {
    lo = boxlo;
    hi = add_few_3(boxlo, prd);
    period = prd;
    coord = x;
  } else {
    lo = boxlo_lambda;
    hi = boxhi_lambda;
    period = prd_lambda;

    // verbatim copy of x2lamda code so that this can remain static.
    // Should be equivalent to x2lamda(x,lambda);
    double delta[3];
    delta[0] = x[0] - boxlo[0];
    delta[1] = x[1] - boxlo[1];
    delta[2] = x[2] - boxlo[2];

    lambda[0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
    lambda[1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
    lambda[2] = h_inv[2]*delta[2];

  }

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
    coord[0] = MAX(coord[0], lo[0]);
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

  if (triclinic) {
    // Verbatim copy of lamda2x so that this can remain static
    // This should be equivalent with lamda2x(coord,x)
    x[0] = h[0]*coord[0] + h[5]*coord[1] + h[4]*coord[2] + boxlo[0];
    x[1] = h[1]*coord[1] + h[3]*coord[2] + boxlo[1];
    x[2] = h[2]*coord[2] + boxlo[2];

    return x;
  } else {
    return coord;
  }
}

}

#endif

/* ERROR/WARNING messages:

E: Illegal simulation box

The lower bound of the simulation box is greater than the upper bound.

*/
