// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Naga Vydyanathan (NVIDIA)
------------------------------------------------------------------------- */

#include "math_special_kokkos.h"

#include <cmath>

#include "meam_kokkos.h"
using namespace MathSpecialKokkos;

//-----------------------------------------------------------------------------
// Compute G(gamma) based on selection flag ibar:
//  0 => G = sqrt(1+gamma)
//  1 => G = exp(gamma/2)
//  2 => not implemented
//  3 => G = 2/(1+exp(-gamma))
//  4 => G = sqrt(1+gamma)
// -5 => G = +-sqrt(abs(1+gamma))
//
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double MEAMKokkos<DeviceType>::G_gam(const double gamma, const int ibar, int &errorflag) const
{
  double gsmooth_switchpoint;

  switch (ibar) {
    case 0:
    case 4:
      gsmooth_switchpoint = -gsmooth_factor / (gsmooth_factor + 1);
      if (gamma < gsmooth_switchpoint) {
        // e.g. gsmooth_factor is 99, {:
        // gsmooth_switchpoint = -0.99
        // G = 0.01*(-0.99/gamma)**99
        double G = 1 / (gsmooth_factor + 1) * pow((gsmooth_switchpoint / gamma), gsmooth_factor);
        return sqrt(G);
      } else {
        return sqrt(1.0 + gamma);
      }
    case 1:
      return MathSpecialKokkos::fm_exp(gamma / 2.0);
    case 3:
      return 2.0 / (1.0 + MathSpecialKokkos::fm_exp(-gamma));
    case -5:
      if ((1.0 + gamma) >= 0) {
        return sqrt(1.0 + gamma);
      } else {
        return -sqrt(-1.0 - gamma);
      }
  }
  errorflag = 1;
  return 0.0;
}

//-----------------------------------------------------------------------------
// Compute G(gamma and dG(gamma) based on selection flag ibar:
//  0 => G = sqrt(1+gamma)
//  1 => G = exp(gamma/2)
//  2 => not implemented
//  3 => G = 2/(1+exp(-gamma))
//  4 => G = sqrt(1+gamma)
// -5 => G = +-sqrt(abs(1+gamma))
//
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double MEAMKokkos<DeviceType>::dG_gam(const double gamma, const int ibar, double& dG) const
{
  double gsmooth_switchpoint;
  double G;

  switch (ibar) {
    case 0:
    case 4:
      gsmooth_switchpoint = -gsmooth_factor / (gsmooth_factor + 1);
      if (gamma < gsmooth_switchpoint) {
        // e.g. gsmooth_factor is 99, {:
        // gsmooth_switchpoint = -0.99
        // G = 0.01*(-0.99/gamma)**99
        G = 1 / (gsmooth_factor + 1) * pow((gsmooth_switchpoint / gamma), gsmooth_factor);
        G = sqrt(G);
        dG = -gsmooth_factor * G / (2.0 * gamma);
        return G;
      } else {
        G = sqrt(1.0 + gamma);
        dG = 1.0 / (2.0 * G);
        return G;
      }
    case 1:
      G = MathSpecialKokkos::fm_exp(gamma / 2.0);
      dG = G / 2.0;
      return G;
    case 3:
      G = 2.0 / (1.0 + MathSpecialKokkos::fm_exp(-gamma));
      dG = G * (2.0 - G) / 2;
      return G;
    case -5:
      if ((1.0 + gamma) >= 0) {
        G = sqrt(1.0 + gamma);
        dG = 1.0 / (2.0 * G);
        return G;
      } else {
        G = -sqrt(-1.0 - gamma);
        dG = -1.0 / (2.0 * G);
        return G;
      }
  }
  dG = 1.0;
  return 0.0;
}

//-----------------------------------------------------------------------------
// Compute ZBL potential
//
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double MEAMKokkos<DeviceType>::zbl(const double r, const int z1, const int z2) const
{
  int i;
  const double c[] = { 0.028171, 0.28022, 0.50986, 0.18175 };
  const double d[] = { 0.20162, 0.40290, 0.94229, 3.1998 };
  const double azero = 0.4685;
  const double cc = 14.3997;
  double a, x;
  // azero = (9pi^2/128)^1/3 (0.529) Angstroms
  a = azero / (pow(z1, 0.23) + pow(z2, 0.23));
  double result = 0.0;
  x = r / a;
  for (i = 0; i <= 3; i++) {
    result = result + c[i] * MathSpecialKokkos::fm_exp(-d[i] * x);
  }
  if (r > 0.0)
    result = result * z1 * z2 / r * cc;
  return result;
}

//-----------------------------------------------------------------------------
// Compute embedding function F(rhobar) and derivative F'(rhobar), eqn I.5
//
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double MEAMKokkos<DeviceType>::embedding(const double A, const double Ec, const double rhobar, double& dF) const
{
  const double AEc = A * Ec;

  if (rhobar > 0.0) {
      const double lrb = log(rhobar);
      dF = AEc * (1.0 + lrb);
      return AEc * rhobar * lrb;
  } else {
    if (emb_lin_neg == 0) {
      dF = 0.0;
      return 0.0;
    } else {
      dF = - AEc;
      return - AEc * rhobar;
    }
  }
}

//-----------------------------------------------------------------------------
// Compute Rose energy function, I.16
//
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double MEAMKokkos<DeviceType>::erose(const double r, const double re, const double alpha, const double Ec, const double repuls,
            const double attrac, const int form) const
{
  double astar, a3;
  double result = 0.0;

  if (r > 0.0) {
    astar = alpha * (r / re - 1.0);
    a3 = 0.0;
    if (astar >= 0)
      a3 = attrac;
    else if (astar < 0)
      a3 = repuls;

    if (form == 1)
      result = -Ec * (1 + astar + (-attrac + repuls / r) * MathSpecialKokkos::cube(astar)) * MathSpecialKokkos::fm_exp(-astar);
    else if (form == 2)
      result = -Ec * (1 + astar + a3 * MathSpecialKokkos::cube(astar)) * MathSpecialKokkos::fm_exp(-astar);
    else
      result = -Ec * (1 + astar + a3 * MathSpecialKokkos::cube(astar) / (r / re)) * MathSpecialKokkos::fm_exp(-astar);
  }
  return result;
}

//-----------------------------------------------------------------------------
// Shape factors for various configurations
//
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void MEAMKokkos<DeviceType>::get_shpfcn(const lattice_t latt, const double sthe, const double cthe, double (&s)[3]) const
{
  switch (latt) {
    case FCC:
    case BCC:
    case B1:
    case B2:
    case SC:
      s[0] = 0.0;
      s[1] = 0.0;
      s[2] = 0.0;
      break;
    case HCP:
      s[0] = 0.0;
      s[1] = 0.0;
      s[2] = 1.0 / 3.0;
      break;
    case CH4: // CH4 actually needs shape factor for diamond for C, dimer for H
    case DIA:
    case DIA3:
      s[0] = 0.0;
      s[1] = 0.0;
      s[2] = 32.0 / 9.0;
      break;
    case DIM:
      s[0] = 1.0;
      s[1] = 2.0 / 3.0;
      // s(4) = 1.d0 // this should be 0.4 unless (1-legendre) is multiplied in the density calc.
      s[2] = 0.40; // this is (1-legendre) where legendre = 0.6 in dynamo is accounted.
      break;
    case LIN: // linear, theta being 180
      s[0] = 0.0;
      s[1] = 8.0 / 3.0; // 4*(co**4 + si**4 - 1.0/3.0) in zig become 4*(1-1/3)
      s[2] = 0.0;
      break;
    case ZIG: //zig-zag
    case TRI: //trimer e.g. H2O
      s[0] = 4.0*pow(cthe,2);
      s[1] = 4.0*(pow(cthe,4) + pow(sthe,4) - 1.0/3.0);
      s[2] = 4.0*(pow(cthe,2) * (3*pow(sthe,4) + pow(cthe,4)));
      s[2] = s[2] - 0.6*s[0]; //legend in dyn, 0.6 is default value.
      break;
    default:
      s[0] = 0.0;
      // call error('Lattice not defined in get_shpfcn.')
  }
}

//-----------------------------------------------------------------------------
// Number of neighbors for the reference structure
//
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int MEAMKokkos<DeviceType>::get_Zij(const lattice_t latt) const
{
  switch (latt) {
    case FCC:
      return 12;
    case BCC:
      return 8;
    case HCP:
      return 12;
    case DIA:
    case DIA3:
      return 4;
    case DIM:
      return 1;
    case B1:
    case SC:
      return 6;
    case C11:
      return 10;
    case L12:
      return 12;
    case B2:
      return 8;
    case CH4: // DYNAMO currently implemented this way while it needs two Z values, 4 and 1
      return 4;
    case LIN:
    case ZIG:
    case TRI:
      return 2;
      // call error('Lattice not defined in get_Zij.')
  }
  return 0;
}
