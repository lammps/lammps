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
   Contributing author: Sebastian Hütter (OvGU)
------------------------------------------------------------------------- */

#include "meam.h"

#include "math_special.h"

#include <cmath>

using namespace LAMMPS_NS;

//-----------------------------------------------------------------------------
// Compute G(gamma) based on selection flag ibar:
//   0 => G = sqrt(1+gamma)
//   1 => G = exp(gamma/2)
//   2 => not implemented
//   3 => G = 2/(1+exp(-gamma))
//   4 => G = sqrt(1+gamma)
//  -5 => G = +-sqrt(abs(1+gamma))
//
double MEAM::G_gam(const double gamma, const int ibar, int &errorflag) const
{
  double gsmooth_switchpoint;

  switch (ibar) {
    case 0:
    case 4:
      gsmooth_switchpoint = -gsmooth_factor / (gsmooth_factor + 1);
      if (gamma < gsmooth_switchpoint) {
        //         e.g. gsmooth_factor is 99, {:
        //         gsmooth_switchpoint = -0.99
        //         G = 0.01*(-0.99/gamma)**99
        double G = 1 / (gsmooth_factor + 1) * pow((gsmooth_switchpoint / gamma), gsmooth_factor);
        return sqrt(G);
      } else {
        return sqrt(1.0 + gamma);
      }
    case 1:
      return MathSpecial::fm_exp(gamma / 2.0);
    case 3:
      return 2.0 / (1.0 + MathSpecial::fm_exp(-gamma));
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
//   0 => G = sqrt(1+gamma)
//   1 => G = exp(gamma/2)
//   2 => not implemented
//   3 => G = 2/(1+exp(-gamma))
//   4 => G = sqrt(1+gamma)
//  -5 => G = +-sqrt(abs(1+gamma))
//
double MEAM::dG_gam(const double gamma, const int ibar, double &dG) const
{
  double gsmooth_switchpoint;
  double G;

  switch (ibar) {
    case 0:
    case 4:
      gsmooth_switchpoint = -gsmooth_factor / (gsmooth_factor + 1);
      if (gamma < gsmooth_switchpoint) {
        //         e.g. gsmooth_factor is 99, {:
        //         gsmooth_switchpoint = -0.99
        //         G = 0.01*(-0.99/gamma)**99
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
      G = MathSpecial::fm_exp(gamma / 2.0);
      dG = G / 2.0;
      return G;
    case 3:
      G = 2.0 / (1.0 + MathSpecial::fm_exp(-gamma));
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
double MEAM::zbl(const double r, const int z1, const int z2)
{
  int i;
  const double c[] = {0.028171, 0.28022, 0.50986, 0.18175};
  const double d[] = {0.20162, 0.40290, 0.94229, 3.1998};
  const double azero = 0.4685;
  const double cc = 14.3997;
  double a, x;
  // azero = (9pi^2/128)^1/3 (0.529) Angstroms
  a = azero / (pow(z1, 0.23) + pow(z2, 0.23));
  double result = 0.0;
  x = r / a;
  for (i = 0; i <= 3; i++) { result = result + c[i] * MathSpecial::fm_exp(-d[i] * x); }
  if (r > 0.0) result = result * z1 * z2 / r * cc;
  return result;
}

//-----------------------------------------------------------------------------
// Compute embedding function F(rhobar) and derivative F'(rhobar), eqn I.5
//
double MEAM::embedding(const double A, const double Ec, const double rhobar, double &dF) const
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
      dF = -AEc;
      return -AEc * rhobar;
    }
  }
}

//-----------------------------------------------------------------------------
// Compute Rose energy function, I.16
//
double MEAM::erose(const double r, const double re, const double alpha, const double Ec,
                   const double repuls, const double attrac, const int form)
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
      result = -Ec * (1 + astar + (-attrac + repuls / r) * MathSpecial::cube(astar)) *
          MathSpecial::fm_exp(-astar);
    else if (form == 2)
      result = -Ec * (1 + astar + a3 * MathSpecial::cube(astar)) * MathSpecial::fm_exp(-astar);
    else
      result = -Ec * (1 + astar + a3 * MathSpecial::cube(astar) / (r / re)) *
          MathSpecial::fm_exp(-astar);
  }
  return result;
}

//-----------------------------------------------------------------------------
// Shape factors for various configurations
//
void MEAM::get_shpfcn(const lattice_t latt, const double sthe, const double cthe, double (&s)[3])
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
    case CH4:    // CH4 actually needs shape factor for diamond for C, dimer for H
    case DIA:
    case DIA3:
      s[0] = 0.0;
      s[1] = 0.0;
      s[2] = 32.0 / 9.0;
      break;
    case DIM:
      s[0] = 1.0;
      s[1] = 2.0 / 3.0;
      //        s(3) = 1.d0 // this should be 0.4 unless (1-legendre) is multiplied in the density calc.
      s[2] = 0.40;    // this is (1-legendre) where legendre = 0.6 in dynamo is accounted.
      break;
    case LIN:    //linear, theta being 180
      s[0] = 0.0;
      s[1] = 8.0 / 3.0;    // 4*(co**4 + si**4 - 1.0/3.0) in zig become 4*(1-1/3)
      s[2] = 0.0;
      break;
    case ZIG:    //zig-zag
    case TRI:    //trimer e.g. H2O
      s[0] = 4.0 * pow(cthe, 2);
      s[1] = 4.0 * (pow(cthe, 4) + pow(sthe, 4) - 1.0 / 3.0);
      s[2] = 4.0 * (pow(cthe, 2) * (3 * pow(sthe, 4) + pow(cthe, 4)));
      s[2] = s[2] - 0.6 * s[0];    //legend in dyn, 0.6 is default value.
      break;
    default:
      s[0] = 0.0;
      //        call error('Lattice not defined in get_shpfcn.')
  }
}

//-----------------------------------------------------------------------------
// Number of first neighbors for reference structure
//
int MEAM::get_Zij(const lattice_t latt)
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
    case CH4:    // DYNAMO currently implemented this way while it needs two Z values, 4 and 1
      return 4;
    case LIN:
    case ZIG:
    case TRI:
      return 2;
      //        call error('Lattice not defined in get_Zij.')
  }
  return 0;
}

//-----------------------------------------------------------------------------
// Number of second neighbors for the reference structure
//   a = distance ratio R1/R2 (a2nn in dynamo)
//   numscr = number of atoms that screen the 2NN bond
//   S = second neighbor screening function (xfac, a part of b2nn in dynamo)
//
int MEAM::get_Zij2(const lattice_t latt, const double cmin, const double cmax, const double stheta,
                   double &a, double &S)
{

  double C, x, sijk;
  int Zij2 = 0, numscr = 0;

  switch (latt) {

    case FCC:
      Zij2 = 6;
      a = sqrt(2.0);
      numscr = 4;
      break;

    case BCC:
      Zij2 = 6;
      a = 2.0 / sqrt(3.0);
      numscr = 4;
      break;

    case HCP:
      Zij2 = 6;
      a = sqrt(2.0);
      numscr = 4;
      break;

    case B1:
    case SC:
      Zij2 = 12;
      a = sqrt(2.0);
      numscr = 2;
      break;

    case DIA:    // 2NN
      Zij2 = 12;
      a = sqrt(8.0 / 3.0);
      numscr = 1;
      if (cmin < 0.500001) {
        //          call error('can not do 2NN MEAM for dia')
      }
      break;

    case DIA3:    // 3NN
      Zij2 = 12;
      a = sqrt(11.0 / 3.0);
      numscr = 4;
      if (cmin < 0.500001) {
        //          call error('can not do 2NN MEAM for dia')
      }
      break;

    case CH4:    //does not have 2nn structure so it returns 0
    case LIN:    //line
    case DIM:
      //        this really shouldn't be allowed; make sure screening is zero
      a = 1.0;
      S = 0.0;
      return 0;

    case TRI:    //TRI
      Zij2 = 1;
      a = 2.0 * stheta;
      numscr = 2;
      break;

    case ZIG:    //zig-zag
      Zij2 = 2;
      a = 2.0 * stheta;
      numscr = 1;
      break;

    case L12:
      Zij2 = 6;
      a = sqrt(2.0);
      numscr = 4;
      break;

    case B2:
      Zij2 = 6;
      a = 2.0 / sqrt(3.0);
      numscr = 4;
      break;
    case C11:
      // unsupported lattice flag C11 in get_Zij
      break;
    default:
      // unknown lattic flag in get Zij
      //        call error('Lattice not defined in get_Zij.')
      break;
  }

  // Compute screening for each first neighbor
  if (latt == DIA3) {    // special case for 3NN diamond structure
    C = 1.0;
  } else {
    C = 4.0 / (a * a) - 1.0;
  }
  x = (C - cmin) / (cmax - cmin);
  sijk = fcut(x);
  // There are numscr first neighbors screening the second neighbors
  S = MathSpecial::powint(sijk, numscr);
  return Zij2;
}

int MEAM::get_Zij2_b2nn(const lattice_t latt, const double cmin, const double cmax, double &S)
{

  double x, sijk, C;
  int numscr = 0, Zij2 = 0;
  switch (latt) {
    case ZIG:    //zig-zag for b11s and b22s term
    case TRI:    //trimer for b11s
      Zij2 = 2;
      numscr = 1;
      break;
    default:
      // unknown lattic flag in get Zij
      //        call error('Lattice not defined in get_Zij.')
      break;
  }
  C = 1.0;
  x = (C - cmin) / (cmax - cmin);
  sijk = fcut(x);
  S = MathSpecial::powint(sijk, numscr);
  return Zij2;
}
