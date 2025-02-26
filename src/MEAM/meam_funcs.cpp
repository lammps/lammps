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
   Contributing author: Sebastian Huetter (OvGU)
------------------------------------------------------------------------- */

#include "meam.h"

#include "math_special.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MEAM_NS;


/* ----------------------------------------------------------------------
   Compute G(gamma) based on selection flag ibar:
     0 => G = sqrt(1+gamma)
     1 => G = exp(gamma/2)
     2 => not implemented
     3 => G = 2/(1+exp(-gamma))
     4 => G = sqrt(1+gamma)
    -5 => G = +-sqrt(abs(1+gamma))
------------------------------------------------------------------------- */

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
        if (gamma == 0.0) return 0.0; // avoid division by zero. For gamma = 0.0 => G = 1 / inf
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

/* ----------------------------------------------------------------------
   Compute G(gamma and dG(gamma) based on selection flag ibar:
     0 => G = sqrt(1+gamma)
     1 => G = exp(gamma/2)
     2 => not implemented
     3 => G = 2/(1+exp(-gamma))
     4 => G = sqrt(1+gamma)
    -5 => G = +-sqrt(abs(1+gamma))

------------------------------------------------------------------------- */

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

/* ----------------------------------------------------------------------
   Compute embedding function F(rhobar) and derivative F'(rhobar), eqn I.5
------------------------------------------------------------------------- */

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

/* ----------------------------------------------------------------------
   Compute ZBL potential
------------------------------------------------------------------------- */

double MEAM_NS::zbl(const double r, const int z1, const int z2)
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

/* ----------------------------------------------------------------------
   Compute Rose energy function, I.16
------------------------------------------------------------------------- */

double MEAM_NS::erose(const double r, const double re, const double alpha, const double Ec,
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
