// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_amoeba.h"

#include "atom.h"
#include "error.h"
#include "math_special.h"
#include "neigh_list.h"

#include <cmath>

using namespace LAMMPS_NS;

using MathSpecial::square;
using MathSpecial::cube;
using MathSpecial::powint;

enum{VDWL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};

/* ----------------------------------------------------------------------
   hal = buffered 14-7 Vdwl interactions
   adapted from Tinker ehal1c() routine
------------------------------------------------------------------------- */

void PairAmoeba::hal()
{
  int i,j,ii,jj,itype,jtype,iclass,jclass,iv,jv;
  int special_which;
  double e,de,eps;
  double rv,rv7;
  double xi,yi,zi;
  double xr,yr,zr;
  double redi,rediv;
  double redj,redjv;
  double dedx,dedy,dedz;
  double rho,tau,tau7;
  double dtau,gtau;
  double taper,dtaper;
  double rik,rik2,rik3;
  double rik4,rik5;
  double rik6,rik7;
  double vxx,vyy,vzz;
  double vyx,vzx,vzy;
  double factor_hal;

  int inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // set cutoffs and taper coeffs

  choose(VDWL);

  // owned atoms

  double **f = atom->f;

  // neigh list

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // find van der Waals energy and derivatives via neighbor list

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = amtype[i];
    iclass = amtype2class[itype];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    redi = kred[iclass];
    rediv = 1.0 - redi;
    xi = xred[i][0];
    yi = xred[i][1];
    zi = xred[i][2];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      special_which = sbmask15(j);
      factor_hal = special_hal[special_which];
      if (factor_hal == 0.0) continue;
      j &= NEIGHMASK15;

      xr = xi - xred[j][0];
      yr = yi - xred[j][1];
      zr = zi - xred[j][2];
      rik2 = xr*xr + yr*yr + zr*zr;

      if (rik2 > off2) continue;

      // compute the energy contribution for this interaction

      jtype = amtype[j];
      jclass = amtype2class[jtype];

      // check for an interaction distance less than the cutoff
      // special_which = 3 is a 1-4 neighbor with its own sigma,epsilon

      rik = sqrt(rik2);
      rv = radmin[iclass][jclass];
      eps = epsilon[iclass][jclass];
      if (special_which == 3) {
        rv = radmin4[iclass][jclass];
        eps = epsilon4[iclass][jclass];
      }
      eps *= factor_hal;

      rv7 = powint(rv,7);
      rik6 = cube(rik2);
      rik7 = rik6 * rik;
      rho = rik7 + ghal*rv7;
      tau = (dhal+1.0) / (rik + dhal*rv);
      tau7 = powint(tau,7);
      dtau = tau / (dhal+1.0);
      gtau = eps*tau7*rik6*(ghal+1.0)*square(rv7/rho);
      e = eps*tau7*rv7*((ghal+1.0)*rv7/rho-2.0);
      de = -7.0 * (dtau*e+gtau);

      // use energy switching if near the cutoff distance

      if (rik2 > cut2) {
        rik3 = rik2 * rik;
        rik4 = rik2 * rik2;
        rik5 = rik2 * rik3;
        taper = c5*rik5 + c4*rik4 + c3*rik3 + c2*rik2 + c1*rik + c0;
        dtaper = 5.0*c5*rik4 + 4.0*c4*rik3 + 3.0*c3*rik2 + 2.0*c2*rik + c1;
        de = e*dtaper + de*taper;
        e *= taper;
      }

      ehal += e;

      // find the chain rule terms for derivative components

      de = de / rik;
      dedx = de * xr;
      dedy = de * yr;
      dedz = de * zr;

      // increment the total van der Waals energy and derivatives
      // if jv < 0, trigger an error, needed H-bond partner is missing

      iv = red2local[i];
      jv = red2local[j];
      if (jv < 0)
        error->one(FLERR,"AMOEBA hal cannot find H bond partner - "
                   "ghost comm is too short");

      if (i == iv) {
        f[i][0] -= dedx;
        f[i][1] -= dedy;
        f[i][2] -= dedz;
      } else {
        f[i][0] -= dedx*redi;
        f[i][1] -= dedy*redi;
        f[i][2] -= dedz*redi;
        f[iv][0] -= dedx*rediv;
        f[iv][1] -= dedy*rediv;
        f[iv][2] -= dedz*rediv;
      }

      if (j == jv) {
        f[j][0] += dedx;
        f[j][1] += dedy;
        f[j][2] += dedz;
      } else {
        redj = kred[jclass];
        redjv = 1.0 - redj;
        f[j][0] += dedx*redj;
        f[j][1] += dedy*redj;
        f[j][2] += dedz*redj;
        f[jv][0] += dedx*redjv;
        f[jv][1] += dedy*redjv;
        f[jv][2] += dedz*redjv;
      }

      // increment the internal virial tensor components

      if (vflag_global) {
        vxx = xr * dedx;
        vyx = yr * dedx;
        vzx = zr * dedx;
        vyy = yr * dedy;
        vzy = zr * dedy;
        vzz = zr * dedz;

        virhal[0] -= vxx;
        virhal[1] -= vyy;
        virhal[2] -= vzz;
        virhal[3] -= vyx;
        virhal[4] -= vzx;
        virhal[5] -= vzy;
      }
    }
  }
}
