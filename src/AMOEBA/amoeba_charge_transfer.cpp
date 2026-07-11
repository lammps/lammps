// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_amoeba.h"

#include "atom.h"
#include "neigh_list.h"

#include <cmath>

using namespace LAMMPS_NS;

enum{VDWL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};

/* ----------------------------------------------------------------------
   charge_transfer = HIPPO charge transfer forces
   adapted from Tinker echgtrn1b() routine
------------------------------------------------------------------------- */

void PairAmoeba::charge_transfer()
{
  int i,j,ii,jj,itype,jtype,iclass,jclass;
  double e,de,felec;
  double rr1,r,r2;
  double r3,r4,r5;
  double xi,yi,zi;
  double xr,yr,zr;
  double chgi,chgj;
  double alphai,alphaj;
  double expi,expj;
  double frcx,frcy,frcz;
  double vxx,vyy,vzz;
  double vxy,vxz,vyz;
  double taper,dtaper;
  double factor_mpole;

  int inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // set cutoffs and taper coeffs

  choose(QFER);

  // owned atoms

  double **x = atom->x;
  double **f = atom->f;

  // neigh list

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // set the energy unit conversion factor

  felec = electric / am_dielectric;

  // find charge transfer energy and derivatives via neighbor list

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = amtype[i];
    iclass = amtype2class[itype];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    xi = x[i][0];
    yi = x[i][1];
    zi = x[i][2];
    chgi = chgct[iclass];
    alphai = dmpct[iclass];
    if (alphai == 0.0) alphai = 100.0;

    // evaluate all sites within the cutoff distance

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_mpole = special_mpole[sbmask15(j)];
      if (factor_mpole == 0.0) continue;
      j &= NEIGHMASK15;

      xr = x[j][0] - xi;
      yr = x[j][1] - yi;
      zr = x[j][2] - zi;
      r2 = xr*xr + yr* yr + zr*zr;
      if (r2 > off2) continue;

      jtype = amtype[j];
      jclass = amtype2class[jtype];

      r = sqrt(r2);
      rr1 = 1.0 / r;
      chgj = chgct[jclass];
      alphaj = dmpct[jclass];
      if (alphaj == 0.0) alphaj = 100.0;

      expi = exp(-alphai*r);
      expj = exp(-alphaj*r);
      e = -chgi*expj - chgj*expi;
      de = chgi*expj*alphaj + chgj*expi*alphai;
      e = felec * e * factor_mpole;
      de = felec * de * factor_mpole;

      // use energy switching if near the cutoff distance

      if (r2 > cut2) {
        r3 = r2 * r;
        r4 = r2 * r2;
        r5 = r2 * r3;
        taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0;
        dtaper = 5.0*c5*r4 + 4.0*c4*r3 + 3.0*c3*r2 + 2.0*c2*r + c1;
        de = e*dtaper + de*taper;
        e *= taper;
      }

      eqxfer += e;

      // compute the force components for this interaction

      frcx = de * xr * rr1;
      frcy = de * yr * rr1;
      frcz = de * zr * rr1;

      // increment the total charge transfer energy and derivatives

      f[i][0] += frcx;
      f[i][1] += frcy;
      f[i][2] += frcz;
      f[j][0] -= frcx;
      f[j][1] -= frcy;
      f[j][2] -= frcz;

      // increment the internal virial tensor components

      if (vflag_global) {
        vxx = xr * frcx;
        vxy = yr * frcx;
        vxz = zr * frcx;
        vyy = yr * frcy;
        vyz = zr * frcy;
        vzz = zr * frcz;

        virqxfer[0] -= vxx;
        virqxfer[1] -= vyy;
        virqxfer[2] -= vzz;
        virqxfer[3] -= vxy;
        virqxfer[4] -= vxz;
        virqxfer[5] -= vyz;
      }
    }
  }
}

