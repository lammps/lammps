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
#include "domain.h"
#include "error.h"
#include "fix_store_peratom.h"
#include "neigh_list.h"

#include <cmath>

using namespace LAMMPS_NS;

enum{NOFRAME,ZONLY,ZTHENX,BISECTOR,ZBISECT,THREEFOLD};

// ----------------------------------------------------------------------
// utility methods used by different parts of the AMOEBA force field
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   kmpole performs one-time assignment of
     xyzaxis multipole neighbors to each owned atom
     any of the values can be 0 if not used
     yaxis can later be negative due to chkpole()
   also sets polaxe and pole[13] multipole for each owned atom
------------------------------------------------------------------------- */

void PairAmoeba::kmpole()
{
  bool path;
  int i,j,k,m,j12,k12,m12,k13,m13,flag;
  int iframe,nframe;
  int itype,jtype,ktype,mtype,xtype,ytype,ztype;
  tagint jneigh,kneigh,mneigh;

  // DEBUG vectors

  tagint bondneigh[12];
  tagint angleneigh[36];

  amtype = atom->ivector[index_amtype];
  int *polaxe = atom->ivector[index_polaxe];
  double **xyzaxis = atom->darray[index_xyzaxis];
  double **pole = fixpole->astore;

  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int nlocal = atom->nlocal;

  int nmissing = 0;

  for (i = 0; i < nlocal; i++) {
    itype = amtype[i];
    nframe = nmultiframe[itype];

    // flag is used to prevent matching multiple times
    // only first match is used

    flag = 0;

    // create a sorted version of bond/angle neighs from special[][]
    // NOTE: this is to try and do it identically to Tinker
    //   b/c I think in Tinker, which case is seen first can depend on atom order

    for (j = 0; j < nspecial[i][0]; j++)
      bondneigh[j] = special[i][j];
    for (m = 0; m < nspecial[i][0]; m++) {
      tagint smallest = MAXTAGINT;
      for (j = m; j < nspecial[i][0]; j++) {
        if (bondneigh[j] < smallest) {
          smallest = bondneigh[j];
          k = j;
          bondneigh[k] = bondneigh[m];
          bondneigh[m] = smallest;
        }
      }
    }

    for (j = nspecial[i][0]; j < nspecial[i][1]; j++)
      angleneigh[j] = special[i][j];
    for (m = nspecial[i][0]; m < nspecial[i][1]; m++) {
      tagint smallest = MAXTAGINT;
      for (j = m; j < nspecial[i][1]; j++) {
        if (angleneigh[j] < smallest) {
          smallest = angleneigh[j];
          k = j;
        }
        angleneigh[k] = angleneigh[m];
        angleneigh[m] = smallest;
      }
    }

    // assign xyz axis and fpole via only 1-2 connected atoms

    for (iframe = 0; iframe < nframe; iframe++) {
      xtype = xpole[itype][iframe];
      ytype = ypole[itype][iframe];
      ztype = zpole[itype][iframe];
      for (j12 = 0; j12 < nspecial[i][0]; j12++) {
        jneigh = bondneigh[j12];
        j = atom->map(jneigh);
        if (j < 0)
          error->one(FLERR,"AMOEBA kmpole() could not find bond partner");
        jtype = amtype[j];
        if (jtype == ztype) {
          for (k12 = 0; k12 < nspecial[i][0]; k12++) {
            if (k12 == j12) continue;
            kneigh = bondneigh[k12];
            k = atom->map(kneigh);
            if (k < 0)
              error->one(FLERR,"AMOEBA kmpole() could not find bond partner");
            ktype = amtype[k];
            if (ktype == xtype) {
              if (ytype == 0 && !flag) {
                flag = 1;
                xyzaxis[i][2] = ubuf(jneigh).d;
                xyzaxis[i][0] = ubuf(kneigh).d;
                xyzaxis[i][1] = 0.0;
                polaxe[i] = mpaxis[itype][iframe];
                for (j = 0; j < 13; j++)
                  pole[i][j] = fpole[itype][iframe][j];
              } else {
                for (m12 = 0; m12 < nspecial[i][0]; m12++) {
                  if (m12 == j12 || m12 == k12) continue;
                  mneigh = bondneigh[m12];
                  m = atom->map(mneigh);
                  if (m < 0)
                    error->one(FLERR,"AMOEBA kmpole() could not find bond partner");
                  mtype = amtype[m];
                  if (mtype == ytype && !flag) {
                    flag = 1;
                    xyzaxis[i][2] = ubuf(jneigh).d;
                    xyzaxis[i][0] = ubuf(kneigh).d;
                    xyzaxis[i][1] = ubuf(mneigh).d;
                    polaxe[i] = mpaxis[itype][iframe];
                    for (j = 0; j < 13; j++)
                      pole[i][j] = fpole[itype][iframe][j];
                  }
                }
              }
            }
          }
        }
      }
    }

    if (flag) continue;

    // assign xyz axis via 1-2 and 1-3 connected atoms

    for (iframe = 0; iframe < nframe; iframe++) {
      xtype = xpole[itype][iframe];
      ytype = ypole[itype][iframe];
      ztype = zpole[itype][iframe];
      for (j12 = 0; j12 < nspecial[i][0]; j12++) {
        jneigh = bondneigh[j12];
        j = atom->map(jneigh);
        if (j < 0) error->one(FLERR,"AMOEBA kmpole() could not find bond partner");
        jtype = amtype[j];
        if (jtype == ztype) {
          for (k13 = nspecial[i][0]; k13 < nspecial[i][1]; k13++) {
            kneigh = angleneigh[k13];
            k = atom->map(kneigh);
            if (k < 0) error->one(FLERR,"AMOEBA kmpole() could not find bond partner");
            ktype = amtype[k];
            path = false;
            for (m12 = 0; m12 < nspecial[k][0]; m12++)
              if (special[k][m12] == jneigh) path = true;
            if (!path) continue;

            if (ktype == xtype) {
              if (ytype == 0 && !flag) {
                flag = 1;
                xyzaxis[i][2] = ubuf(jneigh).d;
                xyzaxis[i][0] = ubuf(kneigh).d;
                xyzaxis[i][1] = 0.0;
                polaxe[i] = mpaxis[itype][iframe];
                for (j = 0; j < 13; j++)
                  pole[i][j] = fpole[itype][iframe][j];
              } else {
                for (m13 = nspecial[i][0]; m13 < nspecial[i][1]; m13++) {
                  if (m13 == k13) continue;
                  mneigh = angleneigh[m13];
                  m = atom->map(mneigh);
                  if (m < 0)
                    error->one(FLERR,"AMOEBA kmpole() could not find bond partner");
                  mtype = amtype[m];
                  path = false;
                  for (m12 = 0; m12 < nspecial[m][0]; m12++)
                    if (special[m][m12] == jneigh) path = true;
                  if (!path) continue;
                  if (mtype == ytype && !flag) {
                    flag = 1;
                    xyzaxis[i][2] = ubuf(jneigh).d;
                    xyzaxis[i][0] = ubuf(kneigh).d;
                    xyzaxis[i][1] = ubuf(mneigh).d;
                    polaxe[i] = mpaxis[itype][iframe];
                    for (j = 0; j < 13; j++)
                      pole[i][j] = fpole[itype][iframe][j];
                  }
                }
              }
            }
          }
        }
      }
    }

    if (flag) continue;

    // assign xyz axis via only a z-defining atom

    for (iframe = 0; iframe < nframe; iframe++) {
      xtype = xpole[itype][iframe];
      ytype = ypole[itype][iframe];
      ztype = zpole[itype][iframe];
      for (j12 = 0; j12 < nspecial[i][0]; j12++) {
        jneigh = bondneigh[j12];
        j = atom->map(jneigh);
        if (j < 0) error->one(FLERR,"AMOEBA kmpole() could not find bond partner");
        jtype = amtype[j];
        if (jtype == ztype) {
          if (xtype == 0 && !flag) {
            flag = 1;
            xyzaxis[i][2] = ubuf(jneigh).d;
            xyzaxis[i][0] = xyzaxis[i][1] = 0.0;
            polaxe[i] = mpaxis[itype][iframe];
            for (j = 0; j < 13; j++)
              pole[i][j] = fpole[itype][iframe][j];
          }
        }
      }
    }

    if (flag) continue;

    // assign xyz axis via no connected atoms

    for (iframe = 0; iframe < nframe; iframe++) {
      xtype = xpole[itype][iframe];
      ytype = ypole[itype][iframe];
      ztype = zpole[itype][iframe];
      if (ztype == 0 && !flag) {
        flag = 1;
        xyzaxis[i][2] = xyzaxis[i][1] = xyzaxis[i][0] = 0.0;
        polaxe[i] = mpaxis[itype][iframe];
        for (j = 0; j < 13; j++)
          pole[i][j] = fpole[itype][iframe][j];
      }
    }

    if (flag) continue;

    // flag error if could not assign xyz axis

    nmissing++;
  }

  // error check on missing settings

  int nmissing_all;
  MPI_Allreduce(&nmissing,&nmissing_all,1,MPI_INT,MPI_SUM,world);
  if (nmissing_all)
    error->all(FLERR, "Pair amoeba: {} multipole settings missing\n", nmissing_all);
}

/* ----------------------------------------------------------------------
   chkpole inverts atomic multipole moments as necessary
   at sites with chiral local reference frame definitions
   called every timestep for each atom I
------------------------------------------------------------------------- */

void PairAmoeba::chkpole(int i)
{
  bool check;
  int ib,ic,id;
  double xad,yad,zad;
  double xbd,ybd,zbd;
  double xcd,ycd,zcd;
  double c1,c2,c3,vol;

  double **pole = fixpole->astore;

  int *polaxe = atom->ivector[index_polaxe];
  double **xyzaxis = atom->darray[index_xyzaxis];
  tagint yaxisID = (tagint) ubuf(xyzaxis[i][1]).i;

  // test for chirality inversion
  // if not, return

  check = true;
  if (polaxe[i] != ZTHENX) check = false;
  if (yaxisID == 0) check = false;
  if (!check) return;

  ib = zaxis2local[i];
  ic = xaxis2local[i];
  id = yaxis2local[i];

  // compute the signed parallelpiped volume at chiral site

  double **x = atom->x;

  xad = x[i][0] - x[id][0];
  yad = x[i][1] - x[id][1];
  zad = x[i][2] - x[id][2];
  xbd = x[ib][0] - x[id][0];
  ybd = x[ib][1] - x[id][1];
  zbd = x[ib][2] - x[id][2];
  xcd = x[ic][0] - x[id][0];
  ycd = x[ic][1] - x[id][1];
  zcd = x[ic][2] - x[id][2];
  c1 = ybd*zcd - zbd*ycd;
  c2 = ycd*zad - zcd*yad;
  c3 = yad*zbd - zad*ybd;
  vol = xad*c1 + xbd*c2 + xcd*c3;

  // invert atomic multipole components involving the y-axis
  // flip sign in permanent yaxis, not yaxis2local

  if ((yaxisID < 0 && vol > 0.0) || (yaxisID > 0 && vol < 0.0)) {
    xyzaxis[i][1] = -ubuf(yaxisID).d;
    pole[i][2] = -pole[i][2];
    pole[i][5] = -pole[i][5];
    pole[i][7] = -pole[i][7];
    pole[i][9] = -pole[i][9];
    pole[i][11] = -pole[i][11];
  }
}

/* ----------------------------------------------------------------------
   rotmat finds the rotation matrix that rotates the local
   coordinate system into the global frame at a multipole site
   called every timestep for each atom I
------------------------------------------------------------------------- */

void PairAmoeba::rotmat(int i)
{
  int ix,iy,iz;
  double r,dot;
  double xi,yi,zi;
  double dx,dy,dz;
  double dx1,dy1,dz1;
  double dx2,dy2,dz2;
  double dx3,dy3,dz3;

  int *polaxe = atom->ivector[index_polaxe];

  // get coordinates and frame definition for the multipole site

  double **x = atom->x;

  xi = x[i][0];
  yi = x[i][1];
  zi = x[i][2];

  iz = zaxis2local[i];
  ix = xaxis2local[i];
  iy = yaxis2local[i];

  // use the identity matrix as the default rotation matrix

  rotate[0][0] = 1.0;
  rotate[0][1] = 0.0;
  rotate[0][2] = 0.0;
  rotate[2][0] = 0.0;
  rotate[2][1] = 0.0;
  rotate[2][2] = 1.0;

  // z-only frame rotation matrix elements for z-axis only

  if (polaxe[i] == ZONLY) {
    dx = x[iz][0] - xi;
    dy = x[iz][1] - yi;
    dz = x[iz][2] - zi;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    rotate[2][0] = dx / r;
    rotate[2][1] = dy / r;
    rotate[2][2] = dz / r;
    dx = 1.0;
    dy = 0.0;
    dz = 0.0;
    dot = rotate[2][0];
    if (fabs(dot) > 0.866) {
      dx = 0.0;
      dy = 1.0;
      dot = rotate[2][1];
    }
    dx = dx - dot*rotate[2][0];
    dy = dy - dot*rotate[2][1];
    dz = dz - dot*rotate[2][2];
    r = sqrt(dx*dx + dy*dy + dz*dz);
    rotate[0][0] = dx / r;
    rotate[0][1] = dy / r;
    rotate[0][2] = dz / r;

  // z-then-x frame rotation matrix elements for z- and x-axes

  } else if (polaxe[i] == ZTHENX) {
    dx = x[iz][0] - xi;
    dy = x[iz][1] - yi;
    dz = x[iz][2] - zi;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    rotate[2][0] = dx / r;
    rotate[2][1] = dy / r;
    rotate[2][2] = dz / r;
    dx = x[ix][0] - xi;
    dy = x[ix][1] - yi;
    dz = x[ix][2] - zi;
    dot = dx*rotate[2][0] + dy*rotate[2][1] + dz*rotate[2][2];
    dx = dx - dot*rotate[2][0];
    dy = dy - dot*rotate[2][1];
    dz = dz - dot*rotate[2][2];
    r = sqrt(dx*dx + dy*dy + dz*dz);
    rotate[0][0] = dx / r;
    rotate[0][1] = dy / r;
    rotate[0][2] = dz / r;

  // bisector frame rotation matrix elements for z- and x-axes

  } else if (polaxe[i] == BISECTOR) {
    dx = x[iz][0] - xi;
    dy = x[iz][1] - yi;
    dz = x[iz][2] - zi;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    dx1 = dx / r;
    dy1 = dy / r;
    dz1 = dz / r;
    dx = x[ix][0] - xi;
    dy = x[ix][1] - yi;
    dz = x[ix][2] - zi;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    dx2 = dx / r;
    dy2 = dy / r;
    dz2 = dz / r;
    dx = dx1 + dx2;
    dy = dy1 + dy2;
    dz = dz1 + dz2;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    rotate[2][0] = dx / r;
    rotate[2][1] = dy / r;
    rotate[2][2] = dz / r;
    dot = dx2*rotate[2][0] + dy2*rotate[2][1] + dz2*rotate[2][2];
    dx = dx2 - dot*rotate[2][0];
    dy = dy2 - dot*rotate[2][1];
    dz = dz2 - dot*rotate[2][2];
    r = sqrt(dx*dx + dy*dy + dz*dz);
    rotate[0][0] = dx / r;
    rotate[0][1] = dy / r;
    rotate[0][2] = dz / r;

  // z-bisect frame rotation matrix elements for z- and x-axes

  } else if (polaxe[i] == ZBISECT) {
    dx = x[iz][0] - xi;
    dy = x[iz][1] - yi;
    dz = x[iz][2] - zi;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    rotate[2][0] = dx / r;
    rotate[2][1] = dy / r;
    rotate[2][2] = dz / r;
    dx = x[ix][0] - xi;
    dy = x[ix][1] - yi;
    dz = x[ix][2] - zi;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    dx1 = dx / r;
    dy1 = dy / r;
    dz1 = dz / r;
    dx = x[iy][0] - xi;
    dy = x[iy][1] - yi;
    dz = x[iy][2] - zi;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    dx2 = dx / r;
    dy2 = dy / r;
    dz2 = dz / r;
    dx = dx1 + dx2;
    dy = dy1 + dy2;
    dz = dz1 + dz2;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    dx = dx / r;
    dy = dy / r;
    dz = dz / r;
    dot = dx*rotate[2][0] + dy*rotate[2][1] + dz*rotate[2][2];
    dx = dx - dot*rotate[2][0];
    dy = dy - dot*rotate[2][1];
    dz = dz - dot*rotate[2][2];
    r = sqrt(dx*dx + dy*dy + dz*dz);
    rotate[0][0] = dx / r;
    rotate[0][1] = dy / r;
    rotate[0][2] = dz / r;

  // 3-fold frame rotation matrix elements for z- and x-axes

  } else if (polaxe[i] == THREEFOLD) {
    dx = x[iz][0] - xi;
    dy = x[iz][1] - yi;
    dz = x[iz][2] - zi;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    dx1 = dx / r;
    dy1 = dy / r;
    dz1 = dz / r;
    dx = x[ix][0] - xi;
    dy = x[ix][1] - yi;
    dz = x[ix][2] - zi;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    dx2 = dx / r;
    dy2 = dy / r;
    dz2 = dz / r;
    dx = x[iy][0] - xi;
    dy = x[iy][1] - yi;
    dz = x[iy][2] - zi;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    dx3 = dx / r;
    dy3 = dy / r;
    dz3 = dz / r;
    dx = dx1 + dx2 + dx3;
    dy = dy1 + dy2 + dy3;
    dz = dz1 + dz2 + dz3;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    rotate[2][0] = dx / r;
    rotate[2][1] = dy / r;
    rotate[2][2] = dz / r;
    dot = dx2*rotate[2][0] + dy2*rotate[2][1] + dz2*rotate[2][2];
    dx = dx2 - dot*rotate[2][0];
    dy = dy2 - dot*rotate[2][1];
    dz = dz2 - dot*rotate[2][2];
    r = sqrt(dx*dx + dy*dy + dz*dz);
    rotate[0][0] = dx / r;
    rotate[0][1] = dy / r;
    rotate[0][2] = dz / r;
  }

  // finally, find rotation matrix elements for the y-axis

  rotate[1][0] = rotate[0][2]*rotate[2][1] - rotate[0][1]*rotate[2][2];
  rotate[1][1] = rotate[0][0]*rotate[2][2] - rotate[0][2]*rotate[2][0];
  rotate[1][2] = rotate[0][1]*rotate[2][0] - rotate[0][0]*rotate[2][1];
}

/* ----------------------------------------------------------------------
   rotsite rotates the local frame atomic multipoles at a
   specified site into the global coordinate frame by applying a rotation matrix
   called every timestep for each atom Isite
------------------------------------------------------------------------- */

void PairAmoeba::rotsite(int isite)
{
  int i,j,k,m;
  double mp[3][3];
  double rp[3][3];

  double **pole = fixpole->astore;

  // monopoles have the same value in any coordinate frame

  rpole[isite][0] = pole[isite][0];

  // rotate the dipoles to the global coordinate frame

  for (i = 1; i <= 3; i++) {
    rpole[isite][i] = 0.0;
    for (j = 1; j <= 3; j++)
      rpole[isite][i] += pole[isite][j] * rotate[j-1][i-1];
  }

  // rotate the quadrupoles to the global coordinate frame

  k = 4;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      mp[j][i] = pole[isite][k];
      rp[j][i] = 0.0;
      k++;
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      if (j < i) rp[j][i] = rp[i][j];
      else {
        for (k = 0; k < 3; k++)
          for (m = 0; m < 3; m++)
            rp[j][i] += rotate[k][i]*rotate[m][j] * mp[m][k];
      }
    }
  }

  k = 4;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      rpole[isite][k] = rp[j][i];
      k++;
    }
  }
}

/* ----------------------------------------------------------------------
   scan standard neighbor list and make it compatible with 1-5 neighbors
   if IJ entry is a 1-2,1-3,1-4 neighbor then adjust offset to SBBITS15
   else scan special15 to see if a 1-5 neighbor and adjust offset to SBBITS15
   else do nothing to IJ entry
------------------------------------------------------------------------- */

void PairAmoeba::add_onefive_neighbors()
{
  int i,j,k,ii,jj,which,inum,jnum,n15;
  tagint jtag;
  int *ilist,*jlist,*numneigh,**firstneigh;
  tagint *list15;

  // test for overflow of reduced allowed size of neighbor list

  if (atom->nlocal + atom->nghost > NEIGHMASK15)
    error->one(FLERR,"Pair amoeba neighbor list overflow");

  // neigh list

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // reset special neighbor flags to include 1-5 neighbors

  tagint *tag = atom->tag;
  int *nspecial15 = atom->nspecial15;
  tagint **special15 = atom->special15;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    n15 = nspecial15[i];
    list15 = special15[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      which = j >> SBBITS & 3;
      j &= NEIGHMASK;
      jtag = tag[j];

      if (!which) {
        for (k = 0; k < n15; k++) {
          if (list15[k] == jtag) {
            which = 4;
            break;
          }
        }
      }

      if (which) jlist[jj] = j ^ (which << SBBITS15);
    }
  }
}

/* ----------------------------------------------------------------------
   update local indices of hydrogen neighbors for owned and ghost atoms
   red2local = used for offset of hydrogen positions in Vdwl term
   called on reneighboring steps, only for AMOEBA
------------------------------------------------------------------------- */

void PairAmoeba::find_hydrogen_neighbors()
{
  int index;
  tagint id;

  // grab current ptr for redID
  // redID[i] = atom ID that atom I is bonded to
  // red2local[i] = local index of that atom
  // for furthest away ghost atoms, bond partner can be missing
  // in that case red2local = -1, but only an error if accessed in hal()

  double *redID = atom->dvector[index_redID];

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  for (int i = 0; i < nall; i++) {
    if (redID[i] == 0.0) red2local[i] = i;
    else {
      id = (tagint) ubuf(redID[i]).i;
      index = atom->map(id);
      if (index >= 0) index = domain->closest_image(i,index);
      red2local[i] = index;
    }
  }
}

/* ----------------------------------------------------------------------
   update local indices of bond topology neighbors for owned atoms
   xyz axis2local = used for multipole orientation
   called on reneighboring steps
------------------------------------------------------------------------- */

void PairAmoeba::find_multipole_neighbors()
{
  int index;
  tagint xaxisID,yaxisID,zaxisID;

  // grab current pts for xaxis,yaxis,zaxis
  // xyzaxis[i] = atom IDs that atom I uses for its multipole orientation
  // can be zero if not used, in which case set local index to self
  // yaxis can be negative, in which case use absolute value

  double **xyzaxis = atom->darray[index_xyzaxis];

  int nlocal = atom->nlocal;
  int nmissing = 0;

  for (int i = 0; i < nlocal; i++) {
    xaxisID = (tagint) ubuf(xyzaxis[i][0]).i;
    yaxisID = (tagint) ubuf(xyzaxis[i][1]).i;
    zaxisID = (tagint) ubuf(xyzaxis[i][2]).i;

    if (xaxisID) {
      index = atom->map(xaxisID);
      if (index == -1) nmissing++;
      else {
        index = domain->closest_image(i,index);
        xaxis2local[i] = index;
      }
    } else xaxis2local[i] = i;

    if (yaxisID) {
      if (xyzaxis[i][1] < 0) yaxisID = -yaxisID;
      index = atom->map(yaxisID);
      if (index == -1) nmissing++;
      else {
        index = domain->closest_image(i,index);
        yaxis2local[i] = index;
      }
    } else yaxis2local[i] = i;

    if (zaxisID) {
      index = atom->map(zaxisID);
      if (index == -1) nmissing++;
      else {
        index = domain->closest_image(i,index);
        zaxis2local[i] = index;
      }
    } else zaxis2local[i] = i;
  }

  // error check on missing neighbors

  int nmissing_all;
  MPI_Allreduce(&nmissing,&nmissing_all,1,MPI_INT,MPI_SUM,world);
  if (nmissing_all)
    error->all(FLERR, "Pair amoeba: {} multipole neighbors missing\n", nmissing_all);
}

/* ----------------------------------------------------------------------
   torque2force takes the torque values on a single site defined by
   a local coordinate frame and converts to Cartesian forces on
   the original site and sites specifying the local frame, also
   gives the x,y,z-force components needed for virial computation

   force distribution for the 3-fold local frame by Chao Lu,
   Ponder Lab, Washington University, July 2016

   literature reference:

   P. L. Popelier and A. J. Stone, "Formulae for the First and
   Second Derivatives of Anisotropic Potentials with Respect to
   Geometrical Parameters", Molecular Physics, 82, 411-425 (1994)

   C. Segui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
   Representation of Electrostatics in Classical Force Fields:
   Efficient Implementation of Multipolar Interactions in
   Biomolecular Simulations", Journal of Chemical Physics, 120,
   73-87 (2004)

   i = single site = local atom index
   trq = torque on that atom
   frc xyz = returned xyz force components
   f = force vector for local atoms, multiple atoms will be updated
------------------------------------------------------------------------- */

void PairAmoeba::torque2force(int i, double *trq,
                              double *frcx, double *frcy, double *frcz,
                              double **f)
{
  int j;
  int ia,ib,ic,id;
  int axetyp;
  double du,dv,dw,dot;
  double usiz,vsiz,wsiz;
  double psiz,rsiz,ssiz;
  double t1siz,t2siz;
  double uvsiz,uwsiz,vwsiz;
  double ursiz,ussiz;
  double vssiz,wssiz;
  double delsiz,dphiddel;
  double uvcos,urcos;
  double vscos,wscos;
  double upcos,vpcos,wpcos;
  double rwcos,rucos,rvcos;
  double ut1cos,ut2cos;
  double uvsin,ursin;
  double vssin,wssin;
  double rwsin,rusin,rvsin;
  double ut1sin,ut2sin;
  double dphidu,dphidv,dphidw;
  double dphidr,dphids;
  double u[3],v[3],w[3];
  double p[3],r[3],s[3];
  double t1[3],t2[3];
  double uv[3],uw[3],vw[3];
  double ur[3],us[3];
  double vs[3],ws[3];
  double del[3],eps[3];

  double **x = atom->x;

  // zero out force components on local frame-defining atoms

  for (j = 0; j < 3; j++) {
    frcz[j] = 0.0;
    frcx[j] = 0.0;
    frcy[j] = 0.0;
  }

  // get the local frame type and the frame-defining atoms

  int *polaxe = atom->ivector[index_polaxe];
  axetyp = polaxe[i];
  if (axetyp == NOFRAME) return;

  ia = zaxis2local[i];
  ib = i;
  ic = xaxis2local[i];
  id = yaxis2local[i];

  // construct the three rotation axes for the local frame

  u[0] = x[ia][0] - x[ib][0];
  u[1] = x[ia][1] - x[ib][1];
  u[2] = x[ia][2] - x[ib][2];
  usiz = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);

  if (axetyp != ZONLY) {
    v[0] = x[ic][0] - x[ib][0];
    v[1] = x[ic][1] - x[ib][1];
    v[2] = x[ic][2] - x[ib][2];
    vsiz = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  } else {
    v[0] = 1.0;
    v[1] = 0.0;
    v[2] = 0.0;
    vsiz = 1.0;
    dot = u[0] / usiz;
    if (fabs(dot) > 0.866) {
      v[0] = 0.0;
      v[1] = 1.0;
    }
  }

  if (axetyp == ZBISECT || axetyp == THREEFOLD) {
    w[0] = x[id][0] - x[ib][0];
    w[1] = x[id][1] - x[ib][1];
    w[2] = x[id][2] - x[ib][2];
  } else {
    w[0] = u[1]*v[2] - u[2]*v[1];
    w[1] = u[2]*v[0] - u[0]*v[2];
    w[2] = u[0]*v[1] - u[1]*v[0];
  }

  wsiz = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

  for (j = 0; j < 3; j++) {
    u[j] /= usiz;
    v[j] /= vsiz;
    w[j] /= wsiz;
  }

  // build some additional axes needed for the Z-Bisect method

  if (axetyp == ZBISECT) {
    r[0] = v[0] + w[0];
    r[1] = v[1] + w[1];
    r[2] = v[2] + w[2];
    rsiz = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    s[0] = u[1]*r[2] - u[2]*r[1];
    s[1] = u[2]*r[0] - u[0]*r[2];
    s[2] = u[0]*r[1] - u[1]*r[0];
    ssiz = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
    for (j = 0; j < 3; j++) {
      r[j] /= rsiz;
      s[j] /= ssiz;
    }
  }

  // find the perpendicular and angle for each pair of axes

  uv[0] = v[1]*u[2] - v[2]*u[1];
  uv[1] = v[2]*u[0] - v[0]*u[2];
  uv[2] = v[0]*u[1] - v[1]*u[0];
  uvsiz = sqrt(uv[0]*uv[0] + uv[1]*uv[1] + uv[2]*uv[2]);
  uw[0] = w[1]*u[2] - w[2]*u[1];
  uw[1] = w[2]*u[0] - w[0]*u[2];
  uw[2] = w[0]*u[1] - w[1]*u[0];
  uwsiz = sqrt(uw[0]*uw[0] + uw[1]*uw[1] + uw[2]*uw[2]);
  vw[0] = w[1]*v[2] - w[2]*v[1];
  vw[1] = w[2]*v[0] - w[0]*v[2];
  vw[2] = w[0]*v[1] - w[1]*v[0];
  vwsiz = sqrt(vw[0]*vw[0] + vw[1]*vw[1] + vw[2]*vw[2]);
  for (j = 0; j < 3; j++) {
    uv[j] /= uvsiz;
    uw[j] /= uwsiz;
    vw[j] /= vwsiz;
  }

  if (axetyp == ZBISECT) {
    ur[0] = r[1]*u[2] - r[2]*u[1];
    ur[1] = r[2]*u[0] - r[0]*u[2];
    ur[2] = r[0]*u[1] - r[1]*u[0];
    ursiz = sqrt(ur[0]*ur[0] + ur[1]*ur[1] + ur[2]*ur[2]);
    us[0] = s[1]*u[2] - s[2]*u[1];
    us[1] = s[2]*u[0] - s[0]*u[2];
    us[2] = s[0]*u[1] - s[1]*u[0];
    ussiz = sqrt(us[0]*us[0] + us[1]*us[1] + us[2]*us[2]);
    vs[0] = s[1]*v[2] - s[2]*v[1];
    vs[1] = s[2]*v[0] - s[0]*v[2];
    vs[2] = s[0]*v[1] - s[1]*v[0];
    vssiz = sqrt(vs[0]*vs[0] + vs[1]*vs[1] + vs[2]*vs[2]);
    ws[0] = s[1]*w[2] - s[2]*w[1];
    ws[1] = s[2]*w[0] - s[0]*w[2];
    ws[2] = s[0]*w[1] - s[1]*w[0];
    wssiz = sqrt(ws[0]*ws[0] + ws[1]*ws[1] + ws[2]*ws[2]);
    for (j = 0; j < 3; j++) {
      ur[j] /= ursiz;
      us[j] /= ussiz;
      vs[j] /= vssiz;
      ws[j] /= wssiz;
    }
  }

  // get sine and cosine of angles between the rotation axes

  uvcos = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
  uvsin = sqrt(1.0 - uvcos*uvcos);
  if (axetyp == ZBISECT) {
    urcos = u[0]*r[0] + u[1]*r[1] + u[2]*r[2];
    ursin = sqrt(1.0 - urcos*urcos);
    vscos = v[0]*s[0] + v[1]*s[1] + v[2]*s[2];
    vssin = sqrt(1.0 - vscos*vscos);
    wscos = w[0]*s[0] + w[1]*s[1] + w[2]*s[2];
    wssin = sqrt(1.0 - wscos*wscos);
  }

  // compute the projection of v and w onto the ru-plane

  if (axetyp == ZBISECT) {
    for (j = 0; j < 3; j++) {
      t1[j] = v[j] - s[j]*vscos;
      t2[j] = w[j] - s[j]*wscos;
    }
    t1siz = sqrt(t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2]);
    t2siz = sqrt(t2[0]*t2[0] + t2[1]*t2[1] + t2[2]*t2[2]);
    for (j = 0; j < 3; j++) {
      t1[j] /= t1siz;
      t2[j] /= t2siz;
    }
    ut1cos = u[0]*t1[0] + u[1]*t1[1] + u[2]*t1[2];
    ut1sin = sqrt(1.0 - ut1cos*ut1cos);
    ut2cos = u[0]*t2[0] + u[1]*t2[1] + u[2]*t2[2];
    ut2sin = sqrt(1.0 - ut2cos*ut2cos);
  }

  // negative of dot product of torque with unit vectors gives
  // result of infinitesimal rotation along these vectors

  dphidu = -trq[0]*u[0] - trq[1]*u[1] - trq[2]*u[2];
  dphidv = -trq[0]*v[0] - trq[1]*v[1] - trq[2]*v[2];
  dphidw = -trq[0]*w[0] - trq[1]*w[1] - trq[2]*w[2];
  if (axetyp == ZBISECT) {
    dphidr = -trq[0]*r[0] - trq[1]*r[1] - trq[2]*r[2];
    dphids = -trq[0]*s[0] - trq[1]*s[1] - trq[2]*s[2];
  }

  // force distribution for the Z-Only local coordinate method

  if (axetyp == ZONLY) {
    for (j = 0; j < 3; j++) {
      du = uv[j]*dphidv/(usiz*uvsin) + uw[j]*dphidw/usiz;
      f[ia][j] -= du;
      f[ib][j] += du;
      frcz[j] += du;
    }

  // force distribution for the Z-then-X local coordinate method

  } else if (axetyp == ZTHENX) {
    for (j = 0; j < 3; j++) {
      du = uv[j]*dphidv/(usiz*uvsin) + uw[j]*dphidw/usiz;
      dv = -uv[j]*dphidu/(vsiz*uvsin);
      f[ia][j] -= du;
      f[ic][j] -= dv;
      f[ib][j] += du + dv;
      frcz[j] += du;
      frcx[j] += dv;
    }

  // force distribution for the Bisector local coordinate method

  } else if (axetyp == BISECTOR) {
    for (j = 0; j < 3; j++) {
      du = uv[j]*dphidv/(usiz*uvsin) + 0.5*uw[j]*dphidw/usiz;
      dv = -uv[j]*dphidu/(vsiz*uvsin) + 0.5*vw[j]*dphidw/vsiz;
      f[ia][j] -= du;
      f[ic][j] -= dv;
      f[ib][j] += du + dv;
      frcz[j] += du;
      frcx[j] += dv;
    }

  // force distribution for the Z-Bisect local coordinate method

  } else if (axetyp == ZBISECT) {
    for (j = 0; j < 3; j++) {
      du = ur[j]*dphidr/(usiz*ursin) + us[j]*dphids/usiz;
      dv = (vssin*s[j]-vscos*t1[j])*dphidu / (vsiz*(ut1sin+ut2sin));
      dw = (wssin*s[j]-wscos*t2[j])*dphidu / (wsiz*(ut1sin+ut2sin));
      f[ia][j] -= du;
      f[ic][j] -= dv;
      f[id][j] -= dw;
      f[ib][j] += du + dv + dw;
      frcz[j] += du;
      frcx[j] += dv;
      frcy[j] += dw;
    }

  // force distribution for the 3-Fold local coordinate method

  } else if (axetyp == THREEFOLD) {
    p[0] = u[0] + v[0] + w[0];
    p[1] = u[1] + v[1] + w[1];
    p[2] = u[2] + v[2] + w[2];
    psiz = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) ;
    for (j = 0; j < 3; j++) p[j] /= psiz;

    wpcos = w[0]*p[0] + w[1]*p[1] + w[2]*p[2];
    upcos = u[0]*p[0] + u[1]*p[1] + u[2]*p[2];
    vpcos = v[0]*p[0] + v[1]*p[1] + v[2]*p[2];
    r[0] = u[0] + v[0];
    r[1] = u[1] + v[1];
    r[2] = u[2] + v[2];
    rsiz = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    for (j = 0; j < 3; j++) r[j] /= rsiz;

    rwcos = r[0]*w[0] + r[1]*w[1] + r[2]*w[2];
    rwsin = sqrt(1.0 - rwcos*rwcos);
    dphidr = -trq[0]*r[0] - trq[1]*r[1] - trq[2]*r[2];
    del[0] = r[1]*w[2] - r[2]*w[1];
    del[1] = r[2]*w[0] - r[0]*w[2] ;
    del[2] = r[0]*w[1] - r[1]*w[0];
    delsiz = sqrt(del[0]*del[0] + del[1]*del[1] + del[2]*del[2]);
    for (j = 0; j < 3; j++) del[j] /= delsiz;

    dphiddel = -trq[0]*del[0] - trq[1]*del[1] - trq[2]*del[2];
    eps[0] = del[1]*w[2] - del[2]*w[1];
    eps[1] = del[2]*w[0] - del[0]*w[2];
    eps[2] = del[0]*w[1] - del[1]*w[0];
    for (j = 0; j < 3; j++) {
      dw = del[j]*dphidr/(wsiz*rwsin) + eps[j]*dphiddel*wpcos/(wsiz*psiz);
      f[id][j] -= dw;
      f[ib][j] += dw;
      frcy[j] += dw;
    }
    r[0] = v[0] + w[0];
    r[1] = v[1] + w[1];
    r[2] = v[2] + w[2];
    rsiz = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    for (j = 0; j < 3; j++) r[j] /= rsiz;

    rucos = r[0]*u[0] + r[1]*u[1] + r[2]*u[2];
    rusin = sqrt(1.0 - rucos*rucos) ;
    dphidr = -trq[0]*r[0] - trq[1]*r[1] - trq[2]*r[2];
    del[0] = r[1]*u[2] - r[2]*u[1];
    del[1] = r[2]*u[0] - r[0]*u[2];
    del[2] = r[0]*u[1] - r[1]*u[0];
    delsiz = sqrt(del[0]*del[0] + del[1]*del[1] + del[2]*del[2]);
    for (j = 0; j < 3; j++) del[j] /= delsiz;

    dphiddel = -trq[0]*del[0] - trq[1]*del[1] - trq[2]*del[2];
    eps[0] = del[1]*u[2] - del[2]*u[1];
    eps[1] = del[2]*u[0] - del[0]*u[2];
    eps[2] = del[0]*u[1] - del[1]*u[0];
    for (j = 0; j < 3; j++) {
      du = del[j]*dphidr/(usiz*rusin) + eps[j]*dphiddel*upcos/(usiz*psiz);
      f[ia][j] -= du;
      f[ib][j] += du;
      frcz[j] += du;
    }
    r[0] = u[0] + w[0];
    r[1] = u[1] + w[1];
    r[2] = u[2] + w[2];
    rsiz = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    for (j = 0; j < 3; j++) r[j] /= rsiz;

    rvcos = r[0]*v[0] + r[1]*v[1] + r[2]*v[2] ;
    rvsin = sqrt(1.0 - rvcos*rvcos);
    dphidr = -trq[0]*r[0] - trq[1]*r[1] - trq[2]*r[2];
    del[0] = r[1]*v[2] - r[2]*v[1];
    del[1] = r[2]*v[0] - r[0]*v[2];
    del[2] = r[0]*v[1] - r[1]*v[0];
    delsiz = sqrt(del[0]*del[0] + del[1]*del[1] + del[2]*del[2]);
    for (j = 0; j < 3; j++) del[j] /= delsiz;

    dphiddel = -trq[0]*del[0] - trq[1]*del[1] - trq[2]*del[2];
    eps[0] = del[1]*v[2] - del[2]*v[1];
    eps[1] = del[2]*v[0] - del[0]*v[2];
    eps[2] = del[0]*v[1] - del[1]*v[0];
    for (j = 0; j < 3; j++) {
      dv = del[j]*dphidr/(vsiz*rvsin) + eps[j]*dphiddel*vpcos/(vsiz*psiz);
      f[ic][j] -= dv;
      f[ib][j] += dv;
      frcx[j] += dv;
    }
  }
}
