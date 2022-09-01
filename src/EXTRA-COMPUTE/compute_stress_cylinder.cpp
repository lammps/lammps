/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_stress_cylinder.h"

#include "atom.h"
#include "citeme.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <mpi.h>

using namespace LAMMPS_NS;
using namespace MathConst;
using MathSpecial::square;

/*-----------------------------------------------------------------------------------
  Contributing authors: Cody K. Addington (North Carolina State University)
                        (Kinetic contribution) : Olav Galteland,
                        (Norwegian University of Science and Technology),
                        olav.galteland@ntnu.no
------------------------------------------------------------------------------------*/

static const char cite_compute_stress_cylinder[] =
    "compute stress/cylinder: doi:10.1063/1.5037054\n\n"
    "@Article{Addington,\n"
    " author = {C. K. Addington and Y. Long and K. E. Gubbins},\n"
    " title = {The Pressure in Interfaces Having Cylindrical Geometry},\n"
    " journal = {J.~Chem.\\ Phys.},\n"
    " year =    2018,\n"
    " volume =  149,\n"
    " number =  8,\n"
    " pages =   {084109}\n"
    "}\n\n";

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Calculate the configurational components of the stress tensor in
  cylindrical geometry, according to the formulation of Addington et al. (2018)
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

ComputeStressCylinder::ComputeStressCylinder(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), Pvr_temp(nullptr), Pvr_all(nullptr), Pvz_temp(nullptr),
    Pvz_all(nullptr), Pvphi_temp(nullptr), Pvphi_all(nullptr), R(nullptr), Rinv(nullptr),
    R2(nullptr), PrAinv(nullptr), PzAinv(nullptr), R2kin(nullptr), density_temp(nullptr),
    invVbin(nullptr), density_all(nullptr), tangent(nullptr), ephi_x(nullptr), ephi_y(nullptr),
    binz(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_compute_stress_cylinder);
  if ((narg != 7) && (narg != 9)) error->all(FLERR, "Illegal compute stress/cylinder command");

  zlo = utils::numeric(FLERR, arg[3], false, lmp);
  zhi = utils::numeric(FLERR, arg[4], false, lmp);
  Rmax = utils::numeric(FLERR, arg[5], false, lmp);
  bin_width = utils::numeric(FLERR, arg[6], false, lmp);

  // Option to include/exclude kinetic contribution. Default is to include
  kinetic_flag = 1;
  int iarg = 7;
  if (narg > iarg) {
    if (strcmp("ke", arg[iarg]) == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Invalid compute stress/cylinder command");
      kinetic_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Unknown compute stress/cylinder command");
  }

  if ((bin_width <= 0.0) || (bin_width > Rmax))
    error->all(FLERR, "Illegal compute stress/cylinder command");
  if ((zhi < zlo) || ((zhi - zlo) < bin_width))
    error->all(FLERR, "Illegal compute stress/cylinder command");
  if ((zhi > domain->boxhi[2]) || (zlo < domain->boxlo[2]))
    error->all(FLERR, "Illegal compute stress/cylinder command");

  nbins = (int) (Rmax / bin_width);
  nzbins = (int) ((zhi - zlo) / bin_width);

  // NOTE: at 2^22 = 4.2M bins, we will be close to exhausting allocatable
  // memory on a 32-bit environment. so we use this as an upper limit.

  if ((nbins < 1) || (nzbins < 1) || (nbins > 2 << 22) || (nzbins > 2 << 22))
    error->all(FLERR, "Illegal compute stress/cylinder command");

  array_flag = 1;
  vector_flag = 0;
  extarray = 0;
  size_array_cols = 5;    // r, number density, pvr, pvphi, pz
  size_array_rows = nbins;

  if (kinetic_flag == 1) {
    size_array_cols = 8;    // r, number density, pkr, pkphi, pkz, pvr, pvphi, pz
    Pkr_temp = new double[nbins];
    Pkr_all = new double[nbins];
    Pkz_temp = new double[nbins];
    Pkz_all = new double[nbins];
    Pkphi_temp = new double[nbins];
    Pkphi_all = new double[nbins];
  }
  Pvr_temp = new double[nbins];
  Pvr_all = new double[nbins];
  Pvz_temp = new double[nbins];
  Pvz_all = new double[nbins];
  Pvphi_temp = new double[nbins];
  Pvphi_all = new double[nbins];
  R = new double[nbins];
  R2 = new double[nbins];
  PrAinv = new double[nbins];
  PzAinv = new double[nbins];
  Rinv = new double[nbins];
  binz = new double[nzbins];

  R2kin = new double[nbins];
  density_temp = new double[nbins];
  invVbin = new double[nbins];
  density_all = new double[nbins];

  nphi = 360;
  tangent = new double[nphi];
  ephi_x = new double[nphi];
  ephi_y = new double[nphi];

  memory->create(array, size_array_rows, size_array_cols, "PN:array");

  nktv2p = force->nktv2p;
}

/* ---------------------------------------------------------------------- */

ComputeStressCylinder::~ComputeStressCylinder()
{
  memory->destroy(array);
  if (kinetic_flag == 1) {
    delete[] Pkr_temp;
    delete[] Pkr_all;
    delete[] Pkz_temp;
    delete[] Pkz_all;
    delete[] Pkphi_temp;
    delete[] Pkphi_all;
  }
  delete[] R;
  delete[] Rinv;
  delete[] R2;
  delete[] R2kin;
  delete[] invVbin;
  delete[] density_temp;
  delete[] density_all;
  delete[] tangent;
  delete[] ephi_x;
  delete[] ephi_y;
  delete[] Pvr_temp;
  delete[] Pvr_all;
  delete[] Pvz_temp;
  delete[] Pvz_all;
  delete[] Pvphi_temp;
  delete[] Pvphi_all;
  delete[] PrAinv;
  delete[] PzAinv;
  delete[] binz;
}

/* ---------------------------------------------------------------------- */

void ComputeStressCylinder::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "No pair style is defined for compute stress/cylinder");
  if (force->pair->single_enable == 0)
    error->all(FLERR, "Pair style does not support compute stress/cylinder");

  double phi;
  for (int iphi = 0; iphi < nphi; iphi++) {
    phi = ((double) iphi) * MY_PI / 180.0;
    tangent[iphi] = tan(phi);
    ephi_x[iphi] = -sin(phi);
    ephi_y[iphi] = cos(phi);
  }
  for (int iq = 0; iq < nbins; iq++) {
    R[iq] = ((double) iq + 0.5) * bin_width;
    Rinv[iq] = 1.0 / R[iq];
    R2[iq] = R[iq] * R[iq];
    R2kin[iq] = (((double) iq) + 1.0) * bin_width;
    R2kin[iq] *= R2kin[iq];
    PrAinv[iq] = 1.0 / (2.0 * MY_PI * (zhi - zlo) * R[iq]);
  }
  PphiAinv = 1.0 / ((zhi - zlo) * bin_width * 2.0 * (double) nphi);

  invVbin[0] = 1.0 / ((zhi - zlo) * MY_PI * R2kin[0]);
  PzAinv[0] = 1.0 / (MY_PI * R2kin[0] * ((double) nzbins));

  for (int jq = 1; jq < nbins; jq++) {
    invVbin[jq] = 1.0 / ((zhi - zlo) * MY_PI * (R2kin[jq] - R2kin[jq - 1]));
    PzAinv[jq] = 1.0 / (MY_PI * (R2kin[jq] - R2kin[jq - 1]) * ((double) nzbins));
  }

  // need an occasional half neighbor list
  neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);

  for (int zzz = 0; zzz < nzbins; zzz++) binz[zzz] = (((double) zzz) + 0.5) * bin_width + zlo;
}

/* ---------------------------------------------------------------------- */

void ComputeStressCylinder::init_list(int /* id */, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   count pairs and compute pair info on this proc
   only count pair once if newton_pair is off
   both atom I,J must be in group
   if flag is set, compute requested info about pair
------------------------------------------------------------------------- */

void ComputeStressCylinder::compute_array()
{
  invoked_array = update->ntimestep;

  int ibin;
  // clear pressures
  for (ibin = 0; ibin < nbins; ibin++) {
    if (kinetic_flag == 1) {
      Pkr_temp[ibin] = 0.0;
      Pkr_all[ibin] = 0.0;
      Pkphi_temp[ibin] = 0.0;
      Pkphi_all[ibin] = 0.0;
      Pkz_temp[ibin] = 0.0;
      Pkz_all[ibin] = 0.0;
    }
    density_temp[ibin] = 0.0;
    density_all[ibin] = 0.0;
    Pvr_temp[ibin] = 0.0;
    Pvr_all[ibin] = 0.0;
    Pvphi_temp[ibin] = 0.0;
    Pvphi_all[ibin] = 0.0;
    Pvz_temp[ibin] = 0.0;
    Pvz_all[ibin] = 0.0;
  }

  // what processor am I?
  int me;
  MPI_Comm_rank(world, &me);

  int i, j, ii, jj, inum, jnum, itype, jtype;
  tagint itag, jtag;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, fpair, factor_coul, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double vr, vp;
  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)
  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // calculate number density and kinetic contribution (by radius)
  double temp_R2;
  for (i = 0; i < nlocal; i++)
    if ((x[i][2] < zhi) && (x[i][2] > zlo)) {
      temp_R2 = x[i][0] * x[i][0] + x[i][1] * x[i][1];
      if (temp_R2 > R2kin[nbins - 1]) continue;    // outside of Rmax

      for (j = 0; j < nbins; j++)
        if (temp_R2 < R2kin[j]) break;

      density_temp[j] += invVbin[j];

      // Check if kinetic option is set to yes
      if (kinetic_flag == 1) {
        if ((temp_R2 != 0.0) && (x[i][0] != 0.0)) {
          // Radial velocity times R
          vr = (x[i][0] * v[i][0] + x[i][1] * v[i][1]);
          // Azimuthal velocity divided by R
          vp = (v[i][1] / x[i][0] - x[i][1] * v[i][0] / (x[i][0] * x[i][0])) /
              (square(x[i][1] / x[i][0]) + 1.0);

          Pkr_temp[j] += mass[type[i]] * vr * vr / temp_R2;
          Pkphi_temp[j] += mass[type[i]] * temp_R2 * vp * vp;
          Pkz_temp[j] += mass[type[i]] * v[i][2] * v[i][2];
        }
      }
    }
  MPI_Allreduce(density_temp, density_all, nbins, MPI_DOUBLE, MPI_SUM, world);
  for (i = 0; i < nbins; i++) array[i][1] = density_all[i];    // NEW

  // loop over neighbors of my atoms
  // skip if I or J are not in group
  // for newton = 0 and J = ghost atom,
  //   need to insure I,J pair is only output by one proc
  //   use same itag,jtag logic as in Neighbor::neigh_half_nsq()
  // for flag = 0, just count pair interactions within force cutoff
  // for flag = 1, calculate requested output fields

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  double r1 = 0.0;
  double r2 = 0.0;
  double risq, rjsq;
  double A, B, C, D;
  double alpha1, alpha2;
  double xi, yi, zi, dx, dy, dz;
  double xR, yR, zR, fn;
  double alpha, xL, yL, zL, L2, ftphi, ftz;
  double sqrtD;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itag = tag[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    r1 = x[i][0] * x[i][0] + x[i][1] * x[i][1];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;

      // itag = jtag is possible for long cutoffs that include images of self
      // do calculation only on appropriate processor
      if (newton_pair == 0 && j >= nlocal) {
        jtag = tag[j];
        if (itag > jtag) {
          if ((itag + jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag + jtag) % 2 == 1) continue;
        } else {
          if (x[j][2] < ztmp) continue;
          if (x[j][2] == ztmp) {
            if (x[j][1] < ytmp) continue;
            if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
          }
        }
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      r2 = x[j][0] * x[j][0] + x[j][1] * x[j][1];

      // ri is smaller of r1 and r2
      if (r2 < r1) {
        risq = r2;
        rjsq = r1;
        xi = x[j][0];
        yi = x[j][1];
        zi = x[j][2];
        dx = x[i][0] - x[j][0];
        dy = x[i][1] - x[j][1];
        dz = x[i][2] - x[j][2];
      } else {
        risq = r1;
        rjsq = r2;
        xi = x[i][0];
        yi = x[i][1];
        zi = x[i][2];
        dx = x[j][0] - x[i][0];
        dy = x[j][1] - x[i][1];
        dz = x[j][2] - x[i][2];
      }

      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      if (rsq >= cutsq[itype][jtype]) continue;

      pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

      A = dx * dx + dy * dy;
      B = 2.0 * (xi * dx + yi * dy);

      // normal pressure contribution P_rhorho
      for (ibin = 0; ibin < nbins; ibin++) {
        // completely inside of R
        if (rjsq < R2[ibin]) continue;

        C = risq - R2[ibin];
        D = B * B - 4.0 * A * C;

        // completely outside of R or zero size bin
        if ((D < 0.0) || (A == 0.0)) continue;

        sqrtD = sqrt(D);
        alpha1 = 0.5 * (-B + sqrtD) / A;
        alpha2 = 0.5 * (-B - sqrtD) / A;

        if ((alpha1 > 0.0) && (alpha1 < 1.0)) {
          zR = zi + alpha1 * dz;
          if ((zR < zhi) && (zR > zlo)) {
            xR = xi + alpha1 * dx;
            yR = yi + alpha1 * dy;
            fn = fpair * fabs(xR * dx + yR * dy);

            Pvr_temp[ibin] += fn;
          }
        }
        if ((alpha2 > 0.0) && (alpha2 < 1.0)) {
          zR = zi + alpha2 * dz;
          if ((zR < zhi) && (zR > zlo)) {
            xR = xi + alpha2 * dx;
            yR = yi + alpha2 * dy;
            fn = fpair * fabs(xR * dx + yR * dy);

            Pvr_temp[ibin] += fn;
          }
        }
      }

      // azimuthal pressure contribution (P_phiphi)
      for (int iphi = 0; iphi < nphi; iphi++) {
        if ((dx * tangent[iphi] - dy) == 0.0) continue;
        alpha = (yi - xi * tangent[iphi]) / (dx * tangent[iphi] - dy);

        // no intersection with phi surface
        if ((alpha >= 1.0) || (alpha <= 0.0)) continue;

        // no contribution (outside of averaging region)
        zL = zi + alpha * dz;
        if ((zL > zhi) || (zL < zlo)) continue;

        xL = xi + alpha * dx;
        yL = yi + alpha * dy;

        L2 = xL * xL + yL * yL;

        // no intersection (outside of Rmax)
        if (L2 > R2kin[nbins - 1]) continue;

        ftphi = fabs(dx * ephi_x[iphi] + dy * ephi_y[iphi]) * fpair;

        // add to appropriate bin
        for (ibin = 0; ibin < nbins; ibin++)
          if (L2 < R2kin[ibin]) {
            Pvphi_temp[ibin] += ftphi;
            break;
          }
      }

      // z pressure contribution (P_zz)
      for (int zbin = 0; zbin < nzbins; zbin++) {
        // check if interaction contributes
        if ((x[i][2] > binz[zbin]) && (x[j][2] > binz[zbin])) continue;
        if ((x[i][2] < binz[zbin]) && (x[j][2] < binz[zbin])) continue;

        alpha = (binz[zbin] - zi) / dz;

        xL = xi + alpha * dx;
        yL = yi + alpha * dy;

        L2 = xL * xL + yL * yL;

        if (L2 > R2kin[nbins - 1]) continue;

        ftz = fabs(dz) * fpair;

        // add to appropriate bin
        for (ibin = 0; ibin < nbins; ibin++)
          if (L2 < R2kin[ibin]) {
            Pvz_temp[ibin] += ftz;
            break;
          }
      }
    }
  }

  // calculate pressure (force over area)
  for (ibin = 0; ibin < nbins; ibin++) {
    if (kinetic_flag == 1) {
      Pkr_temp[ibin] *= invVbin[ibin];
      Pkphi_temp[ibin] *= invVbin[ibin];
      Pkz_temp[ibin] *= invVbin[ibin];
    }
    Pvr_temp[ibin] *= PrAinv[ibin] * Rinv[ibin];
    Pvphi_temp[ibin] *= PphiAinv;
    Pvz_temp[ibin] *= PzAinv[ibin];
  }

  // communicate these values across processors
  if (kinetic_flag == 1) {
    MPI_Allreduce(Pkr_temp, Pkr_all, nbins, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(Pkphi_temp, Pkphi_all, nbins, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(Pkz_temp, Pkz_all, nbins, MPI_DOUBLE, MPI_SUM, world);
  }
  MPI_Allreduce(Pvr_temp, Pvr_all, nbins, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(Pvphi_temp, Pvphi_all, nbins, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(Pvz_temp, Pvz_all, nbins, MPI_DOUBLE, MPI_SUM, world);

  // populate array
  for (ibin = 0; ibin < nbins; ibin++) {
    array[ibin][0] = R[ibin];
    if (kinetic_flag == 1) {
      array[ibin][2] = Pkr_all[ibin] * nktv2p;
      array[ibin][3] = Pkphi_all[ibin] * nktv2p;
      array[ibin][4] = Pkz_all[ibin] * nktv2p;
      array[ibin][5] = Pvr_all[ibin] * nktv2p;
      array[ibin][6] = Pvphi_all[ibin] * nktv2p;
      array[ibin][7] = Pvz_all[ibin] * nktv2p;
    } else {
      array[ibin][2] = Pvr_all[ibin] * nktv2p;
      array[ibin][3] = Pvphi_all[ibin] * nktv2p;
      array[ibin][4] = Pvz_all[ibin] * nktv2p;
    }
  }
}

/* ----------------------------------------------------------------------
memory usage of data
------------------------------------------------------------------------- */
double ComputeStressCylinder::memory_usage()
{
  double bytes =
      (3.0 * (double) nphi + 16.0 * (double) nbins + (5.0 + 3.0 * kinetic_flag) * (double) nbins) *
      sizeof(double);
  return bytes;
}
