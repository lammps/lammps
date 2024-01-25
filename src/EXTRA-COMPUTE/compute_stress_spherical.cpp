/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS dir1ectory.
------------------------------------------------------------------------- */

#include "compute_stress_spherical.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;
using MathSpecial::cube;
using MathSpecial::square;

static constexpr double SMALL = 1.0e-10;

/*-----------------------------------------------------------------------------------
  Contributing author: Olav Galteland (Norwegian University of Science and Technology)
                        olav.galteland@ntnu.no
------------------------------------------------------------------------------------*/

static const char cite_compute_stress_sphere[] =
    "compute stress/spherical: doi:10.48550/arXiv.2201.13060\n\n"
    "@article{galteland2022defining,\n"
    "title={Defining the Pressures of a Fluid in a Nanoporous, Heterogeneous Medium},\n"
    "author={Galteland, Olav and Rauter, Michael T and Varughese, Kevin K and Bedeaux, Dick and "
    "   Kjelstrup, Signe},\n"
    "journal={arXiv preprint arXiv:2201.13060},\n"
    "year={2022}\n"
    "}\n\n";

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

ComputeStressSpherical::ComputeStressSpherical(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), dens(nullptr), pkrr(nullptr), pktt(nullptr), pkpp(nullptr),
    pcrr(nullptr), pctt(nullptr), pcpp(nullptr), tdens(nullptr), tpkrr(nullptr), tpktt(nullptr),
    tpkpp(nullptr), tpcrr(nullptr), tpctt(nullptr), tpcpp(nullptr), list(nullptr)
{

  if (lmp->citeme) lmp->citeme->add(cite_compute_stress_sphere);
  if (narg != 8)
    error->all(FLERR, "Illegal compute stress/spherical command. Illegal number of arguments.");

  x0 = utils::numeric(FLERR, arg[3], false, lmp);
  y0 = utils::numeric(FLERR, arg[4], false, lmp);
  z0 = utils::numeric(FLERR, arg[5], false, lmp);
  bin_width = utils::numeric(FLERR, arg[6], false, lmp);
  Rmax = utils::numeric(FLERR, arg[7], false, lmp);
  nbins = (int) (Rmax / bin_width) + 1;
  double tmp_width = Rmax / nbins;
  if ((fabs(bin_width - tmp_width) > SMALL) && (comm->me == 0))
    utils::logmesg(lmp, "Adjusting bin width for compute {} from {:.6f} to {:.6f}\n", style,
                   bin_width, tmp_width);
  bin_width = tmp_width;

  if (bin_width <= 0.0)
    error->all(FLERR, "Illegal compute stress/spherical command. Bin width must be > 0");

  array_flag = 1;
  vector_flag = 0;
  extarray = 0;
  size_array_cols = 8;    // r, dens, pkrr, pktt, pkpp, pcrr, pctt, pcpp
  size_array_rows = nbins;

  memory->create(invV, nbins, "compute/stress/spherical:invV");
  memory->create(dens, nbins, "compute/stress/spherical:dens");
  memory->create(pkrr, nbins, "compute/stress/spherical:pkrr");
  memory->create(pktt, nbins, "compute/stress/spherical:pktt");
  memory->create(pkpp, nbins, "compute/stress/spherical:pkpp");
  memory->create(pcrr, nbins, "compute/stress/spherical:pcrr");
  memory->create(pctt, nbins, "compute/stress/spherical:pctt");
  memory->create(pcpp, nbins, "compute/stress/spherical:pcpp");
  memory->create(tdens, nbins, "compute/stress/spherical:tdens");
  memory->create(tpkrr, nbins, "compute/stress/spherical:tpkrr");
  memory->create(tpktt, nbins, "compute/stress/spherical:tpktt");
  memory->create(tpkpp, nbins, "compute/stress/spherical:tpkpp");
  memory->create(tpcrr, nbins, "compute/stress/spherical:tpcrr");
  memory->create(tpctt, nbins, "compute/stress/spherical:tpctt");
  memory->create(tpcpp, nbins, "compute/stress/spherical:tpcpp");
  memory->create(array, size_array_rows, size_array_cols, "compute/stress/spherical:array");
}

/* ---------------------------------------------------------------------- */

ComputeStressSpherical::~ComputeStressSpherical()
{
  memory->destroy(invV);
  memory->destroy(dens);
  memory->destroy(pkrr);
  memory->destroy(pktt);
  memory->destroy(pkpp);
  memory->destroy(pcrr);
  memory->destroy(pctt);
  memory->destroy(pcpp);
  memory->destroy(tdens);
  memory->destroy(tpkrr);
  memory->destroy(tpktt);
  memory->destroy(tpkpp);
  memory->destroy(tpcrr);
  memory->destroy(tpctt);
  memory->destroy(tpcpp);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeStressSpherical::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "No pair style is defined for compute stress/spherical");
  if (force->pair->single_enable == 0)
    error->all(FLERR, "Pair style does not support compute stress/spherical");

  // Inverse volume of each spherical shell (bin)
  for (int bin = 0; bin < nbins; bin++)
    invV[bin] = 0.75 / (MY_PI * (cube((bin + 1) * bin_width) - cube(bin * bin_width)));

  neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputeStressSpherical::init_list(int /* id */, NeighList *ptr)
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

void ComputeStressSpherical::compute_array()
{
  invoked_array = update->ntimestep;

  int bin;
  // Zero arrays
  for (int bin = 0; bin < nbins; bin++) {
    tdens[bin] = 0.0;
    tpkrr[bin] = 0.0;
    tpktt[bin] = 0.0;
    tpkpp[bin] = 0.0;
    tpcrr[bin] = 0.0;
    tpctt[bin] = 0.0;
    tpcpp[bin] = 0.0;
  }

  int i, j, ii, jj, inum, jnum, itype, jtype;
  tagint itag, jtag;
  double ri[3], xtmp, ytmp, ztmp, tmp;
  double rsq, fpair, factor_coul, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double r, vr, vt, vp, theta;
  double **x = atom->x;
  double **v = atom->v;
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

  // calculate number density and kinetic contribution to pressure
  for (i = 0; i < nlocal; i++) {
    ri[0] = x[i][0] - x0;
    ri[1] = x[i][1] - y0;
    ri[2] = x[i][2] - z0;
    for (j = 0; j < 3; j++) {
      tmp = domain->boxhi[j] - domain->boxlo[j];
      if (ri[j] > 0.5 * tmp)
        ri[j] -= tmp;
      else if (ri[j] < -0.5 * tmp)
        ri[j] += tmp;
    }
    r = sqrt(ri[0] * ri[0] + ri[1] * ri[1] + ri[2] * ri[2]);
    if (r >= Rmax) continue;
    bin = (int) (r / bin_width);

    // Avoiding division by zero
    if ((r != 0.0) && (ri[0] != 0.0)) {
      theta = acos(ri[2] / r);
      tdens[bin] += 1.0;
      vr = (ri[0] * v[i][0] + ri[1] * v[i][1] + ri[2] * v[i][2]) / r;
      vt = r * sin(theta) / (square(ri[1] / ri[0]) + 1.0) *
          ((ri[0] * v[i][1] - ri[1] * v[i][0]) / (ri[0] * ri[0]));
      vp = (ri[2] * vr - r * v[i][2]) / (r * sqrt(1.0 - square(ri[2] / r)));
      tpkrr[bin] += vr * vr;
      tpktt[bin] += vt * vt;
      tpkpp[bin] += vp * vp;
    }
  }

  // loop over neighbors of my atoms
  // skip if I or J are not in group
  // for newton = 0 and J = ghost atom,
  //   need to ensure I,J pair is only output by one proc
  //   use same itag,jtag logic as in Neighbor::neigh_half_nsq()
  // for flag = 0, just count pair interactions within force cutoff
  // for flag = 1, calculate requested output fields

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  double qi[3], l1, l2, l3, l4, R1, R2, Fa, Fb, l_sum;
  double rij, f, ririj, sqr, la, lb, sql0, lambda0;
  double rsqxy, ririjxy, sqrixy, sqlxy0, A, B, C;
  int end_bin;

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

      ri[0] = x[j][0] - xtmp;
      ri[1] = x[j][1] - ytmp;
      ri[2] = x[j][2] - ztmp;

      rsq = ri[0] * ri[0] + ri[1] * ri[1] + ri[2] * ri[2];
      jtype = type[j];

      // Check if inside cut-off
      if (rsq >= cutsq[itype][jtype]) continue;

      pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);
      qi[0] = xtmp - x0;
      qi[1] = ytmp - y0;
      qi[2] = ztmp - z0;
      for (int k = 0; k < 3; k++) {
        tmp = domain->boxhi[k] - domain->boxlo[k];
        if (qi[k] > 0.5 * tmp)
          qi[k] -= tmp;
        else if (qi[k] < -0.5 * tmp)
          qi[k] += tmp;
      }
      l_sum = 0.0;
      rij = sqrt(rsq);
      f = -rij * fpair;
      ririj = qi[0] * ri[0] + qi[1] * ri[1] + qi[2] * ri[2];
      sqr = qi[0] * qi[0] + qi[1] * qi[1] + qi[2] * qi[2];
      la = 0.0;
      lb = 0.0;
      sql0 = sqr - ririj * ririj / rsq;
      lambda0 = -ririj / rsq;
      rsqxy = ri[0] * ri[0] + ri[1] * ri[1];
      ririjxy = qi[0] * ri[0] + qi[1] * ri[1];
      sqrixy = qi[0] * qi[0] + qi[1] * qi[1];
      sqlxy0 = (rsqxy != 0.0) ? sqrixy - ririjxy * ririjxy / rsqxy : 0.0;
      A = square(qi[0] * ri[1] - qi[1] * ri[0]);
      if (sqlxy0 > SMALL) C = sqrt(rsqxy * sqrixy - ririjxy * ririjxy);
      if (sql0 > SMALL) B = sqrt(rsq * sqr - ririj * ririj);
      end_bin = (int) (sqrt(rsq + 2.0 * ririj + sqr) / bin_width);

      while (lb < 1.0) {
        l1 = lb + SMALL;
        bin = (int) floor(sqrt(rsq * l1 * l1 + 2.0 * ririj * l1 + sqr) / bin_width);
        R1 = bin * bin_width;
        R2 = (bin + 1) * bin_width;

        // we must not take the square root of a negative number or divide by zero.
        if ((R1 * R1 < sql0) || (R2 * R2 < sql0) || (rsq <= 0.0)) {
          lb = 1.0;
        } else {
          l1 = lambda0 + sqrt((R1 * R1 - sql0) / rsq);
          l2 = lambda0 - sqrt((R1 * R1 - sql0) / rsq);
          l3 = lambda0 + sqrt((R2 * R2 - sql0) / rsq);
          l4 = lambda0 - sqrt((R2 * R2 - sql0) / rsq);

          if (l4 >= 0.0 && l4 <= 1.0 && l4 > la)
            lb = l4;
          else if (l3 >= 0.0 && l3 <= 1.0 && l3 > la)
            lb = l3;
          else if (l2 >= 0.0 && l2 <= 1.0 && l2 > la)
            lb = l2;
          else if (l1 >= 0.0 && l1 <= 1.0 && l1 > la)
            lb = l1;
          else
            lb = 1.0;
        }

        if (bin == end_bin) lb = 1.0;
        if (la > lb) error->all(FLERR, "Error: la > lb\n");
        if (bin >= 0 && bin < nbins) {
          if (sql0 > SMALL) {
            Fa = -B * atan2(rsq * la + ririj, B);
            Fb = -B * atan2(rsq * lb + ririj, B);
            tpcrr[bin] -= (f / rij) * (rsq * (lb - la) + Fb - Fa);
            tpcpp[bin] -= (f / rij) * (Fa - Fb) / 2.0;
          } else
            tpcrr[bin] -= f * rij * (lb - la);

          if (sqlxy0 > SMALL)
            tpctt[bin] -= (f / rij) * A / C *
                (atan2(rsqxy * lb + ririjxy, C) - atan2(rsqxy * la + ririjxy, C));
        }
        l_sum += lb - la;
        la = lb;
      }

      // Error check
      if (fabs(l_sum - 1.0) > SMALL)
        error->all(FLERR, "ERROR: The sum of the fractional line segments is not 1.0");
    }
  }

  // normalize pressure
  for (bin = 0; bin < nbins; bin++) {
    tdens[bin] *= invV[bin];
    tpkrr[bin] *= invV[bin];
    tpktt[bin] *= invV[bin];
    tpkpp[bin] *= invV[bin];
    tpcrr[bin] *= invV[bin];
    tpctt[bin] *= invV[bin];
    tpcpp[bin] *= invV[bin];
  }
  // communicate across processors
  MPI_Allreduce(tdens, dens, nbins, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(tpkrr, pkrr, nbins, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(tpktt, pktt, nbins, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(tpkpp, pkpp, nbins, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(tpcrr, pcrr, nbins, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(tpctt, pctt, nbins, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(tpcpp, pcpp, nbins, MPI_DOUBLE, MPI_SUM, world);

  // populate array to output.
  for (bin = 0; bin < nbins; bin++) {
    array[bin][0] = (bin + 0.5) * bin_width;
    array[bin][1] = dens[bin];
    array[bin][2] = pkrr[bin];
    array[bin][3] = pktt[bin];
    array[bin][4] = pkpp[bin];
    array[bin][5] = pcrr[bin];
    array[bin][6] = pctt[bin];
    array[bin][7] = pcpp[bin];
  }
}

double ComputeStressSpherical::memory_usage()
{
  return 15.0 * (double) (nbins + size_array_rows * size_array_cols) * sizeof(double);
}
