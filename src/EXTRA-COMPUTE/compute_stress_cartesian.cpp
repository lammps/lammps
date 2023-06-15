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

#include "compute_stress_cartesian.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define SMALL 1.0e-10
/*-----------------------------------------------------------------------------------
  Contributing author: Olav Galteland (Norwegian University of Science and Technology)
                        olav.galteland@ntnu.no
------------------------------------------------------------------------------------*/

static const char cite_compute_stress_cartesian[] =
    "compute stress/cartesian: doi:10.3390/nano11010165\n\n"
    "@article{galteland2021nanothermodynamic,\n"
    "title={Nanothermodynamic Description and Molecular Simulation of a\n"
    "   Single-Phase Fluid in a Slit Pore},\n"
    "author={Galteland, Olav and Bedeaux, Dick and Kjelstrup, Signe},\n"
    "journal={Nanomaterials},\n"
    "volume={11},\n"
    "number={1},\n"
    "pages={165},\n"
    "year={2021},\n"
    "publisher={Multidisciplinary Digital Publishing Institute}\n"
    "}\n\n";

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

ComputeStressCartesian::ComputeStressCartesian(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), dens(nullptr), pkxx(nullptr), pkyy(nullptr), pkzz(nullptr),
    pcxx(nullptr), pcyy(nullptr), pczz(nullptr), tdens(nullptr), tpkxx(nullptr), tpkyy(nullptr),
    tpkzz(nullptr), tpcxx(nullptr), tpcyy(nullptr), tpczz(nullptr), list(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_compute_stress_cartesian);

  // narg == 5 for one-dimensional and narg == 7 for two-dimensional
  if (narg == 5)
    dims = 1;
  else if (narg == 7)
    dims = 2;
  else
    error->all(FLERR, "Illegal compute stress/cartesian command. Illegal number of arguments.");

  if (strcmp(arg[3], "x") == 0)
    dir1 = 0;
  else if (strcmp(arg[3], "y") == 0)
    dir1 = 1;
  else if (strcmp(arg[3], "z") == 0)
    dir1 = 2;
  else
    error->all(FLERR, "Illegal compute stress/cartesian direction: {}", arg[3]);

  dir2 = 0;
  bin_width1 = utils::numeric(FLERR, arg[4], false, lmp);
  bin_width2 = domain->boxhi[dir2] - domain->boxlo[dir2];
  nbins1 = (int) ((domain->boxhi[dir1] - domain->boxlo[dir1]) / bin_width1);
  nbins2 = 1;

  // no triclinic boxes
  if (domain->triclinic) error->all(FLERR, "Compute stress/cartesian requires an orthogonal box");

  // adjust bin width if not a perfect match
  double tmp_binwidth = (domain->boxhi[dir1] - domain->boxlo[dir1]) / nbins1;
  if ((fabs(tmp_binwidth - bin_width1) > SMALL) && (comm->me == 0))
    utils::logmesg(lmp, "Adjusting first bin width for compute {} from {:.6f} to {:.6f}\n", style,
                   bin_width1, tmp_binwidth);
  bin_width1 = tmp_binwidth;

  if (bin_width1 <= 0.0)
    error->all(FLERR, "Illegal compute stress/cartesian command. First bin width must be > 0");
  else if (bin_width1 > domain->boxhi[dir1] - domain->boxlo[dir1])
    error->all(FLERR, "Illegal compute stress/cartesian command. First bin width > box.");

  invV = bin_width1;
  if (dims == 2) {
    if (strcmp(arg[5], "x") == 0)
      dir2 = 0;
    else if (strcmp(arg[5], "y") == 0)
      dir2 = 1;
    else if (strcmp(arg[5], "z") == 0)
      dir2 = 2;
    else
      error->all(FLERR, "Illegal compute stress/cartesian direction {}", arg[5]);

    bin_width2 = utils::numeric(FLERR, arg[6], false, lmp);
    nbins2 = (int) ((domain->boxhi[dir2] - domain->boxlo[dir2]) / bin_width2);

    // adjust bin width if not a perfect match
    tmp_binwidth = (domain->boxhi[dir2] - domain->boxlo[dir2]) / nbins2;
    if ((fabs(tmp_binwidth - bin_width2) > SMALL) && (comm->me == 0))
      utils::logmesg(lmp, "Adjusting second bin width for compute {} from {:.6f} to {:.6f}\n",
                     style, bin_width2, tmp_binwidth);
    bin_width2 = tmp_binwidth;

    invV *= bin_width2;

    if (bin_width2 <= 0.0)
      error->all(FLERR, "Illegal compute stress/cartesian command. Second bin width must be > 0");
    else if (bin_width2 > domain->boxhi[dir2] - domain->boxlo[dir2])
      error->all(FLERR, "Illegal compute stress/cartesian command. Second bin width > box");
  }

  // check for overflow
  if (((bigint) nbins1 * nbins2) > MAXSMALLINT)
    error->all(FLERR, "Too many bins in compute stress/cartesian");

  // check for variable box dimension
  int box_incompatible = 0;
  for (auto ifix : modify->get_fix_list()) {
    if (((dir1 == 0) && (ifix->box_change & Fix::BOX_CHANGE_X)) ||
        ((dir1 == 1) && (ifix->box_change & Fix::BOX_CHANGE_Y)) ||
        ((dir1 == 2) && (ifix->box_change & Fix::BOX_CHANGE_Z)))
      box_incompatible = 1;
    if (dims == 2) {
      if (((dir2 == 0) && (ifix->box_change & Fix::BOX_CHANGE_X)) ||
          ((dir2 == 1) && (ifix->box_change & Fix::BOX_CHANGE_Y)) ||
          ((dir2 == 2) && (ifix->box_change & Fix::BOX_CHANGE_Z)))
        box_incompatible = 1;
    }
  }
  if (box_incompatible)
    error->all(FLERR, "Must not use compute stress/cartesian on variable box dimension");

  for (int i = 0; i < 3; i++)
    if ((dims == 1 && i != dir1) || (dims == 2 && (i != dir1 && i != dir2)))
      invV *= domain->boxhi[i] - domain->boxlo[i];
  invV = 1.0 / invV;

  array_flag = 1;
  vector_flag = 0;
  extarray = 0;
  size_array_cols = 7 + dims;    // dir1, dir2, number density, pkxx, pkyy, pkzz, pcxx, pcyy, pczz
  size_array_rows = nbins1 * nbins2;

  memory->create(dens, nbins1 * nbins2, "dens");
  memory->create(pkxx, nbins1 * nbins2, "pkxx");
  memory->create(pkyy, nbins1 * nbins2, "pkyy");
  memory->create(pkzz, nbins1 * nbins2, "pkzz");
  memory->create(pcxx, nbins1 * nbins2, "pcxx");
  memory->create(pcyy, nbins1 * nbins2, "pcyy");
  memory->create(pczz, nbins1 * nbins2, "pczz");
  memory->create(tdens, nbins1 * nbins2, "tdens");
  memory->create(tpkxx, nbins1 * nbins2, "tpkxx");
  memory->create(tpkyy, nbins1 * nbins2, "tpkyy");
  memory->create(tpkzz, nbins1 * nbins2, "tpkzz");
  memory->create(tpcxx, nbins1 * nbins2, "tpcxx");
  memory->create(tpcyy, nbins1 * nbins2, "tpcyy");
  memory->create(tpczz, nbins1 * nbins2, "tpczz");
  memory->create(array, size_array_rows, size_array_cols, "stress:cartesian:output");
}

/* ---------------------------------------------------------------------- */

ComputeStressCartesian::~ComputeStressCartesian()
{
  memory->destroy(dens);
  memory->destroy(pkxx);
  memory->destroy(pkyy);
  memory->destroy(pkzz);
  memory->destroy(pcxx);
  memory->destroy(pcyy);
  memory->destroy(pczz);
  memory->destroy(tdens);
  memory->destroy(tpkxx);
  memory->destroy(tpkyy);
  memory->destroy(tpkzz);
  memory->destroy(tpcxx);
  memory->destroy(tpcyy);
  memory->destroy(tpczz);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeStressCartesian::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "No pair style is defined for compute stress/cartesian");
  if (force->pair->single_enable == 0)
    error->all(FLERR, "Pair style does not support compute stress/cartesian");

  // need an occasional half neighbor list.
  neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputeStressCartesian::init_list(int /* id */, NeighList *ptr)
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

void ComputeStressCartesian::compute_array()
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  int bin, bin1, bin2;
  tagint itag, jtag;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, fpair, factor_coul, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

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
  double *boxlo = domain->boxlo;

  // invoke half neighbor list (will copy or build if necessary)
  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Zero arrays
  for (bin = 0; bin < nbins1 * nbins2; bin++) {
    tdens[bin] = 0;
    tpkxx[bin] = 0;
    tpkyy[bin] = 0;
    tpkzz[bin] = 0;
    tpcxx[bin] = 0;
    tpcyy[bin] = 0;
    tpczz[bin] = 0;
  }

  // calculate number density and kinetic contribution to pressure
  for (i = 0; i < nlocal; i++) {
    bin1 = (int) ((x[i][dir1] - boxlo[dir1]) / bin_width1) % nbins1;
    bin2 = 0;
    if (dims == 2) bin2 = (int) ((x[i][dir2] - boxlo[dir2]) / bin_width2) % nbins2;

    // Apply periodic boundary conditions and avoid out of range access
    if (domain->periodicity[dir1] == 1) {
      if (bin1 < 0)
        bin1 = (bin1 + nbins1) % nbins1;
      else if (bin1 >= nbins1)
        bin1 = (bin1 - nbins1) % nbins1;
    } else if (bin1 < 0)
      bin1 = 0;
    else if (bin1 >= nbins1)
      bin1 = nbins1 - 1;

    if (domain->periodicity[dir2] == 1) {
      if (bin2 < 0)
        bin2 = (bin2 + nbins2) % nbins2;
      else if (bin2 >= nbins2)
        bin2 = (bin2 - nbins2) % nbins2;
    } else if (bin2 < 0)
      bin2 = 0;
    else if (bin2 >= nbins2)
      bin2 = nbins2 - 1;

    j = bin1 + bin2 * nbins1;
    tdens[j] += 1;
    tpkxx[j] += mass[type[i]] * v[i][0] * v[i][0];
    tpkyy[j] += mass[type[i]] * v[i][1] * v[i][1];
    tpkzz[j] += mass[type[i]] * v[i][2] * v[i][2];
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

  double xi1, xi2;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    xi1 = x[i][dir1] - boxlo[dir1];
    xi2 = x[i][dir2] - boxlo[dir2];
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
      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;

      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      // Check if inside cut-off
      if (rsq >= cutsq[itype][jtype]) continue;
      pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);
      compute_pressure(fpair, xi1, xi2, delx, dely, delz);
    }
  }

  // normalize pressure
  for (bin = 0; bin < nbins1 * nbins2; bin++) {
    tdens[bin] *= invV;
    tpkxx[bin] *= invV;
    tpkyy[bin] *= invV;
    tpkzz[bin] *= invV;
    tpcxx[bin] *= invV;
    tpcyy[bin] *= invV;
    tpczz[bin] *= invV;
  }

  // communicate across processors
  MPI_Allreduce(tdens, dens, nbins1 * nbins2, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(tpkxx, pkxx, nbins1 * nbins2, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(tpkyy, pkyy, nbins1 * nbins2, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(tpkzz, pkzz, nbins1 * nbins2, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(tpcxx, pcxx, nbins1 * nbins2, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(tpcyy, pcyy, nbins1 * nbins2, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(tpczz, pczz, nbins1 * nbins2, MPI_DOUBLE, MPI_SUM, world);

  // populate array to output.
  for (bin = 0; bin < nbins1 * nbins2; bin++) {
    array[bin][0] = (bin % nbins1 + 0.5) * bin_width1;
    if (dims == 2) array[bin][1] = ((int) (bin / nbins1) + 0.5) * bin_width2;
    array[bin][0 + dims] = dens[bin];
    array[bin][1 + dims] = pkxx[bin];
    array[bin][2 + dims] = pkyy[bin];
    array[bin][3 + dims] = pkzz[bin];
    array[bin][4 + dims] = pcxx[bin];
    array[bin][5 + dims] = pcyy[bin];
    array[bin][6 + dims] = pczz[bin];
  }
}

void ComputeStressCartesian::compute_pressure(double fpair, double xi, double yi, double delx,
                                              double dely, double delz)
{
  int bin1, bin2, next_bin1, next_bin2;
  double la = 0.0, lb = 0.0, l_sum = 0.0;
  double rij[3] = {delx, dely, delz};
  double l1 = 0.0, l2, rij1, rij2;
  rij1 = rij[dir1];
  rij2 = rij[dir2];

  next_bin1 = (int) floor(xi / bin_width1);
  next_bin2 = (int) floor(yi / bin_width2);

  // Integrating along line
  while (lb < 1.0) {
    bin1 = next_bin1;
    bin2 = next_bin2;

    if (rij1 > 0)
      l1 = ((bin1 + 1) * bin_width1 - xi) / rij1;
    else
      l1 = (bin1 * bin_width1 - xi) / rij1;
    if (rij2 > 0)
      l2 = ((bin2 + 1) * bin_width2 - yi) / rij2;
    else
      l2 = (bin2 * bin_width2 - yi) / rij2;

    if ((l1 < l2 || l2 < lb + SMALL) && l1 <= 1.0 && l1 > lb) {
      lb = l1;
      next_bin1 = bin1 + (int) (rij1 / fabs(rij1));
    } else if (l2 <= 1.0 && l2 > lb) {
      lb = l2;
      next_bin2 = bin2 + (int) (rij2 / fabs(rij2));
    } else
      lb = 1.0;

    // Periodic boundary conditions
    if (domain->periodicity[dir1] == 1) {
      if (bin1 < 0)
        bin1 = (bin1 + nbins1) % nbins1;
      else if (bin1 >= nbins1)
        bin1 = (bin1 - nbins1) % nbins1;
    } else if (bin1 < 0)
      bin1 = 0;
    else if (bin1 >= nbins1)
      bin1 = nbins1 - 1;

    if (domain->periodicity[dir2] == 1) {
      if (bin2 < 0)
        bin2 = (bin2 + nbins2) % nbins2;
      else if (bin2 >= nbins2)
        bin2 = (bin2 - nbins2) % nbins2;
    } else if (bin2 < 0)
      bin2 = 0;
    else if (bin2 >= nbins2)
      bin2 = nbins2 - 1;

    if (bin1 + bin2 * nbins1 > nbins1 * nbins2) error->all(FLERR, "Bin outside: lb={:.16g}", lb);

    tpcxx[bin1 + bin2 * nbins1] += (fpair * delx * delx * (lb - la));
    tpcyy[bin1 + bin2 * nbins1] += (fpair * dely * dely * (lb - la));
    tpczz[bin1 + bin2 * nbins1] += (fpair * delz * delz * (lb - la));

    l_sum += lb - la;
    la = lb;
  }

  if (l_sum > 1.0 + SMALL || l_sum < 1.0 - SMALL)
    error->all(FLERR, "Sum of fractional line segments does not equal 1.");
}

/* ----------------------------------------------------------------------
memory usage of data
------------------------------------------------------------------------- */

double ComputeStressCartesian::memory_usage()
{
  return (14.0 + dims + 7) * (double) (nbins1 * nbins2) * sizeof(double);
}
