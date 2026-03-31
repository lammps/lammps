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
   Contributing author: David Immel (d.immel@fz-juelich.de, FZJ, Germany)
------------------------------------------------------------------------- */

#include "fix_lambda_la_csp_apip.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_store_atom.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <algorithm>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixLambdaLACSPAPIP::FixLambdaLACSPAPIP(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), ngh_pairs(nullptr), f_lambda(nullptr), distsq(nullptr), nearest(nullptr),
    list(nullptr), fixstore_pairs(nullptr), fixstore_la_avg(nullptr), fixstore_la_inp(nullptr),
    fixstore_la_norm(nullptr), prefactor1(nullptr), prefactor2(nullptr)
{
  comm_reverse = 2;
  comm_forward = 2;
  comm_forward_flag = FORWARD_INP_LAMBDA;
  restart_global = 1;

  prefactor1_size = prefactor2_size = ngh_pairs_size = 0;
  maxneigh = 0;

  store_stats = false;

  // set defaults
  cut_lo = cut_hi = -1;
  threshold_lo = threshold_hi = -1;
  lambda_non_group = 1;    // fast
  const_ngh_flag = true;
  calculate_forces_flag = true;
  csp_cutsq = 25;

  tags_stored = false;
  counter_changed_csp_nghs = 0;

  if (narg < 8) error->all(FLERR, "fix lambda/la/csp/apip requires five arguments");

  threshold_lo = utils::numeric(FLERR, arg[3], false, lmp);
  threshold_hi = utils::numeric(FLERR, arg[4], false, lmp);
  threshold_width = threshold_hi - threshold_lo;

  cut_lo = utils::numeric(FLERR, arg[5], false, lmp);
  cut_hi = utils::numeric(FLERR, arg[6], false, lmp);
  cut_width = cut_hi - cut_lo;
  cut_hi_sq = cut_hi * cut_hi;

  if (strcmp(arg[7], "fcc") == 0)
    nnn = 12;
  else if (strcmp(arg[7], "bcc") == 0)
    nnn = 14;
  else
    nnn = utils::inumeric(FLERR, arg[7], false, lmp);

  for (int i = 8; i < narg; i++) {
    if (strcmp(arg[i], "csp_mode") == 0) {
      if (i + 1 == narg)
        error->all(FLERR,
                   "the csp_mode option of fix lambda/la/csp/apip requires an additional argument");
      if (strcmp(arg[i + 1], "dynamic") == 0)
        const_ngh_flag = false;
      else if (strcmp(arg[i + 1], "static") == 0)
        const_ngh_flag = true;
      else
        error->all(FLERR, "expected dynamic or static instead of {}", arg[i + 1]);
      i++;
    } else if (strcmp(arg[i], "csp_cut") == 0) {
      if (i + 1 == narg)
        error->all(FLERR,
                   "the csp_cut option of fix lambda/la/csp/apip requires an additional argument");
      csp_cutsq = pow(utils::numeric(FLERR, arg[i + 1], false, lmp), 2);
      i++;
    } else if (strcmp(arg[i], "forces") == 0) {
      if (i + 1 == narg)
        error->all(FLERR,
                   "the forces option of fix lambda/la/csp/apip requires an additional argument");
      calculate_forces_flag = utils::logical(FLERR, arg[i + 1], false, lmp);
      i++;
    } else if (strcmp(arg[i], "lambda_non_group") == 0) {
      if (i + 1 == narg)
        error->all(FLERR,
                   "the lambda_non_group option of fix lambda/la/csp/apip requires an additional "
                   "argument");
      if (strcmp(arg[i + 1], "fast") == 0)
        lambda_non_group = 1;
      else if (strcmp(arg[i + 1], "precise") == 0)
        lambda_non_group = 0;
      else
        lambda_non_group = utils::numeric(FLERR, arg[i + 1], false, lmp);
      i++;
    } else if (strcmp(arg[i], "store_peratom") == 0) {
      if (i + 1 == narg)
        error->all(
            FLERR,
            "the store_stats option of fix lambda/la/csp/apip requires an additional argument");
      peratom_freq = utils::inumeric(FLERR, arg[i + 1], false, lmp);
      if (peratom_freq <= 0) error->all(FLERR, "store_peratom frequency needs to be positive.");
      i++;
    } else {
      error->all(FLERR, "unknown argument {}", arg[i]);
    }
  }

  // verify arguments
  if (cut_lo > cut_hi || cut_lo < 0) error->all(FLERR, "0 <= cut_lo <= cut_hi required");
  if (threshold_lo > threshold_hi || threshold_lo < 0)
    error->all(FLERR, "0 <= threshold_lo <= threshold_hi required");
  if (lambda_non_group < 0 || lambda_non_group > 1)
    error->all(FLERR, "0 <= lambda_non_group <= 1 required");
  if (!const_ngh_flag) {
    scalar_flag = 1;
    extscalar = 1;
  }

  cutsq_combined = csp_cutsq > cut_hi_sq ? csp_cutsq : cut_hi_sq;

  if (calculate_forces_flag) { virial_global_flag = virial_peratom_flag = thermo_virial = 1; }

  if (!atom->apip_lambda_flag) {
    error->all(FLERR, "fix lambda/la/csp/apip requires atomic style with lambda.");
  }

  if (peratom_freq > 0) {
    store_stats = true;
    peratom_flag = 1;
    size_peratom_cols = 5;

    size_f_lambda = atom->nlocal;
    memory->create(f_lambda, size_f_lambda, size_peratom_cols, "pair:lambda:la:csp:apip:f:lambda");
    array_atom = f_lambda;
    for (int i = 0; i < atom->nlocal; i++) {
      for (int j = 0; j < size_peratom_cols; j++) f_lambda[i][j] = 0;
    }
  }
}

/* ---------------------------------------------------------------------- */

FixLambdaLACSPAPIP::~FixLambdaLACSPAPIP()
{
  memory->destroy(ngh_pairs);
  if (store_stats) memory->destroy(f_lambda);
  memory->destroy(distsq);
  memory->destroy(nearest);
  memory->destroy(prefactor1);
  memory->destroy(prefactor2);
  if (fixstore_pairs && modify->nfix) modify->delete_fix(fixstore_pairs->id);
  if (fixstore_la_avg && modify->nfix) modify->delete_fix(fixstore_la_avg->id);
  if (fixstore_la_inp && modify->nfix) modify->delete_fix(fixstore_la_inp->id);
  if (fixstore_la_norm && modify->nfix) modify->delete_fix(fixstore_la_norm->id);
  fixstore_pairs = fixstore_la_avg = fixstore_la_inp = fixstore_la_norm = nullptr;
}

/**
  * allocate per-particle storage for the last ngh_pairs
  * fix could already be allocated if fix lambda/la/csp/apip is re-specified
  */

void FixLambdaLACSPAPIP::post_constructor()
{
  // 1. CSP pairs
  std::string cmd;
  cmd = id;
  cmd += "LAST_CSP_NGH_PAIR";

  // delete existing fix store if existing
  fixstore_pairs = dynamic_cast<FixStoreAtom *>(modify->get_fix_by_id(cmd));
  // check nfix in case all fixes have already been deleted
  if (fixstore_pairs && modify->nfix) modify->delete_fix(fixstore_pairs->id);
  fixstore_pairs = nullptr;

  char str_values[40];
  sprintf(str_values, "%d", nnn);

  // arguments of peratom:
  cmd += " all STORE/ATOM ";
  cmd += str_values;    // n1
  cmd += " 0 0 1";      // n2 gflag rflag
  fixstore_pairs = dynamic_cast<FixStoreAtom *>(modify->add_fix(cmd));

  // do not carry the CSP-pairs with atoms during normal atom migration yet
  // activate after the CSP-pairs are calculated
  fixstore_pairs->disable = 1;

  // 2. local averaging data
  cmd = id;
  cmd += "LA_INP";
  fixstore_la_inp = dynamic_cast<FixStoreAtom *>(modify->get_fix_by_id(cmd));
  if (fixstore_la_inp && modify->nfix) modify->delete_fix(fixstore_la_inp->id);
  fixstore_la_inp = nullptr;
  cmd += " all STORE/ATOM 1 0 1 0";
  fixstore_la_inp = dynamic_cast<FixStoreAtom *>(modify->add_fix(cmd));

  cmd = id;
  cmd += "LA_NORM";
  fixstore_la_norm = dynamic_cast<FixStoreAtom *>(modify->get_fix_by_id(cmd));
  if (fixstore_la_norm && modify->nfix) modify->delete_fix(fixstore_la_norm->id);
  fixstore_la_norm = nullptr;
  cmd += " all STORE/ATOM 1 0 1 0";
  fixstore_la_norm = dynamic_cast<FixStoreAtom *>(modify->add_fix(cmd));

  cmd = id;
  cmd += "LA_AVG";
  fixstore_la_avg = dynamic_cast<FixStoreAtom *>(modify->get_fix_by_id(cmd));
  if (fixstore_la_avg && modify->nfix) modify->delete_fix(fixstore_la_avg->id);
  fixstore_la_avg = nullptr;
  cmd += " all STORE/ATOM 1 0 1 0";
  fixstore_la_avg = dynamic_cast<FixStoreAtom *>(modify->add_fix(cmd));
}

/* ---------------------------------------------------------------------- */

int FixLambdaLACSPAPIP::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_NEIGHBOR;
  mask |= PRE_REVERSE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLambdaLACSPAPIP::init()
{
  // full neighbour list for thermostating
  auto *req = neighbor->add_request(this, NeighConst::REQ_FULL);
  req->set_cutoff(sqrt(cutsq_combined));

  if (atom->tag_enable == 0) error->all(FLERR, "fix lambda/la/csp/apip requires atom IDs");

  // only one fix lambda/la/csp/apip
  int count = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style, "lambda/la/csp/apip") == 0) count++;
  }
  if (count > 1) error->all(FLERR, "More than one fix lambda/la/csp/apip.");

  if (force->pair->cutforce < cut_hi)
    error->all(FLERR, "cutoff of potential ({}) smaller than cutoff of weighting function ({})",
               force->pair->cutforce, cut_hi);
  if (force->pair->cutforce * force->pair->cutforce < csp_cutsq)
    error->all(FLERR, "cutoff of potential ({}) smaller than cutoff of the CSP ({})",
               force->pair->cutforce, sqrt(csp_cutsq));

  if (strcmp(atom->atom_style, "apip"))
    error->all(FLERR, "fix lambda/la/csp/apip requires atom style apip");
}

/**
  *  Init neighbor list.
  */

void FixLambdaLACSPAPIP::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/**
  *  Build atom map if required.
  */

void FixLambdaLACSPAPIP::setup_post_neighbor()
{
  if (const_ngh_flag && atom->map_style == Atom::MAP_NONE) {
    atom->map_init();
    atom->map_set();
  }
}

/**
  *  Rebuild atom map if required.
  */

void FixLambdaLACSPAPIP::post_neighbor()
{
  if (atom->map_style != Atom::MAP_NONE) {
    atom->map_init();
    atom->map_set();
  }
}

/**
  *  Compute lambda, csp, csp_avg, csp_norm for all local atoms.
  */

void FixLambdaLACSPAPIP::setup_pre_force(int /*vflag*/)
{
  if (!const_ngh_flag || !tags_stored)
    pre_force_dyn_pairs();
  else
    pre_force_const_pairs();

  comm_forward_flag = FORWARD_INP_LAMBDA;
  comm->forward_comm(this);
}

void FixLambdaLACSPAPIP::pre_force(int /*vflag*/)
{
  if (const_ngh_flag)
    pre_force_const_pairs();
  else
    pre_force_dyn_pairs();

  comm_forward_flag = FORWARD_INP_LAMBDA;
  comm->forward_comm(this);
}

/**
  *  Compute lambda, csp, csp_avg, csp_norm for all local atoms with dynamic neighbour pairs.
  */

void FixLambdaLACSPAPIP::pre_force_dyn_pairs()
{
  int i, j, k, ii, jj, kk, n_nearest, n_pairs, inum, jnum;
  int *ilist, *jlist, *numneigh, **firstneigh, *mask;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, value, weight;
  double *lambda, *csp, *csp_avg, *csp_norm;

  tagint *tag = atom->tag;

  mask = atom->mask;
  lambda = atom->apip_lambda;
  csp = fixstore_la_inp->vstore;
  csp_avg = fixstore_la_avg->vstore;
  csp_norm = fixstore_la_norm->vstore;

  int nlocal = atom->nlocal;
  int nmax = atom->nmax;

  // grow ngh_pairs array if necessary

  if (atom->nlocal > ngh_pairs_size) {
    memory->destroy(ngh_pairs);
    ngh_pairs_size = atom->nlocal;
    memory->create(ngh_pairs, ngh_pairs_size, nnn, "fix_lambda_la_csp_apip:ngh_pairs");
  }

  double **stored_tags = fixstore_pairs->astore;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // npairs = number of unique pairs

  int nhalf = nnn / 2;
  int npairs = nnn * (nnn - 1) / 2;
  auto *pairs_value = new double[npairs];
  auto *pairs_j = new int[npairs];
  auto *pairs_k = new int[npairs];
  auto *pairs_index = new int[npairs];
  auto *pairs_used_now = new CSPpairAPIP[nhalf];
  auto *pairs_used_prev = new CSPpairAPIP[nhalf];

  // compute centro-symmetry parameter for each atom in group

  double **x = atom->x;

  // zero values
  for (i = 0; i < nlocal; i++) csp[i] = 0;
  for (i = 0; i < nmax; i++) csp_norm[i] = csp_avg[i] = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    // 1. compute csp[i]

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // ensure distsq and nearest arrays are long enough

    if (jnum > maxneigh) {
      memory->destroy(distsq);
      memory->destroy(nearest);
      maxneigh = jnum;
      memory->create(distsq, maxneigh, "lambda:la:csp:apip:distsq");
      memory->create(nearest, maxneigh, "lambda:la:csp:apip:nearest");
    }

    // loop over list of all neighbors within force cutoff
    // distsq[] = distance sq to each
    // nearest[] = atom indices of neighbors

    n_nearest = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < cutsq_combined) {
        distsq[n_nearest] = rsq;
        nearest[n_nearest++] = j;
      }
    }

    if (n_nearest < nnn) {
      error->all(FLERR, "Atom has only {} neighbours, but {} are required to calculate the CSP.",
                 n_nearest, nnn);
    }

    // store nnn nearest neighs in 1st nnn locations of distsq and nearest

    select2(nnn, n_nearest, distsq, nearest);

    // R = Ri + Rj for each of npairs i,j pairs among nnn neighbors
    // pairs = squared length of each R

    n_pairs = 0;
    for (j = 0; j < nnn; j++) {
      jj = nearest[j];
      for (k = j + 1; k < nnn; k++) {
        kk = nearest[k];
        delx = x[jj][0] + x[kk][0] - 2.0 * xtmp;
        dely = x[jj][1] + x[kk][1] - 2.0 * ytmp;
        delz = x[jj][2] + x[kk][2] - 2.0 * ztmp;
        pairs_value[n_pairs] = delx * delx + dely * dely + delz * delz;
        pairs_j[n_pairs] = jj;
        pairs_k[n_pairs] = kk;
        pairs_index[n_pairs++] = n_pairs;
      }
    }

    // store nhalf smallest pair distances in 1st nhalf locations of pairs

    select2(nhalf, npairs, pairs_value, pairs_index);

    // centrosymmetry = sum of nhalf smallest squared values
    // compute csp with detected neighbour pairs
    value = 0.0;
    for (j = 0; j < nhalf; j++) {
      value += pairs_value[j];
      // store used neighbor pairs
      ngh_pairs[i][j] = pairs_j[pairs_index[j]];
      ngh_pairs[i][j + nhalf] = pairs_k[pairs_index[j]];

      // get new tags
      pairs_used_now[j] = CSPpairAPIP(tag[ngh_pairs[i][j]], tag[ngh_pairs[i][j + nhalf]]);
    }

    csp[i] = value;

    // get previous tags
    if (tags_stored) {
      for (j = 0; j < nhalf; j++) {
        pairs_used_prev[j] =
            CSPpairAPIP((tagint) stored_tags[i][j], (tagint) stored_tags[i][j + nhalf]);
      }

      // sort tag arrays
      std::sort(pairs_used_now, pairs_used_now + nhalf);
      std::sort(pairs_used_prev, pairs_used_prev + nhalf);

      // count number of different pairs
      for (j = 0; j < nhalf; j++) {
        if (pairs_used_now[j] != pairs_used_prev[j]) {
          counter_changed_csp_nghs++;
          break;
        }
      }
    }

    // store current tags
    for (j = 0; j < nnn; j++) { stored_tags[i][j] = tag[ngh_pairs[i][j]]; }

    // 2. compute csp avg
    // csp_avg_i += csp_i * w(r_ij)
    // csp_norm_i += w(r_ij)

    // own particle
    weight = weighting_function_poly(0);
    csp_avg[i] += value * weight;
    csp_norm[i] += weight;

    // neighbouring particles
    for (j = 0; j < n_nearest; j++) {
      if (distsq[j] <= cut_hi_sq) {
        weight = weighting_function_poly(sqrt(distsq[j]));
        jj = nearest[j];
        csp_avg[jj] += value * weight;
        csp_norm[jj] += weight;
      }
    }
  }
  tags_stored = true;

  // carry the CSP-pairs with atoms during normal atom migration
  fixstore_pairs->disable = 0;

  // reverse communication of csp_norm and csp_avgs of ghost atoms
  comm->reverse_comm(this);

  for (i = 0; i < nlocal; i++) {
    csp_avg[i] /= csp_norm[i];

    if (mask[i] & groupbit)
      lambda[i] = switching_function_poly(csp_avg[i]);
    else
      lambda[i] = lambda_non_group;
  }

  delete[] pairs_value;
  delete[] pairs_index;
  delete[] pairs_j;
  delete[] pairs_k;
  delete[] pairs_used_now;
  delete[] pairs_used_prev;
}

/**
  *  Compute lambda, csp, csp_avg, csp_norm for all local atoms with constant neighbour pairs.
  */

void FixLambdaLACSPAPIP::pre_force_const_pairs()
{
  int i, j, ii, jj, kk, inum, jnum;
  int *ilist, *jlist, *numneigh, **firstneigh, *mask;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, value, weight, delx_j, dely_j, delz_j, delx_k,
      dely_k, delz_k;
  double *lambda, *csp, *csp_avg, *csp_norm, **x, **stored_tags;

  int nlocal = atom->nlocal;
  int nmax = atom->nmax;
  int nhalf = nnn / 2;

  mask = atom->mask;
  lambda = atom->apip_lambda;
  csp = fixstore_la_inp->vstore;
  csp_avg = fixstore_la_avg->vstore;
  csp_norm = fixstore_la_norm->vstore;
  x = atom->x;
  stored_tags = fixstore_pairs->astore;

  // grow ngh_pairs array if necessary

  if (atom->nlocal > ngh_pairs_size) {
    memory->destroy(ngh_pairs);
    ngh_pairs_size = atom->nlocal;
    memory->create(ngh_pairs, ngh_pairs_size, nnn, "fix_lambda_la_csp_apip:ngh_pairs");
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero values
  for (i = 0; i < nlocal; i++) csp[i] = 0;
  for (i = 0; i < nmax; i++) csp_norm[i] = csp_avg[i] = 0;

  // 1. compute csp[i]

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    value = 0;
    for (j = 0; j < nhalf; j++) {

      jj = atom->map((tagint) stored_tags[i][j]);
      if (jj == -1)
        error->one(FLERR, "atom ID {} csp neighbour with ID {} not present", atom->tag[i],
                   stored_tags[i][j]);
      while (jj >= 0) {
        // calculate distance to jj
        delx_j = x[jj][0] - xtmp;
        dely_j = x[jj][1] - ytmp;
        delz_j = x[jj][2] - ztmp;
        // correct periodic image?
        if (fabs(delx_j) < domain->prd_half[0] && fabs(dely_j) < domain->prd_half[1] &&
            fabs(delz_j) < domain->prd_half[2]) {
          // correct periodic image
          break;
        } else {
          // next periodic image
          jj = atom->sametag[jj];
        }
      }
      if (jj == -1)
        error->one(FLERR, "atom ID {} no correct image of csp neighbour with ID {} found",
                   atom->tag[i], stored_tags[i][j]);
      ngh_pairs[i][j] = jj;

      kk = atom->map((tagint) stored_tags[i][j + nhalf]);
      if (kk == -1)
        error->one(FLERR, "atom ID {} csp neighbour with ID {} not present", atom->tag[i],
                   stored_tags[i][j + nhalf]);
      while (kk >= 0) {
        // calculate distance to kk
        delx_k = x[kk][0] - xtmp;
        dely_k = x[kk][1] - ytmp;
        delz_k = x[kk][2] - ztmp;
        // correct periodic image?
        if (fabs(delx_k) < domain->prd_half[0] && fabs(dely_k) < domain->prd_half[1] &&
            fabs(delz_k) < domain->prd_half[2]) {
          // correct periodic image
          break;
        } else {
          // next periodic image
          kk = atom->sametag[kk];
        }
      }
      if (kk == -1)
        error->one(FLERR, "atom ID {} no correct image of csp neighbour with ID {} found",
                   atom->tag[i], stored_tags[i][j + nhalf]);
      ngh_pairs[i][j + nhalf] = kk;

      delx = delx_j + delx_k;
      dely = dely_j + dely_k;
      delz = delz_j + delz_k;

      value += delx * delx + dely * dely + delz * delz;
    }
    csp[i] = value;
  }

  // 2. compute csp avg
  // csp_avg_i += csp_i * w(r_ij)
  // csp_norm_i += w(r_ij)

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // own particle
    weight = weighting_function_poly(0);
    csp_avg[i] += csp[i] * weight;
    csp_norm[i] += weight;

    // neighbouring particles
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      // this is a full neighbour list
      // avoid double counting
      if (i > j) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq <= cut_hi_sq) {
        weight = weighting_function_poly(sqrt(rsq));

        csp_avg[j] += csp[i] * weight;
        csp_norm[j] += weight;

        if (j >= nlocal) continue;

        csp_norm[i] += weight;
        csp_avg[i] += csp[j] * weight;
      }
    }
  }

  // reverse communication of csp_norm and csp_avgs of ghost atoms
  comm->reverse_comm(this);

  for (i = 0; i < nlocal; i++) {
    csp_avg[i] /= csp_norm[i];

    if (mask[i] & groupbit)
      lambda[i] = switching_function_poly(csp_avg[i]);
    else
      lambda[i] = lambda_non_group;
  }
}

/* --------------------------------------------------------------------- */

void FixLambdaLACSPAPIP::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag, vflag);
}

/**
  * Calculate forces and store per-atom stats.
  */

void FixLambdaLACSPAPIP::pre_reverse(int /*eflag*/, int vflag)
{
  store_la();
  calculate_forces(vflag);
}

/**
  * Calculate derivative of the switching parameter for the forces.
  */

void FixLambdaLACSPAPIP::calculate_forces(int vflag)
{
  if (!calculate_forces_flag) return;

  int i, j, ii, jj, inum, jnum, i_pair, i1, i2, i3;
  int *ilist, *jlist, *numneigh, **firstneigh, *mask;
  double **x, **f, *lambda, *csp, *csp_avg, *csp_norm, *e_fast, *e_precise;
  double xtmp, ytmp, ztmp, fpair, delx, dely, delz, r, rsq, cspavgtmp, prefactortmp,
      delx1, dely1, delz1, delx2, dely2, delz2, tmp, ftmp[3];

  int nlocal = atom->nlocal;
  int nhalf = nnn / 2;

  v_init(vflag);

  x = atom->x;
  f = atom->f;
  lambda = atom->apip_lambda;
  csp = fixstore_la_inp->vstore;
  csp_avg = fixstore_la_avg->vstore;
  csp_norm = fixstore_la_norm->vstore;
  e_fast = atom->apip_e_fast;
  e_precise = atom->apip_e_precise;
  mask = atom->mask;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  if (atom->nmax > prefactor1_size) {
    memory->destroy(prefactor1);
    prefactor1_size = atom->nmax;
    memory->create(prefactor1, prefactor1_size, "pair:la:csp:apip:prefactor1");
  }
  if (nlocal > prefactor2_size) {
    memory->destroy(prefactor2);
    prefactor2_size = nlocal;
    memory->create(prefactor2, prefactor2_size, "pair:la:csp:apip:prefactor2");
  }

  // e_fast and e_precise are known only for own atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    // the derivative must be calculated from the switching function defined in FixLambdaLACSPAPIP::switching_function_poly
    if (lambda[i] != 0 && lambda[i] != 1 && (mask[i] & groupbit)) {
      // calculate derviative of lambda
      tmp = 1 - 2 * (1 + (csp_avg[i] - threshold_hi) / threshold_width);
      prefactor1[i] = -1.875 / threshold_width * (1 - 2 * tmp * tmp + pow(tmp, 4));
    } else {
      prefactor1[i] = 0;
    }

    // e_fast and e_precise are computed only for lambda in (0,1)
    // lambda in (0,1) implies that the derivative of the switching function is non-zero
    if (prefactor1[i] != 0) prefactor1[i] *= (e_fast[i] - e_precise[i]) / csp_norm[i];
    prefactor2[i] = prefactor1[i] * weighting_function_poly(0);
  }

  // communication of prefactor1 and csp
  // The csp was computed in pre_force, but is not known for ghosts yet.
  // send values of own atoms to neighbouring processors
  // get values of ghosts
  comm_forward_flag = FORWARD_PREFACTOR;
  comm->forward_comm(this);

  // store all forces before force computation
  store_f_lambda_before();

  // compute derivative of the radial weight function first
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    prefactortmp = prefactor1[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    cspavgtmp = csp_avg[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq >= cut_hi_sq) continue;

      r = sqrt(rsq);
      prefactor2[i] += prefactor1[j] * weighting_function_poly(r);
      fpair = prefactortmp * (csp[j] - cspavgtmp) * der_weighting_function_poly(r) / r;

      if (fpair == 0) continue;

      f[i][0] -= delx * fpair;
      f[i][1] -= dely * fpair;
      f[i][2] -= delz * fpair;
      f[j][0] += delx * fpair;
      f[j][1] += dely * fpair;
      f[j][2] += delz * fpair;

      if (evflag) ev_tally2(i, j, fpair, delx, dely, delz);
    }
  }

  // compute derivative of the CSP
  for (ii = 0; ii < inum; ii++) {
    i2 = ilist[ii];

    fpair = -2 * prefactor2[i2];

    if (fpair == 0) continue;

    xtmp = x[i2][0];
    ytmp = x[i2][1];
    ztmp = x[i2][2];

    for (i_pair = 0; i_pair < nhalf; i_pair++) {
      i1 = ngh_pairs[i2][i_pair];
      i3 = ngh_pairs[i2][i_pair + nhalf];

      delx1 = x[i1][0] - xtmp;
      dely1 = x[i1][1] - ytmp;
      delz1 = x[i1][2] - ztmp;

      delx2 = x[i3][0] - xtmp;
      dely2 = x[i3][1] - ytmp;
      delz2 = x[i3][2] - ztmp;

      ftmp[0] = (delx1 + delx2) * fpair;
      ftmp[1] = (dely1 + dely2) * fpair;
      ftmp[2] = (delz1 + delz2) * fpair;

      f[i1][0] += ftmp[0];
      f[i1][1] += ftmp[1];
      f[i1][2] += ftmp[2];

      f[i2][0] -= 2 * ftmp[0];
      f[i2][1] -= 2 * ftmp[1];
      f[i2][2] -= 2 * ftmp[2];

      f[i3][0] += ftmp[0];
      f[i3][1] += ftmp[1];
      f[i3][2] += ftmp[2];

      if (evflag) ev_tally3(i1, i2, i3, ftmp, delx1, dely1, delz1, delx2, dely2, delz2);
    }
  }

  // calculate force due to gradient lambda by comparison with before stored force
  store_f_lambda_after();
}

/* ----------------------------------------------------------------------
   select routine from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   and sort auxiliary array at same time
------------------------------------------------------------------------- */

void FixLambdaLACSPAPIP::select2(int k, int n, double *arr, int *iarr)
{
  int i, ir, j, l, mid, ia;
  double a;

  arr--;
  iarr--;
  l = 1;
  ir = n;
  while (true) {
    if (ir <= l + 1) {
      if (ir == l + 1 && arr[ir] < arr[l]) {
        std::swap(arr[l], arr[ir]);
        std::swap(iarr[l], iarr[ir]);
      }
      return;
    } else {
      mid = (l + ir) >> 1;
      std::swap(arr[mid], arr[l + 1]);
      std::swap(iarr[mid], iarr[l + 1]);
      if (arr[l] > arr[ir]) {
        std::swap(arr[l], arr[ir]);
        std::swap(iarr[l], iarr[ir]);
      }
      if (arr[l + 1] > arr[ir]) {
        std::swap(arr[l + 1], arr[ir]);
        std::swap(iarr[l + 1], iarr[ir]);
      }
      if (arr[l] > arr[l + 1]) {
        std::swap(arr[l], arr[l + 1]);
        std::swap(iarr[l], iarr[l + 1]);
      }
      i = l + 1;
      j = ir;
      a = arr[l + 1];
      ia = iarr[l + 1];
      while (true) {
        do i++;
        while (arr[i] < a);
        do j--;
        while (arr[j] > a);
        if (j < i) break;
        std::swap(arr[i], arr[j]);
        std::swap(iarr[i], iarr[j]);
      }
      arr[l + 1] = arr[j];
      arr[j] = a;
      iarr[l + 1] = iarr[j];
      iarr[j] = ia;
      if (j >= k) ir = j - 1;
      if (j <= k) l = i;
    }
  }
}

// helper function
// similar to cutoff_func_poly in ace_radial.cpp
// compare Phys Rev Mat 6, 013804 (2022) APPENDIX C: RADIAL AND CUTOFF FUNCTIONS 2. Cutoff function
// the first two derivatives of the switching function lambda vanishes at the boundaries of the switching region
double FixLambdaLACSPAPIP::switching_function_poly(double input)
{
  // NOTE: the derivative of this switching function is implemented and used in pair_la_csp_apip.cpp
  // calculate lambda
  if (input <= threshold_lo) {
    return 1;
  } else if (input >= threshold_hi) {
    return 0;
  } else {
    double tmp = 1 - 2 * (1 + (input - threshold_hi) / threshold_width);
    return 0.5 + 7.5 / 2. * (tmp / 4. - pow(tmp, 3) / 6. + pow(tmp, 5) / 20.);
  }
}

double FixLambdaLACSPAPIP::weighting_function_poly(double r)
{
  // calculate lambda
  if (r <= cut_lo) {
    return 1;
  } else if (r >= cut_hi) {
    return 0;
  } else {
    double tmp = 1 - 2 * (1 + (r - cut_hi) / cut_width);
    return 0.5 + 7.5 / 2. * (tmp / 4. - pow(tmp, 3) / 6. + pow(tmp, 5) / 20.);
  }
}

double FixLambdaLACSPAPIP::der_weighting_function_poly(double r)
{
  // calculate lambda
  if (r <= cut_lo || r >= cut_hi) {
    return 0;
  } else {
    double tmp = 1 - 2 * (1 + (r - cut_hi) / cut_width);
    return -1.875 / cut_width * (1 - 2 * tmp * tmp + pow(tmp, 4));
  }
}

/**
  * extract lambda(time averaged lambda_input) and lambda_input_history_len
  */

void *FixLambdaLACSPAPIP::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "fix_lambda_la_csp_apip:ngh_pairs") == 0) { return (void *) ngh_pairs; }
  dim = 0;
  if (strcmp(str, "fix_lambda_la_csp_apip:cut_lo") == 0) { return &cut_lo; }
  if (strcmp(str, "fix_lambda_la_csp_apip:cut_hi") == 0) { return &cut_hi; }
  if (strcmp(str, "fix_lambda_la_csp_apip:cut_hi_sq") == 0) { return &cut_hi_sq; }
  if (strcmp(str, "fix_lambda_la_csp_apip:cut_width") == 0) { return &cut_width; }
  if (strcmp(str, "fix_lambda_la_csp_apip:lambda_non_group") == 0) { return &lambda_non_group; }
  if (strcmp(str, "fix_lambda_la_csp_apip:cutsq_combined") == 0) { return &cutsq_combined; }
  if (strcmp(str, "fix_lambda_la_csp_apip:nnn") == 0) { return &nnn; }
  return nullptr;
}

/* ---------------------------------------------------------------------- */
int FixLambdaLACSPAPIP::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;
  double *avg = fixstore_la_avg->vstore;
  double *norm = fixstore_la_norm->vstore;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = avg[i];
    buf[m++] = norm[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixLambdaLACSPAPIP::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;
  double *avg = fixstore_la_avg->vstore;
  double *norm = fixstore_la_norm->vstore;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    avg[j] += buf[m++];
    norm[j] += buf[m++];
  }
}

/**
  * Send lambda/la_inp or prefactor to neighbours.
  */

int FixLambdaLACSPAPIP::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                          int * /*pbc*/)
{
  int i, j, m;
  m = 0;

  if (comm_forward_flag == FORWARD_INP_LAMBDA) {

    double *la_inp = fixstore_la_inp->vstore;
    double *lambda = atom->apip_lambda;

    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = la_inp[j];
      buf[m++] = lambda[j];
    }

  } else if (comm_forward_flag == FORWARD_PREFACTOR) {

    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = prefactor1[j];
    }
  }

  return m;
}

/**
  * Receive lambda/la_inp or prefactor to neighbours.
  */

void FixLambdaLACSPAPIP::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  if (comm_forward_flag == FORWARD_INP_LAMBDA) {

    double *lambda = atom->apip_lambda;
    double *la_inp = fixstore_la_inp->vstore;
    for (i = first; i < last; i++) {
      la_inp[i] = buf[m++];
      lambda[i] = buf[m++];
    }

  } else if (comm_forward_flag == FORWARD_PREFACTOR) {

    for (i = first; i < last; i++) prefactor1[i] = buf[m++];
  }
}

/**
  * Initial f_lambda.
  * This function is called from the force calculation routine.
  */

void FixLambdaLACSPAPIP::store_f_lambda_before()
{
  if ((!store_stats) || update->ntimestep % peratom_freq) { return; }

  int nlocal = atom->nlocal;
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) {
    f_lambda[i][0] = -f[i][0];
    f_lambda[i][1] = -f[i][1];
    f_lambda[i][2] = -f[i][2];
  }
}

/**
  * Calculate f_lambda.
  * This function is called from the force calculation routine.
  */

void FixLambdaLACSPAPIP::store_f_lambda_after()
{
  if ((!store_stats) || update->ntimestep % peratom_freq) { return; }

  int nlocal = atom->nlocal;
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) {
    f_lambda[i][0] += f[i][0];
    f_lambda[i][1] += f[i][1];
    f_lambda[i][2] += f[i][2];
  }
}

/**
  * Store per-atom stats regarding the local averaging and allocate the per-atom array.
  */

void FixLambdaLACSPAPIP::store_la()
{
  if ((!store_stats) || update->ntimestep % peratom_freq) { return; }

  int nlocal = atom->nlocal;
  double *inp = fixstore_la_inp->vstore;
  double *avg = fixstore_la_avg->vstore;

  // allocate more memory if required
  if (nlocal > size_f_lambda) {
    memory->destroy(f_lambda);
    size_f_lambda = nlocal;
    memory->create(f_lambda, size_f_lambda, size_peratom_cols, "pair:lambda:la:csp:apip:f:lambda");
    array_atom = f_lambda;
    // zero forces if required
    if (!calculate_forces_flag)
      for (int i = 0; i < nlocal; i++) f_lambda[i][0] = f_lambda[i][1] = f_lambda[i][2] = 0;
  }

  for (int i = 0; i < nlocal; i++) {
    // store local averaging stats
    f_lambda[i][3] = inp[i];
    f_lambda[i][4] = avg[i];
  }
}

CSPpairAPIP::CSPpairAPIP(tagint i0, tagint i1)
{
  if (i0 < i1) {
    tag_smaller = i0;
    tag_larger = i1;
  } else {
    tag_smaller = i1;
    tag_larger = i0;
  }
}

double FixLambdaLACSPAPIP::compute_scalar()
{
  int counter;
  MPI_Allreduce(&counter_changed_csp_nghs, &counter, 1, MPI_INT, MPI_SUM, world);
  return counter;
}

/**
   * store scalar information regarding the stored CSP-pairs
   */

void FixLambdaLACSPAPIP::write_restart(FILE *fp)
{
  int n = 0;
  double list[4];
  list[n++] = tags_stored;
  list[n++] = counter_changed_csp_nghs;
  list[n++] = const_ngh_flag;
  list[n++] = nnn;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list, sizeof(double), n, fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixLambdaLACSPAPIP::restart(char *buf)
{
  // store arguments passed to the fix for comparison with the ones in the restart file
  int nnn_tmp = nnn;
  bool const_ngh_flag_tmp = const_ngh_flag;

  int n = 0;
  auto *list = (double *) buf;

  tags_stored = static_cast<int>(list[n++]);
  counter_changed_csp_nghs = static_cast<int>(list[n++]);
  const_ngh_flag = static_cast<int>(list[n++]);
  nnn = (static_cast<int>(list[n++]));

  // simple comparisons first
  if (nnn_tmp != nnn)
    error->all(FLERR, "fix lambda/la/csp/apip: nnn = {} != {} = nnn in restart file", nnn_tmp, nnn);
  if (const_ngh_flag_tmp != const_ngh_flag)
    error->all(FLERR,
               "fix lambda/la/csp/apip: const_ngh_flag = {} != {} = const_ngh_flag in restart file",
               const_ngh_flag_tmp, const_ngh_flag);

  if (tags_stored) fixstore_pairs->disable = 0;
}

/**
  * compute the virial similar to Pair::ev_tally
  * There is no potential energy involved.
  * newton_pair is true since there are no double computations on different processors.
  */

void FixLambdaLACSPAPIP::ev_tally2(int i, int j, double fpair, double delx, double dely,
                                   double delz)
{
  double v[6];

  v[0] = delx * delx * fpair;
  v[1] = dely * dely * fpair;
  v[2] = delz * delz * fpair;
  v[3] = delx * dely * fpair;
  v[4] = delx * delz * fpair;
  v[5] = dely * delz * fpair;

  if (vflag_global) {
    virial[0] += v[0];
    virial[1] += v[1];
    virial[2] += v[2];
    virial[3] += v[3];
    virial[4] += v[4];
    virial[5] += v[5];
  }

  if (vflag_atom) {
    vatom[i][0] += 0.5 * v[0];
    vatom[i][1] += 0.5 * v[1];
    vatom[i][2] += 0.5 * v[2];
    vatom[i][3] += 0.5 * v[3];
    vatom[i][4] += 0.5 * v[4];
    vatom[i][5] += 0.5 * v[5];

    vatom[j][0] += 0.5 * v[0];
    vatom[j][1] += 0.5 * v[1];
    vatom[j][2] += 0.5 * v[2];
    vatom[j][3] += 0.5 * v[3];
    vatom[j][4] += 0.5 * v[4];
    vatom[j][5] += 0.5 * v[5];
  }
}

/**
  * compute the virial similar to Angle::ev_tally
  * Differences:
  * 1. There is no potential energy involved.
  * 2. newton_bond is true since there are no double computations
  *    on different processors.
  * 3. f1 = f3
  */

void FixLambdaLACSPAPIP::ev_tally3(int i, int j, int k, double *f, double delx1, double dely1,
                                   double delz1, double delx2, double dely2, double delz2)
{
  double v[6];

  v[0] = (delx1 + delx2) * f[0];
  v[1] = (dely1 + dely2) * f[1];
  v[2] = (delz1 + delz2) * f[2];
  v[3] = (delx1 + delx2) * f[1];
  v[4] = (delx1 + delx2) * f[2];
  v[5] = (dely1 + dely2) * f[2];

  if (vflag_global) {
    virial[0] += v[0];
    virial[1] += v[1];
    virial[2] += v[2];
    virial[3] += v[3];
    virial[4] += v[4];
    virial[5] += v[5];
  }

  if (vflag_atom) {
    vatom[i][0] += THIRD * v[0];
    vatom[i][1] += THIRD * v[1];
    vatom[i][2] += THIRD * v[2];
    vatom[i][3] += THIRD * v[3];
    vatom[i][4] += THIRD * v[4];
    vatom[i][5] += THIRD * v[5];

    vatom[j][0] += THIRD * v[0];
    vatom[j][1] += THIRD * v[1];
    vatom[j][2] += THIRD * v[2];
    vatom[j][3] += THIRD * v[3];
    vatom[j][4] += THIRD * v[4];
    vatom[j][5] += THIRD * v[5];

    vatom[k][0] += THIRD * v[0];
    vatom[k][1] += THIRD * v[1];
    vatom[k][2] += THIRD * v[2];
    vatom[k][3] += THIRD * v[3];
    vatom[k][4] += THIRD * v[4];
    vatom[k][5] += THIRD * v[5];
  }

  // per-atom centroid virial

  if (cvflag_atom) {

    // r0 = (r1+r2+r3)/3
    // rij = ri-rj
    // total virial = r10*f1 + r20*f2 + r30*f3
    // del1: r12
    // del2: r32

    double a1[3];

    // a1 = r10 = (2*r12 -   r32)/3
    a1[0] = THIRD * (2 * delx1 - delx2);
    a1[1] = THIRD * (2 * dely1 - dely2);
    a1[2] = THIRD * (2 * delz1 - delz2);

    cvatom[i][0] += a1[0] * f[0];
    cvatom[i][1] += a1[1] * f[1];
    cvatom[i][2] += a1[2] * f[2];
    cvatom[i][3] += a1[0] * f[1];
    cvatom[i][4] += a1[0] * f[2];
    cvatom[i][5] += a1[1] * f[2];
    cvatom[i][6] += a1[1] * f[0];
    cvatom[i][7] += a1[2] * f[0];
    cvatom[i][8] += a1[2] * f[1];

    double a2[3];
    double f2[3];

    // a2 = r20 = ( -r12 -   r32)/3
    a2[0] = THIRD * (-delx1 - delx2);
    a2[1] = THIRD * (-dely1 - dely2);
    a2[2] = THIRD * (-delz1 - delz2);

    f2[0] = -2 * f[0];
    f2[1] = -2 * f[1];
    f2[2] = -2 * f[2];

    cvatom[j][0] += a2[0] * f2[0];
    cvatom[j][1] += a2[1] * f2[1];
    cvatom[j][2] += a2[2] * f2[2];
    cvatom[j][3] += a2[0] * f2[1];
    cvatom[j][4] += a2[0] * f2[2];
    cvatom[j][5] += a2[1] * f2[2];
    cvatom[j][6] += a2[1] * f2[0];
    cvatom[j][7] += a2[2] * f2[0];
    cvatom[j][8] += a2[2] * f2[1];

    double a3[3];

    // a3 = r30 = ( -r12 + 2*r32)/3
    a3[0] = THIRD * (-delx1 + 2 * delx2);
    a3[1] = THIRD * (-dely1 + 2 * dely2);
    a3[2] = THIRD * (-delz1 + 2 * delz2);

    cvatom[k][0] += a3[0] * f[0];
    cvatom[k][1] += a3[1] * f[1];
    cvatom[k][2] += a3[2] * f[2];
    cvatom[k][3] += a3[0] * f[1];
    cvatom[k][4] += a3[0] * f[2];
    cvatom[k][5] += a3[1] * f[2];
    cvatom[k][6] += a3[1] * f[0];
    cvatom[k][7] += a3[2] * f[0];
    cvatom[k][8] += a3[2] * f[1];
  }
}
