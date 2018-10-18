/* ----------------------------------------------------------------------
         LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
         http://lammps.sandia.gov, Sandia National Laboratories
         Steve Plimpton, sjplimp@sandia.gov

         Copyright (2003) Sandia Corporation.	Under the terms of Contract
         DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
         certain rights in this software.	This software is distributed
under
         the GNU General Public License.

         See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
         Contributing author: PM Larsen (MIT)
------------------------------------------------------------------------- */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "atom.h"
#include "comm.h"
#include "compute_ptm_atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include "ptm_functions.h"

#define MAX_NEIGHBORS 30
#define NUM_COLUMNS 7
#define UNKNOWN 0
#define OTHER 8

using namespace LAMMPS_NS;

static const char cite_user_ptm_package[] =
    "USER-PTM package:\n\n"
    "@Article{larsen2016ptm,\n"
    " author={Larsen, Peter Mahler and Schmidt, S{\\o}ren and Schi{\\o}tz, "
    "Jakob},\n"
    " title={Robust structural identification via polyhedral template "
    "matching},\n"
    " journal={Modelling~Simul.~Mater.~Sci.~Eng.},\n"
    " year={2016},\n"
    " number={5},\n"
    " volume={24},\n"
    " pages={055007},\n"
    " DOI = {10.1088/0965-0393/24/5/055007}"
    "}\n\n";

/* ---------------------------------------------------------------------- */

ComputePTMAtom::ComputePTMAtom(LAMMPS *lmp, int narg, char **arg)
    : Compute(lmp, narg, arg), list(NULL), output(NULL) {
  if (narg != 5)
    error->all(FLERR, "Illegal compute ptm/atom command");

  char *structures = arg[3];
  char *ptr = structures;

  const char *strings[] = {"fcc",  "hcp",  "bcc", "ico",    "sc",
                           "dcub", "dhex", "all", "default"};
  int32_t flags[] = {
      PTM_CHECK_FCC,
      PTM_CHECK_HCP,
      PTM_CHECK_BCC,
      PTM_CHECK_ICO,
      PTM_CHECK_SC,
      PTM_CHECK_DCUB,
      PTM_CHECK_DHEX,
      PTM_CHECK_ALL,
      PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_BCC | PTM_CHECK_ICO};

  input_flags = 0;
  while (*ptr != '\0') {

    bool found = false;
    for (int i = 0; i < 9; i++) {
      int len = strlen(strings[i]);
      if (strncmp(ptr, strings[i], len) == 0) {
        input_flags |= flags[i];
        ptr += len;
        found = true;
        break;
      }
    }

    if (!found)
      error->all(FLERR,
                 "Illegal compute ptm/atom command (invalid structure type)");

    if (*ptr == '\0')
      break;

    if (*ptr != '-')
      error->all(FLERR,
                 "Illegal compute ptm/atom command (invalid structure type)");

    ptr++;
  }

  double threshold = force->numeric(FLERR, arg[4]);
  if (threshold < 0.0)
    error->all(FLERR,
               "Illegal compute ptm/atom command (threshold is negative)");
  rmsd_threshold = threshold;
  if (rmsd_threshold == 0)
    rmsd_threshold = INFINITY;

  peratom_flag = 1;
  size_peratom_cols = NUM_COLUMNS;
  create_attribute = 1;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputePTMAtom::~ComputePTMAtom() { memory->destroy(output); }

/* ---------------------------------------------------------------------- */

void ComputePTMAtom::init() {
  if (force->pair == NULL)
    error->all(FLERR, "Compute ptm/atom requires a pair style be defined");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style, "ptm/atom") == 0)
      count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR, "More than one compute ptm/atom defined");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputePTMAtom::init_list(int id, NeighList *ptr) { list = ptr; }

/* ---------------------------------------------------------------------- */

typedef struct {
  int index;
  double d;
} ptmnbr_t;

static bool sorthelper_compare(ptmnbr_t const &a, ptmnbr_t const &b) {
  return a.d < b.d;
}

static int get_neighbors(double *pos, int jnum, int *jlist, double **x,
                         double (*nbr)[3]) {

  ptmnbr_t *nbr_order = new ptmnbr_t[jnum];

  for (int jj = 0; jj < jnum; jj++) {
    int j = jlist[jj];
    j &= NEIGHMASK;

    double dx = pos[0] - x[j][0];
    double dy = pos[1] - x[j][1];
    double dz = pos[2] - x[j][2];
    double rsq = dx * dx + dy * dy + dz * dz;

    nbr_order[jj].index = j;
    nbr_order[jj].d = rsq;
  }

  std::sort(nbr_order, nbr_order + jnum, &sorthelper_compare);
  int num_nbrs = std::min(MAX_NEIGHBORS, jnum);

  nbr[0][0] = nbr[0][1] = nbr[0][2] = 0;
  for (int jj = 0; jj < num_nbrs; jj++) {

    int j = nbr_order[jj].index;
    nbr[jj + 1][0] = x[j][0] - pos[0];
    nbr[jj + 1][1] = x[j][1] - pos[1];
    nbr[jj + 1][2] = x[j][2] - pos[2];
  }

  delete[] nbr_order;
  return num_nbrs;
}

void ComputePTMAtom::compute_peratom() {
  // PTM global initialization.  If already initialized this function does
  // nothing.
  ptm_initialize_global();

  // initialize PTM local storage
  ptm_local_handle_t local_handle = ptm_initialize_local();

  invoked_peratom = update->ntimestep;

  // grow arrays if necessary
  if (atom->nmax > nmax) {
    memory->destroy(output);
    nmax = atom->nmax;

    memory->create(output, nmax, NUM_COLUMNS, "ptm:ptm_output");
    array_atom = output;
  }

  // invoke full neighbor list (will copy or build if necessary)
  neighbor->build_one(list);

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int ii = 0; ii < inum; ii++) {

    int i = ilist[ii];
    output[i][0] = UNKNOWN;
    if (!(mask[i] & groupbit))
      continue;

    double *pos = x[i];

    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    if (jnum <= 0)
      continue;

    // get neighbours ordered by increasing distance
    double nbr[MAX_NEIGHBORS + 1][3];
    int num_nbrs = get_neighbors(pos, jnum, jlist, x, nbr);

    // check that we have enough neighbours for the desired structure types
    int32_t flags = 0;
    if (num_nbrs >= PTM_NUM_NBRS_SC && (input_flags & PTM_CHECK_SC))
      flags |= PTM_CHECK_SC;
    if (num_nbrs >= PTM_NUM_NBRS_FCC && (input_flags & PTM_CHECK_FCC))
      flags |= PTM_CHECK_FCC;
    if (num_nbrs >= PTM_NUM_NBRS_HCP && (input_flags & PTM_CHECK_HCP))
      flags |= PTM_CHECK_HCP;
    if (num_nbrs >= PTM_NUM_NBRS_ICO && (input_flags & PTM_CHECK_ICO))
      flags |= PTM_CHECK_ICO;
    if (num_nbrs >= PTM_NUM_NBRS_BCC && (input_flags & PTM_CHECK_BCC))
      flags |= PTM_CHECK_BCC;
    if (num_nbrs >= PTM_NUM_NBRS_DCUB && (input_flags & PTM_CHECK_DCUB))
      flags |= PTM_CHECK_DCUB;
    if (num_nbrs >= PTM_NUM_NBRS_DHEX && (input_flags & PTM_CHECK_DHEX))
      flags |= PTM_CHECK_DHEX;

    // now run PTM
    int8_t mapping[MAX_NEIGHBORS + 1];
    int32_t type, alloy_type;
    double scale, rmsd, interatomic_distance, lattice_constant;
    double q[4], F[9], F_res[3], U[9], P[9];
    ptm_index(local_handle, flags, num_nbrs + 1, nbr, NULL, true, &type,
              &alloy_type, &scale, &rmsd, q, F, F_res, U, P, mapping,
              &interatomic_distance, &lattice_constant);

    if (rmsd > rmsd_threshold) {
      type = PTM_MATCH_NONE;
    }

    // printf("%d type=%d rmsd=%f\n", i, type, rmsd);

    if (type == PTM_MATCH_NONE)
      type = OTHER;

    output[i][0] = type;
    output[i][1] = rmsd;
    output[i][2] = interatomic_distance;
    output[i][3] = q[0];
    output[i][4] = q[1];
    output[i][5] = q[2];
    output[i][6] = q[3];
  }

  // printf("finished ptm analysis\n");
  ptm_uninitialize_local(local_handle);
}

/* ----------------------------------------------------------------------
         memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePTMAtom::memory_usage() {
  double bytes = nmax * NUM_COLUMNS * sizeof(double);
  bytes += nmax * sizeof(double);
  return bytes;
}
