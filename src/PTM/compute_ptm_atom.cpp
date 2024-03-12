// clang-format off
/* ----------------------------------------------------------------------
         LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
         https://www.lammps.org/, Sandia National Laboratories
         LAMMPS development team: developers@lammps.org

         Copyright (2003) Sandia Corporation.  Under the terms of Contract
         DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
         certain rights in this software.  This software is distributed
under
         the GNU General Public License.

         See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
         Contributing author: PM Larsen (MIT)
------------------------------------------------------------------------- */

#include "compute_ptm_atom.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>

#include "ptm_constants.h"
#include "ptm_functions.h"
#include "ptm_initialize_data.h"

#define NUM_COLUMNS 7
#define PTM_LAMMPS_UNKNOWN -1
#define PTM_LAMMPS_OTHER 0

using namespace LAMMPS_NS;

static const char cite_user_ptm_package[] =
    "PTM package: doi:10.1088/0965-0393/24/5/055007\n\n"
    "@Article{larsen2016ptm,\n"
    " author={Larsen, Peter Mahler and Schmidt, S{\\o}ren and\n"
    "    Schi{\\o}tz, Jakob},\n"
    " title={Robust Structural Identification via Polyhedral Template\n"
    "    Matching},\n"
    " journal={Model.\\ Simulat.\\ Mater.\\ Sci.\\ Eng.},\n"
    " year={2016},\n"
    " number={5},\n"
    " volume={24},\n"
    " pages={055007},\n"
    " DOI = {10.1088/0965-0393/24/5/055007}"
    "}\n\n";

/* ---------------------------------------------------------------------- */

ComputePTMAtom::ComputePTMAtom(LAMMPS *lmp, int narg, char **arg)
    : Compute(lmp, narg, arg), list(nullptr), output(nullptr) {
  if (narg < 5 || narg > 6)
    error->all(FLERR, "Illegal compute ptm/atom command");

  char *structures = arg[3];
  char *ptr = structures;

  const char *strings[] = {"fcc",  "hcp",  "bcc", "ico",    "sc",
                           "dcub", "dhex", "graphene", "all", "default"};
  int num_strings = sizeof(strings) / sizeof(const char*);

  int32_t flags[] = {
      PTM_CHECK_FCC,
      PTM_CHECK_HCP,
      PTM_CHECK_BCC,
      PTM_CHECK_ICO,
      PTM_CHECK_SC,
      PTM_CHECK_DCUB,
      PTM_CHECK_DHEX,
      PTM_CHECK_GRAPHENE,
      PTM_CHECK_ALL,
      PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_BCC | PTM_CHECK_ICO};

  if (lmp->citeme)
    lmp->citeme->add(cite_user_ptm_package);

  input_flags = 0;
  while (*ptr != '\0') {

    bool found = false;
    for (int i = 0; i < num_strings; i++) {
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

  double threshold = utils::numeric(FLERR, arg[4],false,lmp);
  if (threshold < 0.0)
    error->all(FLERR,
               "Illegal compute ptm/atom command (threshold is negative)");
  rmsd_threshold = threshold;
  if (rmsd_threshold == 0)
    rmsd_threshold = INFINITY;

  auto  group_name = (char *)"all";
  if (narg > 5) {
    group_name = arg[5];
  }
  int igroup = group->find(group_name);
  if (igroup == -1) error->all(FLERR,"Could not find fix group ID");
  group2bit = group->bitmask[igroup];

  peratom_flag = 1;
  size_peratom_cols = NUM_COLUMNS;
  create_attribute = 1;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputePTMAtom::~ComputePTMAtom() { memory->destroy(output); }

/* ---------------------------------------------------------------------- */

void ComputePTMAtom::init() {
  if (force->pair == nullptr)
    error->all(FLERR, "Compute ptm/atom requires a pair style be defined");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style, "ptm/atom") == 0)
      count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR, "More than one compute ptm/atom defined");

  // need an occasional full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputePTMAtom::init_list(int /* id */, NeighList *ptr) { list = ptr; }

/* ---------------------------------------------------------------------- */

typedef struct
{
  double **x;
  int *numneigh;
  int **firstneigh;
  int *ilist;
  int nlocal;
  int *mask;
  int group2bit;

} ptmnbrdata_t;


typedef struct {
  int index;
  double d;
} ptmnbr_t;

static bool sorthelper_compare(ptmnbr_t const &a, ptmnbr_t const &b) {
  return a.d < b.d;
}

static int get_neighbours(void* vdata, size_t central_index, size_t atom_index, int num, size_t* nbr_indices, int32_t* numbers, double (*nbr_pos)[3])
{
  auto  data = (ptmnbrdata_t*)vdata;
  int *mask = data->mask;
  int group2bit = data->group2bit;

  double **x = data->x;
  double *pos = x[atom_index];

  int *jlist = nullptr;
  int jnum = 0;
  if ((int)atom_index < data->nlocal) {
    jlist = data->firstneigh[atom_index];
    jnum = data->numneigh[atom_index];
  }
  else {
    jlist = data->firstneigh[central_index];
    jnum = data->numneigh[central_index];
  }

  std::vector<ptmnbr_t> nbr_order;

  for (int jj = 0; jj < jnum; jj++) {
    int j = jlist[jj];
    if (!(mask[j] & group2bit))
      continue;

    j &= NEIGHMASK;
    if (j == (int)atom_index)
      continue;

    double dx = pos[0] - x[j][0];
    double dy = pos[1] - x[j][1];
    double dz = pos[2] - x[j][2];
    double rsq = dx * dx + dy * dy + dz * dz;

    ptmnbr_t nbr = {j, rsq};
    nbr_order.push_back(nbr);
  }

  std::sort(nbr_order.begin(), nbr_order.end(), &sorthelper_compare);
  int num_nbrs = std::min(num - 1, (int)nbr_order.size());

  nbr_pos[0][0] = nbr_pos[0][1] = nbr_pos[0][2] = 0;
  nbr_indices[0] = atom_index;
  numbers[0] = 0;
  for (int jj = 0; jj < num_nbrs; jj++) {

    int j = nbr_order[jj].index;
    nbr_pos[jj + 1][0] = x[j][0] - pos[0];
    nbr_pos[jj + 1][1] = x[j][1] - pos[1];
    nbr_pos[jj + 1][2] = x[j][2] - pos[2];

    nbr_indices[jj + 1] = j;
    numbers[jj + 1] = 0;
  }

  return num_nbrs + 1;
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
  ptmnbrdata_t nbrlist = {x, numneigh, firstneigh, ilist, atom->nlocal, mask, group2bit};

  // zero output

  memset(&output[0][0],0,nmax*NUM_COLUMNS*sizeof(double));

  for (int ii = 0; ii < inum; ii++) {

    int i = ilist[ii];
    output[i][0] = PTM_LAMMPS_UNKNOWN;
    if (!(mask[i] & groupbit))
      continue;

    int jnum = numneigh[i];
    if (jnum <= 0)
      continue;


    // now run PTM
    int32_t type, alloy_type;
    double scale, rmsd, interatomic_distance;
    double q[4];
    bool standard_orientations = false;

    rmsd = INFINITY;
    interatomic_distance = q[0] = q[1] = q[2] = q[3] = 0.0;

    ptm_index(local_handle, i, get_neighbours, (void*)&nbrlist,
              input_flags, standard_orientations,
              &type, &alloy_type, &scale, &rmsd, q,
              nullptr, nullptr, nullptr, nullptr, &interatomic_distance, nullptr, nullptr);

    if (rmsd > rmsd_threshold) type = PTM_MATCH_NONE;
    if (type == PTM_MATCH_NONE) type = PTM_LAMMPS_OTHER;

    output[i][0] = type;
    output[i][1] = rmsd;
    output[i][2] = interatomic_distance;
    output[i][3] = q[0];
    output[i][4] = q[1];
    output[i][5] = q[2];
    output[i][6] = q[3];
  }

  ptm_uninitialize_local(local_handle);
}

/* ----------------------------------------------------------------------
         memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePTMAtom::memory_usage() {
  double bytes = (double)nmax * NUM_COLUMNS * sizeof(double);
  bytes += (double)nmax * sizeof(double);
  return bytes;
}
