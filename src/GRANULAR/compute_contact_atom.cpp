// clang-format off
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

#include "compute_contact_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeContactAtom::ComputeContactAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), group2(nullptr), contact(nullptr)
{
  if ((narg != 3) && (narg != 4)) error->all(FLERR, "Illegal compute contact/atom command");

  jgroup = group->find("all");
  jgroupbit = group->bitmask[jgroup];
  if (narg == 4) {
    group2 = utils::strdup(arg[3]);
    jgroup = group->find(group2);
    if (jgroup == -1) error->all(FLERR, "Compute contact/atom group2 ID {} does not exist", group2);
    jgroupbit = group->bitmask[jgroup];
  }

  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_reverse = 1;

  nmax = 0;

  // error checks

  if (!atom->sphere_flag) error->all(FLERR, "Compute contact/atom requires atom style sphere");
}

/* ---------------------------------------------------------------------- */

ComputeContactAtom::~ComputeContactAtom()
{
  memory->destroy(contact);
  delete[] group2;
}

/* ---------------------------------------------------------------------- */

void ComputeContactAtom::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute contact/atom requires a pair style be defined");

  if (modify->get_compute_by_style("contact/atom").size() > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute contact/atom");

  // need an occasional neighbor list

  neighbor->add_request(this, NeighConst::REQ_SIZE | NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputeContactAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeContactAtom::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,radsumsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow contact array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(contact);
    nmax = atom->nmax;
    memory->create(contact,nmax,"contact/atom:contact");
    vector_atom = contact;
  }

  // invoke neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute number of contacts for each atom in group
  // contact if distance <= sum of radii
  // tally for both I and J

  double **x = atom->x;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  bool update_i_flag, update_j_flag;

  for (i = 0; i < nall; i++) contact[i] = 0.0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    // Only proceed if i is either part of the compute group or will contribute to contacts
    if (! (mask[i] & groupbit) && ! (mask[i] & jgroupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      // Only tally for atoms in compute group (groupbit) if neighbor is in group2 (jgroupbit)
      update_i_flag = (mask[i] & groupbit) && (mask[j] & jgroupbit);
      update_j_flag = (mask[j] & groupbit) && (mask[i] & jgroupbit);
      if (! update_i_flag && ! update_j_flag) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      radsum = radi + radius[j];
      radsumsq = radsum * radsum;
      if (rsq <= radsumsq) {
        if (update_i_flag) contact[i] += 1.0;
        if (update_j_flag) contact[j] += 1.0;
      }
    }
  }

  // communicate ghost atom counts between neighbor procs if necessary

  if (force->newton_pair) comm->reverse_comm(this);
}

/* ---------------------------------------------------------------------- */

int ComputeContactAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    buf[m++] = contact[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeContactAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    contact[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeContactAtom::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  return bytes;
}
