/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <cstring>
#include "compute_contact_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeContactAtom::ComputeContactAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  contact(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal compute contact/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_reverse = 1;

  nmax = 0;

  // error checks

  if (!atom->sphere_flag)
    error->all(FLERR,"Compute contact/atom requires atom style sphere");
}

/* ---------------------------------------------------------------------- */

ComputeContactAtom::~ComputeContactAtom()
{
  memory->destroy(contact);
}

/* ---------------------------------------------------------------------- */

void ComputeContactAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute contact/atom requires a pair style be defined");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"contact/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute contact/atom");

  // need an occasional neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->size = 1;
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
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

  for (i = 0; i < nall; i++) contact[i] = 0.0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      radi = radius[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radsum = radi + radius[j];
        radsumsq = radsum*radsum;
        if (rsq <= radsumsq) {
          contact[i] += 1.0;
          contact[j] += 1.0;
        }
      }
    }
  }

  // communicate ghost atom counts between neighbor procs if necessary

  if (force->newton_pair) comm->reverse_comm_compute(this);
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
  double bytes = nmax * sizeof(double);
  return bytes;
}
