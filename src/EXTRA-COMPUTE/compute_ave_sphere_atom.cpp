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

#include "compute_ave_sphere_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"
#include "math_const.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;


/* ---------------------------------------------------------------------- */

ComputeAveSphereAtom::ComputeAveSphereAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  result(nullptr)
{
  if (narg < 3 || narg > 5) error->all(FLERR,"Illegal compute ave/sphere/atom command");

  // process optional args

  cutoff = 0.0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute ave/sphere/atom command");
      cutoff = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (cutoff <= 0.0) error->all(FLERR,"Illegal compute ave/sphere/atom command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute ave/sphere/atom command");
  }

  peratom_flag = 1;
  size_peratom_cols = 2;
  comm_forward = 3;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeAveSphereAtom::~ComputeAveSphereAtom()
{
  if (copymode) return;

  memory->destroy(result);
}

/* ---------------------------------------------------------------------- */

void ComputeAveSphereAtom::init()
{
  if (!force->pair && cutoff == 0.0)
    error->all(FLERR,"Compute ave/sphere/atom requires a cutoff be specified "
               "or a pair style be defined");

  double skin = neighbor->skin;
  if (cutoff != 0.0) {
    double cutghost;            // as computed by Neighbor and Comm
    if (force->pair)
      cutghost = MAX(force->pair->cutforce+skin,comm->cutghostuser);
    else
      cutghost = comm->cutghostuser;

    if (cutoff > cutghost)
      error->all(FLERR,"Compute ave/sphere/atom cutoff exceeds ghost atom range - "
                 "use comm_modify cutoff command");
  }

  int cutflag = 1;
  if (force->pair) {
    if (cutoff == 0.0) {
      cutoff = force->pair->cutforce;
    }
    if (cutoff <= force->pair->cutforce+skin) cutflag = 0;
  }

  cutsq = cutoff*cutoff;
  sphere_vol = 4.0/3.0*MY_PI*cutsq*cutoff;

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
  if (cutflag) {
    neighbor->requests[irequest]->cut = 1;
    neighbor->requests[irequest]->cutoff = cutoff;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeAveSphereAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeAveSphereAtom::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int count;
  double vsum[3],vavg[3],vnet[3];

  invoked_peratom = update->ntimestep;

  // grow result array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(result);
    nmax = atom->nmax;
    memory->create(result,nmax,2,"ave/sphere/atom:result");
    array_atom = result;
  }

  // need velocities of ghost atoms

  comm->forward_comm(this);

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute properties for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // i atom contribution

      count = 1;
      vsum[0] = v[i][0];
      vsum[1] = v[i][1];
      vsum[2] = v[i][2];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          count++;
          vsum[0] += v[j][0];
          vsum[1] += v[j][1];
          vsum[2] += v[j][2];
        }
      }

      vavg[0] = vsum[0]/count;
      vavg[1] = vsum[1]/count;
      vavg[2] = vsum[2]/count;

      // i atom contribution

      count = 1;
      vnet[0] = v[i][0] - vavg[0];
      vnet[1] = v[i][1] - vavg[1];
      vnet[2] = v[i][2] - vavg[2];
      double ke_sum = vnet[0]*vnet[0] + vnet[1]*vnet[1] + vnet[2]*vnet[2];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          count++;
          vnet[0] = v[j][0] - vavg[0];
          vnet[1] = v[j][1] - vavg[1];
          vnet[2] = v[j][2] - vavg[2];
          ke_sum += vnet[0]*vnet[0] + vnet[1]*vnet[1] + vnet[2]*vnet[2];
        }
      }
      double density = count/sphere_vol;
      double temp = ke_sum/3.0/count;
      result[i][0] = density;
      result[i][1] = temp;
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeAveSphereAtom::pack_forward_comm(int n, int *list, double *buf,
                                        int /*pbc_flag*/, int * /*pbc*/)
{
  double **v = atom->v;

  int i,m=0;
  for (i = 0; i < n; ++i) {
    buf[m++] = v[list[i]][0];
    buf[m++] = v[list[i]][1];
    buf[m++] = v[list[i]][2];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeAveSphereAtom::unpack_forward_comm(int n, int first, double *buf)
{
  double **v = atom->v;

  int i,last,m=0;
  last = first + n;
  for (i = first; i < last; ++i) {
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeAveSphereAtom::memory_usage()
{
  double bytes = (double)2*nmax * sizeof(double);
  return bytes;
}
