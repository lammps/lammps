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

#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include "fix_hyper_global.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "modify.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define DELTA 16384
#define VECLEN 5

// NOTE: count/output # of timesteps on which bias is non-zero
// NOTE: should there be a virial contribution from boosted bond?
// NOTE: allow newton off?  see Note in pre_reverse()

/* ---------------------------------------------------------------------- */

FixHyperGlobal::FixHyperGlobal(LAMMPS *lmp, int narg, char **arg) :
  FixHyper(lmp, narg, arg), blist(NULL), xold(NULL), tagold(NULL)
{
  if (atom->map_style == 0)
    error->all(FLERR,"Fix hyper/global command requires atom map");

  if (narg != 7) error->all(FLERR,"Illegal fix hyper/global command");

  hyperflag = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 11;
  global_freq = 1;
  extscalar = 0;
  extvector = 0;

  cutbond = force->numeric(FLERR,arg[3]);
  qfactor = force->numeric(FLERR,arg[4]);
  vmax = force->numeric(FLERR,arg[5]);
  tequil = force->numeric(FLERR,arg[6]);

  if (cutbond < 0.0 || qfactor <= 0.0 || vmax < 0.0 || tequil <= 0.0)
    error->all(FLERR,"Illegal fix hyper/global command");

  invqfactorsq = 1.0 / (qfactor*qfactor);
  cutbondsq = cutbond*cutbond;
  beta = 1.0 / (force->boltz * tequil);

  maxbond = 0;
  nblocal = 0;
  blist = NULL;

  maxold = 0;
  xold = NULL;
  tagold = NULL;

  me = comm->me;
  firstflag = 1;
  bcastflag = 0;
  for (int i = 0; i < VECLEN; i++) outvec[i] = 0.0;

  nevent = 0;
  nevent_atom = 0;
  t_hyper = 0.0;
}

/* ---------------------------------------------------------------------- */

FixHyperGlobal::~FixHyperGlobal()
{
  memory->sfree(blist);
  memory->destroy(xold);
  memory->destroy(tagold);
}

/* ---------------------------------------------------------------------- */

int FixHyperGlobal::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= PRE_FORCE;
  mask |= PRE_REVERSE;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHyperGlobal::init_hyper()
{
  maxdriftsq = 0.0;
  maxbondlen = 0.0;
  nobias = 0;
}

/* ---------------------------------------------------------------------- */

void FixHyperGlobal::init()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Hyper global requires newton pair on");

  dt = update->dt;

  // need an occasional half neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void FixHyperGlobal::init_list(int /* id */, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixHyperGlobal::setup_pre_neighbor()
{
  pre_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixHyperGlobal::setup_pre_reverse(int eflag, int vflag)
{
  // no increment in nobias or hyper time when pre-run forces are calculated

  int nobias_hold = nobias;
  double t_hyper_hold = t_hyper;

  pre_reverse(eflag,vflag);

  nobias = nobias_hold;
  t_hyper = t_hyper_hold;
}

/* ---------------------------------------------------------------------- */

void FixHyperGlobal::pre_neighbor()
{
  int m,iold,jold,ilocal,jlocal;
  double distsq;

  // reset local IDs for owned bond atoms, since atoms have migrated
  // uses xold and tagold from when bonds were created
  // must be done after ghost atoms are setup via comm->borders()

  double **x = atom->x;

  int flag = 0;

  for (m = 0; m < nblocal; m++) {
    iold = blist[m].iold;
    jold = blist[m].jold;
    ilocal = atom->map(tagold[iold]);
    jlocal = atom->map(tagold[jold]);
    ilocal = domain->closest_image(xold[iold],ilocal);
    jlocal = domain->closest_image(xold[iold],jlocal);
    blist[m].i = ilocal;
    blist[m].j = jlocal;

    if (ilocal < 0 || jlocal < 0) flag++;
    else {
      distsq = MathExtra::distsq3(x[ilocal],xold[iold]);
      maxdriftsq = MAX(distsq,maxdriftsq);
    }
  }

  if (flag) error->one(FLERR,"Fix hyper/global bond atom not found");
}

/* ---------------------------------------------------------------------- */

void FixHyperGlobal::pre_reverse(int /* eflag */, int /* vflag */)
{
  int i,j,m,imax,jmax;
  double delx,dely,delz;
  double r,r0,estrain,rmax,r0max,emax,dt_boost;
  double vbias,fbias,fbiasr;

  // compute current strain of each owned bond
  // emax = maximum strain of any bond I own
  // imax,jmax = local indices of my 2 atoms in that bond

  double **x = atom->x;
  emax = 0.0;

  for (m = 0; m < nblocal; m++) {
    i = blist[m].i;
    j = blist[m].j;
    delx = x[i][0] - x[j][0];
    dely = x[i][1] - x[j][1];
    delz = x[i][2] - x[j][2];
    r = sqrt(delx*delx + dely*dely + delz*delz);
    maxbondlen = MAX(r,maxbondlen);
    r0 = blist[m].r0;
    estrain = fabs(r-r0) / r0;

    if (estrain > emax) {
      emax = estrain;
      rmax = r;
      r0max = r0;
      imax = i;
      jmax = j;
    }
  }

  // perform Allreduce() on pairme = double/int pair
  // finds max strain and what proc owns it
  // owner = proc that owns that bond

  pairme.value = emax;
  pairme.proc = me;
  MPI_Allreduce(&pairme,&pairall,1,MPI_DOUBLE_INT,MPI_MAXLOC,world);
  owner = pairall.proc;

  bcastflag = 1;

  // all procs acquire t_hyper from owner proc
  // MPI_Bcast call by owner proc is below

  for (int i = 0; i < VECLEN; i++) outvec[i] = 0.0;

  if (owner != me) {
    MPI_Bcast(&t_hyper,1,MPI_DOUBLE,owner,world);
    return;
  }

  // I own the bond with max strain
  // compute Vbias and apply force to atoms imax,jmax
  // NOTE: logic would need to be different for newton off

  double **f = atom->f;

  vbias = fbias = 0.0;
  dt_boost = 1.0;

  if (emax < qfactor) {
    vbias = vmax * (1.0 - emax*emax*invqfactorsq);
    fbias = 2.0 * vmax * emax / (qfactor*qfactor * r0max);
    dt_boost = exp(beta*vbias);

    delx = x[imax][0] - x[jmax][0];
    dely = x[imax][1] - x[jmax][1];
    delz = x[imax][2] - x[jmax][2];

    fbiasr = fbias / rmax;
    f[imax][0] += delx*fbiasr;
    f[imax][1] += dely*fbiasr;
    f[imax][2] += delz*fbiasr;

    f[jmax][0] -= delx*fbiasr;
    f[jmax][1] -= dely*fbiasr;
    f[jmax][2] -= delz*fbiasr;
  } else nobias++;

  // output quantities

  outvec[0] = vbias;
  outvec[1] = dt_boost;
  outvec[2] = emax;
  outvec[3] = atom->tag[imax];
  outvec[4] = atom->tag[jmax];

  t_hyper += dt_boost*dt;
  MPI_Bcast(&t_hyper,1,MPI_DOUBLE,owner,world);
}

/* ---------------------------------------------------------------------- */

void FixHyperGlobal::build_bond_list(int natom)
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  if (natom) {
    nevent++;
    nevent_atom += natom;
  }

  // trigger neighbor list build

  neighbor->build_one(list);

  // identify bonds assigned to each owned atom
  // do not create a bond between two non-group atoms

  double **x = atom->x;
  int *mask = atom->mask;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  nblocal = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      // skip if neither atom I or J are in fix group

      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutbondsq) {
        if (nblocal == maxbond) grow_bond();
        blist[nblocal].i = i;
        blist[nblocal].j = j;
        blist[nblocal].iold = i;
        blist[nblocal].jold = j;
        blist[nblocal].r0 = sqrt(rsq);
        nblocal++;
      }
    }
  }

  // store IDs and coords for owned+ghost atoms at time of bond creation
  // realloc xold and tagold as needed

  if (atom->nmax > maxold) {
    memory->destroy(xold);
    memory->destroy(tagold);
    maxold = atom->nmax;
    memory->create(xold,maxold,3,"hyper/global:xold");
    memory->create(tagold,maxold,"hyper/global:tagold");
  }

  tagint *tag = atom->tag;
  int nall = atom->nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    xold[i][0] = x[i][0];
    xold[i][1] = x[i][1];
    xold[i][2] = x[i][2];
    tagold[i] = tag[i];
  }
}

/* ----------------------------------------------------------------------
   grow bond list by a chunk
------------------------------------------------------------------------- */

void FixHyperGlobal::grow_bond()
{
  // NOTE: could add int arg to do initial large alloc:
  // maxbond = maxbond/DELTA * DELTA; maxbond += DELTA;

  maxbond += DELTA;
  if (maxbond < 0 || maxbond > MAXSMALLINT)
    error->one(FLERR,"Fix hyper/local per-processor bond count is too big");
  blist = (OneBond *)
    memory->srealloc(blist,maxbond*sizeof(OneBond),"hyper/local:blist");
}

/* ---------------------------------------------------------------------- */

double FixHyperGlobal::compute_scalar()
{
  if (bcastflag) {
    MPI_Bcast(outvec,VECLEN,MPI_DOUBLE,owner,world);
    bcastflag = 0;
  }
  return outvec[0];
}

/* ---------------------------------------------------------------------- */

double FixHyperGlobal::compute_vector(int i)
{
  if (bcastflag) {
    MPI_Bcast(outvec,VECLEN,MPI_DOUBLE,owner,world);
    bcastflag = 0;
  }

  // 11 vector outputs returned for i = 0-10

  // i = 0 = boost factor on this step
  // i = 1 = max strain of any bond on this step
  // i = 2 = ID of atom I in max-strain bond on this step
  // i = 3 = ID of atom J in max-strain bond on this step
  // i = 4 = ave bonds/atom on this step

  // i = 5 = fraction of steps with no bias during this run
  // i = 6 = max drift of any atom during this run
  // i = 7 = max bond length during this run

  // i = 8 = cummulative hyper time since fix created
  // i = 9 = cummulative # of event timesteps since fix created
  // i = 10 = cummulative # of atoms in events since fix created

  if (i == 0) return outvec[1];
  if (i == 1) return outvec[2];
  if (i == 2) return outvec[3];
  if (i == 3) return outvec[4];

  if (i == 4) {
    int allbonds;     // NOTE: bigint?
    MPI_Allreduce(&nblocal,&allbonds,1,MPI_INT,MPI_SUM,world);
    return 2.0*allbonds/atom->natoms;
  }

  if (i == 5) {
    if (update->ntimestep == update->firststep) return 0.0;
    int allnobias;
    MPI_Allreduce(&nobias,&allnobias,1,MPI_INT,MPI_SUM,world);
    return 1.0*allnobias / (update->ntimestep - update->firststep);
  }

  if (i == 6) {
    double alldriftsq;
    MPI_Allreduce(&maxdriftsq,&alldriftsq,1,MPI_DOUBLE,MPI_MAX,world);
    return sqrt(alldriftsq);
  }

  if (i == 7) {
    double allbondlen;
    MPI_Allreduce(&maxbondlen,&allbondlen,1,MPI_DOUBLE,MPI_MAX,world);
    return allbondlen;
  }

  if (i == 8) return t_hyper;
  if (i == 9) return (double) nevent;
  if (i == 10) return (double) nevent_atom;

  return 0.0;
}

/* ----------------------------------------------------------------------
   wrapper on compute_vector()
   used by hyper.cpp to call FixHyper
------------------------------------------------------------------------- */

double FixHyperGlobal::query(int i)
{
  if (i == 1) return compute_vector(8);  // cummulative hyper time
  if (i == 2) return compute_vector(9);  // nevent
  if (i == 3) return compute_vector(10); // nevent_atom
  if (i == 4) return compute_vector(4);  // ave bonds/atom
  if (i == 5) return compute_vector(6);  // maxdrift
  if (i == 6) return compute_vector(7);  // maxbondlen
  if (i == 7) return compute_vector(5);  // fraction with zero bias

  error->all(FLERR,"Invalid query to fix hyper/global");

  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixHyperGlobal::memory_usage()
{
  double bytes = maxbond * sizeof(OneBond);    // blist
  return bytes;
}
