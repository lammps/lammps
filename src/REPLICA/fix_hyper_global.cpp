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

#include "fix_hyper_global.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define DELTABOND 16384
#define VECLEN 5

// possible enhancements
//   should there be a virial contribution from boosted bond?
//   allow newton off?

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
  size_vector = 12;
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
  old2now = NULL;

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
  memory->destroy(old2now);
}

/* ---------------------------------------------------------------------- */

int FixHyperGlobal::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
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
  negstrain = 0;
}

/* ---------------------------------------------------------------------- */

void FixHyperGlobal::init()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Hyper global requires newton pair on");

  if (atom->molecular && me == 0)
    error->warning(FLERR,"Hyper global for molecular systems "
                   "requires care in defining hyperdynamic bonds");

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
  // no increment in these quantities when pre-run forces are calculated

  int nobias_hold = nobias;
  int negstrain_hold = negstrain;
  double t_hyper_hold = t_hyper;

  pre_reverse(eflag,vflag);

  nobias = nobias_hold;
  negstrain = negstrain_hold;
  t_hyper = t_hyper_hold;
}

/* ---------------------------------------------------------------------- */

void FixHyperGlobal::pre_neighbor()
{
  int i,m,iold,jold,ilocal,jlocal;
  // double distsq;

  // reset local indices for owned bond atoms, since atoms have migrated
  // must be done after ghost atoms are setup via comm->borders()
  // first time this is done for a particular I or J atom:
  //   use tagold and xold from when bonds were created
  //   atom->map() finds atom ID if it exists, owned index if possible
  //   closest current I or J atoms to old I may now be ghost atoms
  //   closest_image() returns the ghost atom index in that case
  // also compute max drift of any atom in a bond
  //   drift = displacement from quenched coord while event has not yet occured
  // NOTE: drift calc is now done in bond_build(), between 2 quenched states

  for (i = 0; i < nall_old; i++) old2now[i] = -1;

  for (m = 0; m < nblocal; m++) {
    iold = blist[m].iold;
    jold = blist[m].jold;
    ilocal = old2now[iold];
    jlocal = old2now[jold];

    if (ilocal < 0) {
      ilocal = atom->map(tagold[iold]);
      ilocal = domain->closest_image(xold[iold],ilocal);
      if (ilocal < 0)
        error->one(FLERR,"Fix hyper/global bond atom not found");
      old2now[iold] = ilocal;
      //distsq = MathExtra::distsq3(x[ilocal],xold[iold]);
      //maxdriftsq = MAX(distsq,maxdriftsq);
    }
    if (jlocal < 0) {
      jlocal = atom->map(tagold[jold]);
      jlocal = domain->closest_image(xold[iold],jlocal);   // closest to iold
      if (jlocal < 0)
        error->one(FLERR,"Fix hyper/global bond atom not found");
      old2now[jold] = jlocal;
      //distsq = MathExtra::distsq3(x[jlocal],xold[jold]);
      //maxdriftsq = MAX(distsq,maxdriftsq);
    }

    blist[m].i = ilocal;
    blist[m].j = jlocal;
  }
}

/* ---------------------------------------------------------------------- */

void FixHyperGlobal::pre_reverse(int /* eflag */, int /* vflag */)
{
  int i,j,m,imax,jmax;
  double delx,dely,delz;
  double r,r0,estrain,rmax,r0max,dt_boost;
  double ebias,vbias,fbias,fbiasr;

  // compute current strain of each owned bond
  // emax = maximum abs value of strain of any bond I own
  // imax,jmax = local indices of my 2 atoms in that bond
  // rmax,r0max = current and relaxed lengths of that bond

  double **x = atom->x;
  double emax = 0.0;

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

  // I own the bond with max absolute value of strain
  // compute bias force on atoms imax,jmax if strain < q, else zero
  // Ebias = current strain = (r-r0) / r0
  // Vbias = bias potential = Vmax (1 - Ebias^2/q^2)
  // Fbias = bias force as function of strain
  //       = -dVbias/dEbias = 2 Vmax Ebias / q^2
  // Fix = x component of force on atom I
  //     = Fbias dEbias/dr dr/dxi, dEbias/dr = 1/r0, dr/dxi = delx/r
  // dt_boost = time boost factor = exp(Vbias/kT)

  double **f = atom->f;

  vbias = fbias = 0.0;
  dt_boost = 1.0;

  if (emax < qfactor) {
    ebias = (rmax-r0max) / r0max;
    vbias = vmax * (1.0 - ebias*ebias*invqfactorsq);
    fbias = 2.0 * vmax * ebias * invqfactorsq;
    dt_boost = exp(beta*vbias);

    delx = x[imax][0] - x[jmax][0];
    dely = x[imax][1] - x[jmax][1];
    delz = x[imax][2] - x[jmax][2];

    fbiasr = fbias / r0max / rmax;
    f[imax][0] += delx*fbiasr;
    f[imax][1] += dely*fbiasr;
    f[imax][2] += delz*fbiasr;

    f[jmax][0] -= delx*fbiasr;
    f[jmax][1] -= dely*fbiasr;
    f[jmax][2] -= delz*fbiasr;

    if (ebias < 0.0) negstrain++;

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
  int i,j,m,ii,jj,iold,jold,ilocal,jlocal,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,distsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  if (natom) {
    nevent++;
    nevent_atom += natom;
  }

  // compute max distance any bond atom has moved between 2 quenched states
  // xold[iold] = last quenched coord for iold
  // x[ilocal] = current quenched coord for same atom

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  for (m = 0; m < nblocal; m++) {
    iold = blist[m].iold;
    ilocal = atom->map(tagold[iold]);
    ilocal = domain->closest_image(xold[iold],ilocal);
    distsq = MathExtra::distsq3(x[ilocal],xold[iold]);
    maxdriftsq = MAX(distsq,maxdriftsq);
    jold = blist[m].jold;
    jlocal = atom->map(tagold[jold]);
    jlocal = domain->closest_image(xold[iold],jlocal);
    distsq = MathExtra::distsq3(x[jlocal],xold[jold]);
    maxdriftsq = MAX(distsq,maxdriftsq);
  }

  // trigger neighbor list build

  neighbor->build_one(list);

  // identify bonds assigned to each owned atom
  // do not create a bond between two non-group atoms

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

  // store per-atom quantities for owned+ghost atoms at time of bond creation
  // nall_old = value of nall at time bonds are built

  tagint *tag = atom->tag;

  if (nall > maxold) {
    memory->destroy(xold);
    memory->destroy(tagold);
    memory->destroy(old2now);
    maxold = atom->nmax;
    memory->create(xold,maxold,3,"hyper/global:xold");
    memory->create(tagold,maxold,"hyper/global:tagold");
    memory->create(old2now,maxold,"hyper/global:old2now");
  }

  memcpy(&xold[0][0],&x[0][0],3*nall*sizeof(double));
  for (i = 0; i < nall; i++) tagold[i] = tag[i];

  nlocal_old = nlocal;
  nall_old = nall;
}

/* ----------------------------------------------------------------------
   grow bond list by a chunk
------------------------------------------------------------------------- */

void FixHyperGlobal::grow_bond()
{
  if (maxbond + DELTABOND > MAXSMALLINT)
    error->one(FLERR,"Fix hyper/global bond count is too big");
  maxbond += DELTABOND;
  blist = (OneBond *)
    memory->srealloc(blist,maxbond*sizeof(OneBond),"hyper/global:blist");
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

  // 12 vector outputs returned for i = 0-11

  // i = 0 = boost factor on this step
  // i = 1 = max strain of any bond on this step (positive or negative)
  // i = 2 = ID of atom I in max-strain bond on this step
  // i = 3 = ID of atom J in max-strain bond on this step
  // i = 4 = ave bonds/atom on this step

  // i = 5 = fraction of steps where bond has no bias during this run
  // i = 6 = fraction of steps where bond has negative strain during this run
  // i = 7 = max drift distance of any atom during this run
  // i = 8 = max bond length during this run

  // i = 9 = cummulative hyper time since fix created
  // i = 10 = cummulative # of event timesteps since fix created
  // i = 11 = cummulative # of atoms in events since fix created

  if (i == 0) return outvec[1];
  if (i == 1) return outvec[2];
  if (i == 2) return outvec[3];
  if (i == 3) return outvec[4];

  if (i == 4) {
    bigint mybonds = nblocal;
    bigint allbonds;
    MPI_Allreduce(&mybonds,&allbonds,1,MPI_LMP_BIGINT,MPI_SUM,world);
    return 2.0*allbonds/atom->natoms;
  }

  if (i == 5) {
    if (update->ntimestep == update->firststep) return 0.0;
    int allnobias;
    MPI_Allreduce(&nobias,&allnobias,1,MPI_INT,MPI_SUM,world);
    return 1.0*allnobias / (update->ntimestep - update->firststep);
  }

  if (i == 6) {
    if (update->ntimestep == update->firststep) return 0.0;
    int allnegstrain;
    MPI_Allreduce(&negstrain,&allnegstrain,1,MPI_INT,MPI_SUM,world);
    return 1.0*allnegstrain / (update->ntimestep - update->firststep);
  }

  if (i == 7) {
    double alldriftsq;
    MPI_Allreduce(&maxdriftsq,&alldriftsq,1,MPI_DOUBLE,MPI_MAX,world);
    return sqrt(alldriftsq);
  }

  if (i == 8) {
    double allbondlen;
    MPI_Allreduce(&maxbondlen,&allbondlen,1,MPI_DOUBLE,MPI_MAX,world);
    return allbondlen;
  }

  if (i == 9) return t_hyper;
  if (i == 10) return (double) nevent;
  if (i == 11) return (double) nevent_atom;

  return 0.0;
}

/* ----------------------------------------------------------------------
   wrapper on compute_vector()
   used by hyper.cpp to call FixHyper
------------------------------------------------------------------------- */

double FixHyperGlobal::query(int i)
{
  if (i == 1) return compute_vector(9);  // cummulative hyper time
  if (i == 2) return compute_vector(10); // nevent
  if (i == 3) return compute_vector(11); // nevent_atom
  if (i == 4) return compute_vector(4);  // ave bonds/atom
  if (i == 5) return compute_vector(7);  // maxdrift
  if (i == 6) return compute_vector(8);  // maxbondlen
  if (i == 7) return compute_vector(5);  // fraction with zero bias
  if (i == 8) return compute_vector(6);  // fraction with negative strain

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
