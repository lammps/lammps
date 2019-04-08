 /* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU Gene<ral Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include "fix_hyper_local.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
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

#define DELTABOND 16384
#define DELTABIAS 16
#define COEFFINIT 1.0
#define FCCBONDS 12
#define BIG 1.0e20

enum{STRAIN,STRAINDOMAIN,BIASFLAG,BIASCOEFF};
enum{IGNORE,WARN,ERROR};

/* ---------------------------------------------------------------------- */

FixHyperLocal::FixHyperLocal(LAMMPS *lmp, int narg, char **arg) :
  FixHyper(lmp, narg, arg), blist(NULL), biascoeff(NULL), numbond(NULL),
  maxhalf(NULL), eligible(NULL), maxhalfstrain(NULL), old2now(NULL),
  tagold(NULL), xold(NULL), maxstrain(NULL), maxstrain_domain(NULL),
  biasflag(NULL), bias(NULL), cpage(NULL), clist(NULL), numcoeff(NULL)
{
  // error checks

  if (atom->map_style == 0)
    error->all(FLERR,"Fix hyper/local command requires atom map");

  // parse args

  if (narg < 10) error->all(FLERR,"Illegal fix hyper/local command");

  hyperflag = 2;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 21;
  local_flag = 1;
  size_local_rows = 0;
  size_local_cols = 0;
  local_freq = 1;

  global_freq = 1;
  extscalar = 0;
  extvector = 0;

  cutbond = force->numeric(FLERR,arg[3]);
  qfactor = force->numeric(FLERR,arg[4]);
  vmax = force->numeric(FLERR,arg[5]);
  tequil = force->numeric(FLERR,arg[6]);
  dcut = force->numeric(FLERR,arg[7]);
  alpha_user = force->numeric(FLERR,arg[8]);
  boost_target = force->numeric(FLERR,arg[9]);

  if (cutbond < 0.0 || qfactor < 0.0 || vmax < 0.0 ||
      tequil <= 0.0 || dcut <= 0.0 || alpha_user <= 0.0 || boost_target < 1.0)
    error->all(FLERR,"Illegal fix hyper/local command");

  invqfactorsq = 1.0 / (qfactor*qfactor);
  cutbondsq = cutbond*cutbond;
  dcutsq = dcut*dcut;
  beta = 1.0 / (force->boltz * tequil);

  // optional args

  checkghost = 0;
  checkbias = 0;

  int iarg = 10;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"check/ghost") == 0) {
      checkghost = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"check/bias") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix hyper/local command");
      checkbias = 1;
      checkbias_every = force->inumeric(FLERR,arg[iarg+1]);
      if (strcmp(arg[iarg+2],"error") == 0) checkbias_flag = ERROR;
      else if (strcmp(arg[iarg+2],"warn") == 0) checkbias_flag = WARN;
      else if (strcmp(arg[iarg+2],"ignore") == 0) checkbias_flag = IGNORE;
      else error->all(FLERR,"Illegal fix hyper/local command");
      iarg += 3;
    } else error->all(FLERR,"Illegal fix hyper/local command");
  }

  // per-atom data structs

  maxbond = nblocal = 0;
  blist = NULL;
  biascoeff = NULL;
  allbonds = 0;

  maxatom = 0;
  maxstrain = NULL;
  maxstrain_domain = NULL;
  biasflag = NULL;

  maxlocal = nlocal_old = 0;
  numbond = NULL;
  maxhalf = NULL;
  eligible = NULL;
  maxhalfstrain = NULL;

  maxall = nall_old = 0;
  xold = NULL;
  tagold = NULL;
  old2now = NULL;

  nbias = maxbias = 0;
  bias = NULL;

  // data structs for persisting bias coeffs when bond list is reformed
  // maxbondperatom = max # of bonds any atom is part of
  // FCCBONDS = 12 is a good estimate for fcc lattices
  // will be reset in build_bond() if necessary

  maxcoeff = 0;
  maxbondperatom = FCCBONDS;
  numcoeff = NULL;
  clist = NULL;
  cpage = new MyPage<OneCoeff>;
  cpage->init(maxbondperatom,1024*maxbondperatom,1);

  // set comm sizes needed by this fix
  // reverse = 2 is for sending atom index + value, though total likely < 1
  // reverse comm for bias coeffs has variable size, so not tallied here

  comm_forward = 1;
  comm_reverse = 2;

  me = comm->me;
  firstflag = 1;

  sumbiascoeff = 0.0;
  avebiascoeff = 0.0;

  starttime = update->ntimestep;
  nostrainyet = 1;
  nnewbond = 0;
  nevent = 0;
  nevent_atom = 0;
  mybias = 0.0;
}

/* ---------------------------------------------------------------------- */

FixHyperLocal::~FixHyperLocal()
{
  memory->destroy(blist);
  memory->destroy(biascoeff);

  memory->destroy(maxstrain);
  memory->destroy(maxstrain_domain);
  memory->destroy(biasflag);

  memory->destroy(numbond);
  memory->destroy(maxhalf);
  memory->destroy(eligible);
  memory->destroy(maxhalfstrain);

  memory->destroy(xold);
  memory->destroy(tagold);
  memory->destroy(old2now);

  memory->destroy(bias);

  memory->destroy(numcoeff);
  memory->sfree(clist);
  delete cpage;
}

/* ---------------------------------------------------------------------- */

int FixHyperLocal::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= PRE_REVERSE;
  mask |= MIN_PRE_NEIGHBOR;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::init_hyper()
{
  ghost_toofar = 0;
  checkbias_count = 0;
  maxdriftsq = 0.0;
  maxbondlen = 0.0;
  avebiascoeff = 0.0;
  minbiascoeff = BIG;
  maxbiascoeff = 0.0;
  nbias_running = 0;
  nobias_running = 0;
  negstrain_running = 0;
  rmaxever = 0.0;
  rmaxeverbig = 0.0;

  nbondbuild = 0;
  time_bondbuild = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::init()
{
  // for newton off, bond force bias will not be applied correctly
  //   for bonds that straddle 2 procs
  // warn if molecular system, since near-neighbors may not appear in neigh list
  //   user should not be including bonded atoms as hyper "bonds"

  if (force->newton_pair == 0)
    error->all(FLERR,"Hyper local requires newton pair on");

  if (atom->molecular && me == 0)
    error->warning(FLERR,"Hyper local for molecular systems "
                   "requires care in defining hyperdynamic bonds");

  // cutghost = communication cutoff as calculated by Neighbor and Comm
  // error if cutghost is smaller than Dcut
  // warn if no drift distance added to cutghost

  if (firstflag) {
    double cutghost;
    if (force->pair)
      cutghost = MAX(force->pair->cutforce+neighbor->skin,comm->cutghostuser);
    else
      cutghost = comm->cutghostuser;

    if (cutghost < dcut)
      error->all(FLERR,"Fix hyper/local domain cutoff exceeds ghost atom range - "
                 "use comm_modify cutoff command");
    if (cutghost < dcut+cutbond/2.0 && me == 0)
      error->warning(FLERR,"Fix hyper/local ghost atom range "
                     "may not allow for atom drift between events");
  }

  alpha = update->dt / alpha_user;

  // need occasional full neighbor list with cutoff = Dcut
  // used for finding maxstrain of neighbor bonds out to Dcut
  // do not need to include neigh skin in cutoff,
  //   b/c this list will be built every time build_bond() is called
  // NOTE: what if pair style list cutoff > Dcut
  //   or what if neigh skin is huge?

  int irequest_full = neighbor->request(this,instance_me);
  neighbor->requests[irequest_full]->id = 1;
  neighbor->requests[irequest_full]->pair = 0;
  neighbor->requests[irequest_full]->fix = 1;
  neighbor->requests[irequest_full]->half = 0;
  neighbor->requests[irequest_full]->full = 1;
  neighbor->requests[irequest_full]->cut = 1;
  neighbor->requests[irequest_full]->cutoff = dcut;
  neighbor->requests[irequest_full]->occasional = 1;

  // also need occasional half neighbor list derived from pair style
  // used for building local bond list
  // no specified cutoff, should be longer than cutbond
  // this list will also be built (or derived/copied)
  //   every time build_bond() is called

  int irequest_half = neighbor->request(this,instance_me);
  neighbor->requests[irequest_half]->id = 2;
  neighbor->requests[irequest_half]->pair = 0;
  neighbor->requests[irequest_half]->fix = 1;
  neighbor->requests[irequest_half]->occasional = 1;

  // extra timing output

  //timefirst = timesecond = timethird = timefourth = timefifth =
  //  timesixth = timeseventh = timetotal = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::init_list(int id, NeighList *ptr)
{
  if (id == 1) listfull = ptr;
  else if (id == 2) listhalf = ptr;
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::setup_pre_neighbor()
{
  // called for dynamics and minimization

  pre_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::setup_pre_reverse(int eflag, int vflag)
{
  // only called for dynamics, not minimization
  // setupflag prevents boostostat update of bias coeffs in setup
  // also prevents increments of nbias_running, nobias_running,
  //   negstrain_running, sumbiascoeff

  setupflag = 1;
  pre_reverse(eflag,vflag);
  setupflag = 0;
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::pre_neighbor()
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
        error->one(FLERR,"Fix hyper/local bond atom not found");
      old2now[iold] = ilocal;
      //distsq = MathExtra::distsq3(x[ilocal],xold[iold]);
      //maxdriftsq = MAX(distsq,maxdriftsq);
    }
    if (jlocal < 0) {
      jlocal = atom->map(tagold[jold]);
      jlocal = domain->closest_image(xold[iold],jlocal);   // close to I atom
      if (jlocal < 0)
        error->one(FLERR,"Fix hyper/local bond atom not found");
      old2now[jold] = jlocal;
      //distsq = MathExtra::distsq3(x[jlocal],xold[jold]);
      //maxdriftsq = MAX(distsq,maxdriftsq);
    }

    blist[m].i = ilocal;
    blist[m].j = jlocal;
  }

  // set remaining old2now values to point to current local atom indices
  // if old2now >= 0, already set by bond loop above
  // only necessary for tagold entries > 0
  //   because if tagold = 0, atom is not active in Dcut neighbor list
  // must be done after atoms migrate and ghost atoms setup via comm->borders()
  // does not matter which atom (owned or ghost) that atom->map() finds
  //   b/c old2now is only used to access maxstrain() or biasflag()
  //   which will be identical for every copy of the same atom ID

  for (iold = 0; iold < nall_old; iold++) {
    if (old2now[iold] >= 0) continue;
    if (tagold[iold] == 0) continue;
    ilocal = atom->map(tagold[iold]);
    old2now[iold] = ilocal;
    if (ilocal < 0) ghost_toofar++;
  }
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::pre_reverse(int /* eflag */, int /* vflag */)
{
  int i,j,m,ii,jj,inum,jnum,iold,jold,ibond,nbond,ijhalf,ncount;
  double delx,dely,delz;
  double r,r0,estrain,emax,ebias,vbias,fbias,fbiasr;
  double halfstrain,selfstrain;
  int *ilist,*jlist,*numneigh,**firstneigh;

  //double time1,time2,time3,time4,time5,time6,time7,time8;
  //time1 = MPI_Wtime();

  nostrainyet = 0;

  // reallocate per-atom maxstrain and biasflag vectors if necessary

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (maxatom < nall) {
    memory->destroy(maxstrain);
    memory->destroy(maxstrain_domain);
    if (checkbias) memory->destroy(biasflag);
    maxatom = atom->nmax;
    memory->create(maxstrain,maxatom,"hyper/local:maxstrain");
    memory->create(maxstrain_domain,maxatom,"hyper/local:maxstrain_domain");
    if (checkbias) memory->create(biasflag,maxatom,"hyper/local:biasflag");
  }

  // one max strain bond per old owned atom is eligible for biasing

  for (iold = 0; iold < nlocal_old; iold++) eligible[iold] = 1;

  // -------------------------------------------------------------
  // stage 1:
  // maxstrain[i] = max abs value of strain of any bond atom I is part of
  // reverse/forward comm so know it for all current owned and ghost atoms
  // -------------------------------------------------------------

  // compute estrain = current abs value strain of each owned bond
  // blist = bondlist from last event
  // mark atom I ineligible if it has no bonds
  // also store:
  //   maxhalf = which owned bond is maxstrain for each old atom I
  //   maxhalfstrain = abs value strain of that bond for each old atom I

  for (i = 0; i < nall; i++) maxstrain[i] = 0.0;

  double **x = atom->x;

  m = 0;
  for (iold = 0; iold < nlocal_old; iold++) {
    nbond = numbond[iold];
    if (!nbond) {
      eligible[iold] = 0;
      continue;
    }
    halfstrain = 0.0;
    for (ibond = 0; ibond < nbond; ibond++) {
      i = blist[m].i;
      j = blist[m].j;
      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      r = sqrt(delx*delx + dely*dely + delz*delz);
      maxbondlen = MAX(r,maxbondlen);
      r0 = blist[m].r0;
      estrain = fabs(r-r0) / r0;
      maxstrain[i] = MAX(maxstrain[i],estrain);
      maxstrain[j] = MAX(maxstrain[j],estrain);
      if (estrain > halfstrain) {
        halfstrain = estrain;
        ijhalf = m;
      }
      m++;
    }
    maxhalf[iold] = ijhalf;
    maxhalfstrain[iold] = halfstrain;
  }

  //time2 = MPI_Wtime();

  // reverse comm acquires maxstrain of all current owned atoms
  //   needed b/c only saw half the bonds of each atom
  //   also needed b/c bond list may refer to old owned atoms that are now ghost
  // forward comm acquires maxstrain of all current ghost atoms

  commflag = STRAIN;
  comm->reverse_comm_fix(this);
  comm->forward_comm_fix(this);

  //time3 = MPI_Wtime();

  // -------------------------------------------------------------
  // stage 2:
  // maxstrain_domain[i] = maxstrain of atom I and all its J neighs out to Dcut
  // reverse/forward comm so know it for all current owned and ghost atoms
  // -------------------------------------------------------------

  // use full Dcut neighbor list to check maxstrain of all neighbor atoms
  // neighlist is from last event
  //   has old indices for I,J (reneighboring may have occurred)
  //   use old2now[] to convert to current indices
  //   if J is unknown (drifted ghost),
  //     assume it was part of an event and its strain = qfactor
  // mark atom I ineligible for biasing if:
  //   its maxhalfstrain < maxstrain (J atom owns the IJ bond)
  //   its maxstrain < maxstrain_domain
  //   ncount > 1 (break tie by making all atoms with tie value ineligible)
  // if ncount > 1, also flip sign of maxstrain_domain for atom I

  for (i = 0; i < nall; i++) maxstrain_domain[i] = 0.0;

  inum = listfull->inum;
  ilist = listfull->ilist;
  numneigh = listfull->numneigh;
  firstneigh = listfull->firstneigh;

  double rmax = rmaxever;
  double rmaxbig = rmaxeverbig;
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  // first two lines of outer loop should be identical to this:
  // for (iold = 0; iold < nlocal_old; iold++)

  for (ii = 0; ii < inum; ii++) {
    iold = ilist[ii];
    if (eligible[iold] == 0) continue;
    jlist = firstneigh[iold];
    jnum = numneigh[iold];

    // I or J may be ghost atoms
    // will always know I b/c atoms do not drift that far
    // but may no longer know J if hops outside cutghost
    // in that case, assume it performed an event, its strain = qfactor
    // this assumes cutghost is sufficiently longer than Dcut

    i = old2now[iold];
    emax = selfstrain = maxstrain[i];
    ncount = 0;

    for (jj = 0; jj < jnum; jj++) {
      jold = jlist[jj];
      j = old2now[jold];

      // special case for missing (drifted) J atom

      if (j < 0) {
        emax = MAX(emax,qfactor);
        if (selfstrain == qfactor) ncount++;
        continue;
      }

      emax = MAX(emax,maxstrain[j]);
      if (selfstrain == maxstrain[j]) ncount++;

      // optional diagnostic
      // tally largest distance from subbox that a ghost atom is (rmaxbig)
      // and the largest distance if strain < qfactor (rmax)

      if (checkghost) {
        if (j >= nlocal) {
          if (x[j][0] < sublo[0]) rmaxbig = MAX(rmaxbig,sublo[0]-x[j][0]);
          if (x[j][1] < sublo[1]) rmaxbig = MAX(rmaxbig,sublo[1]-x[j][1]);
          if (x[j][2] < sublo[2]) rmaxbig = MAX(rmaxbig,sublo[2]-x[j][2]);
          if (x[j][0] > subhi[0]) rmaxbig = MAX(rmaxbig,x[j][0]-subhi[0]);
          if (x[j][1] > subhi[1]) rmaxbig = MAX(rmaxbig,x[j][1]-subhi[1]);
          if (x[j][2] > subhi[2]) rmaxbig = MAX(rmaxbig,x[j][2]-subhi[2]);
          if (maxstrain[j] < qfactor) {
            if (x[j][0] < sublo[0]) rmax = MAX(rmax,sublo[0]-x[j][0]);
            if (x[j][1] < sublo[1]) rmax = MAX(rmax,sublo[1]-x[j][1]);
            if (x[j][2] < sublo[2]) rmax = MAX(rmax,sublo[2]-x[j][2]);
            if (x[j][0] > subhi[0]) rmax = MAX(rmax,x[j][0]-subhi[0]);
            if (x[j][1] > subhi[1]) rmax = MAX(rmax,x[j][1]-subhi[1]);
            if (x[j][2] > subhi[2]) rmax = MAX(rmax,x[j][2]-subhi[2]);
          }
        }
      }
    }

    if (maxhalfstrain[iold] < selfstrain) eligible[iold] = 0;
    if (selfstrain < emax) eligible[iold] = 0;
    else if (ncount > 1) {
      eligible[iold] = 0;
      emax = -emax;
    }
    maxstrain_domain[i] = emax;
  }

  //time4 = MPI_Wtime();

  // reverse comm to acquire maxstrain_domain from ghost atoms
  //   needed b/c neigh list may refer to old owned atoms that are now ghost
  // forward comm acquires maxstrain_domain of all current ghost atoms

  commflag = STRAINDOMAIN;
  comm->reverse_comm_fix(this);
  comm->forward_comm_fix(this);

  //time5 = MPI_Wtime();

  // -------------------------------------------------------------
  // stage 3:
  // create bias = list of Nbias biased bonds this proc owns
  // -------------------------------------------------------------

  // identify biased bonds and add to bias list
  // bias the I,J maxhalf bond of atom I only if all these conditions hold:
  //   maxstrain[i] = maxstrain_domain[i] (checked in stage 2)
  //   maxstrain[j] = maxstrain_domain[j] (checked here)
  //   I is not part of an I,J bond with > strain owned by some J (checked in 2)
  //   no ties with other maxstrain bonds in atom I's domain (chedcked in 2)

  nbias = 0;
  for (iold = 0; iold < nlocal_old; iold++) {
    if (eligible[iold] == 0) continue;
    j = blist[maxhalf[iold]].j;
    if (maxstrain[j] != maxstrain_domain[j]) continue;
    if (nbias == maxbias) {
      maxbias += DELTABIAS;
      memory->grow(bias,maxbias,"hyper/local:bias");
    }
    bias[nbias++] = maxhalf[iold];
  }

  //time6 = MPI_Wtime();

  // -------------------------------------------------------------
  // stage 4:
  // apply bias force to bonds with locally max strain
  // -------------------------------------------------------------

  double **f = atom->f;

  int nobias = 0;
  int negstrain = 0;
  mybias = 0.0;

  for (int ibias = 0; ibias < nbias; ibias++) {
    m = bias[ibias];
    i = blist[m].i;
    j = blist[m].j;

    if (maxstrain[i] >= qfactor) {
      nobias++;
      continue;
    }

    delx = x[i][0] - x[j][0];
    dely = x[i][1] - x[j][1];
    delz = x[i][2] - x[j][2];
    r = sqrt(delx*delx + dely*dely + delz*delz);
    r0 = blist[m].r0;
    ebias = (r-r0) / r0;
    vbias = biascoeff[m] * vmax * (1.0 - ebias*ebias*invqfactorsq);
    fbias = biascoeff[m] * 2.0 * vmax * ebias * invqfactorsq;

    fbiasr = fbias / r0 / r;
    f[i][0] += delx*fbiasr;
    f[i][1] += dely*fbiasr;
    f[i][2] += delz*fbiasr;

    f[j][0] -= delx*fbiasr;
    f[j][1] -= dely*fbiasr;
    f[j][2] -= delz*fbiasr;

    if (ebias < 0.0) negstrain++;
    mybias += vbias;
  }

  //time7 = MPI_Wtime();

  // -------------------------------------------------------------
  // apply boostostat to bias coeffs of all bonds I own
  // -------------------------------------------------------------

  // no boostostat update when pre_reverse called from setup()
  // nbias_running, nobias_running, negstrain_running only incremented
  //   on run steps

  if (setupflag) return;
  nbias_running += nbias;
  nobias_running += nobias;
  negstrain_running += negstrain;

  // loop over bonds I own to adjust bias coeff
  // delta in boost coeff is function of maxboost_domain vs target boost
  // maxboost_domain is function of two maxstrain_domains for I,J

  double emaxi,emaxj,maxboost_domain,bc;
  double mybiascoeff = 0.0;

  for (m = 0; m < nblocal; m++) {
    i = blist[m].i;
    j = blist[m].j;
    emaxi = maxstrain_domain[i];
    emaxj = maxstrain_domain[j];
    emax = MAX(emaxi,emaxj);
    if (emax < qfactor) vbias = vmax * (1.0 - emax*emax*invqfactorsq);
    else vbias = 0.0;

    maxboost_domain = exp(beta * biascoeff[m]*vbias);
    biascoeff[m] -= alpha * (maxboost_domain-boost_target) / boost_target;

    // stats

    bc = biascoeff[m];
    mybiascoeff += bc;
    minbiascoeff = MIN(minbiascoeff,bc);
    maxbiascoeff = MAX(maxbiascoeff,bc);
  }

  // -------------------------------------------------------------
  // diagnostics, some optional
  // -------------------------------------------------------------

  MPI_Allreduce(&mybiascoeff,&sumbiascoeff,1,MPI_DOUBLE,MPI_SUM,world);
  if (allbonds) avebiascoeff += sumbiascoeff/allbonds;

  // if requested, monitor ghost distance from processor sub-boxes

  if (checkghost) {
    double rmax2[2],rmax2all[2];
    rmax2[0] = rmax;
    rmax2[1] = rmaxbig;
    MPI_Allreduce(&rmax2,&rmax2all,2,MPI_DOUBLE,MPI_MAX,world);
    rmaxever = rmax2all[0];
    rmaxeverbig = rmax2all[1];
  }

  // if requsted, check for any biased bonds that are too close to each other
  // keep a running count for output
  // requires 2 additional local comm operations

  if (checkbias && update->ntimestep % checkbias_every == 0) {

    // mark each atom in a biased bond with ID of partner
    // this may mark some ghost atoms

    for (i = 0; i < nall; i++) biasflag[i] = 0;

    tagint *tag = atom->tag;

    for (int ibias = 0; ibias < nbias; ibias++) {
      m = bias[ibias];
      i = blist[m].i;
      j = blist[m].j;
      biasflag[i] = tag[j];
      biasflag[j] = tag[i];
    }

    // reverse comm to acquire biasflag from ghost atoms
    // forward comm to set biasflag for all ghost atoms

    commflag = BIASFLAG;
    comm->reverse_comm_fix(this);
    comm->forward_comm_fix(this);

    // loop over Dcut full neighbor list
    // I and J may be ghost atoms
    // only continue if I is a biased atom
    // if J is unknown (drifted ghost) just ignore
    // if J is biased and is not bonded to I, then flag as too close

    for (ii = 0; ii < inum; ii++) {
      iold = ilist[ii];
      i = old2now[iold];
      if (biasflag[i] == 0) continue;

      jlist = firstneigh[iold];
      jnum = numneigh[iold];

      for (jj = 0; jj < jnum; jj++) {
        jold = jlist[jj];
        j = old2now[jold];
        if (j < 0) continue;
        if (biasflag[j] && biasflag[j] != tag[i]) checkbias_count++;
      }
    }

    if (checkbias_flag != IGNORE) {
      int allcount;
      MPI_Allreduce(&checkbias_count,&allcount,1,MPI_INT,MPI_SUM,world);
      if (allcount) {
        char str[128];
        sprintf(str,"Fix hyper/local biased bonds too close: "
                "cumulative atom count %d",allcount);
        if (checkbias_flag == WARN) {
          if (me == 0) error->warning(FLERR,str);
        } else error->all(FLERR,str);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::min_pre_neighbor()
{
  pre_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::build_bond_list(int natom)
{
  int i,j,ii,jj,m,n,iold,jold,ilocal,jlocal,inum,jnum,nbond;
  tagint jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,distsq,oldcoeff;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double time1,time2;
  time1 = MPI_Wtime();

  if (natom) {
    nevent++;
    nevent_atom += natom;
  }

  // compute max distance any bond atom has moved between 2 quenched states
  // xold[iold] = last quenched coord for iold
  // x[ilocal] = current quenched coord for same atom
  // use of old2now calculates distsq only once per atom

  double **x = atom->x;

  for (i = 0; i < nall_old; i++) old2now[i] = -1;

  for (m = 0; m < nblocal; m++) {
    iold = blist[m].iold;
    if (old2now[iold] < 0) {
      ilocal = atom->map(tagold[iold]);
      ilocal = domain->closest_image(xold[iold],ilocal);
      if (ilocal < 0) error->one(FLERR,"Fix hyper/local bond atom not found");
      old2now[iold] = ilocal;
      distsq = MathExtra::distsq3(x[ilocal],xold[iold]);
      maxdriftsq = MAX(distsq,maxdriftsq);
    }
    jold = blist[m].jold;
    if (old2now[jold] < 0) {
      jold = blist[m].jold;
      jlocal = atom->map(tagold[jold]);
      jlocal = domain->closest_image(xold[iold],jlocal);  // close to I atom
      if (jlocal < 0) error->one(FLERR,"Fix hyper/local bond atom not found");
      old2now[jold] = jlocal;
      distsq = MathExtra::distsq3(x[jlocal],xold[jold]);
      maxdriftsq = MAX(distsq,maxdriftsq);
    }
  }

  // store old bond coeffs so can persist them in new blist
  // while loop allows growing value of maxbondperatom
  // will loop at most 2 times, stops when maxbondperatom is large enough
  // requires reverse comm, no forward comm:
  //    b/c new coeff list is stored only by current owned atoms

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (maxcoeff < nall) {
    memory->destroy(numcoeff);
    memory->sfree(clist);
    maxcoeff = atom->nmax;
    memory->create(numcoeff,maxcoeff,"hyper/local:numcoeff");
    clist = (OneCoeff **) memory->smalloc(maxcoeff*sizeof(OneCoeff *),
                                         "hyper/local:clist");
  }

  while (1) {
    if (firstflag) break;
    for (i = 0; i < nall; i++) numcoeff[i] = 0;
    for (i = 0; i < nall; i++) clist[i] = NULL;
    cpage->reset();

    for (m = 0; m < nblocal; m++) {
      i = blist[m].i;
      j = blist[m].j;

      if (numcoeff[i] == 0) clist[i] = cpage->get(maxbondperatom);
      if (numcoeff[j] == 0) clist[j] = cpage->get(maxbondperatom);

      if (numcoeff[i] < maxbondperatom) {
        clist[i][numcoeff[i]].biascoeff = biascoeff[m];
        clist[i][numcoeff[i]].tag = tag[j];
      }
      numcoeff[i]++;

      if (numcoeff[j] < maxbondperatom) {
        clist[j][numcoeff[j]].biascoeff = biascoeff[m];
        clist[j][numcoeff[j]].tag = tag[i];
      }
      numcoeff[j]++;
    }

    int mymax = 0;
    for (i = 0; i < nall; i++) mymax = MAX(mymax,numcoeff[i]);
    int maxcoeffall;
    MPI_Allreduce(&mymax,&maxcoeffall,1,MPI_INT,MPI_MAX,world);

    if (maxcoeffall > maxbondperatom) {
      maxbondperatom = maxcoeffall;
      cpage->init(maxbondperatom,1024*maxbondperatom,1);
      continue;
    }

    commflag = BIASCOEFF;
    comm->reverse_comm_fix_variable(this);

    mymax = 0;
    for (i = 0; i < nall; i++) mymax = MAX(mymax,numcoeff[i]);
    MPI_Allreduce(&mymax,&maxcoeffall,1,MPI_INT,MPI_MAX,world);
    if (maxcoeffall <= maxbondperatom) break;

    maxbondperatom = maxcoeffall;
    cpage->init(maxbondperatom,1024*maxbondperatom,1);
  }

  // reallocate vectors that are maxlocal and maxall length if necessary

  if (nlocal > maxlocal) {
    memory->destroy(eligible);
    memory->destroy(numbond);
    memory->destroy(maxhalf);
    memory->destroy(maxhalfstrain);
    maxlocal = nlocal;
    memory->create(eligible,maxlocal,"hyper/local:eligible");
    memory->create(numbond,maxlocal,"hyper/local:numbond");
    memory->create(maxhalf,maxlocal,"hyper/local:maxhalf");
    memory->create(maxhalfstrain,maxlocal,"hyper/local:maxhalfstrain");
  }

  if (nall > maxall) {
    memory->destroy(xold);
    memory->destroy(tagold);
    memory->destroy(old2now);
    maxall = atom->nmax;
    memory->create(xold,maxall,3,"hyper/local:xold");
    memory->create(tagold,maxall,"hyper/local:tagold");
    memory->create(old2now,maxall,"hyper/local:old2now");
  }

  // nlocal_old = value of nlocal at time bonds are built
  // nall_old = value of nall at time bonds are built
  // archive current atom coords in xold
  // tagold will be set to non-zero below for accessed atoms
  // numbond will be set below

  nlocal_old = nlocal;
  nall_old = nall;

  memcpy(&xold[0][0],&x[0][0],3*nall*sizeof(double));
  for (i = 0; i < nall; i++) tagold[i] = 0;
  for (i = 0; i < nlocal; i++) numbond[i] = 0;

  // trigger neighbor list builds for both lists
  // insure the I loops in both are from 1 to nlocal

  neighbor->build_one(listfull);
  neighbor->build_one(listhalf);

  if (listfull->inum != nlocal || listhalf->inum != nlocal)
    error->one(FLERR,"Invalid neighbor list in fix hyper/local bond build");

  // set tagold = 1 for all J atoms used in full neighbor list
  // tagold remains 0 for unused atoms, skipped in pre_neighbor

  inum = listfull->inum;
  ilist = listfull->ilist;
  numneigh = listfull->numneigh;
  firstneigh = listfull->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    tagold[i] = tag[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      tagold[j] = tag[j];
    }
  }

  // identify bonds assigned to each owned atom
  // do not create a bond between two non-group atoms

  int *mask = atom->mask;

  inum = listhalf->inum;
  ilist = listhalf->ilist;
  numneigh = listhalf->numneigh;
  firstneigh = listhalf->firstneigh;

  nblocal = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    tagold[i] = tag[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    nbond = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      jtag = tag[j];
      tagold[j] = jtag;

      // skip if neither atom I or J are in fix group

      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutbondsq) {
        nbond++;

        if (nblocal == maxbond) grow_bond();
        blist[nblocal].i = i;
        blist[nblocal].j = j;
        blist[nblocal].iold = i;
        blist[nblocal].jold = j;
        blist[nblocal].r0 = sqrt(rsq);

        // set biascoeff to old coeff for same I,J pair or to default

        if (firstflag) oldcoeff = 0.0;
        else {
          oldcoeff = 0.0;
          jtag = tag[j];
          n = numcoeff[i];
          for (m = 0; m < n; m++) {
            if (clist[i][m].tag == jtag) {
              oldcoeff = clist[i][m].biascoeff;
              break;
            }
          }
        }

        if (oldcoeff > 0.0) biascoeff[nblocal] = oldcoeff;
        else {
          biascoeff[nblocal] = COEFFINIT;
          nnewbond++;
        }

        nblocal++;
      }
    }

    numbond[i] = nbond;
  }

  // this fix allows access to biascoeffs as local data

  size_local_rows = nblocal;

  // allbonds = total # of bonds in system

  bigint bondcount = nblocal;
  MPI_Allreduce(&bondcount,&allbonds,1,MPI_LMP_BIGINT,MPI_SUM,world);

  time2 = MPI_Wtime();

  if (firstflag) nnewbond = 0;
  else {
    time_bondbuild += time2-time1;
    nbondbuild++;
  }
  firstflag = 0;
}

/* ---------------------------------------------------------------------- */

int FixHyperLocal::pack_forward_comm(int n, int *list, double *buf,
                                     int /* pbc_flag */, int * /* pbc */)
{
  int i,j,m;

  m = 0;

  // STRAIN
  // pack maxstrain vector
  // must send to all ghosts out to Dcut

  if (commflag == STRAIN) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = maxstrain[j];
    }

  // STRAINDOMAIN
  // pack maxstrain_domain vector
  // could just send to nearby ghosts in bonds
  // don't see easy way to determine precisely which atoms that is

  } else if (commflag == STRAINDOMAIN) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = maxstrain_domain[j];
    }

  // BIASFLAG
  // pack biasflag vector
  // must send to all ghosts out to Dcut

  } else if (commflag == BIASFLAG) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(biasflag[j]).d;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  // STRAIN
  // unpack maxstrain vector

  if (commflag == STRAIN) {
    for (i = first; i < last; i++) {
      maxstrain[i] = buf[m++];
    }

  // STRAINREGION
  // unpack maxstrain_domain vector

  } else if (commflag == STRAINDOMAIN) {
    for (i = first; i < last; i++) {
      maxstrain_domain[i] = buf[m++];
    }

  // BIASFLAG
  // unpack biasflag vector

  } else if (commflag == BIASFLAG) {
    for (i = first; i < last; i++) {
      biasflag[i] = (tagint) ubuf(buf[m++]).i;
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixHyperLocal::pack_reverse_comm(int n, int first, double *buf)
{
  int i,j,m,last;

  m = 0;
  last = first + n;

  // STRAIN
  // pack maxstrain vector
  // only pack for nonzero values

  if (commflag == STRAIN) {
    int nonzero = 0;
    m++;                           // placeholder for count of atoms
    for (i = first; i < last; i++) {
      if (maxstrain[i] == 0.0) continue;
      nonzero++;
      buf[m++] = ubuf(i-first).d;  // which atom is next
      buf[m++] = maxstrain[i];     // value
    }
    buf[0] = ubuf(nonzero).d;

  // STRAINDOMAIN
  // pack maxstrain_domain vector
  // only pack for nonzero values

  } else if (commflag == STRAINDOMAIN) {
    int nonzero = 0;
    m++;                           // placeholder for count of atoms
    for (i = first; i < last; i++) {
      if (maxstrain_domain[i] == 0.0) continue;
      nonzero++;
      buf[m++] = ubuf(i-first).d;      // which atom is next
      buf[m++] = maxstrain_domain[i];  // value
    }
    buf[0] = ubuf(nonzero).d;

  // BIASFLAG
  // pack biasflag vector
  // could just pack for nonzero values, like STRAIN and STRAINDOMAIN

  } else if (commflag == BIASFLAG) {
    for (i = first; i < last; i++) {
      buf[m++] = ubuf(biasflag[i]).d;
    }

  // BIASCOEFF
  // pack list of biascoeffs
  // only pack for atoms with nonzero # of bias coeffs
  // this will skip majority of ghost atoms

  } else if (commflag == BIASCOEFF) {
    int ncoeff;
    int nonzero = 0;
    m++;                             // placeholder for count of atoms
    for (i = first; i < last; i++) {
      if (numcoeff[i] == 0) continue;
      nonzero++;
      ncoeff = numcoeff[i];
      buf[m++] = ubuf(i-first).d;   // which atom is next
      buf[m++] = ubuf(ncoeff).d;    // # of bias coeffs
      for (j = 0; j < ncoeff; j++) {
        buf[m++] = clist[i][j].biascoeff;
        buf[m++] = ubuf(clist[i][j].tag).d;
      }
    }
    buf[0] = ubuf(nonzero).d;
  }

  return m;
}

/* ----------------------------------------------------------------------
   callback by comm->reverse_comm_fix_variable() in build_bond()
   same logic as BIASCOEFF option in pack_reverse_comm()
   m = returned size of message
------------------------------------------------------------------------- */

int FixHyperLocal::pack_reverse_comm_size(int n, int first)
{
  int last = first + n;
  int m = 1;
  for (int i = first; i < last; i++) {
    if (numcoeff[i]) m += 2 + 2*numcoeff[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,k,m;

  // return if n = 0
  // b/c if there are no atoms (n = 0), the message will not have
  //   been sent by Comm::reverse_comm_fix() or reverse_comm_fix_variable()
  // so must not read nonzero from first buf location (would be zero anyway)

  if (n == 0) return;

  // STRAIN
  // unpack maxstrain vector
  // nonzero # of entries, each has offset to which atom in receiver's list
  // use MAX, b/c want maximum abs value strain for each atom's bonds

  m = 0;

  if (commflag == STRAIN) {
    int offset;
    int nonzero = (int) ubuf(buf[m++]).i;   // # of atoms with values

    for (int iatom = 0; iatom < nonzero; iatom++) {
      offset = (int) ubuf(buf[m++]).i;      // offset into list for which atom
      j = list[offset];
      maxstrain[j] = MAX(maxstrain[j],buf[m]);
      m++;
    }

  // STRAINDOMAIN
  // unpack maxstrain_domain vector
  // use MAX, b/c want maximum abs value strain for each atom's domain
  // could also use SUM, b/c exactly one ghost or owned value is non-zero

  } else if (commflag == STRAINDOMAIN) {
    int offset;
    int nonzero = (int) ubuf(buf[m++]).i;   // # of atoms with values
    for (int iatom = 0; iatom < nonzero; iatom++) {
      offset = (int) ubuf(buf[m++]).i;      // offset into list for which atom
      j = list[offset];
      maxstrain_domain[j] = MAX(maxstrain_domain[j],buf[m]);
      m++;
    }

  // BIASFLAG
  // unpack biasflag vector

  } else if (commflag == BIASFLAG) {
    for (i = 0; i < n; i++) {
      j = list[i];
      biasflag[j] = (tagint) ubuf(buf[m++]).i;
    }

  // BIASCOEFF
  // unpack list of biascoeffs
  // nonzero # of entries, each has offset to which atom in receiver's list
  // protect against overflow of clist vector
  // if that happens, caller will re-setup cpage and reverse comm again

  } else if (commflag == BIASCOEFF) {
    int offset,ncoeff;
    int nonzero = (int) ubuf(buf[m++]).i;   // # of atoms with coeffs
    for (int iatom = 0; iatom < nonzero; iatom++) {
      offset = (int) ubuf(buf[m++]).i;      // offset into list for which atom
      j = list[offset];
      ncoeff = (int) ubuf(buf[m++]).i;      // # of bias coeffs
      for (k = 0; k < ncoeff; k++) {
        if (numcoeff[j] == 0) clist[j] = cpage->get(maxbondperatom);
        if (numcoeff[j] < maxbondperatom) {
          clist[j][numcoeff[j]].biascoeff = buf[m++];
          clist[j][numcoeff[j]].tag = (tagint) ubuf(buf[m++]).i;
        } else m += 2;
        numcoeff[j]++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   grow bond list and bias coeff vector by a chunk
------------------------------------------------------------------------- */

void FixHyperLocal::grow_bond()
{
  if (maxbond + DELTABOND > MAXSMALLINT)
    error->one(FLERR,"Fix hyper/local bond count is too big");
  maxbond += DELTABOND;
  blist = (OneBond *)
    memory->srealloc(blist,maxbond*sizeof(OneBond),"hyper/local:blist");
  memory->grow(biascoeff,maxbond,"hyper/local:biascoeff");
  vector_local = biascoeff;
}

/* ---------------------------------------------------------------------- */

double FixHyperLocal::compute_scalar()
{
  double allbias;
  MPI_Allreduce(&mybias,&allbias,1,MPI_DOUBLE,MPI_SUM,world);
  return allbias;
}

/* ---------------------------------------------------------------------- */

double FixHyperLocal::compute_vector(int i)
{
  // 21 vector outputs returned for i = 0-20

  // i = 0 = # of biased bonds on this step
  // i = 1 = max strain of any bond on this step
  // i = 2 = average bias coeff for all bonds on this step
  // i = 3 = ave bonds/atom on this step
  // i = 4 = ave neighbor bonds/bond on this step

  // i = 5 = max bond length during this run
  // i = 6 = average # of biased bonds/step during this run
  // i = 7 = fraction of biased bonds with no bias during this run
  // i = 8 = fraction of biased bonds with negative strain during this run
  // i = 9 = average bias coeff for all bonds during this run
  // i = 10 = min bias coeff for any bond during this run
  // i = 11 = max bias coeff for any bond during this run

  // i = 12 = max drift distance of any atom during this run
  // i = 13 = max distance from proc subbox of any ghost atom with
  //          maxstrain < qfactor during this run
  // i = 14 = max distance from proc subbox of any ghost atom with
  //          any maxstrain during this run
  // i = 15 = count of ghost atoms that could not be found
  //          on reneighbor steps during this run
  // i = 16 = count of bias overlaps (< Dcut) found during this run

  // i = 17 = cumulative hyper time since fix created
  // i = 18 = cumulative # of event timesteps since fix created
  // i = 19 = cumulative # of atoms in events since fix created
  // i = 20 = cumulative # of new bonds formed since fix created

  if (i == 0) {
    int nbiasall;
    MPI_Allreduce(&nbias,&nbiasall,1,MPI_INT,MPI_SUM,world);
    return (double) nbiasall;
  }

  if (i == 1) {
    if (nostrainyet) return 0.0;
    int nlocal = atom->nlocal;
    double emax = 0.0;
    for (int j = 0; j < nlocal; j++)
      emax = MAX(emax,maxstrain[j]);
    double eall;
    MPI_Allreduce(&emax,&eall,1,MPI_DOUBLE,MPI_MAX,world);
    return eall;
  }

  if (i == 2) {
    if (allbonds) return sumbiascoeff/allbonds;
    return 1.0;
  }

  if (i == 3) {
    return 2.0*allbonds/atom->natoms;
  }

  if (i == 4) {
    bigint allneigh,thisneigh;
    thisneigh = listfull->ipage->ndatum;
    MPI_Allreduce(&thisneigh,&allneigh,1,MPI_LMP_BIGINT,MPI_SUM,world);
    const double natoms = atom->natoms;
    const double neighsperatom = static_cast<double>(allneigh)/natoms;
    const double bondsperatom = static_cast<double>(allbonds)/natoms;
    return neighsperatom * bondsperatom;
  }

  if (i == 5) {
    double allbondlen;
    MPI_Allreduce(&maxbondlen,&allbondlen,1,MPI_DOUBLE,MPI_MAX,world);
    return allbondlen;
  }

  if (i == 6) {
    if (update->ntimestep == update->firststep) return 0.0;
    int allbias_running;
    MPI_Allreduce(&nbias_running,&allbias_running,1,MPI_INT,MPI_SUM,world);
    return 1.0*allbias_running / (update->ntimestep - update->firststep);
  }

  if (i == 7) {
    int allbias_running,allnobias_running;
    MPI_Allreduce(&nbias_running,&allbias_running,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&nobias_running,&allnobias_running,1,MPI_INT,MPI_SUM,world);
    if (allbias_running) return 1.0*allnobias_running / allbias_running;
    return 0.0;
  }

  if (i == 8) {
    int allbias_running,allnegstrain_running;
    MPI_Allreduce(&nbias_running,&allbias_running,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&negstrain_running,&allnegstrain_running,1,MPI_INT,
                  MPI_SUM,world);
    if (allbias_running) return 1.0*allnegstrain_running / allbias_running;
    return 0.0;
  }

  if (i == 9) {
    if (update->ntimestep == update->firststep) return 0.0;
    return avebiascoeff / (update->ntimestep - update->firststep);
  }

  if (i == 10) {
    double biascoeff;
    MPI_Allreduce(&minbiascoeff,&biascoeff,1,MPI_DOUBLE,MPI_MIN,world);
    return biascoeff;
  }

  if (i == 11) {
    double biascoeff;
    MPI_Allreduce(&maxbiascoeff,&biascoeff,1,MPI_DOUBLE,MPI_MAX,world);
    return biascoeff;
  }

  if (i == 12) {
    double alldriftsq;
    MPI_Allreduce(&maxdriftsq,&alldriftsq,1,MPI_DOUBLE,MPI_MAX,world);
    return (double) sqrt(alldriftsq);
  }

  if (i == 13) return rmaxever;
  if (i == 14) return rmaxeverbig;

  if (i == 15) {
    int allghost_toofar;
    MPI_Allreduce(&ghost_toofar,&allghost_toofar,1,MPI_INT,MPI_SUM,world);
    return 1.0*allghost_toofar;
  }

  if (i == 16) {
    int allclose;
    MPI_Allreduce(&checkbias_count,&allclose,1,MPI_INT,MPI_SUM,world);
    return 1.0*allclose;
  }

  if (i == 17) {
    return boost_target * update->dt * (update->ntimestep - starttime);
  }

  if (i == 18) return (double) nevent;
  if (i == 19) return (double) nevent_atom;

  if (i == 20) {
    int allnewbond;
    MPI_Allreduce(&nnewbond,&allnewbond,1,MPI_INT,MPI_SUM,world);
    return (double) allnewbond;
  }

  return 0.0;
}

/* ----------------------------------------------------------------------
   wrapper on compute_vector()
   used by hyper.cpp to call FixHyper
------------------------------------------------------------------------- */

double FixHyperLocal::query(int i)
{
  if (i == 1) return compute_vector(17);  // cummulative hyper time
  if (i == 2) return compute_vector(18);  // nevent
  if (i == 3) return compute_vector(19);  // nevent_atom
  if (i == 4) return compute_vector(3);   // ave bonds/atom
  if (i == 5) return compute_vector(12);  // maxdrift
  if (i == 6) return compute_vector(5);   // maxbondlen
  if (i == 7) return compute_vector(7);   // fraction with zero bias
  if (i == 8) return compute_vector(8);   // fraction with negative strain

  // unique to local hyper

  if (i == 9) return compute_vector(20);   // number of new bonds
  if (i == 10) return 1.0*maxbondperatom;  // max bonds/atom
  if (i == 11) return compute_vector(6);   // ave # of biased bonds/step
  if (i == 12) return compute_vector(9);   // ave bias coeff over all bonds
  if (i == 13) return compute_vector(10);  // min bias cooef for any bond
  if (i == 14) return compute_vector(11);  // max bias cooef for any bond
  if (i == 15) return compute_vector(4);   // neighbor bonds/bond
  if (i == 16) return compute_vector(2);   // ave bias coeff now
  if (i == 17) return time_bondbuild;      // CPU time for build_bond calls
  if (i == 18) return rmaxever;            // ghost atom distance for < maxstrain
  if (i == 19) return rmaxeverbig;         // ghost atom distance for any strain
  if (i == 20) return compute_vector(15);  // count of ghost atoms not found
  if (i == 21) return compute_vector(16);  // count of bias overlaps

  error->all(FLERR,"Invalid query to fix hyper/local");

  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of per-atom and per-bond data structs
------------------------------------------------------------------------- */

double FixHyperLocal::memory_usage()
{
  double bytes = maxbond * sizeof(OneBond);       // blist
  bytes = maxbond * sizeof(double);               // per-bond bias coeffs
  bytes += 3*maxlocal * sizeof(int);              // numbond,maxhalf,eligible
  bytes += maxlocal * sizeof(double);             // maxhalfstrain
  bytes += maxall * sizeof(int);                  // old2now
  bytes += maxall * sizeof(tagint);               // tagold
  bytes += 3*maxall * sizeof(double);             // xold
  bytes += 2*maxall * sizeof(double);             // maxstrain,maxstrain_domain
  if (checkbias) bytes += maxall * sizeof(tagint);  // biasflag
  bytes += maxcoeff * sizeof(int);                // numcoeff
  bytes += maxcoeff * sizeof(OneCoeff *);         // clist
  bytes += maxlocal*maxbondperatom * sizeof(OneCoeff);  // cpage estimate
  return bytes;
}
