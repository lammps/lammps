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

#define DELTABOOST 16
#define BOOSTINIT 1.0
#define COEFFMAX 1.2
#define BIG 1.0e20

enum{STRAIN,STRAINREGION,BIASFLAG};
enum{IGNORE,WARN,ERROR};

/* ---------------------------------------------------------------------- */

FixHyperLocal::FixHyperLocal(LAMMPS *lmp, int narg, char **arg) :
  FixHyper(lmp, narg, arg), old2now(NULL), xold(NULL), tagold(NULL),
  bonds(NULL), numbond(NULL), maxstrain(NULL), maxstrain_region(NULL),
  maxstrain_bondindex(NULL), biasflag(NULL), boost(NULL),
  histo(NULL), allhisto(NULL)
{
  // error checks
  // solution for tagint != int is to worry about storing
  //   local index vs global ID in same variable
  //   maybe need to declare them all tagint, not int

  if (atom->map_style == 0)
    error->all(FLERR,"Fix hyper/local command requires atom map");

  if (sizeof(tagint) != sizeof(int))
     error->all(FLERR,"Fix hyper/local requires tagint = int");

  // parse args

  if (narg < 10) error->all(FLERR,"Illegal fix hyper/local command");

  hyperflag = 2;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 23;

  global_freq = 1;
  extscalar = 0;
  extvector = 0;

  cutbond = force->numeric(FLERR,arg[3]);
  qfactor = force->numeric(FLERR,arg[4]);
  vmax = force->numeric(FLERR,arg[5]);
  tequil = force->numeric(FLERR,arg[6]);
  dcut = force->numeric(FLERR,arg[7]);
  alpha_user = force->numeric(FLERR,arg[8]);
  boosttarget = force->numeric(FLERR,arg[9]);

  if (cutbond < 0.0 || qfactor < 0.0 || vmax < 0.0 ||
      tequil <= 0.0 || dcut <= 0.0 || alpha_user <= 0.0 || boosttarget < 1.0)
    error->all(FLERR,"Illegal fix hyper/local command");

  invqfactorsq = 1.0 / (qfactor*qfactor);
  cutbondsq = cutbond*cutbond;
  dcutsq = dcut*dcut;
  beta = 1.0 / (force->boltz * tequil);

  // optional args

  histoflag = 0;
  lostbond = IGNORE;
  checkbias = 0;
  checkcoeff = 0;

  int iarg = 10;
  while (iarg < narg) {
    /*  NOTE: do not enable this yet, need to think about it differently
    if (strcmp(arg[iarg],"histo") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix hyper/local command");
      histoflag = 1;
      histo_every = force->inumeric(FLERR,arg[iarg+1]);
      histo_count = force->inumeric(FLERR,arg[iarg+2]);
      histo_delta = force->numeric(FLERR,arg[iarg+3]);
      histo_print = force->inumeric(FLERR,arg[iarg+4]);
      if (histo_every <= 0 || histo_count % 2 ||
          histo_delta <= 0.0 || histo_print <= 0)
        error->all(FLERR,"Illegal fix hyper/local command");
      iarg += 5;
    */

    if (strcmp(arg[iarg],"lostbond") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix hyper/local command");
      if (strcmp(arg[iarg+1],"error") == 0) lostbond = ERROR;
      else if (strcmp(arg[iarg+1],"warn") == 0) lostbond = WARN;
      else if (strcmp(arg[iarg+1],"ignore") == 0) lostbond = IGNORE;
      else error->all(FLERR,"Illegal fix hyper/local command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"check/bias") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix hyper/local command");
      checkbias = 1;
      checkbias_every = force->inumeric(FLERR,arg[iarg+1]);
      if (strcmp(arg[iarg+2],"error") == 0) checkbias_flag = ERROR;
      else if (strcmp(arg[iarg+2],"warn") == 0) checkbias_flag = WARN;
      else if (strcmp(arg[iarg+2],"ignore") == 0) checkbias_flag = IGNORE;
      else error->all(FLERR,"Illegal fix hyper/local command");
      iarg += 3;

    } else if (strcmp(arg[iarg],"check/coeff") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix hyper/local command");
      checkcoeff = 1;
      checkcoeff_every = force->inumeric(FLERR,arg[iarg+1]);
      if (strcmp(arg[iarg+2],"error") == 0) checkcoeff_flag = ERROR;
      else if (strcmp(arg[iarg+2],"warn") == 0) checkcoeff_flag = WARN;
      else if (strcmp(arg[iarg+2],"ignore") == 0) checkcoeff_flag = IGNORE;
      else error->all(FLERR,"Illegal fix hyper/local command");
      iarg += 3;

    } else error->all(FLERR,"Illegal fix hyper/local command");
  }

  // per-atom data structs

  maxbond = 0;
  bonds = NULL;
  numbond = NULL;
  maxstrain = NULL;
  maxstrain_region = NULL;
  maxstrain_bondindex = NULL;
  biasflag = NULL;

  maxold = old_nall = 0;
  old2now = NULL;
  xold = NULL;      // NOTE: don't really need except to monitor drift
  tagold = NULL;

  nboost = maxboost = 0;
  boost = NULL;

  // maxbondperatom = max # of bonds any atom is part of
  // will be reset in bond_build()
  // set comm size needed by this fix

  maxbondperatom = 1;
  comm_forward = 1;
  comm_reverse = 1;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  grow_arrays(atom->nmax);
  atom->add_callback(0);

  me = comm->me;
  firstflag = 1;

  allbonds = 0;
  allboost = 0.0;

  starttime = update->ntimestep;
  nostrainyet = 1;
  nnewbond = 0;
  nevent = 0;
  nevent_atom = 0;
  mybias = 0.0;

  histo = allhisto = NULL;
  if (histoflag) {
    invhisto_delta = 1.0 / histo_delta;
    histo_lo = 1.0 - (histo_count/2 * histo_delta);
    histo = new bigint[histo_count+2];
    allhisto = new bigint[histo_count+2];
  }
}

/* ---------------------------------------------------------------------- */

FixHyperLocal::~FixHyperLocal()
{
  memory->destroy(bonds);
  memory->destroy(numbond);

  atom->delete_callback(id,0);

  memory->destroy(maxstrain);
  memory->destroy(maxstrain_region);
  memory->destroy(maxstrain_bondindex);
  memory->destroy(biasflag);

  memory->destroy(old2now);
  memory->destroy(xold);
  memory->destroy(tagold);
  memory->destroy(boost);
  delete [] histo;
  delete [] allhisto;
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
  lostbond_partner = 0;
  lostbond_coeff = 0.0;
  checkbias_count = 0;
  checkcoeff_count = 0;
  maxdriftsq = 0.0;
  maxbondlen = 0.0;
  maxboostcoeff = 0.0;
  minboostcoeff = BIG;
  sumboostcoeff = 0.0;
  nboost_running = 0;
  nobias_running = 0;
  rmaxever = 0.0;
  rmaxeverbig = 0.0;

  nbondbuild = 0;
  time_bondbuild = 0.0;

  if (histoflag) histo_steps = 0;
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::init()
{
  // for newton off, bond force bias will not be applied correctly
  //   bonds that straddle 2 procs
  // warn if molecular system, since near-neighbors may not appear in neigh list
  //   user should not be including bonded atoms as hyper "bonds"

  if (force->newton_pair == 0)
    error->all(FLERR,"Hyper local requires newton pair on");

  if (atom->molecular && me == 0)
    error->warning(FLERR,"Hyper local for molecular systems "
                   "requires care in defining hyperdynamics bonds");

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
      error->all(FLERR,"Fix hyper/local bond cutoff exceeds ghost atom range - "
                 "use comm_modify cutoff command");
    if (cutghost < dcut+cutbond/2.0 && me == 0)
      error->warning(FLERR,"Fix hyper/local ghost atom range "
                     "may not allow for atom drift between events");
  }


  alpha = update->dt / alpha_user;

  // need an occasional full neighbor list with cutoff = Dcut
  // do not need to include neigh skin in cutoff,
  //   b/c this list will be built every time bond_build() is called
  // NOTE: what if pair style list cutoff > Dcut
  //   or what if neigh skin is huge?

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->cut = 1;
  neighbor->requests[irequest]->cutoff = dcut;
  neighbor->requests[irequest]->occasional = 1;

  // extra timing output

  //timefirst = timesecond = timethird = timefourth = timefifth =
  //  timesixth = timeseventh = timetotal = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::init_list(int /* id */, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::setup_pre_neighbor()
{
  // called for dynamics and minimization   NOTE: check if min is needed?

  pre_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::setup_pre_reverse(int eflag, int vflag)
{
  // only called for dynamics, not minimization
  // setupflag prevents boostostat update of boost coeffs in setup
  // also prevents increments of nboost_running, nbias_running, sumboostcoeff

  setupflag = 1;
  pre_reverse(eflag,vflag);
  setupflag = 0;
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::pre_neighbor()
{
  int i,m,n,ilocal,jlocal;

  // convert global ID bond partners back to local indices
  // need to use closest_image() so can calculate bond lengths
  // error flag should not happen for a well-behaved system
  // b/c are only looking up bond partners inside or near my sub-domain

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int missing = 0;
  double missing_coeff = 0.0;

  for (i = 0; i < nlocal; i++) {
    n = numbond[i];
    for (m = 0; m < n; m++) {
      jlocal = atom->map(bonds[i][m].jtag);
      if (jlocal >= 0) bonds[i][m].j = domain->closest_image(i,jlocal);
      else {
        bonds[i][m].j = -1;
        missing++;
        missing_coeff += bonds[i][m].boostcoeff;
        if (lostbond != IGNORE) {
          char str[128];
          sprintf(str,"Fix hyper/local bond info missing for bond "
                  TAGINT_FORMAT "," TAGINT_FORMAT
                  " with coeff %g at step " BIGINT_FORMAT,
                  atom->tag[i],bonds[i][m].jtag,bonds[i][m].boostcoeff,
                  update->ntimestep);
          if (lostbond == ERROR) error->one(FLERR,str);
          if (lostbond == WARN) error->warning(FLERR,str);
        }
      }
    }
  }

  lostbond_partner += missing;
  lostbond_coeff += missing_coeff;

  // set old2now to point to current local atom indices
  // only necessary for tagold entries > 0
  //   because if tagold = 0, atom is not active in Dcut neighbor list
  // must be done after atoms migrate and ghost atoms setup via comm->borders()
  // does not matter if there are multiple ghost copies of a global ID
  //   since only need the ghost atom strain, not its coordinates
  // NOTE: maybe need not use closest image, b/c old2now only used to access
  //   maxstrain values, which will be same for any copy of ghost atom ??
  //   also need it to compute maxdriftsq correctly when a proc spans a PBC

  double distsq;

  for (i = 0; i < old_nall; i++) {
    if (tagold[i] == 0) continue;
    ilocal = atom->map(tagold[i]);
    ilocal = domain->closest_image(xold[i],ilocal);
    old2now[i] = ilocal;

    if (ilocal >= 0) {
      distsq = MathExtra::distsq3(x[ilocal],xold[i]);
      maxdriftsq = MAX(distsq,maxdriftsq);
    } else ghost_toofar++;
  }
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::pre_reverse(int /* eflag */, int /* vflag */)
{
  int i,j,m,ii,jj,inum,jnum,iold,jold,nbond,bondindex;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double r,r0,estrain,emax,vbias,fbias,fbiasr,boostcoeff;
  int *ilist,*jlist,*numneigh,**firstneigh;

  //double time1,time2,time3,time4,time5,time6,time7,time8;
  //time1 = MPI_Wtime();

  // compute current maxstrain and maxstrain_bond for each owned atom
  // use per-atom full bond list
  // this is double-calculating for IJ and JI bonds
  //   could compute once, but would have to find/store index of JI bond
  // order two I,J atoms consistently for IJ and JI calcs
  //   to insure no round-off issue when comparing maxstrain values of I,J

  double **x = atom->x;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int mybonds = 0;
  nostrainyet = 0;

  for (i = 0; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itag = tag[i];
    emax = 0.0;
    bondindex = -1;
    jtag = 0;

    nbond = numbond[i];
    mybonds += nbond;
    for (m = 0; m < nbond; m++) {
      j = bonds[i][m].j;
      if (j < 0) continue;
      jtag = bonds[i][m].jtag;
      r0 = bonds[i][m].r0;
      if (itag < jtag) {
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
      } else {
        delx = x[j][0] - xtmp;
        dely = x[j][1] - ytmp;
        delz = x[j][2] - ztmp;;
      }
      r = sqrt(delx*delx + dely*dely + delz*delz);
      maxbondlen = MAX(r,maxbondlen);
      estrain = fabs(r-r0) / r0;
      if (estrain > emax) {
        emax = estrain;
        bondindex = m;
      }
    }
    maxstrain[i] = emax;
    maxstrain_bondindex[i] = bondindex;
  }

  //time2 = MPI_Wtime();

  // forward comm to acquire maxstrain of all ghost atoms

  commflag = STRAIN;
  comm->forward_comm_fix(this);

  //time3 = MPI_Wtime();

  // use original Dcut neighbor list to check maxstrain of all neighbor atoms
  // set maxstrain_region of I atoms = maxstrain of I and all J neighs
  // neighbor list has old indices for IJ b/c reneighboring may have occurred
  //   use old2now[] to convert to current indices
  // if neighbor is not currently known (too far away),
  //   then assume it was part of an event and its strain = qfactor
  // this double loop sets maxstrain_region of mostly owned atoms
  //   but possibly some ghost atoms as well

  int nall = nlocal + atom->nghost;
  for (i = 0; i < nall; i++) maxstrain_region[i] = 0.0;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // find largest distance from subbox that a ghost atom is with strain < qfactor

  double rmax = rmaxever;
  double rmaxbig = rmaxeverbig;
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  for (ii = 0; ii < inum; ii++) {
    iold = ilist[ii];
    jlist = firstneigh[iold];
    jnum = numneigh[iold];

    // I and J may be ghost atoms
    // will always know I b/c atoms do not drift that far
    // but may no longer know J if hops outside cutghost
    // in that case, assume it performed an event, its strain = qfactor

    i = old2now[iold];
    emax = maxstrain[i];

    for (jj = 0; jj < jnum; jj++) {
      jold = jlist[jj];
      j = old2now[jold];
      if (j >= 0) emax = MAX(emax,maxstrain[j]);
      else {
        emax = MAX(emax,qfactor);
        continue;
      }

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

    maxstrain_region[i] = emax;
  }

  double rmax2[2],rmax2all[2];
  rmax2[0] = rmax;
  rmax2[1] = rmaxbig;
  MPI_Allreduce(&rmax2,&rmax2all,2,MPI_DOUBLE,MPI_MAX,world);
  rmaxever = rmax2all[0];
  rmaxeverbig = rmax2all[1];

  MPI_Allreduce(&mybonds,&allbonds,1,MPI_INT,MPI_SUM,world);

  //time4 = MPI_Wtime();

  // reverse comm to acquire maxstrain_region from ghost atoms
  // needed b/c neighbor list referred to old owned atoms,
  //   so above loop may set maxstrain_region of ghost atoms
  // forward comm to acquire maxstrain_region of all ghost atoms

  commflag = STRAINREGION;
  comm->reverse_comm_fix(this);
  comm->forward_comm_fix(this);

  //time5 = MPI_Wtime();

  // identify biased bonds and add to boost list
  // for each max-strain bond IJ of atom I:
  //   bias this bond only if all these conditions hold:
  //     itag < jtag, so bond is only biased once
  //     maxstrain[i] = maxstrain_region[i]
  //     maxstrain[j] = maxstrain_region[j]
  // NOTE: also need to check that maxstrain[i] = maxstrain[j] ??
  //       don't think so, b/c maxstrain_region[i] includes maxstrain[i]

  nboost = 0;

  for (i = 0; i < nlocal; i++) {
    if (numbond[i] == 0) continue;
    itag = tag[i];
    j = bonds[i][maxstrain_bondindex[i]].j;
    jtag = tag[j];
    if (itag > jtag) continue;

    if (maxstrain[i] != maxstrain_region[i]) continue;
    if (maxstrain[j] != maxstrain_region[j]) continue;

    if (nboost == maxboost) {
      maxboost += DELTABOOST;
      memory->grow(boost,maxboost,"hyper/local:boost");
    }
    boost[nboost++] = i;
  }

  //time6 = MPI_Wtime();

  // apply boost force to bonds with locally max strain

  double **f = atom->f;

  int nobias = 0;
  mybias = 0.0;

  for (int iboost = 0; iboost < nboost; iboost++) {
    i = boost[iboost];
    emax = maxstrain[i];
    if (emax >= qfactor) {
      nobias++;
      continue;
    }

    m = maxstrain_bondindex[i];
    j = bonds[i][m].j;
    r0 = bonds[i][m].r0;
    boostcoeff = bonds[i][m].boostcoeff;

    vbias = boostcoeff * vmax * (1.0 - emax*emax*invqfactorsq);
    fbias = boostcoeff * 2.0 * vmax * emax / (qfactor*qfactor * r0);

    delx = x[i][0] - x[j][0];
    dely = x[i][1] - x[j][1];
    delz = x[i][2] - x[j][2];
    r = sqrt(delx*delx + dely*dely + delz*delz);
    fbiasr = fbias / r;

    f[i][0] += delx*fbiasr;
    f[i][1] += dely*fbiasr;
    f[i][2] += delz*fbiasr;

    f[j][0] -= delx*fbiasr;
    f[j][1] -= dely*fbiasr;
    f[j][2] -= delz*fbiasr;

    mybias += vbias;
  }

  //time7 = MPI_Wtime();

  // no boostostat update of boost coeffs when pre_reverse called from setup()
  // nboost_running, nobias_running, sumboostcoeff only incremented on run steps
  // NOTE: maybe should also not bias any bonds on firststep of this fix

  if (setupflag) return;
  nboost_running += nboost;
  nobias_running += nobias;

  // apply boostostat to boost coefficients of all bonds of all owned atoms
  // use per-atom full bond list
  // this is double-calculating for IJ and JI bonds
  //   should be identical for both, b/c emax is the same
  //   could compute once, but would have to find/store index of JI bond
  // delta in boost coeff is function of maxboost_region vs target boost
  // maxboost_region is function of two maxstrain_regions for I,J
  // NOTE: if J is lost to I but not vice versa, then biascoeff IJ != JI

  double myboost = 0.0;
  double emaxi,emaxj,maxboost_region;

  for (i = 0; i < nlocal; i++) {
    emaxi = maxstrain_region[i];
    nbond = numbond[i];
    for (m = 0; m < nbond; m++) {
      j = bonds[i][m].j;
      if (j < 0) continue;
      emaxj = maxstrain_region[j];
      emax = MAX(emaxi,emaxj);
      if (emax < qfactor) vbias = vmax * (1.0 - emax*emax*invqfactorsq);
      else vbias = 0.0;
      boostcoeff = bonds[i][m].boostcoeff;
      maxboost_region = exp(beta * boostcoeff*vbias);
      boostcoeff -= alpha * (maxboost_region-boosttarget) / boosttarget;
      // COMMENT OUT for now - need better way to bound boostcoeff
      //boostcoeff = MIN(boostcoeff,COEFFMAX);
      myboost += boostcoeff;
      maxboostcoeff = MAX(maxboostcoeff,boostcoeff);
      minboostcoeff = MIN(minboostcoeff,boostcoeff);
      bonds[i][m].boostcoeff = boostcoeff;
    }
  }

  // running stats

  MPI_Allreduce(&myboost,&allboost,1,MPI_DOUBLE,MPI_SUM,world);
  if (allbonds) sumboostcoeff += allboost/allbonds;

  // histogram the bond coeffs and output as requested
  // do not double count each bond

  if (histoflag && update->ntimestep % histo_every == 0) {
    if (histo_steps == 0)
      for (i = 0; i < histo_count+2; i++) histo[i] = 0;
    histo_steps++;

    int ihisto;
    for (i = 0; i < nlocal; i++) {
      nbond = numbond[i];
      for (m = 0; m < nbond; m++) {
        if (tag[i] > bonds[i][m].jtag) continue;
        boostcoeff = bonds[i][m].boostcoeff;
        if (boostcoeff < histo_lo) ihisto = -1;
        else ihisto = static_cast<int> ((boostcoeff-histo_lo) * invhisto_delta);
        if (ihisto >= histo_count) ihisto = histo_count;
        histo[ihisto+1]++;
      }
    }

    if (update->ntimestep % histo_print == 0) {
      MPI_Allreduce(histo,allhisto,histo_count+2,MPI_LMP_BIGINT,MPI_SUM,world);

      bigint total = 0;
      for (i = 0; i < histo_count+2; i++) total += allhisto[i];

      if (me == 0) {
        if (screen) {
          fprintf(screen,"Histogram of bias coeffs:\n");
          for (i = 0; i < histo_count+2; i++)
            fprintf(screen,"  %g",1.0*allhisto[i]/total);
          fprintf(screen,"\n");
        }
        if (logfile) {
          fprintf(logfile,"Histogram of bias coeffs:\n");
          for (i = 0; i < histo_count+2; i++)
            fprintf(logfile,"  %g",1.0*allhisto[i]/total);
          fprintf(logfile,"\n");
        }
      }
    }
  }

  // check for any biased bonds that are too close to each other
  // keep a running count for output

  if (checkbias && update->ntimestep % checkbias_every == 0) {

    // mark each atom in a biased bond with ID of partner
    // nboost loop will mark some ghost atoms

    for (i = 0; i < nall; i++) biasflag[i] = 0;

    for (int iboost = 0; iboost < nboost; iboost++) {
      i = boost[iboost];
      m = maxstrain_bondindex[i];
      j = bonds[i][m].j;
      biasflag[i] = tag[j];
      biasflag[j] = tag[i];
    }

    // reverse comm to acquire biasflag from ghost atoms
    // needed b/c above loop may set biasflag of ghost atoms
    // forward comm to acquire biasflag of all ghost atoms

    commflag = BIASFLAG;
    comm->reverse_comm_fix(this);
    comm->forward_comm_fix(this);

    // I and J may be ghost atoms
    // only continue if I is a biased atom
    // if J is unknonw (drifted ghost) just ignore
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
  }

  // check for any bond bias coeffcients that do not match
  // cannot check unless both atoms IJ are owned by this proc
  // keep a running count for output

  if (checkcoeff && update->ntimestep % checkcoeff_every == 0) {
    int jb,jbonds;

    for (i = 0; i < nlocal; i++) {
      nbond = numbond[i];
      for (m = 0; m < nbond; m++) {
        if (tag[i] > bonds[i][m].jtag) continue;
        j = bonds[i][m].j;
        if (j < 0) continue;
        if (j >= nlocal) continue;
        itag = tag[i];
        jbonds = numbond[j];
        for (jb = 0; jb < jbonds; jb++)
          if (bonds[j][jb].jtag == itag) break;
        if (jb == jbonds)
          error->one(FLERR,"Fix hyper/local could not find duplicate bond");
        if (bonds[i][m].boostcoeff != bonds[j][jb].boostcoeff)
          checkcoeff_count++;
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
  int i,j,ii,jj,m,n,inum,jnum,nbond;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,oldcoeff;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double time1,time2;
  time1 = MPI_Wtime();

  if (natom) {
    nevent++;
    nevent_atom += natom;
  }

  // trigger Dcut neighbor list build
  // NOTE: turn off special bonds in this Dcut neigh list?

  neighbor->build_one(list);

  // make copy of old bonds to preserve boostcoeffs for bonds that persist
  // allocate new numbond

  OneBond **old_bonds = bonds;
  int *old_numbond = numbond;

  int nmax = atom->nmax;
  memory->create(numbond,nmax,"hyper/local:numbond");

  // old_nall = value of nall at time bonds are built
  // reallocate new xold and tagold if necessary
  // initialize xold to current coords
  // initialize tagold to zero, so atoms not in neighbor list will remain zero

  old_nall = atom->nlocal + atom->nghost;

  if (old_nall > maxold) {
    memory->destroy(xold);
    memory->destroy(tagold);
    memory->destroy(old2now);
    maxold = atom->nmax;
    memory->create(xold,maxold,3,"hyper/local:xold");
    memory->create(tagold,maxold,"hyper/local:tagold");
    memory->create(old2now,maxold,"hyper/local:old2now");
  }

  double **x = atom->x;

  memcpy(&xold[0][0],&x[0][0],3*old_nall*sizeof(double));
  for (i = 0; i < old_nall; i++) tagold[i] = 0;

  // create and populate new bonds data struct
  // while loop allows maxbondperatom to increase once if necessary
  //   don't know new maxbondperatom value until end of loop
  // in practice maxbondperatom will hardly ever increase
  //   since there is a physical max value

  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  while (1) {
    bonds = (OneBond **) memory->create(bonds,nmax,maxbondperatom,
                                        "hyper/local:bonds");
    if (bonds) memset(bonds[0],0,nmax*sizeof(OneBond));
    for (i = 0; i < nlocal; i++) numbond[i] = 0;

    // identify bonds assigned to each owned atom
    // do not create a bond between two non-group atoms
    // set tagold = global ID for all I,J atoms used in neighbor list
    //   tagold remains 0 for unused atoms, skipped in pre_neighbor

    int nbondmax = 0;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itag = tag[i];
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
        // order IJ to insure IJ and JI bonds are stored consistently

        if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

        if (itag < jtag) {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
        } else {
          delx = x[j][0] - xtmp;
          dely = x[j][1] - ytmp;
          delz = x[j][2] - ztmp;
        }

        rsq = delx*delx + dely*dely + delz*delz;

        // NOTE: could create two bonds for IJ both owned from one calc?
        //       have to skip one of 2 bonds in that case

        if (rsq < cutbondsq) {
          if (nbond >= maxbondperatom) {
            nbond++;
            continue;
          }

          bonds[i][nbond].r0 = sqrt(rsq);
          bonds[i][nbond].jtag = tag[j];
          bonds[i][nbond].j = j;

          if (firstflag) oldcoeff = 0.0;
          else {
            oldcoeff = 0.0;
            jtag = tag[j];
            n = old_numbond[i];
            for (m = 0; m < n; m++) {
              if (old_bonds[i][m].jtag == jtag) {
                oldcoeff = old_bonds[i][m].boostcoeff;
                break;
              }
            }
          }

          if (oldcoeff > 0.0) bonds[i][nbond].boostcoeff = oldcoeff;
          else {
            bonds[i][nbond].boostcoeff = BOOSTINIT;
            nnewbond++;
          }
          nbond++;
        }
      }
      numbond[i] = nbond;
      nbondmax = MAX(nbondmax,nbond);
    }

    // maxbondperatom must increase uniformly on all procs
    // since bonds are comunicated when atoms migrate

    int allnbondmax;
    MPI_Allreduce(&nbondmax,&allnbondmax,1,MPI_INT,MPI_MAX,world);
    if (allnbondmax <= maxbondperatom) break;

    maxbondperatom = allnbondmax;
    memory->destroy(bonds);
  }

  // deallocate old_bonds and old_numbond

  memory->destroy(old_bonds);
  memory->destroy(old_numbond);

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

  if (commflag == STRAIN) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = maxstrain[j];
    }

  // STRAINREGION
  // pack maxstrain_region vector

  } else if (commflag == STRAINREGION) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = maxstrain_region[j];
    }

  // BIASFLAG
  // pack biasflag vector

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
  // unpack maxstrain_region vector

  } else if (commflag == STRAINREGION) {
    for (i = first; i < last; i++) {
      maxstrain_region[i] = buf[m++];
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
  int i,m,last;

  m = 0;
  last = first + n;

  // STRAINREGION
  // pack maxstrain_region vector

  if (commflag == STRAINREGION) {
    for (i = first; i < last; i++) {
      buf[m++] = maxstrain_region[i];
    }

  // BIASFLAG
  // pack biasflag vector

  } else if (commflag == BIASFLAG) {
    for (i = first; i < last; i++) {
      buf[m++] = ubuf(biasflag[i]).d;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixHyperLocal::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;

  // STRAINREGION
  // unpack maxstrain_region vector

  if (commflag == STRAINREGION) {
    for (i = 0; i < n; i++) {
      j = list[i];
      maxstrain_region[j] += buf[m++];
    }

  // BIASFLAG
  // unpack biasflag vector

  } else if (commflag == BIASFLAG) {
    for (i = 0; i < n; i++) {
      j = list[i];
      biasflag[j] = (tagint) ubuf(buf[m++]).i;
    }
  }
}

/* ----------------------------------------------------------------------
   allocate atom-based arrays
------------------------------------------------------------------------- */

void FixHyperLocal::grow_arrays(int nmax)
{
  // NOTE: not all of these need to be Nmax in length, could allocate elsewhere

  memory->grow(maxstrain,nmax,"hyper/local:maxstrain");
  memory->grow(maxstrain_bondindex,nmax,"hyper/local:maxstrain_bondindex");
  memory->grow(maxstrain_region,nmax,"hyper/local:maxstrain_region");
  if (checkbias) memory->grow(biasflag,nmax,"hyper/local:biasflag");

  memory->grow(numbond,nmax,"hyper/local:numbond");
  memory->grow(bonds,nmax,maxbondperatom,"hyper/local:bonds");

  // zero so valgrind does not complain about memcpy() in copy()
  // also so loops in pre_neighbor() are OK before
  //   bonds are setup for the first time

  if (bonds) {
    memset(bonds[maxbond],0,(nmax-maxbond)*maxbondperatom*sizeof(OneBond));
    memset(&numbond[maxbond],0,(nmax-maxbond)*sizeof(int));
    maxbond = nmax;
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixHyperLocal::copy_arrays(int i, int j, int /* delflag */)
{
  // avoid valgrind copy-to-self warning

  if (i != j) memcpy(bonds[j],bonds[i],numbond[i]*sizeof(OneBond));
  numbond[j] = numbond[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixHyperLocal::pack_exchange(int i, double *buf)
{
  int m = 1;
  int n = numbond[i];
  buf[m++] = ubuf(n).d;
  for (int j = 0; j < n; j++) {
    buf[m++] = bonds[i][j].r0;
    buf[m++] = bonds[i][j].boostcoeff;
    buf[m++] = ubuf(bonds[i][j].jtag).d;
    buf[m++] = ubuf(bonds[i][j].j).d;
  }

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixHyperLocal::unpack_exchange(int nlocal, double *buf)
{
  int m = 1;
  int n = numbond[nlocal] = (int) ubuf(buf[m++]).i;
  for (int j = 0; j < n; j++) {
    bonds[nlocal][j].r0 = buf[m++];
    bonds[nlocal][j].boostcoeff = buf[m++];
    bonds[nlocal][j].jtag = (tagint) ubuf(buf[m++]).i;
    bonds[nlocal][j].j = (int) ubuf(buf[m++]).i;
  }

  return m;
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
  // 23 vector outputs returned for i = 0-22

  // i = 0 = # of boosted bonds on this step
  // i = 1 = max strain of any bond on this step
  // i = 2 = average bias potential for all bonds on this step
  // i = 3 = ave bonds/atom on this step
  // i = 4 = ave neighbor bonds/bond on this step

  // i = 5 = fraction of steps and bonds with no bias during this run
  // i = 6 = max drift distance of any atom during this run
  // i = 7 = max bond length during this run
  // i = 8 = average # of boosted bonds/step during this run
  // i = 9 = average bias potential for all bonds during this run
  // i = 10 = max bias potential for any bond during this run
  // i = 11 = min bias potential for any bond during this run
  // i = 12 = max dist from my box of any ghost atom with
  //          maxstain < qfactor during this run
  // i = 13 = max dist from my box of any ghost atom with
  //          any maxstrain during this run
  // i = 14 = count of ghost atoms that could not be found
  //          by any proc at any reneighbor step during this run
  // i = 15 = count of lost bond partners during this run
  // i = 16 = average bias coeff for lost bond partners during this run
  // i = 17 = count of bias overlaps found during this run
  // i = 18 = count of non-matching bias coefficients found during this run

  // i = 19 = cummulative hyper time
  // i = 20 = cummulative # of event timesteps since fix created
  // i = 21 = cummulative # of atoms in events since fix created
  // i = 22 = cummulative # of new bonds formed since fix created

  if (i == 0) {
    int nboostall;
    MPI_Allreduce(&nboost,&nboostall,1,MPI_INT,MPI_SUM,world);
    return (double) nboostall;
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
    if (allboost && allbonds) return allboost/allbonds * vmax;
    return 1.0;
  }

  if (i == 3) return 1.0*allbonds/atom->natoms;

  if (i == 4) {
    const int nlocal = atom->nlocal;
    bigint nbonds = 0;
    for (int j = 0; j < nlocal; j++)
      nbonds += numbond[j];
    bigint allbonds;
    MPI_Allreduce(&nbonds,&allbonds,1,MPI_LMP_BIGINT,MPI_SUM,world);
    bigint allneigh,thisneigh;
    thisneigh = list->ipage->ndatum;
    MPI_Allreduce(&thisneigh,&allneigh,1,MPI_LMP_BIGINT,MPI_SUM,world);
    const double natoms = atom->natoms;
    const double neighsperatom = static_cast<double>(allneigh)/natoms;
    const double bondsperatom = 0.5*static_cast<double>(allbonds)/natoms;
    return neighsperatom * bondsperatom;
  }

  if (i == 5) {
    int allboost_running,allnobias_running;
    MPI_Allreduce(&nboost_running,&allboost_running,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&nobias_running,&allnobias_running,1,MPI_INT,MPI_SUM,world);
    if (allboost_running) return 1.0*allnobias_running / allboost_running;
    return 0.0;
  }

  if (i == 6) {
    double alldriftsq;
    MPI_Allreduce(&maxdriftsq,&alldriftsq,1,MPI_DOUBLE,MPI_MAX,world);
    return (double) sqrt(alldriftsq);
  }

  if (i == 7) {
    double allbondlen;
    MPI_Allreduce(&maxbondlen,&allbondlen,1,MPI_DOUBLE,MPI_MAX,world);
    return allbondlen;
  }

  if (i == 8) {
    if (update->ntimestep == update->firststep) return 0.0;
    int allboost_running;
    MPI_Allreduce(&nboost_running,&allboost_running,1,MPI_INT,MPI_SUM,world);
    return 1.0*allboost_running / (update->ntimestep - update->firststep);
  }

  if (i == 9) {
    if (update->ntimestep == update->firststep) return 0.0;
    return sumboostcoeff * vmax / (update->ntimestep - update->firststep);
  }

  if (i == 10) {
    double allboostcoeff;
    MPI_Allreduce(&maxboostcoeff,&allboostcoeff,1,MPI_DOUBLE,MPI_MAX,world);
    return allboostcoeff * vmax;
  }

  if (i == 11) {
    double allboostcoeff;
    MPI_Allreduce(&minboostcoeff,&allboostcoeff,1,MPI_DOUBLE,MPI_MAX,world);
    return allboostcoeff * vmax;
  }

  if (i == 12) return rmaxever;
  if (i == 13) return rmaxeverbig;

  if (i == 14) {
    int allghost_toofar;
    MPI_Allreduce(&ghost_toofar,&allghost_toofar,1,MPI_INT,MPI_SUM,world);
    return 1.0*allghost_toofar;
  }

  if (i == 15) {
    int alllost;
    MPI_Allreduce(&lostbond_partner,&alllost,1,MPI_INT,MPI_SUM,world);
    return 1.0*alllost;
  }

  if (i == 16) {
    int alllost;
    MPI_Allreduce(&lostbond_partner,&alllost,1,MPI_INT,MPI_SUM,world);
    double allcoeff;
    MPI_Allreduce(&lostbond_coeff,&allcoeff,1,MPI_DOUBLE,MPI_SUM,world);
    if (alllost == 0) return 0;
    return allcoeff/alllost;
  }

  if (i == 17) {
    int allclose;
    MPI_Allreduce(&checkbias_count,&allclose,1,MPI_INT,MPI_SUM,world);
    return 1.0*allclose;
  }

  if (i == 18) {
    int allcoeff;
    MPI_Allreduce(&checkcoeff_count,&allcoeff,1,MPI_INT,MPI_SUM,world);
    return 1.0*allcoeff;
  }

  if (i == 19) {
    return boosttarget * update->dt * (update->ntimestep - starttime);
  }

  if (i == 20) return (double) nevent;
  if (i == 21) return (double) nevent_atom;

  if (i == 22) {
    int allnew;
    MPI_Allreduce(&nnewbond,&allnew,1,MPI_INT,MPI_SUM,world);
    return (double) 0.5*allnew;
  }

  return 0.0;
}

/* ----------------------------------------------------------------------
   wrapper on compute_vector()
   used by hyper.cpp to call FixHyper
------------------------------------------------------------------------- */

double FixHyperLocal::query(int i)
{
  if (i == 1) return compute_vector(19);  // cummulative hyper time
  if (i == 2) return compute_vector(20);  // nevent
  if (i == 3) return compute_vector(21);  // nevent_atom
  if (i == 4) return compute_vector(3);   // ave bonds/atom
  if (i == 5) return compute_vector(6);   // maxdrift
  if (i == 6) return compute_vector(7);   // maxbondlen
  if (i == 7) return compute_vector(5);   // fraction with zero bias

  // unique to local hyper

  if (i == 8) return compute_vector(22);   // number of new bonds
  if (i == 9) return 1.0*maxbondperatom;   // max bonds/atom
  if (i == 10) return compute_vector(8);   // ave # of boosted bonds/step
  if (i == 11) return compute_vector(9);   // ave boost coeff over all bonds
  if (i == 12) return compute_vector(10);  // max boost cooef for any bond
  if (i == 13) return compute_vector(11);  // max boost cooef for any bond
  if (i == 14) return compute_vector(4);   // neighbor bonds/bond
  if (i == 15) return compute_vector(2);   // ave boost cooef now
  if (i == 16) return time_bondbuild;      // CPU time for bond_build calls
  if (i == 17) return rmaxever;            // ghost atom distance for < maxstrain
  if (i == 18) return rmaxeverbig;         // ghost atom distance for any strain
  if (i == 19) return compute_vector(14);  // count of ghost atoms not found
  if (i == 20) return compute_vector(15);  // count of lost bond partners
  if (i == 21) return compute_vector(16);  // ave bias coeff of long bonds
  if (i == 22) return compute_vector(17);  // count of bias overlaps
  if (i == 23) return compute_vector(18);  // count of non-matching bias coeffs

  error->all(FLERR,"Invalid query to fix hyper/local");

  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of per-atom and per-bond arrays
------------------------------------------------------------------------- */

double FixHyperLocal::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 2*nmax * sizeof(int);     // numbond, maxstrain_bondindex
  bytes = 2*nmax * sizeof(double);         // maxstrain, maxstrain_region
  bytes += maxbondperatom*nmax * sizeof(OneBond);      // bonds
  bytes += 3*maxold * sizeof(double);      // xold
  bytes += maxold * sizeof(tagint);        // tagold
  bytes += maxold * sizeof(int);           // old2now
  return bytes;
}
