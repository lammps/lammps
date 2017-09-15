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
#include <string.h>
#include <stdio.h>
#include "fix_shear_history.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{DEFAULT,NPARTNER,PERPARTNER};

/* ---------------------------------------------------------------------- */

FixShearHistory::FixShearHistory(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  npartner(NULL), partner(NULL), shearpartner(NULL), pair(NULL), 
  ipage(NULL), dpage(NULL)
{
  if (narg != 4) error->all(FLERR,"Illegal fix SHEAR_HISTORY command");

  restart_peratom = 1;
  create_attribute = 1;

  newton_pair = force->newton_pair;

  dnum = force->inumeric(FLERR,arg[3]);
  dnumbytes = dnum * sizeof(double);

  onesided = 0;
  if (strcmp(id,"LINE_SHEAR_HISTORY") == 0) onesided = 1;
  if (strcmp(id,"TRI_SHEAR_HISTORY") == 0) onesided = 1;

  if (newton_pair) comm_reverse = 1;   // just for single npartner value
                                       // variable-size history communicated via
                                       // reverse_comm_fix_variable()

  // perform initial allocation of atom-based arrays
  // register with atom class

  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  pgsize = oneatom = 0;

  // initialize npartner to 0 so neighbor list creation is OK the 1st time

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) npartner[i] = 0;
  maxtouch = 0;

  nlocal_neigh = nall_neigh = 0;
  commflag = DEFAULT;
}

/* ---------------------------------------------------------------------- */

FixShearHistory::~FixShearHistory()
{
  // unregister this fix so atom class doesn't invoke it any more

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  memory->destroy(npartner);
  memory->sfree(partner);
  memory->sfree(shearpartner);

  // to better detect use-after-delete errors

  pair = NULL;
  npartner = NULL;
  partner = NULL;
  shearpartner = NULL;

  delete [] ipage;
  delete [] dpage;
}

/* ---------------------------------------------------------------------- */

int FixShearHistory::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= MIN_PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixShearHistory::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Granular shear history requires atoms have IDs");

  allocate_pages();
}

/* ----------------------------------------------------------------------
   create pages if first time or if neighbor pgsize/oneatom has changed
   note that latter could cause shear history info to be discarded
------------------------------------------------------------------------- */

void FixShearHistory::allocate_pages()
{
  int create = 0;
  if (ipage == NULL) create = 1;
  if (pgsize != neighbor->pgsize) create = 1;
  if (oneatom != neighbor->oneatom) create = 1;

  if (create) {
    delete [] ipage;
    delete [] dpage;

    pgsize = neighbor->pgsize;
    oneatom = neighbor->oneatom;
    int nmypage = comm->nthreads;
    ipage = new MyPage<tagint>[nmypage];
    dpage = new MyPage<double>[nmypage];
    for (int i = 0; i < nmypage; i++) {
      ipage[i].init(oneatom,pgsize);
      dpage[i].init(dnum*oneatom,dnum*pgsize);
    }
  }
}

/* ----------------------------------------------------------------------
   copy shear partner info from neighbor lists to atom arrays
   should be called whenever neighbor list stores current history info
     and need to store the info with owned atoms
   e.g. so atoms can migrate to new procs or between runs
     when atoms may be added or deleted (neighbor list becomes out-of-date)
   the next granular neigh list build will put this info back into neigh list
   called during run before atom exchanges, including for restart files
   called at end of run via post_run()
   do not call during setup of run (setup_pre_exchange)
     b/c there is no guarantee of a current neigh list (even on continued run)
   if run command does a 2nd run with pre = no, then no neigh list
     will be built, but old neigh list will still have the info
   onesided and newton on and newton off versions
------------------------------------------------------------------------- */

void FixShearHistory::pre_exchange()
{
  if (onesided) pre_exchange_onesided();
  else if (newton_pair) pre_exchange_newton();
  else pre_exchange_no_newton();
}

/* ----------------------------------------------------------------------
   onesided version for sphere contact with line/tri particles
   neighbor list has I = sphere, J = line/tri
   only store history info with spheres
------------------------------------------------------------------------- */

void FixShearHistory::pre_exchange_onesided()
{
  int i,j,ii,jj,m,n,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  // NOTE: all operations until very end are with nlocal_neigh <= current nlocal
  // b/c previous neigh list was built with nlocal_neigh
  // nlocal can be larger if other fixes added atoms at this pre_exchange()

  // zero npartner for owned atoms
  // clear 2 page data structures

  for (i = 0; i < nlocal_neigh; i++) npartner[i] = 0;

  ipage->reset();
  dpage->reset();

  // 1st loop over neighbor list, I = sphere, J = tri
  // only calculate npartner for owned spheres

  tagint *tag = atom->tag;
  NeighList *list = pair->list;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = list->listhistory->firstneigh;
  firstshear = list->listhistory->firstdouble;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    touch = firsttouch[i];

    for (jj = 0; jj < jnum; jj++)
      if (touch[jj]) npartner[i]++;
  }

  // get page chunks to store atom IDs and shear history for owned atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    n = npartner[i];
    partner[i] = ipage->get(n);
    shearpartner[i] = dpage->get(dnum*n);
    if (partner[i] == NULL || shearpartner[i] == NULL)
      error->one(FLERR,"Shear history overflow, boost neigh_modify one");
  }

  // 2nd loop over neighbor list, I = sphere, J = tri
  // store atom IDs and shear history for owned spheres
  // re-zero npartner to use as counter

  for (i = 0; i < nlocal_neigh; i++) npartner[i] = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    allshear = firstshear[i];
    jnum = numneigh[i];
    touch = firsttouch[i];

    for (jj = 0; jj < jnum; jj++) {
      if (touch[jj]) {
        shear = &allshear[dnum*jj];
        j = jlist[jj];
        j &= NEIGHMASK;
        m = npartner[i]++;
        partner[i][m] = tag[j];
        memcpy(&shearpartner[i][dnum*m],shear,dnumbytes);
      }
    }
  }

  // set maxtouch = max # of partners of any owned atom
  // bump up comm->maxexchange_fix if necessary
  
  maxtouch = 0;
  for (i = 0; i < nlocal_neigh; i++) maxtouch = MAX(maxtouch,npartner[i]);
  comm->maxexchange_fix = MAX(comm->maxexchange_fix,(dnum+1)*maxtouch+1);

  // zero npartner values from previous nlocal_neigh to current nlocal

  int nlocal = atom->nlocal;
  for (i = nlocal_neigh; i < nlocal; i++) npartner[i] = 0;
}

/* ----------------------------------------------------------------------
   newton on version, for sphere/sphere contacts
   performs reverse comm to acquire shear partner info from ghost atoms
------------------------------------------------------------------------- */

void FixShearHistory::pre_exchange_newton()
{
  int i,j,ii,jj,m,n,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*shearj,*allshear,**firstshear;

  // NOTE: all operations until very end are with 
  //   nlocal_neigh  <= current nlocal and nall_neigh
  // b/c previous neigh list was built with nlocal_neigh,nghost_neigh
  // nlocal can be larger if other fixes added atoms at this pre_exchange()

  // zero npartner for owned+ghost atoms
  // clear 2 page data structures

  for (i = 0; i < nall_neigh; i++) npartner[i] = 0;

  ipage->reset();
  dpage->reset();

  // 1st loop over neighbor list
  // calculate npartner for owned+ghost atoms

  tagint *tag = atom->tag;
  NeighList *list = pair->list;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = list->listhistory->firstneigh;
  firstshear = list->listhistory->firstdouble;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    touch = firsttouch[i];

    for (jj = 0; jj < jnum; jj++) {
      if (touch[jj]) {
        npartner[i]++;
        j = jlist[jj];
        j &= NEIGHMASK;
        npartner[j]++;
      }
    }
  }

  // perform reverse comm to augment owned npartner counts with ghost counts

  commflag = NPARTNER;
  comm->reverse_comm_fix(this,0);

  // get page chunks to store atom IDs and shear history for owned+ghost atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    n = npartner[i];
    partner[i] = ipage->get(n);
    shearpartner[i] = dpage->get(dnum*n);
    if (partner[i] == NULL || shearpartner[i] == NULL) {
      error->one(FLERR,"Shear history overflow, boost neigh_modify one");
    }
  }

  for (i = nlocal_neigh; i < nall_neigh; i++) {
    n = npartner[i];
    partner[i] = ipage->get(n);
    shearpartner[i] = dpage->get(dnum*n);
    if (partner[i] == NULL || shearpartner[i] == NULL) {
      error->one(FLERR,"Shear history overflow, boost neigh_modify one");
    }
  }

  // 2nd loop over neighbor list
  // store atom IDs and shear history for owned+ghost atoms
  // re-zero npartner to use as counter

  for (i = 0; i < nall_neigh; i++) npartner[i] = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    allshear = firstshear[i];
    jnum = numneigh[i];
    touch = firsttouch[i];

    for (jj = 0; jj < jnum; jj++) {
      if (touch[jj]) {
        shear = &allshear[dnum*jj];
        j = jlist[jj];
        j &= NEIGHMASK;
        m = npartner[i]++;
        partner[i][m] = tag[j];
        memcpy(&shearpartner[i][dnum*m],shear,dnumbytes);
        m = npartner[j]++;
        partner[j][m] = tag[i];
        shearj = &shearpartner[j][dnum*m];
        for (n = 0; n < dnum; n++) shearj[n] = -shear[n];
      }
    }
  }

  // perform reverse comm to augment
  // owned atom partner/shearpartner with ghost info
  // use variable variant b/c size of packed data can be arbitrarily large
  //   if many touching neighbors for large particle

  commflag = PERPARTNER;
  comm->reverse_comm_fix_variable(this);

  // set maxtouch = max # of partners of any owned atom
  // bump up comm->maxexchange_fix if necessary

  maxtouch = 0;
  for (i = 0; i < nlocal_neigh; i++) maxtouch = MAX(maxtouch,npartner[i]);
  comm->maxexchange_fix = MAX(comm->maxexchange_fix,(dnum+1)*maxtouch+1);

  // zero npartner values from previous nlocal_neigh to current nlocal

  int nlocal = atom->nlocal;
  for (i = nlocal_neigh; i < nlocal; i++) npartner[i] = 0;
}

/* ----------------------------------------------------------------------
   newton off version, for sphere/sphere contacts
   newton OFF works with smaller vectors that don't include ghost info
------------------------------------------------------------------------- */

void FixShearHistory::pre_exchange_no_newton()
{
  int i,j,ii,jj,m,n,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*shearj,*allshear,**firstshear;

  // NOTE: all operations until very end are with nlocal_neigh <= current nlocal
  // b/c previous neigh list was built with nlocal_neigh
  // nlocal can be larger if other fixes added atoms at this pre_exchange()

  // zero npartner for owned atoms
  // clear 2 page data structures

  for (i = 0; i < nlocal_neigh; i++) npartner[i] = 0;

  ipage->reset();
  dpage->reset();

  // 1st loop over neighbor list
  // calculate npartner for owned atoms

  tagint *tag = atom->tag;
  NeighList *list = pair->list;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = list->listhistory->firstneigh;
  firstshear = list->listhistory->firstdouble;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    touch = firsttouch[i];

    for (jj = 0; jj < jnum; jj++) {
      if (touch[jj]) {
        npartner[i]++;
        j = jlist[jj];
        j &= NEIGHMASK;
        if (j < nlocal_neigh) npartner[j]++;
      }
    }
  }

  // get page chunks to store atom IDs and shear history for owned atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    n = npartner[i];
    partner[i] = ipage->get(n);
    shearpartner[i] = dpage->get(dnum*n);
    if (partner[i] == NULL || shearpartner[i] == NULL)
      error->one(FLERR,"Shear history overflow, boost neigh_modify one");
  }

  // 2nd loop over neighbor list
  // store atom IDs and shear history for owned atoms
  // re-zero npartner to use as counter

  for (i = 0; i < nlocal_neigh; i++) npartner[i] = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    allshear = firstshear[i];
    jnum = numneigh[i];
    touch = firsttouch[i];

    for (jj = 0; jj < jnum; jj++) {
      if (touch[jj]) {
        shear = &allshear[dnum*jj];
        j = jlist[jj];
        j &= NEIGHMASK;
        m = npartner[i]++;
        partner[i][m] = tag[j];
        memcpy(&shearpartner[i][dnum*m],shear,dnumbytes);
        if (j < nlocal_neigh) {
          m = npartner[j]++;
          partner[j][m] = tag[i];
          shearj = &shearpartner[j][dnum*m];
          for (n = 0; n < dnum; n++) shearj[n] = -shear[n];
        }
      }
    }
  }

  // set maxtouch = max # of partners of any owned atom
  // bump up comm->maxexchange_fix if necessary
  
  maxtouch = 0;
  for (i = 0; i < nlocal_neigh; i++) maxtouch = MAX(maxtouch,npartner[i]);
  comm->maxexchange_fix = MAX(comm->maxexchange_fix,(dnum+1)*maxtouch+1);

  // zero npartner values from previous nlocal_neigh to current nlocal

  int nlocal = atom->nlocal;
  for (i = nlocal_neigh; i < nlocal; i++) npartner[i] = 0;
}

/* ---------------------------------------------------------------------- */

void FixShearHistory::min_pre_exchange()
{
  pre_exchange();
}

/* ---------------------------------------------------------------------- */

void FixShearHistory::post_run()
{
  pre_exchange();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixShearHistory::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax * sizeof(tagint *);
  bytes += nmax * sizeof(double *);

  int nmypage = comm->nthreads;
  for (int i = 0; i < nmypage; i++) {
    bytes += ipage[i].size();
    bytes += dpage[i].size();
  }

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixShearHistory::grow_arrays(int nmax)
{
  memory->grow(npartner,nmax,"shear_history:npartner");
  partner = (tagint **) memory->srealloc(partner,nmax*sizeof(tagint *),
                                         "shear_history:partner");
  shearpartner = (double **) memory->srealloc(shearpartner,
                                              nmax*sizeof(double *),
                                              "shear_history:shearpartner");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixShearHistory::copy_arrays(int i, int j, int delflag)
{
  // just copy pointers for partner and shearpartner
  // b/c can't overwrite chunk allocation inside ipage,dpage
  // incoming atoms in unpack_exchange just grab new chunks
  // so are orphaning chunks for migrating atoms
  // OK, b/c will reset ipage,dpage on next reneighboring

  npartner[j] = npartner[i];
  partner[j] = partner[i];
  shearpartner[j] = shearpartner[i];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixShearHistory::set_arrays(int i)
{
  npartner[i] = 0;
}

/* ----------------------------------------------------------------------
   only called by Comm::reverse_comm_fix_variable for PERPARTNER mode
------------------------------------------------------------------------- */

int FixShearHistory::pack_reverse_comm_size(int n, int first)
{
  int i,last;

  int m = 0;
  last = first + n;

  for (i = first; i < last; i++)
    m += 1 + (dnum+1)*npartner[i];

  return m;
}

/* ----------------------------------------------------------------------
   two modes: NPARTNER and PERPARTNER
------------------------------------------------------------------------- */

int FixShearHistory::pack_reverse_comm(int n, int first, double *buf)
{
  int i,k,last;

  int m = 0;
  last = first + n;

  if (commflag == NPARTNER) {
    for (i = first; i < last; i++) {
      buf[m++] = npartner[i];
    }
  } else if (commflag == PERPARTNER) {
    for (i = first; i < last; i++) {
      buf[m++] = npartner[i];
      for (k = 0; k < npartner[i]; k++) {
        buf[m++] = partner[i][k];
        memcpy(&buf[m],&shearpartner[i][dnum*k],dnumbytes);
        m += dnum;
      }
    }
  } else error->all(FLERR,"Unsupported comm mode in shear history");

  return m;
}

/* ----------------------------------------------------------------------
   two modes: NPARTNER and PERPARTNER
------------------------------------------------------------------------- */

void FixShearHistory::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,k,kk,ncount;

  int m = 0;

  if (commflag == NPARTNER) {
    for (i = 0; i < n; i++) {
      j = list[i];
      npartner[j] += static_cast<int> (buf[m++]);
    }
  } else if (commflag == PERPARTNER) {
    for (i = 0; i < n; i++) {
      j = list[i];
      ncount = static_cast<int> (buf[m++]);
      for (k = 0; k < ncount; k++) {
        kk = npartner[j]++;
        partner[j][kk] = static_cast<tagint> (buf[m++]);
        memcpy(&shearpartner[j][dnum*kk],&buf[m],dnumbytes);
        m += dnum;
      }
    }
  } else error->all(FLERR,"Unsupported comm mode in shear history");
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixShearHistory::pack_exchange(int i, double *buf)
{
  // NOTE: how do I know comm buf is big enough if extreme # of touching neighs
  // Comm::BUFEXTRA may need to be increased

  int m = 0;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    memcpy(&buf[m],&shearpartner[i][dnum*n],dnumbytes);
    m += dnum;
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixShearHistory::unpack_exchange(int nlocal, double *buf)
{
  // allocate new chunks from ipage,dpage for incoming values

  int m = 0;
  npartner[nlocal] = static_cast<int> (buf[m++]);
  maxtouch = MAX(maxtouch,npartner[nlocal]);
  partner[nlocal] = ipage->get(npartner[nlocal]);
  shearpartner[nlocal] = dpage->get(dnum*npartner[nlocal]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<tagint> (buf[m++]);
    memcpy(&shearpartner[nlocal][dnum*n],&buf[m],dnumbytes);
    m += dnum;
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixShearHistory::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    memcpy(&buf[m],&shearpartner[i][dnum*n],dnumbytes);
    m += dnum;
  }
  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixShearHistory::unpack_restart(int nlocal, int nth)
{
  // ipage = NULL if being called from granular pair style init()

  if (ipage == NULL) allocate_pages();

  // skip to Nth set of extra values

  double **extra = atom->extra;

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  // allocate new chunks from ipage,dpage for incoming values

  npartner[nlocal] = static_cast<int> (extra[nlocal][m++]);
  maxtouch = MAX(maxtouch,npartner[nlocal]);
  partner[nlocal] = ipage->get(npartner[nlocal]);
  shearpartner[nlocal] = dpage->get(dnum*npartner[nlocal]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<tagint> (extra[nlocal][m++]);
    memcpy(&shearpartner[nlocal][dnum*n],&extra[nlocal][m],dnumbytes);
    m += dnum;
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixShearHistory::maxsize_restart()
{
  // maxtouch_all = max # of touching partners across all procs

  int maxtouch_all;
  MPI_Allreduce(&maxtouch,&maxtouch_all,1,MPI_INT,MPI_MAX,world);
  return (dnum+1)*maxtouch_all + 2;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixShearHistory::size_restart(int nlocal)
{
  return (dnum+1)*npartner[nlocal] + 2;
}
