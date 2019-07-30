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

#include "fix_neigh_history.h"
#include <mpi.h>
#include <cstring>
#include "my_page.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{DEFAULT,NPARTNER,PERPARTNER}; // also set in fix neigh/history/omp

/* ---------------------------------------------------------------------- */

FixNeighHistory::FixNeighHistory(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  pair(NULL), npartner(NULL), partner(NULL), valuepartner(NULL),
  ipage_atom(NULL), dpage_atom(NULL), ipage_neigh(NULL), dpage_neigh(NULL)
{
  if (narg != 4) error->all(FLERR,"Illegal fix NEIGH_HISTORY command");

  restart_peratom = 1;
  create_attribute = 1;

  newton_pair = force->newton_pair;

  dnum = force->inumeric(FLERR,arg[3]);
  dnumbytes = dnum * sizeof(double);

  zeroes = new double[dnum];
  for (int i = 0; i < dnum; i++) zeroes[i] = 0.0;

  onesided = 0;
  if (strcmp(id,"LINE_NEIGH_HISTORY") == 0) onesided = 1;
  if (strcmp(id,"TRI_NEIGH_HISTORY") == 0) onesided = 1;

  if (newton_pair) comm_reverse = 1;   // just for single npartner value
                                       // variable-size history communicated via
                                       // reverse_comm_fix_variable()

  // perform initial allocation of atom-based arrays
  // register with atom class

  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  pgsize = oneatom = 0;

  // other per-atom vectors

  firstflag = NULL;
  firstvalue = NULL;
  maxatom = 0;

  // per-atom and per-neighbor data structs

  ipage_atom = NULL;
  dpage_atom = NULL;
  ipage_neigh = NULL;
  dpage_neigh = NULL;

  // initialize npartner to 0 so neighbor list creation is OK the 1st time

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) npartner[i] = 0;
  maxpartner = 0;

  nlocal_neigh = nall_neigh = 0;
  commflag = DEFAULT;
}

/* ---------------------------------------------------------------------- */

FixNeighHistory::~FixNeighHistory()
{
  if (copymode) return;

  // unregister this fix so atom class doesn't invoke it any more

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  delete [] zeroes;

  memory->sfree(firstflag);
  memory->sfree(firstvalue);

  memory->destroy(npartner);
  memory->sfree(partner);
  memory->sfree(valuepartner);

  delete [] ipage_atom;
  delete [] dpage_atom;
  delete [] ipage_neigh;
  delete [] dpage_neigh;

  // to better detect use-after-delete errors

  firstflag = NULL;
  firstvalue = NULL;

  pair = NULL;
  npartner = NULL;
  partner = NULL;
  valuepartner = NULL;
}

/* ---------------------------------------------------------------------- */

int FixNeighHistory::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= MIN_PRE_EXCHANGE;
  mask |= POST_NEIGHBOR;
  mask |= MIN_POST_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNeighHistory::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Neighbor history requires atoms have IDs");

  allocate_pages();
}

/* ----------------------------------------------------------------------
   create pages if first time or if neighbor pgsize/oneatom has changed
   note that latter could cause neighbor history info to be discarded
------------------------------------------------------------------------- */

void FixNeighHistory::allocate_pages()
{
  int create = 0;
  if (ipage_atom == NULL) create = 1;
  if (pgsize != neighbor->pgsize) create = 1;
  if (oneatom != neighbor->oneatom) create = 1;

  if (create) {
    delete [] ipage_atom;
    delete [] dpage_atom;
    delete [] ipage_neigh;
    delete [] dpage_neigh;

    pgsize = neighbor->pgsize;
    oneatom = neighbor->oneatom;
    int nmypage = comm->nthreads;
    ipage_atom = new MyPage<tagint>[nmypage];
    dpage_atom = new MyPage<double>[nmypage];
    ipage_neigh = new MyPage<int>[nmypage];
    dpage_neigh = new MyPage<double>[nmypage];
    for (int i = 0; i < nmypage; i++) {
      ipage_atom[i].init(oneatom,pgsize);
      dpage_atom[i].init(dnum*oneatom,dnum*pgsize);
      ipage_neigh[i].init(oneatom,pgsize);
      dpage_neigh[i].init(dnum*oneatom,dnum*pgsize);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNeighHistory::setup_post_neighbor()
{
  post_neighbor();
}

/* ----------------------------------------------------------------------
   copy partner info from neighbor data structs (NDS) to atom arrays
   should be called whenever NDS store current history info
     and need to transfer the info to owned atoms
   e.g. when atoms migrate to new procs, new neigh list built, or between runs
     when atoms may be added or deleted (NDS becomes out-of-date)
   the next post_neighbor() will put this info back into new NDS
   called during run before atom exchanges, including for restart files
   called at end of run via post_run()
   do not call during setup of run (setup_pre_exchange)
     b/c there is no guarantee of a current NDS (even on continued run)
   if run command does a 2nd run with pre = no, then no neigh list
     will be built, but old neigh list will still have the info
   onesided and newton on and newton off versions
------------------------------------------------------------------------- */

void FixNeighHistory::pre_exchange()
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

void FixNeighHistory::pre_exchange_onesided()
{
  int i,j,ii,jj,m,n,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *allflags;
  double *allvalues,*onevalues;

  // NOTE: all operations until very end are with nlocal_neigh <= current nlocal
  // b/c previous neigh list was built with nlocal_neigh
  // nlocal can be larger if other fixes added atoms at this pre_exchange()

  // clear two paged data structures

  ipage_atom->reset();
  dpage_atom->reset();

  // 1st loop over neighbor list, I = sphere, J = tri
  // only calculate npartner for owned spheres

  for (i = 0; i < nlocal_neigh; i++) npartner[i] = 0;

  tagint *tag = atom->tag;
  NeighList *list = pair->list;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    allflags = firstflag[i];

    for (jj = 0; jj < jnum; jj++)
      if (allflags[jj]) npartner[i]++;
  }

  // get page chunks to store partner IDs and values for owned atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    n = npartner[i];
    partner[i] = ipage_atom->get(n);
    valuepartner[i] = dpage_atom->get(dnum*n);
    if (partner[i] == NULL || valuepartner[i] == NULL)
      error->one(FLERR,"Neighbor history overflow, boost neigh_modify one");
  }

  // 2nd loop over neighbor list, I = sphere, J = tri
  // store partner IDs and values for owned+ghost atoms
  // re-zero npartner to use as counter

  for (i = 0; i < nlocal_neigh; i++) npartner[i] = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    allflags = firstflag[i];
    allvalues = firstvalue[i];

    for (jj = 0; jj < jnum; jj++) {
      if (allflags[jj]) {
        onevalues = &allvalues[dnum*jj];
        j = jlist[jj];
        j &= NEIGHMASK;
        m = npartner[i]++;
        partner[i][m] = tag[j];
        memcpy(&valuepartner[i][dnum*m],onevalues,dnumbytes);
      }
    }
  }

  // set maxpartner = max # of partners of any owned atom
  // bump up comm->maxexchange_fix if necessary

  maxpartner = 0;
  for (i = 0; i < nlocal_neigh; i++) maxpartner = MAX(maxpartner,npartner[i]);
  comm->maxexchange_fix = MAX(comm->maxexchange_fix,(dnum+1)*maxpartner+1);

  // zero npartner values from previous nlocal_neigh to current nlocal

  int nlocal = atom->nlocal;
  for (i = nlocal_neigh; i < nlocal; i++) npartner[i] = 0;
}

/* ----------------------------------------------------------------------
   newton ON version
   performs reverse comm to acquire partner values from ghost atoms
------------------------------------------------------------------------- */

void FixNeighHistory::pre_exchange_newton()
{
  int i,j,ii,jj,m,n,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *allflags;
  double *allvalues,*onevalues,*jvalues;

  // NOTE: all operations until very end are with
  //   nlocal_neigh  <= current nlocal and nall_neigh
  // b/c previous neigh list was built with nlocal_neigh & nghost_neigh
  // nlocal can be larger if other fixes added atoms at this pre_exchange()

  // clear two paged data structures

  ipage_atom->reset();
  dpage_atom->reset();

  // 1st loop over neighbor list
  // calculate npartner for owned+ghost atoms

  for (i = 0; i < nall_neigh; i++) npartner[i] = 0;

  tagint *tag = atom->tag;
  NeighList *list = pair->list;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    allflags = firstflag[i];

    for (jj = 0; jj < jnum; jj++) {
      if (allflags[jj]) {
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

  // get page chunks to store partner IDs and values for owned+ghost atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    n = npartner[i];
    partner[i] = ipage_atom->get(n);
    valuepartner[i] = dpage_atom->get(dnum*n);
    if (partner[i] == NULL || valuepartner[i] == NULL) {
      error->one(FLERR,"Neighbor history overflow, boost neigh_modify one");
    }
  }

  for (i = nlocal_neigh; i < nall_neigh; i++) {
    n = npartner[i];
    partner[i] = ipage_atom->get(n);
    valuepartner[i] = dpage_atom->get(dnum*n);
    if (partner[i] == NULL || valuepartner[i] == NULL) {
      error->one(FLERR,"Neighbor history overflow, boost neigh_modify one");
    }
  }

  // 2nd loop over neighbor list
  // store partner IDs and values for owned+ghost atoms
  // re-zero npartner to use as counter

  for (i = 0; i < nall_neigh; i++) npartner[i] = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    allflags = firstflag[i];
    allvalues = firstvalue[i];

    for (jj = 0; jj < jnum; jj++) {
      if (allflags[jj]) {
        onevalues = &allvalues[dnum*jj];
        j = jlist[jj];
        j &= NEIGHMASK;
        m = npartner[i]++;
        partner[i][m] = tag[j];
        memcpy(&valuepartner[i][dnum*m],onevalues,dnumbytes);
        m = npartner[j]++;
        partner[j][m] = tag[i];
        jvalues = &valuepartner[j][dnum*m];
        if (pair->nondefault_history_transfer) 
          pair->transfer_history(onevalues,jvalues);
        else for (n = 0; n < dnum; n++) jvalues[n] = -onevalues[n];
      }
    }
  }

  // perform reverse comm to augment
  // owned atom partner/valuepartner with ghost info
  // use variable variant b/c size of packed data can be arbitrarily large
  //   if many touching neighbors for large particle

  commflag = PERPARTNER;
  comm->reverse_comm_fix_variable(this);

  // set maxpartner = max # of partners of any owned atom
  // bump up comm->maxexchange_fix if necessary

  maxpartner = 0;
  for (i = 0; i < nlocal_neigh; i++) maxpartner = MAX(maxpartner,npartner[i]);
  comm->maxexchange_fix = MAX(comm->maxexchange_fix,(dnum+1)*maxpartner+1);

  // zero npartner values from previous nlocal_neigh to current nlocal

  int nlocal = atom->nlocal;
  for (i = nlocal_neigh; i < nlocal; i++) npartner[i] = 0;
}

/* ----------------------------------------------------------------------
   newton OFF version
   do not need partner values from ghost atoms
   assume J values are negative of I values
------------------------------------------------------------------------- */

void FixNeighHistory::pre_exchange_no_newton()
{
  int i,j,ii,jj,m,n,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *allflags;
  double *allvalues,*onevalues,*jvalues;

  // NOTE: all operations until very end are with nlocal_neigh <= current nlocal
  // b/c previous neigh list was built with nlocal_neigh
  // nlocal can be larger if other fixes added atoms at this pre_exchange()

  // clear two paged data structures

  ipage_atom->reset();
  dpage_atom->reset();

  // 1st loop over neighbor list
  // calculate npartner for owned atoms

  for (i = 0; i < nlocal_neigh; i++) npartner[i] = 0;

  tagint *tag = atom->tag;
  NeighList *list = pair->list;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    allflags = firstflag[i];

    for (jj = 0; jj < jnum; jj++) {
      if (allflags[jj]) {
        npartner[i]++;
        j = jlist[jj];
        j &= NEIGHMASK;
        if (j < nlocal_neigh) npartner[j]++;
      }
    }
  }

  // get page chunks to store partner IDs and values for owned atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    n = npartner[i];
    partner[i] = ipage_atom->get(n);
    valuepartner[i] = dpage_atom->get(dnum*n);
    if (partner[i] == NULL || valuepartner[i] == NULL)
      error->one(FLERR,"Neighbor history overflow, boost neigh_modify one");
  }

  // 2nd loop over neighbor list
  // store partner IDs and values for owned+ghost atoms
  // re-zero npartner to use as counter

  for (i = 0; i < nlocal_neigh; i++) npartner[i] = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    allflags = firstflag[i];
    allvalues = firstvalue[i];

    for (jj = 0; jj < jnum; jj++) {
      if (allflags[jj]) {
        onevalues = &allvalues[dnum*jj];
        j = jlist[jj];
        j &= NEIGHMASK;
        m = npartner[i]++;
        partner[i][m] = tag[j];
        memcpy(&valuepartner[i][dnum*m],onevalues,dnumbytes);
        if (j < nlocal_neigh) {
          m = npartner[j]++;
          partner[j][m] = tag[i];
          jvalues = &valuepartner[j][dnum*m];
          if (pair->nondefault_history_transfer) 
            pair->transfer_history(onevalues, jvalues);
          else for (n = 0; n < dnum; n++) jvalues[n] = -onevalues[n];
        }
      }
    }
  }

  // set maxpartner = max # of partners of any owned atom
  // bump up comm->maxexchange_fix if necessary

  maxpartner = 0;
  for (i = 0; i < nlocal_neigh; i++) maxpartner = MAX(maxpartner,npartner[i]);
  comm->maxexchange_fix = MAX(comm->maxexchange_fix,(dnum+1)*maxpartner+1);

  // zero npartner values from previous nlocal_neigh to current nlocal

  int nlocal = atom->nlocal;
  for (i = nlocal_neigh; i < nlocal; i++) npartner[i] = 0;
}

/* ---------------------------------------------------------------------- */

void FixNeighHistory::min_pre_exchange()
{
  pre_exchange();
}

/* ----------------------------------------------------------------------
   called after neighbor list is build
   recover history info stored temporarily in per-atom partner lists
     and store afresh in per-neighbor firstflag and firstvalue lists
------------------------------------------------------------------------- */

void FixNeighHistory::post_neighbor()
{
  int i,j,m,ii,jj,nn,np,inum,jnum,rflag;
  tagint jtag;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *allflags;
  double *allvalues;

  // store atom counts used for new neighbor list which was just built

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  nlocal_neigh = nlocal;
  nall_neigh = nall;

  // realloc firstflag and firstvalue if needed

  if (maxatom < nlocal) {
    memory->sfree(firstflag);
    memory->sfree(firstvalue);
    maxatom = nall;
    firstflag = (int **)
      memory->smalloc(maxatom*sizeof(int *),"neighbor_history:firstflag");
    firstvalue = (double **)
      memory->smalloc(maxatom*sizeof(double *),"neighbor_history:firstvalue");
  }

  // loop over newly built neighbor list
  // repopulate entire per-neighbor data structs
  //   whether with old-neigh partner info or zeroes

  ipage_neigh->reset();
  dpage_neigh->reset();

  tagint *tag = atom->tag;
  NeighList *list = pair->list;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    firstflag[i] = allflags = ipage_neigh->get(jnum);
    firstvalue[i] = allvalues = dpage_neigh->get(jnum*dnum);
    np = npartner[i];
    nn = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      rflag = sbmask(j) | pair->beyond_contact;
      j &= NEIGHMASK;
      jlist[jj] = j;

      // rflag = 1 if r < radsum in npair_size() method
      // preserve neigh history info if tag[j] is in old-neigh partner list
      // this test could be more geometrically precise for two sphere/line/tri

      if (rflag) {
        jtag = tag[j];
        for (m = 0; m < np; m++)
          if (partner[i][m] == jtag) break;
        if (m < np) {
          allflags[jj] = 1;
          memcpy(&allvalues[nn],&valuepartner[i][dnum*m],dnumbytes);
        } else {
          allflags[jj] = 0;
          memcpy(&allvalues[nn],zeroes,dnumbytes);
        }
      } else {
        allflags[jj] = 0;
        memcpy(&allvalues[nn],zeroes,dnumbytes);
      }
      nn += dnum;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNeighHistory::min_post_neighbor()
{
  post_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixNeighHistory::post_run()
{
  pre_exchange();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixNeighHistory::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);    // npartner
  bytes += nmax * sizeof(tagint *);     // partner
  bytes += nmax * sizeof(double *);     // valuepartner
  bytes += maxatom * sizeof(int *);     // firstflag
  bytes += maxatom * sizeof(double *);  // firstvalue

  int nmypage = comm->nthreads;
  for (int i = 0; i < nmypage; i++) {
    bytes += ipage_atom[i].size();
    bytes += dpage_atom[i].size();
    bytes += ipage_neigh[i].size();
    bytes += dpage_neigh[i].size();
  }

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixNeighHistory::grow_arrays(int nmax)
{
  memory->grow(npartner,nmax,"neighbor_history:npartner");
  partner = (tagint **) memory->srealloc(partner,nmax*sizeof(tagint *),
                                         "neighbor_history:partner");
  valuepartner = (double **) memory->srealloc(valuepartner,
                                              nmax*sizeof(double *),
                                              "neighbor_history:valuepartner");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixNeighHistory::copy_arrays(int i, int j, int /*delflag*/)
{
  // just copy pointers for partner and valuepartner
  // b/c can't overwrite chunk allocation inside ipage_atom,dpage_atom
  // incoming atoms in unpack_exchange just grab new chunks
  // so are orphaning chunks for migrating atoms
  // OK, b/c will reset ipage_atom,dpage_atom on next reneighboring

  npartner[j] = npartner[i];
  partner[j] = partner[i];
  valuepartner[j] = valuepartner[i];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixNeighHistory::set_arrays(int i)
{
  npartner[i] = 0;
}

/* ----------------------------------------------------------------------
   only called by Comm::reverse_comm_fix_variable for PERPARTNER mode
------------------------------------------------------------------------- */

int FixNeighHistory::pack_reverse_comm_size(int n, int first)
{
  int i,last;
  int dnump1 = dnum + 1;

  int m = 0;
  last = first + n;

  for (i = first; i < last; i++)
    m += 1 + dnump1*npartner[i];

  return m;
}

/* ----------------------------------------------------------------------
   two modes: NPARTNER and PERPARTNER
------------------------------------------------------------------------- */

int FixNeighHistory::pack_reverse_comm(int n, int first, double *buf)
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
        memcpy(&buf[m],&valuepartner[i][dnum*k],dnumbytes);
        m += dnum;
      }
    }
  } else error->all(FLERR,"Unsupported comm mode in neighbor history");

  return m;
}

/* ----------------------------------------------------------------------
   two modes: NPARTNER and PERPARTNER
------------------------------------------------------------------------- */

void FixNeighHistory::unpack_reverse_comm(int n, int *list, double *buf)
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
        memcpy(&valuepartner[j][dnum*kk],&buf[m],dnumbytes);
        m += dnum;
      }
    }
  } else error->all(FLERR,"Unsupported comm mode in neighbor history");
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixNeighHistory::pack_exchange(int i, double *buf)
{
  // NOTE: how do I know comm buf is big enough if extreme # of touching neighs
  // Comm::BUFEXTRA may need to be increased

  int m = 0;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    memcpy(&buf[m],&valuepartner[i][dnum*n],dnumbytes);
    m += dnum;
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixNeighHistory::unpack_exchange(int nlocal, double *buf)
{
  // allocate new chunks from ipage_atom,dpage_atom for incoming values

  int m = 0;
  npartner[nlocal] = static_cast<int> (buf[m++]);
  maxpartner = MAX(maxpartner,npartner[nlocal]);
  partner[nlocal] = ipage_atom->get(npartner[nlocal]);
  valuepartner[nlocal] = dpage_atom->get(dnum*npartner[nlocal]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<tagint> (buf[m++]);
    memcpy(&valuepartner[nlocal][dnum*n],&buf[m],dnumbytes);
    m += dnum;
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixNeighHistory::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    memcpy(&buf[m],&valuepartner[i][dnum*n],dnumbytes);
    m += dnum;
  }
  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixNeighHistory::unpack_restart(int nlocal, int nth)
{
  // ipage_atom = NULL if being called from granular pair style init()

  if (ipage_atom == NULL) allocate_pages();

  // skip to Nth set of extra values

  double **extra = atom->extra;

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  // allocate new chunks from ipage_atom,dpage_atom for incoming values

  npartner[nlocal] = static_cast<int> (extra[nlocal][m++]);
  maxpartner = MAX(maxpartner,npartner[nlocal]);
  partner[nlocal] = ipage_atom->get(npartner[nlocal]);
  valuepartner[nlocal] = dpage_atom->get(dnum*npartner[nlocal]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<tagint> (extra[nlocal][m++]);
    memcpy(&valuepartner[nlocal][dnum*n],&extra[nlocal][m],dnumbytes);
    m += dnum;
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixNeighHistory::maxsize_restart()
{
  // maxpartner_all = max # of touching partners across all procs

  int maxpartner_all;
  MPI_Allreduce(&maxpartner,&maxpartner_all,1,MPI_INT,MPI_MAX,world);
  return (dnum+1)*maxpartner_all + 2;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixNeighHistory::size_restart(int nlocal)
{
  return (dnum+1)*npartner[nlocal] + 2;
}
