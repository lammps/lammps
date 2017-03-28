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

#include "neigh_list.h"
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "my_page.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define PGDELTA 1

enum{NSQ,BIN,MULTI};     // also in Neighbor

/* ---------------------------------------------------------------------- */

NeighList::NeighList(LAMMPS *lmp) : Pointers(lmp)
{
  // initializations

  maxatom = 0;

  inum = gnum = 0;
  ilist = NULL;
  numneigh = NULL;
  firstneigh = NULL;
  firstdouble = NULL;

  // defaults, but may be reset by post_constructor()

  occasional = 0;
  ghost = 0;
  ssa = 0;
  copy = 0;
  copymode = 0;
  dnum = 0;

  // ptrs

  iskip = NULL;
  ijskip = NULL;

  listcopy = NULL;
  listskip = NULL;
  listfull = NULL;

  listhistory = NULL;
  fix_history = NULL;

  respamiddle = 0;
  listinner = NULL;
  listmiddle = NULL;

  fix_bond = NULL;

  ipage = NULL;
  dpage = NULL;

  // Kokkos package

  kokkos = 0;
  execution_space = Host;

  // USER-DPD package

  for (int i = 0; i < 8; i++) AIRct_ssa[i] = 0;
  np = NULL;
}

/* ---------------------------------------------------------------------- */

NeighList::~NeighList()
{
  if (copymode) return;
  if (!copy) {
    memory->destroy(ilist);
    memory->destroy(numneigh);
    memory->sfree(firstneigh);
    memory->sfree(firstdouble);

    delete [] ipage;
    delete [] dpage;
  }

  delete [] iskip;
  memory->destroy(ijskip);
}

/* ----------------------------------------------------------------------
   adjust settings to match corresponding NeighRequest
   cannot do this in constructor b/c not all NeighLists are allocated yet
   copy -> set listcopy for list to copy from
   skip -> set listskip for list to skip from, create copy of itype,ijtype
   halffull -> set listfull for full list to derive from
   history -> set LH and FH ptrs in partner list that uses the history info
   respaouter -> set listinner/listmiddle for other rRESPA lists
   bond -> set fix_bond to Fix that made the request
------------------------------------------------------------------------- */

void NeighList::post_constructor(NeighRequest *nq)
{
  // copy request settings used by list itself
  
  occasional = nq->occasional;
  ghost = nq->ghost;
  ssa = nq->ssa;
  copy = nq->copy;
  dnum = nq->dnum;

  if (nq->copy)
    listcopy = neighbor->lists[nq->copylist];

  if (nq->skip) {
    listskip = neighbor->lists[nq->skiplist];
    int ntypes = atom->ntypes;
    iskip = new int[ntypes+1];
    memory->create(ijskip,ntypes+1,ntypes+1,"neigh_list:ijskip");
    int i,j;
    for (i = 1; i <= ntypes; i++) iskip[i] = nq->iskip[i];
    for (i = 1; i <= ntypes; i++)
      for (j = 1; j <= ntypes; j++)
        ijskip[i][j] = nq->ijskip[i][j];
  }

  if (nq->halffull)
    listfull = neighbor->lists[nq->halffulllist];

  if (nq->history) {
    neighbor->lists[nq->historylist]->listhistory = this;
    int tmp;
    neighbor->lists[nq->historylist]->fix_history = 
      (Fix *) ((Pair *) nq->requestor)->extract("history",tmp);
  }
  
  if (nq->respaouter) {
    if (nq->respamiddlelist < 0) {
      respamiddle = 0;
      listinner = neighbor->lists[nq->respainnerlist];
    } else {
      respamiddle = 1;
      listmiddle = neighbor->lists[nq->respamiddlelist];
      listinner = neighbor->lists[nq->respainnerlist];
    }
  }

  if (nq->bond) fix_bond = (Fix *) nq->requestor;
}

/* ---------------------------------------------------------------------- */

void NeighList::setup_pages(int pgsize_caller, int oneatom_caller)
{
  pgsize = pgsize_caller;
  oneatom = oneatom_caller;

  int nmypage = comm->nthreads;
  ipage = new MyPage<int>[nmypage];
  for (int i = 0; i < nmypage; i++)
    ipage[i].init(oneatom,pgsize,PGDELTA);

  if (dnum) {
    dpage = new MyPage<double>[nmypage];
    for (int i = 0; i < nmypage; i++)
      dpage[i].init(dnum*oneatom,dnum*pgsize,PGDELTA);
  } else dpage = NULL;
}

/* ----------------------------------------------------------------------
   grow per-atom data to allow for nlocal/nall atoms
   for parent lists:
     also trigger grow in child list(s) which are not built themselves
     history calls grow() in listhistory
     respaouter calls grow() in respainner, respamiddle
   triggered by neighbor list build
   not called if a copy list
------------------------------------------------------------------------- */

void NeighList::grow(int nlocal, int nall)
{
  // trigger grow() in children before possible return

  if (listhistory) listhistory->grow(nlocal,nall);
  if (listinner) listinner->grow(nlocal,nall);
  if (listmiddle) listmiddle->grow(nlocal,nall);

  // skip if data structs are already big enough
  if (ssa) {
    if ((nlocal * 3) + nall <= maxatom) return;
  } else if (ghost) {
    if (nall <= maxatom) return;
  } else {
    if (nlocal <= maxatom) return;
  }

  if (ssa) maxatom = (nlocal * 3) + nall;
  else maxatom = atom->nmax;

  memory->destroy(ilist);
  memory->destroy(numneigh);
  memory->sfree(firstneigh);
  memory->create(ilist,maxatom,"neighlist:ilist");
  memory->create(numneigh,maxatom,"neighlist:numneigh");
  firstneigh = (int **) memory->smalloc(maxatom*sizeof(int *),
                                        "neighlist:firstneigh");
  if (dnum) {
    memory->sfree(firstdouble);
    firstdouble = (double **) memory->smalloc(maxatom*sizeof(double *),
                                              "neighlist:firstdouble");
  }
}

/* ----------------------------------------------------------------------
   print attributes of this list and associated request
------------------------------------------------------------------------- */

void NeighList::print_attributes()
{
  if (comm->me != 0) return;

  NeighRequest *rq = neighbor->requests[index];

  printf("Neighbor list/request %d:\n",index);
  printf("  %p = requestor ptr (instance %d id %d)\n",
         rq->requestor,rq->requestor_instance,rq->id);
  printf("  %d = pair\n",rq->pair);
  printf("  %d = fix\n",rq->fix);
  printf("  %d = compute\n",rq->compute);
  printf("  %d = command\n",rq->command);
  printf("  %d = neigh\n",rq->neigh);
  printf("\n");
  printf("  %d = half\n",rq->half);
  printf("  %d = full\n",rq->full);
  printf("\n");
  printf("  %d = occasional\n",occasional);
  printf("  %d = newton\n",rq->newton);
  printf("  %d = ghost flag\n",ghost);
  printf("  %d = size\n",rq->size);
  printf("  %d = history\n",rq->history);
  printf("  %d = granonesided\n",rq->granonesided);
  printf("  %d = respainner\n",rq->respainner);
  printf("  %d = respamiddle\n",rq->respamiddle);
  printf("  %d = respaouter\n",rq->respaouter);
  printf("  %d = bond\n",rq->bond);
  printf("  %d = omp\n",rq->omp);
  printf("  %d = intel\n",rq->intel);
  printf("  %d = kokkos host\n",rq->kokkos_host);
  printf("  %d = kokkos device\n",rq->kokkos_device);
  printf("  %d = ssa flag\n",ssa);
  printf("  %d = dnum\n",dnum);
  printf("\n");
  printf("  %d = skip flag\n",rq->skip);
  printf("  %d = off2on\n",rq->off2on);
  printf("  %d = copy flag\n",rq->copy);
  printf("  %d = half/full\n",rq->halffull);
  printf("  %d = history/partner\n",rq->history_partner);
  printf("\n");
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   if growflag = 0, maxatom & maxpage will also be 0
   if stencilflag = 0, maxstencil * maxstencil_multi will also be 0
------------------------------------------------------------------------- */

bigint NeighList::memory_usage()
{
  bigint bytes = 0;
  bytes += memory->usage(ilist,maxatom);
  bytes += memory->usage(numneigh,maxatom);
  bytes += maxatom * sizeof(int *);

  int nmypage = comm->nthreads;

  if (ipage) {
    for (int i = 0; i < nmypage; i++)
      bytes += ipage[i].size();
  }

  if (dnum && dpage) {
    for (int i = 0; i < nmypage; i++) {
      bytes += maxatom * sizeof(double *);
      bytes += dpage[i].size();
    }
  }

  return bytes;
}
