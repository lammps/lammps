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
  dnum = 0;

  // ptrs

  iskip = NULL;
  ijskip = NULL;

  listgranhistory = NULL;
  fix_history = NULL;

  respamiddle = 0;
  listinner = NULL;
  listmiddle = NULL;
  listfull = NULL;
  listcopy = NULL;
  listskip = NULL;

  ipage = NULL;
  dpage = NULL;

  // Kokkos package

  kokkos = 0;
  execution_space = Host;

  // USER-DPD package

  ndxAIR_ssa = NULL;
}

/* ---------------------------------------------------------------------- */

NeighList::~NeighList()
{
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

  if (ssa) {
    memory->sfree(ndxAIR_ssa);
  }
}

/* ----------------------------------------------------------------------
   adjust settings to match corresponding NeighRequest
   cannot do this in constructor b/c not all NeighLists are allocated yet
   copy -> set listcopy for list to copy from
   skip -> set listskip and create local copy of itype,ijtype
   halffull -> preceeding list must be full
   granhistory -> preceeding list must be gran, set its LGH and FH ptrs
   respaouter -> preceeding list(s) must be inner or middle, set LM/LI ptrs
------------------------------------------------------------------------- */

void NeighList::post_constructor(NeighRequest *nq)
{
  occasional = nq->occasional;
  ghost = nq->ghost;
  ssa = nq->ssa;
  copy = nq->copy;
  dnum = nq->dnum;

  if (nq->copy)
    listcopy = neighbor->lists[nq->otherlist];

  if (nq->skip) {
    listskip = neighbor->lists[nq->otherlist];
    int ntypes = atom->ntypes;
    iskip = new int[ntypes+1];
    memory->create(ijskip,ntypes+1,ntypes+1,"neigh_list:ijskip");
    int i,j;
    for (i = 1; i <= ntypes; i++) iskip[i] = nq->iskip[i];
    for (i = 1; i <= ntypes; i++)
      for (j = 1; j <= ntypes; j++)
        ijskip[i][j] = nq->ijskip[i][j];
  }

  if (nq->half_from_full) {
    NeighRequest *oq = neighbor->requests[nq->index-1];
    if (oq->full != 1) 
      error->all(FLERR,"Neighbor half-from-full list does not follow "
                 "full list");
    listfull = neighbor->lists[nq->index-1];
  }

  if (nq->granhistory) {
    NeighRequest *oq = neighbor->requests[nq->index-1];
    if (oq->gran != 1)
      error->all(FLERR,"Neighbor granhistory list does not follow gran list");
    neighbor->lists[nq->index-1]->listgranhistory = this;
    neighbor->lists[nq->index-1]->fix_history = nq->fix_history;
  }
  
  if (nq->respaouter) {
    NeighRequest *oq = neighbor->requests[nq->index-1];
    if (oq->respainner) {
      respamiddle = 0;
      listinner = neighbor->lists[nq->index-1];
    } else if (oq->respamiddle) {
      respamiddle = 1;
      listmiddle = neighbor->lists[nq->index-1];
      oq = neighbor->requests[nq->index-2];
      if (!oq->respainner)
        error->all(FLERR,"Neighbor respa outer list does not follow "
                   "respa list");
      listinner = neighbor->lists[nq->index-2];
    } else
      error->all(FLERR,"Neighbor respa outer list does not follow respa list");
  }
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
     gran calls grow() in granhistory
     respaouter calls grow() in respainner, respamiddle
   triggered by neighbor list build
------------------------------------------------------------------------- */

void NeighList::grow(int nlocal, int nall)
{
  // trigger grow() in children before possible return

  if (listgranhistory) listgranhistory->grow(nlocal,nall);
  if (listinner) listinner->grow(nlocal,nall);
  if (listmiddle) listmiddle->grow(nlocal,nall);

  // skip if data structs are already big enough

  if (ghost) {
    if (nall <= maxatom) return;
  } else {
    if (nlocal <= maxatom) return;
  }

  maxatom = atom->nmax;

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

  if (ssa) {
    if (ndxAIR_ssa) memory->sfree(ndxAIR_ssa);
    ndxAIR_ssa = (uint16_t (*)[8]) memory->smalloc(sizeof(uint16_t)*8*maxatom,
      "neighlist:ndxAIR_ssa");
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
  printf("  %d = occasional\n",occasional);
  printf("  %d = ghost flag\n",ghost);
  printf("  %d = ssa flag\n",ssa);
  printf("  %d = copy flag\n",copy);
  printf("  %d = dnum\n",dnum);
  printf("\n");
  printf("  %d = pair\n",rq->pair);
  printf("  %d = fix\n",rq->fix);
  printf("  %d = compute\n",rq->compute);
  printf("  %d = command\n",rq->command);
  printf("\n");
  printf("  %d = half\n",rq->half);
  printf("  %d = full\n",rq->full);
  printf("  %d = gran\n",rq->gran);
  printf("  %d = granhistory\n",rq->granhistory);
  printf("  %d = respainner\n",rq->respainner);
  printf("  %d = respamiddle\n",rq->respamiddle);
  printf("  %d = respaouter\n",rq->respaouter);
  printf("  %d = half_from_full\n",rq->half_from_full);
  printf("\n");
  printf("  %d = newton\n",rq->newton);
  printf("  %d = granonesided\n",rq->granonesided);
  printf("  %d = omp\n",rq->omp);
  printf("  %d = intel\n",rq->intel);
  printf("  %d = kokkos host\n",rq->kokkos_host);
  printf("  %d = kokkos device\n",rq->kokkos_device);
  printf("  %d = skip\n",rq->skip);
  printf("  %d = otherlist\n",rq->otherlist);
  printf("  %p = listskip\n",(void *)listskip);
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

  if (ndxAIR_ssa) bytes += sizeof(uint16_t) * 8 * maxatom;

  return bytes;
}
