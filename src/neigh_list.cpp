// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "neigh_list.h"
#include "my_page.h"    // IWYU pragma: keep
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "my_page.h"
#include "memory.h"

using namespace LAMMPS_NS;

static constexpr int PGDELTA = 1;

/* ---------------------------------------------------------------------- */

NeighList::NeighList(LAMMPS *lmp) : Pointers(lmp)
{
  // initializations

  maxatom = 0;

  inum = gnum = 0;
  ilist = nullptr;
  numneigh = nullptr;
  firstneigh = nullptr;

  // defaults, but may be reset by post_constructor()

  occasional = 0;
  ghost = 0;
  ssa = 0;
  history = 0;
  respaouter = 0;
  respamiddle = 0;
  respainner = 0;
  copy = 0;
  trim = 0;
  copymode = 0;

  // ptrs

  iskip = nullptr;
  ijskip = nullptr;

  listcopy = nullptr;
  listskip = nullptr;
  listfull = nullptr;

  fix_bond = nullptr;

  ipage = nullptr;

  // extra rRESPA lists

  inum_inner = gnum_inner = 0;
  ilist_inner = nullptr;
  numneigh_inner = nullptr;
  firstneigh_inner = nullptr;

  inum_middle = gnum_middle = 0;
  ilist_middle = nullptr;
  numneigh_middle = nullptr;
  firstneigh_middle = nullptr;

  ipage_inner = nullptr;
  ipage_middle = nullptr;

  // Kokkos package

  kokkos = 0;
  kk2cpu = 0;
  execution_space = Host;

  // DPD-REACT package

  np = nullptr;

  requestor = nullptr;
  requestor_type = NeighList::NONE;
}

/* ---------------------------------------------------------------------- */

NeighList::~NeighList()
{
  if (copymode) return;

  if (!copy || trim || kk2cpu) {
    memory->destroy(ilist);
    memory->destroy(numneigh);
    memory->sfree(firstneigh);
    delete [] ipage;
  }

  if (respainner) {
    memory->destroy(ilist_inner);
    memory->destroy(numneigh_inner);
    memory->sfree(firstneigh_inner);
    delete [] ipage_inner;
  }

  if (respamiddle) {
    memory->destroy(ilist_middle);
    memory->destroy(numneigh_middle);
    memory->sfree(firstneigh_middle);
    delete [] ipage_middle;
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
   respaouter -> set all 3 outer/middle/inner flags
   bond -> set fix_bond to Fix that made the request
------------------------------------------------------------------------- */

void NeighList::post_constructor(NeighRequest *nq)
{
  // copy request settings used by list itself

  occasional = nq->occasional;
  ghost = nq->ghost;
  ssa = nq->ssa;
  history = nq->history;
  respaouter = nq->respaouter;
  respamiddle = nq->respamiddle;
  respainner = nq->respainner;
  copy = nq->copy;
  trim = nq->trim;
  id = nq->id;

  if (nq->copy) {
    listcopy = neighbor->lists[nq->copylist];
    if (listcopy->kokkos && !this->kokkos)
      kk2cpu = 1;
  }

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

  if (respainner) {
    ipage_inner = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++)
      ipage_inner[i].init(oneatom,pgsize,PGDELTA);
  }

  if (respamiddle) {
    ipage_middle = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++)
      ipage_middle[i].init(oneatom,pgsize,PGDELTA);
  }
}

/* ----------------------------------------------------------------------
   grow per-atom data to allow for nlocal/nall atoms
   triggered by neighbor list build
   not called if a copy list
------------------------------------------------------------------------- */

void NeighList::grow(int nlocal, int nall)
{
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

  if (respainner) {
    memory->destroy(ilist_inner);
    memory->destroy(numneigh_inner);
    memory->sfree(firstneigh_inner);
    memory->create(ilist_inner,maxatom,"neighlist:ilist_inner");
    memory->create(numneigh_inner,maxatom,"neighlist:numneigh_inner");
    firstneigh_inner = (int **) memory->smalloc(maxatom*sizeof(int *),
                                                "neighlist:firstneigh_inner");
  }

  if (respamiddle) {
    memory->destroy(ilist_middle);
    memory->destroy(numneigh_middle);
    memory->sfree(firstneigh_middle);
    memory->create(ilist_middle,maxatom,"neighlist:ilist_middle");
    memory->create(numneigh_middle,maxatom,"neighlist:numneigh_middle");
    firstneigh_middle = (int **) memory->smalloc(maxatom*sizeof(int *),
                                                 "neighlist:firstneigh_middle");
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
  printf("  %d = respaouter\n",rq->respaouter);
  printf("  %d = respamiddle\n",rq->respamiddle);
  printf("  %d = respainner\n",rq->respainner);
  printf("  %d = bond\n",rq->bond);
  printf("  %d = omp\n",rq->omp);
  printf("  %d = intel\n",rq->intel);
  printf("  %d = kokkos host\n",rq->kokkos_host);
  printf("  %d = kokkos device\n",rq->kokkos_device);
  printf("  %d = ssa flag\n",ssa);
  printf("\n");
  printf("  %d = skip flag\n",rq->skip);
  printf("  %d = off2on\n",rq->off2on);
  printf("  %d = copy flag\n",rq->copy);
  printf("  %d = trim flag\n",rq->trim);
  printf("  %d = kk2cpu flag\n",kk2cpu);
  printf("  %d = half/full\n",rq->halffull);
  printf("\n");
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   if growflag = 0, maxatom & maxpage will also be 0
   if stencilflag = 0, maxstencil * maxstencil_multi will also be 0
------------------------------------------------------------------------- */

double NeighList::memory_usage()
{
  double bytes = 0;
  bytes += memory->usage(ilist,maxatom);
  bytes += memory->usage(numneigh,maxatom);
  bytes += (double)maxatom * sizeof(int *);

  int nmypage = comm->nthreads;

  if (ipage) {
    for (int i = 0; i < nmypage; i++)
      bytes += ipage[i].size();
  }

  if (respainner) {
    bytes += memory->usage(ilist_inner,maxatom);
    bytes += memory->usage(numneigh_inner,maxatom);
    bytes += (double)maxatom * sizeof(int *);
    if (ipage_inner) {
      for (int i = 0; i < nmypage; i++)
        bytes += ipage_inner[i].size();
    }
  }

  if (respamiddle) {
    bytes += memory->usage(ilist_middle,maxatom);
    bytes += memory->usage(numneigh_middle,maxatom);
    bytes += (double)maxatom * sizeof(int *);
    if (ipage_middle) {
      for (int i = 0; i < nmypage; i++)
        bytes += ipage_middle[i].size();
    }
  }

  return bytes;
}
