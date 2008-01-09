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
#include "neighbor.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define PGDELTA 1

enum{NSQ,BIN,MULTI};     // also in neighbor.cpp

/* ---------------------------------------------------------------------- */

NeighList::NeighList(LAMMPS *lmp, int size) : Pointers(lmp)
{
  maxlocal = 0;
  pgsize = size;

  ilist = NULL;
  numneigh = NULL;
  firstneigh = NULL;
  firstdouble = NULL;

  maxpage = 0;
  pages = NULL;
  dpages = NULL;
  dnum = 0;

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

  maxstencil = 0;
  stencil = NULL;

  maxstencil_multi = 0;
  nstencil_multi = NULL;
  stencil_multi = NULL;
  distsq_multi = NULL;
}

/* ---------------------------------------------------------------------- */

NeighList::~NeighList()
{
  if (!listcopy) {
    memory->sfree(ilist);
    memory->sfree(numneigh);
    memory->sfree(firstneigh);
    memory->sfree(firstdouble);

    for (int i = 0; i < maxpage; i++) memory->sfree(pages[i]);
    memory->sfree(pages);
    if (dnum) {
      for (int i = 0; i < maxpage; i++) memory->sfree(dpages[i]);
      memory->sfree(dpages);
    }
  }

  delete [] iskip;
  memory->destroy_2d_int_array(ijskip);

  if (maxstencil) memory->sfree(stencil);
  if (maxstencil_multi) {
    for (int i = 1; i <= atom->ntypes; i++) {
      memory->sfree(stencil_multi[i]);
      memory->sfree(distsq_multi[i]);
    }
    delete [] nstencil_multi;
    delete [] stencil_multi;
    delete [] distsq_multi;
  }
}

/* ----------------------------------------------------------------------
   grow atom arrays to allow for nmax atoms
   triggered by more atoms on a processor
------------------------------------------------------------------------- */

void NeighList::grow(int nmax)
{
  // skip if grow not needed

  if (nmax <= maxlocal) return;
  maxlocal = nmax;

  memory->sfree(ilist);
  memory->sfree(numneigh);
  memory->sfree(firstneigh);
  ilist = (int *)
    memory->smalloc(maxlocal*sizeof(int),"neighlist:ilist");
  numneigh = (int *)
    memory->smalloc(maxlocal*sizeof(int),"neighlist:numneigh");
  firstneigh = (int **)
    memory->smalloc(maxlocal*sizeof(int *),"neighlist:firstneigh");

  if (dnum) 
    firstdouble = (double **)
      memory->smalloc(maxlocal*sizeof(double *),"neighlist:firstdouble");
}

/* ----------------------------------------------------------------------
   insure stencils are large enough for smax bins
   style = BIN or MULTI
------------------------------------------------------------------------- */

void NeighList::stencil_allocate(int smax, int style)
{
  int i;

  if (style == BIN) {
    if (smax > maxstencil) {
      maxstencil = smax;
      memory->sfree(stencil);
      stencil = (int *) memory->smalloc(smax*sizeof(int),"neighlist:stencil");
    }

  } else {
    int n = atom->ntypes;
    if (maxstencil_multi == 0) {
      nstencil_multi = new int[n+1];
      stencil_multi = new int*[n+1];
      distsq_multi = new double*[n+1];
      for (i = 1; i <= n; i++) {
	nstencil_multi[i] = 0;
	stencil_multi[i] = NULL;
	distsq_multi[i] = NULL;
      }
    }
    if (smax > maxstencil_multi) {
      maxstencil_multi = smax;
      for (i = 1; i <= n; i++) {
	memory->sfree(stencil_multi[i]);
	memory->sfree(distsq_multi[i]);
	stencil_multi[i] = (int *)
	  memory->smalloc(smax*sizeof(int),"neighlist:stencil_multi");
	distsq_multi[i] = (double *) memory->smalloc(smax*sizeof(double),
						     "neighlist:distsq_multi");
      }
    }
  }
}

/* ----------------------------------------------------------------------
   add PGDELTA pages to neighbor list
------------------------------------------------------------------------- */

int **NeighList::add_pages()
{
  int npage = maxpage;
  maxpage += PGDELTA;

  pages = (int **) 
    memory->srealloc(pages,maxpage*sizeof(int *),"neighlist:pages");
  for (int i = npage; i < maxpage; i++)
    pages[i] = (int *) memory->smalloc(pgsize*sizeof(int),
				       "neighlist:pages[i]");

  if (dnum) {
    dpages = (double **) 
      memory->srealloc(dpages,maxpage*sizeof(double *),"neighlist:dpages");
    for (int i = npage; i < maxpage; i++)
      dpages[i] = (double *) memory->smalloc(dnum*pgsize*sizeof(double),
					     "neighlist:dpages[i]");
  }

  return pages;
}

/* ----------------------------------------------------------------------
   copy skip info from request rq into list's iskip,ijskip
------------------------------------------------------------------------- */

void NeighList::copy_skip_info(int *rq_iskip, int **rq_ijskip)
{
  int ntypes = atom->ntypes;
  iskip = new int[ntypes+1];
  ijskip = memory->create_2d_int_array(ntypes+1,ntypes+1,
				       "neigh_list:ijskip");
  int i,j;
  for (i = 1; i <= ntypes; i++) iskip[i] = rq_iskip[i];
  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      ijskip[i][j] = rq_ijskip[i][j];
}

/* ----------------------------------------------------------------------
   print attributes of this list and associated request
------------------------------------------------------------------------- */

void NeighList::print_attributes()
{
  if (comm->me != 0) return;

  NeighRequest *rq = neighbor->requests[index];

  printf("Neighbor list/request %d:\n",index);
  printf("  %d = build flag\n",buildflag);
  printf("  %d = grow flag\n",growflag);
  printf("  %d = stencil flag\n",stencilflag);
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
  printf("  %d = occasional\n",rq->occasional);
  printf("  %d = dnum\n",rq->dnum);
  printf("  %d = copy\n",rq->copy);
  printf("  %d = skip\n",rq->skip);
  printf("  %d = otherlist\n",rq->otherlist);
  printf("  %p = listskip\n",listskip);
  printf("\n");
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   if growflag = 0, maxlocal & maxpage will also be 0
   if stencilflag = 0, maxstencil * maxstencil_multi will also be 0
------------------------------------------------------------------------- */

double NeighList::memory_usage()
{
  double bytes = 0.0;
  bytes += 2 * maxlocal * sizeof(int);
  bytes += maxlocal * sizeof(int *);
  bytes += maxpage*pgsize * sizeof(int);

  if (dnum) {
    bytes += maxlocal * sizeof(double *);
    bytes += dnum * maxpage*pgsize * sizeof(double);
  }

  if (maxstencil) bytes += maxstencil * sizeof(int);
  if (maxstencil_multi) {
    bytes += atom->ntypes * maxstencil_multi * sizeof(int);
    bytes += atom->ntypes * maxstencil_multi * sizeof(double);
  }

  return bytes;
}
