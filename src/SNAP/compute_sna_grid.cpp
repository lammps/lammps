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

#include "compute_grid.h"
#include "compute_sna_grid.h"
#include <cstring>
#include <cstdlib>
#include "sna.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

ComputeSNAGrid::ComputeSNAGrid(LAMMPS *lmp, int narg, char **arg) :
  ComputeGrid(lmp, narg, arg), cutsq(NULL), list(NULL), sna(NULL),
  radelem(NULL), wjelem(NULL)
{
  double rmin0, rfac0;
  int twojmax, switchflag, bzeroflag;
  radelem = NULL;
  wjelem = NULL;

  // skip over arguments used by base class
  // so that argument positions are identical to 
  // regular per-atom compute

  arg += nargbase;
  narg -= nargbase;

  int ntypes = atom->ntypes;
  int nargmin = 6+2*ntypes;

  if (narg < nargmin) error->all(FLERR,"Illegal compute sna/grid command");

  // default values

  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  quadraticflag = 0;

  // offset by 1 to match up with types

  memory->create(radelem,ntypes+1,"sna/grid:radelem");
  memory->create(wjelem,ntypes+1,"sna/grid:wjelem");

  rcutfac = atof(arg[3]);
  rfac0 = atof(arg[4]);
  twojmax = atoi(arg[5]);

  for(int i = 0; i < ntypes; i++)
    radelem[i+1] = atof(arg[6+i]);
  for(int i = 0; i < ntypes; i++)
    wjelem[i+1] = atof(arg[6+ntypes+i]);

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq,ntypes+1,ntypes+1,"sna/grid:cutsq");
  for(int i = 1; i <= ntypes; i++) {
    cut = 2.0*radelem[i]*rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[i][i] = cut*cut;
    for(int j = i+1; j <= ntypes; j++) {
      cut = (radelem[i]+radelem[j])*rcutfac;
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }
  }

  // process optional args

  int iarg = nargmin;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"rmin0") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute sna/grid command");
      rmin0 = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"switchflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute sna/grid command");
      switchflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"bzeroflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute sna/grid command");
      bzeroflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"quadraticflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute sna/grid command");
      quadraticflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal compute sna/grid command");

  }

  snaptr = new SNA(lmp,rfac0,twojmax,
                   rmin0,switchflag,bzeroflag);

  ncoeff = snaptr->ncoeff;
  size_array_cols = ncoeff;
  if (quadraticflag) size_array_cols += (ncoeff*(ncoeff+1))/2;
  array_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeSNAGrid::~ComputeSNAGrid()
{
  memory->destroy(sna);
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);
  delete snaptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSNAGrid::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute sna/grid requires a pair style be defined");

  if (cutmax > force->pair->cutforce)
    error->all(FLERR,"Compute sna/grid cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"sna/grid") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute sna/grid");
  snaptr->init();
}

/* ---------------------------------------------------------------------- */

void ComputeSNAGrid::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSNAGrid::compute_array()
{
  invoked_array = update->ntimestep;

//   // invoke full neighbor list (will copy or build if necessary)

//   neighbor->build_one(list);

//   const int inum = list->inum;
//   const int* const ilist = list->ilist;
//   const int* const numneigh = list->numneigh;
//   int** const firstneigh = list->firstneigh;

  int * const type = atom->type;

  // compute sna for each gridpoint

  double** const x = atom->x;
  const int* const mask = atom->mask;
  const int ntotal = atom->nlocal + atom->nghost;

  // insure rij, inside, and typej are of size jnum
  
  snaptr->grow_rij(ntotal);

  printf("ngridfull = %d\n",ngridfull);
  for (int igrid = 0; igrid < ngridfull; igrid++) {
    double rtmp[3];
    igridfull2x(igrid, rtmp);
    const double xtmp = rtmp[0];
    const double ytmp = rtmp[1];
    const double ztmp = rtmp[2];

    // rij[][3] = displacements between atom I and those neighbors
    // inside = indices of neighbors of I within cutoff
    // typej = types of neighbors of I within cutoff

    int ninside = 0;
    for (int j = 0; j < ntotal; j++) {

      // check that j is in compute group

      if (!(mask[j] & groupbit)) continue;

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];
      if (rsq < cutsq[jtype][jtype] && rsq>1e-20) {
	//        printf("ninside = %d\n",ninside);
        snaptr->rij[ninside][0] = delx;
        snaptr->rij[ninside][1] = dely;
        snaptr->rij[ninside][2] = delz;
        snaptr->inside[ninside] = j;
        snaptr->wj[ninside] = wjelem[jtype];
        snaptr->rcutij[ninside] = 2.0*radelem[jtype]*rcutfac;
        ninside++;
      }
    }

    snaptr->compute_ui(ninside);
    snaptr->compute_zi();
    snaptr->compute_bi();
    for (int icoeff = 0; icoeff < ncoeff; icoeff++)
      sna[igrid][icoeff] = snaptr->blist[icoeff];
    //    printf("igrid = %d %g %g %g %d B0 = %g\n",igrid,xtmp,ytmp,ztmp,ninside,sna[igrid][0]);
    if (quadraticflag) {
      int ncount = ncoeff;
      for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
        double bi = snaptr->blist[icoeff];

        // diagonal element of quadratic matrix

        sna[igrid][ncount++] = 0.5*bi*bi;

        // upper-triangular elements of quadratic matrix

        for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++)
          sna[igrid][ncount++] = bi*snaptr->blist[jcoeff];
      }
    }
  }
  //  gather_global_array();
  copy_local_grid();
}

/* ----------------------------------------------------------------------
   allocate array in base class and then set up pointers
------------------------------------------------------------------------- */

void ComputeSNAGrid::allocate()
{
  ComputeGrid::allocate();
  sna = gridfull;
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeSNAGrid::memory_usage()
{
  double bytes = size_array_rows*size_array_cols * sizeof(double); // grid
  bytes += snaptr->memory_usage();                        // SNA object

  return bytes;
}

