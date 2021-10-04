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

#include "pair_grid.h"
#include "pair_sna_grid.h"
#include "sna.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

PairSNAGrid::PairSNAGrid(LAMMPS *lmp) :
  PairGrid(lmp), cutsq(nullptr),
  radelem(nullptr), wjelem(nullptr)
{
}

/* ---------------------------------------------------------------------- */

PairSNAGrid::~PairSNAGrid()
{
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);
  delete snaptr;

  if (chemflag) memory->destroy(map);
}

/* ---------------------------------------------------------------------- */

void PairSNAGrid::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Pair sna/grid requires a pair style be defined");

  if (cutmax > force->pair->cutforce)
    error->all(FLERR,"Pair sna/grid cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  snaptr->init();
}

/* ---------------------------------------------------------------------- */

void PairSNAGrid::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void PairSNAGrid::compute(int eflag, int vflag)
{

  // compute sna for each gridpoint

  double** const x = atom->x;
  const int* const mask = atom->mask;
  int * const type = atom->type;
  const int ntotal = atom->nlocal + atom->nghost;

  // insure rij, inside, and typej are of size jnum
  
  snaptr->grow_rij(ntotal);

  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
	double xgrid[3];
	grid2x(ix, iy, iz, xgrid);
	const double xtmp = xgrid[0];
	const double ytmp = xgrid[1];
	const double ztmp = xgrid[2];

	// currently, all grid points are type 1
	
	const int itype = 1;
	int ielem = 0;
	if (chemflag)
	  ielem = map[itype];
	const double radi = radelem[itype];

	// rij[][3] = displacements between atom I and those neighbors
	// inside = indices of neighbors of I within cutoff
	// typej = types of neighbors of I within cutoff

	int ninside = 0;
	for (int j = 0; j < ntotal; j++) {

	  const double delx = xtmp - x[j][0];
	  const double dely = ytmp - x[j][1];
	  const double delz = ztmp - x[j][2];
	  const double rsq = delx*delx + dely*dely + delz*delz;
	  int jtype = type[j];
	  int jelem = 0;
	  if (chemflag)
	    jelem = map[jtype];
	  if (rsq < cutsq[jtype][jtype] && rsq > 1e-20) {
	    snaptr->rij[ninside][0] = delx;
	    snaptr->rij[ninside][1] = dely;
	    snaptr->rij[ninside][2] = delz;
	    snaptr->inside[ninside] = j;
	    snaptr->wj[ninside] = wjelem[jtype];
	    snaptr->rcutij[ninside] = 2.0*radelem[jtype]*rcutfac;
	    snaptr->element[ninside] = jelem; // element index for multi-element snap
	    ninside++;
	  }
	}

	snaptr->compute_ui(ninside, ielem);
	snaptr->compute_zi();
	snaptr->compute_bi(ielem);

	// linear contributions

	for (int icoeff = 0; icoeff < ncoeff; icoeff++)
	  gridlocal[ndesc_base+icoeff][iz][iy][ix] = 
	    snaptr->blist[icoeff];

	// quadratic contributions

	if (quadraticflag) {
	  int ncount = ncoeff;
	  for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
	    double bveci = snaptr->blist[icoeff];
	    gridlocal[ndesc_base+ncount++][iz][iy][ix] = 
	      0.5*bveci*bveci;
	    for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++)
	      gridlocal[ndesc_base+ncount++][iz][iy][ix] = 
		bveci*snaptr->blist[jcoeff];
	  }
	}
      }
}


/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSNAGrid::settings(int narg, char ** arg)
{
  double rfac0, rmin0;
  int twojmax, switchflag, bzeroflag, bnormflag, wselfallflag;

  // call base class first

  PairGrid::settings(narg, arg);
    
  // skip over arguments used by base class
  // so that argument positions are identical to
  // regular per-atom compute
  
  arg += nargbase;
  narg -= nargbase;

  int ntypes = atom->ntypes;
  int nargmin = 3+2*ntypes;

  if (narg < nargmin) error->all(FLERR,"Illegal pair sna/grid command");

  // default values

  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  quadraticflag = 0;
  chemflag = 0;
  bnormflag = 0;
  wselfallflag = 0;
  nelements = 1;
  
  // process required arguments

  memory->create(radelem,ntypes+1,"sna/grid/local:radelem"); // offset by 1 to match up with types
  memory->create(wjelem,ntypes+1,"sna/grid/local:wjelem");

  rcutfac = atof(arg[0]);
  rfac0 = atof(arg[1]);
  twojmax = atoi(arg[2]);

  for(int i = 0; i < ntypes; i++)
    radelem[i+1] = atof(arg[3+i]);
  for(int i = 0; i < ntypes; i++)
    wjelem[i+1] = atof(arg[3+ntypes+i]);

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq,ntypes+1,ntypes+1,"sna/grid/local:cutsq");
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
    printf("%d %d %d %s\n",iarg,narg,nargbase,arg[iarg]);
    if (strcmp(arg[iarg],"rmin0") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      rmin0 = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"switchflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      switchflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"bzeroflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      bzeroflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"quadraticflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      quadraticflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"chem") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      chemflag = 1;
      memory->create(map,ntypes+1,"compute_sna_grid_local:map");
      nelements = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      for (int i = 0; i < ntypes; i++) {
        int jelem = utils::inumeric(FLERR,arg[iarg+2+i],false,lmp);
        if (jelem < 0 || jelem >= nelements)
          error->all(FLERR,"Illegal pair sna/grid command");
        map[i+1] = jelem;
      }
      iarg += 2+ntypes;
    } else if (strcmp(arg[iarg],"bnormflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      bnormflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"wselfallflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      wselfallflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal pair sna/grid command");

  }

}
  
/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double PairSNAGrid::memory_usage()
{
  double nbytes = snaptr->memory_usage();    // SNA object
  int n = atom->ntypes+1;
  nbytes += (double)n*sizeof(int);            // map

  return nbytes;
}

