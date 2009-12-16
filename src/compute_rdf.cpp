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

/* ----------------------------------------------------------------------
   Contributing authors: Paul Crozier (SNL), Jeff Greathouse (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "compute_rdf.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeRDF::ComputeRDF(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 8 || (narg-6) % 2) error->all("Illegal compute rdf command");

  array_flag = 1;
  size_array_rows = 1;
  size_array_cols = 1;
  extarray = 1;

  maxbin = atoi(arg[5]);

  npairs = 0;
  rdfpair = memory->create_2d_int_array(atom->ntypes+1,atom->ntypes+1,
					"rdf:rdfpair");

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = 1; j <= atom->ntypes; j++)
      rdfpair[i][j] = 0;

  int itype,jtype;
  for (int i = 6; i < narg; i += 2) {
    itype = atoi(arg[i]);
    jtype = atoi(arg[i+1]);
    if (itype < 1 || jtype < 1 || itype > atom->ntypes || jtype > atom->ntypes)
      error->all("Invalid atom type in compute rdf command");
    npairs++;
    rdfpair[itype][jtype] = npairs;
  }

  hist = memory->create_2d_double_array(maxbin,npairs+1,"rdf:hist");
  array = memory->create_2d_double_array(maxbin,npairs+1,"rdf:array");

  int *nrdfatom = new int[atom->ntypes+1];
  for (int i = 1; i <= atom->ntypes; i++) nrdfatom[i] = 0;

  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) nrdfatom[type[i]]++;

  nrdfatoms = new int[atom->ntypes+1];
  MPI_Allreduce(&nrdfatom[1],&nrdfatoms[1],atom->ntypes,MPI_INT,MPI_SUM,world);
  delete [] nrdfatom;
}

/* ---------------------------------------------------------------------- */

ComputeRDF::~ComputeRDF()
{
  memory->destroy_2d_int_array(rdfpair);
  memory->destroy_2d_double_array(hist);
  delete [] nrdfatoms;
}

/* ---------------------------------------------------------------------- */

void ComputeRDF::init()
{
  if (force->pair) delr = force->pair->cutforce / maxbin;
  else error->all("Compute rdf requires a pair style be defined");
  delrinv = 1.0/delr;

  // need an occasional half neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeRDF::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeRDF::compute_array()
{
  invoked_array = update->ntimestep;

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int nall = atom->nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  int i,j,ii,jj,inum,jnum,itype,jtype,ipair,jpair,bin;
  double xtmp,ytmp,ztmp,delx,dely,delz,r;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // invoke half neighbor list (will copy or build if necessary)

  neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero the histogram counts

  for (int i = 0; i < maxbin; i++)
    for (int j = 0; j < npairs; j++)
      hist[i][j] = 0;

  // tally the RDF
  // both atom i and j must be in fix group
  // itype,jtype must have been specified by user
  // weighting factor must be != 0.0 for this pair
  //   could be 0 and still be in neigh list for long-range Coulombics
  // count the interaction once even if neighbor pair is stored on 2 procs
  // if itype = jtype, count the interaction twice

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];

	if (j >= nall) {
	  if (special_coul[j/nall] == 0.0 && special_lj[j/nall] == 0.0)
	    continue;
	  j %= nall;
	}

        if (mask[j] & groupbit) {
          jtype = type[j];
	  ipair = rdfpair[itype][jtype];
	  jpair = rdfpair[jtype][itype];
	  if (!ipair && !jpair) continue;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          r = sqrt(delx*delx + dely*dely + delz*delz);
          bin = static_cast<int> (r*delrinv);
	  if (bin >= maxbin) continue;

	  if (ipair) hist[bin][ipair-1]++;
	  if (newton_pair || j < nlocal)
	    if (jpair) hist[bin][jpair-1]++;
	}
      }
    }
  }

  // sum histogram across procs

  MPI_Allreduce(hist[0],array[0],maxbin*npairs,MPI_DOUBLE,MPI_SUM,world);
}
