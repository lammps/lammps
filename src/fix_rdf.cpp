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
#include "fix_rdf.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "output.h"
#include "group.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixRDF::FixRDF(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8 || (narg-6) % 2) error->all("Illegal fix rdf command");

  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix rdf command");
  first = 1;

  MPI_Comm_rank(world,&me);
  if (me == 0) {
    fp = fopen(arg[4],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix rdf file %s",arg[4]);
      error->one(str);
    }
  }

  maxbin = atoi(arg[5]);

  n_rdf_pairs = 0;
  rdfpair = memory->create_2d_int_array(atom->ntypes+1,atom->ntypes+1,
					"rdf:rdfpair");

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = 1; j <= atom->ntypes; j++)
      rdfpair[i][j] = 0;

  int itype,jtype;
  for (int i = 6; i < narg; i+=2) {
    itype = atoi(arg[i]);
    jtype = atoi(arg[i+1]);
    if (itype < 1 || jtype < 1 || itype > atom->ntypes || jtype > atom->ntypes)
      error->all("Invalid atom types in fix rdf command");
    n_rdf_pairs++;
    rdfpair[itype][jtype] = n_rdf_pairs;
  }

  hist = memory->create_2d_int_array(n_rdf_pairs,maxbin,"rdf:hist");
  hist_all = memory->create_2d_int_array(n_rdf_pairs,maxbin,"rdf:hist_all");

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

  gr_ave = memory->create_2d_double_array(n_rdf_pairs,maxbin,"rdf:gr_ave");
  ncoord_ave = memory->create_2d_double_array(n_rdf_pairs,maxbin,
					      "rdf:nccord_ave");
}

/* ---------------------------------------------------------------------- */

FixRDF::~FixRDF()
{
  if (me == 0) fclose(fp);

  memory->destroy_2d_int_array(rdfpair);
  memory->destroy_2d_int_array(hist);
  memory->destroy_2d_int_array(hist_all);
  delete [] nrdfatoms;
  memory->destroy_2d_double_array(gr_ave);
  memory->destroy_2d_double_array(ncoord_ave);
}

/* ---------------------------------------------------------------------- */

int FixRDF::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRDF::init()
{
  if (force->pair) delr = force->pair->cutforce / maxbin;
  else error->all("Fix rdf requires a pair style be defined");

  // need an occasional half neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->occasional = 1;
  
  delrinv = 1.0/delr;
  nframes = 0;

  // set running averages to 0.0

  int i,j,bin,irdf;

  for (i = 1; i <= atom->ntypes; i++)
    for (j = 1; j <= atom->ntypes; j++)
      if (rdfpair[i][j]) {
	irdf = rdfpair[i][j] - 1;
	for (bin = 0; bin < maxbin; bin++)
          gr_ave[irdf][bin] = ncoord_ave[irdf][bin] = 0.0;         
      }
}

/* ---------------------------------------------------------------------- */

void FixRDF::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixRDF::setup(int vflag)
{
  if (first) end_of_step();
  first = 0;
}

/* ---------------------------------------------------------------------- */

void FixRDF::end_of_step()
{
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

  for (int i = 0; i < n_rdf_pairs; i++)
    for (int j = 0; j < maxbin; j++)
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
          bin = static_cast<int> (r * delrinv);
	  if (bin >= maxbin) continue;

	  if (ipair) hist[ipair-1][bin]++;
	  if (newton_pair || j < nlocal)
	    if (jpair) hist[jpair-1][bin]++;
	}
      }
    }
  }

  // sum histogram across procs

  MPI_Allreduce(hist[0],hist_all[0],n_rdf_pairs*maxbin,MPI_INT,MPI_SUM,world);
  nframes++;

  if (me == 0) {

    // print RDF for current snapshot

    double rlower,rcenter,rupper,nideal,gr;
    double PI = 4.0*atan(1.0);
    double constant = 4.0*PI / (3.0*domain->xprd*domain->yprd*domain->zprd);
    int irdf;

    fprintf(fp,"# Timestep %d\n",update->ntimestep);
    fprintf(fp,"%s","r");

    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = 1; j <= atom->ntypes; j++)
        if (rdfpair[i][j])
          fprintf(fp,", %s%d,%d%s, %s%d,%d%s",
		  "g(",i,j,",r)","ncoord(",i,j,",r)");

    double *ncoord = new double[n_rdf_pairs];   
    for (int i = 0; i < n_rdf_pairs; i++) ncoord[i] = 0;

    for (bin = 0; bin < maxbin; bin++) {
      rlower = bin*delr;
      rcenter = rlower + 0.5*delr;
      rupper = rlower + delr;
      fprintf(fp,"\n%f ",rcenter); 
      for (int i = 1; i <= atom->ntypes; i++)
        for (int j = 1; j <= atom->ntypes; j++)
          if (rdfpair[i][j]) {
            nideal = constant*nrdfatoms[j] *
              (rupper*rupper*rupper - rlower*rlower*rlower);
            irdf = rdfpair[i][j] - 1;
            if (nrdfatoms[i]*nideal != 0.0) gr = hist_all[irdf][bin] / (nrdfatoms[i]*nideal);
            else gr = 0.0;
            ncoord[irdf] += gr*nideal;
            fprintf(fp,"%f %f ",gr,ncoord[irdf]);  
            gr_ave[irdf][bin] = 
	      ((nframes-1)*gr_ave[irdf][bin] + gr) / nframes;
            ncoord_ave[irdf][bin] = 
	      ((nframes-1)*ncoord_ave[irdf][bin] + ncoord[irdf]) / nframes;
          }
    }

    fprintf(fp,"\n");
    delete [] ncoord;
 
    // if last time in run that RDF is computed, print running averages

    if (update->ntimestep + nevery > update->endstep) {

      fprintf(fp,"# Run Average\n");
      fprintf(fp,"%s","r");

      for (int i = 1; i <= atom->ntypes; i++)
        for (int j = 1; j <= atom->ntypes; j++)
          if (rdfpair[i][j])
            fprintf(fp,", %s%d,%d%s, %s%d,%d%s",
		    "g(",i,j,",r)","ncoord(",i,j,",r)");
	    
      for (bin = 0; bin < maxbin; bin++) {
        rlower = bin*delr;
        rcenter = rlower + 0.5*delr;
        rupper = rlower + delr;
        fprintf(fp,"\n%f ",rcenter); 
        for (int i = 1; i <= atom->ntypes; i++)
          for (int j = 1; j <= atom->ntypes; j++)
            if (rdfpair[i][j]) {
              irdf = rdfpair[i][j] - 1;
              fprintf(fp,"%f %f ",gr_ave[irdf][bin],ncoord_ave[irdf][bin]);  
            }  
      }
      fprintf(fp,"\n");
    }

    fflush(fp);
  }
}
