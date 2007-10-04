/* -----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------ */

#include "math.h"
#include "string.h"
#include "dihedral_hybrid.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EXTRA 1000

/* ---------------------------------------------------------------------- */

DihedralHybrid::DihedralHybrid(LAMMPS *lmp) : Dihedral(lmp)
{
  nstyles = 0;
}

/* ---------------------------------------------------------------------- */

DihedralHybrid::~DihedralHybrid()
{
  if (nstyles) {
    for (int i = 0; i < nstyles; i++) delete styles[i];
    delete [] styles;
    for (int i = 0; i < nstyles; i++) delete [] keywords[i];
    delete [] keywords;
  }

  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(map);
    delete [] ndihedrallist;
    delete [] maxdihedral;
    for (int i = 0; i < nstyles; i++)
      memory->destroy_2d_int_array(dihedrallist[i]);
    delete [] dihedrallist;
  }
}

/* ---------------------------------------------------------------------- */

void DihedralHybrid::compute(int eflag, int vflag)
{
  int i,m,n;

  // save ptrs to original dihedrallist

  int ndihedrallist_orig = neighbor->ndihedrallist;
  int **dihedrallist_orig = neighbor->dihedrallist;

  // if this is re-neighbor step, create sub-style dihedrallists
  // ndihedrallist[] = length of each sub-style list
  // realloc sub-style dihedrallist if necessary
  // load sub-style dihedrallist with 5 values from original dihedrallist

  if (neighbor->ago == 0) {
    for (m = 0; m < nstyles; m++) ndihedrallist[m] = 0;
    for (i = 0; i < ndihedrallist_orig; i++)
      ndihedrallist[map[dihedrallist_orig[i][4]]]++;
    for (m = 0; m < nstyles; m++) {
      if (ndihedrallist[m] > maxdihedral[m]) {
	memory->destroy_2d_int_array(dihedrallist[m]);
	maxdihedral[m] = ndihedrallist[m] + EXTRA;
	dihedrallist[m] = (int **)
	  memory->create_2d_int_array(maxdihedral[m],5,
				      "dihedral_hybrid:dihedrallist");
      }
      ndihedrallist[m] = 0;
    }
    for (i = 0; i < ndihedrallist_orig; i++) {
      m = map[dihedrallist_orig[i][4]];
      n = ndihedrallist[m];
      dihedrallist[m][n][0] = dihedrallist_orig[i][0];
      dihedrallist[m][n][1] = dihedrallist_orig[i][1];
      dihedrallist[m][n][2] = dihedrallist_orig[i][2];
      dihedrallist[m][n][3] = dihedrallist_orig[i][3];
      dihedrallist[m][n][4] = dihedrallist_orig[i][4];
      ndihedrallist[m]++;
    }
  }
  
  // call each sub-style's compute function
  // must set neighbor->dihedrallist to sub-style dihedrallist before call
  // accumulate sub-style energy,virial in hybrid's energy,virial

  energy = 0.0;
  eng_vdwl = eng_coul = 0.0;
  if (vflag) for (n = 0; n < 6; n++) virial[n] = 0.0;

  for (m = 0; m < nstyles; m++) {
    if (styles[m] == NULL) continue;
    neighbor->ndihedrallist = ndihedrallist[m];
    neighbor->dihedrallist = dihedrallist[m];
    styles[m]->compute(eflag,vflag);
    if (eflag) {
      energy += styles[m]->energy;
      eng_vdwl += styles[m]->eng_vdwl;
      eng_coul += styles[m]->eng_coul;
    }
    if (vflag) for (n = 0; n < 6; n++) virial[n] += styles[m]->virial[n];
  }

  // restore ptrs to original dihedrallist

  neighbor->ndihedrallist = ndihedrallist_orig;
  neighbor->dihedrallist = dihedrallist_orig;
}

/* ---------------------------------------------------------------------- */

void DihedralHybrid::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  map = (int *) memory->smalloc((n+1)*sizeof(int),"dihedral:map");
  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;

  ndihedrallist = new int[nstyles];
  maxdihedral = new int[nstyles];
  dihedrallist = new int**[nstyles];
  for (int m = 0; m < nstyles; m++) maxdihedral[m] = 0;
  for (int m = 0; m < nstyles; m++) dihedrallist[m] = NULL;
}

/* ----------------------------------------------------------------------
   create one dihedral style for each arg in list
------------------------------------------------------------------------- */

void DihedralHybrid::settings(int narg, char **arg)
{
  nstyles = narg;
  styles = new Dihedral*[nstyles];
  keywords = new char*[nstyles];

  for (int m = 0; m < nstyles; m++) {
    for (int i = 0; i < m; i++)
      if (strcmp(arg[m],arg[i]) == 0) 
	error->all("Dihedral style hybrid cannot use same dihedral style twice");
    if (strcmp(arg[m],"hybrid") == 0) 
      error->all("Dihedral style hybrid cannot have hybrid as an argument");
    styles[m] = force->new_dihedral(arg[m]);
    keywords[m] = new char[strlen(arg[m])+1];
    strcpy(keywords[m],arg[m]);
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one type
---------------------------------------------------------------------- */

void DihedralHybrid::coeff(int which, int narg, char **arg)
{
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->ndihedraltypes,ilo,ihi);

  // 2nd arg = dihedral style name (harmonic, etc)

  int m;
  for (m = 0; m < nstyles; m++)
    if (strcmp(arg[1],keywords[m]) == 0) break;
  if (m == nstyles) error->all("Dihedral coeff for hybrid has invalid style");

  // set low-level coefficients for each dihedraltype
  // replace 2nd arg with i, call coeff() with no 1st arg
  // if sub-style is NULL for "none", still set setflag

  for (int i = ilo; i <= ihi; i++) {
    sprintf(arg[1],"%d",i);
    map[i] = m;
    if (styles[m]) styles[m]->coeff(which,narg-1,&arg[1]);
    if (styles[m] == NULL) setflag[i] = 1;
    else setflag[i] = styles[m]->setflag[i];
  }
}

/* ---------------------------------------------------------------------- */

void DihedralHybrid::init_style()
{
  for (int m = 0; m < nstyles; m++)
    if (styles[m]) styles[m]->init_style();
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void DihedralHybrid::write_restart(FILE *fp)
{
  fwrite(&nstyles,sizeof(int),1,fp);

  int n;
  for (int m = 0; m < nstyles; m++) {
    n = strlen(keywords[m]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(keywords[m],sizeof(char),n,fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void DihedralHybrid::read_restart(FILE *fp)
{
  allocate();

  int me = comm->me;
  if (me == 0) fread(&nstyles,sizeof(int),1,fp);
  MPI_Bcast(&nstyles,1,MPI_INT,0,world);
  styles = new Dihedral*[nstyles];
  keywords = new char*[nstyles];
  
  int n;
  for (int m = 0; m < nstyles; m++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    keywords[m] = new char[n];
    if (me == 0) fread(keywords[m],sizeof(char),n,fp);
    MPI_Bcast(keywords[m],n,MPI_CHAR,0,world);
    styles[m] = force->new_dihedral(keywords[m]);
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double DihedralHybrid::memory_usage()
{
  double bytes = 0.0;
  for (int m = 0; m < nstyles; m++) bytes += maxdihedral[m]*5 * sizeof(int);
  for (int m = 0; m < nstyles; m++) 
    if (styles[m]) bytes += styles[m]->memory_usage();
  return bytes;
}
