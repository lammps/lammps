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

#include "string.h"
#include "compute_msd.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMSD::ComputeMSD(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute msd command");

  vector_flag = 1;
  size_vector = 4;
  extvector = 0;

  // optional args

  comflag = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"com") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute msd command");
      if (strcmp(arg[iarg+1],"no") == 0) comflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) comflag = 1;
      else error->all(FLERR,"Illegal compute msd command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute msd command");
  }

  // create a new fix store/state style with or without com keyword
  // id = compute-ID + store_state, fix group = compute group

  int n = strlen(id) + strlen("_store_state") + 1;
  id_fix = new char[n];
  strcpy(id_fix,id);
  strcat(id_fix,"_store_state");

  char **newarg = new char*[9];
  newarg[0] = id_fix;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "store/state";
  newarg[3] = (char *) "0";
  newarg[4] = (char *) "xu";
  newarg[5] = (char *) "yu";
  newarg[6] = (char *) "zu";
  newarg[7] = (char *) "com";
  newarg[8] = (char *) "yes";
  if (comflag) modify->add_fix(9,newarg);
  else modify->add_fix(7,newarg);
  delete [] newarg;

  vector = new double[4];
}

/* ---------------------------------------------------------------------- */

ComputeMSD::~ComputeMSD()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);

  delete [] id_fix;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeMSD::init()
{
  // set fix which stores original atom coords

  int ifix = modify->find_fix(id_fix);
  if (ifix < 0) error->all(FLERR,"Could not find compute msd fix ID");
  fix = modify->fix[ifix];

  // nmsd = # of atoms in group

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  nmsd = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) nmsd++;

  int nmsd_all;
  MPI_Allreduce(&nmsd,&nmsd_all,1,MPI_INT,MPI_SUM,world);
  nmsd = nmsd_all;

  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void ComputeMSD::compute_vector()
{
  invoked_vector = update->ntimestep;

  // cm = current center of mass

  double cm[3];
  if (comflag) group->xcm(igroup,masstotal,cm);
  else cm[0] = cm[1] = cm[2] = 0.0;

  // dx,dy,dz = displacement of atom from original position
  // original unwrapped position is stored by fix
  // relative to center of mass if comflag is set
  // for triclinic, need to unwrap current atom coord via h matrix

  double **xoriginal = fix->array_atom;

  double **x = atom->x;
  int *mask = atom->mask;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  
  double dx,dy,dz;
  int xbox,ybox,zbox;

  double msd[4];
  msd[0] = msd[1] = msd[2] = msd[3] = 0.0;

  if (domain->triclinic == 0) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	xbox = (image[i] & 1023) - 512;
	ybox = (image[i] >> 10 & 1023) - 512;
	zbox = (image[i] >> 20) - 512;
	dx = x[i][0] + xbox*xprd - cm[0] - xoriginal[i][0];
	dy = x[i][1] + ybox*yprd - cm[1] - xoriginal[i][1];
	dz = x[i][2] + zbox*zprd - cm[2] - xoriginal[i][2];
	msd[0] += dx*dx;
	msd[1] += dy*dy;
	msd[2] += dz*dz;
	msd[3] += dx*dx + dy*dy + dz*dz;
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	xbox = (image[i] & 1023) - 512;
	ybox = (image[i] >> 10 & 1023) - 512;
	zbox = (image[i] >> 20) - 512;
	dx = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox - 
	  cm[0] - xoriginal[i][0];
	dy = x[i][1] + h[1]*ybox + h[3]*zbox - cm[1] - xoriginal[i][1];
	dz = x[i][2] + h[2]*zbox - cm[2] - xoriginal[i][2];
	msd[0] += dx*dx;
	msd[1] += dy*dy;
	msd[2] += dz*dz;
	msd[3] += dx*dx + dy*dy + dz*dz;
      }
  }

  MPI_Allreduce(msd,vector,4,MPI_DOUBLE,MPI_SUM,world);
  if (nmsd) {
    vector[0] /= nmsd;
    vector[1] /= nmsd;
    vector[2] /= nmsd;
    vector[3] /= nmsd;
  }
}
