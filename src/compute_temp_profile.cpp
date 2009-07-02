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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "compute_temp_profile.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

ComputeTempProfile::ComputeTempProfile(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg < 7) error->all("Illegal compute temp/profile command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;
  tempbias = 1;

  xflag = atoi(arg[3]);
  yflag = atoi(arg[4]);
  zflag = atoi(arg[5]);
  if (zflag && domain->dimension == 2)
    error->all("Compute temp/profile cannot use vz for 2d systemx");

  ncount = 0;
  ivx = ivy = ivz = 0;
  if (xflag) ivx = ncount++;
  if (yflag) ivy = ncount++;
  if (zflag) ivz = ncount++;

  nbinx = nbiny = nbinz = 1;

  if (strcmp(arg[6],"x") == 0) {
    if (narg != 8) error->all("Illegal compute temp/profile command");
    nbinx = atoi(arg[7]);
  } else if (strcmp(arg[6],"y") == 0) {
    if (narg != 8) error->all("Illegal compute temp/profile command");
    nbiny = atoi(arg[7]);
  } else if (strcmp(arg[6],"z") == 0) {
    if (narg != 8) error->all("Illegal compute temp/profile command");
    if (domain->dimension == 2)
      error->all("Compute temp/profile cannot bin z for 2d systems");
    nbinz = atoi(arg[7]);
  } else if (strcmp(arg[6],"xy") == 0) {
    if (narg != 9) error->all("Illegal compute temp/profile command");
    nbinx = atoi(arg[7]);
    nbiny = atoi(arg[8]);
  } else if (strcmp(arg[6],"yz") == 0) {
    if (narg != 9) error->all("Illegal compute temp/profile command");
    if (domain->dimension == 2)
      error->all("Compute temp/profile cannot bin z for 2d systems");
    nbiny = atoi(arg[7]);
    nbinz = atoi(arg[8]);
  } else if (strcmp(arg[6],"xz") == 0) {
    if (narg != 9) error->all("Illegal compute temp/profile command");
    if (domain->dimension == 2)
      error->all("Compute temp/profile cannot bin z for 2d systems");
    nbinx = atoi(arg[7]);
    nbinz = atoi(arg[8]);
  } else if (strcmp(arg[6],"xyz") == 0) {
    if (narg != 10) error->all("Illegal compute temp/profile command");
    if (domain->dimension == 2)
      error->all("Compute temp/profile cannot bin z for 2d systems");
    nbinx = atoi(arg[7]);
    nbiny = atoi(arg[8]);
    nbinz = atoi(arg[9]);
  } else error->all("Illegal compute temp/profile command");

  nbins = nbinx*nbiny*nbinz;
  if (nbins <= 0) error->all("Illegal compute temp/profile command");

  vbin = memory->create_2d_double_array(nbins,ncount+1,"temp/profile:vbin");
  binave = memory->create_2d_double_array(nbins,ncount+1,
					  "temp/profile:binave");
  
  maxatom = 0;
  bin = NULL;

  maxbias = 0;
  vbiasall = NULL;
  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeTempProfile::~ComputeTempProfile()
{
  memory->destroy_2d_double_array(vbin);
  memory->destroy_2d_double_array(binave);
  memory->sfree(bin);
  memory->destroy_2d_double_array(vbiasall);
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempProfile::init()
{
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  dof_compute();

  // ptrs to domain data

  box_change = domain->box_change;
  triclinic = domain->triclinic;
  periodicity = domain->periodicity;

  if (triclinic) {
    boxlo = domain->boxlo_lamda;
    boxhi = domain->boxhi_lamda;
    prd = domain->prd_lamda;
  } else {
    boxlo = domain->boxlo;
    boxhi = domain->boxhi;
    prd = domain->prd;
  }

  if (!box_change) bin_setup();
}

/* ---------------------------------------------------------------------- */

void ComputeTempProfile::dof_compute()
{
  double natoms = group->count(igroup);
  int nper = domain->dimension;
  dof = nper * natoms;
  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempProfile::compute_scalar()
{
  int ibin;
  double vthermal[3];

  invoked_scalar = update->ntimestep;

  bin_average();

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ibin = bin[i];
      if (xflag) vthermal[0] = v[i][0] - binave[ibin][ivx];
      else vthermal[0] = v[i][0];
      if (yflag) vthermal[1] = v[i][1] - binave[ibin][ivy];
      else vthermal[1] = v[i][1];
      if (zflag) vthermal[2] = v[i][2] - binave[ibin][ivz];
      else vthermal[2] = v[i][2];

      if (rmass)
	t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + 
	      vthermal[2]*vthermal[2]) * rmass[i];
      else
	t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + 
	      vthermal[2]*vthermal[2]) * mass[type[i]];
    }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempProfile::compute_vector()
{
  int i,ibin;
  double vthermal[3];

  invoked_vector = update->ntimestep;

  bin_average();

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ibin = bin[i];
      if (xflag) vthermal[0] = v[i][0] - binave[ibin][ivx];
      else vthermal[0] = v[i][0];
      if (yflag) vthermal[1] = v[i][1] - binave[ibin][ivy];
      else vthermal[1] = v[i][1];
      if (zflag) vthermal[2] = v[i][2] - binave[ibin][ivz];
      else vthermal[2] = v[i][2];

      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      t[0] += massone * vthermal[0]*vthermal[0];
      t[1] += massone * vthermal[1]*vthermal[1];
      t[2] += massone * vthermal[2]*vthermal[2];
      t[3] += massone * vthermal[0]*vthermal[1];
      t[4] += massone * vthermal[0]*vthermal[2];
      t[5] += massone * vthermal[1]*vthermal[2];
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempProfile::remove_bias(int i, double *v)
{
  int ibin = bin[i];
  if (xflag) {
    vbias[0] = binave[ibin][ivx];
    v[0] -= vbias[0];
  }
  if (yflag) {
    vbias[1] = binave[ibin][ivy];
    v[1] -= vbias[1];
  }
  if (zflag) {
    vbias[2] = binave[ibin][ivz];
    v[2] -= vbias[2];
  }
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempProfile::remove_bias_all()
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (nlocal > maxbias) {
    memory->destroy_2d_double_array(vbiasall);
    maxbias = atom->nmax;
    vbiasall = memory->create_2d_double_array(maxbias,3,
					      "compute/temp:vbiasall");
  }

  int ibin;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ibin = bin[i];
      if (xflag) {
	vbiasall[i][0] = binave[ibin][ivx];
	v[i][0] -= vbiasall[i][0];
      }
      if (yflag) {
	vbiasall[i][1] = binave[ibin][ivy];
	v[i][1] -= vbiasall[i][1];
      }
      if (zflag) {
	vbiasall[i][2] = binave[ibin][ivz];
	v[i][2] -= vbiasall[i][2];
      }
    }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempProfile::restore_bias(int i, double *v)
{
  if (xflag) v[0] += vbias[0];
  if (yflag) v[1] += vbias[1];
  if (zflag) v[2] += vbias[2];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempProfile::restore_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (xflag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	v[i][0] += vbiasall[i][0];
  }
  if (yflag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	v[i][1] += vbiasall[i][1];
  }
  if (zflag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	v[i][2] += vbiasall[i][2];
  }
}

/* ----------------------------------------------------------------------
   compute average velocity in each bin
------------------------------------------------------------------------- */

void ComputeTempProfile::bin_average()
{
  int i,j,ibin;

  if (box_change) bin_setup();
  bin_assign();

  // clear bins, including particle count

  for (i = 0; i < nbins; i++)
    for (j = 0; j <= ncount; j++)
      vbin[i][j] = 0.0;

  // sum each particle's velocity to appropriate bin

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ibin = bin[i];
      if (xflag) vbin[ibin][ivx] += v[i][0];
      if (yflag) vbin[ibin][ivy] += v[i][1];
      if (zflag) vbin[ibin][ivz] += v[i][2];
      vbin[ibin][ncount] += 1.0;
    }

  // sum bins across processors

  MPI_Allreduce(vbin[0],binave[0],nbins*(ncount+1),MPI_DOUBLE,MPI_SUM,world);

  // compute ave velocity in each bin, checking for no particles

  for (i = 0; i < nbins; i++)
    if (binave[i][ncount] > 0.0)
      for (j = 0; j < ncount; j++)
	binave[i][j] /= binave[i][ncount];
}

/* ----------------------------------------------------------------------
   set bin sizes, redo if box size changes
------------------------------------------------------------------------- */

void ComputeTempProfile::bin_setup()
{
  invdelta[0] = nbinx / prd[0];
  invdelta[1] = nbiny / prd[1];
  invdelta[2] = nbinz / prd[2];
}

/* ----------------------------------------------------------------------
   assign all atoms to bins
------------------------------------------------------------------------- */

void ComputeTempProfile::bin_assign()
{
  // reallocate bin array if necessary

  if (atom->nlocal > maxatom) {
    maxatom = atom->nmax;
    memory->sfree(bin);
    bin = (int *) memory->smalloc(maxatom*sizeof(int),"temp/profile:bin");
  }

  // assign each atom to a bin, accounting for PBC
  // if triclinic, do this in lamda space

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int ibinx,ibiny,ibinz;
  double coord;

  if (triclinic) domain->x2lamda(nlocal);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (nbinx > 1) {
	coord = x[i][0];
	if (periodicity[0]) {
	  if (coord < boxlo[0]) coord += prd[0];
	  if (coord >= boxhi[0]) coord -= prd[0];
	}
	ibinx = static_cast<int> ((coord - boxlo[0]) * invdelta[0]);
	ibinx = MAX(ibinx,0);
	ibinx = MIN(ibinx,nbinx-1);
      } else ibinx = 0;
      
      if (nbiny > 1) {
	coord = x[i][1];
	if (periodicity[1]) {
	  if (coord < boxlo[1]) coord += prd[1];
	  if (coord >= boxhi[1]) coord -= prd[1];
	}
	ibiny = static_cast<int> ((coord - boxlo[1]) * invdelta[1]);
	ibiny = MAX(ibiny,0);
	ibiny = MIN(ibiny,nbiny-1);
      } else ibiny = 0;
      
      if (nbinz > 1) {
	coord = x[i][2];
	if (periodicity[2]) {
	  if (coord < boxlo[2]) coord += prd[2];
	  if (coord >= boxhi[2]) coord -= prd[2];
	}
	ibinz = static_cast<int> ((coord - boxlo[2]) * invdelta[2]);
	ibinz = MAX(ibinz,0);
	ibinz = MIN(ibinz,nbinz-1);
      } else ibinz = 0;
      
      bin[i] = nbinx*nbiny*ibinz + nbinx*ibiny + ibinx;
    }

  if (triclinic) domain->lamda2x(nlocal);
}

/* ---------------------------------------------------------------------- */

double ComputeTempProfile::memory_usage()
{
  double bytes = maxbias * sizeof(double);
  bytes += maxatom * sizeof(int);
  bytes += nbins*(ncount+1) * sizeof(double);
  return bytes;
}
