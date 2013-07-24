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

enum{TENSOR,BIN};

/* ---------------------------------------------------------------------- */

ComputeTempProfile::ComputeTempProfile(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal compute temp/profile command");

  scalar_flag = 1;
  extscalar = 0;
  tempflag = 1;
  tempbias = 1;

  xflag = atoi(arg[3]);
  yflag = atoi(arg[4]);
  zflag = atoi(arg[5]);
  if (zflag && domain->dimension == 2)
    error->all(FLERR,"Compute temp/profile cannot use vz for 2d systemx");

  ncount = 0;
  ivx = ivy = ivz = 0;
  if (xflag) ivx = ncount++;
  if (yflag) ivy = ncount++;
  if (zflag) ivz = ncount++;
  ncount += 2;

  nbinx = nbiny = nbinz = 1;
  int lastarg;

  int iarg = 6;
  if (strcmp(arg[iarg],"x") == 0) {
    if (iarg+2 > narg) error->all(FLERR,"Illegal compute temp/profile command");
    nbinx = atoi(arg[iarg+1]);
    iarg += 2;
  } else if (strcmp(arg[iarg],"y") == 0) {
    if (iarg+2 > narg) error->all(FLERR,"Illegal compute temp/profile command");
    nbiny = atoi(arg[iarg+1]);
    iarg += 2;
  } else if (strcmp(arg[iarg],"z") == 0) {
    if (iarg+2 > narg) error->all(FLERR,"Illegal compute temp/profile command");
    if (domain->dimension == 2)
      error->all(FLERR,"Compute temp/profile cannot bin z for 2d systems");
    nbinz = atoi(arg[iarg+1]);
    iarg += 2;
  } else if (strcmp(arg[iarg],"xy") == 0) {
    if (iarg+3 > narg) error->all(FLERR,"Illegal compute temp/profile command");
    nbinx = atoi(arg[iarg+1]);
    nbiny = atoi(arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"yz") == 0) {
    if (iarg+3 > narg) error->all(FLERR,"Illegal compute temp/profile command");
    if (domain->dimension == 2)
      error->all(FLERR,"Compute temp/profile cannot bin z for 2d systems");
    nbiny = atoi(arg[iarg+1]);
    nbinz = atoi(arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"xz") == 0) {
    if (iarg+3 > narg) error->all(FLERR,"Illegal compute temp/profile command");
    if (domain->dimension == 2)
      error->all(FLERR,"Compute temp/profile cannot bin z for 2d systems");
    nbinx = atoi(arg[iarg+1]);
    nbinz = atoi(arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"xyz") == 0) {
    if (iarg+4 > narg) error->all(FLERR,"Illegal compute temp/profile command");
    if (domain->dimension == 2)
      error->all(FLERR,"Compute temp/profile cannot bin z for 2d systems");
    nbinx = atoi(arg[iarg+1]);
    nbiny = atoi(arg[iarg+2]);
    nbinz = atoi(arg[iarg+3]);
    iarg += 4;
  } else error->all(FLERR,"Illegal compute temp/profile command");

  // optional keywords

  outflag = TENSOR;
  
  while (iarg < narg) {
    if (strcmp(arg[iarg],"out") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute temp/profile command");
      if (strcmp(arg[iarg+1],"tensor") == 0) outflag = TENSOR;
      else if (strcmp(arg[iarg+1],"bin") == 0) outflag = BIN;
      else error->all(FLERR,"Illegal compute temp/profile command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute temp/profile command");
  }

  // setup

  nbins = nbinx*nbiny*nbinz;
  if (nbins <= 0) error->all(FLERR,"Illegal compute temp/profile command");

  memory->create(vbin,nbins,ncount,"temp/profile:vbin");
  memory->create(binave,nbins,ncount,"temp/profile:binave");

  if (outflag == TENSOR) {
    vector_flag = 1;
    size_vector = 6;
    extvector = 1;
    vector = new double[size_vector];
  } else {
    array_flag = 1;
    size_array_rows = nbins;
    size_array_cols = 2;
    extarray = 0;
    memory->create(tbin,nbins,"temp/profile:tbin");
    memory->create(tbinall,nbins,"temp/profile:tbinall");
    memory->create(array,nbins,2,"temp/profile:array");
  }

  maxatom = 0;
  bin = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeTempProfile::~ComputeTempProfile()
{
  memory->destroy(vbin);
  memory->destroy(binave);
  memory->destroy(bin);
  if (outflag == TENSOR) delete [] vector;
  else {
    memory->destroy(tbin);
    memory->destroy(tbinall);
    memory->destroy(array);
  }
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

void ComputeTempProfile::setup()
{
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
  dof_compute();
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

/* ---------------------------------------------------------------------- */

void ComputeTempProfile::compute_array()
{
  int i,ibin;
  double vthermal[3];

  invoked_array = update->ntimestep;

  bin_average();

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < nbins; i++) tbin[i] = 0.0;

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
        tbin[ibin] += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] +
                       vthermal[2]*vthermal[2]) * rmass[i];
      else
        tbin[ibin] += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] +
                       vthermal[2]*vthermal[2]) * mass[type[i]];
    }

  MPI_Allreduce(tbin,tbinall,nbins,MPI_DOUBLE,MPI_SUM,world);

  int nper = domain->dimension;
  for (i = 0; i < nbins; i++) {
    array[i][0] = binave[i][ncount-1];
    if (array[i][0] > 0.0) {
      dof = nper*array[i][0] - extra_dof;
      if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
      else tfactor = 0.0;
      array[i][1] = tfactor*tbinall[i];
    } else array[i][1] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempProfile::remove_bias(int i, double *v)
{
  int ibin = bin[i];
  if (xflag) v[0] -= binave[ibin][ivx];
  if (yflag) v[1] -= binave[ibin][ivy];
  if (zflag) v[2] -= binave[ibin][ivz];
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempProfile::remove_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int ibin;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ibin = bin[i];
      if (xflag) v[i][0] -= binave[ibin][ivx];
      if (yflag) v[i][1] -= binave[ibin][ivy];
      if (zflag) v[i][2] -= binave[ibin][ivz];
    }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempProfile::restore_bias(int i, double *v)
{
  int ibin = bin[i];
  if (xflag) v[0] += binave[ibin][ivx];
  if (yflag) v[1] += binave[ibin][ivy];
  if (zflag) v[2] += binave[ibin][ivz];
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

  int ibin;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ibin = bin[i];
      if (xflag) v[i][0] += binave[ibin][ivx];
      if (yflag) v[i][1] += binave[ibin][ivy];
      if (zflag) v[i][2] += binave[ibin][ivz];
    }
}

/* ----------------------------------------------------------------------
   compute average COM velocity in each bin
------------------------------------------------------------------------- */

void ComputeTempProfile::bin_average()
{
  int i,j,ibin;

  if (box_change) bin_setup();
  bin_assign();

  // clear bins, including particle mass and count

  for (i = 0; i < nbins; i++)
    for (j = 0; j < ncount; j++)
      vbin[i][j] = 0.0;

  // sum each particle's mass-weighted velocity, mass, count to appropriate bin

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  int nc2 = ncount-2;
  int nc1 = ncount-1;

  if (rmass) {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        ibin = bin[i];
        if (xflag) vbin[ibin][ivx] += rmass[i]*v[i][0];
        if (yflag) vbin[ibin][ivy] += rmass[i]*v[i][1];
        if (zflag) vbin[ibin][ivz] += rmass[i]*v[i][2];
        vbin[ibin][nc2] += rmass[i];
        vbin[ibin][nc1] += 1.0;
      }
  } else {
    double onemass;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        ibin = bin[i];
        onemass = mass[type[i]];
        if (xflag) vbin[ibin][ivx] += onemass*v[i][0];
        if (yflag) vbin[ibin][ivy] += onemass*v[i][1];
        if (zflag) vbin[ibin][ivz] += onemass*v[i][2];
        vbin[ibin][nc2] += onemass;
        vbin[ibin][nc1] += 1.0;
      }
  }

  // sum bins across processors

  MPI_Allreduce(vbin[0],binave[0],nbins*ncount,MPI_DOUBLE,MPI_SUM,world);

  // compute ave COM velocity in each bin, checking for no particles

  for (i = 0; i < nbins; i++)
    if (binave[i][nc1] > 0.0)
      for (j = 0; j < nc2; j++)
        binave[i][j] /= binave[i][nc2];
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
    memory->destroy(bin);
    memory->create(bin,maxatom,"temp/profile:bin");
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
  double bytes = maxatom * sizeof(int);
  bytes += nbins*ncount * sizeof(double);
  return bytes;
}
