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
   Contributing author: Andres Jaramillo-Botero (Caltech)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include "compute_temp_deform_eff.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "group.h"
#include "comm.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempDeformEff::ComputeTempDeformEff(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute temp/deform/eff command");

  if (!atom->electron_flag)
    error->all(FLERR,"Compute temp/deform/eff requires atom style electron");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;
  tempbias = 1;

  maxbias = 0;
  vbiasall = NULL;
  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeTempDeformEff::~ComputeTempDeformEff()
{
  memory->destroy(vbiasall);
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeformEff::init()
{
  // check fix deform remap settings

  int i;
  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      if (((FixDeform *) modify->fix[i])->remapflag == Domain::X_REMAP &&
          comm->me == 0)
        error->warning(FLERR,"Using compute temp/deform/eff with inconsistent "
                       "fix deform remap option");
      break;
    }
  if (i == modify->nfix && comm->me == 0)
    error->warning(FLERR,
                   "Using compute temp/deform/eff with no fix deform defined");
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeformEff::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeformEff::dof_compute()
{
  adjust_dof_fix();
  natoms_temp = group->count(igroup);
  dof = domain->dimension * natoms_temp;
  dof -= extra_dof + fix_dof;

  // just include nuclear dof

  int *spin = atom->spin;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int one = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (abs(spin[i]) == 1) one++;
    }
  int nelectrons;
  MPI_Allreduce(&one,&nelectrons,1,MPI_INT,MPI_SUM,world);

  // Assume 3/2 k T per nucleus

  dof -= domain->dimension * nelectrons;

  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempDeformEff::compute_scalar()
{
  double lamda[3],vstream[3],vthermal[3];

  invoked_scalar = update->ntimestep;

  double **x = atom->x;
  double **v = atom->v;
  double *ervel = atom->ervel;
  double *mass = atom->mass;
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double mefactor = domain->dimension/4.0;

  // lamda = 0-1 triclinic lamda coords
  // vstream = streaming velocity = Hrate*lamda + Hratelo
  // vthermal = thermal velocity = v - vstream

  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  double t = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->x2lamda(x[i],lamda);
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] +
        h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      vthermal[0] = v[i][0] - vstream[0];
      vthermal[1] = v[i][1] - vstream[1];
      vthermal[2] = v[i][2] - vstream[2];

      if (mass) {
        t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] +
              vthermal[2]*vthermal[2])* mass[type[i]];
        if (abs(spin[i])==1) t += mefactor*mass[type[i]]*ervel[i]*ervel[i];
      }
    }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeformEff::compute_vector()
{
  double lamda[3],vstream[3],vthermal[3];

  invoked_vector = update->ntimestep;

  double **x = atom->x;
  double **v = atom->v;
  double *ervel = atom->ervel;
  double *mass = atom->mass;
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double mefactor = domain->dimension/4.0;

  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  double massone,t[6];
  for (int i = 0; i < 6; i++) t[i] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->x2lamda(x[i],lamda);
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] +
        h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      vthermal[0] = v[i][0] - vstream[0];
      vthermal[1] = v[i][1] - vstream[1];
      vthermal[2] = v[i][2] - vstream[2];

      massone = mass[type[i]];
      t[0] += massone * vthermal[0]*vthermal[0];
      t[1] += massone * vthermal[1]*vthermal[1];
      t[2] += massone * vthermal[2]*vthermal[2];
      t[3] += massone * vthermal[0]*vthermal[1];
      t[4] += massone * vthermal[0]*vthermal[2];
      t[5] += massone * vthermal[1]*vthermal[2];
      if (abs(spin[i])==1) {
        t[0] += mefactor * massone * ervel[i]*ervel[i];
        t[1] += mefactor * massone * ervel[i]*ervel[i];
        t[2] += mefactor * massone * ervel[i]*ervel[i];
      }
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (int i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
------------------------------------------------------------------------- */

void ComputeTempDeformEff::remove_bias(int i, double *v)
{
  double lamda[3];
  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  domain->x2lamda(atom->x[i],lamda);
  vbias[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] +
    h_rate[4]*lamda[2] + h_ratelo[0];
  vbias[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
  vbias[2] = h_rate[2]*lamda[2] + h_ratelo[2];
  v[0] -= vbias[0];
  v[1] -= vbias[1];
  v[2] -= vbias[2];
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
   NOTE: only removes translational velocity bias from electrons
------------------------------------------------------------------------- */

void ComputeTempDeformEff::remove_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (atom->nmax > maxbias) {
    memory->destroy(vbiasall);
    maxbias = atom->nmax;
    memory->create(vbiasall,maxbias,3,"temp/deform/eff:vbiasall");
  }

  double lamda[3];
  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->x2lamda(atom->x[i],lamda);
      vbiasall[i][0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] +
        h_rate[4]*lamda[2] + h_ratelo[0];
      vbiasall[i][1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vbiasall[i][2] = h_rate[2]*lamda[2] + h_ratelo[2];
      v[i][0] -= vbiasall[i][0];
      v[i][1] -= vbiasall[i][1];
      v[i][2] -= vbiasall[i][2];
    }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempDeformEff::restore_bias(int /*i*/, double *v)
{
  v[0] += vbias[0];
  v[1] += vbias[1];
  v[2] += vbias[2];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempDeformEff::restore_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] += vbiasall[i][0];
      v[i][1] += vbiasall[i][1];
      v[i][2] += vbiasall[i][2];
    }
}

/* ---------------------------------------------------------------------- */

double ComputeTempDeformEff::memory_usage()
{
  double bytes = maxbias * sizeof(double);
  return bytes;
}
