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
   Contributing author: Andres Jaramillo-Botero (Caltech))
------------------------------------------------------------------------- */

#include "mpi.h"
#include "compute_temp_deform_eff.h"
#include "update.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "modify.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NO_REMAP,X_REMAP,V_REMAP};                   // same as fix_deform.cpp

/* ---------------------------------------------------------------------- */

ComputeTempDeformEff::ComputeTempDeformEff(LAMMPS *lmp, int narg, char **arg) :
  ComputeTempDeform(lmp, narg, arg)
{
  // error check

  if (!atom->spin_flag || !atom->ervel_flag) 
    error->all("Compute temp/deform/eff requires atom attributes spin, ervel");
}

/* ---------------------------------------------------------------------- */

void ComputeTempDeformEff::dof_compute()
{
  double natoms = group->count(igroup);
  dof = domain->dimension * natoms;
  dof -= extra_dof + fix_dof;
  
  // just include nuclear dof

  int *spin = atom->spin;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int one = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (spin[i]) one++;
    }
  int nelectrons;
  MPI_Allreduce(&one,&nelectrons,1,MPI_INT,MPI_SUM,world);

  // the -3 recovers an extra_dof taken out because it's used by eradius

  dof -= domain->dimension * nelectrons - 3;
  
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
  double *rmass = atom->rmass;
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

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
      
      if (rmass) {
	t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + 
	      vthermal[2]*vthermal[2]) * rmass[i];
        if (spin[i]) t += 0.75*rmass[i]*ervel[i]*ervel[i];
      } else {
	t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + 
	      vthermal[2]*vthermal[2])* mass[type[i]];
        if (spin[i]) t += 0.75*mass[type[i]]*ervel[i]*ervel[i];
      }
    }
  
  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
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
  double *rmass = atom->rmass;
  int *spin = atom->spin; 
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
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
      
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      t[0] += massone * vthermal[0]*vthermal[0]; 
      t[1] += massone * vthermal[1]*vthermal[1];
      t[2] += massone * vthermal[2]*vthermal[2];
      t[3] += massone * vthermal[0]*vthermal[1];
      t[4] += massone * vthermal[0]*vthermal[2];
      t[5] += massone * vthermal[1]*vthermal[2];
      if (spin[i]) {
        t[0] += 0.75 * massone * ervel[i]*ervel[i];
        t[1] += 0.75 * massone * ervel[i]*ervel[i];
        t[2] += 0.75 * massone * ervel[i]*ervel[i];
      }
    }
  
  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (int i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}
