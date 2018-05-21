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

#include <mpi.h>
#include <cstdlib>
#include <cstring>
#include "compute_temp_drude.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
#include "lattice.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempDrude::ComputeTempDrude(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute temp command");

  vector_flag = 1;
  scalar_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = -1;
  extlist = new int[6];
  extlist[0] = extlist[1] = 0;
  extlist[2] = extlist[3] = extlist[4] = extlist[5] = 1;
  tempflag = 0; // because does not compute a single temperature (scalar and vector)

  vector = new double[6];
  fix_drude = NULL;
  id_temp = NULL;
  temperature = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeTempDrude::~ComputeTempDrude()
{
  delete [] vector;
  delete [] extlist;
  delete [] id_temp;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDrude::init()
{
  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(modify->fix[ifix]->style,"drude") == 0) break;
  if (ifix == modify->nfix) error->all(FLERR, "compute temp/drude requires fix drude");
  fix_drude = (FixDrude *) modify->fix[ifix];

  if (!comm->ghost_velocity)
    error->all(FLERR,"compute temp/drude requires ghost velocities. Use comm_modify vel yes");
}

/* ---------------------------------------------------------------------- */

void ComputeTempDrude::setup()
{
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempDrude::dof_compute()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int dim = domain->dimension;
  int *drudetype = fix_drude->drudetype;

  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);

  bigint dof_core_loc = 0, dof_drude_loc = 0;
  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      if (drudetype[type[i]] == DRUDE_TYPE) // Non-polarizable atom
          dof_drude_loc++;
      else
          dof_core_loc++;
    }
  }
  dof_core_loc *= dim;
  dof_drude_loc *= dim;
  MPI_Allreduce(&dof_core_loc,  &dof_core,  1, MPI_LMP_BIGINT, MPI_SUM, world);
  MPI_Allreduce(&dof_drude_loc, &dof_drude, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  dof_core -= fix_dof;
  vector[2] = dof_core;
  vector[3] = dof_drude;
}

/* ---------------------------------------------------------------------- */

int ComputeTempDrude::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,
                 "Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void ComputeTempDrude::compute_vector()
{
    invoked_vector = update->ntimestep;

    int nlocal = atom->nlocal;
    int *mask = atom->mask;
    int *type = atom->type;
    double *rmass = atom->rmass, *mass = atom->mass;
    double **v = atom->v;
    tagint *drudeid = fix_drude->drudeid;
    int *drudetype = fix_drude->drudetype;
    int dim = domain->dimension;
    double mvv2e = force->mvv2e, kb = force->boltz;

    double mcore, mdrude;
    double ecore, edrude;
    double *vcore, *vdrude;
    double kineng_core_loc = 0., kineng_drude_loc = 0.;
    for (int i=0; i<nlocal; i++){
        if (groupbit & mask[i] && drudetype[type[i]] != DRUDE_TYPE){
            if (drudetype[type[i]] == NOPOL_TYPE) {
                ecore = 0.;
                vcore = v[i];
                if (temperature) temperature->remove_bias(i, vcore);
                for (int k=0; k<dim; k++) ecore += vcore[k]*vcore[k];
                if (temperature) temperature->restore_bias(i, vcore);
                if (rmass) mcore = rmass[i];
                else mcore = mass[type[i]];
                kineng_core_loc += mcore * ecore;
            } else { // CORE_TYPE
                int j = atom->map(drudeid[i]);
                if (rmass) {
                    mcore = rmass[i];
                    mdrude = rmass[j];
                } else {
                    mcore = mass[type[i]];
                    mdrude = mass[type[j]];
                }
                double mtot_inv = 1. / (mcore + mdrude);
                ecore = 0.;
                edrude = 0.;
                vcore = v[i];
                vdrude = v[j];
                if (temperature) {
                    temperature->remove_bias(i, vcore);
                    temperature->remove_bias(j, vdrude);
                }
                for (int k=0; k<dim; k++) {
                    double v1 = mdrude * vdrude[k] + mcore * vcore[k];
                    ecore += v1 * v1;
                    double v2 = vdrude[k] - vcore[k];
                    edrude += v2 * v2;
                }
                if (temperature) {
                    temperature->restore_bias(i, vcore);
                    temperature->restore_bias(j, vdrude);
                }
                kineng_core_loc += mtot_inv * ecore;
                kineng_drude_loc += mtot_inv * mcore * mdrude * edrude;
            }
        }
    }

    if (dynamic) dof_compute();
    kineng_core_loc *= 0.5 * mvv2e;
    kineng_drude_loc *= 0.5 * mvv2e;
    MPI_Allreduce(&kineng_core_loc,&kineng_core,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&kineng_drude_loc,&kineng_drude,1,MPI_DOUBLE,MPI_SUM,world);
    temp_core = 2.0 * kineng_core / (dof_core * kb);
    temp_drude = 2.0 * kineng_drude / (dof_drude * kb);
    vector[0] = temp_core;
    vector[1] = temp_drude;
    vector[4] = kineng_core;
    vector[5] = kineng_drude;
}

double ComputeTempDrude::compute_scalar(){
    compute_vector();
    scalar = vector[0];
    return scalar;
}

