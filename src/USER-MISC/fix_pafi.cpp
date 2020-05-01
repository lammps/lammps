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
/* ------------------------------------------------------------------------
   Contributing authors: Thomas Swinburne (CNRS & CINaM, Marseille, France)

   Please cite the related publication:
   T.D. Swinburne and M.-C. Marinica, Unsupervised calculation of free energy barriers in large crystalline systems, Physical Review Letters 2018
------------------------------------------------------------------------- */


#include "fix_pafi.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include "math_extra.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "compute.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "utils.h"
#include "citeme.h"

using namespace LAMMPS_NS;

static const char cite_fix_pafi_package[] =
  "citation for fix pafi:\n\n"
  "@article{SwinburneMarinica2018,\n"
  "author={T. D. Swinburne and M. C. Marinica},\n"
  "title={Unsupervised calculation of free energy barriers in large "
  "crystalline systems},\n"
  "journal={Physical Review Letters},\n"
  "volume={276},\n"
  "number={1},\n"
  "pages={154-165},\n"
  "year={2018},\n"
  "publisher={APS}\n}\n"
  "Recommended to be coupled with PAFI++ code:\n"
  "https://github.com/tomswinburne/pafi\n";

using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPAFI::FixPAFI(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), idregion(NULL), random(NULL), computename(NULL),
      h(NULL), step_respa(NULL)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_pafi_package);

  if (narg < 11) error->all(FLERR,"Illegal fix pafi command");

  dynamic_group_allow = 0;
  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extvector = 0;
  od_flag = 0;
  com_flag = 1;
  time_integrate = 1;

  int n = strlen(arg[3])+1;
  computename = new char[n];
  strcpy(computename,&arg[3][0]);


  icompute = modify->find_compute(computename);
  if (icompute == -1)
    error->all(FLERR,"Compute ID for fix pafi does not exist");
  PathCompute = modify->compute[icompute];
  if (PathCompute->peratom_flag==0)
    error->all(FLERR,"Compute for fix pafi does not calculate a local array");
  if (PathCompute->size_peratom_cols < 9)
    error->all(FLERR,"Compute for fix pafi must have 9 fields per atom");
");

  if (comm->me==0) {
    if (screen) fprintf(screen,
      "fix pafi compute name,style: %s,%s\n",computename,PathCompute->style);
    if (logfile) fprintf(logfile,
      "fix pafi compute name,style: %s,%s\n",computename,PathCompute->style);
  }


  respa_level_support = 1;
  ilevel_respa = nlevels_respa = 0;

  temperature = force->numeric(FLERR,arg[4]);
  t_period = force->numeric(FLERR,arg[5]);
  seed = force->inumeric(FLERR,arg[6]);
  // TODO UNITS
  gamma = 1. / t_period / force->ftm2v;
  sqrtD = sqrt(1.) * sqrt(24.0*force->boltz/t_period/update->dt/force->mvv2e*temperature) / force->ftm2v;

  // optional args
  iregion = -1;
  idregion = NULL; // not used
  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"overdamped") == 0) {
      if (strcmp(arg[iarg+1],"no") == 0 || strcmp(arg[iarg+1],"0") == 0)
        od_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0 || strcmp(arg[iarg+1],"1") == 0)
        od_flag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"com") == 0) {
      if (strcmp(arg[iarg+1],"no") == 0 || strcmp(arg[iarg+1],"0")
        com_flag = 0;
      else if (strcmp(arg[iarg+1],"yes"  || strcmp(arg[iarg+1],"1") == 0)
        com_flag = 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal fix pafi command");
  }
  force_flag = 0;

  for(int i = 0; i < 10; i++) {
    c_v[i] = 0.;
    c_v_all[i] = 0.;
  }
  for(int i=0; i<5; i++) {
    proj[i] = 0.0;
    proj_all[i] = 0.0;
  }
  for(int i=0; i<4; i++) {
    results[i] = 0.0;
    results_all[i] = 0.0;
  }
  maxatom = 1;
  memory->create(h,maxatom,3,"FixPAFI:h");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);

}

/* ---------------------------------------------------------------------- */

FixPAFI::~FixPAFI()
{
  if (copymode) return;
  delete random;
  delete [] idregion;
  delete [] computename;
  memory->destroy(h);
}

/* ---------------------------------------------------------------------- */

int FixPAFI::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  // nve
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPAFI::init()
{
  // set index and check validity of region
  // nve
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix pafi does not exist");
  }


  icompute = modify->find_compute(computename);
  if (icompute == -1)
    error->all(FLERR,"Compute ID for fix pafi does not exist");
  PathCompute = modify->compute[icompute];
  if (PathCompute->peratom_flag==0)
    error->all(FLERR,"Compute for fix pafi does not calculate a local array");
  if (PathCompute->size_peratom_cols < 9)
    error->all(FLERR,"Compute for fix pafi must have 9 fields per atom");



  if (strstr(update->integrate_style,"respa")) {
    step_respa = ((Respa *) update->integrate)->step; // nve
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,nlevels_respa-1);
    else ilevel_respa = nlevels_respa-1;
  }

}

void FixPAFI::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
}

void FixPAFI::min_setup(int vflag)
{
  if( utils::strmatch(update->minimize_style,"^fire")==0 &&
          utils::strmatch(update->minimize_style,"^quickmin")==0 )
    error->all(FLERR,"fix pafi requires a damped dynamics minimizer");
  min_post_force(vflag);
}


void FixPAFI::post_force(int vflag)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }
  // reallocate norm array if necessary
  if (atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(h);
    memory->create(h,maxatom,3,"FixPAFI:h");
  }

  PathCompute->compute_peratom();
  double **path = PathCompute->array_atom;

  double xum=0.;

  // proj 0,1,2 = f.n, v.n, h.n
  // proj 3,4,5 = psi, f.n**2, f*(1-psi)
  // c_v 0,1,2 = fxcom, fycom, fzcom etc
  for(int i = 0; i < 10; i++) {
    c_v[i] = 0.;
    c_v_all[i] = 0.;
  }
  for(int i = 0; i < 5; i++) {
    proj[i] = 0.;
    proj_all[i] = 0.;
  }

  double deviation[3] = {0.,0.,0.};

  double fn;

  force_flag=0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;

      h[i][0] = random->uniform() - 0.5;
      h[i][1] = random->uniform() - 0.5;
      h[i][2] = random->uniform() - 0.5;

      proj[0] += f[i][0] * path[i][3]; // f.n
      proj[0] += f[i][1] * path[i][4]; // f.n
      proj[0] += f[i][2] * path[i][5]; // f.n

      proj[1] += v[i][0] * path[i][3]; // v.n
      proj[1] += v[i][1] * path[i][4]; // v.n
      proj[1] += v[i][2] * path[i][5]; // v.n

      proj[2] += h[i][0] * path[i][3]; // h.n
      proj[2] += h[i][1] * path[i][4]; // h.n
      proj[2] += h[i][2] * path[i][5]; // h.n

      deviation[0] = x[i][0]-path[i][0]; // x-path
      deviation[1] = x[i][1]-path[i][1]; // x-path
      deviation[2] = x[i][2]-path[i][2]; // x-path
      domain->minimum_image(deviation);

      proj[3] += path[i][6]*deviation[0]; // (x-path).dn/nn = psi
      proj[3] += path[i][7]*deviation[1]; // (x-path).dn/nn = psi
      proj[3] += path[i][8]*deviation[2]; // (x-path).dn/nn = psi

      proj[4] += path[i][3]*deviation[0]; // (x-path).n
      proj[4] += path[i][4]*deviation[1]; // (x-path).n
      proj[4] += path[i][5]*deviation[2]; // (x-path).n

    }
  }

  if(com_flag == 0){
    c_v[9] += 1.0;
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;

        c_v[0] += f[i][0];
        c_v[1] += f[i][1];
        c_v[2] += f[i][2];

        c_v[3] += v[i][0];
        c_v[4] += v[i][1];
        c_v[5] += v[i][2];

        c_v[6] += h[i][0];
        c_v[7] += h[i][1];
        c_v[8] += h[i][2];

        c_v[9] += 1.0;
      }
  }
  MPI_Allreduce(proj,proj_all,5,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(c_v,c_v_all,10,MPI_DOUBLE,MPI_SUM,world);


  // results - f.n*(1-psi), (f.n)^2*(1-psi)^2, 1-psi, dX.n
  if(comm->me ==0) {
    results_all[0] = proj_all[0] * (1.-proj_all[3]);
    results_all[1] = results_all[0] * results_all[0];
    results_all[2] = 1.-proj_all[3];
    results_all[3] = proj_all[4];
  }
  MPI_Bcast(results_all,4,MPI_DOUBLE,0,world);
  force_flag = 1;

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;

      f[i][0] -= proj_all[0] * path[i][3] + c_v_all[0]/c_v_all[9];
      f[i][1] -= proj_all[0] * path[i][4] + c_v_all[1]/c_v_all[9];
      f[i][2] -= proj_all[0] * path[i][5] + c_v_all[2]/c_v_all[9];

      v[i][0] -= proj_all[1] * path[i][3] + c_v_all[3]/c_v_all[9];
      v[i][1] -= proj_all[1] * path[i][4] + c_v_all[4]/c_v_all[9];
      v[i][2] -= proj_all[1] * path[i][5] + c_v_all[5]/c_v_all[9];

      h[i][0] -= proj_all[2] * path[i][3] + c_v_all[6]/c_v_all[9];
      h[i][1] -= proj_all[2] * path[i][4] + c_v_all[7]/c_v_all[9];
      h[i][2] -= proj_all[2] * path[i][5] + c_v_all[8]/c_v_all[9];
    }
  }

  if (od_flag == 0) {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        if(rmass) mass_f = sqrt(rmass[i]);
        else mass_f = sqrt(mass[type[i]]);

        f[i][0] += -gamma * mass_f * mass_f * v[i][0];
        f[i][1] += -gamma * mass_f * mass_f * v[i][1];
        f[i][2] += -gamma * mass_f * mass_f * v[i][2];

        f[i][0] += sqrtD * mass_f * h[i][0];
        f[i][1] += sqrtD * mass_f * h[i][1];
        f[i][2] += sqrtD * mass_f * h[i][2];
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;

        if(rmass) mass_f = sqrt(rmass[i]);
        else mass_f = sqrt(mass[type[i]]);

        f[i][0] += sqrtD * h[i][0] * mass_f;
        f[i][1] += sqrtD * h[i][1] * mass_f;
        f[i][2] += sqrtD * h[i][2] * mass_f;

        f[i][0] /=  gamma * mass_f * mass_f;
        f[i][1] /=  gamma * mass_f * mass_f;
        f[i][2] /=  gamma * mass_f * mass_f;

      }
    }
  }

};

void FixPAFI::post_force_respa(int vflag, int ilevel, int iloop)
{
  // set force to desired value on requested level, 0.0 on other levels

  if (ilevel == ilevel_respa) post_force(vflag);
  else {
    Region *region = NULL;
    if (iregion >= 0) {
      region = domain->regions[iregion];
      region->prematch();
    }
    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        for (int k = 0; k < 3; k++) f[i][k] = 0.0;
      }
  }
};

void FixPAFI::min_post_force(int vflag)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  PathCompute->compute_peratom();
  double **path = PathCompute->array_atom;

  double xum=0.;

  // proj 0,1,2 = f.n, v.n, h.n
  // proj 3,4,5 = psi, f.n**2, f*(1-psi)
  // c_v 0,1,2 = fxcom, fycom, fzcom etc
  for(int i = 0; i < 10; i++) {
    c_v[i] = 0.;
    c_v_all[i] = 0.;
  }
  for(int i = 0; i < 5; i++) {
    proj[i] = 0.;
    proj_all[i] = 0.;
  }

  double deviation[3] = {0.,0.,0.};

  double fn;

  force_flag=0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;

      proj[0] += f[i][0] * path[i][3]; // f.n
      proj[0] += f[i][1] * path[i][4]; // f.n
      proj[0] += f[i][2] * path[i][5]; // f.n

      proj[1] += v[i][0] * path[i][3]; // v.n
      proj[1] += v[i][1] * path[i][4]; // v.n
      proj[1] += v[i][2] * path[i][5]; // v.n

      proj[2] += h[i][0] * path[i][3]; // h.n
      proj[2] += h[i][1] * path[i][4]; // h.n
      proj[2] += h[i][2] * path[i][5]; // h.n

      deviation[0] = x[i][0]-path[i][0]; // x-path
      deviation[1] = x[i][1]-path[i][1]; // x-path
      deviation[2] = x[i][2]-path[i][2]; // x-path
      domain->minimum_image(deviation);

      proj[3] += path[i][6]*deviation[0]; // (x-path).dn/nn = psi
      proj[3] += path[i][7]*deviation[1]; // (x-path).dn/nn = psi
      proj[3] += path[i][8]*deviation[2]; // (x-path).dn/nn = psi

      proj[4] += path[i][3]*deviation[0]; // (x-path).n
      proj[4] += path[i][4]*deviation[1]; // (x-path).n
      proj[4] += path[i][5]*deviation[2]; // (x-path).n

    }
  }

  if(com_flag == 0){
    c_v[9] += 1.0;
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;

        c_v[0] += f[i][0];
        c_v[1] += f[i][1];
        c_v[2] += f[i][2];

        c_v[3] += v[i][0];
        c_v[4] += v[i][1];
        c_v[5] += v[i][2];

        c_v[6] += h[i][0];
        c_v[7] += h[i][1];
        c_v[8] += h[i][2];

        c_v[9] += 1.0;
      }
  }
  MPI_Allreduce(proj,proj_all,5,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(c_v,c_v_all,10,MPI_DOUBLE,MPI_SUM,world);

  results_all[0] = proj_all[0] * (1.-proj_all[3]); // f.n * psi
  results_all[1] = results_all[0] * results_all[0]; // (f.n * psi)^2
  results_all[2] = 1.-proj_all[3]; // psi
  results_all[3] = proj_all[4]; // dX.n

  MPI_Bcast(results_all,4,MPI_DOUBLE,0,world);
  force_flag = 1;

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;

      f[i][0] -= proj_all[0] * path[i][3] + c_v_all[0]/c_v_all[9];
      f[i][1] -= proj_all[0] * path[i][4] + c_v_all[1]/c_v_all[9];
      f[i][2] -= proj_all[0] * path[i][5] + c_v_all[2]/c_v_all[9];

      v[i][0] -= proj_all[1] * path[i][3] + c_v_all[3]/c_v_all[9];
      v[i][1] -= proj_all[1] * path[i][4] + c_v_all[4]/c_v_all[9];
      v[i][2] -= proj_all[1] * path[i][5] + c_v_all[5]/c_v_all[9];

    }
  }
};

double FixPAFI::compute_vector(int n)
{
  return results_all[n];
};



void FixPAFI::initial_integrate(int vflag)
{
  double dtfm;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;

  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  PathCompute->compute_peratom();
  double **path = PathCompute->array_atom;

  for(int i = 0; i < 10; i++) {
    c_v[i] = 0.;
    c_v_all[i] = 0.;
  }
  for(int i = 0; i < 5; i++) {
    proj[i] = 0.;
    proj_all[i] = 0.;
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      proj[0] += f[i][0] * path[i][3]; // f.n
      proj[0] += f[i][1] * path[i][4]; // f.n
      proj[0] += f[i][2] * path[i][5]; // f.n

      proj[1] += v[i][0] * path[i][3]; // v.n
      proj[1] += v[i][1] * path[i][4]; // v.n
      proj[1] += v[i][2] * path[i][5]; // v.n
    }
  }
  if(com_flag == 0){
    c_v[9] += 1.0;
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {

        c_v[0] += v[i][0];
        c_v[1] += v[i][1];
        c_v[2] += v[i][2];

        c_v[3] += f[i][0];
        c_v[4] += f[i][1];
        c_v[5] += f[i][2];

        c_v[9] += 1.0;
      }
  }


  MPI_Allreduce(proj,proj_all,5,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(c_v,c_v_all,10,MPI_DOUBLE,MPI_SUM,world);

  if (od_flag == 0){
    if (rmass) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * (f[i][0]-path[i][3]*proj_all[0] - c_v_all[3]/c_v_all[9]);
          v[i][1] += dtfm * (f[i][1]-path[i][4]*proj_all[0] - c_v_all[4]/c_v_all[9]);
          v[i][2] += dtfm * (f[i][2]-path[i][5]*proj_all[0] - c_v_all[5]/c_v_all[9]);
          x[i][0] += dtv * (v[i][0]-path[i][3]*proj_all[1] - c_v_all[0]/c_v_all[9]);
          x[i][1] += dtv * (v[i][1]-path[i][4]*proj_all[1] - c_v_all[1]/c_v_all[9]);
          x[i][2] += dtv * (v[i][2]-path[i][5]*proj_all[1] - c_v_all[2]/c_v_all[9]);
        }
    } else {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * (f[i][0]-path[i][3]*proj_all[0] - c_v_all[3]/c_v_all[9]);
          v[i][1] += dtfm * (f[i][1]-path[i][4]*proj_all[0] - c_v_all[4]/c_v_all[9]);
          v[i][2] += dtfm * (f[i][2]-path[i][5]*proj_all[0] - c_v_all[5]/c_v_all[9]);
          x[i][0] += dtv * (v[i][0]-path[i][3]*proj_all[1] - c_v_all[0]/c_v_all[9]);
          x[i][1] += dtv * (v[i][1]-path[i][4]*proj_all[1] - c_v_all[1]/c_v_all[9]);
          x[i][2] += dtv * (v[i][2]-path[i][5]*proj_all[1] - c_v_all[2]/c_v_all[9]);
        }
    }
  } else {
    if (rmass) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          dtfm = dtf / rmass[i];
          v[i][0] = 0.;
          v[i][1] = 0.;
          v[i][2] = 0.;
          x[i][0] += dtv * (f[i][0]-path[i][3]*proj_all[0] - c_v_all[3]/c_v_all[9]);
          x[i][1] += dtv * (f[i][1]-path[i][4]*proj_all[0] - c_v_all[4]/c_v_all[9]);
          x[i][2] += dtv * (f[i][2]-path[i][5]*proj_all[0] - c_v_all[5]/c_v_all[9]);
        }
    } else {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          dtfm = dtf / mass[type[i]];
          v[i][0] = 0.;
          v[i][1] = 0.;
          v[i][2] = 0.;
          x[i][0] += dtv * (f[i][0]-path[i][3]*proj_all[0] - c_v_all[3]/c_v_all[9]);
          x[i][1] += dtv * (f[i][1]-path[i][4]*proj_all[0] - c_v_all[4]/c_v_all[9]);
          x[i][2] += dtv * (f[i][2]-path[i][5]*proj_all[0] - c_v_all[5]/c_v_all[9]);
        }
    }
  }
};

/* ---------------------------------------------------------------------- */

void FixPAFI::final_integrate()
{
  double dtfm;

  // update v of atoms in group
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  PathCompute->compute_peratom();
  double **path = PathCompute->array_atom;

  for(int i = 0; i < 10; i++) {
    c_v[i] = 0.;
    c_v_all[i] = 0.;
  }
  for(int i = 0; i < 5; i++) {
    proj[i] = 0.;
    proj_all[i] = 0.;
  }
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      proj[0] += f[i][0] * path[i][3]; // f.n
      proj[0] += f[i][1] * path[i][4]; // f.n
      proj[0] += f[i][2] * path[i][5]; // f.n
    }
  if(com_flag == 0){
    c_v[9] += 1.0;
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        c_v[3] += f[i][0];
        c_v[4] += f[i][1];
        c_v[5] += f[i][2];
        c_v[9] += 1.0;
      }
  }

  MPI_Allreduce(proj,proj_all,5,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(c_v,c_v_all,10,MPI_DOUBLE,MPI_SUM,world);

  if (od_flag == 0){
    if (rmass) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * (f[i][0]-path[i][3]*proj_all[0] - c_v_all[3]/c_v_all[9]);
          v[i][1] += dtfm * (f[i][1]-path[i][4]*proj_all[0] - c_v_all[4]/c_v_all[9]);
          v[i][2] += dtfm * (f[i][2]-path[i][5]*proj_all[0] - c_v_all[5]/c_v_all[9]);
        }
    } else {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * (f[i][0]-path[i][3]*proj_all[0] - c_v_all[3]/c_v_all[9]);
          v[i][1] += dtfm * (f[i][1]-path[i][4]*proj_all[0] - c_v_all[4]/c_v_all[9]);
          v[i][2] += dtfm * (f[i][2]-path[i][5]*proj_all[0] - c_v_all[5]/c_v_all[9]);
        }
    }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        v[i][0] = 0.;
        v[i][1] = 0.;
        v[i][2] = 0.;
      }
  }
};

/* ---------------------------------------------------------------------- */

void FixPAFI::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
};

/* ---------------------------------------------------------------------- */

void FixPAFI::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
};

/* ---------------------------------------------------------------------- */

void FixPAFI::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
};


/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixPAFI::memory_usage()
{
  double bytes = 0.0;
  bytes = maxatom*3 * sizeof(double);
  return bytes;
};
