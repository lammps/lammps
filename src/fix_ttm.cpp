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
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_ttm.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define MAXLINE 1024

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixTTM::FixTTM(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 15) error->all("Illegal fix TTM command");

  int seed = atoi(arg[3]);
  electronic_specific_heat = atof(arg[4]);
  electronic_density = atof(arg[5]);
  electronic_thermal_conductivity = atof(arg[6]);
  gamma_p = atof(arg[7]);
  gamma_s = atof(arg[8]);
  v_0 = atof(arg[9]);
  nxnodes = atoi(arg[10]);
  nynodes = atoi(arg[11]);
  nznodes = atoi(arg[12]);

  fpr = fopen(arg[13],"r");
  if (fpr == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",arg[13]);
    error->one(str);
  }

  nfileevery = atoi(arg[14]);

  if (nfileevery) {
    if (narg != 16) error->all("Illegal fix ttm command");
    MPI_Comm_rank(world,&me);
    if (me == 0) {
      fp = fopen(arg[15],"w");
      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open fix ttm file %s",arg[15]);
        error->one(str);
      }
    }
  }

  if (seed <= 0) error->all("Fix TTM seed must be >= 0");
  if (electronic_specific_heat <= 0.0) error->all("Fix TTM electronic_specific_heat must be > 0.0");
  if (electronic_density <= 0.0) error->all("Fix TTM electronic_density must be > 0.0");
  if (electronic_thermal_conductivity < 0.0) error->all("Fix TTM electronic_thermal_conductivity must be >= 0.0");
  if (gamma_p <= 0.0) error->all("Fix TTM gamma_p must be > 0.0");
  if (gamma_s < 0.0) error->all("Fix TTM gamma_s must be >= 0.0");
  if (v_0 < 0.0) error->all("Fix TTM v_0 must be >= 0.0");
  v_0_sq = v_0*v_0;
  if ((nxnodes <= 0) or (nynodes <= 0) or (nznodes <= 0)) error->all("Fix TTM nnodes must be > 0");

  nsum = nsum_prime = nsum_all = nsum_prime_all = T_initial_set = NULL;
  sum_vsq = sum_vsq_prime = sum_mass_vsq = sum_mass_vsq_prime = sum_vsq_all = sum_vsq_prime_all = NULL;
  sum_mass_vsq_all = sum_mass_vsq_prime_all = T_electron_old = T_electron = T_a = T_a_prime = g_s = g_p = NULL;

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  // allocate per-type arrays for force prefactors

  gfactor1 = new double[atom->ntypes+1];
  gfactor2 = new double[atom->ntypes+1];

  total_nnodes = nxnodes*nynodes*nznodes;

  // allocate memory for 3d vectors
  nsum = memory->create_3d_int_array(nxnodes,nynodes,nznodes,"TTM:nsum");
  nsum_prime = memory->create_3d_int_array(nxnodes,nynodes,nznodes,"TTM:nsum_prime");
  nsum_all = memory->create_3d_int_array(nxnodes,nynodes,nznodes,"TTM:nsum_all");
  nsum_prime_all = memory->create_3d_int_array(nxnodes,nynodes,nznodes,"TTM:nsum_prime_all");
  T_initial_set = memory->create_3d_int_array(nxnodes,nynodes,nznodes,"TTM:T_initial_set");
  sum_vsq = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:sum_vsq");
  sum_vsq_prime = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:sum_vsq_prime");
  sum_mass_vsq = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:sum_mass_vsq");
  sum_mass_vsq_prime = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:sum_mass_vsq_prime");
  sum_vsq_all = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:sum_vsq_all");
  sum_vsq_prime_all = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:sum_vsq_prime_all");
  sum_mass_vsq_all = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:sum_mass_vsq_all");
  sum_mass_vsq_prime_all = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:sum_mass_vsq_prime_all");
  T_electron_old = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:T_electron_old");
  T_electron = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:T_electron");
  T_a = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:T_a");
  T_a_prime = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:T_a_prime");
  g_s = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:g_s");
  g_p = memory->create_3d_double_array(nxnodes,nynodes,nznodes,"TTM:g_p");

  // set initial electron temperatures from user-supplied file

  if (me == 0) read_initial_electron_temperatures();
  MPI_Bcast(&T_electron[0][0][0],total_nnodes,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

FixTTM::~FixTTM()
{
  if (nfileevery && me == 0) fclose(fp);

  delete random;
  delete [] gfactor1;
  delete [] gfactor2;

  memory->destroy_3d_int_array(nsum);
  memory->destroy_3d_int_array(nsum_prime);
  memory->destroy_3d_int_array(nsum_all);
  memory->destroy_3d_int_array(nsum_prime_all);
  memory->destroy_3d_int_array(T_initial_set);
  memory->destroy_3d_double_array(sum_vsq);
  memory->destroy_3d_double_array(sum_vsq_prime); 
  memory->destroy_3d_double_array(sum_mass_vsq);
  memory->destroy_3d_double_array(sum_mass_vsq_prime);
  memory->destroy_3d_double_array(sum_vsq_all);
  memory->destroy_3d_double_array(sum_vsq_prime_all);
  memory->destroy_3d_double_array(sum_mass_vsq_all);
  memory->destroy_3d_double_array(sum_mass_vsq_prime_all);
  memory->destroy_3d_double_array(T_electron_old);
  memory->destroy_3d_double_array(T_electron);
  memory->destroy_3d_double_array(T_a);
  memory->destroy_3d_double_array(T_a_prime);
  memory->destroy_3d_double_array(g_s);
  memory->destroy_3d_double_array(g_p);
}

/* ---------------------------------------------------------------------- */

int FixTTM::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTTM::init()
{
  if (atom->mass == NULL)
    error->all("Cannot use fix TTM without per-type mass defined");

  // set force prefactors

  for (int i = 1; i <= atom->ntypes; i++) {
    gfactor1[i] = - gamma_p / force->ftm2v;
    gfactor2[i] = sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
  }

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ----------------------------------------------------------------------
   read in initial electron temperatures from a user-specified file
   only called by proc 0
------------------------------------------------------------------------- */

void FixTTM::read_initial_electron_temperatures()
{
  char line[MAXLINE];

  for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
    for (int iynode = 0; iynode < nynodes; iynode++) {
      for (int iznode = 0; iznode < nznodes; iznode++) {
        T_initial_set[ixnode][iynode][iznode] = 0;
      } 
    } 
  }

  // read initial electron temperature values from file

  int ixnode,iynode,iznode;
  double T_tmp;
  while (1) {
    if (fgets(line,MAXLINE,fpr) == NULL) break;
    sscanf(line,"%d %d %d %lg",&ixnode,&iynode,&iznode,&T_tmp);
    if (T_tmp < 0.0) error->all("Fix TTM electron temperatures must be > 0.0");
    T_electron[ixnode][iynode][iznode] = T_tmp;
    T_initial_set[ixnode][iynode][iznode] = 1;
  }

  for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
    for (int iynode = 0; iynode < nynodes; iynode++) {
      for (int iznode = 0; iznode < nznodes; iznode++) {
        if (T_initial_set[ixnode][iynode][iznode] == 0)
          error->all("All initial temperatures have not been set in the TTM fix.");
      } 
    } 
  }

  // close file

  fclose(fpr);
}

/* ---------------------------------------------------------------------- */

void FixTTM::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force(int vflag)
{
  update_electron_temperatures();

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double gamma1,gamma2;

  // apply damping and thermostat to all atoms in fix group

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      double xscale = (x[i][0] - domain->boxlo[0])/domain->xprd;
      double yscale = (x[i][1] - domain->boxlo[1])/domain->yprd;
      double zscale = (x[i][2] - domain->boxlo[2])/domain->zprd;
      int ixnode = static_cast<int>(xscale*nxnodes);
      int iynode = static_cast<int>(yscale*nynodes);
      int iznode = static_cast<int>(zscale*nznodes);
      while (ixnode > nxnodes-1) ixnode -= nxnodes;
      while (iynode > nynodes-1) iynode -= nynodes;
      while (iznode > nznodes-1) iznode -= nznodes;
      while (ixnode < 0) ixnode += nxnodes;
      while (iynode < 0) iynode += nynodes;
      while (iznode < 0) iznode += nznodes;

      double tsqrt = sqrt(T_electron[ixnode][iynode][iznode]);

      gamma1 = gfactor1[type[i]];
      double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      if (vsq > v_0_sq) gamma1 *= (gamma_p + gamma_s)/gamma_p;
      gamma2 = gfactor2[type[i]] * tsqrt;

      f[i][0] += gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
      f[i][1] += gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
      f[i][2] += gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTM::reset_dt()
{
  for (int i = 1; i <= atom->ntypes; i++) {
    gfactor2[i] = sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::update_electron_temperatures()
{
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // compute atomic Ta, and Ta' (for high v atoms) for each node

  for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
    for (int iynode = 0; iynode < nynodes; iynode++) {
      for (int iznode = 0; iznode < nznodes; iznode++) {
        nsum[ixnode][iynode][iznode] = 0;
        nsum_prime[ixnode][iynode][iznode] = 0;
        nsum_all[ixnode][iynode][iznode] = 0;
        nsum_prime_all[ixnode][iynode][iznode] = 0;
        sum_vsq[ixnode][iynode][iznode] = 0.0; 
        sum_vsq_prime[ixnode][iynode][iznode] = 0.0;
        sum_mass_vsq[ixnode][iynode][iznode] = 0.0;
        sum_mass_vsq_prime[ixnode][iynode][iznode] = 0.0;
        sum_vsq_all[ixnode][iynode][iznode] = 0.0;
        sum_vsq_prime_all[ixnode][iynode][iznode] = 0.0;
        sum_mass_vsq_all[ixnode][iynode][iznode] = 0.0;
        sum_mass_vsq_prime_all[ixnode][iynode][iznode] = 0.0;
      } 
    } 
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double mass = atom->mass[type[i]];
      double xscale = (x[i][0] - domain->boxlo[0])/domain->xprd;
      double yscale = (x[i][1] - domain->boxlo[1])/domain->yprd;
      double zscale = (x[i][2] - domain->boxlo[2])/domain->zprd;
      int ixnode = static_cast<int>(xscale*nxnodes);
      int iynode = static_cast<int>(yscale*nynodes);
      int iznode = static_cast<int>(zscale*nznodes);
      while (ixnode > nxnodes-1) ixnode -= nxnodes;
      while (iynode > nynodes-1) iynode -= nynodes;
      while (iznode > nznodes-1) iznode -= nznodes;
      while (ixnode < 0) ixnode += nxnodes;
      while (iynode < 0) iynode += nynodes;
      while (iznode < 0) iznode += nznodes;
      double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      nsum[ixnode][iynode][iznode] += 1;
      sum_vsq[ixnode][iynode][iznode] += vsq;
      sum_mass_vsq[ixnode][iynode][iznode] += mass*vsq;
      if (vsq > v_0_sq) {
        nsum_prime[ixnode][iynode][iznode] += 1;
        sum_vsq_prime[ixnode][iynode][iznode] += vsq;
        sum_mass_vsq_prime[ixnode][iynode][iznode] += mass*vsq;
      }
    }
  }

  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nznodes;
  double del_vol = dx*dy*dz;

  MPI_Allreduce(&nsum[0][0][0],&nsum_all[0][0][0],total_nnodes,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&nsum_prime[0][0][0],&nsum_prime_all[0][0][0],total_nnodes,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&sum_vsq[0][0][0],&sum_vsq_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&sum_vsq_prime[0][0][0],&sum_vsq_prime_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&sum_mass_vsq[0][0][0],&sum_mass_vsq_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&sum_mass_vsq_prime[0][0][0],&sum_mass_vsq_prime_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);

  double max_g_p = 0.0;
  for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
    for (int iynode = 0; iynode < nynodes; iynode++) {
      for (int iznode = 0; iznode < nznodes; iznode++) {
        if (nsum_all[ixnode][iynode][iznode] > 0) {
          T_a[ixnode][iynode][iznode] = sum_mass_vsq_all[ixnode][iynode][iznode]/
                                        (3*force->boltz*nsum_all[ixnode][iynode][iznode]/force->mvv2e);
          double g_p_tmp = gamma_p*sum_vsq_all[ixnode][iynode][iznode]/
                                        T_a[ixnode][iynode][iznode]/del_vol;
          max_g_p = MAX(max_g_p,g_p_tmp);
          g_p[ixnode][iynode][iznode] = g_p_tmp;
        } else { 
          T_a[ixnode][iynode][iznode] = 0;
          g_p[ixnode][iynode][iznode] = 0;
        }
        if (nsum_prime_all[ixnode][iynode][iznode] > 0) {
          T_a_prime[ixnode][iynode][iznode] = sum_mass_vsq_prime_all[ixnode][iynode][iznode]/
                                              (3*force->boltz*nsum_prime_all[ixnode][iynode][iznode]/force->mvv2e);
          g_s[ixnode][iynode][iznode] = gamma_s*sum_vsq_prime_all[ixnode][iynode][iznode]/
                                        T_a_prime[ixnode][iynode][iznode]/del_vol;
        } else {
          T_a_prime[ixnode][iynode][iznode] = 0;
          g_s[ixnode][iynode][iznode] = 0;
        }
      }
    }
  }

  // figure out how many inner steps (thermal solves) are required this MD timestep in order to maintain a stable explicit solve
  int num_inner_timesteps = 1;
  double inner_dt = update->dt;
  double stability_criterion = 1.0 - 2.0*inner_dt/(electronic_specific_heat*electronic_density)* 
   (electronic_thermal_conductivity*(1/dx/dx + 1/dy/dy + 1/dz/dz) + 0.5*max_g_p);
  if (stability_criterion < 0.0) {
    inner_dt = 0.5*(electronic_specific_heat*electronic_density)/(electronic_thermal_conductivity*(1/dx/dx + 1/dy/dy + 1/dz/dz) + 0.5*max_g_p);
    num_inner_timesteps = static_cast<int>(update->dt/inner_dt) + 1;
    inner_dt = update->dt/double(num_inner_timesteps);
    if (num_inner_timesteps > 1e6) error->warning("In the TTM fix, trying to do more than a million inner timesteps!");
  }

  for (int ith_inner_timestep = 0; ith_inner_timestep < num_inner_timesteps; ith_inner_timestep++) {
    
    for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
      for (int iynode = 0; iynode < nynodes; iynode++) {
        for (int iznode = 0; iznode < nznodes; iznode++) {
          T_electron_old[ixnode][iynode][iznode] = T_electron[ixnode][iynode][iznode];
        } 
      } 
    }
    
    // compute new electron T profile
    for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
      for (int iynode = 0; iynode < nynodes; iynode++) {
        for (int iznode = 0; iznode < nznodes; iznode++) {
          int right_xnode = ixnode + 1;
          int right_ynode = iynode + 1;
          int right_znode = iznode + 1;
          if (right_xnode == nxnodes) right_xnode = 0;
          if (right_ynode == nynodes) right_ynode = 0;
          if (right_znode == nznodes) right_znode = 0;
          int left_xnode = ixnode - 1;
          int left_ynode = iynode - 1;
          int left_znode = iznode - 1;
          if (left_xnode == -1) left_xnode = nxnodes - 1;
          if (left_ynode == -1) left_ynode = nynodes - 1;
          if (left_znode == -1) left_znode = nznodes - 1;
          T_electron[ixnode][iynode][iznode] = T_electron_old[ixnode][iynode][iznode] + 
           inner_dt/(electronic_specific_heat*electronic_density)*(electronic_thermal_conductivity*
           ((T_electron_old[right_xnode][iynode][iznode] + T_electron_old[left_xnode][iynode][iznode] - 2*T_electron_old[ixnode][iynode][iznode])/dx/dx +
            (T_electron_old[ixnode][right_ynode][iznode] + T_electron_old[ixnode][left_ynode][iznode] - 2*T_electron_old[ixnode][iynode][iznode])/dy/dy +
            (T_electron_old[ixnode][iynode][right_znode] + T_electron_old[ixnode][iynode][left_znode] - 2*T_electron_old[ixnode][iynode][iznode])/dz/dz) +
           g_p[ixnode][iynode][iznode]*(T_a[ixnode][iynode][iznode] - T_electron_old[ixnode][iynode][iznode]) + 
           g_s[ixnode][iynode][iznode]*T_a_prime[ixnode][iynode][iznode]);   
        }
      }
    }
  }

  if (!(update->ntimestep % nfileevery) && (me == 0)) { 
    fprintf(fp,"%d ",update->ntimestep);
    // print nodal temperatures for current timestep
    for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
      for (int iynode = 0; iynode < nynodes; iynode++) {
        for (int iznode = 0; iznode < nznodes; iznode++) {
          fprintf(fp,"%f ",T_a[ixnode][iynode][iznode]);
        }
      }
    }
    fprintf(fp,"\t");
    for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
      for (int iynode = 0; iynode < nynodes; iynode++) {
        for (int iznode = 0; iznode < nznodes; iznode++) {
          fprintf(fp,"%f ",T_electron[ixnode][iynode][iznode]);
        }
      }
    }
    fprintf(fp,"\n");
  }
}
