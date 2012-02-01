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
   Contributing authors: Paul Crozier (SNL)
                         Carolyn Phillips (University of Michigan)
------------------------------------------------------------------------- */

#include "lmptype.h"
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

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

FixTTM::FixTTM(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 15) error->all(FLERR,"Illegal fix ttm command");

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 1;
  nevery = 1;
  restart_peratom = 1;
  restart_global = 1;

  seed = atoi(arg[3]);
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
    error->one(FLERR,str);
  }

  nfileevery = atoi(arg[14]);

  if (nfileevery) {
    if (narg != 16) error->all(FLERR,"Illegal fix ttm command");
    MPI_Comm_rank(world,&me);
    if (me == 0) {
      fp = fopen(arg[15],"w");
      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open fix ttm file %s",arg[15]);
        error->one(FLERR,str);
      }
    }
  }

  // error check

  if (seed <= 0) error->all(FLERR,"Invalid random number seed in fix ttm command");
  if (electronic_specific_heat <= 0.0) 
    error->all(FLERR,"Fix ttm electronic_specific_heat must be > 0.0");
  if (electronic_density <= 0.0) 
    error->all(FLERR,"Fix ttm electronic_density must be > 0.0");
  if (electronic_thermal_conductivity < 0.0)
    error->all(FLERR,"Fix ttm electronic_thermal_conductivity must be >= 0.0");
  if (gamma_p <= 0.0) error->all(FLERR,"Fix ttm gamma_p must be > 0.0");
  if (gamma_s < 0.0) error->all(FLERR,"Fix ttm gamma_s must be >= 0.0");
  if (v_0 < 0.0) error->all(FLERR,"Fix ttm v_0 must be >= 0.0");
  if (nxnodes <= 0 || nynodes <= 0 || nznodes <= 0)
    error->all(FLERR,"Fix ttm number of nodes must be > 0");

  v_0_sq = v_0*v_0;

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  // allocate per-type arrays for force prefactors

  gfactor1 = new double[atom->ntypes+1];
  gfactor2 = new double[atom->ntypes+1];

  // allocate 3d grid variables

  total_nnodes = nxnodes*nynodes*nznodes;

  memory->create(nsum,nxnodes,nynodes,nznodes,"ttm:nsum");
  memory->create(nsum_all,nxnodes,nynodes,nznodes,"ttm:nsum_all");
  memory->create(T_initial_set,nxnodes,nynodes,nznodes,"ttm:T_initial_set");
  memory->create(sum_vsq,nxnodes,nynodes,nznodes,"ttm:sum_vsq");
  memory->create(sum_mass_vsq,nxnodes,nynodes,nznodes,"ttm:sum_mass_vsq");
  memory->create(sum_vsq_all,nxnodes,nynodes,nznodes,"ttm:sum_vsq_all");
  memory->create(sum_mass_vsq_all,nxnodes,nynodes,nznodes,
		 "ttm:sum_mass_vsq_all");
  memory->create(T_electron_old,nxnodes,nynodes,nznodes,"ttm:T_electron_old");
  memory->create(T_electron,nxnodes,nynodes,nznodes,"ttm:T_electron");
  memory->create(net_energy_transfer,nxnodes,nynodes,nznodes,
		 "TTM:net_energy_transfer");
  memory->create(net_energy_transfer_all,nxnodes,nynodes,nznodes,
		 "TTM:net_energy_transfer_all");
 
  flangevin = NULL;
  grow_arrays(atom->nmax);

  // zero out the flangevin array

  for (int i = 0; i < atom->nmax; i++) {
    flangevin[i][0] = 0;
    flangevin[i][1] = 0;
    flangevin[i][2] = 0;
  }

  atom->add_callback(0);
  atom->add_callback(1);

  // set initial electron temperatures from user input file

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

  memory->destroy(nsum);
  memory->destroy(nsum_all);
  memory->destroy(T_initial_set);
  memory->destroy(sum_vsq);
  memory->destroy(sum_mass_vsq);
  memory->destroy(sum_vsq_all);
  memory->destroy(sum_mass_vsq_all);
  memory->destroy(T_electron_old);
  memory->destroy(T_electron);
  memory->destroy(flangevin);
  memory->destroy(net_energy_transfer);
  memory->destroy(net_energy_transfer_all); 
}

/* ---------------------------------------------------------------------- */

int FixTTM::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTTM::init()
{
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix ttm with 2d simulation");
  if (domain->nonperiodic != 0)
    error->all(FLERR,"Cannot use nonperiodic boundares with fix ttm");
  if (domain->triclinic)
    error->all(FLERR,"Cannot use fix ttm with triclinic box");

  // set force prefactors

  for (int i = 1; i <= atom->ntypes; i++) {
    gfactor1[i] = - gamma_p / force->ftm2v;
    gfactor2[i] = 
      sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
  }

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        net_energy_transfer_all[ixnode][iynode][iznode] = 0;

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixTTM::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force_setup(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa_setup(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force(int vflag)
{
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

      if (T_electron[ixnode][iynode][iznode] < 0) 
        error->all(FLERR,"Electronic temperature dropped below zero");

      double tsqrt = sqrt(T_electron[ixnode][iynode][iznode]);

      gamma1 = gfactor1[type[i]];
      double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      if (vsq > v_0_sq) gamma1 *= (gamma_p + gamma_s)/gamma_p;
      gamma2 = gfactor2[type[i]] * tsqrt;

      flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
      flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
      flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);

      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_setup(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // apply langevin forces that have been stored from previous run

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_respa_setup(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force_setup(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTM::reset_dt()
{
  for (int i = 1; i <= atom->ntypes; i++)
    gfactor2[i] = 
      sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
}

/* ----------------------------------------------------------------------
   read in initial electron temperatures from a user-specified file
   only called by proc 0
------------------------------------------------------------------------- */

void FixTTM::read_initial_electron_temperatures()
{
  char line[MAXLINE];

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        T_initial_set[ixnode][iynode][iznode] = 0;

  // read initial electron temperature values from file

  int ixnode,iynode,iznode;
  double T_tmp;
  while (1) {
    if (fgets(line,MAXLINE,fpr) == NULL) break;
    sscanf(line,"%d %d %d %lg",&ixnode,&iynode,&iznode,&T_tmp);
    if (T_tmp < 0.0) error->one(FLERR,"Fix ttm electron temperatures must be > 0.0");
    T_electron[ixnode][iynode][iznode] = T_tmp;
    T_initial_set[ixnode][iynode][iznode] = 1;
  }

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        if (T_initial_set[ixnode][iynode][iznode] == 0)
          error->one(FLERR,"Initial temperatures not all set in fix ttm");

  // close file

  fclose(fpr);
}

/* ---------------------------------------------------------------------- */

void FixTTM::end_of_step()
{
  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
 
  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        net_energy_transfer[ixnode][iynode][iznode] = 0;
      
  for (int i = 0; i < nlocal; i++)
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
      net_energy_transfer[ixnode][iynode][iznode] += 
	(flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] + 
	 flangevin[i][2]*v[i][2]);
    }
  
  MPI_Allreduce(&net_energy_transfer[0][0][0],
		&net_energy_transfer_all[0][0][0],
		total_nnodes,MPI_DOUBLE,MPI_SUM,world);

  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nznodes;
  double del_vol = dx*dy*dz;

  // num_inner_timesteps = # of inner steps (thermal solves)
  // required this MD step to maintain a stable explicit solve

  int num_inner_timesteps = 1;
  double inner_dt = update->dt;
  double stability_criterion = 1.0 - 
    2.0*inner_dt/(electronic_specific_heat*electronic_density) * 
    (electronic_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
  if (stability_criterion < 0.0) {
    inner_dt = 0.5*(electronic_specific_heat*electronic_density) / 
      (electronic_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
    num_inner_timesteps = static_cast<int>(update->dt/inner_dt) + 1;
    inner_dt = update->dt/double(num_inner_timesteps);
    if (num_inner_timesteps > 1000000) 
      error->warning(FLERR,"Too many inner timesteps in fix ttm",0);
  }

  for (int ith_inner_timestep = 0; ith_inner_timestep < num_inner_timesteps; 
       ith_inner_timestep++) {

    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++)
          T_electron_old[ixnode][iynode][iznode] =
            T_electron[ixnode][iynode][iznode];
    
    // compute new electron T profile

    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
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
          T_electron[ixnode][iynode][iznode] = 
            T_electron_old[ixnode][iynode][iznode] + 
            inner_dt/(electronic_specific_heat*electronic_density) * 
            (electronic_thermal_conductivity *
             ((T_electron_old[right_xnode][iynode][iznode] + 
               T_electron_old[left_xnode][iynode][iznode] - 
               2*T_electron_old[ixnode][iynode][iznode])/dx/dx +
              (T_electron_old[ixnode][right_ynode][iznode] + 
               T_electron_old[ixnode][left_ynode][iznode] - 
               2*T_electron_old[ixnode][iynode][iznode])/dy/dy +
              (T_electron_old[ixnode][iynode][right_znode] + 
               T_electron_old[ixnode][iynode][left_znode] - 
               2*T_electron_old[ixnode][iynode][iznode])/dz/dz) -
              (net_energy_transfer_all[ixnode][iynode][iznode])/del_vol);   
        }
  }

  // output nodal temperatures for current timestep

  if ((nfileevery) && !(update->ntimestep % nfileevery)) { 

    // compute atomic Ta for each grid point
  
    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++) {
          nsum[ixnode][iynode][iznode] = 0;
          nsum_all[ixnode][iynode][iznode] = 0;
          sum_vsq[ixnode][iynode][iznode] = 0.0; 
          sum_mass_vsq[ixnode][iynode][iznode] = 0.0;
          sum_vsq_all[ixnode][iynode][iznode] = 0.0;
          sum_mass_vsq_all[ixnode][iynode][iznode] = 0.0;
        } 

    double massone;  
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
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
        sum_mass_vsq[ixnode][iynode][iznode] += massone*vsq;
      }
  
    MPI_Allreduce(&nsum[0][0][0],&nsum_all[0][0][0],total_nnodes,
                  MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&sum_vsq[0][0][0],&sum_vsq_all[0][0][0],total_nnodes,
                  MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&sum_mass_vsq[0][0][0],&sum_mass_vsq_all[0][0][0],
                  total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  
    if (me == 0) {
      fprintf(fp,BIGINT_FORMAT,update->ntimestep);

      double T_a;
      for (int ixnode = 0; ixnode < nxnodes; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            T_a = 0;
            if (nsum_all[ixnode][iynode][iznode] > 0)
              T_a = sum_mass_vsq_all[ixnode][iynode][iznode]/
                (3.0*force->boltz*nsum_all[ixnode][iynode][iznode]/force->mvv2e);
            fprintf(fp," %f",T_a);
          }
            
      fprintf(fp,"\t");
      for (int ixnode = 0; ixnode < nxnodes; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++)
            fprintf(fp,"%f ",T_electron[ixnode][iynode][iznode]);
      fprintf(fp,"\n");
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of 3d grid
------------------------------------------------------------------------- */

double FixTTM::memory_usage()
{
  double bytes = 0.0;
  bytes += 5*total_nnodes * sizeof(int);
  bytes += 14*total_nnodes * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixTTM::grow_arrays(int ngrow) 
{

 memory->grow(flangevin,ngrow,3,"TTM:flangevin");  

}

/* ----------------------------------------------------------------------
  return the energy of the electronic subsystem or the net_energy transfer
   between the subsystems
------------------------------------------------------------------------- */

double FixTTM::compute_vector(int n)
{
  double e_energy = 0.0;
  double transfer_energy = 0.0;

  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nznodes;
  double del_vol = dx*dy*dz;

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        e_energy += 
	  T_electron[ixnode][iynode][iznode]*electronic_specific_heat* 
	  electronic_density*del_vol;
        transfer_energy += 
	  net_energy_transfer_all[ixnode][iynode][iznode]*update->dt;
  }

  if (n == 0) return e_energy;
  if (n == 1) return transfer_energy;
  return 0.0;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write 
------------------------------------------------------------------------- */

void FixTTM::write_restart(FILE *fp)
{
  double *rlist;
  memory->create(rlist,nxnodes*nynodes*nznodes+1,"TTM:rlist");

  int n = 0;
  rlist[n++] = seed;
 
  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        rlist[n++] =  T_electron[ixnode][iynode][iznode];
  
  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(rlist,sizeof(double),n,fp);
  }

  memory->destroy(rlist);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix 
------------------------------------------------------------------------- */

void FixTTM::restart(char *buf)
{
  int n = 0;
  double *rlist = (double *) buf;
  
  // the seed must be changed from the initial seed

  seed = static_cast<int> (0.5*rlist[n++]);
  
  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        T_electron[ixnode][iynode][iznode] = rlist[n++];
  
  delete random;
  random = new RanMars(lmp,seed+comm->me);
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixTTM::pack_restart(int i, double *buf)
{
  buf[0] = 4;
  buf[1] = flangevin[i][0];
  buf[2] = flangevin[i][1];
  buf[3] = flangevin[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixTTM::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  flangevin[nlocal][0] = extra[nlocal][m++];
  flangevin[nlocal][1] = extra[nlocal][m++];
  flangevin[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixTTM::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixTTM::size_restart(int nlocal)
{
  return 4;
}
