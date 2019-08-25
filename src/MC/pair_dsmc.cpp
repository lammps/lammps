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
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include "pair_dsmc.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "update.h"
#include "random_mars.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairDSMC::PairDSMC(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;

  total_number_of_collisions = 0;
  max_particles = max_particle_list = 0;
  next_particle = NULL;
  random = NULL;
}

/* ---------------------------------------------------------------------- */

PairDSMC::~PairDSMC()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(sigma);
    memory->destroy(cut);
    memory->destroy(V_sigma_max);
    memory->destroy(particle_list);
    memory->destroy(first);
    memory->destroy(number);
  }

  delete [] next_particle;
  delete random;
}

/* ---------------------------------------------------------------------- */

void PairDSMC::compute(int /*eflag*/, int /*vflag*/)
{
  double **x = atom->x;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 1; i <= atom->ntypes; ++i)
    for (int j = 0; j < total_ncells; ++j) {
      first[i][j] = -1;
      number[i][j] = 0;
    }

  if (atom->nmax > max_particles) {
    delete [] next_particle;
    max_particles = atom->nmax;
    next_particle = new int[max_particles];
  }

  // find each particle's cell and sort by type
  // assume a constant volume and shape simulation domain
  // skip particle if outside processor domain

  for (int i = 0; i < nlocal; ++i) {
    int xcell = static_cast<int>((x[i][0] - domain->boxlo[0])/cellx);
    int ycell = static_cast<int>((x[i][1] - domain->boxlo[1])/celly);
    int zcell = static_cast<int>((x[i][2] - domain->boxlo[2])/cellz);

    if ((xcell < 0) || (xcell > ncellsx-1) ||
        (ycell < 0) || (ycell > ncellsy-1) ||
        (zcell < 0) || (zcell > ncellsz-1)) continue;

    int icell = xcell + ycell*ncellsx + zcell*ncellsx*ncellsy;
    itype = type[i];
    next_particle[i] = first[itype][icell];
    first[itype][icell] = i;
    number[itype][icell]++;
  }

  for (int icell = 0; icell < total_ncells; ++icell) {

    for (itype = 1; itype <= atom->ntypes; ++itype) {
      number_of_A = number[itype][icell];
      if (number_of_A > max_particle_list) {
        max_particle_list = number_of_A;
        memory->grow(particle_list,atom->ntypes+1,max_particle_list,
                     "pair:particle_list");
      }

      int m = first[itype][icell];
      for (int k = 0; k < number_of_A; k++) {
        particle_list[itype][k] = m;
        m = next_particle[m];
      }
    }

    for (itype = 1; itype <= atom->ntypes; ++itype) {
      imass = mass[itype];
      number_of_A = number[itype][icell];

      for (jtype = itype; jtype <= atom->ntypes; ++jtype) {
        jmass = mass[jtype];
        number_of_B = number[jtype][icell];

        reduced_mass = imass*jmass/(imass + jmass);
        total_mass = imass + jmass;
        jmass_tmass = jmass/total_mass;
        imass_tmass = imass/total_mass;

        // if necessary, recompute V_sigma_max values

        if (recompute_vsigmamax_stride &&
            (update->ntimestep % recompute_vsigmamax_stride == 0))
          recompute_V_sigma_max(icell);

        // # of collisions to perform for itype-jtype pairs

        double &Vs_max = V_sigma_max[itype][jtype];
        double num_of_collisions_double = number_of_A * number_of_B *
          weighting * Vs_max * update->dt / vol;

        if ((itype == jtype) && number_of_B)
          num_of_collisions_double *=
            0.5 * double(number_of_B - 1) / double(number_of_B);

        int num_of_collisions =
          convert_double_to_equivalent_int(num_of_collisions_double);

        if (num_of_collisions > number_of_A)
          error->warning(FLERR,"Pair dsmc: num_of_collisions > number_of_A",0);
        if (num_of_collisions > number_of_B)
          error->warning(FLERR,"Pair dsmc: num_of_collisions > number_of_B",0);

        // perform collisions on pairs of particles in icell

        for (int k = 0; k < num_of_collisions; k++) {
          if ((number_of_A < 1) || (number_of_B < 1)) break;
          if ((itype == jtype) && (number_of_A < 2)) break;
          int ith_A = static_cast<int>(random->uniform()*number_of_A);
          int jth_B = static_cast<int>(random->uniform()*number_of_B);
          int i = particle_list[itype][ith_A];
          int j = particle_list[jtype][jth_B];
          if (i == j) {
            k--;
            continue;
          }
          double probability = V_sigma(i,j)/Vs_max;
          if (probability > random->uniform()) scatter_random(i,j,icell);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairDSMC::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(V_sigma_max,n+1,n+1,"pair:V_sigma_max");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairDSMC::settings(int narg, char **arg)
{
  if (narg != 6) error->all(FLERR,"Illegal pair_style command");

  cut_global = 0.0;
  max_cell_size = force->numeric(FLERR,arg[0]);
  seed = force->inumeric(FLERR,arg[1]);
  weighting = force->numeric(FLERR,arg[2]);
  T_ref = force->numeric(FLERR,arg[3]);
  recompute_vsigmamax_stride = force->inumeric(FLERR,arg[4]);
  vsigmamax_samples = force->inumeric(FLERR,arg[5]);

  // initialize Marsaglia RNG with processor-unique seed

  if (max_cell_size <= 0.0) error->all(FLERR,"Illegal pair_style command");
  if (seed <= 0) error->all(FLERR,"Illegal pair_style command");
  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);

  kT_ref = force->boltz*T_ref;

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDSMC::coeff(int narg, char **arg)
{
  if (narg < 3 || narg > 4) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double sigma_one = force->numeric(FLERR,arg[2]);

  double cut_one = cut_global;
  if (narg == 4) cut_one = force->numeric(FLERR,arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDSMC::init_style()
{
  ncellsx = ncellsy = ncellsz = 1;
  while (((domain->boxhi[0] - domain->boxlo[0])/ncellsx) > max_cell_size)
    ncellsx++;
  while (((domain->boxhi[1] - domain->boxlo[1])/ncellsy) > max_cell_size)
    ncellsy++;
  while (((domain->boxhi[2] - domain->boxlo[2])/ncellsz) > max_cell_size)
    ncellsz++;

  cellx = (domain->boxhi[0] - domain->boxlo[0])/ncellsx;
  celly = (domain->boxhi[1] - domain->boxlo[1])/ncellsy;
  cellz = (domain->boxhi[2] - domain->boxlo[2])/ncellsz;

  if (comm->me == 0) {
    if (screen) fprintf(screen,"DSMC cell size = %g x %g x %g\n",
                        cellx,celly,cellz);
    if (logfile) fprintf(logfile,"DSMC cell size = %g x %g x %g\n",
                         cellx,celly,cellz);
  }

  total_ncells = ncellsx*ncellsy*ncellsz;
  vol = cellx*celly*cellz;

  memory->create(particle_list,atom->ntypes+1,0,"pair:particle_list");
  memory->create(first,atom->ntypes+1,total_ncells,"pair:first");
  memory->create(number,atom->ntypes+1,total_ncells,"pair:number");

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = 1; j <= atom->ntypes; j++)
      V_sigma_max[i][j] = 0.0;

  two_pi = 8.0*atan(1.0);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDSMC::init_one(int i, int j)
{
  if (setflag[i][j] == 0) cut[i][j] = 0.0;
  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDSMC::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDSMC::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDSMC::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&max_cell_size,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDSMC::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&max_cell_size,sizeof(double),1,fp);
    fread(&seed,sizeof(int),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }

  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&max_cell_size,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/*-------------------------------------------------------------------------
  rezero and recompute the V_sigma_max values this timestep for use during
  the next nrezero timesteps
-------------------------------------------------------------------------*/

void PairDSMC::recompute_V_sigma_max(int /*icell*/)
{
  int i,j,k;
  double Vsigma_max = 0;

  if (number_of_A && number_of_B) {
    for (k = 0; k < vsigmamax_samples; k++) {
      i = particle_list[itype]
        [static_cast<int>(random->uniform()*number_of_A)];
      j = particle_list[jtype]
        [static_cast<int>(random->uniform()*number_of_B)];
      if (i == j) continue;
      Vsigma_max = MAX(Vsigma_max,V_sigma(i,j));
    }
  }
  V_sigma_max[itype][jtype] = Vsigma_max;
}

/*-------------------------------------------------------------------------
  VHS model
  compute the velocity vector difference between i and j and multiply by
  their combined collision cross section, sigma, for neutral-neutral
  collisions using the Variable Hard Sphere model
-------------------------------------------------------------------------*/

double PairDSMC::V_sigma(int i, int j)
{
  double relative_velocity_sq,relative_velocity,pair_sigma;
  double delv[3];
  double *vi = atom->v[i];
  double *vj = atom->v[j];

  subtract3d(vi,vj,delv);
  relative_velocity_sq = dot3d(delv,delv);
  relative_velocity = sqrt(relative_velocity_sq);

  // from Bird eq 4.63, and omega=0.67
  // (omega - 0.5) = 0.17
  // 1/GAMMA(2.5 - omega) = 1.06418029298371

  if (relative_velocity_sq != 0.0)
    pair_sigma = sigma[itype][jtype]*
      pow(kT_ref/(0.5*reduced_mass*relative_velocity_sq),0.17) *
      1.06418029298371;
  else
    pair_sigma = 0.0;

  return relative_velocity*pair_sigma;
}

/*-------------------------------------------------------------------------
  generate new velocities for collided particles
-------------------------------------------------------------------------*/

void PairDSMC::scatter_random(int i, int j, int /*icell*/)
{
  double mag_delv,cos_phi,cos_squared,r,theta;
  double delv[3],vcm[3];
  double *vi = atom->v[i];
  double *vj = atom->v[j];

  subtract3d(vi,vj,delv);
  if (itype == jtype) mag_delv = sqrt(dot3d(delv,delv))*0.5;
  else mag_delv = sqrt(dot3d(delv,delv));

  cos_phi = 1.0 - (2.0*random->uniform());
  cos_squared = MIN(1.0,cos_phi*cos_phi);
  r = sqrt(1.0 - cos_squared);
  delv[0] = cos_phi*mag_delv;
  theta = two_pi*random->uniform();
  delv[1] = r*mag_delv*cos(theta);
  delv[2] = r*mag_delv*sin(theta);

  if (itype == jtype) {
    vcm[0] = (vi[0]+vj[0])*0.5;
    vcm[1] = (vi[1]+vj[1])*0.5;
    vcm[2] = (vi[2]+vj[2])*0.5;
    vi[0] = vcm[0] + delv[0];
    vi[1] = vcm[1] + delv[1];
    vi[2] = vcm[2] + delv[2];
    vj[0] = vcm[0] - delv[0];
    vj[1] = vcm[1] - delv[1];
    vj[2] = vcm[2] - delv[2];
  } else {
    vcm[0] = vi[0]*imass_tmass + vj[0]*jmass_tmass;
    vcm[1] = vi[1]*imass_tmass + vj[1]*jmass_tmass;
    vcm[2] = vi[2]*imass_tmass + vj[2]*jmass_tmass;
    vi[0] = vcm[0] + delv[0]*jmass_tmass;
    vi[1] = vcm[1] + delv[1]*jmass_tmass;
    vi[2] = vcm[2] + delv[2]*jmass_tmass;
    vj[0] = vcm[0] - delv[0]*imass_tmass;
    vj[1] = vcm[1] - delv[1]*imass_tmass;
    vj[2] = vcm[2] - delv[2]*imass_tmass;
  }

  total_number_of_collisions++;
}

/* ----------------------------------------------------------------------
   This method converts the double supplied by the calling function into
   an int, which is returned. By adding a random number between 0 and 1
   to the double before converting it to an int, we ensure that,
   statistically, we round down with probability identical to the
   remainder and up the rest of the time. So even though we're using an
   integer, we're statistically matching the exact expression represented
   by the double.
------------------------------------------------------------------------- */

int PairDSMC::convert_double_to_equivalent_int(double input_double)
{
  if (input_double > INT_MAX)
    error->all(FLERR,"Tried to convert a double to int, but input_double > INT_MAX");

  int output_int = static_cast<int>(input_double + random->uniform());
  return output_int;
}
