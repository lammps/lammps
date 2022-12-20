// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Shen,Yuan, Qi,Tingting, and Reed,Evan
   Implementation of the colored thermostat for quantum nuclear effects
------------------------------------------------------------------------- */

#include "fix_qtb.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ----------------------------------------------------------------------
   read parameters
------------------------------------------------------------------------- */
FixQTB::FixQTB(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix qtb command");

  // default parameters

  t_target = 300.0;
  t_period = 1.0;
  fric_coef = 1/t_period;
  seed = 880302;
  f_max = 200.0;
  N_f = 100;

  // reading parameters
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qtb command");
      t_target = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (t_target < 0.0) error->all(FLERR,"Fix qtb temp must be >= 0.0");
      iarg += 2;
    } else if (strcmp(arg[iarg],"damp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qtb command");
      t_period = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (t_period <= 0.0) error->all(FLERR,"Fix qtb damp must be > 0.0");
      fric_coef = 1/t_period;
      iarg += 2;
    } else if (strcmp(arg[iarg],"seed") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qtb command");
      seed = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (seed <= 0) error->all(FLERR,"Illegal fix qtb command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"f_max") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qtb command");
      f_max = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (f_max <= 0) error->all(FLERR,"Illegal fix qtb command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"N_f") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qtb command");
      N_f = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (N_f <= 0) error->all(FLERR,"Illegal fix qtb command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix qtb command");
  }

  maxexchange = 6*N_f+3;

  // allocate qtb
  gfactor1 = nullptr;
  gfactor3 = nullptr;
  omega_H = nullptr;
  time_H = nullptr;
  random_array_0 = nullptr;
  random_array_1 = nullptr;
  random_array_2 = nullptr;
  fran = nullptr;
  id_temp = nullptr;
  temperature = nullptr;

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);

  // allocate per-type arrays for force prefactors
  gfactor1 = new double[atom->ntypes+1];
  gfactor3 = new double[atom->ntypes+1];

  // allocate random-arrays and fran
  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  // allocate omega_H and time_H
  memory->create(omega_H,2*N_f,"qtb:omega_H");
  memory->create(time_H,2*N_f,"qtb:time_H");
}

/* ----------------------------------------------------------------------
   release memories
------------------------------------------------------------------------- */
FixQTB::~FixQTB()
{
  delete random;
  delete [] gfactor1;
  delete [] gfactor3;
  delete [] id_temp;
  memory->destroy(fran);
  memory->destroy(random_array_0);
  memory->destroy(random_array_1);
  memory->destroy(random_array_2);
  memory->destroy(omega_H);
  memory->destroy(time_H);
  atom->delete_callback(id,Atom::GROW);
}

/* ----------------------------------------------------------------------
   setmask
------------------------------------------------------------------------- */
int FixQTB::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ----------------------------------------------------------------------
   fix initiation
------------------------------------------------------------------------- */
void FixQTB::init()
{
  // copy parameters from other classes
  double dtv = update->dt;
  if (atom->mass == nullptr)
    error->all(FLERR,"Cannot use fix msst without per-type mass defined");

  //initiate the counter \mu
  counter_mu=0;

  //set up the h time step for updating the random force \delta{}h=\frac{\pi}{\Omega_{max}}
  if (int(1.0/(2*f_max*dtv)) == 0) {
    if (comm->me == 0) error->warning(FLERR,"Either f_max is too high or the time step "
                                      "is too big, setting f_max to be 1/timestep!\n");
    h_timestep=dtv;
    alpha=1;
  } else {
    alpha=int(1.0/(2*f_max*dtv));
    h_timestep=alpha*dtv;
  }
  if (comm->me == 0 && screen)
    fmt::print(screen,"The effective maximum frequency is now {} inverse time unit "
               "with alpha value as {}!\n", 0.5/h_timestep, alpha);

  // set force prefactors
  if (!atom->rmass) {
    for (int i = 1; i <= atom->ntypes; i++) {
      //gfactor1 is the friction force \gamma{}m_{i}\frac{dv}{dt}
      gfactor1[i] = (atom->mass[i]*fric_coef) / force->ftm2v;
      //gfactor3 is the random force \sqrt{\frac{2\gamma{}m_{i}}{\alpha*\delta{}t}}, \sqrt{12} makes the random array variance equal to unit.
      gfactor3[i] = sqrt(2*fric_coef*atom->mass[i])*sqrt(force->mvv2e)*sqrt(12/h_timestep);  //this still leaves a square energy term from the power spectrum H.
    }
  }

  // generate random number array with zero mean and variance equal 1/12.
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    for (int m = 0; m < 2*N_f; m++) {
      random_array_0[i][m] = random->uniform()-0.5;
      random_array_1[i][m] = random->uniform()-0.5;
      random_array_2[i][m] = random->uniform()-0.5;
    }
  }

  // load omega_H with calculated spectrum at a specific temperature (corrected spectrum), omega_H is the Fourier transformation of time_H
  for (int k = 0; k < 2*N_f; k++) {
    double f_k=(k-N_f)/(2*N_f*h_timestep);  //\omega_k=\frac{2\pi}{\delta{}h}\frac{k}{2N_f} for k from -N_f to N_f-1
    if (k == N_f) {
      omega_H[k]=sqrt(force->boltz * t_target);
    } else {
      double energy_k= force->hplanck * fabs(f_k);
      omega_H[k]=sqrt( energy_k * (0.5+1.0/( exp(energy_k/(force->boltz * t_target)) - 1.0 )) );
      omega_H[k]*=alpha*sin((k-N_f)*MY_PI/(2*alpha*N_f))/sin((k-N_f)*MY_PI/(2*N_f));
    }
  }

  // construct the signal filter H, filter has the unit of of sqrt(energy) \sqrt{2N_f}^{-1}H\left(t_n\right)
  for (int n = 0; n < 2*N_f; n++) {
    time_H[n] = 0;
    double t_n=(n-N_f);
    for (int k = 0; k < 2*N_f; k++) {
      double omega_k=(k-N_f)*MY_PI/N_f;
      time_H[n] += omega_H[k]*(cos(omega_k*t_n));
    }
    time_H[n]/=(2.0*N_f);
   }

  // respa
  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
}

/* ----------------------------------------------------------------------
   no MD, so setup returns post force
------------------------------------------------------------------------- */
void FixQTB::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(nlevels_respa-1);
  }
}

/* ----------------------------------------------------------------------
   post_force
------------------------------------------------------------------------- */
void FixQTB::post_force(int /*vflag*/)
{
  double gamma1,gamma3;

  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  bigint nlocal = atom->nlocal;
  bigint ntotal = atom->natoms;

  //update the colored random force every alpha MD steps
  if (counter_mu == alpha) {
    //propagate h_timestep ahead
    for (int j = 0; j < nlocal; j++) {

      //update random array
      for (int m = 0; m < 2*N_f-1; m++) {
            random_array_0[j][m] = random_array_0[j][m+1];
            random_array_1[j][m] = random_array_1[j][m+1];
            random_array_2[j][m] = random_array_2[j][m+1];
      }
      random_array_0[j][2*N_f-1] = random->uniform()-0.5;
      random_array_1[j][2*N_f-1] = random->uniform()-0.5;
      random_array_2[j][2*N_f-1] = random->uniform()-0.5;
    }

    //reset counter \mu
    counter_mu=0;
  }

  if (counter_mu == 0) {
    for (int j = 0; j < nlocal; j++) {
      fran[j][0] = 0.0;
      fran[j][1] = 0.0;
      fran[j][2] = 0.0;

      //reset random force
      if (mask[j] & groupbit) {
        gamma3 = gfactor3[type[j]];

        for (int m = 0; m < 2*N_f; m++) {
          fran[j][0] += time_H[m] * random_array_0[j][2*N_f-m-1];
          fran[j][1] += time_H[m] * random_array_1[j][2*N_f-m-1];
          fran[j][2] += time_H[m] * random_array_2[j][2*N_f-m-1];
        }
        fran[j][0] = fran[j][0]*gamma3;
        fran[j][1] = fran[j][1]*gamma3;
        fran[j][2] = fran[j][2]*gamma3;
      }
    }
  }

  //reset all the force sums
  fsum[0]=0.0; fsumall[0]=0.0;
  fsum[1]=0.0; fsumall[1]=0.0;
  fsum[2]=0.0; fsumall[2]=0.0;

  for (int j = 0; j < nlocal; j++) {
    //sum over each atom
    if (mask[j] & groupbit) {
      gamma1 = gfactor1[type[j]];

      fsum[0]+=fran[j][0]-gamma1*v[j][0];
      fsum[1]+=fran[j][1]-gamma1*v[j][1];
      fsum[2]+=fran[j][2]-gamma1*v[j][2];
    }
  }

  //compute force sums
  MPI_Allreduce(fsum,fsumall,3,MPI_DOUBLE,MPI_SUM,world);

  //implement random forces
  for (int j = 0; j < nlocal; j++) {
    //make sure there is no net force on the system
    f[j][0] -= fsumall[0]/ntotal;
    f[j][1] -= fsumall[1]/ntotal;
    f[j][2] -= fsumall[2]/ntotal;

    if (mask[j] & groupbit) {
      gamma1 = gfactor1[type[j]];

      f[j][0]+=fran[j][0]-gamma1*v[j][0];
      f[j][1]+=fran[j][1]-gamma1*v[j][1];
      f[j][2]+=fran[j][2]-gamma1*v[j][2];
    }
  }

  //move 1 step forward
  counter_mu++;
}

/* ----------------------------------------------------------------------
   post_force_respa
------------------------------------------------------------------------- */
void FixQTB::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   modifications of fix qtb
------------------------------------------------------------------------- */
int FixQTB::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    delete [] id_temp;
    id_temp = utils::strdup(arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,"Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   memory usage of fix qtb
------------------------------------------------------------------------- */
double FixQTB::memory_usage()
{
  double bytes = 0.0;
  // random_arrays memory usage
  bytes += (double)(atom->nmax* 6*N_f * sizeof(double));
  // fran memory usage
  bytes += (double)(atom->nmax* 3 * sizeof(double));
  bytes += (double)(4*N_f * sizeof(double));
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array for fran and random_array
------------------------------------------------------------------------- */
void FixQTB::grow_arrays(int nmax)
{
  memory->grow(random_array_0,nmax,2*N_f,"qtb:random_array_0");
  memory->grow(random_array_1,nmax,2*N_f,"qtb:random_array_1");
  memory->grow(random_array_2,nmax,2*N_f,"qtb:random_array_2");
  memory->grow(fran,nmax,3,"qtb:fran");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */
void FixQTB::copy_arrays(int i, int j, int /*delflag*/)
{
  for (int m = 0; m < 2*N_f; m++) {
    random_array_0[j][m] = random_array_0[i][m];
    random_array_1[j][m] = random_array_1[i][m];
    random_array_2[j][m] = random_array_2[i][m];
  }

  for (int m = 0; m < 3; m++)
    fran[j][m] = fran[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */
int FixQTB::pack_exchange(int i, double *buf)
{
  int n = 0;
  for (int m = 0; m < 2*N_f; m++) buf[n++] = random_array_0[i][m];
  for (int m = 0; m < 2*N_f; m++) buf[n++] = random_array_1[i][m];
  for (int m = 0; m < 2*N_f; m++) buf[n++] = random_array_2[i][m];
  for (int m = 0; m < 3; m++) buf[n++] = fran[i][m];
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */
int FixQTB::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  for (int m = 0; m < 2*N_f; m++) random_array_0[nlocal][m] = buf[n++];
  for (int m = 0; m < 2*N_f; m++) random_array_1[nlocal][m] = buf[n++];
  for (int m = 0; m < 2*N_f; m++) random_array_2[nlocal][m] = buf[n++];
  for (int m = 0; m < 3; m++) fran[nlocal][m] = buf[n++];
  return n;
}
