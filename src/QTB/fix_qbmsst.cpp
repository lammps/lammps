// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Shen,Yuan, Qi,Tingting, and Reed,Evan
   Implementation of the Multi-Scale Shock Method with quantum nuclear effects
------------------------------------------------------------------------- */

#include "fix_qbmsst.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ----------------------------------------------------------------------
   read parameters
------------------------------------------------------------------------- */

FixQBMSST::FixQBMSST(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix qbmsst command");

  if (strcmp(arg[3],"x") == 0) {
    direction = 0;
    box_change |= BOX_CHANGE_X;
  } else if (strcmp(arg[3],"y") == 0) {
    direction = 1;
    box_change |= BOX_CHANGE_Y;
  } else if (strcmp(arg[3],"z") == 0) {
    direction = 2;
    box_change |= BOX_CHANGE_Z;
  } else {
    error->all(FLERR,"Illegal fix qbmsst command");
  }
  velocity = atof(arg[4]);
  if (velocity < 0)
    error->all(FLERR,"Illegal fix qbmsst command");

  // default parameters

  global_freq = 1;
  extscalar = 1;
  extvector = 0;
  nevery = 1;
  restart_global = 1;
  time_integrate = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 5;
  ecouple_flag = 1;

  qmass = 1.0e1;
  mu = 0.0;
  p0 = 0.0;
  v0 = 1.0;
  e0 = 0.0;
  p0_set = 0;
  v0_set = 0;
  e0_set = 0;
  tscale = 0.01;
  t_period = 1.0;
  fric_coef = 1/t_period;
  seed = 880302;
  f_max = 200.0;
  N_f = 100;
  eta = 1.0;
  beta = 100;
  t_init = 300.0;
  qtb_set = 0;

  // reading parameters

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"q") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qbmsst command");
      qmass = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (qmass < 0.0) error->all(FLERR,"Fix qbmsst qmass must be >= 0.0");
      iarg += 2;
    } else if (strcmp(arg[iarg],"mu") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qbmsst command");
      mu = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (mu < 0.0) error->all(FLERR,"Fix qbmsst mu must be >= 0.0");
      iarg += 2;
    } else if (strcmp(arg[iarg],"p0") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qbmsst command");
      p0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (p0 < 0.0) error->all(FLERR,"Fix qbmsst p0 must be >= 0.0");
      p0_set = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"v0") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qbmsst command");
      v0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (v0 < 0.0) error->all(FLERR,"Fix qbmsst v0 must be >= 0.0");
      v0_set = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"e0") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qbmsst command");
      e0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      e0_set = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"tscale") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qbmsst command");
      tscale = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (tscale < 0.0 || tscale > 1.0) error->all(FLERR,"Fix qbmsst tscale must satisfy 0 <= tscale < 1");
      iarg += 2;
    } else if (strcmp(arg[iarg],"damp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qbmsst command");
      t_period = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (t_period <= 0.0) error->all(FLERR,"Fix qbmsst damp must be > 0.0");
      fric_coef = 1/t_period;
      iarg += 2;
    } else if (strcmp(arg[iarg],"seed") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qbmsst command");
      seed = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (seed <= 0) error->all(FLERR,"Fix qbmsst seed must be a positive integer");
      iarg += 2;
    } else if (strcmp(arg[iarg],"f_max") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qbmsst command");
      f_max = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (f_max <= 0) error->all(FLERR,"Fix qbmsst f_max must be > 0.0");
      iarg += 2;
    } else if (strcmp(arg[iarg],"N_f") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qbmsst command");
      N_f = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (N_f <= 0) error->all(FLERR,"Fix qbmsst N_f must be a positive integer");
      iarg += 2;
    } else if (strcmp(arg[iarg],"eta") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qbmsst command");
      eta = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (eta <= 0) error->all(FLERR,"Fix qbmsst eta must be >= 0.0");
      iarg += 2;
    } else if (strcmp(arg[iarg],"beta") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qbmsst command");
      beta = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (beta <= 0) error->all(FLERR,"Fix qbmsst beta must be a positive integer");
      iarg += 2;
    } else if (strcmp(arg[iarg],"T_init") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qbmsst command");
      t_init = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (t_init <= 0) error->all(FLERR,"Fix qbmsst T_init must be >= 0.0");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix qbmsst command");
  }

  // check for periodicity in controlled dimensions
  if (domain->nonperiodic) error->all(FLERR,"Fix qbmsst requires a periodic box");

  maxexchange = 6*N_f+3;

  // comments
  if (comm->me == 0) {
    std::string msg = "QBMSST parameters:\n";

    if (direction == 0)      msg += "  Shock in x direction\n";
    else if (direction == 1) msg += "  Shock in y direction\n";
    else if (direction == 2) msg += "  Shock in z direction\n";

    msg += fmt::format("  Cell mass-like parameter qmass "
                       "(units of mass^2/length^4) = {:12.5e}\n", qmass);
    msg += fmt::format("  Shock velocity = {:12.5e}\n", velocity);
    msg += fmt::format("  Artificial viscosity (units of "
                       "mass/length/time) = {:12.5e}\n", mu);

    if (p0_set)
      msg += fmt::format("  Initial pressure specified to be {:12.5e}\n", p0);
    else msg += "  Initial pressure calculated on first step\n";
    if (v0_set)
      msg += fmt::format("  Initial volume specified to be {:12.5e}\n", v0);
    else msg += "  Initial volume calculated on first step\n";
    if (e0_set)
      msg += fmt::format("  Initial energy specified to be {:12.5e}\n", e0);
    else msg += "  Initial energy calculated on first step\n";
    utils::logmesg(lmp,msg);
  }

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  //   and thus its KE/temperature contribution should use group all

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} all temp",id_temp));
  tflag = 1;

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  id_press = utils::strdup(std::string(id) + "_press");
  modify->add_compute(fmt::format("{} all pressure {}",id_press, id_temp));
  pflag = 1;

  // create a new compute potential energy compute

  id_pe = utils::strdup(std::string(id) + "_pe");
  modify->add_compute(fmt::format("{} all pe",id_pe));
  peflag = 1;

  // allocate qbmsst
  temperature = nullptr;
  pressure = nullptr;
  pe = nullptr;
  old_velocity = nullptr;
  rfix = nullptr;
  gfactor = nullptr;
  random = nullptr;
  omega_H = nullptr;
  time_H = nullptr;
  random_array_0 = nullptr;
  random_array_1 = nullptr;
  random_array_2 = nullptr;
  fran = nullptr;

  // initialize Marsagxlia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);

  // allocate per-type arrays for force prefactors
  gfactor = new double[atom->ntypes+1];

  // allocate random-arrays and fran
  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  // allocate omega_H and time_H
  memory->create(omega_H,2*N_f,"qbmsst:omega_H");
  memory->create(time_H,2*N_f,"qbmsst:time_H");

  // initiate velocity record array
  memory->create(old_velocity,atom->nlocal,3,"qbmsst:old_velocity");
  atoms_allocated = atom->nlocal;
}

/* ----------------------------------------------------------------------
   release memories
------------------------------------------------------------------------- */

FixQBMSST::~FixQBMSST()
{
  delete [] rfix;
  delete [] gfactor;
  delete random;

  // delete temperature and pressure if fix created them
  if (tflag) modify->delete_compute(id_temp);
  if (pflag) modify->delete_compute(id_press);
  if (peflag) modify->delete_compute(id_pe);
  delete [] id_temp;
  delete [] id_press;
  delete [] id_pe;

  memory->destroy(old_velocity);
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

int FixQBMSST::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ----------------------------------------------------------------------
   fix initiation
------------------------------------------------------------------------- */

void FixQBMSST::init()
{
  // copy parameters from other classes
  dtv = update->dt;
  dthalf = 0.5 * update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  ntotal = atom->natoms;
  boltz = force->boltz;
  nktv2p = force->nktv2p;
  mvv2e = force->mvv2e;
  if (atom->mass == nullptr)
    error->all(FLERR,"Cannot use fix qbmsst without per-type mass defined");

  // set compute ptrs
  int itemp = modify->find_compute(id_temp);
  int ipress = modify->find_compute(id_press);
  int ipe = modify->find_compute(id_pe);
  if (itemp < 0 || ipress < 0|| ipe < 0)
    error->all(FLERR,"Could not find fix qbmsst compute ID");
  if (modify->compute[itemp]->tempflag == 0)
    error->all(FLERR,"Fix qbmsst compute ID does not compute temperature");
  if (modify->compute[ipress]->pressflag == 0)
    error->all(FLERR,"Fix qbmsst compute ID does not compute pressure");
  if (modify->compute[ipe]->peflag == 0)
    error->all(FLERR,"Fix qbmsst compute ID does not compute potential energy");
  temperature = modify->compute[itemp];
  pressure = modify->compute[ipress];
  pe = modify->compute[ipe];

  // initiate the counter l and \mu
  counter_l=0;
  counter_mu=0;

  // initiate qtb temperature
  if (!qtb_set) {
        t_current = t_init; qtb_set=1;
  }
  old_eavg = e0;

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

  //gfactor is the random force \sqrt{\frac{2\gamma{}m_{i}}{\alpha*\delta{}t}}, \sqrt{12} makes the random array variance equal to unit.
  for (int i = 1; i <= atom->ntypes; i++) {
    gfactor[i] = sqrt(2*fric_coef*atom->mass[i])*sqrt(force->mvv2e)*sqrt(12/h_timestep);//this still leaves a square energy term from the power spectrum H.
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

  // initiate dilations
  dilation[0] = dilation[1] = dilation[2] = 1.0;

  // initialize the time derivative of the volume.
  omega[0] = omega[1] = omega[2] = 0.0;

  // compute total mass
  double mass = 0.0;
  for (int i = 0; i < atom->nlocal; i++) mass += atom->mass[atom->type[i]];
  MPI_Allreduce(&mass,&total_mass,1,MPI_DOUBLE,MPI_SUM,world);

  // enable kspace summation of long range forces
  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  // detect if any fix rigid exist so rigid bodies move when box is dilated
  // rfix[] = indices to each fix rigid
  nrigid = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (utils::strmatch(modify->fix[i]->style,"^rigid") ||
        (strcmp(modify->fix[i]->style,"poems") == 0)) nrigid++;
  if (nrigid > 0) {
    rfix = new int[nrigid];
    nrigid = 0;
    for (int i = 0; i < modify->nfix; i++)
      if (utils::strmatch(modify->fix[i]->style,"^rigid") ||
          (strcmp(modify->fix[i]->style,"poems") == 0)) rfix[nrigid++] = i;
  }
}

/* ----------------------------------------------------------------------
   compute T,P before integrator starts
------------------------------------------------------------------------- */
void FixQBMSST::setup(int /*vflag*/)
{
  lagrangian_position = 0.0;

  temperature->compute_vector();
  pressure->compute_vector();
  couple();
  velocity_sum = compute_vsum();

  if (v0_set == 0) {
    v0 = compute_vol();
    v0_set = 1;
    if (comm->me == 0)
      utils::logmesg(lmp,"Fix QBMSST v0 = {:12.5e}\n", v0);
  }

  if (p0_set == 0) {
    p0 = p_current[direction];
    p0_set = 1;

    if (comm->me == 0)
      utils::logmesg(lmp,"Fix QBMSST p0 = {:12.5e}\n", p0);
  }

  if (e0_set == 0) {
    e0 = compute_etotal();
    e0_set = 1;
    old_eavg = e0;

    if (comm->me == 0)
      utils::logmesg(lmp,"Fix QBMSST e0 = to be {:12.5e}\n",e0);
  }

  temperature->compute_vector();
  double *ke_tensor = temperature->vector;
  double ke_temp = ke_tensor[0]+ke_tensor[1]+ke_tensor[2];
  if (ke_temp > 0.0 && tscale > 0.0) {

    // transfer energy from atom velocities to cell volume motion
    // to bias initial compression
    double **v = atom->v;
    int *mask = atom->mask;
    double sqrt_initial_temperature_scaling = sqrt(1.0-tscale);

    double fac1 =  tscale*total_mass/qmass*ke_temp/force->mvv2e;

    omega[direction]=-1*sqrt(fac1);
    double fac2 = omega[direction]/v0;

    if ( comm->me == 0 && tscale != 1.0)
      utils::logmesg(lmp,"Fix QBMSST initial strain rate of {:12.5e} "
                     "established by reducing temperature by "
                     "factor of {:12.5e}\n",fac2,tscale);
    for (int i = 0; i < atom->nlocal; i++) {
      if (mask[i] & groupbit) {
        for (int k = 0; k < 3; k++) {
          v[i][k]*=sqrt_initial_temperature_scaling;
        }
      }
    }
  }

  // trigger virial computation on next timestep
  pressure->addstep(update->ntimestep+1);
  pe->addstep(update->ntimestep+1);
}

/* ----------------------------------------------------------------------
   1st half of Verlet update
------------------------------------------------------------------------- */
void FixQBMSST::initial_integrate(int /*vflag*/)
{
  int sd;
  sd = direction;
  double p_qbmsst;// QBMSST driving pressure.
  int i, k;
  double vol;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  double **x = atom->x;
  double boltz = force->boltz;

  // check to see if old_velocity is correctly allocated
  check_alloc(nlocal);

  // compute new pressure and volume.
  temperature->compute_vector();
  pressure->compute_vector();
  couple();
  vol = compute_vol();

  // decide if the qtb temperature need to be updated or not
  if (counter_l == 0) {
    t_current -= dtv*fric_coef*eta*beta*(old_eavg-e0)/(3*ntotal*boltz);
    if (t_current > 0.0) {
      old_eavg = 0;//clear old energy average

      // load omega_H with calculated spectrum at a specific temperature (corrected spectrum), omega_H is the Fourier transformation of time_H
      for (int k = 0; k < 2*N_f; k++) {
        double f_k=(k-N_f)/(2*N_f*h_timestep);  //\omega_k=\frac{2\pi}{\delta{}h}\frac{k}{2N_f} for k from -N_f to N_f-1
        if (k == N_f) {
          omega_H[k]=sqrt(force->boltz * t_current);
        } else {
          double energy_k= force->hplanck * fabs(f_k);
          omega_H[k]=sqrt( energy_k * (0.5+1.0/( exp(energy_k/(force->boltz * t_current)) - 1.0 )) );
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
    }
  }

  //update the colored random force every alpha MD steps
  if (counter_mu == 0) {
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

      //reset random force
      fran[j][0] = 0.0;
      fran[j][1] = 0.0;
      fran[j][2] = 0.0;
      if (mask[j] & groupbit) {
        double gamma3 = gfactor[type[j]];

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

  //  estimate old energy average in this step
  old_eavg = old_eavg + compute_egrand()/beta;
  counter_l = (counter_l + 1) % beta;
  counter_mu = (counter_mu + 1) % alpha;

  // propagate the time derivative of
  // the volume 1/2 step at fixed vol, r, rdot.
  p_qbmsst = nktv2p * mvv2e * velocity * velocity * total_mass *
    ( v0 - vol)/( v0 * v0);
  double A = total_mass * ( p_current[sd] - p0 - p_qbmsst ) /
    (qmass * nktv2p * mvv2e);
  double B = total_mass * mu / ( qmass * vol );

  // prevent blow-up of the volume.
  if (vol > v0 && A > 0.0) {
    A = -A;
  }

  // use taylor expansion to avoid singularity at B == 0.
  if (B * dthalf > 1.0e-06) {
    omega[sd] = ( omega[sd] + A * ( exp(B * dthalf) - 1.0 ) / B )
      * exp(-B * dthalf);
  } else {
    omega[sd] = omega[sd] + (A - B * omega[sd]) * dthalf +
      0.5 * (B * B * omega[sd] - A * B ) * dthalf * dthalf;
  }

  // propagate velocity sum 1/2 step by
  // temporarily propagating the velocities.
  velocity_sum = compute_vsum();
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      for (k = 0; k < 3; k++) {
        double C = (f[i][k] + fran[i][k])* force->ftm2v / mass[type[i]];//  this term now has a random force part
        double D = mu * omega[sd] * omega[sd] /
          (velocity_sum * mass[type[i]] * vol ) - fric_coef;
        old_velocity[i][k] = v[i][k];
        if (k == direction) {
          D = D - 2.0 * omega[sd] / vol;
        }
        if (fabs(dthalf * D) > 1.0e-06) {
          double expd = exp(D * dthalf);
          v[i][k] = expd * ( C + D * v[i][k] - C / expd ) / D;
        } else {
          v[i][k] = v[i][k] + ( C + D * v[i][k] ) * dthalf +
            0.5 * (D * D * v[i][k] + C * D ) * dthalf * dthalf;
        }
      }
    }
  }

  velocity_sum = compute_vsum();

  // reset the velocities.
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      for (k = 0; k < 3; k++) {
        v[i][k] = old_velocity[i][k];
      }
    }
  }

  // propagate velocities 1/2 step using the new velocity sum.
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      for (k = 0; k < 3; k++) {
        double C = (f[i][k] + fran[i][k])* force->ftm2v / mass[type[i]];//  this term now has a random force part
        double D = mu * omega[sd] * omega[sd] /
          (velocity_sum * mass[type[i]] * vol ) - fric_coef;

        if (k == direction) {
          D = D - 2.0 * omega[sd] / vol;
        }
        if (fabs(dthalf * D) > 1.0e-06) {
          double expd = exp(D * dthalf);
          v[i][k] = expd * ( C + D * v[i][k] - C / expd ) / D;
        } else {
          v[i][k] = v[i][k] + ( C + D * v[i][k] ) * dthalf +
            0.5 * (D * D * v[i][k] + C * D ) * dthalf * dthalf;
        }
      }
    }
  }

  // propagate the volume 1/2 step.
  double vol1 = vol + omega[sd] * dthalf;

  // rescale positions and change box size.
  dilation[sd] = vol1/vol;
  remap(0);

  // propagate particle positions 1 time step.
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }

  // propagate the volume 1/2 step.
  double vol2 = vol1 + omega[sd] * dthalf;

  // rescale positions and change box size.
  dilation[sd] = vol2/vol1;
  remap(0);

  if (kspace_flag) force->kspace->setup();
}

/* ----------------------------------------------------------------------
   2nd half of Verlet update
------------------------------------------------------------------------- */
void FixQBMSST::final_integrate()
{
  int i;

  // v update only for atoms in QBMSST group

  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double vol = compute_vol();
  double p_qbmsst;
  int sd = direction;

  // propagate particle velocities 1/2 step.

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      for (int k = 0; k < 3; k++) {
        double C = (f[i][k] + fran[i][k]) * force->ftm2v / mass[type[i]];//  this term now has a random force part
        double D = mu * omega[sd] * omega[sd] /
          (velocity_sum * mass[type[i]] * vol ) - fric_coef;

        if (k == direction) {
          D = D - 2.0 * omega[sd] / vol;
        }
        if (fabs(dthalf * D) > 1.0e-06) {
          double expd = exp(D * dthalf);
          v[i][k] = expd * ( C + D * v[i][k] - C / expd ) / D;
        } else {
          v[i][k] = v[i][k] + ( C + D * v[i][k] ) * dthalf +
            0.5 * (D * D * v[i][k] + C * D ) * dthalf * dthalf;
        }
      }
    }
  }

  // compute new pressure and volume.

  temperature->compute_vector();

  pressure->compute_vector();
  couple();
  velocity_sum = compute_vsum();
  vol = compute_vol();

  // propagate the time derivative of the volume 1/2 step at fixed V, r, rdot.

  p_qbmsst = nktv2p * mvv2e * velocity * velocity * total_mass *
    ( v0 - vol )/( v0 * v0 );
  double A = total_mass * ( p_current[sd] - p0 - p_qbmsst ) /
    ( qmass * nktv2p * mvv2e );
  double B = total_mass * mu  / ( qmass * vol );

  // prevent blow-up of the volume.

  if (vol > v0 && A > 0.0) {
    A = -A;
  }

  // use taylor expansion to avoid singularity at B == 0.

  if (B * dthalf > 1.0e-06) {
    omega[sd] = ( omega[sd] + A *
                  ( exp(B * dthalf) - 1.0 ) / B ) * exp(-B * dthalf);
  } else {
    omega[sd] = omega[sd] + (A - B * omega[sd]) * dthalf +
      0.5 * (B * B * omega[sd] - A * B ) * dthalf * dthalf;
  }

  // calculate Lagrangian position of computational cell

  lagrangian_position -= velocity*vol/v0*update->dt;

  // trigger virial computation on next timestep

  pressure->addstep(update->ntimestep+1);
  pe->addstep(update->ntimestep+1);
}

/* ----------------------------------------------------------------------
   couple
------------------------------------------------------------------------- */
void FixQBMSST::couple()
{
  double *tensor = pressure->vector;

  p_current[0] = tensor[0];
  p_current[1] = tensor[1];
  p_current[2] = tensor[2];
}

/* ----------------------------------------------------------------------
   change box size
   remap owned or owned+ghost atoms depending on flag
   if rigid bodies exist, scale rigid body centers-of-mass
------------------------------------------------------------------------- */
void FixQBMSST::remap(int flag)
{
  int i,n;
  double oldlo,oldhi,ctr;

  double **v = atom->v;
  if (flag) n = atom->nlocal + atom->nghost;
  else n = atom->nlocal;

  // convert pertinent atoms and rigid bodies to lamda coords

  domain->x2lamda(n);

  if (nrigid)
    for (i = 0; i < nrigid; i++)
      modify->fix[rfix[i]]->deform(0);

  // reset global and local box to new size/shape

  for (i = 0; i < 3; i++) {
    if (direction == i) {
      oldlo = domain->boxlo[i];
      oldhi = domain->boxhi[i];
      ctr = 0.5 * (oldlo + oldhi);
      domain->boxlo[i] = (oldlo-ctr)*dilation[i] + ctr;
      domain->boxhi[i] = (oldhi-ctr)*dilation[i] + ctr;
    }
  }

  domain->set_global_box();
  domain->set_local_box();

  // convert pertinent atoms and rigid bodies back to box coords

  domain->lamda2x(n);

  if (nrigid)
    for (i = 0; i < nrigid; i++)
      modify->fix[rfix[i]]->deform(1);

  for (i = 0; i < n; i++) {
    v[i][direction] = v[i][direction] *
      dilation[direction];
  }
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */
void FixQBMSST::write_restart(FILE *fp)
{
  int n = 0;
  double list[5];
  list[n++] = omega[direction];
  list[n++] = e0;
  list[n++] = v0;
  list[n++] = p0;
  list[n++] = t_current;
  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(&list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */
void FixQBMSST::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  omega[direction] = list[n++];
  e0 = list[n++];
  v0 = list[n++];
  p0 = list[n++];
  t_current = list[n++];
  e0_set = 1;
  v0_set = 1;
  p0_set = 1;
  qtb_set = 1;
}

/* ----------------------------------------------------------------------
   modify parameters
------------------------------------------------------------------------- */
int FixQBMSST::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    id_temp = utils::strdup(arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,"Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning(FLERR,"Temperature for QBMSST is not for group all");

    return 2;

  } else if (strcmp(arg[0],"press") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (pflag) {
      modify->delete_compute(id_press);
      pflag = 0;
    }
    delete [] id_press;
    id_press = utils::strdup(arg[1]);

    int icompute = modify->find_compute(id_press);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify pressure ID");
    pressure = modify->compute[icompute];

    if (pressure->pressflag == 0)
      error->all(FLERR,"Fix_modify pressure ID does not compute pressure");
    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   compute scalar
------------------------------------------------------------------------- */

double FixQBMSST::compute_scalar()
{
  // compute new pressure and volume.

  temperature->compute_vector();
  pressure->compute_vector();
  couple();

  double volume = compute_vol();

  double energy = 0.0;
  int i;

  i = direction;
  energy = qmass * omega[i] * omega[i] / (2.0 * total_mass) * mvv2e;
  energy -= 0.5 * total_mass * velocity * velocity *
    (1.0 - volume/ v0) *
    (1.0 - volume/ v0) * mvv2e;
  energy -= p0 * ( v0 - volume ) / nktv2p;

  return energy;
}

/* ----------------------------------------------------------------------
   return a single element from the following vector,
   [dhug,dray,lgr_vel,lgr_pos,T_qm]
------------------------------------------------------------------------- */
double FixQBMSST::compute_vector(int n)
{
  if (n == 0) {
    return compute_hugoniot();
  } else if (n == 1) {
    return compute_rayleigh();
  } else if (n == 2) {
    return compute_lagrangian_speed();
  } else if (n == 3) {
    return compute_lagrangian_position();
  } else if (n == 4) {
    return t_current;
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   Computes the deviation of the current point
   from the Hugoniot in Kelvin for the QBMSST.
------------------------------------------------------------------------- */
double FixQBMSST::compute_hugoniot()
{
  double v, e, p;
  double dhugo;

  e = compute_etotal();

  temperature->compute_vector();
  pressure->compute_vector();
  p = pressure->vector[direction];

  v = compute_vol();

  dhugo = (0.5 * (p + p0 ) * ( v0 - v)) /
    force->nktv2p + e0 - e;
  dhugo /= temperature->dof * force->boltz;

  return dhugo;
}

/* ----------------------------------------------------------------------
   Computes the deviation of the current point from the Rayleigh
   in pressure units for the QBMSST.
------------------------------------------------------------------------- */
double FixQBMSST::compute_rayleigh()
{
  double v, p;
  double drayleigh;

  temperature->compute_vector();
  pressure->compute_vector();
  p = pressure->vector[direction];

  v = compute_vol();

  drayleigh = p - p0 -
    total_mass * velocity * velocity * force->mvv2e *
    (1.0 - v / v0 ) * force->nktv2p / v0;

  return drayleigh;
}

/* ----------------------------------------------------------------------
   Computes the speed of the QBMSST computational cell in the
   unshocked material rest-frame
------------------------------------------------------------------------- */
double FixQBMSST::compute_lagrangian_speed()
{
  double v = compute_vol();
  return velocity*(1.0-v/v0);
}

 /* ----------------------------------------------------------------------
    Computes the distance behind the
    shock front of the QBMSST computational cell.
 ------------------------------------------------------------------------- */
double FixQBMSST::compute_lagrangian_position()
{
   return lagrangian_position;
}

/* ----------------------------------------------------------------------
   Computes the atomic kinetic + atomic potential energy.  This excludes the QBMSST
   external potential terms in the QBMSST Lagrangian.
------------------------------------------------------------------------- */
double FixQBMSST::compute_etotal()
{
  double epot,ekin,etot;
  epot = pe->compute_scalar();
  ekin = temperature->compute_scalar();
  ekin *= 0.5 * temperature->dof * force->boltz;
  etot = epot+ekin;
  return etot;
}

/* ----------------------------------------------------------------------
   Computes the atomic kinetic + atomic potential energy + QBMSST external potential.
------------------------------------------------------------------------- */
double FixQBMSST::compute_egrand()
{
  double epot,ekin,ecouple,etot;
  epot = pe->compute_scalar();
  ekin = temperature->compute_scalar();
  ekin *= 0.5 * temperature->dof * force->boltz;
  ecouple = compute_scalar();
  etot = epot + ekin + ecouple;
  return etot;
}

/* ----------------------------------------------------------------------
   Computes the atomic kinetic + atomic potential energy.  This excludes the QBMSST
   external potential terms in the QBMSST Lagrangian.
------------------------------------------------------------------------- */
double FixQBMSST::compute_vol()
{
  if (domain->dimension == 3)
    return domain->xprd * domain->yprd * domain->zprd;
  else
    return domain->xprd * domain->yprd;
}

/* ----------------------------------------------------------------------
   Checks to see if the allocated size of old_velocity is >= n
   The number of local atoms can change during a parallel run.
------------------------------------------------------------------------- */
void FixQBMSST::check_alloc(int n)
{
  if (atoms_allocated < n) {
    memory->destroy(old_velocity);
    memory->create(old_velocity,n,3,"qbmsst:old_velocity");
    atoms_allocated = n;
  }
}

/* ----------------------------------------------------------------------
   Computes the atomic kinetic + atomic potential energy.  This excludes the QBMSST
   external potential terms in the QBMSST Lagrangian.
------------------------------------------------------------------------- */
double FixQBMSST::compute_vsum()
{
  double vsum;

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) ;
    }
  }

  MPI_Allreduce(&t,&vsum,1,MPI_DOUBLE,MPI_SUM,world);
  return vsum;
}

/* ----------------------------------------------------------------------
   memory usage of fix qbmsst
------------------------------------------------------------------------- */
double FixQBMSST::memory_usage()
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
void FixQBMSST::grow_arrays(int nmax)
{
  memory->grow(random_array_0,nmax,2*N_f,"qbmsst:random_array_0");
  memory->grow(random_array_1,nmax,2*N_f,"qbmsst:random_array_1");
  memory->grow(random_array_2,nmax,2*N_f,"qbmsst:random_array_2");
  memory->grow(fran,nmax,3,"qbmsst:fran");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */
void FixQBMSST::copy_arrays(int i, int j, int /*delflag*/)
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
int FixQBMSST::pack_exchange(int i, double *buf)
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
int FixQBMSST::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  for (int m = 0; m < 2*N_f; m++) random_array_0[nlocal][m] = buf[n++];
  for (int m = 0; m < 2*N_f; m++) random_array_1[nlocal][m] = buf[n++];
  for (int m = 0; m < 2*N_f; m++) random_array_2[nlocal][m] = buf[n++];
  for (int m = 0; m < 3; m++) fran[nlocal][m] = buf[n++];
  return n;
}
