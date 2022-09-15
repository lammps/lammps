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
   Contributing authors: Laurence Fried (LLNL), Evan Reed (LLNL, Stanford)
     implementation of the Multi-Scale Shock Method
   see Reed, Fried, Joannopoulos, Phys Rev Lett, 90, 235503 (2003)
------------------------------------------------------------------------- */

#include "fix_msst.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_external.h"
#include "force.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMSST::FixMSST(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), old_velocity(nullptr), rfix(nullptr),
  id_temp(nullptr), id_press(nullptr), id_pe(nullptr), temperature(nullptr),
  pressure(nullptr), pe(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal fix msst command");

  restart_global = 1;
  time_integrate = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extscalar = 1;
  extvector = 0;
  ecouple_flag = 1;

  // set defaults

  velocity = 0.0;
  dilation[0] = dilation[1] = dilation[2] = 1.0;
  p0 = 0.0;
  v0 = 1.0;
  e0 = 0.0;
  TS_int = 0;
  T0S0 = 0.0;
  S_elec = 0.0;
  S_elec_1 = 0.0;
  S_elec_2 = 0.0;

  qmass = 1.0e1;
  mu = 0.0;
  direction = 2;
  p0_set = 0;
  v0_set = 0;
  e0_set = 0;
  tscale = 0.01;
  dftb = 0;
  beta = 0.0;

  if (strcmp(arg[3],"x") == 0) {
    direction = 0;
    box_change |= BOX_CHANGE_X;
  } else if (strcmp(arg[3],"y") == 0) {
    direction = 1;
    box_change |= BOX_CHANGE_Y;
  } else if (strcmp(arg[3],"z") == 0) {
    direction = 2;
    box_change |= BOX_CHANGE_Z;
  } else error->all(FLERR,"Illegal fix msst command");

  velocity = utils::numeric(FLERR,arg[4],false,lmp);
  if (velocity < 0) error->all(FLERR,"Illegal fix msst command");

  // optional args

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"q") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix msst command");
      qmass = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"mu") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix msst command");
      mu = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"p0") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix msst command");
      p0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p0_set = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"v0") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix msst command");
      v0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      v0_set = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"e0") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix msst command");
      e0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      e0_set = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"tscale") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix msst command");
      tscale = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (tscale < 0.0 || tscale > 1.0)
        error->all(FLERR,"Fix msst tscale must satisfy 0 <= tscale < 1");
      iarg += 2;
    } else if (strcmp(arg[iarg],"dftb") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix msst command");
      dftb = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"beta") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix msst command");
      beta = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (beta < 0.0 || beta > 1.0)
        error->all(FLERR,"Illegal fix msst command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix msst command");
  }

  // output MSST info

  if (comm->me == 0) {
    std::string mesg = "MSST parameters:\n";
    if (direction == 0) mesg += "  Shock in x direction\n";
    else if (direction == 1) mesg += "  Shock in y direction\n";
    else if (direction == 2) mesg += "  Shock in z direction\n";
    mesg += fmt::format("  Cell mass-like parameter qmass "
                        "(units of mass^2/length^4) = {:.8g}\n", qmass);
    mesg += fmt::format("  Shock velocity = {:.8g}\n", velocity);
    mesg += fmt::format("  Artificial viscosity (units of mass/length/time) = {:.8g}\n", mu);

    if (p0_set)
      mesg += fmt::format("  Initial pressure specified to be {:.8g}\n", p0);
    else mesg += "  Initial pressure calculated on first step\n";

    if (v0_set)
      mesg += fmt::format("  Initial volume specified to be {:.8g}\n", v0);
    else mesg += "  Initial volume calculated on first step\n";

    if (e0_set)
      mesg += fmt::format("  Initial energy specified to be {:.8g}\n", e0);
    else mesg += "  Initial energy calculated on first step\n";
    utils::logmesg(lmp,mesg);
  }

  // check for periodicity in controlled dimensions

  if (domain->nonperiodic) error->all(FLERR,"Fix msst requires a periodic box");

  // create a new temperature compute
  // id = fix-ID + "MSST_temp"
  // compute group = all since pressure is always global (group all)
  //   and thus its KE/temperature contribution should use group all

  std::string fixcmd = std::string(id) + "MSST_temp";
  id_temp = utils::strdup(fixcmd);
  modify->add_compute(fixcmd + " all temp");
  tflag = 1;

  // create a new pressure compute
  // id = fix-ID + "MSST_press", compute group = all
  // pass id_temp as 4th arg to pressure constructor

  fixcmd = std::string(id) + "MSST_press";
  id_press = utils::strdup(fixcmd);
  modify->add_compute(fixcmd + " all pressure " + std::string(id_temp));
  pflag = 1;

  // create a new potential energy compute
  // id = fix-ID + "MSST_pe", compute group = all

  fixcmd = std::string(id) + "MSST_pe";
  id_pe = utils::strdup(fixcmd);
  modify->add_compute(fixcmd + " all pe");
  peflag = 1;

  // initialize the time derivative of the volume

  omega[0] = omega[1] = omega[2] = 0.0;
  nrigid = 0;
  rfix = nullptr;

  maxold = -1;
  old_velocity = nullptr;
}

/* ---------------------------------------------------------------------- */

FixMSST::~FixMSST()
{
  delete [] rfix;

  // delete temperature and pressure if fix created them

  if (tflag) modify->delete_compute(id_temp);
  if (pflag) modify->delete_compute(id_press);
  if (peflag) modify->delete_compute(id_pe);

  delete [] id_temp;
  delete [] id_press;
  delete [] id_pe;

  memory->destroy(old_velocity);
}

/* ---------------------------------------------------------------------- */

int FixMSST::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMSST::init()
{
  if (atom->mass == nullptr)
    error->all(FLERR,"Cannot use fix msst without per-type mass defined");

  // set compute ptrs

  int itemp = modify->find_compute(id_temp);
  int ipress = modify->find_compute(id_press);
  int ipe = modify->find_compute(id_pe);
  if (itemp < 0 || ipress < 0|| ipe < 0)
    error->all(FLERR,"Could not find fix msst compute ID");
  if (modify->compute[itemp]->tempflag == 0)
    error->all(FLERR,"Fix msst compute ID does not compute temperature");
  if (modify->compute[ipress]->pressflag == 0)
    error->all(FLERR,"Fix msst compute ID does not compute pressure");
  if (modify->compute[ipe]->peflag == 0)
    error->all(FLERR,"Fix msst compute ID does not compute potential energy");

  temperature = modify->compute[itemp];
  pressure = modify->compute[ipress];
  pe = modify->compute[ipe];

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dthalf = 0.5 * update->dt;

  boltz = force->boltz;
  nktv2p = force->nktv2p;
  mvv2e = force->mvv2e;

  double mass = 0.0;
  for (int i = 0; i < atom->nlocal; i++) mass += atom->mass[atom->type[i]];
  MPI_Allreduce(&mass,&total_mass,1,MPI_DOUBLE,MPI_SUM,world);

  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  // detect if any fix rigid exist so rigid bodies move when box is dilated
  // rfix[] = indices to each fix rigid

  delete [] rfix;
  nrigid = 0;
  rfix = nullptr;

  for (int i = 0; i < modify->nfix; i++)
    if (utils::strmatch(modify->fix[i]->style,"^rigid")  ||
        utils::strmatch(modify->fix[i]->style,"^poems$")) nrigid++;
  if (nrigid) {
    rfix = new int[nrigid];
    nrigid = 0;
    for (int i = 0; i < modify->nfix; i++)
      if (utils::strmatch(modify->fix[i]->style,"^rigid") ||
          utils::strmatch(modify->fix[i]->style,"^poems$")) rfix[nrigid++] = i;
  }

  // find fix external being used to drive LAMMPS from DFTB+

  if (dftb) {
    for (int i = 0; i < modify->nfix; i++)
      if (utils::strmatch(modify->fix[i]->style,"^external$"))
        fix_external = dynamic_cast<FixExternal *>(modify->fix[i]);
    if (fix_external == nullptr)
      error->all(FLERR,"Fix msst dftb cannot be used w/out fix external");
  }
}

/* ----------------------------------------------------------------------
   compute T,P before integrator starts
------------------------------------------------------------------------- */

void FixMSST::setup(int /*vflag*/)
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
      utils::logmesg(lmp,"Fix MSST v0 = {:.8g}\n", v0);
  }

  if (p0_set == 0) {
    p0 = p_current[direction];
    p0_set = 1;

    if (comm->me == 0)
      utils::logmesg(lmp,"Fix MSST p0 = {:.8g}\n", p0);
  }

  if (e0_set == 0) {
    e0 = compute_etotal();
    e0_set = 1;

    if (comm->me == 0)
      utils::logmesg(lmp,"Fix MSST e0 = {:.8g}\n", e0);
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
      utils::logmesg(lmp,"Fix MSST initial strain rate of {:.8g} "
                     "established by reducing temperature by factor "
                     "of {:.8g}\n",fac2,tscale);
    for (int i = 0; i < atom->nlocal; i++) {
      if (mask[i] & groupbit) {
        for (int k = 0; k < 3; k++) {
          v[i][k]*=sqrt_initial_temperature_scaling;
        }
      }
    }
  }

  // trigger virial computation on next timestep

  pe->addstep(update->ntimestep+1);
  pressure->addstep(update->ntimestep+1);
}

/* ----------------------------------------------------------------------
   1st half of Verlet update
------------------------------------------------------------------------- */

void FixMSST::initial_integrate(int /*vflag*/)
{
  int i,k;
  double p_msst;                       // MSST driving pressure
  double vol;

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  double **x = atom->x;
  int sd = direction;

  // realloc old_velocity if necessary

  if (nlocal > maxold) {
    memory->destroy(old_velocity);
    maxold = atom->nmax;
    memory->create(old_velocity,maxold,3,"msst:old_velocity");
  }

  // for DFTB, extract TS_dftb from fix external
  // must convert energy to mv^2 units

  if (dftb) {
    const double TS_dftb = fix_external->compute_vector(0);
    const double TS = force->ftm2v*TS_dftb;
    // update S_elec terms and compute TS_dot via finite differences
    S_elec_2 = S_elec_1;
    S_elec_1 = S_elec;
    const double Temp = temperature->compute_scalar();
    S_elec = TS/Temp;
    TS_dot = Temp*(3.0*S_elec-4.0*S_elec_1+S_elec_2)/(2.0*update->dt);
    TS_int += (update->dt*TS_dot);
    if (update->ntimestep == 1) T0S0 = TS;
  }

  // compute new pressure and volume

  temperature->compute_vector();
  pressure->compute_vector();
  couple();
  vol = compute_vol();

  // compute etot + extra terms for conserved quantity

  double e_scale = compute_etotal() + compute_scalar();

  // propagate the time derivative of
  // the volume 1/2 step at fixed vol, r, rdot

  p_msst = nktv2p * mvv2e * velocity * velocity * total_mass *
    ( v0 - vol)/( v0 * v0);
  double A = total_mass * ( p_current[sd] - p0 - p_msst ) /
    (qmass * nktv2p * mvv2e);
  double B = total_mass * mu / ( qmass * vol );

  // prevent blow-up of the volume

  if (vol > v0 && A > 0.0) A = -A;

  // use Taylor expansion to avoid singularity at B = 0

  if (B * dthalf > 1.0e-06) {
    omega[sd] = ( omega[sd] + A * ( exp(B * dthalf) - 1.0 ) / B )
      * exp(-B * dthalf);
  } else {
    omega[sd] = omega[sd] + (A - B * omega[sd]) * dthalf +
      0.5 * (B * B * omega[sd] - A * B ) * dthalf * dthalf;
  }

  // propagate velocity sum 1/2 step by
  // temporarily propagating the velocities

  velocity_sum = compute_vsum();

  if (dftb) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        for (k = 0; k < 3; k++) {
          const double C = f[i][k] * force->ftm2v / mass[type[i]];
          const double TS_term = TS_dot/(mass[type[i]]*velocity_sum);
          const double escale_term = force->ftm2v*beta*(e0-e_scale) /
            (mass[type[i]]*velocity_sum);
          double D = mu * omega[sd] * omega[sd] /
            (velocity_sum * mass[type[i]] * vol );
          D += escale_term - TS_term;
          old_velocity[i][k] = v[i][k];
          if (k == direction) D -= 2.0 * omega[sd] / vol;
          if (fabs(dthalf * D) > 1.0e-06) {
            const double expd = exp(D * dthalf);
            v[i][k] = expd * ( C + D * v[i][k] - C / expd ) / D;
          } else {
            v[i][k] = v[i][k] + ( C + D * v[i][k] ) * dthalf +
              0.5 * (D * D * v[i][k] + C * D ) * dthalf * dthalf;
          }
        }
      }
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        for (k = 0; k < 3; k++) {
          const double C = f[i][k] * force->ftm2v / mass[type[i]];
          double D = mu * omega[sd] * omega[sd] /
            (velocity_sum * mass[type[i]] * vol );
          old_velocity[i][k] = v[i][k];
          if (k == direction) {
            D -= 2.0 * omega[sd] / vol;
          }
          if (fabs(dthalf * D) > 1.0e-06) {
            const double expd = exp(D * dthalf);
            v[i][k] = expd * ( C + D * v[i][k] - C / expd ) / D;
          } else {
            v[i][k] = v[i][k] + ( C + D * v[i][k] ) * dthalf +
              0.5 * (D * D * v[i][k] + C * D ) * dthalf * dthalf;
          }
        }
      }
    }
  }

  velocity_sum = compute_vsum();

  // reset the velocities

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      v[i][0] = old_velocity[i][0];
      v[i][1] = old_velocity[i][1];
      v[i][2] = old_velocity[i][2];
    }
  }

  // propagate velocities 1/2 step using the new velocity sum

  if (dftb) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        for (k = 0; k < 3; k++) {
          const double C = f[i][k] * force->ftm2v / mass[type[i]];
          const double TS_term = TS_dot/(mass[type[i]]*velocity_sum);
          const double escale_term = force->ftm2v*beta*(e0-e_scale) /
            (mass[type[i]]*velocity_sum);
          double D = mu * omega[sd] * omega[sd] /
            (velocity_sum * mass[type[i]] * vol );
          D += escale_term - TS_term;
          if (k == direction) D -= 2.0 * omega[sd] / vol;
          if (fabs(dthalf * D) > 1.0e-06) {
            const double expd = exp(D * dthalf);
            v[i][k] = expd * ( C + D * v[i][k] - C / expd ) / D;
          } else {
            v[i][k] = v[i][k] + ( C + D * v[i][k] ) * dthalf +
              0.5 * (D * D * v[i][k] + C * D ) * dthalf * dthalf;
          }
        }
      }
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        for (k = 0; k < 3; k++) {
          const double C = f[i][k] * force->ftm2v / mass[type[i]];
          double D = mu * omega[sd] * omega[sd] /
            (velocity_sum * mass[type[i]] * vol );
          if (k == direction) {
            D -= 2.0 * omega[sd] / vol;
          }
          if (fabs(dthalf * D) > 1.0e-06) {
            const double expd = exp(D * dthalf);
            v[i][k] = expd * ( C + D * v[i][k] - C / expd ) / D;
          } else {
            v[i][k] = v[i][k] + ( C + D * v[i][k] ) * dthalf +
              0.5 * (D * D * v[i][k] + C * D ) * dthalf * dthalf;
          }
        }
      }
    }
  }

  // propagate the volume 1/2 step

  double vol1 = vol + omega[sd] * dthalf;

  // rescale positions and change box size

  dilation[sd] = vol1/vol;
  remap(0);

  // propagate particle positions 1 time step

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }

  // propagate the volume 1/2 step

  double vol2 = vol1 + omega[sd] * dthalf;

  // rescale positions and change box size

  dilation[sd] = vol2/vol1;
  remap(0);

  if (kspace_flag) force->kspace->setup();
}

/* ----------------------------------------------------------------------
   2nd half of Verlet update
------------------------------------------------------------------------- */

void FixMSST::final_integrate()
{
  int i;
  double p_msst;                  // MSST driving pressure

  // v update only for atoms in MSST group

  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double vol = compute_vol();

  int sd = direction;

  // compute etot + extra terms for conserved quantity

  double e_scale = compute_etotal() + compute_scalar();

  // for DFTB, extract TS_dftb from fix external
  // must convert energy to mv^2 units

  if (dftb) {
    const double TS_dftb = fix_external->compute_vector(0);
    const double TS = force->ftm2v*TS_dftb;
    S_elec_2 = S_elec_1;
    S_elec_1 = S_elec;
    const double Temp = temperature->compute_scalar();
    // update S_elec terms and compute TS_dot via finite differences
    S_elec = TS/Temp;
    TS_dot = Temp*(3.0*S_elec-4.0*S_elec_1+S_elec_2)/(2.0*update->dt);
    TS_int += (update->dt*TS_dot);
    if (update->ntimestep == 1) T0S0 = TS;
  }

  // propagate particle velocities 1/2 step

  if (dftb) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        for (int k = 0; k < 3; k++) {
          const double C = f[i][k] * force->ftm2v / mass[type[i]];
          const double TS_term = TS_dot/(mass[type[i]]*velocity_sum);
          const double escale_term = force->ftm2v*beta*(e0-e_scale) /
            (mass[type[i]]*velocity_sum);
          double D = mu * omega[sd] * omega[sd] /
            (velocity_sum * mass[type[i]] * vol );
          D += escale_term - TS_term;
          if (k == direction) D -= 2.0 * omega[sd] / vol;
          if (fabs(dthalf * D) > 1.0e-06) {
            const double expd = exp(D * dthalf);
            v[i][k] = expd * ( C + D * v[i][k] - C / expd ) / D;
          } else {
            v[i][k] = v[i][k] + ( C + D * v[i][k] ) * dthalf +
              0.5 * (D * D * v[i][k] + C * D ) * dthalf * dthalf;
          }
        }
      }
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        for (int k = 0; k < 3; k++) {
          const double C = f[i][k] * force->ftm2v / mass[type[i]];
          double D = mu * omega[sd] * omega[sd] /
            (velocity_sum * mass[type[i]] * vol );
          if (k == direction) {
            D -= 2.0 * omega[sd] / vol;
          }
          if (fabs(dthalf * D) > 1.0e-06) {
            const double expd = exp(D * dthalf);
            v[i][k] = expd * ( C + D * v[i][k] - C / expd ) / D;
          } else {
            v[i][k] = v[i][k] + ( C + D * v[i][k] ) * dthalf +
              0.5 * (D * D * v[i][k] + C * D ) * dthalf * dthalf;
          }
        }
      }
    }
  }

  // compute new pressure and volume

  temperature->compute_vector();

  pressure->compute_vector();
  couple();
  velocity_sum = compute_vsum();
  vol = compute_vol();

  // propagate the time derivative of the volume 1/2 step at fixed V, r, rdot

  p_msst = nktv2p * mvv2e * velocity * velocity * total_mass *
    ( v0 - vol )/( v0 * v0 );
  double A = total_mass * ( p_current[sd] - p0 - p_msst ) /
    ( qmass * nktv2p * mvv2e );
  const double B = total_mass * mu  / ( qmass * vol );

  // prevent blow-up of the volume

  if (vol > v0 && A > 0.0) A = -A;

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

  // trigger energy and virial computation on next timestep

  pe->addstep(update->ntimestep+1);
  pressure->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixMSST::couple()
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

void FixMSST::remap(int flag)
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

void FixMSST::write_restart(FILE *fp)
{
  int n = 0;
  double list[5];
  list[n++] = omega[direction];
  list[n++] = e0;
  list[n++] = v0;
  list[n++] = p0;
  list[n++] = TS_int;
  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(&list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixMSST::restart(char *buf)
{
  int n = 0;
  auto list = (double *) buf;
  omega[direction] = list[n++];
  e0 = list[n++];
  v0 = list[n++];
  p0 = list[n++];
  TS_int = list[n++];
  tscale = 0.0;           // set tscale to zero for restart
  p0_set = 1;
  v0_set = 1;
  e0_set = 1;
}

/* ---------------------------------------------------------------------- */

int FixMSST::modify_param(int narg, char **arg)
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
    if (icompute < 0)
      error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,"Fix_modify temperature ID does not "
                 "compute temperature");
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning(FLERR,"Temperature for MSST is not for group all");

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

/* ---------------------------------------------------------------------- */

double FixMSST::compute_scalar()
{
  // compute new pressure and volume

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

  // subtract off precomputed TS_int integral value
  // TS_int = 0 for non DFTB calculations

  if (dftb) energy -= TS_int;

  return energy;
}

/* ----------------------------------------------------------------------
   return a single element from the following vector,
   [dhug,dray,lgr_vel,lgr_pos]
------------------------------------------------------------------------- */

double FixMSST::compute_vector(int n)
{
  if (n == 0) {
    return compute_hugoniot();
  } else if (n == 1) {
    return compute_rayleigh();
  } else if (n == 2) {
    return compute_lagrangian_speed();
  } else if (n == 3) {
    return compute_lagrangian_position();
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   Computes the deviation of the current point
   from the Hugoniot in Kelvin for the MSST
------------------------------------------------------------------------- */

double FixMSST::compute_hugoniot()
{
  if (!temperature) return 0.0;
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
   in pressure units for the MSST
------------------------------------------------------------------------- */

double FixMSST::compute_rayleigh()
{
  if (!temperature) return 0.0;

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
   Computes the speed of the MSST computational cell in the
   unshocked material rest-frame
------------------------------------------------------------------------- */

double FixMSST::compute_lagrangian_speed()
{
  double v = compute_vol();
  return velocity*(1.0-v/v0);
}

/* ----------------------------------------------------------------------
   Computes the distance behind the
   shock front of the MSST computational cell.
------------------------------------------------------------------------- */

double FixMSST::compute_lagrangian_position()
{
   return lagrangian_position;
}

/* ---------------------------------------------------------------------- */

double FixMSST::compute_etotal()
{
  if (!pe) return 0.0;

  double epot,ekin,etot;
  epot = pe->compute_scalar();
  ekin = temperature->compute_scalar();
  ekin *= 0.5 * temperature->dof * force->boltz;
  etot = epot+ekin;
  return etot;
}

/* ---------------------------------------------------------------------- */

double FixMSST::compute_vol()
{
  if (domain->dimension == 3)
    return domain->xprd * domain->yprd * domain->zprd;
  else
    return domain->xprd * domain->yprd;
}

/* ---------------------------------------------------------------------- */

double FixMSST::compute_vsum()
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
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixMSST::memory_usage()
{
  double bytes = 3*atom->nmax * sizeof(double);
  return bytes;
}
