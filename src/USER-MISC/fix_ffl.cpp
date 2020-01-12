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
   Fast-forward Langevin thermostat, see
   M. Hijazi, D. M. Wilkins, M. Ceriotti, J. Chem. Phys. 148, 184109 (2018)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Lionel Constantin (EPFL), David M. Wilkins (EPFL),
                         Michele Ceriotti (EPFL)
------------------------------------------------------------------------- */

#include "fix_ffl.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "comm.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum {NOBIAS,BIAS};
enum {CONSTANT,EQUAL,ATOM};
enum {NO_FLIP, FLIP_RESCALE, FLIP_HARD, FLIP_SOFT};
//#define FFL_DEBUG 1

#define MAXLINE 1024

/* syntax for fix_ffl:
 * fix nfix id-group ffl tau Tstart Tstop seed [flip_type]
 *
 *                                                                        */

/* ---------------------------------------------------------------------- */


FixFFL::FixFFL(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {


  if (narg < 7)
    error->all(FLERR,"Illegal fix ffl command. Expecting: fix <fix-ID>"
               " <group-ID> ffl <tau> <Tstart> <Tstop> <seed>  ");

  restart_peratom = 1;
  time_integrate = 1;
  scalar_flag = 1;

  //gamma = 1/ time constant(tau)
  if (force->numeric(FLERR,arg[3]) <= 0)
    error->all(FLERR,"Illegal fix ffl tau value, should be greater than 0");
  gamma = 1.0/force->numeric(FLERR,arg[3]);
  ffl_every=1;
  ffl_step=0;

  // start temperature (t ramp)
  t_start = force->numeric(FLERR,arg[4]);

  // final temperature (t ramp)
  t_stop = force->numeric(FLERR,arg[5]);

  // PRNG seed
  int seed = force->inumeric(FLERR,arg[6]);

  // Flip type used, uses rescale if no flip is given
  if (narg == 8) {
    if (strcmp(arg[7],"no_flip") == 0) {
      flip_int = NO_FLIP;
    } else if (strcmp(arg[7],"rescale") == 0) {
      flip_int = FLIP_RESCALE;
    } else if (strcmp(arg[7],"hard") == 0) {
      flip_int = FLIP_HARD;
    } else if (strcmp(arg[7],"soft") == 0) {
      flip_int = FLIP_SOFT;
    } else {
      error->all(FLERR,"Illegal fix ffl flip type, only accepts : rescale - hard - soft - no_flip");
    }
  } else {
    flip_int = FLIP_RESCALE;
  }

  t_target=t_start;

  // initialize Marsaglia RNG with processor-unique seed
  // NB: this means runs will not be the same with different numbers of processors
  if (seed <= 0) error->all(FLERR,"Illegal fix ffl command");
  random = new RanMars(lmp,seed + comm->me);

  // allocate per-type arrays for mass-scaling
  sqrt_m=NULL;
  memory->grow(sqrt_m, atom->ntypes+1,"ffl:sqrt_m");

  // allocates space for temporaries
  ffl_tmp1=ffl_tmp2=NULL;

  grow_arrays(atom->nmax);

  // add callbacks to enable restarts
  atom->add_callback(0);
  atom->add_callback(1);

  energy = 0.0;
}


/* --- Frees up memory used by temporaries and buffers ------------------ */

FixFFL::~FixFFL() {
  delete random;

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  memory->destroy(sqrt_m);
  memory->destroy(ffl_tmp1);
  memory->destroy(ffl_tmp2);
}

/* ---------------------------------------------------------------------- */

int FixFFL::setmask() {
  int mask = 0;

  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  mask |= THERMO_ENERGY;


  return mask;
}

/* ------- Initializes one-time quantities for FFL ---------------------- */

void FixFFL::init() {
  doffl = 1;
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  // set force prefactors
  if (!atom->rmass) {
    for (int i = 1; i <= atom->ntypes; i++) {
      sqrt_m[i] = sqrt(atom->mass[i]);
    }
  }

  if (strstr(update->integrate_style,"respa")) {
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
    step_respa = ((Respa *) update->integrate)->step;
  }

  init_ffl();
}

/* ------- Initializes constants for FFL (change with T and dt) ------- */

void FixFFL::init_ffl() {
  const double kT = t_target * force->boltz / force->mvv2e;

  // compute constants for FFL

  c1 = exp ( - gamma * 0.5 * dtv );
  c2 = sqrt( (1.0 - c1*c1)* kT ); //without the mass term


}



/* ---------------------------------------------------------------------- */

void FixFFL::setup(int vflag) {
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

void FixFFL::ffl_integrate() {
  double **v = atom->v;
  double *rmass = atom->rmass, smi, ismi;
  double factor;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // loads momentum data (mass-scaled) into the temporary vectors for the propagation
  int nk=0;
  double deltae=0.0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass) smi = sqrt(rmass[i]);
      else smi = sqrt_m[type[i]];

      for (int k = 0; k<3; k++) {
        // first loads velocities and accumulates conserved quantity
        ffl_tmp2[nk] = v[i][k] * smi;
        deltae += ffl_tmp2[nk] * ffl_tmp2[nk];
        nk++;
      }
    }
  }

  //fills up a vector of random numbers
  for (int i = 0; i < nk; i++) ffl_tmp1[i] = random->gaussian();


  // unloads momentum data (mass-scaled) from the temporary vectors
  nk=0;
  for (int i = 0; i < nlocal; i++) if (mask[i] & groupbit) {
      if (rmass) ismi = 1.0 / sqrt(rmass[i]);
      else ismi = 1.0/ sqrt_m[type[i]];

      for (int k = 0; k<3; k++) {
        // fetches new velocities and completes computation of the conserved quantity change
        v[i][k]= c1*v[i][k] + c2*ffl_tmp1[nk]*ismi;

        deltae-= v[i][k]*v[i][k] /ismi /ismi;

        //flips the sign of the momentum (HARD FLIP)
        if ( flip_int == FLIP_HARD) {
          if (v[i][k]*ffl_tmp2[nk] < 0.0) v[i][k] = -v[i][k];
        }

        nk++;
      }
    }

  //rescale operation (RESCALE FLIP)
  if (flip_int == FLIP_RESCALE) {
    nk=0;
    for (int i = 0; i < nlocal; i++) if (mask[i] & groupbit) {
      factor = sqrt ((v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) /
                       (ffl_tmp2[nk]*ffl_tmp2[nk] + ffl_tmp2[nk+1]*ffl_tmp2[nk+1]
                        + ffl_tmp2[nk+2]*ffl_tmp2[nk+2]));

      for (int k = 0; k<3; k++) {
        v[i][k]= factor * ffl_tmp2[nk];
        nk++;
      }
    }
  }


  //soft flip operation (SOFT FLIP)
  if (flip_int == FLIP_SOFT) {
    nk=0;
    for (int i = 0; i < nlocal; i++) if (mask[i] & groupbit) {
      factor = v[i][0]*ffl_tmp2[nk] + v[i][1]*ffl_tmp2[nk+1] + v[i][2]*ffl_tmp2[nk+2];
      if (factor < 0) {
        factor =  factor / (ffl_tmp2[nk]*ffl_tmp2[nk] + ffl_tmp2[nk+1]*ffl_tmp2[nk+1]
                            + ffl_tmp2[nk+2]*ffl_tmp2[nk+2]);

        for (int k = 0; k<3; k++) {
          v[i][k] -= 2.0 * factor * ffl_tmp2[nk];
          nk++;
        }
      } else {
        nk += 3;
      }
    }

  }

  energy += deltae*0.5*force->mvv2e;

}

void FixFFL::initial_integrate(int /* vflag */) {
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

  ffl_step--;
  if (doffl && ffl_step<1) ffl_integrate();

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }
}

void FixFFL::final_integrate() {
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

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
  }

  if (doffl && ffl_step<1) {
    ffl_integrate();
    ffl_step = ffl_every;
  }

  // Change the temperature for the next step
  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  t_target = t_start + delta * (t_stop - t_start);
  if (t_stop != t_start) {
    // only updates if it is really necessary
    init_ffl();
  }

}
/* ---------------------------------------------------------------------- */

void FixFFL::initial_integrate_respa(int vflag, int ilevel, int /* iloop */) {
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel==nlevels_respa-1) ffl_integrate();
  doffl=0;
  if (ilevel == 0) initial_integrate(vflag);
  else {
    final_integrate();
  }
}

void FixFFL::final_integrate_respa(int ilevel, int /* iloop */) {

  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  doffl=0;
  final_integrate();
  if (ilevel==nlevels_respa-1) ffl_integrate();
}


double FixFFL::compute_scalar() {

  double energy_me = energy;
  double energy_all;
  MPI_Allreduce(&energy_me,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);

  return energy_all;
}

/* ----------------------------------------------------------------------
   extract thermostat properties
------------------------------------------------------------------------- */

void *FixFFL::extract(const char *str, int &dim) {
  dim = 0;
  if (strcmp(str,"t_target") == 0) {
    return &t_target;
  }
  return NULL;
}


/* ----------------------------------------------------------------------
   Called when a change to the target temperature is requested mid-run
------------------------------------------------------------------------- */

void FixFFL::reset_target(double t_new) {

  t_target = t_start = t_stop = t_new;
}

/* ----------------------------------------------------------------------
   Called when a change to the timestep is requested mid-run
------------------------------------------------------------------------- */

void FixFFL::reset_dt() {
  // set the time integration constants
  dtv = update->dt;
  dtf = 0.5 * update->dt * (force->ftm2v);
  init_ffl();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixFFL::memory_usage() {
  double bytes = atom->nmax*(3*2)*sizeof(double);
  return bytes;
}


/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixFFL::grow_arrays(int nmax) {
  memory->grow(ffl_tmp1, nmax*3,"ffl:tmp1");
  memory->grow(ffl_tmp2, nmax*3,"ffl:tmp2");
  //zeroes out temporary buffers
  for (int i=0; i< nmax*3; ++i) ffl_tmp1[i] = 0.0;
  for (int i=0; i< nmax*3; ++i) ffl_tmp2[i] = 0.0;
}


