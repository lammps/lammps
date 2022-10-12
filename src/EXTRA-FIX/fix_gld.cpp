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
   Contributing authors: Stephen Bond (SNL) and
                         Andrew Baczewski (Michigan State/SNL)
------------------------------------------------------------------------- */

#include "fix_gld.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "random_mars.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

#define GLD_UNIFORM_DISTRO

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   Parses parameters passed to the method, allocates some memory
------------------------------------------------------------------------- */

FixGLD::FixGLD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  step_respa(nullptr), prony_c(nullptr), prony_tau(nullptr), s_gld(nullptr), random(nullptr)
{
  int narg_min = 8;
  // Check to make sure we have the minimal number of inputs
  if (narg < narg_min) error->all(FLERR,"Illegal fix gld command");

  time_integrate = 1;
  restart_peratom = 1;

  // Parse the first set of required input arguments
  // 0 = Fix ID           (e.g., 1)
  // 1 = Group ID         (e.g., all)
  // 2 = gld              (name of this fix)
  // 3 = t_start          (Starting target temperature)
  t_start = utils::numeric(FLERR,arg[3],false,lmp);

  // 4 = t_stop           (Stopping target temperature)
  t_stop = utils::numeric(FLERR,arg[4],false,lmp);

  // 5 = prony_terms      (number of terms in Prony series)
  prony_terms = utils::inumeric(FLERR,arg[5],false,lmp);

  // 6 = seed             (random seed)
  int seed    = utils::inumeric(FLERR,arg[6],false,lmp);

  // 7 = series type
  if (strcmp(arg[7],"pprony") == 0) {
     series_type = 1;   // series type 1 is 'positive Prony series'
  } else {
     error->all(FLERR,"Fix gld series type must be pprony for now");
  }

  // Error checking for the first set of required input arguments
  if (seed <= 0) error->all(FLERR,"Illegal fix gld command");
  if (prony_terms <= 0)
    error->all(FLERR,"Fix gld prony terms must be > 0");
  if (t_start < 0)
    error->all(FLERR,"Fix gld start temperature must be >= 0");
  if (t_stop < 0)
    error->all(FLERR,"Fix gld stop temperature must be >= 0");
  if (narg - narg_min < 2*(prony_terms) )
    error->all(FLERR,"Fix gld needs more prony series coefficients");

  // allocate memory for Prony series force coefficients
  memory->create(prony_c, prony_terms, "gld:prony_c");
  // allocate memory for Prony series timescale coefficients
  memory->create(prony_tau, prony_terms, "gld:prony_tau");
  // allocate memory for Prony series extended variables
  s_gld = nullptr;
  FixGLD::grow_arrays(atom->nmax);
  // add callbacks to enable restarts
  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  // read in the Prony series coefficients
  int iarg = narg_min;
  int icoeff = 0;
  while (iarg < narg && icoeff < prony_terms) {
    double pc = utils::numeric(FLERR,arg[iarg],false,lmp);
    double ptau = utils::numeric(FLERR,arg[iarg+1],false,lmp);

    if (pc < 0)
      error->all(FLERR,"Fix gld c coefficients must be >= 0");
    if (ptau <= 0)
      error->all(FLERR,"Fix gld tau coefficients must be > 0");

    // All atom types to have the same Prony series
    prony_c[icoeff] = pc;
    prony_tau[icoeff] = ptau;

    icoeff += 1;
    iarg += 2;
  }

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);

  // initialize the extended variables
  init_s_gld();

  // optional arguments
  freezeflag = 0;
  zeroflag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"zero") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix gld command");
      zeroflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"frozen") == 0) {
       if (iarg+2 > narg) error->all(FLERR, "Illegal fix gld command");
       freezeflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
       if (freezeflag) {
         for (int i = 0; i < atom->nlocal; i++) {
           if (atom->mask[i] & groupbit) {
             for (int k = 0; k < 3*prony_terms; k=k+3)
             {
               s_gld[i][k] = 0.0;
               s_gld[i][k+1] = 0.0;
               s_gld[i][k+2] = 0.0;
             }
           }
         }
       }
       iarg += 2;
    }
    else error->all(FLERR,"Illegal fix gld command");
  }

  // Initialize the target temperature
  t_target = t_start;
}

/* ----------------------------------------------------------------------
   Destroys memory allocated by the method
------------------------------------------------------------------------- */

FixGLD::~FixGLD()
{
  delete random;
  memory->destroy(prony_c);
  memory->destroy(prony_tau);
  memory->destroy(s_gld);

  // remove callbacks to fix, so atom class stops calling it
  atom->delete_callback(id,Atom::GROW);
  atom->delete_callback(id,Atom::RESTART);
}

/* ----------------------------------------------------------------------
   Specifies when the fix is called during the timestep
------------------------------------------------------------------------- */

int FixGLD::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ----------------------------------------------------------------------
   Initialize the method parameters before a run
------------------------------------------------------------------------- */

void FixGLD::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (utils::strmatch(update->integrate_style,"^respa"))
    step_respa = (dynamic_cast<Respa *>(update->integrate))->step;
}

/* ----------------------------------------------------------------------
   First half of a timestep (V^{n} -> V^{n+1/2}; X^{n} -> X^{n+1})
------------------------------------------------------------------------- */

void FixGLD::initial_integrate(int /*vflag*/)
{
  double dtfm;
  double ftm2v = force->ftm2v;

  double fran[3], fsum[3], fsumall[3];
  bigint count;

  int icoeff;

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

  // set kT to the temperature in mvvv units
  double kT = (force->boltz)*t_target/(force->mvv2e);

  // zero an accumulator for the total random force
  fsum[0] = fsum[1] = fsum[2] = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        // Advance V by dt/2
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        for (int k = 0; k < 3*prony_terms; k=k+3) {
          v[i][0] += dtfm * s_gld[i][k];
          v[i][1] += dtfm * s_gld[i][k+1];
          v[i][2] += dtfm * s_gld[i][k+2];
        }

        // Advance X by dt
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];

        // Advance S by dt
        icoeff = 0;
        for (int k = 0; k < 3*prony_terms; k=k+3) {
          double theta = exp(-dtv/prony_tau[icoeff]);
          double ck = prony_c[icoeff];
          double vmult = (theta-1.)*ck/ftm2v;
          double rmult = sqrt(2.0*kT*ck/dtv)*(1.-theta)/ftm2v;

          // random force
#ifdef GLD_GAUSSIAN_DISTRO
          fran[0] = rmult*random->gaussian();
          fran[1] = rmult*random->gaussian();
          fran[2] = rmult*random->gaussian();
#endif

#ifdef GLD_UNIFORM_DISTRO
          rmult *= sqrt(12.0); // correct variance of uniform distribution
          fran[0] = rmult*(random->uniform() - 0.5);
          fran[1] = rmult*(random->uniform() - 0.5);
          fran[2] = rmult*(random->uniform() - 0.5);
#endif

          // sum of random forces
          fsum[0] += fran[0];
          fsum[1] += fran[1];
          fsum[2] += fran[2];

          s_gld[i][k]   *= theta;
          s_gld[i][k+1] *= theta;
          s_gld[i][k+2] *= theta;
          s_gld[i][k]   += vmult*v[i][0];
          s_gld[i][k+1] += vmult*v[i][1];
          s_gld[i][k+2] += vmult*v[i][2];
          s_gld[i][k]   += fran[0];
          s_gld[i][k+1] += fran[1];
          s_gld[i][k+2] += fran[2];

          icoeff += 1;
        }
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        // Advance V by dt/2
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        for (int k = 0; k < 3*prony_terms; k=k+3) {
          v[i][0] += dtfm * s_gld[i][k];
          v[i][1] += dtfm * s_gld[i][k+1];
          v[i][2] += dtfm * s_gld[i][k+2];
        }

        // Advance X by dt
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];

        // Advance S by dt
        icoeff = 0;
        for (int k = 0; k < 3*prony_terms; k=k+3) {
          double theta = exp(-dtv/prony_tau[icoeff]);
          double ck = prony_c[icoeff];
          double vmult = (theta-1.)*ck/ftm2v;
          double rmult = sqrt(2.0*kT*ck/dtv)*(1.-theta)/ftm2v;

          // random force
#ifdef GLD_GAUSSIAN_DISTRO
          fran[0] = rmult*random->gaussian();
          fran[1] = rmult*random->gaussian();
          fran[2] = rmult*random->gaussian();
#endif

#ifdef GLD_UNIFORM_DISTRO
          rmult *= sqrt(12.0); // correct variance of uniform distribution
          fran[0] = rmult*(random->uniform() - 0.5);
          fran[1] = rmult*(random->uniform() - 0.5);
          fran[2] = rmult*(random->uniform() - 0.5);
#endif

          // sum of random forces
          fsum[0] += fran[0];
          fsum[1] += fran[1];
          fsum[2] += fran[2];

          s_gld[i][k]   *= theta;
          s_gld[i][k+1] *= theta;
          s_gld[i][k+2] *= theta;
          s_gld[i][k]   += vmult*v[i][0];
          s_gld[i][k+1] += vmult*v[i][1];
          s_gld[i][k+2] += vmult*v[i][2];
          s_gld[i][k]   += fran[0];
          s_gld[i][k+1] += fran[1];
          s_gld[i][k+2] += fran[2];

          icoeff += 1;

        }
      }
  }

  // correct the random force, if zeroflag is set
  if (zeroflag) {
    count = group->count(igroup);
    if (count == 0) error->all(FLERR,"Cannot zero gld force for zero atoms");

    MPI_Allreduce(fsum,fsumall,3,MPI_DOUBLE,MPI_SUM,world);
    fsumall[0] /= (count*prony_terms);
    fsumall[1] /= (count*prony_terms);
    fsumall[2] /= (count*prony_terms);
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        for (int k = 0; k < 3*prony_terms; k=k+3) {
          s_gld[i][k]   -= fsumall[0];
          s_gld[i][k+1] -= fsumall[1];
          s_gld[i][k+2] -= fsumall[2];
        }
      }
    }
  }

}

/* ----------------------------------------------------------------------
   Second half of a timestep (V^{n+1/2} -> V^{n+1})
------------------------------------------------------------------------- */

void FixGLD::final_integrate()
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

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        for (int k = 0; k < 3*prony_terms; k=k+3) {
          v[i][0] += dtfm * s_gld[i][k];
          v[i][1] += dtfm * s_gld[i][k+1];
          v[i][2] += dtfm * s_gld[i][k+2];
        }
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        for (int k = 0; k < 3*prony_terms; k=k+3) {
          v[i][0] += dtfm * s_gld[i][k];
          v[i][1] += dtfm * s_gld[i][k+1];
          v[i][2] += dtfm * s_gld[i][k+2];
        }
      }
  }

  // Change the temperature for the next step
  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  t_target = t_start + delta * (t_stop - t_start);
}

/* ---------------------------------------------------------------------- */

void FixGLD::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * (force->ftm2v);

  // innermost level - GLD update of v and x
  // all other levels - GLD update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixGLD::final_integrate_respa(int ilevel, int /*iloop*/)
{
  dtf = 0.5 * step_respa[ilevel] * (force->ftm2v);
  final_integrate();
}

/* ----------------------------------------------------------------------
   Called when a change to the target temperature is requested mid-run
------------------------------------------------------------------------- */

void FixGLD::reset_target(double t_new)
{
  t_target = t_start = t_stop = t_new;
}

/* ----------------------------------------------------------------------
   Called when a change to the timestep is requested mid-run
------------------------------------------------------------------------- */

void FixGLD::reset_dt()
{
  // set the time integration constants
  dtv = update->dt;
  dtf = 0.5 * update->dt * (force->ftm2v);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixGLD::memory_usage()
{
  double bytes = (double)atom->nmax*3*prony_terms*sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixGLD::grow_arrays(int nmax)
{
  memory->grow(s_gld, nmax, 3*prony_terms,"gld:s_gld");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixGLD::copy_arrays(int i, int j, int /*delflag*/)
{
  for (int k = 0; k < 3*prony_terms; k++) {
    s_gld[j][k] = s_gld[i][k];
  }
}

/* ----------------------------------------------------------------------
   Pack extended variables assoc. w/ atom i into buffer for exchange
   with another processor
------------------------------------------------------------------------- */

int FixGLD::pack_exchange(int i, double *buf)
{
  int m = 0;
  for (int k = 0; k < 3*prony_terms; k++) {
    buf[m++] = s_gld[i][k];
  }
  return m;
}

/* ----------------------------------------------------------------------
   Unpack extended variables from exchange with another processor
------------------------------------------------------------------------- */

int FixGLD::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  for (int k = 0; k < 3*prony_terms; k++) {
    s_gld[nlocal][k] = buf[m++];
  }
  return m;
}


/* ----------------------------------------------------------------------
   Pack extended variables assoc. w/ atom i into buffer for
   writing to a restart file
------------------------------------------------------------------------- */

int FixGLD::pack_restart(int i, double *buf)
{
  int m = 0;
  // pack buf[0] this way because other fixes unpack it
  buf[m++] = 3*prony_terms + 1;
  for (int k = 0; k < 3*prony_terms; k=k+3)
  {
    buf[m++] = s_gld[i][k];
    buf[m++] = s_gld[i][k+1];
    buf[m++] = s_gld[i][k+2];
  }
  return m;
}

/* ----------------------------------------------------------------------
   Unpack extended variables to restart the fix from a restart file
------------------------------------------------------------------------- */

void FixGLD::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to the nth set of extended variables
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i< nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  for (int k = 0; k < 3*prony_terms; k=k+3)
  {
    s_gld[nlocal][k] = extra[nlocal][m++];
    s_gld[nlocal][k+1] = extra[nlocal][m++];
    s_gld[nlocal][k+2] = extra[nlocal][m++];
  }
}

/* ----------------------------------------------------------------------
   Returns the number of items in atomic restart data associated with
   local atom nlocal.  Used in determining the total extra data stored by
   fixes on a given processor.
------------------------------------------------------------------------- */

int FixGLD::size_restart(int /*nlocal*/)
{
  return 3*prony_terms+1;
}

/* ----------------------------------------------------------------------
   Returns the maximum number of items in atomic restart data
   Called in Modify::restart for peratom restart.
------------------------------------------------------------------------- */

int FixGLD::maxsize_restart()
{
  return 3*prony_terms+1;
}

/* ----------------------------------------------------------------------
   Initializes the extended variables to equilibrium distribution
   at t_start.
------------------------------------------------------------------------- */

void FixGLD::init_s_gld()
{
  int icoeff;
  double eq_sdev=0.0;

  // set kT to the temperature in mvvv units
  double kT = (force->boltz)*t_start/(force->mvv2e);

#ifdef GLD_GAUSSIAN_DISTRO
  double scale = sqrt(kT)/(force->ftm2v);
#endif

#ifdef GLD_UNIFORM_DISTRO
  double scale = sqrt(12.0*kT)/(force->ftm2v);
#endif

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      icoeff = 0;
      for (int k = 0; k < 3*prony_terms; k=k+3) {
        eq_sdev = scale*sqrt(prony_c[icoeff]/prony_tau[icoeff]);
#ifdef GLD_GAUSSIAN_DISTRO
        s_gld[i][k] = eq_sdev*random->gaussian();
        s_gld[i][k+1] = eq_sdev*random->gaussian();
        s_gld[i][k+2] = eq_sdev*random->gaussian();
#endif

#ifdef GLD_UNIFORM_DISTRO
        s_gld[i][k] = eq_sdev*(random->uniform()-0.5);
        s_gld[i][k+1] = eq_sdev*(random->uniform()-0.5);
        s_gld[i][k+2] = eq_sdev*(random->uniform()-0.5);
#endif
        icoeff += 1;
      }
    }
  }
}
