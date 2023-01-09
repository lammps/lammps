// clang-format off
/* -------------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
   Contributing authors:
             Rodrigo Freitas (UC Berkeley) - rodrigof@berkeley.edu
             Mark Asta (UC Berkeley) - mdasta@berkeley.edu
             Maurice de Koning (Unicamp/Brazil) - dekoning@ifi.unicamp.br
------------------------------------------------------------------------- */

#include "fix_ti_spring.h"

#include "atom.h"
#include "citeme.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "respa.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_fix_ti_spring[] =
  "ti/spring command: doi:10.1016/j.commatsci.2015.10.050\n\n"
  "@article{freitas2016,\n"
  "  author={Freitas, Rodrigo and Asta, Mark and de Koning, Maurice},\n"
  "  title={Nonequilibrium Free-Energy Calculation of Solids Using {LAMMPS}},\n"
  "  journal={Computational Materials Science},\n"
  "  volume={112},\n"
  "  pages={333--341},\n"
  "  year={2016},\n"
  "  publisher={Elsevier}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixTISpring::FixTISpring(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_ti_spring);

  if (narg < 6 || narg > 8)
    error->all(FLERR,"Illegal fix ti/spring command");

  // Flags.
  restart_peratom = 1;
  scalar_flag = 1;
  global_freq = 1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  energy_global_flag = 1;

  // disallow resetting the time step, while this fix is defined
  time_depend = 1;

  // Spring constant.
  k = utils::numeric(FLERR,arg[3],false,lmp);
  if (k <= 0.0) error->all(FLERR,"Illegal fix ti/spring command");

  // Perform initial allocation of atom-based array
  // Register with Atom class
  xoriginal = nullptr;
  FixTISpring::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  // xoriginal = initial unwrapped positions of atoms

  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) domain->unmap(x[i],image[i],xoriginal[i]);
    else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }

  // Time variables.
  t0 = update->ntimestep;  // timestep of original/starting coordinates
  t_switch = utils::bnumeric(FLERR,arg[4],false,lmp); // Number of steps for switching
  t_equil  = utils::bnumeric(FLERR,arg[5],false,lmp); // Number of steps for equilibration
  if ((t_switch <= 0) || (t_equil < 0))
    error->all(FLERR,"Illegal fix ti/spring command");

  // Coupling parameter initialization
  sf = 1;
  if (narg > 6) {
    if (strcmp(arg[6], "function") == 0) sf = utils::inumeric(FLERR,arg[7],false,lmp);
    else error->all(FLERR,"Illegal fix ti/spring switching function");
    if ((sf!=1) && (sf!=2))
      error->all(FLERR,"Illegal fix ti/spring switching function");
  }
  lambda  =  switch_func(0);
  dlambda = dswitch_func(0);

  espring = 0.0;
}

/* ---------------------------------------------------------------------- */

FixTISpring::~FixTISpring()
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,Atom::GROW);
  atom->delete_callback(id,Atom::RESTART);

  // delete locally stored array
  memory->destroy(xoriginal);
}

/* ---------------------------------------------------------------------- */

int FixTISpring::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTISpring::init()
{
  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixTISpring::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixTISpring::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTISpring::post_force(int /*vflag*/)
{
  // do not calculate forces during equilibration
  if ((update->ntimestep - t0) < t_equil) return;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double dx, dy, dz;
  double unwrap[3];

  espring = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - xoriginal[i][0];
      dy = unwrap[1] - xoriginal[i][1];
      dz = unwrap[2] - xoriginal[i][2];
      f[i][0] = (1-lambda) * f[i][0] + lambda * (-k*dx);
      f[i][1] = (1-lambda) * f[i][1] + lambda * (-k*dy);
      f[i][2] = (1-lambda) * f[i][2] + lambda * (-k*dz);
      espring += k * (dx*dx + dy*dy + dz*dz);
    }

  espring *= 0.5;
}

/* ---------------------------------------------------------------------- */

void FixTISpring::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTISpring::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTISpring::initial_integrate(int /*vflag*/)
{
  // Update the coupling parameter value if needed
  if ((update->ntimestep - t0) < t_equil) return;

  const bigint t = update->ntimestep - (t0+t_equil);
  const double r_switch = 1.0/t_switch;

  if ((t >= 0) && (t <= t_switch)) {
    lambda  =  switch_func(t*r_switch);
    dlambda = dswitch_func(t*r_switch);
  }

  if ((t >= t_equil+t_switch) && (t <= (t_equil+2*t_switch))) {
    lambda  =    switch_func(1.0 - (t - t_switch - t_equil)*r_switch);
    dlambda = - dswitch_func(1.0 - (t - t_switch - t_equil)*r_switch);
  }
}

/* ----------------------------------------------------------------------
   energy of stretched springs
------------------------------------------------------------------------- */

double FixTISpring::compute_scalar()
{
  double all;
  MPI_Allreduce(&espring,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

/* ----------------------------------------------------------------------
   information about coupling parameter
------------------------------------------------------------------------- */

double FixTISpring::compute_vector(int n)
{
  linfo[0] = lambda;
  linfo[1] = dlambda;
  return linfo[n];
}

/* ----------------------------------------------------------------------
     memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixTISpring::memory_usage()
{
  double bytes = (double)atom->nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
     allocate atom-based array
------------------------------------------------------------------------- */

void FixTISpring::grow_arrays(int nmax)
{
  memory->grow(xoriginal,nmax,3,"fix_ti/spring:xoriginal");
}

/* ----------------------------------------------------------------------
     copy values within local atom-based array
------------------------------------------------------------------------- */

void FixTISpring::copy_arrays(int i, int j, int /*delflag*/)
{
  xoriginal[j][0] = xoriginal[i][0];
  xoriginal[j][1] = xoriginal[i][1];
  xoriginal[j][2] = xoriginal[i][2];
}

/* ----------------------------------------------------------------------
    pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixTISpring::pack_exchange(int i, double *buf)
{
  buf[0] = xoriginal[i][0];
  buf[1] = xoriginal[i][1];
  buf[2] = xoriginal[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
    unpack values in local atom-based array from exchange with another proc
 ------------------------------------------------------------------------- */

int FixTISpring::unpack_exchange(int nlocal, double *buf)
{
  xoriginal[nlocal][0] = buf[0];
  xoriginal[nlocal][1] = buf[1];
  xoriginal[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
    pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixTISpring::pack_restart(int i, double *buf)
{
  // pack buf[0] this way because other fixes unpack it
  buf[0] = 4;
  buf[1] = xoriginal[i][0];
  buf[2] = xoriginal[i][1];
  buf[3] = xoriginal[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
    unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixTISpring::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  xoriginal[nlocal][0] = extra[nlocal][m++];
  xoriginal[nlocal][1] = extra[nlocal][m++];
  xoriginal[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
     maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixTISpring::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
     size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixTISpring::size_restart(int /*nlocal*/)
{
  return 4;
}

/* ----------------------------------------------------------------------
     Switching function
------------------------------------------------------------------------- */

double FixTISpring::switch_func(double t)
{
  if (sf == 1) return t;

  double t2 = t*t;
  double t5 = t2*t2*t;
  return ((70.0*t2*t2 - 315.0*t2*t + 540.0*t2 - 420.0*t + 126.0)*t5);
}

/* ----------------------------------------------------------------------
     Switching function derivative
------------------------------------------------------------------------- */

double FixTISpring::dswitch_func(double t)
{
  if (sf == 1) return 1.0/t_switch;

  double t2 = t*t;
  double t4 = t2*t2;
  return ((630*t2*t2 - 2520*t2*t + 3780*t2 - 2520*t + 630)*t4) / t_switch;
}
