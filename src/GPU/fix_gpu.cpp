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

#include "string.h"
#include "stdlib.h"
#include "fix_gpu.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "pair_hybrid_overlay.h"
#include "respa.h"
#include "input.h"
#include "timer.h"
#include "modify.h"
#include "update.h"
#include "domain.h"
#include "universe.h"
#include "gpu_extra.h"
#include "neighbor.h"
#include "citeme.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH};

extern int lmp_init_device(MPI_Comm world, MPI_Comm replica,
                           const int first_gpu, const int last_gpu,
                           const int gpu_mode, const double particle_split,
                           const int nthreads, const int t_per_atom,
                           const double cell_size, char *opencl_flags);
extern void lmp_clear_device();
extern double lmp_gpu_forces(double **f, double **tor, double *eatom,
                             double **vatom, double *virial, double &ecoul);

static const char cite_gpu_package[] =
  "GPU package (short-range, long-range and three-body potentials):\n\n"
  "@Article{Brown11,\n"
  " author = {W. M. Brown, P. Wang, S. J. Plimpton, A. N. Tharrington},\n"
  " title = {Implementing Molecular Dynamics on Hybrid High Performance Computers - Short Range Forces},\n"
  " journal = {Comp.~Phys.~Comm.},\n"
  " year =    2011,\n"
  " volume =  182,\n"
  " pages =   {898--911}\n"
  "}\n\n"
  "@Article{Brown12,\n"
  " author = {W. M. Brown, A. Kohlmeyer, S. J. Plimpton, A. N. Tharrington},\n"
  " title = {Implementing Molecular Dynamics on Hybrid High Performance Computers - Particle-Particle Particle-Mesh},\n"
  " journal = {Comp.~Phys.~Comm.},\n"
  " year =    2012,\n"
  " volume =  183,\n"
  " pages =   {449--459}\n"
  "}\n\n"
  "@Article{Brown13,\n"
  " author = {W. M. Brown, Y. Masako},\n"
  " title = {Implementing Molecular Dynamics on Hybrid High Performance Computers – Three-Body Potentials},\n"
  " journal = {Comp.~Phys.~Comm.},\n"
  " year =    2013,\n"
  " volume =  184,\n"
  " pages =   {2785--2793}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixGPU::FixGPU(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_gpu_package);

  if (lmp->cuda)
    error->all(FLERR,"Cannot use GPU package with USER-CUDA package enabled");

  if (narg < 4) error->all(FLERR,"Illegal package gpu command");

  int ngpu = atoi(arg[3]);
  if (ngpu <= 0) error->all(FLERR,"Illegal package gpu command");
  int first_gpu = 0;
  int last_gpu = ngpu-1;
  
  // options

  _gpu_mode = GPU_NEIGH;
  _particle_split = 1.0;
  int nthreads = 1;
  int newtonflag = 0;
  int threads_per_atom = -1;
  double binsize = -1;
  char *opencl_flags = NULL;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"neigh") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package gpu command");
      if (strcmp(arg[iarg]+1,"yes") == 0) _gpu_mode = GPU_NEIGH;
      else if (strcmp(arg[iarg]+1,"no") == 0) _gpu_mode = GPU_FORCE;
      else if (strcmp(arg[iarg]+1,"hybrid") == 0) _gpu_mode = GPU_HYB_NEIGH;
      else error->all(FLERR,"Illegal package gpu command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"split") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package gpu command");
      _particle_split = force->numeric(FLERR,arg[iarg+1]);
      if (_particle_split == 0.0 || _particle_split > 1.0)
        error->all(FLERR,"Illegal package GPU command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"newton") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package gpu command");
      if (strcmp(arg[iarg]+1,"off") == 0) newtonflag = 0;
      else if (strcmp(arg[iarg]+1,"on") == 0) newtonflag = 1;
      else error->all(FLERR,"Illegal package gpu command");
    } else if (strcmp(arg[iarg],"gpuID") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal package gpu command");
      first_gpu = force->inumeric(FLERR,arg[iarg+1]);
      last_gpu = force->inumeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"tpa") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package gpu command");
      threads_per_atom = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"nthreads") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package gpu command");
      nthreads = force->inumeric(FLERR,arg[iarg+1]);
      if (nthreads < 1) error->all(FLERR,"Illegal fix GPU command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"binsize") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package gpu command");
      binsize = force->numeric(FLERR,arg[iarg+1]);
      if (binsize <= 0.0) error->all(FLERR,"Illegal fix GPU command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"device") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package gpu command");
      opencl_flags = arg[iarg+1];
      iarg += 2;
    } else error->all(FLERR,"Illegal package gpu command");
  }

  // error check

  if ((_gpu_mode == GPU_NEIGH || _gpu_mode == GPU_HYB_NEIGH) && 
      domain->triclinic)
    error->all(FLERR,"Cannot use package gpu neigh yes with triclinic box");

  #ifndef _OPENMP
  if (nthreads > 1)
    error->all(FLERR,"No OpenMP support compiled in");
  #endif

  // set newton pair flag
  // require newtonflag = 0 since currently required by all GPU pair styles

  if (newtonflag == 1) error->all(FLERR,"Illegal package gpu command");

  force->newton_pair = newtonflag;
  if (force->newton_pair || force->newton_bond) force->newton = 1;
  else force->newton = 0;

  // pass params to GPU library

  int gpu_flag = lmp_init_device(universe->uworld, world, first_gpu, last_gpu,
                                 _gpu_mode, _particle_split, nthreads,
                                 threads_per_atom, binsize, opencl_flags);
  GPU_EXTRA::check_flag(gpu_flag,error,world);
}

/* ---------------------------------------------------------------------- */

FixGPU::~FixGPU()
{
  lmp_clear_device();
}

/* ---------------------------------------------------------------------- */

int FixGPU::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGPU::init()
{
  // GPU package cannot be used with atom_style template
  
  if (atom->molecular == 2) 
    error->all(FLERR,"GPU package does not (yet) work with "
               "atom_style template");

  // hybrid cannot be used with force/neigh option

  if (_gpu_mode == GPU_NEIGH || _gpu_mode == GPU_HYB_NEIGH)
    if (force->pair_match("hybrid",1) != NULL ||
        force->pair_match("hybrid/overlay",1) != NULL)
      error->all(FLERR,"Cannot use pair hybrid with GPU neighbor list builds");
  if (_particle_split < 0)
    if (force->pair_match("hybrid",1) != NULL ||
        force->pair_match("hybrid/overlay",1) != NULL)
      error->all(FLERR,"GPU split param must be positive "
                 "for hybrid pair styles");

  // make sure fdotr virial is not accumulated multiple times
  
  if (force->pair_match("hybrid",1) != NULL) {
    PairHybrid *hybrid = (PairHybrid *) force->pair;
    for (int i = 0; i < hybrid->nstyles; i++)
      if (strstr(hybrid->keywords[i],"/gpu")==NULL)
        force->pair->no_virial_fdotr_compute = 1;
  } else if (force->pair_match("hybrid/overlay",1) != NULL) {
    PairHybridOverlay *hybrid = (PairHybridOverlay *) force->pair;
    for (int i = 0; i < hybrid->nstyles; i++)
      if (strstr(hybrid->keywords[i],"/gpu")==NULL)
        force->pair->no_virial_fdotr_compute = 1;
  }

  // rRESPA support

  if (strstr(update->integrate_style,"respa"))
    _nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixGPU::setup(int vflag)
{
  if (_gpu_mode == GPU_NEIGH || _gpu_mode == GPU_HYB_NEIGH)
    if (neighbor->exclude_setting()!=0)
      error->all(FLERR,
                 "Cannot use neigh_modify exclude with GPU neighbor builds");

  if (strstr(update->integrate_style,"verlet")) post_force(vflag);
  else {
    // in setup only, all forces calculated on GPU are put in the outer level
    ((Respa *) update->integrate)->copy_flevel_f(_nlevels_respa-1);
    post_force(vflag);
    ((Respa *) update->integrate)->copy_f_flevel(_nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixGPU::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGPU::post_force(int vflag)
{
  timer->stamp();
  double lvirial[6];
  for (int i = 0; i < 6; i++) lvirial[i] = 0.0;
  double my_eng = lmp_gpu_forces(atom->f, atom->torque, force->pair->eatom,
                                 force->pair->vatom, lvirial,
                                 force->pair->eng_coul);

  force->pair->eng_vdwl += my_eng;
  force->pair->virial[0] += lvirial[0];
  force->pair->virial[1] += lvirial[1];
  force->pair->virial[2] += lvirial[2];
  force->pair->virial[3] += lvirial[3];
  force->pair->virial[4] += lvirial[4];
  force->pair->virial[5] += lvirial[5];

  if (force->pair->vflag_fdotr) force->pair->virial_fdotr_compute();
  timer->stamp(TIME_PAIR);
}

/* ---------------------------------------------------------------------- */

void FixGPU::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGPU::post_force_respa(int vflag, int ilevel, int iloop)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixGPU::memory_usage()
{
  double bytes = 0.0;
  // memory usage currently returned by pair routine
  return bytes;
}

