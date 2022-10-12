// clang-format off
/* ---------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation. Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software. This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------ */

/* ---------------------------------------------------------------------
   Contributing authors: Stephen M. Foiles (SNL)
                         James A. Stewart (SNL)
   ------------------------------------------------------------------ */

#include "fix_electron_stopping_fit.h"

#include "atom.h"
#include "citeme.h"
#include "error.h"
#include "force.h"
#include "math_special.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

// ---------------------------------------------------------------------

static const char cite_fix_electron_stopping_fit_c[] =
  "fix electron/stopping/fit command: doi:10.1063/1.5022471, doi:10.1103/PhysRevB.102.024107\n\n"
  "@Article{Stewart2018,\n"
  " author  = {J. A. Stewart and G. Brookman and P. Price and M. Franco and W. Ji and K. Hattar and R. Dingreville},\n"
  " title   = {Characterizing Single Isolated Radiation-Damage Events from Molecular Dynamics via Virtual Diffraction Methods},\n"
  " journal = {Journal of Applied Physics},\n"
  " year    = {2018},\n"
  " volume  = {123},\n"
  " number  = {16},\n"
  " pages   = {165902}\n"
  "}\n\n"
  "@Article{Lee2020,\n"
  " author  = {C. W. Lee and J. A. Stewart and S. M. Foiles and R. Dingreville and A. Schleife },\n"
  " title   = {Multiscale Simulations of Electron and Ion Dynamics in Self-Irradiated Silicon},\n"
  " journal = {Physical Review~B},\n"
  " year    = {2020},\n"
  " volume  = {102},\n"
  " number  = {2},\n"
  " pages   = {024107}\n"
  "}\n\n";

// ---------------------------------------------------------------------

FixElectronStoppingFit::FixElectronStoppingFit(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp,narg,arg), energy_coh_in(nullptr), v_min_sq(nullptr), v_max_sq(nullptr),
  drag_fac_in_1(nullptr), drag_fac_in_2(nullptr),
  drag_fac_1(nullptr), drag_fac_2(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_electron_stopping_fit_c);

  if (narg < 3 + 3*atom->ntypes) {
     error->all(FLERR,"Incorrect number of fix electron/stopping/fit arguments");
  }

  scalar_flag = 1;
  global_freq = 1;

  energy_coh_in = new double[atom->ntypes+1];

  drag_fac_in_1 = new double[atom->ntypes+1];
  drag_fac_in_2 = new double[atom->ntypes+1];

  for (int i = 1; i <= atom->ntypes; i++) {
     energy_coh_in[i] = utils::numeric(FLERR,arg[3*i],false,lmp);
     drag_fac_in_1[i] = utils::numeric(FLERR,arg[3*i+1],false,lmp);
     drag_fac_in_2[i] = utils::numeric(FLERR,arg[3*i+2],false,lmp);
  };

  v_min_sq = new double[atom->ntypes+1];
  v_max_sq = new double[atom->ntypes+1];

  drag_fac_1 = new double[atom->ntypes+1];
  drag_fac_2 = new double[atom->ntypes+1];

  for (int i = 1; i <= atom->ntypes; i++) {
     double mvv;
     mvv = 2.0*energy_coh_in[i]/force->mvv2e;
     v_min_sq[i] = 1.0*mvv/atom->mass[i];
     v_max_sq[i] = 2.0*mvv/atom->mass[i];
     drag_fac_1[i] = drag_fac_in_1[i];
     drag_fac_2[i] = drag_fac_in_2[i];
  };
};

// ---------------------------------------------------------------------

FixElectronStoppingFit::~FixElectronStoppingFit()
{
  delete [] energy_coh_in;
  delete [] drag_fac_in_1;
  delete [] drag_fac_in_2;
  delete [] drag_fac_1;
  delete [] drag_fac_2;
  delete [] v_min_sq;
  delete [] v_max_sq;
};

// ---------------------------------------------------------------------

int FixElectronStoppingFit::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
};

// ---------------------------------------------------------------------

void FixElectronStoppingFit::init()
{
  electronic_loss_this_node = 0.;
  electronic_loss = 0.;
  f_dot_v_prior = 0.;
  f_dot_v_current = 0.;
  last_step = update->ntimestep;
};

// ---------------------------------------------------------------------

void FixElectronStoppingFit::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
     post_force(vflag);
  else {
     (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(nlevels_respa-1);
     post_force_respa(vflag,nlevels_respa-1,0);
     (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(nlevels_respa-1);
  };
};

// ---------------------------------------------------------------------

void FixElectronStoppingFit::post_force(int /*vflag*/)
{
  double **v = atom->v;
  double **f = atom->f;
  int *type  = atom->type;
  int nlocal = atom->nlocal;

  f_dot_v_current = 0.0;
  for (int i = 0; i < nlocal; i++) {
     double vv = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
     if (vv > v_min_sq[type[i]]) {
        double gamma_x;
        double gamma_y;
        double gamma_z;
        double v_mag = sqrt(vv);
        if (vv < v_max_sq[type[i]]) {
           double frac = (vv - v_min_sq[type[i]])/(v_max_sq[type[i]] - v_min_sq[type[i]]);
           gamma_x = frac*(drag_fac_2[type[i]]*v[i][0] + drag_fac_1[type[i]]);
           gamma_y = frac*(drag_fac_2[type[i]]*v[i][1] + drag_fac_1[type[i]]);
           gamma_z = frac*(drag_fac_2[type[i]]*v[i][2] + drag_fac_1[type[i]]);
        } else {
           gamma_x = drag_fac_2[type[i]]*v[i][0] + drag_fac_1[type[i]];
           gamma_y = drag_fac_2[type[i]]*v[i][1] + drag_fac_1[type[i]];
           gamma_z = drag_fac_2[type[i]]*v[i][2] + drag_fac_1[type[i]];
        };
        f[i][0] -= gamma_x*v[i][0];
        f[i][1] -= gamma_y*v[i][1];
        f[i][2] -= gamma_z*v[i][2];
        f_dot_v_current += v_mag*sqrt( MathSpecial::square(gamma_x*v[i][0])
                                     + MathSpecial::square(gamma_y*v[i][1])
                                     + MathSpecial::square(gamma_z*v[i][2]) );
     };
  };
  this_step = update->ntimestep;
  electronic_loss_this_node += (this_step - last_step)*update->dt*0.5*(f_dot_v_prior + f_dot_v_current);
  last_step = this_step;
  f_dot_v_prior = f_dot_v_current;
};

// ---------------------------------------------------------------------

void FixElectronStoppingFit::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
};

// ---------------------------------------------------------------------

double FixElectronStoppingFit::compute_scalar()
{
  MPI_Allreduce(&electronic_loss_this_node,&electronic_loss,1,MPI_DOUBLE,MPI_SUM,world);
  return electronic_loss;
};

// ---------------------------------------------------------------------
