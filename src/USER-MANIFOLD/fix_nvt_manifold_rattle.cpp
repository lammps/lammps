/* ----------------------------------------------------------------------
   Lammps - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   -----------------------------------------------------------------------

   This file is a part of the USER-MANIFOLD package.

   Copyright (2013-2014) Stefan Paquay, Eindhoven University of Technology.
   License: GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This file is part of the user-manifold package written by
   Stefan Paquay at the Eindhoven University of Technology.
   This module makes it possible to do MD with particles constrained
   to pretty arbitrary manifolds characterised by some constraint function
   g(x,y,z) = 0 and its normal grad(g). The number of manifolds available
   right now is limited but can be extended straightforwardly by making
   a new class that inherits from manifold and implements all pure virtual
   methods.

   Thanks to Remy Kusters for beta-testing!

------------------------------------------------------------------------- */


#include "fix_nvt_manifold_rattle.h"
#include <cstring>
#include <cmath>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "group.h"
#include "citeme.h"
#include "modify.h"
#include "compute.h"

#include "manifold.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace user_manifold;

enum {CONSTANT,EQUAL};
enum {NOBIAS,BIAS};




static const char* cite_fix_nvt_manifold_rattle =
  "fix nvt/manifold/rattle command:\n\n"
  "@article{paquay-2016,\n"
  "   author        = {Paquay, Stefan and Kusters, Remy},\n"
  "   doi           = {10.1016/j.bpj.2016.02.017},\n"
  "   issn          = {0006-3495},\n"
  "   journal       = {Biophysical Journal},\n"
  "   month         = apr,\n"
  "   number        = {6},\n"
  "   pages         = {1226--1233},\n"
  "   title         = {{A Method for Molecular Dynamics on Curved Surfaces}},\n"
  "   volume        = {110},\n"
  "   year          = {2016}\n"
  "}\n\n";
/* ---------------------------------------------------------------------- */

FixNVTManifoldRattle::FixNVTManifoldRattle(LAMMPS *lmp, int narg, char **arg,
                                           int error_on_unknown_keyword )
  : FixNVEManifoldRattle(lmp,narg,arg, 0)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_nvt_manifold_rattle);

  if (narg < 6 ) error->all(FLERR,"Illegal fix nvt/manifold/rattle command");

  // Set all bits/settings:
  dof_flag = 1;
  dthalf = dt4 = dt8 = 0;

  t_start = t_stop = t_period = t_current = t_target = ke_target = 0.0;
  t_freq = drag = tdrag_factor = 0;

  boltz = force->boltz, nktv2p = force->nktv2p;
  tdof = 0;
  mtchain = 3;
  factor_eta = 0.0;
  which = got_temp = 0;

  int argi = 6 + ptr_m->nparams();
  while( argi < narg )
  {
    if (strcmp( arg[argi], "temp") == 0) {
      if (argi+3 >= narg)
        error->all(FLERR,"Keyword 'temp' needs 3 arguments");

      t_start  = force->numeric(FLERR, arg[argi+1]);
      t_stop   = force->numeric(FLERR, arg[argi+2]);
      t_period = force->numeric(FLERR, arg[argi+3]);
      t_target = t_start;
      got_temp = 1;

      argi += 4;
    } else if (strcmp( arg[argi], "tchain" ) == 0) {
      if (argi+1 >= narg)
        error->all(FLERR,"Keyword 'tchain' needs 1 argument");

      mtchain = force->inumeric(FLERR, arg[argi+1]);
      argi += 2;
    } else if (error_on_unknown_keyword) {
      char msg[2048];
      sprintf(msg,"Error parsing arg \"%s\".\n", arg[argi]);
      error->all(FLERR, msg);
    } else {
      argi += 1;
    }
  }

  reset_dt();

  if (!got_temp ) error->all(FLERR,"Fix nvt/manifold/rattle needs 'temp'!");

  if (t_period < 0.0) {
    error->all(FLERR,"Fix nvt/manifold/rattle damping parameter must be > 0.0");
  }

  // Create temperature compute:
  const char *fix_id = arg[1];
  int n = strlen(fix_id)+6;
  id_temp = new char[n];
  strcpy(id_temp,fix_id);
  strcat(id_temp,"_temp");
  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char*) "temp";


  modify->add_compute(3,newarg);
  delete [] newarg;
  int icompute = modify->find_compute(id_temp);
  if (icompute < 0) {
    error->all(FLERR,"Temperature ID for fix nvt/manifold/rattle "
               "does not exist");
  }
  temperature = modify->compute[icompute];
  if (temperature->tempbias) which = BIAS;
  else                        which = NOBIAS;

  // Set t_freq from t_period
  t_freq = 1.0 / t_period;

  // Init Nosé-Hoover chain:
  eta        = new double[mtchain];
  eta_dot    = new double[mtchain+1];
  eta_dotdot = new double[mtchain];
  eta_mass   = new double[mtchain];
  eta_dot[mtchain] = 0.0;

  eta_dot[mtchain] = 0.0;
  for( int ich = 0; ich < mtchain; ++ich ){
    eta[ich] = eta_dot[ich] = eta_dotdot[ich] = 0.0;
  }


}

/* ---------------------------------------------------------------------- */

FixNVTManifoldRattle::~FixNVTManifoldRattle()
{
  // Deallocate heap-allocated objects.
  if (eta)        delete[] eta;
  if (eta_dot)    delete[] eta_dot;
  if (eta_dotdot) delete[] eta_dotdot;
  if (eta_mass)   delete[] eta_mass;

  modify->delete_compute(id_temp);
  if (id_temp)    delete[] id_temp;
}




int FixNVTManifoldRattle::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  if (nevery > 0) mask |= END_OF_STEP;

  return mask;
}


/* --------------------------------------------------------------------------
   Check that force modification happens before position and velocity update.
   Make sure respa is not used.
------------------------------------------------------------------------- */
void FixNVTManifoldRattle::init()
{
  // Makes sure the manifold params are set initially.
  update_var_params();

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0) {
    error->all(FLERR,"Temperature ID for fix nvt/manifold/rattle "
               "does not exist");
  }
  temperature = modify->compute[icompute];
  if (temperature->tempbias) which = BIAS;
  else                        which = NOBIAS;

}



void FixNVTManifoldRattle::setup(int /*vflag*/)
{
  compute_temp_target();

  t_current = temperature->compute_scalar();
  tdof = temperature->dof;

  // Compute/set eta-masses:
  double inv_t_freq2 = 1.0 / (t_freq*t_freq);
  eta_mass[0] = tdof * boltz * t_target * inv_t_freq2;
  for( int ich = 1; ich < mtchain; ++ich ){
    eta_mass[ich] = boltz * t_target * inv_t_freq2;
  }

  for( int ich = 1; ich < mtchain; ++ich ){
    eta_dotdot[ich] = (eta_mass[ich-1]*eta_dot[ich-1]*eta_dot[ich-1] -
                       boltz * t_target ) / eta_mass[ich];
  }
}

void FixNVTManifoldRattle::compute_temp_target()
{

  t_current = temperature->compute_scalar();
  tdof      = temperature->dof;

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0){
    delta /= update->endstep - update->beginstep;
  }

  tdof = temperature->dof;
  t_target = t_start + delta * (t_stop-t_start);
  ke_target = tdof * boltz * t_target;
}

void FixNVTManifoldRattle::nhc_temp_integrate()
{
  int ich;
  // t_current = temperature->compute_scalar();
  // tdof = temperature->dof;
  compute_temp_target();

  double expfac, kecurrent = tdof * boltz * t_current;
  double inv_t_freq2 = 1.0 / (t_freq*t_freq);
  eta_mass[0] = tdof * boltz * t_target * inv_t_freq2;
  for( int ich = 1; ich < mtchain; ++ich ){
    eta_mass[ich] = boltz * t_target * inv_t_freq2;
  }

  if (eta_mass[0] > 0.0) {
    eta_dotdot[0] = (kecurrent - ke_target)/eta_mass[0];
  } else {
    eta_dotdot[0] = 0;
  }

  for( ich = mtchain-1; ich > 0; --ich ){
    expfac = exp(-dt8*eta_dot[ich+1]);
    eta_dot[ich] *= expfac;
    eta_dot[ich] += eta_dotdot[ich] * dt4;
    eta_dot[ich] *= tdrag_factor * expfac;

  }

  expfac = exp(-dt8*eta_dot[1]);
  eta_dot[0] *= expfac;
  eta_dot[0] += eta_dotdot[0] * dt4;
  eta_dot[0] *= tdrag_factor * expfac;

  factor_eta = exp(-dthalf*eta_dot[0]);

  if (factor_eta == 0) {
    char msg[2048];
    sprintf(msg, "WTF, factor_eta is 0! dthalf = %f, eta_dot[0] = %f",
            dthalf, eta_dot[0]);
    error->all(FLERR,msg);
  }

  nh_v_temp();

  t_current *= factor_eta*factor_eta;
  kecurrent = tdof * boltz * t_current;

  if (eta_mass[0] > 0.0) {
    eta_dotdot[0] = (kecurrent - ke_target) / eta_mass[0];
  } else {
    eta_dotdot[0] = 0.0;
  }

  for( int ich = 1; ich < mtchain; ++ich ){
    eta[ich] += dthalf*eta_dot[ich];
  }
  eta_dot[0] *= expfac;
  eta_dot[0] += eta_dotdot[0]*dt4;
  eta_dot[0] *= expfac;

  for( int ich = 1; ich < mtchain; ++ich ){
    expfac = exp(-dt8*eta_dot[ich+1]);
    eta_dot[ich] *= expfac;
    eta_dotdot[ich] = (eta_mass[ich-1]*eta_dot[ich-1]*eta_dot[ich-1]
                       - boltz*t_target) / eta_mass[ich];
    eta_dot[ich] *= eta_dotdot[ich] * dt4;
    eta_dot[ich] *= expfac;
  }

}

void FixNVTManifoldRattle::nh_v_temp()
{
  double **v = atom->v;
  int *mask  = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;




  if (which == NOBIAS) {
    for( int i = 0; i < nlocal; ++i ){
      if (mask[i] & groupbit) {
        v[i][0] *= factor_eta;
        v[i][1] *= factor_eta;
        v[i][2] *= factor_eta;
      }
    }
  } else if (which == BIAS) {
    for( int i = 0; i < nlocal; ++i ){
      if (mask[i] & groupbit) {
        temperature->remove_bias(i,v[i]);
        v[i][0] *= factor_eta;
        v[i][1] *= factor_eta;
        v[i][2] *= factor_eta;
        temperature->restore_bias(i,v[i]);
      }
    }
  }
}




// Most of this logic is based on fix_nh:
void FixNVTManifoldRattle::initial_integrate(int /*vflag*/)
{

  update_var_params();

  compute_temp_target();
  nhc_temp_integrate();

  nve_x_rattle(igroup, groupbit);
}

void FixNVTManifoldRattle::final_integrate()
{
  nve_v_rattle(igroup, groupbit);

  nhc_temp_integrate();
}



/* ---------------------------------------------------------------------- */
void FixNVTManifoldRattle::reset_dt()
{
  FixNVEManifoldRattle::reset_dt();

  dthalf = 0.5 * update->dt;
  dt4 = 0.25 * update->dt;
  dt8 = 0.125 * update->dt;
  tdrag_factor = 1.0 - (update->dt * t_freq * drag);

}





double FixNVTManifoldRattle::memory_usage()
{
  double bytes = FixNVEManifoldRattle::memory_usage();
  bytes += (4*mtchain+1)*sizeof(double);

  return bytes;
}
