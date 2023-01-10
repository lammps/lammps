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

#include "fix_manifoldforce.h"  // For stuff

#include "manifold.h"
#include "manifold_factory.h"   // For constructing manifold

#include "atom.h"
#include "error.h"
#include "respa.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace user_manifold;


// Helper functions for parameters/equal style variables in input script
static bool was_var( const char *arg )
{
  return strstr( arg, "v_" ) == arg;
}

static bool str_eq( const char *str1, const char *str2 )
{
  return strcmp(str1,str2) == 0;
}

/* ---------------------------------------------------------------------- */

FixManifoldForce::FixManifoldForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  int me = -1;
  MPI_Comm_rank(world,&me);



  // Command is given as
  // fix <name> <group> manifoldforce manifold_name manifold_args
  if (narg < 5) {
    error->all(FLERR,"Illegal fix manifoldforce! No manifold given");
  }
  const char *m_name = arg[3];
  ptr_m = create_manifold(m_name,lmp,narg,arg);

  // Construct manifold from factory:
  if (!ptr_m)
    error->all(FLERR,"Manifold pointer for manifold '{}' "
                                 "was NULL for some reason", arg[3]);

  // After constructing the manifold, you can safely make
  // room for the parameters
  nvars = ptr_m->nparams();
  if (narg < nvars+4)
    error->all(FLERR,"Manifold {} needs at least {} "
                                 "argument(s)!", m_name, nvars);

  ptr_m->params = new double[nvars];
  if (ptr_m->params == nullptr) {
    error->all(FLERR,"Parameter pointer was NULL!");
  }

  // This part here stores the names/text of each argument,
  // determines which params are equal-style variables,
  // and sets the values of those arguments that were _not_
  // equal style vars (so that they are not overwritten each time step).

  double *params = ptr_m->params;
  for (int i = 0; i < nvars; ++i) {
    if (was_var( arg[i+4] ))
      error->all(FLERR,"Equal-style variables not allowed with fix manifoldforce");

    // Use force->numeric to trigger an error if arg is not a number.
    params[i] = utils::numeric(FLERR,arg[i+4],false,lmp);
  }


  // Perform any further initialization for the manifold that depends on params:
  ptr_m->post_param_init();
}

/* ---------------------------------------------------------------------- */

int FixManifoldForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixManifoldForce::init()
{
  // Check the min-style:
  const bool is_good_min_style = str_eq(update->minimize_style,"hftn")
                                || str_eq(update->minimize_style,"quickmin");
  if (!is_good_min_style) {
    error->all(FLERR,"Fix manifoldforce requires min_style hftn or quickmin");
  }
}

/* ---------------------------------------------------------------------- */

void FixManifoldForce::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    int nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixManifoldForce::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixManifoldForce::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double n[3];
  double invn2;
  double dot;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      // Determine normal of particle:
      ptr_m->n(x[i],n);

      invn2 = 1.0 / ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );
      dot = f[i][0]*n[0] + f[i][1]*n[1] + f[i][2]*n[2];

      f[i][0] -= dot*n[0] * invn2;
      f[i][1] -= dot*n[1] * invn2;
      f[i][2] -= dot*n[2] * invn2;

    }
  }
}

/* ---------------------------------------------------------------------- */

void FixManifoldForce::post_force_respa(int vflag, int /*ilevel*/, int /*iloop*/)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixManifoldForce::min_post_force(int vflag)
{
  post_force(vflag);
}
