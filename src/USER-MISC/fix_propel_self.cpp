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

/* -----------------------------------------------------------------------
   Contributed by Stefan Paquay @ Brandeis University

   Thanks to Liesbeth Janssen @ Eindhoven University for useful discussions!
----------------------------------------------------------------------- */

#include "fix_propel_self.h"

#include <cstdio>
#include <cstring>
#include <string>

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "math.h"
#include "math_const.h"
#include "math_extra.h"
#include "math_vector.h"
#include "modify.h"
#include "random_mars.h"
#include "respa.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define PRINT_DEBUG_OUTPUT 0

/* ---------------------------------------------------------------------- */

FixPropelSelf::FixPropelSelf( LAMMPS *lmp, int narg, char **argv )
  : Fix(lmp, narg, argv), magnitude(0.0),
    mode(VELOCITY), n_types_filter(0), apply_to_type(NULL)
{
  if (narg < 5) error->all(FLERR, "Illegal fix propel/self command");

  // The fix is to support the following cases:
  // 1. Simple atoms, in which case the force points along the velocity
  // 2. Aspherical particles with an orientation.
  // The first argument (mode) is used to differentiate between these.

  // args: fix ID all propel/self mode magnitude
  // Optional args are

  const char *mode_str = argv[3];

  if (strncmp(mode_str, "velocity", 8) == 0) {
    mode = VELOCITY;

  } else if (strncmp(mode_str, "quat", 4) == 0) {

    // This mode should only be supported if the atom style has
    // a quaternion (and if all atoms in the group have it)

    if (!atoms_have_quaternion()) {
      error->all(FLERR, "All fix atoms need to be extended particles");
    }
    mode = QUATERNION;

  } else {
    char msg[2048];
    sprintf(msg, "Illegal mode \"%s\" for fix propel/self", mode_str);
    error->all(FLERR, msg);
  }

  magnitude = force->numeric( FLERR, argv[4] );

  // Handle rest of args:

  int iarg = 5;
  while (iarg < narg) {

    if (strcmp(argv[iarg],"types") == 0) {

      apply_to_type = new int[atom->ntypes+1];
      memset(apply_to_type, 0, atom->ntypes * sizeof(int));

      // consume all following numerical arguments as types

      iarg++;
      int flag=0;
      while (iarg < narg) {
        if (isdigit(argv[iarg][0])) {
          int thistype = force->inumeric(FLERR,argv[iarg]);
          if ((thistype < 1) || (thistype > atom->ntypes))
            error->all(FLERR,"Illegal atom type to types keyword");
          apply_to_type[thistype] = 1;
          flag = 1;
          iarg++;
        } else break;
      }
      if (!flag)
        error->all(FLERR,"'types' keyword requires at least one type");
      else
        n_types_filter = 1;

    } else {
      error->all(FLERR,"Illegal fix propel/self command.");
    }
  }
}

/* ---------------------------------------------------------------------- */

FixPropelSelf::~FixPropelSelf()
{
  delete[] apply_to_type;
}
/* ---------------------------------------------------------------------- */

int FixPropelSelf::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;

  return mask;
}

/* ---------------------------------------------------------------------- */

double FixPropelSelf::memory_usage()
{
  // magnitude + thermostat_orient + mode + n_types_filter + apply_to_type
  double bytes = sizeof(double) + 3*sizeof(int) + sizeof(int*);
  bytes += sizeof(int)*atom->ntypes*n_types_filter;

  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixPropelSelf::post_force(int vflag )
{
  switch(mode) {
  case QUATERNION:
    if (n_types_filter) post_force_quaternion<1>(vflag);
    else                post_force_quaternion<0>(vflag);
    break;
  case VELOCITY:
    if (n_types_filter) post_force_velocity<1>(vflag);
    else                post_force_velocity<0>(vflag);
    break;
  default:
    ;
  }
}

/* ---------------------------------------------------------------------- */

template <int filter_by_type>
void FixPropelSelf::post_force_quaternion(int /* vflag */ )
{
  double **f = atom->f;
  AtomVecEllipsoid *av = static_cast<AtomVecEllipsoid*>(atom->avec);

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int* ellipsoid = atom->ellipsoid;

  AtomVecEllipsoid::Bonus *bonus = av->bonus;

  // Add the active force to the atom force:

  for( int i = 0; i < nlocal; ++i ){
    if( mask[i] & groupbit ){
      if (filter_by_type && !apply_to_type[type[i]]) {
        continue;
      }

      double f_act[3] = { 1.0, 0.0, 0.0 };
      double f_rot[3];

      double *quat  = bonus[ellipsoid[i]].quat;

      double Q[3][3];
      MathExtra::quat_to_mat( quat, Q );
      MathExtra::matvec( Q, f_act, f_rot );

      f[i][0] += magnitude * f_rot[0];
      f[i][1] += magnitude * f_rot[1];
      f[i][2] += magnitude * f_rot[2];
    }
  }
}

/* ---------------------------------------------------------------------- */

template <int filter_by_type>
void FixPropelSelf::post_force_velocity(int /*vflag*/)
{
  double **f = atom->f;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  // Add the active force to the atom force:

  for(int i = 0; i < nlocal; ++i) {
    if( mask[i] & groupbit ){
      if (filter_by_type && !apply_to_type[type[i]]) {
        continue;
      }

      const double *vi = v[i];
      double f_act[3] = { vi[0], vi[1], vi[2] };
      double nv2 = vi[0]*vi[0] + vi[1]*vi[1] + vi[2]*vi[2];
      double fnorm = 0.0;
      const double TOL = 1e-14;

      if (nv2 > TOL) {

        // Without this check you can run into numerical
        // issues because fnorm will blow up.

        fnorm = magnitude / sqrt(nv2);
      }

      f[i][0] += fnorm * f_act[0];
      f[i][1] += fnorm * f_act[1];
      f[i][2] += fnorm * f_act[2];
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixPropelSelf::atoms_have_quaternion()
{
  if (!atom->ellipsoid_flag) {
    error->all(FLERR, "Mode 'quat' requires atom style ellipsoid");
    return 0;
  }

  int *mask = atom->mask;
  int flag=0,flagall=0;

  // Make sure all atoms have ellipsoid data:

  for (int i = 0; i < atom->nlocal; ++i)
    if (mask[i] & groupbit)
      if (atom->ellipsoid[i] < 0) ++flag;

  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall > 0) return 0;

  return 1;
}
