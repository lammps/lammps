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

/*  -----------------------------------------------------------------------
   Contributed by Stefan Paquay @ Brandeis University

   Thanks to Liesbeth Janssen @ Eindhoven University for useful discussions!
   ----------------------------------------------------------------------- */


#include <stdio.h>
#include <string.h>

#include "fix_propel_self.h"

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


static constexpr const bool debug_out = false;


FixPropelSelf::FixPropelSelf( LAMMPS *lmp, int narg, char **argv )
	: Fix(lmp, narg, argv), magnitude(1.0), thermostat_orient(0),
	  mode(VELOCITY), n_types_filter(0), apply_to_type(NULL)
{
  if( narg < 5 ) error->all(FLERR, "Illegal fix propel/self command");

  // The fix is to support the following cases:
  // 1. Simple atoms, in which case the force points along the velocity
  // 2. Cases where the atoms have an internal ellipsoid. Currently those
  //    styles are only body and aspherical particles.
  // The first argument (mode) is used to differentiate between these.
  
  // args: fix ID all propel/self mode magnitude
  // Optional args are
  const char *mode_str = argv[3];
  if (strncmp(mode_str, "velocity", 8) == 0) {
    mode = VELOCITY;
  } else if (strncmp(mode_str, "quaternion", 10) == 0) {
    // This mode should only be supported if the atom style has
    // a quaternion (and if all atoms in the group have it)
    if (verify_atoms_have_quaternion()) {
      error->all(FLERR, "All atoms need a quaternion");
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
    const char *arg = argv[iarg];
    if (strncmp(arg, "types", 5) == 0) {
      // In this case we need to allocate the type list.
      // First find the largest given type:
      int max_type_in_list = 0;
      ++iarg;
      while ( (iarg + n_types_filter < narg) &&
              (!isalpha(argv[iarg + n_types_filter][0]))) {
        int curr_type = force->numeric(FLERR, argv[iarg + n_types_filter]);
        if (curr_type > max_type_in_list) max_type_in_list = curr_type;
        ++n_types_filter;
      }
      if (n_types_filter == 0) {
        error->all(FLERR, "\"types\" option requires at least one type");
      }
      
      // sanity check:
      if (max_type_in_list < n_types_filter) {
        error->warning(FLERR, "Found duplicate type in \"types\" option");
      }
      
      apply_to_type = new int[max_type_in_list+1];
      if (!apply_to_type) error->all(FLERR, "Failed to allocate type storage!");
        
      for (int i = 0; i < max_type_in_list+1; ++i) {
        apply_to_type[i] = 0;
      }
      
      for (int i = 0; i < n_types_filter; ++i) {
        int curr_type = force->numeric(FLERR, argv[iarg + i]);
        apply_to_type[curr_type] = 1;
      }
      // Done handling types argument.
      iarg += n_types_filter; // -1 because we incremented by 1 previously.
    } else {
      char msg[2048];
      sprintf(msg, "Illegal fix propel/self command: "
              "Unrecognized argument %s", arg);
      error->all(FLERR, msg);
    }
  }

}


FixPropelSelf::~FixPropelSelf()
{}


int FixPropelSelf::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;

  return mask;
}


double FixPropelSelf::memory_usage()
{
  // magnitude + thermostat_orient + mode + n_types_filter + apply_to_type
  double bytes = sizeof(double) + 3*sizeof(int) + sizeof(int*);
  bytes += sizeof(int)*n_types_filter;
  
  return bytes;
}


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
    error->all(FLERR, "reached statement that should be unreachable");
  }
}



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
      
      double f_act[3] = { 0.0, 0.0, 1.0 };
      double f_rot[3];

      double *quat  = bonus[ellipsoid[i]].quat;
      tagint *tag = atom->tag;

      double Q[3][3];
      MathExtra::quat_to_mat( quat, Q );
      MathExtra::matvec( Q, f_act, f_rot );

      if (debug_out && comm->me == 0) {
        // Magical reference particle:
        if (tag[i] == 12) {
          fprintf(screen, "rotation quaternion: (%f %f %f %f)\n",
                  quat[0], quat[1], quat[2], quat[3]);
          fprintf(screen, "rotation matrix: / %f %f %f \\\n",
                  Q[0][0] ,Q[0][1], Q[0][2]);
          fprintf(screen, "                 | %f %f %f |\n",
                  Q[1][0] ,Q[1][1], Q[1][2]);
          fprintf(screen, "                 \\ %f %f %f /\n",
                  Q[2][0] ,Q[2][1], Q[2][2]);

          fprintf(screen, "Active force on atom %d: (%f %f %f)\n",
                  tag[i], f_rot[0], f_rot[1], f_rot[2]);
          fprintf(screen, "  Total force before: (%f %f %f)\n",
                  f[i][0], f[i][1], f[i][2]);
        }
      }

      f[i][0] += magnitude * f_rot[0];
      f[i][1] += magnitude * f_rot[1];
      f[i][2] += magnitude * f_rot[2];
    }
  }
}


template <int filter_by_type> 
void FixPropelSelf::post_force_velocity(int /*vflag*/ )
{
  double **f = atom->f;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  int *type = atom->type;

  // Add the active force to the atom force:
  for( int i = 0; i < nlocal; ++i ){
    if( mask[i] & groupbit ){
      if (filter_by_type && !apply_to_type[type[i]]) {
        continue;
      }
	    
      const double *vi = v[i];
      double f_act[3] = { vi[0], vi[1], vi[2] };
      double nv2 = vi[0]*vi[0] + vi[1]*vi[1] + vi[2]*vi[2];
      double fnorm = 0.0;
      constexpr const double TOL = 1e-14;
      if (nv2 > TOL) {
        // Without this check you can run into numerical issues
        // because fnorm will blow up.
        fnorm = magnitude / sqrt(nv2);
      }
      
      if (debug_out && comm->me == 0) {
        // Magical reference particle:
        if (tag[i] == 12) {
          fprintf(screen, "Active force on atom %d: (%f %f %f)\n",
                  tag[i], f_act[0], f_act[1], f_act[2]);
          fprintf(screen, "  Total force before: (%f %f %f)\n",
                  f[i][0], f[i][1], f[i][2]);
        }
      }

      f[i][0] += fnorm * f_act[0];
      f[i][1] += fnorm * f_act[1];
      f[i][2] += fnorm * f_act[2];
    }
  }
}


int FixPropelSelf::verify_atoms_have_quaternion()
{
	int ellipsoid_flag = atom->ellipsoid_flag;
	int body_flag = atom->body_flag;
	int *mask = atom->mask;
  if (! (ellipsoid_flag || body_flag) ){
    error->all(FLERR, "mode quaternion requires body or ellipsoid atoms");
  }
  
  // Make sure all atoms have ellipsoid or body set:
  for (int i = 0; i < atom->nlocal; ++i) {
    if (mask[i] & groupbit) {
      if (ellipsoid_flag && atom->ellipsoid[i] < 0) {
        error->all(FLERR, "Got atom without ellipsoid set");
        // Kind-of pointless return but silences compiler warnings:
        return 1;
      }
      if (body_flag && atom->body[i] < 0) {
        error->all(FLERR, "Got atom without body set");
        // Kind-of pointless return silences compiler warnings:
        return 1;
      }
    }
  }
  
  return 0;
}

