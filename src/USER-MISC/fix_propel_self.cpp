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
  : Fix(lmp, narg, argv)
{
  // if( lmp->citeme) lmp->citeme->add(cite_fix_active);
  if( narg < 4 ) error->all(FLERR, "Illegal fix propel/self command");
  
  AtomVecEllipsoid *av = static_cast<AtomVecEllipsoid*>(atom->avec);
  if (!av) error->all(FLERR, "FixPropelSelf requires atom_style ellipsoid");

  if( debug_out && comm->me == 0 ){
    fprintf(screen, "This is fix active, narg is %d\n", narg );
    fprintf(screen, "args:");
    for( int i = 0; i < narg; ++i ){
      fprintf(screen, " %s", argv[i]);
    }
    fprintf(screen, "\n");
  }

  // args: fix ID all active magnitude prop1 prop2 prop3
  // Optional args are
  magnitude = force->numeric( FLERR, argv[3] );
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
  double bytes = 0.0;
  return bytes;
}



void FixPropelSelf::post_force(int /* vflag */ )
{
  // Then do the rest:
  double **f = atom->f;

  AtomVecEllipsoid *av = static_cast<AtomVecEllipsoid*>(atom->avec);

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  AtomVecEllipsoid::Bonus *bonus = av->bonus;
  // Add the active force to the atom force:
  for( int i = 0; i < nlocal; ++i ){
    if( mask[i] & groupbit ){
      double f_act[3] = { 0.0, 0.0, 1.0 };
      double f_rot[3];

      int* ellipsoid = atom->ellipsoid;
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
