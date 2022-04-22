/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   -----------------------------------------------------------------------

   This file is a part of the MANIFOLD package.

   Copyright (2013-2014) Stefan Paquay, Eindhoven University of Technology.
   License: GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This file is part of the user-manifold package written by
   Stefan Paquay at the Eindhoven University of Technology.
   This module makes it possible to do MD with particles constrained
   to pretty arbitrary manifolds characterized by some constraint function
   g(x,y,z) = 0 and its normal grad(g). The number of manifolds available
   right now is limited but can be extended straightforwardly by making
   a new class that inherits from manifold and implements all pure virtual
   methods.

   Thanks to Remy Kusters for beta-testing!

------------------------------------------------------------------------- */

#include "manifold_factory.h"

#include "manifold_cylinder.h"
#include "manifold_cylinder_dent.h"
#include "manifold_dumbbell.h"
#include "manifold_ellipsoid.h"
#include "manifold_gaussian_bump.h"
#include "manifold_plane.h"
#include "manifold_plane_wiggle.h"
#include "manifold_sphere.h"
#include "manifold_spine.h"
#include "manifold_supersphere.h"
#include "manifold_thylakoid.h"
#include "manifold_torus.h"

#include <cstring>

namespace LAMMPS_NS {
namespace user_manifold {

  template <typename m_type>
  void make_manifold_if(manifold **man_ptr, const char *name, LAMMPS *lmp, int narg, char **arg)
  {
    if (strcmp(m_type::type(), name) == 0) {
      if (*man_ptr == nullptr) { *man_ptr = new m_type(lmp, narg, arg); }
    }
  }

  manifold *create_manifold(const char *mname, LAMMPS *lmp, int narg, char **arg)
  {
    manifold *man = nullptr;
    make_manifold_if<manifold_cylinder>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_cylinder_dent>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_dumbbell>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_ellipsoid>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_gaussian_bump>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_plane>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_plane_wiggle>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_sphere>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_supersphere>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_spine>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_spine_two>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_thylakoid>(&man, mname, lmp, narg, arg);
    make_manifold_if<manifold_torus>(&man, mname, lmp, narg, arg);

    return man;
  }
}    // namespace user_manifold

}    // namespace LAMMPS_NS
