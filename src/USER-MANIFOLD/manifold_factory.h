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


#ifndef LMP_MANIFOLD_FACTORY_H
#define LMP_MANIFOLD_FACTORY_H


#include "manifold.h"
#include "pointers.h"

#include <string>
#include <map>


#include "manifold_cylinder.h"
#include "manifold_cylinder_dent.h"
#include "manifold_dumbbell.h"
#include "manifold_ellipsoid.h"
#include "manifold_plane.h"
#include "manifold_plane_wiggle.h"
#include "manifold_sphere.h"
#include "manifold_supersphere.h"
#include "manifold_spine.h"
#include "manifold_spine_opt.h"
#include "manifold_thylakoid.h"
#include "manifold_torus.h"


inline int str_eql(const char *s1, const char *s2)
{
  return strcmp(s1,s2) == 0;
}

/*
 * Defining USE_PHONY_LAMMPS makes sure that none of the LAMMPS classes are
 * included/compiled. This is done in order to allow other programs to use
 * the manifold_factory without compiling all of LAMMPS itself. The relevant
 * classes/functions are replaced with dummy ones defined in this #ifdef-block:
 */
#ifdef USE_PHONY_LAMMPS
#  ifdef __GNUG__
#  warning Not compiling actual LAMMPS classes!
#  endif


struct Error {
  void all(const char *fname, int line, const char* msg)
  {
    fprintf(stderr,"ERROR: %s (%s:%d)",msg,fname,line);
    std::terminate();
  }

  void one(const char *fname, int line, const char* msg)
  {
    fprintf(stderr,"ERROR: %s (%s:%d)",msg,fname,line);
    std::terminate();
  }
};

struct LAMMPS { };

struct Pointers
{
  Pointers(LAMMPS *) : error( &e ){}
  Error e;
  Error *error;
};

static FILE *screen = fopen("/dev/stdout","w");

#define FLERR __FILE__,__LINE__    // Equivalent to definition in pointers.h
#endif // USE_PHONY_LAMMPS



// Macro to simplify the fix_impl creation:
#define RETURN_MANIFOLD_IF(MNAME, MANIFOLD_TYPE) \
  do { \
    using namespace LAMMPS_NS; \
    if( strcmp( MNAME, MANIFOLD_TYPE::type() ) == 0 ) {\
      return new MANIFOLD_TYPE(lmp,narg,arg); \
    } \
  }while(0)

namespace LAMMPS_NS {

  template <typename m_type>
  void make_manifold_if( manifold **man_ptr, const char *name,
                         LAMMPS *lmp, int narg, char **arg )
  {
    if( strcmp( m_type::type(), name ) == 0 ){
      if( *man_ptr == NULL ){
        *man_ptr = new m_type(lmp, narg, arg);
      }
    }
  }
  
  inline manifold* create_manifold(const char *mname, LAMMPS *lmp,
                                   int narg, char **arg )
  {
    using namespace LAMMPS_NS;
    manifold *man = NULL;
    make_manifold_if<manifold_cylinder>     ( &man, mname, lmp, narg, arg );
    make_manifold_if<manifold_cylinder_dent>( &man, mname, lmp, narg, arg );
    make_manifold_if<manifold_dumbbell>     ( &man, mname, lmp, narg, arg );
    make_manifold_if<manifold_ellipsoid>    ( &man, mname, lmp, narg, arg );
    make_manifold_if<manifold_plane>        ( &man, mname, lmp, narg, arg );
    make_manifold_if<manifold_plane_wiggle> ( &man, mname, lmp, narg, arg );
    make_manifold_if<manifold_sphere>       ( &man, mname, lmp, narg, arg );
    make_manifold_if<manifold_supersphere>  ( &man, mname, lmp, narg, arg );
    make_manifold_if<manifold_spine>        ( &man, mname, lmp, narg, arg );
    // make_manifold_if<manifold_spine_opt>    ( &man, mname, lmp, narg, arg );
    make_manifold_if<manifold_thylakoid>    ( &man, mname, lmp, narg, arg );
    make_manifold_if<manifold_torus>        ( &man, mname, lmp, narg, arg );

    
    return man;
  }

} // namespace LAMMPS_NS


#endif // LMP_MANIFOLD_FACTORY_H
