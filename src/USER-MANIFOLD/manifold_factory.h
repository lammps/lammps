/* -*- c++ -*- ----------------------------------------------------------
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
#include <cstring>

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



/* Here the actual implementation of LAMMPS-related functions begins. */

namespace LAMMPS_NS {

namespace user_manifold {

  // Templated, so needs to be in header.
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

  manifold* create_manifold(const char *, LAMMPS *,
                            int , char ** );

} // namespace user_manifold

} // namespace LAMMPS_NS


#endif // LMP_MANIFOLD_FACTORY_H
