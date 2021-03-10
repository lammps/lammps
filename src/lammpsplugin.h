/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_LAMMPSPLUGIN_H
#define LMP_LAMMPSPLUGIN_H

// C style API and data structures required for dynamic loading

extern "C" {

  typedef void *(lammpsplugin_factory)(void *);
  typedef void (*lammpsplugin_initfunc)(void *);

  typedef struct {
    const char *version;
    const char *style;
    const char *name;
    const char *info;
    lammpsplugin_factory *creator;
  } lammpsplugin_t;

  typedef struct {
    const char *style;
    const char *name;
    const void *handle;
  } lammpsplugin_entry_t;

  // prototype for initializer function required
  // to load a plugin; uses C bindings

  void lammpsplugin_init(void *);
}

namespace LAMMPS_NS
{
  extern void lammpsplugin_load(const char *, void *);
  extern void lammpsplugin_register(lammpsplugin_t *, void *);
  extern int lammpsplugin_get_num_plugins();
  extern const lammpsplugin_entry_t *lammpsplugin_info(int);
  extern int lammpsplugin_find(const char *, const char *);
  extern void lammpsplugin_unload(const char *, const char *, void *);
}

#endif
