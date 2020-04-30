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

/*
   Variant of the C style library interface to LAMMPS
   that uses a shared library and dynamically opens it,
   so this can be used as a prototype code to integrate
   a LAMMPS plugin to some other software.
*/

#include "library.h"
#include "liblammpsplugin.h"
#include <stdlib.h>
#include <dlfcn.h>

liblammpsplugin_t *liblammpsplugin_load(const char *lib)
{
  liblammpsplugin_t *lmp;
  void *handle;

  if (lib == NULL) return NULL;
  handle = dlopen(lib,RTLD_NOW|RTLD_GLOBAL);
  if (handle == NULL) return NULL;
  
  lmp = (liblammpsplugin_t *) malloc(sizeof(liblammpsplugin_t));
  lmp->handle = handle;

#define ADDSYM(symbol) lmp->symbol = dlsym(handle,"lammps_" #symbol)
  ADDSYM(open);
  ADDSYM(open_no_mpi);
  ADDSYM(close);
  ADDSYM(version);
  ADDSYM(file);
  ADDSYM(command);
  ADDSYM(commands_list);
  ADDSYM(commands_string);
  ADDSYM(free);
  ADDSYM(extract_setting);
  ADDSYM(extract_global);
  ADDSYM(extract_box);
  ADDSYM(extract_atom);
  ADDSYM(extract_compute);
  ADDSYM(extract_fix);
  ADDSYM(extract_variable);

  ADDSYM(get_thermo);
  ADDSYM(get_natoms);

  ADDSYM(set_variable);
  ADDSYM(reset_box);

  ADDSYM(gather_atoms);
  ADDSYM(gather_atoms_concat);
  ADDSYM(gather_atoms_subset);
  ADDSYM(scatter_atoms);
  ADDSYM(scatter_atoms_subset);

  ADDSYM(set_fix_external_callback);

  ADDSYM(config_has_package);
  ADDSYM(config_package_count);
  ADDSYM(config_package_name);
  ADDSYM(config_has_gzip_support);
  ADDSYM(config_has_png_support);
  ADDSYM(config_has_jpeg_support);
  ADDSYM(config_has_ffmpeg_support);
  ADDSYM(config_has_exceptions);
  ADDSYM(create_atoms);
#ifdef LAMMPS_EXCEPTIONS
  lmp->has_exceptions = 1;
  ADDSYM(has_error);
  ADDSYM(get_last_error_message);
#else
  lmp->has_exceptions = 0;
  lmp->has_error = NULL;
  lmp->get_last_error_message = NULL;
#endif
  return lmp;
}

int liblammpsplugin_release(liblammpsplugin_t *lmp)
{
  if (lmp == NULL) return 1;
  if (lmp->handle == NULL) return 2;

  dlclose(lmp->handle);
  free((void *)lmp);
  return 0;
}
