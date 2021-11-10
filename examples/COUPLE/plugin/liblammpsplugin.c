/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
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
  ADDSYM(open_fortran);
  ADDSYM(close);

  ADDSYM(mpi_init);
  ADDSYM(mpi_finalize);
  ADDSYM(kokkos_finalize);
  ADDSYM(python_finalize);

  ADDSYM(file);
  ADDSYM(command);
  ADDSYM(commands_list);
  ADDSYM(commands_string);

  ADDSYM(get_natoms);
  ADDSYM(get_thermo);

  ADDSYM(extract_box);
  ADDSYM(reset_box);

  ADDSYM(memory_usage);
  ADDSYM(get_mpi_comm);

  ADDSYM(extract_setting);
  ADDSYM(extract_global_datatype);
  ADDSYM(extract_global);

  ADDSYM(extract_atom_datatype);
  ADDSYM(extract_atom);

  ADDSYM(extract_compute);
  ADDSYM(extract_fix);
  ADDSYM(extract_variable);
  ADDSYM(set_variable);

  ADDSYM(gather_atoms);
  ADDSYM(gather_atoms_concat);
  ADDSYM(gather_atoms_subset);
  ADDSYM(scatter_atoms);
  ADDSYM(scatter_atoms_subset);
  ADDSYM(gather_bonds);

  ADDSYM(create_atoms);

  ADDSYM(find_pair_neighlist);
  ADDSYM(find_fix_neighlist);
  ADDSYM(find_compute_neighlist);
  ADDSYM(neighlist_num_elements);
  ADDSYM(neighlist_element_neighbors);

  ADDSYM(version);
  ADDSYM(get_os_info);

  ADDSYM(config_has_mpi_support);
  ADDSYM(config_has_gzip_support);
  ADDSYM(config_has_png_support);
  ADDSYM(config_has_jpeg_support);
  ADDSYM(config_has_ffmpeg_support);
  ADDSYM(config_has_exceptions);

  ADDSYM(config_has_package);
  ADDSYM(config_package_count);
  ADDSYM(config_package_name);

  ADDSYM(config_accelerator);
  ADDSYM(has_gpu_device);
  ADDSYM(get_gpu_device_info);

  ADDSYM(has_style);
  ADDSYM(style_count);
  ADDSYM(style_name);

  ADDSYM(has_id);
  ADDSYM(id_count);
  ADDSYM(id_name);

  ADDSYM(plugin_count);
  ADDSYM(plugin_name);

  ADDSYM(set_fix_external_callback);
  ADDSYM(fix_external_get_force);
  ADDSYM(fix_external_set_energy_global);
  ADDSYM(fix_external_set_energy_peratom);
  ADDSYM(fix_external_set_virial_global);
  ADDSYM(fix_external_set_virial_peratom);
  ADDSYM(fix_external_set_vector_length);
  ADDSYM(fix_external_set_vector);

  ADDSYM(free);

  ADDSYM(is_running);
  ADDSYM(force_timeout);

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
