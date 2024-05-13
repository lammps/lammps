/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

#include "liblammpsplugin.h"

#if defined(_WIN32)

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#if defined(_WIN32_WINNT)
#undef _WIN32_WINNT
#endif

// target Windows version is windows 7 and later
#define _WIN32_WINNT _WIN32_WINNT_WIN7
#define PSAPI_VERSION 2

#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include <stdlib.h>


liblammpsplugin_t *liblammpsplugin_load(const char *lib)
{
  liblammpsplugin_t *lmp;
  void *handle;

  if (lib == NULL) return NULL;

#ifdef _WIN32
  handle = (void *) LoadLibrary(lib);
#else
  handle = dlopen(lib,RTLD_NOW|RTLD_GLOBAL);
#endif
  if (handle == NULL) return NULL;

  lmp = (liblammpsplugin_t *) calloc(1, sizeof(liblammpsplugin_t));
  lmp->abiversion = LAMMPSPLUGIN_ABI_VERSION;
  lmp->handle = handle;

#ifdef _WIN32
#define ADDSYM(symbol) *(void **) (&lmp->symbol) = (void *) GetProcAddress((HINSTANCE) handle, "lammps_" #symbol)
#else
#define ADDSYM(symbol) *(void **) (&lmp->symbol) = dlsym(handle,"lammps_" #symbol)
#endif

#if defined(LAMMPS_LIB_MPI)
  ADDSYM(open);
#else
  lmp->open = NULL;
#endif

  ADDSYM(open_no_mpi);
  ADDSYM(open_fortran);
  ADDSYM(close);

  ADDSYM(mpi_init);
  ADDSYM(mpi_finalize);
  ADDSYM(kokkos_finalize);
  ADDSYM(python_finalize);

  ADDSYM(error);

  ADDSYM(file);
  ADDSYM(command);
  ADDSYM(commands_list);
  ADDSYM(commands_string);

  ADDSYM(get_natoms);
  ADDSYM(get_thermo);
  ADDSYM(last_thermo);

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
  ADDSYM(extract_variable_datatype);
  ADDSYM(set_variable);
  ADDSYM(set_string_variable);
  ADDSYM(set_internal_variable);
  ADDSYM(variable_info);

  ADDSYM(gather_atoms);
  ADDSYM(gather_atoms_concat);
  ADDSYM(gather_atoms_subset);
  ADDSYM(scatter_atoms);
  ADDSYM(scatter_atoms_subset);

  ADDSYM(gather_bonds);
  ADDSYM(gather_angles);
  ADDSYM(gather_dihedrals);
  ADDSYM(gather_impropers);

  ADDSYM(gather);
  ADDSYM(gather_concat);
  ADDSYM(gather_subset);
  ADDSYM(scatter);
  ADDSYM(scatter_subset);

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

  ADDSYM(encode_image_flags);
  ADDSYM(decode_image_flags);

  ADDSYM(set_fix_external_callback);
  ADDSYM(fix_external_get_force);
  ADDSYM(fix_external_set_energy_global);
  ADDSYM(fix_external_set_energy_peratom);
  ADDSYM(fix_external_set_virial_global);
  ADDSYM(fix_external_set_virial_peratom);
  ADDSYM(fix_external_set_vector_length);
  ADDSYM(fix_external_set_vector);

  ADDSYM(flush_buffers);

  ADDSYM(free);

  ADDSYM(is_running);
  ADDSYM(force_timeout);

  lmp->has_exceptions = lmp->config_has_exceptions();
  if (lmp->has_exceptions) {
    ADDSYM(has_error);
    ADDSYM(get_last_error_message);
  }

  ADDSYM(python_api_version);
  return lmp;
}

int liblammpsplugin_release(liblammpsplugin_t *lmp)
{
  if (lmp == NULL) return 1;
  if (lmp->handle == NULL) return 2;

#ifdef _WIN32
  FreeLibrary((HINSTANCE) lmp->handle);
#else
  dlclose(lmp->handle);
#endif
  free((void *)lmp);
  return 0;
}
