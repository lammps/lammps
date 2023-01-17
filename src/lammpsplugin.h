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

#ifndef LMP_LAMMPSPLUGIN_H
#define LMP_LAMMPSPLUGIN_H

// C style API and data structure required for dynamic loading
#ifdef __cplusplus
extern "C" {
#endif

typedef void *(lammpsplugin_factory1) (void *);
typedef void *(lammpsplugin_factory2) (void *, int, char **);

typedef struct {
  const char *version;
  const char *style;
  const char *name;
  const char *info;
  const char *author;
  union {
    lammpsplugin_factory1 *v1;
    lammpsplugin_factory2 *v2;
  } creator;
  void *handle;
} lammpsplugin_t;

typedef void (*lammpsplugin_regfunc)(lammpsplugin_t *, void *);
typedef void (*lammpsplugin_initfunc)(void *, void *, void *);

// prototype for initializer function required
// to load a plugin; uses C bindings

void lammpsplugin_init(void *, void *, void *);

#ifdef __cplusplus
}
#endif

#endif
