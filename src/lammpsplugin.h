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

#ifdef __cplusplus
extern "C" {
#endif

  typedef void *(lammpsplugin_factory)(void *);
  typedef void (*lammpsplugin_initfunc)(void *);

  typedef struct {
    const char *version;
    const char *style;
    const char *name;
    const char *info;
    lammpsplugin_factory *creator;
  } lammpsplugin_t;

  void lammpsplugin_init(void *);

#ifdef __cplusplus
}
#endif

namespace LAMMPS_NS
{
  extern  void lammpsplugin_load(const char *, void *);
  extern  void lammpsplugin_register(lammpsplugin_t *, void *);
}

#endif
