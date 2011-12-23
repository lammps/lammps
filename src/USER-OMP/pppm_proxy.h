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

#ifdef KSPACE_CLASS

KSpaceStyle(pppm/proxy,PPPMProxy)

#else

#ifndef LMP_PPPM_PROXY_H
#define LMP_PPPM_PROXY_H

#include "pppm.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

  class PPPMProxy : public PPPM, public ThrOMP {
 public:
  PPPMProxy(class LAMMPS *, int, char **);
  virtual ~PPPMProxy () {};

  virtual void compute(int, int);
  virtual void compute_proxy(int eflag, int vflag);

  // setup is delayed until the compute proxy is called
  virtual void setup() { need_setup=1; };
  virtual void setup_proxy();

 protected:
  int need_setup;
};

}

#endif
#endif
