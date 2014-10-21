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

#ifdef FIX_CLASS

FixStyle(smd,FixSMD)

#else

#ifndef LMP_FIX_SMD_H
#define LMP_FIX_SMD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSMD : public Fix {
 public:
  FixSMD(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  double compute_vector(int);

  void write_restart(FILE *);
  void restart(char *);

 private:
  double xc,yc,zc,xn,yn,zn,r0;
  double k_smd,f_smd,v_smd;
  int xflag,yflag,zflag;
  int styleflag;
  double r_old,r_now,pmf;

  int igroup2,group2bit;
  double masstotal,masstotal2;
  int nlevels_respa;
  double ftotal[3],ftotal_all[7];
  int force_flag;

  void smd_tether();
  void smd_couple();
};

}

#endif
#endif
