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

KSpaceStyle(pppm/stagger,PPPMStagger)

#else

#ifndef LMP_PPPM_STAGGER_H
#define LMP_PPPM_STAGGER_H

#include "pppm.h"

namespace LAMMPS_NS {

class PPPMStagger : public PPPM {
 public:
  PPPMStagger(class LAMMPS *, int, char **);
  virtual ~PPPMStagger();
  virtual void init();
  virtual void compute(int, int);
  virtual int timing_1d(int, double &);
  virtual int timing_3d(int, double &);

 protected:
  int nstagger;
  double stagger;
  double **gf_b2;

  virtual double compute_qopt();
  double compute_qopt_ad();
  virtual void compute_gf_denom();
  virtual void compute_gf_ik();
  virtual void compute_gf_ad();
  
  virtual void particle_map();
  virtual void make_rho();
  virtual void fieldforce_ik();
  virtual void fieldforce_ad();
  virtual void fieldforce_peratom();


  inline double gf_denom2(const double &x, const double &y,
                         const double &z) const {
    double sx,sy,sz;
    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double xl = x;
    double yl = y;
    double zl = z;
    sx = sy = sz = 0.0;
    for (int l = 0; l < order; l++) {
      sx += gf_b2[order][l]*xl;
      sy += gf_b2[order][l]*yl;
      sz += gf_b2[order][l]*zl;
      xl *= x2;
      yl *= y2;
      zl *= z2;
    }
    double s = sx*sy*sz;
    return s*s;
  };
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot (yet) use kspace_style pppm/stagger with triclinic systems

This feature is not yet supported.

E: Out of range atoms - cannot compute PPPM

One or more atoms are attempting to map their charge to a PPPM grid
point that is not owned by a processor.  This is likely for one of two
reasons, both of them bad.  First, it may mean that an atom near the
boundary of a processor's sub-domain has moved more than 1/2 the
"neighbor skin distance"_neighbor.html without neighbor lists being
rebuilt and atoms being migrated to new processors.  This also means
you may be missing pairwise interactions that need to be computed.
The solution is to change the re-neighboring criteria via the
"neigh_modify"_neigh_modify command.  The safest settings are "delay 0
every 1 check yes".  Second, it may mean that an atom has moved far
outside a processor's sub-domain or even the entire simulation box.
This indicates bad physics, e.g. due to highly overlapping atoms, too
large a timestep, etc.

*/
