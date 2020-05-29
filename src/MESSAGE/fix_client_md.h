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

FixStyle(client/md,FixClientMD)

#else

#ifndef LMP_FIX_CLIENT_MD_H
#define LMP_FIX_CLIENT_MD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixClientMD : public Fix {
 public:
  FixClientMD(class LAMMPS *, int, char **);
  ~FixClientMD();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void min_post_force(int);
  double compute_scalar();

 private:
  int maxatom,units,server_error;
  double eng;
  double inv_nprocs;
  double fconvert,econvert,pconvert;
  double box[3][3];
  double *xpbc;

  void pack_coords();
  void pack_box();
  void receive_fev(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
