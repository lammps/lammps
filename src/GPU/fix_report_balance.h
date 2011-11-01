/* ----------------------------------------------------------------------
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

FixStyle(report/balance,FixReportBalance)

#else

#ifndef LMP_FIX_REPORT_BALANCE_H
#define LMP_FIX_REPORT_BALANCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixReportBalance : public Fix {
 public:
  FixReportBalance(class LAMMPS *, int, char **);
  ~FixReportBalance();
  int setmask();
  void post_run();

 private:
  void report_time(double, int, int, const char *, double *);
};

}

#endif
#endif
