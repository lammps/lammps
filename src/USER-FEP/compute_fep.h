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

/* ----------------------------------------------------------------------
   Contributing author: Agilio Padua (ICCF,UBP,CNRS)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(fep,ComputeFEP)

#else

#ifndef COMPUTE_FEP_H
#define COMPUTE_FEP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeFEP : public Compute {
 public:
  ComputeFEP(class LAMMPS *, int, char **);
  ~ComputeFEP();
  void init();
  void compute_vector();

 private:
  int npert;
  int pairflag;
  int chgflag;
  int tailflag, volumeflag;
  int fepinitflag;
  double temp_fep;

  int nmax;
  double *q_orig;
  double **f_orig;
  double eng_vdwl_orig,eng_coul_orig;
  double pvirial_orig[6];
  double *peatom_orig,**pvatom_orig;
  double energy_orig;
  double kvirial_orig[6];
  double *keatom_orig,**kvatom_orig;

  struct Perturb {
    int which,ivar;
    char *var;
    char *pstyle,*pparam;
    int ilo,ihi,jlo,jhi;
    int pdim;
    double **array,**array_orig;
    int aparam;
  };

  Perturb *perturb;

  double compute_epair();
  void perturb_params();
  void backup_params();
  void restore_params();
  void allocate_storage();
  void deallocate_storage();
  void backup_qfev();
  void restore_qfev();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for compute fep does not exist

Self-explanatory.

E: Variable for compute fep is invalid style

Self-explanatory.

E: Compute fep pair style does not exist

Self-explanatory.

E: Energy was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied energy, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
