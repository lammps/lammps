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

//FixStyle(neb,FixNEB)
FixStyle(neb/spin,FixNEB_spin)

#else

#ifndef LMP_FIX_NEB_SPIN_H
#define LMP_FIX_NEB_SPIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNEB_spin : public Fix {
 public:
  double veng,plen,nlen,dotpath,dottangrad,gradlen,dotgrad;
  int rclimber;

  FixNEB_spin(class LAMMPS *, int, char **);
  ~FixNEB_spin();
  int setmask();
  void init();
  void min_setup(int);
  void min_post_force(int);
  //void advance_spins(double);
  //double evaluate_dt();

 private:
  int me,nprocs,nprocs_universe;
  double kspring,kspringIni,kspringFinal,kspringPerp,EIniIni,EFinalIni;
  bool StandardNEB,NEBLongRange,PerpSpring,FreeEndIni,FreeEndFinal;
  bool FreeEndFinalWithRespToEIni,FinalAndInterWithRespToEIni;
  bool SpinLattice;
  int ireplica,nreplica;
  int procnext,procprev;
  int cmode;
  MPI_Comm uworld;
  MPI_Comm rootworld;


  char *id_pe;
  class Compute *pe;

  int nebatoms;
  int ntotal;                  // total # of atoms, NEB or not
  int maxlocal;                // size of xprev,xnext,tangent arrays
  double *nlenall;
  double **xprev,**xnext,**fnext;
  // spin quantities
  double **spprev,**spnext,**fmnext;
  double **springF;
  double **tangent;
  double **xsend,**xrecv;      // coords to send/recv to/from other replica
  double **fsend,**frecv;      // coords to send/recv to/from other replica
  // spin quantities
  double **spsend,**sprecv;      // sp to send/recv to/from other replica
  double **fmsend,**fmrecv;      // fm to send/recv to/from other replica
  tagint *tagsend,*tagrecv;    // ditto for atom IDs

                                 // info gathered from all procs in my replica
  double **xsendall,**xrecvall;    // coords to send/recv to/from other replica
  double **fsendall,**frecvall;    // force to send/recv to/from other replica
  // spin quantities
  double **spsendall,**sprecvall;    // sp to send/recv to/from other replica
  double **fmsendall,**fmrecvall;    // fm to send/recv to/from other replica
  tagint *tagsendall,*tagrecvall;  // ditto for atom IDs

  int *counts,*displacements;   // used for MPI_Gather

  double geodesic_distance(double *, double *);
  void inter_replica_comm();
  void reallocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Potential energy ID for fix neb does not exist

Self-explanatory.

E: Too many active NEB atoms

UNDOCUMENTED

E: Too many atoms for NEB

UNDOCUMENTED

U: Atom count changed in fix neb

This is not allowed in a NEB calculation.

*/
