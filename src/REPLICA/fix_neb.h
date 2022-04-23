/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(neb,FixNEB);
// clang-format on
#else

#ifndef LMP_FIX_NEB_H
#define LMP_FIX_NEB_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNEB : public Fix {
 public:
  double veng, plen, nlen, dotpath, dottangrad, gradlen, dotgrad;
  int rclimber;

  FixNEB(class LAMMPS *, int, char **);
  ~FixNEB() override;
  int setmask() override;
  void init() override;
  void min_setup(int) override;
  void min_post_force(int) override;

 private:
  int me, nprocs, nprocs_universe;
  double kspring, kspringIni, kspringFinal, kspringPerp, EIniIni, EFinalIni;
  bool StandardNEB, NEBLongRange, PerpSpring, FreeEndIni, FreeEndFinal;
  bool FreeEndFinalWithRespToEIni, FinalAndInterWithRespToEIni;
  int ireplica, nreplica;
  int procnext, procprev;
  int cmode;
  MPI_Comm uworld;
  MPI_Comm rootworld;

  char *id_pe;
  class Compute *pe;

  int nebatoms;
  int ntotal;      // total # of atoms, NEB or not
  int maxlocal;    // size of xprev,xnext,tangent arrays
  double *nlenall;
  double **xprev, **xnext, **fnext, **springF;
  double **tangent;
  double **xsend, **xrecv;      // coords to send/recv to/from other replica
  double **fsend, **frecv;      // coords to send/recv to/from other replica
  tagint *tagsend, *tagrecv;    // ditto for atom IDs

  // info gathered from all procs in my replica
  double **xsendall, **xrecvall;      // coords to send/recv to/from other replica
  double **fsendall, **frecvall;      // force to send/recv to/from other replica
  tagint *tagsendall, *tagrecvall;    // ditto for atom IDs

  int *counts, *displacements;    // used for MPI_Gather

  void inter_replica_comm();
  void reallocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
