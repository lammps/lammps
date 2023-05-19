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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(dcci,DCCI);
// clang-format on
#else

#ifndef LMP_DCCI_H
#define LMP_DCCI_H

#include "command.h"

namespace LAMMPS_NS {

class DCCI : public Command {
 public:
  DCCI(class LAMMPS *);
  ~DCCI() override;
  void command(int, char **) override;

 private:
  int me, me_universe;    // my proc ID in world and universe
  int iworld, nworlds;    // world info
  double nktv2p;
  MPI_Comm roots;     // MPI comm with 1 root proc from each world
  bigint t_sc;        // Total scaling steps
  int whichfix;       // index of temperature fix to use
  int fixstyle;       // what kind of temperature fix is used
  int *world2root;    // world2root[i] = root proc of world

  double Tcoex, Pcoex;    // temperature and pressure coxistence
  double T_start, T_end;
  double P_start, P_end;
  double lambda, Pcoex_rs;
  int dcci_flag;
  int atomic_flag;
  int *NATOMS;    // numbers of atoms of each phase
  double *PE;     // potential energy of each phase
  double *VOL;    // volume of each phase
  void print_status();
  double lambda_initial;
  double lambda_final;
  int sf;    // scaling function option
  double scaling_function(double, double, double);

  class FixAdaptDCCI *fix_adapt_dcci;
};

}    // namespace LAMMPS_NS

#endif
#endif
