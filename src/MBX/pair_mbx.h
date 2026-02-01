/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(mbx, PairMBX)

#else

#ifndef LMP_PAIR_MBX_H
#define LMP_PAIR_MBX_H

#include "pair.h"


namespace LAMMPS_NS {
  class FixMBX;

class PairMBX : public Pair {
  friend FixMBX;    // accesses cut_global

 public:
  PairMBX(class LAMMPS *);
  virtual ~PairMBX();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);

 protected:
  double cut_global;

  int me;

  FixMBX *fix_MBX;    // owner of MBX objects

  int nmolecule;    // # of molecules in system (would break if number of molecules can change)

  double mbx_e1b;
  double mbx_e2b;
  double mbx_e3b;
  double mbx_e4b;
  double mbx_disp;
  double mbx_buck;
  double mbx_ele;
  double mbx_total_energy;

  double mbx_virial[6];

  virtual void allocate();

  void accumulate_f(bool);
  void accumulate_f_all(bool);    // local + ghost
  void accumulate_f_local(bool);
};

}    // namespace LAMMPS_NS

#endif
#endif
