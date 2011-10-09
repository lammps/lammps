/* -*- c++ -*- -------------------------------------------------------------
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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifndef LMP_THR_DATA_H
#define LMP_THR_DATA_H

namespace LAMMPS_NS {

// per thread data accumulators
class ThrData {
  friend class FixOMP;
  friend class ThrOMP;

 public:
  ThrData(int tid);
  ~ThrData();

  void clear(int);  // erase contents
  void grow_arrays(int,int,int); // grow per atom arrays

  double get_vdwl() const { return eng_vdwl; };
  double get_coul() const { return eng_coul; };
  double get_bond() const { return eng_bond; };
  const double *get_virial() const { return virial; };

  void add_vdwl(const double &val)  { eng_vdwl += val; };
  void add_coul(const double &val)  { eng_coul += val; };
  void add_bond(const double &val)  { eng_bond += val; };
  void add_virial(const double v[6]) 
    { for (int i=0; i<6; ++i) virial[i] += v[i]; };
  
 protected:
  double eng_vdwl;   // non-bonded non-coulomb energy
  double eng_coul;   // non-bonded coulomb energy
  double eng_bond;   // bonded energy
  double virial[6];  // virial

  double *eatom;     // per atom total energy
  double *vatom;     // per atom virial

 private:
  int _maxeatom;     // size of eatom array
  int _maxvatom;     // size of vatom array
  int _tid;          // my thread id

 public:
//  double memory_usage();

 // disabled default methods
 private:
  ThrData() {};
};

}
#endif
