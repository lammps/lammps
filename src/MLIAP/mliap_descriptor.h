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

#ifndef LMP_MLIAP_DESCRIPTOR_H
#define LMP_MLIAP_DESCRIPTOR_H

#include "pointers.h"

namespace LAMMPS_NS {

class MLIAPDescriptor : protected Pointers  {
public:
  MLIAPDescriptor(LAMMPS*);
  ~MLIAPDescriptor();
  virtual void forward(int*, class NeighList*, double**)=0;
  virtual void backward(class PairMLIAP*, class NeighList*, double**, int)=0;
  virtual void init()=0;
  virtual double get_cutoff(int, int)=0;
  virtual double get_cutmax()=0;
  virtual double memory_usage()=0;

  int ndescriptors;              // number of descriptors
  int nelements;                 // # of unique elements
  char **elements;               // names of unique elements

protected:

};

}

#endif

