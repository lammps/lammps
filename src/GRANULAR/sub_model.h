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

#ifndef GRANULAR_SUB_MODEL_H_
#define GRANULAR_SUB_MODEL_H_

#include "contact.h"

using namespace LAMMPS_NS;

namespace Contact{

class SubModel : Pointers{
  SubModel(){};
  virtual ~SubModel();
public:
  int num_coeffs;
  double *coeffs;
  virtual double calculate_forces() = 0;
  void read_restart();
  virtual void parse_coeffs(char **, int);
  void write_restart(FILE*);
  void read_restart(FILE*);
  void read_restart(FILE*, int);
  virtual void coeffs_to_local();
  void allocate_coeffs();
  std::string model_name;
private:
  ContactModel &contact;
  int allocated;
  int material_prop_flag = 0;
  int size_history;

};

}

#endif /* GRANULAR_SUB_MODEL_H_ */
