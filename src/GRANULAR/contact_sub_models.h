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

#ifndef CONTACT_SUB_MODEL_H_
#define CONTACT_SUB_MODEL_H_

#include "contact.h"
#include "pointers.h"

namespace LAMMPS_NS {
namespace Contact {

class SubModel : protected Pointers {
 public:
  SubModel();
  virtual ~SubModel();

  int num_coeffs;
  double *coeffs;
  void read_restart();
  void parse_coeffs(char **, int);
  virtual void mix_coeffs(SubModel*, SubModel*) {};
  virtual void coeffs_to_local() {};
  void allocate_coeffs();
  std::string name;

  int size_history;
  int nondefault_history_transfer;
  double *transfer_history_factor;

  int history_index;
  int beyond_contact;
  int allow_limit_damping;

  ContactModel *contact;

 protected:
  int allocated;

  double mix_stiffnessE(double, double, double, double);
  double mix_stiffnessG(double, double, double, double);
  double mix_geom(double, double);
};

}    // namespace Contact
}    // namespace LAMMPS_NS

#endif /* CONTACT_SUB_MODEL_H_ */
