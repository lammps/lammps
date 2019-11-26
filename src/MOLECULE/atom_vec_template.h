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

#ifdef ATOM_CLASS

AtomStyle(template,AtomVecTemplate)

#else

#ifndef LMP_ATOM_VEC_TEMPLATE_H
#define LMP_ATOM_VEC_TEMPLATE_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecTemplate : public AtomVec {
 public:
  AtomVecTemplate(class LAMMPS *);
  void process_args(int, char **);
  void create_atom(int, double *);
  void data_atom(double *, imageint, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Molecule template ID for atom_style template does not exist

Self-explanatory.

E: Atom style template molecule must have atom types

The defined molecule(s) does not specify atom types.

E: Invalid template index in Atoms section of data file

The template indices must be between 1 to N, where N is the number of
molecules in the template.

E: Invalid template atom in Atoms section of data file

The atom indices must be between 1 to N, where N is the number of
atoms in the template molecule the atom belongs to.

*/
