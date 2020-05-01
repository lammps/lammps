/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the USER-SMD package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

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

#ifdef ATOM_CLASS

AtomStyle(smd,AtomVecSMD)

#else

#ifndef LMP_ATOM_VEC_SMD_H
#define LMP_ATOM_VEC_SMD_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecSMD : public AtomVec {
 public:
  AtomVecSMD(class LAMMPS *);

  void grow_pointers();
  void force_clear(int, size_t);
  void create_atom_post(int);
  void data_atom_post(int);

 private:
  tagint *molecule;
  double *esph,*desph,*vfrac,*rmass,*radius,*contact_radius;
  double *eff_plastic_strain,*eff_plastic_strain_rate,*damage;
  double **x0,**smd_data_9,**smd_stress,**vest;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

E: Invalid radius in Atoms section of data file

Radius must be >= 0.0.

E: Invalid density in Atoms section of data file

Density value cannot be <= 0.0.

*/
