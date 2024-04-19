/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the MACHDYN package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(smd,AtomVecSMD);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_SMD_H
#define LMP_ATOM_VEC_SMD_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecSMD : virtual public AtomVec {
 public:
  AtomVecSMD(class LAMMPS *);

  void grow_pointers() override;
  void force_clear(int, size_t) override;
  void create_atom_post(int) override;
  void data_atom_post(int) override;
  void write_data_restricted_to_general() override;
  void write_data_restore_restricted() override;

 private:
  tagint *molecule;
  double *esph, *desph, *vfrac, *rmass, *radius, *contact_radius;
  double *eff_plastic_strain, *eff_plastic_strain_rate, *damage;
  double **x0, **smd_data_9, **smd_stress, **vest;

  double **x0_hold;
};

}    // namespace LAMMPS_NS

#endif
#endif
