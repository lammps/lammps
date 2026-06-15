/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

   Heavily modified and adapted for LAMMPS by the LAMMPS developers.
------------------------------------------------------------------------- */

#ifndef LMP_REAXFF_INLINE_H
#define LMP_REAXFF_INLINE_H

#include "accelerator_kokkos.h"    // for LAMMPS_INLINE

namespace ReaxFF {
struct LR_data {
  double H;
  double e_vdW, CEvd;
  double e_ele, CEclmb;

  LR_data() = default;
  ~LR_data() = default;
  LR_data(LR_data &&) = delete;
  LR_data(const LR_data &) = delete;
  LR_data &operator=(const LR_data &) = delete;
  LR_data &operator=(LR_data &&) = delete;
  void swap(LR_data &) = delete;
};

struct cubic_spline_coef {
  double a, b, c, d;

  cubic_spline_coef() = default;
  cubic_spline_coef(cubic_spline_coef &&) = default;
  cubic_spline_coef &operator=(cubic_spline_coef &&) = default;
  ~cubic_spline_coef() = default;

  LAMMPS_INLINE
  cubic_spline_coef(const cubic_spline_coef &_c)
  {
    a = _c.a;
    b = _c.b;
    c = _c.c;
    d = _c.d;
  }

  LAMMPS_INLINE
  cubic_spline_coef &operator=(const cubic_spline_coef &rhs) // NOLINT
  {
    a = rhs.a;
    b = rhs.b;
    c = rhs.c;
    d = rhs.d;
    return *this;
  }
};
}    // namespace ReaxFF

#endif
