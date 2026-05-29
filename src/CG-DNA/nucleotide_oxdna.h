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

#ifndef NUCLEOTIDE_OXDNA_H
#define NUCLEOTIDE_OXDNA_H

#include "pointers.h"
#include "constants_oxdna.h"

namespace LAMMPS_NS {

template <class Derived>
class NucleotideOxdna {
 public:
  NucleotideOxdna() = default;
  virtual ~NucleotideOxdna() = default;
  void backbone_site_interface(double e1[3], double e2[3], double e3[3], double rbk[3]) {
    static_cast<Derived*>(this)->backbone_site(e1, e2, e3, rbk);
  }
  void stacking_site_interface(double e1[3], double e2[3], double e3[3], double rstk[3]) {
    static_cast<Derived*>(this)->stacking_site(e1, e2, e3, rstk);
  }
  template <int N>
  void base_site_interface(double e1[3], double e2[3], double e3[3], double rbs[3]) {
    static_cast<Derived*>(this)->template base_site<N>(e1, e2, e3, rbs);
  }
};

/* ----------------------------------------------------------------------
   oxDNA1 nucleotide
------------------------------------------------------------------------- */
class NucleotideOxdna1 : public NucleotideOxdna<NucleotideOxdna1> {
 public:
  inline void backbone_site(double e1[3], double /*e2*/[3], double /*e3*/[3], double rbk[3]) {
    double dx_cbk_oxdna1 = ConstantsOxdna::get_dx_cbk_oxdna1();
    rbk[0] = dx_cbk_oxdna1 * e1[0];
    rbk[1] = dx_cbk_oxdna1 * e1[1];
    rbk[2] = dx_cbk_oxdna1 * e1[2];
  }
  inline void stacking_site(double e1[3], double /*e2*/[3], double /*e3*/[3], double rstk[3]) {
    double dx_cstk_oxdna1 = ConstantsOxdna::get_dx_cstk_oxdna1();
    rstk[0] = dx_cstk_oxdna1 * e1[0];
    rstk[1] = dx_cstk_oxdna1 * e1[1];
    rstk[2] = dx_cstk_oxdna1 * e1[2];
  }
  template <int N>
  inline void base_site(double e1[3], double /*e2*/[3], double /*e3*/[3], double rbs[3]);
};

template <>
inline void NucleotideOxdna1::base_site<0>(double e1[3], double /*e2*/[3],
  double /*e3*/[3], double rbs[3]) {
  double dx_cbs_oxdna1 = ConstantsOxdna::get_dx_cbs_oxdna1();
  rbs[0] = dx_cbs_oxdna1*e1[0];
  rbs[1] = dx_cbs_oxdna1*e1[1];
  rbs[2] = dx_cbs_oxdna1*e1[2];
}

/* ----------------------------------------------------------------------
   oxDNA2 nucleotide
------------------------------------------------------------------------- */
class NucleotideOxdna2 : public NucleotideOxdna<NucleotideOxdna2> {
 public:
  inline void backbone_site(double e1[3], double e2[3], double /*e3*/[3], double rbk[3]) {
    double dx_cbk_oxdna2 = ConstantsOxdna::get_dx_cbk_oxdna2();
    double dy_cbk_oxdna2 = ConstantsOxdna::get_dy_cbk_oxdna2();
    rbk[0] = dx_cbk_oxdna2 * e1[0] + dy_cbk_oxdna2 * e2[0];
    rbk[1] = dx_cbk_oxdna2 * e1[1] + dy_cbk_oxdna2 * e2[1];
    rbk[2] = dx_cbk_oxdna2 * e1[2] + dy_cbk_oxdna2 * e2[2];
  }
};

/* ----------------------------------------------------------------------
   oxDNA3 nucleotide
------------------------------------------------------------------------- */
class NucleotideOxdna3 : public NucleotideOxdna<NucleotideOxdna3> {
 public:
  inline void stacking_site(double e1[3], double /*e2*/[3], double /*e3*/[3], double rstk[3]) {
    double dx_cstk_oxdna3 = ConstantsOxdna::get_dx_cstk_oxdna3();
    rstk[0] = dx_cstk_oxdna3 * e1[0];
    rstk[1] = dx_cstk_oxdna3 * e1[1];
    rstk[2] = dx_cstk_oxdna3 * e1[2];
  }
  template <int N>
  inline void base_site(double e1[3], double /*e2*/[3], double /*e3*/[3], double rbs[3]);
};

template <>
inline void NucleotideOxdna3::base_site<0>(double e1[3], double /*e2*/[3],
  double /*e3*/[3], double rbs[3]) {
  double dx_cbs_pyr_oxdna3 = ConstantsOxdna::get_dx_cbs_pyr_oxdna3();
  rbs[0] = dx_cbs_pyr_oxdna3*e1[0];
  rbs[1] = dx_cbs_pyr_oxdna3*e1[1];
  rbs[2] = dx_cbs_pyr_oxdna3*e1[2];
}
template <>
inline void NucleotideOxdna3::base_site<1>(double e1[3], double /*e2*/[3],
  double /*e3*/[3], double rbs[3]) {
  double dx_cbs_pur_oxdna3 = ConstantsOxdna::get_dx_cbs_pur_oxdna3();
  rbs[0] = dx_cbs_pur_oxdna3*e1[0];
  rbs[1] = dx_cbs_pur_oxdna3*e1[1];
  rbs[2] = dx_cbs_pur_oxdna3*e1[2];
}
template <>
inline void NucleotideOxdna3::base_site<2>(double e1[3], double /*e2*/[3],
  double /*e3*/[3], double rbs[3]) {
  double dx_cbs_pyr_oxdna3 = ConstantsOxdna::get_dx_cbs_pyr_oxdna3();
  rbs[0] = dx_cbs_pyr_oxdna3*e1[0];
  rbs[1] = dx_cbs_pyr_oxdna3*e1[1];
  rbs[2] = dx_cbs_pyr_oxdna3*e1[2];
}
template <>
inline void NucleotideOxdna3::base_site<3>(double e1[3], double /*e2*/[3],
  double /*e3*/[3], double rbs[3]) {
  double dx_cbs_pur_oxdna3 = ConstantsOxdna::get_dx_cbs_pur_oxdna3();
  rbs[0] = dx_cbs_pur_oxdna3*e1[0];
  rbs[1] = dx_cbs_pur_oxdna3*e1[1];
  rbs[2] = dx_cbs_pur_oxdna3*e1[2];
}

/* ----------------------------------------------------------------------
   oxRNA2 nucleotide
------------------------------------------------------------------------- */
class NucleotideOxrna2 : public NucleotideOxdna<NucleotideOxrna2> {
 public:
  inline void backbone_site(double e1[3], double /*e2*/[3], double e3[3], double rbk[3]) {
    double dx_cbk_oxrna2 = ConstantsOxdna::get_dx_cbk_oxrna2();
    double dz_cbk_oxrna2 = ConstantsOxdna::get_dz_cbk_oxrna2();
    rbk[0] = dx_cbk_oxrna2 * e1[0] + dz_cbk_oxrna2 * e3[0];
    rbk[1] = dx_cbk_oxrna2 * e1[1] + dz_cbk_oxrna2 * e3[1];
    rbk[2] = dx_cbk_oxrna2 * e1[2] + dz_cbk_oxrna2 * e3[2];
  }
};

}
#endif
