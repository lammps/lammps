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

#ifdef FIX_CLASS
// clang-format off
FixStyle(eos/table/rx/kk,FixEOStableRXKokkos<LMPDeviceType>);
FixStyle(eos/table/rx/kk/device,FixEOStableRXKokkos<LMPDeviceType>);
FixStyle(eos/table/rx/kk/host,FixEOStableRXKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_EOS_TABLE_RX_KOKKOS_H
#define LMP_FIX_EOS_TABLE_RX_KOKKOS_H

#include "fix_eos_table_rx.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixEOStableRXInit{};
struct TagFixEOStableRXSetup{};
struct TagFixEOStableRXTemperatureLookup{};
struct TagFixEOStableRXTemperatureLookup2{};

template<class DeviceType>
class FixEOStableRXKokkos : public FixEOStableRX {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  FixEOStableRXKokkos(class LAMMPS *, int, char **);
  virtual ~FixEOStableRXKokkos();
  void setup(int);
  void init();
  void post_integrate();
  void end_of_step();

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEOStableRXInit, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEOStableRXSetup, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEOStableRXTemperatureLookup, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEOStableRXTemperatureLookup2, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void energy_lookup(int, double, double &) const;

  KOKKOS_INLINE_FUNCTION
  void temperature_lookup(int, double, double &) const;

 protected:
  //struct Table {
  //  int ninput;
  //  double lo,hi;
  //  double *rfile,*efile;
  //  double *e2file;
  //  double delta,invdelta,deltasq6;
  //  double *r,*e,*de,*e2;
  //};
  //Table *tables, *tables2;

  /*struct TableDeviceConst {
    typename ArrayTypes<DeviceType>::t_int_1d_randomread lo,hi;
    typename ArrayTypes<DeviceType>::t_ffloat_1d_randomread invdelta;
    typename ArrayTypes<DeviceType>::t_ffloat_2d_randomread r,e,de;
  };*/
 //Its faster not to use texture fetch if the number of tables is less than 32!
  struct TableDeviceConst {
    typename ArrayTypes<DeviceType>::t_int_1d lo,hi;
    typename ArrayTypes<DeviceType>::t_ffloat_1d invdelta;
    typename ArrayTypes<DeviceType>::t_ffloat_2d_randomread r,e,de;
  };

  struct TableDevice {
    typename ArrayTypes<DeviceType>::t_int_1d lo,hi;
    typename ArrayTypes<DeviceType>::t_ffloat_1d invdelta;
    typename ArrayTypes<DeviceType>::t_ffloat_2d r,e,de;
  };

  struct TableHost {
    typename ArrayTypes<LMPHostType>::t_int_1d lo,hi;
    typename ArrayTypes<LMPHostType>::t_ffloat_1d invdelta;
    typename ArrayTypes<LMPHostType>::t_ffloat_2d r,e,de;
  };

  TableDeviceConst d_table_const;
  TableDevice* d_table;
  TableHost* h_table;

  int **tabindex;

  double boltz;

  void allocate();
  void error_check();
  int update_table;
  void create_kokkos_tables();

  DAT::tdual_float_1d k_dHf,k_energyCorr,k_tempCorrCoeff,k_moleculeCorrCoeff;
  typename AT::t_float_1d d_dHf,d_energyCorr,d_tempCorrCoeff,d_moleculeCorrCoeff;

  typename AT::t_int_1d mask;
  typename AT::t_efloat_1d uCond,uMech,uChem,uCG,uCGnew,rho,dpdTheta,duChem;
  typename AT::t_float_2d dvector;

  DAT::tdual_int_scalar k_error_flag;
  DAT::tdual_int_scalar k_warning_flag;

  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_forward_comm(int , int *, double *, int, int *);
  void unpack_forward_comm(int , int , double *);

  };
}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: FixEOStableRXKokkos requires a fix rx command.

The fix rx command must come before the pair style command in the input file

E:  There are no rx species specified

There must be at least one species specified through the fix rx command

E:  Invalid eos/table/rx length

The eos/table/rx table must have more than one entry.

E:  eos/table/rx values are not increasing

The equation-of-state must an increasing function

E:  FixEOStableRX requires atom_style with internal temperature and energies (e.g. dpd)

Self-explanatory.

E:  Internal temperature <= zero.

Self-explanatory.

E:  Cannot open eos table/rx potential file %s

Self-explanatory.

E:  Incorrect format in eos table/rx file

Self-explanatory.

E:  Cannot open file %s

Self-explanatory.

E:  Did not find keyword in table file

Self-explanatory.

E:  Illegal fix eos/table/rx command

Incorrect number of arguments specified for the fix eos/table/rx command.

E:  Invalid keyword in fix eos/table/rx parameters

Self-explanatory.

E:  The number of columns in fix eos/table/rx does not match the number of species.

Self-explanatory.  Check format for fix eos/table/rx file.

E:  fix eos/table/rx parameters did not set N

The number of table entries was not set in the eos/table/rx file

W:  Secant solver did not converge because table bounds were exceeded

The secant solver failed to converge, resulting in the lower or upper table bound temperature to be returned

E: NaN detected in secant solver.

Self-explanatory.

E: Maxit exceeded in secant solver

The maximum number of iterations was exceeded in the secant solver

*/
