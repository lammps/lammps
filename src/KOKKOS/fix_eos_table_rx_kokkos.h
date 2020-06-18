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

#ifdef FIX_CLASS

FixStyle(eos/table/rx/kk,FixEOStableRXKokkos<Device>)
FixStyle(eos/table/rx/kk/device,FixEOStableRXKokkos<Device>)
FixStyle(eos/table/rx/kk/host,FixEOStableRXKokkos<Host>)

#else

#ifndef LMP_FIX_EOS_TABLE_RX_KOKKOS_H
#define LMP_FIX_EOS_TABLE_RX_KOKKOS_H

#include "fix_eos_table_rx.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixEOStableRXInit{};
struct TagFixEOStableRXSetup{};
struct TagFixEOStableRXTemperatureLookup{};
struct TagFixEOStableRXTemperatureLookup2{};

template<ExecutionSpace Space>
class FixEOStableRXKokkos : public FixEOStableRX {
 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  typedef EV_FLOAT value_type;
  typedef typename GetFloatType<Space>::type SPACE_FLOAT;

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
  void energy_lookup(int, SPACE_FLOAT, SPACE_FLOAT &) const;

  KOKKOS_INLINE_FUNCTION
  void temperature_lookup(int, SPACE_FLOAT, SPACE_FLOAT &) const;

 protected:
  //struct Table {
  //  int ninput;
  //  KK_FLOAT lo,hi;
  //  KK_FLOAT *rfile,*efile;
  //  KK_FLOAT *e2file;
  //  KK_FLOAT delta,invdelta,deltasq6;
  //  KK_FLOAT *r,*e,*de,*e2;
  //};
  //Table *tables, *tables2;

  struct TableDual {
    DAT::tdual_int_1d k_lo,k_hi;
    DAT::tdual_float_1d k_invdelta;
    DAT::tdual_float_2d k_r,k_e,k_de;
  };

  /*struct TableDeviceConst {
    typename AT::t_int_1d_randomread lo,hi;
    typename AT::t_float_1d_randomread invdelta;
    typename AT::t_float_2d_randomread r,e,de;
  };*/
 //Its faster not to use texture fetch if the number of tables is less than 32!
  struct TableDeviceConst {
    typename AT::t_int_1d lo,hi;
    typename AT::t_float_1d invdelta;
    typename AT::t_float_2d_randomread r,e,de;
  };

  struct TableDevice {
    typename AT::t_int_1d lo,hi;
    typename AT::t_float_1d invdelta;
    typename AT::t_float_2d r,e,de;

    void copy(const TableDual *rhs) {
      lo = DualViewHelper<Space>::view(rhs->k_lo);
      hi = DualViewHelper<Space>::view(rhs->k_hi);
      invdelta = DualViewHelper<Space>::view(rhs->k_invdelta);
      r = DualViewHelper<Space>::view(rhs->k_r);
      e = DualViewHelper<Space>::view(rhs->k_e);
      de = DualViewHelper<Space>::view(rhs->k_de);
    }
  };

  struct TableHost {
    HAT::t_int_1d lo,hi;
    HAT::t_float_1d invdelta;
    HAT::t_float_2d r,e,de;

    void copy(const TableDual *rhs) {
      lo = rhs->k_lo.h_view;
      hi = rhs->k_hi.h_view;
      invdelta = rhs->k_invdelta.h_view;
      r = rhs->k_r.h_view;
      e = rhs->k_e.h_view;
      de = rhs->k_de.h_view;
    }
  };

  TableDual* k_table;
  TableDeviceConst d_table_const;
  TableDevice* d_table;
  TableHost* h_table;

  int **tabindex;

  KK_FLOAT boltz;

  void allocate();
  void error_check();
  int update_table;
  void create_kokkos_tables();

  DAT::tdual_float_1d k_dHf,k_energyCorr,k_tempCorrCoeff,k_moleculeCorrCoeff;
  typename AT::t_float_1d d_dHf,d_energyCorr,d_tempCorrCoeff,d_moleculeCorrCoeff;

  typename AT::t_int_1d mask;
  typename AT::t_float_1d uCond,uMech,uChem,uCG,uCGnew,rho,dpdTheta,duChem;
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
