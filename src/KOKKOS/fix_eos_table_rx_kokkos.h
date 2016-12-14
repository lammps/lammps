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

FixStyle(eos/table/rx/kk,FixEOStableRXKokkos<LMPDeviceType>)
FixStyle(eos/table/rx/kk/device,FixEOStableRXKokkos<LMPDeviceType>)
FixStyle(eos/table/rx/kk/host,FixEOStableRXKokkos<LMPHostType>)

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

  void allocate();
  void error_check();

  //double *dHf;

  typename AT::t_int_1d mask;
  typename AT::t_efloat_1d uCond,uMech,uChem,uCG,uCGnew,rho,dpdTheta,duChem;

  DAT::tdual_int_scalar k_error_flag;
  DAT::tdual_int_scalar k_warning_flag;

  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_forward_comm(int , int *, double *, int, int *);
  void unpack_forward_comm(int , int , double *);

  //int *eosSpecies;
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

The maximum number of interations was exceeded in the secant solver

*/
