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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(table/kk,PairTableKokkos<LMPDeviceType>);
PairStyle(table/kk/device,PairTableKokkos<LMPDeviceType>);
PairStyle(table/kk/host,PairTableKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_TABLE_KOKKOS_H
#define LMP_PAIR_TABLE_KOKKOS_H

#include "pair_table.h"
#include "pair_kokkos.h"
#include "neigh_list_kokkos.h"
#include "atom_kokkos.h"

namespace LAMMPS_NS {

template<class Device,int TABSTYLE>
struct S_TableCompute {
  static constexpr int TabStyle = TABSTYLE;
};

template <class DeviceType, int NEIGHFLAG, int TABSTYLE>
struct PairTableComputeFunctor;

template<class DeviceType>
class PairTableKokkos : public PairTable {
 public:

  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairTableKokkos(class LAMMPS *);
  virtual ~PairTableKokkos();

  virtual void compute(int, int);

  template<int TABSTYLE>
  void compute_style(int, int);

  void settings(int, char **);
  double init_one(int, int);

  void init_style();


 protected:

  /*struct TableDeviceConst {
    typename AT::t_ffloat_2d_randomread cutsq;
    typename AT::t_int_2d_randomread tabindex;
    typename AT::t_int_1d_randomread nshiftbits,nmask;
    typename AT::t_ffloat_1d_randomread innersq,invdelta,deltasq6;
    typename AT::t_ffloat_2d_randomread rsq,drsq,e,de,f,df,e2,f2;
  };*/
 //Its faster not to use texture fetch if the number of tables is less than 32!
  struct TableDeviceConst {
    typename AT::t_ffloat_2d cutsq;
    typename AT::t_int_2d tabindex;
    typename AT::t_int_1d nshiftbits,nmask;
    typename AT::t_ffloat_1d innersq,invdelta,deltasq6;
    typename AT::t_ffloat_2d_randomread rsq,drsq,e,de,f,df,e2,f2;
  };

  struct TableDevice {
    typename AT::t_ffloat_2d cutsq;
    typename AT::t_int_2d tabindex;
    typename AT::t_int_1d nshiftbits,nmask;
    typename AT::t_ffloat_1d innersq,invdelta,deltasq6;
    typename AT::t_ffloat_2d rsq,drsq,e,de,f,df,e2,f2;
  };

  struct TableHost {
    typename ArrayTypes<LMPHostType>::t_ffloat_2d cutsq;
    typename ArrayTypes<LMPHostType>::t_int_2d tabindex;
    typename ArrayTypes<LMPHostType>::t_int_1d nshiftbits,nmask;
    typename ArrayTypes<LMPHostType>::t_ffloat_1d innersq,invdelta,deltasq6;
    typename ArrayTypes<LMPHostType>::t_ffloat_2d rsq,drsq,e,de,f,df,e2,f2;
  };

  TableDeviceConst d_table_const;
  TableDevice* d_table;
  TableHost* h_table;

  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  typename AT::t_ffloat_2d d_cutsq;

  virtual void allocate();
  void compute_table(Table *);

  typename AT::t_x_array_randomread x;
  typename AT::t_x_array_const c_x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

 protected:
  int nlocal,nall,eflag,vflag,neighflag,newton_pair;

  int update_table;
  void create_kokkos_tables();

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_ecoul(const F_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                        const int& /*itype*/, const int& /*jtype*/) const { return 0; }

  friend struct PairComputeFunctor<PairTableKokkos,FULL,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend struct PairComputeFunctor<PairTableKokkos,FULL,false,S_TableCompute<DeviceType,LOOKUP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,false,S_TableCompute<DeviceType,LOOKUP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,LOOKUP> >;

  friend struct PairComputeFunctor<PairTableKokkos,FULL,true,S_TableCompute<DeviceType,LINEAR> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,true,S_TableCompute<DeviceType,LINEAR> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,LINEAR> >;
  friend struct PairComputeFunctor<PairTableKokkos,FULL,false,S_TableCompute<DeviceType,LINEAR> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,false,S_TableCompute<DeviceType,LINEAR> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,LINEAR> >;

  friend struct PairComputeFunctor<PairTableKokkos,FULL,true,S_TableCompute<DeviceType,SPLINE> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,true,S_TableCompute<DeviceType,SPLINE> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,SPLINE> >;
  friend struct PairComputeFunctor<PairTableKokkos,FULL,false,S_TableCompute<DeviceType,SPLINE> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,false,S_TableCompute<DeviceType,SPLINE> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,SPLINE> >;

  friend struct PairComputeFunctor<PairTableKokkos,FULL,true,S_TableCompute<DeviceType,BITMAP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,true,S_TableCompute<DeviceType,BITMAP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,BITMAP> >;
  friend struct PairComputeFunctor<PairTableKokkos,FULL,false,S_TableCompute<DeviceType,BITMAP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,false,S_TableCompute<DeviceType,BITMAP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,BITMAP> >;

  friend void pair_virial_fdotr_compute<PairTableKokkos>(PairTableKokkos*);
};






}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Unknown table style in pair_style command

Style of table is invalid for use with pair_style table command.

E: Illegal number of pair table entries

There must be at least 2 table entries.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot use chosen neighbor list style with lj/cut/kk

That style is not supported by Kokkos.

U: Pair distance < table inner cutoff

Two atoms are closer together than the pairwise table allows.

U: Pair distance > table outer cutoff

Two atoms are further apart than the pairwise table allows.

U: Invalid pair table length

Length of read-in pair table is invalid

U: Invalid pair table cutoff

Cutoffs in pair_coeff command are not valid with read-in pair table.

U: Bitmapped table in file does not match requested table

Setting for bitmapped table in pair_coeff command must match table
in file exactly.

U: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

U: Did not find keyword in table file

Keyword used in pair_coeff command was not found in table file.

U: Bitmapped table is incorrect length in table file

Number of table entries is not a correct power of 2.

U: Invalid keyword in pair table parameters

Keyword used in list of table parameters is not recognized.

U: Pair table parameters did not set N

List of pair table parameters must include N setting.

U: Pair table cutoffs must all be equal to use with KSpace

When using pair style table with a long-range KSpace solver, the
cutoffs for all atom type pairs must all be the same, since the
long-range solver starts at that cutoff.

*/
