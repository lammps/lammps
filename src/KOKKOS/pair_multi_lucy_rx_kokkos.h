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

#ifdef PAIR_CLASS

PairStyle(multi/lucy/rx/kk,PairMultiLucyRXKokkos<Device>)
PairStyle(multi/lucy/rx/kk/device,PairMultiLucyRXKokkos<Device>)
PairStyle(multi/lucy/rx/kk/host,PairMultiLucyRXKokkos<Host>)

#else

#ifndef LMP_PAIR_MULTI_LUCY_RX_KOKKOS_H
#define LMP_PAIR_MULTI_LUCY_RX_KOKKOS_H


#include "pair_multi_lucy_rx.h"
#include "pair_kokkos.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagPairMultiLucyRXPackForwardComm{};
struct TagPairMultiLucyRXUnpackForwardComm{};

struct TagPairMultiLucyRXgetMixingWeights{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int TABSTYLE>
struct TagPairMultiLucyRXCompute{};

struct TagPairMultiLucyRXZero{};

template<int NEIGHFLAG, int NEWTON_PAIR, bool ONE_TYPE>
struct TagPairMultiLucyRXComputeLocalDensity{};

template<ExecutionSpace Space>
class PairMultiLucyRXKokkos : public PairMultiLucyRX, public KokkosBase {
 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  typedef typename GetFloatType<Space>::type SPACE_FLOAT;
  typedef EV_FLOAT value_type;

  PairMultiLucyRXKokkos(class LAMMPS *);
  virtual ~PairMultiLucyRXKokkos();

  void compute(int, int);
  void settings(int, char **);

  template<int TABSTYLE>
  void compute_style(int, int);

  void init_style();
  int pack_forward_comm_kokkos(int, DAT::tdual_int_2d, int, DAT::tdual_float_1d&,
                               int, int *);
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_float_1d&);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  void computeLocalDensity();

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXPackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXUnpackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXgetMixingWeights, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int TABSTYLE>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG,TABSTYLE>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int TABSTYLE>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG,TABSTYLE>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXZero, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, bool ONE_TYPE>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXComputeLocalDensity<NEIGHFLAG,NEWTON_PAIR,ONE_TYPE>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

 private:
  int nlocal;
  int neighflag;
  int eflag,vflag;

  KK_FLOAT cutsq_type11;
  KK_FLOAT rcut_type11;
  KK_FLOAT factor_type11;

  enum{LOOKUP,LINEAR,SPLINE,BITMAP};

  struct TableDual {
    DAT::tdual_int_2d k_tabindex;
    DAT::tdual_float_1d k_innersq,k_invdelta;
    DAT::tdual_float_2d k_rsq,k_e,k_de,k_f,k_df;
  };

  //struct Table {
  //  int ninput,rflag,fpflag,match;
  //  KK_FLOAT rlo,rhi,fplo,fphi,cut;
  //  KK_FLOAT *rfile,*efile,*ffile;
  //  KK_FLOAT *e2file,*f2file;
  //  KK_FLOAT innersq,delta,invdelta,deltasq6;
  //  KK_FLOAT *rsq,*drsq,*e,*de,*f,*df,*e2,*f2;
  //};

  /*struct TableDeviceConst {
    typename AT::t_int_2d_randomread tabindex;
    typename AT::t_float_1d_randomread innersq,invdelta;
    typename AT::t_float_2d_randomread rsq,e,de,f,df;
  };*/
 //Its faster not to use texture fetch if the number of tables is less than 32!
  struct TableDeviceConst {
    typename AT::t_int_2d tabindex;
    typename AT::t_float_1d innersq,invdelta;
    typename AT::t_float_2d_randomread rsq,e,de,f,df;
  };

  struct TableDevice {
    typename AT::t_int_2d tabindex;
    typename AT::t_float_1d innersq,invdelta;
    typename AT::t_float_2d rsq,e,de,f,df;

    TableDevice(const TableDual *rhs) {
      tabindex = DualViewHelper<Space>::view(rhs->k_tabindex);
      innersq = DualViewHelper<Space>::view(rhs->k_innersq);
      invdelta = DualViewHelper<Space>::view(rhs->k_invdelta);
      rsq = DualViewHelper<Space>::view(rhs->k_rsq);
      e = DualViewHelper<Space>::view(rhs->k_e);
      de = DualViewHelper<Space>::view(rhs->k_de);
      f = DualViewHelper<Space>::view(rhs->k_f);
      df = DualViewHelper<Space>::view(rhs->k_df);
    }
  };

  struct TableHost {
    HAT::t_int_2d tabindex;
    HAT::t_float_1d innersq,invdelta;
    HAT::t_float_2d rsq,e,de,f,df;

    TableHost(const TableDual *rhs) {
      tabindex = rhs->k_tabindex.h_view;
      innersq = rhs->k_innersq.h_view;
      invdelta = rhs->k_invdelta.h_view;
      rsq = rhs->k_rsq.h_view;
      e = rhs->k_e.h_view;
      de = rhs->k_de.h_view;
      f = rhs->k_f.h_view;
      df = rhs->k_df.h_view;
    }
  };

  TableDual* k_table;
  TableDeviceConst d_table_const;
  TableDevice* d_table;
  TableHost* h_table;

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  void allocate();
  int update_table;
  void create_kokkos_tables();

  KOKKOS_INLINE_FUNCTION
  void getMixingWeights(int, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &) const;

  typename AT::t_float_1d d_mixWtSite1old,d_mixWtSite2old,d_mixWtSite1,d_mixWtSite2;

  typename AT::t_float_1d_3_randomread x;
  typename AT::t_float_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_float_1d rho;
  typename HAT::t_float_1d h_rho;
  typename AT::t_float_1d uCG, uCGnew;
  typename AT::t_float_2d dvector;

  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  typename AT::t_float_1d d_eatom;
  typename AT::t_float_1d_6 d_vatom;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::tdual_int_scalar k_error_flag;

  DAT::tdual_float_2d k_cutsq;
  typename AT::t_float_2d d_cutsq;

  int iswap;
  int first;
  typename AT::t_int_2d d_sendlist;
  typename AT::t_float_1d_um v_buf;

  friend void pair_virial_fdotr_compute<Space,PairMultiLucyRXKokkos>(PairMultiLucyRXKokkos*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Pair multi/lucy/rx command requires atom_style with density (e.g. dpd, meso)

Self-explanatory

E: Density < table inner cutoff

The local density inner is smaller than the inner cutoff

E: Density > table inner cutoff

The local density inner is greater than the inner cutoff

E: Only LOOKUP and LINEAR table styles have been implemented for pair multi/lucy/rx

Self-explanatory

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E:  Unknown table style in pair_style command

Self-explanatory

E: Illegal number of pair table entries

There must be at least 2 table entries.

E: Illegal pair_coeff command

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: PairMultiLucyRXKokkos requires a fix rx command

The fix rx command must come before the pair style command in the input file

E:  There are no rx species specified

There must be at least one species specified through the fix rx command

E: Invalid pair table length

Length of read-in pair table is invalid

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Did not find keyword in table file

Keyword used in pair_coeff command was not found in table file.

E: Invalid keyword in pair table parameters

Keyword used in list of table parameters is not recognized.

E: Pair table parameters did not set N

List of pair table parameters must include N setting.

*/
