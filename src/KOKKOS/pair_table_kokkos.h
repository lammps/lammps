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

PairStyle(table/kk,PairTableKokkos<LMPDeviceType>)
PairStyle(table/kk/device,PairTableKokkos<LMPDeviceType>)
PairStyle(table/kk/host,PairTableKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_TABLE_KOKKOS_H
#define LMP_PAIR_TABLE_KOKKOS_H

#include "pair.h"
#include "pair_kokkos.h"
#include "neigh_list_kokkos.h"
#include "atom_kokkos.h"

namespace LAMMPS_NS {

template<class Device,int TABSTYLE>
struct S_TableCompute {
  enum {TabStyle = TABSTYLE};
};

template <class DeviceType, int NEIGHFLAG, int TABSTYLE>
class PairTableComputeFunctor;

template<class DeviceType>
class PairTableKokkos : public Pair {
 public:

  enum {COUL_FLAG=0};
  typedef DeviceType device_type;

  PairTableKokkos(class LAMMPS *);
  virtual ~PairTableKokkos();

  virtual void compute(int, int);
  
  template<int TABSTYLE> 
  void compute_style(int, int);

  /*template<int EVFLAG, int NEIGHFLAG, int NEWTON_PAIR, int TABSTYLE>
  KOKKOS_FUNCTION
  EV_FLOAT compute_item(const int& i,
                        const NeighListKokkos<DeviceType> &list) const;
*/
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

  void init_style();


 protected:
  enum{LOOKUP,LINEAR,SPLINE,BITMAP};

  int tabstyle,tablength;
  /*struct TableDeviceConst {
    typename ArrayTypes<DeviceType>::t_ffloat_2d_randomread cutsq;
    typename ArrayTypes<DeviceType>::t_int_2d_randomread tabindex;
    typename ArrayTypes<DeviceType>::t_int_1d_randomread nshiftbits,nmask;
    typename ArrayTypes<DeviceType>::t_ffloat_1d_randomread innersq,invdelta,deltasq6;
    typename ArrayTypes<DeviceType>::t_ffloat_2d_randomread rsq,drsq,e,de,f,df,e2,f2;
  };*/
 //Its faster not to use texture fetch if the number of tables is less than 32! 
  struct TableDeviceConst {
    typename ArrayTypes<DeviceType>::t_ffloat_2d cutsq;
    typename ArrayTypes<DeviceType>::t_int_2d tabindex;
    typename ArrayTypes<DeviceType>::t_int_1d nshiftbits,nmask;
    typename ArrayTypes<DeviceType>::t_ffloat_1d innersq,invdelta,deltasq6;
    typename ArrayTypes<DeviceType>::t_ffloat_2d_randomread rsq,drsq,e,de,f,df,e2,f2;
  };

  struct TableDevice {
    typename ArrayTypes<DeviceType>::t_ffloat_2d cutsq;
    typename ArrayTypes<DeviceType>::t_int_2d tabindex;
    typename ArrayTypes<DeviceType>::t_int_1d nshiftbits,nmask;
    typename ArrayTypes<DeviceType>::t_ffloat_1d innersq,invdelta,deltasq6;
    typename ArrayTypes<DeviceType>::t_ffloat_2d rsq,drsq,e,de,f,df,e2,f2;
  };

  struct TableHost {
    typename ArrayTypes<LMPHostType>::t_ffloat_2d cutsq;
    typename ArrayTypes<LMPHostType>::t_int_2d tabindex;
    typename ArrayTypes<LMPHostType>::t_int_1d nshiftbits,nmask;
    typename ArrayTypes<LMPHostType>::t_ffloat_1d innersq,invdelta,deltasq6;
    typename ArrayTypes<LMPHostType>::t_ffloat_2d rsq,drsq,e,de,f,df,e2,f2;
  };

  struct Table {
    int ninput,rflag,fpflag,match,ntablebits;
    int nshiftbits,nmask;
    double rlo,rhi,fplo,fphi,cut;
    double *rfile,*efile,*ffile;
    double *e2file,*f2file;
    double innersq,delta,invdelta,deltasq6;
    double *rsq,*drsq,*e,*de,*f,*df,*e2,*f2;
  };
  int ntables;
  Table *tables;
  TableDeviceConst d_table_const;
  TableDevice* d_table;
  TableHost* h_table;

  int **tabindex;
  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  typename ArrayTypes<DeviceType>::t_ffloat_2d d_cutsq;

  void allocate();
  void read_table(Table *, char *, char *);
  void param_extract(Table *, char *);
  void bcast_table(Table *);
  void spline_table(Table *);
  void compute_table(Table *);
  void null_table(Table *);
  void free_table(Table *);
  void spline(double *, double *, int, double, double, double *);
  double splint(double *, double *, double *, int, double);

  typename ArrayTypes<DeviceType>::t_x_array_randomread x;
  typename ArrayTypes<DeviceType>::t_x_array_const c_x;
  typename ArrayTypes<DeviceType>::t_f_array f;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread type;
  typename ArrayTypes<DeviceType>::t_efloat_1d d_eatom;
  typename ArrayTypes<DeviceType>::t_virial_array d_vatom;

 protected:
  int nlocal,nall,eflag,vflag,neighflag,newton_pair;
  class AtomKokkos *atomKK;
  int update_table;
  void create_kokkos_tables();
  void cleanup_copy();

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_ecoul(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
    return 0;
  }

  friend class PairComputeFunctor<PairTableKokkos,FULL,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableKokkos,HALF,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableKokkos,N2,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableKokkos,FULLCLUSTER,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableKokkos,FULL,false,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableKokkos,HALF,false,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableKokkos,N2,false,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableKokkos,FULLCLUSTER,false,S_TableCompute<DeviceType,LOOKUP> >;

  friend class PairComputeFunctor<PairTableKokkos,FULL,true,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableKokkos,HALF,true,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableKokkos,N2,true,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableKokkos,FULLCLUSTER,true,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableKokkos,FULL,false,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableKokkos,HALF,false,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableKokkos,N2,false,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableKokkos,FULLCLUSTER,false,S_TableCompute<DeviceType,LINEAR> >;

  friend class PairComputeFunctor<PairTableKokkos,FULL,true,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableKokkos,HALF,true,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableKokkos,N2,true,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableKokkos,FULLCLUSTER,true,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableKokkos,FULL,false,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableKokkos,HALF,false,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableKokkos,N2,false,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableKokkos,FULLCLUSTER,false,S_TableCompute<DeviceType,SPLINE> >;

  friend class PairComputeFunctor<PairTableKokkos,FULL,true,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableKokkos,HALF,true,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableKokkos,N2,true,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableKokkos,FULLCLUSTER,true,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableKokkos,FULL,false,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableKokkos,HALF,false,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableKokkos,N2,false,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableKokkos,FULLCLUSTER,false,S_TableCompute<DeviceType,BITMAP> >;
/*template<int FULL_NEIGH>
  KOKKOS_INLINE_FUNCTION
    void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
                  const F_FLOAT &fpair, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const;
*/
};
/*
template <class DeviceType, int NEIGHFLAG, int TABSTYLE>
struct PairTableComputeFunctor  {
  typedef DeviceType device_type ;
  typedef EV_FLOAT value_type;

  PairTableKokkos<DeviceType> c;
  NeighListKokkos<DeviceType> list;

  PairTableComputeFunctor(PairTableKokkos<DeviceType>* c_ptr,
                          NeighListKokkos<DeviceType>* list_ptr):
  c(*c_ptr),list(*list_ptr) {};
  ~PairTableComputeFunctor() {c.cleanup_copy();list.clean_copy();};

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    if (c.newton_pair) c.template compute_item<0,NEIGHFLAG,1,TABSTYLE>(i,list);
    else c.template compute_item<0,NEIGHFLAG,0,TABSTYLE>(i,list);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type &energy_virial) const {
    if (c.newton_pair)
      energy_virial += c.template compute_item<1,NEIGHFLAG,1,TABSTYLE>(i,list);
    else
      energy_virial += c.template compute_item<1,NEIGHFLAG,0,TABSTYLE>(i,list);
  }

  KOKKOS_INLINE_FUNCTION
  static void init(volatile value_type &update) {
    update.evdwl = 0;
    update.ecoul = 0;
    update.v[0] = 0;
    update.v[1] = 0;
    update.v[2] = 0;
    update.v[3] = 0;
    update.v[4] = 0;
    update.v[5] = 0;
  }
  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type &update,
                   const volatile value_type &source) {
    update.evdwl += source.evdwl;
    update.ecoul += source.ecoul;
    update.v[0] += source.v[0];
    update.v[1] += source.v[1];
    update.v[2] += source.v[2];
    update.v[3] += source.v[3];
    update.v[4] += source.v[4];
    update.v[5] += source.v[5];
  }
};

*/





}

#endif
#endif

/* ERROR/WARNING messages:

E: Pair distance < table inner cutoff

Two atoms are closer together than the pairwise table allows.

E: Pair distance > table outer cutoff

Two atoms are further apart than the pairwise table allows.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Unknown table style in pair_style command

Style of table is invalid for use with pair_style table command.

E: Illegal number of pair table entries

There must be at least 2 table entries.

E: Invalid pair table length

Length of read-in pair table is invalid

E: Invalid pair table cutoff

Cutoffs in pair_coeff command are not valid with read-in pair table.

E: Bitmapped table in file does not match requested table

Setting for bitmapped table in pair_coeff command must match table
in file exactly.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

E: Did not find keyword in table file

Keyword used in pair_coeff command was not found in table file.

E: Bitmapped table is incorrect length in table file

Number of table entries is not a correct power of 2.

E: Invalid keyword in pair table parameters

Keyword used in list of table parameters is not recognized.

E: Pair table parameters did not set N

List of pair table parameters must include N setting.

E: Pair table cutoffs must all be equal to use with KSpace

When using pair style table with a long-range KSpace solver, the
cutoffs for all atom type pairs must all be the same, since the
long-range solver starts at that cutoff.

*/
