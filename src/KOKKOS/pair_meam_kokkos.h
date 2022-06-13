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

#ifdef PAIR_CLASS

PairStyle(meam/c/kk,PairMEAMKokkos<LMPDeviceType>)
PairStyle(meam/c/kk/device,PairMEAMKokkos<LMPDeviceType>)
PairStyle(meam/c/kk/host,PairMEAMKokkos<LMPHostType>)
PairStyle(meam/kk,PairMEAMKokkos<LMPDeviceType>)
PairStyle(meam/kk/device,PairMEAMKokkos<LMPDeviceType>)
PairStyle(meam/kk/host,PairMEAMKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_MEAMC_KOKKOS_H
#define LMP_PAIR_MEAMC_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_meam.h"
#include "neigh_list_kokkos.h"
#include "meam_kokkos.h"


namespace LAMMPS_NS {
struct TagPairMEAMKernelNeighStrip{};
struct TagPairMEAMKernelA{};
struct TagPairMEAMPackForwardComm{};
struct TagPairMEAMUnpackForwardComm{};

template<class DeviceType>
class MEAMKokkos;

template<class DeviceType>
class PairMEAMKokkos : public PairMEAM, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  //typedef EV_FLOAT value_type;

  PairMEAMKokkos(class LAMMPS *);
  virtual ~PairMEAMKokkos();
  void compute(int, int);
  void coeff(int, char **);
  void init_style();

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMEAMPackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMEAMUnpackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMEAMKernelNeighStrip,  const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMEAMKernelA,  const int, int&) const;

  int pack_forward_comm_kokkos(int, DAT::tdual_int_2d, int, DAT::tdual_xfloat_1d&,
                       int, int *);
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_xfloat_1d&);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 protected:
  class MEAMKokkos<DeviceType> *meam_inst_kk;
  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;
  
  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename ArrayTypes<DeviceType>::t_efloat_1d d_eatom;
  typename ArrayTypes<DeviceType>::t_virial_array d_vatom;

  DAT::tdual_int_1d k_offset; 
  HAT::t_int_1d h_offset;
  typename AT::t_int_1d_randomread d_offset;
  
  DAT::tdual_int_1d k_map; 
  typename AT::t_int_1d_randomread d_map;
  typename AT::t_int_1d_randomread d_ilist_half;
  typename AT::t_int_1d_randomread d_numneigh_half;
  typename AT::t_neighbors_2d d_neighbors_half;
  typename AT::t_int_1d_randomread d_numneigh_full;
  typename AT::t_neighbors_2d d_neighbors_full;
  typename AT::t_int_2d d_sendlist;
  typename AT::t_xfloat_1d_um v_buf;
   
  int neighflag,nlocal,nall,eflag,vflag;
  int iswap; 
  int first;

  friend void pair_virial_fdotr_compute<PairMEAMKokkos>(PairMEAMKokkos*);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: MEAM library error %d

A call to the MEAM Fortran library returned an error.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style MEAM requires newton pair on

See the newton command.  This is a restriction to use the MEAM
potential.

E: Cannot open MEAM potential file %s

The specified MEAM potential file cannot be opened.  Check that the
path and name are correct.

E: Incorrect format in MEAM potential file

Incorrect number of words per line in the potential file.

E: Unrecognized lattice type in MEAM file 1

The lattice type in an entry of the MEAM library file is not
valid.

E: Did not find all elements in MEAM library file

The requested elements were not found in the MEAM file.

E: Keyword %s in MEAM parameter file not recognized

Self-explanatory.

E: Unrecognized lattice type in MEAM file 2

The lattice type in an entry of the MEAM parameter file is not
valid.

*/
