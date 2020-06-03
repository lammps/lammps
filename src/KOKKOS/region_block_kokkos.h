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

#ifdef REGION_CLASS

RegionStyle(block/kk,RegBlockKokkos<Device>)
RegionStyle(block/kk/device,RegBlockKokkos<Device>)
RegionStyle(block/kk/host,RegBlockKokkos<Host>)

#else

#ifndef LMP_REGION_BLOCK_KOKKOS_H
#define LMP_REGION_BLOCK_KOKKOS_H

#include "region_block.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagRegBlockMatchAll{};

template<ExecutionSpace Space>
class RegBlockKokkos : public RegBlock, public KokkosBase {
  friend class FixPour;

 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  RegBlockKokkos(class LAMMPS *, int, char **);
  ~RegBlockKokkos();
  void match_all_kokkos(int, DAT::tdual_int_1d);

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRegBlockMatchAll, const int&) const;

 private:
  int groupbit;
  typename AT::t_int_1d d_match;

  typename AT::t_float_1d_3_randomread x;
  typename AT::t_int_1d_randomread mask;

  KOKKOS_INLINE_FUNCTION
  int k_inside(KK_FLOAT, KK_FLOAT, KK_FLOAT) const;
  KOKKOS_INLINE_FUNCTION
  int match(KK_FLOAT, KK_FLOAT, KK_FLOAT) const;
  KOKKOS_INLINE_FUNCTION
  void inverse_transform(KK_FLOAT &, KK_FLOAT &, KK_FLOAT &) const;
  KOKKOS_INLINE_FUNCTION
  void rotate(KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT) const;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot use region INF or EDGE when box does not exist

Regions that extend to the box boundaries can only be used after the
create_box command has been used.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
