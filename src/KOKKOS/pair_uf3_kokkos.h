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

/* ----------------------------------------------------------------------
   Contributing author: Ajinkya Hire (Univ. of Florida),
                        Hendrik Kraß (Univ. of Constance),
                        Matthias Rupp (Luxembourg Institute of Science and Technology),
                        Richard Hennig (Univ of Florida)
---------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(uf3/kk,PairUF3Kokkos<LMPDeviceType>)
PairStyle(uf3/kk/device,PairUF3Kokkos<LMPDeviceType>)
// clang-format on
#else

#ifndef LMP_PAIR_UF3_KOKKOS_H
#define LMP_PAIR_UF3_KOKKOS_H

#include "kokkos.h"
#include "pair_kokkos.h"
#include "pair_uf3.h"

template <int NEIGHFLAG, int EVFLAG> struct TagPairUF3ComputeFullA {};
struct TagPairUF3ComputeShortNeigh {};

namespace LAMMPS_NS {

template <class DeviceType> class PairUF3Kokkos : public PairUF3 {
 public:
  PairUF3Kokkos(class LAMMPS *);
  ~PairUF3Kokkos() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void allocate();
  void init_style() override;
  void init_list(int, class NeighList *) override;    // needed for ptr to full neigh list
  double init_one(int, int) override;                 // needed for cutoff radius for neighbour list
  double single(int, int, int, int, double, double, double, double &) override;

  template <typename T, typename V> void copy_2d(V &d, T **h, int m, int n);
  template <typename T, typename V> void copy_3d(V &d, T ***h, int m, int n, int o);

  template <int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION void operator()(TagPairUF3ComputeFullA<NEIGHFLAG, EVFLAG>, const int &,
                                         EV_FLOAT &) const;

  template <int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION void operator()(TagPairUF3ComputeFullA<NEIGHFLAG, EVFLAG>,
                                         const int &) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairUF3ComputeShortNeigh, const int &) const;

  enum { EnabledNeighFlags = FULL };
  enum { COUL_FLAG = 0 };
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

 protected:
  typename AT::tdual_ffloat_2d k_cutsq;//Create a DualView, defination of tdual_ffloat_2d in kokkos_type.h
  typename AT::t_ffloat_2d d_cutsq; //t_ffloat_2d = t_dev ==> Creates a new View d_cutsq
  //the type of d_cutsq is decided by the Device(not host) type for the DualView k_cutsq
  //Meaning the memory location of d_cutsq is the same as the Device(not host) memory location of
  //k_cutsq
  typedef Kokkos::DualView<F_FLOAT***, Kokkos::LayoutRight, DeviceType> tdual_ffloat_3d;
  typedef Kokkos::DualView<F_FLOAT****, Kokkos::LayoutRight, DeviceType> tdual_ffloat_4d;
  tdual_ffloat_3d k_cut_3b;
  tdual_ffloat_4d k_min_cut_3b;
  typename tdual_ffloat_3d::t_dev d_cut_3b;
  typename tdual_ffloat_4d::t_dev d_min_cut_3b;
  template <typename TYPE> void destroy_3d(TYPE data, typename TYPE::value_type*** &array);
  template <typename TYPE> void destroy_4d(TYPE data, typename TYPE::value_type**** &array);
  Kokkos::View<F_FLOAT **, LMPDeviceType::array_layout, LMPDeviceType> /*d_cutsq,*/ d_cut_3b_list;
  //Kokkos::View<F_FLOAT ***, LMPDeviceType::array_layout, LMPDeviceType> d_cut_3b;

  Kokkos::View<F_FLOAT **, LMPDeviceType::array_layout, LMPDeviceType> d_coefficients_2b;
  Kokkos::View<F_FLOAT **, LMPDeviceType::array_layout, LMPDeviceType> d_dncoefficients_2b;
  Kokkos::View<F_FLOAT **, LMPDeviceType::array_layout, LMPDeviceType> d_n2b_knot;
  Kokkos::View<F_FLOAT *, LMPDeviceType::array_layout, LMPDeviceType> d_n2b_knot_spacings;
  Kokkos::View<int **, LMPDeviceType::array_layout, LMPDeviceType> map2b;
  Kokkos::View<F_FLOAT[4][4], LMPDeviceType::array_layout, LMPDeviceType> constants;
  Kokkos::View<F_FLOAT[3][3], LMPDeviceType::array_layout, LMPDeviceType> dnconstants;
  Kokkos::View<F_FLOAT ***, LMPDeviceType::array_layout, LMPDeviceType> d_n3b_knot_matrix;
  Kokkos::View<F_FLOAT ****, LMPDeviceType::array_layout, LMPDeviceType> d_coefficients_3b;
  Kokkos::View<F_FLOAT *****, LMPDeviceType::array_layout, LMPDeviceType> d_dncoefficients_3b;
  Kokkos::View<F_FLOAT **, LMPDeviceType::array_layout, LMPDeviceType> d_n3b_knot_spacings;
  Kokkos::View<F_FLOAT **, LMPDeviceType::array_layout, LMPDeviceType> d_n3b_knot_matrix_spacings;
  Kokkos::View<int ***, LMPDeviceType::array_layout, LMPDeviceType> map3b;

  Kokkos::View<F_FLOAT **[16], LMPDeviceType::array_layout, LMPDeviceType> constants_2b;
  Kokkos::View<F_FLOAT **[9], LMPDeviceType::array_layout, LMPDeviceType> dnconstants_2b;
  Kokkos::View<F_FLOAT ***[16], LMPDeviceType::array_layout, LMPDeviceType> constants_3b;
  Kokkos::View<F_FLOAT ***[9], LMPDeviceType::array_layout, LMPDeviceType> dnconstants_3b;

  std::vector<F_FLOAT> get_constants(double *knots, double coefficient);
  std::vector<F_FLOAT> get_dnconstants(double *knots, double coefficient);

  int coefficients_created = 0;
  void create_coefficients();
  void create_3b_coefficients();
  void create_2b_coefficients();
  std::vector<F_FLOAT> get_coefficients(const double *knots, const double coefficient) const;
  std::vector<F_FLOAT> get_dncoefficients(const double *knots, const double coefficient) const;

  template <int EVFLAG>
  void twobody(const int itype, const int jtype, const F_FLOAT r, F_FLOAT &evdwl,
               F_FLOAT &fpair) const;
  template <int EVFLAG>
  void threebody(const int itype, const int jtype, const int ktype, const F_FLOAT value_rij,
                 const F_FLOAT value_rik, const F_FLOAT value_rjk, F_FLOAT &evdwl3,
                 F_FLOAT (&fforce)[3]) const;

  template <int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION void
  ev_tally(EV_FLOAT &ev, const int &i, const int &j, const F_FLOAT &epair, const F_FLOAT &fpair,
           const F_FLOAT &delx, const F_FLOAT &dely, const F_FLOAT &delz) const;

  template <int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION void ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
                                        const F_FLOAT &evdwl, const F_FLOAT &ecoul, F_FLOAT *fj,
                                        F_FLOAT *fk, F_FLOAT *drji, F_FLOAT *drki) const;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_tagint_1d tag;
  typename AT::t_int_1d_randomread type;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  using ScatterFType = Kokkos::Experimental::ScatterView<F_FLOAT *[3], Kokkos::LayoutRight,
                                                         typename DeviceType::memory_space>;
  ScatterFType fscatter;
  using ScatterVType = Kokkos::Experimental::ScatterView<F_FLOAT *[6], Kokkos::LayoutRight,
                                                         typename DeviceType::memory_space>;
  ScatterVType vscatter;
  using ScatterCVType = Kokkos::Experimental::ScatterView<F_FLOAT *[9], Kokkos::LayoutRight,
                                                          typename DeviceType::memory_space>;
  ScatterCVType cvscatter;
  using ScatterEType = Kokkos::Experimental::ScatterView<E_FLOAT *, LMPDeviceType::array_layout,
                                                         typename DeviceType::memory_space>;
  ScatterEType escatter;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  int neighflag, newton_pair;
  int nlocal, nall, eflag, vflag;

  int inum;
  Kokkos::View<int **, DeviceType> d_neighbors_short;
  Kokkos::View<int *, DeviceType> d_numneigh_short;

  friend void pair_virial_fdotr_compute<PairUF3Kokkos>(PairUF3Kokkos *);
};

KOKKOS_INLINE_FUNCTION int min(int i, int j)
{
  return i < j ? i : j;
}
KOKKOS_INLINE_FUNCTION int max(int i, int j)
{
  return i > j ? i : j;
}

}    // namespace LAMMPS_NS

#endif
#endif

