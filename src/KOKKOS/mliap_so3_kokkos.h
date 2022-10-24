/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Matt Bettencourt (NVIDIA)
 ------------------------------------------------------------------------- */

#ifndef LMP_MLIAP_SO3_KOKKOS_H
#define LMP_MLIAP_SO3_KOKKOS_H

#include "kokkos_type.h"
#include "pointers.h"

namespace LAMMPS_NS {

template <class DeviceType> class MLIAP_SO3Kokkos : protected Pointers {
 public:
  MLIAP_SO3Kokkos(LAMMPS *, double vrcut, int vlmax, int vnmax, double valpha);
  MLIAP_SO3Kokkos(LAMMPS *_lmp) : Pointers(_lmp){};
  ~MLIAP_SO3Kokkos() override;

  void init();
  double memory_usage();

  using MemoryDeviceType = typename KKDevice<DeviceType>::value;
  using float_1d = Kokkos::View<double *, Kokkos::LayoutRight, MemoryDeviceType>;
  using float_2d = Kokkos::View<double **, Kokkos::LayoutRight, MemoryDeviceType>;
  using float_3d = Kokkos::View<double ***, Kokkos::LayoutRight, MemoryDeviceType>;
  using float_4d = Kokkos::View<double ****, Kokkos::LayoutRight, MemoryDeviceType>;
  using int_1d = Kokkos::View<int *, MemoryDeviceType>;

  int ncoeff;
  float_2d m_plist_r;
  float_3d k_dplist_r;

 private:
  double alloc_init, alloc_arrays;
  int_1d m_ellpl1, m_ellm1;
  float_1d m_pfac, m_Ylms;
  int m_pfac_l1, m_pfac_l2;
  float_1d m_dfac0, m_dfac1, m_dfac2, m_dfac3, m_dfac4, m_dfac5;
  int m_dfac_l1, m_dfac_l2;
  double m_rcut, m_alpha;
  int m_lmax, m_nmax, m_Nmax;
  float_1d m_g_array, m_w;
  float_1d m_rootpq;
  int_1d m_idxu_block, m_idxylm;
  int m_idxu_count, m_idxy_count;
  int m_numYlms;

  float_3d m_clisttot_r, m_clisttot_i;
  float_1d m_rip_array, m_rip_darray;
  float_1d m_sbes_array, m_sbes_darray;
  int m_init_arrays;

  float_2d m_ulist_r, m_ulist_i;

  float_3d m_dYlm_r, m_dYlm_i;
  float_4d m_dclist;

  // arguments for operators.
  int_1d t_numneighs;
  int_1d t_jelems;
  float_1d t_wjelem;
  float_2d t_rij;
  int t_nmax;
  int t_lmax;
  double t_rcut;
  double t_alpha;
  int_1d t_ij;

  static constexpr bigint m_temp_memory_size = 512 * 1024 * 1024;
  int m_chunk_size;

 public:
  void spectrum(int nlocal, DAT::tdual_int_1d numneighs, DAT::tdual_int_1d jelems,
                DAT::tdual_float_1d wjelem, DAT::tdual_float_2d rij, DAT::tdual_int_1d k_ij,
                int nmax, int lmax, double rcut, double alpha, int totaln, int ncoefs);
  struct MLIAPSO3SpectrumTag {};
  KOKKOS_FUNCTION
  void operator()(const MLIAPSO3SpectrumTag &, int ii) const;

  void spectrum_dxdr(int nlocal, DAT::tdual_int_1d numneighs, DAT::tdual_int_1d jelems,
                     DAT::tdual_float_1d wjelem, DAT::tdual_float_2d rij, DAT::tdual_int_1d k_ij,
                     int nmax, int lmax, double rcut, double alpha, bigint npairs, int ncoefs);
  struct MLIAPSO3SpectrumDXDRTag {};
  KOKKOS_FUNCTION
  void operator()(const MLIAPSO3SpectrumDXDRTag &, int ii) const;

  KOKKOS_FUNCTION
  double Cosine(double Rij, double Rc) const;
  KOKKOS_FUNCTION
  double CosinePrime(double Rij, double Rc) const;
  KOKKOS_FUNCTION
  double compute_sfac(double r, double rcut) const;
  KOKKOS_FUNCTION
  double compute_dsfac(double r, double rcut) const;

  struct MLIAPSO3GetSBESArrayTag {};
  KOKKOS_FUNCTION
  void operator()(const MLIAPSO3GetSBESArrayTag &, int ii) const;

  struct MLIAPSO3GetRipArrayTag {};
  KOKKOS_FUNCTION
  void operator()(const MLIAPSO3GetRipArrayTag &, int ii) const;

  void init_arrays(int nlocal, int ncoefs);
  void init_garray(int nmax, int lmax, double rcut, double alpha, double *w, int lw1,
                   double *g_array, int lg2);

  template <typename UlistView>
  KOKKOS_FUNCTION void compute_uarray_recursive(double x, double y, double z, double r, int twol,
                                                UlistView ulist_r, UlistView ulist_i,
                                                int_1d idxu_block, float_1d rootpqarray) const;
  void compute_ncoeff();

  KOKKOS_FUNCTION
  static int get_sum(int istart, int iend, int id, int imult);

  double compute_g(double r, int n, int nmax, double rcut, double *w, int lw1);
  double phi(double r, int alpha, double rcut);

  template <typename ViewType>
  KOKKOS_FUNCTION void compute_pi(int nmax, int lmax, ViewType clisttot_r, ViewType clisttot_i,
                                  int lcl2, float_2d plist_r, int indpl) const;

  void compute_W(int nmax, double *arr);
};

}    // namespace LAMMPS_NS

#endif
