#ifndef LMP_MEAMKOKKOS_H
#define LMP_MEAMKOKKOS_H

#include "meam.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "kokkos.h"
#include <cmath>
#include <cstdlib>

namespace LAMMPS_NS {

struct TagMEAMDensFinal{};
template<int NEIGHFLAG>
struct TagMEAMDensInit{};
struct TagMEAMInitialize{};
template<int NEIGHFLAG>
struct TagMEAMforce{};

template<class DeviceType>
class MEAMKokkos : public MEAM
{
public:
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  MEAMKokkos(Memory* mem);
  ~MEAMKokkos();

  KOKKOS_INLINE_FUNCTION
  void operator()(TagMEAMDensFinal, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagMEAMDensInit<NEIGHFLAG>, const int&, EV_FLOAT&) const;
  
  KOKKOS_INLINE_FUNCTION
  void operator()(TagMEAMInitialize, const int&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagMEAMforce<NEIGHFLAG>, const int&, EV_FLOAT&) const;
protected:
//Parameters to meam_dens_init - is there a better way to do this?
  int ntype;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread d_offset;
  typename AT::t_int_1d_randomread fmap;
  typename AT::t_x_array_randomread x;
  typename AT::t_int_1d_randomread d_numneigh_half;
  typename AT::t_int_1d_randomread d_numneigh_full;
  typename AT::t_neighbors_2d d_neighbors_half;
  typename AT::t_neighbors_2d d_neighbors_full;
  typename AT::t_int_1d_randomread d_ilist_half;
  typename AT::t_f_array f;
  typename ArrayTypes<DeviceType>::t_virial_array d_vatom;
//Parameters to meam_dens_final - is there a better way to do this?
  int eflag_either, eflag_global, eflag_atom, vflag_atom;
  double *eng_vdwl;
  typename ArrayTypes<DeviceType>::t_efloat_1d d_eatom;

public:
  void meam_dens_setup(int, int, int);
  void meam_setup_done(void);
  void meam_dens_init(int , int , typename AT::t_int_1d_randomread , typename AT::t_int_1d_randomread, typename AT::t_x_array_randomread, typename AT::t_int_1d_randomread, 
                      typename AT::t_int_1d_randomread , int* , typename AT::t_int_1d_randomread, typename AT::t_neighbors_2d,typename AT::t_neighbors_2d,typename AT::t_int_1d_randomread, int );
  void meam_dens_final(int , int , int , int , double* ,
                       typename ArrayTypes<DeviceType>::t_efloat_1d , int , typename AT::t_int_1d_randomread , typename AT::t_int_1d_randomread , int& );
  void meam_force(int , int , int , int , int , double* ,
                  typename ArrayTypes<DeviceType>::t_efloat_1d , int , typename AT::t_int_1d_randomread , typename AT::t_int_1d_randomread , typename AT::t_x_array_randomread , typename AT::t_int_1d_randomread , 
                  typename AT::t_int_1d_randomread , typename AT::t_f_array , typename ArrayTypes<DeviceType>::t_virial_array ,typename AT::t_int_1d_randomread , typename AT::t_int_1d_randomread, typename AT::t_neighbors_2d, typename AT::t_neighbors_2d, int);
  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void getscreen(int , int, typename AT::t_x_array_randomread , typename AT::t_int_1d_randomread,
                 typename AT::t_int_1d_randomread, int , typename AT::t_int_1d_randomread , typename AT::t_int_1d_randomread ) const;
  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void calc_rho1(int , int , typename AT::t_int_1d_randomread , typename AT::t_int_1d_randomread , typename AT::t_x_array_randomread , typename AT::t_int_1d_randomread, int ) const;
  KOKKOS_INLINE_FUNCTION
  double fcut(const double xi) const;
  KOKKOS_INLINE_FUNCTION
  double dfcut(const double xi, double& dfc) const;
  KOKKOS_INLINE_FUNCTION
  double dCfunc(const double, const double, const double) const;
  KOKKOS_INLINE_FUNCTION
  void dCfunc2(const double, const double, const double, double&, double&) const;
  KOKKOS_INLINE_FUNCTION
  double G_gam(const double, const int, int&) const;
  KOKKOS_INLINE_FUNCTION
  double dG_gam(const double, const int, double&) const;
  KOKKOS_INLINE_FUNCTION
  double zbl(const double, const int, const int) const; 
  KOKKOS_INLINE_FUNCTION
  double erose(const double, const double, const double, const double, const double, const double, const int) const;
  KOKKOS_INLINE_FUNCTION
  void get_shpfcn(const lattice_t latt, double (&s)[3]) const;
  KOKKOS_INLINE_FUNCTION
  int get_Zij(const lattice_t ) const;
  KOKKOS_INLINE_FUNCTION
  int get_Zij2(const lattice_t, const double, const double, double&, double&) const; 
public:
  DAT::tdual_ffloat_1d k_rho, k_rho0, k_rho1, k_rho2, k_rho3, k_frhop;
  typename ArrayTypes<DeviceType>::t_ffloat_1d d_rho, d_rho0,d_rho1, d_rho2, d_rho3, d_frhop;
  HAT::t_ffloat_1d h_rho, h_rho0, h_rho1, h_rho2, h_rho3, h_frhop;
  DAT::tdual_ffloat_1d k_gamma, k_dgamma1, k_dgamma2, k_dgamma3, k_arho2b;
  typename ArrayTypes<DeviceType>::t_ffloat_1d d_gamma, d_dgamma1, d_dgamma2, d_dgamma3, d_arho2b;
  HAT::t_ffloat_1d h_gamma, h_dgamma1, h_dgamma2, h_dgamma3, h_arho2b;
  DAT::tdual_ffloat_2d k_arho1, k_arho2, k_arho3, k_arho3b, k_t_ave, k_tsq_ave;
  typename ArrayTypes<DeviceType>::t_ffloat_2d d_arho1, d_arho2, d_arho3, d_arho3b, d_t_ave, d_tsq_ave;
  HAT::t_ffloat_2d h_arho1, h_arho2, h_arho3, h_arho3b, h_t_ave, h_tsq_ave;
  DAT::tdual_ffloat_2d k_phir, k_phirar, k_phirar1, k_phirar2, k_phirar3, k_phirar4, k_phirar5, k_phirar6;
  typename ArrayTypes<DeviceType>::t_ffloat_2d d_phir, d_phirar, d_phirar1, d_phirar2, d_phirar3, d_phirar4, d_phirar5, d_phirar6;
  HAT::t_ffloat_2d h_phir, h_phirar, h_phirar1, h_phirar2, h_phirar3, h_phirar4, h_phirar5, h_phirar6;
  DAT::tdual_ffloat_1d k_scrfcn, k_dscrfcn, k_fcpair;
  typename ArrayTypes<DeviceType>::t_ffloat_1d d_scrfcn, d_dscrfcn, d_fcpair;
  HAT::t_ffloat_1d h_scrfcn, h_dscrfcn, h_fcpair;

  
};
KOKKOS_INLINE_FUNCTION
static bool iszero_kk(const double f) {
  return fabs(f) < 1e-20;
}


KOKKOS_INLINE_FUNCTION
static double fdiv_zero_kk(const double n, const double d) {
  if (iszero_kk(d))
    return 0.0;
  return n / d;
}

// Functions we need for compat

}
#include "meam_impl_kokkos.h"
//#include "meam_setup_done_kokkos.h"
//#include "meam_funcs_kokkos.h"
//#include "meam_dens_init_kokkos.h"
//#include "meam_dens_final_kokkos.h"
//#include "meam_force_kokkos.h"
#endif
