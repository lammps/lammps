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

#ifdef NPAIR_CLASS
// clang-format off
typedef NPairSSAKokkos<LMPHostType> NPairSSAKokkosHost;
NPairStyle(half/bin/newton/ssa/kk/host,
           NPairSSAKokkosHost,
           NP_HALF | NP_BIN | NP_NEWTON | NP_ORTHO | NP_SSA | NP_GHOST | NP_KOKKOS_HOST);

typedef NPairSSAKokkos<LMPDeviceType> NPairSSAKokkosDevice;
NPairStyle(half/bin/newton/ssa/kk/device,
           NPairSSAKokkosDevice,
           NP_HALF | NP_BIN | NP_NEWTON | NP_ORTHO | NP_SSA | NP_GHOST | NP_KOKKOS_DEVICE);
// clang-format on
#else

// clang-format off
#ifndef LMP_NPAIR_SSA_KOKKOS_H
#define LMP_NPAIR_SSA_KOKKOS_H

#include "npair.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class NPairSSAKokkos : public NPair {
 public:
  typedef ArrayTypes<DeviceType> AT;

  // SSA Work plan data structures
  int ssa_phaseCt;
  DAT::tdual_int_1d k_ssa_phaseLen;
  DAT::tdual_int_1d_3 k_ssa_phaseOff;
  DAT::tdual_int_2d k_ssa_itemLoc;
  DAT::tdual_int_2d k_ssa_itemLen;
  typename AT::t_int_1d ssa_phaseLen;
  typename AT::t_int_1d_3 ssa_phaseOff;
  typename AT::t_int_2d ssa_itemLoc;
  typename AT::t_int_2d ssa_itemLen;

  const int ssa_gphaseCt;
  DAT::tdual_int_1d k_ssa_gphaseLen;
  DAT::tdual_int_2d k_ssa_gitemLoc;
  DAT::tdual_int_2d k_ssa_gitemLen;
  typename AT::t_int_1d ssa_gphaseLen;
  typename AT::t_int_2d ssa_gitemLoc;
  typename AT::t_int_2d ssa_gitemLen;

  NPairSSAKokkos(class LAMMPS *);
  void copy_neighbor_info() override;
  void copy_bin_info() override;
  void copy_stencil_info() override;
  void build(class NeighList *) override;
 private:
  // data from Neighbor class

  DAT::tdual_xfloat_2d k_cutneighsq;

  // exclusion data from Neighbor class

  DAT::tdual_int_1d k_ex1_type,k_ex2_type;
  DAT::tdual_int_2d k_ex_type;
  DAT::tdual_int_1d k_ex1_group,k_ex2_group;
  DAT::tdual_int_1d k_ex1_bit,k_ex2_bit;
  DAT::tdual_int_1d k_ex_mol_group;
  DAT::tdual_int_1d k_ex_mol_bit;
  DAT::tdual_int_1d k_ex_mol_intra;

  // data from NBinSSA class

  int atoms_per_bin;
  DAT::tdual_int_1d k_bincount;
  DAT::tdual_int_2d k_bins;
  int ghosts_per_gbin;
  DAT::tdual_int_1d k_gbincount;
  DAT::tdual_int_2d k_gbins;
  int lbinxlo, lbinxhi, lbinylo, lbinyhi, lbinzlo, lbinzhi;

  // data from NStencilSSA class

  int nstencil;
  DAT::tdual_int_1d k_stencil;  // # of J neighs for each I
  DAT::tdual_int_1d_3 k_stencilxyz;
  DAT::tdual_int_1d k_nstencil_ssa;
  int sx1, sy1, sz1;
};

template<class DeviceType>
class NPairSSAKokkosExecute
{
  typedef ArrayTypes<DeviceType> AT;

 public:
  NeighListKokkos<DeviceType> neigh_list;

  // data from Neighbor class

  const typename AT::t_xfloat_2d_randomread cutneighsq;

  // exclusion data from Neighbor class

  const int exclude;

  const int nex_type;
  const typename AT::t_int_1d_const ex1_type,ex2_type;
  const typename AT::t_int_2d_const ex_type;

  const int nex_group;
  const typename AT::t_int_1d_const ex1_group,ex2_group;
  const typename AT::t_int_1d_const ex1_bit,ex2_bit;

  const int nex_mol;
  const typename AT::t_int_1d_const ex_mol_group;
  const typename AT::t_int_1d_const ex_mol_bit;
  const typename AT::t_int_1d_const ex_mol_intra;

  // data from NBinSSA class

  const typename AT::t_int_1d bincount;
  const typename AT::t_int_1d_const c_bincount;
  typename AT::t_int_2d bins;
  typename AT::t_int_2d_const c_bins;
  const typename AT::t_int_1d gbincount;
  const typename AT::t_int_1d_const c_gbincount;
  typename AT::t_int_2d gbins;
  typename AT::t_int_2d_const c_gbins;
  const int lbinxlo, lbinxhi, lbinylo, lbinyhi, lbinzlo, lbinzhi;


  // data from NStencil class

  const int nstencil;
  const int sx1, sy1, sz1;
  typename AT::t_int_1d d_stencil;  // # of J neighs for each I
  typename AT::t_int_1d_3 d_stencilxyz;
  typename AT::t_int_1d d_nstencil_ssa;

  // data from Atom class

  const typename AT::t_x_array_randomread x;
  const typename AT::t_int_1d_const type,mask;
  const typename AT::t_tagint_1d_const molecule;
  const typename AT::t_tagint_1d_const tag;
  const typename AT::t_tagint_2d_const special;
  const typename AT::t_int_2d_const nspecial;
  const int molecular;
  int moltemplate;

  int special_flag[4];

  const int nbinx,nbiny,nbinz;
  const int mbinx,mbiny,mbinz;
  const int mbinxlo,mbinylo,mbinzlo;
  const X_FLOAT bininvx,bininvy,bininvz;
  X_FLOAT bboxhi[3],bboxlo[3];

  const int nlocal;

  typename AT::t_int_scalar resize;
  typename AT::t_int_scalar new_maxneighs;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_resize;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_new_maxneighs;

  const int xperiodic, yperiodic, zperiodic;
  const int xprd_half, yprd_half, zprd_half;

  // SSA Work plan data structures
  int ssa_phaseCt;
  typename AT::t_int_1d d_ssa_phaseLen;
  typename AT::t_int_1d_3_const d_ssa_phaseOff;
  typename AT::t_int_2d d_ssa_itemLoc;
  typename AT::t_int_2d d_ssa_itemLen;
  int ssa_gphaseCt;
  typename AT::t_int_1d d_ssa_gphaseLen;
  typename AT::t_int_2d d_ssa_gitemLoc;
  typename AT::t_int_2d d_ssa_gitemLen;

  NPairSSAKokkosExecute(
        const NeighListKokkos<DeviceType> &_neigh_list,
        const typename AT::t_xfloat_2d_randomread &_cutneighsq,
        const typename AT::t_int_1d &_bincount,
        const typename AT::t_int_2d &_bins,
        const typename AT::t_int_1d &_gbincount,
        const typename AT::t_int_2d &_gbins,
        const int _lbinxlo, const int _lbinxhi,
        const int _lbinylo, const int _lbinyhi,
        const int _lbinzlo, const int _lbinzhi,
        const int _nstencil, const int _sx1, const int _sy1, const int _sz1,
        const typename AT::t_int_1d &_d_stencil,
        const typename AT::t_int_1d_3 &_d_stencilxyz,
        const typename AT::t_int_1d &_d_nstencil_ssa,
        const int _ssa_phaseCt,
        const typename AT::t_int_1d &_d_ssa_phaseLen,
        const typename AT::t_int_1d_3 &_d_ssa_phaseOff,
        const typename AT::t_int_2d &_d_ssa_itemLoc,
        const typename AT::t_int_2d &_d_ssa_itemLen,
        const int _ssa_gphaseCt,
        const typename AT::t_int_1d &_d_ssa_gphaseLen,
        const typename AT::t_int_2d &_d_ssa_gitemLoc,
        const typename AT::t_int_2d &_d_ssa_gitemLen,
        const int _nlocal,
        const typename AT::t_x_array_randomread &_x,
        const typename AT::t_int_1d_const &_type,
        const typename AT::t_int_1d_const &_mask,
        const typename AT::t_tagint_1d_const &_molecule,
        const typename AT::t_tagint_1d_const &_tag,
        const typename AT::t_tagint_2d_const &_special,
        const typename AT::t_int_2d_const &_nspecial,
        const int &_molecular,
        const int & _nbinx,const int & _nbiny,const int & _nbinz,
        const int & _mbinx,const int & _mbiny,const int & _mbinz,
        const int & _mbinxlo,const int & _mbinylo,const int & _mbinzlo,
        const X_FLOAT &_bininvx,const X_FLOAT &_bininvy,const X_FLOAT &_bininvz,
        const int & _exclude,const int & _nex_type,
        const typename AT::t_int_1d_const & _ex1_type,
        const typename AT::t_int_1d_const & _ex2_type,
        const typename AT::t_int_2d_const & _ex_type,
        const int & _nex_group,
        const typename AT::t_int_1d_const & _ex1_group,
        const typename AT::t_int_1d_const & _ex2_group,
        const typename AT::t_int_1d_const & _ex1_bit,
        const typename AT::t_int_1d_const & _ex2_bit,
        const int & _nex_mol,
        const typename AT::t_int_1d_const & _ex_mol_group,
        const typename AT::t_int_1d_const & _ex_mol_bit,
        const typename AT::t_int_1d_const & _ex_mol_intra,
        const X_FLOAT *_bboxhi, const X_FLOAT* _bboxlo,
        const int & _xperiodic, const int & _yperiodic, const int & _zperiodic,
        const int & _xprd_half, const int & _yprd_half, const int & _zprd_half):
    neigh_list(_neigh_list), cutneighsq(_cutneighsq),
    exclude(_exclude),nex_type(_nex_type),
    ex1_type(_ex1_type),ex2_type(_ex2_type),ex_type(_ex_type),
    nex_group(_nex_group),
    ex1_group(_ex1_group),ex2_group(_ex2_group),
    ex1_bit(_ex1_bit),ex2_bit(_ex2_bit),nex_mol(_nex_mol),
    ex_mol_group(_ex_mol_group),ex_mol_bit(_ex_mol_bit),
    ex_mol_intra(_ex_mol_intra),
    bincount(_bincount),c_bincount(_bincount),bins(_bins),c_bins(_bins),
    gbincount(_gbincount),c_gbincount(_gbincount),gbins(_gbins),c_gbins(_gbins),
    lbinxlo(_lbinxlo),lbinxhi(_lbinxhi),
    lbinylo(_lbinylo),lbinyhi(_lbinyhi),
    lbinzlo(_lbinzlo),lbinzhi(_lbinzhi),
    nstencil(_nstencil),sx1(_sx1),sy1(_sy1),sz1(_sz1),
    d_stencil(_d_stencil),d_stencilxyz(_d_stencilxyz),d_nstencil_ssa(_d_nstencil_ssa),
    x(_x),type(_type),mask(_mask),molecule(_molecule),
    tag(_tag),special(_special),nspecial(_nspecial),molecular(_molecular),
    nbinx(_nbinx),nbiny(_nbiny),nbinz(_nbinz),
    mbinx(_mbinx),mbiny(_mbiny),mbinz(_mbinz),
    mbinxlo(_mbinxlo),mbinylo(_mbinylo),mbinzlo(_mbinzlo),
    bininvx(_bininvx),bininvy(_bininvy),bininvz(_bininvz),
    nlocal(_nlocal),
    xperiodic(_xperiodic),yperiodic(_yperiodic),zperiodic(_zperiodic),
    xprd_half(_xprd_half),yprd_half(_yprd_half),zprd_half(_zprd_half),
    ssa_phaseCt(_ssa_phaseCt),
    d_ssa_phaseLen(_d_ssa_phaseLen),
    d_ssa_phaseOff(_d_ssa_phaseOff),
    d_ssa_itemLoc(_d_ssa_itemLoc),
    d_ssa_itemLen(_d_ssa_itemLen),
    ssa_gphaseCt(_ssa_gphaseCt),
    d_ssa_gphaseLen(_d_ssa_gphaseLen),
    d_ssa_gitemLoc(_d_ssa_gitemLoc),
    d_ssa_gitemLen(_d_ssa_gitemLen)
    {

    if (molecular == 2) moltemplate = 1;
    else moltemplate = 0;

    bboxlo[0] = _bboxlo[0]; bboxlo[1] = _bboxlo[1]; bboxlo[2] = _bboxlo[2];
    bboxhi[0] = _bboxhi[0]; bboxhi[1] = _bboxhi[1]; bboxhi[2] = _bboxhi[2];

    resize = typename AT::t_int_scalar("NPairSSAKokkosExecute::resize");
    h_resize = Kokkos::create_mirror_view(resize);
    h_resize() = 1;
    new_maxneighs = typename AT::
      t_int_scalar("NPairSSAKokkosExecute::new_maxneighs");
    h_new_maxneighs = Kokkos::create_mirror_view(new_maxneighs);
    h_new_maxneighs() = neigh_list.maxneighs;
  };

  ~NPairSSAKokkosExecute() {neigh_list.copymode = 1;};

  KOKKOS_FUNCTION
  void build_locals_onePhase(const bool firstTry, int workPhase) const;

  KOKKOS_FUNCTION
  void build_ghosts_onePhase(int workPhase) const;

  KOKKOS_INLINE_FUNCTION
  int coord2bin(const X_FLOAT & x,const X_FLOAT & y,const X_FLOAT & z, int* i) const
  {
    int ix,iy,iz;

    if (x >= bboxhi[0])
      ix = static_cast<int> ((x-bboxhi[0])*bininvx) + nbinx;
    else if (x >= bboxlo[0]) {
      ix = static_cast<int> ((x-bboxlo[0])*bininvx);
      ix = MIN(ix,nbinx-1);
    } else
      ix = static_cast<int> ((x-bboxlo[0])*bininvx) - 1;

    if (y >= bboxhi[1])
      iy = static_cast<int> ((y-bboxhi[1])*bininvy) + nbiny;
    else if (y >= bboxlo[1]) {
      iy = static_cast<int> ((y-bboxlo[1])*bininvy);
      iy = MIN(iy,nbiny-1);
    } else
      iy = static_cast<int> ((y-bboxlo[1])*bininvy) - 1;

    if (z >= bboxhi[2])
      iz = static_cast<int> ((z-bboxhi[2])*bininvz) + nbinz;
    else if (z >= bboxlo[2]) {
      iz = static_cast<int> ((z-bboxlo[2])*bininvz);
      iz = MIN(iz,nbinz-1);
    } else
      iz = static_cast<int> ((z-bboxlo[2])*bininvz) - 1;

    i[0] = ix - mbinxlo;
    i[1] = iy - mbinylo;
    i[2] = iz - mbinzlo;

    return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
  }

  KOKKOS_INLINE_FUNCTION
  int exclusion(const int &i,const int &j, const int &itype,const int &jtype) const;

  KOKKOS_INLINE_FUNCTION
  int find_special(const int &i, const int &j) const;

  KOKKOS_INLINE_FUNCTION
  int minimum_image_check(double dx, double dy, double dz) const {
    if (xperiodic && fabs(dx) > xprd_half) return 1;
    if (yperiodic && fabs(dy) > yprd_half) return 1;
    if (zperiodic && fabs(dz) > zprd_half) return 1;
    return 0;
  }

};

}

#endif
#endif

