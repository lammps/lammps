// clang-format off
/* ----------------------------------------------------------------------
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
   Contributing authors:  Stan Moore (SNL)

   Tomas Oppelstrup (LLNL): Optimization which reduces the number
   of iterations in the L,m1,m2 loops (by a factor of up to 10), and
   avoids evaluation of Ylm functions of negative m
------------------------------------------------------------------------- */

#include "compute_orientorder_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "math_const.h"
#include "math_special.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;
using MathSpecial::factorial;

#ifdef DBL_EPSILON
static constexpr double MY_EPSILON = (10.0 * DBL_EPSILON);
#else
static constexpr double MY_EPSILON = (10.0 * 2.220446049250313e-16);
#endif

static constexpr double QEPSILON = 1.0e-6;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeOrientOrderAtomKokkos<DeviceType>::ComputeOrientOrderAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeOrientOrderAtom(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  host_flag = (execution_space == Host);

  d_qnormfac = t_sna_1d("orientorder/atom:qnormfac",nqlist);
  d_qnormfac2 = t_sna_1d("orientorder/atom:qnormfac2",nqlist);

  auto h_qnormfac = Kokkos::create_mirror_view(d_qnormfac);
  auto h_qnormfac2 = Kokkos::create_mirror_view(d_qnormfac2);

  for (int il = 0; il < nqlist; il++) {
    h_qnormfac[il] = qnormfac[il];
    h_qnormfac2[il] = qnormfac2[il];
  }

  Kokkos::deep_copy(d_qnormfac,h_qnormfac);
  Kokkos::deep_copy(d_qnormfac2,h_qnormfac2);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeOrientOrderAtomKokkos<DeviceType>::~ComputeOrientOrderAtomKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_qnarray,qnarray);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeOrientOrderAtomKokkos<DeviceType>::init()
{
  ComputeOrientOrderAtom::init();

  d_qlist = t_sna_1i("orientorder/atom:qlist",nqlist);
  auto h_qlist = Kokkos::create_mirror_view(d_qlist);
  for (int i = 0; i < nqlist; i++)
    h_qlist(i) = qlist[i];
  Kokkos::deep_copy(d_qlist,h_qlist);

  // need an occasional full neighbor list

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct FindMaxNumNeighs {
  typedef DeviceType device_type;
  NeighListKokkos<DeviceType> k_list;

  FindMaxNumNeighs(NeighListKokkos<DeviceType>* nl): k_list(*nl) {}
  ~FindMaxNumNeighs() {k_list.copymode = 1;}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& ii, int& maxneigh) const {
    const int i = k_list.d_ilist[ii];
    const int num_neighs = k_list.d_numneigh[i];
    if (maxneigh < num_neighs) maxneigh = num_neighs;
  }
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeOrientOrderAtomKokkos<DeviceType>::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  // grow order parameter array if necessary

  if (atom->nmax > nmax) {
    memoryKK->destroy_kokkos(k_qnarray,qnarray);
    nmax = atom->nmax;
    memoryKK->create_kokkos(k_qnarray,qnarray,nmax,ncol,"orientorder/atom:qnarray");
    array_atom = qnarray;
    d_qnarray = k_qnarray.template view<DeviceType>();
  }

  chunk_size = MIN(chunksize,inum); // "chunksize" variable is set by user
  chunk_offset = 0;

  if (chunk_size > (int)d_ncount.extent(0)) {
    d_qnm = t_sna_3c("orientorder/atom:qnm",chunk_size,nqlist,qmax+1);
    d_ncount = t_sna_1i("orientorder/atom:ncount",chunk_size);
  }

  copymode = 1;

  // ensure distsq and nearest arrays are long enough

  maxneigh = 0;
  Kokkos::parallel_reduce("ComputeOrientOrderAtomKokkos::find_max_neighs",inum, FindMaxNumNeighs<DeviceType>(k_list), Kokkos::Max<int>(maxneigh));

  if (chunk_size > (int)d_distsq.extent(0) || maxneigh > (int)d_distsq.extent(1)) {
    d_distsq = t_sna_2d_lr("orientorder/atom:distsq",chunk_size,maxneigh);
    d_nearest = t_sna_2i_lr("orientorder/atom:nearest",chunk_size,maxneigh);
    d_rlist = t_sna_3d_lr("orientorder/atom:rlist",chunk_size,maxneigh,3);

    d_distsq_um = d_distsq;
    d_rlist_um = d_rlist;
    d_nearest_um = d_nearest;
  }

  // compute order parameter for each atom in group
  // use full neighbor list to count atoms less than cutoff

  atomKK->sync(execution_space,X_MASK|MASK_MASK);
  x = atomKK->k_x.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int vector_length_default = 1;
  int team_size_default = 1;
  if (!host_flag)
    team_size_default = 32;//max_neighs;

  while (chunk_offset < inum) { // chunk up loop to prevent running out of memory

    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    Kokkos::deep_copy(d_qnm,{0.0,0.0});

    //Neigh
    {
      int vector_length = vector_length_default;
      int team_size = team_size_default;
      check_team_size_for<TagComputeOrientOrderAtomNeigh>(chunk_size,team_size,vector_length);
      typename Kokkos::TeamPolicy<DeviceType, TagComputeOrientOrderAtomNeigh> policy_neigh(chunk_size,team_size,vector_length);
      Kokkos::parallel_for("ComputeOrientOrderAtomNeigh",policy_neigh,*this);
    }

    //Select3
    typename Kokkos::RangePolicy<DeviceType, TagComputeOrientOrderAtomSelect3> policy_select3(0,chunk_size);
    Kokkos::parallel_for("ComputeOrientOrderAtomSelect3",policy_select3,*this);

    //BOOP1
    {
      int vector_length = vector_length_default;
      int team_size = team_size_default;
      check_team_size_for<TagComputeOrientOrderAtomBOOP1>(chunk_size,team_size,vector_length);
      typename Kokkos::TeamPolicy<DeviceType, TagComputeOrientOrderAtomBOOP1> policy_boop1(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      Kokkos::parallel_for("ComputeOrientOrderAtomBOOP1",policy_boop1,*this);
    }

    //BOOP2
    typename Kokkos::RangePolicy<DeviceType, TagComputeOrientOrderAtomBOOP2> policy_boop2(0,chunk_size);
    Kokkos::parallel_for("ComputeOrientOrderAtomBOOP2",policy_boop2,*this);

    chunk_offset += chunk_size;
  } // end while

  copymode = 0;

  k_qnarray.template modify<DeviceType>();
  k_qnarray.template sync<LMPHostType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeOrientOrderAtomKokkos<DeviceType>::operator() (TagComputeOrientOrderAtomNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagComputeOrientOrderAtomNeigh>::member_type& team) const
{
  const int ii = team.league_rank();
  const int i = d_ilist[ii + chunk_offset];
  if (mask[i] & groupbit) {
    const X_FLOAT xtmp = x(i,0);
    const X_FLOAT ytmp = x(i,1);
    const X_FLOAT ztmp = x(i,2);
    const int jnum = d_numneigh[i];

    // loop over list of all neighbors within force cutoff
    // distsq[] = distance sq to each
    // rlist[] = distance vector to each
    // nearest[] = atom indices of neighbors

    int ncount = 0;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,jnum),
        [&] (const int jj, int& count) {
      Kokkos::single(Kokkos::PerThread(team), [&] () {
        int j = d_neighbors(i,jj);
        j &= NEIGHMASK;
        const F_FLOAT delx = x(j,0) - xtmp;
        const F_FLOAT dely = x(j,1) - ytmp;
        const F_FLOAT delz = x(j,2) - ztmp;
        const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq)
         count++;
      });
    },ncount);

    d_ncount(ii) = ncount;

    if (team.team_rank() == 0)
    Kokkos::parallel_scan(Kokkos::ThreadVectorRange(team,jnum),
        [&] (const int jj, int& offset, bool final) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      const F_FLOAT delx = x(j,0) - xtmp;
      const F_FLOAT dely = x(j,1) - ytmp;
      const F_FLOAT delz = x(j,2) - ztmp;
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq) {
        if (final) {
          d_distsq(ii,offset) = rsq;
          d_rlist(ii,offset,0) = delx;
          d_rlist(ii,offset,1) = dely;
          d_rlist(ii,offset,2) = delz;
          d_nearest(ii,offset) = j;
        }
        offset++;
      }
    });
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeOrientOrderAtomKokkos<DeviceType>::operator() (TagComputeOrientOrderAtomSelect3,const int& ii) const {

  const int i = d_ilist[ii + chunk_offset];
  const int ncount = d_ncount(ii);

  // if not nnn neighbors, order parameter = 0

  if ((ncount == 0) || (ncount < nnn)) {
    for (int jj = 0; jj < ncol; jj++)
      d_qnarray(i,jj) = 0.0;
    return;
  }

  // if nnn > 0, use only nearest nnn neighbors

  if (nnn > 0) {
    select3(nnn, ncount, ii);
    d_ncount(ii) = nnn;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeOrientOrderAtomKokkos<DeviceType>::operator() (TagComputeOrientOrderAtomBOOP1,const typename Kokkos::TeamPolicy<DeviceType, TagComputeOrientOrderAtomBOOP1>::member_type& team) const {

  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ncount = d_ncount(ii);
  if (jj >= ncount) return;

  // if not nnn neighbors, order parameter = 0

  if ((ncount == 0) || (ncount < nnn))
    return;

  calc_boop1(ncount, ii, jj);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeOrientOrderAtomKokkos<DeviceType>::operator() (TagComputeOrientOrderAtomBOOP2,const int& ii) const {
  const int ncount = d_ncount(ii);

  // if not nnn neighbors, order parameter = 0

  if ((ncount == 0) || (ncount < nnn))
    return;

  calc_boop2(ncount, ii);
}

/* ----------------------------------------------------------------------
   select3 routine from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   sort auxiliary arrays at same time
------------------------------------------------------------------------- */

// Use no-op do while to create single statement

#define SWAP(view,i,j) do {       \
    tmp = view(i); view(i) = view(j); view(j) = tmp; \
  } while (0)

#define ISWAP(view,i,j) do {        \
    itmp = view(i); view(i) = view(j); view(j) = itmp; \
  } while (0)

#define SWAP3(view,i,j) do {                  \
    tmp = view(i,0); view(i,0) = view(j,0); view(j,0) = tmp; \
    tmp = view(i,1); view(i,1) = view(j,1); view(j,1) = tmp; \
    tmp = view(i,2); view(i,2) = view(j,2); view(j,2) = tmp; \
  } while (0)

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeOrientOrderAtomKokkos<DeviceType>::select3(int k, int n, int ii) const
{
  int i,ir,j,l,mid,ia,itmp;
  double a,tmp,a3[3];

  auto arr = Kokkos::subview(d_distsq_um, ii, Kokkos::ALL);
  auto iarr = Kokkos::subview(d_nearest_um, ii, Kokkos::ALL);
  auto arr3 = Kokkos::subview(d_rlist_um, ii, Kokkos::ALL, Kokkos::ALL);

  l = 0;
  ir = n-1;
  while (true) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
        SWAP(arr,l,ir);
        ISWAP(iarr,l,ir);
        SWAP3(arr3,l,ir);
      }
      return;
    } else {
      mid=((l+ir+2) >> 1) - 1;
      SWAP(arr,mid,l+1);
      ISWAP(iarr,mid,l+1);
      SWAP3(arr3,mid,l+1);
      if (arr[l] > arr[ir]) {
        SWAP(arr,l,ir);
        ISWAP(iarr,l,ir);
        SWAP3(arr3,l,ir);
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr,l+1,ir);
        ISWAP(iarr,l+1,ir);
        SWAP3(arr3,l+1,ir);
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr,l,l+1);
        ISWAP(iarr,l,l+1);
        SWAP3(arr3,l,l+1);
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      ia = iarr[l+1];
      a3[0] = arr3(l+1,0);
      a3[1] = arr3(l+1,1);
      a3[2] = arr3(l+1,2);
      while (true) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        SWAP(arr,i,j);
        ISWAP(iarr,i,j);
        SWAP3(arr3,i,j);
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      iarr[l+1] = iarr[j];
      iarr[j] = ia;
      arr3(l+1,0) = arr3(j,0);
      arr3(l+1,1) = arr3(j,1);
      arr3(l+1,2) = arr3(j,2);
      arr3(j,0) = a3[0];
      arr3(j,1) = a3[1];
      arr3(j,2) = a3[2];
      if (j+1 >= k) ir = j-1;
      if (j+1 <= k) l = i;
    }
  }
}

/* ----------------------------------------------------------------------
   calculate the bond orientational order parameters
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeOrientOrderAtomKokkos<DeviceType>::calc_boop1(int /*ncount*/, int ii, int ineigh) const
{
  const double r0 = d_rlist(ii,ineigh,0);
  const double r1 = d_rlist(ii,ineigh,1);
  const double r2 = d_rlist(ii,ineigh,2);
  const double rmag = sqrt(r0*r0 + r1*r1 + r2*r2);
  if (rmag <= MY_EPSILON) {
    return;
  }

  const double costheta = r2 / rmag;
  SNAcomplex expphi = {r0,r1};
  const double rxymag = sqrt(expphi.re*expphi.re+expphi.im*expphi.im);
  if (rxymag <= MY_EPSILON) {
    expphi.re = 1.0;
    expphi.im = 0.0;
  } else {
    const double rxymaginv = 1.0/rxymag;
    expphi.re *= rxymaginv;
    expphi.im *= rxymaginv;
  }

  for (int il = 0; il < nqlist; il++) {
    const int l = d_qlist[il];

    // calculate spherical harmonics
    // Ylm, -l <= m <= l
    // sign convention: sign(Yll(0,0)) = (-1)^l

    //d_qnm(ii,il,l).re += polar_prefactor(l, 0, costheta);
    const double polar_pf = polar_prefactor(l, 0, costheta);
    Kokkos::atomic_add(&(d_qnm(ii,il,0).re), polar_pf);
    SNAcomplex expphim = {expphi.re,expphi.im};
    for (int m = 1; m <= +l; m++) {
      const double prefactor = polar_prefactor(l, m, costheta);
      SNAcomplex ylm = {prefactor * expphim.re, prefactor * expphim.im};
      Kokkos::atomic_add(&(d_qnm(ii,il,m).re), ylm.re);
      Kokkos::atomic_add(&(d_qnm(ii,il,m).im), ylm.im);
      // Skip calculation of qnm for m<0 due to symmetry
      SNAcomplex tmp;
      tmp.re = expphim.re*expphi.re - expphim.im*expphi.im;
      tmp.im = expphim.re*expphi.im + expphim.im*expphi.re;
      expphim.re = tmp.re;
      expphim.im = tmp.im;
    }
  }
}

/* ----------------------------------------------------------------------
   calculate the bond orientational order parameters
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeOrientOrderAtomKokkos<DeviceType>::calc_boop2(int ncount, int ii) const
{
  const int i = d_ilist[ii + chunk_offset];

  // convert sums to averages

  double facn = 1.0 / ncount;
  for (int il = 0; il < nqlist; il++) {
    int l = d_qlist[il];
    for (int m = 0; m < l+1; m++) {
      d_qnm(ii,il,m).re *= facn;
      d_qnm(ii,il,m).im *= facn;
    }
  }

  // calculate Q_l
  // NOTE: optional W_l_hat and components of Q_qlcomp use these stored Q_l values

  int jj = 0;
  for (int il = 0; il < nqlist; il++) {
    int l = d_qlist[il];
    double qm_sum = d_qnm(ii,il,0).re*d_qnm(ii,il,0).re;
    for (int m = 1; m < l+1; m++)
      qm_sum += 2.0*(d_qnm(ii,il,m).re*d_qnm(ii,il,m).re + d_qnm(ii,il,m).im*d_qnm(ii,il,m).im);
    d_qnarray(i,jj++) = d_qnormfac(il) * sqrt(qm_sum);
  }

  // calculate W_l

  int nterms = 0;
  int widx_count = 0;
  if (wlflag || wlhatflag) {
    for (int il = 0; il < nqlist; il++) {
      int l = d_qlist[il];
      double wlsum = 0.0;
      for (int m1 = -l; m1 <= 0; m1++) {
        const int sgn = 1 - 2*(m1&1); // sgn = (-1)^m1
        for (int m2 = 0; m2 <= ((-m1)>>1); m2++) {
          const int m3 = -(m1 + m2);
          // Loop enforces -L<=m1<=0<=m2<=m3<=L, and m1+m2+m3=0

          // For even L, W3j is invariant under permutation of
          // (m1,m2,m3) and (m1,m2,m3)->(-m1,-m2,-m3). The loop
          // structure enforces visiting only one member of each
          // such symmetry (invariance) group.

          // m1 <= 0, and Qlm[-m] = (-1)^m*conjg(Qlm[m])
          SNAcomplex Q1Q2;
          Q1Q2.re = (d_qnm(ii,il,-m1).re*d_qnm(ii,il,m2).re + d_qnm(ii,il,-m1).im*d_qnm(ii,il,m2).im)*sgn;
          Q1Q2.im = (d_qnm(ii,il,-m1).re*d_qnm(ii,il,m2).im - d_qnm(ii,il,-m1).im*d_qnm(ii,il,m2).re)*sgn;
          const double Q1Q2Q3 = Q1Q2.re*d_qnm(ii,il,m3).re - Q1Q2.im*d_qnm(ii,il,m3).im;
          const double c = d_w3jlist[widx_count++];
          wlsum += Q1Q2Q3*c;

        }
      }
      d_qnarray(i,jj++) = wlsum/d_qnormfac2(il);
      nterms++;
    }
  }

  // calculate W_l_hat

  if (wlhatflag) {
    const int jptr = jj-nterms;
    if (!wlflag) jj = jptr;
    for (int il = 0; il < nqlist; il++) {
      if (d_qnarray(i,il) < QEPSILON)
        d_qnarray(i,jj++) = 0.0;
      else {
        const double qnfac = d_qnormfac(il)/d_qnarray(i,il);
        d_qnarray(i,jj++) = d_qnarray(i,jptr+il) * (qnfac*qnfac*qnfac) * d_qnormfac2(il);
      }
    }
  }

  // Calculate components of Q_l/|Q_l|, for l=qlcomp

  if (qlcompflag) {
    const int il = iqlcomp;
    const int l = qlcomp;
    if (d_qnarray(i,il) < QEPSILON)
      for (int m = 0; m < 2*l+1; m++) {
        d_qnarray(i,jj++) = 0.0;
        d_qnarray(i,jj++) = 0.0;
      }
    else {
      const double qnfac = d_qnormfac(il)/d_qnarray(i,il);
      for (int m = -l; m < 0; m++) {
        // Computed only qnm for m>=0.
        // qnm[-m] = (-1)^m * conjg(qnm[m])
        const int sgn = 1 - 2*(m&1); // sgn = (-1)^m
        d_qnarray(i,jj++) =  d_qnm(ii,il,-m).re * qnfac * sgn;
        d_qnarray(i,jj++) = -d_qnm(ii,il,-m).im * qnfac * sgn;
      }
      for (int m = 0; m < l+1; m++) {
        d_qnarray(i,jj++) = d_qnm(ii,il,m).re * qnfac;
        d_qnarray(i,jj++) = d_qnm(ii,il,m).im * qnfac;
      }
    }
  }

}

/* ----------------------------------------------------------------------
   polar prefactor for spherical harmonic Y_l^m, where
   Y_l^m (theta, phi) = prefactor(l, m, cos(theta)) * exp(i*m*phi)
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double ComputeOrientOrderAtomKokkos<DeviceType>::polar_prefactor(int l, int m, double costheta) const
{
  const int mabs = abs(m);

  double prefactor = 1.0;
  for (int i=l-mabs+1; i < l+mabs+1; ++i)
    prefactor *= static_cast<double>(i);

  prefactor = sqrt(static_cast<double>(2*l+1)/(MY_4PI*prefactor))
    * associated_legendre(l,mabs,costheta);

  if ((m < 0) && (m % 2)) prefactor = -prefactor;

  return prefactor;
}

/* ----------------------------------------------------------------------
   associated legendre polynomial
   sign convention: P(l,l) = (2l-1)!!(-sqrt(1-x^2))^l
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double ComputeOrientOrderAtomKokkos<DeviceType>::associated_legendre(int l, int m, double x) const
{
  if (l < m) return 0.0;

  double p(1.0), pm1(0.0), pm2(0.0);

  if (m != 0) {
    const double msqx = -sqrt(1.0-x*x);
    for (int i=1; i < m+1; ++i)
      p *= static_cast<double>(2*i-1) * msqx;
  }

  for (int i=m+1; i < l+1; ++i) {
    pm2 = pm1;
    pm1 = p;
    p = (static_cast<double>(2*i-1)*x*pm1
         - static_cast<double>(i+m-1)*pm2) / static_cast<double>(i-m);
  }

  return p;
}

/* ----------------------------------------------------------------------
  Initialize table of Wigner 3j symbols
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeOrientOrderAtomKokkos<DeviceType>::init_wigner3j()
{
  ComputeOrientOrderAtom::init_wigner3j();

  d_w3jlist = t_sna_1d("computeorientorderatom:w3jlist",widx_max);
  auto h_w3jlist = Kokkos::create_mirror_view(d_w3jlist);

  for (int i = 0; i< widx_max; i++)
    h_w3jlist(i) = w3jlist[i];

  Kokkos::deep_copy(d_w3jlist,h_w3jlist);
}

/* ----------------------------------------------------------------------
   check max team size
------------------------------------------------------------------------- */

template<class DeviceType>
template<class TagStyle>
void ComputeOrientOrderAtomKokkos<DeviceType>::check_team_size_for(int inum, int &team_size, int vector_length) {
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelForTag());

  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

namespace LAMMPS_NS {
template class ComputeOrientOrderAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeOrientOrderAtomKokkos<LMPHostType>;
#endif
}
