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

#include "atom_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType, int HALF_NEIGH>
void NeighborKokkos::full_bin_kokkos(NeighListKokkos<DeviceType> *list)
{
  const int nall = includegroup?atom->nfirst:atom->nlocal;
  list->grow(nall);

  NeighborKokkosExecute<DeviceType> 
    data(*list,
         k_cutneighsq.view<DeviceType>(),
         k_bincount.view<DeviceType>(),
         k_bins.view<DeviceType>(),nall,
         atomKK->k_x.view<DeviceType>(),
         atomKK->k_type.view<DeviceType>(),
         atomKK->k_mask.view<DeviceType>(),
         atomKK->k_molecule.view<DeviceType>(),
         nbinx,nbiny,nbinz,mbinx,mbiny,mbinz,mbinxlo,mbinylo,mbinzlo,
         bininvx,bininvy,bininvz,
         bboxhi,bboxlo);

  k_cutneighsq.sync<DeviceType>();
  atomKK->sync(Device,X_MASK|TYPE_MASK|MASK_MASK);
  Kokkos::deep_copy(list->d_stencil,list->h_stencil);

  while(data.h_resize() > 0) {
    data.h_resize() = 0;
    deep_copy(data.resize, data.h_resize);

    MemsetZeroFunctor<DeviceType> f_zero;
    f_zero.ptr = (void*) k_bincount.view<DeviceType>().ptr_on_device();
    Kokkos::parallel_for(mbins, f_zero);
    DeviceType::fence();

    NeighborKokkosBinAtomsFunctor<DeviceType> f(data);

    Kokkos::parallel_for(atom->nlocal+atom->nghost, f);
    DeviceType::fence();

    deep_copy(data.h_resize, data.resize);
    if(data.h_resize()) {

      atoms_per_bin += 16;
      k_bins = DAT::tdual_int_2d("bins", mbins, atoms_per_bin);
      data.bins = k_bins.view<DeviceType>();
      data.c_bins = data.bins;
    }
  }

  if(list->d_neighbors.dimension_0()<nall) {
    list->d_neighbors = typename ArrayTypes<DeviceType>::t_neighbors_2d("neighbors", nall*1.1, list->maxneighs);
    list->d_numneigh = typename ArrayTypes<DeviceType>::t_int_1d("numneigh", nall*1.1);
    data.neigh_list.d_neighbors = list->d_neighbors;
    data.neigh_list.d_numneigh = list->d_numneigh;
  }
  data.h_resize()=1;
  while(data.h_resize()) {
    data.h_new_maxneighs() = list->maxneighs;
  data.h_resize() = 0;

  Kokkos::deep_copy(data.resize, data.h_resize);
  Kokkos::deep_copy(data.new_maxneighs, data.h_new_maxneighs);
#if DEVICE==2
    #define BINS_PER_BLOCK 2
    const int factor = atoms_per_bin<64?2:1;
    Kokkos::ParallelWorkRequest config((mbins+factor-1)/factor,atoms_per_bin*factor);
#else
    const int factor = 1;
#endif

if(newton_pair) {
  NeighborKokkosBuildFunctor<DeviceType,HALF_NEIGH,1> f(data,atoms_per_bin * 5 * sizeof(X_FLOAT) * factor);
#if DEVICE==2
  Kokkos::parallel_for(config, f);
#else
  Kokkos::parallel_for(nall, f);
#endif
} else {
  NeighborKokkosBuildFunctor<DeviceType,HALF_NEIGH,0> f(data,atoms_per_bin * 5 * sizeof(X_FLOAT) * factor);
#if DEVICE==2
  Kokkos::parallel_for(config, f);
#else
  Kokkos::parallel_for(nall, f);
#endif
}
  DeviceType::fence();
    deep_copy(data.h_resize, data.resize);

    if(data.h_resize()) {
      deep_copy(data.h_new_maxneighs, data.new_maxneighs);
      list->maxneighs = data.h_new_maxneighs() * 1.2;
      list->d_neighbors = typename ArrayTypes<DeviceType>::t_neighbors_2d("neighbors", list->d_neighbors.dimension_0(), list->maxneighs);
      data.neigh_list.d_neighbors = list->d_neighbors;
      data.neigh_list.maxneighs = list->maxneighs;
    }
  }

  list->inum = nall;
  list->gnum = 0;

}

/* ---------------------------------------------------------------------- */

template<class Device>
KOKKOS_INLINE_FUNCTION
void NeighborKokkosExecute<Device>::binatomsItem(const int &i) const
{
  const int ibin = coord2bin(x(i, 0), x(i, 1), x(i, 2));

  const int ac = Kokkos::atomic_fetch_add(&bincount[ibin], (int)1);
  if(ac < bins.dimension_1()) {
    bins(ibin, ac) = i;
  } else {
    resize() = 1;
  }
}

/* ---------------------------------------------------------------------- */

template<class Device> template<int HalfNeigh,int GhostNewton>
void NeighborKokkosExecute<Device>::
   build_Item(const int &i) const
{
  /* if necessary, goto next page and add pages */
  int n = 0;

  // get subview of neighbors of i

  const AtomNeighbors neighbors_i = neigh_list.get_neighbors(i);
  const X_FLOAT xtmp = x(i, 0);
  const X_FLOAT ytmp = x(i, 1);
  const X_FLOAT ztmp = x(i, 2);
  const int itype = type(i);

  const int ibin = coord2bin(xtmp, ytmp, ztmp);

  const int nstencil = neigh_list.nstencil;
  const typename ArrayTypes<Device>::t_int_1d_const_um stencil
    = neigh_list.d_stencil;

  // loop over all bins in neighborhood (includes ibin)
  if(HalfNeigh)
  for(int m = 0; m < c_bincount(ibin); m++) {
    const int j = c_bins(ibin,m);
  // printf("%i %i %i\n",i,ibin,m,c_bincount(ibin),j);
    const int jtype = type(j);
    //for same bin as atom i skip j if i==j and skip atoms "below and to the left" if using HalfNeighborlists
    if((j == i) || (HalfNeigh && !GhostNewton && (j < i))  ||
        (HalfNeigh && GhostNewton && ((j < i) || ((j >= nlocal) &&
                                       ((x(j, 2) < ztmp) || (x(j, 2) == ztmp && x(j, 1) < ytmp) ||
                                        (x(j, 2) == ztmp && x(j, 1)  == ytmp && x(j, 0) < xtmp)))))
      ) continue;
    //if(Exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;


    const X_FLOAT delx = xtmp - x(j, 0);
    const X_FLOAT dely = ytmp - x(j, 1);
    const X_FLOAT delz = ztmp - x(j, 2);
    const X_FLOAT rsq = delx * delx + dely * dely + delz * delz;
    if(rsq <= cutneighsq(itype,jtype)) {
      if(n<neigh_list.maxneighs) neighbors_i(n) = j;
      n++;
    }
  }

  for(int k = 0; k < nstencil; k++) {
    const int jbin = ibin + stencil[k];
    // get subview of jbin
    if(!GhostNewton&&HalfNeigh&&(ibin==jbin)) continue;
    //const ArrayTypes<Device>::t_int_1d_const_um =Kokkos::subview<t_int_1d_const_um>(bins,jbin,ALL);
      for(int m = 0; m < c_bincount(jbin); m++) {
        const int j = c_bins(jbin,m);
        //if(i==0)
        //printf("%i %i %i %i %i %i %i\n",i,jbin,m,c_bincount(jbin),j,k,stencil[k]);
        const int jtype = type(j);

        if(HalfNeigh && !GhostNewton && (j < i)) continue;
        if(!HalfNeigh && j==i) continue;
        //if(Exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

        const X_FLOAT delx = xtmp - x(j, 0);
        const X_FLOAT dely = ytmp - x(j, 1);
        const X_FLOAT delz = ztmp - x(j, 2);
        const X_FLOAT rsq = delx * delx + dely * dely + delz * delz;
      //if(i==0)
        //printf("%i %i %lf %lf NEIGHS\n",i,j,rsq,cutneighsq(itype,jtype));

        if(rsq <= cutneighsq(itype,jtype)) {
          if(n<neigh_list.maxneighs) neighbors_i(n) = j;
          n++;
        }

      }
  }

  neigh_list.d_numneigh(i) = n;

  if(n >= neigh_list.maxneighs) {
    resize() = 1;

    if(n >= new_maxneighs()) new_maxneighs() = n;
  }
  neigh_list.d_ilist(i) = i;
}

#if DEVICE==2
extern __shared__ X_FLOAT sharedmem[];

/* ---------------------------------------------------------------------- */

template<class DeviceType> template<int HalfNeigh>
__device__ inline
void NeighborKokkosExecute<DeviceType>::build_ItemCuda(DeviceType dev) const
{
  /* loop over atoms in i's bin,
  */
  const int atoms_per_bin = c_bins.dimension_1();
  const int BINS_PER_TEAM = blockDim.x/atoms_per_bin;
  const int MY_BIN = threadIdx.x/atoms_per_bin;
  const int MY_II = threadIdx.x%atoms_per_bin;

  const int ibin = (blockIdx.x)*BINS_PER_TEAM+MY_BIN;

  if(ibin >=c_bincount.dimension_0()) return;
  X_FLOAT* other_x = sharedmem;
  other_x = other_x + 5*atoms_per_bin*MY_BIN;

  int* other_id = (int*) &other_x[4 * atoms_per_bin];

  int bincount_current = c_bincount[ibin];

  const int i = MY_II < bincount_current ? c_bins(ibin, MY_II) : -1;
  /* if necessary, goto next page and add pages */

  int n = 0;

  X_FLOAT xtmp;
  X_FLOAT ytmp;
  X_FLOAT ztmp;
  int itype;
  const AtomNeighbors neighbors_i = neigh_list.get_neighbors((i>=0&&i<nlocal)?i:0);

  if(i >= 0) {
    xtmp = x(i, 0);
    ytmp = x(i, 1);
    ztmp = x(i, 2);
    itype = type(i);
    other_x[MY_II] = xtmp;
    other_x[MY_II + atoms_per_bin] = ytmp;
    other_x[MY_II + 2 * atoms_per_bin] = ztmp;
    other_x[MY_II + 3 * atoms_per_bin] = itype;
  }
  other_id[MY_II] = i;
  int test = (__syncthreads_count(i >= 0 && i <= nlocal) == 0);

  if(test) return;

  if(i >= 0 && i < nlocal) {
    #pragma unroll 4
    for(int m = 0; m < bincount_current; m++) {
      int j = other_id[m];

      //for same bin as atom i skip j if i==j and skip atoms "below and to the left" if using halfneighborlists
      //if(j==i) continue;
      if((j == i) || (HalfNeigh && (j < i)))  continue;

      const X_FLOAT delx = xtmp - other_x[m];
      const X_FLOAT dely = ytmp - other_x[m + atoms_per_bin];
      const X_FLOAT delz = ztmp - other_x[m + 2 * atoms_per_bin];
      const int jtype = other_x[m + 3 * atoms_per_bin];
      const X_FLOAT rsq = delx * delx + dely * dely + delz * delz;
      if((rsq <= cutneighsq(itype,jtype)) && (n < neigh_list.maxneighs)) neighbors_i(n++) = j;
    }
  }
  __syncthreads();

  const int nstencil = neigh_list.nstencil;
  const typename ArrayTypes<DeviceType>::t_int_1d_const_um stencil
    = neigh_list.d_stencil;
  for(int k = 0; k < nstencil; k++) {
    const int jbin = ibin + stencil[k];

    if(ibin == jbin) continue;

    bincount_current = c_bincount[jbin];
    int j = MY_II < bincount_current ? c_bins(jbin, MY_II) : -1;

    if(j >= 0) {
      other_x[MY_II] = x(j, 0);
      other_x[MY_II + atoms_per_bin] = x(j, 1);
      other_x[MY_II + 2 * atoms_per_bin] = x(j, 2);
      other_x[MY_II + 3 * atoms_per_bin] = type(j);
     }

    other_id[MY_II] = j;

    __syncthreads();

    if(i >= 0 && i < nlocal) {
      #pragma unroll 8
      for(int m = 0; m < bincount_current; m++) {
        const int j = other_id[m];

        if(HalfNeigh && (j < i))  continue;

        const X_FLOAT delx = xtmp - other_x[m];
        const X_FLOAT dely = ytmp - other_x[m + atoms_per_bin];
        const X_FLOAT delz = ztmp - other_x[m + 2 * atoms_per_bin];
        const int jtype = other_x[m + 3 * atoms_per_bin];
        const X_FLOAT rsq = delx * delx + dely * dely + delz * delz;
        if((rsq <= cutneighsq(itype,jtype)) && (n < neigh_list.maxneighs)) neighbors_i(n++) = j;
      }
    }
    __syncthreads();
  }

  if(i >= 0 && i < nlocal) {
    neigh_list.d_numneigh(i) = n;
    neigh_list.d_ilist(i) = i;
  }

  if(n >= neigh_list.maxneighs) {
    resize() = 1;

    if(n >= new_maxneighs()) new_maxneighs() = n;
  }
}
#endif

template<class DeviceType>
void NeighborKokkos::full_bin_cluster_kokkos(NeighListKokkos<DeviceType> *list)
{
  const int nall = includegroup?atom->nfirst:atom->nlocal;
  list->grow(nall);

  NeighborKokkosExecute<DeviceType>
    data(*list,
         k_cutneighsq.view<DeviceType>(),
         k_bincount.view<DeviceType>(),
         k_bins.view<DeviceType>(),nall,
         atomKK->k_x.view<DeviceType>(),
         atomKK->k_type.view<DeviceType>(),
         atomKK->k_mask.view<DeviceType>(),
         atomKK->k_molecule.view<DeviceType>(),
         nbinx,nbiny,nbinz,mbinx,mbiny,mbinz,mbinxlo,mbinylo,mbinzlo,
         bininvx,bininvy,bininvz,
         bboxhi,bboxlo);

  k_cutneighsq.sync<DeviceType>();
  atomKK->sync(Device,X_MASK|TYPE_MASK|MASK_MASK);
  Kokkos::deep_copy(list->d_stencil,list->h_stencil);
  DeviceType::fence();

  while(data.h_resize() > 0) {
    data.h_resize() = 0;
    deep_copy(data.resize, data.h_resize);

    MemsetZeroFunctor<DeviceType> f_zero;
    f_zero.ptr = (void*) k_bincount.view<DeviceType>().ptr_on_device();
    Kokkos::parallel_for(mbins, f_zero);
    DeviceType::fence();

    NeighborKokkosBinAtomsFunctor<DeviceType> f(data);

    Kokkos::parallel_for(atom->nlocal+atom->nghost, f);
    DeviceType::fence();

    deep_copy(data.h_resize, data.resize);
    if(data.h_resize()) {

      atoms_per_bin += 16;
      k_bins = DAT::tdual_int_2d("bins", mbins, atoms_per_bin);
      data.bins = k_bins.view<DeviceType>();
      data.c_bins = data.bins;
    }
  }

  if(list->d_neighbors.dimension_0()<nall) {
    list->d_neighbors = typename ArrayTypes<DeviceType>::t_neighbors_2d("neighbors", nall*1.1, list->maxneighs);
    list->d_numneigh = typename ArrayTypes<DeviceType>::t_int_1d("numneigh", nall*1.1);
    data.neigh_list.d_neighbors = list->d_neighbors;
    data.neigh_list.d_numneigh = list->d_numneigh;
  }
  data.h_resize()=1;
  while(data.h_resize()) {
    data.h_new_maxneighs() = list->maxneighs;
  data.h_resize() = 0;

  Kokkos::deep_copy(data.resize, data.h_resize);
  Kokkos::deep_copy(data.new_maxneighs, data.h_new_maxneighs);
#if DEVICE==2
    #define BINS_PER_BLOCK 2
    const int factor = atoms_per_bin<64?2:1;
    Kokkos::ParallelWorkRequest config((mbins+factor-1)/factor,atoms_per_bin*factor);
#else
    const int factor = 1;
#endif

if(newton_pair) {
  NeighborClusterKokkosBuildFunctor<DeviceType,NeighClusterSize> f(data,atoms_per_bin * 5 * sizeof(X_FLOAT) * factor);
//#if DEVICE==2
//  Kokkos::parallel_for(config, f);
//#else
  Kokkos::parallel_for(nall, f);
//#endif
} else {
  NeighborClusterKokkosBuildFunctor<DeviceType,NeighClusterSize> f(data,atoms_per_bin * 5 * sizeof(X_FLOAT) * factor);
//#if DEVICE==2
//  Kokkos::parallel_for(config, f);
//#else
  Kokkos::parallel_for(nall, f);
//#endif
}
  DeviceType::fence();
    deep_copy(data.h_resize, data.resize);

    if(data.h_resize()) {
      deep_copy(data.h_new_maxneighs, data.new_maxneighs);
      list->maxneighs = data.h_new_maxneighs() * 1.2;
      list->d_neighbors = typename ArrayTypes<DeviceType>::t_neighbors_2d("neighbors", list->d_neighbors.dimension_0(), list->maxneighs);
      data.neigh_list.d_neighbors = list->d_neighbors;
      data.neigh_list.maxneighs = list->maxneighs;
    }
  }

  list->inum = nall;
  list->gnum = 0;

}

/* ---------------------------------------------------------------------- */

template<class Device> template<int ClusterSize>
void NeighborKokkosExecute<Device>::
   build_cluster_Item(const int &i) const
{
  /* if necessary, goto next page and add pages */
  int n = 0;

  // get subview of neighbors of i

  const AtomNeighbors neighbors_i = neigh_list.get_neighbors(i);
  const X_FLOAT xtmp = x(i, 0);
  const X_FLOAT ytmp = x(i, 1);
  const X_FLOAT ztmp = x(i, 2);
  const int itype = type(i);

  const int ibin = coord2bin(xtmp, ytmp, ztmp);

  const int nstencil = neigh_list.nstencil;
  const typename ArrayTypes<Device>::t_int_1d_const_um stencil
    = neigh_list.d_stencil;

  for(int k = 0; k < nstencil; k++) {
    const int jbin = ibin + stencil[k];
      for(int m = 0; m < c_bincount(jbin); m++) {
        const int j = c_bins(jbin,m);
        bool skip = i == j;
        for(int k = 0; k< (n<neigh_list.maxneighs?n:neigh_list.maxneighs); k++)
          if((j-(j%ClusterSize)) == neighbors_i(k)) {skip=true;};//{m += ClusterSize - j&(ClusterSize-1)-1; skip=true;}

        if(!skip) {
          const int jtype = type(j);

          const X_FLOAT delx = xtmp - x(j, 0);
          const X_FLOAT dely = ytmp - x(j, 1);
          const X_FLOAT delz = ztmp - x(j, 2);
          const X_FLOAT rsq = delx * delx + dely * dely + delz * delz;

          if(rsq <= cutneighsq(itype,jtype)) {
            if(n<neigh_list.maxneighs) neighbors_i(n) = (j-(j%ClusterSize));
            n++;
            //m += ClusterSize - j&(ClusterSize-1)-1;
          }
        }

      }
  }

  neigh_list.d_numneigh(i) = n;

  if(n >= neigh_list.maxneighs) {
    resize() = 1;

    if(n >= new_maxneighs()) new_maxneighs() = n;
  }
  neigh_list.d_ilist(i) = i;
}
