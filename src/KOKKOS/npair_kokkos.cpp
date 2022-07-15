// clang-format off
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

#include "npair_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "update.h"
#include "neighbor_kokkos.h"
#include "nbin_kokkos.h"
#include "nstencil.h"
#include "force.h"
#include "kokkos.h"
#include "transpose_helper_kokkos.h"

namespace LAMMPS_NS {

/* ---------------------------------------------------------------------- */

template<class DeviceType, int HALF_NEIGH, int GHOST, int TRI, int SIZE>
NPairKokkos<DeviceType,HALF_NEIGH,GHOST,TRI,SIZE>::NPairKokkos(LAMMPS *lmp) : NPair(lmp) {

  last_stencil_old = -1;

  // use 1D view for scalars to reduce GPU memory operations

  d_scalars = typename AT::t_int_1d("neighbor:scalars",2);
  h_scalars = HAT::t_int_1d("neighbor:scalars_mirror",2);

  d_resize = Kokkos::subview(d_scalars,0);
  d_new_maxneighs = Kokkos::subview(d_scalars,1);

  h_resize = Kokkos::subview(h_scalars,0);
  h_new_maxneighs = Kokkos::subview(h_scalars,1);
}

/* ----------------------------------------------------------------------
   copy needed info from Neighbor class to this build class
   ------------------------------------------------------------------------- */

template<class DeviceType, int HALF_NEIGH, int GHOST, int TRI, int SIZE>
void NPairKokkos<DeviceType,HALF_NEIGH,GHOST,TRI,SIZE>::copy_neighbor_info()
{
  NPair::copy_neighbor_info();

  NeighborKokkos* neighborKK = (NeighborKokkos*) neighbor;

  // general params

  newton_pair = force->newton_pair;
  k_cutneighsq = neighborKK->k_cutneighsq;

  // overwrite per-type Neighbor cutoffs with custom value set by requestor
  // only works for style = BIN (checked by Neighbor class)

  if (cutoff_custom > 0.0) {
    int n = atom->ntypes;
    auto k_mycutneighsq = DAT::tdual_xfloat_2d("neigh:cutneighsq,",n+1,n+1);
    for (int i = 1; i <= n; i++)
      for (int j = 1; j <= n; j++)
        k_mycutneighsq.h_view(i,j) = cutoff_custom * cutoff_custom;
    k_cutneighsq = k_mycutneighsq;
  }

  k_cutneighsq.modify<LMPHostType>();

  // exclusion info

  k_ex1_type = neighborKK->k_ex1_type;
  k_ex2_type = neighborKK->k_ex2_type;
  k_ex_type = neighborKK->k_ex_type;
  k_ex1_group = neighborKK->k_ex1_group;
  k_ex2_group = neighborKK->k_ex2_group;
  k_ex1_bit = neighborKK->k_ex1_bit;
  k_ex2_bit = neighborKK->k_ex2_bit;
  k_ex_mol_group = neighborKK->k_ex_mol_group;
  k_ex_mol_bit = neighborKK->k_ex_mol_bit;
  k_ex_mol_intra = neighborKK->k_ex_mol_intra;
}

/* ----------------------------------------------------------------------
 copy per-atom and per-bin vectors from NBin class to this build class
 ------------------------------------------------------------------------- */

template<class DeviceType, int HALF_NEIGH, int GHOST, int TRI, int SIZE>
void NPairKokkos<DeviceType,HALF_NEIGH,GHOST,TRI,SIZE>::copy_bin_info()
{
  NPair::copy_bin_info();

  NBinKokkos<DeviceType>* nbKK = (NBinKokkos<DeviceType>*) nb;

  atoms_per_bin = nbKK->atoms_per_bin;
  k_bincount = nbKK->k_bincount;
  k_bins = nbKK->k_bins;
  k_atom2bin = nbKK->k_atom2bin;
}

/* ----------------------------------------------------------------------
 copy needed info from NStencil class to this build class
 ------------------------------------------------------------------------- */

template<class DeviceType, int HALF_NEIGH, int GHOST, int TRI, int SIZE>
void NPairKokkos<DeviceType,HALF_NEIGH,GHOST,TRI,SIZE>::copy_stencil_info()
{
  NPair::copy_stencil_info();
  nstencil = ns->nstencil;

  if (ns->last_stencil != last_stencil_old) {
    // copy stencil to device as it may have changed

    last_stencil_old = ns->last_stencil;

    int maxstencil = ns->get_maxstencil();

    if (maxstencil > (int)k_stencil.extent(0))
      k_stencil = DAT::tdual_int_1d("neighlist:stencil",maxstencil);
    for (int k = 0; k < maxstencil; k++)
      k_stencil.h_view(k) = ns->stencil[k];
    k_stencil.modify<LMPHostType>();
    k_stencil.sync<DeviceType>();
    if (GHOST) {
      if (maxstencil > (int)k_stencilxyz.extent(0))
        k_stencilxyz = DAT::tdual_int_1d_3("neighlist:stencilxyz",maxstencil);
      for (int k = 0; k < maxstencil; k++) {
        k_stencilxyz.h_view(k,0) = ns->stencilxyz[k][0];
        k_stencilxyz.h_view(k,1) = ns->stencilxyz[k][1];
        k_stencilxyz.h_view(k,2) = ns->stencilxyz[k][2];
      }
      k_stencilxyz.modify<LMPHostType>();
      k_stencilxyz.sync<DeviceType>();
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, int HALF_NEIGH, int GHOST, int TRI, int SIZE>
void NPairKokkos<DeviceType,HALF_NEIGH,GHOST,TRI,SIZE>::build(NeighList *list_)
{
  NeighListKokkos<DeviceType>* list = (NeighListKokkos<DeviceType>*) list_;
  const int nlocal = includegroup?atom->nfirst:atom->nlocal;
  int nall = nlocal;
  if (GHOST)
    nall += atom->nghost;

  if (nall == 0) return;

  list->grow(nall);

  NeighborKokkosExecute<DeviceType>
    data(*list,
         k_cutneighsq.view<DeviceType>(),
         k_bincount.view<DeviceType>(),
         k_bins.view<DeviceType>(),
         k_atom2bin.view<DeviceType>(),
         mbins,nstencil,
         k_stencil.view<DeviceType>(),
         k_stencilxyz.view<DeviceType>(),
         nlocal,nall,lmp->kokkos->neigh_transpose,
         atomKK->k_x.view<DeviceType>(),
         atomKK->k_radius.view<DeviceType>(),
         atomKK->k_type.view<DeviceType>(),
         atomKK->k_mask.view<DeviceType>(),
         atomKK->k_molecule.view<DeviceType>(),
         atomKK->k_tag.view<DeviceType>(),
         atomKK->k_special.view<DeviceType>(),
         atomKK->k_nspecial.view<DeviceType>(),
         atomKK->molecular,
         nbinx,nbiny,nbinz,mbinx,mbiny,mbinz,mbinxlo,mbinylo,mbinzlo,
         bininvx,bininvy,bininvz,
         exclude, nex_type,
         k_ex1_type.view<DeviceType>(),
         k_ex2_type.view<DeviceType>(),
         k_ex_type.view<DeviceType>(),
         nex_group,
         k_ex1_group.view<DeviceType>(),
         k_ex2_group.view<DeviceType>(),
         k_ex1_bit.view<DeviceType>(),
         k_ex2_bit.view<DeviceType>(),
         nex_mol,
         k_ex_mol_group.view<DeviceType>(),
         k_ex_mol_bit.view<DeviceType>(),
         k_ex_mol_intra.view<DeviceType>(),
         bboxhi,bboxlo,
         domain->xperiodic,domain->yperiodic,domain->zperiodic,
         domain->xprd_half,domain->yprd_half,domain->zprd_half,
         skin,d_resize,h_resize,d_new_maxneighs,h_new_maxneighs);

  k_cutneighsq.sync<DeviceType>();
  k_ex1_type.sync<DeviceType>();
  k_ex2_type.sync<DeviceType>();
  k_ex_type.sync<DeviceType>();
  k_ex1_group.sync<DeviceType>();
  k_ex2_group.sync<DeviceType>();
  k_ex1_bit.sync<DeviceType>();
  k_ex2_bit.sync<DeviceType>();
  k_ex_mol_group.sync<DeviceType>();
  k_ex_mol_bit.sync<DeviceType>();
  k_ex_mol_intra.sync<DeviceType>();
  k_bincount.sync<DeviceType>();
  k_bins.sync<DeviceType>();
  k_atom2bin.sync<DeviceType>();

  if (atom->molecular != Atom::ATOMIC) {
    if (exclude)
      atomKK->sync(Device,X_MASK|RADIUS_MASK|TYPE_MASK|MASK_MASK|MOLECULE_MASK|TAG_MASK|SPECIAL_MASK);
    else
      atomKK->sync(Device,X_MASK|RADIUS_MASK|TYPE_MASK|TAG_MASK|SPECIAL_MASK);
  } else {
    if (exclude)
      atomKK->sync(Device,X_MASK|RADIUS_MASK|TYPE_MASK|MASK_MASK);
    else
      atomKK->sync(Device,X_MASK|RADIUS_MASK|TYPE_MASK);
  }

  data.special_flag[0] = special_flag[0];
  data.special_flag[1] = special_flag[1];
  data.special_flag[2] = special_flag[2];
  data.special_flag[3] = special_flag[3];

  data.h_resize()=1;
  while (data.h_resize()) {
    data.h_new_maxneighs() = list->maxneighs;
    data.h_resize() = 0;

    Kokkos::deep_copy(d_scalars, h_scalars);

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    #define BINS_PER_BLOCK 2
    const int factor = atoms_per_bin <64?2:1;
#else
    const int factor = 1;
#endif

    if (GHOST) {
      NPairKokkosBuildFunctorGhost<DeviceType,HALF_NEIGH> f(data,atoms_per_bin * 5 * sizeof(X_FLOAT) * factor);
#ifdef LMP_KOKKOS_GPU
      if (ExecutionSpaceFromDevice<DeviceType>::space == Device) {
        int team_size = atoms_per_bin*factor;
        int team_size_max = Kokkos::TeamPolicy<DeviceType>(team_size,Kokkos::AUTO).team_size_max(f,Kokkos::ParallelForTag());
        if (team_size <= team_size_max) {
          Kokkos::TeamPolicy<DeviceType> config((mbins+factor-1)/factor,team_size);
          Kokkos::parallel_for(config, f);
        } else { // fall back to flat method
          f.sharedsize = 0;
          Kokkos::parallel_for(nall, f);
        }
      } else
        Kokkos::parallel_for(nall, f);
#else
      Kokkos::parallel_for(nall, f);
#endif
    } else {
      if (newton_pair) {
        if (SIZE) {
          NPairKokkosBuildFunctorSize<DeviceType,TRI?0:HALF_NEIGH,1,TRI> f(data,atoms_per_bin * 6 * sizeof(X_FLOAT) * factor);
#ifdef LMP_KOKKOS_GPU
          if (ExecutionSpaceFromDevice<DeviceType>::space == Device) {
            int team_size = atoms_per_bin*factor;
            int team_size_max = Kokkos::TeamPolicy<DeviceType>(team_size,Kokkos::AUTO).team_size_max(f,Kokkos::ParallelForTag());
            if (team_size <= team_size_max) {
              Kokkos::TeamPolicy<DeviceType> config((mbins+factor-1)/factor,team_size);
              Kokkos::parallel_for(config, f);
            } else { // fall back to flat method
              f.sharedsize = 0;
              Kokkos::parallel_for(nall, f);
            }
          } else
            Kokkos::parallel_for(nall, f);
#else
          Kokkos::parallel_for(nall, f);
#endif
        } else {
          NPairKokkosBuildFunctor<DeviceType,TRI?0:HALF_NEIGH,1,TRI> f(data,atoms_per_bin * 5 * sizeof(X_FLOAT) * factor);
#ifdef LMP_KOKKOS_GPU
          if (ExecutionSpaceFromDevice<DeviceType>::space == Device) {
            int team_size = atoms_per_bin*factor;
            int team_size_max = Kokkos::TeamPolicy<DeviceType>(team_size,Kokkos::AUTO).team_size_max(f,Kokkos::ParallelForTag());
            if (team_size <= team_size_max) {
              Kokkos::TeamPolicy<DeviceType> config((mbins+factor-1)/factor,team_size);
              Kokkos::parallel_for(config, f);
            } else { // fall back to flat method
              f.sharedsize = 0;
              Kokkos::parallel_for(nall, f);
            }
          } else
            Kokkos::parallel_for(nall, f);
#else
          Kokkos::parallel_for(nall, f);
#endif
        }
      } else {
        if (SIZE) {
          NPairKokkosBuildFunctorSize<DeviceType,HALF_NEIGH,0,0> f(data,atoms_per_bin * 5 * sizeof(X_FLOAT) * factor);
#ifdef LMP_KOKKOS_GPU
          if (ExecutionSpaceFromDevice<DeviceType>::space == Device) {
            int team_size = atoms_per_bin*factor;
            int team_size_max = Kokkos::TeamPolicy<DeviceType>(team_size,Kokkos::AUTO).team_size_max(f,Kokkos::ParallelForTag());
            if (team_size <= team_size_max) {
              Kokkos::TeamPolicy<DeviceType> config((mbins+factor-1)/factor,team_size);
              Kokkos::parallel_for(config, f);
            } else { // fall back to flat method
              f.sharedsize = 0;
              Kokkos::parallel_for(nall, f);
            }
          } else
            Kokkos::parallel_for(nall, f);
#else
          Kokkos::parallel_for(nall, f);
#endif
        } else {
          NPairKokkosBuildFunctor<DeviceType,HALF_NEIGH,0,0> f(data,atoms_per_bin * 5 * sizeof(X_FLOAT) * factor);
#ifdef LMP_KOKKOS_GPU
          if (ExecutionSpaceFromDevice<DeviceType>::space == Device) {
            int team_size = atoms_per_bin*factor;
            int team_size_max = Kokkos::TeamPolicy<DeviceType>(team_size,Kokkos::AUTO).team_size_max(f,Kokkos::ParallelForTag());
            if (team_size <= team_size_max) {
              Kokkos::TeamPolicy<DeviceType> config((mbins+factor-1)/factor,team_size);
              Kokkos::parallel_for(config, f);
            } else { // fall back to flat method
              f.sharedsize = 0;
              Kokkos::parallel_for(nall, f);
            }
          } else
            Kokkos::parallel_for(nall, f);
#else
          Kokkos::parallel_for(nall, f);
#endif
        }
      }
    }
    Kokkos::deep_copy(h_scalars, d_scalars);

    if (data.h_resize()) {
      list->maxneighs = data.h_new_maxneighs() * 1.2;
      int maxatoms = list->d_neighbors.extent(0);
      data.neigh_list.d_neighbors = typename AT::t_neighbors_2d();
      list->d_neighbors = typename AT::t_neighbors_2d();
      list->d_neighbors = typename AT::t_neighbors_2d(Kokkos::NoInit("neighlist:neighbors"), maxatoms, list->maxneighs);
      data.neigh_list.d_neighbors = list->d_neighbors;
      data.neigh_list.maxneighs = list->maxneighs;

      if (lmp->kokkos->neigh_transpose) {
        data.neigh_list.d_neighbors_transpose = typename AT::t_neighbors_2d_lr();
        list->d_neighbors_transpose = typename AT::t_neighbors_2d_lr();
        list->d_neighbors_transpose = typename AT::t_neighbors_2d_lr(Kokkos::NoInit("neighlist:neighbors"), maxatoms, list->maxneighs);
        data.neigh_list.d_neighbors_transpose = list->d_neighbors_transpose;
      }
    }
  }

  if (GHOST) {
    list->inum = atom->nlocal;
    list->gnum = nall - atom->nlocal;
  } else {
    list->inum = nall;
    list->gnum = 0;
  }

  list->k_ilist.template modify<DeviceType>();

  if (lmp->kokkos->neigh_transpose)
    TransposeHelperKokkos<DeviceType, typename AT::t_neighbors_2d,
      typename AT::t_neighbors_2d_lr>(list->d_neighbors, list->d_neighbors_transpose);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int NeighborKokkosExecute<DeviceType>::find_special(const int &i, const int &j) const
{
  const int n1 = nspecial(i,0);
  const int n2 = nspecial(i,1);
  const int n3 = nspecial(i,2);

  for (int k = 0; k < n3; k++) {
    if (special(i,k) == tag(j)) {
      if (k < n1) {
        if (special_flag[1] == 0) return -1;
        else if (special_flag[1] == 1) return 0;
        else return 1;
      } else if (k < n2) {
        if (special_flag[2] == 0) return -1;
        else if (special_flag[2] == 1) return 0;
        else return 2;
      } else {
        if (special_flag[3] == 0) return -1;
        else if (special_flag[3] == 1) return 0;
        else return 3;
      }
    }
  }
  return 0;
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int NeighborKokkosExecute<DeviceType>::exclusion(const int &i,const int &j,
                                             const int &itype,const int &jtype) const
{
  int m;

  if (nex_type && ex_type(itype,jtype)) return 1;

  if (nex_group) {
    for (m = 0; m < nex_group; m++) {
      if (mask(i) & ex1_bit(m) && mask(j) & ex2_bit(m)) return 1;
      if (mask(i) & ex2_bit(m) && mask(j) & ex1_bit(m)) return 1;
    }
  }

  if (nex_mol) {
    for (m = 0; m < nex_mol; m++)
      if (ex_mol_intra[m]) { // intra-chain: exclude i-j pair if on same molecule
        if (mask[i] & ex_mol_bit[m] && mask[j] & ex_mol_bit[m] &&
            molecule[i] == molecule[j]) return 1;
      } else                 // exclude i-j pair if on different molecules
        if (mask[i] & ex_mol_bit[m] && mask[j] & ex_mol_bit[m] &&
            molecule[i] != molecule[j]) return 1;
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType> template<int HalfNeigh,int Newton,int Tri>
KOKKOS_FUNCTION
void NeighborKokkosExecute<DeviceType>::
   build_Item(const int &i) const
{
  /* if necessary, goto next page and add pages */
  int n = 0;
  int which = 0;
  int moltemplate;
  if (molecular == Atom::TEMPLATE) moltemplate = 1;
  else moltemplate = 0;
  // get subview of neighbors of i

  const AtomNeighbors neighbors_i = neigh_transpose ?
    neigh_list.get_neighbors_transpose(i) : neigh_list.get_neighbors(i);
  const X_FLOAT xtmp = x(i, 0);
  const X_FLOAT ytmp = x(i, 1);
  const X_FLOAT ztmp = x(i, 2);
  const int itype = type(i);

  const int ibin = c_atom2bin(i);

  const typename ArrayTypes<DeviceType>::t_int_1d_const_um stencil
    = d_stencil;

  // loop over all bins in neighborhood (includes ibin)
  if (HalfNeigh)
  for (int m = 0; m < c_bincount(ibin); m++) {
    const int j = c_bins(ibin,m);
    const int jtype = type(j);

    //for same bin as atom i skip j if i==j and skip atoms "below and to the left" if using HalfNeighborlists
    if ((j == i) || (HalfNeigh && !Newton && (j < i))  ||
        (HalfNeigh && Newton && ((j < i) || ((j >= nlocal) &&
                                       ((x(j, 2) < ztmp) || (x(j, 2) == ztmp && x(j, 1) < ytmp) ||
                                        (x(j, 2) == ztmp && x(j, 1)  == ytmp && x(j, 0) < xtmp)))))
      ) continue;
    if (exclude && exclusion(i,j,itype,jtype)) continue;

    const X_FLOAT delx = xtmp - x(j, 0);
    const X_FLOAT dely = ytmp - x(j, 1);
    const X_FLOAT delz = ztmp - x(j, 2);
    const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    if (rsq <= cutneighsq(itype,jtype)) {
      if (molecular != Atom::ATOMIC) {
        if (!moltemplate)
          which = find_special(i,j);
            /* else if (imol >= 0) */
            /*   which = find_special(onemols[imol]->special[iatom], */
            /*                        onemols[imol]->nspecial[iatom], */
            /*                        tag[j]-tagprev); */
            /* else which = 0; */
        if (which == 0) {
          if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
          else n++;
        } else if (minimum_image_check(delx,dely,delz)) {
          if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
          else n++;
        }
        else if (which > 0) {
          if (n < neigh_list.maxneighs) neighbors_i(n++) = j ^ (which << SBBITS);
          else n++;
        }
      } else {
        if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
        else n++;
      }
    }
  }

  for (int k = 0; k < nstencil; k++) {
    const int jbin = ibin + stencil[k];

    // get subview of jbin
    if (HalfNeigh && (ibin==jbin)) continue;
    //const ArrayTypes<DeviceType>::t_int_1d_const_um =Kokkos::subview<t_int_1d_const_um>(bins,jbin,ALL);
      for (int m = 0; m < c_bincount(jbin); m++) {

        const int j = c_bins(jbin,m);
        const int jtype = type(j);

        if (HalfNeigh && !Newton && (j < i)) continue;
        if (!HalfNeigh && j==i) continue;
        if (Tri) {
          if (x(j,2) < ztmp) continue;
          if (x(j,2) == ztmp) {
            if (x(j,1) < ytmp) continue;
            if (x(j,1) == ytmp) {
              if (x(j,0) < xtmp) continue;
              if (x(j,0) == xtmp && j <= i) continue;
            }
          }
        }
        if (exclude && exclusion(i,j,itype,jtype)) continue;

        const X_FLOAT delx = xtmp - x(j, 0);
        const X_FLOAT dely = ytmp - x(j, 1);
        const X_FLOAT delz = ztmp - x(j, 2);
        const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighsq(itype,jtype)) {
          if (molecular != Atom::ATOMIC) {
            if (!moltemplate)
              which = NeighborKokkosExecute<DeviceType>::find_special(i,j);
            /* else if (imol >= 0) */
            /*   which = find_special(onemols[imol]->special[iatom], */
            /*                        onemols[imol]->nspecial[iatom], */
            /*                        tag[j]-tagprev); */
            /* else which = 0; */
            if (which == 0) {
              if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
              else n++;
            } else if (minimum_image_check(delx,dely,delz)) {
              if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
              else n++;
            }
            else if (which > 0) {
              if (n < neigh_list.maxneighs) neighbors_i(n++) = j ^ (which << SBBITS);
              else n++;
            }
          } else {
            if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
            else n++;
          }
        }

      }
  }

  neigh_list.d_numneigh(i) = n;

  if (n > neigh_list.maxneighs) {
    resize() = 1;

    if (n > new_maxneighs()) new_maxneighs() = n; // avoid atomics, safe because in while loop
  }

  neigh_list.d_ilist(i) = i;
}

/* ---------------------------------------------------------------------- */

#ifdef KOKKOS_ENABLE_HIP
#include <hip/hip_version.h>
#if HIP_VERSION_MAJOR < 3 || (HIP_VERSION_MAJOR == 3 && HIP_VERSION_MINOR < 7)
// ROCm versions < 3.7 are missing __syncthreads_count, so we define a functional
// but (probably) not performant workaround
__device__ __forceinline__ int __syncthreads_count(int predicate) {
  __shared__ int test_block[1];
  if (!(threadIdx.x || threadIdx.y || threadIdx.z))
    test_block[0] = 0;
  __syncthreads();
  atomicAdd(test_block, predicate);
  __threadfence_block();
  return test_block[0];
}
#endif
#endif

#ifdef LMP_KOKKOS_GPU
template<class DeviceType> template<int HalfNeigh,int Newton,int Tri>
LAMMPS_DEVICE_FUNCTION inline
void NeighborKokkosExecute<DeviceType>::build_ItemGPU(typename Kokkos::TeamPolicy<DeviceType>::member_type dev,
                                                      size_t sharedsize) const
{
  auto* sharedmem = static_cast<X_FLOAT *>(dev.team_shmem().get_shmem(sharedsize));
  /* loop over atoms in i's bin,
  */
  const int atoms_per_bin = c_bins.extent(1);
  const int BINS_PER_TEAM = dev.team_size()/atoms_per_bin <1?1:dev.team_size()/atoms_per_bin;
  const int TEAMS_PER_BIN = atoms_per_bin/dev.team_size()<1?1:atoms_per_bin/dev.team_size();
  const int MY_BIN = dev.team_rank()/atoms_per_bin;

  const int ibin = dev.league_rank()*BINS_PER_TEAM+MY_BIN;

  if (ibin >= mbins) return;

  X_FLOAT* other_x = sharedmem + 5*atoms_per_bin*MY_BIN;
  int* other_id = (int*) &other_x[4 * atoms_per_bin];

  int bincount_current = c_bincount[ibin];

  for (int kk = 0; kk < TEAMS_PER_BIN; kk++) {
    const int MY_II = dev.team_rank()%atoms_per_bin+kk*dev.team_size();
    const int i = MY_II < bincount_current ? c_bins(ibin, MY_II) : -1;
    /* if necessary, goto next page and add pages */

    int n = 0;

    X_FLOAT xtmp;
    X_FLOAT ytmp;
    X_FLOAT ztmp;
    int itype;
    const int index = (i >= 0 && i < nlocal) ? i : 0;
    const AtomNeighbors neighbors_i = neigh_transpose ?
    neigh_list.get_neighbors_transpose(index) : neigh_list.get_neighbors(index);

    if (i >= 0) {
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

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    int test = (__syncthreads_count(i >= 0 && i < nlocal) == 0);
    if (test) return;
#elif defined(KOKKOS_ENABLE_SYCL)
    int not_done = (i >= 0 && i < nlocal);
    dev.team_reduce(Kokkos::Max<int>(not_done));
    if(not_done == 0) return;
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
    dev.team_barrier();
#endif

    if (i >= 0 && i < nlocal) {
      #pragma unroll 4
      for (int m = 0; m < bincount_current; m++) {
        int j = other_id[m];
        const int jtype = other_x[m + 3 * atoms_per_bin];

        //for same bin as atom i skip j if i==j and skip atoms "below and to the left" if using halfneighborlists
        if ((j == i) ||
           (HalfNeigh && !Newton && (j < i))  ||
           (HalfNeigh && Newton &&
              ((j < i) ||
              ((j >= nlocal) && ((x(j, 2) < ztmp) || (x(j, 2) == ztmp && x(j, 1) < ytmp) ||
                (x(j, 2) == ztmp && x(j, 1)  == ytmp && x(j, 0) < xtmp)))))
          ) continue;
          if (Tri) {
            if (x(j,2) < ztmp) continue;
            if (x(j,2) == ztmp) {
              if (x(j,1) < ytmp) continue;
              if (x(j,1) == ytmp) {
                if (x(j,0) < xtmp) continue;
                if (x(j,0) == xtmp && j <= i) continue;
              }
            }
          }
        if (exclude && exclusion(i,j,itype,jtype)) continue;
        const X_FLOAT delx = xtmp - other_x[m];
        const X_FLOAT dely = ytmp - other_x[m + atoms_per_bin];
        const X_FLOAT delz = ztmp - other_x[m + 2 * atoms_per_bin];
        const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighsq(itype,jtype)) {
          if (molecular != Atom::ATOMIC) {
            int which = 0;
            if (!moltemplate)
              which = NeighborKokkosExecute<DeviceType>::find_special(i,j);
            /* else if (imol >= 0) */
            /*   which = find_special(onemols[imol]->special[iatom], */
            /*                        onemols[imol]->nspecial[iatom], */
            /*                        tag[j]-tagprev); */
            /* else which = 0; */
            if (which == 0) {
              if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
              else n++;
            } else if (minimum_image_check(delx,dely,delz)) {
              if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
              else n++;
            }
            else if (which > 0) {
              if (n < neigh_list.maxneighs) neighbors_i(n++) = j ^ (which << SBBITS);
              else n++;
            }
          } else {
            if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
            else n++;
          }
        }

      }
    }
    dev.team_barrier();

    const typename ArrayTypes<DeviceType>::t_int_1d_const_um stencil
      = d_stencil;
    for (int k = 0; k < nstencil; k++) {
      const int jbin = ibin + stencil[k];

      if (ibin == jbin) continue;

      bincount_current = c_bincount[jbin];
      int j = MY_II < bincount_current ? c_bins(jbin, MY_II) : -1;

      if (j >= 0) {
        other_x[MY_II] = x(j, 0);
        other_x[MY_II + atoms_per_bin] = x(j, 1);
        other_x[MY_II + 2 * atoms_per_bin] = x(j, 2);
        other_x[MY_II + 3 * atoms_per_bin] = type(j);
      }

      other_id[MY_II] = j;

      dev.team_barrier();

      if (i >= 0 && i < nlocal) {
        #pragma unroll 8
        for (int m = 0; m < bincount_current; m++) {
          const int j = other_id[m];
          const int jtype = other_x[m + 3 * atoms_per_bin];

          //if(HalfNeigh && (j < i))  continue;
          if (HalfNeigh && !Newton && (j < i)) continue;
          if (!HalfNeigh && j==i) continue;
          if (Tri) {
            if (x(j,2) < ztmp) continue;
            if (x(j,2) == ztmp) {
              if (x(j,1) < ytmp) continue;
              if (x(j,1) == ytmp) {
                if (x(j,0) < xtmp) continue;
                if (x(j,0) == xtmp && j <= i) continue;
              }
            }
          }
          if (exclude && exclusion(i,j,itype,jtype)) continue;

          const X_FLOAT delx = xtmp - other_x[m];
          const X_FLOAT dely = ytmp - other_x[m + atoms_per_bin];
          const X_FLOAT delz = ztmp - other_x[m + 2 * atoms_per_bin];
          const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;

          if (rsq <= cutneighsq(itype,jtype)) {
            if (molecular != Atom::ATOMIC) {
              int which = 0;
              if (!moltemplate)
                which = NeighborKokkosExecute<DeviceType>::find_special(i,j);
              /* else if (imol >= 0) */
              /*   which = find_special(onemols[imol]->special[iatom], */
              /*                        onemols[imol]->nspecial[iatom], */
              /*                        tag[j]-tagprev); */
              /* else which = 0; */
              if (which == 0) {
                if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
                else n++;
              } else if (minimum_image_check(delx,dely,delz)) {
                if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
                else n++;
              }
              else if (which > 0) {
                if (n < neigh_list.maxneighs) neighbors_i(n++) = j ^ (which << SBBITS);
                else n++;
              }
            } else {
              if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
              else n++;
            }
          }

        }
      }
      dev.team_barrier();
    }

    if (i >= 0 && i < nlocal) {
      neigh_list.d_numneigh(i) = n;
      neigh_list.d_ilist(i) = i;
    }

    if (n > neigh_list.maxneighs) {
      resize() = 1;

      if (n > new_maxneighs()) new_maxneighs() = n; // avoid atomics, safe because in while loop
    }
  }
}
#endif

/* ---------------------------------------------------------------------- */

template<class DeviceType>  template<int HalfNeigh>
KOKKOS_FUNCTION
void NeighborKokkosExecute<DeviceType>::
   build_ItemGhost(const int &i) const
{
  /* if necessary, goto next page and add pages */
  int n = 0;
  int which = 0;
  int moltemplate;
  if (molecular == Atom::TEMPLATE) moltemplate = 1;
  else moltemplate = 0;
  // get subview of neighbors of i

  const AtomNeighbors neighbors_i = neigh_transpose ?
    neigh_list.get_neighbors_transpose(i) : neigh_list.get_neighbors(i);
  const X_FLOAT xtmp = x(i, 0);
  const X_FLOAT ytmp = x(i, 1);
  const X_FLOAT ztmp = x(i, 2);
  const int itype = type(i);

  const typename ArrayTypes<DeviceType>::t_int_1d_const_um stencil
    = d_stencil;
  const typename ArrayTypes<DeviceType>::t_int_1d_3_const_um stencilxyz
    = d_stencilxyz;

  // loop over all atoms in surrounding bins in stencil including self
  // when i is a ghost atom, must check if stencil bin is out of bounds
  // skip i = j
  // no molecular test when i = ghost atom

  if (i < nlocal) {
    const int ibin = c_atom2bin(i);
    for (int k = 0; k < nstencil; k++) {
      const int jbin = ibin + stencil[k];
      for (int m = 0; m < c_bincount(jbin); m++) {
        const int j = c_bins(jbin,m);

        if (HalfNeigh && j <= i) continue;
        else if (j == i) continue;

        const int jtype = type[j];
        if (exclude && exclusion(i,j,itype,jtype)) continue;

        const X_FLOAT delx = xtmp - x(j,0);
        const X_FLOAT dely = ytmp - x(j,1);
        const X_FLOAT delz = ztmp - x(j,2);
        const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighsq(itype,jtype)) {
          if (molecular != Atom::ATOMIC) {
            if (!moltemplate)
              which = find_special(i,j);
            /* else if (imol >= 0) */
            /*   which = find_special(onemols[imol]->special[iatom], */
            /*                        onemols[imol]->nspecial[iatom], */
            /*                        tag[j]-tagprev); */
            /* else which = 0; */
            if (which == 0) {
              if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
              else n++;
            } else if (minimum_image_check(delx,dely,delz)) {
              if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
              else n++;
            }
            else if (which > 0) {
              if (n < neigh_list.maxneighs) neighbors_i(n++) = j ^ (which << SBBITS);
              else n++;
            }
          } else {
            if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
            else n++;
          }
        }
      }
    }
  } else {
    int binxyz[3];
    const int ibin = coord2bin(xtmp, ytmp, ztmp, binxyz);
    const int xbin = binxyz[0];
    const int ybin = binxyz[1];
    const int zbin = binxyz[2];
    for (int k = 0; k < nstencil; k++) {
      const int xbin2 = xbin + stencilxyz(k,0);
      const int ybin2 = ybin + stencilxyz(k,1);
      const int zbin2 = zbin + stencilxyz(k,2);
      if (xbin2 < 0 || xbin2 >= mbinx ||
          ybin2 < 0 || ybin2 >= mbiny ||
          zbin2 < 0 || zbin2 >= mbinz) continue;
      const int jbin = ibin + stencil[k];
      for (int m = 0; m < c_bincount(jbin); m++) {
        const int j = c_bins(jbin,m);

        if (HalfNeigh && j <= i) continue;
        else if (j == i) continue;

        const int jtype = type[j];
        if (exclude && exclusion(i,j,itype,jtype)) continue;

        const X_FLOAT delx = xtmp - x(j,0);
        const X_FLOAT dely = ytmp - x(j,1);
        const X_FLOAT delz = ztmp - x(j,2);
        const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighsq(itype,jtype)) {
          if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
          else n++;
        }
      }
    }
  }

  neigh_list.d_numneigh(i) = n;

  if (n > neigh_list.maxneighs) {
    resize() = 1;

    if (n > new_maxneighs()) new_maxneighs() = n; // avoid atomics, safe because in while loop
  }
  neigh_list.d_ilist(i) = i;
}

/* ---------------------------------------------------------------------- */

#ifdef LMP_KOKKOS_GPU
template<class DeviceType> template<int HalfNeigh>
LAMMPS_DEVICE_FUNCTION inline
void NeighborKokkosExecute<DeviceType>::build_ItemGhostGPU(typename Kokkos::TeamPolicy<DeviceType>::member_type dev,
                                                      size_t sharedsize) const
{
  auto* sharedmem = static_cast<X_FLOAT *>(dev.team_shmem().get_shmem(sharedsize));
  // loop over atoms in i's bin

  const int atoms_per_bin = c_bins.extent(1);
  const int BINS_PER_TEAM = dev.team_size()/atoms_per_bin <1?1:dev.team_size()/atoms_per_bin;
  const int TEAMS_PER_BIN = atoms_per_bin/dev.team_size()<1?1:atoms_per_bin/dev.team_size();
  const int MY_BIN = dev.team_rank()/atoms_per_bin;

  const int ibin = dev.league_rank()*BINS_PER_TEAM+MY_BIN;

  if (ibin >= mbins) return;

  X_FLOAT* other_x = sharedmem + 5*atoms_per_bin*MY_BIN;
  int* other_id = (int*) &other_x[4 * atoms_per_bin];

  int bincount_current = c_bincount[ibin];

  for (int kk = 0; kk < TEAMS_PER_BIN; kk++) {
    const int MY_II = dev.team_rank()%atoms_per_bin+kk*dev.team_size();
    const int i = MY_II < bincount_current ? c_bins(ibin, MY_II) : -1;

    int n = 0;

    X_FLOAT xtmp;
    X_FLOAT ytmp;
    X_FLOAT ztmp;
    int itype;
    const int index = (i >= 0 && i < nall) ? i : 0;
    const AtomNeighbors neighbors_i = neigh_transpose ?
    neigh_list.get_neighbors_transpose(index) : neigh_list.get_neighbors(index);

    if (i >= 0) {
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
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    int test = (__syncthreads_count(i >= 0 && i < nall) == 0);
    if (test) return;
#elif defined(KOKKOS_ENABLE_SYCL)
    int not_done = (i >= 0 && i < nall);
    dev.team_reduce(Kokkos::Max<int>(not_done));
    if (not_done == 0) return;
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
    dev.team_barrier();
#endif

    int which = 0;
    int moltemplate;
    if (molecular == Atom::TEMPLATE) moltemplate = 1;
    else moltemplate = 0;

    const typename ArrayTypes<DeviceType>::t_int_1d_const_um stencil
      = d_stencil;
    const typename ArrayTypes<DeviceType>::t_int_1d_3_const_um stencilxyz
      = d_stencilxyz;

    // loop over all atoms in surrounding bins in stencil including self
    // when i is a ghost atom, must check if stencil bin is out of bounds
    // skip i = j
    // no molecular test when i = ghost atom

    int ghost = (i >= nlocal && i < nall);
    int binxyz[3];
    if (ghost)
      coord2bin(xtmp, ytmp, ztmp, binxyz);
    const int xbin = binxyz[0];
    const int ybin = binxyz[1];
    const int zbin = binxyz[2];
    for (int k = 0; k < nstencil; k++) {
      int active = 1;
      if (ghost) {
        const int xbin2 = xbin + stencilxyz(k,0);
        const int ybin2 = ybin + stencilxyz(k,1);
        const int zbin2 = zbin + stencilxyz(k,2);
        if (xbin2 < 0 || xbin2 >= mbinx ||
            ybin2 < 0 || ybin2 >= mbiny ||
            zbin2 < 0 || zbin2 >= mbinz) active = 0;
      }

      const int jbin = ibin + stencil[k];
      bincount_current = c_bincount[jbin];
      int j = MY_II < bincount_current ? c_bins(jbin, MY_II) : -1;

      if (j >= 0) {
        other_x[MY_II] = x(j, 0);
        other_x[MY_II + atoms_per_bin] = x(j, 1);
        other_x[MY_II + 2 * atoms_per_bin] = x(j, 2);
        other_x[MY_II + 3 * atoms_per_bin] = type(j);
      }

      other_id[MY_II] = j;

      dev.team_barrier();

      if (active && i >= 0 && i < nall) {
        #pragma unroll 4
        for (int m = 0; m < bincount_current; m++) {
          const int j = other_id[m];

          if (HalfNeigh && j <= i) continue;
          else if (j == i) continue;

          const int jtype = other_x[m + 3 * atoms_per_bin];
          if (exclude && exclusion(i,j,itype,jtype)) continue;

          const X_FLOAT delx = xtmp - other_x[m];
          const X_FLOAT dely = ytmp - other_x[m + atoms_per_bin];
          const X_FLOAT delz = ztmp - other_x[m + 2 * atoms_per_bin];
          const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;

          if (rsq <= cutneighsq(itype,jtype)) {
            if (molecular != Atom::ATOMIC && !ghost) {
              if (!moltemplate)
                which = NeighborKokkosExecute<DeviceType>::find_special(i,j);
              /* else if (imol >= 0) */
              /*   which = find_special(onemols[imol]->special[iatom], */
              /*                        onemols[imol]->nspecial[iatom], */
              /*                        tag[j]-tagprev); */
              /* else which = 0; */
              if (which == 0) {
                if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
                else n++;
              } else if (minimum_image_check(delx,dely,delz)) {
                if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
                else n++;
              }
              else if (which > 0) {
                if (n < neigh_list.maxneighs) neighbors_i(n++) = j ^ (which << SBBITS);
                else n++;
              }
            } else {
              if (n < neigh_list.maxneighs) neighbors_i(n++) = j;
              else n++;
            }
          }
        }
      }
      dev.team_barrier();
    }

    if (i >= 0 && i < nall) {
      neigh_list.d_numneigh(i) = n;
      neigh_list.d_ilist(i) = i;
    }

    if (n > neigh_list.maxneighs) {
      resize() = 1;

      if (n > new_maxneighs()) new_maxneighs() = n; // avoid atomics, safe because in while loop
    }
  }
}
#endif

/* ---------------------------------------------------------------------- */

template<class DeviceType> template<int HalfNeigh,int Newton,int Tri>
KOKKOS_FUNCTION
void NeighborKokkosExecute<DeviceType>::
   build_ItemSize(const int &i) const
{
  /* if necessary, goto next page and add pages */
  int n = 0;

  // get subview of neighbors of i

  const AtomNeighbors neighbors_i = neigh_transpose ?
    neigh_list.get_neighbors_transpose(i) : neigh_list.get_neighbors(i);
  const X_FLOAT xtmp = x(i, 0);
  const X_FLOAT ytmp = x(i, 1);
  const X_FLOAT ztmp = x(i, 2);
  const X_FLOAT radi = radius(i);
  const int itype = type(i);

  const int ibin = c_atom2bin(i);

  const typename ArrayTypes<DeviceType>::t_int_1d_const_um stencil
    = d_stencil;

  const int mask_history = 3 << SBBITS;

  // loop over all bins in neighborhood (includes ibin)
  if (HalfNeigh)
  for (int m = 0; m < c_bincount(ibin); m++) {
    const int j = c_bins(ibin,m);
    const int jtype = type(j);

    //for same bin as atom i skip j if i==j and skip atoms "below and to the left" if using HalfNeighborlists
    if ((j == i) || (HalfNeigh && !Newton && (j < i))  ||
        (HalfNeigh && Newton && ((j < i) || ((j >= nlocal) &&
                                       ((x(j, 2) < ztmp) || (x(j, 2) == ztmp && x(j, 1) < ytmp) ||
                                        (x(j, 2) == ztmp && x(j, 1)  == ytmp && x(j, 0) < xtmp)))))
      ) continue;
    if (exclude && exclusion(i,j,itype,jtype)) continue;

    const X_FLOAT delx = xtmp - x(j, 0);
    const X_FLOAT dely = ytmp - x(j, 1);
    const X_FLOAT delz = ztmp - x(j, 2);
    const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    const X_FLOAT radsum = radi + radius(j);
    const X_FLOAT cutsq = (radsum + skin) * (radsum + skin);

    if (rsq <= cutsq) {
      if (n < neigh_list.maxneighs) {
        if (neigh_list.history && rsq < radsum*radsum) neighbors_i(n++) = j ^ mask_history;
        else neighbors_i(n++) = j;
      }
      else n++;
    }
  }

  for (int k = 0; k < nstencil; k++) {
    const int jbin = ibin + stencil[k];

    // get subview of jbin
    if (HalfNeigh && (ibin==jbin)) continue;
    //const ArrayTypes<DeviceType>::t_int_1d_const_um =Kokkos::subview<t_int_1d_const_um>(bins,jbin,ALL);
    for (int m = 0; m < c_bincount(jbin); m++) {

      const int j = c_bins(jbin,m);
      const int jtype = type(j);

      if (HalfNeigh && !Newton && (j < i)) continue;
      if (!HalfNeigh && j==i) continue;
      if (Tri) {
        if (x(j,2) < ztmp) continue;
        if (x(j,2) == ztmp) {
          if (x(j,1) < ytmp) continue;
          if (x(j,1) == ytmp) {
            if (x(j,0) < xtmp) continue;
            if (x(j,0) == xtmp && j <= i) continue;
          }
        }
      }
      if (exclude && exclusion(i,j,itype,jtype)) continue;

      const X_FLOAT delx = xtmp - x(j, 0);
      const X_FLOAT dely = ytmp - x(j, 1);
      const X_FLOAT delz = ztmp - x(j, 2);
      const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      const X_FLOAT radsum = radi + radius(j);
      const X_FLOAT cutsq = (radsum + skin) * (radsum + skin);

      if (rsq <= cutsq) {
        if (n < neigh_list.maxneighs) {
          if (neigh_list.history && rsq < radsum*radsum) neighbors_i(n++) = j ^ mask_history;
          else neighbors_i(n++) = j;
        }
        else n++;
      }
    }
  }

  neigh_list.d_numneigh(i) = n;

  if (n > neigh_list.maxneighs) {
    resize() = 1;

    if (n > new_maxneighs()) new_maxneighs() = n; // avoid atomics, safe because in while loop
  }

  neigh_list.d_ilist(i) = i;
}

/* ---------------------------------------------------------------------- */

#ifdef LMP_KOKKOS_GPU
template<class DeviceType> template<int HalfNeigh,int Newton,int Tri>
LAMMPS_DEVICE_FUNCTION inline
void NeighborKokkosExecute<DeviceType>::build_ItemSizeGPU(typename Kokkos::TeamPolicy<DeviceType>::member_type dev,
                                                          size_t sharedsize) const
{
  auto* sharedmem = static_cast<X_FLOAT *>(dev.team_shmem().get_shmem(sharedsize));
  /* loop over atoms in i's bin,
   */
  const int atoms_per_bin = c_bins.extent(1);
  const int BINS_PER_TEAM = dev.team_size()/atoms_per_bin <1?1:dev.team_size()/atoms_per_bin;
  const int TEAMS_PER_BIN = atoms_per_bin/dev.team_size()<1?1:atoms_per_bin/dev.team_size();
  const int MY_BIN = dev.team_rank()/atoms_per_bin;

  const int ibin = dev.league_rank()*BINS_PER_TEAM+MY_BIN;

  if (ibin >= mbins) return;

  X_FLOAT* other_x = sharedmem + 6*atoms_per_bin*MY_BIN;
  int* other_id = (int*) &other_x[5 * atoms_per_bin];

  int bincount_current = c_bincount[ibin];

  for (int kk = 0; kk < TEAMS_PER_BIN; kk++) {
    const int MY_II = dev.team_rank()%atoms_per_bin+kk*dev.team_size();
    const int i = MY_II < bincount_current ? c_bins(ibin, MY_II) : -1;
    /* if necessary, goto next page and add pages */

    int n = 0;

    X_FLOAT xtmp;
    X_FLOAT ytmp;
    X_FLOAT ztmp;
    X_FLOAT radi;
    int itype;
    const int index = (i >= 0 && i < nlocal) ? i : 0;
    const AtomNeighbors neighbors_i = neigh_transpose ?
    neigh_list.get_neighbors_transpose(index) : neigh_list.get_neighbors(index);
    const int mask_history = 3 << SBBITS;

    if (i >= 0) {
      xtmp = x(i, 0);
      ytmp = x(i, 1);
      ztmp = x(i, 2);
      radi = radius(i);
      itype = type(i);
      other_x[MY_II] = xtmp;
      other_x[MY_II + atoms_per_bin] = ytmp;
      other_x[MY_II + 2 * atoms_per_bin] = ztmp;
      other_x[MY_II + 3 * atoms_per_bin] = itype;
      other_x[MY_II + 4 * atoms_per_bin] = radi;
    }
    other_id[MY_II] = i;
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    int test = (__syncthreads_count(i >= 0 && i < nlocal) == 0);
    if (test) return;
#elif defined(KOKKOS_ENABLE_SYCL)
    int not_done = (i >= 0 && i < nlocal);
    dev.team_reduce(Kokkos::Max<int>(not_done));
    if (not_done == 0) return;
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
    dev.team_barrier();
#endif

    if (i >= 0 && i < nlocal) {
      #pragma unroll 4
      for (int m = 0; m < bincount_current; m++) {
        int j = other_id[m];
        const int jtype = other_x[m + 3 * atoms_per_bin];

        //for same bin as atom i skip j if i==j and skip atoms "below and to the left" if using halfneighborlists
        if ((j == i) ||
           (HalfNeigh && !Newton && (j < i))  ||
           (HalfNeigh && Newton &&
            ((j < i) ||
             ((j >= nlocal) && ((x(j, 2) < ztmp) || (x(j, 2) == ztmp && x(j, 1) < ytmp) ||
                                (x(j, 2) == ztmp && x(j, 1)  == ytmp && x(j, 0) < xtmp)))))
           ) continue;
        if (Tri) {
          if (x(j,2) < ztmp) continue;
          if (x(j,2) == ztmp) {
            if (x(j,1) < ytmp) continue;
            if (x(j,1) == ytmp) {
              if (x(j,0) < xtmp) continue;
              if (x(j,0) == xtmp && j <= i) continue;
            }
          }
        }
        if (exclude && exclusion(i,j,itype,jtype)) continue;
        const X_FLOAT delx = xtmp - other_x[m];
        const X_FLOAT dely = ytmp - other_x[m + atoms_per_bin];
        const X_FLOAT delz = ztmp - other_x[m + 2 * atoms_per_bin];
        const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;
        const X_FLOAT radsum = radi + other_x[m + 4 * atoms_per_bin];
        const X_FLOAT cutsq = (radsum + skin) * (radsum + skin);

        if (rsq <= cutsq) {
          if (n < neigh_list.maxneighs) {
            if (neigh_list.history && rsq < radsum*radsum) neighbors_i(n++) = j ^ mask_history;
            else neighbors_i(n++) = j;
          }
          else n++;
        }
      }
    }
    dev.team_barrier();

    const typename ArrayTypes<DeviceType>::t_int_1d_const_um stencil
      = d_stencil;
    for (int k = 0; k < nstencil; k++) {
      const int jbin = ibin + stencil[k];

      if (ibin == jbin) continue;

      bincount_current = c_bincount[jbin];
      int j = MY_II < bincount_current ? c_bins(jbin, MY_II) : -1;

      if (j >= 0) {
        other_x[MY_II] = x(j, 0);
        other_x[MY_II + atoms_per_bin] = x(j, 1);
        other_x[MY_II + 2 * atoms_per_bin] = x(j, 2);
        other_x[MY_II + 3 * atoms_per_bin] = type(j);
        other_x[MY_II + 4 * atoms_per_bin] = radius(j);
      }

      other_id[MY_II] = j;

      dev.team_barrier();

      if (i >= 0 && i < nlocal) {
        #pragma unroll 8
        for (int m = 0; m < bincount_current; m++) {
          const int j = other_id[m];
          const int jtype = other_x[m + 3 * atoms_per_bin];

          if (HalfNeigh && (j < i))  continue;
          if (HalfNeigh && !Newton && (j < i)) continue;
          if (!HalfNeigh && j==i) continue;
          if (Tri) {
            if (x(j,2) < ztmp) continue;
            if (x(j,2) == ztmp) {
              if (x(j,1) < ytmp) continue;
              if (x(j,1) == ytmp) {
                if (x(j,0) < xtmp) continue;
                if (x(j,0) == xtmp && j <= i) continue;
              }
            }
          }
          if (exclude && exclusion(i,j,itype,jtype)) continue;

          const X_FLOAT delx = xtmp - other_x[m];
          const X_FLOAT dely = ytmp - other_x[m + atoms_per_bin];
          const X_FLOAT delz = ztmp - other_x[m + 2 * atoms_per_bin];
          const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;
          const X_FLOAT radsum = radi + other_x[m + 4 * atoms_per_bin];
          const X_FLOAT cutsq = (radsum + skin) * (radsum + skin);

          if (rsq <= cutsq) {
            if (n < neigh_list.maxneighs) {
              if (neigh_list.history && rsq < radsum*radsum) neighbors_i(n++) = j ^ mask_history;
              else neighbors_i(n++) = j;
            }
            else n++;
          }
        }
      }
      dev.team_barrier();
    }

    if (i >= 0 && i < nlocal) {
      neigh_list.d_numneigh(i) = n;
      neigh_list.d_ilist(i) = i;
    }

    if (n > neigh_list.maxneighs) {
      resize() = 1;

      if (n > new_maxneighs()) new_maxneighs() = n; // avoid atomics, safe because in while loop
    }
  }
}
#endif

}

namespace LAMMPS_NS {
template class NPairKokkos<LMPDeviceType,0,0,0,0>;
template class NPairKokkos<LMPDeviceType,0,1,0,0>;
template class NPairKokkos<LMPDeviceType,1,0,0,0>;
template class NPairKokkos<LMPDeviceType,1,1,0,0>;
template class NPairKokkos<LMPDeviceType,1,0,1,0>;
template class NPairKokkos<LMPDeviceType,1,0,0,1>;
template class NPairKokkos<LMPDeviceType,1,0,1,1>;
#ifdef LMP_KOKKOS_GPU
template class NPairKokkos<LMPHostType,0,0,0,0>;
template class NPairKokkos<LMPHostType,0,1,0,0>;
template class NPairKokkos<LMPHostType,1,0,0,0>;
template class NPairKokkos<LMPHostType,1,1,0,0>;
template class NPairKokkos<LMPHostType,1,0,1,0>;
template class NPairKokkos<LMPHostType,1,0,0,1>;
template class NPairKokkos<LMPHostType,1,0,1,1>;
#endif
}
