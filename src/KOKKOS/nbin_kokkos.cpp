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

#include "nbin_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "memory_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

static constexpr double SMALL = 1.0e-6;
#define CUT2BIN_RATIO 100

/* ---------------------------------------------------------------------- */

template<class DeviceType>
NBinKokkos<DeviceType>::NBinKokkos(LAMMPS *lmp) : NBinStandard(lmp) {
  atoms_per_bin = 16;

  d_resize = typename AT::t_int_scalar("NeighborKokkosFunctor::resize");
  h_resize = Kokkos::create_mirror_view(d_resize);
  h_resize() = 1;

  kokkos = 1;
}

/* ----------------------------------------------------------------------
   setup neighbor binning geometry
   bin numbering in each dimension is global:
     0 = 0.0 to binsize, 1 = binsize to 2*binsize, etc
     nbin-1,nbin,etc = bbox-binsize to bbox, bbox to bbox+binsize, etc
     -1,-2,etc = -binsize to 0.0, -2*binsize to -binsize, etc
   code will work for any binsize
     since next(xyz) and stencil extend as far as necessary
     binsize = 1/2 of cutoff is roughly optimal
   for orthogonal boxes:
     a dim must be filled exactly by integer # of bins
     in periodic, procs on both sides of PBC must see same bin boundary
     in non-periodic, coord2bin() still assumes this by use of nbin xyz
   for triclinic boxes:
     tilted simulation box cannot contain integer # of bins
     stencil & neigh list built differently to account for this
   mbinlo = lowest global bin any of my ghost atoms could fall into
   mbinhi = highest global bin any of my ghost atoms could fall into
   mbin = number of bins I need in a dimension
------------------------------------------------------------------------- */

template<class DeviceType>
void NBinKokkos<DeviceType>::bin_atoms_setup(int nall)
{
  if (mbins > (int)k_bins.d_view.extent(0)) {
    MemoryKokkos::realloc_kokkos(k_bins,"Neighbor::d_bins",mbins,atoms_per_bin);
    bins = k_bins.view<DeviceType>();

    MemoryKokkos::realloc_kokkos(k_bincount,"Neighbor::d_bincount",mbins);
    bincount = k_bincount.view<DeviceType>();
  }
  if (nall > (int)k_atom2bin.d_view.extent(0)) {
    MemoryKokkos::realloc_kokkos(k_atom2bin,"Neighbor::d_atom2bin",nall);
    atom2bin = k_atom2bin.view<DeviceType>();
  }
}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms
------------------------------------------------------------------------- */

template<class DeviceType>
void NBinKokkos<DeviceType>::bin_atoms()
{
  last_bin = update->ntimestep;

  k_bins.template sync<DeviceType>();
  k_bincount.template sync<DeviceType>();
  k_atom2bin.template sync<DeviceType>();

  h_resize() = 1;

  while (h_resize() > 0) {
    h_resize() = 0;
    Kokkos::deep_copy(d_resize, h_resize);

    MemsetZeroFunctor<DeviceType> f_zero;
    f_zero.ptr = (void*) k_bincount.view<DeviceType>().data();
    Kokkos::parallel_for(mbins, f_zero);

    atomKK->sync(ExecutionSpaceFromDevice<DeviceType>::space,X_MASK);
    x = atomKK->k_x.view<DeviceType>();

    bboxlo_[0] = bboxlo[0]; bboxlo_[1] = bboxlo[1]; bboxlo_[2] = bboxlo[2];
    bboxhi_[0] = bboxhi[0]; bboxhi_[1] = bboxhi[1]; bboxhi_[2] = bboxhi[2];

    NPairKokkosBinAtomsFunctor<DeviceType> f(*this);

    Kokkos::parallel_for(atom->nlocal+atom->nghost, f);

    Kokkos::deep_copy(h_resize, d_resize);
    if (h_resize()) {

      atoms_per_bin += 16;
      k_bins = DAT::tdual_int_2d("bins", mbins, atoms_per_bin);
      bins = k_bins.view<DeviceType>();
      c_bins = bins;
    }
  }

  k_bins.template modify<DeviceType>();
  k_bincount.template modify<DeviceType>();
  k_atom2bin.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NBinKokkos<DeviceType>::binatomsItem(const int &i) const
{
  const int ibin = coord2bin(x(i, 0), x(i, 1), x(i, 2));

  atom2bin(i) = ibin;
  const int ac = Kokkos::atomic_fetch_add(&bincount[ibin], (int)1);
  if (ac < (int)bins.extent(1)) {
    bins(ibin, ac) = i;
  } else {
    d_resize() = 1;
  }
}

namespace LAMMPS_NS {
template class NBinKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class NBinKokkos<LMPHostType>;
#endif
}
