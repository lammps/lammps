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

#include "neigh_list_kokkos.h"
#include "atom.h"
#include "memory.h"

using namespace LAMMPS_NS;

enum{NSQ,BIN,MULTI};

/* ---------------------------------------------------------------------- */

template<class Device>
void NeighListKokkos<Device>::clean_copy()
{
  ilist = NULL;
  numneigh = NULL;
  firstneigh = NULL;
  firstdouble = NULL;
  dnum = 0;
  iskip = NULL;
  ijskip = NULL;
  
  ipage = NULL;
  dpage = NULL;
  maxstencil = 0;
  ghostflag = 0;
  maxstencil_multi = 0;
}

/* ---------------------------------------------------------------------- */

template<class Device>
void NeighListKokkos<Device>::grow(int nmax)
{
  // skip if this list is already long enough to store nmax atoms

  if (nmax <= maxatoms) return;
  maxatoms = nmax;

  d_ilist = 
    typename ArrayTypes<Device>::t_int_1d("neighlist:ilist",maxatoms);
  d_numneigh = 
    typename ArrayTypes<Device>::t_int_1d("neighlist:numneigh",maxatoms);
  d_neighbors = 
    typename ArrayTypes<Device>::t_neighbors_2d("neighlist:neighbors",
                                                maxatoms,maxneighs);

  memory->sfree(firstneigh);
  memory->sfree(firstdouble);

  firstneigh = (int **) memory->smalloc(maxatoms*sizeof(int *),
                                        "neighlist:firstneigh");
  if (dnum)
    firstdouble = (double **) memory->smalloc(maxatoms*sizeof(double *),
                                              "neighlist:firstdouble");
}

/* ---------------------------------------------------------------------- */

template<class Device>
void NeighListKokkos<Device>::stencil_allocate(int smax, int style)
{
  int i;

  if (style == BIN) {
    if (smax > maxstencil) {
      maxstencil = smax;
      d_stencil = 
        memory->create_kokkos(d_stencil,h_stencil,stencil,maxstencil,
                              "neighlist:stencil");
      if (ghostflag) {
        memory->destroy(stencilxyz);
        memory->create(stencilxyz,maxstencil,3,"neighlist:stencilxyz");
      }
    }

  } else {
    int n = atom->ntypes;
    if (maxstencil_multi == 0) {
      nstencil_multi = new int[n+1];
      stencil_multi = new int*[n+1];
      distsq_multi = new double*[n+1];
      for (i = 1; i <= n; i++) {
        nstencil_multi[i] = 0;
        stencil_multi[i] = NULL;
        distsq_multi[i] = NULL;
      }
    }
    if (smax > maxstencil_multi) {
      maxstencil_multi = smax;
      for (i = 1; i <= n; i++) {
        memory->destroy(stencil_multi[i]);
        memory->destroy(distsq_multi[i]);
        memory->create(stencil_multi[i],maxstencil_multi,
                       "neighlist:stencil_multi");
        memory->create(distsq_multi[i],maxstencil_multi,
                       "neighlist:distsq_multi");
      }
    }
  }
}

template class NeighListKokkos<LMPDeviceType>;
#if DEVICE==2
template class NeighListKokkos<LMPHostType>;
#endif
