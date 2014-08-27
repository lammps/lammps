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

#ifndef LMP_ATOM_VEC_KOKKOS_H
#define LMP_ATOM_VEC_KOKKOS_H

#include "atom_vec.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class AtomVecKokkos : public AtomVec {
 public:
  AtomVecKokkos(class LAMMPS *);
  virtual ~AtomVecKokkos() {}

  virtual void sync(ExecutionSpace space, unsigned int mask) = 0;
  virtual void modified(ExecutionSpace space, unsigned int mask) = 0;

  virtual int 
    pack_comm_self(const int &n, const DAT::tdual_int_2d &list, 
                   const int & iswap, const int nfirst, 
                   const int &pbc_flag, const int pbc[]) = 0;
  //{return 0;}
  virtual int 
    pack_comm_kokkos(const int &n, const DAT::tdual_int_2d &list, 
                     const int & iswap, const DAT::tdual_xfloat_2d &buf,
                     const int &pbc_flag, const int pbc[]) = 0;
  //{return 0;}
  virtual void 
    unpack_comm_kokkos(const int &n, const int &nfirst, 
                       const DAT::tdual_xfloat_2d &buf) = 0;
  virtual int 
    pack_border_kokkos(int n, DAT::tdual_int_2d k_sendlist, 
                       DAT::tdual_xfloat_2d buf,int iswap,
                       int pbc_flag, int *pbc, ExecutionSpace space) = 0;
  //{return 0;};
  virtual void 
    unpack_border_kokkos(const int &n, const int &nfirst, 
                         const DAT::tdual_xfloat_2d &buf, 
                         ExecutionSpace space) = 0;

  virtual int 
    pack_exchange_kokkos(const int &nsend, DAT::tdual_xfloat_2d &buf, 
                         DAT::tdual_int_1d k_sendlist,
                         DAT::tdual_int_1d k_copylist,
                         ExecutionSpace space, int dim, X_FLOAT lo, X_FLOAT hi) = 0;
  //{return 0;};
  virtual int 
    unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf, int nrecv,
                           int nlocal, int dim, X_FLOAT lo, X_FLOAT hi,
                           ExecutionSpace space) = 0;
  //{return 0;};

 protected:
  class AtomKokkos *atomKK;
  class CommKokkos *commKK;
};

}

#endif

/* ERROR/WARNING messages:

*/
