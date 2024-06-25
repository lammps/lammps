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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(pod/kk,PairPODKokkos<LMPDeviceType>);
PairStyle(pod/kk/device,PairPODKokkos<LMPDeviceType>);
PairStyle(pod/kk/host,PairPODKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_POD_KOKKOS_H
#define LMP_PAIR_POD_KOKKOS_H

#include "pair_pod.h"
#include "kokkos_type.h"
#include "pair_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairPODKokkos : public PairPOD {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairPODKokkos(class LAMMPS *);
  ~PairPODKokkos() override;

  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  int inum, maxneigh;
  int host_flag;

  int eflag, vflag;
  int neighflag;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;

  typedef Kokkos::DualView<F_FLOAT**, DeviceType> tdual_fparams;
  tdual_fparams k_cutsq, k_scale;
  typedef Kokkos::View<F_FLOAT**, DeviceType> t_fparams;
  t_fparams d_cutsq, d_scale;
  typename AT::t_int_1d d_map;

  friend void pair_virial_fdotr_compute<PairPODKokkos>(PairPODKokkos*);

  void grow(int, int);
  void copy_from_pod_class(EAPOD *podptr);
  void divideInterval(int *intervals, int N, int M);
  int calculateNumberOfIntervals(int N, int intervalSize);
  void grow_atoms(int Ni);
  void grow_pairs(int Nij);

  void allocate() override;
  double memory_usage() override;

  typedef Kokkos::View<int*, DeviceType> t_pod_1i;
  typedef Kokkos::View<double*, DeviceType> t_pod_1d;
//   typedef Kokkos::View<int**, DeviceType> t_pod_2i;
//   typedef Kokkos::View<double**, DeviceType> t_pod_2d;
//   typedef Kokkos::View<double**[3], DeviceType> t_pod_3d3;

  int atomBlockSize;        // size of each atom block
  int nAtomBlocks;          // number of atoms blocks
  int atomBlocks[101];      // atom blocks
  double comptime[100];
  int timing;

  int ni;            // number of atoms i in the current atom block
  int nij;           // number of pairs (i,j) in the current atom block
  int nimax;         // maximum number of atoms i
  int nijmax;        // maximum number of pairs (i,j)

  int nelements; // number of elements
  int onebody;   // one-body descriptors
  int besseldegree; // degree of Bessel functions
  int inversedegree; // degree of inverse functions
  int nbesselpars;  // number of Bessel parameters
  int nCoeffPerElement; // number of coefficients per element = (nl1 + Mdesc*nClusters)
  int ns;      // number of snapshots for radial basis functions
  int nl1, nl2, nl3, nl4, nl23, nl33, nl34, nl44, nl;   // number of local descriptors
  int nrbf2, nrbf3, nrbf4, nrbfmax;            // number of radial basis functions
  int nabf3, nabf4;                            // number of angular basis functions
  int K3, K4, Q4;                                  // number of monomials

  // environmental variables
  int nClusters; // number of environment clusters
  int nComponents; // number of principal components
  int Mdesc; // number of base descriptors

  double rin;  // inner cut-off radius
  double rcut; // outer cut-off radius
  double rmax; // rcut - rin
  double rcutsq;

  t_pod_1d rij;         // (xj - xi) for all pairs (I, J)
  t_pod_1d fij;         // force for all pairs (I, J)
  t_pod_1d ei;          // energy for each atom I
  t_pod_1i typeai;         // types of atoms I only
  t_pod_1i numij;          // number of pairs (I, J) for each atom I
  t_pod_1i idxi;           // storing linear indices of atom I for all pairs (I, J)
  t_pod_1i ai;             // IDs of atoms I for all pairs (I, J)
  t_pod_1i aj;             // IDs of atoms J for all pairs (I, J)
  t_pod_1i ti;             // types of atoms I for all pairs (I, J)
  t_pod_1i tj;             // types of atoms J for all pairs (I, J)

  t_pod_1d besselparams;
  t_pod_1d Phi;  // eigenvectors matrix ns x ns
  t_pod_1d rbf;  // radial basis functions nij x nrbfmax
  t_pod_1d rbfx; // x-derivatives of radial basis functions nij x nrbfmax
  t_pod_1d rbfy; // y-derivatives of radial basis functions nij x nrbfmax
  t_pod_1d rbfz; // z-derivatives of radial basis functions nij x nrbfmax
  t_pod_1d abf;  // angular basis functions nij x K3
  t_pod_1d abfx; // x-derivatives of angular basis functions nij x K3
  t_pod_1d abfy; // y-derivatives of angular basis functions nij x K3
  t_pod_1d abfz; // z-derivatives of angular basis functions nij x K3
  t_pod_1d sumU; // sum of radial basis functions ni x K3 x nrbfmax x nelements
  t_pod_1d forcecoeff; // force coefficients ni x K3 x nrbfmax x nelements
  t_pod_1d Proj; // PCA Projection matrix
  t_pod_1d Centroids; // centroids of the clusters
  t_pod_1d bd;   // base descriptors ni x Mdesc
  t_pod_1d cb;   // force coefficients for base descriptors ni x Mdesc
  t_pod_1d pd;   // environment probability descriptors ni x (1 + nComponents + 3*nClusters)
  t_pod_1d coefficients; // coefficients nCoeffPerElement x nelements
  t_pod_1i pq3, pn3, pc3; // arrays to compute 3-body angular basis functions
  t_pod_1i pa4, pb4, pc4; // arrays to compute 4-body angular basis functions
  t_pod_1i ind33l, ind33r; // nl33
  t_pod_1i ind34l, ind34r; // nl34
  t_pod_1i ind44l, ind44r; // nl44
  t_pod_1i elemindex;

  void set_array_to_zero(t_pod_1d a, int N);

  int NeighborCount(t_pod_1i, double, int, int);

  void NeighborList(t_pod_1d l_rij, t_pod_1i l_numij,  t_pod_1i l_typeai, t_pod_1i l_idxi,
    t_pod_1i l_ai, t_pod_1i l_aj, t_pod_1i l_ti, t_pod_1i l_tj, double l_rcutsq, int gi1, int Ni);

  void radialbasis(t_pod_1d rbft, t_pod_1d rbftx, t_pod_1d rbfty, t_pod_1d rbftz,
    t_pod_1d rij, t_pod_1d l_besselparams, double l_rin, double l_rmax, int l_besseldegree,
    int l_inversedegree, int l_nbesselpars, int Nij);

  void matrixMultiply(t_pod_1d a, t_pod_1d b, t_pod_1d c, int r1, int c1, int c2);

  void angularbasis(t_pod_1d l_abf, t_pod_1d l_abfx, t_pod_1d l_abfy, t_pod_1d l_abfz,
        t_pod_1d l_rij, t_pod_1i l_pq3, int l_K3, int N);

  void radialangularsum(t_pod_1d l_sumU, t_pod_1d l_rbf, t_pod_1d l_abf, t_pod_1i l_tj,
    t_pod_1i l_numij, const int l_nelements, const int l_nrbf3, const int l_K3, const int Ni, const int Nij);

  void twobodydesc(t_pod_1d d2, t_pod_1d l_rbf, t_pod_1i l_idxi, t_pod_1i l_tj, int l_nrbf2, const int Ni, const int Nij);

  void threebodydesc(t_pod_1d d3, t_pod_1d l_sumU, t_pod_1i l_pc3, t_pod_1i l_pn3,
        int l_nelements, int l_nrbf3, int l_nabf3, int l_K3, const int Ni);

  void fourbodydesc(t_pod_1d d4,  t_pod_1d l_sumU, t_pod_1i l_pa4, t_pod_1i l_pb4, t_pod_1i l_pc4,
      int l_nelements, int l_nrbf3, int l_nrbf4, int l_nabf4, int l_K3, int l_Q4, int Ni);

  void crossdesc(t_pod_1d d12, t_pod_1d d1, t_pod_1d d2, t_pod_1i ind1, t_pod_1i ind2, int n12, int Ni);

  void crossdesc_reduction(t_pod_1d cb1, t_pod_1d cb2, t_pod_1d c12, t_pod_1d d1,
        t_pod_1d d2, t_pod_1i ind1, t_pod_1i ind2, int n12, int Ni);
  void blockatom_base_descriptors(t_pod_1d bd, int Ni, int Nij);
  void blockatom_base_coefficients(t_pod_1d ei, t_pod_1d cb, t_pod_1d B, int Ni);
  void blockatom_environment_descriptors(t_pod_1d ei, t_pod_1d cb, t_pod_1d B, int Ni);

  void twobody_forces(t_pod_1d fij, t_pod_1d cb2, t_pod_1d l_rbfx, t_pod_1d l_rbfy, t_pod_1d l_rbfz,
          t_pod_1i l_idxi, t_pod_1i l_tj, int l_nrbf2, const int Ni, const int Nij);
  void threebody_forces(t_pod_1d fij, t_pod_1d cb3, t_pod_1d l_rbf, t_pod_1d l_rbfx,
    t_pod_1d l_rbfy, t_pod_1d l_rbfz, t_pod_1d l_abf, t_pod_1d l_abfx, t_pod_1d l_abfy, t_pod_1d l_abfz,
    t_pod_1d l_sumU, t_pod_1i l_idxi, t_pod_1i l_tj, t_pod_1i l_pc3, t_pod_1i l_pn3, t_pod_1i l_elemindex,
    int l_nelements, int l_nrbf3, int l_nabf3, int l_K3, int Ni, int Nij);
  void fourbody_forces(t_pod_1d fij, t_pod_1d cb4, t_pod_1d l_rbf, t_pod_1d l_rbfx, t_pod_1d l_rbfy,
    t_pod_1d l_rbfz, t_pod_1d l_abf, t_pod_1d l_abfx, t_pod_1d l_abfy, t_pod_1d l_abfz, t_pod_1d l_sumU,
    t_pod_1i l_idxi, t_pod_1i l_tj, t_pod_1i l_pa4, t_pod_1i l_pb4, t_pod_1i l_pc4,
    int l_nelements, int l_nrbf3, int l_nrbf4, int l_nabf4, int l_K3, int l_Q4, int Ni, int Nij);

  void threebody_forcecoeff(t_pod_1d fb3, t_pod_1d cb3, t_pod_1d l_sumU, t_pod_1i l_pc3,
    t_pod_1i l_pn3, t_pod_1i l_elemindex, int l_nelements, int l_nrbf3, int l_nabf3, int l_K3, int Ni);

  void fourbody_forcecoeff(t_pod_1d fb4, t_pod_1d cb4, t_pod_1d l_sumU, t_pod_1i l_pa4,
    t_pod_1i l_pb4, t_pod_1i l_pc4, int l_nelements, int l_nrbf3, int l_nrbf4, int l_nabf4, int l_K3, int l_Q4, int Ni);

  void allbody_forces(t_pod_1d fij, t_pod_1d l_forcecoeff, t_pod_1d l_rbf, t_pod_1d l_rbfx,
    t_pod_1d l_rbfy, t_pod_1d l_rbfz, t_pod_1d l_abf, t_pod_1d l_abfx, t_pod_1d l_abfy, t_pod_1d l_abfz,
    t_pod_1i l_idxi, t_pod_1i l_tj, int l_nelements, int l_nrbf3, int l_K3, int Nij);

  void blockatom_energyforce(t_pod_1d l_ei, t_pod_1d l_fij, int Ni, int Nij);
  void tallyenergy(t_pod_1d l_ei, int istart, int Ni);
  void tallyforce(t_pod_1d l_fij, t_pod_1i l_ai, t_pod_1i l_aj, int Nij);
  void tallystress(t_pod_1d l_fij, t_pod_1d l_rij, t_pod_1i l_ai, t_pod_1i l_aj, int Nij);

  void savematrix2binfile(std::string filename, t_pod_1d d_A, int nrows, int ncols);
  void saveintmatrix2binfile(std::string filename, t_pod_1i d_A, int nrows, int ncols);
  void savedatafordebugging();
};
}    // namespace LAMMPS_NS

#endif
#endif
