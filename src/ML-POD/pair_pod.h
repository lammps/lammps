/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(pod,PairPOD);
// clang-format on
#else

#ifndef LMP_PAIR_POD_H
#define LMP_PAIR_POD_H

#include "pair.h"

namespace LAMMPS_NS {

class PairPOD : public Pair {
 public:
  PairPOD(class LAMMPS *);
  ~PairPOD() override;
  void compute(int, int) override;

  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;

  void lammpsNeighborList(double *rij1, int *ai1, int *aj1, int *ti1, int *tj1, double **x,
                          int **firstneigh, int *atomtype, int *map, int *numneigh, double rcutsq,
                          int i);
  void NeighborCount(double **x, int **firstneigh, int *ilist, int *numneigh, double rcutsq,
                     int i1);
  void NeighborList(double **x, int **firstneigh, int *atomtype, int *map, int *ilist,
                    int *numneigh, double rcutsq, int i1);
  void tallyenergy(double *ei, int istart, int Ni);
  void tallystress(double *fij, double *rij, int *ai, int *aj, int nlocal, int N);
  void tallyforce(double **force, double *fij, int *ai, int *aj, int N);
  void divideInterval(int *intervals, int N, int M);
  int calculateNumberOfIntervals(int N, int intervalSize);
  int numberOfNeighbors();

  void copy_data_from_pod_class();
  void radialbasis(double *rbft, double *rbftx, double *rbfty, double *rbftz, double *rij, int Nij);
  void orthogonalradialbasis(int Nij);
  void angularbasis(double *tm, double *tmu, double *tmv, double *tmw, int N);
  void radialangularsum(int Ni, int Nij);
  void radialangularsum2(int Ni);
  void twobodydesc(double *d2, int Ni, int Nij);
  void twobodydescderiv(double *dd2, int Nij);
  void twobodydescderiv(double *d2, double *dd2, int Ni, int Nij);
  void threebodydesc(double *d3, int Ni);
  void threebodydescderiv(double *dd3, int Nij);
  void extractsumU(int Ni);
  void fourbodydesc(double *d4, int Ni);
  void fourbodydescderiv(double *dd4, int Nij);
  void crossdesc(double *d12, double *d1, double *d2, int *ind1, int *ind2, int n12, int Ni);
  void crossdescderiv(double *dd12, double *d1, double *d2, double *dd1, double *dd2, int *ind1,
                      int *ind2, int *idxi, int n12, int Ni, int Nij);
  void blockatombase_descriptors(double *bd1, double *bdd1, int Ni, int Nij);
  void blockatomenergyforce(double *ei, double *fij, int Ni, int Nij);

  void crossdesc_reduction(double *cb1, double *cb2, double *c12, double *d1, double *d2, int *ind1,
                           int *ind2, int n12, int Ni);
  void blockatom_base_descriptors(double *bd1, int Ni, int Nij);
  void blockatom_base_coefficients(double *ei, double *cb, double *B, int Ni);
  void blockatom_environment_descriptors(double *ei, double *cb, double *B, int Ni);
  void blockatom_energyforce(double *ei, double *fij, int Ni, int Nij);
  void blockatom_energies(double *ei, int Ni, int Nij);
  void blockatom_forces(double *fij, int Ni, int Nij);

  void twobody_forces(double *fij, double *cb2, int Ni, int Nij);
  void threebody_forces(double *fij, double *cb3, int Ni, int Nij);
  void fourbody_forces(double *fij, double *cb4, int Ni, int Nij);

  void threebody_forcecoeff(double *fb3, double *cb3, int Ni);
  void fourbody_forcecoeff(double *fb4, double *cb4, int Ni);
  void allbody_forces(double *fij, double *forcecoeff, int Nij);

  void savematrix2binfile(std::string filename, double *A, int nrows, int ncols);
  void saveintmatrix2binfile(std::string filename, int *A, int nrows, int ncols);
  void savedatafordebugging();

 protected:
  class EAPOD *fastpodptr;
  virtual void allocate();
  void grow_atoms(int Ni);
  void grow_pairs(int Nij);

  int atomBlockSize;      // size of each atom block
  int nAtomBlocks;        // number of atoms blocks
  int atomBlocks[101];    // atom blocks

  int ni;        // total number of atoms i
  int nij;       // total number of pairs (i,j)
  int nimax;     // maximum number of atoms i
  int nijmax;    // maximum number of pairs (i,j)

  int nelements;           // number of elements
  int onebody;             // one-body descriptors
  int besseldegree;        // degree of Bessel functions
  int inversedegree;       // degree of inverse functions
  int nbesselpars;         // number of Bessel parameters
  int nCoeffPerElement;    // number of coefficients per element = (nl1 + Mdesc*nClusters)
  int ns;                  // number of snapshots for radial basis functions
  int nl1, nl2, nl3, nl4, nl23, nl33, nl34, nl44, nl;    // number of local descriptors
  int nrbf2, nrbf3, nrbf4, nrbfmax;                      // number of radial basis functions
  int nabf3, nabf4;                                      // number of angular basis functions
  int K3, K4, Q4;                                        // number of monomials

  // environmental variables
  int nClusters;      // number of environment clusters
  int nComponents;    // number of principal components
  int Mdesc;          // number of base descriptors

  double rin;     // inner cut-off radius
  double rcut;    // outer cut-off radius
  double rmax;    // rcut - rin

  double *rij;    // (xj - xi) for all pairs (I, J)
  double *fij;    // force for all pairs (I, J)
  double *ei;     // energy for each atom I
  int *typeai;    // types of atoms I only
  int *numij;     // number of pairs (I, J) for each atom I
  int *idxi;      // storing linear indices of atom I for all pairs (I, J)
  int *ai;        // IDs of atoms I for all pairs (I, J)
  int *aj;        // IDs of atoms J for all pairs (I, J)
  int *ti;        // types of atoms I for all pairs (I, J)
  int *tj;        // types of atoms J  for all pairs (I, J)

  double besselparams[3];
  double *Phi;             // eigenvectors matrix ns x ns
  double *rbf;             // radial basis functions nij x nrbfmax
  double *rbfx;            // x-derivatives of radial basis functions nij x nrbfmax
  double *rbfy;            // y-derivatives of radial basis functions nij x nrbfmax
  double *rbfz;            // z-derivatives of radial basis functions nij x nrbfmax
  double *abf;             // angular basis functions nij x K3
  double *abfx;            // x-derivatives of angular basis functions nij x K3
  double *abfy;            // y-derivatives of angular basis functions nij x K3
  double *abfz;            // z-derivatives of angular basis functions nij x K3
  double *abftm;           // angular basis functions 4 x K3
  double *sumU;            // sum of radial basis functions ni x K3 x nrbfmax x nelements
  double *forcecoeff;      // force coefficients ni x K3 x nrbfmax x nelements
  double *Proj;            // PCA Projection matrix
  double *Centroids;       // centroids of the clusters
  double *bd;              // base descriptors ni x Mdesc
  double *cb;              // force coefficients for base descriptors ni x Mdesc
  double *pd;              // environment probability descriptors ni x nClusters
  double *bdd;             // base descriptors derivatives 3 x nij x Mdesc
  double *pdd;             // environment probability descriptors derivatives 3 x nij x nClusters
  double *coefficients;    // coefficients nCoeffPerElement x nelements
  int *pq3, *pn3, *pc3;    // arrays to compute 3-body angular basis functions
  int *pa4, *pb4, *pc4;    // arrays to compute 4-body angular basis functions
  int *ind33l, *ind33r;    // nl33
  int *ind34l, *ind34r;    // nl34
  int *ind44l, *ind44r;    // nl44
  int *elemindex;

  bool peratom_warn;    // print warning about missing per-atom energies or stresses
};

}    // namespace LAMMPS_NS

#endif
#endif
