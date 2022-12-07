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

#ifndef LMP_FASTPOD_H
#define LMP_FASTPOD_H

#include "pointers.h"

#define DGEMM dgemm_
#define DSYEV dsyev_

extern "C" {
void DGEMM(char *, char *, int *, int *, int *, double *, double *, int *, double *, int *,
           double *, double *, int *);
void DSYEV(char *, char *, int *, double *, int *, double *, double *, int *, int *);
}

namespace LAMMPS_NS {

class FASTPOD : protected Pointers {
private:
  int indexmap3(int *indx, int n1, int n2, int n3, int N1, int N2);

  void init3bodyarray(int *np, int *pq, int *pc, int Pa3);

  void init4bodyarray(int *pa4, int *pb4, int *pc4, int Pa4);

  void init2body();

  void init3body(int Pa3);

  void init4body(int Pa4);

  void snapshots(double *rbf, double *xij, int N);

  void eigenvaluedecomposition(double *Phi, double *Lambda, int N);

  void myneighbors(double *rij, double *x, int *ai, int *aj, int *ti, int *tj,
        int *jlist, int *pairnumsum, int *atomtype, int *alist, int i);

  void twobodycoeff(double *newcoeff2, double *coeff2);

  double threebodycoeff(double *cU, double *coeff3, double *sumU, int N);

  double fourbodycoeff(double *cU, double *sumU, double *coeff4, int N);

  void radialfunctions(double *rbf, double *rij, double *besselparams, double rin,
        double rmax, int besseldegree, int inversedegree, int nbesselpars, int N);

  void radialbasis(double *rbf, double *rbfx, double *rbfy, double *rbfz, double *rij, double *besselparams, double rin,
        double rmax, int besseldegree, int inversedegree, int nbesselpars, int N);

  void orthogonalradialbasis(double *orthorbf, double *rij, double *Phi, double *besselparams,
        double rin, double rmax, int besseldegree, int inversedegree, int nbesselpars, int nrbf2, int N);

  void angularfunctions(double *abf, double *rij, double *tm, int *pq, int N, int K);

  void angularbasis(double *abf, double *abfx, double *abfy, double *abfz, double *rij, double *tm, int *pq, int N, int K);

  void radialangularfunctions(double *U, double *rbf, double *abf, int N, int K, int M);

  void radialangularbasis(double *U, double *Ux, double *Uy, double *Uz,
        double *rbf, double *rbfx, double *rbfy, double *rbfz, double *abf,
        double *abfx, double *abfy, double *abfz, int N, int K, int M);

  void sumradialangularfunctions(double *sumU, double *U, int *atomtype, int N, int K, int M, int Ne);

  void unifiedbasis(double *U, double *Ux, double *Uy, double *Uz, double *sumU, double *rij,
        double *Phi, double *besselparams, double *tmpmem, double rin, double rcut, int *pdegree,
        int *tj, int *pq, int nbesselpars, int nrbf, int K, int nelements, int Nj);

  void tallytwobodyglobdesc(double *gd, double *d, int *elemindex, int nrbf, int nelements, int ti);

  void tallytwobodyglobdescderiv(double *gdd, double *dd, int *ai, int *aj, int *ti, int *tj,
        int *elemindex, int nrbf, int nelements, int natom, int N);

  void tallyglobdesc(double *gd, double *d, int ndesc, int ti);

  void tallyglobdescderiv(double *gdd, double *dd,  int *ai, int *aj, int natom, int N, int ndesc, int ti);

  double tallytwobodylocalforce(double *fij, double *coeff2,  double *rbf, double *rbfx,
        double *rbfy, double *rbfz, int *tj, int nbf, int N);

  void tallylocalforce(double *fij, double *cU, double *Ux, double *Uy, double *Uz,
        int *atomtype, int N, int K, int M, int Ne);

  void MatMul(double *c, double *a, double *b, int r1, int c1, int c2);

  void scalarproduct(double *d, double c, int N);

  double dotproduct(double *c, double *d, int ndesc);

  void mvproduct(double *fij, double *c, double *dd, int N, int ndesc);

public:
  std::vector<std::string> species;

  double rin;
  double rcut;

  int nelements;
  int pbc[3];
  int *elemindex ;

  int onebody;   // one-body descriptors
  int besseldegree;
  int inversedegree;
  int pdegree[2];
  int nbesselpars;
  double besselparams[3];
  double *Phi ;    // eigenvectors
  double *Lambda ; // eigenvalues
  double *coeff;  // coefficients
  double *newcoeff ;  // coefficients
  double *tmpmem;

  int Njmax;
  int ncoeff;  // number of coefficients in the input file
  int ns;      // number of snapshots for radial basis functions
  int nd1, nd2, nd3, nd4, nd5, nd6, nd7, nd;   // number of global descriptors
  int nl1, nl2, nl3, nl4, nl5, nl6, nl7, nl;   // number of local descriptors
  int nrbf2, nrbf3, nrbf4;                     // number of radial basis functions
  int nabf3, nabf4;                            // number of angular basis functions
  int P3, P4;                                  // angular polynomial degrees
  int K3, K4, Q4;                              // number of monomials
  int *pn3, *pq3, *pc3;          // arrays to compute 3-body angular basis functions
  int *pq4, *pa4, *pb4, *pc4;// arrays to compute 3-body angular basis functions
  int *tmpint;
  int nintmem; // number of integers in tmpint array
  int ndblmem; // number of doubles in tmpmem array

  // four-body descriptors
  int *ind23, *ind32, nrbf23, nabf23, P23, n23, n32, nl23, nd23;

  // five-body descriptors
  int *ind33, nrbf33, nabf33, P33, n33, nl33, nd33;

  // six-body descriptors
  int *ind34, *ind43, nrbf34, nabf34, nabf43, P34, n34, n43, nl34, nd34;

  // seven-body descriptors
  int *ind44, nrbf44, nabf44, P44, n44, nl44, nd44;


  FASTPOD(LAMMPS *, const std::string &pod_file, const std::string &coeff_file);

  FASTPOD(LAMMPS *lmp) : Pointers(lmp){};
  ~FASTPOD() override;

  void print_matrix(const char* desc, int m, int n, int* a, int lda );
  void print_matrix(const char* desc, int m, int n, double* a, int lda );

  void read_pod_file(std::string pod_file);
  int read_coeff_file(std::string coeff_file);

  int estimate_memory(int Nj);

  void mknewcoeff();

  void mknewcoeff(double *c);

  void onebodydescriptors(double *gd1, double *gdd1, int *ti, int natom, int i);

  double onebodyenergy(double *coeff1, int *ti);

  void twobodydescderiv(double *d2, double *dd2, double *rbf, double *rbfx,
        double *rbfy, double *rbfz, int *tj, int N);

  void twobodydescriptors(double *d2, double *dd2, double *rij, double *tempmem, int *tj, int Nj);

  void twobodydescriptors(double *gd2, double *gdd2, double *d2, double *dd2, double *rij,
        double *tempmem, int *ai, int *aj, int *ti, int *tj, int Nj, int natom);

  double twobodyenergyforce(double *fij, double *rij, double *coeff2, double *tempmem, int *tj, int Nj);

  void threebodydesc(double *d3, double *sumU, int N);

  void threebodydescderiv(double *dd3, double *sumU, double *Ux, double *Uy, double *Uz,
        int *atomtype, int N);

  void threebodydescriptors(double *d3, double *dd3, double *rij, double *tempmem, int *tj, int Nj);

  void threebodydescriptors(double *gd3, double *gdd3, double *d3, double *dd3, double *rij,
        double *tempmem, int *ai, int *aj, int *ti, int *tj, int Nj, int natom);

  double threebodyenergyforce(double *fij, double *rij, double *coeff3, double *tempmem, int *tj, int Nj);

  void fourbodydescderiv(double *d4, double *dd4, double *sumU, double *Ux, double *Uy, double *Uz,
      int *atomtype, int N);

  void fourbodydescriptors(double *d4, double *dd4, double *rij, double *tempmem, int *tj, int Nj);

  void fourbodydescriptors(double *gd4, double *gdd4, double *d4, double *dd4, double *rij,
        double *tempmem, int *ai, int *aj, int *ti, int *tj, int Nj, int natom);

  double fourbodyenergyforce(double *fij, double *rij, double *coeff4, double *tempmem,
         int *tj, int Nj);

  void descriptors(double *gd, double *gdd, double *x, int *atomtype, int *alist,
          int *jlist, int *pairnumsum, int natom);

  double localenergyforce(double *fij, double *rij, double *tempmem, int *ti, int *tj, int Nj);

  double atomicenergyforce(double *fij, double *rij, double *tempmem, int *ti, int *tj, int Nj);

  double peratomenergyforce(double *fij, double *rij, double *tempmem, int *ti, int *tj, int Nj);

  double energyforce(double *force, double *x, int *atomtype, int *alist,
          int *jlist, int *pairnumsum, int natom);

  void tallyforce(double *force, double *fij,  int *ai, int *aj, int N);
  void tallyforce(double **force, double *fij,  int *ai, int *aj, int N);

  void fourbodydesc23(double* d23, double* d2, double *d3);
  void fourbodydescderiv23(double* dd23, double* d2, double *d3, double* dd2, double *dd3, int N);
  void fivebodydesc33(double* d33, double *d3);
  void fivebodydescderiv33(double* dd33, double *d3, double *dd3, int N);
  void sixbodydesc34(double* d34, double* d3, double *d4);
  void sixbodydescderiv34(double* dd34, double* d3, double *d4, double* dd3, double *dd4, int N);
  void sevenbodydesc44(double* d44, double *d4);
  void sevenbodydescderiv44(double* dd44, double *d4, double *dd4, int N);
  void fourbodyfij23(double *fij, double *cf, double *coeff23, double *d2, double *d3, double *dd2, double *dd3, int N);
  void fivebodyfij33(double *fij, double *cf, double *coeff33, double *d3, double *dd3, int N);
  void sixbodyfij34(double *fij, double *cf, double *coeff34, double *d3, double *d4, double *dd3, double *dd4, int N);
  void sevenbodyfij44(double *fij, double *cf, double *coeff44, double *d4, double *dd4, int N);

  void fourbodydescriptors23(double *gd23, double *gdd23, double *d23, double *dd23,
        double* d2, double *d3, double* dd2, double *dd3, int *ai, int *aj, int *ti, int *tj,
        int Nj, int natom);
  void fivebodydescriptors33(double *gd33, double *gdd33, double *d33, double *dd33,
        double *d3, double *dd3, int *ai, int *aj, int *ti, int *tj, int Nj, int natom);
  void sixbodydescriptors34(double *gd34, double *gdd34, double *d34, double *dd34,
        double* d3, double *d4, double* dd3, double *dd4, int *ai, int *aj, int *ti, int *tj,
        int Nj, int natom);
  void sevenbodydescriptors44(double *gd44, double *gdd44, double *d44, double *dd44,
        double *d4, double *dd4, int *ai, int *aj, int *ti, int *tj, int Nj, int natom);
};

}    // namespace LAMMPS_NS

#endif

