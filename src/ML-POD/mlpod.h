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

#ifndef LMP_MLPOD_H
#define LMP_MLPOD_H

#include "pointers.h"

#define DDOT ddot_
#define DGEMV dgemv_
#define DGEMM dgemm_
#define DGETRF dgetrf_
#define DGETRI dgetri_
#define DSYEV dsyev_
#define DPOSV dposv_

extern "C" {
double DDOT(int *, double *, int *, double *, int *);
void DGEMV(char *, int *, int *, double *, double *, int *, double *, int *, double *, double *,
           int *);
void DGEMM(char *, char *, int *, int *, int *, double *, double *, int *, double *, int *,
           double *, double *, int *);
void DGETRF(int *, int *, double *, int *, int *, int *);
void DGETRI(int *, double *, int *, int *, double *, int *, int *);
void DSYEV(char *, char *, int *, double *, int *, double *, double *, int *, int *);
void DPOSV(char *, int *, int *, double *, int *, double *, int *, int *);
}

namespace LAMMPS_NS {

class MLPOD : protected Pointers {

 private:
  // functions for reading input files

  void read_pod(const std::string &pod_file);
  void read_coeff_file(const std::string &coeff_file);

  // functions for calculating/collating POD descriptors/coefficients for energies

  void podradialbasis(double *rbf, double *drbf, double *xij, double *besselparams, double rin,
                      double rmax, int besseldegree, int inversedegree, int nbesselpars, int N);
  void pod1body(double *eatom, double *fatom, int *atomtype, int nelements, int natom);
  void podtally2b(double *eatom, double *fatom, double *eij, double *fij, int *ai, int *aj, int *ti,
                  int *tj, int *elemindex, int nelements, int nbf, int natom, int N);
  void pod3body(double *eatom, double *fatom, double *rij, double *e2ij, double *f2ij,
                double *tmpmem, int *elemindex, int *pairnumsum, int *ai, int *aj, int *ti, int *tj,
                int nrbf, int nabf, int nelements, int natom, int Nij);
  void poddesc(double *eatom1, double *fatom1, double *eatom2, double *fatom2, double *eatom3,
               double *fatom3, double *rij, double *Phi, double *besselparams, double *tmpmem,
               double rin, double rcut, int *pairnumsum, int *atomtype, int *ai, int *aj, int *ti,
               int *tj, int *elemindex, int *pdegree, int nbesselpars, int nrbf2, int nrbf3,
               int nabf, int nelements, int Nij, int natom);
  double quadratic_coefficients(double *c2, double *c3, double *d2, double *d3, double *coeff23,
                                int *quadratic, int nc2, int nc3);
  double quadratic_coefficients(double *c3, double *d3, double *coeff33, int *quadratic, int nc3);
  double cubic_coefficients(double *c2, double *c3, double *c4, double *d2, double *d3, double *d4,
                            double *coeff234, int *cubic, int nc2, int nc3, int nc4);
  double cubic_coefficients(double *c3, double *d3, double *coeff333, int *cubic, int nc3);
  double quadratic_coefficients(double *ce2, double *ce3, double *c2, double *c3, double *d2,
                                double *d3, double *coeff23, int *quadratic, int nc2, int nc3);
  double quadratic_coefficients(double *ce3, double *c3, double *d3, double *coeff33,
                                int *quadratic, int nc3);
  double cubic_coefficients(double *ce2, double *ce3, double *ce4, double *c2, double *c3,
                            double *c4, double *d2, double *d3, double *d4, double *coeff234,
                            int *cubic, int nc2, int nc3, int nc4);
  double cubic_coefficients(double *ce3, double *c3, double *d3, double *coeff333, int *cubic,
                            int nc3);

  // functions for calculating/collating SNAP descriptors/coefficients for energies

  void snapSetup(int twojmax, int ntypes);
  void InitSnap();
  void snapComputeUlist(double *Sr, double *Si, double *dSr, double *dSi, double *rootpqarray,
                        double *rij, double *wjelem, double *radelem, double rmin0, double rfac0,
                        double rcutfac, int *idxu_block, int *ti, int *tj, int twojmax,
                        int idxu_max, int ijnum, int switch_flag);
  void snapZeroUarraytot2(double *Stotr, double *Stoti, double wself, int *idxu_block, int *type,
                          int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max,
                          int nelements, int twojmax, int inum);
  void snapAddUarraytot(double *Stotr, double *Stoti, double *Sr, double *Si, int *map, int *ai,
                        int *tj, int idxu_max, int inum, int ijnum, int chemflag);
  void snapComputeZi2(double *zlist_r, double *zlist_i, double *Stotr, double *Stoti,
                      double *cglist, int *idxz, int *idxu_block, int *idxcg_block, int twojmax,
                      int idxu_max, int idxz_max, int nelements, int bnorm_flag, int inum);
  void snapComputeBi1(double *blist, double *zlist_r, double *zlist_i, double *Stotr, double *Stoti,
                      int *idxb, int *idxu_block, int *idxz_block, int twojmax, int idxb_max,
                      int idxu_max, int idxz_max, int nelements, int inum);
  void snapComputeDbidrj(double *dblist, double *zlist_r, double *zlist_i, double *dulist_r,
                         double *dulist_i, int *idxb, int *idxu_block, int *idxz_block, int *map,
                         int *ai, int *tj, int twojmax, int idxb_max, int idxu_max, int idxz_max,
                         int nelements, int bnorm_flag, int chemflag, int inum, int ijnum);
  void snapdesc(double *blist, double *bd, double *rij, double *tmpmem, int *atomtype, int *ai,
                int *aj, int *ti, int *tj, int natom, int Nij);

  // functions for calculating/collating POD descriptors/coefficients for forces

  void podradialbasis(double *rbf, double *xij, double *besselparams, double rin, double rmax,
                      int besseldegree, int inversedegree, int nbesselpars, int N);
  void pod1body(double *eatom, int *atomtype, int nelements, int natom);
  void podtally2b(double *eatom, double *eij, int *ai, int *ti, int *tj, int *elemindex,
                  int nelements, int nbf, int natom, int N);
  void pod3body(double *eatom, double *yij, double *e2ij, double *tmpmem, int *elemindex,
                int *pairnumsum, int *ai, int *ti, int *tj, int nrbf, int nabf, int nelements,
                int natom, int Nij);
  void poddesc_ij(double *eatom1, double *eatom2, double *eatom3, double *rij, double *Phi,
                  double *besselparams, double *tmpmem, double rin, double rcut, int *pairnumsum,
                  int *atomtype, int *ai, int *ti, int *tj, int *elemindex, int *pdegree,
                  int nbesselpars, int nrbf2, int nrbf3, int nabf, int nelements, int Nij,
                  int natom);
  void snapComputeUij(double *Sr, double *Si, double *rootpqarray, double *rij, double *wjelem,
                      double *radelem, double rmin0, double rfac0, double rcutfac, int *idxu_block,
                      int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag);
  void snapdesc_ij(double *blist, double *rij, double *tmpmem, int *atomtype, int *ai, int *ti,
                   int *tj, int natom, int Nij);
  void pod2body_force(double *force, double *fij, double *coeff2, int *ai, int *aj, int *ti,
                      int *tj, int *elemindex, int nelements, int nbf, int natom, int Nij);
  void pod3body_force(double *force, double *yij, double *e2ij, double *f2ij, double *coeff3,
                      double *tmpmem, int *elemindex, int *pairnumsum, int *ai, int *aj, int *ti,
                      int *tj, int nrbf, int nabf, int nelements, int natom, int Nij);
  void snapTallyForce(double *force, double *dbdr, double *coeff4, int *ai, int *aj, int *ti,
                      int ijnum, int ncoeff, int ntype);
  void pod4body_force(double *force, double *rij, double *coeff4, double *tmpmem, int *atomtype,
                      int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij);
  void pod2body_force(double **force, double *fij, double *coeff2, int *ai, int *aj, int *ti,
                      int *tj, int *elemindex, int nelements, int nbf, int natom, int Nij);
  void pod3body_force(double **force, double *yij, double *e2ij, double *f2ij, double *coeff3,
                      double *tmpmem, int *elemindex, int *pairnumsum, int *ai, int *aj, int *ti,
                      int *tj, int nrbf, int nabf, int nelements, int natom, int Nij);
  void snapTallyForce(double **force, double *dbdr, double *coeff4, int *ai, int *aj, int *ti,
                      int ijnum, int ncoeff, int ntype);
  void pod4body_force(double **force, double *rij, double *coeff4, double *tmpmem, int *atomtype,
                      int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij);

  // eigenproblem functions

  void podeigenvaluedecomposition(double *Phi, double *Lambda, double *besselparams, double rin,
                                  double rcut, int besseldegree, int inversedegree, int nbesselpars,
                                  int N);

 public:
  MLPOD(LAMMPS *, const std::string &pod_file, const std::string &coeff_file);

  MLPOD(LAMMPS *lmp) : Pointers(lmp){};
  ~MLPOD() override;

  struct podstruct {
    podstruct();
    virtual ~podstruct();

    std::vector<std::string> species;
    int twobody[3];
    int threebody[4];
    int fourbody[4];
    int *pbc;
    int *elemindex;

    int quadratic22[2];
    int quadratic23[2];
    int quadratic24[2];
    int quadratic33[2];
    int quadratic34[2];
    int quadratic44[2];
    int cubic234[3];
    int cubic333[3];
    int cubic444[3];
    int nelements;
    int onebody;
    int besseldegree;
    int inversedegree;

    int quadraticpod;

    double rin;
    double rcut;
    double *besselparams;
    double *coeff;
    double *Phi2, *Phi3, *Phi4, *Lambda2, *Lambda3, *Lambda4;

    // variables declaring number of snapshots, descriptors, and combinations

    int nbesselpars = 3;
    int ns2, ns3,
        ns4;    // number of snapshots for radial basis functions for linear POD potentials
    int nc2, nc3, nc4;             // number of chemical  combinations for linear POD potentials
    int nbf1, nbf2, nbf3, nbf4;    // number of basis functions for linear POD potentials
    int nd1, nd2, nd3, nd4;        // number of descriptors for linear POD potentials
    int nd22, nd23, nd24, nd33, nd34, nd44;    // number of descriptors for quadratic POD potentials
    int nd234, nd333, nd444;                   // number of descriptors for cubic POD potentials
    int nrbf3, nabf3, nrbf4, nabf4;
    int nd, nd1234;

    int snaptwojmax;    // also used to tell if SNAP is used when allocating/deallocating
    int snapchemflag;
    double snaprfac0;
    double snapelementradius[10];
    double snapelementweight[10];
  };

  struct snastruct {
    int twojmax;
    int ncoeff;
    int idxb_max;
    int idxu_max;
    int idxz_max;
    int idxcg_max;
    int ntypes;
    int nelements;
    int ndoubles;    // number of multi-element pairs
    int ntriples;    // number of multi-element triplets
    int bnormflag;
    int chemflag;
    int switchflag;
    int bzeroflag;
    int wselfallflag;

    double wself;
    double rmin0;
    double rfac0;
    double rcutfac;
    double rcutmax;

    int *map;    // map types to [0,nelements)
    int *idx_max;
    int *idxz;
    int *idxz_block;
    int *idxb;
    int *idxb_block;
    int *idxu_block;
    int *idxcg_block;

    double *rcutsq;
    double *radelem;
    double *wjelem;
    double *bzero;
    double *fac;
    double *rootpqarray;
    double *cglist;
  };

  podstruct pod;
  snastruct sna;

  // functions for collecting/collating arrays

  void podMatMul(double *c, double *a, double *b, int r1, int c1, int c2);
  void podArraySetValue(double *y, double a, int n);
  void podArrayCopy(double *y, double *x, int n);
  void podArrayFill(int *output, int start, int length);

  // functions for calculating energy and force descriptors

  void podNeighPairs(double *xij, double *x, int *ai, int *aj, int *ti, int *tj, int *pairlist,
                     int *pairnumsum, int *atomtype, int *alist, int inum, int dim);
  void linear_descriptors(double *gd, double *efatom, double *y, double *tmpmem, int *atomtype,
                          int *alist, int *pairlist, int *pairnum, int *pairnumsum, int *tmpint,
                          int natom, int Nij);
  void quadratic_descriptors(double *d23, double *dd23, double *d2, double *d3, double *dd2,
                             double *dd3, int M2, int M3, int N);
  void quadratic_descriptors(double *d33, double *dd33, double *d3, double *dd3, int M3, int N);
  void cubic_descriptors(double *d234, double *dd234, double *d2, double *d3, double *d4,
                         double *dd2, double *dd3, double *dd4, int M2, int M3, int M4, int N);
  void cubic_descriptors(double *d333, double *Dd333, double *d3, double *Dd3, int M3, int N);
  double calculate_energyforce(double *force, double *gd, double *gdd, double *coeff, double *tmp,
                               int natom);
  double energyforce_calculation(double *f, double *gd, double *gdd, double *coeff, double *y,
                                 int *atomtype, int *alist, int *pairlist, int *pairnum,
                                 int *pairnumsum, int *tmpint, int natom, int Nij);

  // functions for calculating energies and forces

  void podNeighPairs(double *rij, double *x, int *idxi, int *ai, int *aj, int *ti, int *tj,
                     int *pairnumsum, int *atomtype, int *jlist, int *alist, int inum);
  int lammpsNeighPairs(double *rij, double **x, double rcutsq, int *idxi, int *ai, int *aj, int *ti,
                       int *tj, int *pairnumsum, int *atomtype, int *numneigh, int *ilist,
                       int **jlist, int inum);
  void linear_descriptors_ij(double *gd, double *eatom, double *rij, double *tmpmem,
                             int *pairnumsum, int *atomtype, int *ai, int *ti, int *tj, int natom,
                             int Nij);
  double calculate_energy(double *effectivecoeff, double *gd, double *coeff);
  double calculate_energy(double *energycoeff, double *forcecoeff, double *gd, double *gdall,
                          double *coeff);
  void calculate_force(double *force, double *effectivecoeff, double *rij, double *tmpmem,
                       int *pairnumsum, int *atomtype, int *idxi, int *ai, int *aj, int *ti,
                       int *tj, int natom, int Nij);
  void calculate_force(double **force, double *effectivecoeff, double *rij, double *tmpmem,
                       int *pairnumsum, int *atomtype, int *idxi, int *ai, int *aj, int *ti,
                       int *tj, int natom, int Nij);
  double energyforce_calculation(double *force, double *podcoeff, double *effectivecoeff,
                                 double *gd, double *rij, double *tmpmem, int *pairnumsum,
                                 int *atomtype, int *idxi, int *ai, int *aj, int *ti, int *tj,
                                 int natom, int Nij);

};

}    // namespace LAMMPS_NS

#endif
