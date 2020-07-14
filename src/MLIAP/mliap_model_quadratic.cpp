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

#include "mliap_model_quadratic.h"
#include "pair_mliap.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define MAXWORD 3

/* ---------------------------------------------------------------------- */

MLIAPModelQuadratic::MLIAPModelQuadratic(LAMMPS* lmp, char* coefffilename) :
  MLIAPModel(lmp, coefffilename)
{
  nonlinearflag = 1;
  ndescriptors = sqrt(2*nparams)-1;
}

/* ---------------------------------------------------------------------- */

MLIAPModelQuadratic::MLIAPModelQuadratic(LAMMPS* lmp, int nelements_in, int nparams_in) : 
  MLIAPModel(lmp, nelements_in, nparams_in)
{
  nonlinearflag = 1;
  ndescriptors = sqrt(2*nparams)-1;
}

/* ---------------------------------------------------------------------- */

MLIAPModelQuadratic::~MLIAPModelQuadratic(){}

/* ----------------------------------------------------------------------
   Calculate model gradients w.r.t descriptors for each atom dE(B_i)/dB_i
   ---------------------------------------------------------------------- */

void MLIAPModelQuadratic::gradient(PairMLIAP* pairmliap, NeighList* list, double **descriptors, double **beta, int eflag)
{
  int i;
  int *type = atom->type;

  for (int ii = 0; ii < list->inum; ii++) {
    i = list->ilist[ii];
    const int itype = type[i];
    const int ielem = pairmliap->map[itype];
    double* coeffi = coeffelem[ielem];

    for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
      beta[ii][icoeff] = coeffi[icoeff+1];

    int k = ndescriptors+1;
    for (int icoeff = 0; icoeff < ndescriptors; icoeff++) {
      double bveci = descriptors[ii][icoeff];
      beta[ii][icoeff] += coeffi[k]*bveci;
      k++;
      for (int jcoeff = icoeff+1; jcoeff < ndescriptors; jcoeff++) {
        double bvecj = descriptors[ii][jcoeff];
        beta[ii][icoeff] += coeffi[k]*bvecj;
        beta[ii][jcoeff] += coeffi[k]*bveci;
        k++;
      }
    }

    // add in contributions to global and per-atom energy
    // this is optional and has no effect on force calculation

    if (eflag) {

      // energy of atom I

      double* coeffi = coeffelem[ielem];
      double etmp = coeffi[0];

      // E_i = beta.B_i + 0.5*B_i^t.alpha.B_i

      for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
        etmp += coeffi[icoeff+1]*descriptors[ii][icoeff];

      // quadratic contributions

      int k = ndescriptors+1;
      for (int icoeff = 0; icoeff < ndescriptors; icoeff++) {
        double bveci = descriptors[ii][icoeff];
        etmp += 0.5*coeffi[k++]*bveci*bveci;
        for (int jcoeff = icoeff+1; jcoeff < ndescriptors; jcoeff++) {
          double bvecj = descriptors[ii][jcoeff];
          etmp += coeffi[k++]*bveci*bvecj;
        }
      }
      pairmliap->e_tally(i,etmp);
    }
  }
}


/* ----------------------------------------------------------------------
   Calculate model double gradients w.r.t descriptors and parameters
   for each atom energy gamma_lk = d2E(B)/dB_k/dsigma_l, 
   where sigma_l is a parameter, B_k a descriptor, 
   and atom subscript i is omitted

   gamma is in CSR format:
      nnz = number of non-zero values
      gamma_row_index[inz] = l indices, 0 <= l < nparams 
      gamma_col_indexiinz] = k indices, 0 <= k < ndescriptors
      gamma[i][inz] = non-zero values, 0 <= inz < nnz

   egradient is derivative of energy w.r.t. parameters
   ---------------------------------------------------------------------- */

void MLIAPModelQuadratic::param_gradient(int *map, NeighList* list, 
                                         double **descriptors, 
                                         int **gamma_row_index, int **gamma_col_index, 
                                         double **gamma, double *egradient)
{
  int i;
  int *type = atom->type;

  // zero out energy gradients

  for (int l = 0; l < nelements*nparams; l++)
    egradient[l] = 0.0;
    
  for (int ii = 0; ii < list->inum; ii++) {

    i = list->ilist[ii];
    const int itype = type[i];
    const int ielem = map[itype];
    const int elemoffset = nparams*ielem;

    // linear contributions

    int l = elemoffset+1;
    for (int icoeff = 0; icoeff < ndescriptors; icoeff++) {
      gamma[ii][icoeff] = 1.0;
      gamma_row_index[ii][icoeff] = l++;
      gamma_col_index[ii][icoeff] = icoeff;
    }

    // quadratic contributions

    int inz = ndescriptors;
    for (int icoeff = 0; icoeff < ndescriptors; icoeff++) {
      double bveci = descriptors[ii][icoeff];
      gamma[ii][inz] = bveci;
      gamma_row_index[ii][inz] = l++;
      gamma_col_index[ii][inz] = icoeff;
      inz++;
      for (int jcoeff = icoeff+1; jcoeff < ndescriptors; jcoeff++) {
        double bvecj = descriptors[ii][jcoeff];
        gamma[ii][inz] = bvecj; // derivative w.r.t. B[icoeff]
        gamma_row_index[ii][inz] = l;
        gamma_col_index[ii][inz] = icoeff;
        inz++;
        gamma[ii][inz] = bveci; // derivative w.r.t. B[jcoeff]
        gamma_row_index[ii][inz] = l;
        gamma_col_index[ii][inz] = jcoeff;
        inz++;
        l++;
      }
    }

    // gradient of energy of atom I w.r.t. parameters
    
    l = elemoffset;
    egradient[l++] += 1.0;
    for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
      egradient[l++] += descriptors[ii][icoeff];
    
    // quadratic contributions
    
    for (int icoeff = 0; icoeff < ndescriptors; icoeff++) {
      double bveci = descriptors[ii][icoeff];
      egradient[l++] += 0.5*bveci*bveci;
      for (int jcoeff = icoeff+1; jcoeff < ndescriptors; jcoeff++) {
        double bvecj = descriptors[ii][jcoeff];
        egradient[l++] += bveci*bvecj;
      }
    }
  }

}

/* ----------------------------------------------------------------------
   count the number of non-zero entries in gamma matrix
   ---------------------------------------------------------------------- */

int MLIAPModelQuadratic::get_gamma_nnz()
{
  int inz = ndescriptors;
  for (int icoeff = 0; icoeff < ndescriptors; icoeff++) {
    inz++;
    for (int jcoeff = icoeff+1; jcoeff < ndescriptors; jcoeff++) {
        inz++;
        inz++;
    }
  }

  return inz;
}

