// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "mliap_model_quadratic.h"

#include "mliap_data.h"
#include "error.h"
#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MLIAPModelQuadratic::MLIAPModelQuadratic(LAMMPS* lmp, char* coefffilename) :
  MLIAPModelSimple(lmp, coefffilename)
{
  if (coefffilename) read_coeffs(coefffilename);
  if (nparams > 0) ndescriptors = sqrt(2*nparams)-1;
  nonlinearflag = 1;
}

/* ---------------------------------------------------------------------- */

MLIAPModelQuadratic::~MLIAPModelQuadratic() {}

/* ----------------------------------------------------------------------
   get number of parameters
   ---------------------------------------------------------------------- */

int MLIAPModelQuadratic::get_nparams()
{
  if (nparams == 0) {
    if (ndescriptors == 0) error->all(FLERR,"ndescriptors not defined");
    else nparams = ndescriptors + 1 + (ndescriptors*(ndescriptors+1))/2;
  }

  return nparams;
}

/* ----------------------------------------------------------------------
   Calculate model gradients w.r.t descriptors for each atom dE(B_i)/dB_i
   ---------------------------------------------------------------------- */

void MLIAPModelQuadratic::compute_gradients(MLIAPData* data)
{
  data->energy = 0.0;

  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int ielem = data->ielems[ii];

    double* coeffi = coeffelem[ielem];
    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++)
      data->betas[ii][icoeff] = coeffi[icoeff+1];

    int k = ndescriptors+1;
    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
      double bveci = data->descriptors[ii][icoeff];
      data->betas[ii][icoeff] += coeffi[k]*bveci;
      k++;
      for (int jcoeff = icoeff+1; jcoeff < data->ndescriptors; jcoeff++) {
        double bvecj = data->descriptors[ii][jcoeff];
        data->betas[ii][icoeff] += coeffi[k]*bvecj;
        data->betas[ii][jcoeff] += coeffi[k]*bveci;
        k++;
      }
    }

    // add in contributions to global and per-atom energy
    // this is optional and has no effect on force calculation

    if (data->eflag) {

      // energy of atom I

      double* coeffi = coeffelem[ielem];
      double etmp = coeffi[0];

      // E_i = beta.B_i + 0.5*B_i^t.alpha.B_i

      for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++)
        etmp += coeffi[icoeff+1]*data->descriptors[ii][icoeff];

      // quadratic contributions

      int k = ndescriptors+1;
      for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
        double bveci = data->descriptors[ii][icoeff];
        etmp += 0.5*coeffi[k++]*bveci*bveci;
        for (int jcoeff = icoeff+1; jcoeff < data->ndescriptors; jcoeff++) {
          double bvecj = data->descriptors[ii][jcoeff];
          etmp += coeffi[k++]*bveci*bvecj;
        }
      }
      data->energy += etmp;
      data->eatoms[ii] = etmp;
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

void MLIAPModelQuadratic::compute_gradgrads(class MLIAPData* data)
{
  // zero out energy gradients

  for (int l = 0; l < data->nelements*data->nparams; l++)
    data->egradient[l] = 0.0;

  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int ielem = data->ielems[ii];
    const int elemoffset = data->nparams*ielem;

    // linear contributions

    int l = elemoffset+1;
    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
      data->gamma[ii][icoeff] = 1.0;
      data->gamma_row_index[ii][icoeff] = l++;
      data->gamma_col_index[ii][icoeff] = icoeff;
    }

    // quadratic contributions

    int inz = data->ndescriptors;
    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
      double bveci = data->descriptors[ii][icoeff];
      data->gamma[ii][inz] = bveci;
      data->gamma_row_index[ii][inz] = l++;
      data->gamma_col_index[ii][inz] = icoeff;
      inz++;
      for (int jcoeff = icoeff+1; jcoeff < data->ndescriptors; jcoeff++) {
        double bvecj = data->descriptors[ii][jcoeff];
        data->gamma[ii][inz] = bvecj; // derivative w.r.t. B[icoeff]
        data->gamma_row_index[ii][inz] = l;
        data->gamma_col_index[ii][inz] = icoeff;
        inz++;
        data->gamma[ii][inz] = bveci; // derivative w.r.t. B[jcoeff]
        data->gamma_row_index[ii][inz] = l;
        data->gamma_col_index[ii][inz] = jcoeff;
        inz++;
        l++;
      }
    }

    // gradient of energy of atom I w.r.t. parameters

    l = elemoffset;
    data->egradient[l++] += 1.0;
    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++)
      data->egradient[l++] += data->descriptors[ii][icoeff];

    // quadratic contributions

    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
      double bveci = data->descriptors[ii][icoeff];
      data->egradient[l++] += 0.5*bveci*bveci;
      for (int jcoeff = icoeff+1; jcoeff < ndescriptors; jcoeff++) {
        double bvecj = data->descriptors[ii][jcoeff];
        data->egradient[l++] += bveci*bvecj;
      }
    }
  }

}

/* ----------------------------------------------------------------------
   count the number of non-zero entries in gamma matrix
   ---------------------------------------------------------------------- */

int MLIAPModelQuadratic::get_gamma_nnz(class MLIAPData* data)
{
  int inz = ndescriptors;
  for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
    inz++;
    for (int jcoeff = icoeff+1; jcoeff < data->ndescriptors; jcoeff++) {
        inz++;
        inz++;
    }
  }

  return inz;
}

void MLIAPModelQuadratic::compute_force_gradients(class MLIAPData* data) {
  // zero out energy gradients

  for (int l = 0; l < data->nelements*data->nparams; l++)
    data->egradient[l] = 0.0;

  int ij = 0;
  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int i = data->iatoms[ii];
    const int ielem = data->ielems[ii];
    const int elemoffset = data->nparams*ielem;

    for (int jj = 0; jj < data->numneighs[ii]; jj++) {
      const int j = data->jatoms[ij];
      int l = elemoffset+1;
      for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
        data->gradforce[i][l]               += data->graddesc[ij][icoeff][0];
        data->gradforce[i][l+data->yoffset] += data->graddesc[ij][icoeff][1];
        data->gradforce[i][l+data->zoffset] += data->graddesc[ij][icoeff][2];
        data->gradforce[j][l]               -= data->graddesc[ij][icoeff][0];
        data->gradforce[j][l+data->yoffset] -= data->graddesc[ij][icoeff][1];
        data->gradforce[j][l+data->zoffset] -= data->graddesc[ij][icoeff][2];
        l++;
      }

      // quadratic contributions

      for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
        double bveci = data->descriptors[ii][icoeff];
        data->gradforce[i][l]               += data->graddesc[ij][icoeff][0]*bveci;
        data->gradforce[i][l+data->yoffset] += data->graddesc[ij][icoeff][1]*bveci;
        data->gradforce[i][l+data->zoffset] += data->graddesc[ij][icoeff][2]*bveci;
        data->gradforce[j][l]               -= data->graddesc[ij][icoeff][0]*bveci;
        data->gradforce[j][l+data->yoffset] -= data->graddesc[ij][icoeff][1]*bveci;
        data->gradforce[j][l+data->zoffset] -= data->graddesc[ij][icoeff][2]*bveci;
        l++;
        for (int jcoeff = icoeff+1; jcoeff < data->ndescriptors; jcoeff++) {
          double bvecj = data->descriptors[ii][jcoeff];
          data->gradforce[i][l]               += data->graddesc[ij][icoeff][0]*bvecj + data->graddesc[ij][jcoeff][0]*bveci;
          data->gradforce[i][l+data->yoffset] += data->graddesc[ij][icoeff][1]*bvecj + data->graddesc[ij][jcoeff][1]*bveci;
          data->gradforce[i][l+data->zoffset] += data->graddesc[ij][icoeff][2]*bvecj + data->graddesc[ij][jcoeff][2]*bveci;
          data->gradforce[j][l]               -= data->graddesc[ij][icoeff][0]*bvecj + data->graddesc[ij][jcoeff][0]*bveci;
          data->gradforce[j][l+data->yoffset] -= data->graddesc[ij][icoeff][1]*bvecj + data->graddesc[ij][jcoeff][1]*bveci;
          data->gradforce[j][l+data->zoffset] -= data->graddesc[ij][icoeff][2]*bvecj + data->graddesc[ij][jcoeff][2]*bveci;
          l++;
        }
      }
      ij++;
    }

    // gradient of energy of atom I w.r.t. parameters

    int l = elemoffset;
    data->egradient[l++] += 1.0;
    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++)
      data->egradient[l++] += data->descriptors[ii][icoeff];

    // quadratic contributions

    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
      double bveci = data->descriptors[ii][icoeff];
      data->egradient[l++] += 0.5*bveci*bveci;
      for (int jcoeff = icoeff+1; jcoeff < data->ndescriptors; jcoeff++) {
        double bvecj = data->descriptors[ii][jcoeff];
        data->egradient[l++] += bveci*bvecj;
      }
    }

  }

}
