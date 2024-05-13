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

/* ----------------------------------------------------------------------
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "mliap_model_linear.h"

#include "mliap_data.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MLIAPModelLinear::MLIAPModelLinear(LAMMPS* lmp, char* coefffilename) :
  MLIAPModelSimple(lmp, coefffilename)
{
  if (nparams > 0) ndescriptors = nparams - 1;
}

/* ----------------------------------------------------------------------
   get number of parameters
   ---------------------------------------------------------------------- */

int MLIAPModelLinear::get_nparams()
{
  if (nparams == 0) {
    if (ndescriptors == 0) error->all(FLERR,"ndescriptors not defined");
    else nparams = ndescriptors + 1;
  }

  return nparams;
}

/* ----------------------------------------------------------------------
   Calculate model gradients w.r.t descriptors
   for each atom beta_i = dE(B_i)/dB_i
   ---------------------------------------------------------------------- */

void MLIAPModelLinear::compute_gradients(MLIAPData* data)
{
  data->energy = 0.0;

  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int ielem = data->ielems[ii];

    double const* coeffi = coeffelem[ielem];
    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++)
      data->betas[ii][icoeff] = coeffi[icoeff+1];

    // add in contributions to global and per-atom energy
    // this is optional and has no effect on force calculation
    if (data->eflag) {
      // energy of atom I
      double etmp = coeffi[0];

      // E_i = beta.B_i
      for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++)
        etmp += coeffi[icoeff+1]*data->descriptors[ii][icoeff];

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

void MLIAPModelLinear::compute_gradgrads(class MLIAPData* data)
{
  // zero out energy gradients

  for (int l = 0; l < data->nelements*data->nparams; l++)
    data->egradient[l] = 0.0;

  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int ielem = data->ielems[ii];
    const int elemoffset = data->nparams*ielem;

    int l = elemoffset+1;
    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
      data->gamma[ii][icoeff] = 1.0;
      data->gamma_row_index[ii][icoeff] = l++;
      data->gamma_col_index[ii][icoeff] = icoeff;
    }

    // gradient of energy of atom I w.r.t. parameters

    l = elemoffset;
    data->egradient[l++] += 1.0;
    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++)
      data->egradient[l++] += data->descriptors[ii][icoeff];

  }

}

/* ----------------------------------------------------------------------
   calculate gradients of forces w.r.t. parameters
   egradient is derivative of energy w.r.t. parameters
   ---------------------------------------------------------------------- */

void MLIAPModelLinear::compute_force_gradients(class MLIAPData* data)
{

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
        data->gradforce[i][l]         += data->graddesc[ij][icoeff][0];
        data->gradforce[i][l+data->yoffset] += data->graddesc[ij][icoeff][1];
        data->gradforce[i][l+data->zoffset] += data->graddesc[ij][icoeff][2];
        data->gradforce[j][l]         -= data->graddesc[ij][icoeff][0];
        data->gradforce[j][l+data->yoffset] -= data->graddesc[ij][icoeff][1];
        data->gradforce[j][l+data->zoffset] -= data->graddesc[ij][icoeff][2];
        l++;
      }
      ij++;
    }

    // gradient of energy of atom I w.r.t. parameters

    int l = elemoffset;
    data->egradient[l++] += 1.0;
    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++)
      data->egradient[l++] += data->descriptors[ii][icoeff];

  }

}

/* ----------------------------------------------------------------------
   count the number of non-zero entries in gamma matrix
   ---------------------------------------------------------------------- */

int MLIAPModelLinear::get_gamma_nnz(class MLIAPData* data)
{
  int inz = data->ndescriptors;
  return inz;
}

