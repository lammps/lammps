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

#include "mliap_model_linear.h"
#include "pair_mliap.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MLIAPModelLinear::MLIAPModelLinear(LAMMPS* lmp, char* coefffilename) :
  MLIAPModel(lmp, coefffilename)
{
  ndescriptors = nparams - 1;
}

/* ---------------------------------------------------------------------- */

MLIAPModelLinear::MLIAPModelLinear(LAMMPS* lmp, int nelements_in, int nparams_in) : 
  MLIAPModel(lmp, nelements_in, nparams_in)
{
  ndescriptors = nparams - 1;
}

/* ---------------------------------------------------------------------- */

MLIAPModelLinear::~MLIAPModelLinear(){}

/* ----------------------------------------------------------------------
   Calculate model gradients w.r.t descriptors 
   for each atom beta_i = dE(B_i)/dB_i
   ---------------------------------------------------------------------- */

void MLIAPModelLinear::gradient(int natomdesc, int *iatomdesc, int *ielemdesc, 
                                double **descriptors, double **beta, PairMLIAP* pairmliap, int eflag)
{
  for (int ii = 0; ii < natomdesc; ii++) {
    const int i = iatomdesc[ii];
    const int ielem = ielemdesc[ii];

    double* coeffi = coeffelem[ielem];
    for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
      beta[ii][icoeff] = coeffi[icoeff+1];

    // add in contributions to global and per-atom energy
    // this is optional and has no effect on force calculation

    if (eflag) {

      // energy of atom I

      double* coeffi = coeffelem[ielem];
      double etmp = coeffi[0];

      // E_i = beta.B_i

      for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
        etmp += coeffi[icoeff+1]*descriptors[ii][icoeff];

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

void MLIAPModelLinear::param_gradient(int natommliap, int *iatommliap, int *ielemmliap,
                                         double **descriptors, 
                                         int **gamma_row_index, int **gamma_col_index, 
                                         double **gamma, double *egradient)
{
  // zero out energy gradients

  for (int l = 0; l < nelements*nparams; l++)
    egradient[l] = 0.0;
    
  for (int ii = 0; ii < natommliap; ii++) {
    const int i = iatommliap[ii];
    const int ielem = ielemmliap[ii];
    const int elemoffset = nparams*ielem;

    int l = elemoffset+1;
    for (int icoeff = 0; icoeff < ndescriptors; icoeff++) {
      gamma[ii][icoeff] = 1.0;
      gamma_row_index[ii][icoeff] = l++;
      gamma_col_index[ii][icoeff] = icoeff;
    }

    // gradient of energy of atom I w.r.t. parameters
    
    l = elemoffset;
    egradient[l++] += 1.0;
    for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
      egradient[l++] += descriptors[ii][icoeff];
    
  }

}

/* ----------------------------------------------------------------------
   count the number of non-zero entries in gamma matrix
   ---------------------------------------------------------------------- */

int MLIAPModelLinear::get_gamma_nnz()
{
  int inz = ndescriptors;
  return inz;
}

void MLIAPModelLinear::compute_force_gradients(double **descriptors, int numlistdesc, 
                                               int *iatomdesc, int *ielemdesc, int *numneighdesc, 
                                               int *jatomdesc, int *jelemdesc, double ***graddesc, 
                                               int yoffset, int zoffset, double **gradforce, 
                                               double *egradient) {


  // zero out energy gradients

  for (int l = 0; l < nelements*nparams; l++)
    egradient[l] = 0.0;
    
  int ij = 0;
  for (int ii = 0; ii < numlistdesc; ii++) {
    const int i = iatomdesc[ii];
    const int ielem = ielemdesc[ii];
    const int elemoffset = nparams*ielem;

    for (int jj = 0; jj < numneighdesc[ii]; jj++) {
      const int j = jatomdesc[ij];
      const int jelem = ielemdesc[ij];

      int l = elemoffset+1;
      for (int icoeff = 0; icoeff < ndescriptors; icoeff++) {
        gradforce[i][l]         += graddesc[ij][icoeff][0];
        gradforce[i][l+yoffset] += graddesc[ij][icoeff][1];
        gradforce[i][l+zoffset] += graddesc[ij][icoeff][2];
        gradforce[j][l]         -= graddesc[ij][icoeff][0];
        gradforce[j][l+yoffset] -= graddesc[ij][icoeff][1];
        gradforce[j][l+zoffset] -= graddesc[ij][icoeff][2];
        l++;
      }
      ij++;
    }

    // gradient of energy of atom I w.r.t. parameters
    
    int l = elemoffset;
    egradient[l++] += 1.0;
    for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
      egradient[l++] += descriptors[ii][icoeff];
    
  }
  
}
