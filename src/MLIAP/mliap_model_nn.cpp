/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "mliap_model_nn.h"
#include "pair_mliap.h"
#include "mliap_data.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MLIAPModelNN::MLIAPModelNN(LAMMPS* lmp, char* coefffilename) :
  MLIAPModelSimple(lmp, coefffilename)
{
  ndescriptors = ndescriptors;
}

/* ---------------------------------------------------------------------- */

MLIAPModelNN::~MLIAPModelNN(){}

/* ----------------------------------------------------------------------
   get number of parameters
   ---------------------------------------------------------------------- */

int MLIAPModelNN::get_nparams()
{
  if (nparams == 0) {
    if (ndescriptors == 0) error->all(FLERR,"ndescriptors not defined");
  }
  return nparams;
}

/*  ----------------------------------------------------------------------
   Calculate model gradients w.r.t descriptors
   for each atom beta_i = dE(B_i)/dB_i
   ---------------------------------------------------------------------- */

void MLIAPModelNN::compute_gradients(MLIAPData* data)
{
  data->energy = 0.0;

  for (int ii = 0; ii < data->natoms; ii++) {
      const int ielem = data->ielems[ii];
      const int nl = data->nlayers;
    
      double* coeffi = coeffelem[ielem];
      double** scalei = scale[ielem];
      double **nodes, **dnodes, **bnodes;
      
      nodes = new double*[nl];
      dnodes = new double*[nl];
      bnodes = new double*[nl];
      for (int l=0; l<nl; ++l) {
        nodes[l] = new double[nnodes[l]];
        dnodes[l] = new double[nnodes[l]];
        bnodes[l] = new double[nnodes[l]];
      }
    
      // forwardprop
      // input - hidden1
      for (int n=0; n < nnodes[0]; n++) {
        nodes[0][n] = 0;
        for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
          nodes[0][n] += coeffi[n*((data->ndescriptors)+1)+icoeff+1] * (data->descriptors[ii][icoeff] - scalei[0][icoeff]) / scalei[1][icoeff];
        }
        if (activation[0] == 1) {
          nodes[0][n] = sigm(nodes[0][n] + coeffi[n*((data->ndescriptors)+1)], dnodes[0][n]);
        }
        else if (activation[0] == 2) {
          nodes[0][n] = tanh(nodes[0][n] + coeffi[n*((data->ndescriptors)+1)], dnodes[0][n]);
        }
        else if (activation[0] == 3) {
          nodes[0][n] = relu(nodes[0][n] + coeffi[n*((data->ndescriptors)+1)], dnodes[0][n]);
        }
        else {
          nodes[0][n] += coeffi[n*((data->ndescriptors)+1)];
          dnodes[0][n] = 1;
        }
      }

      // hidden~output
      int k = 0;
      if (nl > 1) {
        k += ((data->ndescriptors)+1)*nnodes[0];
        for (int l=1; l < nl; l++) {
          for (int n=0; n < nnodes[l]; n++) {
            nodes[l][n] = 0;
            for (int j=0; j < nnodes[l-1]; j++) {
              nodes[l][n] += coeffi[k+n*(nnodes[l-1]+1)+j+1] * nodes[l-1][j];
            }
            if (activation[l] == 1) {
              nodes[l][n] = sigm(nodes[l][n] + coeffi[k+n*(nnodes[l-1]+1)], dnodes[l][n]);
            }
            else if (activation[l] == 2) {
              nodes[l][n] = tanh(nodes[l][n] + coeffi[k+n*(nnodes[l-1]+1)], dnodes[l][n]);
            }
            else if (activation[l] == 3) {
              nodes[l][n] = relu(nodes[l][n] + coeffi[k+n*(nnodes[l-1]+1)], dnodes[l][n]);
            }
            else {
              nodes[l][n] += coeffi[k+n*(nnodes[l-1]+1)];
              dnodes[l][n] = 1;
            }
          }
	  k += (nnodes[l-1]+1)*nnodes[l];
        }
      }

      // backwardprop
      // output layer dnode initialized to 1.
      
      for (int n=0; n<nnodes[nl-1]; n++) {
        if (activation[nl-1] == 0) {
            bnodes[nl-1][n] = 1;
        }
        else {
            bnodes[nl-1][n] = dnodes[nl-1][n];
        }
      }

      if (nl > 1) {
        for (int l=nl-1; l>0; l--) {
          k -= (nnodes[l-1]+1)*nnodes[l];
          for (int n=0; n<nnodes[l-1]; n++) {
            bnodes[l-1][n] = 0;
            for (int j=0; j<nnodes[l]; j++) {
              bnodes[l-1][n] += coeffi[k+j*(nnodes[l-1]+1)+n+1] * bnodes[l][j];
            }
            if (activation[l-1] >= 1) {
              bnodes[l-1][n] *= dnodes[l-1][n];
            }
          }
        }
      }
      
      for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
        data->betas[ii][icoeff] = 0;
        for (int j=0; j<nnodes[0]; j++) {
          data->betas[ii][icoeff] += coeffi[j*((data->ndescriptors)+1)+icoeff+1] * bnodes[0][j];
        }
	data->betas[ii][icoeff] = data->betas[ii][icoeff]/scalei[1][icoeff];
      }
      
      if (data->eflag) {

      // energy of atom I (E_i)
      
        double etmp = nodes[nl-1][0];
          
        data->energy += etmp;
        data->eatoms[ii] = etmp;
      }
      // Deleting the variables
      
      for (int n=0; n<nl; n++) {
        delete nodes[n];
        delete dnodes[n];
        delete bnodes[n];
      }
          
      delete nodes;
      delete dnodes;
      delete bnodes;
      
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

void MLIAPModelNN::compute_gradgrads(class MLIAPData* data)
{
  error->all(FLERR,"compute_gradgrads not implemented");
}

/* ----------------------------------------------------------------------
   calculate gradients of forces w.r.t. parameters
   egradient is derivative of energy w.r.t. parameters
   ---------------------------------------------------------------------- */

void MLIAPModelNN::compute_force_gradients(class MLIAPData* data)
{
  error->all(FLERR,"compute_force_gradients not implemented");
}

/* ----------------------------------------------------------------------
   count the number of non-zero entries in gamma matrix
   ---------------------------------------------------------------------- */

int MLIAPModelNN::get_gamma_nnz(class MLIAPData* data)
{
  // todo: get_gamma_nnz
  return 0;
}

