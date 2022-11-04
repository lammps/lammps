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
   Contributing author: Pedro Antonio Santos Flórez (UNLV)
------------------------------------------------------------------------- */

#include "mliap_model_nn.h"

#include "mliap_data.h"

#include "comm.h"
#include "error.h"
#include "memory.h"
#include "tokenizer.h"

#include <cstring>

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

MLIAPModelNN::MLIAPModelNN(LAMMPS *_lmp, char *coefffilename) : MLIAPModel(_lmp, coefffilename)
{
  nnodes = nullptr;
  activation = nullptr;
  scale = nullptr;
  if (coefffilename) MLIAPModelNN::read_coeffs(coefffilename);
  nonlinearflag = 1;
}

/* ---------------------------------------------------------------------- */

MLIAPModelNN::~MLIAPModelNN()
{
  memory->destroy(nnodes);
  memory->destroy(activation);
  memory->destroy(scale);
}

/* ----------------------------------------------------------------------
   get number of parameters
   ---------------------------------------------------------------------- */

int MLIAPModelNN::get_nparams()
{
  if (nparams == 0)
    if (ndescriptors == 0) error->all(FLERR, "ndescriptors not defined");

  return nparams;
}

void MLIAPModelNN::read_coeffs(char *coefffilename)
{

  // open coefficient file on proc 0

  FILE *fpcoeff;
  if (comm->me == 0) {
    fpcoeff = utils::open_potential(coefffilename, lmp, nullptr);
    if (fpcoeff == nullptr)
      error->one(FLERR, "Cannot open MLIAPModel coeff file {}: {}", coefffilename,
                 utils::getsyserror());
  }

  char line[MAXLINE], *ptr;
  int n, eof = 0, nwords = 0;
  while (nwords == 0) {
    if (comm->me == 0) {
      ptr = fgets(line, MAXLINE, fpcoeff);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fpcoeff);
      } else
        n = strlen(line) + 1;
    }
    MPI_Bcast(&eof, 1, MPI_INT, 0, world);
    if (eof) break;
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    MPI_Bcast(line, n, MPI_CHAR, 0, world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line, '#'))) *ptr = '\0';
    nwords = utils::count_words(line);
  }
  if (nwords != 2) error->all(FLERR, "Incorrect format in MLIAPModel coefficient file");

  // words = ptrs to all words in line
  // strip single and double quotes from words

  try {
    ValueTokenizer coeffs(line);
    nelements = coeffs.next_int();
    nparams = coeffs.next_int();
  } catch (TokenizerException &e) {
    error->all(FLERR, "Incorrect format in MLIAPModel coefficient file: {}", e.what());
  }

  // set up coeff lists

  memory->destroy(coeffelem);
  memory->create(coeffelem, nelements, nparams, "mliap_snap_model:coeffelem");

  int stats = 0;
  int ielem = 0;
  int l = 0;

  while (true) {
    if (comm->me == 0) {
      ptr = fgets(line, MAXLINE, fpcoeff);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fpcoeff);
      } else
        n = strlen(line) + 1;
    }

    MPI_Bcast(&eof, 1, MPI_INT, 0, world);
    if (eof) break;
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    MPI_Bcast(line, n, MPI_CHAR, 0, world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line, '#'))) *ptr = '\0';

    nwords = utils::trim_and_count_words(line);
    if (nwords == 0) continue;

    ValueTokenizer values(line, "\"' \t\n\t\f");
    if (stats == 0) {    // Header NET
      auto tstr = values.next_string();
      if (tstr.substr(0, 3) != "NET")
        error->all(FLERR, "Incorrect format in MLIAPModel coefficient file");

      ndescriptors = values.next_int();
      nlayers = values.next_int();

      memory->create(activation, nlayers, "mliap_model:activation");
      memory->create(nnodes, nlayers, "mliap_model:nnodes");
      memory->create(scale, nelements, 2, ndescriptors, "mliap_model:scale");

      for (int ilayer = 0; ilayer < nlayers; ilayer++) {
        tstr = values.next_string();
        nnodes[ilayer] = values.next_int();

        if (tstr == "linear")
          activation[ilayer] = 0;
        else if (tstr == "sigmoid")
          activation[ilayer] = 1;
        else if (tstr == "tanh")
          activation[ilayer] = 2;
        else if (tstr == "relu")
          activation[ilayer] = 3;
        else
          activation[ilayer] = 4;
      }

      stats = 1;

    } else if (stats == 1) {
      scale[ielem][0][l] = values.next_double();
      for (int icoeff = 1; icoeff < nwords; icoeff++) {
        scale[ielem][0][l + icoeff] = values.next_double();
      }
      l += nwords;
      if (l == ndescriptors) {
        stats = 2;
        l = 0;
      }

    } else if (stats == 2) {
      scale[ielem][1][l] = values.next_double();
      for (int icoeff = 1; icoeff < nwords; icoeff++) {
        scale[ielem][1][l + icoeff] = values.next_double();
      }
      l += nwords;
      if (l == ndescriptors) {
        stats = 3;
        l = 0;
      }

      // set up coeff lists

    } else if (stats == 3) {
      if (nwords > 30) error->all(FLERR, "Incorrect format in MLIAPModel coefficient file");

      coeffelem[ielem][l] = values.next_double();
      for (int icoeff = 1; icoeff < nwords; icoeff++) {
        coeffelem[ielem][l + icoeff] = values.next_double();
      }
      l += nwords;
      if (l == nparams) {
        stats = 1;
        l = 0;
        ielem++;
      }
    }
  }
}

/*  ----------------------------------------------------------------------
   Calculate model gradients w.r.t descriptors
   for each atom beta_i = dE(B_i)/dB_i
   ---------------------------------------------------------------------- */

void MLIAPModelNN::compute_gradients(MLIAPData *data)
{
  data->energy = 0.0;

  for (int ii = 0; ii < data->nlistatoms; ii++) {
    const int ielem = data->ielems[ii];
    const int nl = nlayers;

    double *coeffi = coeffelem[ielem];
    double **scalei = scale[ielem];
    double **nodes, **dnodes, **bnodes;

    nodes = new double *[nl];
    dnodes = new double *[nl];
    bnodes = new double *[nl];
    for (int l = 0; l < nl; ++l) {
      nodes[l] = new double[nnodes[l]];
      dnodes[l] = new double[nnodes[l]];
      bnodes[l] = new double[nnodes[l]];
    }

    // forwardprop
    // input - hidden1

    for (int n = 0; n < nnodes[0]; n++) {
      nodes[0][n] = 0;
      for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
        nodes[0][n] += coeffi[n * ((data->ndescriptors) + 1) + icoeff + 1] *
            (data->descriptors[ii][icoeff] - scalei[0][icoeff]) / scalei[1][icoeff];
      }
      if (activation[0] == 1) {
        nodes[0][n] = sigm(nodes[0][n] + coeffi[n * ((data->ndescriptors) + 1)], dnodes[0][n]);
      } else if (activation[0] == 2) {
        nodes[0][n] = tanh(nodes[0][n] + coeffi[n * ((data->ndescriptors) + 1)], dnodes[0][n]);
      } else if (activation[0] == 3) {
        nodes[0][n] = relu(nodes[0][n] + coeffi[n * ((data->ndescriptors) + 1)], dnodes[0][n]);
      } else {
        nodes[0][n] += coeffi[n * ((data->ndescriptors) + 1)];
        dnodes[0][n] = 1;
      }
    }

    // hidden~output

    int k = 0;
    if (nl > 1) {
      k += ((data->ndescriptors) + 1) * nnodes[0];
      for (int l = 1; l < nl; l++) {
        for (int n = 0; n < nnodes[l]; n++) {
          nodes[l][n] = 0;
          for (int j = 0; j < nnodes[l - 1]; j++) {
            nodes[l][n] += coeffi[k + n * (nnodes[l - 1] + 1) + j + 1] * nodes[l - 1][j];
          }
          if (activation[l] == 1) {
            nodes[l][n] = sigm(nodes[l][n] + coeffi[k + n * (nnodes[l - 1] + 1)], dnodes[l][n]);
          } else if (activation[l] == 2) {
            nodes[l][n] = tanh(nodes[l][n] + coeffi[k + n * (nnodes[l - 1] + 1)], dnodes[l][n]);
          } else if (activation[l] == 3) {
            nodes[l][n] = relu(nodes[l][n] + coeffi[k + n * (nnodes[l - 1] + 1)], dnodes[l][n]);
          } else {
            nodes[l][n] += coeffi[k + n * (nnodes[l - 1] + 1)];
            dnodes[l][n] = 1;
          }
        }
        k += (nnodes[l - 1] + 1) * nnodes[l];
      }
    }

    // backwardprop
    // output layer dnode initialized to 1.

    for (int n = 0; n < nnodes[nl - 1]; n++) {
      if (activation[nl - 1] == 0) {
        bnodes[nl - 1][n] = 1;
      } else {
        bnodes[nl - 1][n] = dnodes[nl - 1][n];
      }
    }

    if (nl > 1) {
      for (int l = nl - 1; l > 0; l--) {
        k -= (nnodes[l - 1] + 1) * nnodes[l];
        for (int n = 0; n < nnodes[l - 1]; n++) {
          bnodes[l - 1][n] = 0;
          for (int j = 0; j < nnodes[l]; j++) {
            bnodes[l - 1][n] += coeffi[k + j * (nnodes[l - 1] + 1) + n + 1] * bnodes[l][j];
          }
          if (activation[l - 1] >= 1) { bnodes[l - 1][n] *= dnodes[l - 1][n]; }
        }
      }
    }

    for (int icoeff = 0; icoeff < data->ndescriptors; icoeff++) {
      data->betas[ii][icoeff] = 0;
      for (int j = 0; j < nnodes[0]; j++) {
        data->betas[ii][icoeff] +=
            coeffi[j * ((data->ndescriptors) + 1) + icoeff + 1] * bnodes[0][j];
      }
      data->betas[ii][icoeff] = data->betas[ii][icoeff] / scalei[1][icoeff];
    }

    if (data->eflag) {

      // energy of atom I (E_i)

      double etmp = nodes[nl - 1][0];

      data->energy += etmp;
      data->eatoms[ii] = etmp;
    }
    // Deleting the variables

    for (int n = 0; n < nl; n++) {
      delete[] nodes[n];
      delete[] dnodes[n];
      delete[] bnodes[n];
    }

    delete[] nodes;
    delete[] dnodes;
    delete[] bnodes;
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

void MLIAPModelNN::compute_gradgrads(class MLIAPData * /*data*/)
{
  error->all(FLERR, "compute_gradgrads not implemented");
}

/* ----------------------------------------------------------------------
   calculate gradients of forces w.r.t. parameters
   egradient is derivative of energy w.r.t. parameters
   ---------------------------------------------------------------------- */

void MLIAPModelNN::compute_force_gradients(class MLIAPData * /*data*/)
{
  error->all(FLERR, "compute_force_gradients not implemented");
}

/* ----------------------------------------------------------------------
   count the number of non-zero entries in gamma matrix
   ---------------------------------------------------------------------- */

int MLIAPModelNN::get_gamma_nnz(class MLIAPData * /*data*/)
{
  // todo: get_gamma_nnz
  return 0;
}

double MLIAPModelNN::memory_usage()
{
  double bytes = 0;

  bytes += (double) nelements * nparams * sizeof(double);             // coeffelem
  bytes += (double) nelements * 2 * ndescriptors * sizeof(double);    // scale
  bytes += (int) nlayers * sizeof(int);                               // nnodes
  bytes += (int) nlayers * sizeof(int);                               // activation
  return bytes;
}
