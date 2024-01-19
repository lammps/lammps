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

#include "mliap_model.h"

#include "comm.h"
#include "error.h"
#include "memory.h"
#include "tokenizer.h"

#include <cstring>

using namespace LAMMPS_NS;

static constexpr int MAXLINE = 1024;
static constexpr int MAXWORD = 3;

/* ---------------------------------------------------------------------- */

MLIAPModel::MLIAPModel(LAMMPS *lmp, char *) : Pointers(lmp), coeffelem(nullptr)
{
  nparams = 0;
  nelements = 0;
  ndescriptors = 0;
  nonlinearflag = 0;
}

/* ---------------------------------------------------------------------- */

MLIAPModel::~MLIAPModel()
{
  memory->destroy(coeffelem);
}

/* ----------------------------------------------------------------------
   placeholder
------------------------------------------------------------------------- */

void MLIAPModel::init() {}

/* ----------------------------------------------------------------------
   set number of elements
   ---------------------------------------------------------------------- */

void MLIAPModel::set_nelements(int nelements_in)
{
  nelements = nelements_in;
}

/* ----------------------------------------------------------------------
   set number of descriptors
   ---------------------------------------------------------------------- */

void MLIAPModel::set_ndescriptors(int ndescriptors_in)
{
  ndescriptors = ndescriptors_in;
}

/* ---------------------------------------------------------------------- */

MLIAPModelSimple::MLIAPModelSimple(LAMMPS *lmp, char *coefffilename) :
    MLIAPModel(lmp, coefffilename)
{
  if (coefffilename) MLIAPModelSimple::read_coeffs(coefffilename);
}

/* ---------------------------------------------------------------------- */

void MLIAPModelSimple::read_coeffs(char *coefffilename)
{

  // open coefficient file on proc 0

  FILE *fpcoeff;
  if (comm->me == 0) {
    fpcoeff = utils::open_potential(coefffilename, lmp, nullptr);
    if (fpcoeff == nullptr)
      error->one(FLERR, "Cannot open MLIAPModel coeff file {}: {}", coefffilename,
                 utils::getsyserror());
  }

  char line[MAXLINE] = {'\0'};
  char *ptr;
  int eof = 0;

  int n;
  int nwords = 0;
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
    error->all(FLERR,
               "Incorrect format in MLIAPModel coefficient "
               "file: {}",
               e.what());
  }

  // set up coeff lists

  memory->destroy(coeffelem);
  memory->create(coeffelem, nelements, nparams, "mliap_snap_model:coeffelem");

  // Loop over nelements blocks in the coefficient file

  for (int ielem = 0; ielem < nelements; ielem++) {
    for (int icoeff = 0; icoeff < nparams; icoeff++) {
      if (comm->me == 0) {
        ptr = fgets(line, MAXLINE, fpcoeff);
        if (ptr == nullptr) {
          eof = 1;
          fclose(fpcoeff);
        } else
          n = strlen(line) + 1;
      }

      MPI_Bcast(&eof, 1, MPI_INT, 0, world);
      if (eof) error->all(FLERR, "Incorrect format in MLIAPModel coefficient file");
      MPI_Bcast(&n, 1, MPI_INT, 0, world);
      MPI_Bcast(line, n, MPI_CHAR, 0, world);

      try {
        ValueTokenizer coeffs(utils::trim_comment(line));
        if (coeffs.count() != 1) throw TokenizerException("Wrong number of items", "");
        coeffelem[ielem][icoeff] = coeffs.next_double();
      } catch (TokenizerException &e) {
        error->all(FLERR,
                   "Incorrect format in MLIAPModel "
                   "coefficient file: {}",
                   e.what());
      }
    }
  }
  if (comm->me == 0) fclose(fpcoeff);
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double MLIAPModelSimple::memory_usage()
{
  double bytes = 0;

  bytes += (double) nelements * nparams * sizeof(double);    // coeffelem
  return bytes;
}
