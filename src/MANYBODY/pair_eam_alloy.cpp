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
   Contributing authors: Stephen Foiles (SNL), Murray Daw (SNL)
------------------------------------------------------------------------- */

#include "pair_eam_alloy.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "potential_file_reader.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairEAMAlloy::PairEAMAlloy(LAMMPS *lmp) : PairEAM(lmp)
{
  one_coeff = 1;
  manybody_flag = 1;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read DYNAMO setfl file
------------------------------------------------------------------------- */

void PairEAMAlloy::coeff(int narg, char **arg)
{
  int i, j;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR, "Number of element to type mappings does not match number of atom types");

  // read EAM setfl file

  if (setfl) {
    for (i = 0; i < setfl->nelements; i++) delete[] setfl->elements[i];
    delete[] setfl->elements;
    memory->destroy(setfl->mass);
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->z2r);
    delete setfl;
  }
  setfl = new Setfl();
  read_file(arg[2]);

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if "NULL"

  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i], "NULL") == 0) {
      map[i - 2] = -1;
      continue;
    }
    for (j = 0; j < setfl->nelements; j++)
      if (strcmp(arg[i], setfl->elements[j]) == 0) break;
    if (j < setfl->nelements)
      map[i - 2] = j;
    else
      error->all(FLERR, "No matching element in EAM potential file");
  }

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++) setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass of atom type if i = j

  int count = 0;
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(FLERR, i, setfl->mass[map[i]]);
        count++;
      }
      scale[i][j] = 1.0;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   read a multi-element DYNAMO setfl file
------------------------------------------------------------------------- */

void PairEAMAlloy::read_file(char *filename)
{
  Setfl *file = setfl;

  // read potential file
  if (comm->me == 0) {
    PotentialFileReader reader(lmp, filename, "eam/alloy", unit_convert_flag);

    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY, unit_convert);
    try {
      reader.skip_line();
      reader.skip_line();
      reader.skip_line();

      // extract element names from nelements line
      ValueTokenizer values = reader.next_values(1);
      file->nelements = values.next_int();

      if ((int) values.count() != file->nelements + 1)
        error->one(FLERR, "Incorrect element names in EAM potential file");

      file->elements = new char *[file->nelements];
      for (int i = 0; i < file->nelements; i++)
        file->elements[i] = utils::strdup(values.next_string());

      values = reader.next_values(5);
      file->nrho = values.next_int();
      file->drho = values.next_double();
      file->nr = values.next_int();
      file->dr = values.next_double();
      file->cut = values.next_double();

      if ((file->nrho <= 0) || (file->nr <= 0) || (file->dr <= 0.0))
        error->one(FLERR, "Invalid EAM potential file");

      memory->create(file->mass, file->nelements, "pair:mass");
      memory->create(file->frho, file->nelements, file->nrho + 1, "pair:frho");
      memory->create(file->rhor, file->nelements, file->nr + 1, "pair:rhor");
      memory->create(file->z2r, file->nelements, file->nelements, file->nr + 1, "pair:z2r");

      for (int i = 0; i < file->nelements; i++) {
        values = reader.next_values(2);
        values.next_int();    // ignore
        file->mass[i] = values.next_double();

        reader.next_dvector(&file->frho[i][1], file->nrho);
        reader.next_dvector(&file->rhor[i][1], file->nr);
        if (unit_convert) {
          for (int j = 1; j < file->nrho; ++j) file->frho[i][j] *= conversion_factor;
        }
      }

      for (int i = 0; i < file->nelements; i++) {
        for (int j = 0; j <= i; j++) {
          reader.next_dvector(&file->z2r[i][j][1], file->nr);
          if (unit_convert) {
            for (int k = 1; k < file->nr; ++k) file->z2r[i][j][k] *= conversion_factor;
          }
        }
      }
    } catch (TokenizerException &e) {
      error->one(FLERR, e.what());
    }
  }

  // broadcast potential information
  MPI_Bcast(&file->nelements, 1, MPI_INT, 0, world);

  MPI_Bcast(&file->nrho, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->drho, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->nr, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->dr, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->cut, 1, MPI_DOUBLE, 0, world);

  // allocate memory on other procs
  if (comm->me != 0) {
    file->elements = new char *[file->nelements];
    for (int i = 0; i < file->nelements; i++) file->elements[i] = nullptr;
    memory->create(file->mass, file->nelements, "pair:mass");
    memory->create(file->frho, file->nelements, file->nrho + 1, "pair:frho");
    memory->create(file->rhor, file->nelements, file->nr + 1, "pair:rhor");
    memory->create(file->z2r, file->nelements, file->nelements, file->nr + 1, "pair:z2r");
  }

  // broadcast file->elements string array
  for (int i = 0; i < file->nelements; i++) {
    int n;
    if (comm->me == 0) n = strlen(file->elements[i]) + 1;
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    if (comm->me != 0) file->elements[i] = new char[n];
    MPI_Bcast(file->elements[i], n, MPI_CHAR, 0, world);
  }

  // broadcast file->mass, frho, rhor
  for (int i = 0; i < file->nelements; i++) {
    MPI_Bcast(&file->mass[i], 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&file->frho[i][1], file->nrho, MPI_DOUBLE, 0, world);
    MPI_Bcast(&file->rhor[i][1], file->nr, MPI_DOUBLE, 0, world);
  }

  // broadcast file->z2r
  for (int i = 0; i < file->nelements; i++) {
    for (int j = 0; j <= i; j++) { MPI_Bcast(&file->z2r[i][j][1], file->nr, MPI_DOUBLE, 0, world); }
  }
}

/* ----------------------------------------------------------------------
   copy read-in setfl potential to standard array format
------------------------------------------------------------------------- */

void PairEAMAlloy::file2array()
{
  int i, j, m, n;
  int ntypes = atom->ntypes;

  // set function params directly from setfl file

  nrho = setfl->nrho;
  nr = setfl->nr;
  drho = setfl->drho;
  dr = setfl->dr;
  rhomax = (nrho - 1) * drho;

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of setfl elements + 1 for zero array

  nfrho = setfl->nelements + 1;
  memory->destroy(frho);
  memory->create(frho, nfrho, nrho + 1, "pair:frho");

  // copy each element's frho to global frho

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nrho; m++) frho[i][m] = setfl->frho[i][m];

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho - 1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to element (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0)
      type2frho[i] = map[i];
    else
      type2frho[i] = nfrho - 1;

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of setfl elements

  nrhor = setfl->nelements;
  memory->destroy(rhor);
  memory->create(rhor, nrhor, nr + 1, "pair:rhor");

  // copy each element's rhor to global rhor

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nr; m++) rhor[i][m] = setfl->rhor[i][m];

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for setfl files, I,J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++) type2rhor[i][j] = map[i];

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of setfl elements

  nz2r = setfl->nelements * (setfl->nelements + 1) / 2;
  memory->destroy(z2r);
  memory->create(z2r, nz2r, nr + 1, "pair:z2r");

  // copy each element pair z2r to global z2r, only for I >= J

  n = 0;
  for (i = 0; i < setfl->nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= nr; m++) z2r[n][m] = setfl->z2r[i][j][m];
      n++;
    }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow, icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}
