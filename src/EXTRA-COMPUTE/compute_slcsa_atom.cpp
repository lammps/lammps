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
   Contributing author: Paul Lafourcade (CEA-DAM-DIF, Arpajon, France)
------------------------------------------------------------------------- */

#include "compute_slcsa_atom.h"

#include "arg_info.h"
#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "potential_file_reader.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;

static const char cite_compute_slcsa_atom_c[] =
    "compute slcsa/atom command: doi:10.1088/0965-0393/21/5/055020\n\n"
    "@Article{Lafourcade2023,\n"
    " author = {P. Lafourcade and J.-B. Maillet and C. Denoual and E. Duval and A. Allera and A. "
    "M. Goryaeva and M.-C. Marinica},\n"
    " title = {Robust crystal structure identification at extreme conditions using a "
    "density-independent spectral descriptor and supervised learning},\n"
    " journal = {Computational Materials Science},\n"
    " year =    2023,\n"
    " volume =  XX,\n"
    " pages =   {XXXXXX}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */
ComputeSLCSAAtom::ComputeSLCSAAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), list(nullptr), lda_scalings(nullptr),
    database_mean_descriptor(nullptr), lr_bias(nullptr), lr_decision(nullptr), icov_list(nullptr),
    mean_projected_descriptors(nullptr), maha_thresholds(nullptr), full_descriptor(nullptr),
    projected_descriptor(nullptr), scores(nullptr), probas(nullptr), prodright(nullptr),
    dmaha(nullptr), classification(nullptr)
{
  // command : compute c1 all slcsa/atom jmax nclasses parameters_file.dat
  // example : compute c1 all slcsa/atom 8 4 slcsa_parameters.dat
  // example : compute c1 all slcsa/atom 8 4 database_mean_descriptor.dat lda_scalings.dat lr_decision.dat lr_bias.dat mahalanobis_data.dat c_b1[*]
  // Steps :
  // 1. bs=bs-xbar
  // 2. dred=coefs_lda*bs
  // 3. scores=decision_lr*dred + lr_bias
  // 4. probas=exp(scores)/sum(exp(scores))
  // 5. cs=argmax(probas)

  // Read the parameters file in one bloc
  // File structure :
  // # database mean descriptor
  // vector with bso4dim rows x 1 col
  // # LDA dimension reduction matrix
  // matrix bso4dim rows x nclasses-1 cols
  // # LR decision matrix
  // matrix with nclasses rows x nclasses-1 cols
  // # LR bias vector
  // vector with 1 row x nclasses cols

  if (narg != 11) utils::missing_cmd_args(FLERR, "compute slcsa/atom", error);

  int twojmax = utils::inumeric(FLERR, arg[3], false, lmp);
  if (twojmax < 0)
    error->all(FLERR, "Illegal compute slcsa/atom command: twojmax must be a non-negative integer");
  ncomps = compute_ncomps(twojmax);

  nclasses = utils::inumeric(FLERR, arg[4], false, lmp);
  if (nclasses < 2)
    error->all(FLERR, "Illegal compute slcsa/atom command: nclasses must be greater than 1");

  database_mean_descriptor_file = arg[5];
  lda_scalings_file = arg[6];
  lr_decision_file = arg[7];
  lr_bias_file = arg[8];
  maha_file = arg[9];

  if (comm->me == 0) {
    auto mesg = fmt::format(
        "Files used:\n  {:24}: {}\n  {:24}: {}\n  {:24}: {}\n  {:24}: {}\n  {:24}: {}\n",
        "database mean descriptor", database_mean_descriptor_file, "lda scalings",
        lda_scalings_file, "lr decision", lr_decision_file, "lr bias", lr_bias_file, "maha stats",
        maha_file);
    utils::logmesg(lmp, mesg);
  }

  int expand = 0;
  char **earg;
  int nvalues = utils::expand_args(FLERR, narg - 10, &arg[10], 1, earg, lmp);
  if (earg != &arg[10]) expand = 1;
  arg = earg;

  ArgInfo argi(arg[0]);
  value_t val;
  val.id = "";
  val.val.c = nullptr;
  val.which = argi.get_type();
  val.argindex = argi.get_index1();
  val.id = argi.get_name();
  if ((val.which == ArgInfo::FIX) || (val.which == ArgInfo::VARIABLE) ||
      (val.which == ArgInfo::UNKNOWN) || (val.which == ArgInfo::NONE) || (argi.get_dim() > 1))
    error->all(FLERR, "Invalid compute slcsa/atom argument: {}", arg[0]);

  // if wildcard expansion occurred, free earg memory from exapnd_args()

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete[] earg[i];
    memory->sfree(earg);
  }

  val.val.c = modify->get_compute_by_id(val.id);
  if (!val.val.c) error->all(FLERR, "Compute ID {} for fix slcsa/atom does not exist", val.id);
  if (val.val.c->peratom_flag == 0)
    error->all(FLERR, "Compute slcsa/atom compute {} does not calculate per-atom values", val.id);
  if (val.argindex == 0 && val.val.c->size_peratom_cols != 0)
    error->all(FLERR, "Compute slcsa/atom compute {} does not calculate a per-atom vector", val.id);
  if (val.argindex && val.val.c->size_peratom_cols == 0)
    error->all(FLERR, "Compute slcsa/atom compute {} does not calculate a per-atom array", val.id);
  if (val.argindex && val.argindex > val.val.c->size_peratom_cols)
    error->all(FLERR, "Compute slcsa/atom compute {} array is accessed out-of-range", val.id);
  descriptorval = val;
  memory->create(database_mean_descriptor, ncomps, "slcsa/atom:database_mean_descriptor");
  memory->create(lda_scalings, ncomps, nclasses - 1, "slcsa/atom:lda_scalings");
  memory->create(lr_decision, nclasses, nclasses - 1, "slcsa/atom:lr_decision");
  memory->create(lr_bias, nclasses, "slcsa/atom:lr_bias");
  memory->create(maha_thresholds, nclasses, "slcsa/atom:maha_thresholds");
  memory->create(icov_list, nclasses, nclasses - 1, nclasses - 1, "slcsa/atom:icov_list");
  memory->create(mean_projected_descriptors, nclasses, nclasses - 1,
                 "slcsa/atom:mean_projected_descriptors");

  if (comm->me == 0) {

    if (strcmp(database_mean_descriptor_file, "NULL") == 0) {
      error->one(FLERR,
                 "Cannot open database mean descriptor file {}: ", database_mean_descriptor_file,
                 utils::getsyserror());
    } else {
      PotentialFileReader reader(lmp, database_mean_descriptor_file,
                                 "database mean descriptor file");
      int nread = 0;
      while (nread < ncomps) {
        auto values = reader.next_values(0);
        database_mean_descriptor[nread] = values.next_double();
        nread++;
      }
    }

    if (strcmp(lda_scalings_file, "NULL") == 0) {
      error->one(FLERR, "Cannot open database linear discriminant analysis scalings file {}: ",
                 lda_scalings_file, utils::getsyserror());
    } else {
      PotentialFileReader reader(lmp, lda_scalings_file, "lda scalings file");
      int nread = 0;
      while (nread < ncomps) {
        auto values = reader.next_values(nclasses - 1);
        lda_scalings[nread][0] = values.next_double();
        lda_scalings[nread][1] = values.next_double();
        lda_scalings[nread][2] = values.next_double();
        nread++;
      }
    }

    if (strcmp(lr_decision_file, "NULL") == 0) {
      error->one(FLERR, "Cannot open logistic regression decision file {}: ", lr_decision_file,
                 utils::getsyserror());
    } else {
      PotentialFileReader reader(lmp, lr_decision_file, "lr decision file");
      int nread = 0;
      while (nread < nclasses) {
        auto values = reader.next_values(nclasses - 1);
        lr_decision[nread][0] = values.next_double();
        lr_decision[nread][1] = values.next_double();
        lr_decision[nread][2] = values.next_double();
        nread++;
      }
    }

    if (strcmp(lr_bias_file, "NULL") == 0) {
      error->one(FLERR, "Cannot open logistic regression bias file {}: ", lr_bias_file,
                 utils::getsyserror());
    } else {
      PotentialFileReader reader(lmp, lr_bias_file, "lr bias file");
      auto values = reader.next_values(nclasses);
      lr_bias[0] = values.next_double();
      lr_bias[1] = values.next_double();
      lr_bias[2] = values.next_double();
      lr_bias[3] = values.next_double();
    }

    if (strcmp(maha_file, "NULL") == 0) {
      error->one(FLERR, "Cannot open mahalanobis stats file {}: ", maha_file, utils::getsyserror());
    } else {
      PotentialFileReader reader(lmp, maha_file, "mahalanobis stats file");
      int nvalues = nclasses * ((nclasses - 1) * (nclasses - 1) + nclasses);
      auto values = reader.next_values(nvalues);

      for (int i = 0; i < nclasses; i++) {
        maha_thresholds[i] = values.next_double();
        for (int j = 0; j < nclasses - 1; j++)
          mean_projected_descriptors[i][j] = values.next_double();
        for (int k = 0; k < nclasses - 1; k++)
          for (int l = 0; l < nclasses - 1; l++) icov_list[i][k][l] = values.next_double();
      }

      for (int i = 0; i < nclasses; i++) {
        auto mesg = fmt::format("For class {}  maha threshold = {:.6}\n", i, maha_thresholds[i]);
        mesg += "  mean B:\n";
        for (int j = 0; j < nclasses - 1; j++)
          mesg += fmt::format("   {:11.6}\n", mean_projected_descriptors[i][j]);
        mesg += "  icov:\n";
        for (int j = 0; j < nclasses - 1; j++) {
          mesg += fmt::format("   {:11.6} {:11.6} {:11.6}\n", icov_list[i][j][0],
                              icov_list[i][j][1], icov_list[i][j][2]);
        }
        utils::logmesg(lmp, mesg);
      }
    }
  }

  MPI_Bcast(&database_mean_descriptor[0], ncomps, MPI_DOUBLE, 0, world);
  MPI_Bcast(&lda_scalings[0][0], ncomps * (nclasses - 1), MPI_DOUBLE, 0, world);
  MPI_Bcast(&lr_decision[0][0], nclasses * (nclasses - 1), MPI_DOUBLE, 0, world);
  MPI_Bcast(&lr_bias[0], nclasses, MPI_DOUBLE, 0, world);
  MPI_Bcast(&maha_thresholds[0], nclasses, MPI_DOUBLE, 0, world);
  MPI_Bcast(&mean_projected_descriptors[0][0], nclasses * (nclasses - 1), MPI_DOUBLE, 0, world);
  MPI_Bcast(&icov_list[0][0][0], nclasses * (nclasses - 1) * (nclasses - 1), MPI_DOUBLE, 0, world);

  peratom_flag = 1;
  size_peratom_cols = nclasses + 1;
  ncols = nclasses + 1;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeSLCSAAtom::~ComputeSLCSAAtom()
{
  memory->destroy(classification);
  memory->destroy(database_mean_descriptor);
  memory->destroy(lda_scalings);
  memory->destroy(lr_decision);
  memory->destroy(lr_bias);
  memory->destroy(maha_thresholds);
  memory->destroy(mean_projected_descriptors);
  memory->destroy(icov_list);
  memory->destroy(full_descriptor);
  memory->destroy(projected_descriptor);
  memory->destroy(scores);
  memory->destroy(probas);
  memory->destroy(prodright);
  memory->destroy(dmaha);
}

/* ---------------------------------------------------------------------- */

void ComputeSLCSAAtom::init()
{

  if (modify->get_compute_by_style(style).size() > 1)
    if (comm->me == 0) error->warning(FLERR, "More than one compute {}", style);
}

/* ---------------------------------------------------------------------- */

void ComputeSLCSAAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSLCSAAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow per-atom if necessary

  if (atom->nmax > nmax) {
    memory->destroy(classification);
    nmax = atom->nmax;
    memory->create(classification, nmax, ncols, "slcsa/atom:classification");
    array_atom = classification;
  }

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (descriptorval.which == ArgInfo::COMPUTE) {
    if (!(descriptorval.val.c->invoked_flag & Compute::INVOKED_PERATOM)) {
      descriptorval.val.c->compute_peratom();
      descriptorval.val.c->invoked_flag |= Compute::INVOKED_PERATOM;
    }
    double **compute_array = descriptorval.val.c->array_atom;

    memory->create(full_descriptor, ncomps, "slcsa/atom:local descriptor");
    memory->create(projected_descriptor, nclasses - 1, "slcsa/atom:reduced descriptor");
    memory->create(scores, nclasses, "slcsa/atom:scores");
    memory->create(probas, nclasses, "slcsa/atom:probas");
    memory->create(prodright, nclasses - 1, "slcsa/atom:prodright");
    memory->create(dmaha, nclasses, "slcsa/atom:prodright");

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        for (int j = 0; j < ncomps; j++) full_descriptor[j] = compute_array[i][j];
        // Here comes the LDA + LR process
        // 1st step : Retrieve mean database descriptor
        for (int j = 0; j < ncomps; j++) full_descriptor[j] -= database_mean_descriptor[j];
        // 2nd step : Matrix multiplication to go from ncompsx1 -> (nclasses-1)*1
        for (int j = 0; j < nclasses - 1; j++) {
          projected_descriptor[j] = 0.0;
          for (int k = 0; k < ncomps; k++) {
            projected_descriptor[j] += full_descriptor[k] * lda_scalings[k][j];
          }
        }
        // 3rd step : Matrix multiplication
        for (int j = 0; j < nclasses; j++) {
          scores[j] = lr_bias[j];
          for (int k = 0; k < nclasses - 1; k++) {
            scores[j] += lr_decision[j][k] * projected_descriptor[k];
          }
        }
        // 4th step : Matrix multiplication
        double sumexpscores = 0.0;
        for (int j = 0; j < nclasses; j++) sumexpscores += exp(scores[j]);
        for (int j = 0; j < nclasses; j++) probas[j] = exp(scores[j]) / sumexpscores;

        classification[i][nclasses] = argmax(probas, nclasses);

        // 5th step : Mahalanobis distance
        for (int j = 0; j < nclasses; j++) {
          prodright[0] = 0.0;
          prodright[1] = 0.0;
          prodright[2] = 0.0;
          for (int k = 0; k < nclasses - 1; k++) {
            for (int l = 0; l < nclasses - 1; l++) {
              prodright[k] += (icov_list[j][k][l] *
                               (projected_descriptor[k] - mean_projected_descriptors[j][k]));
            }
          }
          double prodleft = 0.0;
          for (int k = 0; k < nclasses - 1; k++) {
            prodleft +=
                (prodright[k] * (projected_descriptor[k] - mean_projected_descriptors[j][k]));
          }
          classification[i][j] = sqrt(prodleft);
        }
        // 6th step : Sanity check
        int locclass = classification[i][nclasses];

        if (classification[i][locclass] > maha_thresholds[locclass]) {
          classification[i][nclasses] = -1.0;
        }

      } else {
        for (int j = 0; j < ncols; j++) classification[i][j] = -1.0;
      }
    }
    memory->destroy(full_descriptor);
    memory->destroy(projected_descriptor);
    memory->destroy(scores);
    memory->destroy(probas);
    memory->destroy(prodright);
    memory->destroy(dmaha);
  }
}

int ComputeSLCSAAtom::compute_ncomps(int twojmax)
{
  int ncount;

  ncount = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) ncount++;

  return ncount;
}

int ComputeSLCSAAtom::argmax(double arr[], int size)
{
  int maxIndex = 0;            // Initialize the index of the maximum value to the first element.
  double maxValue = arr[0];    // Initialize the maximum value to the first element.

  for (int i = 1; i < size; ++i) {
    if (arr[i] > maxValue) {
      // If a greater value is found, update the maxIndex and maxValue.
      maxIndex = i;
      maxValue = arr[i];
    }
  }

  return maxIndex;
}
