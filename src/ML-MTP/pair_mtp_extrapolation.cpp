/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

//
// Contributing author, Richard Meng, Queen's University at Kingston, 10.02.25, contact@richardzjm.com
//

#include "pair_mtp_extrapolation.h"
#include "mtp_radial_basis.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <fstream>

using namespace LAMMPS_NS;

PairMTPExtrapolation::PairMTPExtrapolation(LAMMPS *lmp) : PairMTP(lmp)
{
  nextra = 1;                      // Number of extra coefficients (1 for extrapolation)
  pvector = new double[nextra];    // Pointer directly to the max extrapolation grade
  pvector[0] = 0.0;
};

/* ---------------------------------------------------------------------- */

PairMTPExtrapolation::~PairMTPExtrapolation()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(active_set);
    memory->destroy(inverse_active_set);
    memory->destroy(radial_jacobian);
    memory->destroy(energy_ders_wrt_coeffs);
    if (!configuration_mode) memory->destroy(nbh_extrapolation_grades);
    if (mlip3_style) memory->destroy(write_buffer);
  }
}

/* ----------------------------------------------------------------------
   Straightfoward MTP implementation based on MLIP3
   ---------------------------------------------------------------------- */

void PairMTPExtrapolation::compute(int eflag, int vflag)
{
  // If we are not extrapolating per fix pair and not extrapolating continously, we can just call the base class compute
  if (!extrapolation_flag && !mlip3_style) {
    PairMTP::compute(eflag, vflag);
    return;
  }

  max_grade = 0;

  ev_setup(eflag, vflag);

  double **x = atom->x;      // atomic positons
  double **f = atom->f;      // atomic forces
  int *type = atom->type;    //atomic types

  int inum = list->inum;             // The number of central atoms (neigbhourhoods)
  int *ilist = list->ilist;          // List of central atom ids
  int *numneigh = list->numneigh;    // List of the number of neighbours for each central atom
  int **firstneigh =
      list->firstneigh;    //List  (head of array) of neighbours for a given central atom

  // Resize the nbh extrapolation grades if needed. No need to initialize, first access is write
  if (!configuration_mode && nbh_count < inum) {
    memory->grow(nbh_extrapolation_grades, inum, "nbh_extrapolation_grades");
    nbh_count = inum;
  }

  // If are in configuration, we need to reset the working array once per compute call / config
  if (configuration_mode)
    std::fill(&energy_ders_wrt_coeffs[0], &energy_ders_wrt_coeffs[0] + coeff_count, 0.0);

  // Loop over all provided neighbourhoods
  for (int ii = 0; ii < inum; ii++) {
    int valid_count = 0;
    const int i = ilist[ii];    // Set central atom index
    const int itype = map[type[i]];
    int jnum = numneigh[i];    // Set number of neighbours
    double nbh_energy = 0;
    const double xi[3] = {x[i][0], x[i][1],
                          x[i][2]};    // Cache the position of the central atom for efficiency

    // Resize the jacobian and cutoff if needed. No need to initialize, first access is write
    if (jac_size < jnum) {
      memory->grow(moment_jacobian, jnum, alpha_index_basic_count, 3, "moment_jacobian");
      memory->grow(valid_j, jnum, "valid_j");
      jac_size = jnum;
    }

    // Reset the working arrays
    std::fill(&moment_tensor_vals[0], &moment_tensor_vals[0] + alpha_moment_count, 0.0);
    std::fill(&nbh_energy_ders_wrt_moments[0], &nbh_energy_ders_wrt_moments[0] + alpha_moment_count,
              0.0);
    std::fill(&radial_moment_ders[0], &radial_moment_ders[0] + alpha_moment_count, 0.0);
    std::fill(&radial_jacobian[0][0][0],
              &radial_jacobian[0][0][0] +
                  (alpha_index_basic_count * species_count * radial_coeff_count_per_pair),
              0.0);

    if (!configuration_mode)
      std::fill(&energy_ders_wrt_coeffs[0], &energy_ders_wrt_coeffs[0] + coeff_count, 0.0);

    // ------------ Begin Alpha Basic Calc ------------
    // Loop over all neighbours
    for (int jj = 0; jj < jnum; jj++) {
      int j = firstneigh[i][jj];    //List of neighbours
      j &= NEIGHMASK;
      const int jtype = map[type[j]];
      const double r[3] = {x[j][0] - xi[0], x[j][1] - xi[1], x[j][2] - xi[2]};
      const double rsq = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];

      if (rsq > max_cutoff_sq) continue;

      valid_j[valid_count] = j;

      const double dist = std::sqrt(rsq);
      radial_basis->calc_radial_basis_ders(dist);    // Calculate radial basis

      // Precompute the coord and distance power
      for (int k = 1; k < max_alpha_index_basic; k++) {
        dist_powers[k] = dist_powers[k - 1] * dist;
        for (int a = 0; a < 3; a++) coord_powers[k][a] = coord_powers[k - 1][a] * r[a];
      }

      // Compute the radial basis values
      int pair_offset = itype * species_count + jtype;
      for (int mu = 0; mu < radial_func_count; mu++) {
        double val = 0;
        double der = 0;
        int offset = (pair_offset * radial_coeff_count_per_pair) + mu * radial_basis_size;

        for (int ri = 0; ri < radial_basis_size; ri++) {
          val += radial_basis_coeffs[offset + ri] * radial_basis->radial_basis_vals[ri];
          der += radial_basis_coeffs[offset + ri] * radial_basis->radial_basis_ders[ri];
        }
        radial_vals[mu] = val;
        radial_ders[mu] = der;
      }

      //Calculate the alpha basics
      for (int k = 0; k < alpha_index_basic_count; k++) {
        double val = 0;
        double der = 0;
        int mu = alpha_index_basic[k][0];

        val = radial_vals[mu];
        der = radial_ders[mu];

        // Normalize by the rank of alpha's coresponding tensor
        int norm_rank = alpha_index_basic[k][1] + alpha_index_basic[k][2] + alpha_index_basic[k][3];
        double norm_fac = 1.0 / dist_powers[norm_rank];

        double pow0 = coord_powers[alpha_index_basic[k][1]][0];
        double pow1 = coord_powers[alpha_index_basic[k][2]][1];
        double pow2 = coord_powers[alpha_index_basic[k][3]][2];
        double pow = pow0 * pow1 * pow2;

        // Calculate the radial jacobian
        int mu_offset = mu * radial_basis_size;
        for (int ri = 0; ri < radial_basis_size; ri++) {
          radial_jacobian[k][jtype][mu_offset + ri] +=
              radial_basis->radial_basis_vals[ri] * norm_fac * pow;
        }

        val *= norm_fac;
        der = der * norm_fac - norm_rank * val / dist;
        moment_tensor_vals[k] += val * pow;

        // Get the component's derivatives too
        pow *= der / dist;
        moment_jacobian[valid_count][k][0] = pow * r[0];
        moment_jacobian[valid_count][k][1] = pow * r[1];
        moment_jacobian[valid_count][k][2] = pow * r[2];

        if (alpha_index_basic[k][1] != 0) {
          moment_jacobian[valid_count][k][0] += val * alpha_index_basic[k][1] *
              coord_powers[alpha_index_basic[k][1] - 1][0] * pow1 * pow2;
        }    //Chain rule for nonzero rank
        if (alpha_index_basic[k][2] != 0) {
          moment_jacobian[valid_count][k][1] += val * alpha_index_basic[k][2] * pow0 *
              coord_powers[alpha_index_basic[k][2] - 1][1] * pow2;
        }    //Chain rule for nonzero rank
        if (alpha_index_basic[k][3] != 0) {
          moment_jacobian[valid_count][k][2] += val * alpha_index_basic[k][3] * pow0 * pow1 *
              coord_powers[alpha_index_basic[k][3] - 1][2];
        }    //Chain rule for nonzero rank
      }
      valid_count++;
    }

    // ------------ Contruct Other Alphas  ------------
    for (int k = 0; k < alpha_index_times_count; k++) {
      double val0 = moment_tensor_vals[alpha_index_times[k][0]];
      double val1 = moment_tensor_vals[alpha_index_times[k][1]];
      int val2 = alpha_index_times[k][2];
      moment_tensor_vals[alpha_index_times[k][3]] += val2 * val0 * val1;
    }

    // ------------ Convolve Basis Set From Alpha Map ------------
    // Calculate the energies only if needed. We always need the basis member for extrapolation.
    int linear_basis_offset = radial_coeff_count + species_count;
    if (eflag_either) {
      nbh_energy = species_coeffs[itype];    // Essentially the reference point energy per species
      for (int k = 0; k < alpha_scalar_count; k++) {
        double basis_member = moment_tensor_vals[alpha_moment_mapping[k]];
        energy_ders_wrt_coeffs[linear_basis_offset + k] += basis_member;
        nbh_energy += linear_coeffs[k] * basis_member;
      }
      // Tally energies per flags if needed
      if (eflag_atom) eatom[i] = nbh_energy;
      if (eflag_global) eng_vdwl += nbh_energy;
    } else
      for (int k = 0; k < alpha_scalar_count; k++)
        energy_ders_wrt_coeffs[linear_basis_offset + k] +=
            moment_tensor_vals[alpha_moment_mapping[k]];

    // ------------ Also add the species coefficient ------------
    energy_ders_wrt_coeffs[radial_coeff_count + itype] += 1;

    // =========== Begin Backpropogation ===========

    //------------ Step 1: NBH energy derivative is the corresponding linear combination------------
    for (int k = 0; k < alpha_scalar_count; k++)
      nbh_energy_ders_wrt_moments[alpha_moment_mapping[k]] = linear_coeffs[k];

    //------------ Step 2: Propogate chain rule through the alpha times to the alpha basics ------------
    for (int k = alpha_index_times_count - 1; k >= 0; k--) {
      int a0 = alpha_index_times[k][0];
      int a1 = alpha_index_times[k][1];
      int multipiler = alpha_index_times[k][2];
      int a3 = alpha_index_times[k][3];

      double val0 = moment_tensor_vals[a0];
      double val1 = moment_tensor_vals[a1];
      double val3 = nbh_energy_ders_wrt_moments[a3];

      nbh_energy_ders_wrt_moments[a1] += val3 * multipiler * val0;
      nbh_energy_ders_wrt_moments[a0] += val3 * multipiler * val1;
    }

    //------------ Step 3: Multiply energy ders wrt moment by the moment jacobian to get forces ------------
    for (int jj = 0; jj < valid_count; jj++) {
      int j = valid_j[jj];

      double temp_force[3] = {0, 0, 0};
      for (int k = 0; k < alpha_index_basic_count; k++) {
        // Backprop relative to positions for forces
        for (int a = 0; a < 3; a++) {
          //Calculate forces
          temp_force[a] += nbh_energy_ders_wrt_moments[k] * moment_jacobian[jj][k][a];
        }
      }

      f[i][0] += temp_force[0];
      f[i][1] += temp_force[1];
      f[i][2] += temp_force[2];

      f[j][0] -= temp_force[0];
      f[j][1] -= temp_force[1];
      f[j][2] -= temp_force[2];

      //Calculate virial stress
      if (vflag) {
        // We only need to calculate rel pos again if stress are needed
        const double r[3] = {x[j][0] - xi[0], x[j][1] - xi[1], x[j][2] - xi[2]};
        virial[0] -= temp_force[0] * r[0];    //xx
        virial[1] -= temp_force[1] * r[1];    //yy
        virial[2] -= temp_force[2] * r[2];    //zz

        virial[3] -= (temp_force[0] * r[1] + temp_force[1] * r[0]) / 2;    //xy
        virial[4] -= (temp_force[0] * r[2] + temp_force[2] * r[0]) / 2;    //xz
        virial[5] -= (temp_force[1] * r[2] + temp_force[2] * r[1]) / 2;    //yz

        if (vflag_atom) {
          vatom[i][0] -= temp_force[0] * r[0];    //xx
          vatom[i][1] -= temp_force[1] * r[1];    //yy
          vatom[i][2] -= temp_force[2] * r[2];    //zz

          vatom[i][3] -= (temp_force[0] * r[1] + temp_force[1] * r[0]) / 2;    //xy
          vatom[i][4] -= (temp_force[0] * r[2] + temp_force[2] * r[0]) / 2;    //xz
          vatom[i][5] -= (temp_force[1] * r[2] + temp_force[2] * r[1]) / 2;    //yz
        }
      }
    }

    //------------ Step 3.5: Multiply energy ders wrt moment by the radial jacobian to get rad ders ------------
    for (int k = 0; k < alpha_index_basic_count; k++)
      for (int jjtype = 0; jjtype < species_count; jjtype++) {
        int offset = (itype * species_count + jjtype) * radial_coeff_count_per_pair;
        for (int ri = 0; ri < radial_coeff_count_per_pair; ri++)
          energy_ders_wrt_coeffs[offset + ri] +=
              nbh_energy_ders_wrt_moments[k] * radial_jacobian[k][jjtype][ri];
      }

    // Directly calculate extraplation grade for neighbourhood mode
    if (!configuration_mode) {
      double grade = calculate_extrapolation_grade();
      max_grade = std::max(grade, max_grade);
      nbh_extrapolation_grades[i] = grade;
    }
  }

  compile_grades();

  if (mlip3_style) evaluate_grades();    // Evaluate grades per MLIP-3 two-threshold style
}

/* ----------------------------------------------------------------------
   Extrapolation Calculation Function
------------------------------------------------------------------------- */
double PairMTPExtrapolation::calculate_extrapolation_grade()
{
  double max_grade = 0;
  for (int i = 0; i < coeff_count; i++) {
    double current_grade = 0;
    for (int j = 0; j < coeff_count; j++) {
      current_grade += energy_ders_wrt_coeffs[j] * inverse_active_set[i][j];
    }
    max_grade = std::max(std::abs(current_grade), max_grade);
  }
  return max_grade;
}

/* ----------------------------------------------------------------------
   Collective Reduction Operation 
------------------------------------------------------------------------- */
void PairMTPExtrapolation::compile_grades()
{
  // MPI reduce operations based on selection mode
  if (configuration_mode) {    // Configuration mode

    // Sum contributions across all processes
    MPI_Allreduce(MPI_IN_PLACE, energy_ders_wrt_coeffs, coeff_count, MPI_DOUBLE, MPI_SUM, world);

    max_grade = calculate_extrapolation_grade();

    if (atom->natoms > 0)
      max_grade /= atom->natoms;    // Normalize
    else
      max_grade = 0.0;

  } else {    // Neighbourhood mode
    MPI_Allreduce(MPI_IN_PLACE, &max_grade, 1, MPI_DOUBLE, MPI_MAX, world);
  }
  if (comm->me == 0) pvector[0] = max_grade;    // Expose the max grade
}

/* ----------------------------------------------------------------------
   Evaluate Thresholds
------------------------------------------------------------------------- */
void PairMTPExtrapolation::evaluate_grades()
{
  if (max_grade >= select_threshold) write_config();
  if (max_grade >= break_threshold && comm->me == 0) {
    std::fflush(preselected_file);    // Ensure the writing buffers are flushed before breaking.
    std::fclose(preselected_file);
    error->one(FLERR, "Exceeded Break Threshold: {:.5f}. Terminating simulation.\n", max_grade);
  }
}
/* ----------------------------------------------------------------------
   Write current config to file
------------------------------------------------------------------------- */
void PairMTPExtrapolation::write_config()
{
  int inum = list->inum;
  int *type = atom->type;
  double **x = atom->x;
  int index_offset = 0;

  MPI_Scan(&inum, &index_offset, 1, MPI_INT, MPI_SUM, world);
  index_offset -= inum;

  const int MAX_LINE = 128;
  bigint needed_size = (bigint) inum * MAX_LINE;

  if (needed_size > write_buffer_size) {
    write_buffer_size = needed_size;
    memory->grow(write_buffer, write_buffer_size, "write_buffer");
  }

  bigint local_buffer_size = 0;
  for (int ii = 0; ii < inum; ii++) {
    const int i = ii;
    const int itype = map[type[i]];
    const double xi[3] = {x[i][0], x[i][1], x[i][2]};
    const int global_i = i + index_offset + 1;

    int n = 0;
    // snprintf returns the EXACT number of characters it wants to write
    if (!configuration_mode) {
      n = snprintf(write_buffer + local_buffer_size, MAX_LINE, "%d\t%d\t%.6f\t%.6f\t%.6f\t%.5f\n",
                   global_i, itype, xi[0], xi[1], xi[2], nbh_extrapolation_grades[i]);
    } else {
      n = snprintf(write_buffer + local_buffer_size, MAX_LINE, "%d\t%d\t%.6f\t%.6f\t%.6f\n",
                   global_i, itype, xi[0], xi[1], xi[2]);
    }

    // Memory usage check
    if (n >= MAX_LINE) n = MAX_LINE - 1;
    local_buffer_size += n;
  }

  // Find max memory needed and ensure Rank 0 can handle it
  bigint max_bytes_any_proc;
  MPI_Reduce(&local_buffer_size, &max_bytes_any_proc, 1, MPI_LMP_BIGINT, MPI_MAX, 0, world);
  if (comm->me == 0) {
    if (max_bytes_any_proc > write_buffer_size) {
      write_buffer_size = max_bytes_any_proc;
      memory->grow(write_buffer, write_buffer_size, "write_buffer");
    }

    // Standard MLIP header prints
    std::fprintf(preselected_file, "BEGIN_CFG\n");
    std::fprintf(preselected_file, "Size\n%ld\n", (long) atom->natoms);
    std::fprintf(preselected_file, "Supercell\n");
    std::fprintf(preselected_file, "%.6f %.6f %.6f\n", domain->xprd, 0.0, 0.0);
    std::fprintf(preselected_file, "%.6f %.6f %.6f\n", domain->xy, domain->yprd, 0.0);
    std::fprintf(preselected_file, "%.6f %.6f %.6f\n", domain->xz, domain->yz, domain->zprd);

    const char *header = (!configuration_mode)
        ? "AtomData:  id type       cartes_x      cartes_y      cartes_z       nbh_grades\n"
        : "AtomData:  id type       cartes_x      cartes_y      cartes_z\n";
    std::fputs(header, preselected_file);

    // Write Rank 0's data
    std::fwrite(write_buffer, 1, local_buffer_size, preselected_file);

    // Receive and write from others
    for (int i = 1; i < comm->nprocs; i++) {
      MPI_Status status;
      MPI_Recv(write_buffer, write_buffer_size, MPI_CHAR, i, 0, world, &status);
      int n_received;
      MPI_Get_count(&status, MPI_CHAR, &n_received);
      std::fwrite(write_buffer, 1, n_received, preselected_file);
    }

    std::fprintf(preselected_file, "Feature   MV_grade\t%.6f\n", max_grade);
    std::fprintf(preselected_file, "END_CFG\n\n");
  } else {
    MPI_Send(write_buffer, local_buffer_size, MPI_CHAR, 0, 0, world);
  }
}
/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMTPExtrapolation::settings(int narg, char **arg)
{
  if ((narg == 2 && LAMMPS_NS::utils::lowercase(arg[0]) == "chunksize") ||
      (narg == 5 && LAMMPS_NS::utils::lowercase(arg[3]) == "chunksize")) {
    if (comm->me == 0) utils::logmesg(lmp, "Ignoring chunksize settings!\n");
    narg -= 2;    // Ignore the chunksize settings
  } else if (narg != 0 && narg != 3)
    error->all(FLERR,
               "Pair mtp/extrapolation accepts no arguments."
               "Or exactly 4 arguments: {output_file}. {selection_threshold} "
               "{break_threshold}.");

  if (narg == 3) {
    mlip3_style = true;
    select_threshold = utils::numeric(FLERR, arg[1], true, lmp);
    break_threshold = utils::numeric(FLERR, arg[2], true, lmp);
  }

  if (comm->me == 0)
    if (mlip3_style)
      utils::logmesg(lmp,
                     "Extrapolation Scheme: {} mode, with a selection threshold of {} "
                     "and break threshold of {}.\n",
                     (configuration_mode ? "Configuration" : "Neighborhood"), select_threshold,
                     break_threshold);
    else
      utils::logmesg(lmp, "Extrapolation Mode: {} mode.\n",
                     (configuration_mode ? "Configuration" : "Neighborhood"));

  if (mlip3_style)
    if (comm->me == 0) preselected_file = std::fopen(arg[0], "w");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */
void PairMTPExtrapolation::coeff(int narg, char **arg)
{
  const int n = atom->ntypes;
  if (narg != 3 + n)
    error->all(FLERR,
               "Incorrect args for pair coefficients. Pair_coeff mtp/extrapolation only accepts "
               "the MTP potential file and the species mappping.");

  // Read in MTP and allocate memory
  FILE *mtp_file = utils::open_potential(arg[2], lmp, nullptr);
  PairMTPExtrapolation::read_file(mtp_file);
  fclose(mtp_file);

  PairMTP::prepare_map(narg - 3, arg + 3);
}

/* ----------------------------------------------------------------------
   MTP file parsing helper function. Includes memory allocation. Excludes some radial basis hyperparameters (in radial basis constructor instead).
------------------------------------------------------------------------- */
void PairMTPExtrapolation::read_file(FILE *mtp_file)
{
  PairMTP::read_file(mtp_file);

  // Some size calcs
  coeff_count = radial_coeff_count + species_count + alpha_scalar_count;
  int num_doubles = coeff_count * coeff_count;

  // Now we allocate memory for the additional memory needed for calculations
  memory->create(active_set, coeff_count, coeff_count, "active_set");
  memory->create(inverse_active_set, coeff_count, coeff_count, "inverse_active_set");
  memory->create(radial_moment_ders, alpha_moment_count, "radial_moment_ders");
  memory->create(energy_ders_wrt_coeffs, coeff_count, "energy_ders_wrt_coeffs");
  memory->create(radial_jacobian, alpha_index_basic_count, species_count,
                 radial_coeff_count_per_pair, "radial_jacobian");
  // We initialize the extrapolation grades during compute since its size depend on problem size.

  if (comm->me == 0) {
    const std::string new_separators = "=, ";
    const std::string separators = TOKENIZER_DEFAULT_SEPARATORS + new_separators;
    TextFileReader tfr(mtp_file, "ml-mtp");
    tfr.ignore_comments = false;

    char *line = tfr.next_line();
    if (line == nullptr) {
      error->one(
          FLERR,
          "No selection state found! Consider training/retraining or disabling extrapolation!\n");
    }

    ValueTokenizer line_tokens = ValueTokenizer(std::string(line), separators);
    std::string keyword = line_tokens.next_string();
    if (keyword != "#MVS_v1.1")
      lmp->error->one(
          FLERR,
          "Error in reading MTP file selection state. Please verify MVS version is #MVS_v1.1!");
    tfr.ignore_comments = true;    // Accept comments after reading the version which is a comment

    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
    if (keyword != "energy_weight")
      lmp->error->one(FLERR, "Error in reading MTP file, energy_weight");
    const double energy_weight = line_tokens.next_double();

    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
    if (keyword != "force_weight")
      lmp->error->one(FLERR, "Error in reading MTP file, force_weight");

    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
    if (keyword != "stress_weight")
      lmp->error->one(FLERR, "Error in reading MTP file, stress_weight");

    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
    if (keyword != "site_en_weight")
      lmp->error->one(FLERR, "Error in reading MTP file, site_en_weight");
    const double site_en_weight = line_tokens.next_double();

    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
    if (keyword != "weight_scaling")
      lmp->error->one(FLERR, "Error in reading MTP file, weight_scaling");

    if (energy_weight + site_en_weight > 1)
      lmp->error->one(FLERR,
                      "Error, the MTP currently only supports configuration mode (energy_weight=1) "
                      "or neighbourhood mode (site_en_weight=1). "
                      "Please retrain the MTP with the correct modes!");

    configuration_mode = (energy_weight == 1);

    fgetc(mtp_file);    // We need to skip foward 1 character. There is a # before the binary data.
    utils::sfread(FLERR, &active_set[0][0], sizeof(double), num_doubles, mtp_file, nullptr,
                  lmp->error);
    utils::sfread(FLERR, &inverse_active_set[0][0], sizeof(double), num_doubles, mtp_file, nullptr,
                  lmp->error);
  }

  //Broadcast active set to others
  MPI_Bcast(&configuration_mode, 1, MPI_INT, 0, world);
  MPI_Bcast(&active_set[0][0], num_doubles, MPI_DOUBLE, 0, world);
  MPI_Bcast(&inverse_active_set[0][0], num_doubles, MPI_DOUBLE, 0, world);
  allocated = 1;
}

/* ----------------------------------------------------------------------
  Flag to indicate if we are computing extrapolation grades on this iteration
 ---------------------------------------------------------------------- */
void *PairMTPExtrapolation::extract(const char *str, int &dim)
{
  dim = 0;
  //check if str=="gamma_flag" then compute extrapolation grades on this iteration
  if (strcmp(str, "extrapolation_flag") == 0) return (void *) &extrapolation_flag;

  return nullptr;
}

/* ----------------------------------------------------------------------
   peratom requests from FixPair
   return ptr to requested data
   also return ncol = # of quantites per atom
     0 = per-atom vector
     1 or more = # of columns in per-atom array
   return NULL if str is not recognized
---------------------------------------------------------------------- */
void *PairMTPExtrapolation::extract_peratom(const char *str, int &ncol)
{
  if (strcmp(str, "extrapolation") == 0) {
    if (configuration_mode)
      error->one(FLERR, "Please use the MLIP-3 style extrapolation for configuration mode MTPs!");

    ncol = 0;
    return (void *) nbh_extrapolation_grades;
  }

  return nullptr;
}
