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

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "mtp_radial_basis.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "text_file_reader.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

using namespace LAMMPS_NS;

PairMTPExtrapolation::PairMTPExtrapolation(LAMMPS *lmp) : PairMTP(lmp)
{
  nextra = 1;                      // Number of extra coefficients (1 for extrapolation)
  pvector = new double[nextra];    // Pointer directly to the max extrapolation grade
  pvector[0] = 0.0;

  active_set = nullptr;
  inverse_active_set = nullptr;
  radial_basis_cache = nullptr;
  energy_ders_wrt_coeffs = nullptr;
  nbh_extrapolation_grades = nullptr;
  write_buffer = nullptr;
  preselected_file = nullptr;

  coeff_count = 0;
  radial_basis_cache_size = 0;
  extrapolation_flag = 0;
  mlip3_style = false;
  configuration_mode = 0;
  weight_scaling = 2;
  select_threshold = 0.0;
  break_threshold = 0.0;
  max_grade = 0.0;
  nbh_count = 0;
  write_buffer_size = 0;
};

/* ---------------------------------------------------------------------- */

PairMTPExtrapolation::~PairMTPExtrapolation()
{
  if (copymode) return;

  delete[] pvector;

  if (allocated) {
    memory->destroy(active_set);
    memory->destroy(inverse_active_set);
    memory->destroy(radial_basis_cache);
    memory->destroy(energy_ders_wrt_coeffs);
    memory->destroy(nbh_extrapolation_grades);
    memory->destroy(write_buffer);
  }

  // settings() may have opened this before coeff() ran, so close it unconditionally.
  if (comm->me == 0 && preselected_file) {
    std::fclose(preselected_file);
    preselected_file = nullptr;
  }
}

/* ----------------------------------------------------------------------
   Straightforward MTP implementation based on MLIP3
   ---------------------------------------------------------------------- */
void PairMTPExtrapolation::compute(int eflag, int vflag)
{
  // If we are not extrapolating per fix pair and not extrapolating continously, we can just call the base class compute
  if (!extrapolation_flag && !mlip3_style) {
    PairMTP::compute(eflag, vflag);
    return;
  }

  max_grade = 0;

  ev_init(eflag, vflag);

  // The candidate vector needs every scalar moment, so the pruned force-only
  // contraction list the base style uses does not apply here.
  const int *times_data = alpha_index_times_count ? alpha_index_times[0] : nullptr;

  double **x = atom->x;      // atomic positions
  double **f = atom->f;      // atomic forces
  int *type = atom->type;    //atomic types

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int inum = list->inum;             // The number of central atoms (neighbourhoods)
  int *ilist = list->ilist;          // List of central atom ids
  int *numneigh = list->numneigh;    // List of the number of neighbours for each central atom
  int **firstneigh =
      list->firstneigh;    //List  (head of array) of neighbours for a given central atom

  // Resize the nbh extrapolation grades if needed.
  if (!configuration_mode) {
    if (nbh_count < atom->nmax) {
      memory->grow(nbh_extrapolation_grades, atom->nmax, "nbh_extrapolation_grades");
      nbh_count = atom->nmax;
    }
    if (inum < nlocal) std::fill(nbh_extrapolation_grades, nbh_extrapolation_grades + nlocal, 0.0);
  }

  // If are in configuration, we need to reset the working array once per compute call / config
  if (configuration_mode)
    std::fill(energy_ders_wrt_coeffs, energy_ders_wrt_coeffs + coeff_count, 0.0);

  // Loop over all provided neighbourhoods
  for (int ii = 0; ii < inum; ii++) {
    int valid_count = 0;
    const int i = ilist[ii];
    const int itype = map[type[i]];
    const int jnum = numneigh[i];
    double nbh_energy = 0;
    const double xi[3] = {x[i][0], x[i][1], x[i][2]};

    // Resize per neighbor arrays
    if (cache_size < jnum) {
      memory->grow(neighbor_cache, jnum, 1 + 2 * radial_func_count, "neighbor_cache");
      memory->grow(cached_j, jnum, "cached_j");
      cache_size = jnum;
    }
    if (radial_basis_cache_size < jnum) {
      memory->grow(radial_basis_cache, jnum, radial_basis_size, "radial_basis_cache");
      radial_basis_cache_size = jnum;
    }

    // Reset the working arrays
    std::fill(moment_tensor_vals, moment_tensor_vals + alpha_moment_count, 0.0);
    std::fill(nbh_energy_ders_wrt_moments, nbh_energy_ders_wrt_moments + alpha_moment_count, 0.0);

    if (!configuration_mode)
      std::fill(energy_ders_wrt_coeffs, energy_ders_wrt_coeffs + coeff_count, 0.0);

    // ------------ Calculate Basic Moments ------------
    for (int jj = 0; jj < jnum; jj++) {
      int j = firstneigh[i][jj];
      j &= NEIGHMASK;
      const int jtype = map[type[j]];
      const double r[3] = {x[j][0] - xi[0], x[j][1] - xi[1], x[j][2] - xi[2]};
      const double rsq = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];

      if (rsq > max_cutoff_sq) continue;
      cached_j[valid_count] = j;

      const double dist = std::sqrt(rsq);
      const double inv_dist = 1.0 / dist;
      const double u[3] = {r[0] * inv_dist, r[1] * inv_dist, r[2] * inv_dist};
      neighbor_cache[valid_count][0] = inv_dist;
      double *vals = neighbor_cache[valid_count] + 1;
      double *ders = vals + radial_func_count;
      radial_basis->calc_radial_basis_ders(dist);
      const double *basis_vals = radial_basis->radial_basis_vals;
      const double *basis_ders = radial_basis->radial_basis_ders;
      std::copy(basis_vals, basis_vals + radial_basis_size, radial_basis_cache[valid_count]);

      // Evaluate each shared angular monomial once.
      for (int k = 1; k < angular_count; k++)
        angular_vals[k] = angular_vals[angular_parent[k]] * u[angular_axis[k]];

      // Compute the radial basis values and derivatives
      const int pair_offset = itype * species_count + jtype;
      for (int mu = 0; mu < radial_func_count; mu++) {
        double val = 0;
        double der = 0;
        const int offset = (pair_offset * radial_coeff_count_per_pair) + mu * radial_basis_size;

        for (int ri = 0; ri < radial_basis_size; ri++) {
          val += radial_basis_coeffs[offset + ri] * basis_vals[ri];
          der += radial_basis_coeffs[offset + ri] * basis_ders[ri];
        }
        vals[mu] = val;
        ders[mu] = der;

        // Accumulate into the basic moment elements
        for (int t = mu_offsets[mu]; t < mu_offsets[mu + 1]; t++) {
          const int k = basic_by_mu[t];
          const double ang = angular_vals[angular_by_mu[t]];
          moment_tensor_vals[k] += val * ang;
        }
      }
      valid_count++;
    }

    // ------------ Construct Composite Moment Values  ------------
    for (int k = 0; k < alpha_index_times_count; k++) {
      const int *term = times_data + 4 * k;
      moment_tensor_vals[term[3]] +=
          term[2] * moment_tensor_vals[term[0]] * moment_tensor_vals[term[1]];
    }

    // ------------ If Energies Are Needed Compute Basis Set From Alpha Map ------------
    const int linear_basis_offset = radial_coeff_count + species_count;
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

    energy_ders_wrt_coeffs[radial_coeff_count + itype] += 1;

    // =========== Begin Backpropagation ===========
    //------------ NBH energy derivative is the corresponding linear combination------------
    for (int k = 0; k < alpha_scalar_count; k++)
      nbh_energy_ders_wrt_moments[alpha_moment_mapping[k]] = linear_coeffs[k];

    //------------ Propagate chain rule through the composite moment elements times to the basics ------------
    for (int k = alpha_index_times_count - 1; k >= 0; k--) {
      const int *term = times_data + 4 * k;
      const int a0 = term[0];
      const int a1 = term[1];

      const double w = term[2] * nbh_energy_ders_wrt_moments[term[3]];

      nbh_energy_ders_wrt_moments[a1] += w * moment_tensor_vals[a0];
      nbh_energy_ders_wrt_moments[a0] += w * moment_tensor_vals[a1];
    }

    for (int t = 0; t < alpha_index_basic_count; t++)
      basic_ders_by_mu[t] = nbh_energy_ders_wrt_moments[basic_by_mu[t]];

    //------------ Compute forces from basic moment derivatives ------------
    double fi[3] = {f[i][0], f[i][1], f[i][2]};
    for (int jj = 0; jj < valid_count; jj++) {
      int j = cached_j[jj];
      const double inv_dist = neighbor_cache[jj][0];
      const double u[3] = {(x[j][0] - xi[0]) * inv_dist, (x[j][1] - xi[1]) * inv_dist,
                           (x[j][2] - xi[2]) * inv_dist};
      const double *vals = neighbor_cache[jj] + 1;
      const double *ders = vals + radial_func_count;
      const double *basis_vals = radial_basis_cache[jj];
      double *pair_ders = energy_ders_wrt_coeffs +
          (itype * species_count + map[type[j]]) * radial_coeff_count_per_pair;
      for (int k = 1; k < angular_count; k++)
        angular_vals[k] = angular_vals[angular_parent[k]] * u[angular_axis[k]];
      std::fill(angular_ders, angular_ders + angular_count, 0.0);

      double temp_force[3] = {0, 0, 0};
      double radial_force = 0;
      for (int mu = 0; mu < radial_func_count; mu++) {
        const int end = mu_offsets[mu + 1];
        if (mu_offsets[mu] == end) continue;
        const double val = vals[mu];
        double radial_sum = 0;
        for (int t = mu_offsets[mu]; t < end; t++) {
          const int angular = angular_by_mu[t];
          const double adj = basic_ders_by_mu[t];
          radial_sum += adj * angular_vals[angular];
          angular_ders[angular] += val * adj;
        }
        radial_force += ders[mu] * radial_sum;
        double *coeff_ders = pair_ders + mu * radial_basis_size;
        for (int ri = 0; ri < radial_basis_size; ri++)
          coeff_ders[ri] += radial_sum * basis_vals[ri];
      }
      // Reverse the shared angular products; no division by components of u.
      for (int k = angular_count - 1; k > 0; k--) {
        const int parent = angular_parent[k];
        const int axis = angular_axis[k];
        const double adj = angular_ders[k];
        temp_force[axis] += adj * angular_vals[parent];
        angular_ders[parent] += adj * u[axis];
      }
      radial_force -=
          inv_dist * (temp_force[0] * u[0] + temp_force[1] * u[1] + temp_force[2] * u[2]);
      for (int a = 0; a < 3; a++) temp_force[a] = inv_dist * temp_force[a] + radial_force * u[a];

      fi[0] += temp_force[0];
      fi[1] += temp_force[1];
      fi[2] += temp_force[2];

      f[j][0] -= temp_force[0];
      f[j][1] -= temp_force[1];
      f[j][2] -= temp_force[2];

      // Accumulate virial stress only if requested
      if (vflag_either) {
        const double del[3] = {xi[0] - x[j][0], xi[1] - x[j][1], xi[2] - x[j][2]};
        ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, temp_force[0], temp_force[1],
                     temp_force[2], del[0], del[1], del[2]);
      }
    }
    f[i][0] = fi[0];
    f[i][1] = fi[1];
    f[i][2] = fi[2];

    // Directly calculate extrapolation grade for neighbourhood mode
    if (!configuration_mode) {
      double grade = calculate_extrapolation_grade(itype);
      max_grade = std::max(grade, max_grade);
      nbh_extrapolation_grades[i] = grade;
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();

  compile_grades();
  if (mlip3_style) evaluate_grades();    // Evaluate grades per MLIP-3 two-threshold style
}

/* ----------------------------------------------------------------------
   Extrapolation Calculation Function
------------------------------------------------------------------------- */
double PairMTPExtrapolation::calculate_extrapolation_grade(int itype)
{
  const int begin = itype < 0 ? 0 : itype * species_count * radial_coeff_count_per_pair;
  const int end = itype < 0 ? coeff_count : begin + species_count * radial_coeff_count_per_pair;
  const int linear_offset = radial_coeff_count + species_count;
  double max_grade = 0;
  for (int i = 0; i < coeff_count; i++) {
    double current_grade = 0;
    const double *row = inverse_active_set[i];
    for (int j = begin; j < end; j++) { current_grade += energy_ders_wrt_coeffs[j] * row[j]; }
    if (itype >= 0) {
      const int species_offset = radial_coeff_count + itype;
      current_grade += energy_ders_wrt_coeffs[species_offset] * row[species_offset];
      for (int j = linear_offset; j < coeff_count; j++)
        current_grade += energy_ders_wrt_coeffs[j] * row[j];
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
  if (configuration_mode) {    // Configuration mode

    // Sum contributions across all processes
    MPI_Allreduce(MPI_IN_PLACE, energy_ders_wrt_coeffs, coeff_count, MPI_DOUBLE, MPI_SUM, world);
    max_grade = calculate_extrapolation_grade();

    if (atom->natoms > 0)
      max_grade /= std::pow((double) atom->natoms, 0.5 * weight_scaling);    // Normalize
    else
      max_grade = 0.0;

  } else {    // Neighbourhood mode
    MPI_Allreduce(MPI_IN_PLACE, &max_grade, 1, MPI_DOUBLE, MPI_MAX, world);
  }

  // ComputePair MPI_SUMs pvector, so only rank 0 sets it; the rest stay 0.0 and
  // the sum yields the already-Allreduced max. Can be improved.
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
    preselected_file = nullptr;
    error->one(FLERR, "Exceeded Break Threshold: {:.5f}. Terminating simulation.\n", max_grade);
  }
}
/* ----------------------------------------------------------------------
   Write current config to file
------------------------------------------------------------------------- */

// This function will likely be remove and configurations should be written with dump instead.
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
    // ii instead of ilist[ii] Avoids mirroring the device ilist to the host on the KOKKOS path.
    // May be affected by lost atoms and sorting issues.
    // The whole function will likely be removed in the future.
    const int i = ii;
    const int itype = map[type[i]];
    const double xi[3] = {x[i][0], x[i][1], x[i][2]};
    const int global_i = i + index_offset + 1;

    int n = 0;
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
  mlip3_style = false;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "mlip3_style") == 0) {
      if (iarg + 3 >= narg)
        utils::missing_cmd_args(FLERR, "pair_style mtp/extrapolation mlip3_style", error);

      mlip3_style = true;
      if (comm->me == 0) {
        if (preselected_file) std::fclose(preselected_file);
        preselected_file = std::fopen(arg[iarg + 1], "w");
        if (!preselected_file)
          error->one(FLERR, "Cannot open mtp/extrapolation output file {}: {}", arg[iarg + 1],
                     utils::getsyserror());
      }
      select_threshold = utils::numeric(FLERR, arg[iarg + 2], true, lmp);
      break_threshold = utils::numeric(FLERR, arg[iarg + 3], true, lmp);
      iarg += 4;
    } else {
      error->all(FLERR, "Unknown pair_style mtp/extrapolation keyword: {}", arg[iarg]);
    }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */
void PairMTPExtrapolation::coeff(int narg, char **arg)
{
  const int n = atom->ntypes;
  if (narg != 3 + n) error->all(FLERR, "Incorrect args for pair coefficients.");

  // Read in MTP and allocate memory
  FILE *mtp_file = nullptr;
  if (comm->me == 0) {
    mtp_file = utils::open_potential(arg[2], lmp, nullptr);
    if (mtp_file == nullptr)
      error->one(FLERR, "Cannot open MTP potential file {}: {}", arg[2], utils::getsyserror());
  }
  PairMTPExtrapolation::read_file(mtp_file);
  if (mtp_file) fclose(mtp_file);

  if (comm->me == 0) {
    if (mlip3_style)
      utils::logmesg(lmp,
                     "Extrapolation Scheme: {} mode, with a selection threshold of {} "
                     "and break threshold of {}.\n",
                     (configuration_mode ? "Configuration" : "Neighborhood"), select_threshold,
                     break_threshold);
    else
      utils::logmesg(lmp, "Extrapolation Mode: {} mode.\n",
                     (configuration_mode ? "Configuration" : "Neighborhood"));
  }

  PairMTP::prepare_map(narg - 3, arg + 3);
}

/* ----------------------------------------------------------------------
   MTP file parsing helper function. Includes memory allocation. Excludes some radial basis hyperparameters (in radial basis constructor instead).
------------------------------------------------------------------------- */
void PairMTPExtrapolation::read_file(FILE *mtp_file)
{
  PairMTP::read_file(mtp_file);

  coeff_count = radial_coeff_count + species_count + alpha_scalar_count;
  int num_doubles = coeff_count * coeff_count;
  radial_basis_cache_size = 0;

  // Now we allocate memory for the additional memory needed for calculations
  memory->create(active_set, coeff_count, coeff_count, "active_set");
  memory->create(inverse_active_set, coeff_count, coeff_count, "inverse_active_set");
  memory->create(energy_ders_wrt_coeffs, coeff_count, "energy_ders_wrt_coeffs");

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
      error->one(
          FLERR,
          "Error in reading MTP file selection state. Please verify MVS version is #MVS_v1.1!");
    tfr.ignore_comments = true;    // Accept comments after reading the version which is a comment

    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
    if (keyword != "energy_weight") error->one(FLERR, "Error in reading MTP file, energy_weight");
    const double energy_weight = line_tokens.next_double();

    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
    if (keyword != "force_weight") error->one(FLERR, "Error in reading MTP file, force_weight");
    const double force_weight = line_tokens.next_double();

    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
    if (keyword != "stress_weight") error->one(FLERR, "Error in reading MTP file, stress_weight");
    const double stress_weight = line_tokens.next_double();

    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
    if (keyword != "site_en_weight") error->one(FLERR, "Error in reading MTP file, site_en_weight");
    const double site_en_weight = line_tokens.next_double();

    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
    if (keyword != "weight_scaling") error->one(FLERR, "Error in reading MTP file, weight_scaling");
    weight_scaling = line_tokens.next_int();

    const bool cfg_mode = (energy_weight == 1 && site_en_weight == 0);
    const bool nbh_mode = (energy_weight == 0 && site_en_weight == 1);

    if ((!cfg_mode && !nbh_mode) || force_weight != 0 || stress_weight != 0)
      error->one(FLERR,
                 "Error, the MTP currently only supports configuration mode "
                 "(energy_weight=1) or neighbourhood mode (site_en_weight=1), "
                 "with force_weight=0 and stress_weight=0. Got energy_weight={}, "
                 "force_weight={}, stress_weight={}, site_en_weight={}. "
                 "Please retrain the MTP with the correct modes!",
                 energy_weight, force_weight, stress_weight, site_en_weight);

    configuration_mode = cfg_mode;

    fgetc(mtp_file);    // We need to skip foward 1 character. There is a # before the binary data.
    utils::sfread(FLERR, &active_set[0][0], sizeof(double), num_doubles, mtp_file, nullptr, error);
    utils::sfread(FLERR, &inverse_active_set[0][0], sizeof(double), num_doubles, mtp_file, nullptr,
                  error);
  }

  //Broadcast active set to others
  MPI_Bcast(&configuration_mode, 1, MPI_INT, 0, world);
  MPI_Bcast(&weight_scaling, 1, MPI_INT, 0, world);
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
