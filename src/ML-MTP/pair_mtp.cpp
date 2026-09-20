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
// Contributing author, Richard Meng, Queen's University at Kingston, 22.11.24, contact@richardzjm.com
//

#include "pair_mtp.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "info.h"
#include "memory.h"
#include "mtp_radial_basis.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "text_file_reader.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <vector>

using namespace LAMMPS_NS;

PairMTP::PairMTP(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;

  radial_basis = nullptr;
  radial_basis_coeffs = nullptr;
  linear_coeffs = nullptr;
  species_coeffs = nullptr;
  alpha_index_basic = nullptr;
  alpha_index_times = nullptr;
  alpha_moment_mapping = nullptr;
  nbh_energy_ders_wrt_moments = nullptr;
  basic_to_angular = nullptr;
  basic_by_mu = nullptr;
  angular_by_mu = nullptr;
  mu_offsets = nullptr;
  angular_parent = nullptr;
  angular_axis = nullptr;
  angular_vals = nullptr;
  angular_ders = nullptr;
  basic_ders_by_mu = nullptr;
  moment_tensor_vals = nullptr;
  cached_j = nullptr;
  neighbor_cache = nullptr;
  map = nullptr;

  cache_size = 0;
  angular_count = 0;
  force_index_times_count = 0;
  scaling = 1.0;
  potential_name = "Untitled";
  potential_tag = "";
}

PairMTP::~PairMTP()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(moment_tensor_vals);
    memory->destroy(radial_basis_coeffs);
    memory->destroy(linear_coeffs);
    memory->destroy(species_coeffs);
    memory->destroy(alpha_index_basic);
    memory->destroy(alpha_index_times);
    memory->destroy(alpha_moment_mapping);
    memory->destroy(nbh_energy_ders_wrt_moments);
    memory->destroy(cached_j);
    memory->destroy(neighbor_cache);

    memory->destroy(basic_to_angular);
    memory->destroy(basic_by_mu);
    memory->destroy(angular_by_mu);
    memory->destroy(mu_offsets);
    memory->destroy(angular_parent);
    memory->destroy(angular_axis);
    memory->destroy(angular_vals);
    memory->destroy(angular_ders);
    memory->destroy(basic_ders_by_mu);

    delete radial_basis;
    radial_basis = nullptr;

    delete[] map;
    map = nullptr;
  }

  if (elements) {
    for (int i = 0; i < nelements; i++) delete[] elements[i];
    delete[] elements;
    elements = nullptr;
  }
}

/* ----------------------------------------------------------------------
   Straightforward MTP implementation based on MLIP3
   ---------------------------------------------------------------------- */
void PairMTP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  const bool need_energy = eflag_atom || eflag_global;
  const int *times_data = alpha_index_times_count ? alpha_index_times[0] : nullptr;
  const int forward_count = need_energy ? alpha_index_times_count : force_index_times_count;

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

    // Clear moment and derivative arrays
    std::fill(moment_tensor_vals, moment_tensor_vals + alpha_moment_count, 0.0);
    std::fill(nbh_energy_ders_wrt_moments, nbh_energy_ders_wrt_moments + alpha_moment_count, 0.0);

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
          moment_tensor_vals[k] += val * angular_vals[angular_by_mu[t]];
        }
      }
      valid_count++;
    }

    // ------------ Construct Composite Moment Values  ------------
    for (int k = 0; k < forward_count; k++) {
      const int *term = times_data + 4 * k;
      moment_tensor_vals[term[3]] +=
          term[2] * moment_tensor_vals[term[0]] * moment_tensor_vals[term[1]];
    }

    // ------------ If Energies Are Needed Compute Basis Set From Alpha Map ------------
    if (eflag_atom || eflag_global) {
      nbh_energy = species_coeffs[itype];    // Essentially the reference point energy per species
      for (int k = 0; k < alpha_scalar_count; k++)
        nbh_energy += linear_coeffs[k] * moment_tensor_vals[alpha_moment_mapping[k]];

      if (eflag_atom) eatom[i] = nbh_energy;
      if (eflag_global) eng_vdwl += nbh_energy;
    }

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
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */
void PairMTP::settings(int narg, char **arg)
{
  if (comm->me == 0) {
    if (narg > 1)
      utils::logmesg(lmp,
                     "pair_style mtp does not accept arguments. Ignoring excessive "
                     "arguments!\n");
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */
void PairMTP::coeff(int narg, char **arg)
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
  read_file(mtp_file);
  if (mtp_file) fclose(mtp_file);

  prepare_map(narg - 3, arg + 3);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */
void PairMTP::init_style()
{
  if (force->newton_pair != 1) error->all(FLERR, "Pair style MTP requires Newton Pair on");

  // Request a full neighbourhood list which is needed for MTP
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMTP::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR, Error::NOLASTLINE,
               "All pair coeffs are not set. Status\n" + Info::get_pair_coeff_status(lmp));

  return radial_basis->max_cutoff;
}

/* ----------------------------------------------------------------------
   MTP file parsing helper function. Includes memory allocation. Excludes some radial basis hyperparameters (in radial basis constructor instead).
------------------------------------------------------------------------- */
void PairMTP::read_file(FILE *mtp_file)
{
  //Open the MTP file on proc 0
  if (comm->me == 0) {
    TextFileReader tfr(mtp_file, "ml-mtp");
    tfr.ignore_comments = true;
    const std::string new_separators = "=, ";
    const std::string separators = TOKENIZER_DEFAULT_SEPARATORS + new_separators;

    ValueTokenizer line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    std::string keyword = line_tokens.next_string();
    if (keyword != "MTP")    // Files checking
      error->one(FLERR, "Only MTP potential files are accepted.");

    // Version checking, tokenized so trailing whitespace and CRLF are tolerated
    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
    if (keyword != "version" || line_tokens.next_string() != "1.1.0")
      error->one(FLERR, "MTP file must have version \"1.1.0\"");

    // Read the potential name (optional field)
    line_tokens = ValueTokenizer(tfr.next_line(), separators);
    keyword = line_tokens.next_string();
    if (keyword == "potential_name") {
      try {
        potential_name = line_tokens.next_string();
      } catch (TokenizerException const &) {
        potential_name = "";
      }
      line_tokens = ValueTokenizer(tfr.next_line(), separators);
      keyword = line_tokens.next_string();
    }

    //Check the scaling
    if (keyword == "scaling") {
      scaling = line_tokens.next_double();
      line_tokens = ValueTokenizer(tfr.next_line(), separators);
      keyword = line_tokens.next_string();
    } else {
      scaling = 1;
    }

    // Read the species count
    if (keyword != "species_count")
      error->one(FLERR, "Error reading MTP file. Species count not found.");
    species_count = line_tokens.next_int();
    if (species_count < 1) error->one(FLERR, "MTP species count must be positive.");
    utils::logmesg(lmp, "There are {} species in this MTP.\n", species_count);

    // Read the potential tag (also optional field)
    line_tokens = ValueTokenizer(tfr.next_line(), separators);
    keyword = line_tokens.next_string();
    if (keyword == "potential_tag") {
      try {
        potential_tag = line_tokens.next_string();
      } catch (TokenizerException const &) {
        potential_tag = "";
      }
      line_tokens = ValueTokenizer(tfr.next_line(), separators);
      keyword = line_tokens.next_string();
    }

    // Read the radial basis type
    if (keyword != "radial_basis_type")
      error->one(FLERR, "Error reading MTP file. No radial basis set type is specified.");
    std::string radial_basis_type = line_tokens.next_string();

    // Set the type of radial basis.
    if (radial_basis_type == "RBChebyshev") {
      radial_basis = new RBChebyshev(tfr, lmp);
      radial_basis_type_index = 1;
    } else
      error->one(FLERR,
                 "Error reading MTP file. The specified radial basis set type, {}, was not "
                 "found/available.",
                 radial_basis_type);

    radial_basis->scaling = scaling;
    radial_basis_size = radial_basis->size;

    // Read the basis function count
    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
    if (keyword != "radial_funcs_count")
      error->one(FLERR, "Error in reading MTP file. Cannot read radial function count.");
    radial_func_count = line_tokens.next_int();
    if (radial_func_count < 1) error->one(FLERR, "MTP radial function count must be positive.");

    // Check for magnetic basis which is currently unsupported.
    line_tokens = ValueTokenizer(tfr.next_line(), separators);
    keyword = line_tokens.next_string();
    if (keyword != "radial_coeffs") {
      if (keyword == "magnetic_basis_type")
        error->one(FLERR, "Magnetic basis is currently not supported.");
      else
        error->one(FLERR, "Error in reading MTP file. Cannot read radial coeffs.");
    }

    // Allocate memory for radial basis
    int pairs_count = species_count * species_count;
    int coeffs_per_pair = radial_basis_size * radial_func_count;
    memory->create(radial_basis_coeffs, pairs_count * coeffs_per_pair, "radial_basis_coeffs");

    // Read the radial basis coeffs
    for (int i = 0; i < pairs_count; i++) {

      //Read which pairs are being allocated
      line_tokens = ValueTokenizer(tfr.next_line(), separators + "-");
      int type1 = line_tokens.next_int();
      int type2 = line_tokens.next_int();
      if (type1 < 0 || type1 >= species_count || type2 < 0 || type2 >= species_count)
        error->one(FLERR, "Invalid species pair {}-{} in MTP radial coefficients.", type1, type2);

      // Read the coeffs for the pair with offset in the array pointer.
      int pair_offset = (type1 * species_count + type2) * coeffs_per_pair;
      for (int j = 0; j < radial_func_count; j++) {
        line_tokens = ValueTokenizer(tfr.next_line(), separators + "{,}");
        for (int k = 0; k < radial_basis_size; k++) {
          radial_basis_coeffs[pair_offset + (j * radial_basis_size) + k] =
              line_tokens.next_double();
        }
      }
    }

    // Get the total alpha count
    line_tokens = ValueTokenizer(tfr.next_line(), separators);
    keyword = line_tokens.next_string();
    if (keyword != "alpha_moments_count")
      error->one(FLERR, "Error reading MTP file. Alpha moment count not found.");
    alpha_moment_count = line_tokens.next_int();
    if (alpha_moment_count < 1) error->one(FLERR, "MTP alpha moment count must be positive.");
    memory->create(moment_tensor_vals, alpha_moment_count, "moment_tensor_vals");
    memory->create(nbh_energy_ders_wrt_moments, alpha_moment_count, "nbh_energy_ders_wrt_moments");

    // Get the basic alpha count
    line_tokens = ValueTokenizer(tfr.next_line(), separators);
    keyword = line_tokens.next_string();
    if (keyword != "alpha_index_basic_count")
      error->one(FLERR, "Error reading MTP file. Alpha moment count not found.");
    alpha_index_basic_count = line_tokens.next_int();
    if (alpha_index_basic_count < 1 || alpha_index_basic_count > alpha_moment_count)
      error->one(FLERR, "MTP alpha index basic count is out of range.");

    // Read the basic alphas
    int radial_func_max = 0;
    tfr.set_bufsize((alpha_index_basic_count * 20 + 20) * sizeof(char));

    line_tokens = ValueTokenizer(tfr.next_line(), separators + "{},");
    keyword = line_tokens.next_string();
    if (keyword != "alpha_index_basic")
      error->one(FLERR, "Error reading MTP file. Alpha index basic not found.");
    memory->create(alpha_index_basic, alpha_index_basic_count, 4, "alpha_index_basic");
    for (int i = 0; i < alpha_index_basic_count; i++) {
      for (int j = 0; j < 4; j++) {
        int index = line_tokens.next_int();
        alpha_index_basic[i][j] = index;
      }
      if (alpha_index_basic[i][0] < 0 || alpha_index_basic[i][0] >= radial_func_count)
        error->one(FLERR, "Radial function index out of range in alpha_index_basic.");
      if (alpha_index_basic[i][1] < 0 || alpha_index_basic[i][2] < 0 || alpha_index_basic[i][3] < 0)
        error->one(FLERR, "Negative angular exponent in alpha_index_basic.");
      if (alpha_index_basic[i][0] > radial_func_max) radial_func_max = alpha_index_basic[i][0];
    }
    if (radial_func_max != radial_func_count - 1)    //Index validity check
      error->one(FLERR, "Wrong number of radial functions specified!");

    //Find the maximum alpha basic index
    max_alpha_index_basic = 0;
    for (int i = 0; i < alpha_index_basic_count; i++)
      max_alpha_index_basic =
          std::max(max_alpha_index_basic,
                   alpha_index_basic[i][1] + alpha_index_basic[i][2] + alpha_index_basic[i][3]);
    max_alpha_index_basic++;    // Add 1 to account for zeroth order indicies

    // Get the alpha times count
    line_tokens = ValueTokenizer(tfr.next_line(), separators);
    keyword = line_tokens.next_string();
    if (keyword != "alpha_index_times_count")
      error->one(FLERR, "Error reading MTP file. Alpha index times count not found.");
    alpha_index_times_count = line_tokens.next_int();
    if (alpha_index_times_count < 0) error->one(FLERR, "MTP alpha index times count is negative.");

    // Read the alphas times
    tfr.set_bufsize((alpha_index_times_count * 40 + 20) * sizeof(char));
    line_tokens = ValueTokenizer(tfr.next_line(), separators + "{},");
    keyword = line_tokens.next_string();
    if (keyword != "alpha_index_times")
      error->one(FLERR, "Error reading MTP file. Alpha index times not found.");
    memory->create(alpha_index_times, alpha_index_times_count, 4, "alpha_index_times");
    for (int i = 0; i < alpha_index_times_count; i++) {
      for (int j = 0; j < 4; j++) { alpha_index_times[i][j] = line_tokens.next_int(); }
      for (int j = 0; j < 4; j++)
        if (j != 2 &&
            (alpha_index_times[i][j] < 0 || alpha_index_times[i][j] >= alpha_moment_count))
          error->one(FLERR, "Moment index out of range in alpha_index_times.");
    }

    // Get the alpha scalar count
    line_tokens = ValueTokenizer(tfr.next_line(), separators);
    keyword = line_tokens.next_string();
    if (keyword != "alpha_scalar_moments")
      error->one(FLERR, "Error reading MTP file. Alpha scalar moment count not found.");
    alpha_scalar_count = line_tokens.next_int();
    if (alpha_scalar_count < 0) error->one(FLERR, "MTP alpha scalar moment count is negative.");

    //Read the alpha moment mappings
    line_tokens = ValueTokenizer(tfr.next_line(), separators + "{},");
    keyword = line_tokens.next_string();
    if (keyword != "alpha_moment_mapping")
      error->one(FLERR, "Error reading MTP file. Alpha moment mappings not found.");
    memory->create(alpha_moment_mapping, alpha_scalar_count, "alpha_moment_mapping");
    for (int i = 0; i < alpha_scalar_count; i++) {
      alpha_moment_mapping[i] = line_tokens.next_int();
      if (alpha_moment_mapping[i] < 0 || alpha_moment_mapping[i] >= alpha_moment_count)
        error->one(FLERR, "Moment index out of range in alpha_moment_mapping.");
    }

    //Read the species coefficients
    line_tokens = ValueTokenizer(tfr.next_line(), separators + "{},");
    keyword = line_tokens.next_string();
    if (keyword != "species_coeffs")
      error->one(FLERR, "Error reading MTP file. Species coefficients not found.");
    memory->create(species_coeffs, species_count, "species_coeffs");
    for (int i = 0; i < species_count; i++) { species_coeffs[i] = line_tokens.next_double(); }

    //Read the linear MTP basis coefficients
    line_tokens = ValueTokenizer(tfr.next_line(), separators + "{},");
    keyword = line_tokens.next_string();
    if (keyword != "moment_coeffs")
      error->one(FLERR, "Error reading MTP file. Moment coefficients not found.");
    memory->create(linear_coeffs, alpha_scalar_count, "moment_coeffs");
    for (int i = 0; i < alpha_scalar_count; i++) { linear_coeffs[i] = line_tokens.next_double(); }
  }    // Proc 0

  // ---------- Now broadcast to all the other procs ----------
  //Radial Basis Set Type First
  MPI_Bcast(&radial_basis_type_index, 1, MPI_INT, 0,
            world);    //index of the radial basis

  //Then Single Values
  MPI_Bcast(&scaling, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&species_count, 1, MPI_INT, 0, world);
  MPI_Bcast(&radial_basis_size, 1, MPI_INT, 0, world);
  MPI_Bcast(&radial_func_count, 1, MPI_INT, 0, world);
  MPI_Bcast(&alpha_moment_count, 1, MPI_INT, 0, world);
  MPI_Bcast(&alpha_index_basic_count, 1, MPI_INT, 0, world);
  MPI_Bcast(&max_alpha_index_basic, 1, MPI_INT, 0, world);
  MPI_Bcast(&alpha_index_times_count, 1, MPI_INT, 0, world);
  MPI_Bcast(&alpha_scalar_count, 1, MPI_INT, 0, world);

  const int pairs_count = species_count * species_count;
  radial_coeff_count_per_pair = radial_basis_size * radial_func_count;
  radial_coeff_count = pairs_count * radial_coeff_count_per_pair;

  // Now we allocate memory for all the arrays.
  if (comm->me != 0) {    // Non-zero proc
    //First we reconstruct the radial basis set
    if (radial_basis_type_index == 1) {
      radial_basis = new RBChebyshev(radial_basis_size, lmp);
      radial_basis->scaling = scaling;
    }

    //Alpha index
    memory->create(alpha_index_basic, alpha_index_basic_count, 4, "alpha_index_basic");
    memory->create(alpha_index_times, alpha_index_times_count, 4, "alpha_index_times");
    memory->create(alpha_moment_mapping, alpha_scalar_count, "alpha_moment_mapping");

    //Working buffers
    memory->create(moment_tensor_vals, alpha_moment_count, "moment_tensor_vals");
    memory->create(nbh_energy_ders_wrt_moments, alpha_moment_count, "nbh_energy_ders_wrt_moments");

    //Coefficients
    memory->create(radial_basis_coeffs, radial_coeff_count, "radial_basis_coeffs");
    memory->create(linear_coeffs, alpha_scalar_count, "linear_coeffs");
    memory->create(species_coeffs, species_count, "species_coeffs");
  }

  //We can then populate the cutoffs
  MPI_Bcast(&radial_basis->min_cutoff, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&radial_basis->max_cutoff, 1, MPI_DOUBLE, 0, world);
  min_cutoff = radial_basis->min_cutoff;
  max_cutoff = radial_basis->max_cutoff;
  max_cutoff_sq = max_cutoff * max_cutoff;

  // Now we B Cast arrays
  // Alphas
  MPI_Bcast(&alpha_index_basic[0][0], alpha_index_basic_count * 4, MPI_INT, 0, world);
  if (alpha_index_times_count)
    MPI_Bcast(&alpha_index_times[0][0], alpha_index_times_count * 4, MPI_INT, 0, world);
  MPI_Bcast(alpha_moment_mapping, alpha_scalar_count, MPI_INT, 0, world);

  // Coefficients
  MPI_Bcast(radial_basis_coeffs, radial_coeff_count, MPI_DOUBLE, 0, world);
  MPI_Bcast(linear_coeffs, alpha_scalar_count, MPI_DOUBLE, 0, world);
  MPI_Bcast(species_coeffs, species_count, MPI_DOUBLE, 0, world);

  // Set working buffers
  prepare_angular();
  memory->create(mu_offsets, radial_func_count + 1, "mu_offsets");
  memory->create(basic_by_mu, alpha_index_basic_count, "basic_by_mu");
  memory->create(angular_by_mu, alpha_index_basic_count, "angular_by_mu");
  memory->create(basic_ders_by_mu, alpha_index_basic_count, "basic_ders_by_mu");
  std::fill(mu_offsets, mu_offsets + radial_func_count + 1, 0);
  for (int k = 0; k < alpha_index_basic_count; k++) mu_offsets[alpha_index_basic[k][0] + 1]++;
  for (int mu = 0; mu < radial_func_count; mu++) mu_offsets[mu + 1] += mu_offsets[mu];
  std::vector<int> radial_next(mu_offsets, mu_offsets + radial_func_count);
  for (int k = 0; k < alpha_index_basic_count; k++)
    basic_by_mu[radial_next[alpha_index_basic[k][0]]++] = k;
  for (int t = 0; t < alpha_index_basic_count; t++)
    angular_by_mu[t] = basic_to_angular[basic_by_mu[t]];

  // Sanity check the contraction graph: basic moments own [0, alpha_index_basic_count)
  // and no term may read a moment that has not been produced yet.
  std::vector<char> produced(alpha_moment_count, 0);
  std::fill(produced.begin(), produced.begin() + alpha_index_basic_count, 1);
  for (int k = 0; k < alpha_index_times_count; k++) {
    if (alpha_index_times[k][3] < alpha_index_basic_count)
      error->all(FLERR, "MTP contraction {} overwrites a basic moment", k);
    if (!produced[alpha_index_times[k][0]] || !produced[alpha_index_times[k][1]])
      error->all(FLERR, "MTP contraction {} reads a moment that is not yet computed", k);
    produced[alpha_index_times[k][3]] = 1;
  }

  // Contractions whose product is never consumed are dead when only forces are wanted.
  // Sink them below the live terms so both passes stream one table: the forward pass
  // stops at force_index_times_count, the reverse pass walks the whole list. A dead
  // output feeds nothing, so the partitioned order stays topologically valid.
  std::vector<char> moment_used(alpha_moment_count, 0);
  for (int k = 0; k < alpha_index_times_count; k++) {
    moment_used[alpha_index_times[k][0]] = 1;
    moment_used[alpha_index_times[k][1]] = 1;
  }

  std::vector<std::array<int, 4>> dead;
  int live = 0;
  for (int k = 0; k < alpha_index_times_count; k++) {
    const int *term = alpha_index_times[k];
    if (moment_used[term[3]]) {
      if (live != k) std::copy(term, term + 4, alpha_index_times[live]);
      live++;
    } else {
      dead.push_back({{term[0], term[1], term[2], term[3]}});
    }
  }
  force_index_times_count = live;
  for (const auto &term : dead) std::copy(term.begin(), term.end(), alpha_index_times[live++]);

  cache_size = 0;

  allocated = 1;
}

/* ----------------------------------------------------------------------
   Share angular monomials across radial channels. Build on every MPI rank
   after alpha_index_basic is broadcast; the supplied contraction graph is
   untouched. Lexicographic exponent order puts every parent before its child.
------------------------------------------------------------------------- */
void PairMTP::prepare_angular()
{
  using Powers = std::array<int, 3>;
  std::map<Powers, int> indices;
  indices.emplace(Powers{{0, 0, 0}}, 0);

  for (int k = 0; k < alpha_index_basic_count; k++) {
    Powers powers{{alpha_index_basic[k][1], alpha_index_basic[k][2], alpha_index_basic[k][3]}};
    if (powers[0] < 0 || powers[1] < 0 || powers[2] < 0)
      error->all(FLERR, "Invalid negative MTP angular exponent");
    while (indices.emplace(powers, 0).second) {
      const int axis = powers[2] ? 2 : (powers[1] ? 1 : 0);
      --powers[axis];
    }
  }

  angular_count = 0;
  for (auto &entry : indices) entry.second = angular_count++;
  memory->create(basic_to_angular, alpha_index_basic_count, "basic_to_angular");
  memory->create(angular_parent, angular_count, "angular_parent");
  memory->create(angular_axis, angular_count, "angular_axis");
  memory->create(angular_vals, angular_count, "angular_vals");
  memory->create(angular_ders, angular_count, "angular_ders");
  angular_parent[0] = angular_axis[0] = 0;
  angular_vals[0] = 1.0;

  for (const auto &entry : indices) {
    const int index = entry.second;
    if (index == 0) continue;
    Powers parent = entry.first;
    const int axis = parent[2] ? 2 : (parent[1] ? 1 : 0);
    --parent[axis];
    angular_parent[index] = indices.at(parent);
    angular_axis[index] = axis;
  }
  for (int k = 0; k < alpha_index_basic_count; k++)
    basic_to_angular[k] = indices.at(
        Powers{{alpha_index_basic[k][1], alpha_index_basic[k][2], alpha_index_basic[k][3]}});
}

/* ----------------------------------------------------------------------
   Sets up the MTP element mapping and allocates memory for cutoffs and setflag
------------------------------------------------------------------------- */

void PairMTP::prepare_map(int narg, char **arg)
{
  const int n = atom->ntypes;
  const int np1 = n + 1;

  // Set up basic flags
  memory->create(setflag, np1, np1, "pair:setflag");
  memory->create(cutsq, np1, np1, "pair:cutsq");
  map = new int[np1];
  map_element2type(narg, arg);

  // Readjust Map
  for (int i = 1; i < np1; i++) {
    std::string entry = arg[i - 1];
    if (entry == "NULL") {
      map[i] = -1;
    } else {
      int type = utils::inumeric(FLERR, entry, false, lmp);
      if (type >= species_count || type < 0)
        error->all(FLERR, "The provided MTP does not support atom type {}!", type);
      map[i] = type;
    }
  }
  // MTP uses a single global cutoff
  std::fill(&cutsq[0][0], &cutsq[0][0] + np1 * np1, max_cutoff_sq);
}
