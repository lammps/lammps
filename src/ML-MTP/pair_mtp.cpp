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

#include "mtp_radial_basis.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;

PairMTP::PairMTP(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  no_virial_fdotr_compute = 0;
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
    memory->destroy(moment_jacobian);
    memory->destroy(nbh_energy_ders_wrt_moments);
    memory->destroy(within_cutoff);

    delete radial_basis;
    radial_basis = nullptr;
  }
}

/* ----------------------------------------------------------------------
   Straightfoward MTP implementation based on MLIP3
   ---------------------------------------------------------------------- */
void PairMTP::compute(int eflag, int vflag)
{
  ev_setup(eflag, vflag);

  double **x = atom->x;      // atomic positons
  double **f = atom->f;      // atomic forces
  int *type = atom->type;    //atomic types

  int inum = list->inum;             // The number of central atoms (neigbhourhoods)
  int *ilist = list->ilist;          // List of central atom ids
  int *numneigh = list->numneigh;    // List of the number of neighbours for each central atom
  int **firstneigh =
      list->firstneigh;    //List  (head of array) of neighbours for a given central atom

  // Loop over all provided neighbourhoods
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    const int itype = type[i] - 1;    // Set central atom type. Convert back to zero indexing.
    if (itype >= species_count) error->one(FLERR, "Too few species count in the MTP potential!");
    const int jnum = numneigh[i];
    double nbh_energy = 0;
    const double xi[3] = {x[i][0], x[i][1], x[i][2]};

    // Resize per neighbor arrays
    if (jac_size < jnum) {
      memory->grow(moment_jacobian, jnum, alpha_index_basic_count, 3, "moment_jacobian");
      memory->grow(within_cutoff, jnum, "within_cutoff");
      jac_size = jnum;
    }

    // Clear moment and derivative arrays
    std::fill(&moment_tensor_vals[0], &moment_tensor_vals[0] + alpha_moment_count, 0.0);
    std::fill(&nbh_energy_ders_wrt_moments[0], &nbh_energy_ders_wrt_moments[0] + alpha_moment_count,
              0.0);

    // ------------ Calculate Basic Moments ------------
    for (int jj = 0; jj < jnum; jj++) {
      int j = firstneigh[i][jj];
      j &= NEIGHMASK;
      const int jtype = type[j] - 1;    // Convert to zero indexing
      if (jtype >= species_count) error->one(FLERR, "Too few species count in the MTP potential!");

      const double r[3] = {x[j][0] - xi[0], x[j][1] - xi[1], x[j][2] - xi[2]};
      const double rsq = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];

      // Check cutoff and store for backwards pass
      if (rsq > cutsq[itype + 1][jtype + 1]) {    //1 indexing
        within_cutoff[jj] = false;
        continue;
      }
      within_cutoff[jj] = true;

      const double dist = std::sqrt(rsq);
      radial_basis->calc_radial_basis_ders(dist);

      // Precompute the coord and distance powers
      for (int k = 1; k < max_alpha_index_basic; k++) {
        dist_powers[k] = dist_powers[k - 1] * dist;
        for (int a = 0; a < 3; a++) coord_powers[k][a] = coord_powers[k - 1][a] * r[a];
      }

      // Compute the radial basis values and derivatives
      for (int mu = 0; mu < radial_func_count; mu++) {
        double val = 0;
        double der = 0;
        const int pair_offset = itype * species_count + jtype;
        const int offset = (pair_offset * radial_coeff_count_per_pair) + mu * radial_basis_size;

        for (int ri = 0; ri < radial_basis_size; ri++) {
          val += radial_basis_coeffs[offset + ri] * radial_basis->radial_basis_vals[ri];
          der += radial_basis_coeffs[offset + ri] * radial_basis->radial_basis_ders[ri];
        }
        radial_vals[mu] = val;
        radial_ders[mu] = der;
      }

      // Accumulate into the basic moment elements
      for (int k = 0; k < alpha_index_basic_count; k++) {
        double val = 0;
        double der = 0;
        int mu = alpha_index_basic[k][0];

        val = radial_vals[mu];
        der = radial_ders[mu];

        // Normalize by the rank of alpha's coresponding tensor
        int norm_rank = alpha_index_basic[k][1] + alpha_index_basic[k][2] + alpha_index_basic[k][3];
        double norm_fac = 1.0 / dist_powers[norm_rank];
        val *= norm_fac;
        der = der * norm_fac - norm_rank * val / dist;
        double pow0 = coord_powers[alpha_index_basic[k][1]][0];
        double pow1 = coord_powers[alpha_index_basic[k][2]][1];
        double pow2 = coord_powers[alpha_index_basic[k][3]][2];
        double pow = pow0 * pow1 * pow2;
        moment_tensor_vals[k] += val * pow;

        // Calculate the Jacobian from derivatives
        pow *= der / dist;
        moment_jacobian[jj][k][0] = pow * r[0];
        moment_jacobian[jj][k][1] = pow * r[1];
        moment_jacobian[jj][k][2] = pow * r[2];
        if (alpha_index_basic[k][1] != 0) {
          moment_jacobian[jj][k][0] += val * alpha_index_basic[k][1] *
              coord_powers[alpha_index_basic[k][1] - 1][0] * pow1 * pow2;
        }    //Chain rule for nonzero rank
        if (alpha_index_basic[k][2] != 0) {
          moment_jacobian[jj][k][1] += val * alpha_index_basic[k][2] * pow0 *
              coord_powers[alpha_index_basic[k][2] - 1][1] * pow2;
        }    //Chain rule for nonzero rank
        if (alpha_index_basic[k][3] != 0) {
          moment_jacobian[jj][k][2] += val * alpha_index_basic[k][3] * pow0 * pow1 *
              coord_powers[alpha_index_basic[k][3] - 1][2];
        }    //Chain rule for nonzero rank
      }
    }

    // ------------ Contruct Composite Moment Values  ------------
    for (int k = 0; k < alpha_index_times_count; k++) {
      double val0 = moment_tensor_vals[alpha_index_times[k][0]];
      double val1 = moment_tensor_vals[alpha_index_times[k][1]];
      int val2 = alpha_index_times[k][2];
      moment_tensor_vals[alpha_index_times[k][3]] += val2 * val0 * val1;
    }

    // ------------ If Energies Are Needed Compute Basis Set From Alpha Map ------------
    if (eflag_atom || eflag_global) {
      nbh_energy = species_coeffs[itype];    // Essentially the reference point energy per species
      for (int k = 0; k < alpha_scalar_count; k++)
        nbh_energy += linear_coeffs[k] * moment_tensor_vals[alpha_moment_mapping[k]];

      if (eflag_atom) eatom[i] = nbh_energy;
      if (eflag_global) eng_vdwl += nbh_energy;
    }

    // =========== Begin Backpropogation ===========
    //------------ NBH energy derivative is the corresponding linear combination------------
    for (int k = 0; k < alpha_scalar_count; k++)
      nbh_energy_ders_wrt_moments[alpha_moment_mapping[k]] = linear_coeffs[k];

    //------------ Propogate chain rule through the composite moment elements times to the basics ------------
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

    //------------  Multiply energy ders wrt basic moments by the Jacobian to get forces ------------
    for (int jj = 0; jj < jnum; jj++) {
      int j = firstneigh[i][jj];
      j &= NEIGHMASK;
      if (!within_cutoff[jj]) continue;

      double temp_force[3] = {0, 0, 0};
      for (int k = 0; k < alpha_index_basic_count; k++)
        for (int a = 0; a < 3; a++)
          temp_force[a] += nbh_energy_ders_wrt_moments[k] * moment_jacobian[jj][k][a];

      f[i][0] += temp_force[0];
      f[i][1] += temp_force[1];
      f[i][2] += temp_force[2];

      f[j][0] -= temp_force[0];
      f[j][1] -= temp_force[1];
      f[j][2] -= temp_force[2];

      // Accumulate virial stress only if requested
      if (vflag) {
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
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */
void PairMTP::settings(int narg, char **arg)
{
  if (comm->me == 0) {
    if (narg < 1)
      error->one(FLERR, "pair_style mtp only accepts 1 argument, the MTP potential file");
    if (narg > 1)
      utils::logmesg(
          lmp,
          "pair_style mtp only accepts 1 argument, the MTP potential file. Ignoring excessive "
          "arguments!\n");
  }
  FILE *mtp_file = utils::open_potential(arg[0], lmp, nullptr);
  read_file(mtp_file);
  fclose(mtp_file);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */
void PairMTP::coeff(int narg, char **arg)
{
  // The potential file is specified in the setting function instead.
  if (narg != 2) error->all(FLERR, "Only \"pair_coeff * *\" is permitted");
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
  if (setflag[i][j] == 0) error->all(FLERR, "Not all pair coeffs are set. See types {}-{}.", i, j);

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
    std::string version_line = std::string(tfr.next_line());
    if (version_line != "version = 1.1.0\n")    // Version checking
      error->one(FLERR, "MTP file must have version \"1.1.0\"");

    // Read the potential name (optional field)
    line_tokens = ValueTokenizer(tfr.next_line(), separators);
    keyword = line_tokens.next_string();
    if (keyword == "potential_name") {
      try {
        potential_name = line_tokens.next_string();
      } catch (TokenizerException e) {
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
    utils::logmesg(lmp, "There are {} species.\n", species_count);

    int np1 = species_count + 1;    // Lammps is 1 indexed instead of MLIP which is 0 indexed
    memory->create(setflag, np1, np1, "pair:setflag");
    memory->create(cutsq, np1, np1, "pair:cutsq");

    // Read the potential tag (also optional field)
    line_tokens = ValueTokenizer(tfr.next_line(), separators);
    keyword = line_tokens.next_string();
    if (keyword == "potential_tag") {
      try {
        potential_tag = line_tokens.next_string();
      } catch (TokenizerException e) {
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
                 "Error reading MTP file. The specified radial basis set type, {}, was not found..",
                 radial_basis_type);

    radial_basis->scaling = scaling;
    radial_basis_size = radial_basis->size;

    // Read the basis function count
    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
    if (keyword != "radial_funcs_count")
      lmp->error->one(FLERR, "Error in reading MTP file. Cannot read radial function count.");
    radial_func_count = line_tokens.next_int();

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
    int radial_coeff_count_per_pair = radial_basis_size * radial_func_count;
    memory->create(radial_basis_coeffs, pairs_count * radial_coeff_count_per_pair,
                   "radial_basis_coeffs");

    // Read the radial basis coeffs
    double rcutmaxsq = radial_basis->max_cutoff * radial_basis->max_cutoff;
    for (int i = 0; i < pairs_count; i++) {

      //Read which pairs are being allocated
      line_tokens = ValueTokenizer(tfr.next_line(), separators + "-");
      int type1 = line_tokens.next_int();
      int type2 = line_tokens.next_int();
      setflag[type1 + 1][type2 + 1] = 1;
      cutsq[type1 + 1][type2 + 1] = rcutmaxsq;

      // Read the coeffs for the pair with offset in the array pointer.
      int pair_offset = (type1 * species_count + type2) * radial_coeff_count_per_pair;
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
    memory->create(moment_tensor_vals, alpha_moment_count, "moment_tensor_vals");
    memory->create(nbh_energy_ders_wrt_moments, alpha_moment_count, "nbh_energy_ders_wrt_moments");

    // Get the basic alpha count
    line_tokens = ValueTokenizer(tfr.next_line(), separators);
    keyword = line_tokens.next_string();
    if (keyword != "alpha_index_basic_count")
      error->one(FLERR, "Error reading MTP file. Alpha moment count not found.");
    alpha_index_basic_count = line_tokens.next_int();

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

    // Read the alphas times
    tfr.set_bufsize((alpha_index_times_count * 40 + 20) * sizeof(char));
    line_tokens = ValueTokenizer(tfr.next_line(), separators + "{},");
    keyword = line_tokens.next_string();
    if (keyword != "alpha_index_times")
      error->one(FLERR, "Error reading MTP file. Alpha index times not found.");
    memory->create(alpha_index_times, alpha_index_times_count, 4, "alpha_index_times");
    for (int i = 0; i < alpha_index_times_count; i++) {
      for (int j = 0; j < 4; j++) { alpha_index_times[i][j] = line_tokens.next_int(); }
    }

    // Get the alpha scalar count
    line_tokens = ValueTokenizer(tfr.next_line(), separators);
    keyword = line_tokens.next_string();
    if (keyword != "alpha_scalar_moments")
      error->one(FLERR, "Error reading MTP file. Alpha scalar moment count not found.");
    alpha_scalar_count = line_tokens.next_int();

    //Read the alpha moment mappings
    line_tokens = ValueTokenizer(tfr.next_line(), separators + "{},");
    keyword = line_tokens.next_string();
    if (keyword != "alpha_moment_mapping")
      error->one(FLERR, "Error reading MTP file. Alpha moment mappings not found.");
    memory->create(alpha_moment_mapping, alpha_scalar_count, "alpha_moment_mapping");
    for (int i = 0; i < alpha_scalar_count; i++) {
      alpha_moment_mapping[i] = line_tokens.next_int();
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
  const int np1 = (species_count + 1);

  //Working buffers
  memory->create(dist_powers, max_alpha_index_basic, "dist_powers");
  memory->create(coord_powers, max_alpha_index_basic, 3, "coord_powers");
  memory->create(radial_vals, radial_func_count, "radial_vals");
  memory->create(radial_ders, radial_func_count, "radial_ders");

  // Now we allocate memory for all the arrays.
  if (comm->me != 0) {    // Non-zero proc
    //First we reconstruct the radial basis set
    if (radial_basis_type_index == 1) {
      radial_basis = new RBChebyshev(radial_basis_size, lmp);
      radial_basis->scaling = scaling;
    }

    //Flags
    memory->create(cutsq, np1, np1, "pair:cutsq");
    memory->create(setflag, np1, np1, "pair:setflag");

    //Alpha index
    memory->create(alpha_index_basic, alpha_index_basic_count, 4, "alpha_index_basic");
    memory->create(alpha_index_times, alpha_index_times_count, 4, "alpha_index_times");
    memory->create(alpha_moment_mapping, alpha_scalar_count, "alpha_moment_mapping");

    //Working buffers
    memory->create(moment_tensor_vals, alpha_moment_count, "moment_tensor_vals");
    memory->create(nbh_energy_ders_wrt_moments, alpha_moment_count, "nbh_energy_ders_wrt_moments");
    //Jacobian and within_cutoff will be first created with memory->grow during calculation.

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

  //Now we B Cast arrays
  //Flags
  MPI_Bcast(&cutsq[0][0], np1 * np1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&setflag[0][0], np1 * np1, MPI_INT, 0, world);

  // //Alphas
  MPI_Bcast(&alpha_index_basic[0][0], alpha_index_basic_count * 4, MPI_INT, 0, world);
  MPI_Bcast(&alpha_index_times[0][0], alpha_index_times_count * 4, MPI_INT, 0, world);
  MPI_Bcast(alpha_moment_mapping, alpha_scalar_count, MPI_INT, 0, world);

  // Coefficients
  MPI_Bcast(radial_basis_coeffs, radial_coeff_count, MPI_DOUBLE, 0, world);
  MPI_Bcast(linear_coeffs, alpha_scalar_count, MPI_DOUBLE, 0, world);
  MPI_Bcast(species_coeffs, species_count, MPI_DOUBLE, 0, world);

  // Set working buffers
  dist_powers[0] = coord_powers[0][0] = coord_powers[0][1] = coord_powers[0][2] = 1;

  allocated = 1;
}