/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#include "fix_electrode_conp.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "electrode_accel_interface.h"
#include "electrode_math.h"
#include "electrode_matrix.h"
#include "electrode_vector.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "pointers.h"
#include "text_file_reader.h"
#include "variable.h"

#include <cassert>
#include <numeric>

using namespace LAMMPS_NS;
using namespace MathConst;

extern "C" {
void dgetrf_(const int *M, const int *N, double *A, const int *lda, int *ipiv, int *info);
void dgetri_(const int *N, double *A, const int *lda, const int *ipiv, double *work,
             const int *lwork, int *info);
}

enum { CONST, EQUAL };
enum {
  V,      // voltage
  QSB,    // q_sb
  MC,     // macro_capacitance
  ME      // macro_elastance
};

//     0        1      2              3    4
// fix fxupdate group1 electrode/conp pot1 eta couple group2 pot2
FixElectrodeConp::FixElectrodeConp(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), f_inv(nullptr), f_mat(nullptr), f_vec(nullptr), array_compute(nullptr),
    ele_vector(nullptr), capacitance(nullptr), pair(nullptr), mat_neighlist(nullptr),
    vec_neighlist(nullptr), recvcounts(nullptr), displs(nullptr), iele_gathered(nullptr),
    buf_gathered(nullptr), potential_i(nullptr), potential_iele(nullptr), charge_iele(nullptr)
{
  read_inv = read_mat = false;
  symm = false;
  ffield = false;
  thermo_time = 0.;
  scalar_flag = 1;
  vector_flag = 1;
  top_group = 0;
  intelflag = false;
  tfflag = false;
  timer_flag = false;

  update_time = 0;
  mult_time = 0;

  // read fix command
  fixname = std::string(arg[0]);
  groups = std::vector<int>(1, igroup);
  group_bits = std::vector<int>(1, groupbit);
  group_psi_var_names = std::vector<std::string>(1);
  group_psi_var_styles = std::vector<int>(1, CONST);
  group_psi = std::vector<double>(1);
  etypes_neighlists = false;
  if (strstr(arg[3], "v_") == arg[3]) {
    std::string vname = arg[3];
    group_psi_var_names[0] = vname.substr(2);
    group_psi_var_styles[0] = EQUAL;
  } else
    group_psi[0] = utils::numeric(FLERR, arg[3], false, lmp);
  char *eta_str = arg[4];
  eta = utils::numeric(FLERR, eta_str, false, lmp);
  int iarg = 5;
  while (iarg < narg) {
    if ((strcmp(arg[iarg], "couple") == 0)) {
      if (iarg + 3 > narg) error->all(FLERR, "Need two arguments after couple keyword");
      int id = group->find(arg[++iarg]);
      if (id < 0) error->all(FLERR, "Group does not exist");
      groups.push_back(id);
      group_bits.push_back(group->bitmask[id]);
      ++iarg;
      if (strstr(arg[iarg], "v_") == arg[iarg]) {
        std::string vname = arg[iarg];
        group_psi_var_names.push_back(vname.substr(2));
        group_psi_var_styles.push_back(EQUAL);
        group_psi.push_back(0.);
      } else {
        std::string null;
        group_psi_var_names.push_back(null);
        group_psi_var_styles.push_back(CONST);
        group_psi.push_back(utils::numeric(FLERR, arg[iarg], false, lmp));
      }
    } else if ((strncmp(arg[iarg], "symm", 4) == 0)) {
      if (iarg + 2 > narg) error->all(FLERR, "Need yes/no command after symm keyword");
      char *symm_arg = arg[++iarg];
      if ((strcmp(symm_arg, "yes") == 0) || (strcmp(symm_arg, "on") == 0)) {
        symm = true;
      } else if ((strcmp(symm_arg, "no") == 0) || (strcmp(symm_arg, "off") == 0)) {
        symm = false;
      } else {
        error->all(FLERR, "Invalid argument after symm keyword");
      }
    } else if ((strncmp(arg[iarg], "write", 5) == 0)) {
      if (iarg + 2 > narg) error->all(FLERR, "Need one argument after write command");
      if (comm->me == 0) {
        if ((strcmp(arg[iarg], "write_inv") == 0)) {    // capacitance matrix
          f_inv = fopen(arg[++iarg], "w");
          if (f_inv == nullptr)
            error->one(FLERR,
                       fmt::format("Cannot open capacitance matrix file {}: {}", arg[iarg],
                                   utils::getsyserror()));
        } else if ((strcmp(arg[iarg], "write_mat") == 0)) {    // b vector
          f_mat = fopen(arg[++iarg], "w");
          if (f_mat == nullptr)
            error->one(FLERR,
                       fmt::format("Cannot open elastance matrix file {}: {}", arg[iarg],
                                   utils::getsyserror()));
        } else if ((strcmp(arg[iarg], "write_vec") == 0)) {    // b vector
          f_vec = fopen(arg[++iarg], "w");
          if (f_vec == nullptr)
            error->one(
                FLERR,
                fmt::format("Cannot open vector file {}: {}", arg[iarg], utils::getsyserror()));
        } else {
          error->all(FLERR, "Illegal fix electrode/conp command with write");
        }
      } else {
        iarg++;
      }
    } else if ((strncmp(arg[iarg], "read", 4) == 0)) {
      if (iarg + 2 > narg) error->all(FLERR, "Need one argument after read command");
      if ((strcmp(arg[iarg], "read_inv") == 0)) {
        read_inv = true;
        input_file_inv = arg[++iarg];
      } else if ((strcmp(arg[iarg], "read_mat") == 0)) {
        read_mat = true;
        input_file_mat = arg[++iarg];
      } else {
        error->all(FLERR, "Illegal fix electrode/conp command with read");
      }
    } else if ((strcmp(arg[iarg], "etypes") == 0)) {
      if (iarg + 2 > narg) error->all(FLERR, "Need one argument after etypes command");
      int ilo, ihi;
      utils::bounds(FLERR, arg[++iarg], 1, atom->ntypes, ilo, ihi, error);
      for (int i = ilo; i <= ihi; ++i) etypes.push_back(i);
      etypes_neighlists = true;
    } else if ((strcmp(arg[iarg], "temp") == 0)) {
      if (iarg + 4 > narg) error->all(FLERR, "Need three arguments after temp command");
      if (strcmp(this->style, "electrode/thermo") != 0)
        error->all(FLERR, "temp keyword not available for this electrode fix");
      thermo_temp = force->boltz / force->qe2f * utils::numeric(FLERR, arg[++iarg], false, lmp);
      thermo_time = utils::numeric(FLERR, arg[++iarg], false, lmp);
      thermo_init = utils::inumeric(FLERR, arg[++iarg], false, lmp);
    } else if ((strcmp(arg[iarg], "ffield") == 0)) {
      if (iarg + 2 > narg) error->all(FLERR, "Need yes/no command after ffield keyword");
      char *ffield_arg = arg[++iarg];
      if ((strcmp(ffield_arg, "yes") == 0) || (strcmp(ffield_arg, "on") == 0)) {
        ffield = true;
      } else if ((strcmp(ffield_arg, "no") == 0) || (strcmp(ffield_arg, "off") == 0)) {
        ffield = false;
      } else {
        error->all(FLERR, "Invalid argument after ffield keyword");
      }
    } else {
      error->all(FLERR, "Illegal fix electrode/conp command");
    }
    iarg++;
  }

  // union of all coupled groups
  std::string union_group = "conp_group";
  std::string group_cmd = union_group + " union";
  for (int g : groups) {
    group_cmd += " ";
    group_cmd += group->names[g];
  }
  group->assign(group_cmd);
  igroup = group->find(union_group);
  if (igroup < 0) error->all(FLERR, "Failed to create union of groups");
  // construct computes
  ele_vector = new ElectrodeVector(lmp, igroup, igroup, eta, true);
  if (!(read_inv || read_mat)) { array_compute = new ElectrodeMatrix(lmp, igroup, eta); }

  // error checks
  assert(groups.size() == group_bits.size());
  assert(groups.size() == group_psi.size());
  assert(groups.size() == group_psi_var_styles.size());
  assert(groups.size() == group_psi_var_names.size());
  assert(igroup == ele_vector->igroup);
  if (!(read_mat || read_inv)) assert(igroup == array_compute->igroup);
  if (read_inv && read_mat) error->all(FLERR, "Cannot read matrix from two files");
  if (f_mat && read_inv)
    error->all(FLERR,
               "Cannot write elastance matrix if reading capacitance matrix "
               "from file");
  num_of_groups = static_cast<int>(groups.size());
  size_vector = num_of_groups;

  // check groups are consistent
  int *mask = atom->mask;
  int groups_overlap = 0;
  for (int i = 0; i < atom->nlocal; i++) {
    int m = mask[i];
    int matches = 0;
    for (int bit : group_bits)
      if (m & bit) matches++;
    if (matches > 1) {
      groups_overlap++;
    } else {
      assert(!matches == !(m & group->bitmask[igroup]));
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &groups_overlap, 1, MPI_INT, MPI_SUM, world);
  if (groups_overlap) error->all(FLERR, "Groups may not overlap");
  groupbit = group->bitmask[igroup];
  ngroup = group->count(igroup);

  accel_interface = new ElectrodeAccelInterface(lmp);

  memory->create(iele_gathered, ngroup, "FixElectrode:iele_gathered");
  memory->create(buf_gathered, ngroup, "FixElectrode:buf_gathered");
  memory->create(potential_iele, ngroup, "FixElectrode:potential_iele");
  memory->create(charge_iele, ngroup, "FixElectrode:charge_iele");

  atom->add_callback(Atom::GROW);    // atomvec track local electrode atoms
  comm_reverse = 1;

  nlocalele = 0;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

//                     0   1   2       3 return 4 if successful
// fix_modify fxupdate set v   [group] [var]
// fix_modify fxupdate set qsb [group] [var]
//                     0   1   2        3        4 return 5 if successful
// fix_modify fxupdate set mc  [group1] [group2] [var]
// fix_modify fxupdate set me  [group1] [group2] [var]
int FixElectrodeConp::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "tf") == 0) {
    if (narg < 4) error->all(FLERR, fmt::format("Incorrect number of arguments for fix_modify tf"));
    tfflag = true;
    // read atom type, Thomas-Fermi length, and voronoi volume (reciprocal
    // number density)
    int const type = utils::inumeric(FLERR, arg[1], false, lmp);
    double const len = utils::numeric(FLERR, arg[2], false, lmp);
    double const voronoi = utils::numeric(FLERR, arg[3], false, lmp);
    // check type exists and is completely in electrode
    int not_in_ele = 0;
    int in_ele = 0;
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->type[i] == type) {
        if (atom->mask[i] & groupbit)
          in_ele++;
        else
          not_in_ele++;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &in_ele, 1, MPI_INT, MPI_SUM, world);
    if (in_ele == 0) error->all(FLERR, "No atoms of type in electrode");
    MPI_Allreduce(MPI_IN_PLACE, &not_in_ele, 1, MPI_INT, MPI_SUM, world);
    if (not_in_ele)
      error->warning(FLERR,
                     "Not all atoms of type in electrode; Thomas-Fermi parameters will be ignored "
                     "for electrolyte");
    // insert into map, replace if already exists
    auto entry = tf_types.find(type);
    if (entry != end(tf_types)) tf_types.erase(entry);
    tf_types.insert(std::pair<int, double>(type, MY_4PI * len * len / voronoi));
    return 4;

  } else if (strcmp(arg[0], "timer") == 0) {
    if (narg < 2)
      error->all(FLERR, fmt::format("Incorrect number of arguments for fix_modify timer"));
    timer_flag = utils::logical(FLERR, arg[1], false, lmp);
    return 2;

  } else if (strcmp(arg[0], "set") == 0) {
    if (strcmp(arg[1], "v") == 0 || strcmp(arg[1], "qsb") == 0) {
      if (narg != 4)
        error->all(FLERR,
                   fmt::format("Incorrect number of arguments for fix_modify set {}", arg[1]));
      int outputgrp = groupnum_from_name(arg[2]);
      int outputvar = input->variable->find(arg[3]);
      if (outputvar < 0)
        error->all(FLERR,
                   fmt::format("Variable {} for fix_modify electrode does not exist", arg[3]));
      if (!input->variable->internalstyle(outputvar))
        error->all(FLERR,
                   fmt::format("Variable {} for fix_modify electrode is invalid style", arg[3]));
      if (strcmp(arg[1], "v") == 0)
        setvars_types.push_back(V);
      else
        setvars_types.push_back(QSB);
      setvars_groups.push_back(outputgrp);
      setvars_vars.push_back(outputvar);
      return 4;
    } else if (strcmp(arg[1], "mc") == 0 || strcmp(arg[1], "me") == 0) {
      if (narg != 5)
        error->all(FLERR,
                   fmt::format("Incorrect number of arguments for fix_modify set {}", arg[1]));
      int outputgrpi = groupnum_from_name(arg[2]);
      int outputgrpj = groupnum_from_name(arg[3]);
      int outputgrpij = outputgrpi * num_of_groups + outputgrpj;
      int outputvar = input->variable->find(arg[4]);
      if (outputvar < 0)
        error->all(FLERR,
                   fmt::format("Variable {} for fix_modify electrode does not exist", arg[4]));
      if (!input->variable->internalstyle(outputvar))
        error->all(FLERR,
                   fmt::format("Variable {} for fix_modify electrode is invalid style", arg[4]));
      if (strcmp(arg[1], "mc") == 0)
        setvars_types.push_back(MC);
      else
        setvars_types.push_back(ME);
      setvars_groups.push_back(outputgrpij);
      setvars_vars.push_back(outputvar);
      return 5;
    } else
      error->all(FLERR, "Invalid set option for fix_modify electrode");
  } else
    error->all(FLERR, "Invalid argument for fix_modify electrode");
  return 0;
}

/* ---------------------------------------------------------------------- */

int FixElectrodeConp::modify_param(const std::string &param_str)
{
  auto args = utils::split_words(param_str);
  char **newarg = new char *[args.size()];
  int i = 0;
  for (const auto &arg : args) { newarg[i++] = (char *) arg.c_str(); }
  int tmp = modify_param(args.size(), newarg);
  delete[] newarg;
  return tmp;
}

/* ---------------------------------------------------------------------- */

int FixElectrodeConp::groupnum_from_name(char *groupname)
{
  int id = group->find(groupname);
  if (id < 0) error->all(FLERR, fmt::format("Group {} does not exist", groupname));
  for (int g = 0; g < num_of_groups; g++) {
    if (groups[g] == id) return g;
  }
  error->all(FLERR, fmt::format("Group {} is not coupled by fix electrode", groupname));
  return -1;    // dummy return value
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::init()
{
  pair = nullptr;    // not sure if needed -- remove if unnecessary
  pair = (Pair *) force->pair_match("coul", 0);
  if (pair == nullptr) {    // couldn't find a pair with name coul -- maybe hybrid
    // return 1st hybrid substyle containing 'coul'
    pair = (Pair *) force->pair_match("coul", 0, 1);
  }
  if (pair == nullptr) error->all(FLERR, "Fix electrode couldn't find a Coulombic pair style");

  // check for package intel
  accel_interface->intel_find_fix();
  if (etypes_neighlists)
    request_etypes_neighlists();
  else {
    auto Req = neighbor->add_request(this);
    if (intelflag) Req->enable_intel();
  }
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::init_list(int id, NeighList *ptr)
{
  if (etypes_neighlists) {
    if (id == 1)
      mat_neighlist = ptr;
    else if (id == 2)
      vec_neighlist = ptr;
  } else
    mat_neighlist = vec_neighlist = ptr;
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::post_constructor()
{
  if (!ffield) return;
  // ffield: test conditions and set up efield
  if (num_of_groups != 2) error->all(FLERR, "Number of electrodes must be two with ffield yes");
  if (!symm) error->all(FLERR, "Keyword symm off not allowed with ffield yes");
  if (domain->zperiodic == 0 || domain->boundary[2][0] != 0 || domain->boundary[2][1] != 0)
    error->all(FLERR, "Periodic z boundaries required with ffield yes");

  top_group = get_top_group();
  // assign variable names:
  std::string var_vtop = fixname + "_ffield_vtop";
  std::string var_vbot = fixname + "_ffield_vbot";
  std::string var_efield = fixname + "_ffield_zfield";
  // set variables:
  input->variable->set(fmt::format("{} equal f_{}[{}]", var_vbot, fixname, 1 + 1 - top_group));
  input->variable->set(fmt::format("{} equal f_{}[{}]", var_vtop, fixname, 1 + top_group));
  input->variable->set(fmt::format("{} equal (v_{}-v_{})/lz", var_efield, var_vbot, var_vtop));
  // check for other efields and warn if found
  if (modify->get_fix_by_style("efield").size() > 0)
    error->warning(FLERR, "Other efield fixes found -- please make sure this is intended!");
  // call fix command:
  // fix [varstem]_efield all efield 0.0 0.0 [var_vdiff]/lz
  std::string efield_call = fixname + "_efield all efield 0.0 0.0 v_" + var_efield;
  modify->add_fix(efield_call, 1);
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::setup_post_neighbor()
{
  int const nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *tag = atom->tag;

  // setvars asserts:
  assert(setvars_groups.size() == setvars_vars.size());
  assert(setvars_groups.size() == setvars_types.size());

  // if Thomas-Fermi, make sure all electrode atoms have parameters
  if (tfflag) {
    int unset_tf = 0;
    int *type = atom->type;
    for (int i = 0; i < nlocal; i++) {
      if ((groupbit & mask[i]) && (tf_types.count(type[i]) == 0)) unset_tf++;
    }
    MPI_Allreduce(MPI_IN_PLACE, &unset_tf, 1, MPI_INT, MPI_SUM, world);
    if (unset_tf)
      error->all(FLERR, fmt::format("Thomas-Fermi parameters not set for all types in electrode"));
  }

  // get equal-style variable ids:
  group_psi_var_ids = std::vector<int>(num_of_groups, -1);
  for (int g = 0; g < num_of_groups; g++) {
    if (group_psi_var_styles[g] == CONST) continue;
    const char *var_name = group_psi_var_names[g].c_str();
    int var_id = input->variable->find(var_name);
    if (var_id < 0)
      error->all(FLERR, fmt::format("Variable '{}' for fix electrode does not exist", var_name));
    if (!input->variable->equalstyle(var_id))
      error->all(FLERR,
                 fmt::format("Variable '{}' for fix electrode is not equal-style", var_name));
    group_psi_var_ids[g] = var_id;
  }

  // pair and list setups:

  evscale = force->qe2f / force->qqrd2e;
  ele_vector->setup(pair, vec_neighlist, timer_flag);
  if (!(read_mat || read_inv)) {
    if (etypes_neighlists) neighbor->build_one(mat_neighlist, 0);
    array_compute->setup(tag_to_iele, pair, mat_neighlist);
    if (tfflag) { array_compute->setup_tf(tf_types); }
  }
  // setup psi with target potentials
  iele_to_group = std::vector<int>(ngroup, -1);
  sd_vectors = std::vector<std::vector<double>>(num_of_groups, std::vector<double>(ngroup));
  sb_charges = std::vector<double>(num_of_groups);
  for (int i = 0; i < nlocal; i++) {
    for (int g = 0; g < num_of_groups; g++) {
      if (mask[i] & group_bits[g]) { iele_to_group[tag_to_iele[tag[i]]] = g; }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &iele_to_group.front(), ngroup, MPI_INT, MPI_MAX, world);

  auto const order_matrix = [](std::vector<tagint> order, double **mat) {
    size_t n = order.size();
    std::vector<std::vector<double>> ordered_mat(n, std::vector<double>(n));
    for (size_t i = 0; i < n; i++) {
      bigint const gi = order[i];
      for (size_t j = 0; j < n; j++) { ordered_mat[gi][order[j]] = mat[i][j]; }
    }
    return ordered_mat;
  };

  // capacitance matrix
  memory->create(capacitance, ngroup, ngroup, "fix_electrode:capacitance");
  if (read_inv) {
    read_from_file(input_file_inv, capacitance, "capacitance");
  } else {
    // temporarily hold elastance in "capacitance"
    if (read_mat)
      read_from_file(input_file_mat, capacitance, "elastance");
    else
      array_compute->compute_array(capacitance, timer_flag);
    if (f_mat && !(read_inv))
      write_to_file(f_mat, taglist_bygroup, order_matrix(group_idx, capacitance));
    invert();    // TODO  uncommented lots of stuff here
  }
  if (symm) symmetrize();

  // build sd vectors and macro matrices
  MPI_Barrier(world);
  double start = MPI_Wtime();
  if (ffield) {
    compute_sd_vectors_ffield();
  } else {
    compute_sd_vectors();
  }
  compute_macro_matrices();
  MPI_Barrier(world);
  if (timer_flag && (comm->me == 0))
    utils::logmesg(
        lmp, fmt::format("SD-vector and macro matrices time: {:.4g} s\n", MPI_Wtime() - start));

  // initial charges and b vector
  update_charges();

  // write to files, ordered by group
  memset(potential_i, 0, atom->nmax * sizeof(double));
  ele_vector->compute_vector(potential_i);
  if (force->newton_pair) comm->reverse_comm(this);
  buffer_and_gather(potential_i, potential_iele);
  if (comm->me == 0) {
    if (f_vec) {
      std::vector<std::vector<double>> vec(ngroup, std::vector<double>(1));
      for (int i = 0; i < ngroup; i++) { vec[group_idx[i]][0] = potential_iele[i]; }
      write_to_file(f_vec, taglist_bygroup, vec);
    }
    if (f_inv) { write_to_file(f_inv, taglist_bygroup, order_matrix(group_idx, capacitance)); }
  }
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::setup_pre_reverse(int eflag, int /*vflag*/)
{
  // correct forces for initial timestep
  gausscorr(eflag, true);
  self_energy(eflag);
  potential_energy(eflag);
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::invert()
{
  MPI_Barrier(world);
  double invert_time = MPI_Wtime();
  if (timer_flag && (comm->me == 0)) utils::logmesg(lmp, "CONP inverting matrix\n");
  int m = ngroup, n = ngroup, lda = ngroup;
  std::vector<int> ipiv(ngroup);
  int const lwork = ngroup * ngroup;
  std::vector<double> work(lwork);

  int info_rf, info_ri;
  dgetrf_(&m, &n, &capacitance[0][0], &lda, &ipiv.front(), &info_rf);
  dgetri_(&n, &capacitance[0][0], &lda, &ipiv.front(), &work.front(), &lwork, &info_ri);
  if (info_rf != 0 || info_ri != 0) error->all(FLERR, "CONP matrix inversion failed!");
  MPI_Barrier(world);
  if (timer_flag && (comm->me == 0))
    utils::logmesg(lmp, fmt::format("Invert time: {:.4g} s\n", MPI_Wtime() - invert_time));
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::symmetrize()
{
  // S matrix to enforce charge neutrality constraint
  std::vector<double> AinvE(ngroup, 0.);
  double EAinvE = 0.0;
  for (int i = 0; i < ngroup; i++) {
    double AinvEtmp = 0.0;
    for (int j = 0; j < ngroup; j++) { AinvEtmp += capacitance[i][j]; }
    AinvE[i] = AinvEtmp;    // use temp accumulator to enable vectorization
    EAinvE += AinvE[i];
  }
  for (int i = 0; i < ngroup; i++) {
    double iAinvE = AinvE[i];
    for (int j = 0; j < ngroup; j++) { capacitance[i][j] -= AinvE[j] * iAinvE / EAinvE; }
  }
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::setup_pre_exchange()    // create_taglist
{
  int *mask = atom->mask;
  int const nlocal = atom->nlocal;
  int const nprocs = comm->nprocs;
  tagint *tag = atom->tag;

  recvcounts = new int[nprocs];
  displs = new int[nprocs];

  // assign a tag to each matrix index sorted by group and by tag
  taglist_bygroup = std::vector<tagint>();
  nlocalele = 0;
  for (int gbit : group_bits) {
    std::vector<tagint> taglist_local;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & gbit) {
        taglist_local.push_back(tag[i]);
        nlocalele++;
      }
    }
    // gather from all cpus for this group
    int gnum_local = taglist_local.size();
    MPI_Allgather(&gnum_local, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
    displs[0] = 0;
    for (int i = 1; i < nprocs; i++) { displs[i] = displs[i - 1] + recvcounts[i - 1]; }
    int const gnum = displs[nprocs - 1] + recvcounts[nprocs - 1];
    std::vector<tagint> taglist_all(gnum);
    MPI_Allgatherv(&taglist_local.front(), gnum_local, MPI_LMP_TAGINT, &taglist_all.front(),
                   recvcounts, displs, MPI_LMP_TAGINT, world);
    std::sort(taglist_all.begin(), taglist_all.end());
    for (tagint t : taglist_all) taglist_bygroup.push_back(t);
  }

  // taglist only sorted by tag not group, same order as in computes
  taglist = taglist_bygroup;
  std::sort(taglist.begin(), taglist.end());

  tag_to_iele = std::unordered_map<tagint, int>();
  tag_to_iele.reserve(taglist.size());
  for (size_t i = 0; i < taglist.size(); i++) {
    tag_to_iele.insert(std::pair<tagint, int>(taglist[i], i));
  }

  nlocalele_outdated = 1;    // force regather

  // group_idx allows mapping a vector that is sorted by taglist to being
  // ordered by taglist_bygroup
  group_idx = std::vector<tagint>(taglist_bygroup.size());
  for (std::size_t i{0}; i < taglist_bygroup.size(); i++) {
    group_idx[i] = (tagint) tag_to_iele[taglist_bygroup[i]];
  }

  // if memory_usage > 0.5 GiB, warn with expected usage
  double mem_needed = memory_usage();
  mem_needed /= (1024 * 1024 * 1024);    // convert to GiB
  if (mem_needed > 0.5)
    error->warning(FLERR,
                   fmt::format("Please ensure there is sufficient memory for fix electrode "
                               "(anticipated usage is {:.1f} per proc)",
                               mem_needed));
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::pre_force(int)
{
  update_charges();
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::pre_reverse(int eflag, int /*vflag*/)
{
  gausscorr(eflag, true);
  self_energy(eflag);
  potential_energy(eflag);
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::compute_sd_vectors()
{
  for (int g = 0; g < num_of_groups; g++) {
    for (int j = 0; j < ngroup; j++) {
      if (iele_to_group[j] == g) {
        for (int k = 0; k < ngroup; k++) { sd_vectors[g][k] += capacitance[k][j] * evscale; }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::compute_sd_vectors_ffield()
{
  double **x = atom->x;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  double zprd = domain->prd[2];
  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      int const i_iele = tag_to_iele[tag[i]];
      double const zprd_offset = (mask[i] & group_bits[top_group]) ? 0.0 : 1.0;
      double const evscale_elez = evscale * (x[i][2] / zprd + zprd_offset);
      for (int g = 0; g < num_of_groups; g++) {
        double gmult = (g == top_group) ? -1.0 : 1.0;
        for (int k = 0; k < ngroup; k++) {
          sd_vectors[g][k] += gmult * capacitance[k][i_iele] * evscale_elez;
        }
      }
    }
  }
  for (int g = 0; g < num_of_groups; g++) {
    MPI_Allreduce(MPI_IN_PLACE, &sd_vectors[g].front(), ngroup, MPI_DOUBLE, MPI_SUM, world);
  }
}

/* ---------------------------------------------------------------------- */

int FixElectrodeConp::get_top_group()
{
  double *zmax = new double[num_of_groups];
  double **x = atom->x;
  for (int g = 0; g < num_of_groups; g++) { zmax[g] = domain->boxlo[2]; }
  int *mask = atom->mask;
  for (int i = 0; i < atom->nlocal; i++) {
    for (int g = 0; g < num_of_groups; g++) {
      if (mask[i] & group_bits[g]) {
        if (x[i][2] > zmax[g]) zmax[g] = x[i][2];
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, zmax, num_of_groups, MPI_DOUBLE, MPI_MAX, world);
  int gmax = 0;
  for (int g = 0; g < num_of_groups; g++) { gmax = (zmax[g] > zmax[gmax]) ? g : gmax; }
  delete[] zmax;
  return gmax;
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::update_charges()
{
  MPI_Barrier(world);
  double start = MPI_Wtime();
  if (atom->nmax > nmax) {
    memory->destroy(potential_i);
    nmax = atom->nmax;
    memory->create(potential_i, nmax, "FixElectrode:potential_i");
  }

  std::fill(sb_charges.begin(), sb_charges.end(), 0.);
  double *q = atom->q;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int const nlocal = atom->nlocal;
  int const nall = nlocal + atom->nghost;
  memset(potential_i, 0, atom->nmax * sizeof(double));
  ele_vector->compute_vector(potential_i);
  if (force->newton_pair) comm->reverse_comm(this);
  buffer_and_gather(potential_i, potential_iele);
  MPI_Barrier(world);
  pre_update();
  MPI_Barrier(world);
  double mult_start = MPI_Wtime();
  for (int i_iele = 0; i_iele < nlocalele; i_iele++) {
    double q_tmp = 0;
    int const iele = list_iele[i_iele];
    double *_noalias caprow = capacitance[iele];
    for (int j = 0; j < ngroup; j++) { q_tmp -= caprow[j] * potential_iele[j]; }
    buf_iele[i_iele] = q_tmp;
    sb_charges[iele_to_group[iele]] += q_tmp;
  }
  gather_elevec(charge_iele);
  MPI_Allreduce(MPI_IN_PLACE, &sb_charges.front(), num_of_groups, MPI_DOUBLE, MPI_SUM, world);

  update_setvars(QSB);
  update_psi();         // use for equal-style and conq
  update_setvars(V);    // push psi into 'fix_modify set' vars

  for (int g = 0; g < num_of_groups; g++)
    for (int j = 0; j < ngroup; j++) { charge_iele[j] += sd_vectors[g][j] * group_psi[g]; }

  for (int i = 0; i < nall; i++) {
    if (!(groupbit & mask[i])) continue;
    int const iele = tag_to_iele[tag[i]];
    q[i] = charge_iele[iele];
  }

  MPI_Barrier(world);
  mult_time += MPI_Wtime() - mult_start;
  update_time += MPI_Wtime() - start;
  accel_interface->intel_pack_buffers();
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::update_psi()
{
  for (int g = 0; g < num_of_groups; g++) {
    if (group_psi_var_styles[g] == CONST) continue;
    group_psi[g] = input->variable->compute_equal(group_psi_var_ids[g]);
  }
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::update_setvars(int vtype)
{
  int num_of_vars = static_cast<int>(setvars_groups.size());
  for (int v = 0; v < num_of_vars; v++) {
    if (vtype != setvars_types[v]) continue;
    switch (vtype) {
      case V:
        input->variable->internal_set(setvars_vars[v], group_psi[setvars_groups[v]]);
        break;
      case QSB:
        input->variable->internal_set(setvars_vars[v], sb_charges[setvars_groups[v]]);
        break;
      case MC: {
        int groupi = setvars_groups[v] / num_of_groups;
        int groupj = setvars_groups[v] % num_of_groups;
        input->variable->internal_set(setvars_vars[v], macro_capacitance[groupi][groupj]);
      } break;
      case ME: {
        int groupi = setvars_groups[v] / num_of_groups;
        int groupj = setvars_groups[v] % num_of_groups;
        input->variable->internal_set(setvars_vars[v], macro_elastance[groupi][groupj]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::compute_macro_matrices()
{
  macro_capacitance =
      std::vector<std::vector<double>>(num_of_groups, std::vector<double>(num_of_groups));
  for (int g = 0; g < num_of_groups; g++) {
    for (int k = 0; k < ngroup; k++) { macro_capacitance[iele_to_group[k]][g] += sd_vectors[g][k]; }
  }

  if (symm) {
    // scaling with C[0][0] improves numerical stability
    double scalar = macro_capacitance[0][0];
    macro_capacitance.back() = std::vector<double>(num_of_groups, scalar);
  }

  macro_elastance =
      std::vector<std::vector<double>>(num_of_groups, std::vector<double>(num_of_groups));

  // TODO: does determinant ever = 0 if !symm? what to do then?
  if (num_of_groups == 1) {
    macro_elastance[0][0] = 1 / macro_capacitance[0][0];
  } else if (num_of_groups == 2) {
    double detinv = 1 /
        (macro_capacitance[0][0] * macro_capacitance[1][1] -
         macro_capacitance[0][1] * macro_capacitance[1][0]);
    macro_elastance[0][0] = macro_capacitance[1][1] * detinv;
    macro_elastance[1][1] = macro_capacitance[0][0] * detinv;
    macro_elastance[0][1] = -macro_capacitance[0][1] * detinv;
    macro_elastance[1][0] = -macro_capacitance[1][0] * detinv;
  } else {
    int m = num_of_groups;
    int n = m, lda = m;
    std::vector<int> ipiv(m);
    int const lwork = m * m;
    std::vector<double> work(lwork);
    std::vector<double> tmp(lwork);

    for (int i = 0; i < num_of_groups; i++) {
      for (int j = 0; j < num_of_groups; j++) {
        int idx = i * num_of_groups + j;
        tmp[idx] = macro_capacitance[i][j];
      }
    }

    int info_rf, info_ri;
    dgetrf_(&m, &n, &tmp.front(), &lda, &ipiv.front(), &info_rf);
    dgetri_(&n, &tmp.front(), &lda, &ipiv.front(), &work.front(), &lwork, &info_ri);
    if (info_rf != 0 || info_ri != 0) error->all(FLERR, "CONP macro matrix inversion failed!");
    for (int i = 0; i < num_of_groups; i++) {
      for (int j = 0; j < num_of_groups; j++) {
        int idx = i * num_of_groups + j;
        macro_elastance[i][j] = tmp[idx];
      }
    }
  }

  update_setvars(MC);
  update_setvars(ME);
}

/* ---------------------------------------------------------------------- */

double FixElectrodeConp::compute_scalar()
{
  return potential_energy(0);
}
/* ---------------------------------------------------------------------- */

double FixElectrodeConp::compute_vector(int i)
{
  return group_psi[i];
}

/* ---------------------------------------------------------------------- */

double FixElectrodeConp::potential_energy(int eflag)
{
  // corrections to energy due to potential psi
  double const qqrd2e = force->qqrd2e;
  int const nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  double *q = atom->q;
  double energy = 0;
  for (int i = 0; i < nlocal; i++) {
    if (groupbit & mask[i]) {
      int const i_iele = tag_to_iele[tag[i]];
      double e = -qqrd2e * q[i] * group_psi[iele_to_group[i_iele]] * evscale;
      energy += e;
      if (eflag) {
        force->pair->ev_tally(i, i, nlocal, force->newton_pair, 0., e, 0, 0, 0,
                              0);    // 0 evdwl, 0 fpair, 0 delxyz
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPI_DOUBLE, MPI_SUM, world);
  return energy;
}
/* ---------------------------------------------------------------------- */

double FixElectrodeConp::self_energy(int eflag)
{
  // corrections to energy due to self interaction
  double const qqrd2e = force->qqrd2e;
  int const nlocal = atom->nlocal;
  double const pre = eta / sqrt(MY_2PI) * qqrd2e;
  int *mask = atom->mask;
  int *type = atom->type;
  double *q = atom->q;
  double energy = 0;
  for (int i = 0; i < nlocal; i++) {
    if (groupbit & mask[i]) {
      double const q2 = q[i] * q[i];
      double e = pre * q2;
      if (tfflag && (groupbit & mask[i])) e += 0.5 * qqrd2e * q2 * tf_types[type[i]];
      energy += e;
      if (eflag) {
        force->pair->ev_tally(i, i, nlocal, force->newton_pair, 0., e, 0, 0, 0,
                              0);    // 0 evdwl, 0 fpair, 0 delxyz
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPI_DOUBLE, MPI_SUM, world);
  return energy;
}

/* ---------------------------------------------------------------------- */

double FixElectrodeConp::gausscorr(int eflag, bool fflag)
{
  // correction to short range interaction due to eta

  int evflag = pair->evflag;
  double const qqrd2e = force->qqrd2e;
  int const nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int newton_pair = force->newton_pair;
  int inum = vec_neighlist->inum;
  int *ilist = vec_neighlist->ilist;
  int *numneigh = vec_neighlist->numneigh;
  int **firstneigh = vec_neighlist->firstneigh;
  double energy_sr = 0.;
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    bool i_in_ele = groupbit & mask[i];
    double qtmp = q[i];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];
    int itype = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int const j = jlist[jj] & NEIGHMASK;
      bool j_in_ele = groupbit & mask[j];
      if (!(i_in_ele || j_in_ele)) continue;
      double eta_ij = (i_in_ele && j_in_ele) ? eta / MY_SQRT2 : eta;

      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];
      double rsq = delx * delx + dely * dely + delz * delz;
      int jtype = type[j];

      if (rsq < force->pair->cutsq[itype][jtype]) {
        double r2inv = 1.0 / rsq;
        double r = sqrt(rsq);
        double erfc_etar = 0.;
        double derfcr = ElectrodeMath::safe_derfcr(eta_ij * r, erfc_etar);
        double prefactor = qqrd2e * qtmp * q[j] / r;
        energy_sr -= prefactor * erfc_etar;

        double fpair = prefactor * derfcr * r2inv;
        if (fflag) {
          f[i][0] += delx * fpair;
          f[i][1] += dely * fpair;
          f[i][2] += delz * fpair;
          if (newton_pair || j < nlocal) {
            f[j][0] -= delx * fpair;
            f[j][1] -= dely * fpair;
            f[j][2] -= delz * fpair;
          }
        }

        double ecoul = 0.;
        if (eflag) ecoul = -prefactor * erfc_etar;

        if (evflag) {
          force->pair->ev_tally(i, j, nlocal, newton_pair, 0., ecoul, fpair, delx, dely, delz);
        }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &energy_sr, 1, MPI_DOUBLE, MPI_SUM, world);
  return energy_sr;
}

/* ---------------------------------------------------------------------- */

FixElectrodeConp::~FixElectrodeConp()
{
  if (timer_flag && (comm->me == 0)) {
    try {
      utils::logmesg(lmp, fmt::format("Multiplication time: {:.4g} s\n", mult_time));
      utils::logmesg(lmp, fmt::format("Update time: {:.4g} s\n", update_time));
    } catch (std::exception &) {
    }
  }
  if (!modify->get_fix_by_id(id))             // avoid segfault if derived fixes' ctor throws err
    atom->delete_callback(id, Atom::GROW);    // atomvec track local electrode atoms
  delete[] recvcounts;
  delete[] displs;
  memory->destroy(iele_gathered);
  memory->destroy(buf_gathered);
  memory->destroy(potential_i);
  memory->destroy(potential_iele);
  memory->destroy(charge_iele);

  if (!(read_mat || read_inv)) delete array_compute;
  delete ele_vector;
  memory->destroy(capacitance);
  delete accel_interface;
  if (f_inv) fclose(f_inv);
  if (f_mat) fclose(f_mat);
  if (f_vec) fclose(f_vec);
}

/* ---------------------------------------------------------------------- */

int FixElectrodeConp::setmask()
{
  int mask = 0;
  mask |= FixConst::PRE_EXCHANGE;
  mask |= FixConst::POST_NEIGHBOR;
  mask |= FixConst::PRE_FORCE;
  mask |= FixConst::PRE_REVERSE;
  //mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::write_to_file(FILE *file, const std::vector<tagint> &tags,
                                     const std::vector<std::vector<double>> &mat)
{
  for (const auto &t : tags) fmt::print(file, "{:20}", t);
  fputs("\n", file);
  for (const auto &vec : mat) {
    for (const auto &x : vec) fmt::print(file, "{:20.11e}", x);
    fputs("\n", file);
  }
}

/*----------------------------------------------------------------------- */

void FixElectrodeConp::read_from_file(const std::string &input_file, double **array,
                                      const std::string &filetype)
{
  if (comm->me == 0) {
    std::vector<std::vector<double>> matrix;
    std::vector<tagint> tags;
    try {
      TextFileReader reader(input_file, filetype);
      reader.set_bufsize(ngroup * 20 + 4);

      // get line with tags
      auto values = reader.next_values(ngroup);
      for (int i = 0; i < ngroup; ++i) tags.push_back(values.next_tagint());

      std::vector<double> a_line;
      for (int i = 0; i < ngroup; ++i) {
        a_line.clear();
        values = reader.next_values(ngroup);
        for (int j = 0; j < ngroup; ++j) a_line.push_back(values.next_double());
        matrix.push_back(a_line);
      }
    } catch (std::exception &e) {
      error->one(FLERR, "Error parsing {} file: {}", filetype, e.what());
    }

    std::vector<tagint> idx;
    for (const auto &t : taglist) {
      for (std::size_t i = 0; i < tags.size(); i++) {
        if (t == tags[i]) {
          idx.push_back(i);
          break;
        }
      }
    }
    if ((bigint) idx.size() != ngroup)
      error->all(FLERR, fmt::format("Read tags do not match taglist of electrode/conp"));
    for (bigint i = 0; i < ngroup; i++) {
      bigint const ii = idx[i];
      for (bigint j = 0; j < ngroup; j++) array[i][j] = matrix[ii][idx[j]];
    }
  }
  MPI_Bcast(&array[0][0], ngroup * ngroup, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::request_etypes_neighlists()
{
  int const ntypes = atom->ntypes;
  int *iskip_mat = new int[ntypes + 1];
  int *iskip_vec = new int[ntypes + 1];
  int **ijskip_mat;
  memory->create(ijskip_mat, ntypes + 1, ntypes + 1, "fixelectrode:ijskip_mat");
  int **ijskip_vec;
  memory->create(ijskip_vec, ntypes + 1, ntypes + 1, "fixelectrode:ijskip_vec");
  for (int itype = 1; itype <= ntypes; ++itype) {
    // itype is 1-indexed -- follow LAMMPS convention
    iskip_mat[itype] = 1;    // alist skips all except etypes by default
    iskip_vec[itype] = 0;
    for (int jtype = 1; jtype <= ntypes; ++jtype) { ijskip_mat[itype][jtype] = 1; }
  }
  for (int etype : etypes) {
    iskip_mat[etype] = 0;
    ijskip_mat[etype][etype] = 0;
  }
  // now, iskip_mat[itype] == 0 iff etype
  // set ijskip_vec[itype][jtype] == 0 if (i is etype XOR j is etype)
  for (int itype = 1; itype <= ntypes; ++itype) {
    for (int jtype = 1; jtype <= ntypes; ++jtype) {
      bool ele_and_sol = (iskip_mat[itype] != iskip_mat[jtype]);
      ijskip_vec[itype][jtype] = (ele_and_sol) ? 0 : 1;
    }
  }

  if (!(read_inv || read_mat)) {
    auto matReq = neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
    matReq->set_skip(iskip_mat, ijskip_mat);
    matReq->set_id(1);
    if (intelflag) matReq->enable_intel();
  } else {
    delete[] iskip_mat;
    memory->destroy(ijskip_mat);
  }

  auto vecReq = neighbor->add_request(this);
  vecReq->set_skip(iskip_vec, ijskip_vec);
  vecReq->set_id(2);
  if (intelflag) vecReq->enable_intel();
}

int FixElectrodeConp::pack_exchange(int i, double * /* buf */)
{
  if (atom->mask[i] & groupbit) {
    nlocalele_outdated = 1;
    nlocalele--;    // decrement nlocalele if we are packing away a particle
  }
  return 0;
}

int FixElectrodeConp::unpack_exchange(int nlocal, double * /* buf */)
{
  if (atom->mask[nlocal] & groupbit) {    // this should work
    nlocalele_outdated = 1;
    nlocalele++;    // increment nlocalele if we are unpacking a particle
  }
  return 0;
}

void FixElectrodeConp::gather_list_iele()
{
  MPI_Allreduce(MPI_IN_PLACE, &nlocalele_outdated, 1, MPI_INT, MPI_SUM, world);
  if (nlocalele_outdated == 0) return;

  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int const nlocal = atom->nlocal;
  list_iele.clear();
  list_iele.reserve(nlocalele);
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) list_iele.push_back(tag_to_iele[tag[i]]);
  }

  nlocalele = static_cast<int>(list_iele.size());    // just for safety

  MPI_Allgather(&nlocalele, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
  displs[0] = 0;
  int const nprocs = comm->nprocs;
  for (int i = 1; i < nprocs; i++) { displs[i] = displs[i - 1] + recvcounts[i - 1]; }

  MPI_Allgatherv(&list_iele[0], nlocalele, MPI_INT, iele_gathered, recvcounts, displs, MPI_INT,
                 world);
  nlocalele_outdated = 0;
}

void FixElectrodeConp::gather_elevec(double *elevec)
{
  MPI_Allgatherv(&buf_iele[0], nlocalele, MPI_DOUBLE, buf_gathered, recvcounts, displs, MPI_DOUBLE,
                 world);

  for (int i = 0; i < ngroup; i++) { elevec[iele_gathered[i]] = buf_gathered[i]; }
}

void FixElectrodeConp::buffer_and_gather(double *ivec, double *elevec)
{
  gather_list_iele();

  buf_iele.reserve(nlocalele);    // avoid unexpected reallocs
  for (int i_iele = 0; i_iele < nlocalele; i_iele++) {
    buf_iele[i_iele] = ivec[atom->map(taglist[list_iele[i_iele]])];
  }

  gather_elevec(elevec);
}

double FixElectrodeConp::memory_usage()
{
  int const nprocs = comm->nprocs;
  int const nmax = atom->nmax;
  double bytes = 0.0;
  bytes += nmax * (sizeof(double));    // potential_i
  bytes += ngroup *
      (sizeof(int) + 3 * sizeof(double));       // iele_gathered, buf_gathered, pot / charge_iele
  bytes += ngroup * ngroup * sizeof(double);    // capacitance
  bytes += nprocs * (2 * sizeof(int));          // displs, recvcounts
  bytes += list_iele.capacity() * sizeof(int);
  bytes += buf_iele.capacity() * sizeof(double);

  return bytes;
}

int FixElectrodeConp::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) { buf[m++] = potential_i[i]; }

  return m;
}

void FixElectrodeConp::unpack_reverse_comm(int n, int *list, double *buf)
{
  for (int i = 0; i < n; i++) { potential_i[list[i]] += buf[i]; }
}
