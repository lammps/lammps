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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (GU), Kamila Savvidi (TUHH), Robert Meissner (Hereon, TUHH)
------------------------------------------------------------------------- */

#include "fix_electrode_conp.h"

#include "atom.h"
#include "charge_solver.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "electrode_cg.h"
#include "electrode_inv.h"
#include "electrode_mat_cg.h"
#include "electrode_math.h"
#include "electrode_matrix.h"
#include "electrode_taglist.h"
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
#include "variable.h"

#include <cassert>
#include <cmath>
#include <cstring>
#include <exception>
#include <memory>
#include <utility>

using namespace LAMMPS_NS;
using namespace MathConst;

static const char cite_fix_electrode[] =
    "fix electrode command: https://doi.org/10.1063/5.0099239\n\n"
    "@article{Ahrens2022\n"
    "author = {Ahrens-Iwers, Ludwig J.V. and Janssen, Mahijs and Tee, Shern R. and Mei{\\ss}ner, "
    "Robert H.},\n"
    "doi = {10.1063/5.0099239},\n"
    "title = {{ELECTRODE: An electrochemistry package for LAMMPS}},\n"
    "journal = {The Journal of Chemical Physics},\n"
    "year = {2022}\n"
    "volume = {157},\n"
    "pages = {084801},\n"
    "}\n";

//     0        1      2              3
// fix fxupdate group1 electrode/conp pot1 couple group2 pot2
FixElectrodeConp::FixElectrodeConp(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), charge_solver(nullptr), potential_i(nullptr), elyt_vector(nullptr),
    elec_vector(nullptr), matrix(nullptr), pair(nullptr), mat_neighlist(nullptr),
    vec_neighlist(nullptr), electrode_taglist(nullptr)

{
  if (lmp->citeme) lmp->citeme->add(cite_fix_electrode);
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Fix {} requires an atom map, see atom_modify", style);

  // fix.h output flags
  scalar_flag = 1;
  vector_flag = 1;
  extscalar = 1;
  extvector = 0;
  extarray = 0;

  virial_global_flag = 1;    // use virials of this fix
  thermo_virial = 1;         // set vflags for v_tally

  bool default_algo = true;
  algo = Algo::MATRIX_INV;
  matrix_algo = true;
  cg_threshold = 0.;
  write_inv = write_mat = write_vec = read_inv = read_mat = false;
  bool symm = false;
  ffield = false;
  taglist_constructed = false;
  thermo_time = 0.;

  top_group = 0;
  intelflag = false;
  tfflag = false;
  etapropflag = enflag = hardnessflag = false;
  pairflag = false;
  timer_flag = false;

  update_time = 0;
  mult_time = 0;

  qtotal = 0.;
  qtotal_var_style = VarStyle::UNSET;

  // read fix command
  fixname = std::string(arg[0]);
  groups = std::vector<int>(1, igroup);
  group_bits = std::vector<int>(1, groupbit);
  group_psi_var_names = std::vector<std::string>(1);
  group_psi_var_styles = std::vector<VarStyle>(1, VarStyle::CONST);
  group_psi_const = std::vector<double>(1);
  etypes_neighlists = false;
  if (strstr(arg[3], "v_") == arg[3]) {
    std::string vname = arg[3];
    group_psi_var_names[0] = vname.substr(2);
    group_psi_var_styles[0] = VarStyle::EQUAL;
  } else
    group_psi_const[0] = utils::numeric(FLERR, arg[3], false, lmp);
  bool etaflag = false;
  bool deprecated_single_eta = false;
  int iarg = 4;
  while (iarg < narg) {
    if ((strcmp(arg[iarg], "couple") == 0)) {
      if (iarg + 3 > narg) error->all(FLERR, "Need two arguments after couple keyword");
      int id = group->find(arg[++iarg]);
      if (id < 0) error->all(FLERR, "Group {} does not exist", arg[iarg]);
      groups.push_back(id);
      group_bits.push_back(group->bitmask[id]);
      ++iarg;
      if (strstr(arg[iarg], "v_") == arg[iarg]) {
        std::string vname = arg[iarg];
        group_psi_var_names.push_back(vname.substr(2));
        group_psi_var_styles.push_back(VarStyle::EQUAL);
        group_psi_const.push_back(0.);
      } else {
        std::string null;
        group_psi_var_names.push_back(null);
        group_psi_var_styles.push_back(VarStyle::CONST);
        group_psi_const.push_back(utils::numeric(FLERR, arg[iarg], false, lmp));
      }
    } else if ((strcmp(arg[iarg], "algo") == 0)) {
      if (!default_algo) error->one(FLERR, "Algorithm can be set only once");
      default_algo = false;
      if (iarg + 2 > narg) error->all(FLERR, "Need at least one argument after algo command");
      char *algo_arg = arg[++iarg];
      bool cg_algo = false;
      if ((strcmp(algo_arg, "mat_inv") == 0)) {
        algo = Algo::MATRIX_INV;
        matrix_algo = true;
      } else if ((strcmp(algo_arg, "mat_cg") == 0)) {
        algo = Algo::MATRIX_CG;
        matrix_algo = true;
        cg_algo = true;
      } else if ((strcmp(algo_arg, "cg") == 0)) {
        algo = Algo::CG;
        matrix_algo = false;
        cg_algo = true;
      } else {
        error->all(FLERR, "Unknown algo keyword {}", algo_arg);
      }
      if (cg_algo) {
        if (iarg + 2 > narg) error->all(FLERR, "Need one argument after algo cg command");
        cg_threshold = utils::numeric(FLERR, arg[++iarg], false, lmp);
      }
    } else if ((strncmp(arg[iarg], "write", 5) == 0)) {
      if (iarg + 2 > narg) error->all(FLERR, "Need one argument after write command");
      if ((strcmp(arg[iarg], "write_inv") == 0)) {    // capacitance matrix
        write_inv = true;
        output_file_inv = arg[++iarg];
      } else if ((strcmp(arg[iarg], "write_mat") == 0)) {    // elastance matrix
        write_mat = true;
        output_file_mat = arg[++iarg];
      } else if ((strcmp(arg[iarg], "write_vec") == 0)) {    // b vector
        write_vec = true;
        output_file_vec = arg[++iarg];
      } else {
        error->all(FLERR, "Illegal fix {} command with write", style);
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
        error->all(FLERR, "Illegal fix {} command with read", style);
      }
    } else if ((strcmp(arg[iarg], "temp") == 0)) {
      if (iarg + 4 > narg) error->all(FLERR, "Need three arguments after temp command");
      if (!utils::strmatch(style, "electrode/thermo"))
        error->all(FLERR, "temp keyword not available for fix {}", style);
      thermo_temp = force->boltz / force->qe2f * utils::numeric(FLERR, arg[++iarg], false, lmp);
      thermo_time = utils::numeric(FLERR, arg[++iarg], false, lmp);
      thermo_init = utils::inumeric(FLERR, arg[++iarg], false, lmp);
    } else if ((strcmp(arg[iarg], "qtotal") == 0)) {
      if (iarg + 2 > narg) error->all(FLERR, "Need one argument after qtotal keyword");
      ++iarg;
      if (strstr(arg[iarg], "v_") == arg[iarg]) {
        std::string vname = arg[iarg];
        qtotal_var_name = vname.substr(2);
        qtotal_var_style = VarStyle::EQUAL;
      } else {
        qtotal = utils::numeric(FLERR, arg[iarg], false, lmp);
        qtotal_var_style = VarStyle::CONST;
      }
    } else if ((strcmp(arg[iarg], "eta") == 0)) {
      if (etaflag || pairflag)
        error->all(FLERR,
                   "eta keyword cannot be used if eta or pair keyword has already been used");
      if (iarg + 2 > narg) error->all(FLERR, "Need one argument after eta command");
      etaflag = true;
      char *eta_str = arg[++iarg];
      int is_double, cols, ghost;
      if (utils::is_double(std::string(eta_str))) {
        eta = utils::numeric(FLERR, eta_str, false, lmp);
      } else if ((eta_index = atom->find_custom_ghost(eta_str + 2, is_double, cols, ghost)) != -1) {
        etapropflag = true;
        if (!is_double)
          error->all(FLERR, "eta keyword requires double-valued property/atom vector");
        if (cols != 0) error->all(FLERR, "eta keyword requires property/atom vector not an array");
        if (!ghost) error->all(FLERR, "eta keyword requires property/atom fix with ghost on");
      } else {
        error->all(FLERR,
                   "eta keyword requires a value or the name of previously defined property");
      }
    } else if ((strcmp(arg[iarg], "pair") == 0)) {
      if (etaflag || pairflag)
        error->all(FLERR,
                   "pair keyword cannot be used if eta or pair coammnd has already been used");
      if (iarg + 2 > narg) error->all(FLERR, "Need one argument after pair command");
      pairflag = true;
      pair_str = arg[++iarg];
    } else if ((strcmp(arg[iarg], "electronegativity") == 0) ||
               (strcmp(arg[iarg], "hardness") == 0)) {
      char *keyword = arg[iarg];
      if (iarg + 2 > narg) error->all(FLERR, "Need one argument after {} command", keyword);
      int is_double, cols;
      int index = atom->find_custom(arg[++iarg] + 2, is_double, cols);
      if (index == -1)
        error->all(FLERR, "{} keyword requires name of previously defined property", keyword);
      if (!is_double)
        error->all(FLERR, "{} keyword requires double-valued property/atom vector", keyword);
      if (cols != 0)
        error->all(FLERR, "{} keyword requires property/atom vector not an array", keyword);
      if (strcmp(keyword, "electronegativity") == 0) {
        en_index = index;
        enflag = true;
      } else if (strcmp(keyword, "hardness") == 0) {
        hardness_index = index;
        hardnessflag = true;
      }
    }
    // toggle parameters
    else if ((strcmp(arg[iarg], "etypes") == 0)) {
      etypes_neighlists = utils::logical(FLERR, arg[++iarg], false, lmp);
    } else if ((strncmp(arg[iarg], "symm", 4) == 0)) {
      symm = utils::logical(FLERR, arg[++iarg], false, lmp);
    } else if ((strcmp(arg[iarg], "ffield") == 0)) {
      ffield = utils::logical(FLERR, arg[++iarg], false, lmp);
    } else if (iarg == 4) {    // deprecated option to specify eta as fourth argument
      char *eta_str = arg[iarg];
      bool etanull = (strcmp(eta_str, "NULL") == 0);
      if (!etanull) {
        if (utils::is_double(std::string(eta_str))) {
          eta = utils::numeric(FLERR, eta_str, false, lmp);
          deprecated_single_eta = true;
        } else {
          error->all(FLERR, "Unknown keyword {} for fix {} command", arg[iarg], style);
        }
      }
      if (comm->me == 0)
        error->warning(FLERR,
                       "Setting eta as the fourth argument is deprecated and will be removed in "
                       "the future. Use the eta command.");
    } else {
      error->all(FLERR, "Unknown keyword {} for fix {} command", arg[iarg], style);
    }
    iarg++;
  }

  if (symm) {
    if (qtotal_var_style != VarStyle::UNSET) {
      error->all(FLERR, "{} cannot use qtotal keyword with symm on", this->style);
    }
    if (comm->me == 0) {
      if (ffield)
        error->warning(FLERR,
                       "The symm keyword is deprecated and will be removed in the future. "
                       "Symmetrization is automatically enabled when using ffield.");
      else
        error->warning(FLERR,
                       "The symm keyword is deprecated and will be removed in the future. Use "
                       "'qtotal 0' instead.");
    }
    qtotal_var_style = VarStyle::CONST;
    qtotal = 0.;
  }
  if (ffield) {
    if (algo != Algo::MATRIX_INV)
      error->all(FLERR, "ffield field is only implemented for matrix inversion");
    // TODO compatibility with qtotal?
  }
  if (!(etaflag || pairflag || deprecated_single_eta))
    error->all(FLERR, "The eta or pair keyword must be used");
  if (comm->me == 0) {
    if (deprecated_single_eta && (pairflag || etaflag))
      error->warning(FLERR,
                     "The eta parameter has been set as fourth argument but will be ignored "
                     "because the eta or pair command has been used");
  }
  if (pairflag && etypes_neighlists)
    error->all(FLERR, "The etypes and pair keyword are not compatible");

  // computatonal potential
  group_psi = std::vector<double>(groups.size());
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
  // construct computes and charge solver
  need_array_compute = !(read_inv || read_mat) && matrix_algo;
  need_elec_vector = algo == Algo::CG;
  // Might work with the plan to create "compute potential/atom"
  elyt_vector = new ElectrodeVector(lmp, 0, arg, igroup, igroup, eta, true);
  if (need_elec_vector) {
    elec_vector = new ElectrodeVector(lmp, 0, arg, igroup, igroup, eta, false);
  }
  assert(groups.size() == group_bits.size());
  assert(groups.size() == group_psi.size());
  assert(groups.size() == group_psi_const.size());
  assert(groups.size() == group_psi_var_styles.size());
  assert(groups.size() == group_psi_var_names.size());
  assert(igroup == elyt_vector->igroup);
  if (need_elec_vector) assert(igroup == elec_vector->igroup);
  if (algo != Algo::MATRIX_INV) {
    if (read_inv || write_inv)
      error->all(
          FLERR,
          "Selected algorithm does not use inverted matrix. Cannot read/write inverted matrix.");
  }
  if (!matrix_algo && (read_mat || write_mat || write_vec)) {
    error->all(FLERR,
               "Selected algorithm does not use matrix. Cannot read/write matrix or vector.");
  }
  if (read_inv && read_mat) error->all(FLERR, "Cannot read matrix from two files");
  if (write_mat && read_inv)
    error->all(FLERR, "Cannot write elastance matrix if reading capacitance matrix from file");
  num_of_groups = static_cast<int>(groups.size());
  size_vector = num_of_groups;
  array_flag = !!(algo == Algo::MATRIX_INV);
  if (array_flag) {
    size_array_rows = num_of_groups;
    size_array_cols = 2 + 2 * num_of_groups;
  }

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

  atom->add_callback(Atom::GROW);    // atomvec track local electrode atoms
  comm_forward = 1;

  nlocalele = 0;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

int FixElectrodeConp::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "tf") == 0) {
    if (narg < 4) error->all(FLERR, "Incorrect number of arguments for fix_modify {}", style);
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
    if (not_in_ele && (comm->me == 0))
      error->warning(FLERR,
                     "Not all atoms of type in electrode; Thomas-Fermi parameters will be ignored "
                     "for electrolyte");
    // insert into map, replace if already exists
    auto entry = tf_types.find(type);
    if (entry != end(tf_types)) tf_types.erase(entry);
    tf_types.insert(std::pair<int, double>(type, MY_4PI * len * len / voronoi));
    return 4;

  } else if (strcmp(arg[0], "timer") == 0) {
    if (narg < 2) error->all(FLERR, "Incorrect number of arguments for fix_modify {} timer", style);
    timer_flag = utils::logical(FLERR, arg[1], false, lmp);
    return 2;

  } else
    error->all(FLERR, "Unknown argument {} for fix_modify {}", arg[0], style);
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

void FixElectrodeConp::init()
{
  pair = nullptr;    // not sure if needed -- remove if unnecessary
  if (pairflag) {
    pair = force->pair_match(pair_str, 1);
    if (pair == nullptr)
      error->all(FLERR, "Fix electrode couldn't find the pair style {}", pair_str);
  } else {
    pair = force->pair_match("coul", 0);
    if (pair ==
        nullptr) {    // couldn't find a pair with name coul -- maybe hybrid return 1st hybrid substyle containing 'coul'
      pair = force->pair_match("coul", 0, 1);
    }
  }
  if (pair == nullptr) error->all(FLERR, "Fix electrode couldn't find a Coulombic pair style");

  // error if more than one fix electrode/*
  if (modify->get_fix_by_style("^electrode").size() > 1)
    error->all(FLERR, "More than one fix electrode");

  // make sure electrode atoms are not integrated if a matrix is used for electrode-electrode interaction
  int const nlocal = atom->nlocal;
  int *mask = atom->mask;
  if (matrix_algo) {
    std::vector<Fix *> integrate_fixes;
    for (auto fix : modify->get_fix_list()) {
      if (fix->time_integrate == 0) continue;
      int electrode_mover = 0;
      int fix_groupbit = fix->groupbit;
      for (int j = 0; j < nlocal; j++)
        if ((mask[j] & fix_groupbit) && (mask[j] & groupbit)) electrode_mover = 1;
      MPI_Allreduce(MPI_IN_PLACE, &electrode_mover, 1, MPI_INT, MPI_SUM, world);
      if (electrode_mover && comm->me == 0) integrate_fixes.push_back(fix);
    }
    if (comm->me == 0)
      for (const auto fix : integrate_fixes)
        error->warning(FLERR,
                       "Electrode atoms are integrated by fix {} {}, but fix electrode is using a "
                       "matrix method. For mobile electrodes use the conjugate gradient algorithm "
                       "without matrix ('algo cg').",
                       fix->id, fix->style);
  }

  // check for package intel
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

void FixElectrodeConp::post_constructor()    // TODO move to solver?
{
  if (!ffield) return;
  // ffield: test conditions and set up efield
  if (num_of_groups != 2) error->all(FLERR, "Number of electrodes must be two with ffield yes");
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
  if ((modify->get_fix_by_style("^efield").size() > 0) && (comm->me == 0))
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

  // if Thomas-Fermi, make sure all electrode atoms have parameters
  if (tfflag) {
    int unset_tf = 0;
    int *type = atom->type;
    for (int i = 0; i < nlocal; i++) {
      if ((groupbit & mask[i]) && (tf_types.count(type[i]) == 0)) unset_tf++;
    }
    MPI_Allreduce(MPI_IN_PLACE, &unset_tf, 1, MPI_INT, MPI_SUM, world);
    if (unset_tf)
      error->all(FLERR, "Thomas-Fermi parameters not set for all types in fix {}", style);
  }

  // get equal-style variable ids:
  group_psi_var_ids = std::vector<int>(num_of_groups, -1);
  for (int g = 0; g < num_of_groups; g++) {
    assert(group_psi_var_styles[g] != VarStyle::UNSET);
    if (group_psi_var_styles[g] == VarStyle::CONST) continue;
    const char *var_name = group_psi_var_names[g].c_str();
    int var_id = input->variable->find(var_name);
    if (var_id < 0) error->all(FLERR, "Variable '{}' for fix {} does not exist", var_name, style);
    if (!input->variable->equalstyle(var_id))
      error->all(FLERR, "Variable '{}' for fix {} is not equal-style", var_name, style);
    group_psi_var_ids[g] = var_id;
  }
  if (qtotal_var_style == VarStyle::EQUAL) {
    const char *var_name = qtotal_var_name.c_str();
    int var_id = input->variable->find(var_name);
    if (var_id < 0) error->all(FLERR, "Variable '{}' for fix electrode does not exist", var_name);
    if (!input->variable->equalstyle(var_id))
      error->all(FLERR, "Variable '{}' for fix electrode is not equal-style", var_name);
    qtotal_var_id = var_id;
  }

  // pair and list setups:

  evscale = force->qe2f / force->qqrd2e;
  elyt_vector->setup_general(pair, vec_neighlist, pairflag, timer_flag);
  if (etapropflag) elyt_vector->setup_eta(eta_index);
  if (need_elec_vector) {
    elec_vector->setup_general(pair, mat_neighlist, pairflag, timer_flag);
    if (etapropflag) elec_vector->setup_eta(eta_index);
    if (tfflag) elec_vector->setup_tf(tf_types);
    if (hardnessflag) elec_vector->setup_hardness(hardness_index);
  }

  if (matrix_algo) {
    assert(taglist_constructed);
    memory->destroy(matrix);
    memory->create(matrix, ngroup, ngroup, "fix_electrode:matrix");
    if (read_mat)
      electrode_taglist->read_from_file(input_file_mat, matrix, "elastance");
    else if (!read_inv) {
      if (etypes_neighlists) neighbor->build_one(mat_neighlist);
      auto array_compute = std::unique_ptr<ElectrodeMatrix>(new ElectrodeMatrix(lmp, igroup, eta));
      array_compute->setup(electrode_taglist->get_tag_to_iele(), pair, mat_neighlist, pairflag);
      if (etapropflag) array_compute->setup_eta(eta_index);
      if (tfflag) array_compute->setup_tf(tf_types);
      if (hardnessflag) array_compute->setup_hardness(hardness_index);
      array_compute->compute_array(matrix, timer_flag);
    } else
      assert(algo == Algo::MATRIX_INV);
    // write_mat before proceeding
    if (write_mat) electrode_taglist->write_to_file(output_file_mat, matrix);
  }
  // construct charge solver
  switch (algo) {
    case Algo::MATRIX_INV: {
      assert(taglist_constructed);
      ElectrodeInv *inv = new ElectrodeInv(lmp);
      if (read_inv) {
        if (comm->me == 0 && ffield)
          error->warning(FLERR,
                         "Symmetrizing matrix from file. Make sure the provided matrix has not "
                         "been symmetrized yet.");
        electrode_taglist->read_from_file(input_file_inv, matrix, "capacitance");
        inv->set_capacitance(ngroup, matrix);
      } else {
        inv->set_elastance(ngroup, matrix);
      }
      assert(taglist_constructed);
      inv->setup_solver(groupbit, electrode_taglist->get_tag_to_iele(), group_bits, ffield);
      charge_solver = inv;
      break;
    }
    case Algo::MATRIX_CG: {
      ElectrodeMatCG *mat_cg = new ElectrodeMatCG(lmp);
      mat_cg->set_elastance(ngroup, matrix);
      mat_cg->setup_solver(cg_threshold, electrode_taglist->get_tag_to_iele());
      charge_solver = mat_cg;
      break;
    }
    case Algo::CG: {
      ElectrodeCG *cg = new ElectrodeCG(lmp);
      cg->setup_solver(cg_threshold, elec_vector);
      charge_solver = cg;
      break;
    }
    default:
      error->all(FLERR, "This algorithm is not implemented, yet");
  }
  if (qtotal_var_style == VarStyle::CONST) charge_solver->set_constraint(qtotal);
  // initial charges and b vector
  update_charges();

  // write to files, ordered by group
  if (write_vec) {
    memset(potential_i, 0, atom->nmax * sizeof(double));
    elyt_vector->compute_pot(potential_i);
    double *potential_iele;
    memory->create(potential_iele, ngroup, "FixElectrode:potential_iele");
    charge_solver->buffer_and_gather(potential_i, potential_iele);
    electrode_taglist->write_to_file(output_file_vec, potential_iele);
    memory->destroy(potential_iele);
  }
  if (write_inv) electrode_taglist->write_to_file(output_file_inv, matrix);
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::setup_pre_reverse(int eflag, int vflag)
{
  if (pair->did_tally_callback() && (comm->me == 0))
    error->warning(FLERR, "Computation of virials in fix {} is incompatible with TALLY package",
                   style);
  // correct forces for initial timestep
  ev_init(eflag, vflag);
  gausscorr(eflag, vflag, true);
  self_energy(eflag);
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::setup_pre_exchange()
{
  nlocalele_outdated = 1;    // force regather

  if (!matrix_algo) return;

  if (!taglist_constructed) {
    electrode_taglist = new ElectrodeTaglist(lmp, group_bits);
    taglist_constructed = true;
  }

  //if memory_usage > 0.5 GiB, warn with expected usage
  double mem_needed = memory_usage();
  mem_needed /= (1024 * 1024 * 1024);    // convert to GiB
  if ((mem_needed > 0.5) && (comm->me == 0))
    error->warning(FLERR,
                   "Please ensure there is sufficient memory for fix {} "
                   "(anticipated usage is at least {:.1f} GiB per proc)",
                   style, mem_needed);
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::pre_force(int)
{
  update_charges();
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::pre_reverse(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  gausscorr(eflag, vflag, true);
  self_energy(eflag);
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

/* ----------------------------------------------------------------------
    Solve the equation of the constant potential method.
    I.e., minimize the energy as function of electrode atom charges.
------------------------------------------------------------------------- */

void FixElectrodeConp::update_charges()
{
  MPI_Barrier(world);
  double start = MPI_Wtime();
  if (atom->nmax > nmax) {
    memory->destroy(potential_i);
    nmax = atom->nmax;
    memory->create(potential_i, nmax, "FixElectrode:potential_i");
  }
  gather_list_iele();
  memset(potential_i, 0., atom->nmax * sizeof(double));
  elyt_vector->compute_pot(potential_i);
  if (enflag) add_electronegativity(potential_i);
  charge_solver->set_elyt_pot(potential_i);
  update_psi_set_constraint();
  set_charges(charge_solver->solve(group_psi));
  MPI_Barrier(world);
  update_time += MPI_Wtime() - start;
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::set_charges(std::vector<double> q_local)
{
  assert((int) q_local.size() == nlocalele);
  double *q = atom->q;
  for (int i = 0; i < nlocalele; i++) q[atom->map(taglist_local[i])] = q_local[i];
  comm->forward_comm(this);
  intel_pack_buffers();
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::update_psi_set_constraint()
{
  for (int g = 0; g < num_of_groups; g++) {
    if (group_psi_var_styles[g] == VarStyle::CONST)
      group_psi[g] = group_psi_const[g];
    else
      group_psi[g] = input->variable->compute_equal(group_psi_var_ids[g]);
  }
  if (qtotal_var_style != VarStyle::UNSET) {
    if (qtotal_var_style == VarStyle::EQUAL) qtotal = input->variable->compute_equal(qtotal_var_id);
    charge_solver->set_constraint(qtotal);
  }
}

/* ---------------------------------------------------------------------- */

double FixElectrodeConp::compute_scalar()
{
  return potential_energy();
}

/* ---------------------------------------------------------------------- */

double FixElectrodeConp::compute_vector(int i)
{
  return charge_solver->get_potential(i);
}

/* ---------------------------------------------------------------------- */

double FixElectrodeConp::compute_array(int i, int j)
{
  if (j == 0)
    return charge_solver->get_sb_charges(i);
  else if (j <= num_of_groups)
    return charge_solver->get_macro_capacitance(i, j - 1);
  else if (j <= 2 * num_of_groups)
    return charge_solver->get_macro_elastance(i, j - num_of_groups - 1);
  else
    return 0.;    // avoid -Wreturn-type warning
}

/* ---------------------------------------------------------------------- */

double FixElectrodeConp::potential_energy()
{
  // corrections to energy due to potential psi
  double const qqrd2e = force->qqrd2e;
  int const nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *q = atom->q;
  double energy = 0;
  for (int i = 0, iele = 0; i < nlocal; i++) {
    if (groupbit & mask[i]) {
      energy -= qqrd2e * q[i] * group_psi[iele_to_group_local[iele]] * evscale;
      iele++;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPI_DOUBLE, MPI_SUM, world);
  return energy;
}

/* ---------------------------------------------------------------------- */

double FixElectrodeConp::self_energy(int eflag)
{
  // corrections to energy due to self interaction
  double energy = 0.;
  double const qqrd2e = force->qqrd2e;
  int const nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  double *q = atom->q;
  if (tfflag) {
    for (int i = 0; i < nlocal; i++) {
      if (groupbit & mask[i]) {
        double const e = 0.5 * qqrd2e * q[i] * q[i] * tf_types[type[i]];
        energy += e;
        if (eflag) force->pair->ev_tally(i, i, nlocal, force->newton_pair, 0., e, 0, 0, 0, 0);
      }
    }
  }
  if (hardnessflag) {
    double *d_hardness = atom->dvector[hardness_index];
    for (int i = 0; i < nlocal; i++) {
      if (groupbit & mask[i]) {
        double const e = 0.5 * q[i] * q[i] * d_hardness[i];
        energy += e;
        if (eflag) force->pair->ev_tally(i, i, nlocal, force->newton_pair, 0., e, 0, 0, 0, 0);
      }
    }
  }
  if (!pairflag) {
    double const pre = 1. / sqrt(MY_2PI) * qqrd2e;
    for (int i = 0; i < nlocal; i++) {
      if (groupbit & mask[i]) {
        double ieta = etapropflag ? atom->dvector[eta_index][i] : eta;
        double e = ieta * pre * q[i] * q[i];
        energy += e;
        if (eflag) { force->pair->ev_tally(i, i, nlocal, force->newton_pair, 0., e, 0, 0, 0, 0); }
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPI_DOUBLE, MPI_SUM, world);
  return energy;
}

/* ---------------------------------------------------------------------- */

double FixElectrodeConp::gausscorr(int eflag, int vflag, bool fflag)
{
  // correction to short range interaction due to eta
  if (pairflag) return 0.;
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
    double const eta_i = etapropflag ? atom->dvector[eta_index][i] : eta;
    int itype = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int const j = jlist[jj] & NEIGHMASK;
      bool j_in_ele = groupbit & mask[j];
      if (!(i_in_ele || j_in_ele)) continue;

      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];
      double rsq = delx * delx + dely * dely + delz * delz;
      int jtype = type[j];

      if (rsq < force->pair->cutsq[itype][jtype]) {
        double const eta_j = etapropflag ? atom->dvector[eta_index][j] : eta;
        double eta_ij;
        if (i_in_ele && j_in_ele)
          eta_ij = eta_i * eta_j / sqrt(eta_i * eta_i + eta_j * eta_j);
        else if (i_in_ele)
          eta_ij = eta_i;
        else {
          assert(j_in_ele);
          eta_ij = eta_j;
        }
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
        if (eflag) {
          double ecoul = -prefactor * erfc_etar;
          force->pair->ev_tally(i, j, nlocal, newton_pair, 0., ecoul, 0., 0., 0., 0.);
        }
        if (vflag) v_tally(i, j, nlocal, newton_pair, fpair, delx, dely, delz);
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &energy_sr, 1, MPI_DOUBLE, MPI_SUM, world);
  return energy_sr;
}

/* ---------------------------------------------------------------------- */

FixElectrodeConp::~FixElectrodeConp()
{
  if (comm->me == 0) {
    try {
      if (timer_flag) {
        if (charge_solver != nullptr)
          utils::logmesg(lmp, "Multiplication time: {:.4g} s\n", charge_solver->get_mult_time());
        utils::logmesg(lmp, "Update time: {:.4g} s\n", update_time);
      }
    } catch (std::exception &) {
    }
  }

  if (modify->get_fix_by_id(id)) atom->delete_callback(id, Atom::GROW);

  memory->destroy(potential_i);

  delete elyt_vector;
  memory->destroy(matrix);
  if (need_elec_vector) delete elec_vector;
  if (charge_solver != nullptr) delete charge_solver;
  if (taglist_constructed) delete electrode_taglist;
}

/* ---------------------------------------------------------------------- */

int FixElectrodeConp::setmask()
{
  int mask = 0;
  mask |= FixConst::PRE_EXCHANGE;
  mask |= FixConst::POST_NEIGHBOR;
  mask |= FixConst::PRE_FORCE;
  mask |= FixConst::PRE_REVERSE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::request_etypes_neighlists()
{
  int const ntypes = atom->ntypes;
  // construct etypes
  int *mask = atom->mask;
  int *type = atom->type;
  auto elec = std::vector<int>(ntypes, 0);
  auto elyt = std::vector<int>(ntypes, 0);
  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit)
      elec[type[i] - 1] += 1;
    else
      elyt[type[i] - 1] += 1;
  }
  MPI_Allreduce(MPI_IN_PLACE, elec.data(), ntypes, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, elyt.data(), ntypes, MPI_INT, MPI_SUM, world);
  etypes.clear();
  for (int i = 0; i < ntypes; i++) {
    if (!elec[i] == !elyt[i]) error->all(FLERR, "Types overlap, cannot use etypes keyword");
    if (elec[i]) etypes.push_back(i + 1);
  }
  // construct skip arrays
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

  if (need_array_compute) {
    auto matReq = neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
    matReq->set_skip(iskip_mat, ijskip_mat);
    matReq->set_id(1);
    if (intelflag) matReq->enable_intel();
  } else if (need_elec_vector) {
    auto matReq = neighbor->add_request(this);
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

/* ---------------------------------------------------------------------- */

int FixElectrodeConp::pack_exchange(int i, double * /* buf */)
{
  if (atom->mask[i] & groupbit) {
    nlocalele_outdated = 1;
    nlocalele--;    // decrement nlocalele if we are packing away a particle
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

int FixElectrodeConp::unpack_exchange(int nlocal, double * /* buf */)
{
  if (atom->mask[nlocal] & groupbit) {    // this should work
    nlocalele_outdated = 1;
    nlocalele++;    // increment nlocalele if we are unpacking a particle
  }
  return 0;
}

/* ----------------------------------------------------------------------
    Update taglist_local and iele_to_group_local, when necessary
------------------------------------------------------------------------- */

void FixElectrodeConp::gather_list_iele()
{
  MPI_Allreduce(MPI_IN_PLACE, &nlocalele_outdated, 1, MPI_INT, MPI_SUM, world);
  if (nlocalele_outdated == 0) return;

  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int const nlocal = atom->nlocal;
  taglist_local.clear();
  iele_to_group_local.clear();
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      tagint const t = tag[i];
      taglist_local.push_back(t);
      for (int g = 0; g < num_of_groups; g++)
        if (mask[i] & group_bits[g]) iele_to_group_local.push_back(g);
    }
  }
  nlocalele = static_cast<int>(taglist_local.size());    // just for safety
  assert((int) iele_to_group_local.size() == nlocalele);

  charge_solver->update_solver(taglist_local, iele_to_group_local);
  nlocalele_outdated = 0;
}

/* ---------------------------------------------------------------------- */

double FixElectrodeConp::memory_usage()
{
  int const nmax = atom->nmax;
  double bytes = 0.0;
  if (taglist_constructed) bytes += electrode_taglist->memory_usage();
  if (charge_solver != nullptr) bytes += charge_solver->memory_use();
  bytes += nmax * (sizeof(double));    // potential_i
  bytes += taglist_local.capacity() * sizeof(tagint);
  bytes += iele_to_group_local.capacity() * sizeof(int);
  return bytes;
}

/* ---------------------------------------------------------------------- */

int FixElectrodeConp::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                        int * /*pbc*/)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int const j = list[i];
    buf[m++] = atom->q[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::unpack_forward_comm(int n, int first, double *buf)
{
  int const last = first + n;
  for (int i = first, m = 0; i < last; i++) atom->q[i] = buf[m++];
}

/* ----------------------------------------------------------------------
   Tally virial of pair interactions in pre_reverse. This cannot be done with pair->ev_tally()
   because compute_fdotr is called before pre_reverse, i.e. Virials need to be tallied even if fdotr
   is used.
------------------------------------------------------------------------- */

void FixElectrodeConp::v_tally(int i, int j, int nlocal, int newton_pair, double fpair, double delx,
                               double dely, double delz)
{
  double v[6];
  if (vflag_either) {
    v[0] = delx * delx * fpair;
    v[1] = dely * dely * fpair;
    v[2] = delz * delz * fpair;
    v[3] = delx * dely * fpair;
    v[4] = delx * delz * fpair;
    v[5] = dely * delz * fpair;

    if (vflag_global) {
      if (newton_pair) {
        virial[0] += v[0];
        virial[1] += v[1];
        virial[2] += v[2];
        virial[3] += v[3];
        virial[4] += v[4];
        virial[5] += v[5];
      } else {
        if (i < nlocal) {
          virial[0] += 0.5 * v[0];
          virial[1] += 0.5 * v[1];
          virial[2] += 0.5 * v[2];
          virial[3] += 0.5 * v[3];
          virial[4] += 0.5 * v[4];
          virial[5] += 0.5 * v[5];
        }
        if (j < nlocal) {
          virial[0] += 0.5 * v[0];
          virial[1] += 0.5 * v[1];
          virial[2] += 0.5 * v[2];
          virial[3] += 0.5 * v[3];
          virial[4] += 0.5 * v[4];
          virial[5] += 0.5 * v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_pair || i < nlocal) {
        vatom[i][0] += 0.5 * v[0];
        vatom[i][1] += 0.5 * v[1];
        vatom[i][2] += 0.5 * v[2];
        vatom[i][3] += 0.5 * v[3];
        vatom[i][4] += 0.5 * v[4];
        vatom[i][5] += 0.5 * v[5];
      }
      if (newton_pair || j < nlocal) {
        vatom[j][0] += 0.5 * v[0];
        vatom[j][1] += 0.5 * v[1];
        vatom[j][2] += 0.5 * v[2];
        vatom[j][3] += 0.5 * v[3];
        vatom[j][4] += 0.5 * v[4];
        vatom[j][5] += 0.5 * v[5];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixElectrodeConp::add_electronegativity(double *b)
{
  assert(enflag);
  int *mask = atom->mask;
  double *d_en = atom->dvector[en_index];
  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) b[i] += d_en[i] / force->qqrd2e;
  }
}

/* ---------------------------------------------------------------------- */
