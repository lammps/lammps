/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_mbx.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "pair_mbx.h"
#include "respa.h"
#include "safe_pointers.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <mpi.h>

#include "bblock/system.h"

#define _MAX_SIZE_MOL_NAME 16
// Subject to further increase _MAX_SIZE_MOL_NAME
#define _MAX_ATOMS_PER_MONOMER 8
#define SMALL 1.0e-4


namespace LAMMPS_NS {
//PImpl idiom to hide MBX implementation details
struct MBXImpl {
  MBXImpl() : ptr_mbx(nullptr), ptr_mbx_local(nullptr) {}
  ~MBXImpl()
  {
    delete ptr_mbx;
    delete ptr_mbx_local;
  }
  bblock::System *ptr_mbx;
  bblock::System *ptr_mbx_local;
};
} // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace FixConst;


static constexpr double MBX_CUTOFF_TOL = 1e-9;

std::string FixMBX::cite_pair_mbx = std::string(
    "pair mbx command:\n\n" \
    "@article{10.1063/5.0156036,\n" \
    " author = {Riera, Marc and Knight, Christopher and Bull-Vulpe, Ethan F. and Zhu, Xuanyu and " \
    "Agnew, Henry and Smith, Daniel G. A. and Simmonett, Andrew C. and Paesani, Francesco},\n" \
    " title = \"{MBX: A many-body energy and force calculator for data-driven many-body " \
    "simulations}\",\n" \
    " journal = {The Journal of Chemical Physics},\n" \
    " volume = {159},\n" \
    " number = {5},\n" \
    " pages = {054802},\n" \
    " year = {2023},\n" \
    " doi = {10.1063/5.0156036},\n" \
    " version = {") + bblock::System::get_mbx_version() + "}\n"  \
    "}\n\n";

/* ---------------------------------------------------------------------- */

// Parse atom ID or range for dp1 monomer
// Input can be a single integer (e.g., "5") or a range (e.g., "1*11")
std::pair<int, int> FixMBX::parse_dp1_range(const std::string &dp1_str)
{
  try {
    size_t pos = 0;
    int val = std::stoi(dp1_str, &pos);
    if (pos == dp1_str.size()) {
      // Single integer
      if (val < 1)
        throw std::invalid_argument("[MBX] Atom ID " + dp1_str + " in dp1 must be positive");
      return {val, val};
    }
  } catch (...) {
    // Not a single integer, try to parse as range
  }

  // Try to parse as a range (e.g., "1*11")
  size_t star = dp1_str.find('*');
  if (star == std::string::npos)
    throw std::invalid_argument("[MBX] Invalid dp1 format: " + dp1_str);

  std::string start_str = dp1_str.substr(0, star);
  std::string end_str = dp1_str.substr(star + 1);
  int start = 0, end = 0;
  try {
    size_t spos = 0, epos = 0;
    start = std::stoi(start_str, &spos);
    end = std::stoi(end_str, &epos);

    if (spos != start_str.size() || epos != end_str.size())
      throw std::invalid_argument("Non-integer characters found in range");
  } catch (...) {
    throw std::invalid_argument("[MBX] Invalid integers in dp1 range: " + dp1_str);
  }
  if (start < 1 || end < start)
    throw std::invalid_argument("[MBX] Invalid range values for dp1: " + dp1_str);

  return {start, end};
}

// Validate that the input arguments to fix mbx are correct
// throws specific errors if the arguments are malformed
bool FixMBX::validateMBXFixParameters(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, ("[MBX] Input line too short"));

  int num_monomers = 0;
  try {
    num_monomers = std::stoi(arg[0]);
  } catch (...) {
    error->all(FLERR, ("[MBX] num_monomers is not a valid integer: " + std::string(arg[0])));
  }
  if (num_monomers < 1) error->all(FLERR, ("[MBX] num_monomers must be positive"));

  std::map<int, std::string> mbx_atom_id_mapping;    // atom ID mapping for all monomers

  int input_validation_index = 1;    // part of arg currently being validated

  // Lambda to check dp1 monomer syntax
  auto check_external_dp1 = [&](int n_atoms, char **current_monomer_atoms) -> bool {
    if (n_atoms != 1)
      error->all(
          FLERR,
          ("[MBX] Wrong number of arguments for dp1: expected 1, got " + std::to_string(n_atoms)));
    const std::string atom_id_str = current_monomer_atoms[0];

    int start_index, end_index;
    try {
      std::tie(start_index, end_index) = parse_dp1_range(atom_id_str);
    } catch (const std::exception &e) {
      error->all(FLERR, e.what());
    }
    if (start_index < 1 || end_index < start_index)
      error->all(FLERR, ("[MBX] Invalid range values for dp1: " + atom_id_str));

    for (int at = start_index; at <= end_index; ++at) {
      if (mbx_atom_id_mapping.count(at))
        error->all(FLERR, ("[MBX] Already defined atom IDs found in dp1: " + std::to_string(at)));
      mbx_atom_id_mapping[at] = "X";
    }
    return true;
  };

  // Lambda to check syntax of a monomer with n_atoms and atom IDs in current_monomer
  auto check_monomer_syntax = [&](int n_atoms, char **current_monomer) -> bool {
    std::string current_monomer_name = current_monomer[0];
    std::vector<std::string> current_monomer_atoms;

    for (int i = 1; i <= n_atoms; ++i) current_monomer_atoms.push_back(current_monomer[i]);
    std::vector<std::string> expected_monomer_atom_ids;

    // special handling for dp1 monomer
    if (current_monomer_name == "dp1") return check_external_dp1(n_atoms, &current_monomer[1]);

    try {
      add_monomer_atom_types(const_cast<char *>(current_monomer_name.c_str()),
                             expected_monomer_atom_ids);
    } catch (const std::exception &e) {
      error->all(FLERR, ("[MBX] Invalid monomer name " + current_monomer_name));
    }

    if (expected_monomer_atom_ids.size() != n_atoms)
      error->all(FLERR,
                 ("[MBX] Wrong number of atoms: expected " + std::to_string(expected_monomer_atom_ids.size()) + ", got " +
                  std::to_string(n_atoms)));

    // validate that atom IDs are positive integers
    std::vector<int> atom_ids;
    for (size_t i = 0; i < expected_monomer_atom_ids.size(); ++i) {
      int at = 0;
      try {
        at = std::stoi(current_monomer_atoms[i]);
      } catch (...) {
        error->all(FLERR,
                   ("[MBX] Atom ID " + current_monomer_atoms[i] + " is not a valid integer"));
      }
      if (at < 1)
        error->all(FLERR, ("[MBX] Atom ID " + current_monomer_atoms[i] + " must be positive"));
      atom_ids.push_back(at);
    }

    // check that number of unique atom IDs is correct
    // this verifies that the monomer has the expected number of different elements
    std::set<int> unique_atom_ids(atom_ids.begin(), atom_ids.end());
    std::set<std::string> unique_monomer_atom_ids(expected_monomer_atom_ids.begin(),
                                                  expected_monomer_atom_ids.end());
    if (unique_atom_ids.size() < unique_monomer_atom_ids.size())
      error->all(FLERR,
                 ("[MBX] Wrong number of unique atom IDs in " + current_monomer_name +
                  ": expected " + std::to_string(unique_monomer_atom_ids.size()) + ", got " +
                  std::to_string(unique_atom_ids.size())));
    std::map<int, std::string> atom_mapping;

    // check that atom ID mapping is consistent
    // atom ID mapping must match expected, such as OHH for h2o
    for (size_t i = 0; i < atom_ids.size(); ++i) {
      int at = atom_ids[i];
      if (!atom_mapping.count(at))    // first time seeing this atom ID
        atom_mapping[at] = expected_monomer_atom_ids[i];
      else if (atom_mapping[at] != expected_monomer_atom_ids[i]) {    // inconsistent mapping detected
        // construct error message
        std::string expected_monomer_atom_ids_string = "";
        for (const auto &mat : expected_monomer_atom_ids) { expected_monomer_atom_ids_string += mat + " "; }
        std::string atom_ids_string = "";
        for (const auto &at2 : atom_ids) { atom_ids_string += std::to_string(at2) + " "; }
        error->all(FLERR,
                   ("[MBX] Incorrect atom ID mapping in " + current_monomer_name + ". Expected " +
                    expected_monomer_atom_ids_string + "but got " + atom_ids_string));
      }
      if (mbx_atom_id_mapping.count(at))    // atom ID already defined in another monomer
        error->all(FLERR,
                   ("[MBX] Already defined atom IDs found in " + current_monomer_name + ": " +
                    std::to_string(at)));
    }

    // check that atom IDs are contiguous
    int minimum_index = *std::min_element(atom_ids.begin(), atom_ids.end());
    int maximum_index = *std::max_element(atom_ids.begin(), atom_ids.end());
    for (int i = minimum_index; i <= maximum_index; ++i) {
      if (std::find(atom_ids.begin(), atom_ids.end(), i) == atom_ids.end())
        error->all(FLERR,
                   ("[MBX] Atom IDs must be contiguous in " + current_monomer_name + ". Missing " +
                    std::to_string(i)));
    }
    for (const auto &kv : atom_mapping)    // transfer monomer mapping to global mapping
      mbx_atom_id_mapping[kv.first] = kv.second;
    return true;
  };

  // validate each monomer
  for (int monomer_index = 0; monomer_index < num_monomers; ++monomer_index) {
    if (input_validation_index >= narg)
      error->all(FLERR, ("[MBX] Not enough arguments to read a monomer name"));
    std::string monomer_name = arg[input_validation_index];
    int num_atoms = 0;

    // ensure monomer name is valid and get number of atoms
    try {
      bool is_ext = false;
      num_atoms = get_num_atoms_per_monomer(const_cast<char *>(monomer_name.c_str()), is_ext);
    } catch (const std::exception &e) {
      error->all(FLERR, ("[MBX] Invalid monomer name " + monomer_name));
    }
    if (input_validation_index + num_atoms >= narg)
      error->all(FLERR, ("[MBX] Not enough arguments to read monomer atoms for " + monomer_name));

    check_monomer_syntax(num_atoms, &arg[input_validation_index]);
    input_validation_index += num_atoms + 1;
  }

  // process remaining optional keywords
  while (input_validation_index < narg) {
    if (strcmp(arg[input_validation_index], "json") == 0) {
      if (input_validation_index + 1 >= narg)
        error->all(FLERR, ("[MBX] Not enough arguments to read json filename"));
      input_validation_index += 2;
    } else if (strcmp(arg[input_validation_index], "print/verbose") == 0) {
      input_validation_index += 1;
    } else if (strcmp(arg[input_validation_index], "print/dipoles") ==
               0) {    // dipoles are now always printed by default
      input_validation_index += 1;
    } else if (strcmp(arg[input_validation_index], "aspc/reset") == 0) {
      if (input_validation_index + 1 >= narg)
        error->all(FLERR, ("[MBX] Not enough arguments to read aspc/reset value"));
      int aspc_reset = 0;
      try {
        aspc_reset = std::stoi(arg[input_validation_index + 1]);
      } catch (...) {
        error->all(FLERR, ("[MBX] aspc/reset value is not a valid integer"));
      }
      if (aspc_reset < 1) error->all(FLERR, ("[MBX] aspc/reset value must be positive"));
      input_validation_index += 2;
    } else {
      error->all(FLERR, ("[MBX] Unknown keyword: " + std::string(arg[input_validation_index])));
    }
  }

  return true;
}

FixMBX::FixMBX(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  // Constructor for fix mbx, called interally by pair mbx.
  // Expected arguments:
  // _FIX_MBX_INTERNAL all MBX num_mol_types mon_name atom_mapping <mon_name2> <atom_mapping2> ... json mbx.json
  // num_mol_types = number of monomer types in the system
  // mon_name = name of the monomer type (e.g. h2o, ch4, etc)
  // atom mapping = list of LAMMPS atom types that correspond to the atoms in the monomer
  // json arg = specifies the name of the MBX json configuration file, such as mbx.json

  if (lmp->citeme) lmp->citeme->add(cite_pair_mbx);

  me = comm->me;
  nprocs = comm->nprocs;
  mbx_impl = new MBXImpl; //PImpl idiom to hide MBX implementation details

  // // validate input arguments
  bool validation_result = validateMBXFixParameters(narg - 3, &arg[3]);

  if (narg < 6) error->all(FLERR, "[MBX] Illegal fix mbx command");

  num_mol_types = utils::inumeric(FLERR, arg[3], false, lmp);

  if (num_mol_types < 1) error->all(FLERR, "[MBX] Illegal fix mbx command");

  num_atoms_per_mol = NULL;
  mol_names = NULL;
  lower_atom_type_index_in_mol = NULL;
  higher_atom_type_index_in_mol = NULL;

  memory->create(num_atoms_per_mol, num_mol_types, "fixmbx:num_atoms_per_mol");
  memory->create(mol_names, num_mol_types, _MAX_SIZE_MOL_NAME, "fixmbx:mol_names");
  memory->create(lower_atom_type_index_in_mol, num_mol_types,
                 "fixmbx:lower_atom_type_index_in_mol");
  memory->create(higher_atom_type_index_in_mol, num_mol_types,
                 "fixmbx:higher_atom_type_index_in_mol");
  memory->create(order_in_mol, num_mol_types, _MAX_ATOMS_PER_MONOMER, "fixmbx:order_in_mol");

  // Extract information about min and max indexes
  int iarg = 4;
  bool is_ext = false;
  for (int i = 0; i < num_mol_types; ++i) {
    std::string current_monomer_name = arg[iarg++];

    if (strlen(current_monomer_name.c_str()) >= _MAX_SIZE_MOL_NAME)
      error->all(FLERR,
                 "[MBX] Monomer name too long: did developer correctly add support for monomer?");

    strcpy(mol_names[i], current_monomer_name.c_str());

    int current_lower_index = _MAX_ATOMS_PER_MONOMER + 1;
    int current_higher_index = -1;

    int current_n_atoms =
        get_num_atoms_per_monomer(const_cast<char *>(current_monomer_name.c_str()), is_ext);

    if (current_n_atoms > _MAX_ATOMS_PER_MONOMER)
      error->all(FLERR,
                 "[MBX] num_atoms_per_mol > _MAX_ATOMS_PER_MONOMER : did developer correctly add "
                 "support for monomer?");

    if (current_monomer_name == "dp1") {
      // dp1 can either by a single atom ID or a range of IDs in the form "1*11"
      int start, end;
      try {

        std::tie(start, end) = parse_dp1_range(std::string(arg[iarg++]));
      } catch (const std::exception &e) {
        error->all(FLERR, e.what());
      }

      lower_atom_type_index_in_mol[i] = start;
      higher_atom_type_index_in_mol[i] = end;

      order_in_mol[i][0] = start;    // only one atom in dp1

    } else {    // handling all other monomers
      // find min and max atom type index in the mapping
      for (int j = 0; j < current_n_atoms; j++) {
        int current_index = utils::inumeric(FLERR, arg[iarg + j], false, lmp);
        if (current_index < current_lower_index) current_lower_index = current_index;
        if (current_index > current_higher_index) current_higher_index = current_index;

        lower_atom_type_index_in_mol[i] = current_lower_index;
        higher_atom_type_index_in_mol[i] = current_higher_index;
      }

      for (int j = 0; j < current_n_atoms; j++) {
        order_in_mol[i][j] = utils::inumeric(FLERR, arg[iarg++], false, lmp);
      }
    }
  }

  // process remaining optional keywords

  use_json = 0;
  std::string json_file;
  print_verbose = 0;
  print_dipoles = 1;    // dipoles are now always printed by default
  aspc_step_reset = 1000;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "json") == 0) {
      json_file = std::string(arg[++iarg]);
      use_json = 1;
    } else if (strcmp(arg[iarg], "print/verbose") == 0) {
      print_verbose = 1;
    } else if (strcmp(arg[iarg], "print/dipoles") == 0) {
      print_dipoles = 1;    // dipoles are now always printed by default
    } else if (strcmp(arg[iarg], "aspc/reset") == 0) {
      aspc_step_reset = atoi(arg[++iarg]);
    } else {
      error->all(FLERR, "[MBX] Unknown keyword in fix mbx command: " + std::string(arg[iarg]));
    }

    iarg++;
  }

  mbx_aspc_enabled = false;

  // verify that this fix is not called directly
  // fix mbx should only be called indirectly by pair mbx
  pair_mbx = nullptr;
  pair_mbx = (PairMBX *) force->pair_match("^mbx", 0);
  if (!pair_mbx) error->all(FLERR, "[MBX] Pair mbx is missing");

  array_atom = NULL;
  mbx_dip = NULL;
  mol_type = NULL;
  mol_anchor = NULL;
  mol_local = NULL;

  grow_arrays(atom->nmax);

  // MRR Call function to fill up arrays
  mbx_fill_system_information_from_atom();

  atom->add_callback(0);
  atom->add_callback(1);

  mbx_num_atoms = 0;
  mbx_num_ext = 0;


  // instance of MBX with just local monomers

  mbx_num_atoms_local = 0;
  mbx_num_ext_local = 0;

  // check that LAMMPS proc mapping matches PME solver
  // MPI ranks must be mapped in "xyz" order for PME
  if (comm->style != 0) error->all(FLERR, "[MBX] Fix mbx must be used with comm_style brick");

  if (comm->layout != Comm::LAYOUT_UNIFORM)
    error->all(FLERR, "[MBX] Fix mbx must be used with comm layout of equal-sized bricks");

  {
    // compute expected proc coordinates for this MPI rank
    int proc_x = me % comm->procgrid[0];
    int proc_y = (me % (comm->procgrid[0] * comm->procgrid[1])) / comm->procgrid[0];
    int proc_z = me / (comm->procgrid[0] * comm->procgrid[1]);

    int e = 0;
    // compare actual location to expected
    if ((proc_x != comm->myloc[0]) || (proc_y != comm->myloc[1]) || (proc_z != comm->myloc[2]))
      e = 1;
    int err = 0;
    MPI_Allreduce(&e, &err, 1, MPI_INT, MPI_SUM, world);
    if (err)
      error->all(
          FLERR,
          "[MBX] Inconsistent proc mapping: 'processors * * * map xyz' required for PME solver");
  }

  // setup json, if requested

  if (use_json) {
    int size = 0;
    if (me == 0) {
      // Test if file present
      SafeFilePtr fp = fopen(json_file.c_str(), "r");
      if (fp == NULL) {
        error->one(FLERR, "Cannot open file " + json_file);
      }

      std::ifstream t(json_file);
      t.seekg(0, std::ios::end);
      size = t.tellg();
      json_settings.resize(size);
      t.seekg(0);
      t.read(&json_settings[0], size);
    }

    MPI_Bcast(&size, 1, MPI_INT, 0, world);
    if (me) json_settings.resize(size);

    MPI_Bcast(&json_settings[0], size + 1, MPI_CHAR, 0, world);
  }

  memory->create(mbxt_count, FixMBX::MBXT_LABELS::NUM_TIMERS, "fixmbx:mbxt_count");
  memory->create(mbxt_time, FixMBX::MBXT_LABELS::NUM_TIMERS, "fixmbx:mbxt_time");
  memory->create(mbxt_time_start, FixMBX::MBXT_LABELS::NUM_TIMERS, "fixmbx:mbxt_time_start");

  for (int i = 0; i < FixMBX::MBXT_LABELS::NUM_TIMERS; ++i) {
    mbxt_time[i] = 0.0;
    mbxt_count[i] = 0;
  }

  first_step = true;

  // MBX currently requires 4-byte tags

  if (sizeof(tagint) != sizeof(int)) error->all(FLERR, "[MBX] Tagints required to be of type int.");

  aspc_dip_hist = NULL;
  aspc_order = 6;    // hard-coded in MBX
  aspc_max_num_hist = aspc_order + 2;
  aspc_per_atom_size = aspc_max_num_hist * 3;    // (# of histories) * (# of dimensions)
  aspc_num_hist = 0;
  aspc_step = 0;

  peratom_flag = 1;
  size_peratom_cols = 9;
  peratom_freq = 1;

  comm_forward = aspc_per_atom_size;

  mbxt_initial_time = MPI_Wtime();
}

/* ---------------------------------------------------------------------- */

FixMBX::~FixMBX()
{
  if (mbx_aspc_enabled) memory->destroy(aspc_dip_hist);

  if (print_dipoles) memory->destroy(mbx_dip);

  memory->destroy(num_atoms_per_mol);
  memory->destroy(mol_names);

  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id, 0);
  atom->delete_callback(id, 1);

  memory->destroy(mol_local);
  memory->destroy(mol_type);
  memory->destroy(mol_anchor);
  memory->destroy(lower_atom_type_index_in_mol);
  memory->destroy(higher_atom_type_index_in_mol);
  memory->destroy(order_in_mol);

  if (mbx_impl->ptr_mbx_local) {
    // accumulate timing info from pme electrostatics

    std::vector<size_t> tmpi = mbx_impl->ptr_mbx_local->GetInfoElectrostaticsCounts();
    std::vector<double> tmpd = mbx_impl->ptr_mbx_local->GetInfoElectrostaticsTimings();

    for (int i = 0; i < tmpi.size(); ++i) {
      mbxt_count[FixMBX::MBXT_LABELS::ELE_PERMDIP_REAL + i] += tmpi[i];
      mbxt_time[FixMBX::MBXT_LABELS::ELE_PERMDIP_REAL + i] += tmpd[i];
    }

    // accumulate timing info from dispersion pme

    std::vector<size_t> tmpi_d = mbx_impl->ptr_mbx_local->GetInfoDispersionCounts();
    std::vector<double> tmpd_d = mbx_impl->ptr_mbx_local->GetInfoDispersionTimings();

    for (int i = 0; i < tmpi_d.size(); ++i) {
      mbxt_count[FixMBX::MBXT_LABELS::DISP_PME_SETUP + i] += tmpi_d[i];
      mbxt_time[FixMBX::MBXT_LABELS::DISP_PME_SETUP + i] += tmpd_d[i];
    }

  }
  delete mbx_impl;

  if (print_verbose)
      mbxt_write_summary();    // this and collecting times should be gated by 'timer full' request

  memory->destroy(mbxt_count);
  memory->destroy(mbxt_time);
  memory->destroy(mbxt_time_start);
}

/* ---------------------------------------------------------------------- */

void FixMBX::post_constructor() {}

/* ---------------------------------------------------------------------- */

int FixMBX::setmask()
{
  int mask = 0;
  mask |= POST_NEIGHBOR;
  mask |= MIN_POST_NEIGHBOR;
  mask |= PRE_FORCE;
  //  mask |= PRE_FORCE_RESPA;
  mask |= MIN_PRE_FORCE;

  mask |= PRE_EXCHANGE;    // only needs to be set when using ASPC integrator

  mask |= POST_FORCE;    // only needs to be set when printing dipoles

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMBX::init()
{
  if (!atom->q_flag) error->all(FLERR, "[MBX] Fix mbx requires atom attribute q");

  if (!atom->molecule_flag) error->all(FLERR, "[MBX] Fix mbx requires atom attribute molecule");

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR, "[MBX] Fix mbx group has no atoms");

}

/* ---------------------------------------------------------------------- */

// Fill mol_type and mol_anchor arrays from atom data
void FixMBX::mbx_fill_system_information_from_atom()
{
  const int nlocal = atom->nlocal;
  const int nghost = atom->nghost;
  const int nall = nlocal + nghost;


  bigint natoms = atom->natoms;

  tagint *tag = atom->tag;
  int *molecule = atom->molecule;
  double **x = atom->x;

  int mtype = -1;
  for (int i = 0; i < nall; ++i) {
    // Assign mol_type
    for (int j = 0; j < num_mol_types; ++j)
      if (atom->type[i] >= lower_atom_type_index_in_mol[j] and
          atom->type[i] <= higher_atom_type_index_in_mol[j]) {
        mol_type[i] = j;
        mtype = j;
        break;
        // If j is max and no type has been found, types in mbx fix do not match types in data file
      } else if (j == num_mol_types - 1) {
        error->all(FLERR,
                   "[MBX] The atom types in fix mbx do not match the atom types in the data file");
      }
  }

  // Reset anchors
  int *last_anchor = mol_anchor + nall;
  memset(mol_anchor, 0, sizeof(int)*nall);

  for (int i = 0; i < nall; ++i) {
    // Assign anchor TODO careful, not necessarily true
    // Create another peratom property -> index within molecules
    tagint itag = tag[i];
    int mtype = mol_type[i];
    bool is_ext = false;
    int na = get_num_atoms_per_monomer(mol_names[mtype], is_ext);
    if (is_ext) {
      mol_anchor[i] = 1;
      continue;
    }

    bool isanchor = true;
    for (int j = 0; j < na; j++) {
      int idx = atom->map(tag[i] + j);
      if (idx < 0) {
        isanchor = false;
        break;
      }

      isanchor = isanchor && atom->type[idx] == order_in_mol[mtype][j];
      if (!isanchor) break;
    }

    if (isanchor) mol_anchor[i] = 1;
  }
}

/* ---------------------------------------------------------------------- */

void FixMBX::setup_post_neighbor()
{

  // Figure out if there is a gcmc fix somewhere
  has_gcmc = false;
  int ifix = -1;
  for (int i = 0; i < modify->nfix; ++i)
    if (strcmp(modify->fix[i]->style, "gcmc") == 0) {
      if (ifix == -1)
        ifix = i;
      else
        error->all(FLERR, "[MBX] Only one GCMC fix instance allowed to be active");
    }
  if (ifix != -1) has_gcmc = true;

  grow_arrays(atom->nmax);

  // MRR Call function to fill up arrays
  mbx_fill_system_information_from_atom();
  post_neighbor();

  first_step = false;
}

/* ---------------------------------------------------------------------- */

void FixMBX::post_neighbor()
{
  // setup after neighbor build

  const int nlocal = atom->nlocal;
  const int nghost = atom->nghost;
  const int nall = nlocal + nghost;

  tagint *tag = atom->tag;
  int *molecule = atom->molecule;
  double **x = atom->x;

  mbx_fill_system_information_from_atom();

  // tear down existing MBX objects

  if (mbx_impl->ptr_mbx) delete mbx_impl->ptr_mbx;

  if (mbx_impl->ptr_mbx_local) {
    // accumulate timing info from pme electrostatics

    std::vector<size_t> tmpi = mbx_impl->ptr_mbx_local->GetInfoElectrostaticsCounts();
    std::vector<double> tmpd = mbx_impl->ptr_mbx_local->GetInfoElectrostaticsTimings();

    for (int i = 0; i < tmpi.size(); ++i) {
      mbxt_count[MBXT_LABELS::ELE_PERMDIP_REAL + i] += tmpi[i];
      mbxt_time[MBXT_LABELS::ELE_PERMDIP_REAL + i] += tmpd[i];
    }

    // accumulate timing info from dispersion pme

    std::vector<size_t> tmpi_d = mbx_impl->ptr_mbx_local->GetInfoDispersionCounts();
    std::vector<double> tmpd_d = mbx_impl->ptr_mbx_local->GetInfoDispersionTimings();

    for (int i = 0; i < tmpi_d.size(); ++i) {
      mbxt_count[MBXT_LABELS::DISP_PME_SETUP + i] += tmpi_d[i];
      mbxt_time[MBXT_LABELS::DISP_PME_SETUP + i] += tmpd_d[i];
    }

    delete mbx_impl->ptr_mbx_local;
  }

  // recreate main instance of MBX object

  mbx_impl->ptr_mbx = new bblock::System();
  mbx_impl->ptr_mbx_local = new bblock::System();

  // initialize all MBX instances

  if (aspc_step == aspc_step_reset) {
    mbx_impl->ptr_mbx_local->ResetDipoleHistory();
    aspc_num_hist = 0;
    aspc_step = 0;
  }

  mbx_init();
  mbx_init_local();

  if (mbx_aspc_enabled) aspc_step++;
}

void FixMBX::min_post_neighbor()
{
  post_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixMBX::setup(int vflag)
{
  mbx_get_dipoles_local();
}

/* ---------------------------------------------------------------------- */

void FixMBX::min_setup(int vflag)
{
  setup(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMBX::setup_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMBX::min_setup_pre_force(int vflag)
{
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMBX::init_storage() {}

/* ---------------------------------------------------------------------- */

void FixMBX::pre_force(int vflag)
{
  // update coordinates in MBX objects

  if (has_gcmc) { post_neighbor(); }
  mbx_update_xyz();
  mbx_update_xyz_local();
}

/* ---------------------------------------------------------------------- */

void FixMBX::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMBX::setup_pre_exchange()
{
  pre_exchange();
}

/* ---------------------------------------------------------------------- */

void FixMBX::pre_exchange()
{
  if (!mbx_aspc_enabled) return;

  // save copy of dipole history

  aspc_num_hist = mbx_impl->ptr_mbx_local->GetNumDipoleHistory();

  //  printf("# of histories= %i\n",aspc_num_hist);

  if (aspc_num_hist > aspc_max_num_hist)
    error->all(FLERR, "[MBX] Inconsistent # of ASPC histories");

  // induced dipole history does not include additional sites (e.g. water's M-site)

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  tagint *tag = atom->tag;
  double **x = atom->x;

  if (mbx_num_atoms_local == 0) { return; }

  for (int h = 0; h < aspc_num_hist; ++h) {
    std::vector<double> mbx_dip_history = mbx_impl->ptr_mbx_local->GetDipoleHistory(h);


    int indx = 0;
    for (int i = 0; i < nall; ++i) {
      if (mol_anchor[i] && mol_local[i]) {
        const int mtype = mol_type[i];

        // to be replaced with integer comparison

        bool include_monomer = true;
        tagint anchor = atom->tag[i];

        // this will save history for both local and ghost particles
        // comm->exchange() will sync ghost histories w/ local particles in new decomposition

        bool is_ext = false;
        int na = get_include_monomer(mol_names[mtype], anchor, include_monomer, is_ext);

        // add info

        if (include_monomer) {
          for (int j = 0; j < na; ++j) {
            const int ii = atom->map(anchor + j);
            aspc_dip_hist[ii][h * 3] = mbx_dip_history[indx++];
            aspc_dip_hist[ii][h * 3 + 1] = mbx_dip_history[indx++];
            aspc_dip_hist[ii][h * 3 + 2] = mbx_dip_history[indx++];
          }
        }
      }    // if(anchor)

    }    // for(nall)

  }    // for(num_hist)

  // pack dipole history into arrays for exchange
}

/* ---------------------------------------------------------------------- */

void FixMBX::post_force(int vflag)
{
  if (!print_dipoles) return;

  mbx_get_dipoles_local();
}

void FixMBX::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMBX::mbx_get_dipoles_local()
{
  if (!print_dipoles) return;

  // only need to grab dipoles on step where output is written

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  tagint *tag = atom->tag;

  // conversion factor for e*Anstrom --> Debye
  // const double qe_Debye = 1.0 / 0.2081943;

  // zero dipole array

  for (int i = 0; i < atom->nmax; ++i)
    for (int j = 0; j < 9; ++j) mbx_dip[i][j] = 0.0;


  {
    std::vector<double> mu_perm;
    std::vector<double> mu_ind;
    std::vector<double> mu_tot;

    mbx_impl->ptr_mbx_local->GetMolecularDipoles(mu_perm, mu_ind);

    int indx = 0;
    for (int i = 0; i < nall; ++i) {
      if (mol_anchor[i] && mol_local[i]) {
        const int mtype = mol_type[i];

        // to be replaced with integer comparison

        bool include_monomer = true;
        tagint anchor = atom->tag[i];

        // this will save history for both local and ghost particles
        // comm->exchange() will sync ghost histories w/ local particles in new decomposition

        bool is_ext = false;
        int na = get_include_monomer(mol_names[mtype], anchor, include_monomer, is_ext);

        // add info

        if (include_monomer) {
          for (int j = 0; j < 1; ++j) {
            const int ii = atom->map(anchor + j);
            mbx_dip[ii][0] = mu_perm[indx * 3];
            mbx_dip[ii][1] = mu_perm[indx * 3 + 1];
            mbx_dip[ii][2] = mu_perm[indx * 3 + 2];

            mbx_dip[ii][3] = mu_ind[indx * 3];
            mbx_dip[ii][4] = mu_ind[indx * 3 + 1];
            mbx_dip[ii][5] = mu_ind[indx * 3 + 2];

            mbx_dip[ii][6] = mbx_dip[ii][0] + mbx_dip[ii][3];
            mbx_dip[ii][7] = mbx_dip[ii][1] + mbx_dip[ii][4];
            mbx_dip[ii][8] = mbx_dip[ii][2] + mbx_dip[ii][5];

            indx++;
          }
        }
      }    // if(anchor)

    }    // for(nall)
  }

  array_atom = mbx_dip;
}

/* ---------------------------------------------------------------------- */

int FixMBX::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < aspc_per_atom_size; ++j) buf[m++] = aspc_dip_hist[list[i]][j];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixMBX::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  for (int i = first; i < first + n; ++i) {
    for (int j = 0; j < aspc_per_atom_size; ++j) aspc_dip_hist[i][j] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

void FixMBX::grow_arrays(int nmax)
{
  memory->grow(mol_type, nmax, "fixmbx:mol_type");
  memory->grow(mol_anchor, nmax, "fixmbx:mol_anchor");
  memory->grow(mol_local, nmax, "fixmbx:mol_local");

  if (mbx_aspc_enabled)
    memory->grow(aspc_dip_hist, nmax, aspc_per_atom_size, "fixmbx:mbx_dip_hist");

  if (print_dipoles) memory->grow(mbx_dip, nmax, 9, "fixmbx:mbx_dip");
}

/* ----------------------------------------------------------------------
   copy per-atom values
------------------------------------------------------------------------- */

void FixMBX::copy_arrays(int i, int j, int /*delflag*/)
{
  mol_type[j] = mol_type[i];
  mol_anchor[j] = mol_anchor[i];
  mol_local[j] = mol_local[i];

  if (mbx_aspc_enabled)
    for (int k = 0; k < aspc_per_atom_size; ++k) aspc_dip_hist[j][k] = aspc_dip_hist[i][k];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixMBX::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = mol_type[i];
  buf[n++] = mol_anchor[i];
  buf[n++] = mol_local[i];

  if (mbx_aspc_enabled)
    for (int j = 0; j < aspc_per_atom_size; ++j) buf[n++] = aspc_dip_hist[i][j];

  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixMBX::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  mol_type[nlocal] = buf[n++];
  mol_anchor[nlocal] = buf[n++];
  mol_local[nlocal] = buf[n++];

  if (mbx_aspc_enabled)
    for (int j = 0; j < aspc_per_atom_size; ++j) aspc_dip_hist[nlocal][j] = buf[n++];

  return n;
}

/* ----------------------------------------------------------------------
   Initialize new MBX instance with all molecules
------------------------------------------------------------------------- */

void FixMBX::mbx_init()
{
  mbxt_start(MBXT_LABELS::INIT);

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  tagint *tag = atom->tag;
  double **x = atom->x;
  double *q = atom->q;

  mbx_num_atoms = 0;
  mbx_num_ext = 0;

  std::vector<size_t> molec;

  double ximage[3];

  int nm = 0;

  const double xlo = domain->boxlo[0];
  const double ylo = domain->boxlo[1];
  const double zlo = domain->boxlo[2];

  std::vector<double> xyz_ext;
  std::vector<double> chg_ext;
  std::vector<size_t> islocal_ext;
  std::vector<int> tag_ext;

  // loop over all atoms on proc (local + ghost)

  for (int i = 0; i < nall; ++i) {
    // if anchor-atom, then add to MBX (currently assume molecule is intact)

    if (mol_anchor[i]) {
      std::vector<std::string> names;
      std::vector<double> xyz;

      const int mtype = mol_type[i];

      int is_local = (i < nlocal);

      bool is_ext = false;
      int na = get_num_atoms_per_monomer(mol_names[mtype], is_ext);

      // ids of particles in molecule on proc

      tagint anchor = tag[i];

      int amap[_MAX_ATOMS_PER_MONOMER];
      bool add_monomer = true;
      for (int j = 1; j < na; ++j) {
        amap[j] = atom->map(anchor + j);
        if (amap[j] == -1) add_monomer = false;
      }

      // test if external charged particle
      if (strcmp("dp1", mol_names[mtype]) == 0) {
        add_monomer = false;

        xyz_ext.push_back(x[i][0] - xlo);
        xyz_ext.push_back(x[i][1] - ylo);
        xyz_ext.push_back(x[i][2] - zlo);
        chg_ext.push_back(q[i]);
        islocal_ext.push_back(is_local);
        tag_ext.push_back(anchor);

        mbx_num_ext++;

        // add info for monomer

      } else if (add_monomer) {
        // add coordinates

        xyz.push_back(x[i][0] - xlo);
        xyz.push_back(x[i][1] - ylo);
        xyz.push_back(x[i][2] - zlo);

        for (int j = 1; j < na; ++j) {
          domain->closest_image(x[i], x[amap[j]], ximage);
          xyz.push_back(ximage[0] - xlo);
          xyz.push_back(ximage[1] - ylo);
          xyz.push_back(ximage[2] - zlo);
        }

        // add types
        add_monomer_atom_types(mol_names[mtype], names);

        // add monomer

        molec.push_back(nm);
        nm++;

        mbx_impl->ptr_mbx->AddMonomer(xyz, names, mol_names[mtype], is_local, anchor);
        mbx_impl->ptr_mbx->AddMolecule(molec);

        mbx_num_atoms += na;
      }

    }    // if(mol_anchor)

  }    // for(i<nall)

  if (mbx_num_atoms + mbx_num_ext == 0) {
    mbxt_stop(MBXT_LABELS::INIT);

    return;
  }

  int *pg = comm->procgrid;
  mbx_impl->ptr_mbx->SetMPI(world, pg[0], pg[1], pg[2]);

  // set MBX solvers

  if (use_json) {
    mbx_impl->ptr_mbx->SetUpFromJson(json_settings);

    // make sure cutoffs are consistent

    double mbx_cut = mbx_impl->ptr_mbx->GetRealspaceCutoff();
    double diff_sq = (mbx_cut - pair_mbx->cut_global) * (mbx_cut - pair_mbx->cut_global);
    if (diff_sq > MBX_CUTOFF_TOL) error->one(FLERR, "[MBX] cutoff not consistent with LAMMPS");
    double mbx_2b_cut = mbx_impl->ptr_mbx->Get2bCutoff();
    if (mbx_2b_cut > mbx_cut)
      error->one(FLERR,
                 "[MBX] 2-body PIP cutoff must be less than or equal to realspace cutoff. (This "
                 "may be changed in a future release.)");
    double mbx_3b_cut = mbx_impl->ptr_mbx->Get3bCutoff();
    if (mbx_3b_cut > mbx_cut)
      error->one(FLERR,
                 "[MBX] 3-body PIP cutoff must be less than or equal to realspace cutoff. (This "
                 "may be changed in a future release.)");
    double mbx_4b_cut = mbx_impl->ptr_mbx->Get4bCutoff();
    if (mbx_4b_cut > mbx_cut)
      error->one(FLERR,
                 "[MBX] 4-body PIP cutoff must be less than or equal to realspace cutoff. (This "
                 "may be changed in a future release.)");
  } else {
    mbx_impl->ptr_mbx->SetRealspaceCutoff(pair_mbx->cut_global);
    mbx_impl->ptr_mbx->SetUpFromJson();
  }

  // load external charged particles
  if (mbx_num_ext > 0) {
    mbx_impl->ptr_mbx->SetExternalChargesAndPositions(chg_ext, xyz_ext, islocal_ext, tag_ext);
  }

  // setup MBX solver(s); need to keep pbc turned off, which currently disables electrostatic solver

  std::vector<double> box;
  mbx_impl->ptr_mbx->SetPBC(box);


  // check for incompatible pair styles
  // electrostatics should be handled entirely by MBX
  Pair *pairstyles_coullong = force->pair_match(".*coul/long.*", 0);
  Pair *pairstyles_coulcut = force->pair_match(".*coul/cut.*", 0);
  Pair *pairstyles_coulexclude = force->pair_match("coul/exclude", 0);

  if (!pairstyles_coulexclude && mbx_num_ext > 0)
    error->warning(FLERR,
                   "[MBX] dp1 monomers present, but coul/exclude pair style not found. If using "
                   "special_bonds, please include coul/exclude: ");
  if (pairstyles_coulcut) {
    error->all(FLERR,
                   "[MBX] Incompatible coul/cut pair style: coulombic interactions should be "
                   "handled internally by MBX: ");
  }
  if (pairstyles_coullong) {
    error->all(FLERR,
                   "[MBX] Incompatible  coul/long pair style: coulombic interactions should be "
                   "handled internally by MBX: ");
  }

  mbxt_stop(MBXT_LABELS::INIT);
}

/* ----------------------------------------------------------------------
   Initialize new MBX instance with all molecules that have local atom
------------------------------------------------------------------------- */

void FixMBX::mbx_init_local()
{
  mbxt_start(MBXT_LABELS::INIT_LOCAL);

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  tagint *tag = atom->tag;
  double **x = atom->x;
  double *q = atom->q;

  mbx_num_atoms_local = 0;
  mbx_num_ext_local = 0;

  for (int i = 0; i < nall; ++i) {
    mol_local[i] = 0;
    if (mol_anchor[i]) mol_local[i] = 1;
  }

  // remove ghost monomers that are periodic images of local monomer
  // -- just an artifact for small systems and PBC

  for (int i = nlocal; i < nall; ++i) {
    if (mol_anchor[i]) {
      const int indx = atom->map(atom->tag[i]);
      if (indx < nlocal && mol_local[indx]) mol_local[i] = 0;
    }
  }

  // remove ghost monomers that are periodic images of ghost monomer
  // -- just an artifact for small systems and PBC

  // if two ghost have the same tag, then they are periodic
  // -- images of each other, so remove the second one

  for (int i = nlocal; i < nall - 1; ++i) {
    if (mol_anchor[i] && mol_local[i]) {
      for (int j = i + 1; j < nall; ++j) {
        if (mol_anchor[j] && mol_local[j]) {
          if (atom->tag[i] == atom->tag[j]) {
            mol_local[j] = 0;
          }
        }
      }
    }
  }

  std::vector<size_t> molec;

  double ximage[3];
  const double xlo = domain->boxlo[0];
  const double ylo = domain->boxlo[1];
  const double zlo = domain->boxlo[2];

  std::vector<double> xyz_ext;
  std::vector<double> chg_ext;
  std::vector<size_t> islocal_ext;
  std::vector<int> tag_ext;

  // loop over all atoms on proc

  int nm = 0;

  for (int i = 0; i < nall; ++i) {
    // if anchor-atom, then add to MBX (currently assume molecule is intact)

    if (mol_anchor[i] && mol_local[i]) {
      std::vector<std::string> names;
      std::vector<double> xyz;

      const int mtype = mol_type[i];

      int is_local = (i < nlocal);

      bool is_ext = false;
      int na = get_num_atoms_per_monomer(mol_names[mtype], is_ext);

      // ids of particles in molecule on proc

      tagint anchor = tag[i];

      // TODO fix this if monomer becomes too large
      int amap[_MAX_ATOMS_PER_MONOMER];
      bool add_monomer = true;
      for (int j = 1; j < na; ++j) {
        amap[j] = atom->map(anchor + j);
        if (amap[j] == -1) add_monomer = false;
      }

      // test if external charged particle
      if (strcmp("dp1", mol_names[mtype]) == 0) {
        add_monomer = false;

        xyz_ext.push_back(x[i][0] - xlo);
        xyz_ext.push_back(x[i][1] - ylo);
        xyz_ext.push_back(x[i][2] - zlo);
        chg_ext.push_back(q[i]);
        islocal_ext.push_back(is_local);
        tag_ext.push_back(anchor);

        mbx_num_ext_local++;

        // add info for monomer

      } else if (add_monomer) {
        // add coordinates

        xyz.push_back(x[i][0] - xlo);
        xyz.push_back(x[i][1] - ylo);
        xyz.push_back(x[i][2] - zlo);

        for (int j = 1; j < na; ++j) {
          domain->closest_image(x[i], x[amap[j]], ximage);
          xyz.push_back(ximage[0] - xlo);
          xyz.push_back(ximage[1] - ylo);
          xyz.push_back(ximage[2] - zlo);
        }

        add_monomer_atom_types(mol_names[mtype], names);

        molec.push_back(nm++);

        mbx_impl->ptr_mbx_local->AddMonomer(xyz, names, mol_names[mtype], is_local, anchor);
        mbx_impl->ptr_mbx_local->AddMolecule(molec);

        mbx_num_atoms_local += na;
      }

    }    // if(mol_anchor)

  }    // for(i<nall)

  // setup MPI in MBX solver

  int *pg = comm->procgrid;
  mbx_impl->ptr_mbx_local->SetMPI(world, pg[0], pg[1], pg[2]);

  // set MBX solvers

  if (use_json) {
    mbx_impl->ptr_mbx_local->SetUpFromJson(json_settings);

    // make sure cutoffs are consistent

    double mbx_cut = mbx_impl->ptr_mbx_local->GetRealspaceCutoff();
    double diff_sq = (mbx_cut - pair_mbx->cut_global) * (mbx_cut - pair_mbx->cut_global);
    if (diff_sq > MBX_CUTOFF_TOL) error->one(FLERR, "[MBX] cutoff not consistent with LAMMPS");
    double mbx_2b_cut = mbx_impl->ptr_mbx_local->Get2bCutoff();
    if (mbx_2b_cut > mbx_cut)
      error->one(FLERR,
                 "[MBX] 2-body PIP cutoff must be less than or equal to realspace cutoff. (This "
                 "may be changed in a future release.)");
    double mbx_3b_cut = mbx_impl->ptr_mbx_local->Get3bCutoff();
    if (mbx_3b_cut > mbx_cut)
      error->one(FLERR,
                 "[MBX] 3-body PIP cutoff must be less than or equal to realspace cutoff. (This "
                 "may be changed in a future release.)");
    double mbx_4b_cut = mbx_impl->ptr_mbx_local->Get4bCutoff();
    if (mbx_4b_cut > mbx_cut)
      error->one(FLERR,
                 "[MBX] 4-body PIP cutoff must be less than or equal to realspace cutoff. (This "
                 "may be changed in a future release.)");
  } else {
    mbx_impl->ptr_mbx_local->SetRealspaceCutoff(pair_mbx->cut_global);
    mbx_impl->ptr_mbx_local->SetUpFromJson();
  }

  // load external charged particles
  if (mbx_num_ext_local > 0) {
    mbx_impl->ptr_mbx_local->SetExternalChargesAndPositions(chg_ext, xyz_ext, islocal_ext, tag_ext);
  }

  std::vector<double> box;

  if (domain->nonperiodic && (domain->xperiodic || domain->yperiodic || domain->zperiodic))
    error->all(FLERR, "[MBX] System must be fully periodic or non-periodic with MBX");

  // get ewald values from MBX and verify they make sense relative to LAMMPS periodicity
  double elec_alpha, elec_grid, disp_alpha, disp_grid;
  size_t elec_spline, disp_spline;

  mbx_impl->ptr_mbx_local->GetEwaldParamsElectrostatics(elec_alpha, elec_grid, elec_spline);
  mbx_impl->ptr_mbx_local->GetEwaldParamsDispersion(disp_alpha, disp_grid, disp_spline);

  if ((elec_alpha > 0.0) && (!domain->xperiodic && !domain->yperiodic && !domain->zperiodic))
    error->all(FLERR,
               "[MBX] Electrostatic Ewald parameters set (alpha_ewald_elec = " +
                   std::to_string(elec_alpha) + "), but system is not periodic");
  if ((disp_alpha > 0.0) && (!domain->xperiodic && !domain->yperiodic && !domain->zperiodic))
    error->all(FLERR,
               "[MBX] Dispersion Ewald parameters set (alpha_ewald_disp = " +
                   std::to_string(disp_alpha) + "), but system is not periodic");
  if ((elec_alpha == 0.0 || disp_alpha == 0.0) &&
      (domain->xperiodic || domain->yperiodic || domain->zperiodic))
    error->warning(FLERR, "[MBX] System is periodic, but Ewald alpha parameters not set");

  box = std::vector<double>(9, 0.0);

  box[0] = domain->xprd;

  box[3] = domain->xy;
  box[4] = domain->yprd;

  box[6] = domain->xz;
  box[7] = domain->yz;
  box[8] = domain->zprd;

  mbx_impl->ptr_mbx_local->SetPBC(box);
  mbx_impl->ptr_mbx_local->SetBoxPMElocal(box);

  mbx_impl->ptr_mbx_local->SetPeriodicity(!domain->nonperiodic);

  std::vector<int> egrid = mbx_impl->ptr_mbx_local->GetFFTDimensionElectrostatics(1);
  std::vector<int> dgrid =
      mbx_impl->ptr_mbx_local->GetFFTDimensionDispersion(1);    // will return mesh even for gas-phase

  if (print_verbose && first_step  && comm->me == 0) {
    std::string mbx_settings_ = mbx_impl->ptr_mbx_local->GetCurrentSystemConfig();
    if (screen) {
      fprintf(screen, "\n[MBX] 'Local' Settings\n%s\n", mbx_settings_.c_str());
      fprintf(screen, "[MBX] LOCAL electrostatics FFT grid= %i %i %i\n", egrid[0], egrid[1],
              egrid[2]);
      fprintf(screen, "[MBX] LOCAL dispersion FFT grid= %i %i %i\n", dgrid[0], dgrid[1], dgrid[2]);
    }
    if (logfile) {
      fprintf(logfile, "\n[MBX] 'Local' Settings\n%s\n", mbx_settings_.c_str());
      fprintf(logfile, "[MBX] LOCAL electrostatics FFT grid= %i %i %i\n", egrid[0], egrid[1],
              egrid[2]);
      fprintf(logfile, "[MBX] LOCAL dispersion FFT grid= %i %i %i\n", dgrid[0], dgrid[1], dgrid[2]);
    }
  }

  // check if using cg or aspc integrator for MBX dipoles

  if (first_step) {
    std::string dip_method = mbx_impl->ptr_mbx->GetDipoleMethod();

    if (dip_method == "aspc") {
      mbx_aspc_enabled = true;

      memory->create(aspc_dip_hist, atom->nmax, aspc_per_atom_size, "fixmbx::aspc_dip_hist");
    } else if (!(dip_method == "cg")) {
      error->one(FLERR, "[MBX] requested dip_method not supported with LAMMPS");
    }
  }

  if (mbx_aspc_enabled) mbx_init_dipole_history_local();

  mbxt_stop(MBXT_LABELS::INIT_LOCAL);
}


/* ----------------------------------------------------------------------
   Update MBX instance for all molecules
------------------------------------------------------------------------- */

void FixMBX::mbx_update_xyz()
{
  mbxt_start(MBXT_LABELS::UPDATE_XYZ);

  // update coordinates

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  tagint *tag = atom->tag;
  double **x = atom->x;
  double *q = atom->q;

  if (mbx_num_atoms + mbx_num_ext == 0) {
    mbxt_stop(MBXT_LABELS::UPDATE_XYZ);
    return;
  }

  const double xlo = domain->boxlo[0];
  const double ylo = domain->boxlo[1];
  const double zlo = domain->boxlo[2];

  double ximage[3];

  std::vector<double> xyz(mbx_num_atoms * 3);

  std::vector<double> xyz_ext(mbx_num_ext * 3);
  std::vector<double> chg_ext(mbx_num_ext);

  int indx = 0;
  int indx_ext = 0;
  for (int i = 0; i < nall; ++i) {
    if (mol_anchor[i]) {
      const int mtype = mol_type[i];

      bool is_ext = false;
      int na = get_num_atoms_per_monomer(mol_names[mtype], is_ext);

      // ids of particles in molecule on proc

      tagint anchor = tag[i];

      int amap[_MAX_ATOMS_PER_MONOMER];
      bool add_monomer = true;
      for (int j = 1; j < na; ++j) {
        amap[j] = atom->map(anchor + j);
        if (amap[j] == -1) add_monomer = false;
      }

      // test if external charged particle
      if (strcmp("dp1", mol_names[mtype]) == 0) {
        add_monomer = false;

        xyz_ext[indx_ext * 3] = x[i][0] - xlo;
        xyz_ext[indx_ext * 3 + 1] = x[i][1] - ylo;
        xyz_ext[indx_ext * 3 + 2] = x[i][2] - zlo;
        chg_ext[indx_ext] = q[i];

        indx_ext++;

        // add info for monomer

      } else if (add_monomer) {
        // add coordinates

        xyz[indx * 3] = x[i][0] - xlo;
        xyz[indx * 3 + 1] = x[i][1] - ylo;
        xyz[indx * 3 + 2] = x[i][2] - zlo;

        for (int j = 1; j < na; ++j) {
          domain->closest_image(x[i], x[amap[j]], ximage);
          xyz[(indx + j) * 3] = ximage[0] - xlo;
          xyz[(indx + j) * 3 + 1] = ximage[1] - ylo;
          xyz[(indx + j) * 3 + 2] = ximage[2] - zlo;
        }

        indx += na;
      }

    }    // if(mol_anchor)

  }    // for(i<nall)

  if (xyz.size() != indx * 3) error->one(FLERR, "Inconsistent # of atoms");
  mbx_impl->ptr_mbx->SetRealXyz(xyz);

  if (xyz_ext.size() != indx_ext * 3) error->one(FLERR, "Inconsistent # of external charges");
  if (mbx_num_ext > 0) { mbx_impl->ptr_mbx->SetExternalChargesAndPositions(chg_ext, xyz_ext); }

  mbxt_stop(MBXT_LABELS::UPDATE_XYZ);
}

/* ----------------------------------------------------------------------
   Update MBX instance for local molecules + plus halo
------------------------------------------------------------------------- */

void FixMBX::mbx_update_xyz_local()
{
  mbxt_start(MBXT_LABELS::UPDATE_XYZ_LOCAL);

  // update if box changes
  // need to update box passed to PME solver

  if (domain->box_change) {
    std::vector<double> box;

    if (!domain->nonperiodic) {
      box = std::vector<double>(9, 0.0);

      box[0] = domain->xprd;

      box[3] = domain->xy;
      box[4] = domain->yprd;

      box[6] = domain->xz;
      box[7] = domain->yz;
      box[8] = domain->zprd;

    } else if (domain->xperiodic || domain->yperiodic || domain->zperiodic)
      error->all(FLERR, "[MBX] System must be fully periodic or non-periodic with MBX");

    // get ewald values from MBX and verify they make sense relative to LAMMPS periodicity
    double elec_alpha, elec_grid, disp_alpha, disp_grid;
    size_t elec_spline, disp_spline;

    mbx_impl->ptr_mbx_local->GetEwaldParamsElectrostatics(elec_alpha, elec_grid, elec_spline);
    mbx_impl->ptr_mbx_local->GetEwaldParamsDispersion(disp_alpha, disp_grid, disp_spline);

    if ((elec_alpha > 0.0) && (!domain->xperiodic || !domain->yperiodic || !domain->zperiodic))
      error->all(FLERR,
                 "[MBX] Electrostatic Ewald parameters set (alpha_ewald_elec = " +
                     std::to_string(elec_alpha) + "), but system is not periodic");
    if ((disp_alpha > 0.0) && (!domain->xperiodic || !domain->yperiodic || !domain->zperiodic))
      error->all(FLERR,
                 "[MBX] Dispersion Ewald parameters set (alpha_ewald_disp = " +
                     std::to_string(disp_alpha) + "), but system is not periodic");
    if ((elec_alpha == 0.0 || disp_alpha == 0.0) &&
        (!domain->xperiodic || !domain->yperiodic || !domain->zperiodic))
      error->warning(FLERR, "[MBX] System is periodic, but Ewald alpha parameters not set");

    mbx_impl->ptr_mbx_local->SetPBC(box);
    mbx_impl->ptr_mbx_local->SetBoxPMElocal(box);
  }

  // update coordinates

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  tagint *tag = atom->tag;
  double **x = atom->x;
  double *q = atom->q;

  if (mbx_num_atoms_local + mbx_num_ext_local == 0) {
    mbxt_stop(MBXT_LABELS::UPDATE_XYZ);
    return;
  }

  const double xlo = domain->boxlo[0];
  const double ylo = domain->boxlo[1];
  const double zlo = domain->boxlo[2];

  double ximage[3];

  std::vector<double> xyz(mbx_num_atoms_local * 3);

  std::vector<double> xyz_ext(mbx_num_ext_local * 3);
  std::vector<double> chg_ext(mbx_num_ext_local);

  int indx = 0;
  int indx_ext = 0;
  for (int i = 0; i < nall; ++i) {
    if (mol_anchor[i] && mol_local[i]) {
      const int mtype = mol_type[i];

      bool is_ext = false;
      int na = get_num_atoms_per_monomer(mol_names[mtype], is_ext);

      // ids of particles in molecule on proc

      tagint anchor = tag[i];

      int amap[_MAX_ATOMS_PER_MONOMER];
      bool add_monomer = true;
      for (int j = 1; j < na; ++j) {
        amap[j] = atom->map(anchor + j);
        if (amap[j] == -1) add_monomer = false;
      }

      // test if external charged particle
      if (strcmp("dp1", mol_names[mtype]) == 0) {
        add_monomer = false;

        xyz_ext[indx_ext * 3] = x[i][0] - xlo;
        xyz_ext[indx_ext * 3 + 1] = x[i][1] - ylo;
        xyz_ext[indx_ext * 3 + 2] = x[i][2] - zlo;
        chg_ext[indx_ext] = q[i];

        indx_ext++;

        // add info for monomer

      } else if (add_monomer) {
        // add coordinates

        xyz[indx * 3] = x[i][0] - xlo;
        xyz[indx * 3 + 1] = x[i][1] - ylo;
        xyz[indx * 3 + 2] = x[i][2] - zlo;

        for (int j = 1; j < na; ++j) {
          domain->closest_image(x[i], x[amap[j]], ximage);
          xyz[(indx + j) * 3] = ximage[0] - xlo;
          xyz[(indx + j) * 3 + 1] = ximage[1] - ylo;
          xyz[(indx + j) * 3 + 2] = ximage[2] - zlo;
        }

        indx += na;
      }

    }    // if(mol_anchor)

  }    // for(i<nall)

  if (xyz.size() != indx * 3) error->one(FLERR, "Inconsistent # of atoms");
  mbx_impl->ptr_mbx_local->SetRealXyz(xyz);

  if (xyz_ext.size() != indx_ext * 3) error->one(FLERR, "Inconsistent # of external charges");
  if (mbx_num_ext_local > 0) { mbx_impl->ptr_mbx_local->SetExternalChargesAndPositions(chg_ext, xyz_ext); }

  mbxt_stop(MBXT_LABELS::UPDATE_XYZ_LOCAL);
}


/* ----------------------------------------------------------------------
   Initialize dipole history for local molecules + plus halo
------------------------------------------------------------------------- */

void FixMBX::mbx_init_dipole_history_local()
{

  if (aspc_num_hist == 0) return;

  // sync dipole histories of ghost particles

  comm->forward_comm(this);

  // update coordinates

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  tagint *tag = atom->tag;
  double **x = atom->x;

  if (mbx_num_atoms_local == 0) {
    return;
  }

  const double xlo = domain->boxlo[0];
  const double ylo = domain->boxlo[1];
  const double zlo = domain->boxlo[2];

  double ximage[3];

  mbx_impl->ptr_mbx_local->SetNumDipoleHistory(aspc_num_hist);

  std::vector<double> mbx_dip_history = std::vector<double>(mbx_num_atoms_local * 3);

  // following debug only works if all ranks contribute

  for (int h = 0; h < aspc_num_hist; ++h) {

    int indx = 0;
    for (int i = 0; i < nall; ++i) {
      if (mol_anchor[i] && mol_local[i]) {
        const int mtype = mol_type[i];

        bool is_ext = false;
        int na = get_num_atoms_per_monomer(mol_names[mtype], is_ext);

        // ids of particles in molecule on proc

        tagint anchor = tag[i];

        int amap[_MAX_ATOMS_PER_MONOMER];
        bool add_monomer = true;
        for (int j = 1; j < na; ++j) {
          amap[j] = atom->map(anchor + j);
          if (amap[j] == -1) add_monomer = false;
        }

        // add info

        if (add_monomer) {
          // add coordinates

          for (int j = 0; j < na; ++j) {
            const int ii = atom->map(anchor + j);
            mbx_dip_history[indx++] = aspc_dip_hist[ii][h * 3];
            mbx_dip_history[indx++] = aspc_dip_hist[ii][h * 3 + 1];
            mbx_dip_history[indx++] = aspc_dip_hist[ii][h * 3 + 2];

          }    // for(na)

        }    // if(add_monomer)

      }    // if(mol_anchor)

    }    // for(i<nall)

    if (mbx_num_atoms_local * 3 != indx) error->one(FLERR, "Inconsistent # of atoms");
    //      printf("calling SetDipoleHistory");
    mbx_impl->ptr_mbx_local->SetDipoleHistory(h, mbx_dip_history);

  }    // for(hist)

}

/* ----------------------------------------------------------------------
   Helper functions for timing
------------------------------------------------------------------------- */

void FixMBX::mbxt_start(int T)
{
  mbxt_time_start[T] = MPI_Wtime();
}

void FixMBX::mbxt_stop(int T)
{
  mbxt_time[T] += MPI_Wtime() - mbxt_time_start[T];
  mbxt_count[T]++;
}

void FixMBX::mbxt_print_time(const char *name, int T, double *d)
{
  double tavg = d[T];
  double tmin = d[MBXT_LABELS::NUM_TIMERS + T];
  double tmax = d[MBXT_LABELS::NUM_TIMERS * 2 + T];

  double p = tmax / d[MBXT_LABELS::NUM_TIMERS * 3] * 100.0;

  if (screen)
    fprintf(screen, "[MBX] %-20s:  %12.5g  %12.5g  %12.5g  %8i %8.2f%%\n", name, tmin, tavg, tmax,
            mbxt_count[T], p);

  if (logfile)
    fprintf(logfile, "[MBX] %-20s:  %12.5g  %12.5g  %12.5g  %8i %8.2f%%\n", name, tmin, tavg, tmax,
            mbxt_count[T], p);
}

void FixMBX::mbxt_write_summary()
{
  double t[MBXT_LABELS::NUM_TIMERS * 3 + 1];
  double *tavg = &t[0];
  double *tmin = &t[MBXT_LABELS::NUM_TIMERS];
  double *tmax = &t[MBXT_LABELS::NUM_TIMERS * 2];

  // total runtime since fix created

  t[MBXT_LABELS::NUM_TIMERS * 3] = MPI_Wtime() - mbxt_initial_time;

  MPI_Reduce(mbxt_time, tavg, MBXT_LABELS::NUM_TIMERS, MPI_DOUBLE, MPI_SUM, 0, world);
  MPI_Reduce(mbxt_time, tmin, MBXT_LABELS::NUM_TIMERS, MPI_DOUBLE, MPI_MIN, 0, world);
  MPI_Reduce(mbxt_time, tmax, MBXT_LABELS::NUM_TIMERS, MPI_DOUBLE, MPI_MAX, 0, world);

  if (me) return;

  for (int i = 0; i < MBXT_LABELS::NUM_TIMERS; ++i) tavg[i] /= (double) nprocs;

  if (screen) {
    fprintf(screen, "\n[MBX] Total MBX fix/pair time= %f seconds\n", t[MBXT_LABELS::NUM_TIMERS * 3]);
    fprintf(screen, "[MBX] Timing Summary\n");
    fprintf(screen,
            "[MBX] kernel                      tmin          tavg          tmax         count   "
            "%%total\n");
    fprintf(
        screen,
        "[MBX] "
        "-----------------------------------------------------------------------------------\n");
  }
  if (logfile) {
    fprintf(logfile, "\n[MBX] Total MBX fix/pair time= %f seconds\n", t[MBXT_LABELS::NUM_TIMERS * 3]);
    fprintf(logfile, "[MBX] Timing Summary\n");
    fprintf(logfile,
            "[MBX] kernel                      tmin          tavg          tmax         count   "
            "%%total\n");
    fprintf(
        logfile,
        "[MBX] "
        "-----------------------------------------------------------------------------------\n");
  }

  mbxt_print_time("INIT", MBXT_LABELS::INIT, t);
  mbxt_print_time("UPDATE_XYZ", MBXT_LABELS::UPDATE_XYZ, t);
  mbxt_print_time("ACCUMULATE_F", MBXT_LABELS::ACCUMULATE_F, t);

  mbxt_print_time("E1B", MBXT_LABELS::E1B, t);
  mbxt_print_time("E2B", MBXT_LABELS::E2B_GHOST, t);
  mbxt_print_time("E3B", MBXT_LABELS::E3B_GHOST, t);
  mbxt_print_time("E4B", MBXT_LABELS::E4B_GHOST, t);
  mbxt_print_time("DISP", MBXT_LABELS::DISP, t);
  mbxt_print_time("DISP_PME", MBXT_LABELS::DISP_PME, t);
  mbxt_print_time("BUCK", MBXT_LABELS::BUCK, t);
  mbxt_print_time("ELE", MBXT_LABELS::ELE, t);
  mbxt_print_time("INIT_LOCAL", MBXT_LABELS::INIT_LOCAL, t);
  mbxt_print_time("UPDATE_XYZ_LOCAL", MBXT_LABELS::UPDATE_XYZ_LOCAL, t);
  mbxt_print_time("ACCUMULATE_F_LOCAL", MBXT_LABELS::ACCUMULATE_F_LOCAL, t);

  if (screen) {
    fprintf(screen, "\n\n[MBX] Electrostatics Summary\n");
    fprintf(screen,
            "[MBX] kernel                      tmin          tavg          tmax         count   "
            "%%total\n");
    fprintf(
        screen,
        "[MBX] "
        "-----------------------------------------------------------------------------------\n");
  }
  if (logfile) {
    fprintf(logfile, "\n\n[MBX] Electrostatics Summary\n");
    fprintf(logfile,
            "[MBX] kernel                      tmin          tavg          tmax         count   "
            "%%total\n");
    fprintf(
        logfile,
        "[MBX] "
        "-----------------------------------------------------------------------------------\n");
  }

  mbxt_print_time("ELE_PERMDIP_REAL", MBXT_LABELS::ELE_PERMDIP_REAL, t);
  mbxt_print_time("ELE_PERMDIP_PME", MBXT_LABELS::ELE_PERMDIP_PME, t);

  mbxt_print_time("ELE_DIPFIELD_REAL", MBXT_LABELS::ELE_DIPFIELD_REAL, t);
  mbxt_print_time("ELE_DIPFIELD_PME", MBXT_LABELS::ELE_DIPFIELD_PME, t);

  mbxt_print_time("ELE_GRAD_REAL", MBXT_LABELS::ELE_GRAD_REAL, t);
  mbxt_print_time("ELE_GRAD_PME", MBXT_LABELS::ELE_GRAD_PME, t);
  mbxt_print_time("ELE_GRAD_FIN", MBXT_LABELS::ELE_GRAD_FIN, t);

  mbxt_print_time("ELE_PME_SETUP", MBXT_LABELS::ELE_PME_SETUP, t);
  mbxt_print_time("ELE_PME_C", MBXT_LABELS::ELE_PME_C, t);
  mbxt_print_time("ELE_PME_D", MBXT_LABELS::ELE_PME_D, t);
  mbxt_print_time("ELE_PME_E", MBXT_LABELS::ELE_PME_E, t);

  mbxt_print_time("DISP_PME_SETUP", MBXT_LABELS::DISP_PME_SETUP, t);
  mbxt_print_time("DISP_PME_E", MBXT_LABELS::DISP_PME_E, t);

  mbxt_print_time("ELE_COMM_REVFOR", MBXT_LABELS::ELE_COMM_REVFOR, t);
  mbxt_print_time("ELE_COMM_REVSET", MBXT_LABELS::ELE_COMM_REVSET, t);
  mbxt_print_time("ELE_COMM_REV", MBXT_LABELS::ELE_COMM_REV, t);
  mbxt_print_time("ELE_COMM_FORSET", MBXT_LABELS::ELE_COMM_FORSET, t);
  mbxt_print_time("ELE_COMM_FOR", MBXT_LABELS::ELE_COMM_FOR, t);
}

/* ----------------------------------------------------------------------
   Helper functions for monomers
------------------------------------------------------------------------- */

int FixMBX::get_num_atoms_per_monomer(char *name, bool &inc_e)
{
  int na;
  inc_e = false;

  if (strcmp("h2o", name) == 0)
    na = 3;
  else if (strcmp("li+", name) == 0)
    na = 1;
  else if (strcmp("na+", name) == 0)
    na = 1;
  else if (strcmp("k+", name) == 0)
    na = 1;
  else if (strcmp("rb+", name) == 0)
    na = 1;
  else if (strcmp("cs+", name) == 0)
    na = 1;
  else if (strcmp("dp1", name) == 0) {
    na = 1;
    inc_e = true;
  } else if (strcmp("f-", name) == 0)
    na = 1;
  else if (strcmp("cl-", name) == 0)
    na = 1;
  else if (strcmp("br-", name) == 0)
    na = 1;
  else if (strcmp("i-", name) == 0)
    na = 1;
  else if (strcmp("co2", name) == 0)
    na = 3;
  else if (strcmp("ch4", name) == 0)
    na = 5;
  else if (strcmp("he", name) == 0)
    na = 1;
  else if (strcmp("ar", name) == 0)
    na = 1;
  else if (strcmp("h2", name) == 0)
    na = 2;
  else if (strcmp("n2o5", name) == 0)
    na = 7;
  else if (strcmp("so4a", name) == 0)
    na = 5;
  else if (strcmp("co3a", name) == 0)
    na = 4;
  else if (strcmp("no3a", name) == 0)
    na = 4;
  else if (strcmp("dp2", name) == 0)
    na = 2;
  else
    error->one(FLERR, "Unsupported molecule type in MBX");

  return na;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixMBX::get_include_monomer(char *name, int anchor, bool &inc, bool &inc_e)
{
  inc = true;
  inc_e = false;
  int na = get_num_atoms_per_monomer(name, inc_e);

  for (int ii = 1; ii < na; ii++) {
    if (atom->map(anchor + ii) < 0) {
      inc = false;
      break;
    }
  }

  if (strcmp("dp1", name) == 0) {
    inc = false;
    inc_e = true;
  }

  return na;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixMBX::add_monomer_atom_types(char *name, std::vector<std::string> &n)
{
  if (strcmp("h2o", name) == 0) {
    n.push_back("O");
    n.push_back("H");
    n.push_back("H");
  } else if (strcmp("li+", name) == 0) {
    n.push_back("Li");
  } else if (strcmp("na+", name) == 0) {
    n.push_back("Na");
  } else if (strcmp("k+", name) == 0) {
    n.push_back("K");
  } else if (strcmp("rb+", name) == 0) {
    n.push_back("Rb");
  } else if (strcmp("cs+", name) == 0) {
    n.push_back("Cs");
  } else if (strcmp("f-", name) == 0) {
    n.push_back("F");
  } else if (strcmp("cl-", name) == 0) {
    n.push_back("Cl");
  } else if (strcmp("br-", name) == 0) {
    n.push_back("Br");
  } else if (strcmp("i-", name) == 0) {
    n.push_back("I");
  } else if (strcmp("he", name) == 0) {
    n.push_back("He");
  } else if (strcmp("co2", name) == 0) {
    n.push_back("C");
    n.push_back("O");
    n.push_back("O");
  } else if (strcmp("ch4", name) == 0) {
    n.push_back("C");
    n.push_back("H");
    n.push_back("H");
    n.push_back("H");
    n.push_back("H");
  } else if (strcmp("ar", name) == 0) {
    n.push_back("Ar");
  } else if (strcmp("h2", name) == 0) {
    n.push_back("H");
    n.push_back("H");
  } else if (strcmp("n2o5", name) == 0) {
    n.push_back("O");
    n.push_back("N");
    n.push_back("N");
    n.push_back("O");
    n.push_back("O");
    n.push_back("O");
    n.push_back("O");
  } else if (strcmp("so4a", name) == 0) {
    n.push_back("S");
    n.push_back("O");
    n.push_back("O");
    n.push_back("O");
    n.push_back("O");
  } else if (strcmp("co3a", name) == 0) {
    n.push_back("C");
    n.push_back("O");
    n.push_back("O");
    n.push_back("O");
  } else if (strcmp("no3a", name) == 0) {
    n.push_back("N");
    n.push_back("O");
    n.push_back("O");
    n.push_back("O");
  } else if (strcmp("dp2", name) == 0) {
    n.push_back("X");
    n.push_back("X");
  }
  else
    error->one(FLERR, "Unsupported molecule type in MBX");
}
