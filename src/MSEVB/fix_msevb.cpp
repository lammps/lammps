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
   fix msevb: Multi-State Empirical Valence Bond method.
   See doc/src/fix_msevb.rst for full documentation.

   Implementation split across files:
     fix_msevb.cpp             - framework (constructor, hooks, MPI, thermo)
     fix_msevb_topology.cpp    - reference topology, site detection, state changes
     fix_msevb_eigensystem.cpp - Hamiltonian, coupling, eigensolve, force mixing
     fix_msevb_transfer.cpp    - PE gathering, permanent transfer
     fix_msevb_superimpose.cpp - BFS graph-isomorphism for template matching

   Contributing Authors:
     Blake Armstrong (blake.armstrong@curtin.edu.au)
     Peter Spackman (peter.spackman@curtin.edu.au)
------------------------------------------------------------------------- */

#include "fix_msevb.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "json.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "universe.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <cstdio>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMSEVB::FixMSEVB(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), samerank(MPI_COMM_NULL), commbuf(nullptr), pe(nullptr), epot(nullptr),
    my_epot(nullptr), epot_request(MPI_REQUEST_NULL), hamiltonian(nullptr), H_work(nullptr),
    eigenvalues(nullptr), eigenvectors(nullptr), amplitudes(nullptr), fd_occ(nullptr),
    eigensys_nmax(0), my_nlocal(nullptr), all_nlocal(nullptr), weights(nullptr), weights_nmax(0),
    saved_forces(nullptr), excess_forces(nullptr), excess_forces_nmax(0),
    excess_forces_nmax_serial(0), excess_force_mode(0), reaction_enabled(0), nsites(0),
    nsites_parallel(0), nsites_serial(0), nstates(1), sites(nullptr), sites_nmax(0),
    chain_H_flat(nullptr), chain_X_flat(nullptr), chain_Y_flat(nullptr), chain_rxn_flat(nullptr),
    chain_glove_flat(nullptr), glove_flat(nullptr), glove_nmax(0), epot_nmax(0), ref_type(nullptr),
    ref_charge(nullptr), ref_molecule(nullptr), ref_num_bond(nullptr), ref_bond_type_flat(nullptr),
    ref_bond_atom_flat(nullptr), ref_num_angle(nullptr), ref_angle_type_flat(nullptr),
    ref_angle_atom1_flat(nullptr), ref_angle_atom2_flat(nullptr), ref_angle_atom3_flat(nullptr),
    ref_num_dihedral(nullptr), ref_dihedral_type_flat(nullptr), ref_dihedral_atom1_flat(nullptr),
    ref_dihedral_atom2_flat(nullptr), ref_dihedral_atom3_flat(nullptr),
    ref_dihedral_atom4_flat(nullptr), ref_num_improper(nullptr), ref_improper_type_flat(nullptr),
    ref_improper_atom1_flat(nullptr), ref_improper_atom2_flat(nullptr),
    ref_improper_atom3_flat(nullptr), ref_improper_atom4_flat(nullptr), ref_nspecial_flat(nullptr),
    ref_special_flat(nullptr), ref_maxspecial(0), ref_maxtag(0), ref_bond_per_atom(0),
    ref_angle_per_atom(0), ref_dihedral_per_atom(0), ref_improper_per_atom(0),
    ref_snapshot_valid(0), list(nullptr), coupling_type(COUPLING_NONE), coupling_lambda(0.0),
    coupling_zeta(0.0), coupling_v12(0.0), coupling_alpha(0.0), coupling_gamma_v(0.0),
    coupling_a(0.0), coupling_b(0.0), coupling_taper(0.0), coupling_enabled(0),
    temp_compute(nullptr), press_compute(nullptr), enumerate_product_states(0),
    fermi_dirac_enabled(0), fd_temperature(0.0), fd_RT(0.0), max_shells(1), file_flag(false),
    file_every(0), fpout(nullptr), json_init(0), reactive_group_bit(0), scf_topology(false),
    scf_max_iter(10), min_terminate(false), max_states(0)
{
  force_reneighbor = 1;

  npartitions = universe->nworlds;
  ipartition = universe->iworld;

  if (npartitions < 2) {
    auto msg = fmt::format("WARNING: fix msevb will run very slowly with only one partition. "
                           "\nWARNING run with more for better performance.");
    utils::logmesg(lmp, msg);
  }

  // ---- Parse arguments ------------------------------------------------
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "coupling") == 0) {
      iarg++;
      if (iarg >= narg)
        error->universe_all(FLERR, "Fix msevb: missing coupling style after 'coupling'");
      coupling_enabled = 1;
      // The coupling style keyword is optional and defaults to energy/gaussian.
      // Style names contain a '/' (or are the literal 'none'); anything else is
      // treated as the first coupling parameter for the default energy/gaussian
      // style and left for the parameter loop below to consume.
      if (strcmp(arg[iarg], "energy/gaussian") == 0) {
        coupling_type = COUPLING_GRIMME2015;
        iarg++;
      } else if (strcmp(arg[iarg], "geometry/gaussian") == 0) {
        coupling_type = COUPLING_RAITERI2011;
        iarg++;
      } else if (strcmp(arg[iarg], "geometry/exp") == 0) {
        coupling_type = COUPLING_VUILLEUMIER1998;
        iarg++;
      } else if (strcmp(arg[iarg], "none") == 0) {
        coupling_type = COUPLING_NONE;
        iarg++;
      } else if (strchr(arg[iarg], '/') != nullptr) {
        error->universe_all(FLERR,
                            fmt::format("Fix msevb: unknown coupling style '{}'", arg[iarg]));
      } else {
        coupling_type = COUPLING_GRIMME2015;    // default style: energy/gaussian
      }
      // Greedily consume coupling params that follow the style keyword
      while (iarg < narg) {
        if (strcmp(arg[iarg], "lambda") == 0) {
          if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for lambda");
          coupling_lambda = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
          iarg += 2;
        } else if (strcmp(arg[iarg], "zeta") == 0) {
          if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for zeta");
          coupling_zeta = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
          iarg += 2;
        } else if (strcmp(arg[iarg], "v12") == 0) {
          if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for v12");
          coupling_v12 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
          iarg += 2;
        } else if (strcmp(arg[iarg], "alpha") == 0) {
          if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for alpha");
          coupling_alpha = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
          iarg += 2;
        } else if (strcmp(arg[iarg], "gamma") == 0) {
          if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for gamma");
          coupling_gamma_v = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
          iarg += 2;
        } else if (strcmp(arg[iarg], "a") == 0) {
          if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for a");
          coupling_a = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
          iarg += 2;
        } else if (strcmp(arg[iarg], "b") == 0) {
          if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for b");
          coupling_b = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
          iarg += 2;
        } else if (strcmp(arg[iarg], "taper") == 0) {
          if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for taper");
          coupling_taper = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
          iarg += 2;
        } else {
          break;    // not a coupling param, return to outer loop
        }
      }
    } else if (strcmp(arg[iarg], "reaction") == 0) {
      // New syntax: reaction pre.mol post.mol react.map cutoff
      if (iarg + 4 >= narg)
        error->universe_all(FLERR,
                            "Fix msevb: reaction requires 4 args: "
                            "pre_mol post_mol mapfile cutoff");
      ReactionDef rd;
      rd.pre_mol_id = arg[iarg + 1];
      rd.post_mol_id = arg[iarg + 2];
      rd.map_file_path = arg[iarg + 3];
      double cutoff = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      rd.cutoff_sq = cutoff * cutoff;
      // type_H, type_Y, glove_n, ibonding etc. filled in init()
      rxndefs.push_back(rd);
      reaction_enabled = 1;
      iarg += 5;
      // Optional per-reaction keywords: coupling, taper, offset, shells
      // (accepted in any order)
      while (iarg < narg) {
        if (strcmp(arg[iarg], "coupling") == 0) {
          iarg++;
          auto &rr = rxndefs.back();
          rr.coupling_set = true;
          if (iarg >= narg)
            error->universe_all(FLERR,
                                "Fix msevb: missing coupling style after reaction 'coupling'");
          // Style keyword is optional and defaults to energy/gaussian (see global
          // coupling parse above for the disambiguation rule).
          if (strcmp(arg[iarg], "energy/gaussian") == 0) {
            rr.coupling_type = COUPLING_GRIMME2015;
            iarg++;
          } else if (strcmp(arg[iarg], "geometry/gaussian") == 0) {
            rr.coupling_type = COUPLING_RAITERI2011;
            iarg++;
          } else if (strcmp(arg[iarg], "geometry/exp") == 0) {
            rr.coupling_type = COUPLING_VUILLEUMIER1998;
            iarg++;
          } else if (strcmp(arg[iarg], "none") == 0) {
            rr.coupling_type = COUPLING_NONE;
            iarg++;
          } else if (strchr(arg[iarg], '/') != nullptr) {
            error->universe_all(
                FLERR,
                fmt::format("Fix msevb: unknown per-reaction coupling style '{}'", arg[iarg]));
          } else {
            rr.coupling_type = COUPLING_GRIMME2015;    // default style: energy/gaussian
          }
          // Parse coupling parameters for this reaction
          while (iarg < narg) {
            if (strcmp(arg[iarg], "lambda") == 0) {
              if (iarg + 1 >= narg)
                error->universe_all(FLERR, "Fix msevb: missing value for lambda");
              rr.coupling_lambda = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
              iarg += 2;
            } else if (strcmp(arg[iarg], "zeta") == 0) {
              if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for zeta");
              rr.coupling_zeta = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
              iarg += 2;
            } else if (strcmp(arg[iarg], "v12") == 0) {
              if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for v12");
              rr.coupling_v12 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
              iarg += 2;
            } else if (strcmp(arg[iarg], "alpha") == 0) {
              if (iarg + 1 >= narg)
                error->universe_all(FLERR, "Fix msevb: missing value for alpha");
              rr.coupling_alpha = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
              iarg += 2;
            } else if (strcmp(arg[iarg], "gamma") == 0) {
              if (iarg + 1 >= narg)
                error->universe_all(FLERR, "Fix msevb: missing value for gamma");
              rr.coupling_gamma_v = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
              iarg += 2;
            } else if (strcmp(arg[iarg], "a") == 0) {
              if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for a");
              rr.coupling_a = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
              iarg += 2;
            } else if (strcmp(arg[iarg], "b") == 0) {
              if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for b");
              rr.coupling_b = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
              iarg += 2;
            } else if (strcmp(arg[iarg], "taper") == 0) {
              if (iarg + 1 >= narg)
                error->universe_all(FLERR, "Fix msevb: missing value for taper");
              rr.coupling_taper = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
              iarg += 2;
            } else {
              break;    // not a coupling param, back to per-reaction loop
            }
          }
        } else if (strcmp(arg[iarg], "taper") == 0) {
          if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for taper");
          rxndefs.back().coupling_taper = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
          iarg += 2;
        } else if (strcmp(arg[iarg], "shells") == 0) {
          if (iarg + 1 >= narg)
            error->universe_all(FLERR, "Fix msevb: missing value for per-reaction shells");
          rxndefs.back().shells = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
          if (rxndefs.back().shells < 1)
            error->universe_all(FLERR, "Fix msevb: per-reaction shells must be >= 1");
          iarg += 2;
        } else {
          break;    // not a per-reaction keyword, back to outer loop
        }
      }
    } else if (strcmp(arg[iarg], "taper") == 0) {
      if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for taper");
      coupling_taper = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "product_states") == 0) {
      enumerate_product_states = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg], "max_states") == 0) {
      if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for max_states");
      max_states = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (max_states < 0) error->universe_all(FLERR, "Fix msevb: max_states must be >= 0");
      iarg += 2;
    } else if (strcmp(arg[iarg], "shells") == 0) {
      if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for shells");
      max_shells = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (max_shells < 1) error->universe_all(FLERR, "Fix msevb: shells must be >= 1");
      iarg += 2;
    } else if (strcmp(arg[iarg], "fermi_dirac") == 0) {
      if (iarg + 1 >= narg)
        error->universe_all(FLERR, "Fix msevb: missing temperature for fermi_dirac");
      fd_temperature = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      fermi_dirac_enabled = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "file") == 0) {
      if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: file requires a filename");
      file_name = arg[iarg + 1];
      file_flag = true;
      iarg += 2;
      // Optional 'every N': also write a record every N steps (in addition to
      // the records written whenever a reaction occurs).
      if (iarg < narg && strcmp(arg[iarg], "every") == 0) {
        if (iarg + 1 >= narg)
          error->universe_all(FLERR, "Fix msevb: file 'every' requires a value");
        file_every = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
        if (file_every < 1) error->universe_all(FLERR, "Fix msevb: file 'every' must be >= 1");
        iarg += 2;
      }
      // Open the JSON file and write its header on the universe root only.
      if (ipartition == 0 && comm->me == 0) {
        fpout = fopen(file_name.c_str(), "w");
        if (!fpout)
          error->one(
              FLERR,
              fmt::format("Fix msevb: cannot open file '{}': {}", file_name, utils::getsyserror()));
        utils::print(fpout,
                     "{\n"
                     "  \"application\": \"LAMMPS\",\n"
                     "  \"format\": \"dump\",\n"
                     "  \"style\": \"evb_states\",\n"
                     "  \"title\": \"fix msevb\",\n"
                     "  \"revision\": 1,\n"
                     "  \"timesteps\": [\n");
        fflush(fpout);
      }
    } else if (strcmp(arg[iarg], "scf_topology") == 0) {
      if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for scf_topology");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        scf_topology = true;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        scf_topology = false;
      else
        error->universe_all(FLERR, "Fix msevb: scf_topology must be 'yes' or 'no'");
      iarg += 2;
    } else if (strcmp(arg[iarg], "scf_max_iter") == 0) {
      if (iarg + 1 >= narg) error->universe_all(FLERR, "Fix msevb: missing value for scf_max_iter");
      scf_max_iter = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (scf_max_iter < 1) error->universe_all(FLERR, "Fix msevb: scf_max_iter must be >= 1");
      iarg += 2;
    } else {
      error->universe_all(FLERR, fmt::format("Fix msevb: unknown keyword '{}'", arg[iarg]));
    }
  }

  if (!reaction_enabled) error->universe_all(FLERR, "Fix msevb: requires 'reaction' args");

  // Propagate global coupling to reactions that have no per-reaction override.
  // If no global coupling was given, every reaction must supply its own.
  for (auto &rd : rxndefs) {
    if (!rd.coupling_set) {
      if (!coupling_enabled)
        error->universe_all(FLERR,
                            fmt::format("Fix msevb: reaction '{}' has no per-reaction coupling "
                                        "and no global coupling was specified; either add a "
                                        "global 'coupling' keyword or supply 'coupling' inside "
                                        "every reaction block",
                                        rd.pre_mol_id));
      rd.coupling_type = coupling_type;
      rd.coupling_lambda = coupling_lambda;
      rd.coupling_zeta = coupling_zeta;
      rd.coupling_v12 = coupling_v12;
      rd.coupling_alpha = coupling_alpha;
      rd.coupling_gamma_v = coupling_gamma_v;
      rd.coupling_a = coupling_a;
      rd.coupling_b = coupling_b;
      rd.coupling_taper = coupling_taper;
    }
  }

  // Propagate global shells default to reactions, then recompute max_shells
  for (auto &rd : rxndefs) {
    if (rd.shells < 0) rd.shells = max_shells;    // inherit global (defaults to 1)
  }
  {
    int effective_max = 1;
    for (const auto &rd : rxndefs) effective_max = MAX(effective_max, rd.shells);
    max_shells = effective_max;
  }

  // Validate coupling parameters for each reaction
  for (size_t r = 0; r < rxndefs.size(); r++) {
    const auto &rd = rxndefs[r];
    if (rd.coupling_type == COUPLING_RAITERI2011 &&
        (rd.coupling_lambda == 0.0 || rd.coupling_zeta == 0.0))
      error->universe_all(
          FLERR,
          fmt::format("Fix msevb: reaction {} geometry/gaussian requires lambda and zeta", r));
    else if (rd.coupling_type == COUPLING_VUILLEUMIER1998 && rd.coupling_v12 == 0.0)
      error->universe_all(FLERR,
                          fmt::format("Fix msevb: reaction {} geometry/exp requires v12", r));
    else if (rd.coupling_type == COUPLING_GRIMME2015 &&
             (rd.coupling_a == 0.0 || rd.coupling_b == 0.0))
      error->universe_all(
          FLERR, fmt::format("Fix msevb: reaction {} energy/gaussian requires a and b", r));
  }

  if (fermi_dirac_enabled) {
    fd_RT = fd_temperature * force->boltz;
    if (fd_RT <= 0.0)
      error->universe_all(FLERR, "Fix msevb: fermi_dirac temperature must be positive");
  }

  // ---- Fix flags ------------------------------------------------------
  no_change_box = 1;
  time_depend = 1;
  scalar_flag = 1;
  extscalar = 1;
  energy_global_flag = 1;    // compute_scalar contributes to compute pe
  virial_global_flag = 1;    // virial[6] contributes to compute pressure
  vector_flag = 1;
  // The global vector holds one entry per EVB state, so its length tracks the
  // number of states detected each step (nstates).  It is initialised to the
  // reference-only state count here and updated every step in
  // detect_reactive_sites() once nstates is known.
  size_vector = nstates;
  extvector = 1;
  comm_forward = 2;    // type + charge per atom
  nmax = 6;

  // ---- Fixed-size per-partition arrays --------------------------------
  int np = npartitions;
  epot_nmax = np;
  memory->create(epot, np, "msevb:epot");
  memory->create(my_epot, np, "msevb:my_epot");
  memory->create(my_nlocal, np, "msevb:my_nlocal");
  memory->create(all_nlocal, np, "msevb:all_nlocal");
  for (int i = 0; i < np; i++) {
    epot[i] = my_epot[i] = 0.0;
    my_nlocal[i] = all_nlocal[i] = 0;
  }
  epot_ground = 0.0;

  // ---- Inter-partition communicator -----------------------------------
  // init() runs at the start of every run; free the previous communicator
  // before re-splitting so repeated runs do not leak MPI handles.
  if (samerank != MPI_COMM_NULL) MPI_Comm_free(&samerank);
  int color = comm->me;
  int key = ipartition;
  MPI_Comm_split(universe->uworld, color, key, &samerank);

  // Verify identical domain decomposition across partitions
  for (int i = 0; i < np; i++) {
    my_nlocal[i] = 0;
    all_nlocal[i] = 0;
  }
  my_nlocal[ipartition] = atom->nlocal;
  MPI_Allreduce(my_nlocal, all_nlocal, np, MPI_INT, MPI_SUM, samerank);
  {
    int fail = 0;
    for (int i = 1; i < np; i++)
      if (all_nlocal[i] != all_nlocal[0]) fail = 1;
    int allfail = 0;
    MPI_Allreduce(&fail, &allfail, 1, MPI_INT, MPI_MAX, universe->uworld);
    if (allfail)
      error->universe_all(FLERR,
                          "Fix msevb: domain decomposition must be "
                          "identical across all partitions");
  }

  // ---- Computes -------------------------------------------------------
  // init() runs at the start of every run; create the computes only once and
  // reuse them, otherwise add_compute() errors on a duplicate ID on rerun.
  if (!pe) {
    id_pe = std::string(id) + "_pe";
    // Exclude fix energies to avoid circular dependency: our compute_scalar()
    // correction would contaminate the diagonal PE used for the Hamiltonian.
    // The standard thermo "pe" compute includes fix energies → shows mixed PE.
    pe = modify->add_compute(id_pe + " all pe pair bond angle dihedral improper kspace");
  }
  pe->addstep(update->ntimestep);

  if (!temp_compute) {
    id_temp = std::string(id) + "_temp";
    temp_compute = modify->add_compute(id_temp + " all temp");
  }
  temp_compute->addstep(update->ntimestep);

  if (!press_compute) {
    id_press = std::string(id) + "_press";
    press_compute = modify->add_compute(id_press + " all pressure " + id_temp);
  }
  press_compute->addstep(update->ntimestep);

  for (int i = 0; i < 6; i++) pressure[i] = 0.0;

  // Create the reactive-atom tracking group in the constructor so it
  // exists as soon as the fix is defined — before init() runs and before
  // any dump commands that reference the group name are processed.
  // The group starts empty and is populated by update_reactive_group()
  // after the first detect_reactive_sites() call.
  if (reaction_enabled) {
    reactive_group_name = std::string(id) + "_atoms";
    int gid = group->find(reactive_group_name.c_str());
    if (gid < 0) {
      char *gargv[2];
      gargv[0] = const_cast<char *>(reactive_group_name.c_str());
      gargv[1] = const_cast<char *>("empty");
      group->assign(2, gargv);
      gid = group->find(reactive_group_name.c_str());
    }
    if (gid < 0)
      error->all(
          FLERR,
          fmt::format("Fix msevb: failed to create reactive-atom group '{}'", reactive_group_name));
    reactive_group_bit = group->bitmask[gid];
    if (comm->me == 0)
      utils::logmesg(lmp,
                     fmt::format("Fix msevb: reactive-atom group '{}' created "
                                 "(use in dump commands to write MSEVB atoms)\n",
                                 reactive_group_name));
  }
}

/* ---------------------------------------------------------------------- */

FixMSEVB::~FixMSEVB()
{
  if (samerank != MPI_COMM_NULL) MPI_Comm_free(&samerank);
  // Computes are created in init(); guard deletes in case the fix is
  // destroyed before init() ran (e.g. an input error after the fix command).
  if (pe) modify->delete_compute(id_pe);
  if (temp_compute) modify->delete_compute(id_temp);
  if (press_compute) modify->delete_compute(id_press);
  memory->destroy(commbuf);
  memory->destroy(epot);
  memory->destroy(my_epot);
  memory->destroy(hamiltonian);
  memory->destroy(H_work);
  memory->destroy(eigenvalues);
  memory->destroy(eigenvectors);
  memory->destroy(amplitudes);
  memory->destroy(fd_occ);
  memory->destroy(my_nlocal);
  memory->destroy(all_nlocal);
  memory->destroy(weights);
  memory->destroy(saved_forces);
  delete[] excess_forces;
  memory->destroy(sites);
  memory->destroy(chain_H_flat);
  memory->destroy(chain_X_flat);
  memory->destroy(chain_Y_flat);
  memory->destroy(chain_rxn_flat);
  memory->destroy(chain_glove_flat);
  memory->destroy(glove_flat);
  destroy_ref_topology();
  if (fpout) {
    // Close the "timesteps" array and the root object.  Works whether or not
    // any records were written (json_init tracks that an entry exists, but the
    // closing brackets are valid for an empty array too).
    utils::print(fpout, "\n  ]\n}\n");
    fclose(fpout);
    fpout = nullptr;
  }
  // Note: do NOT call group->assign("delete") here.  By the time the
  // destructor runs during LAMMPS teardown, lmp->group may already be
  // destroyed.  LAMMPS cleans up all groups itself; we just clear our
  // cached bit so it is obviously invalid.
  reactive_group_bit = 0;
}

/* ---------------------------------------------------------------------- */

int FixMSEVB::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= POST_RUN;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMSEVB::init()
{
  int onenmax = MAX(nmax, 3 * atom->nmax);
  MPI_Allreduce(&onenmax, &nmax, 1, MPI_INT, MPI_MAX, universe->uworld);
  memory->destroy(commbuf);
  memory->create(commbuf, nmax, "msevb:commbuf");

  if (modify->get_fix_by_style("^balance").size() > 0)
    error->universe_all(FLERR, "Fix msevb is not compatible with load balancing");
  if (modify->get_fix_by_style("^msevb").size() > 1)
    error->universe_all(FLERR, "There may be only one fix msevb at a time");
  if (utils::strmatch(update->integrate_style, "^respa"))
    error->universe_all(FLERR, "Must not use run style respa with fix msevb");

  MPI_Bcast(&domain->boxlo[0], 3, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(&domain->boxhi[0], 3, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(&domain->yz, 1, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(&domain->xz, 1, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(&domain->xy, 1, MPI_DOUBLE, 0, samerank);
  domain->set_global_box();
  domain->set_local_box();

  // ---- Template reaction initialisation --------------------------------
  // Load molecules (must have been declared with 'molecule' command before
  // this fix), parse map files, compute topology diffs, and determine
  // glove_nmax.
  glove_nmax = 0;
  for (auto &rd : rxndefs) {
    // Find pre-reaction molecule
    int pre_idx = atom->find_molecule(rd.pre_mol_id.c_str());
    if (pre_idx < 0)
      error->universe_all(FLERR,
                          fmt::format("Fix msevb: reaction pre-mol '{}' not found; "
                                      "use 'molecule' command before fix msevb",
                                      rd.pre_mol_id));
    rd.pre_mol = atom->molecules[pre_idx];

    // Find post-reaction molecule
    int post_idx = atom->find_molecule(rd.post_mol_id.c_str());
    if (post_idx < 0)
      error->universe_all(
          FLERR, fmt::format("Fix msevb: reaction post-mol '{}' not found", rd.post_mol_id));
    rd.post_mol = atom->molecules[post_idx];

    // Parse the map file (fills ibonding, jbonding, is_edge, pre_to_post)
    parse_template_mapfile(rd);

    // Auto-detect which of ibonding/jbonding is H (the transferring atom).
    // H is the atom that breaks a bond to its pre-template neighbor X (the
    // donor) and forms a new bond to Y (the acceptor).  The react.map
    // InitiatorIDs may list them in either order (H first or Y first), so
    // we try ibonding first; if no broken bond is found, try jbonding and
    // swap the two if jbonding is actually H.
    rd.ix_bonding = -1;
    rd.type_X = -1;
    {
      const Molecule *pre = rd.pre_mol;
      const Molecule *post = rd.post_mol;

      // Helper: check if 'candidate' loses a bond to some neighbor (not
      // 'other') between pre and post.  Returns the pre-index of the lost
      // neighbor, or -1.
      auto find_broken_bond = [&](int candidate, int other) -> int {
        int post_cand = rd.pre_to_post[candidate];
        if (post_cand < 0) return -1;
        for (int j = 0; j < pre->nspecial[candidate][0]; j++) {
          int ni = (int) (pre->special[candidate][j] - 1);
          if (ni == other) continue;
          int ni_post = (rd.pre_to_post[ni] >= 0) ? rd.pre_to_post[ni] : -1;
          if (ni_post < 0) continue;
          bool bonded_in_post = false;
          for (int k = 0; k < post->nspecial[post_cand][0]; k++) {
            if ((int) (post->special[post_cand][k] - 1) == ni_post) {
              bonded_in_post = true;
              break;
            }
          }
          if (!bonded_in_post) return ni;
        }
        return -1;
      };

      // Try ibonding as H (the common case: ibonding = H, jbonding = Y)
      int x_idx = find_broken_bond(rd.ibonding, rd.jbonding);
      if (x_idx >= 0) {
        rd.ix_bonding = x_idx;
        rd.type_X = pre->type[x_idx];
      } else {
        // ibonding has no broken bonds — try jbonding as H instead
        x_idx = find_broken_bond(rd.jbonding, rd.ibonding);
        if (x_idx >= 0) {
          // jbonding is actually H; swap so ibonding = H, jbonding = Y
          std::swap(rd.ibonding, rd.jbonding);
          rd.ix_bonding = x_idx;
          rd.type_X = pre->type[x_idx];
        }
      }
    }

    // Set type_H / type_Y from pre-template (after possible swap above)
    rd.type_H = rd.pre_mol->type[rd.ibonding];
    rd.type_Y = rd.pre_mol->type[rd.jbonding];

    // Compute bond/type topology diff
    compute_topology_diff(rd);

    // Warn if the post-reaction template would exceed the per-atom topology
    // capacities allocated at read_data time (extra/*/per/atom).
    check_reaction_topology_limits(rd);

    rd.glove_n = rd.pre_mol->natoms;
    glove_nmax = MAX(glove_nmax, rd.glove_n);
  }

  if (reaction_enabled && ipartition == 0) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

// Parse REACTER-format map file for a template reaction.
// Fills: rd.ibonding, rd.jbonding, rd.is_edge[], rd.pre_to_post[].
// All atom IDs in the file are 1-based; stored internally as 0-based.
void FixMSEVB::parse_template_mapfile(ReactionDef &rd)
{
  const int natoms_pre = rd.pre_mol->natoms;
  const int natoms_post = rd.post_mol->natoms;

  // fix msevb applies a topology *difference* to a fixed set of atoms; it
  // cannot insert or delete atoms.  The pre and post templates must therefore
  // describe the same atoms.  Reject a change in atom count up front.
  if (natoms_pre != natoms_post)
    error->universe_all(
        FLERR,
        fmt::format("Fix msevb: pre template '{}' has {} atoms but post template '{}' has {}; "
                    "fix msevb requires atom-conserving reactions and cannot insert or delete "
                    "atoms (the pre- and post-reaction templates must contain the same number "
                    "of atoms)",
                    rd.pre_mol_id, natoms_pre, rd.post_mol_id, natoms_post));

  rd.glove_n = natoms_pre;
  rd.is_edge.assign(natoms_pre, 0);
  rd.pre_to_post.assign(natoms_pre, -1);
  rd.ibonding = -1;
  rd.jbonding = -1;

  FILE *fp = fopen(rd.map_file_path.c_str(), "r");
  if (!fp)
    error->universe_all(FLERR,
                        fmt::format("Fix msevb: cannot open map file '{}'", rd.map_file_path));

  char line[1024];
  int section = 0;    // 0=none 1=InitiatorIDs 2=EdgeIDs 3=Equivalences

  auto read_line = [&]() -> bool {
    while (fgets(line, sizeof(line), fp)) {
      char *nl = strchr(line, '\n');
      if (nl) *nl = '\0';
      char *p = line;
      while (*p == ' ' || *p == '\t') p++;
      if (*p == '\0' || *p == '#') continue;
      return true;
    }
    return false;
  };

  while (read_line()) {
    char *p = line;
    while (*p == ' ' || *p == '\t') p++;

    if (strncmp(p, "InitiatorIDs", 12) == 0) {
      section = 1;
    } else if (strncmp(p, "EdgeIDs", 7) == 0) {
      section = 2;
    } else if (strncmp(p, "Equivalences", 12) == 0) {
      section = 3;
    } else {
      switch (section) {
        case 1: {
          int a, b;
          if (sscanf(p, "%d %d", &a, &b) == 2) {
            rd.ibonding = a - 1;
            rd.jbonding = b - 1;
          }
          section = 0;
          break;
        }
        case 2: {
          int id;
          char *tok = strtok(p, " \t");
          while (tok) {
            if (sscanf(tok, "%d", &id) == 1) {
              if (id >= 1 && id <= natoms_pre) rd.is_edge[id - 1] = 1;
            }
            tok = strtok(nullptr, " \t");
          }
          break;
        }
        case 3: {
          int pre_id, post_id;
          if (sscanf(p, "%d %d", &pre_id, &post_id) == 2) {
            if (pre_id < 1 || pre_id > natoms_pre || post_id < 1 || post_id > natoms_post)
              error->universe_all(
                  FLERR,
                  fmt::format("Fix msevb: Equivalences entry '{} {}' in map file '{}' is out of "
                              "range for templates with {} (pre) and {} (post) atoms",
                              pre_id, post_id, rd.map_file_path, natoms_pre, natoms_post));
            // fix msevb keeps each atom's identity across the reaction, so the
            // map must pair every atom with itself.  A pre/post mismatch means
            // the map is renumbering atoms, which is not supported.
            if (pre_id != post_id)
              error->universe_all(
                  FLERR,
                  fmt::format("Fix msevb: Equivalences entry maps atom {} to atom {} in map file "
                              "'{}'; fix msevb does not support changing atom IDs (the pre- and "
                              "post-reaction templates must use identical atom numbering, i.e. "
                              "every Equivalences pair must map an atom to itself)",
                              pre_id, post_id, rd.map_file_path));
            rd.pre_to_post[pre_id - 1] = post_id - 1;
          }
          break;
        }
      }
    }
  }
  fclose(fp);

  if (rd.ibonding < 0 || rd.jbonding < 0)
    error->universe_all(
        FLERR,
        fmt::format("Fix msevb: map file '{}' missing InitiatorIDs section", rd.map_file_path));

  // Every atom must appear in the Equivalences section.  A missing entry would
  // drop an atom from the reaction, i.e. effectively change the atom count.
  for (int i = 0; i < natoms_pre; i++)
    if (rd.pre_to_post[i] < 0)
      error->universe_all(
          FLERR,
          fmt::format("Fix msevb: map file '{}' Equivalences section does not map atom {}; every "
                      "atom in the pre-reaction template must be mapped one-to-one to the "
                      "post-reaction template",
                      rd.map_file_path, i + 1));
}

/* ---------------------------------------------------------------------- */

// Compute the topology diff between pre and post templates.
// Fills: rd.type_changes, rd.bond_breaks, rd.bond_creates, rd.bond_retypes.
void FixMSEVB::compute_topology_diff(ReactionDef &rd)
{
  const Molecule *pre = rd.pre_mol;
  const Molecule *post = rd.post_mol;
  const int natoms_pre = pre->natoms;
  const int natoms_post = post->natoms;

  rd.type_changes.clear();
  rd.bond_breaks.clear();
  rd.bond_creates.clear();
  rd.bond_retypes.clear();

  // Build inverse map: post atom 0-based idx → pre atom 0-based idx.
  std::vector<int> post_to_pre(natoms_post, -1);
  for (int i = 0; i < natoms_pre; i++)
    if (rd.pre_to_post[i] >= 0) post_to_pre[rd.pre_to_post[i]] = i;

  // 1. Type changes: pre atom i whose type changes in post.
  for (int i = 0; i < natoms_pre; i++) {
    int post_i = rd.pre_to_post[i];
    if (post_i < 0) continue;
    if (pre->type[i] != post->type[post_i]) rd.type_changes.push_back({i, post->type[post_i]});
  }

  // Helper: check if atoms (a, b) are 1-2 bonded in a molecule (0-based).
  auto bonded_in_mol = [](const Molecule *mol, int a, int b) -> bool {
    for (int j = 0; j < mol->nspecial[a][0]; j++)
      if ((int) (mol->special[a][j] - 1) == b) return true;
    return false;
  };

  // 2. Bond breaks: bonds present in pre but absent in post.
  for (int i = 0; i < natoms_pre; i++) {
    for (int j = 0; j < pre->nspecial[i][0]; j++) {
      int ni = (int) (pre->special[i][j] - 1);
      if (ni <= i) continue;
      int post_i = rd.pre_to_post[i];
      int post_ni = rd.pre_to_post[ni];
      if (post_i < 0 || post_ni < 0) continue;
      if (!bonded_in_mol(post, post_i, post_ni)) rd.bond_breaks.push_back({i, ni});
    }
  }

  // Helper: find the bond type between atoms a and b in a molecule (0-based).
  // Returns 0 if not found.
  auto bond_type_in_mol = [](const Molecule *mol, int a, int b) -> int {
    for (int k = 0; k < mol->num_bond[a]; k++)
      if ((int) (mol->bond_atom[a][k] - 1) == b) return mol->bond_type[a][k];
    // bond may be stored on b
    for (int k = 0; k < mol->num_bond[b]; k++)
      if ((int) (mol->bond_atom[b][k] - 1) == a) return mol->bond_type[b][k];
    return 0;
  };

  // 3. Bond creates: bonds present in post but absent in pre.
  for (int p = 0; p < natoms_post; p++) {
    for (int j = 0; j < post->nspecial[p][0]; j++) {
      int q = (int) (post->special[p][j] - 1);
      if (q <= p) continue;
      int pre_p = post_to_pre[p];
      int pre_q = post_to_pre[q];
      if (pre_p < 0 || pre_q < 0) continue;
      if (!bonded_in_mol(pre, pre_p, pre_q)) {
        int bt = bond_type_in_mol(post, p, q);
        rd.bond_creates.push_back({pre_p, pre_q, bt});
      }
    }
  }

  // 4. Bond retypes: bonds present in both pre and post but with different
  // type.
  for (int i = 0; i < natoms_pre; i++) {
    for (int j = 0; j < pre->nspecial[i][0]; j++) {
      int ni = (int) (pre->special[i][j] - 1);
      if (ni <= i) continue;
      int post_i = rd.pre_to_post[i];
      int post_ni = rd.pre_to_post[ni];
      if (post_i < 0 || post_ni < 0) continue;
      if (!bonded_in_mol(post, post_i, post_ni)) continue;    // broken, not a retype
      int pre_bt = bond_type_in_mol(pre, i, ni);
      int post_bt = bond_type_in_mol(post, post_i, post_ni);
      if (pre_bt != post_bt && post_bt > 0) rd.bond_retypes.push_back({i, ni, post_bt});
    }
  }
}

/* ----------------------------------------------------------------------
   Warn if a reaction's post-reaction template needs more bonds/angles/
   dihedrals/impropers/special neighbors per atom than the system was
   allocated to hold.

   The per-atom topology capacities (atom->bond_per_atom, angle_per_atom,
   dihedral_per_atom, improper_per_atom and atom->maxspecial) are fixed when
   read_data runs, optionally enlarged with the extra/<kind>/per/atom keywords.
   fix msevb cannot grow these arrays, so a reaction that raises an atom's
   connectivity beyond them will overflow per-atom storage when it fires.

   The check is based on the templates alone.  Because a matched atom may carry
   additional topology outside the template, the post-template maxima are a
   lower bound on the real post-reaction counts: a failing check guarantees an
   overflow, while a clean check does not guarantee its absence.  Hence a
   warning rather than a hard error.
------------------------------------------------------------------------- */

void FixMSEVB::check_reaction_topology_limits(const ReactionDef &rd)
{
  // per-atom capacities are global, so emit the warning from one rank only
  if (universe->me != 0) return;

  const Molecule *pre = rd.pre_mol;
  const Molecule *post = rd.post_mol;

  auto check = [&](const char *kind, int pre_n, int post_n, int cap, const char *flag) {
    if (post_n <= cap) return;
    error->warning(
        FLERR,
        fmt::format("Fix msevb: reaction '{}' -> '{}' needs up to {} {}/atom "
                    "(pre-reaction template: {}/atom), but the system is allocated for only "
                    "{} {}/atom; this reaction may overflow per-atom storage when it fires. "
                    "Increase the capacity with '{}' on the read_data command",
                    rd.pre_mol_id, rd.post_mol_id, post_n, kind, pre_n, cap, kind, flag));
  };

  check("bonds", pre->bond_per_atom, post->bond_per_atom, atom->bond_per_atom,
        "extra/bond/per/atom");
  check("angles", pre->angle_per_atom, post->angle_per_atom, atom->angle_per_atom,
        "extra/angle/per/atom");
  check("dihedrals", pre->dihedral_per_atom, post->dihedral_per_atom, atom->dihedral_per_atom,
        "extra/dihedral/per/atom");
  check("impropers", pre->improper_per_atom, post->improper_per_atom, atom->improper_per_atom,
        "extra/improper/per/atom");
  check("special neighbors", pre->maxspecial, post->maxspecial, atom->maxspecial,
        "extra/special/per/atom");
}

/* ---------------------------------------------------------------------- */

void FixMSEVB::init_list(int /*id*/, NeighList *ptr)
{ list = ptr; }

/* ---------------------------------------------------------------------- */

void FixMSEVB::setup(int vflag)
{
  if (universe->me == 0) {
    std::string msg;
    msg += fmt::format("\nFix msevb: {} partition(s)\n", npartitions);
    msg += "----------------------------------------------------------------\n";

    for (size_t r = 0; r < rxndefs.size(); r++) {
      const auto &rd = rxndefs[r];
      msg += fmt::format("\n  Reaction {}:  {} -> {}  (map: {})\n", r, rd.pre_mol_id,
                         rd.post_mol_id, rd.map_file_path);
      msg += fmt::format("    {:<17}{:.4f} A\n", "cutoff:", std::sqrt(rd.cutoff_sq));
      msg += fmt::format("    {:<17}H={}  X={}  Y={}\n", "atom types:", rd.type_H, rd.type_X,
                         rd.type_Y);
      msg += fmt::format("    {:<17}{}\n", "shells:", rd.shells);

      // Coupling
      const char *src = rd.coupling_set ? "per-reaction" : "global";
      std::string cpl;
      switch (rd.coupling_type) {
        case COUPLING_NONE:
          cpl = fmt::format("none  [{}]", src);
          break;
        case COUPLING_RAITERI2011:
          cpl = fmt::format("geometry/gaussian  lambda={:.4f}  zeta={:.4f}", rd.coupling_lambda,
                            rd.coupling_zeta);
          if (rd.coupling_taper > 0.0) cpl += fmt::format("  taper={:.4f}", rd.coupling_taper);
          cpl += fmt::format("  [{}]", src);
          break;
        case COUPLING_VUILLEUMIER1998:
          cpl = fmt::format("geometry/exp  v12={:.4f}  alpha={:.4f}"
                            "  gamma={:.4f}",
                            rd.coupling_v12, rd.coupling_alpha, rd.coupling_gamma_v);
          if (rd.coupling_taper > 0.0) cpl += fmt::format("  taper={:.4f}", rd.coupling_taper);
          cpl += fmt::format("  [{}]", src);
          break;
        case COUPLING_GRIMME2015:
          cpl = fmt::format("energy/gaussian  a={:.4f}  b={:.4f}", rd.coupling_a, rd.coupling_b);
          if (rd.coupling_taper > 0.0) cpl += fmt::format("  taper={:.4f}", rd.coupling_taper);
          cpl += fmt::format("  [{}]", src);
          break;
      }
      msg += fmt::format("    {:<17}{}\n", "coupling:", cpl);

      // Topology diff
      msg += fmt::format("    {:<17}{} type-change  {} bond-break"
                         "  {} bond-create  {} bond-retype\n",
                         "topology:", rd.type_changes.size(), rd.bond_breaks.size(),
                         rd.bond_creates.size(), rd.bond_retypes.size());
    }

    msg += "\n  Global options:\n";
    msg += fmt::format("    {:<17}{}\n", "max shells:", max_shells);
    if (fermi_dirac_enabled)
      msg +=
          fmt::format("    {:<17}T={:.2f} K  (RT={:.6f})\n", "fermi-dirac:", fd_temperature, fd_RT);
    if (enumerate_product_states) msg += fmt::format("    {:<17}enabled\n", "product states:");
    if (file_flag) {
      if (file_every > 0)
        msg += fmt::format("    {:<17}{}  (on reaction + every {} steps)\n", "file:", file_name,
                           file_every);
      else
        msg += fmt::format("    {:<17}{}  (on reaction)\n", "file:", file_name);
    }

    msg += "----------------------------------------------------------------\n\n";

    utils::logmesg(lmp, msg);
  }
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMSEVB::min_setup(int vflag)
{
  // Called by modify->setup() during minimization (whichflag==2) instead of
  // setup().  Forces have already been computed by the minimizer setup; run
  // the EVB post-force pass so that epot_ground is set, forces are EVB-mixed,
  // and the step-0 thermo output shows the correct ground-state energy.
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

int FixMSEVB::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    buf[m++] = ubuf(atom->type[j]).d;
    buf[m++] = atom->q[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixMSEVB::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0, last = first + n;
  for (int i = first; i < last; i++) {
    atom->type[i] = (int) ubuf(buf[m++]).i;
    atom->q[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void FixMSEVB::sync_positions()
{
  const int nall = atom->nlocal;
  if (nall > 0) MPI_Bcast(&atom->x[0][0], 3 * nall, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(&domain->boxlo[0], 3, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(&domain->boxhi[0], 3, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(&domain->yz, 1, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(&domain->xz, 1, MPI_DOUBLE, 0, samerank);
  MPI_Bcast(&domain->xy, 1, MPI_DOUBLE, 0, samerank);
  domain->set_global_box();
  domain->set_local_box();
}

/* ---------------------------------------------------------------------- */

void FixMSEVB::post_integrate()
{
  check_consistency_atoms();
  sync_positions();
  restore_reference_topology();
  detect_reactive_sites();
  update_reactive_group();
  apply_per_partition_state_changes();

  if (nsites > 0) next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixMSEVB::setup_pre_force(int vflag)
{
  if (reaction_enabled && !ref_snapshot_valid) snapshot_reference_topology();

  if (reaction_enabled) {
    restore_reference_topology();
    detect_reactive_sites();
    update_reactive_group();
    apply_per_partition_state_changes();

    if (nsites > 0) {
      if (domain->triclinic) domain->x2lamda(atom->nlocal);
      domain->pbc();
      comm->exchange();
      comm->borders();
      if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
      sync_before_neighbor_build();
      neighbor->build(1);
    }
  }
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMSEVB::pre_force(int /*vflag*/)
{
  if (!reaction_enabled || nsites == 0) return;
  int istate = (ipartition < nstates) ? ipartition : 0;
  if (istate == 0) return;
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

void FixMSEVB::min_pre_force(int vflag)
{
  // If we already received a stop signal, skip everything — forces remain
  // zeroed (force_clear() ran before this hook) and pair/bond will overwrite
  // them, but min_post_force will zero them again so CG sees zero gradient.
  if (min_terminate) return;

  // Symmetric stop handshake: each running partition sends 1, a finished
  // partition sends 0 from post_run(). MPI_MIN == 0 means at least one
  // partition is done — all should stop. This works regardless of which
  // partition finishes first.
  int local = 1, result = 1;
  MPI_Allreduce(&local, &result, 1, MPI_INT, MPI_MIN, samerank);
  if (result == 0) {
    min_terminate = true;
    return;
  }

  // During minimization there is no post_integrate, so site detection and
  // per-partition topology changes are handled here instead.
  check_consistency_atoms();
  sync_positions();
  restore_reference_topology();
  detect_reactive_sites();
  update_reactive_group();
  apply_per_partition_state_changes();

  if (nsites > 0) {
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    comm->borders();
    if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
    sync_before_neighbor_build();
    neighbor->build(1);
  }

  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMSEVB::min_post_force(int vflag)
{
  if (min_terminate) {
    // pair/bond set forces on stale topology — zero them so CG sees zero
    // gradient and exits cleanly via FTOL on the next iteration check.
    double **f = atom->f;
    const int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) f[i][0] = f[i][1] = f[i][2] = 0.0;
    return;
  }
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMSEVB::post_run()
{
  // During minimization (whichflag==2), partitions may finish at different
  // steps.  The first to finish sends 0 here to match the Allreduce in
  // min_pre_force on any still-running partition, causing it to terminate.
  // If min_terminate is already set, this partition was stopped by another
  // partition's signal — skip the Allreduce (no partner is waiting for us).
  if (update->whichflag == 2 && !min_terminate) {
    int local = 0, result = 0;
    MPI_Allreduce(&local, &result, 1, MPI_INT, MPI_MIN, samerank);
  }
  min_terminate = false;
}

/* ---------------------------------------------------------------------- */

void FixMSEVB::post_force(int vflag)
{
  // Initialize virial tallying (zeros virial[6] when thermo_virial is set)
  v_init(vflag);

  // Grow commbuf if needed
  if (3 * atom->nmax > nmax) {
    nmax = 3 * atom->nmax;
    memory->grow(commbuf, nmax, "msevb:commbuf");
  }

  grow_epot_arrays();
  grow_eigensystem_arrays();

  // Check if any reaction uses non-NONE coupling
  bool any_non_none = false;
  for (const auto &rd : rxndefs) {
    if (rd.coupling_type != COUPLING_NONE) any_non_none = true;
  }

  // Need force mixing when coupling is active or FD is enabled
  const bool need_forces = any_non_none || fermi_dirac_enabled;

  // 1. Gather potential energies (non-blocking)
  gather_potential_energies();

  // 2. Wait for PE, build Hamiltonian, compute coupling values
  MPI_Wait(&epot_request, MPI_STATUS_IGNORE);
  build_hamiltonian();
  if (any_non_none) compute_coupling_values();

  if (excess_force_mode == 0) {
    // ---- Approach A: save/restore + excess_forces buffer ----

    // 3. Excess states (saves/restores atom->f, fills excess_forces)
    if (nsites_serial > 0) compute_excess_states();

    // 4. Solve eigensystem
    solve_eigensystem();

    // 5. Force mixing
    if (need_forces) {
      compute_mixing_weights();
      weight_based_hellmann_feynman_forces();
    }

  } else {
    // ---- Approach B: split energy-only + deferred force application ----

    // 3. Excess state energies + coupling values (no force storage)
    if (nsites_serial > 0) compute_excess_energies();

    // 4. Solve eigensystem
    solve_eigensystem();

    // 5. Force mixing (parallel states only)
    if (need_forces) {
      compute_mixing_weights();
      weight_based_hellmann_feynman_forces();
    }

    // 6. Apply excess forces with known weights
    if (nsites_serial > 0 && need_forces) apply_excess_forces();
  }

  // --- Permanent transfer + optional SCF topology loop -----------------
  //
  // If scf_topology is enabled and a permanent transfer occurs, re-evaluate
  // the EVB system (topology detection → force recomputation → EVB solve)
  // until no further transfer occurs or scf_max_iter is reached.  Positions
  // are never updated; only the reference topology evolves.
  //
  // Transfer output is written BEFORE the SCF loop begins so that the
  // nsites/nstates recorded in the output file reflect the state at the
  // moment of transfer, not the (potentially zero) value after the SCF
  // re-detection sees the new reference topology.
  int dominant_state = 0;
  double dominant_amp = 0.0;
  bool any_transferred = false;
  bool transfer_output_written = false;

  {
    bool transferred = do_permanent_transfer(dominant_state, dominant_amp);
    any_transferred = transferred;

    if (scf_topology && transferred) {
      // Write output now, while nsites/nstates still reflect the topology
      // that was present when the transfer was detected.
      if (file_flag && ipartition == 0 && comm->me == 0) {
        write_msevb_json(update->ntimestep, dominant_state, dominant_amp);
        transfer_output_written = true;
      }

      int scf_iter = 0;
      while (transferred && scf_iter < scf_max_iter) {
        scf_iter++;

        // Re-detect reactive sites.  do_permanent_transfer() already restored
        // the new reference topology on all partitions, so no position sync or
        // restore_reference_topology() is needed here.
        detect_reactive_sites();
        update_reactive_group();
        apply_per_partition_state_changes();

        // Rebuild neighbor list — specials changed, neighbor list must be
        // rebuilt before force evaluation.  Mirrors setup_pre_force() and
        // the batch loop inside compute_excess_states().
        if (nsites > 0) {
          if (domain->triclinic) domain->x2lamda(atom->nlocal);
          domain->pbc();
          comm->exchange();
          comm->borders();
          if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
          sync_before_neighbor_build();
          neighbor->build(1);
        }

        // Broadcast ghost type/charge for non-reference partitions (pre_force).
        if (reaction_enabled && nsites > 0) {
          int istate = (ipartition < nstates) ? ipartition : 0;
          if (istate != 0) comm->forward_comm(this);
        }

        // Recompute forces for this partition's current topology.
        {
          double **f = atom->f;
          const int cur_nlocal = atom->nlocal;
          for (int i = 0; i < cur_nlocal; i++) f[i][0] = f[i][1] = f[i][2] = 0.0;
          modified_forces_on_host();
          sync_before_force_compute();
          if (force->pair && force->pair->compute_flag) force->pair->compute(1, vflag);
          if (atom->molecular != Atom::ATOMIC) {
            if (force->bond) force->bond->compute(1, vflag);
            if (force->angle) force->angle->compute(1, vflag);
            if (force->dihedral) force->dihedral->compute(1, vflag);
            if (force->improper) force->improper->compute(1, vflag);
          }
          if (force->kspace && force->kspace->compute_flag) force->kspace->compute(1, vflag);
          if (force->newton) comm->reverse_comm();
          modified_after_force_compute();
          sync_forces_to_host();
        }

        // Gather PE directly from force-style energy accumulators.
        // Cannot use pe->compute_scalar() here because it caches the result
        // per timestep; instead we sum force-object energies directly,
        // mirroring the approach used in compute_excess_states().
        {
          grow_epot_arrays();
          const double scalefac = 1.0 / comm->nprocs;
          double eng = 0.0;
          if (force->pair) eng += force->pair->eng_vdwl + force->pair->eng_coul;
          if (atom->molecular != Atom::ATOMIC) {
            if (force->bond) eng += force->bond->energy;
            if (force->angle) eng += force->angle->energy;
            if (force->dihedral) eng += force->dihedral->energy;
            if (force->improper) eng += force->improper->energy;
          }
          MPI_Allreduce(MPI_IN_PLACE, &eng, 1, MPI_DOUBLE, MPI_SUM, world);
          if (force->kspace) eng += force->kspace->energy;

          for (int i = 0; i < npartitions; i++) my_epot[i] = 0.0;
          my_epot[ipartition] = scalefac * eng;
          MPI_Allreduce(my_epot, epot, npartitions, MPI_DOUBLE, MPI_SUM, universe->uworld);
        }

        // Re-run EVB evaluation with new topology.
        grow_eigensystem_arrays();
        build_hamiltonian();
        if (any_non_none) compute_coupling_values();

        if (excess_force_mode == 0) {
          if (nsites_serial > 0) compute_excess_states();
          solve_eigensystem();
          if (need_forces) {
            compute_mixing_weights();
            weight_based_hellmann_feynman_forces();
          }
        } else {
          if (nsites_serial > 0) compute_excess_energies();
          solve_eigensystem();
          if (need_forces) {
            compute_mixing_weights();
            weight_based_hellmann_feynman_forces();
          }
          if (nsites_serial > 0 && need_forces) apply_excess_forces();
        }

        // Check for another transfer in the new topology.
        transferred = do_permanent_transfer(dominant_state, dominant_amp);
        any_transferred |= transferred;

        // Write output for each subsequent transfer while nsites/nstates
        // still reflect the state at which this transfer was detected.
        if (transferred && file_flag && ipartition == 0 && comm->me == 0)
          write_msevb_json(update->ntimestep, dominant_state, dominant_amp);
      }

      if (transferred && universe->me == 0)
        utils::logmesg(lmp,
                       fmt::format("WARNING: MSEVB SCF topology did not converge in {} "
                                   "iterations at step {}\n",
                                   scf_max_iter, update->ntimestep));
    }
  }

  // Compute EVB-mixed virial tensor and pressure (using final converged state).
  //
  // Strategy: gather per-partition virials from force objects, amplitude-weight
  // them across partitions via samerank MPI, add coupling_virial, then compute
  // the correction so that LAMMPS thermo sees the mixed pressure.
  //
  // Force object virials (pair, bond, angle, ...) are per-rank and need
  // allreduce within the partition.  Kspace virial is already globally reduced.
  {
    // 1. Gather per-rank virial from force objects (excluding kspace)
    double my_vir[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    if (force->pair)
      for (int k = 0; k < 6; k++) my_vir[k] += force->pair->virial[k];
    if (atom->molecular != Atom::ATOMIC) {
      if (force->bond)
        for (int k = 0; k < 6; k++) my_vir[k] += force->bond->virial[k];
      if (force->angle)
        for (int k = 0; k < 6; k++) my_vir[k] += force->angle->virial[k];
      if (force->dihedral)
        for (int k = 0; k < 6; k++) my_vir[k] += force->dihedral->virial[k];
      if (force->improper)
        for (int k = 0; k < 6; k++) my_vir[k] += force->improper->virial[k];
    }

    // Allreduce within partition to get partition virial (pair+bond+... part)
    double partition_vir[6];
    MPI_Allreduce(my_vir, partition_vir, 6, MPI_DOUBLE, MPI_SUM, world);

    // Add kspace virial (already globally reduced within partition)
    if (force->kspace)
      for (int k = 0; k < 6; k++) partition_vir[k] += force->kspace->virial[k];

    // 2. Amplitude-weight and allreduce across partitions via samerank
    double amp = (ipartition < nstates) ? amplitudes[ipartition] : 0.0;
    const double scalefac = 1.0 / comm->nprocs;
    double weighted_vir[6];
    for (int k = 0; k < 6; k++) weighted_vir[k] = amp * scalefac * partition_vir[k];
    double mixed_vir[6];
    MPI_Allreduce(weighted_vir, mixed_vir, 6, MPI_DOUBLE, MPI_SUM, universe->uworld);

    // Add coupling virial (already in energy units)
    for (int k = 0; k < 6; k++) mixed_vir[k] += coupling_virial[k];

    // 3. Compute mixed pressure (for extract / JSON output)
    double volume = domain->xprd * domain->yprd * domain->zprd;
    double vir_to_press = force->nktv2p / volume;
    press_compute->compute_vector();
    for (int k = 0; k < 6; k++)
      pressure[k] = press_compute->vector[k] + (mixed_vir[k] - partition_vir[k]) * vir_to_press;
    press_compute->addstep(update->ntimestep + 1);

    // 4. Store virial correction for LAMMPS thermo (fix_modify virial yes).
    // compute_pressure allreduces fix->virial across ranks in the partition,
    // so store per-rank share = total_correction / nprocs.
    if (vflag_global) {
      double inv_nprocs = 1.0 / comm->nprocs;
      for (int k = 0; k < 6; k++) virial[k] = (mixed_vir[k] - partition_vir[k]) * inv_nprocs;
    }
  }

  // Transfer output was already written before the SCF loop (if SCF ran);
  // skip it here to avoid double-writing.  For SCF-disabled runs, write now.
  if (file_flag && !transfer_output_written && any_transferred && ipartition == 0 &&
      comm->me == 0) {
    write_msevb_json(update->ntimestep, dominant_state, dominant_amp);
    return;
  }
  if (transfer_output_written) return;

  // Periodic record: only when 'every N' was given and no reaction wrote this step.
  if (file_flag && file_every > 0 && (update->ntimestep % file_every == 0) && ipartition == 0 &&
      comm->me == 0) {
    write_msevb_json(update->ntimestep, dominant_state, dominant_amp);
    return;
  }
}

/* ----------------------------------------------------------------------
   Update the reactive-atom tracking group to contain exactly the atoms
   involved in the current set of detected reactive sites.

   Included atoms: every non-zero entry in chain_glove_flat for every
   non-product-state site across all chain depths.  The glove is the
   full BFS match of the pre-reaction molecule template, so for a
   hydronium/water reaction this includes the transferring H, the donor
   oxygen and its remaining hydrogens, the acceptor oxygen and all of
   its hydrogens — the complete reactive neighbourhood, not just H/X/Y.

   Product-state sites (n_components > 0) are skipped because their
   component sites already carry the full glove individually.

   The group bitmask is cleared on all local atoms first, then set for
   every matched atom that is owned by this rank.  A forward_comm is not
   needed: dump commands operate on local atoms and LAMMPS propagates
   group membership via atom->mask during the normal comm cycle.
---------------------------------------------------------------------- */

void FixMSEVB::update_reactive_group()
{
  if (!reaction_enabled || reactive_group_bit == 0) return;

  const int nlocal = atom->nlocal;
  int *mask = atom->mask;
  const int bit = reactive_group_bit;
  const int notbit = ~bit;

  // Clear the group bit on every local atom.
  for (int i = 0; i < nlocal; i++) mask[i] &= notbit;

  if (nsites == 0 || glove_nmax == 0) return;

  // Set the group bit for every atom that appears in the matched glove
  // of any active reactive site.  The glove contains every real atom
  // matched by the BFS template (e.g. for H3O+/H2O: the transferring H,
  // the donor O, the donor's other H atoms, the acceptor O, the
  // acceptor's other H atoms — the full pre-template neighbourhood).
  // Edge atoms are stored as 0 in the glove and are skipped.
  //
  // For multi-shell sites, chain_glove_flat stores a separate glove for
  // each depth, so all atoms in every step of the proton-hopping chain
  // are included.
  //
  // Product-state sites (n_components > 0) are skipped here because
  // their component sites are already processed individually and carry
  // the full glove data.
  for (int k = 0; k < nsites; k++) {
    if (sites[k].n_components > 0) continue;

    const int clen = sites[k].chain_len;
    const int rxn = sites[k].rxn_idx;
    const int gn = rxndefs[rxn].glove_n;    // number of atoms in this rxn template

    for (int d = 0; d < clen; d++) {
      // chain_glove_flat[(k * max_shells + d) * glove_nmax + g] is the
      // real global tag of pre-template atom g at chain depth d.
      const tagint *gptr = chain_glove_flat + (k * max_shells + d) * glove_nmax;

      for (int g = 0; g < gn; g++) {
        tagint tag = gptr[g];
        if (tag == 0) continue;    // edge atom or unmatched slot
        int idx = atom->map(tag);
        if (idx >= 0 && idx < nlocal) mask[idx] |= bit;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMSEVB::write_msevb_json(bigint timestep, int max_state, double max_amp)
{
  if (!fpout) return;    // only the universe root holds an open file

  const int ns = nstates;

  json rec;
  rec["timestep"] = timestep;
  rec["nstates"] = ns;
  rec["nsites_parallel"] = nsites_parallel;
  rec["nsites_serial"] = nsites_serial;

  // Helper: build the transfer chain (one link per depth) of site `si`.
  auto build_chain = [&](int si) {
    json chain = json::array();
    for (int d = 0; d < sites[si].chain_len; d++) {
      json link;
      link["X"] = chain_X_flat[si * max_shells + d];
      link["H"] = chain_H_flat[si * max_shells + d];
      link["Y"] = chain_Y_flat[si * max_shells + d];
      link["reaction"] = chain_rxn_flat[si * max_shells + d];
      chain.push_back(link);
    }
    return chain;
  };

  // Reactive states detected this step (one entry per site).
  json states = json::array();
  for (int k = 0; k < nsites; k++) {
    const ReactiveSite &s = sites[k];
    const bool is_product = s.n_components > 0;
    json st;
    st["index"] = k + 1;
    st["tag_X"] = s.tag_X;    // 0 means no bond-breaking partner
    st["tag_H"] = s.tag_H;
    st["tag_Y"] = s.tag_Y;
    st["parent_state"] = s.parent_state;
    // Product (combination) states apply several disjoint transfer chains at
    // once; mark them so they are distinguishable from single-transfer states,
    // and expose their constituent component states (whose own entries carry
    // the full per-chain data).  For a single-transfer state, "chain" holds its
    // own transfer chain as before.
    if (is_product) {
      st["kind"] = "product";
      st["n_components"] = s.n_components;
      json comps = json::array();
      json chain = json::array();    // concatenation of all component chains
      for (int ci = 0; ci < s.n_components; ci++) {
        const int comp = s.components[ci];
        json cj;
        cj["state"] = comp + 1;
        cj["chain"] = build_chain(comp);
        comps.push_back(cj);
        for (const auto &link : cj["chain"]) chain.push_back(link);
      }
      st["components"] = comps;
      st["chain"] = chain;
    } else {
      st["kind"] = (k < nsites_parallel) ? "parallel" : "serial";
      st["chain"] = build_chain(k);
    }
    states.push_back(st);
  }
  rec["states"] = states;

  // Effective Hamiltonian (row-major ns x ns).
  json hmat = json::array();
  for (int i = 0; i < ns; i++) {
    json row = json::array();
    for (int j = 0; j < ns; j++) row.push_back(hamiltonian[i * ns + j]);
    hmat.push_back(row);
  }
  rec["hamiltonian"] = hmat;

  // Eigenvalues and ground-state amplitudes (density-matrix diagonal under
  // Fermi-Dirac mixing, |c_n|^2 otherwise).
  json evals = json::array();
  json amps = json::array();
  for (int i = 0; i < ns; i++) {
    evals.push_back(eigenvalues[i]);
    amps.push_back(amplitudes[i]);
  }
  rec["eigenvalues"] = evals;
  rec["amplitudes"] = amps;

  if (fermi_dirac_enabled) {
    json fd;
    fd["temperature"] = fd_temperature;
    json occ = json::array();
    for (int k = 0; k < ns; k++) occ.push_back(fd_occ[k]);
    fd["occupations"] = occ;
    fd["E_mixed"] = epot_ground;
    rec["fermi_dirac"] = fd;
  }

  rec["E_ground"] = epot_ground;
  rec["dominant_state"] = max_state;
  rec["dominant_amplitude"] = max_amp;
  // Null when no reference-topology update occurs this step.
  if (max_state > 0)
    rec["reference_update"] = max_state;
  else
    rec["reference_update"] = nullptr;

  // Stream the record into the open "timesteps" array, separated by commas.
  if (json_init) utils::print(fpout, ",\n");
  json_init = 1;
  const std::string text = rec.dump(2);
  fwrite(text.data(), 1, text.size(), fpout);
  fflush(fpout);
}

/* ---------------------------------------------------------------------- */

// Kokkos sync hooks — no-ops in base class, overridden in fix_msevb/kk.
void FixMSEVB::sync_before_force_compute() {}
void FixMSEVB::modified_after_force_compute() {}
void FixMSEVB::sync_forces_to_host() {}
void FixMSEVB::modified_forces_on_host() {}
void FixMSEVB::modified_topology_on_host() {}
void FixMSEVB::sync_before_neighbor_build() {}

/* ---------------------------------------------------------------------- */

double FixMSEVB::compute_scalar()
{
  if (min_terminate) return 0.0;
  // Return PE correction: E_mixed - E_partition0
  // When fix_modify energy yes is set, compute_pe adds this to the raw
  // partition PE, giving the correct EVB-mixed potential energy.
  // Without fix_modify energy yes, f_evb shows the correction directly.
  return epot_ground - epot[ipartition];
}

double FixMSEVB::compute_vector(int n)
{
  // The vector is sized to the number of EVB states detected this step
  // (size_vector == nstates), so every state's energy is exposed, including
  // serial states beyond the partition count.  epot[] is grown to nstates
  // each step in grow_epot_arrays().
  if (n < 0 || n >= nstates) return 0.0;
  return epot[n];
}

void *FixMSEVB::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "pe") == 0) return &epot_ground;
  dim = 1;
  if (strcmp(str, "pressure") == 0) return pressure;
  return nullptr;
}
