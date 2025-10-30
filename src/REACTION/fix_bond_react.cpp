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
Contributing Author: Jacob Gissinger (jgissing@stevens.edu)
------------------------------------------------------------------------- */

#include "fix_bond_react.h"

#include "atom.h"
#include "atom_vec.h"
#include "citeme.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "random_mars.h"
#include "reset_atoms_mol.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include "superpose3d.h"

#include <cctype>
#include <cmath>
#include <cstring>

#include <algorithm>
#include <random>
#include <utility>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

static const char cite_fix_bond_react[] =
    "fix bond/react: reacter.org doi:10.1016/j.polymer.2017.09.038, "
    "doi:10.1021/acs.macromol.0c02012, doi:10.1016/j.cpc.2024.109287\n\n"
    "@Article{Gissinger17,\n"
    " author = {J. R. Gissinger and B. D. Jensen and K. E. Wise},\n"
    " title = {Modeling Chemical Reactions in Classical Molecular Dynamics Simulations},\n"
    " journal = {Polymer},\n"
    " year =    2017,\n"
    " volume =  128,\n"
    " pages =   {211--217}\n"
    "}\n\n"
    "@Article{Gissinger20,\n"
    " author = {J. R. Gissinger, B. D. Jensen, K. E. Wise},\n"
    " title = {{REACTER}: A Heuristic Method for Reactive Molecular Dynamics},\n"
    " journal = {Macromolecules},\n"
    " year =    2020,\n"
    " volume =  53,\n"
    " number =  22,\n"
    " pages =   {9953--9961}\n"
    "}\n\n"
    "@Article{Gissinger24,\n"
    " author = {J. R. Gissinger, B. D. Jensen, K. E. Wise},\n"
    " title = {Molecular Modeling of Reactive Systems with REACTER},\n"
    " journal = {Computer Physics Communications},\n"
    " year =    2024,\n"
    " volume =  304,\n"
    " number =  109287\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */
// clang-format off

FixBondReact::FixBondReact(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_bond_react);

  fix1 = nullptr;
  fix2 = nullptr;
  fix3 = nullptr;
  reset_mol_ids = nullptr;

  if (narg < 8) utils::missing_cmd_args(FLERR,"fix bond/react", error);

  newton_bond = force->newton_bond;

  restart_global = 1;
  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  global_freq = 1;
  extvector = 0;
  cuff = 1;
  status = Status::PROCEED;

  // reaction functions used by 'custom' constraint
  nrxnfunction = 3;
  rxnfunclist.resize(nrxnfunction);
  peratomflag.resize(nrxnfunction);
  rxnfunclist[0] = "rxnsum";
  peratomflag[0] = 1;
  rxnfunclist[1] = "rxnave";
  peratomflag[1] = 1;
  rxnfunclist[2] = "rxnbond";
  peratomflag[2] = 0;
  nvvec = 0;
  ncustomvars = 0;
  vvec = nullptr;

  nxspecial = nullptr;
  xspecial = nullptr;

  // these group names are reserved for use exclusively by bond/react
  master_group = "bond_react_MASTER_group";

  // by using fixed group names, only one instance of fix bond/react is allowed.
  if (modify->get_fix_by_style("^bond/react").size() != 0)
    error->all(FLERR, Error::NOLASTLINE, "Only one instance of fix bond/react allowed at a time");

  // let's find number of reactions specified
  int nrxns = 0;
  for (int i = 3; i < narg; i++) {
    if (strcmp(arg[i],"react") == 0) {
      nrxns++;
      i = i + 6; // skip past mandatory arguments
      if (i > narg) utils::missing_cmd_args(FLERR,"fix bond/react react", error);
    }
  }

  if (nrxns == 0)
    error->all(FLERR, Error::NOLASTLINE, "Fix bond/react is missing mandatory 'react' keyword");

  size_vector = nrxns;

  int iarg = 3;
  stabilization_flag = 0;
  molid_mode = Reset_Mol_IDs::YES;
  shuffle_seed = 0;
  int hang_catch = 0;
  int num_common_keywords = 50; // generous arbitrary limit
  while (true) {
    if (++hang_catch > num_common_keywords) error->all(FLERR, iarg, "Incorrect fix bond/react command syntax");
    if (strcmp(arg[iarg],"stabilization") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"fix bond/react stabilization", error);
      stabilization_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      if (stabilization_flag) {
        if (iarg+4 > narg)
          utils::missing_cmd_args(FLERR, "fix bond/react stabilization yes", error);
        exclude_group = arg[iarg+2];
        nve_limit_xmax = arg[iarg+3];
        iarg += 4;
      } else iarg += 2;
    } else if (strcmp(arg[iarg],"reset_mol_ids") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"fix bond/react reset_mol_ids", error);
      std::string str = arg[iarg+1];
      if (str == "yes") molid_mode = Reset_Mol_IDs::YES;
      else if (str == "no") molid_mode = Reset_Mol_IDs::NO;
      else if (str == "molmap") molid_mode = Reset_Mol_IDs::MOLMAP;
      else error->all(FLERR, iarg+1, "Unknown option {} for 'reset_mol_ids' keyword", str);
      iarg += 2;
    } else if (strcmp(arg[iarg],"max_rxn") == 0) {
      if (iarg+1 > narg) utils::missing_cmd_args(FLERR,"fix bond/react max_rxn", error);
      struct MaxRxnLimit maxlimit;
      maxlimit.Nrxns = 0;
      int j = iarg+1;
      while (isalpha(arg[j++][0])) {
        maxlimit.Nrxns++;
        if (j > narg+2) utils::missing_cmd_args(FLERR,"fix bond/react rate_limit", error);
      }
      if (maxlimit.Nrxns == 0) error->all(FLERR, iarg, "Illegal fix bond/react command: "
                                         "at least one rxn-ID should be listed directly after the 'max_rxn' keyword");
      if (iarg+maxlimit.Nrxns+2 > narg) utils::missing_cmd_args(FLERR,"fix bond/react max_rxn", error);
      for (int i = 0; i < maxlimit.Nrxns; i++) {
        std::string tmpstr = arg[iarg+1+i];
        maxlimit.rxn_names.push_back(tmpstr);
      }
      maxlimit.max_rxn = utils::inumeric(FLERR,arg[iarg+maxlimit.Nrxns+1],false,lmp);
      if (maxlimit.max_rxn < 0) error->all(FLERR, iarg, "Illegal fix bond/react command: "
                                         "'max_rxn' cannot be negative");
      max_rxn_limits.push_back(maxlimit);
      iarg += maxlimit.Nrxns+2;
    } else if (strcmp(arg[iarg],"rate_limit") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"fix bond/react rate_limit", error);
      struct RateLimit rlm;
      rlm.Nrxns = 0;
      int j = iarg+1;
      while (isalpha(arg[j++][0])) {
        rlm.Nrxns++;
        if (j > narg+2) utils::missing_cmd_args(FLERR,"fix bond/react rate_limit", error);
      }
      if (rlm.Nrxns == 0) error->all(FLERR, iarg, "Illegal fix bond/react command: "
                                         "at least one rxn-ID should be listed directly after the 'rate_limit' keyword");
      if (iarg+rlm.Nrxns+4 > narg) utils::missing_cmd_args(FLERR,"fix bond/react rate_limit", error);
      for (int i = 0; i < rlm.Nrxns; i++) {
        std::string tmpstr = arg[iarg+1+i];
        rlm.rxn_names.push_back(tmpstr);
      }
      char *myarg = arg[iarg+rlm.Nrxns+1]; // Nlimit
      if (strncmp(myarg,"v_",2) == 0) {
        rlm.var_flag = 1;
        rlm.var_id = input->variable->find(myarg);
        if (rlm.var_id < 0)
          error->all(FLERR,"Fix bond/react: Variable name {} for rate_limit does not exist",myarg);
        if (!input->variable->equalstyle(rlm.var_id))
          error->all(FLERR,"Fix bond/react: Variable {} for rate_limit is not equal-style",myarg);
      } else rlm.Nlimit = utils::inumeric(FLERR,myarg,false,lmp);
      rlm.Nsteps = utils::inumeric(FLERR,arg[iarg+rlm.Nrxns+2],false,lmp);
      rate_limits.push_back(rlm);
      iarg += rlm.Nrxns+3;
    } else if (strcmp(arg[iarg],"shuffle_seed") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"fix bond/react seed", error);
      shuffle_seed = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"react") == 0) {
      break;
    } else error->all(FLERR, iarg, "Unknown fix bond/react command keyword {}", arg[iarg]);
  }

  if (molid_mode == Reset_Mol_IDs::YES) {
    delete reset_mol_ids;
    reset_mol_ids = new ResetAtomsMol(lmp);
    reset_mol_ids->create_computes(id,group->names[igroup]);
  }

  // set up common variables as vectors of length 'nrxns'
  rxns.resize(nrxns);

  rescale_charges_anyflag = 0;
  int id = 0;
  for (auto &rxn : rxns) {
    rxn.ID = id++;
    rxn.fraction = 1.0;
    rxn.seed = 12345;
    rxn.stabilize_steps_flag = 0;
    rxn.custom_charges_fragid = -1;
    rxn.rescale_charges_flag = 0;
    rxn.create_atoms_flag = 0;
    rxn.modify_create_fragid = -1;
    rxn.overlapsq = 0.0;
    rxn.mol_total_charge = 0.0;
    rxn.molecule_keyword = Molecule_Keys::OFF;
    rxn.limit_duration = 60;
    rxn.reaction_count = 0;
    rxn.local_rxn_count = 0;
    rxn.ghostly_rxn_count = 0;
    rxn.reaction_count_total = 0;
    rxn.v_rmin = rxn.v_rmax = -1;
    rxn.v_nevery = rxn.v_prob = -1;
  }

  for (auto &rxn : rxns) {

    if (strcmp(arg[iarg],"react") != 0) error->all(FLERR,"Illegal fix bond/react command: "
                                                   "'react' or 'stabilization' has incorrect arguments");
    iarg++;

    rxn.name = arg[iarg++];
    if (rxn.name.size()+1 > MAXNAME) error->all(FLERR,"Reaction name (react-ID) is too long (limit: 255 characters)");

    int groupid = group->find(arg[iarg++]);
    if (groupid == -1) error->all(FLERR,"Could not find fix group ID");
    rxn.groupbits = group->bitmask[groupid];

    if (strncmp(arg[iarg],"v_",2) == 0) {
      rxn.v_nevery = input->variable->find(&arg[iarg][2]);
      validate_variable_keyword(&arg[iarg][2], rxn.v_nevery);
    } else {
      rxn.nevery = utils::inumeric(FLERR,arg[iarg],false,lmp);
      if (rxn.nevery <= 0) error->all(FLERR,"Illegal fix bond/react command: "
                                       "'Nevery' must be a positive integer");
    }
    iarg++;

    double cutoff;
    if (strncmp(arg[iarg],"v_",2) == 0) {
      rxn.v_rmin = input->variable->find(&arg[iarg][2]);
      validate_variable_keyword(&arg[iarg][2], rxn.v_rmin);
      cutoff = input->variable->compute_equal(rxn.v_rmin);
    } else cutoff = utils::numeric(FLERR,arg[iarg],false,lmp);
      if (cutoff < 0.0) error->all(FLERR,"Illegal fix bond/react command: "
                                   "'Rmin' cannot be negative");
      rxn.rminsq = cutoff*cutoff;
    iarg++;

    if (strncmp(arg[iarg],"v_",2) == 0) {
      rxn.v_rmax = input->variable->find(&arg[iarg][2]);
      validate_variable_keyword(&arg[iarg][2], rxn.v_rmax);
      cutoff = input->variable->compute_equal(rxn.v_rmax);
    } else cutoff = utils::numeric(FLERR,arg[iarg],false,lmp);
      if (cutoff < 0.0) error->all(FLERR,"Illegal fix bond/react command:"
                                   "'Rmax' cannot be negative");
      rxn.rmaxsq = cutoff*cutoff;
    iarg++;

    int mol_idx = atom->find_molecule(arg[iarg++]);
    if (mol_idx == -1) error->all(FLERR,"Pre-reaction molecule template ID for "
                                             "fix bond/react does not exist");
    rxn.reactant = atom->molecules[mol_idx];
    mol_idx = atom->find_molecule(arg[iarg++]);
    if (mol_idx == -1) error->all(FLERR,"Post-reaction molecule template ID for "
                                           "fix bond/react does not exist");
    rxn.product = atom->molecules[mol_idx];

    //read map file
    rxn.mapfilename = arg[iarg];
    iarg++;

    while (iarg < narg && strcmp(arg[iarg],"react") != 0) {
      if (strcmp(arg[iarg],"prob") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'prob' keyword has too few arguments");
        // check if probability is a variable
        if (strncmp(arg[iarg+1],"v_",2) == 0) {
          rxn.v_prob = input->variable->find(&arg[iarg][2]);
          validate_variable_keyword(&arg[iarg+1][2], rxn.v_prob);
          rxn.fraction = input->variable->compute_equal(rxn.v_prob);
        } else {
          // otherwise probability should be a number
          rxn.fraction = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        }
        rxn.seed = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
        if (rxn.fraction < 0.0 || rxn.fraction > 1.0)
          error->all(FLERR,"Illegal fix bond/react command: "
                     "probability fraction must between 0 and 1, inclusive");
        if (rxn.seed <= 0) error->all(FLERR,"Illegal fix bond/react command: "
                                       "probability seed must be positive");
        iarg += 3;
      } else if (strcmp(arg[iarg],"stabilize_steps") == 0) {
        if (stabilization_flag == 0) error->all(FLERR,"Stabilize_steps keyword "
                                                "used without stabilization keyword");
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'stabilize_steps' has too few arguments");
        rxn.limit_duration = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
        rxn.stabilize_steps_flag = 1;
        iarg += 2;
      } else if (strcmp(arg[iarg],"custom_charges") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'custom_charges' has too few arguments");
        if (strcmp(arg[iarg+1],"no") == 0) rxn.custom_charges_fragid = -1; //default
        else {
          rxn.custom_charges_fragid = rxn.reactant->findfragment(arg[iarg+1]);
          if (rxn.custom_charges_fragid < 0) error->one(FLERR,"Fix bond/react: Molecule fragment for "
                                                         "'custom_charges' keyword does not exist");
        }
        iarg += 2;
      } else if (strcmp(arg[iarg],"rescale_charges") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'rescale_charges' has too few arguments");
        if (strcmp(arg[iarg+1],"no") == 0) rxn.rescale_charges_flag = 0; //default
        else if (strcmp(arg[iarg+1],"yes") == 0) {
          if (!atom->q_flag) error->all(FLERR,"Illegal fix bond/react command: cannot use "
                                      "'rescale_charges' without atomic charges enabled");
          if (!rxn.product->qflag) error->all(FLERR,"Illegal fix bond/react command: cannot use "
                                      "'rescale_charges' without Charges section in post-reaction template");
          rxn.rescale_charges_flag = 1; // overloaded below to also indicate number of atoms to update
          rescale_charges_anyflag = 1;
          cuff = 2; // index shift for extra values carried around by mega_gloves
        } else error->one(FLERR,"Bond/react: Illegal option for 'rescale_charges' keyword");
        iarg += 2;
      } else if (strcmp(arg[iarg],"molecule") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'molecule' has too few arguments");
        if (strcmp(arg[iarg+1],"off") == 0) rxn.molecule_keyword = Molecule_Keys::OFF; //default
        else if (strcmp(arg[iarg+1],"inter") == 0) rxn.molecule_keyword = Molecule_Keys::INTER;
        else if (strcmp(arg[iarg+1],"intra") == 0) rxn.molecule_keyword = Molecule_Keys::INTRA;
        else error->one(FLERR,"Fix bond/react: Illegal option for 'molecule' keyword");
        iarg += 2;
      } else if (strcmp(arg[iarg],"modify_create") == 0) {
        if (iarg++ > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'modify_create' has too few arguments");
        while (iarg < narg && strcmp(arg[iarg],"react") != 0) {
          if (strcmp(arg[iarg],"fit") == 0) {
            if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                          "'modify_create' has too few arguments");
            if (strcmp(arg[iarg+1],"all") == 0) rxn.modify_create_fragid = -1; //default
            else {
              rxn.modify_create_fragid = rxn.product->findfragment(arg[iarg+1]);
              if (rxn.modify_create_fragid < 0) error->one(FLERR,"Fix bond/react: Molecule fragment for "
                                                             "'modify_create' keyword does not exist");
            }
            iarg += 2;
          } else if (strcmp(arg[iarg],"overlap") == 0) {
            if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                          "'modify_create' has too few arguments");
            rxn.overlapsq = utils::numeric(FLERR,arg[iarg+1],false,lmp);
            rxn.overlapsq *= rxn.overlapsq;
            iarg += 2;
          } else break;
        }
      } else if (strcmp(arg[iarg],"rate_limit") == 0) {
        error->all(FLERR,"Fix bond/react: 'rate_limit' as an 'individual keyword' has been deprecated. "
                         "Please use the 'rate_limit' common keyword instead, which can be applied to one or more reactions.");
      } else if (strcmp(arg[iarg],"max_rxn") == 0) {
        error->all(FLERR,"Fix bond/react: 'max_rxn' as an 'individual keyword' has been deprecated. "
                         "Please use the 'max_rxn' common keyword instead, which can be applied to one or more reactions.");
      } else error->all(FLERR,"Illegal fix bond/react command: unknown keyword");
    }
  }

  for (auto &rlm : rate_limits) {
    for (int i = 0; i < rlm.Nrxns; i++) {
      int existflag = 0;
      for (auto &rxn : rxns) {
        if (rlm.rxn_names[i] == rxn.name) {
          rlm.rxnIDs.push_back(rxn.ID);
          existflag = 1;
          break;
        }
      }
      if (existflag == 0) error->all(FLERR, "Fix bond/react: Invalid reaction name {} listed for rate_limit", rlm.rxn_names[i]);
    }
    rlm.store_rxn_counts.assign(rlm.Nsteps,-1);
  }
  for (auto &maxlimit : max_rxn_limits) {
    for (int i = 0; i < maxlimit.Nrxns; i++) {
      int existflag = 0;
      for (auto &rxn : rxns) {
        if (maxlimit.rxn_names[i] == rxn.name) {
          maxlimit.rxnIDs.push_back(rxn.ID);
          existflag = 1;
          break;
        }
      }
      if (existflag == 0) error->all(FLERR, "Fix bond/react: Invalid reaction name {} listed for max_rxn", maxlimit.rxn_names[i]);
    }
  }

  max_natoms = 0; // the number of atoms in largest molecule template
  for (auto &rxn : rxns) max_natoms = MAX(max_natoms,rxn.product->natoms);

  for (auto &rxn : rxns) {
    rxn.nnewmolids = 0;
    rxn.atoms.resize(max_natoms);
    int idx = 1;
    for (auto &atm : rxn.atoms) {
      atm.edge = 0;
      atm.recharged = 1; // update all partial charges by default
      atm.deleted = 0;
      atm.created = 0;
      atm.newmolid = 0;
      atm.chiral.fill(0);
      // default amap to their own molecule template atom ID
      // all but created atoms will be updated
      atm.amap.fill(idx++);
    }
  }

  if (molid_mode == Reset_Mol_IDs::MOLMAP) {
    for (auto &rxn : rxns) {
      if (!rxn.reactant->moleculeflag || !rxn.product->moleculeflag) {
        if (comm->me == 0)
          error->warning(FLERR,"Fix bond/react ('reset_mol_ids molmap' option): Pre- and post-reaction templates must "
                               "both contain a 'Molecules' section for molecule IDs to be updated for a given reaction");
        break;
      }
    }
    // 'new' mol IDs are ones that exist in post-reaction but not in pre-reaction
    // let's condense these and shift to be indexed from 1
    for (auto &rxn : rxns) {
      if (rxn.reactant->moleculeflag && rxn.product->moleculeflag) {
        for (int j = 0; j < rxn.product->natoms; j++) {
          if (rxn.atoms[j].newmolid != 0) continue;
          int molid_isnew = 1;
          for (int k = 0; k < rxn.reactant->natoms; k++) {
            if (rxn.product->molecule[j] == rxn.reactant->molecule[k]) {
              molid_isnew = 0;
              break;
            }
          }
          if (molid_isnew == 1) {
            rxn.nnewmolids++;
            for (int k = j; k < rxn.product->natoms; k++) {
              if (rxn.product->molecule[k] == rxn.product->molecule[j])
                rxn.atoms[k].newmolid = rxn.nnewmolids;
            }
          }
        }
      }
    }
  }

  // read all map files afterward
  for (auto &rxn : rxns) {
    fp = fopen(rxn.mapfilename.c_str(),"r");
    if (fp == nullptr) error->one(FLERR, "Fix bond/react: Cannot open map file {}", rxn.mapfilename);
    rxn.reactant->check_attributes();
    rxn.product->check_attributes();
    read_map_file(rxn);
    fclose(fp);
    rxn.iatomtype = rxn.reactant->type[rxn.ibonding-1];
    rxn.jatomtype = rxn.reactant->type[rxn.jbonding-1];
    find_landlocked_atoms(rxn);
    if (rxn.custom_charges_fragid >= 0) CustomCharges(rxn.custom_charges_fragid,rxn);
  }

  // charge rescaling values must be calculated after calling CustomCharges
  for (auto &rxn : rxns) {
    if (rxn.rescale_charges_flag) {
      rxn.rescale_charges_flag = 0; // will now store number of updated atoms
      for (int j = 0; j < rxn.product->natoms; j++) {
        int jj = rxn.atoms[j].amap[1]-1;
        if (rxn.atoms[jj].recharged == 1 && rxn.atoms[jj].deleted == 0) {
          rxn.mol_total_charge += rxn.product->q[j];
          rxn.rescale_charges_flag++;
        }
      }
    }
  }

  // get the names of per-atom variables needed by 'rxn' functions of custom constraint
  customvarnames();

  // initialize Marsaglia RNG with processor-unique seed (Arrhenius prob)
  for (auto &rxn : rxns)
    for (auto &constraint : rxn.constraints)
      if (constraint.type == Reaction::Constraint::Type::ARRHENIUS)
        constraint.arrhenius.rrhandom = new RanMars(lmp, (int) constraint.arrhenius.seed + comm->me);

  if (atom->molecular != Atom::MOLECULAR)
    error->all(FLERR,"Fix bond/react: Cannot use fix bond/react with non-molecular systems");

  // check if bonding atoms are 1-2, 1-3, or 1-4 bonded neighbors
  // if so, we don't need non-bonded neighbor list
  for (auto &rxn : rxns) {
    rxn.closeneigh = -1; // indicates will search non-bonded neighbors
    for (int k = 0; k < rxn.reactant->nspecial[rxn.ibonding-1][2]; k++) {
      if (rxn.reactant->special[rxn.ibonding-1][k] == rxn.jbonding) {
        rxn.closeneigh = 2; // index for 1-4 neighbor
        if (k < rxn.reactant->nspecial[rxn.ibonding-1][1])
          rxn.closeneigh = 1; // index for 1-3 neighbor
        if (k < rxn.reactant->nspecial[rxn.ibonding-1][0])
          rxn.closeneigh = 0; // index for 1-2 neighbor
        break;
      }
    }
  }

  // initialize Marsaglia RNG with processor-unique seed ('prob' keyword)

  random = new RanMars*[rxns.size()];
  for (auto &rxn : rxns) {
    random[rxn.ID] = new RanMars(lmp, rxn.seed + comm->me);
  }

  // set comm sizes needed by this fix
  // forward is big due to comm of broken bonds and 1-2 neighbors

  comm_forward = MAX(2,2+atom->maxspecial);
  comm_reverse = 2;

  // allocate arrays local to this fix
  nmax = 0;
  partner = finalpartner = nullptr;
  distsq = nullptr;
  allnattempt = 0;
  my_num_mega = 0;
  local_num_mega = 0;
  ghostly_num_mega = 0;
  global_megasize = 0;
  global_mega_glove = nullptr;
  custom_exclude_flag = 0;

  // used to store restart info
  set = new Set[rxns.size()];
  memset(set,0,rxns.size()*sizeof(Set));
}

/* ---------------------------------------------------------------------- */

FixBondReact::~FixBondReact()
{
  for (auto &rxn : rxns)
    for (auto &constraint : rxn.constraints)
      if (constraint.type == Reaction::Constraint::Type::ARRHENIUS)
        delete constraint.arrhenius.rrhandom;

  for (int i = 0; i < rxns.size(); i++) delete random[i];
  delete[] random;

  delete reset_mol_ids;

  memory->destroy(partner);
  memory->destroy(finalpartner);
  memory->destroy(distsq);
  if (vvec != nullptr) memory->destroy(vvec);
  memory->destroy(global_mega_glove);

  if (stabilization_flag == 1) {
    // delete fixes if not already deleted
    if (!id_fix1.empty() && modify->get_fix_by_id(id_fix1)) modify->delete_fix(id_fix1);
    if (!id_fix3.empty() && modify->get_fix_by_id(id_fix3)) modify->delete_fix(id_fix3);
  }
  if (!id_fix2.empty() && modify->get_fix_by_id(id_fix2)) modify->delete_fix(id_fix2);

  delete[] set;

  if (group) {
    group->assign(master_group + " delete");
    if (stabilization_flag == 1) group->assign(exclude_group + " delete");
  }
}

/* ---------------------------------------------------------------------- */

int FixBondReact::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  mask |= POST_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
let's add an internal nve/limit fix for relaxation of reaction sites
also let's add our per-atom property fix here!
this per-atom property will state the timestep an atom was 'limited'
it will have the name 'i_limit_tags' and will be intitialized to 0 (not in group)
------------------------------------------------------------------------- */

void FixBondReact::post_constructor()
{
  // let's add the limit_tags per-atom property fix
  id_fix2 = "bond_react_props_internal";
  if (!modify->get_fix_by_id(id_fix2))
    fix2 = modify->add_fix(id_fix2 + " all property/atom i_limit_tags i_react_tags ghost yes");

  // create master_group if not already existing
  // NOTE: limit_tags and react_tags automaticaly intitialized to zero (unless read from restart)
  group->find_or_create(master_group.c_str());
  std::string cmd = fmt::format("{} dynamic all property limit_tags",master_group);
  group->assign(cmd);

  if (stabilization_flag == 1) {
    int groupid = group->find(exclude_group);
    // create exclude_group if not already existing, or use as parent group if static
    if (groupid == -1 || group->dynamic[groupid] == 0) {

      // create stabilization per-atom property
      id_fix3 = "bond_react_stabilization_internal";
      if (!modify->get_fix_by_id(id_fix3))
        fix3 = modify->add_fix(id_fix3 + " all property/atom i_statted_tags ghost yes");

      statted_id = "statted_tags";

      // if static group exists, use as parent group
      // also, rename dynamic exclude_group by appending '_REACT'
      std::string exclude_PARENT_group = exclude_group;
      exclude_group = exclude_PARENT_group + "_REACT";

      group->find_or_create(exclude_group.c_str());
      if (groupid == -1)
        cmd = fmt::format("{} dynamic all property statted_tags", exclude_group);
      else
        cmd = fmt::format("{} dynamic {} property statted_tags", exclude_group, exclude_PARENT_group);
      group->assign(cmd);

      // on to statted_tags (system-wide thermostat)
      // initialize per-atom statted_flags to 1
      // (only if not already initialized by restart)
      if (fix3 && fix3->restart_reset != 1) {
        int flag,cols;
        int index = atom->find_custom("statted_tags",flag,cols);
        int *i_statted_tags = atom->ivector[index];

        for (int i = 0; i < atom->nlocal; i++)
          i_statted_tags[i] = 1;
      }
    } else {
      // sleeping code, for future capabilities
      custom_exclude_flag = 1;
      // first we have to find correct fix group reference
      Fix *fix = modify->get_fix_by_id("GROUP_" + exclude_group);

      // this returns names of corresponding property
      int unused;
      char *idprop;
      idprop = (char *) fix->extract("property",unused);
      if (idprop == nullptr)
        error->all(FLERR,"Exclude group must be a per-atom property group");
      statted_id = idprop;

      // initialize per-atom statted_tags to 1
      // need to correct for smooth restarts
      //int flag,cols;
      //int index = atom->find_custom(statted_id,flag,cols);
      //int *i_statted_tags = atom->ivector[index];
      //for (int i = 0; i < atom->nlocal; i++)
      //  i_statted_tags[i] = 1;
    }

    // let's create a new nve/limit fix to limit newly reacted atoms
    id_fix1 = "bond_react_MASTER_nve_limit";
    if (!modify->get_fix_by_id(id_fix1))
      fix1 = modify->add_fix(fmt::format("{} {} nve/limit {}", id_fix1, master_group, nve_limit_xmax));
  }
}

/* ---------------------------------------------------------------------- */

void FixBondReact::init()
{

  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;

  // check cutoff for iatomtype,jatomtype
  if (!utils::strmatch(force->pair_style,"^hybrid"))
    for (auto &rxn : rxns)
      if (force->pair == nullptr || (rxn.closeneigh < 0 && rxn.rmaxsq > force->pair->cutsq[rxn.iatomtype][rxn.jatomtype]))
        error->all(FLERR,"Fix bond/react: Fix bond/react cutoff is longer than pairwise cutoff");

  // need a half neighbor list, built every Nevery steps
  neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);

  lastcheck = -1;
}

/* ---------------------------------------------------------------------- */

void FixBondReact::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
  Identify all pairs of potentially reactive atoms for this time step.
  This function is modified from LAMMPS' fix bond/create.
---------------------------------------------------------------------- */

void FixBondReact::post_integrate()
{
  // update store_rxn_count on every step
  for (auto &rlm : rate_limits) {
    int rxn_count_sum = 0;
    for (auto i : rlm.rxnIDs) rxn_count_sum += rxns[i].reaction_count_total;
    rlm.store_rxn_counts.push_front(rxn_count_sum);
    rlm.store_rxn_counts.pop_back();
  }

  // check if any reactions could occur on this timestep
  int nevery_check = 1;
  for (auto &rxn : rxns) {
    if (rxn.v_nevery > -1)
      rxn.nevery = ceil(input->variable->compute_equal(rxn.v_nevery)); // NOLINT
    if (rxn.nevery <= 0)
      error->all(FLERR,"Illegal fix bond/react command: "
                 "'Nevery' must be a positive integer");
    if (!(update->ntimestep % rxn.nevery)) {
      nevery_check = 0;
      break;
    }
  }

  for (auto &rxn : rxns) {
    rxn.reaction_count = 0;
    rxn.local_rxn_count = 0;
    rxn.ghostly_rxn_count = 0;
    rxn.nlocalkeep = INT_MAX;
    rxn.nghostlykeep = INT_MAX;
    if (rxn.v_prob > -1) rxn.fraction = input->variable->compute_equal(rxn.v_prob);
  }

  if (nevery_check) {
    unlimit_bond();
    return;
  }

  // acquire updated ghost atom positions
  // necessary b/c are calling this after integrate, but before Verlet comm

  comm->forward_comm();

  // resize bond partner list and initialize it
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(partner);
    memory->destroy(finalpartner);
    memory->destroy(distsq);
    nmax = atom->nmax;
    memory->create(partner,nmax,"bond/react:partner");
    memory->create(finalpartner,nmax,"bond/react:finalpartner");
    memory->create(distsq,nmax,2,"bond/react:distsq");
  }

  // reset 'rxn_attempt' counts
  for (auto &rxn : rxns) rxn.attempts.clear();
  // reset per-bond compute map flag
  atoms2bondflag = 0;

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  // loop over neighbors of my atoms
  // each atom sets one closest eligible partner atom ID to bond with

  tagint *tag = atom->tag;
  int *type = atom->type;

  neighbor->build_one(list);

  // here we define a full special list
  // may need correction for unusual special bond settings
  nxspecial = atom->nspecial;
  xspecial = atom->special;

  // check if we are over rate_limits limits
  std::vector<int> rxn_limit_flag(rxns.size(), 1);
  for (auto &rlm : rate_limits) {
    int myrxn_count = rlm.store_rxn_counts[rlm.Nsteps-1];
    int nrxns_delta, my_nrate;
    if (myrxn_count != -1) {
      int rxn_count_sum = 0;
      for (auto i : rlm.rxnIDs) rxn_count_sum += rxns[i].reaction_count_total;
      nrxns_delta = rxn_count_sum - myrxn_count;
      if (rlm.var_flag == 1) {
        my_nrate = input->variable->compute_equal(rlm.var_id); // NOLINT
      } else my_nrate = rlm.Nlimit;
    }
    if (myrxn_count == -1 || nrxns_delta >= my_nrate)
      for (auto i : rlm.rxnIDs) rxn_limit_flag[i] = 0;
  }
  for (auto &maxlimit : max_rxn_limits) {
    int rxn_count_sum = 0;
    for (auto i : maxlimit.rxnIDs) rxn_count_sum += rxns[i].reaction_count_total;
    if (rxn_count_sum >= maxlimit.max_rxn)
      for (auto i : maxlimit.rxnIDs) rxn_limit_flag[i] = 0;
  }

  int j;
  for (auto &rxn : rxns) {
    if ((update->ntimestep % rxn.nevery) ||
        (rxn_limit_flag[rxn.ID] == 0)) continue;

    for (int ii = 0; ii < nall; ii++) {
      partner[ii] = 0;
      finalpartner[ii] = 0;
      distsq[ii][0] = 0.0;
      distsq[ii][1] = BIG;
    }

    // fork between far and close_partner here
    rxnptr = &rxn; // for reverse_comm
    if (rxn.closeneigh < 0) {
      far_partner(rxn);
      // reverse comm of distsq and partner
      // not needed if newton_pair off since I,J pair was seen by both procs
      commflag = 2;
      if (force->newton_pair) comm->reverse_comm(this);
    } else {
      close_partner(rxn);
      commflag = 2;
      comm->reverse_comm(this);
    }

    // each atom now knows its winning partner
    // forward comm of partner, so ghosts have it

    commflag = 2;
    comm->forward_comm(this,1);

    // consider for reaction:
    // only if both atoms list each other as winning bond partner
    // if other atom is owned by another proc, it should do same thing

    int temp_nattempt = 0;
    for (int i = 0; i < nlocal; i++) {
      if (partner[i] == 0) {
        continue;
      }

      j = atom->map(partner[i]);
      if (partner[j] != tag[i]) {
        continue;
      }

      // store final bond partners and count the rxn possibility once

      finalpartner[i] = tag[j];
      finalpartner[j] = tag[i];

      if (tag[i] < tag[j]) temp_nattempt++;
    }

    // cycle loop if no even eligible bonding atoms were found (on any proc)
    int some_chance;
    MPI_Allreduce(&temp_nattempt,&some_chance,1,MPI_INT,MPI_SUM,world);
    if (!some_chance) continue;

    // communicate final partner

    commflag = 3;
    comm->forward_comm(this);

    // add instance to 'attempt' only if this processor
    // owns the atoms with smaller global ID
    // NOTE: we no longer care about ghost-ghost instances as bond/create did
    // this is because we take care of updating topology later (and differently)
    for (int i = 0; i < nlocal; i++) {

      if (finalpartner[i] == 0) continue;

      j = atom->map(finalpartner[i]);
      if (tag[i] < tag[j]) {
        // to ensure types remain in same order
        if (rxn.iatomtype == type[i]) {
          rxn.attempts.push_back({tag[i], finalpartner[i]});
          // add another attempt if initiator atoms are same type
          if (rxn.iatomtype == rxn.jatomtype) rxn.attempts.push_back({finalpartner[i], tag[i]});
        } else {
          rxn.attempts.push_back({finalpartner[i], tag[i]});
        }
      }
    }
  }

  // break loop if no even eligible bonding atoms were found (on any proc)
  int some_chance;

  allnattempt = 0;
  for (auto &rxn : rxns)
    allnattempt += rxn.attempts.size();

  MPI_Allreduce(&allnattempt,&some_chance,1,MPI_INT,MPI_SUM,world);
  if (!some_chance) {
    unlimit_bond();
    return;
  }

  // evaluate custom constraint variable values here and forward_comm
  get_customvars();
  commflag = 1;
  comm->forward_comm(this,ncustomvars);

  // run through the superimpose algorithm
  // this checks if simulation topology matches unreacted mol template
  superimpose_algorithm();
  // free atoms that have been limited after reacting
  unlimit_bond();
}

/* ----------------------------------------------------------------------
  Search non-bonded neighbor lists if bonding atoms are not in special list
------------------------------------------------------------------------- */

void FixBondReact::far_partner(Reaction &rxn)
{
  int inum,jnum,itype,jtype,possible;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // loop over neighbors of my atoms
  // each atom sets one closest eligible partner atom ID to bond with

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int *type = atom->type;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // per-atom property indicating if in bond/react master group
  int flag,cols;
  int index1 = atom->find_custom("limit_tags",flag,cols);
  int *i_limit_tags = atom->ivector[index1];

  int i,j;

  for (int ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & rxn.groupbits)) continue;
    if (i_limit_tags[i] != 0) continue;
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      if (!(mask[j] & rxn.groupbits)) {
        continue;
      }

      if (i_limit_tags[j] != 0) {
        continue;
      }

      if (rxn.molecule_keyword == Molecule_Keys::INTER) {
        if (atom->molecule[i] == atom->molecule[j]) continue;
      } else if (rxn.molecule_keyword == Molecule_Keys::INTRA) {
        if (atom->molecule[i] != atom->molecule[j]) continue;
      }

      jtype = type[j];
      possible = 0;
      if (itype == rxn.iatomtype && jtype == rxn.jatomtype) {
        possible = 1;
      } else if (itype == rxn.jatomtype && jtype == rxn.iatomtype) {
        possible = 1;
      }

      if (possible == 0) continue;

      // do not allow bonding atoms within special list
      for (int k = 0; k < nxspecial[i][2]; k++)
        if (xspecial[i][k] == tag[j]) possible = 0;
      if (!possible) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      domain->minimum_image(FLERR, delx,dely,delz); // ghost location fix
      rsq = delx*delx + dely*dely + delz*delz;

      if (rxn.v_rmin > -1) {
        double cutoff = input->variable->compute_equal(rxn.v_rmin);
        rxn.rminsq = cutoff*cutoff;
      }
      if (rxn.v_rmax > -1) {
        double cutoff = input->variable->compute_equal(rxn.v_rmax);
        rxn.rmaxsq = cutoff*cutoff;
      }
      if (rsq >= rxn.rmaxsq || rsq <= rxn.rminsq) {
        continue;
      }
      if (rsq < distsq[i][1]) {
        partner[i] = tag[j];
        distsq[i][1] = rsq;
      }
      if (rsq < distsq[j][1]) {
        partner[j] = tag[i];
        distsq[j][1] = rsq;
      }
    }
  }
}

/* ----------------------------------------------------------------------
  Slightly simpler to find bonding partner when a close neighbor
------------------------------------------------------------------------- */

void FixBondReact::close_partner(Reaction &rxn)
{
  int n,i1,i2,itype,jtype;
  double delx,dely,delz,rsq;

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;

  // per-atom property indicating if in bond/react master group
  int flag,cols;
  int index1 = atom->find_custom("limit_tags",flag,cols);
  int *i_limit_tags = atom->ivector[index1];

  // loop over special list
  for (int ii = 0; ii < atom->nlocal; ii++) {
    itype = type[ii];
    n = 0;
    if (rxn.closeneigh != 0)
      n = nxspecial[ii][rxn.closeneigh-1];
    for (; n < nxspecial[ii][rxn.closeneigh]; n++) {
      i1 = ii;
      i2 = atom->map(xspecial[ii][n]);
      jtype = type[i2];
      if (!(mask[i1] & rxn.groupbits)) continue;
      if (!(mask[i2] & rxn.groupbits)) continue;
      if (i_limit_tags[i1] != 0) continue;
      if (i_limit_tags[i2] != 0) continue;
      if (itype != rxn.iatomtype || jtype != rxn.jatomtype) continue;

      if (rxn.molecule_keyword == Molecule_Keys::INTER) {
        if (atom->molecule[i1] == atom->molecule[i2]) continue;
      } else if (rxn.molecule_keyword == Molecule_Keys::INTRA) {
        if (atom->molecule[i1] != atom->molecule[i2]) continue;
      }

      delx = x[i1][0] - x[i2][0];
      dely = x[i1][1] - x[i2][1];
      delz = x[i1][2] - x[i2][2];
      domain->minimum_image(FLERR, delx,dely,delz); // ghost location fix
      rsq = delx*delx + dely*dely + delz*delz;

      if (rxn.v_rmin > -1) {
        double cutoff = input->variable->compute_equal(rxn.v_rmin);
        rxn.rminsq = cutoff*cutoff;
      }
      if (rxn.v_rmax > -1) {
        double cutoff = input->variable->compute_equal(rxn.v_rmax);
        rxn.rmaxsq = cutoff*cutoff;
      }
      if (rsq >= rxn.rmaxsq || rsq <= rxn.rminsq) continue;

      if (rxn.closeneigh == 0) {
        if (rsq > distsq[i1][0]) {
          partner[i1] = tag[i2];
          distsq[i1][0] = rsq;
        }
        if (rsq > distsq[i2][0]) {
          partner[i2] = tag[i1];
          distsq[i2][0] = rsq;
        }
      } else {
        if (rsq < distsq[i1][1]) {
          partner[i1] = tag[i2];
          distsq[i1][1] = rsq;
        }
        if (rsq < distsq[i2][1]) {
          partner[i2] = tag[i1];
          distsq[i2][1] = rsq;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
  Set up global variables. Loop through all pairs; loop through Pioneers
  until Superimpose Algorithm is completed for each pair.
  various statuses of superimpose algorithm:
  ACCEPT: site successfully matched to pre-reacted template
  REJECT: site does not match pre-reacted template
  PROCEED: normal execution (non-guessing mode)
  CONTINUE: a neighbor has been assigned, skip to next neighbor
  GUESSFAIL: a guess has failed (if no more restore points, status = 'REJECT')
  RESTORE: restore mode, load most recent restore point
------------------------------------------------------------------------- */

void FixBondReact::superimpose_algorithm()
{
  const int nprocs = comm->nprocs;
  my_num_mega = 0;
  local_num_mega = 0;
  ghostly_num_mega = 0;

  // indicates local ghosts of other procs
  int tmp;
  localsendlist = (int *) comm->extract("localsendlist",tmp);

  // quick description of important global indices you'll see floating about:
  // pion: the pioneer loop index
  // neigh: in the first neighbor index
  // trace: retraces the first nieghbors. once you choose a first neighbor, you then check for other nieghbors of same type
  // pioneers: during Superimpose Algorithm, atoms which have been assigned, but whose first neighbors haven't
  // glove: global IDs. index indicates is pre-reaction ID-1, value is mapped sim atom ID
  // glove_counter: used to determine when to terminate Superimpose Algorithm

  Superimpose super;
  Superimpose::StatePoint &sp = super.sp;
  int &avail_guesses = super.avail_guesses;
  std::vector<int> &guess_branch = super.guess_branch;
  guess_branch.resize(MAXGUESS, 0);

  sp.glove.resize(max_natoms);
  sp.pioneers.resize(max_natoms);
  sp.pioneer_count.resize(max_natoms);

  restore_pts.resize(MAXGUESS);
  for (auto &restore_pt : restore_pts) {
    restore_pt.glove.resize(max_natoms);
    restore_pt.pioneers.resize(max_natoms);
    restore_pt.pioneer_count.resize(max_natoms);
  }
  my_mega_glove.resize(max_natoms+cuff); // mega_glove indexing seems inside-out
  for (auto &vec : my_mega_glove) vec.resize(allnattempt, 0);

  // let's finally begin the superimpose loop
  for (auto &rxn : rxns) {
    for (auto &rxn_attempt : rxn.attempts) {

      status = Status::PROCEED;

      sp.pion = sp.neigh = sp.trace = sp.glove_counter = 0;
      std::fill(sp.glove.begin(), sp.glove.end(), 0);
      std::fill(guess_branch.begin(), guess_branch.end(), 0);

      sp.glove[rxn.ibonding-1] = rxn_attempt[0];
      sp.glove_counter++;
      sp.glove[rxn.jbonding-1] = rxn_attempt[1];
      sp.glove_counter++;

      // special case, only two atoms in reaction templates
      // then: bonding reactant nspecials guaranteed to be equal, and either 0 or 1
      if (sp.glove_counter == rxn.reactant->natoms) {
        tagint local_atom1 = atom->map(sp.glove[rxn.ibonding-1]);
        tagint local_atom2 = atom->map(sp.glove[rxn.jbonding-1]);
        if ( (nxspecial[local_atom1][0] == rxn.reactant->nspecial[rxn.ibonding-1][0] &&
              nxspecial[local_atom2][0] == nxspecial[local_atom1][0]) &&
             (nxspecial[local_atom1][0] == 0 ||
              xspecial[local_atom1][0] == atom->tag[local_atom2]) &&
             check_constraints(rxn, sp.glove)) {
          if (rxn.fraction < 1.0 &&
              random[rxn.ID]->uniform() >= rxn.fraction) {
            status = Status::REJECT;
          } else {
            status = Status::ACCEPT;
            my_mega_glove[0][my_num_mega] = (double) rxn.ID;
            if (rxn.rescale_charges_flag) my_mega_glove[1][my_num_mega] = get_totalcharge(rxn, sp.glove);
            for (int i = 0; i < rxn.reactant->natoms; i++) {
              my_mega_glove[i+cuff][my_num_mega] = (double) sp.glove[i];
            }
            my_num_mega++;
          }
        } else status = Status::REJECT;
      }

      avail_guesses = 0;

      std::fill(sp.pioneer_count.begin(), sp.pioneer_count.end(), 0);

      for (int i = 0; i < rxn.reactant->nspecial[rxn.ibonding-1][0]; i++)
        sp.pioneer_count[rxn.reactant->special[rxn.ibonding-1][i]-1]++;

      for (int i = 0; i < rxn.reactant->nspecial[rxn.jbonding-1][0]; i++)
        sp.pioneer_count[rxn.reactant->special[rxn.jbonding-1][i]-1]++;


      int hang_catch = 0;
      while (status != Status::ACCEPT && status != Status::REJECT) {

        for (int i = 0; i < max_natoms; i++) sp.pioneers[i] = 0;

        for (int i = 0; i < rxn.reactant->natoms; i++) {
          if (sp.glove[i] != 0 && sp.pioneer_count[i] < rxn.reactant->nspecial[i][0] && rxn.atoms[i].edge == 0) {
            sp.pioneers[i] = 1;
          }
        }

        // run through the pioneers
        // due to use of restore points, 'pion' index can change in loop
        for (sp.pion = 0; sp.pion < rxn.reactant->natoms; sp.pion++) {
          if (sp.pioneers[sp.pion] || status == Status::GUESSFAIL) {
            make_a_guess(super, rxn);
            if (status == Status::ACCEPT || status == Status::REJECT) break;
          }
        }

        // reaction site found successfully!
        if (status == Status::ACCEPT) {
          if (rxn.fraction < 1.0 &&
              random[rxn.ID]->uniform() >= rxn.fraction) status = Status::REJECT;
          else {
            my_mega_glove[0][my_num_mega] = (double) rxn.ID;
            if (rxn.rescale_charges_flag) my_mega_glove[1][my_num_mega] = get_totalcharge(rxn, sp.glove);
            for (int i = 0; i < rxn.reactant->natoms; i++) {
              my_mega_glove[i+cuff][my_num_mega] = (double) sp.glove[i];
            }
            my_num_mega++;
          }
        }
        hang_catch++;
        // let's go ahead and catch the simplest of hangs
        //if (hang_catch > rxn.reactant->natoms*4)
        if (hang_catch > atom->nlocal*30) {
          error->one(FLERR,"Fix bond/react: Excessive iteration of superimpose algorithm. "
              "Please check that all pre-reaction template atoms are linked to an initiator atom, "
              "via at least one path that does not involve edge atoms.");
        }
      }
    }
  }

  global_megasize = 0;

  local_mega_glove.resize(max_natoms+cuff); // mega_glove indexing seems inside-out
  for (auto &vec : local_mega_glove) vec.resize(my_num_mega, 0);
  ghostly_mega_glove.resize(max_natoms+cuff); // mega_glove indexing seems inside-out
  for (auto &vec : ghostly_mega_glove) vec.resize(my_num_mega, 0);

  dedup_mega_gloves(Dedup_Modes::LOCAL); // make sure atoms aren't added to more than one reaction
  glove_ghostcheck(); // split into 'local' and 'global'
  ghost_glovecast(); // consolidate all mega_gloves to all processors

  std::vector<int> mpi_send(rxns.size()), mpi_recv(rxns.size());
  for (auto &rxn : rxns) mpi_send[rxn.ID] = rxn.local_rxn_count;
  MPI_Allreduce(mpi_send.data(), mpi_recv.data(), rxns.size(), MPI_INT, MPI_SUM, world);
  for (auto &rxn : rxns) rxn.reaction_count = mpi_send[rxn.ID];

  int rxnflag = 0;
  int *delta_rxn;
  memory->create(delta_rxn, rxns.size(), "bond/react:delta_rxn");
  if (comm->me == 0)
    for (auto &rxn : rxns) {
      delta_rxn[rxn.ID] = rxn.reaction_count + rxn.ghostly_rxn_count;
      rxnflag += delta_rxn[rxn.ID];
    }

  MPI_Bcast(&delta_rxn[0], rxns.size(), MPI_INT, 0, world);
  MPI_Bcast(&rxnflag, 1, MPI_INT, 0, world);

  if (!rxnflag) return;

  // C++11 and later compatible version of Park pRNG
  std::minstd_rand park_rng;
  if (shuffle_seed == 0) {
    std::random_device rnd;
    park_rng.seed(rnd());
  } else {
    park_rng.seed(shuffle_seed);
  }

  std::vector<int> oversteps(rxns.size(), 0);
  if (comm->me == 0) {
    // check if we overstepped our reaction limit, via either max_rxn or rate_limit
    for (auto maxlimit : max_rxn_limits) {
      int rxn_count_sum = 0;
      int delta_rxn_sum = 0;
      for (auto i : maxlimit.rxnIDs) {
        rxn_count_sum += rxns[i].reaction_count_total;
        delta_rxn_sum += delta_rxn[i];
      }
      int max_limit_overstep_sum = rxn_count_sum + delta_rxn_sum - maxlimit.max_rxn;
      if (max_limit_overstep_sum > 0) {
        if (maxlimit.Nrxns == 1) {
          oversteps[maxlimit.rxnIDs[0]] = MAX(oversteps[maxlimit.rxnIDs[0]], max_limit_overstep_sum);
        } else {
          std::vector<int> dummy_list;
          for (auto i : maxlimit.rxnIDs)
            for (int j = 0; j < delta_rxn[i]; j++)
              dummy_list.push_back(i);
          std::shuffle(dummy_list.begin(), dummy_list.end(), park_rng);
          std::vector<int> max_limit_overstep(rxns.size(),0);
          for (int i = 0; i < max_limit_overstep_sum; i++)
            max_limit_overstep[dummy_list[i]]++;
          for (auto &rxn : rxns)
            oversteps[rxn.ID] = MAX(oversteps[rxn.ID], max_limit_overstep[rxn.ID]);
        }
      }
    }

    for (auto rlm : rate_limits) {
      int myrxn_count = rlm.store_rxn_counts[rlm.Nsteps-1];
      if (myrxn_count != -1) {
        int rxn_count_sum = 0;
        int delta_rxn_sum = 0;
        for (auto i : rlm.rxnIDs) {
          rxn_count_sum += rxns[i].reaction_count_total;
          delta_rxn_sum += delta_rxn[i];
        }
        int nrxn_delta = rxn_count_sum + delta_rxn_sum - myrxn_count;
        int my_nrate;
        if (rlm.var_flag == 1) {
          my_nrate = input->variable->compute_equal(rlm.var_id); // NOLINT
        } else my_nrate = rlm.Nlimit;
        int rate_limit_overstep_sum = nrxn_delta - my_nrate;
        if (rate_limit_overstep_sum > 0) {
          if (rlm.Nrxns == 1) {
            oversteps[rlm.rxnIDs[0]] = MAX(oversteps[rlm.rxnIDs[0]], rate_limit_overstep_sum);
          } else {
            std::vector<int> dummy_list;
            for (auto i : rlm.rxnIDs)
              for (int j = 0; j < delta_rxn[i]; j++)
                dummy_list.push_back(i);
            std::shuffle(dummy_list.begin(), dummy_list.end(), park_rng);
            std::vector<int> rate_limit_overstep(rxns.size(),0);
            for (int i = 0; i < rate_limit_overstep_sum; i++)
              rate_limit_overstep[dummy_list[i]]++;
            for (auto &rxn : rxns)
              oversteps[rxn.ID] = MAX(oversteps[rxn.ID], rate_limit_overstep[rxn.ID]);
          }
        }
      }
    }
  }
  MPI_Bcast(oversteps.data(),rxns.size(),MPI_INT,0,world);

  for (auto &rxn : rxns) {
    if (oversteps[rxn.ID] > 0) {
      // let's randomly choose rxns to skip, unbiasedly from local and ghostly
      int *local_rxncounts;
      int *all_localkeep;
      memory->create(local_rxncounts,nprocs,"bond/react:local_rxncounts");
      memory->create(all_localkeep,nprocs,"bond/react:all_localkeep");
      MPI_Gather(&rxn.local_rxn_count,1,MPI_INT,local_rxncounts,1,MPI_INT,0,world);
      if (comm->me == 0) {
        // when using variable input for rate_limit, rate_limit_overstep could be > delta_rxn (below)
        // we need to limit overstep to the number of reactions on this timestep
        // essentially skipping all reactions, would be more efficient to use a skip_all flag
        if (oversteps[rxn.ID] > delta_rxn[rxn.ID]) oversteps[rxn.ID] = delta_rxn[rxn.ID];
        int nkeep = delta_rxn[rxn.ID] - oversteps[rxn.ID];
        int *rxn_by_proc;
        memory->create(rxn_by_proc,delta_rxn[rxn.ID],"bond/react:rxn_by_proc");
        for (int j = 0; j < delta_rxn[rxn.ID]; j++)
          rxn_by_proc[j] = -1; // corresponds to ghostly
        int itemp = 0;
        for (int j = 0; j < nprocs; j++)
          for (int k = 0; k < local_rxncounts[j]; k++)
            rxn_by_proc[itemp++] = j;
        std::shuffle(&rxn_by_proc[0],&rxn_by_proc[delta_rxn[rxn.ID]], park_rng);
        for (int j = 0; j < nprocs; j++)
          all_localkeep[j] = 0;
        rxn.nghostlykeep = 0;
        for (int j = 0; j < nkeep; j++) {
          if (rxn_by_proc[j] == -1) rxn.nghostlykeep++;
          else all_localkeep[rxn_by_proc[j]]++;
        }
        memory->destroy(rxn_by_proc);
      }
      MPI_Scatter(&all_localkeep[0],1,MPI_INT,&rxn.nlocalkeep,1,MPI_INT,0,world);
      MPI_Bcast(&rxn.nghostlykeep,1,MPI_INT,0,world);
      memory->destroy(local_rxncounts);
      memory->destroy(all_localkeep);
    }
  }
  memory->destroy(delta_rxn);

  // this updates topology next step
  next_reneighbor = update->ntimestep;

  update_everything(); // change topology
}

/* ----------------------------------------------------------------------
  Screen for obvious algorithm fails. This is the return point when a guess
  has failed: check for available restore points.
------------------------------------------------------------------------- */

void FixBondReact::make_a_guess(Superimpose &super, Reaction &rxn)
{
  Superimpose::StatePoint &sp = super.sp;
  int &avail_guesses = super.avail_guesses;

  int *type = atom->type;
  int nfirst_neighs = rxn.reactant->nspecial[sp.pion][0];

  // per-atom property indicating if in bond/react master group
  int flag,cols;
  int index1 = atom->find_custom("limit_tags",flag,cols);
  int *i_limit_tags = atom->ivector[index1];

  if (status == Status::GUESSFAIL && avail_guesses == 0) {
    status = Status::REJECT;
    return;
  }

  if (status == Status::GUESSFAIL && avail_guesses > 0) {
    // load restore point
    for (int i = 0; i < rxn.reactant->natoms; i++) {
      sp.glove[i] = restore_pts[avail_guesses-1].glove[i];
      sp.pioneer_count[i] = restore_pts[avail_guesses-1].pioneer_count[i];
      sp.pioneers[i] = restore_pts[avail_guesses-1].pioneers[i];
    }
    sp.pion = restore_pts[avail_guesses-1].pion;
    sp.neigh = restore_pts[avail_guesses-1].neigh;
    sp.trace = restore_pts[avail_guesses-1].trace;
    sp.glove_counter = restore_pts[avail_guesses-1].glove_counter;
    status = Status::RESTORE;
    neighbor_loop(super, rxn);
    if (status != Status::PROCEED) return;
  }

  nfirst_neighs = rxn.reactant->nspecial[sp.pion][0];

  //  check if any of first neighbors are in bond_react_MASTER_group
  //  if so, this constitutes a fail
  //  because still undergoing a previous reaction!
  //  could technically fail unnecessarily during a wrong guess if near edge atoms
  //  we accept this temporary and infrequent decrease in reaction occurrences

  for (int i = 0; i < nxspecial[atom->map(sp.glove[sp.pion])][0]; i++) {
    if (atom->map(xspecial[atom->map(sp.glove[sp.pion])][i]) < 0) {
      error->one(FLERR,"Fix bond/react: Fix bond/react needs ghost atoms from further away"); // parallel issues.
    }
    if (i_limit_tags[(int)atom->map(xspecial[atom->map(sp.glove[sp.pion])][i])] != 0) {
      status = Status::GUESSFAIL;
      return;
    }
  }

  // check for same number of neighbors between unreacted mol and simulation
  if (nfirst_neighs != nxspecial[atom->map(sp.glove[sp.pion])][0]) {
    status = Status::GUESSFAIL;
    return;
  }

  // make sure all neighbors aren't already assigned
  // an issue discovered for coarse-grained example
  int assigned_count = 0;
  for (int i = 0; i < nfirst_neighs; i++)
    for (int j = 0; j < rxn.reactant->natoms; j++)
      if (xspecial[atom->map(sp.glove[sp.pion])][i] == sp.glove[j]) {
        assigned_count++;
        break;
      }

  if (assigned_count == nfirst_neighs) status = Status::GUESSFAIL;

  // check if all neigh atom types are the same between simulation and unreacted mol
  int *mol_ntypes = new int[atom->ntypes];
  int *lcl_ntypes = new int[atom->ntypes];

  for (int i = 0; i < atom->ntypes; i++) {
    mol_ntypes[i] = 0;
    lcl_ntypes[i] = 0;
  }

  for (int i = 0; i < nfirst_neighs; i++) {
    mol_ntypes[(int)rxn.reactant->type[(int)rxn.reactant->special[sp.pion][i]-1]-1]++;
    lcl_ntypes[(int)type[(int)atom->map(xspecial[atom->map(sp.glove[sp.pion])][i])]-1]++; //added -1
  }

  for (int i = 0; i < atom->ntypes; i++) {
    if (mol_ntypes[i] != lcl_ntypes[i]) {
      status = Status::GUESSFAIL;
      delete[] mol_ntypes;
      delete[] lcl_ntypes;
      return;
    }
  }

  delete[] mol_ntypes;
  delete[] lcl_ntypes;

  // okay everything seems to be in order. let's assign some ID pairs!!!
  neighbor_loop(super, rxn);
}

/* ----------------------------------------------------------------------
  Loop through all First Bonded Neighbors of the current Pioneer.
  Prepare appropriately if we are in Restore Mode.
------------------------------------------------------------------------- */

void FixBondReact::neighbor_loop(Superimpose &super, Reaction &rxn)
{
  Superimpose::StatePoint &sp = super.sp;

  int nfirst_neighs = rxn.reactant->nspecial[sp.pion][0];

  if (status == Status::RESTORE) {
    check_a_neighbor(super, rxn);
    return;
  }

  for (sp.neigh = 0; sp.neigh < nfirst_neighs; sp.neigh++) {
    if (sp.glove[(int)rxn.reactant->special[sp.pion][sp.neigh]-1] == 0) {
      check_a_neighbor(super, rxn);
    }
  }
  // status should still = PROCEED
}

/* ----------------------------------------------------------------------
  Check if we can assign this First Neighbor to pre-reacted template
  without guessing. If so, do it! If not, call crosscheck_the_nieghbor().
------------------------------------------------------------------------- */

void FixBondReact::check_a_neighbor(Superimpose &super, Reaction &rxn)
{
  Superimpose::StatePoint &sp = super.sp;

  int *type = atom->type;
  int nfirst_neighs = rxn.reactant->nspecial[sp.pion][0];

  if (status != Status::RESTORE) {
    // special consideration for hydrogen atoms (and all first neighbors bonded to no other atoms) (and aren't edge atoms)
    if (rxn.reactant->nspecial[(int)rxn.reactant->special[sp.pion][sp.neigh]-1][0] == 1 && rxn.atoms[(int)rxn.reactant->special[sp.pion][sp.neigh]-1].edge == 0) {

      for (int i = 0; i < nfirst_neighs; i++) {

        if (type[(int)atom->map(xspecial[(int)atom->map(sp.glove[sp.pion])][i])] == rxn.reactant->type[(int)rxn.reactant->special[sp.pion][sp.neigh]-1] &&
            nxspecial[(int)atom->map(xspecial[(int)atom->map(sp.glove[sp.pion])][i])][0] == 1) {

          int already_assigned = 0;
          for (int j = 0; j < rxn.reactant->natoms; j++) {
            if (sp.glove[j] == xspecial[atom->map(sp.glove[sp.pion])][i]) {
              already_assigned = 1;
              break;
            }
          }

          if (already_assigned == 0) {
            sp.glove[(int)rxn.reactant->special[sp.pion][sp.neigh]-1] = xspecial[(int)atom->map(sp.glove[sp.pion])][i];

            //another check for ghost atoms. perhaps remove the one in make_a_guess
            if (atom->map(sp.glove[(int)rxn.reactant->special[sp.pion][sp.neigh]-1]) < 0) {
              error->one(FLERR,"Fix bond/react: Fix bond/react needs ghost atoms from further away");
            }

            for (int j = 0; j < rxn.reactant->nspecial[rxn.reactant->special[sp.pion][sp.neigh]-1][0]; j++) {
              sp.pioneer_count[rxn.reactant->special[rxn.reactant->special[sp.pion][sp.neigh]-1][j]-1]++;
            }

            sp.glove_counter++;
            if (sp.glove_counter == rxn.reactant->natoms) {
              if (ring_check(rxn, sp.glove) && check_constraints(rxn, sp.glove)) status = Status::ACCEPT;
              else status = Status::GUESSFAIL;
              return;
            }
            // status should still == PROCEED
            return;
          }
        }
      }
      // we are here if no matching atom found
      status = Status::GUESSFAIL;
      return;
    }
  }

  crosscheck_the_neighbor(super, rxn);
  if (status != Status::PROCEED) {
    if (status == Status::CONTINUE)
      status = Status::PROCEED;
    return;
  }

  // finally ready to match non-duplicate, non-edge atom IDs!!

  for (int i = 0; i < nfirst_neighs; i++) {

    if (type[atom->map((int)xspecial[(int)atom->map(sp.glove[sp.pion])][i])] == rxn.reactant->type[(int)rxn.reactant->special[sp.pion][sp.neigh]-1]) {
      int already_assigned = 0;

      //check if a first neighbor of the pioneer is already assigned to pre-reacted template
      for (int j = 0; j < rxn.reactant->natoms; j++) {
        if (sp.glove[j] == xspecial[atom->map(sp.glove[sp.pion])][i]) {
          already_assigned = 1;
          break;
        }
      }

      if (already_assigned == 0) {
        sp.glove[(int)rxn.reactant->special[sp.pion][sp.neigh]-1] = xspecial[(int)atom->map(sp.glove[sp.pion])][i];

        //another check for ghost atoms. perhaps remove the one in make_a_guess
        if (atom->map(sp.glove[(int)rxn.reactant->special[sp.pion][sp.neigh]-1]) < 0) {
          error->one(FLERR,"Fix bond/react: Fix bond/react needs ghost atoms from further away");
        }

        for (int ii = 0; ii < rxn.reactant->nspecial[rxn.reactant->special[sp.pion][sp.neigh]-1][0]; ii++) {
          sp.pioneer_count[rxn.reactant->special[rxn.reactant->special[sp.pion][sp.neigh]-1][ii]-1]++;
        }

        sp.glove_counter++;
        if (sp.glove_counter == rxn.reactant->natoms) {
          if (ring_check(rxn, sp.glove) && check_constraints(rxn, sp.glove)) status = Status::ACCEPT;
          else status = Status::GUESSFAIL;
          return;
          // will never complete here when there are edge atoms
          // ...actually that could be wrong if people get creative...shouldn't affect anything
        }
        // status should still = PROCEED
        return;
      }
    }
  }
  // status is still 'PROCEED' if we are here!
}

/* ----------------------------------------------------------------------
  Check if there a viable guess to be made. If so, prepare to make a
  guess by recording a restore point.
------------------------------------------------------------------------- */

void FixBondReact::crosscheck_the_neighbor(Superimpose &super, Reaction &rxn)
{
  Superimpose::StatePoint &sp = super.sp;
  int &avail_guesses = super.avail_guesses;

  int nfirst_neighs = rxn.reactant->nspecial[sp.pion][0];

  if (status == Status::RESTORE) {
    inner_crosscheck_loop(super, rxn);
    return;
  }

  for (sp.trace = 0; sp.trace < nfirst_neighs; sp.trace++) {
    if (sp.neigh != sp.trace && rxn.reactant->type[(int)rxn.reactant->special[sp.pion][sp.neigh]-1] == rxn.reactant->type[(int)rxn.reactant->special[sp.pion][sp.trace]-1] &&
        sp.glove[rxn.reactant->special[sp.pion][sp.trace]-1] == 0) {

      if (avail_guesses == MAXGUESS) {
        error->warning(FLERR,"Fix bond/react: Fix bond/react failed because MAXGUESS set too small. ask developer for info");
        status = Status::GUESSFAIL;
        return;
      }
      avail_guesses++;
      for (int i = 0; i < rxn.reactant->natoms; i++) {
        restore_pts[avail_guesses-1].glove[i] = sp.glove[i];
        restore_pts[avail_guesses-1].pioneer_count[i] = sp.pioneer_count[i];
        restore_pts[avail_guesses-1].pioneers[i] = sp.pioneers[i];
      }
      restore_pts[avail_guesses-1].pion = sp.pion;
      restore_pts[avail_guesses-1].neigh = sp.neigh;
      restore_pts[avail_guesses-1].trace = sp.trace;
      restore_pts[avail_guesses-1].glove_counter = sp.glove_counter;

      inner_crosscheck_loop(super, rxn);
      return;
    }
  }
  // status is still 'PROCEED' if we are here!
}

/* ----------------------------------------------------------------------
  We are ready to make a guess. If there are multiple possible choices
  for this guess, keep track of these.
------------------------------------------------------------------------- */

void FixBondReact::inner_crosscheck_loop(Superimpose &super, Reaction &rxn)
{
  Superimpose::StatePoint &sp = super.sp;
  int &avail_guesses = super.avail_guesses;
  std::vector<int> &guess_branch = super.guess_branch;

  int *type = atom->type;
  // arbitrarily limited to 5 identical first neighbors
  tagint tag_choices[5];
  int nfirst_neighs = rxn.reactant->nspecial[sp.pion][0];

  int num_choices = 0;
  for (int i = 0; i < nfirst_neighs; i++) {
    if (type[(int)atom->map(xspecial[atom->map(sp.glove[sp.pion])][i])] == rxn.reactant->type[(int)rxn.reactant->special[sp.pion][sp.neigh]-1]) {
      if (num_choices == 5) { // here failed because too many identical first neighbors. but really no limit if situation arises
        status = Status::GUESSFAIL;
        return;
      }
      tag_choices[num_choices++] = xspecial[atom->map(sp.glove[sp.pion])][i];
    }
  }

  // guess branch is for when multiple identical neighbors. then we guess each one in turn
  // guess_branch must work even when avail_guesses = 0. so index accordingly!
  // ...actually, avail_guesses should never be zero here anyway
  if (guess_branch[avail_guesses-1] == 0) guess_branch[avail_guesses-1] = num_choices;

  for (int i=1; i < num_choices; ++i) {
    tagint hold = tag_choices[i];
    int j = i - 1;
    while ((j >= 0) && (tag_choices[j] > hold)) {
      tag_choices[j+1] = tag_choices[j];
      --j;
    }
    tag_choices[j+1] = hold;
  }

  for (int i = guess_branch[avail_guesses-1]-1; i >= 0; i--) {
    int already_assigned = 0;
    for (int j = 0; j < rxn.reactant->natoms; j++) {
      if (sp.glove[j] == tag_choices[i]) {
        already_assigned = 1;
        break;
      }
    }
    if (already_assigned == 1) {
      guess_branch[avail_guesses-1]--;
      if (guess_branch[avail_guesses-1] == 0) {
        status = Status::REJECT;
        return;
      }
    } else {
      sp.glove[rxn.reactant->special[sp.pion][sp.neigh]-1] = tag_choices[i];
      guess_branch[avail_guesses-1]--;
      break;
    }
  }

  //another check for ghost atoms. perhaps remove the one in make_a_guess
  if (atom->map(sp.glove[(int)rxn.reactant->special[sp.pion][sp.neigh]-1]) < 0) {
    error->one(FLERR,"Fix bond/react: Fix bond/react needs ghost atoms from further away");
  }

  if (guess_branch[avail_guesses-1] == 0) avail_guesses--;

  for (int i = 0; i < rxn.reactant->nspecial[rxn.reactant->special[sp.pion][sp.neigh]-1][0]; i++) {
    sp.pioneer_count[rxn.reactant->special[rxn.reactant->special[sp.pion][sp.neigh]-1][i]-1]++;
  }
  sp.glove_counter++;
  if (sp.glove_counter == rxn.reactant->natoms) {
    if (ring_check(rxn, sp.glove) && check_constraints(rxn, sp.glove)) status = Status::ACCEPT;
    else status = Status::GUESSFAIL;
    return;
  }
  status = Status::CONTINUE;
}

/* ----------------------------------------------------------------------
  Check that newly assigned atoms have correct bonds
  Necessary for certain ringed structures
------------------------------------------------------------------------- */

int FixBondReact::ring_check(Reaction &rxn, std::vector<tagint> &glove)
{
  // ring_check can be made more efficient by re-introducing 'frozen' atoms
  // 'frozen' atoms have been assigned and also are no longer pioneers

  // double check the number of neighbors match for all non-edge atoms
  // otherwise, atoms at 'end' of symmetric ring can behave like edge atoms
  for (int i = 0; i < rxn.reactant->natoms; i++)
    if (rxn.atoms[i].edge == 0 &&
        rxn.reactant->nspecial[i][0] != nxspecial[atom->map(glove[i])][0])
      return 0;

  for (int i = 0; i < rxn.reactant->natoms; i++) {
    for (int j = 0; j < rxn.reactant->nspecial[i][0]; j++) {
      int ring_fail = 1;
      int ispecial = rxn.reactant->special[i][j];
      for (int k = 0; k < nxspecial[atom->map(glove[i])][0]; k++) {
        if (xspecial[atom->map(glove[i])][k] == glove[ispecial-1]) {
          ring_fail = 0;
          break;
        }
      }
      if (ring_fail == 1) return 0;
    }
  }
  return 1;
}

/* ----------------------------------------------------------------------
evaluate constraints: return 0 if any aren't satisfied
------------------------------------------------------------------------- */

int FixBondReact::check_constraints(Reaction &rxn, std::vector<tagint> &glove)
{
  double x1[3],x2[3],x3[3],x4[3];
  double delx,dely,delz,rsq;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double rsq1,rsq2,r1,r2,c,t,prrhob;
  // for computation of dihedrals
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,ra2inv,rb2inv,rabinv;
  double s,phi;
  int ANDgate;

  tagint atom1,atom2;
  double **x = atom->x;

  for (auto &constraint : rxn.constraints) constraint.satisfied = true;

  for (auto &constraint : rxn.constraints) {
    if (constraint.type == Reaction::Constraint::Type::DISTANCE) {
      get_IDcoords(constraint.idtypes[0], constraint.ids[0], x1, rxn.reactant, glove);
      get_IDcoords(constraint.idtypes[1], constraint.ids[1], x2, rxn.reactant, glove);
      delx = x1[0] - x2[0];
      dely = x1[1] - x2[1];
      delz = x1[2] - x2[2];
      domain->minimum_image(FLERR, delx,dely,delz); // ghost location fix
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < constraint.distance.rminsq || rsq > constraint.distance.rmaxsq) constraint.satisfied = false;
    } else if (constraint.type == Reaction::Constraint::Type::ANGLE) {
      get_IDcoords(constraint.idtypes[0], constraint.ids[0], x1, rxn.reactant, glove);
      get_IDcoords(constraint.idtypes[1], constraint.ids[1], x2, rxn.reactant, glove);
      get_IDcoords(constraint.idtypes[2], constraint.ids[2], x3, rxn.reactant, glove);

      // 1st bond
      delx1 = x1[0] - x2[0];
      dely1 = x1[1] - x2[1];
      delz1 = x1[2] - x2[2];
      domain->minimum_image(FLERR, delx1,dely1,delz1); // ghost location fix
      rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
      r1 = sqrt(rsq1);

      // 2nd bond
      delx2 = x3[0] - x2[0];
      dely2 = x3[1] - x2[1];
      delz2 = x3[2] - x2[2];
      domain->minimum_image(FLERR, delx2,dely2,delz2); // ghost location fix
      rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
      r2 = sqrt(rsq2);

      // angle (cos and sin)
      c = delx1*delx2 + dely1*dely2 + delz1*delz2;
      c /= r1*r2;
      if (c > 1.0) c = 1.0;
      if (c < -1.0) c = -1.0;
      if (acos(c) < constraint.angle.amin || acos(c) > constraint.angle.amax) constraint.satisfied = false;
    } else if (constraint.type == Reaction::Constraint::Type::DIHEDRAL) {
      // phi calculation from dihedral style harmonic
      get_IDcoords(constraint.idtypes[0], constraint.ids[0], x1, rxn.reactant, glove);
      get_IDcoords(constraint.idtypes[1], constraint.ids[1], x2, rxn.reactant, glove);
      get_IDcoords(constraint.idtypes[2], constraint.ids[2], x3, rxn.reactant, glove);
      get_IDcoords(constraint.idtypes[3], constraint.ids[3], x4, rxn.reactant, glove);

      vb1x = x1[0] - x2[0];
      vb1y = x1[1] - x2[1];
      vb1z = x1[2] - x2[2];
      domain->minimum_image(FLERR, vb1x,vb1y,vb1z);

      vb2x = x3[0] - x2[0];
      vb2y = x3[1] - x2[1];
      vb2z = x3[2] - x2[2];
      domain->minimum_image(FLERR, vb2x,vb2y,vb2z);

      vb2xm = -vb2x;
      vb2ym = -vb2y;
      vb2zm = -vb2z;
      domain->minimum_image(FLERR, vb2xm,vb2ym,vb2zm);

      vb3x = x4[0] - x3[0];
      vb3y = x4[1] - x3[1];
      vb3z = x4[2] - x3[2];
      domain->minimum_image(FLERR, vb3x,vb3y,vb3z);

      ax = vb1y*vb2zm - vb1z*vb2ym;
      ay = vb1z*vb2xm - vb1x*vb2zm;
      az = vb1x*vb2ym - vb1y*vb2xm;
      bx = vb3y*vb2zm - vb3z*vb2ym;
      by = vb3z*vb2xm - vb3x*vb2zm;
      bz = vb3x*vb2ym - vb3y*vb2xm;

      rasq = ax*ax + ay*ay + az*az;
      rbsq = bx*bx + by*by + bz*bz;
      rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
      rg = sqrt(rgsq);

      ra2inv = rb2inv = 0.0;
      if (rasq > 0) ra2inv = 1.0/rasq;
      if (rbsq > 0) rb2inv = 1.0/rbsq;
      rabinv = sqrt(ra2inv*rb2inv);

      c = (ax*bx + ay*by + az*bz)*rabinv;
      s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

      if (c > 1.0) c = 1.0;
      if (c < -1.0) c = -1.0;
      phi = atan2(s,c);

      ANDgate = 0;
      if (constraint.dihedral.amin < constraint.dihedral.amax) {
        if (phi > constraint.dihedral.amin && phi < constraint.dihedral.amax) ANDgate = 1;
      } else {
        if (phi > constraint.dihedral.amin || phi < constraint.dihedral.amax) ANDgate = 1;
      }
      if (constraint.dihedral.amin2 < constraint.dihedral.amax2) {
        if (phi > constraint.dihedral.amin2 && phi < constraint.dihedral.amax2) ANDgate = 1;
      } else {
        if (phi > constraint.dihedral.amin2 || phi < constraint.dihedral.amax2) ANDgate = 1;
      }
      if (ANDgate != 1) constraint.satisfied = false;
    } else if (constraint.type == Reaction::Constraint::Type::ARRHENIUS) {
      std::vector<tagint> myglove(glove.begin(), glove.begin() + rxn.reactant->natoms);
      t = get_temperature(myglove);
      prrhob = constraint.arrhenius.A*pow(t,constraint.arrhenius.n)*
        exp(-constraint.arrhenius.E_a/(force->boltz*t));
      if (prrhob < constraint.arrhenius.rrhandom->uniform()) constraint.satisfied = false;
    } else if (constraint.type == Reaction::Constraint::Type::RMSD) {
      // call superpose
      int iatom;
      int iref = -1; // choose first atom as reference
      int n2superpose = 0;
      double **xfrozen; // coordinates for the "frozen" target molecule
      double **xmobile; // coordinates for the "mobile" molecule
      int ifragment = constraint.ids[0];
      if (ifragment >= 0) {
        for (int j = 0; j < rxn.reactant->natoms; j++)
          if (rxn.reactant->fragmentmask[ifragment][j]) n2superpose++;
        memory->create(xfrozen,n2superpose,3,"bond/react:xfrozen");
        memory->create(xmobile,n2superpose,3,"bond/react:xmobile");
        int myincr = 0;
        for (int j = 0; j < rxn.reactant->natoms; j++) {
          if (rxn.reactant->fragmentmask[ifragment][j]) {
            iatom = atom->map(glove[j]);
            if (iref == -1) iref = iatom;
            iatom = domain->closest_image(iref,iatom);
            for (int k = 0; k < 3; k++) {
              xfrozen[myincr][k] = x[iatom][k];
              xmobile[myincr][k] = rxn.reactant->x[j][k];
            }
            myincr++;
          }
        }
      } else {
        int iatom;
        int iref = -1; // choose first atom as reference
        n2superpose = rxn.reactant->natoms;
        memory->create(xfrozen,n2superpose,3,"bond/react:xfrozen");
        memory->create(xmobile,n2superpose,3,"bond/react:xmobile");
        for (int j = 0; j < n2superpose; j++) {
          iatom = atom->map(glove[j]);
          if (iref == -1) iref = iatom;
          iatom = domain->closest_image(iref,iatom);
          for (int k = 0; k < 3; k++) {
            xfrozen[j][k] = x[iatom][k];
            xmobile[j][k] = rxn.reactant->x[j][k];
          }
        }
      }
      Superpose3D<double, double **> superposer(n2superpose);
      double rmsd = superposer.Superpose(xfrozen, xmobile);
      memory->destroy(xfrozen);
      memory->destroy(xmobile);
      if (rmsd > constraint.rmsd.rmsdmax) constraint.satisfied = false;
    } else if (constraint.type == Reaction::Constraint::Type::CUSTOM) {
      constraint.satisfied = custom_constraint(constraint.custom.str, rxn, glove); // NOLINT
    }
  }

  if (!rxn.constraints.empty()) {
    std::string evalstr = rxn.constraintstr;
    for (auto &constraint : rxn.constraints) {
      evalstr.replace(evalstr.find('C'), 1, constraint.satisfied ? "1" : "0");
    }
    std::vector<char> buffer(evalstr.begin(), evalstr.end());
    buffer.push_back('\0');
    double verdict = input->variable->evaluate_boolean(buffer.data());
    if (verdict == 0.0) return 0;
  }

  // let's also check chirality within 'check_constraint'
  for (int i = 0; i < rxn.reactant->natoms; i++) {
    if (rxn.atoms[i].chiral[0] == 1) {
      double my4coords[12];
      // already ensured, by transitive property, that chiral simulation atom has four neighs
      for (int j = 0; j < 4; j++) {
        atom1 = atom->map(glove[i]);
        // loop over known types involved in chiral center
        for (int jj = 0; jj < 4; jj++) {
          if (atom->type[atom->map(xspecial[atom1][j])] == rxn.atoms[i].chiral[jj+2]) {
            atom2 = atom->map(xspecial[atom1][j]);
            atom2 = domain->closest_image(atom1,atom2);
            for (int k = 0; k < 3; k++) {
              my4coords[3*jj+k] = x[atom2][k];
            }
            break;
          }
        }
      }
      if (get_chirality(my4coords) != rxn.atoms[i].chiral[1]) return 0;
    }
  }

  return 1;
}

/* ----------------------------------------------------------------------
return pre-reaction atom or fragment location
fragment: given pre-reacted molID (reactant) and fragID,
          return geometric center (of mapped simulation atoms)
------------------------------------------------------------------------- */

void FixBondReact::get_IDcoords(Reaction::Constraint::IDType idtype, int myID,
                                double *center, Molecule *mol, std::vector<tagint> &glove)
{
  double **x = atom->x;
  if (idtype == Reaction::Constraint::IDType::ATOM) {
    int iatom = atom->map(glove[myID-1]);
    for (int i = 0; i < 3; i++)
      center[i] = x[iatom][i];
  } else {
    int iref = -1; // choose first atom as reference
    int iatom;
    int nfragatoms = 0;
    for (int i = 0; i < 3; i++)
      center[i] = 0;

    for (int i = 0; i < mol->natoms; i++) {
      if (mol->fragmentmask[myID][i]) {
        if (iref == -1) iref = atom->map(glove[i]);
        iatom = atom->map(glove[i]);
        iatom = domain->closest_image(iref,iatom);
        for (int j = 0; j < 3; j++)
          center[j] += x[iatom][j];
        nfragatoms++;
      }
    }
    if (nfragatoms > 0)
      for (int i = 0; i < 3; i++) center[i] /= nfragatoms;
  }
}

/* ----------------------------------------------------------------------
compute local temperature: average over all atoms in reaction template
------------------------------------------------------------------------- */

double FixBondReact::get_temperature(std::vector<tagint> &glove)
{
  int i,ilocal;
  double adof = domain->dimension;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;

  double t = 0.0;

  if (rmass) {
    for (i = 0; i < glove.size(); i++) {
      ilocal = atom->map(glove[i]);
      t += (v[ilocal][0]*v[ilocal][0] + v[ilocal][1]*v[ilocal][1] +
            v[ilocal][2]*v[ilocal][2]) * rmass[ilocal];
    }
  } else {
    for (i = 0; i < glove.size(); i++) {
      ilocal = atom->map(glove[i]);
      t += (v[ilocal][0]*v[ilocal][0] + v[ilocal][1]*v[ilocal][1] +
            v[ilocal][2]*v[ilocal][2]) * mass[type[ilocal]];
    }
  }

  // final temperature
  double dof = adof*glove.size();
  double tfactor = force->mvv2e / (dof * force->boltz);
  t *= tfactor;
  return t;
}

/* ----------------------------------------------------------------------
compute sum of partial charges in rxn site, for updated atoms
------------------------------------------------------------------------- */

double FixBondReact::get_totalcharge(Reaction &rxn, std::vector<tagint> &glove)
{
  int j,jj;
  double *q = atom->q;
  double sim_total_charge = 0.0;
  for (j = 0; j < rxn.reactant->natoms; j++) {
    jj = rxn.atoms[j].amap[1]-1;
    if (rxn.atoms[jj].recharged == 1)
      sim_total_charge += q[atom->map(glove[jj])];
  }
  return sim_total_charge;
}

/* ----------------------------------------------------------------------
get per-atom variable names used by custom constraint
------------------------------------------------------------------------- */

void FixBondReact::customvarnames()
{
  std::size_t pos,pos1,pos2,pos3;
  int prev3;
  std::string varstr,argstr,varid;

  // search all constraints' varstr for special 'rxn' functions
  //   add variable names to customvarstrs
  //   add values to customvars

  for (auto &rxn : rxns) {
    for (auto &constraint : rxn.constraints) {
      if (constraint.type == Reaction::Constraint::Type::CUSTOM) {
        varstr = constraint.custom.str;
        prev3 = -1;
        while (true) {
          // find next reaction special function occurrence
          pos1 = std::string::npos;
          for (int i = 0; i < nrxnfunction; i++) {
            if (peratomflag[i] == 0) continue;
            pos = varstr.find(rxnfunclist[i],prev3+1);
            if (pos == std::string::npos) continue;
            if (pos < pos1) pos1 = pos;
          }
          if (pos1 == std::string::npos) break;

          pos2 = varstr.find("(",pos1);
          pos3 = varstr.find(")",pos2);
          if (pos2 == std::string::npos || pos3 == std::string::npos)
            error->all(FLERR,"Fix bond/react: Illegal rxn function syntax\n");
          prev3 = (int)pos3;
          argstr = varstr.substr(pos2+1,pos3-pos2-1);
          argstr.erase(remove_if(argstr.begin(), argstr.end(), isspace), argstr.end()); // remove whitespace
          pos2 = argstr.find(",");
          if (pos2 != std::string::npos) varid = argstr.substr(0,pos2);
          else varid = argstr;
          // check if we already know about this variable
          int varidflag = 0;
          for (int j = 0; j < ncustomvars; j++) {
            if (customvarstrs[j] == varid) {
              varidflag = 1;
              break;
            }
          }
          if (!varidflag) {
            customvarstrs.resize(ncustomvars+1);
            customvarstrs[ncustomvars++] = varid;
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
evaluate per-atom variables needed for custom constraint
------------------------------------------------------------------------- */

void FixBondReact::get_customvars()
{
  double *tempvvec;
  std::string varid;
  int nall = atom->nlocal + atom->nghost;

  memory->create(tempvvec,nall,"bond/react:tempvvec");
  if (vvec == nullptr) {
    memory->create(vvec,nall,ncustomvars,"bond/react:vvec");
    nvvec = nall;
  }
  if (nvvec < nall) {
    memory->grow(vvec,nall,ncustomvars,"bond/react:vvec");
    nvvec = nall;
  }
  for (int i = 0; i < ncustomvars; i++) {
    varid = customvarstrs[i];
    if (varid.substr(0,2) != "v_") error->all(FLERR,"Fix bond/react: Reaction special function variable "
                                     "name should begin with 'v_'");
    varid = varid.substr(2);
    int ivar = input->variable->find(varid.c_str());
    if (ivar < 0)
      error->all(FLERR,"Fix bond/react: Reaction special function variable "
                                   "name does not exist");
    if (!input->variable->atomstyle(ivar))
      error->all(FLERR,"Fix bond/react: Reaction special function must "
                                   "reference an atom-style variable");

    input->variable->compute_atom(ivar,igroup,tempvvec,1,0);
    for (int j = 0; j < nall; j++) vvec[j][i] = tempvvec[j];
  }
  memory->destroy(tempvvec);
}

/* ----------------------------------------------------------------------
evaulate expression for variable constraint
------------------------------------------------------------------------- */

bool FixBondReact::custom_constraint(const std::string &varstr, Reaction &rxn, std::vector<tagint> &glove)
{
  std::size_t pos,pos1,pos2,pos3;
  int irxnfunc;
  int prev3 = -1;
  std::string argstr,varid,fragid,evlcat;
  std::vector<std::string> evlstr;

  // search varstr for special 'rxn' functions
  while (true) {
    // find next reaction special function occurrence
    pos1 = std::string::npos;
    for (int i = 0; i < nrxnfunction; i++) {
      pos = varstr.find(rxnfunclist[i],prev3+1);
      if (pos == std::string::npos) continue;
      if (pos < pos1) {
        pos1 = pos;
        irxnfunc = i;
      }
    }
    if (pos1 == std::string::npos) break;

    fragid = "all"; // operate over entire reaction site by default
    pos2 = varstr.find("(",pos1);
    pos3 = varstr.find(")",pos2);
    if (pos2 == std::string::npos || pos3 == std::string::npos)
      error->one(FLERR,"Fix bond/react: Illegal rxn function syntax\n");
    evlstr.push_back(varstr.substr(prev3+1,pos1-(prev3+1)));
    prev3 = pos3;
    argstr = varstr.substr(pos2+1,pos3-pos2-1);
    argstr.erase(remove_if(argstr.begin(), argstr.end(), isspace), argstr.end()); // remove whitespace
    pos2 = argstr.find(",");
    if (pos2 != std::string::npos) {
      varid = argstr.substr(0,pos2);
      fragid = argstr.substr(pos2+1);
    } else varid = argstr;
    evlstr.push_back(std::to_string(rxnfunction(rxnfunclist[irxnfunc], varid, fragid, rxn.reactant, glove)));
  }
  evlstr.push_back(varstr.substr(prev3+1));

  for (auto & evl : evlstr) evlcat += evl;
  return static_cast<bool>(input->variable->compute_equal(evlcat));
}

/* ----------------------------------------------------------------------
currently three 'rxn' functions: rxnsum, rxnave, and rxnbond
------------------------------------------------------------------------- */

double FixBondReact::rxnfunction(const std::string& rxnfunc, const std::string& varid,
                                 const std::string& fragid, Molecule *mol, std::vector<tagint> &glove)
{
  int ifrag = -1;
  if (fragid != "all") {
    ifrag = mol->findfragment(fragid.c_str());
    if (ifrag < 0) error->one(FLERR,"Bond/react: Molecule fragment "
                              "in reaction special function does not exist");
  }

  // start with 'rxnbond' per-bond function
  // for 'rxnbond', varid corresponds to 'compute bond/local' name,
  //                and fragid is a pre-reaction fragment containing the two atoms in the bond
  if (rxnfunc == "rxnbond") {
    int icompute,ibond;
    double perbondval;
    std::set<tagint> aset;
    std::string computeid = varid;

    if (computeid.substr(0,2) != "c_") error->one(FLERR,"Bond/react: Reaction special function compute "
                                         "name should begin with 'c_'");
    computeid = computeid.substr(2);
    icompute = modify->find_compute(computeid);
    if (icompute < 0) error->one(FLERR,"Bond/react: Reaction special function compute name does not exist");
    cperbond = modify->compute[icompute];
    std::string compute_style = cperbond->style;
    if (compute_style != "bond/local") error->one(FLERR,"Bond/react: Compute used by reaction "
                                         "special function 'rxnbond' must be of style 'bond/local'");
    if (cperbond->size_local_cols > 0) error->one(FLERR,"Bond/react: 'Compute bond/local' used by reaction "
                                         "special function 'rxnbond' must compute one value");

    if (atoms2bondflag == 0) {
      atoms2bondflag = 1;
      get_atoms2bond(cperbond->groupbit);
    }

    for (int i = 0; i < mol->natoms; i++) {
      if (mol->fragmentmask[ifrag][i]) {
        aset.insert(glove[i]);
      }
    }
    if (aset.size() != 2) error->one(FLERR,"Bond/react: Molecule fragment of reaction special function 'rxnbond' "
                     "must contain exactly two atoms");

    if (cperbond->invoked_local != lmp->update->ntimestep)
      cperbond->compute_local();

    auto it = atoms2bond.find(aset);
    if (it == atoms2bond.end()) error->one(FLERR,"Bond/react: Unable to locate bond referenced by "
                                            "reaction special function 'rxnbond'");
    ibond = it->second;
    perbondval = cperbond->vector_local[ibond];
    return perbondval;
  }

  int ivar = -1;
  for (int i = 0; i < ncustomvars; i++) {
    if (varid == customvarstrs[i]) {
      ivar = i;
      break;
    }
  }
  // variable name should always be found, at this point
  // however, let's double check for completeness
  if (ivar < 0)
    error->one(FLERR,"Fix bond/react: Reaction special function variable "
                                 "name does not exist");

  int iatom;
  int nsum = 0;
  double sumvvec = 0;
  if (rxnfunc == "rxnsum" || rxnfunc == "rxnave") {
    if (fragid == "all") {
      for (int i = 0; i < mol->natoms; i++) {
        iatom = atom->map(glove[i]);
        sumvvec += vvec[iatom][ivar];
      }
      nsum = mol->natoms;
    } else {
      for (int i = 0; i < mol->natoms; i++) {
        if (mol->fragmentmask[ifrag][i]) {
          iatom = atom->map(glove[i]);
          sumvvec += vvec[iatom][ivar];
          nsum++;
        }
      }
    }
  }

  if (rxnfunc == "rxnsum") return sumvvec;
  if (rxnfunc == "rxnave") return sumvvec/nsum;
  return 0.0;
}

/* ----------------------------------------------------------------------
populate map to get bond index from atom IDs
------------------------------------------------------------------------- */

void FixBondReact::get_atoms2bond(int cgroupbit)
{
  int i,m,atom1,atom2,btype,nb;
  std::set<tagint> aset;

  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *mask = atom->mask;

  m = 0;
  atoms2bond.clear();
  for (atom1 = 0; atom1 < nlocal; atom1++) {
    if (!(mask[atom1] & cgroupbit)) continue;
    nb = num_bond[atom1];
    for (i = 0; i < nb; i++) {
      btype = bond_type[atom1][i];
      atom2 = atom->map(bond_atom[atom1][i]);
      if (atom2 < 0 || !(mask[atom2] & cgroupbit)) continue;
      if (newton_bond == 0 && tag[atom1] > tag[atom2]) continue;
      if (btype == 0) continue;
      aset = {tag[atom1], tag[atom2]};
      atoms2bond.insert(std::make_pair(aset,m++));
    }
  }
}

/* ----------------------------------------------------------------------
return handedness (1 or -1) of a chiral center, given ordered set of coordinates
------------------------------------------------------------------------- */

int FixBondReact::get_chirality(double four_coords[12])
{
  // define oriented plane with first three coordinates
  double vec1[3],vec2[3],vec3[3],vec4[3],mean3[3],dot;

  for (int i = 0; i < 3; i++) {
    vec1[i] = four_coords[i]-four_coords[i+3];
    vec2[i] = four_coords[i+3]-four_coords[i+6];
  }

  MathExtra::cross3(vec1,vec2,vec3);

  for (int i = 0; i < 3; i++) {
    mean3[i] = (four_coords[i] + four_coords[i+3] +
                four_coords[i+6])/3;
    vec4[i] = four_coords[i+9] - mean3[i];
  }

  dot = MathExtra::dot3(vec3,vec4);
  dot = dot/fabs(dot);
  return (int) dot;
}

/* ----------------------------------------------------------------------
  Determine which pre-reacted template atoms are at least three bonds
  away from edge atoms.
------------------------------------------------------------------------- */

void FixBondReact::find_landlocked_atoms(Reaction &rxn)
{
  // landlocked_atoms are atoms for which all topology is contained in reacted template
  // if dihedrals/impropers exist: this means that edge atoms are not in their 1-3 neighbor list
  //   note: due to various usage/definitions of impropers, treated same as dihedrals
  // if angles exist: this means edge atoms not in their 1-2 neighbors list
  // if just bonds: this just means that edge atoms are not landlocked
  // Note: landlocked defined in terms of reacted template
  // if no edge atoms (small reacting molecule), all atoms are landlocked
  // we can delete all current topology of landlocked atoms and replace

  // always remove edge atoms from landlocked list
  for (int i = 0; i < rxn.product->natoms; i++) {
    if (rxn.atoms[i].created == 0 && rxn.atoms[rxn.atoms[i].amap[1]-1].edge == 1)
      rxn.atoms[i].landlocked = 0;
    else rxn.atoms[i].landlocked = 1;
  }
  int nspecial_limit = -1;
  if (force->angle && rxn.product->angleflag) nspecial_limit = 0;

  if ((force->dihedral && rxn.product->dihedralflag) ||
      (force->improper && rxn.product->improperflag)) nspecial_limit = 1;

  if (nspecial_limit != -1) {
    for (int i = 0; i < rxn.product->natoms; i++) {
      for (int j = 0; j < rxn.product->nspecial[i][nspecial_limit]; j++) {
        for (int k = 0; k < rxn.reactant->natoms; k++) {
          if (rxn.atoms[rxn.product->special[i][j]-1].amap[1] == k+1 && rxn.atoms[k].edge == 1) {
            rxn.atoms[i].landlocked = 0;
          }
        }
      }
    }
  }

  // bad molecule templates check
  // if atoms change types, but aren't landlocked, that's bad
  for (int i = 0; i < rxn.product->natoms; i++) {
    if ((rxn.atoms[i].created == 0) &&
        (rxn.product->type[i] != rxn.reactant->type[rxn.atoms[i].amap[1]-1]) &&
        (rxn.atoms[i].landlocked == 0))
      error->all(FLERR, "Fix bond/react: Atom type affected by reaction {} is too close "
                 "to template edge", rxn.name);
  }

  // additionally, if a bond changes type, but neither involved atom is landlocked, bad
  // would someone want to change an angle type but not bond or atom types? (etc.) ...hopefully not yet
  for (int i = 0; i < rxn.product->natoms; i++) {
    if (rxn.atoms[i].created == 0) {
      if (rxn.atoms[i].landlocked == 0) {
        for (int j = 0; j < rxn.product->num_bond[i]; j++) {
          int product_atomj = rxn.product->bond_atom[i][j];
          if (rxn.atoms[product_atomj-1].landlocked == 0) {
            int onemol_atomi = rxn.atoms[i].amap[1];
            int onemol_batom;
            for (int m = 0; m < rxn.reactant->num_bond[onemol_atomi-1]; m++) {
              onemol_batom = rxn.reactant->bond_atom[onemol_atomi-1][m];
              if ((onemol_batom == rxn.atoms[product_atomj-1].amap[1]) &&
                  (rxn.product->bond_type[i][j] != rxn.reactant->bond_type[onemol_atomi-1][m]))
                error->all(FLERR, "Fix bond/react: Bond type affected by reaction {} is "
                           "too close to template edge",rxn.name);
            }
            if (newton_bond) {
              int onemol_atomj = rxn.atoms[product_atomj-1].amap[1];
              for (int m = 0; m < rxn.reactant->num_bond[onemol_atomj-1]; m++) {
                onemol_batom = rxn.reactant->bond_atom[onemol_atomj-1][m];
                if ((onemol_batom == rxn.atoms[i].amap[1]) &&
                    (rxn.product->bond_type[i][j] != rxn.reactant->bond_type[onemol_atomj-1][m]))
                  error->all(FLERR, "Fix bond/react: Bond type affected by reaction {} is "
                             "too close to template edge",rxn.name);
              }
            }
          }
        }
      }
    }
  }

  // additionally, if a deleted atom is bonded to an atom that is not deleted, bad
  for (int i = 0; i < rxn.reactant->natoms; i++) {
    if (rxn.atoms[i].deleted == 1) {
      int ii = rxn.atoms[i].ramap[1] - 1;
      for (int j = 0; j < rxn.product->nspecial[ii][0]; j++) {
        if (rxn.atoms[rxn.atoms[rxn.product->special[ii][j]-1].amap[1]-1].deleted == 0) {
          error->all(FLERR,"Fix bond/react: A deleted atom cannot remain bonded to an atom that is not deleted");
        }
      }
    }
  }

  // also, if atoms change number of bonds, but aren't landlocked, that could be bad
  int warnflag = 0;
  if (comm->me == 0)
    for (int i = 0; i < rxn.product->natoms; i++) {
      if ((rxn.atoms[i].created == 0) &&
          (rxn.product->nspecial[i][0] != rxn.reactant->nspecial[rxn.atoms[i].amap[1]-1][0]) &&
          (rxn.atoms[i].landlocked == 0)) {
        warnflag = 1;
        break;
      }
    }

  // also, if an atom changes any of its bonds, but is not landlocked, that could be bad
  int thereflag;
  if (comm->me == 0)
    for (int i = 0; i < rxn.product->natoms; i++) {
      if (rxn.atoms[i].landlocked == 1) continue;
      for (int j = 0; j < rxn.product->nspecial[i][0]; j++) {
        int oneneighID = rxn.atoms[rxn.product->special[i][j]-1].amap[1];
        int ii = rxn.atoms[i].amap[1] - 1;
        thereflag = 0;
        for (int k = 0; k < rxn.reactant->nspecial[ii][0]; k++) {
          if (oneneighID == rxn.reactant->special[ii][k]) {
            thereflag = 1;
            break;
          }
        }
        if (thereflag == 0) {
          warnflag = 1;
          break;
        }
      }
      if (warnflag == 1) break;
    }

  if (comm->me == 0 && warnflag == 1) error->warning(FLERR, "Fix bond/react: Atom affected "
                       "by reaction {} is too close to template edge",rxn.name);

  // finally, if a created atom is not landlocked, bad!
  for (int i = 0; i < rxn.product->natoms; i++) {
    if (rxn.atoms[i].created == 1 && rxn.atoms[i].landlocked == 0) {
      error->one(FLERR,"Fix bond/react: Created atom too close to template edge");
    }
  }
}

/* ----------------------------------------------------------------------
let's dedup global_mega_glove
allows for same site undergoing different pathways, in parallel
------------------------------------------------------------------------- */

void FixBondReact::dedup_mega_gloves(Dedup_Modes dedup_mode)
{
  // dedup_mode == LOCAL for local_dedup
  // dedup_mode == GLOBAL for global_mega_glove

  if (dedup_mode == Dedup_Modes::GLOBAL)
    for (auto &rxn : rxns)
      rxn.ghostly_rxn_count = 0;

  int dedup_size = 0;
  if (dedup_mode == Dedup_Modes::LOCAL) {
    dedup_size = my_num_mega;
  } else if (dedup_mode == Dedup_Modes::GLOBAL) {
    dedup_size = global_megasize;
  }

  double **dedup_glove;
  memory->create(dedup_glove,max_natoms+cuff,dedup_size,"bond/react:dedup_glove");

  if (dedup_mode == Dedup_Modes::LOCAL) {
    for (int i = 0; i < dedup_size; i++) {
      for (int j = 0; j < max_natoms+cuff; j++) {
        dedup_glove[j][i] = my_mega_glove[j][i];
      }
    }
  } else if (dedup_mode == Dedup_Modes::GLOBAL) {
    for (int i = 0; i < dedup_size; i++) {
      for (int j = 0; j < max_natoms+cuff; j++) {
        dedup_glove[j][i] = global_mega_glove[j][i];
      }
    }
  }

  // dedup_mask is size dedup_size and filters reactions that have been deleted
  // a value of 1 means this reaction instance has been deleted
  int *dedup_mask = new int[dedup_size];
  for (int i = 0; i < dedup_size; i++) {
    dedup_mask[i] = 0;
  }

  // let's randomly mix up our reaction instances first
  // then we can feel okay about ignoring ones we've already deleted (or accepted)
  // based off std::shuffle
  auto *temp_rxn = new double[max_natoms+cuff];
  for (int i = dedup_size-1; i > 0; --i) { //dedup_size
    // choose random entry to swap current one with
    int k = floor(random[0]->uniform()*(i+1)); // NOLINT

    // swap entries
    for (int j = 0; j < max_natoms+cuff; j++)
      temp_rxn[j] = dedup_glove[j][i];

    for (int j = 0; j < max_natoms+cuff; j++) {
      dedup_glove[j][i] = dedup_glove[j][k];
      dedup_glove[j][k] = temp_rxn[j];
    }
  }
  delete[] temp_rxn;

  for (int i = 0; i < dedup_size; i++) {
    if (dedup_mask[i] == 0) {
      int myrxnid1 = dedup_glove[0][i]; // NOLINT
      for (int j = 0; j < rxns[myrxnid1].reactant->natoms; j++) {
        int check1 = dedup_glove[j+cuff][i]; // NOLINT
        for (int ii = i + 1; ii < dedup_size; ii++) {
          if (dedup_mask[ii] == 0) {
            int myrxnid2 = dedup_glove[0][ii]; // NOLINT
            for (int jj = 0; jj < rxns[myrxnid2].reactant->natoms; jj++) {
              int check2 = dedup_glove[jj+cuff][ii]; // NOLINT
              if (check2 == check1) {
                dedup_mask[ii] = 1;
                break;
              }
            }
          }
        }
      }
    }
  }

  // we must update local_mega_glove and local_megasize
  // we can simply overwrite local_mega_glove column by column
  if (dedup_mode == Dedup_Modes::LOCAL) {
    int my_new_megasize = 0;
    for (int i = 0; i < my_num_mega; i++) {
      if (dedup_mask[i] == 0) {
        for (int j = 0; j < max_natoms+cuff; j++) {
          my_mega_glove[j][my_new_megasize] = dedup_glove[j][i];
        }
        my_new_megasize++;
      }
    }
    my_num_mega = my_new_megasize;
  }

  // we must update global_mega_glove and global_megasize
  // we can simply overwrite global_mega_glove column by column
  if (dedup_mode == Dedup_Modes::GLOBAL) {
    int new_global_megasize = 0;
    for (int i = 0; i < global_megasize; i++) {
      if (dedup_mask[i] == 0) {
        rxns[(int) dedup_glove[0][i]].ghostly_rxn_count++;
        for (int j = 0; j < max_natoms + cuff; j++) {
          global_mega_glove[j][new_global_megasize] = dedup_glove[j][i];
        }
        new_global_megasize++;
      }
    }
    global_megasize = new_global_megasize;
  }

  memory->destroy(dedup_glove);
  delete[] dedup_mask;
}

/* ----------------------------------------------------------------------
let's unlimit movement of newly bonded atoms after n timesteps.
we give them back to the system thermostat
------------------------------------------------------------------------- */

void FixBondReact::unlimit_bond()
{
  // let's now unlimit in terms of i_limit_tags
  // we just run through all nlocal, looking for > limit_duration
  // then we return i_limit_tag to 0 (which removes from dynamic group)
  int flag, cols;
  int index1 = atom->find_custom("limit_tags",flag,cols);
  int *i_limit_tags = atom->ivector[index1];

  int *i_statted_tags;
  if (stabilization_flag == 1) {
    int index2 = atom->find_custom(statted_id.c_str(),flag,cols);
    i_statted_tags = atom->ivector[index2];
  }

  int index3 = atom->find_custom("react_tags",flag,cols);
  int *i_react_tags = atom->ivector[index3];

  int unlimitflag = 0;
  for (int i = 0; i < atom->nlocal; i++) {
    // unlimit atoms for next step! this resolves # of procs disparity, mostly
    // first '1': indexing offset, second '1': for next step
    if (i_limit_tags[i] != 0 && (update->ntimestep + 1 - i_limit_tags[i]) > rxns[i_react_tags[i]].limit_duration) {
      unlimitflag = 1;
      i_limit_tags[i] = 0;
      if (stabilization_flag == 1) i_statted_tags[i] = 1;
      i_react_tags[i] = 0;
    }
  }

  // really should only communicate this per-atom property, not entire reneighboring
  MPI_Allreduce(MPI_IN_PLACE,&unlimitflag,1,MPI_INT,MPI_MAX,world);
  if (unlimitflag) next_reneighbor = update->ntimestep;
}

/* ----------------------------------------------------------------------
check mega_glove for ghosts
if so, flag for broadcasting for perusal by all processors
------------------------------------------------------------------------- */

void FixBondReact::glove_ghostcheck()
{
  // here we add glove to either local_mega_glove or ghostly_mega_glove
  // ghostly_mega_glove includes atoms that are ghosts, either of this proc or another
  // 'ghosts of another' indication taken from comm->sendlist
  // also includes local gloves that overlap with ghostly gloves, to get dedup right

  for (auto &rxn : rxns) rxn.local_rxn_count = 0;

  for (int i = 0; i < my_num_mega; i++) {
    int rxnID = (int) my_mega_glove[0][i];
    auto &rxn = rxns[rxnID];
    int ghostly = 0;
  #if !defined(MPI_STUBS)
    if (comm->style == Comm::BRICK) {
      if (rxn.create_atoms_flag == 1) {
        ghostly = 1;
      } else {
        for (int j = 0; j < rxn.reactant->natoms; j++) {
          int ilocal = atom->map((tagint) my_mega_glove[j+cuff][i]);
          if (ilocal >= atom->nlocal || localsendlist[ilocal] == 1) {
            ghostly = 1;
            break;
          }
        }
      }
    } else {
      ghostly = 1;
    }
  #endif

    if (ghostly == 1) {
      for (int j = 0; j < rxn.reactant->natoms+cuff; j++) {
        ghostly_mega_glove[j][ghostly_num_mega] = my_mega_glove[j][i];
      }
      ghostly_num_mega++;
    } else {
      rxn.local_rxn_count++;
      for (int j = 0; j < rxn.reactant->natoms+cuff; j++) {
        local_mega_glove[j][local_num_mega] = my_mega_glove[j][i];
      }
      local_num_mega++;
    }
  }
}

/* ----------------------------------------------------------------------
broadcast entries of mega_glove which contain nonlocal atoms for perusal by all processors
------------------------------------------------------------------------- */

void FixBondReact::ghost_glovecast()
{
#if !defined(MPI_STUBS)
  const int nprocs = comm->nprocs;

  global_megasize = 0;

  int *allncols = new int[nprocs];
  for (int i = 0; i < nprocs; i++)
    allncols[i] = 0;
  MPI_Allgather(&ghostly_num_mega, 1, MPI_INT, allncols, 1, MPI_INT, world);
  for (int i = 0; i < nprocs; i++)
    global_megasize = global_megasize + allncols[i];

  if (global_megasize == 0) {
    delete[] allncols;
    return;
  }

  int *allstarts = new int[nprocs];

  int start = 0;
  for (int i = 0; i < comm->me; i++) {
    start += allncols[i];
  }
  MPI_Allgather(&start, 1, MPI_INT, allstarts, 1, MPI_INT, world);
  MPI_Datatype columnunsized, column;
  int sizes[2]    = {max_natoms+cuff, global_megasize};
  int subsizes[2] = {max_natoms+cuff, 1};
  int starts[2]   = {0,0};
  MPI_Type_create_subarray (2, sizes, subsizes, starts, MPI_ORDER_C,
                            MPI_DOUBLE, &columnunsized);
  MPI_Type_create_resized (columnunsized, 0, sizeof(double), &column);
  MPI_Type_commit(&column);

  memory->destroy(global_mega_glove);
  memory->create(global_mega_glove,max_natoms+cuff,global_megasize,"bond/react:global_mega_glove");

  for (int i = 0; i < max_natoms+cuff; i++)
    for (int j = 0; j < global_megasize; j++)
      global_mega_glove[i][j] = 0;

  if (ghostly_num_mega > 0) {
    for (int i = 0; i < max_natoms+cuff; i++) {
      for (int j = 0; j < ghostly_num_mega; j++) {
        global_mega_glove[i][j+start] = ghostly_mega_glove[i][j];
      }
    }
  }
  // let's send to root, dedup, then broadcast
  if (comm->me == 0) {
    MPI_Gatherv(MPI_IN_PLACE, ghostly_num_mega, column, // Note: some values ignored for MPI_IN_PLACE
                &(global_mega_glove[0][0]), allncols, allstarts,
                column, 0, world);
  } else {
    MPI_Gatherv(&(global_mega_glove[0][start]), ghostly_num_mega, column,
                &(global_mega_glove[0][0]), allncols, allstarts,
                column, 0, world);
  }

  if (comm->me == 0) dedup_mega_gloves(Dedup_Modes::GLOBAL); // global_mega_glove mode
  MPI_Bcast(&global_megasize,1,MPI_INT,0,world);
  MPI_Bcast(&(global_mega_glove[0][0]), global_megasize, column, 0, world);

  delete[] allstarts;
  delete[] allncols;

  MPI_Type_free(&column);
  MPI_Type_free(&columnunsized);
#endif
}

/* ----------------------------------------------------------------------
update molecule IDs, charges, types, special lists and all topology
------------------------------------------------------------------------- */

void FixBondReact::update_everything()
{
  int nlocal = atom->nlocal; // must be redefined after create atoms
  int *type = atom->type;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  tagint *tag = atom->tag;
  AtomVec *avec = atom->avec;

  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;

  // used when deleting atoms
  int ndel,ndelone;
  int *mark;
  int nmark = nlocal;
  memory->create(mark,nmark,"bond/react:mark");
  for (int i = 0; i < nmark; i++) mark[i] = 0;

  // used when creating atoms
  addatomtag = 0;
  for (int i = 0; i < nlocal; i++) addatomtag = MAX(addatomtag,tag[i]);
  MPI_Allreduce(MPI_IN_PLACE,&addatomtag,1,MPI_LMP_TAGINT,MPI_MAX,world);
  addatoms.clear();

  // flag used to delete special interactions
  int *delflag;
  memory->create(delflag,atom->maxspecial,"bond/react:delflag");

  // used when creating atoms
  int inserted_atoms_flag = 0;

  // update atom->nbonds, etc.
  // TODO: correctly tally with 'newton off'
  int delta_bonds = 0;
  int delta_angle = 0;
  int delta_dihed = 0;
  int delta_imprp = 0;

  // use the following per-atom arrays to keep track of reacting atoms

  int flag,cols;
  int index1 = atom->find_custom("limit_tags",flag,cols);
  int *i_limit_tags = atom->ivector[index1];

  int *i_statted_tags;
  if (stabilization_flag == 1) {
    int index2 = atom->find_custom(statted_id.c_str(),flag,cols);
    i_statted_tags = atom->ivector[index2];
  }

  int index3 = atom->find_custom("react_tags",flag,cols);
  int *i_react_tags = atom->ivector[index3];

  // pass through twice
  // redefining 'update_num_mega' and 'update_mega_glove' each time
  //  first pass: when glove is all local atoms
  //  second pass: search for local atoms in global_mega_glove
  // add check for local atoms as well

  int update_num_mega;
  tagint **update_mega_glove;
  // for now, keeping rxnID in update_mega_glove, but not rest of cuff in update_mega_glove
  int maxmega = MAX(local_num_mega,global_megasize);
  memory->create(update_mega_glove,max_natoms+1,maxmega,"bond/react:update_mega_glove");

  double *sim_total_charges;
  if (rescale_charges_anyflag) memory->create(sim_total_charges,maxmega,"bond/react:sim_total_charges");

  for (int pass = 0; pass < 2; pass++) {
    update_num_mega = 0;
    int *noccur = new int[rxns.size()];
    for (int i = 0; i < rxns.size(); i++) noccur[i] = 0;
    if (pass == 0) {
      for (int i = 0; i < local_num_mega; i++) {
        auto &rxn = rxns[(int) local_mega_glove[0][i]];
        // reactions already shuffled from dedup procedure, so can skip first N
        // wait, this check needs to be after add atoms, because they can also be 'skipped' due to overlap
        if (noccur[rxn.ID] >= rxn.nlocalkeep) continue;

        // this will be overwritten if reaction skipped by create_atoms below
        update_mega_glove[0][update_num_mega] = (tagint) local_mega_glove[0][i];
        for (int j = 0; j < max_natoms; j++)
          update_mega_glove[j+1][update_num_mega] = (tagint) local_mega_glove[j+cuff][i];

        // atoms inserted here for serial MPI_STUBS build only
        if (rxn.create_atoms_flag == 1) {
          if (insert_atoms_setup(update_mega_glove,update_num_mega,rxn)) inserted_atoms_flag = 1;
          else continue;
        }
        noccur[rxn.ID]++;

        if (rxn.rescale_charges_flag) sim_total_charges[update_num_mega] = local_mega_glove[1][i];
        update_num_mega++;
      }
      MPI_Allreduce(MPI_IN_PLACE, &noccur[0], rxns.size(), MPI_INT, MPI_SUM, world);
      for (auto &rxn : rxns) rxn.reaction_count_total += noccur[rxn.ID];
    } else if (pass == 1) {
      for (int i = 0; i < global_megasize; i++) {
        auto &rxn = rxns[(int) global_mega_glove[0][i]];
        // reactions already shuffled from dedup procedure, so can skip first N
        if (noccur[rxn.ID] >= rxn.nghostlykeep) continue;

        // this will be overwritten if reaction skipped by create_atoms below
        update_mega_glove[0][update_num_mega] = (tagint) global_mega_glove[0][i];
        for (int j = 0; j < max_natoms; j++)
          update_mega_glove[j+1][update_num_mega] = (tagint) global_mega_glove[j+cuff][i];

        // we can insert atoms here, now that reactions are finalized
        // can't do it any earlier, due to skipped reactions (max_rxn)
        // for MPI build, reactions that create atoms are always treated as 'global'

        if (rxn.create_atoms_flag == 1) {
          if (insert_atoms_setup(update_mega_glove,update_num_mega,rxn)) inserted_atoms_flag = 1;
          else continue;
        }
        noccur[rxn.ID]++;
        rxn.reaction_count_total++;

        if (rxn.rescale_charges_flag) sim_total_charges[update_num_mega] = global_mega_glove[1][i];
        update_num_mega++;
      }
    }
    delete[] noccur;

    // find current max molecule ID and shift for each proc
    tagint moloffset = 0;
    if (molid_mode == Reset_Mol_IDs::MOLMAP) {
      tagint maxmol_all = 0;
      for (int i = 0; i < atom->nlocal; i++) maxmol_all = MAX(maxmol_all, atom->molecule[i]);
      MPI_Allreduce(MPI_IN_PLACE, &maxmol_all, 1, MPI_LMP_TAGINT, MPI_MAX, world);
      // find number of new molids needed for each proc
      if (pass == 0) {
        tagint molcreate = 0;
        for (int i = 0; i < update_num_mega; i++) {
          auto &rxn = rxns[(int) update_mega_glove[0][i]];
          molcreate += rxn.nnewmolids;
        }
        MPI_Scan(&molcreate, &moloffset, 1, MPI_LMP_TAGINT, MPI_SUM, world);
        moloffset = moloffset - molcreate + maxmol_all;
      }
      if (pass == 1) moloffset = maxmol_all;
    }

    if (update_num_mega == 0) continue;

    // for 'reset_mol_ids molmap', update molecule IDs
    // assumes consistent molecule IDs in pre- and post-reaction template
    // NOTE: all procs assumed to have same update_mega_glove for second pass
    // NOTE: must be done before add atoms, because add_atoms deletes ghost info
    if (molid_mode == Reset_Mol_IDs::MOLMAP) {
      for (int i = 0; i < update_num_mega; i++) {
        auto &rxn = rxns[(int) update_mega_glove[0][i]];
        if (!rxn.reactant->moleculeflag || !rxn.product->moleculeflag) continue;
        tagint molmapid = -1;
        for (int j = 0; j < rxn.product->natoms; j++) {
          int neednewid = 0;
          tagint *thismolid;
          if (rxn.atoms[j].created == 1) {
            for (auto & myaddatom : addatoms) {
              if (myaddatom.tag == update_mega_glove[j+1][i]) {
                thismolid = &(myaddatom.molecule);
                neednewid = 1;
                break;
              }
            }
          } else {
            int jj = rxn.atoms[j].amap[1]-1;
            int jlocal = atom->map(update_mega_glove[jj+1][i]);
            if (jlocal < nlocal && jlocal >= 0) {
              thismolid = &(atom->molecule[jlocal]);
              neednewid = 1;
            }
          }
          if (neednewid == 1) {
            if (rxn.atoms[j].newmolid != 0) {
              molmapid = moloffset + rxn.atoms[j].newmolid;
            } else {
              for (int k = 0; k < rxn.reactant->natoms; k++) {
                if (rxn.product->molecule[j] == rxn.reactant->molecule[k]) {
                  int klocal = atom->map(update_mega_glove[k+1][i]);
                  if (klocal >= 0) {
                    molmapid = atom->molecule[klocal];
                    break;
                  }
                }
              }
            }
            if (molmapid != -1) {
              *thismolid = molmapid;
            } else {
              error->one(FLERR,"Fix bond/react: unable to assign molecule ID using 'molmap' option. "
                                  "Need ghost atoms from further away");
            }
            thismolid = nullptr;
          }
        }
        moloffset += rxn.nnewmolids;
      }
    }

    // insert all atoms for all rxns here
    if (inserted_atoms_flag == 1) {
      // clear to-be-overwritten ghost info
      atom->nghost = 0;
      atom->avec->clear_bonus();

      for (auto & myaddatom : addatoms) {
        atom->avec->create_atom(myaddatom.type,myaddatom.x);
        int n = atom->nlocal - 1;
        atom->tag[n] = myaddatom.tag;
        atom->molecule[n] = myaddatom.molecule;
        atom->mask[n] = myaddatom.mask;
        atom->image[n] = myaddatom.image;
        atom->v[n][0] = myaddatom.v[0];
        atom->v[n][1] = myaddatom.v[1];
        atom->v[n][2] = myaddatom.v[2];
        if (atom->rmass) atom->rmass[n]= myaddatom.rmass;
        modify->create_attribute(n);
      }

      // reset atom->map
      if (atom->map_style != Atom::MAP_NONE) {
        atom->map_init();
        atom->map_set();
      }
    }

    // mark to-delete atoms
    nlocal = atom->nlocal;
    if (nlocal > nmark) {
      memory->grow(mark,nlocal,"bond/react:mark");
      for (int i = nmark; i < nlocal; i++) mark[i] = 0;
      nmark = nlocal;
    }
    for (int i = 0; i < update_num_mega; i++) {
      auto &rxn = rxns[(int) update_mega_glove[0][i]];
      for (int j = 0; j < rxn.reactant->natoms; j++) {
        int iatom = atom->map(update_mega_glove[j+1][i]);
        if (rxn.atoms[j].deleted == 1 && iatom >= 0 && iatom < nlocal) {
          mark[iatom] = 1;
        }
      }
    }

    // update charges and types of landlocked atoms
    // also keep track of 'stabilization' groups here
    int n_custom_charge;
    double charge_rescale_addend;
    for (int i = 0; i < update_num_mega; i++) {
      charge_rescale_addend = 0;
      auto &rxn = rxns[(int) update_mega_glove[0][i]];
      if (rxn.rescale_charges_flag) {
        n_custom_charge = rxn.rescale_charges_flag;
        charge_rescale_addend = (sim_total_charges[i]-rxn.mol_total_charge)/n_custom_charge;
      }
      for (int j = 0; j < rxn.product->natoms; j++) {
        int jj = rxn.atoms[j].amap[1]-1;
        int ilocal = atom->map(update_mega_glove[jj+1][i]);
        if (ilocal >= 0 && ilocal < nlocal) {

          // update->ntimestep could be 0. so add 1 throughout
          i_limit_tags[ilocal] = update->ntimestep + 1;
          if (stabilization_flag == 1) i_statted_tags[ilocal] = 0;
          i_react_tags[ilocal] = rxn.ID;

          if (rxn.atoms[j].landlocked == 1)
            type[ilocal] = rxn.product->type[j];
          if (rxn.product->qflag && atom->q_flag && rxn.atoms[jj].recharged == 1) {
            double *q = atom->q;
            q[ilocal] = rxn.product->q[j]+charge_rescale_addend;
          }
        }
      }
    }

    int insert_num;
    // very nice and easy to completely overwrite special bond info for landlocked atoms
    for (int i = 0; i < update_num_mega; i++) {
      auto &rxn = rxns[(int) update_mega_glove[0][i]];
      for (int j = 0; j < rxn.product->natoms; j++) {
        int jj = rxn.atoms[j].amap[1]-1;
        int ilocal = atom->map(update_mega_glove[jj+1][i]);
        if (ilocal < nlocal && ilocal >= 0) {
          if (rxn.atoms[j].landlocked == 1) {
            for (int k = 0; k < 3; k++) {
              nspecial[ilocal][k] = rxn.product->nspecial[j][k];
            }
            for (int p = 0; p < rxn.product->nspecial[j][2]; p++) {
              special[ilocal][p] = update_mega_glove[rxn.atoms[rxn.product->special[j][p]-1].amap[1]][i];
            }
          }
          // now delete and replace landlocked atoms from non-landlocked atoms' special info
          // delete 1-2, 1-3, 1-4 specials individually. only delete if special exists in pre-reaction template
          if (rxn.atoms[j].landlocked == 0) {
            int ispec, fspec, imolspec, fmolspec, nspecdel[3];
            for (int k = 0; k < 3; k++) nspecdel[k] = 0;
            for (int k = 0; k < atom->maxspecial; k++) delflag[k] = 0;
            for (int specn = 0; specn < 3; specn++) {
              if (specn == 0) {
                imolspec = 0;
                ispec = 0;
              } else {
                imolspec = rxn.reactant->nspecial[jj][specn-1];
                ispec = nspecial[ilocal][specn-1];
              }
              fmolspec = rxn.reactant->nspecial[jj][specn];
              fspec = nspecial[ilocal][specn];
              for (int k = ispec; k < fspec; k++) {
                for (int p = imolspec; p < fmolspec; p++) {
                  if (update_mega_glove[rxn.reactant->special[jj][p]][i] == special[ilocal][k]) {
                    delflag[k] = 1;
                    for (int m = 2; m >= specn; m--) nspecdel[m]++;
                    break;
                  }
                }
              }
            }
            int incr = 0;
            for (int k = 0; k < nspecial[ilocal][2]; k++)
              if (delflag[k] == 0) special[ilocal][incr++] = special[ilocal][k];
            for (int m = 0; m < 3; m++) nspecial[ilocal][m] -= nspecdel[m];
            // now reassign from reacted template
            for (int k = 0; k < rxn.product->nspecial[j][2]; k++) {
              if (k > rxn.product->nspecial[j][1] - 1) {
                insert_num = nspecial[ilocal][2]++;
              } else if (k > rxn.product->nspecial[j][0] - 1) {
                insert_num = nspecial[ilocal][1]++;
                nspecial[ilocal][2]++;
              } else {
                insert_num = nspecial[ilocal][0]++;
                nspecial[ilocal][1]++;
                nspecial[ilocal][2]++;
              }
              if (nspecial[ilocal][2] > atom->maxspecial)
                error->one(FLERR,"Fix bond/react special bond generation overflow");
              for (int n = nspecial[ilocal][2]-1; n > insert_num; n--) {
                special[ilocal][n] = special[ilocal][n-1];
              }
              special[ilocal][insert_num] = update_mega_glove[rxn.atoms[rxn.product->special[j][k]-1].amap[1]][i];
            }
          }
        }
      }
    }

    // next let's update bond info
    // cool thing is, newton_bond issues are already taken care of in templates
    // same with class2 improper issues, which is why this fix started in the first place
    // also need to find any instances of bond history to update histories
    auto histories = modify->get_fix_by_style("BOND_HISTORY");
    int n_histories = histories.size();

    for (int i = 0; i < update_num_mega; i++) {
      auto &rxn = rxns[(int) update_mega_glove[0][i]];
      // let's first delete all bond info about landlocked atoms
      for (int j = 0; j < rxn.product->natoms; j++) {
        int jj = rxn.atoms[j].amap[1]-1;
        if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
          if (rxn.atoms[j].landlocked == 1) {
            delta_bonds -= num_bond[atom->map(update_mega_glove[jj+1][i])];
            // If deleting all bonds, first cache then remove all histories
            if (n_histories > 0)
              for (auto &ihistory: histories) {
                for (int n = 0; n < num_bond[atom->map(update_mega_glove[jj+1][i])]; n++)
                  dynamic_cast<FixBondHistory *>(ihistory)->cache_history(atom->map(update_mega_glove[jj+1][i]), n);
                for (int n = 0; n < num_bond[atom->map(update_mega_glove[jj+1][i])]; n++)
                  dynamic_cast<FixBondHistory *>(ihistory)->delete_history(atom->map(update_mega_glove[jj+1][i]), 0);
              }
            num_bond[atom->map(update_mega_glove[jj+1][i])] = 0;
          }
          if (rxn.atoms[j].landlocked == 0) {
            for (int p = num_bond[atom->map(update_mega_glove[jj+1][i])]-1; p > -1 ; p--) {
              for (int n = 0; n < rxn.product->natoms; n++) {
                int nn = rxn.atoms[n].amap[1]-1;
                if (n!=j && bond_atom[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i] && rxn.atoms[n].landlocked == 1) {
                  // Cache history information, shift history, then delete final element
                  if (n_histories > 0)
                    for (auto &ihistory: histories)
                      dynamic_cast<FixBondHistory *>(ihistory)->cache_history(atom->map(update_mega_glove[jj+1][i]), p);
                  for (int m = p; m < num_bond[atom->map(update_mega_glove[jj+1][i])]-1; m++) {
                    bond_type[atom->map(update_mega_glove[jj+1][i])][m] = bond_type[atom->map(update_mega_glove[jj+1][i])][m+1];
                    bond_atom[atom->map(update_mega_glove[jj+1][i])][m] = bond_atom[atom->map(update_mega_glove[jj+1][i])][m+1];
                    if (n_histories > 0)
                      for (auto &ihistory: histories)
                        dynamic_cast<FixBondHistory *>(ihistory)->shift_history(atom->map(update_mega_glove[jj+1][i]),m,m+1);
                  }
                  if (n_histories > 0)
                    for (auto &ihistory: histories)
                      dynamic_cast<FixBondHistory *>(ihistory)->delete_history(atom->map(update_mega_glove[jj+1][i]),
                                                                 num_bond[atom->map(update_mega_glove[jj+1][i])]-1);
                  num_bond[atom->map(update_mega_glove[jj+1][i])]--;
                  delta_bonds--;
                }
              }
            }
          }
        }
      }
      // now let's add the new bond info.
      for (int j = 0; j < rxn.product->natoms; j++) {
        int jj = rxn.atoms[j].amap[1]-1;
        if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
          if (rxn.atoms[j].landlocked == 1)  {
            num_bond[atom->map(update_mega_glove[jj+1][i])] = rxn.product->num_bond[j];
            delta_bonds += rxn.product->num_bond[j];
            for (int p = 0; p < rxn.product->num_bond[j]; p++) {
              bond_type[atom->map(update_mega_glove[jj+1][i])][p] = rxn.product->bond_type[j][p];
              bond_atom[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[rxn.atoms[rxn.product->bond_atom[j][p]-1].amap[1]][i];
              // Check cached history data to see if bond regenerated
              if (n_histories > 0)
                for (auto &ihistory: histories)
                  dynamic_cast<FixBondHistory *>(ihistory)->check_cache(atom->map(update_mega_glove[jj+1][i]), p);
            }
          }
          if (rxn.atoms[j].landlocked == 0) {
            for (int p = 0; p < rxn.product->num_bond[j]; p++) {
              if (rxn.atoms[rxn.product->bond_atom[j][p]-1].landlocked == 1) {
                insert_num = num_bond[atom->map(update_mega_glove[jj+1][i])];
                bond_type[atom->map(update_mega_glove[jj+1][i])][insert_num] = rxn.product->bond_type[j][p];
                bond_atom[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[rxn.atoms[rxn.product->bond_atom[j][p]-1].amap[1]][i];
                // Check cached history data to see if bond regenerated
                if (n_histories > 0)
                  for (auto &ihistory: histories)
                    dynamic_cast<FixBondHistory *>(ihistory)->check_cache(atom->map(update_mega_glove[jj+1][i]), insert_num);
                num_bond[atom->map(update_mega_glove[jj+1][i])]++;
                if (num_bond[atom->map(update_mega_glove[jj+1][i])] > atom->bond_per_atom)
                  error->one(FLERR,"Fix bond/react topology/atom exceed system topology/atom");
                delta_bonds++;
              }
            }
          }
        }
      }
    }

    if (n_histories > 0)
      for (auto &ihistory: histories)
        dynamic_cast<FixBondHistory *>(ihistory)->clear_cache();

    // Angles! First let's delete all angle info:
    if (force->angle) {
      int *num_angle = atom->num_angle;
      int **angle_type = atom->angle_type;
      tagint **angle_atom1 = atom->angle_atom1;
      tagint **angle_atom2 = atom->angle_atom2;
      tagint **angle_atom3 = atom->angle_atom3;

      for (int i = 0; i < update_num_mega; i++) {
        auto &rxn = rxns[(int) update_mega_glove[0][i]];
        for (int j = 0; j < rxn.product->natoms; j++) {
          int jj = rxn.atoms[j].amap[1]-1;
          if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
            if (rxn.atoms[j].landlocked == 1) {
              delta_angle -= num_angle[atom->map(update_mega_glove[jj+1][i])];
              num_angle[atom->map(update_mega_glove[jj+1][i])] = 0;
            }
            if (rxn.atoms[j].landlocked == 0) {
              for (int p = num_angle[atom->map(update_mega_glove[jj+1][i])]-1; p > -1; p--) {
                for (int n = 0; n < rxn.product->natoms; n++) {
                  int nn = rxn.atoms[n].amap[1]-1;
                  if (n!=j && rxn.atoms[n].landlocked == 1 &&
                      (angle_atom1[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i] ||
                       angle_atom2[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i] ||
                       angle_atom3[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i])) {
                    for (int m = p; m < num_angle[atom->map(update_mega_glove[jj+1][i])]-1; m++) {
                      angle_type[atom->map(update_mega_glove[jj+1][i])][m] = angle_type[atom->map(update_mega_glove[jj+1][i])][m+1];
                      angle_atom1[atom->map(update_mega_glove[jj+1][i])][m] = angle_atom1[atom->map(update_mega_glove[jj+1][i])][m+1];
                      angle_atom2[atom->map(update_mega_glove[jj+1][i])][m] = angle_atom2[atom->map(update_mega_glove[jj+1][i])][m+1];
                      angle_atom3[atom->map(update_mega_glove[jj+1][i])][m] = angle_atom3[atom->map(update_mega_glove[jj+1][i])][m+1];
                    }
                    num_angle[atom->map(update_mega_glove[jj+1][i])]--;
                    delta_angle--;
                    break;
                  }
                }
              }
            }
          }
        }
        // now let's add the new angle info.
        if (rxn.product->angleflag) {
          for (int j = 0; j < rxn.product->natoms; j++) {
            int jj = rxn.atoms[j].amap[1]-1;
            if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
              if (rxn.atoms[j].landlocked == 1) {
                num_angle[atom->map(update_mega_glove[jj+1][i])] = rxn.product->num_angle[j];
                delta_angle += rxn.product->num_angle[j];
                for (int p = 0; p < rxn.product->num_angle[j]; p++) {
                  angle_type[atom->map(update_mega_glove[jj+1][i])][p] = rxn.product->angle_type[j][p];
                  angle_atom1[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[rxn.atoms[rxn.product->angle_atom1[j][p]-1].amap[1]][i];
                  angle_atom2[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[rxn.atoms[rxn.product->angle_atom2[j][p]-1].amap[1]][i];
                  angle_atom3[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[rxn.atoms[rxn.product->angle_atom3[j][p]-1].amap[1]][i];
                }
              }
              if (rxn.atoms[j].landlocked == 0) {
                for (int p = 0; p < rxn.product->num_angle[j]; p++) {
                  if (rxn.atoms[rxn.product->angle_atom1[j][p]-1].landlocked == 1 ||
                      rxn.atoms[rxn.product->angle_atom2[j][p]-1].landlocked == 1 ||
                      rxn.atoms[rxn.product->angle_atom3[j][p]-1].landlocked == 1) {
                    insert_num = num_angle[atom->map(update_mega_glove[jj+1][i])];
                    angle_type[atom->map(update_mega_glove[jj+1][i])][insert_num] = rxn.product->angle_type[j][p];
                    angle_atom1[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[rxn.atoms[rxn.product->angle_atom1[j][p]-1].amap[1]][i];
                    angle_atom2[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[rxn.atoms[rxn.product->angle_atom2[j][p]-1].amap[1]][i];
                    angle_atom3[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[rxn.atoms[rxn.product->angle_atom3[j][p]-1].amap[1]][i];
                    num_angle[atom->map(update_mega_glove[jj+1][i])]++;
                    if (num_angle[atom->map(update_mega_glove[jj+1][i])] > atom->angle_per_atom)
                      error->one(FLERR,"Fix bond/react topology/atom exceed system topology/atom");
                    delta_angle++;
                  }
                }
              }
            }
          }
        }
      }
    }

    // Dihedrals! first let's delete all dihedral info for landlocked atoms
    if (force->dihedral) {
      int *num_dihedral = atom->num_dihedral;
      int **dihedral_type = atom->dihedral_type;
      tagint **dihedral_atom1 = atom->dihedral_atom1;
      tagint **dihedral_atom2 = atom->dihedral_atom2;
      tagint **dihedral_atom3 = atom->dihedral_atom3;
      tagint **dihedral_atom4 = atom->dihedral_atom4;

      for (int i = 0; i < update_num_mega; i++) {
        auto &rxn = rxns[(int) update_mega_glove[0][i]];
        for (int j = 0; j < rxn.product->natoms; j++) {
          int jj = rxn.atoms[j].amap[1]-1;
          if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
            if (rxn.atoms[j].landlocked == 1) {
              delta_dihed -= num_dihedral[atom->map(update_mega_glove[jj+1][i])];
              num_dihedral[atom->map(update_mega_glove[jj+1][i])] = 0;
            }
            if (rxn.atoms[j].landlocked == 0) {
              for (int p = num_dihedral[atom->map(update_mega_glove[jj+1][i])]-1; p > -1; p--) {
                for (int n = 0; n < rxn.product->natoms; n++) {
                  int nn = rxn.atoms[n].amap[1]-1;
                  if (n!=j && rxn.atoms[n].landlocked == 1 &&
                      (dihedral_atom1[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i] ||
                       dihedral_atom2[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i] ||
                       dihedral_atom3[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i] ||
                       dihedral_atom4[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i])) {
                    for (int m = p; m < num_dihedral[atom->map(update_mega_glove[jj+1][i])]-1; m++) {
                      dihedral_type[atom->map(update_mega_glove[jj+1][i])][m] = dihedral_type[atom->map(update_mega_glove[jj+1][i])][m+1];
                      dihedral_atom1[atom->map(update_mega_glove[jj+1][i])][m] = dihedral_atom1[atom->map(update_mega_glove[jj+1][i])][m+1];
                      dihedral_atom2[atom->map(update_mega_glove[jj+1][i])][m] = dihedral_atom2[atom->map(update_mega_glove[jj+1][i])][m+1];
                      dihedral_atom3[atom->map(update_mega_glove[jj+1][i])][m] = dihedral_atom3[atom->map(update_mega_glove[jj+1][i])][m+1];
                      dihedral_atom4[atom->map(update_mega_glove[jj+1][i])][m] = dihedral_atom4[atom->map(update_mega_glove[jj+1][i])][m+1];
                    }
                    num_dihedral[atom->map(update_mega_glove[jj+1][i])]--;
                    delta_dihed--;
                    break;
                  }
                }
              }
            }
          }
        }
        // now let's add new dihedral info
        if (rxn.product->dihedralflag) {
          for (int j = 0; j < rxn.product->natoms; j++) {
            int jj = rxn.atoms[j].amap[1]-1;
            if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
              if (rxn.atoms[j].landlocked == 1) {
                num_dihedral[atom->map(update_mega_glove[jj+1][i])] = rxn.product->num_dihedral[j];
                delta_dihed += rxn.product->num_dihedral[j];
                for (int p = 0; p < rxn.product->num_dihedral[j]; p++) {
                  dihedral_type[atom->map(update_mega_glove[jj+1][i])][p] = rxn.product->dihedral_type[j][p];
                  dihedral_atom1[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[rxn.atoms[rxn.product->dihedral_atom1[j][p]-1].amap[1]][i];
                  dihedral_atom2[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[rxn.atoms[rxn.product->dihedral_atom2[j][p]-1].amap[1]][i];
                  dihedral_atom3[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[rxn.atoms[rxn.product->dihedral_atom3[j][p]-1].amap[1]][i];
                  dihedral_atom4[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[rxn.atoms[rxn.product->dihedral_atom4[j][p]-1].amap[1]][i];
                }
              }
              if (rxn.atoms[j].landlocked == 0) {
                for (int p = 0; p < rxn.product->num_dihedral[j]; p++) {
                  if (rxn.atoms[rxn.product->dihedral_atom1[j][p]-1].landlocked == 1 ||
                      rxn.atoms[rxn.product->dihedral_atom2[j][p]-1].landlocked == 1 ||
                      rxn.atoms[rxn.product->dihedral_atom3[j][p]-1].landlocked == 1 ||
                      rxn.atoms[rxn.product->dihedral_atom4[j][p]-1].landlocked == 1) {
                    insert_num = num_dihedral[atom->map(update_mega_glove[jj+1][i])];
                    dihedral_type[atom->map(update_mega_glove[jj+1][i])][insert_num] = rxn.product->dihedral_type[j][p];
                    dihedral_atom1[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[rxn.atoms[rxn.product->dihedral_atom1[j][p]-1].amap[1]][i];
                    dihedral_atom2[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[rxn.atoms[rxn.product->dihedral_atom2[j][p]-1].amap[1]][i];
                    dihedral_atom3[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[rxn.atoms[rxn.product->dihedral_atom3[j][p]-1].amap[1]][i];
                    dihedral_atom4[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[rxn.atoms[rxn.product->dihedral_atom4[j][p]-1].amap[1]][i];
                    num_dihedral[atom->map(update_mega_glove[jj+1][i])]++;
                    if (num_dihedral[atom->map(update_mega_glove[jj+1][i])] > atom->dihedral_per_atom)
                      error->one(FLERR,"Fix bond/react topology/atom exceed system topology/atom");
                    delta_dihed++;
                  }
                }
              }
            }
          }
        }
      }
    }

    // finally IMPROPERS!!!! first let's delete all improper info for landlocked atoms
    if (force->improper) {
      int *num_improper = atom->num_improper;
      int **improper_type = atom->improper_type;
      tagint **improper_atom1 = atom->improper_atom1;
      tagint **improper_atom2 = atom->improper_atom2;
      tagint **improper_atom3 = atom->improper_atom3;
      tagint **improper_atom4 = atom->improper_atom4;

      for (int i = 0; i < update_num_mega; i++) {
        auto &rxn = rxns[(int) update_mega_glove[0][i]];
        for (int j = 0; j < rxn.product->natoms; j++) {
          int jj = rxn.atoms[j].amap[1]-1;
          if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
            if (rxn.atoms[j].landlocked == 1) {
              delta_imprp -= num_improper[atom->map(update_mega_glove[jj+1][i])];
              num_improper[atom->map(update_mega_glove[jj+1][i])] = 0;
            }
            if (rxn.atoms[j].landlocked == 0) {
              for (int p = num_improper[atom->map(update_mega_glove[jj+1][i])]-1; p > -1; p--) {
                for (int n = 0; n < rxn.product->natoms; n++) {
                  int nn = rxn.atoms[n].amap[1]-1;
                  if (n!=j && rxn.atoms[n].landlocked == 1 &&
                      (improper_atom1[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i] ||
                       improper_atom2[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i] ||
                       improper_atom3[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i] ||
                       improper_atom4[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i])) {
                    for (int m = p; m < num_improper[atom->map(update_mega_glove[jj+1][i])]-1; m++) {
                      improper_type[atom->map(update_mega_glove[jj+1][i])][m] = improper_type[atom->map(update_mega_glove[jj+1][i])][m+1];
                      improper_atom1[atom->map(update_mega_glove[jj+1][i])][m] = improper_atom1[atom->map(update_mega_glove[jj+1][i])][m+1];
                      improper_atom2[atom->map(update_mega_glove[jj+1][i])][m] = improper_atom2[atom->map(update_mega_glove[jj+1][i])][m+1];
                      improper_atom3[atom->map(update_mega_glove[jj+1][i])][m] = improper_atom3[atom->map(update_mega_glove[jj+1][i])][m+1];
                      improper_atom4[atom->map(update_mega_glove[jj+1][i])][m] = improper_atom4[atom->map(update_mega_glove[jj+1][i])][m+1];
                    }
                    num_improper[atom->map(update_mega_glove[jj+1][i])]--;
                    delta_imprp--;
                    break;
                  }
                }
              }
            }
          }
        }
        // now let's add new improper info
        if (rxn.product->improperflag) {
          for (int j = 0; j < rxn.product->natoms; j++) {
            int jj = rxn.atoms[j].amap[1]-1;
            if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
              if (rxn.atoms[j].landlocked == 1) {
                num_improper[atom->map(update_mega_glove[jj+1][i])] = rxn.product->num_improper[j];
                delta_imprp += rxn.product->num_improper[j];
                for (int p = 0; p < rxn.product->num_improper[j]; p++) {
                  improper_type[atom->map(update_mega_glove[jj+1][i])][p] = rxn.product->improper_type[j][p];
                  improper_atom1[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[rxn.atoms[rxn.product->improper_atom1[j][p]-1].amap[1]][i];
                  improper_atom2[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[rxn.atoms[rxn.product->improper_atom2[j][p]-1].amap[1]][i];
                  improper_atom3[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[rxn.atoms[rxn.product->improper_atom3[j][p]-1].amap[1]][i];
                  improper_atom4[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[rxn.atoms[rxn.product->improper_atom4[j][p]-1].amap[1]][i];
                }
              }
              if (rxn.atoms[j].landlocked == 0) {
                for (int p = 0; p < rxn.product->num_improper[j]; p++) {
                  if (rxn.atoms[rxn.product->improper_atom1[j][p]-1].landlocked == 1 ||
                      rxn.atoms[rxn.product->improper_atom2[j][p]-1].landlocked == 1 ||
                      rxn.atoms[rxn.product->improper_atom3[j][p]-1].landlocked == 1 ||
                      rxn.atoms[rxn.product->improper_atom4[j][p]-1].landlocked == 1) {
                    insert_num = num_improper[atom->map(update_mega_glove[jj+1][i])];
                    improper_type[atom->map(update_mega_glove[jj+1][i])][insert_num] = rxn.product->improper_type[j][p];
                    improper_atom1[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[rxn.atoms[rxn.product->improper_atom1[j][p]-1].amap[1]][i];
                    improper_atom2[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[rxn.atoms[rxn.product->improper_atom2[j][p]-1].amap[1]][i];
                    improper_atom3[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[rxn.atoms[rxn.product->improper_atom3[j][p]-1].amap[1]][i];
                    improper_atom4[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[rxn.atoms[rxn.product->improper_atom4[j][p]-1].amap[1]][i];
                    num_improper[atom->map(update_mega_glove[jj+1][i])]++;
                    if (num_improper[atom->map(update_mega_glove[jj+1][i])] > atom->improper_per_atom)
                      error->one(FLERR,"Fix bond/react topology/atom exceed system topology/atom");
                    delta_imprp++;
                  }
                }
              }
            }
          }
        }
      }
    }

  }

  memory->destroy(update_mega_glove);
  if (rescale_charges_anyflag) memory->destroy(sim_total_charges);

  // delete atoms. taken from fix_evaporate. but don't think it needs to be in pre_exchange
  // loop in reverse order to avoid copying marked atoms
  ndel = ndelone = 0;
  for (int i = atom->nlocal-1; i >= 0; i--) {
    if (mark[i] == 1) {
      avec->copy(atom->nlocal-1,i,1);
      atom->nlocal--;
      ndelone++;

      if (atom->avec->bonds_allow) {
        if (force->newton_bond) delta_bonds += atom->num_bond[i];
        else {
          for (int j = 0; j < atom->num_bond[i]; j++) {
            if (tag[i] < atom->bond_atom[i][j]) delta_bonds++;
          }
        }
      }
      if (atom->avec->angles_allow) {
        if (force->newton_bond) delta_angle += atom->num_angle[i];
        else {
          for (int j = 0; j < atom->num_angle[i]; j++) {
            int m = atom->map(atom->angle_atom2[i][j]);
            if (m >= 0 && m < nlocal) delta_angle++;
          }
        }
      }
      if (atom->avec->dihedrals_allow) {
        if (force->newton_bond) delta_dihed += atom->num_dihedral[i];
        else {
          for (int j = 0; j < atom->num_dihedral[i]; j++) {
            int m = atom->map(atom->dihedral_atom2[i][j]);
            if (m >= 0 && m < nlocal) delta_dihed++;
          }
        }
      }
      if (atom->avec->impropers_allow) {
        if (force->newton_bond) delta_imprp += atom->num_improper[i];
        else {
          for (int j = 0; j < atom->num_improper[i]; j++) {
            int m = atom->map(atom->improper_atom2[i][j]);
            if (m >= 0 && m < nlocal) delta_imprp++;
          }
        }
      }
    }
  }
  memory->destroy(mark);
  memory->destroy(delflag);

  MPI_Allreduce(&ndelone,&ndel,1,MPI_INT,MPI_SUM,world);

  atom->natoms -= ndel;
  // done deleting atoms
  // something to think about: this could done much more concisely if
  // all atom-level info (bond,angles, etc...) were kinda inherited from a common data struct --JG

  int Tdelta_bonds;
  MPI_Allreduce(&delta_bonds,&Tdelta_bonds,1,MPI_INT,MPI_SUM,world);
  atom->nbonds += Tdelta_bonds;

  int Tdelta_angle;
  MPI_Allreduce(&delta_angle,&Tdelta_angle,1,MPI_INT,MPI_SUM,world);
  atom->nangles += Tdelta_angle;

  int Tdelta_dihed;
  MPI_Allreduce(&delta_dihed,&Tdelta_dihed,1,MPI_INT,MPI_SUM,world);
  atom->ndihedrals += Tdelta_dihed;

  int Tdelta_imprp;
  MPI_Allreduce(&delta_imprp,&Tdelta_imprp,1,MPI_INT,MPI_SUM,world);
  atom->nimpropers += Tdelta_imprp;

  if (ndel && (atom->map_style != Atom::MAP_NONE)) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }
}

/* ----------------------------------------------------------------------
setup for inserting created atoms
atoms for all rxns are actually created all at once in update_everything
------------------------------------------------------------------------- */

int FixBondReact::insert_atoms_setup(tagint **my_update_mega_glove, int iupdate, Reaction &rxn)
{
  // inserting atoms based off fix_deposit->pre_exchange
  int flag;
  imageint *imageflags;
  double **coords,lamda[3],rotmat[3][3];
  double *newcoord;
  double t,delx,dely,delz,rsq;

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int dimension = domain->dimension;

  memory->create(coords,rxn.product->natoms,3,"bond/react:coords");
  memory->create(imageflags,rxn.product->natoms,"bond/react:imageflags");

  double *sublo,*subhi;
  if (domain->triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // only proc that owns reacting atom (use ibonding),
  // fits post-reaction template to reaction site, for creating atoms
  int n2superpose = 0;
  for (int j = 0; j < rxn.product->natoms; j++) {
    if (rxn.modify_create_fragid >= 0)
      if (!rxn.product->fragmentmask[rxn.modify_create_fragid][j]) continue;
    if (!rxn.atoms[j].created && !rxn.atoms[rxn.atoms[j].amap[1]].deleted)
      n2superpose++;
  }

  int ifit = atom->map(my_update_mega_glove[rxn.ibonding+1][iupdate]); // use this local ID to find fitting proc
  Superpose3D<double, double **> superposer(n2superpose);
  int fitroot = 0;
  if (ifit >= 0 && ifit < atom->nlocal) {
    fitroot = comm->me;

    // get 'temperatere' averaged over site, used for created atoms' vels
    // note: row_offset for my_update_mega_glove is unity, not 'cuff'
    std::vector<tagint> myglove(rxn.reactant->natoms);
    for (int i = 0; i < rxn.reactant->natoms; i++) myglove[i] = my_update_mega_glove[i+1][iupdate];
    t = get_temperature(myglove);

    double **xfrozen; // coordinates for the "frozen" target molecule
    double **xmobile; // coordinates for the "mobile" molecule
    memory->create(xfrozen,n2superpose,3,"bond/react:xfrozen");
    memory->create(xmobile,n2superpose,3,"bond/react:xmobile");
    tagint iatom;
    tagint iref = -1; // choose first atom as reference
    int fit_incr = 0;
    for (int j = 0; j < rxn.product->natoms; j++) {
      if (rxn.modify_create_fragid >= 0)
        if (!rxn.product->fragmentmask[rxn.modify_create_fragid][j]) continue;
      int ipre = rxn.atoms[j].amap[1]-1; // equiv pre-reaction template index
      if (!rxn.atoms[j].created && !rxn.atoms[ipre].deleted) {
        if (atom->map(my_update_mega_glove[ipre+1][iupdate]) < 0) {
          error->warning(FLERR," eligible atoms skipped for created-atoms fit on rank {}\n",
                         comm->me);
          continue;
        }
        iatom = atom->map(my_update_mega_glove[ipre+1][iupdate]);
        if (iref == -1) iref = iatom;
        iatom = domain->closest_image(iref,iatom);
        for (int k = 0; k < 3; k++) {
          xfrozen[fit_incr][k] = x[iatom][k];
          xmobile[fit_incr][k] = rxn.product->x[j][k];
        }
        fit_incr++;
      }
    }
    superposer.Superpose(xfrozen, xmobile);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        rotmat[i][j] = superposer.R[i][j];
    memory->destroy(xfrozen);
    memory->destroy(xmobile);
  }
  MPI_Allreduce(MPI_IN_PLACE,&fitroot,1,MPI_INT,MPI_SUM,world);
  MPI_Bcast(&t,1,MPI_DOUBLE,fitroot,world);

  // get coordinates and image flags
  for (int m = 0; m < rxn.product->natoms; m++) {
    if (rxn.atoms[m].created == 1) {
      // apply optimal rotation/translation for created atom coords
      // also map coords back into simulation box
      if (fitroot == comm->me) {
        MathExtra::matvec(rotmat,rxn.product->x[m],coords[m]);
        for (int i = 0; i < 3; i++) coords[m][i] += superposer.T[i];
        imageflags[m] = atom->image[ifit];
        domain->remap(coords[m],imageflags[m]);
      }
      MPI_Bcast(&imageflags[m],1,MPI_LMP_IMAGEINT,fitroot,world);
      MPI_Bcast(coords[m],3,MPI_DOUBLE,fitroot,world);
    }
  }

  // check distance between any existing atom and inserted atom
  // if less than near, abort
  if (rxn.overlapsq > 0) {
    int abortflag = 0;
    for (int m = 0; m < rxn.product->natoms; m++) {
      if (rxn.atoms[m].created == 1) {
        for (int i = 0; i < nlocal; i++) {
          delx = coords[m][0] - x[i][0];
          dely = coords[m][1] - x[i][1];
          delz = coords[m][2] - x[i][2];
          domain->minimum_image(FLERR, delx,dely,delz);
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < rxn.overlapsq) {
            abortflag = 1;
            break;
          }
        }
        if (abortflag) break;
      }
    }
    // also check against previous to-be-added atoms
    if (!abortflag) {
      for (auto & myaddatom : addatoms) {
        for (int m = 0; m < rxn.product->natoms; m++) {
          if (rxn.atoms[m].created == 1) {
            delx = coords[m][0] - myaddatom.x[0];
            dely = coords[m][1] - myaddatom.x[1];
            delz = coords[m][2] - myaddatom.x[2];
            domain->minimum_image(FLERR, delx,dely,delz);
            rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < rxn.overlapsq) {
              abortflag = 1;
              break;
            }
          }
        }
        if (abortflag) break;
      }
    }

    MPI_Allreduce(MPI_IN_PLACE,&abortflag,1,MPI_INT,MPI_MAX,world);
    if (abortflag) {
      memory->destroy(coords);
      memory->destroy(imageflags);
      return 0;
    }
  }

  // check if new atoms are in my sub-box or above it if I am highest proc
  // if so, add atom to my list via create_atom()
  // initialize additional info about the atoms
  // set group mask to "all" plus fix group
  int preID; // new amap index
  int add_count = 0;
  for (int m = 0; m < rxn.product->natoms; m++) {
    if (rxn.atoms[m].created == 1) {
      // increase atom count
      add_count++;
      preID = rxn.reactant->natoms+add_count;

      if (domain->triclinic) {
        domain->x2lamda(coords[m],lamda);
        newcoord = lamda;
      } else newcoord = coords[m];

      flag = 0;
      if (newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
          newcoord[1] >= sublo[1] && newcoord[1] < subhi[1] &&
          newcoord[2] >= sublo[2] && newcoord[2] < subhi[2]) flag = 1;
      else if (dimension == 3 && newcoord[2] >= domain->boxhi[2]) {
        if (comm->layout != Comm::LAYOUT_TILED) {
          if (comm->myloc[2] == comm->procgrid[2]-1 &&
              newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
              newcoord[1] >= sublo[1] && newcoord[1] < subhi[1]) flag = 1;
        } else {
          if (comm->mysplit[2][1] == 1.0 &&
              newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
              newcoord[1] >= sublo[1] && newcoord[1] < subhi[1]) flag = 1;
        }
      } else if (dimension == 2 && newcoord[1] >= domain->boxhi[1]) {
        if (comm->layout != Comm::LAYOUT_TILED) {
          if (comm->myloc[1] == comm->procgrid[1]-1 &&
              newcoord[0] >= sublo[0] && newcoord[0] < subhi[0]) flag = 1;
        } else {
          if (comm->mysplit[1][1] == 1.0 &&
              newcoord[0] >= sublo[0] && newcoord[0] < subhi[0]) flag = 1;
        }
      }

      int root = 0;
      addatomtag++;
      if (flag) {
        struct AddAtom myaddatom;
        root = comm->me;

        myaddatom.type = rxn.product->type[m];
        myaddatom.x[0] = coords[m][0];
        myaddatom.x[1] = coords[m][1];
        myaddatom.x[2] = coords[m][2];
        myaddatom.tag = addatomtag;

        // locally update mega_glove
        my_update_mega_glove[preID][iupdate] = myaddatom.tag;

        myaddatom.mask = 1 | groupbit;
        myaddatom.image = imageflags[m];
        if (atom->molecule_flag) myaddatom.molecule = 0;

        // guess a somewhat reasonable initial velocity based on reaction site
        // further control is possible using bond_react_MASTER_group
        // compute |velocity| corresponding to a given temperature t, using specific atom's mass
        myaddatom.rmass = atom->rmass ? rxn.product->rmass[m] : atom->mass[rxn.product->type[m]];
        double vtnorm = sqrt(t / (force->mvv2e / (dimension * force->boltz)) / myaddatom.rmass);
        double myv[3];
        myv[0] = random[rxn.ID]->uniform();
        myv[1] = random[rxn.ID]->uniform();
        myv[2] = random[rxn.ID]->uniform();
        double vnorm = sqrt(myv[0]*myv[0] + myv[1]*myv[1] + myv[2]*myv[2]);
        myaddatom.v[0] = myv[0]/vnorm*vtnorm;
        myaddatom.v[1] = myv[1]/vnorm*vtnorm;
        myaddatom.v[2] = myv[2]/vnorm*vtnorm;
        addatoms.push_back(myaddatom);
      }
      // globally update mega_glove and amap
      MPI_Allreduce(MPI_IN_PLACE,&root,1,MPI_INT,MPI_SUM,world);
      MPI_Bcast(&my_update_mega_glove[preID][iupdate],1,MPI_LMP_TAGINT,root,world);
      rxn.atoms[m].amap[0] = m+1;
      rxn.atoms[m].amap[1] = preID;
      rxn.atoms[preID-1].ramap[0] = preID;
      rxn.atoms[preID-1].ramap[1] = m+1;
    }
  }

  // reset global natoms here
  // reset atom map elsewhere, after all calls to 'insert_atoms_setup'
  atom->natoms += add_count;
  if (atom->natoms < 0)
    error->all(FLERR,"Too many total atoms");
  if (addatomtag >= MAXTAGINT)
    error->all(FLERR,"New atom IDs exceed maximum allowed ID");
  // atom creation successful
  memory->destroy(coords);
  memory->destroy(imageflags);
  return 1;
}

/* ----------------------------------------------------------------------
add equal-style variable to keyword argument list
------------------------------------------------------------------------- */

void FixBondReact::validate_variable_keyword(const char *myarg, int var_id)
{
  if (var_id < 0)
    error->all(FLERR,"Fix bond/react: Variable name {} does not exist",myarg);
  if (!input->variable->equalstyle(var_id))
    error->all(FLERR,"Fix bond/react: Variable {} is not equal-style",myarg);
}

/* ----------------------------------------------------------------------
read map file
------------------------------------------------------------------------- */

void FixBondReact::read_map_file(Reaction &rxn)
{
  int rv, nedge, nequivalent, nchiral, ndelete, ncreate = 0;
  char line[MAXLINE] = {'\0'};
  char keyword[MAXLINE] = {'\0'};
  char *eof,*ptr;

  // skip 1st line of file
  eof = fgets(line,MAXLINE,fp);
  if (eof == nullptr) error->one(FLERR,"Fix bond/react: Unexpected end of superimpose file");

  // read header lines
  // skip blank lines or lines that start with "#"
  // stop when read an unrecognized line

  while (true) {

    readline(line);

    // trim anything from '#' onward
    // if line is blank, continue

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    if (strstr(line,"edgeIDs")) sscanf(line,"%d",&nedge);
    else if (strstr(line,"equivalences")) {
      rv = sscanf(line,"%d",&nequivalent);
      if (rv != 1) error->one(FLERR, "Map file header is incorrectly formatted");
      if (nequivalent != rxn.reactant->natoms)
        error->one(FLERR,"Fix bond/react: Number of equivalences in map file must "
                   "equal number of atoms in reaction templates");
    }
    else if (strstr(line,"deleteIDs")) {
      rv = sscanf(line,"%d",&ndelete);
      if (rv != 1) error->one(FLERR, "Map file header is incorrectly formatted");
    } else if (strstr(line,"createIDs")) {
      rv = sscanf(line,"%d",&ncreate);
      if (rv != 1) error->one(FLERR, "Map file header is incorrectly formatted");
    } else if (strstr(line,"chiralIDs")) {
      rv = sscanf(line,"%d",&nchiral);
      if (rv != 1) error->one(FLERR, "Map file header is incorrectly formatted");
    } else if (strstr(line,"constraints")) {
      int nconstraints;
      rv = sscanf(line,"%d",&nconstraints);
      if (rv != 1) error->one(FLERR, "Map file header is incorrectly formatted");
      rxn.constraints.resize(nconstraints);
      for (int i = 0; i < nconstraints; i++) rxn.constraints[i].ID = i;
    } else break;
  }

  if (ncreate == 0 && rxn.reactant->natoms != rxn.product->natoms)
    error->all(FLERR,"Fix bond/react: Reaction templates must contain the same number of atoms");
  else if (ncreate > 0 && rxn.reactant->natoms + ncreate != rxn.product->natoms)
    error->all(FLERR,"Fix bond/react: Incorrect number of created atoms");

  // grab keyword and skip next line

  parse_keyword(0,line,keyword);
  readline(line);

  // loop over sections of superimpose file

  int equivflag = 0, bondflag = 0;
  while (strlen(keyword)) {
    if (strcmp(keyword,"InitiatorIDs") == 0 || strcmp(keyword,"BondingIDs") == 0) {
      if (strcmp(keyword,"BondingIDs") == 0)
        if (comm->me == 0) error->warning(FLERR,"Fix bond/react: The BondingIDs section title has been deprecated. Please use InitiatorIDs instead.");
      bondflag = 1;
      readline(line);
      rv = sscanf(line,"%d",&rxn.ibonding);
      if (rv != 1) error->one(FLERR, "InitiatorIDs section is incorrectly formatted");
      if (rxn.ibonding > rxn.reactant->natoms)
        error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
      readline(line);
      rv = sscanf(line,"%d",&rxn.jbonding);
      if (rv != 1) error->one(FLERR, "InitiatorIDs section is incorrectly formatted");
      if (rxn.jbonding > rxn.reactant->natoms)
        error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    } else if (strcmp(keyword,"EdgeIDs") == 0) {
      EdgeIDs(line, rxn, nedge);
    } else if (strcmp(keyword,"Equivalences") == 0) {
      equivflag = 1;
      Equivalences(line, rxn, nequivalent);
    } else if (strcmp(keyword,"DeleteIDs") == 0) {
      DeleteAtoms(line, rxn, ndelete);
    } else if (strcmp(keyword,"CreateIDs") == 0) {
      CreateAtoms(line, rxn, ncreate);
    } else if (strcmp(keyword,"ChiralIDs") == 0) {
      ChiralCenters(line, rxn, nchiral);
    } else if (strcmp(keyword,"Constraints") == 0) {
      ReadConstraints(line, rxn);
    } else error->one(FLERR,"Fix bond/react: Unknown section in map file");

    parse_keyword(1,line,keyword);

  }

  // error check
  for (int i = 0; i < rxn.reactant->natoms; i++) {
    int my_equiv = rxn.atoms[i].ramap[1];
    if (rxn.atoms[my_equiv-1].created == 1)
      error->all(FLERR,"Fix bond/react: Created atoms cannot also be listed in Equivalences section\n");
  }

  // error check
  if (bondflag == 0 || equivflag == 0)
    error->all(FLERR,"Fix bond/react: Map file missing InitiatorIDs or Equivalences section\n");
}

void FixBondReact::EdgeIDs(char *line, Reaction &rxn, int nedge)
{
  // puts a 1 at edge(edgeID)

  int tmp,rv;
  for (int i = 0; i < nedge; i++) {
    readline(line);
    rv = sscanf(line,"%d",&tmp);
    if (rv != 1) error->one(FLERR, "EdgeIDs section is incorrectly formatted");
    if (tmp > rxn.reactant->natoms)
      error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    rxn.atoms[tmp-1].edge = 1;
  }
}

void FixBondReact::Equivalences(char *line, Reaction &rxn, int nequivalent)
{
  int tmp1,tmp2,rv;
  for (int i = 0; i < nequivalent; i++) {
    readline(line);
    rv = sscanf(line,"%d %d",&tmp1,&tmp2);
    if (rv != 2) error->one(FLERR, "Equivalences section is incorrectly formatted");
    if (tmp1 > rxn.reactant->natoms || tmp2 > rxn.product->natoms)
      error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    //equivalences is-> clmn 1: post-reacted, clmn 2: pre-reacted
    rxn.atoms[tmp2-1].amap[0] = tmp2;
    rxn.atoms[tmp2-1].amap[1] = tmp1;
    //reverse_equiv is-> clmn 1: pre-reacted, clmn 2: post-reacted
    rxn.atoms[tmp1-1].ramap[0] = tmp1;
    rxn.atoms[tmp1-1].ramap[1] = tmp2;
  }
  // sanity check for one-to-one mapping for equivalences
  for (int i = 0; i < rxn.product->natoms; i++) {
    if (rxn.atoms[i].created == 1) continue;
    for (int j = i+1; j < rxn.product->natoms; j++) {
      if (rxn.atoms[j].created == 1) continue;
      if (rxn.atoms[i].amap[0] == rxn.atoms[j].amap[0] ||
          rxn.atoms[i].amap[1] == rxn.atoms[j].amap[1]) {
        error->one(FLERR,"Fix bond/react: Repeated atoms IDs in Equivalences section");
      }
    }
  }
}

void FixBondReact::DeleteAtoms(char *line, Reaction &rxn, int ndelete)
{
  int tmp,rv;
  for (int i = 0; i < ndelete; i++) {
    readline(line);
    rv = sscanf(line,"%d",&tmp);
    if (rv != 1) error->one(FLERR, "DeleteIDs section is incorrectly formatted");
    if (tmp > rxn.reactant->natoms)
      error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    rxn.atoms[tmp-1].deleted = 1;
  }
}

void FixBondReact::CreateAtoms(char *line, Reaction &rxn, int ncreate)
{
  rxn.create_atoms_flag = 1;
  int tmp,rv;
  for (int i = 0; i < ncreate; i++) {
    readline(line);
    rv = sscanf(line,"%d",&tmp);
    if (rv != 1) error->one(FLERR, "CreateIDs section is incorrectly formatted");
    if (tmp > rxn.product->natoms)
      error->one(FLERR,"Fix bond/react: Invalid atom ID in CreateIDs section of map file");
    rxn.atoms[tmp-1].created = 1;
  }
  if (rxn.product->xflag == 0)
    error->one(FLERR,"Fix bond/react: 'Coords' section required in post-reaction template when creating new atoms");
  if (atom->rmass_flag && !rxn.product->rmassflag)
    error->one(FLERR, "Fix bond/react: 'Masses' section required in post-reaction template when creating new atoms if per-atom masses are defined.");
}

void FixBondReact::CustomCharges(int ifragment, Reaction &rxn)
{
  for (int i = 0; i < rxn.reactant->natoms; i++)
    if (rxn.reactant->fragmentmask[ifragment][i])
      rxn.atoms[i].recharged = 1;
    else
      rxn.atoms[i].recharged = 0;
}

void FixBondReact::ChiralCenters(char *line, Reaction &rxn, int nchiral)
{
  int tmp,rv;
  for (int i = 0; i < nchiral; i++) {
    readline(line);
    rv = sscanf(line,"%d",&tmp);
    if (rv != 1) error->one(FLERR, "ChiralIDs section is incorrectly formatted");
    if (tmp > rxn.reactant->natoms)
      error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    rxn.atoms[tmp-1].chiral[0] = 1;
    if (rxn.reactant->xflag == 0)
      error->one(FLERR,"Fix bond/react: Molecule template 'Coords' section required for chiralIDs keyword");
    if ((int) rxn.reactant->nspecial[tmp-1][0] != 4)
      error->one(FLERR,"Fix bond/react: Chiral atoms must have exactly four first neighbors");
    for (int j = 0; j < 4; j++) {
      for (int k = j+1; k < 4; k++) {
        if (rxn.reactant->type[rxn.reactant->special[tmp-1][j]-1] ==
            rxn.reactant->type[rxn.reactant->special[tmp-1][k]-1])
          error->one(FLERR,"Fix bond/react: First neighbors of chiral atoms must be of mutually different types");
      }
    }
    // record order of atom types, and coords
    double my4coords[12];
    for (int j = 0; j < 4; j++) {
      rxn.atoms[tmp-1].chiral[j+2] = rxn.reactant->type[rxn.reactant->special[tmp-1][j]-1];
      for (int k = 0; k < 3; k++) {
        my4coords[3*j+k] = rxn.reactant->x[rxn.reactant->special[tmp-1][j]-1][k];
      }
    }
    // get orientation
    rxn.atoms[tmp-1].chiral[1] = get_chirality(my4coords);
  }
}

void FixBondReact::ReadConstraints(char *line, Reaction &rxn)
{
  int rv;
  double tmp[MAXCONARGS];
  char **strargs,*ptr,*lptr;
  memory->create(strargs,MAXCONARGS,MAXLINE,"bond/react:strargs");
  auto *constraint_type = new char[MAXLINE];
  rxn.constraintstr = "("; // string for boolean constraint logic
  for (auto &constraint : rxn.constraints) {
    readline(line);
    // find left parentheses, add to constraintstr, and update line
    for (int j = 0; j < (int)strlen(line); j++) {
      if (line[j] == '(') rxn.constraintstr += "(";
      if (isalpha(line[j])) {
        line = line + j;
        break;
      }
    }
    // 'C' indicates where to sub in next constraint
    rxn.constraintstr += "C";
    // special consideration for 'custom' constraint
    // find final double quote, or skip two words
    lptr = line;
    if ((ptr = strrchr(lptr,'\"'))) lptr = ptr+1;
    else {
      while (lptr[0] != ' ') lptr++; // skip first 'word'
      while (lptr[0] == ' ' || lptr[0] == '\t') lptr++; // skip blanks
      while (lptr[0] != ' ') lptr++; // skip second 'word'
    }
    // find right parentheses
    for (int j = 0; j < (int)strlen(lptr); j++)
      if (lptr[j] == ')') rxn.constraintstr += ")";
    // find logic symbols, and trim line via ptr
    if ((ptr = strstr(lptr,"&&"))) {
      rxn.constraintstr += "&&";
      *ptr = '\0';
    } else if ((ptr = strstr(lptr,"||"))) {
      rxn.constraintstr += "||";
      *ptr = '\0';
    } else if (constraint.ID+1 < rxn.constraints.size()) {
      rxn.constraintstr += "&&";
    }
    if ((ptr = strchr(lptr,')')))
      *ptr = '\0';
    rv = sscanf(line,"%s",constraint_type);
    if (rv != 1) error->one(FLERR, "Constraints section is incorrectly formatted");
    if (strcmp(constraint_type,"distance") == 0) {
      constraint.type = Reaction::Constraint::Type::DISTANCE;
      rv = sscanf(line,"%*s %s %s %lg %lg",strargs[0],strargs[1],&tmp[0],&tmp[1]);
      if (rv != 4) error->one(FLERR, "Distance constraint is incorrectly formatted");
      readID(strargs[0], constraint, rxn, 0);
      readID(strargs[1], constraint, rxn, 1);
      // cutoffs
      constraint.distance.rminsq = tmp[0]*tmp[0]; // using square of distance
      constraint.distance.rmaxsq = tmp[1]*tmp[1];
    } else if (strcmp(constraint_type,"angle") == 0) {
      constraint.type = Reaction::Constraint::Type::ANGLE;
      rv = sscanf(line,"%*s %s %s %s %lg %lg",strargs[0],strargs[1],strargs[2],&tmp[0],&tmp[1]);
      if (rv != 5) error->one(FLERR, "Angle constraint is incorrectly formatted");
      readID(strargs[0], constraint, rxn, 0);
      readID(strargs[1], constraint, rxn, 1);
      readID(strargs[2], constraint, rxn, 2);
      constraint.angle.amin = tmp[0]/180.0 * MY_PI;
      constraint.angle.amax = tmp[1]/180.0 * MY_PI;
    } else if (strcmp(constraint_type,"dihedral") == 0) {
      constraint.type = Reaction::Constraint::Type::DIHEDRAL;
      tmp[2] = 181.0; // impossible range
      tmp[3] = 182.0;
      rv = sscanf(line,"%*s %s %s %s %s %lg %lg %lg %lg",strargs[0],strargs[1],
             strargs[2],strargs[3],&tmp[0],&tmp[1],&tmp[2],&tmp[3]);
      if (rv != 6 && rv != 8) error->one(FLERR, "Dihedral constraint is incorrectly formatted");
      readID(strargs[0], constraint, rxn, 0);
      readID(strargs[1], constraint, rxn, 1);
      readID(strargs[2], constraint, rxn, 2);
      readID(strargs[3], constraint, rxn, 3);
      constraint.dihedral.amin = tmp[0]/180.0 * MY_PI;
      constraint.dihedral.amax = tmp[1]/180.0 * MY_PI;
      constraint.dihedral.amin2 = tmp[2]/180.0 * MY_PI;
      constraint.dihedral.amax2 = tmp[3]/180.0 * MY_PI;
    } else if (strcmp(constraint_type,"arrhenius") == 0) {
      constraint.type = Reaction::Constraint::Type::ARRHENIUS;
      rv = sscanf(line,"%*s %lg %lg %lg %lg",&tmp[0],&tmp[1],&tmp[2],&tmp[3]);
      if (rv != 4) error->one(FLERR, "Arrhenius constraint is incorrectly formatted");
      constraint.arrhenius.A = tmp[0];
      constraint.arrhenius.n = tmp[1];
      constraint.arrhenius.E_a = tmp[2];
      constraint.arrhenius.seed = tmp[3];
    } else if (strcmp(constraint_type,"rmsd") == 0) {
      constraint.type = Reaction::Constraint::Type::RMSD;
      strcpy(strargs[0],"0");
      rv = sscanf(line,"%*s %lg %s",&tmp[0],strargs[0]);
      if (rv != 1 && rv != 2) error->one(FLERR, "RMSD constraint is incorrectly formatted");
      constraint.rmsd.rmsdmax = tmp[0];
      constraint.ids[0] = -1; // optional molecule fragment
      if (isalpha(strargs[0][0])) {
        int ifragment = rxn.reactant->findfragment(strargs[0]);
        if (ifragment < 0) error->one(FLERR,"Fix bond/react: Molecule fragment does not exist");
        else constraint.ids[0] = ifragment;
      }
    } else if (strcmp(constraint_type,"custom") == 0) {
      constraint.type = Reaction::Constraint::Type::CUSTOM;
      std::vector<std::string> args = utils::split_words(line);
      constraint.custom.str = args[1];
    } else error->one(FLERR,"Fix bond/react: Illegal constraint type in 'Constraints' section of map file");
  }
  rxn.constraintstr += ")"; // close boolean constraint logic string
  delete[] constraint_type;
  memory->destroy(strargs);
}

/* ----------------------------------------------------------------------
if ID starts with character, assume it is a pre-reaction molecule fragment ID
otherwise, it is a pre-reaction atom ID
---------------------------------------------------------------------- */

void FixBondReact::readID(char *strarg, Reaction::Constraint &constraint, Reaction &rxn, int i)
{
  if (isalpha(strarg[0])) {
    constraint.idtypes[i] = Reaction::Constraint::IDType::FRAG; // fragment vs. atom ID flag
    int ifragment = rxn.reactant->findfragment(strarg);
    if (ifragment < 0)
      error->one(FLERR,"Fix bond/react: Molecule fragment {} does not exist", strarg);
    constraint.ids[i] = ifragment;
  } else {
    constraint.idtypes[i] = Reaction::Constraint::IDType::ATOM; // fragment vs. atom ID flag
    int iatom = utils::inumeric(FLERR, strarg, true, lmp);
    if (iatom > rxn.reactant->natoms)
      error->one(FLERR,"Fix bond/react: Invalid template atom ID {} in map file", strarg);
    constraint.ids[i] = iatom;
  }
}

void FixBondReact::readline(char *line)
{
  int n;
  if (comm->me == 0) {
    if (fgets(line,MAXLINE,fp) == nullptr) n = 0;
    else n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  if (n == 0) error->all(FLERR,"Fix bond/react: Unexpected end of map file");
  MPI_Bcast(line,n,MPI_CHAR,0,world);
}

void FixBondReact::parse_keyword(int flag, char *line, char *keyword)
{
  if (flag) {

    // read upto non-blank line plus 1 following line
    // eof is set to 1 if any read hits end-of-file

    int eof = 0;
    if (comm->me == 0) {
      if (fgets(line,MAXLINE,fp) == nullptr) eof = 1;
      while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
        if (fgets(line,MAXLINE,fp) == nullptr) eof = 1;
      }
      if (fgets(keyword,MAXLINE,fp) == nullptr) eof = 1;
    }

    // if eof, set keyword empty and return

    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) {
      keyword[0] = '\0';
      return;
    }

    // bcast keyword line to all procs

    int n;
    if (comm->me == 0) n = strlen(line) + 1;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);
  }

  // copy non-whitespace portion of line into keyword

  int start = strspn(line," \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t'
         || line[stop] == '\n' || line[stop] == '\r') stop--;
  line[stop+1] = '\0';
  strcpy(keyword,&line[start]);
}

/* ---------------------------------------------------------------------- */

double FixBondReact::compute_vector(int n)
{
  // now we print just the totals for each reaction instance
  return (double) rxns[n].reaction_count_total;

}

/* ---------------------------------------------------------------------- */

void FixBondReact::post_integrate_respa(int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_integrate();
}

/* ---------------------------------------------------------------------- */

void FixBondReact::post_force(int /*vflag*/)
{
  if (molid_mode == Reset_Mol_IDs::YES) reset_mol_ids->reset();
}

/* ---------------------------------------------------------------------- */

int FixBondReact::pack_forward_comm(int n, int *list, double *buf,
                                    int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,k,m,ns;

  m = 0;

  if (commflag == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      for (k = 0; k < ncustomvars; k++)
        buf[m++] = vvec[j][k];
    }
    return m;
  }

  if (commflag == 2) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(partner[j]).d;
    }
    return m;
  }

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = ubuf(finalpartner[j]).d;
    ns = nxspecial[j][0];
    buf[m++] = ubuf(ns).d;
    for (k = 0; k < ns; k++)
      buf[m++] = ubuf(xspecial[j][k]).d;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixBondReact::unpack_forward_comm(int n, int first, double *buf)
{
  int i,j,k,m,ns,last;

  m = 0;
  last = first + n;

  if (commflag == 1) {
    for (i = first; i < last; i++)
      for (k = 0; k < ncustomvars; k++)
        vvec[i][k] = buf[m++];
  } else if (commflag == 2) {
    for (i = first; i < last; i++)
      partner[i] = (tagint) ubuf(buf[m++]).i;
  } else {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      finalpartner[i] = (tagint) ubuf(buf[m++]).i;
      ns = (int) ubuf(buf[m++]).i;
      nxspecial[i][0] = ns;
      for (j = 0; j < ns; j++)
        xspecial[i][j] = (tagint) ubuf(buf[m++]).i;
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixBondReact::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  for (i = first; i < last; i++) {
    buf[m++] = ubuf(partner[i]).d;
    if (rxnptr->closeneigh != 0)
      buf[m++] = distsq[i][1];
    else
      buf[m++] = distsq[i][0];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixBondReact::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    if (rxnptr->closeneigh != 0) {
      if (buf[m+1] < distsq[j][1]) {
        partner[j] = (tagint) ubuf(buf[m++]).i;
        distsq[j][1] = buf[m++];
      } else m += 2;
    } else {
      if (buf[m+1] > distsq[j][0]) {
        partner[j] = (tagint) ubuf(buf[m++]).i;
        distsq[j][0] = buf[m++];
      } else m += 2;
    }
  }
}

/* ----------------------------------------------------------------------
   write Set data to restart file
------------------------------------------------------------------------- */

void FixBondReact::write_restart(FILE *fp)
{
  int revision = 2;
  set[0].nrxns = rxns.size();
  set[0].nratelimits = rate_limits.size();

  for (int i = 0; i < rxns.size(); i++) {
    set[i].reaction_count_total = rxns[i].reaction_count_total;

    strncpy(set[i].rxn_name,rxns[i].name.c_str(),MAXNAME-1);
    set[i].rxn_name[MAXNAME-1] = '\0';
  }

  // to store, for each RateLimit: Nrxns rxn_IDs[Nrxns] NSteps store_rxn_counts[Nsteps]
  // NOTE: rxn_IDs only valid in reference to this restart file's reaction list
  int rbufcount = rate_limits.size()*2;
  for (auto rlm : rate_limits)
    rbufcount += rlm.Nsteps + rlm.Nrxns;

  int ii = 0;
  int *rbuf;
  if (rbufcount) {
    memory->create(rbuf,rbufcount,"bond/react:rbuf");
    for (auto &rlm : rate_limits) {
      rbuf[ii++] = rlm.Nrxns; // need memcpy?
      for (auto myrxnID : rlm.rxnIDs) rbuf[ii++] = myrxnID;
      rbuf[ii++] = rlm.Nsteps;
      for (auto mycount : rlm.store_rxn_counts) rbuf[ii++] = mycount;
    }
  }

  if (comm->me == 0) {
    int size = rxns.size()*sizeof(Set)+(rbufcount+2)*sizeof(int);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(&revision,sizeof(int),1,fp);
    fwrite(set,sizeof(Set),rxns.size(),fp);
    fwrite(&rbufcount,sizeof(int),1,fp);
    if (rbufcount) fwrite(rbuf,sizeof(int),rbufcount,fp);
  }
  if (rbufcount) memory->destroy(rbuf);
}

/* ----------------------------------------------------------------------
   use selected state info from restart file to restart the Fix
   bond/react restart revisions numbers added after LAMMPS version 3 Nov 2022
------------------------------------------------------------------------- */

void FixBondReact::restart(char *buf)
{
  int revision, r_nrxns, r_nratelimits, ibufcount;
  int iptr = 0;
  int *ibuf;

  if (lmp->restart_ver > utils::date2num("3 Nov 2022")) {
    revision = buf[iptr];
    iptr += sizeof(int);
  } else revision = 0;

  Set *set_restart = (Set *) &buf[iptr];
  r_nrxns = set_restart[0].nrxns;
  iptr += sizeof(Set)*r_nrxns;

  for (int i = 0; i < r_nrxns; i++)
    for (int j = 0; j < rxns.size(); j++)
      if (strcmp(set_restart[i].rxn_name,rxns[j].name.c_str()) == 0)
        rxns[j].reaction_count_total = set_restart[i].reaction_count_total;

  if (revision > 1) {
    std::vector<RateLimit> restart_rate_limits;
    r_nratelimits = set_restart[0].nratelimits;
    ibufcount = buf[iptr];
    iptr += sizeof(int);
    if (ibufcount > 0) {
      memory->create(ibuf,ibufcount,"bond/react:ibuf");
      memcpy(&ibuf[0],&buf[iptr],sizeof(int)*ibufcount);
    }
    int ii = 0;
    for (int i = 0; i < r_nratelimits; i++) {
      struct RateLimit r_rlm;
      r_rlm.Nrxns = ibuf[ii++];
      for (int i = 0; i < r_rlm.Nrxns; i++) {
        r_rlm.rxnIDs.push_back(ibuf[ii++]);
        std::string myrxn_name = set_restart[r_rlm.rxnIDs[i]].rxn_name;
        r_rlm.rxn_names.push_back(myrxn_name);
      }
      r_rlm.Nsteps = ibuf[ii++];
      for (int i = 0; i < r_rlm.Nsteps; i++) r_rlm.store_rxn_counts.push_back(ibuf[ii++]);
      restart_rate_limits.push_back(r_rlm);
    }
    // restore rate_limits store_rxn_counts if all rxn_names match
    // assumes there are no repeats - perhaps should error-check this?
    for (auto &rlm : rate_limits) {
      for (auto r_rlm : restart_rate_limits) {
        if (rlm.Nrxns != r_rlm.Nrxns) continue;
        int nmatch = 0;
        for (int i = 0; i < rlm.Nrxns; i++)
          for (int j = 0; j < r_rlm.Nrxns; j++)
            if (rlm.rxn_names[i] == r_rlm.rxn_names[j]) nmatch++;
        if (nmatch == rlm.Nrxns)
          std::copy(r_rlm.store_rxn_counts.begin(), r_rlm.store_rxn_counts.end(), rlm.store_rxn_counts.begin());
      }
    }
    if (ibufcount > 0) memory->destroy(ibuf);
  }
}

/* ----------------------------------------------------------------------
memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixBondReact::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = (double)nmax * sizeof(int);
  bytes = 2*nmax * sizeof(tagint);
  bytes += (double)nmax * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixBondReact::print_bb()
{
#if 0
  //fix bond/create cargo code. eg nbonds needs to be added

for (int i = 0; i < atom->nlocal; i++) {
  // printf("TAG " TAGINT_FORMAT ": %d nbonds: ",atom->tag[i],atom->num_bond[i]);
  for (int j = 0; j < atom->num_bond[i]; j++) {
  // printf(" " TAGINT_FORMAT,atom->bond_atom[i][j]);
  }
  // printf("\n");
  // printf("TAG " TAGINT_FORMAT ": %d nangles: ",atom->tag[i],atom->num_angle[i]);
  for (int j = 0; j < atom->num_angle[i]; j++) {
  // printf(" " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT ",",
      atom->angle_atom1[i][j], atom->angle_atom2[i][j],
      atom->angle_atom3[i][j]);
  }
  // printf("\n");
  // printf("TAG " TAGINT_FORMAT ": %d ndihedrals: ",atom->tag[i],atom->num_dihedral[i]);
  for (int j = 0; j < atom->num_dihedral[i]; j++) {
  // printf(" " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " "
      TAGINT_FORMAT ",", atom->dihedral_atom1[i][j],
    atom->dihedral_atom2[i][j],atom->dihedral_atom3[i][j],
    atom->dihedral_atom4[i][j]);
  }
  // printf("\n");
  // printf("TAG " TAGINT_FORMAT ": %d nimpropers: ",atom->tag[i],atom->num_improper[i]);
  for (int j = 0; j < atom->num_improper[i]; j++) {
  // printf(" " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " "
      TAGINT_FORMAT ",",atom->improper_atom1[i][j],
    atom->improper_atom2[i][j],atom->improper_atom3[i][j],
    atom->improper_atom4[i][j]);
  }
  // printf("\n");
  // printf("TAG " TAGINT_FORMAT ": %d %d %d nspecial: ",atom->tag[i],
  atom->nspecial[i][0],atom->nspecial[i][1],atom->nspecial[i][2]);
  for (int j = 0; j < atom->nspecial[i][2]; j++) {
    printf(" " TAGINT_FORMAT,atom->special[i][j]);
  }
  // printf("\n");
}
#endif
}
