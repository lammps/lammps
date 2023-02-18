// clang-format off
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
Contributing Author: Jacob Gissinger (jacob.r.gissinger@gmail.com)
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
  "fix bond/react: reacter.org doi:10.1016/j.polymer.2017.09.038, doi:10.1021/acs.macromol.0c02012\n\n"
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
  "}\n\n";

#define BIG 1.0e20
#define DELTA 16
#define MAXGUESS 20 // max # of guesses allowed by superimpose algorithm
#define MAXCONARGS 14 // max # of arguments for any type of constraint + rxnID
#define NUMVARVALS 5 // max # of keyword values that have variables as input

// various statuses of superimpose algorithm:
// ACCEPT: site successfully matched to pre-reacted template
// REJECT: site does not match pre-reacted template
// PROCEED: normal execution (non-guessing mode)
// CONTINUE: a neighbor has been assigned, skip to next neighbor
// GUESSFAIL: a guess has failed (if no more restore points, status = 'REJECT')
// RESTORE: restore mode, load most recent restore point
enum{ACCEPT,REJECT,PROCEED,CONTINUE,GUESSFAIL,RESTORE};

// types of available reaction constraints
enum{DISTANCE,ANGLE,DIHEDRAL,ARRHENIUS,RMSD,CUSTOM};

// ID type used by constraint
enum{ATOM,FRAG};

// keyword values that accept variables as input
enum{NEVERY,RMIN,RMAX,PROB,NRATE};

// flag for one-proc vs shared reaction sites
enum{LOCAL,GLOBAL};

// values for molecule_keyword
enum{OFF,INTER,INTRA};

/* ---------------------------------------------------------------------- */

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
  attempted_rxn = 0;
  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  global_freq = 1;
  extvector = 0;
  rxnID = 0;
  maxnconstraints = 0;
  narrhenius = 0;
  status = PROCEED;

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
  onemol_nxspecial = nullptr;
  twomol_nxspecial = nullptr;
  xspecial = nullptr;
  onemol_xspecial = nullptr;
  twomol_xspecial = nullptr;

  // these group names are reserved for use exclusively by bond/react
  master_group = (char *) "bond_react_MASTER_group";

  // by using fixed group names, only one instance of fix bond/react is allowed.
  if (modify->get_fix_by_style("^bond/react").size() != 0)
    error->all(FLERR,"Only one instance of fix bond/react allowed at a time");

  // let's find number of reactions specified
  nreacts = 0;
  for (int i = 3; i < narg; i++) {
    if (strcmp(arg[i],"react") == 0) {
      nreacts++;
      i = i + 6; // skip past mandatory arguments
      if (i > narg) error->all(FLERR,"Illegal fix bond/react command: "
                               "'react' has too few arguments");
    }
  }

  if (nreacts == 0) error->all(FLERR,"Illegal fix bond/react command: "
                               "missing mandatory 'react' argument");

  size_vector = nreacts;

  int iarg = 3;
  stabilization_flag = 0;
  reset_mol_ids_flag = 1;
  int num_common_keywords = 2;
  for (int m = 0; m < num_common_keywords; m++) {
    if (strcmp(arg[iarg],"stabilization") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'stabilization' keyword has too few arguments");
      stabilization_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      if (stabilization_flag) {
        if (iarg+4 > narg) error->all(FLERR,"Illegal fix bond/react command:"
                                      "'stabilization' keyword has too few arguments");
        exclude_group = utils::strdup(arg[iarg+2]);
        nve_limit_xmax = arg[iarg+3];
        iarg += 4;
      } else iarg += 2;
    } else if (strcmp(arg[iarg],"reset_mol_ids") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'reset_mol_ids' keyword has too few arguments");
      reset_mol_ids_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"react") == 0) {
      break;
    } else error->all(FLERR,"Illegal fix bond/react command: unknown keyword");
  }

  if (reset_mol_ids_flag) {
    delete reset_mol_ids;
    reset_mol_ids = new ResetAtomsMol(lmp);
    reset_mol_ids->create_computes(id,group->names[igroup]);
  }

  // set up common variables as vectors of length 'nreacts'
  // nevery, cutoff, onemol, twomol, superimpose file

  // this looks excessive
  // the price of vectorization (all reactions in one command)?
  memory->create(rxn_name,nreacts,MAXNAME,"bond/react:rxn_name");
  memory->create(nevery,nreacts,"bond/react:nevery");
  memory->create(cutsq,nreacts,2,"bond/react:cutsq");
  memory->create(unreacted_mol,nreacts,"bond/react:unreacted_mol");
  memory->create(reacted_mol,nreacts,"bond/react:reacted_mol");
  memory->create(fraction,nreacts,"bond/react:fraction");
  memory->create(max_rxn,nreacts,"bond/react:max_rxn");
  memory->create(nlocalskips,nreacts,"bond/react:nlocalskips");
  memory->create(nghostlyskips,nreacts,"bond/react:nghostlyskips");
  memory->create(seed,nreacts,"bond/react:seed");
  memory->create(limit_duration,nreacts,"bond/react:limit_duration");
  memory->create(rate_limit,3,nreacts,"bond/react:rate_limit");
  memory->create(stabilize_steps_flag,nreacts,"bond/react:stabilize_steps_flag");
  memory->create(custom_charges_fragid,nreacts,"bond/react:custom_charges_fragid");
  memory->create(rescale_charges_flag,nreacts,"bond/react:rescale_charges_flag");
  memory->create(create_atoms_flag,nreacts,"bond/react:create_atoms_flag");
  memory->create(modify_create_fragid,nreacts,"bond/react:modify_create_fragid");
  memory->create(overlapsq,nreacts,"bond/react:overlapsq");
  memory->create(molecule_keyword,nreacts,"bond/react:molecule_keyword");
  memory->create(nconstraints,nreacts,"bond/react:nconstraints");
  memory->create(constraintstr,nreacts,MAXLINE,"bond/react:constraintstr");
  memory->create(var_flag,NUMVARVALS,nreacts,"bond/react:var_flag");
  memory->create(var_id,NUMVARVALS,nreacts,"bond/react:var_id");
  memory->create(iatomtype,nreacts,"bond/react:iatomtype");
  memory->create(jatomtype,nreacts,"bond/react:jatomtype");
  memory->create(ibonding,nreacts,"bond/react:ibonding");
  memory->create(jbonding,nreacts,"bond/react:jbonding");
  memory->create(closeneigh,nreacts,"bond/react:closeneigh");
  memory->create(groupbits,nreacts,"bond/react:groupbits");
  memory->create(reaction_count,nreacts,"bond/react:reaction_count");
  memory->create(local_rxn_count,nreacts,"bond/react:local_rxn_count");
  memory->create(ghostly_rxn_count,nreacts,"bond/react:ghostly_rxn_count");
  memory->create(reaction_count_total,nreacts,"bond/react:reaction_count_total");

  for (int i = 0; i < nreacts; i++) {
    fraction[i] = 1;
    seed[i] = 12345;
    max_rxn[i] = INT_MAX;
    for (int j = 0; j < 3; j++)
      rate_limit[j][i] = 0;
    stabilize_steps_flag[i] = 0;
    custom_charges_fragid[i] = -1;
    rescale_charges_flag[i] = 0;
    create_atoms_flag[i] = 0;
    modify_create_fragid[i] = -1;
    overlapsq[i] = 0;
    molecule_keyword[i] = OFF;
    nconstraints[i] = 0;
    // set default limit duration to 60 timesteps
    limit_duration[i] = 60;
    reaction_count[i] = 0;
    local_rxn_count[i] = 0;
    ghostly_rxn_count[i] = 0;
    reaction_count_total[i] = 0;
    for (int j = 0; j < NUMVARVALS; j++) {
      var_flag[j][i] = 0;
      var_id[j][i] = 0;
    }
  }

  char **files;
  files = new char*[nreacts];

  for (int rxn = 0; rxn < nreacts; rxn++) {

    if (strcmp(arg[iarg],"react") != 0) error->all(FLERR,"Illegal fix bond/react command: "
                                                   "'react' or 'stabilization' has incorrect arguments");

    iarg++;

    int n = strlen(arg[iarg]) + 1;
    if (n > MAXNAME) error->all(FLERR,"Reaction name (react-ID) is too long (limit: 256 characters)");
    strcpy(rxn_name[rxn],arg[iarg++]);

    int groupid = group->find(arg[iarg++]);
    if (groupid == -1) error->all(FLERR,"Could not find fix group ID");
    groupbits[rxn] = group->bitmask[groupid];

    if (strncmp(arg[iarg],"v_",2) == 0) read_variable_keyword(&arg[iarg][2],NEVERY,rxn);
    else {
      nevery[rxn] = utils::inumeric(FLERR,arg[iarg],false,lmp);
      if (nevery[rxn] <= 0) error->all(FLERR,"Illegal fix bond/react command: "
                                       "'Nevery' must be a positive integer");
    }
    iarg++;

    double cutoff;
    if (strncmp(arg[iarg],"v_",2) == 0) {
      read_variable_keyword(&arg[iarg][2],RMIN,rxn);
      cutoff = input->variable->compute_equal(var_id[RMIN][rxn]);
    } else cutoff = utils::numeric(FLERR,arg[iarg],false,lmp);
      if (cutoff < 0.0) error->all(FLERR,"Illegal fix bond/react command: "
                                   "'Rmin' cannot be negative");
      cutsq[rxn][0] = cutoff*cutoff;
    iarg++;

    if (strncmp(arg[iarg],"v_",2) == 0) {
      read_variable_keyword(&arg[iarg][2],RMAX,rxn);
      cutoff = input->variable->compute_equal(var_id[RMAX][rxn]);
    } else cutoff = utils::numeric(FLERR,arg[iarg],false,lmp);
      if (cutoff < 0.0) error->all(FLERR,"Illegal fix bond/react command:"
                                   "'Rmax' cannot be negative");
      cutsq[rxn][1] = cutoff*cutoff;
    iarg++;

    unreacted_mol[rxn] = atom->find_molecule(arg[iarg++]);
    if (unreacted_mol[rxn] == -1) error->all(FLERR,"Unreacted molecule template ID for "
                                             "fix bond/react does not exist");
    reacted_mol[rxn] = atom->find_molecule(arg[iarg++]);
    if (reacted_mol[rxn] == -1) error->all(FLERR,"Reacted molecule template ID for "
                                           "fix bond/react does not exist");

    //read map file
    files[rxn] = utils::strdup(arg[iarg]);
    iarg++;

    while (iarg < narg && strcmp(arg[iarg],"react") != 0) {
      if (strcmp(arg[iarg],"prob") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'prob' keyword has too few arguments");
        // check if probability is a variable
        if (strncmp(arg[iarg+1],"v_",2) == 0) {
          read_variable_keyword(&arg[iarg+1][2],PROB,rxn);
          fraction[rxn] = input->variable->compute_equal(var_id[PROB][rxn]);
        } else {
          // otherwise probability should be a number
          fraction[rxn] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        }
        seed[rxn] = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
        if (fraction[rxn] < 0.0 || fraction[rxn] > 1.0)
          error->all(FLERR,"Illegal fix bond/react command: "
                     "probability fraction must between 0 and 1, inclusive");
        if (seed[rxn] <= 0) error->all(FLERR,"Illegal fix bond/react command: "
                                       "probability seed must be positive");
        iarg += 3;
      } else if (strcmp(arg[iarg],"max_rxn") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'max_rxn' has too few arguments");
        max_rxn[rxn] = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
        if (max_rxn[rxn] < 0) error->all(FLERR,"Illegal fix bond/react command: "
                                         "'max_rxn' cannot be negative");
        iarg += 2;
      } else if (strcmp(arg[iarg],"rate_limit") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'rate_limit' has too few arguments");
        rate_limit[0][rxn] = 1; // serves as flag for rate_limit keyword
        if (strncmp(arg[iarg+1],"v_",2) == 0) read_variable_keyword(&arg[iarg+1][2],NRATE,rxn);
        else rate_limit[1][rxn] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        rate_limit[2][rxn] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        iarg += 3;
      } else if (strcmp(arg[iarg],"stabilize_steps") == 0) {
        if (stabilization_flag == 0) error->all(FLERR,"Stabilize_steps keyword "
                                                "used without stabilization keyword");
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'stabilize_steps' has too few arguments");
        limit_duration[rxn] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        stabilize_steps_flag[rxn] = 1;
        iarg += 2;
      } else if (strcmp(arg[iarg],"custom_charges") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'custom_charges' has too few arguments");
        if (strcmp(arg[iarg+1],"no") == 0) custom_charges_fragid[rxn] = -1; //default
        else {
          custom_charges_fragid[rxn] = atom->molecules[unreacted_mol[rxn]]->findfragment(arg[iarg+1]);
          if (custom_charges_fragid[rxn] < 0) error->one(FLERR,"Fix bond/react: Molecule fragment for "
                                                         "'custom_charges' keyword does not exist");
        }
        iarg += 2;
      } else if (strcmp(arg[iarg],"rescale_charges") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'rescale_charges' has too few arguments");
        if (strcmp(arg[iarg+1],"no") == 0) rescale_charges_flag[rxn] = 0; //default
        else if (strcmp(arg[iarg+1],"yes") == 0) rescale_charges_flag[rxn] = 1;
        else error->one(FLERR,"Bond/react: Illegal option for 'rescale_charges' keyword");
        iarg += 2;
      } else if (strcmp(arg[iarg],"molecule") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'molecule' has too few arguments");
        if (strcmp(arg[iarg+1],"off") == 0) molecule_keyword[rxn] = OFF; //default
        else if (strcmp(arg[iarg+1],"inter") == 0) molecule_keyword[rxn] = INTER;
        else if (strcmp(arg[iarg+1],"intra") == 0) molecule_keyword[rxn] = INTRA;
        else error->one(FLERR,"Fix bond/react: Illegal option for 'molecule' keyword");
        iarg += 2;
      } else if (strcmp(arg[iarg],"modify_create") == 0) {
        if (iarg++ > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'modify_create' has too few arguments");
        while (iarg < narg && strcmp(arg[iarg],"react") != 0 ) {
          if (strcmp(arg[iarg],"fit") == 0) {
            if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                          "'modify_create' has too few arguments");
            if (strcmp(arg[iarg+1],"all") == 0) modify_create_fragid[rxn] = -1; //default
            else {
              modify_create_fragid[rxn] = atom->molecules[reacted_mol[rxn]]->findfragment(arg[iarg+1]);
              if (modify_create_fragid[rxn] < 0) error->one(FLERR,"Fix bond/react: Molecule fragment for "
                                                             "'modify_create' keyword does not exist");
            }
            iarg += 2;
          } else if (strcmp(arg[iarg],"overlap") == 0) {
            if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                          "'modify_create' has too few arguments");
            overlapsq[rxn] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
            overlapsq[rxn] *= overlapsq[rxn];
            iarg += 2;
          } else break;
        }
      } else error->all(FLERR,"Illegal fix bond/react command: unknown keyword");
    }
  }

  max_natoms = 0; // the number of atoms in largest molecule template
  max_rate_limit_steps = 0;
  for (int myrxn = 0; myrxn < nreacts; myrxn++) {
    twomol = atom->molecules[reacted_mol[myrxn]];
    max_natoms = MAX(max_natoms,twomol->natoms);
    max_rate_limit_steps = MAX(max_rate_limit_steps,rate_limit[2][myrxn]);
  }

  memory->create(equivalences,max_natoms,2,nreacts,"bond/react:equivalences");
  memory->create(reverse_equiv,max_natoms,2,nreacts,"bond/react:reverse_equiv");
  memory->create(edge,max_natoms,nreacts,"bond/react:edge");
  memory->create(landlocked_atoms,max_natoms,nreacts,"bond/react:landlocked_atoms");
  memory->create(store_rxn_count,max_rate_limit_steps,nreacts,"bond/react:store_rxn_count");
  memory->create(custom_charges,max_natoms,nreacts,"bond/react:custom_charges");
  memory->create(delete_atoms,max_natoms,nreacts,"bond/react:delete_atoms");
  memory->create(create_atoms,max_natoms,nreacts,"bond/react:create_atoms");
  memory->create(chiral_atoms,max_natoms,6,nreacts,"bond/react:chiral_atoms");

  for (int j = 0; j < nreacts; j++) {
    for (int i = 0; i < max_natoms; i++) {
      edge[i][j] = 0;
      custom_charges[i][j] = 1; // update all partial charges by default
      delete_atoms[i][j] = 0;
      create_atoms[i][j] = 0;
      for (int k = 0; k < 6; k++) {
        chiral_atoms[i][k][j] = 0;
      }
      // default equivalences to their own mol index
      // all but created atoms will be updated
      for (int m = 0; m < 2; m++) {
        equivalences[i][m][j] = i+1;
      }
    }
    for (int i = 0; i < max_rate_limit_steps; i++) {
      store_rxn_count[i][j] = -1;
    }
  }

  // read all map files afterward
  for (int i = 0; i < nreacts; i++) {
    open(files[i]);
    onemol = atom->molecules[unreacted_mol[i]];
    twomol = atom->molecules[reacted_mol[i]];
    onemol->check_attributes();
    twomol->check_attributes();
    get_molxspecials();
    read_map_file(i);
    fclose(fp);
    if (ncreate == 0 && onemol->natoms != twomol->natoms)
      error->all(FLERR,"Fix bond/react: Reaction templates must contain the same number of atoms");
    else if (ncreate > 0 && onemol->natoms + ncreate != twomol->natoms)
      error->all(FLERR,"Fix bond/react: Incorrect number of created atoms");
    iatomtype[i] = onemol->type[ibonding[i]-1];
    jatomtype[i] = onemol->type[jbonding[i]-1];
    find_landlocked_atoms(i);
    if (custom_charges_fragid[i] >= 0) CustomCharges(custom_charges_fragid[i],i);
  }

  // get the names of per-atom variables needed by 'rxn' functions of custom constraint
  customvarnames();

  // initialize Marsaglia RNG with processor-unique seed (Arrhenius prob)

  rrhandom = new RanMars*[narrhenius];
  int tmp = 0;
  for (int i = 0; i < nreacts; i++) {
    for (int j = 0; j < nconstraints[i]; j++) {
      if (constraints[j][i].type == ARRHENIUS) {
        rrhandom[tmp++] = new RanMars(lmp,(int) constraints[j][i].par[4] + comm->me);
      }
    }
  }

  for (int i = 0; i < nreacts; i++) {
    delete [] files[i];
  }
  delete [] files;

  if (atom->molecular != Atom::MOLECULAR)
    error->all(FLERR,"Fix bond/react: Cannot use fix bond/react with non-molecular systems");

  // check if bonding atoms are 1-2, 1-3, or 1-4 bonded neighbors
  // if so, we don't need non-bonded neighbor list
  for (int myrxn = 0; myrxn < nreacts; myrxn++) {
    closeneigh[myrxn] = -1; // indicates will search non-bonded neighbors
    onemol = atom->molecules[unreacted_mol[myrxn]];
    get_molxspecials();
    for (int k = 0; k < onemol_nxspecial[ibonding[myrxn]-1][2]; k++) {
      if (onemol_xspecial[ibonding[myrxn]-1][k] == jbonding[myrxn]) {
        closeneigh[myrxn] = 2; // index for 1-4 neighbor
        if (k < onemol_nxspecial[ibonding[myrxn]-1][1])
          closeneigh[myrxn] = 1; // index for 1-3 neighbor
        if (k < onemol_nxspecial[ibonding[myrxn]-1][0])
          closeneigh[myrxn] = 0; // index for 1-2 neighbor
        break;
      }
    }
  }

  // initialize Marsaglia RNG with processor-unique seed ('prob' keyword)

  random = new RanMars*[nreacts];
  for (int i = 0; i < nreacts; i++) {
    random[i] = new RanMars(lmp,seed[i] + comm->me);
  }

  // set comm sizes needed by this fix
  // forward is big due to comm of broken bonds and 1-2 neighbors

  comm_forward = MAX(2,2+atom->maxspecial);
  comm_reverse = 2;

  // allocate arrays local to this fix
  nmax = 0;
  partner = finalpartner = nullptr;
  distsq = nullptr;
  maxattempt = 0;
  attempt = nullptr;
  nattempt = nullptr;
  allnattempt = 0;
  my_num_mega = 0;
  local_num_mega = 0;
  ghostly_num_mega = 0;
  restore =  nullptr;

  // zero out stats
  global_megasize = 0;
  avail_guesses = 0;
  glove_counter = 0;
  guess_branch = new int[MAXGUESS]();
  pioneer_count = new int[max_natoms];
  my_mega_glove = nullptr;
  local_mega_glove = nullptr;
  ghostly_mega_glove = nullptr;
  global_mega_glove = nullptr;

  // these are merely loop indices that became important
  pion = neigh = trace = 0;

  id_fix1 = nullptr;
  id_fix2 = nullptr;
  id_fix3 = nullptr;
  statted_id = nullptr;
  custom_exclude_flag = 0;

  // used to store restart info
  set = new Set[nreacts];
  memset(set,0,nreacts*sizeof(Set));
}

/* ---------------------------------------------------------------------- */

FixBondReact::~FixBondReact()
{
  for (int i = 0; i < narrhenius; i++) {
    delete rrhandom[i];
  }
  delete [] rrhandom;

  for (int i = 0; i < nreacts; i++) {
    delete random[i];
  }
  delete [] random;

  delete reset_mol_ids;

  memory->destroy(partner);
  memory->destroy(finalpartner);
  memory->destroy(nattempt);
  memory->destroy(distsq);
  memory->destroy(attempt);
  memory->destroy(edge);
  memory->destroy(equivalences);
  memory->destroy(reverse_equiv);
  memory->destroy(landlocked_atoms);
  memory->destroy(store_rxn_count);
  memory->destroy(custom_charges);
  memory->destroy(delete_atoms);
  memory->destroy(create_atoms);
  memory->destroy(chiral_atoms);
  if (vvec != nullptr) memory->destroy(vvec);

  memory->destroy(rxn_name);
  memory->destroy(nevery);
  memory->destroy(cutsq);
  memory->destroy(unreacted_mol);
  memory->destroy(reacted_mol);
  memory->destroy(fraction);
  memory->destroy(seed);
  memory->destroy(max_rxn);
  memory->destroy(nlocalskips);
  memory->destroy(nghostlyskips);
  memory->destroy(limit_duration);
  memory->destroy(var_flag);
  memory->destroy(var_id);
  memory->destroy(rate_limit);
  memory->destroy(stabilize_steps_flag);
  memory->destroy(custom_charges_fragid);
  memory->destroy(rescale_charges_flag);
  memory->destroy(molecule_keyword);
  memory->destroy(nconstraints);
  memory->destroy(constraintstr);
  memory->destroy(create_atoms_flag);
  memory->destroy(modify_create_fragid);
  memory->destroy(overlapsq);

  memory->destroy(iatomtype);
  memory->destroy(jatomtype);
  memory->destroy(ibonding);
  memory->destroy(jbonding);
  memory->destroy(closeneigh);
  memory->destroy(groupbits);
  memory->destroy(reaction_count);
  memory->destroy(local_rxn_count);
  memory->destroy(ghostly_rxn_count);
  memory->destroy(reaction_count_total);

  if (newton_bond == 0) {
    memory->destroy(xspecial);
    memory->destroy(nxspecial);
    memory->destroy(onemol_xspecial);
    memory->destroy(onemol_nxspecial);
    memory->destroy(twomol_xspecial);
    memory->destroy(twomol_nxspecial);
  }

  if (attempted_rxn == 1) {
    memory->destroy(restore_pt);
    memory->destroy(restore);
    memory->destroy(glove);
    memory->destroy(pioneers);
    memory->destroy(my_mega_glove);
    memory->destroy(local_mega_glove);
    memory->destroy(ghostly_mega_glove);
  }

  memory->destroy(global_mega_glove);

  if (stabilization_flag == 1) {
    // delete fixes if not already deleted
    if (id_fix1 && modify->get_fix_by_id(id_fix1)) modify->delete_fix(id_fix1);
    delete[] id_fix1;

    if (id_fix3 && modify->get_fix_by_id(id_fix3)) modify->delete_fix(id_fix3);
    delete[] id_fix3;
  }

  if (id_fix2 && modify->get_fix_by_id(id_fix2)) modify->delete_fix(id_fix2);
  delete[] id_fix2;

  delete[] statted_id;
  delete[] guess_branch;
  delete[] pioneer_count;
  delete[] set;

  if (group) {
    group->assign(std::string(master_group) + " delete");
    if (stabilization_flag == 1) {
      group->assign(std::string(exclude_group) + " delete");
      delete[] exclude_group;
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixBondReact::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
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
  id_fix2 = utils::strdup("bond_react_props_internal");
  if (!modify->get_fix_by_id(id_fix2))
    fix2 = modify->add_fix(std::string(id_fix2) +
                           " all property/atom i_limit_tags i_react_tags ghost yes");

  // create master_group if not already existing
  // NOTE: limit_tags and react_tags automaticaly intitialized to zero (unless read from restart)
  group->find_or_create(master_group);
  std::string cmd = fmt::format("{} dynamic all property limit_tags",master_group);
  group->assign(cmd);

  if (stabilization_flag == 1) {
    int groupid = group->find(exclude_group);
    // create exclude_group if not already existing, or use as parent group if static
    if (groupid == -1 || group->dynamic[groupid] == 0) {

      // create stabilization per-atom property
      id_fix3 = utils::strdup("bond_react_stabilization_internal");
      if (!modify->get_fix_by_id(id_fix3))
        fix3 = modify->add_fix(std::string(id_fix3) +
                               " all property/atom i_statted_tags ghost yes");

      statted_id = utils::strdup("statted_tags");

      // if static group exists, use as parent group
      // also, rename dynamic exclude_group by appending '_REACT'
      char *exclude_PARENT_group;
      exclude_PARENT_group = utils::strdup(exclude_group);
      delete[] exclude_group;
      exclude_group = utils::strdup(std::string(exclude_PARENT_group)+"_REACT");

      group->find_or_create(exclude_group);
      if (groupid == -1)
        cmd = fmt::format("{} dynamic all property statted_tags",exclude_group);
      else
        cmd = fmt::format("{} dynamic {} property statted_tags",exclude_group,exclude_PARENT_group);
      group->assign(cmd);
      delete[] exclude_PARENT_group;

      // on to statted_tags (system-wide thermostat)
      // initialize per-atom statted_flags to 1
      // (only if not already initialized by restart)
      if (fix3->restart_reset != 1) {
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
      Fix *fix = modify->get_fix_by_id(std::string("GROUP_")+exclude_group);

      // this returns names of corresponding property
      int unused;
      char *idprop;
      idprop = (char *) fix->extract("property",unused);
      if (idprop == nullptr)
        error->all(FLERR,"Exclude group must be a per-atom property group");
      statted_id = utils::strdup(idprop);

      // initialize per-atom statted_tags to 1
      // need to correct for smooth restarts
      //int flag,cols;
      //int index = atom->find_custom(statted_id,flag,cols);
      //int *i_statted_tags = atom->ivector[index];
      //for (int i = 0; i < atom->nlocal; i++)
      //  i_statted_tags[i] = 1;
    }

    // let's create a new nve/limit fix to limit newly reacted atoms
    id_fix1 = utils::strdup("bond_react_MASTER_nve_limit");
    if (!modify->get_fix_by_id(id_fix1))
      fix1 = modify->add_fix(fmt::format("{} {} nve/limit  {}",
                                         id_fix1,master_group,nve_limit_xmax));
  }
}

/* ---------------------------------------------------------------------- */

void FixBondReact::init()
{

  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;

  // check cutoff for iatomtype,jatomtype
  for (int i = 0; i < nreacts; i++) {
    if (!utils::strmatch(force->pair_style,"^hybrid"))
      if (force->pair == nullptr || cutsq[i][1] > force->pair->cutsq[iatomtype[i]][jatomtype[i]])
        error->all(FLERR,"Fix bond/react: Fix bond/react cutoff is longer than pairwise cutoff");
  }

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
  This function is modified from LAMMPS’ fix bond/create.
---------------------------------------------------------------------- */

void FixBondReact::post_integrate()
{
  // update store_rxn_count on every step
  for (int myrxn = 0; myrxn < nreacts; myrxn++) {
    if (rate_limit[0][myrxn] == 1) {
      for (int i = rate_limit[2][myrxn]-1; i > 0; i--) {
        store_rxn_count[i][myrxn] = store_rxn_count[i-1][myrxn];
      }
      store_rxn_count[0][myrxn] = reaction_count_total[myrxn];
    }
  }

  // check if any reactions could occur on this timestep
  int nevery_check = 1;
  for (int i = 0; i < nreacts; i++) {
    if (var_flag[NEVERY][i])
      nevery[i] = ceil(input->variable->compute_equal(var_id[NEVERY][i]));
    if (nevery[i] <= 0)
      error->all(FLERR,"Illegal fix bond/react command: "
                 "'Nevery' must be a positive integer");
    if (!(update->ntimestep % nevery[i])) {
      nevery_check = 0;
      break;
    }
  }

  for (int i = 0; i < nreacts; i++) {
    reaction_count[i] = 0;
    local_rxn_count[i] = 0;
    ghostly_rxn_count[i] = 0;
    nlocalskips[i] = 0;
    nghostlyskips[i] = 0;
    // update reaction probability
    if (var_flag[PROB][i])
      fraction[i] = input->variable->compute_equal(var_id[PROB][i]);
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
    memory->destroy(nattempt);
    nmax = atom->nmax;
    memory->create(partner,nmax,"bond/react:partner");
    memory->create(finalpartner,nmax,"bond/react:finalpartner");
    memory->create(distsq,nmax,2,"bond/react:distsq");
    memory->create(nattempt,nreacts,"bond/react:nattempt");
  }

  // reset 'attempt' counts
  for (int i = 0; i < nreacts; i++) {
    nattempt[i] = 0;
  }
  // reset per-bond compute map flag
  atoms2bondflag = 0;

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  // loop over neighbors of my atoms
  // each atom sets one closest eligible partner atom ID to bond with

  tagint *tag = atom->tag;
  int *type = atom->type;

  neighbor->build_one(list,1);

  // here we define a full special list, independent of Newton setting
  if (newton_bond == 1) {
    nxspecial = atom->nspecial;
    xspecial = atom->special;
  } else {
    int nall = atom->nlocal + atom->nghost;
    memory->destroy(nxspecial);
    memory->destroy(xspecial);
    memory->create(nxspecial,nall,3,"bond/react:nxspecial");
    memory->create(xspecial,nall,atom->maxspecial,"bond/react:xspecial");
    for (int i = 0; i < atom->nlocal; i++) {
      nxspecial[i][0] = atom->num_bond[i];
      for (int j = 0; j < nxspecial[i][0]; j++) {
        xspecial[i][j] = atom->bond_atom[i][j];
      }
      nxspecial[i][1] = atom->nspecial[i][1];
      nxspecial[i][2] = atom->nspecial[i][2];
      int joffset = nxspecial[i][0] - atom->nspecial[i][0];
      for (int j = nxspecial[i][0]; j < nxspecial[i][2]; j++) {
        xspecial[i][j+joffset] = atom->special[i][j];
      }
    }
  }

  int j;
  for (rxnID = 0; rxnID < nreacts; rxnID++) {
    int rate_limit_flag = 1;
    if (rate_limit[0][rxnID] == 1) {
      int myrxn_count = store_rxn_count[rate_limit[2][rxnID]-1][rxnID];
      if (myrxn_count == -1) rate_limit_flag = 0;
      else {
        int nrxns_delta = reaction_count_total[rxnID] - myrxn_count;
        int my_nrate;
        if (var_flag[NRATE][rxnID] == 1) {
          my_nrate = input->variable->compute_equal(var_id[NRATE][rxnID]);
        } else my_nrate = rate_limit[1][rxnID];
        if (nrxns_delta >= my_nrate) rate_limit_flag = 0;
      }
    }
    if ((update->ntimestep % nevery[rxnID]) ||
        (max_rxn[rxnID] <= reaction_count_total[rxnID]) ||
        (rate_limit_flag == 0)) continue;
    for (int ii = 0; ii < nall; ii++) {
      partner[ii] = 0;
      finalpartner[ii] = 0;
      distsq[ii][0] = 0.0;
      distsq[ii][1] = BIG;
    }

    // fork between far and close_partner here
    if (closeneigh[rxnID] < 0) {
      far_partner();
      // reverse comm of distsq and partner
      // not needed if newton_pair off since I,J pair was seen by both procs
      commflag = 2;
      if (force->newton_pair) comm->reverse_comm(this);
    } else {
      close_partner();
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
        if (nattempt[rxnID] > maxattempt-2) {
          maxattempt += DELTA;
          // third dim of 'attempt': bond/react integer ID
          memory->grow(attempt,maxattempt,2,nreacts,"bond/react:attempt");
        }
        // to ensure types remain in same order
        if (iatomtype[rxnID] == type[i]) {
          attempt[nattempt[rxnID]][0][rxnID] = tag[i];
          attempt[nattempt[rxnID]][1][rxnID] = finalpartner[i];
          nattempt[rxnID]++;
          // add another attempt if initiator atoms are same type
          if (iatomtype[rxnID] == jatomtype[rxnID]) {
            attempt[nattempt[rxnID]][0][rxnID] = finalpartner[i];
            attempt[nattempt[rxnID]][1][rxnID] = tag[i];
            nattempt[rxnID]++;
          }
        } else {
          attempt[nattempt[rxnID]][0][rxnID] = finalpartner[i];
          attempt[nattempt[rxnID]][1][rxnID] = tag[i];
          nattempt[rxnID]++;
        }
      }
    }
  }

  // break loop if no even eligible bonding atoms were found (on any proc)
  int some_chance;

  allnattempt = 0;
  for (int i = 0; i < nreacts; i++)
    allnattempt += nattempt[i];

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

void FixBondReact::far_partner()
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
    if (!(mask[i] & groupbits[rxnID])) continue;
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
      if (!(mask[j] & groupbits[rxnID])) {
        continue;
      }

      if (i_limit_tags[j] != 0) {
        continue;
      }

      if (molecule_keyword[rxnID] == INTER) {
        if (atom->molecule[i] == atom->molecule[j]) continue;
      } else if (molecule_keyword[rxnID] == INTRA) {
        if (atom->molecule[i] != atom->molecule[j]) continue;
      }

      jtype = type[j];
      possible = 0;
      if (itype == iatomtype[rxnID] && jtype == jatomtype[rxnID]) {
        possible = 1;
      } else if (itype == jatomtype[rxnID] && jtype == iatomtype[rxnID]) {
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
      domain->minimum_image(delx,dely,delz); // ghost location fix
      rsq = delx*delx + dely*dely + delz*delz;

      if (var_flag[RMIN][rxnID]) {
        double cutoff = input->variable->compute_equal(var_id[RMIN][rxnID]);
        cutsq[rxnID][0] = cutoff*cutoff;
      }
      if (var_flag[RMAX][rxnID]) {
        double cutoff = input->variable->compute_equal(var_id[RMAX][rxnID]);
        cutsq[rxnID][1] = cutoff*cutoff;
      }
      if (rsq >= cutsq[rxnID][1] || rsq <= cutsq[rxnID][0]) {
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

void FixBondReact::close_partner()
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
    if (closeneigh[rxnID] != 0)
      n = nxspecial[ii][closeneigh[rxnID]-1];
    for (; n < nxspecial[ii][closeneigh[rxnID]]; n++) {
      i1 = ii;
      i2 = atom->map(xspecial[ii][n]);
      jtype = type[i2];
      if (!(mask[i1] & groupbits[rxnID])) continue;
      if (!(mask[i2] & groupbits[rxnID])) continue;
      if (i_limit_tags[i1] != 0) continue;
      if (i_limit_tags[i2] != 0) continue;
      if (itype != iatomtype[rxnID] || jtype != jatomtype[rxnID]) continue;

      if (molecule_keyword[rxnID] == INTER) {
        if (atom->molecule[i1] == atom->molecule[i2]) continue;
      } else if (molecule_keyword[rxnID] == INTRA) {
        if (atom->molecule[i1] != atom->molecule[i2]) continue;
      }

      delx = x[i1][0] - x[i2][0];
      dely = x[i1][1] - x[i2][1];
      delz = x[i1][2] - x[i2][2];
      domain->minimum_image(delx,dely,delz); // ghost location fix
      rsq = delx*delx + dely*dely + delz*delz;

      if (var_flag[RMIN][rxnID]) {
        double cutoff = input->variable->compute_equal(var_id[RMIN][rxnID]);
        cutsq[rxnID][0] = cutoff*cutoff;
      }
      if (var_flag[RMAX][rxnID]) {
        double cutoff = input->variable->compute_equal(var_id[RMAX][rxnID]);
        cutsq[rxnID][1] = cutoff*cutoff;
      }
      if (rsq >= cutsq[rxnID][1] || rsq <= cutsq[rxnID][0]) continue;

      if (closeneigh[rxnID] == 0) {
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
  // 'pion' is the pioneer loop index
  // 'neigh' in the first neighbor index
  // 'trace' retraces the first nieghbors
  // trace: once you choose a first neighbor, you then check for other nieghbors of same type

  if (attempted_rxn == 1) {
    memory->destroy(restore_pt);
    memory->destroy(restore);
    memory->destroy(glove);
    memory->destroy(pioneers);
    memory->destroy(my_mega_glove);
    memory->destroy(local_mega_glove);
    memory->destroy(ghostly_mega_glove);
  }

  memory->create(glove,max_natoms,2,"bond/react:glove");
  memory->create(restore_pt,MAXGUESS,4,"bond/react:restore_pt");
  memory->create(pioneers,max_natoms,"bond/react:pioneers");
  memory->create(restore,max_natoms,MAXGUESS*4,"bond/react:restore");
  memory->create(my_mega_glove,max_natoms+1,allnattempt,"bond/react:local_mega_glove");

  for (int i = 0; i < max_natoms+1; i++)
    for (int j = 0; j < allnattempt; j++)
      my_mega_glove[i][j] = 0;

  attempted_rxn = 1;

  // let's finally begin the superimpose loop
  for (rxnID = 0; rxnID < nreacts; rxnID++) {
    for (lcl_inst = 0; lcl_inst < nattempt[rxnID]; lcl_inst++) {

      onemol = atom->molecules[unreacted_mol[rxnID]];
      twomol = atom->molecules[reacted_mol[rxnID]];
      get_molxspecials();

      status = PROCEED;

      glove_counter = 0;
      for (int i = 0; i < max_natoms; i++) {
        for (int j = 0; j < 2; j++) {
          glove[i][j] = 0;
        }
      }

      for (int i = 0; i < MAXGUESS; i++) {
        guess_branch[i] = 0;
      }

      int myibonding = ibonding[rxnID];
      int myjbonding = jbonding[rxnID];

      glove[myibonding-1][0] = myibonding;
      glove[myibonding-1][1] = attempt[lcl_inst][0][rxnID];
      glove_counter++;
      glove[myjbonding-1][0] = myjbonding;
      glove[myjbonding-1][1] = attempt[lcl_inst][1][rxnID];
      glove_counter++;

      // special case, only two atoms in reaction templates
      // then: bonding onemol_nxspecials guaranteed to be equal, and either 0 or 1
      if (glove_counter == onemol->natoms) {
        tagint local_atom1 = atom->map(glove[myibonding-1][1]);
        tagint local_atom2 = atom->map(glove[myjbonding-1][1]);
        if ( (nxspecial[local_atom1][0] == onemol_nxspecial[myibonding-1][0] &&
              nxspecial[local_atom2][0] == nxspecial[local_atom1][0]) &&
             (nxspecial[local_atom1][0] == 0 ||
              xspecial[local_atom1][0] == atom->tag[local_atom2]) &&
             check_constraints()) {
          if (fraction[rxnID] < 1.0 &&
              random[rxnID]->uniform() >= fraction[rxnID]) {
            status = REJECT;
          } else {
            status = ACCEPT;
            my_mega_glove[0][my_num_mega] = rxnID;
            for (int i = 0; i < onemol->natoms; i++) {
              my_mega_glove[i+1][my_num_mega] = glove[i][1];
            }
            my_num_mega++;
          }
        } else status = REJECT;
      }

      avail_guesses = 0;

      for (int i = 0; i < max_natoms; i++)
        pioneer_count[i] = 0;

      for (int i = 0; i < onemol_nxspecial[myibonding-1][0]; i++)
        pioneer_count[onemol_xspecial[myibonding-1][i]-1]++;

      for (int i = 0; i < onemol_nxspecial[myjbonding-1][0]; i++)
        pioneer_count[onemol_xspecial[myjbonding-1][i]-1]++;


      int hang_catch = 0;
      while (!(status == ACCEPT || status == REJECT)) {

        for (int i = 0; i < max_natoms; i++) {
          pioneers[i] = 0;
        }

        for (int i = 0; i < onemol->natoms; i++) {
          if (glove[i][0] !=0 && pioneer_count[i] < onemol_nxspecial[i][0] && edge[i][rxnID] == 0) {
            pioneers[i] = 1;
          }
        }

        // run through the pioneers
        // due to use of restore points, 'pion' index can change in loop
        for (pion = 0; pion < onemol->natoms; pion++) {
          if (pioneers[pion] || status == GUESSFAIL) {
            make_a_guess();
            if (status == ACCEPT || status == REJECT) break;
          }
        }

        // reaction site found successfully!
        if (status == ACCEPT) {
          if (fraction[rxnID] < 1.0 &&
              random[rxnID]->uniform() >= fraction[rxnID]) status = REJECT;
          else {
            my_mega_glove[0][my_num_mega] = rxnID;
            for (int i = 0; i < onemol->natoms; i++) {
              my_mega_glove[i+1][my_num_mega] = glove[i][1];
            }
            my_num_mega++;
          }
        }
        hang_catch++;
        // let's go ahead and catch the simplest of hangs
        //if (hang_catch > onemol->natoms*4)
        if (hang_catch > atom->nlocal*30) {
          error->one(FLERR,"Fix bond/react: Excessive iteration of superimpose algorithm. "
              "Please check that all pre-reaction template atoms are linked to an initiator atom, "
              "via at least one path that does not involve edge atoms.");
        }
      }
    }
  }

  global_megasize = 0;

  memory->create(local_mega_glove,max_natoms+1,my_num_mega,"bond/react:local_mega_glove");
  memory->create(ghostly_mega_glove,max_natoms+1,my_num_mega,"bond/react:ghostly_mega_glove");

  for (int i = 0; i < max_natoms+1; i++) {
    for (int j = 0; j < my_num_mega; j++) {
      local_mega_glove[i][j] = 0;
      ghostly_mega_glove[i][j] = 0;
    }
  }

  dedup_mega_gloves(LOCAL); // make sure atoms aren't added to more than one reaction
  glove_ghostcheck(); // split into 'local' and 'global'
  ghost_glovecast(); // consolidate all mega_gloves to all processors

  MPI_Allreduce(&local_rxn_count[0],&reaction_count[0],nreacts,MPI_INT,MPI_SUM,world);

  int rxnflag = 0;
  if (comm->me == 0)
    for (int i = 0; i < nreacts; i++) {
      reaction_count_total[i] += reaction_count[i] + ghostly_rxn_count[i];
      rxnflag += reaction_count[i] + ghostly_rxn_count[i];
    }

  MPI_Bcast(&reaction_count_total[0], nreacts, MPI_INT, 0, world);
  MPI_Bcast(&rxnflag, 1, MPI_INT, 0, world);

  if (!rxnflag) return;

  // C++11 and later compatible version of Park pRNG
  std::random_device rnd;
  std::minstd_rand park_rng(rnd());

  // check if we overstepped our reaction limit, via either max_rxn or rate_limit
  for (int i = 0; i < nreacts; i++) {
    int overstep = 0;
    int max_rxn_overstep = reaction_count_total[i] - max_rxn[i];
    overstep = MAX(overstep,max_rxn_overstep);
    if (rate_limit[0][i] == 1) {
      int myrxn_count = store_rxn_count[rate_limit[2][i]-1][i];
      if (myrxn_count != -1) {
        int nrxn_delta = reaction_count_total[i] - myrxn_count;
        int my_nrate;
        if (var_flag[NRATE][i] == 1) {
          my_nrate = input->variable->compute_equal(var_id[NRATE][i]);
        } else my_nrate = rate_limit[1][i];
        int rate_limit_overstep = nrxn_delta - my_nrate;
        overstep = MAX(overstep,rate_limit_overstep);
      }
    }

    if (overstep > 0) {
      // let's randomly choose rxns to skip, unbiasedly from local and ghostly
      int *local_rxncounts;
      int *all_localskips;
      memory->create(local_rxncounts,nprocs,"bond/react:local_rxncounts");
      memory->create(all_localskips,nprocs,"bond/react:all_localskips");
      MPI_Gather(&local_rxn_count[i],1,MPI_INT,local_rxncounts,1,MPI_INT,0,world);
      if (comm->me == 0) {
        int delta_rxn = reaction_count[i] + ghostly_rxn_count[i];
        // when using variable input for rate_limit, rate_limit_overstep could be > delta_rxn (below)
        // we need to limit overstep to the number of reactions on this timestep
        // essentially skipping all reactions, would be more efficient to use a skip_all flag
        if (overstep > delta_rxn) overstep = delta_rxn;
        int *rxn_by_proc;
        memory->create(rxn_by_proc,delta_rxn,"bond/react:rxn_by_proc");
        for (int j = 0; j < delta_rxn; j++)
          rxn_by_proc[j] = -1; // corresponds to ghostly
        int itemp = 0;
        for (int j = 0; j < nprocs; j++)
          for (int k = 0; k < local_rxncounts[j]; k++)
            rxn_by_proc[itemp++] = j;
        std::shuffle(&rxn_by_proc[0],&rxn_by_proc[delta_rxn], park_rng);
        for (int j = 0; j < nprocs; j++)
          all_localskips[j] = 0;
        nghostlyskips[i] = 0;
        for (int j = 0; j < overstep; j++) {
          if (rxn_by_proc[j] == -1) nghostlyskips[i]++;
          else all_localskips[rxn_by_proc[j]]++;
        }
        memory->destroy(rxn_by_proc);
        reaction_count_total[i] -= overstep;
      }
      MPI_Scatter(&all_localskips[0],1,MPI_INT,&nlocalskips[i],1,MPI_INT,0,world);
      MPI_Bcast(&nghostlyskips[i],1,MPI_INT,0,world);
      memory->destroy(local_rxncounts);
      memory->destroy(all_localskips);
    }
  }
  MPI_Bcast(&reaction_count_total[0], nreacts, MPI_INT, 0, world);

  // this updates topology next step
  next_reneighbor = update->ntimestep;

  update_everything(); // change topology
}

/* ----------------------------------------------------------------------
  Screen for obvious algorithm fails. This is the return point when a guess
  has failed: check for available restore points.
------------------------------------------------------------------------- */

void FixBondReact::make_a_guess()
{
  int *type = atom->type;
  int nfirst_neighs = onemol_nxspecial[pion][0];

  // per-atom property indicating if in bond/react master group
  int flag,cols;
  int index1 = atom->find_custom("limit_tags",flag,cols);
  int *i_limit_tags = atom->ivector[index1];

  if (status == GUESSFAIL && avail_guesses == 0) {
    status = REJECT;
    return;
  }

  if (status == GUESSFAIL && avail_guesses > 0) {
    // load restore point
    for (int i = 0; i < onemol->natoms; i++) {
      glove[i][0] = restore[i][(avail_guesses*4)-4];
      glove[i][1] = restore[i][(avail_guesses*4)-3];
      pioneer_count[i] = restore[i][(avail_guesses*4)-2];
      pioneers[i] = restore[i][(avail_guesses*4)-1];
    }
    pion = restore_pt[avail_guesses-1][0];
    neigh = restore_pt[avail_guesses-1][1];
    trace = restore_pt[avail_guesses-1][2];
    glove_counter = restore_pt[avail_guesses-1][3];
    status = RESTORE;
    neighbor_loop();
    if (status != PROCEED) return;
  }

  nfirst_neighs = onemol_nxspecial[pion][0];

  //  check if any of first neighbors are in bond_react_MASTER_group
  //  if so, this constitutes a fail
  //  because still undergoing a previous reaction!
  //  could technically fail unnecessarily during a wrong guess if near edge atoms
  //  we accept this temporary and infrequent decrease in reaction occurrences

  for (int i = 0; i < nxspecial[atom->map(glove[pion][1])][0]; i++) {
    if (atom->map(xspecial[atom->map(glove[pion][1])][i]) < 0) {
      error->one(FLERR,"Fix bond/react: Fix bond/react needs ghost atoms from further away"); // parallel issues.
    }
    if (i_limit_tags[(int)atom->map(xspecial[atom->map(glove[pion][1])][i])] != 0) {
      status = GUESSFAIL;
      return;
    }
  }

  // check for same number of neighbors between unreacted mol and simulation
  if (nfirst_neighs != nxspecial[atom->map(glove[pion][1])][0]) {
    status = GUESSFAIL;
    return;
  }

  // make sure all neighbors aren't already assigned
  // an issue discovered for coarse-grained example
  int assigned_count = 0;
  for (int i = 0; i < nfirst_neighs; i++)
    for (int j = 0; j < onemol->natoms; j++)
      if (xspecial[atom->map(glove[pion][1])][i] == glove[j][1]) {
        assigned_count++;
        break;
      }

  if (assigned_count == nfirst_neighs) status = GUESSFAIL;

  // check if all neigh atom types are the same between simulation and unreacted mol
  int *mol_ntypes = new int[atom->ntypes];
  int *lcl_ntypes = new int[atom->ntypes];

  for (int i = 0; i < atom->ntypes; i++) {
    mol_ntypes[i] = 0;
    lcl_ntypes[i] = 0;
  }

  for (int i = 0; i < nfirst_neighs; i++) {
    mol_ntypes[(int)onemol->type[(int)onemol_xspecial[pion][i]-1]-1]++;
    lcl_ntypes[(int)type[(int)atom->map(xspecial[atom->map(glove[pion][1])][i])]-1]++; //added -1
  }

  for (int i = 0; i < atom->ntypes; i++) {
    if (mol_ntypes[i] != lcl_ntypes[i]) {
      status = GUESSFAIL;
      delete [] mol_ntypes;
      delete [] lcl_ntypes;
      return;
    }
  }

  delete [] mol_ntypes;
  delete [] lcl_ntypes;

  // okay everything seems to be in order. let's assign some ID pairs!!!
  neighbor_loop();
}

/* ----------------------------------------------------------------------
  Loop through all First Bonded Neighbors of the current Pioneer.
  Prepare appropriately if we are in Restore Mode.
------------------------------------------------------------------------- */

void FixBondReact::neighbor_loop()
{
  int nfirst_neighs = onemol_nxspecial[pion][0];

  if (status == RESTORE) {
    check_a_neighbor();
    return;
  }

  for (neigh = 0; neigh < nfirst_neighs; neigh++) {
    if (glove[(int)onemol_xspecial[pion][neigh]-1][0] == 0) {
      check_a_neighbor();
    }
  }
  // status should still = PROCEED
}

/* ----------------------------------------------------------------------
  Check if we can assign this First Neighbor to pre-reacted template
  without guessing. If so, do it! If not, call crosscheck_the_nieghbor().
------------------------------------------------------------------------- */

void FixBondReact::check_a_neighbor()
{
  int *type = atom->type;
  int nfirst_neighs = onemol_nxspecial[pion][0];

  if (status != RESTORE) {
    // special consideration for hydrogen atoms (and all first neighbors bonded to no other atoms) (and aren't edge atoms)
    if (onemol_nxspecial[(int)onemol_xspecial[pion][neigh]-1][0] == 1 && edge[(int)onemol_xspecial[pion][neigh]-1][rxnID] == 0) {

      for (int i = 0; i < nfirst_neighs; i++) {

        if (type[(int)atom->map(xspecial[(int)atom->map(glove[pion][1])][i])] == onemol->type[(int)onemol_xspecial[pion][neigh]-1] &&
            nxspecial[(int)atom->map(xspecial[(int)atom->map(glove[pion][1])][i])][0] == 1) {

          int already_assigned = 0;
          for (int j = 0; j < onemol->natoms; j++) {
            if (glove[j][1] == xspecial[atom->map(glove[pion][1])][i]) {
              already_assigned = 1;
              break;
            }
          }

          if (already_assigned == 0) {
            glove[(int)onemol_xspecial[pion][neigh]-1][0] = onemol_xspecial[pion][neigh];
            glove[(int)onemol_xspecial[pion][neigh]-1][1] = xspecial[(int)atom->map(glove[pion][1])][i];

            //another check for ghost atoms. perhaps remove the one in make_a_guess
            if (atom->map(glove[(int)onemol_xspecial[pion][neigh]-1][1]) < 0) {
              error->one(FLERR,"Fix bond/react: Fix bond/react needs ghost atoms from further away");
            }

            for (int j = 0; j < onemol_nxspecial[onemol_xspecial[pion][neigh]-1][0]; j++) {
              pioneer_count[onemol_xspecial[onemol_xspecial[pion][neigh]-1][j]-1]++;
            }

            glove_counter++;
            if (glove_counter == onemol->natoms) {
              if (ring_check() && check_constraints()) status = ACCEPT;
              else status = GUESSFAIL;
              return;
            }
            // status should still == PROCEED
            return;
          }
        }
      }
      // we are here if no matching atom found
      status = GUESSFAIL;
      return;
    }
  }

  crosscheck_the_neighbor();
  if (status != PROCEED) {
    if (status == CONTINUE)
      status = PROCEED;
    return;
  }

  // finally ready to match non-duplicate, non-edge atom IDs!!

  for (int i = 0; i < nfirst_neighs; i++) {

    if (type[atom->map((int)xspecial[(int)atom->map(glove[pion][1])][i])] == onemol->type[(int)onemol_xspecial[pion][neigh]-1]) {
      int already_assigned = 0;

      //check if a first neighbor of the pioneer is already assigned to pre-reacted template
      for (int j = 0; j < onemol->natoms; j++) {
        if (glove[j][1] == xspecial[atom->map(glove[pion][1])][i]) {
          already_assigned = 1;
          break;
        }
      }

      if (already_assigned == 0) {
        glove[(int)onemol_xspecial[pion][neigh]-1][0] = onemol_xspecial[pion][neigh];
        glove[(int)onemol_xspecial[pion][neigh]-1][1] = xspecial[(int)atom->map(glove[pion][1])][i];

        //another check for ghost atoms. perhaps remove the one in make_a_guess
        if (atom->map(glove[(int)onemol_xspecial[pion][neigh]-1][1]) < 0) {
          error->one(FLERR,"Fix bond/react: Fix bond/react needs ghost atoms from further away");
        }

        for (int ii = 0; ii < onemol_nxspecial[onemol_xspecial[pion][neigh]-1][0]; ii++) {
          pioneer_count[onemol_xspecial[onemol_xspecial[pion][neigh]-1][ii]-1]++;
        }

        glove_counter++;
        if (glove_counter == onemol->natoms) {
          if (ring_check() && check_constraints()) status = ACCEPT;
          else status = GUESSFAIL;
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

void FixBondReact::crosscheck_the_neighbor()
{
  int nfirst_neighs = onemol_nxspecial[pion][0];

  if (status == RESTORE) {
    inner_crosscheck_loop();
    return;
  }

  for (trace = 0; trace < nfirst_neighs; trace++) {
    if (neigh!=trace && onemol->type[(int)onemol_xspecial[pion][neigh]-1] == onemol->type[(int)onemol_xspecial[pion][trace]-1] &&
        glove[onemol_xspecial[pion][trace]-1][0] == 0) {

      if (avail_guesses == MAXGUESS) {
        error->warning(FLERR,"Fix bond/react: Fix bond/react failed because MAXGUESS set too small. ask developer for info");
        status = GUESSFAIL;
        return;
      }
      avail_guesses++;
      for (int i = 0; i < onemol->natoms; i++) {
        restore[i][(avail_guesses*4)-4] = glove[i][0];
        restore[i][(avail_guesses*4)-3] = glove[i][1];
        restore[i][(avail_guesses*4)-2] = pioneer_count[i];
        restore[i][(avail_guesses*4)-1] = pioneers[i];
        restore_pt[avail_guesses-1][0] = pion;
        restore_pt[avail_guesses-1][1] = neigh;
        restore_pt[avail_guesses-1][2] = trace;
        restore_pt[avail_guesses-1][3] = glove_counter;
      }

      inner_crosscheck_loop();
      return;
    }
  }
  // status is still 'PROCEED' if we are here!
}

/* ----------------------------------------------------------------------
  We are ready to make a guess. If there are multiple possible choices
  for this guess, keep track of these.
------------------------------------------------------------------------- */

void FixBondReact::inner_crosscheck_loop()
{
  int *type = atom->type;
  // arbitrarily limited to 5 identical first neighbors
  tagint tag_choices[5];
  int nfirst_neighs = onemol_nxspecial[pion][0];

  int num_choices = 0;
  for (int i = 0; i < nfirst_neighs; i++) {
    if (type[(int)atom->map(xspecial[atom->map(glove[pion][1])][i])] == onemol->type[(int)onemol_xspecial[pion][neigh]-1]) {
      if (num_choices == 5) { // here failed because too many identical first neighbors. but really no limit if situation arises
        status = GUESSFAIL;
        return;
      }
      tag_choices[num_choices++] = xspecial[atom->map(glove[pion][1])][i];
    }
  }

  // guess branch is for when multiple identical neighbors. then we guess each one in turn
  // guess_branch must work even when avail_guesses = 0. so index accordingly!
  // ...actually, avail_guesses should never be zero here anyway
  if (guess_branch[avail_guesses-1] == 0) guess_branch[avail_guesses-1] = num_choices;

  for (int i=1; i < num_choices; ++i) {
    tagint hold = tag_choices[i];
    int j = i - 1;
    while ((j >=0) && (tag_choices[j] > hold)) {
      tag_choices[j+1] = tag_choices[j];
      --j;
    }
    tag_choices[j+1] = hold;
  }

  for (int i = guess_branch[avail_guesses-1]-1; i >= 0; i--) {
    int already_assigned = 0;
    for (int j = 0; j < onemol->natoms; j++) {
      if (glove[j][1] == tag_choices[i]) {
        already_assigned = 1;
        break;
      }
    }
    if (already_assigned == 1) {
      guess_branch[avail_guesses-1]--;
      if (guess_branch[avail_guesses-1] == 0) {
        status = REJECT;
        return;
      }
    } else {
      glove[onemol_xspecial[pion][neigh]-1][0] = onemol_xspecial[pion][neigh];
      glove[onemol_xspecial[pion][neigh]-1][1] = tag_choices[i];
      guess_branch[avail_guesses-1]--;
      break;
    }
  }

  //another check for ghost atoms. perhaps remove the one in make_a_guess
  if (atom->map(glove[(int)onemol_xspecial[pion][neigh]-1][1]) < 0) {
    error->one(FLERR,"Fix bond/react: Fix bond/react needs ghost atoms from further away");
  }

  if (guess_branch[avail_guesses-1] == 0) avail_guesses--;

  for (int i = 0; i < onemol_nxspecial[onemol_xspecial[pion][neigh]-1][0]; i++) {
    pioneer_count[onemol_xspecial[onemol_xspecial[pion][neigh]-1][i]-1]++;
  }
  glove_counter++;
  if (glove_counter == onemol->natoms) {
    if (ring_check() && check_constraints()) status = ACCEPT;
    else status = GUESSFAIL;
    return;
  }
  status = CONTINUE;
}

/* ----------------------------------------------------------------------
  Check that newly assigned atoms have correct bonds
  Necessary for certain ringed structures
------------------------------------------------------------------------- */

int FixBondReact::ring_check()
{
  // ring_check can be made more efficient by re-introducing 'frozen' atoms
  // 'frozen' atoms have been assigned and also are no longer pioneers

  // double check the number of neighbors match for all non-edge atoms
  // otherwise, atoms at 'end' of symmetric ring can behave like edge atoms
  for (int i = 0; i < onemol->natoms; i++)
    if (edge[i][rxnID] == 0 &&
        onemol_nxspecial[i][0] != nxspecial[atom->map(glove[i][1])][0])
      return 0;

  for (int i = 0; i < onemol->natoms; i++) {
    for (int j = 0; j < onemol_nxspecial[i][0]; j++) {
      int ring_fail = 1;
      int ispecial = onemol_xspecial[i][j];
      for (int k = 0; k < nxspecial[atom->map(glove[i][1])][0]; k++) {
        if (xspecial[atom->map(glove[i][1])][k] == glove[ispecial-1][1]) {
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

int FixBondReact::check_constraints()
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

  int *satisfied;
  memory->create(satisfied,nconstraints[rxnID],"bond/react:satisfied");
  for (int i = 0; i < nconstraints[rxnID]; i++)
    satisfied[i] = 1;

  for (int i = 0; i < nconstraints[rxnID]; i++) {
    if (constraints[i][rxnID].type == DISTANCE) {
      get_IDcoords(constraints[i][rxnID].idtype[0], constraints[i][rxnID].id[0], x1);
      get_IDcoords(constraints[i][rxnID].idtype[1], constraints[i][rxnID].id[1], x2);
      delx = x1[0] - x2[0];
      dely = x1[1] - x2[1];
      delz = x1[2] - x2[2];
      domain->minimum_image(delx,dely,delz); // ghost location fix
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < constraints[i][rxnID].par[0] || rsq > constraints[i][rxnID].par[1]) satisfied[i] = 0;
    } else if (constraints[i][rxnID].type == ANGLE) {
      get_IDcoords(constraints[i][rxnID].idtype[0], constraints[i][rxnID].id[0], x1);
      get_IDcoords(constraints[i][rxnID].idtype[1], constraints[i][rxnID].id[1], x2);
      get_IDcoords(constraints[i][rxnID].idtype[2], constraints[i][rxnID].id[2], x3);

      // 1st bond
      delx1 = x1[0] - x2[0];
      dely1 = x1[1] - x2[1];
      delz1 = x1[2] - x2[2];
      domain->minimum_image(delx1,dely1,delz1); // ghost location fix
      rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
      r1 = sqrt(rsq1);

      // 2nd bond
      delx2 = x3[0] - x2[0];
      dely2 = x3[1] - x2[1];
      delz2 = x3[2] - x2[2];
      domain->minimum_image(delx2,dely2,delz2); // ghost location fix
      rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
      r2 = sqrt(rsq2);

      // angle (cos and sin)
      c = delx1*delx2 + dely1*dely2 + delz1*delz2;
      c /= r1*r2;
      if (c > 1.0) c = 1.0;
      if (c < -1.0) c = -1.0;
      if (acos(c) < constraints[i][rxnID].par[0] || acos(c) > constraints[i][rxnID].par[1]) satisfied[i] = 0;
    } else if (constraints[i][rxnID].type == DIHEDRAL) {
      // phi calculation from dihedral style harmonic
      get_IDcoords(constraints[i][rxnID].idtype[0], constraints[i][rxnID].id[0], x1);
      get_IDcoords(constraints[i][rxnID].idtype[1], constraints[i][rxnID].id[1], x2);
      get_IDcoords(constraints[i][rxnID].idtype[2], constraints[i][rxnID].id[2], x3);
      get_IDcoords(constraints[i][rxnID].idtype[3], constraints[i][rxnID].id[3], x4);

      vb1x = x1[0] - x2[0];
      vb1y = x1[1] - x2[1];
      vb1z = x1[2] - x2[2];
      domain->minimum_image(vb1x,vb1y,vb1z);

      vb2x = x3[0] - x2[0];
      vb2y = x3[1] - x2[1];
      vb2z = x3[2] - x2[2];
      domain->minimum_image(vb2x,vb2y,vb2z);

      vb2xm = -vb2x;
      vb2ym = -vb2y;
      vb2zm = -vb2z;
      domain->minimum_image(vb2xm,vb2ym,vb2zm);

      vb3x = x4[0] - x3[0];
      vb3y = x4[1] - x3[1];
      vb3z = x4[2] - x3[2];
      domain->minimum_image(vb3x,vb3y,vb3z);

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
      if (constraints[i][rxnID].par[0] < constraints[i][rxnID].par[1]) {
        if (phi > constraints[i][rxnID].par[0] && phi < constraints[i][rxnID].par[1]) ANDgate = 1;
      } else {
        if (phi > constraints[i][rxnID].par[0] || phi < constraints[i][rxnID].par[1]) ANDgate = 1;
      }
      if (constraints[i][rxnID].par[2] < constraints[i][rxnID].par[3]) {
        if (phi > constraints[i][rxnID].par[2] && phi < constraints[i][rxnID].par[3]) ANDgate = 1;
      } else {
        if (phi > constraints[i][rxnID].par[2] || phi < constraints[i][rxnID].par[3]) ANDgate = 1;
      }
      if (ANDgate != 1) satisfied[i] = 0;
    } else if (constraints[i][rxnID].type == ARRHENIUS) {
      t = get_temperature(glove,0,1);
      prrhob = constraints[i][rxnID].par[1]*pow(t,constraints[i][rxnID].par[2])*
        exp(-constraints[i][rxnID].par[3]/(force->boltz*t));
      if (prrhob < rrhandom[(int) constraints[i][rxnID].par[0]]->uniform()) satisfied[i] = 0;
    } else if (constraints[i][rxnID].type == RMSD) {
      // call superpose
      int iatom;
      int iref = -1; // choose first atom as reference
      int n2superpose = 0;
      double **xfrozen; // coordinates for the "frozen" target molecule
      double **xmobile; // coordinates for the "mobile" molecule
      int ifragment = constraints[i][rxnID].id[0];
      if (ifragment >= 0) {
        for (int j = 0; j < onemol->natoms; j++)
          if (onemol->fragmentmask[ifragment][j]) n2superpose++;
        memory->create(xfrozen,n2superpose,3,"bond/react:xfrozen");
        memory->create(xmobile,n2superpose,3,"bond/react:xmobile");
        int myincr = 0;
        for (int j = 0; j < onemol->natoms; j++) {
          if (onemol->fragmentmask[ifragment][j]) {
            iatom = atom->map(glove[j][1]);
            if (iref == -1) iref = iatom;
            iatom = domain->closest_image(iref,iatom);
            for (int k = 0; k < 3; k++) {
              xfrozen[myincr][k] = x[iatom][k];
              xmobile[myincr][k] = onemol->x[j][k];
            }
            myincr++;
          }
        }
      } else {
        int iatom;
        int iref = -1; // choose first atom as reference
        n2superpose = onemol->natoms;
        memory->create(xfrozen,n2superpose,3,"bond/react:xfrozen");
        memory->create(xmobile,n2superpose,3,"bond/react:xmobile");
        for (int j = 0; j < n2superpose; j++) {
          iatom = atom->map(glove[j][1]);
          if (iref == -1) iref = iatom;
          iatom = domain->closest_image(iref,iatom);
          for (int k = 0; k < 3; k++) {
            xfrozen[j][k] = x[iatom][k];
            xmobile[j][k] = onemol->x[j][k];
          }
        }
      }
      Superpose3D<double, double **> superposer(n2superpose);
      double rmsd = superposer.Superpose(xfrozen, xmobile);
      memory->destroy(xfrozen);
      memory->destroy(xmobile);
      if (rmsd > constraints[i][rxnID].par[0]) satisfied[i] = 0;
    } else if (constraints[i][rxnID].type == CUSTOM) {
      satisfied[i] = custom_constraint(constraints[i][rxnID].str);
    }
  }

  if (nconstraints[rxnID] > 0) {
    char evalstr[MAXLINE],*ptr;
    strcpy(evalstr,constraintstr[rxnID]);
    for (int i = 0; i < nconstraints[rxnID]; i++) {
      ptr = strchr(evalstr,'C');
      *ptr = satisfied[i] ? '1' : '0';
    }
    double verdict = input->variable->evaluate_boolean(evalstr);
    if (verdict == 0.0) {
      memory->destroy(satisfied);
      return 0;
    }
  }

  // let's also check chirality within 'check_constraint'
  for (int i = 0; i < onemol->natoms; i++) {
    if (chiral_atoms[i][0][rxnID] == 1) {
      double my4coords[12];
      // already ensured, by transitive property, that chiral simulation atom has four neighs
      for (int j = 0; j < 4; j++) {
        atom1 = atom->map(glove[i][1]);
        // loop over known types involved in chiral center
        for (int jj = 0; jj < 4; jj++) {
          if (atom->type[atom->map(xspecial[atom1][j])] == chiral_atoms[i][jj+2][rxnID]) {
            atom2 = atom->map(xspecial[atom1][j]);
            atom2 = domain->closest_image(atom1,atom2);
            for (int k = 0; k < 3; k++) {
              my4coords[3*jj+k] = x[atom2][k];
            }
            break;
          }
        }
      }
      if (get_chirality(my4coords) != chiral_atoms[i][1][rxnID]) {
        memory->destroy(satisfied);
        return 0;
      }
    }
  }

  memory->destroy(satisfied);
  return 1;
}

/* ----------------------------------------------------------------------
return pre-reaction atom or fragment location
fragment: given pre-reacted molID (onemol) and fragID,
          return geometric center (of mapped simulation atoms)
------------------------------------------------------------------------- */

void FixBondReact::get_IDcoords(int mode, int myID, double *center)
{
  double **x = atom->x;
  if (mode == ATOM) {
    int iatom = atom->map(glove[myID-1][1]);
    for (int i = 0; i < 3; i++)
      center[i] = x[iatom][i];
  } else {
    int iref = -1; // choose first atom as reference
    int iatom;
    int nfragatoms = 0;
    for (int i = 0; i < 3; i++)
      center[i] = 0;

    for (int i = 0; i < onemol->natoms; i++) {
      if (onemol->fragmentmask[myID][i]) {
        if (iref == -1)
          iref = atom->map(glove[i][1]);
        iatom = atom->map(glove[i][1]);
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

double FixBondReact::get_temperature(tagint **myglove, int row_offset, int col)
{
  int i,ilocal;
  double adof = domain->dimension;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;

  double t = 0.0;

  if (rmass) {
    for (i = 0; i < onemol->natoms; i++) {
      ilocal = atom->map(myglove[i+row_offset][col]);
      t += (v[ilocal][0]*v[ilocal][0] + v[ilocal][1]*v[ilocal][1] +
            v[ilocal][2]*v[ilocal][2]) * rmass[ilocal];
    }
  } else {
    for (i = 0; i < onemol->natoms; i++) {
      ilocal = atom->map(myglove[i+row_offset][col]);
      t += (v[ilocal][0]*v[ilocal][0] + v[ilocal][1]*v[ilocal][1] +
            v[ilocal][2]*v[ilocal][2]) * mass[type[ilocal]];
    }
  }

  // final temperature
  double dof = adof*onemol->natoms;
  double tfactor = force->mvv2e / (dof * force->boltz);
  t *= tfactor;
  return t;
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

  for (rxnID = 0; rxnID < nreacts; rxnID++) {
    for (int i = 0; i < nconstraints[rxnID]; i++) {
      if (constraints[i][rxnID].type == CUSTOM) {
        varstr = constraints[i][rxnID].str;
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

double FixBondReact::custom_constraint(const std::string& varstr)
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
    evlstr.push_back(std::to_string(rxnfunction(rxnfunclist[irxnfunc], varid, fragid)));
  }
  evlstr.push_back(varstr.substr(prev3+1));

  for (auto & evl : evlstr) evlcat += evl;
  return input->variable->compute_equal(evlcat);
}

/* ----------------------------------------------------------------------
currently three 'rxn' functions: rxnsum, rxnave, and rxnbond
------------------------------------------------------------------------- */

double FixBondReact::rxnfunction(const std::string& rxnfunc, const std::string& varid,
                                 const std::string& fragid)
{
  int ifrag = -1;
  if (fragid != "all") {
    ifrag = onemol->findfragment(fragid.c_str());
    if (ifrag < 0) error->one(FLERR,"Bond/react: Molecule fragment "
                              "in reaction special function does not exist");
  }

  // start with 'rxnbond' per-bond function
  // for 'rxnbond', varid corresponds to 'compute bond/local' name,
  //                and fragid is a pre-reaction fragment containing the two atoms in the bond
  if (rxnfunc == "rxnbond") {
    int icompute,ibond,nsum;
    double perbondval;
    std::set<tagint> aset;
    std::string computeid = varid;
    std::map<std::set<tagint>,int>::iterator it;

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

    nsum = 0;
    for (int i = 0; i < onemol->natoms; i++) {
      if (onemol->fragmentmask[ifrag][i]) {
        aset.insert(glove[i][1]);
        nsum++;
      }
    }
    if (nsum != 2) error->one(FLERR,"Bond/react: Molecule fragment of reaction special function 'rxnbond' "
                     "must contain exactly two atoms");

    if (cperbond->invoked_local != lmp->update->ntimestep)
      cperbond->compute_local();

    it = atoms2bond.find(aset);
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
      for (int i = 0; i < onemol->natoms; i++) {
        iatom = atom->map(glove[i][1]);
        sumvvec += vvec[iatom][ivar];
      }
      nsum = onemol->natoms;
    } else {
      for (int i = 0; i < onemol->natoms; i++) {
        if (onemol->fragmentmask[ifrag][i]) {
          iatom = atom->map(glove[i][1]);
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
  Get xspecials for current molecule templates
------------------------------------------------------------------------- */

void FixBondReact::get_molxspecials()
{
  if (newton_bond == 1) {
    onemol_nxspecial = onemol->nspecial;
    onemol_xspecial = onemol->special;
    twomol_nxspecial = twomol->nspecial;
    twomol_xspecial = twomol->special;
  } else {
    memory->destroy(onemol_nxspecial);
    memory->destroy(onemol_xspecial);
    memory->create(onemol_nxspecial,onemol->natoms,3,"bond/react:onemol_nxspecial");
    memory->create(onemol_xspecial,onemol->natoms,atom->maxspecial,"bond/react:onemol_xspecial");
    for (int i = 0; i < onemol->natoms; i++) {
      onemol_nxspecial[i][0] = onemol->num_bond[i];
      for (int j = 0; j < onemol_nxspecial[i][0]; j++) {
        onemol_xspecial[i][j] = onemol->bond_atom[i][j];
      }
      onemol_nxspecial[i][1] = onemol->nspecial[i][1];
      onemol_nxspecial[i][2] = onemol->nspecial[i][2];
      int joffset = onemol_nxspecial[i][0] - onemol->nspecial[i][0];
      for (int j = onemol_nxspecial[i][0]; j < onemol_nxspecial[i][2]; j++) {
        onemol_xspecial[i][j+joffset] = onemol->special[i][j];
      }
    }
    memory->destroy(twomol_nxspecial);
    memory->destroy(twomol_xspecial);
    memory->create(twomol_nxspecial,twomol->natoms,3,"bond/react:twomol_nxspecial");
    memory->create(twomol_xspecial,twomol->natoms,atom->maxspecial,"bond/react:twomol_xspecial");
    for (int i = 0; i < twomol->natoms; i++) {
      twomol_nxspecial[i][0] = twomol->num_bond[i];
      for (int j = 0; j < twomol_nxspecial[i][0]; j++) {
        twomol_xspecial[i][j] = twomol->bond_atom[i][j];
      }
      twomol_nxspecial[i][1] = twomol->nspecial[i][1];
      twomol_nxspecial[i][2] = twomol->nspecial[i][2];
      int joffset = twomol_nxspecial[i][0] - twomol->nspecial[i][0];
      for (int j = twomol_nxspecial[i][0]; j < twomol_nxspecial[i][2]; j++) {
        twomol_xspecial[i][j+joffset] = twomol->special[i][j];
      }
    }
  }
}

/* ----------------------------------------------------------------------
  Determine which pre-reacted template atoms are at least three bonds
  away from edge atoms.
------------------------------------------------------------------------- */

void FixBondReact::find_landlocked_atoms(int myrxn)
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
  for (int i = 0; i < twomol->natoms; i++) {
    if (create_atoms[i][myrxn] == 0 && edge[equivalences[i][1][myrxn]-1][myrxn] == 1)
      landlocked_atoms[i][myrxn] = 0;
    else landlocked_atoms[i][myrxn] = 1;
  }
  int nspecial_limit = -1;
  if (force->angle && twomol->angleflag) nspecial_limit = 0;

  if ((force->dihedral && twomol->dihedralflag) ||
      (force->improper && twomol->improperflag)) nspecial_limit = 1;

  if (nspecial_limit != -1) {
    for (int i = 0; i < twomol->natoms; i++) {
      for (int j = 0; j < twomol_nxspecial[i][nspecial_limit]; j++) {
        for (int k = 0; k < onemol->natoms; k++) {
          if (equivalences[twomol_xspecial[i][j]-1][1][myrxn] == k+1 && edge[k][myrxn] == 1) {
            landlocked_atoms[i][myrxn] = 0;
          }
        }
      }
    }
  }

  // bad molecule templates check
  // if atoms change types, but aren't landlocked, that's bad
  for (int i = 0; i < twomol->natoms; i++) {
    if ((create_atoms[i][myrxn] == 0) &&
        (twomol->type[i] != onemol->type[equivalences[i][1][myrxn]-1]) &&
        (landlocked_atoms[i][myrxn] == 0))
      error->all(FLERR, "Fix bond/react: Atom type affected by reaction {} is too close "
                 "to template edge", rxn_name[myrxn]);
  }

  // additionally, if a bond changes type, but neither involved atom is landlocked, bad
  // would someone want to change an angle type but not bond or atom types? (etc.) ...hopefully not yet
  for (int i = 0; i < twomol->natoms; i++) {
    if (create_atoms[i][myrxn] == 0) {
      if (landlocked_atoms[i][myrxn] == 0) {
        for (int j = 0; j < twomol->num_bond[i]; j++) {
          int twomol_atomj = twomol->bond_atom[i][j];
          if (landlocked_atoms[twomol_atomj-1][myrxn] == 0) {
            int onemol_atomi = equivalences[i][1][myrxn];
            int onemol_batom;
            for (int m = 0; m < onemol->num_bond[onemol_atomi-1]; m++) {
              onemol_batom = onemol->bond_atom[onemol_atomi-1][m];
              if ((onemol_batom == equivalences[twomol_atomj-1][1][myrxn]) &&
                  (twomol->bond_type[i][j] != onemol->bond_type[onemol_atomi-1][m]))
                error->all(FLERR, "Fix bond/react: Bond type affected by reaction {} is "
                           "too close to template edge",rxn_name[myrxn]);
            }
            if (newton_bond) {
              int onemol_atomj = equivalences[twomol_atomj-1][1][myrxn];
              for (int m = 0; m < onemol->num_bond[onemol_atomj-1]; m++) {
                onemol_batom = onemol->bond_atom[onemol_atomj-1][m];
                if ((onemol_batom == equivalences[i][1][myrxn]) &&
                    (twomol->bond_type[i][j] != onemol->bond_type[onemol_atomj-1][m]))
                  error->all(FLERR, "Fix bond/react: Bond type affected by reaction {} is "
                             "too close to template edge",rxn_name[myrxn]);
              }
            }
          }
        }
      }
    }
  }

  // additionally, if a deleted atom is bonded to an atom that is not deleted, bad
  for (int i = 0; i < onemol->natoms; i++) {
    if (delete_atoms[i][myrxn] == 1) {
      int ii = reverse_equiv[i][1][myrxn] - 1;
      for (int j = 0; j < twomol_nxspecial[ii][0]; j++) {
        if (delete_atoms[equivalences[twomol_xspecial[ii][j]-1][1][myrxn]-1][myrxn] == 0) {
          error->all(FLERR,"Fix bond/react: A deleted atom cannot remain bonded to an atom that is not deleted");
        }
      }
    }
  }

  // also, if atoms change number of bonds, but aren't landlocked, that could be bad
  if (comm->me == 0)
    for (int i = 0; i < twomol->natoms; i++) {
      if ((create_atoms[i][myrxn] == 0) &&
          (twomol_nxspecial[i][0] != onemol_nxspecial[equivalences[i][1][myrxn]-1][0]) &&
          (landlocked_atoms[i][myrxn] == 0))
        error->warning(FLERR, "Fix bond/react: Atom affected by reaction {} is too close "
                       "to template edge",rxn_name[myrxn]);
          break;
    }

  // finally, if a created atom is not landlocked, bad!
  for (int i = 0; i < twomol->natoms; i++) {
    if (create_atoms[i][myrxn] == 1 && landlocked_atoms[i][myrxn] == 0) {
      error->one(FLERR,"Fix bond/react: Created atom too close to template edge");
    }
  }
}

/* ----------------------------------------------------------------------
let's dedup global_mega_glove
allows for same site undergoing different pathways, in parallel
------------------------------------------------------------------------- */

void FixBondReact::dedup_mega_gloves(int dedup_mode)
{
  // dedup_mode == LOCAL for local_dedup
  // dedup_mode == GLOBAL for global_mega_glove

  if (dedup_mode == GLOBAL)
    for (int i = 0; i < nreacts; i++)
      ghostly_rxn_count[i] = 0;

  int dedup_size = 0;
  if (dedup_mode == LOCAL) {
    dedup_size = my_num_mega;
  } else if (dedup_mode == GLOBAL) {
    dedup_size = global_megasize;
  }

  tagint **dedup_glove;
  memory->create(dedup_glove,max_natoms+1,dedup_size,"bond/react:dedup_glove");

  if (dedup_mode == LOCAL) {
    for (int i = 0; i < dedup_size; i++) {
      for (int j = 0; j < max_natoms+1; j++) {
        dedup_glove[j][i] = my_mega_glove[j][i];
      }
    }
  } else if (dedup_mode == GLOBAL) {
    for (int i = 0; i < dedup_size; i++) {
      for (int j = 0; j < max_natoms+1; j++) {
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
  int *temp_rxn = new int[max_natoms+1];
  for (int i = dedup_size-1; i > 0; --i) { //dedup_size
    // choose random entry to swap current one with
    int k = floor(random[0]->uniform()*(i+1));

    // swap entries
    for (int j = 0; j < max_natoms+1; j++)
      temp_rxn[j] = dedup_glove[j][i];

    for (int j = 0; j < max_natoms+1; j++) {
      dedup_glove[j][i] = dedup_glove[j][k];
      dedup_glove[j][k] = temp_rxn[j];
    }
  }
  delete [] temp_rxn;

  for (int i = 0; i < dedup_size; i++) {
    if (dedup_mask[i] == 0) {
      int myrxnid1 = dedup_glove[0][i];
      onemol = atom->molecules[unreacted_mol[myrxnid1]];
      for (int j = 0; j < onemol->natoms; j++) {
        int check1 = dedup_glove[j+1][i];
        for (int ii = i + 1; ii < dedup_size; ii++) {
          if (dedup_mask[ii] == 0) {
            int myrxnid2 = dedup_glove[0][ii];
            twomol = atom->molecules[unreacted_mol[myrxnid2]];
            for (int jj = 0; jj < twomol->natoms; jj++) {
              int check2 = dedup_glove[jj+1][ii];
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
  if (dedup_mode == LOCAL) {
    int my_new_megasize = 0;
    for (int i = 0; i < my_num_mega; i++) {
      if (dedup_mask[i] == 0) {
        for (int j = 0; j < max_natoms+1; j++) {
          my_mega_glove[j][my_new_megasize] = dedup_glove[j][i];
        }
        my_new_megasize++;
      }
    }
    my_num_mega = my_new_megasize;
  }

  // we must update global_mega_glove and global_megasize
  // we can simply overwrite global_mega_glove column by column
  if (dedup_mode == GLOBAL) {
    int new_global_megasize = 0;
    for (int i = 0; i < global_megasize; i++) {
      if (dedup_mask[i] == 0) {
        ghostly_rxn_count[dedup_glove[0][i]]++;
        for (int j = 0; j < max_natoms + 1; j++) {
          global_mega_glove[j][new_global_megasize] = dedup_glove[j][i];
        }
        new_global_megasize++;
      }
    }
    global_megasize = new_global_megasize;
  }

  memory->destroy(dedup_glove);
  delete [] dedup_mask;
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
    int index2 = atom->find_custom(statted_id,flag,cols);
    i_statted_tags = atom->ivector[index2];
  }

  int index3 = atom->find_custom("react_tags",flag,cols);
  int *i_react_tags = atom->ivector[index3];

  int unlimitflag = 0;
  for (int i = 0; i < atom->nlocal; i++) {
    // unlimit atoms for next step! this resolves # of procs disparity, mostly
    // first '1': indexing offset, second '1': for next step
    if (i_limit_tags[i] != 0 && (update->ntimestep + 1 - i_limit_tags[i]) > limit_duration[i_react_tags[i]]) {
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

  for (int i = 0; i < nreacts; i++)
    local_rxn_count[i] = 0;

  for (int i = 0; i < my_num_mega; i++) {
    rxnID = my_mega_glove[0][i];
    onemol = atom->molecules[unreacted_mol[rxnID]];
    int ghostly = 0;
  #if !defined(MPI_STUBS)
    if (comm->style == Comm::BRICK) {
      if (create_atoms_flag[rxnID] == 1) {
        ghostly = 1;
      } else {
        for (int j = 0; j < onemol->natoms; j++) {
          int ilocal = atom->map(my_mega_glove[j+1][i]);
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
      ghostly_mega_glove[0][ghostly_num_mega] = rxnID;
      for (int j = 0; j < onemol->natoms+1; j++) {
        ghostly_mega_glove[j][ghostly_num_mega] = my_mega_glove[j][i];
      }
      ghostly_num_mega++;
    } else {
      local_mega_glove[0][local_num_mega] = rxnID;
      local_rxn_count[rxnID]++;
      for (int j = 0; j < onemol->natoms+1; j++) {
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
    delete [] allncols;
    return;
  }

  int *allstarts = new int[nprocs];

  int start = 0;
  for (int i = 0; i < comm->me; i++) {
    start += allncols[i];
  }
  MPI_Allgather(&start, 1, MPI_INT, allstarts, 1, MPI_INT, world);
  MPI_Datatype columnunsized, column;
  int sizes[2]    = {max_natoms+1, global_megasize};
  int subsizes[2] = {max_natoms+1, 1};
  int starts[2]   = {0,0};
  MPI_Type_create_subarray (2, sizes, subsizes, starts, MPI_ORDER_C,
                            MPI_LMP_TAGINT, &columnunsized);
  MPI_Type_create_resized (columnunsized, 0, sizeof(tagint), &column);
  MPI_Type_commit(&column);

  memory->destroy(global_mega_glove);
  memory->create(global_mega_glove,max_natoms+1,global_megasize,"bond/react:global_mega_glove");

  for (int i = 0; i < max_natoms+1; i++)
    for (int j = 0; j < global_megasize; j++)
      global_mega_glove[i][j] = 0;

  if (ghostly_num_mega > 0) {
    for (int i = 0; i < max_natoms+1; i++) {
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

  if (comm->me == 0) dedup_mega_gloves(GLOBAL); // global_mega_glove mode
  MPI_Bcast(&global_megasize,1,MPI_INT,0,world);
  MPI_Bcast(&(global_mega_glove[0][0]), global_megasize, column, 0, world);

  delete [] allstarts;
  delete [] allncols;

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

  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;

  // used when deleting atoms
  int ndel,ndelone;
  int *mark;
  int nmark = nlocal;
  memory->create(mark,nmark,"bond/react:mark");
  for (int i = 0; i < nmark; i++) mark[i] = 0;

  // flag used to delete special interactions
  int *delflag;
  memory->create(delflag,atom->maxspecial,"bond/react:delflag");

  tagint *tag = atom->tag;
  AtomVec *avec = atom->avec;

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
    int index2 = atom->find_custom(statted_id,flag,cols);
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
  memory->create(update_mega_glove,max_natoms+1,MAX(local_num_mega,global_megasize),"bond/react:update_mega_glove");

  for (int pass = 0; pass < 2; pass++) {
    update_num_mega = 0;
    int *iskip = new int[nreacts];
    for (int i = 0; i < nreacts; i++) iskip[i] = 0;
    if (pass == 0) {
      for (int i = 0; i < local_num_mega; i++) {
        rxnID = local_mega_glove[0][i];
        // reactions already shuffled from dedup procedure, so can skip first N
        if (iskip[rxnID]++ < nlocalskips[rxnID]) continue;

        // atoms inserted here for serial MPI_STUBS build only
        if (create_atoms_flag[rxnID] == 1) {
          onemol = atom->molecules[unreacted_mol[rxnID]];
          twomol = atom->molecules[reacted_mol[rxnID]];
          if (insert_atoms(local_mega_glove,i)) {
            inserted_atoms_flag = 1;
          } else { // create aborted
            reaction_count_total[rxnID]--;
            continue;
          }
        }

        for (int j = 0; j < max_natoms+1; j++)
          update_mega_glove[j][update_num_mega] = local_mega_glove[j][i];
        update_num_mega++;
      }
    } else if (pass == 1) {
      for (int i = 0; i < global_megasize; i++) {
        rxnID = global_mega_glove[0][i];
        // reactions already shuffled from dedup procedure, so can skip first N
        if (iskip[rxnID]++ < nghostlyskips[rxnID]) continue;

        // we can insert atoms here, now that reactions are finalized
        // can't do it any earlier, due to skipped reactions (max_rxn)
        // for MPI build, reactions that create atoms are always treated as 'global'
        if (create_atoms_flag[rxnID] == 1) {
          onemol = atom->molecules[unreacted_mol[rxnID]];
          twomol = atom->molecules[reacted_mol[rxnID]];
          if (insert_atoms(global_mega_glove,i)) {
            inserted_atoms_flag = 1;
          } else { // create aborted
            reaction_count_total[rxnID]--;
            continue;
          }
        }

        for (int j = 0; j < max_natoms+1; j++)
          update_mega_glove[j][update_num_mega] = global_mega_glove[j][i];
        update_num_mega++;
      }
    }
    delete [] iskip;

    if (update_num_mega == 0) continue;

    // if inserted atoms and global map exists, reset map now instead
    //   of waiting for comm since other pre-exchange fixes may use it
    // invoke map_init() b/c atom count has grown
    // do this once after all atom insertions
    if (inserted_atoms_flag == 1 && atom->map_style != Atom::MAP_NONE) {
      atom->map_init();
      atom->map_set();
    }

    // mark to-delete atoms
    nlocal = atom->nlocal;
    if (nlocal > nmark) {
      memory->grow(mark,nlocal,"bond/react:mark");
      for (int i = nmark; i < nlocal; i++) mark[i] = 0;
      nmark = nlocal;
    }
    for (int i = 0; i < update_num_mega; i++) {
      rxnID = update_mega_glove[0][i];
      onemol = atom->molecules[unreacted_mol[rxnID]];
      for (int j = 0; j < onemol->natoms; j++) {
        int iatom = atom->map(update_mega_glove[j+1][i]);
        if (delete_atoms[j][rxnID] == 1 && iatom >= 0 && iatom < nlocal) {
          mark[iatom] = 1;
        }
      }
    }

    // get charge rescale delta
    double charge_rescale_addend = 0;
    if (rescale_charges_flag[rxnID] == 1) {
      double sim_total_charge = 0;
      double mol_total_charge = 0;
      int n_custom_charge = 0;
      for (int i = 0; i < update_num_mega; i++) {
        rxnID = update_mega_glove[0][i];
        twomol = atom->molecules[reacted_mol[rxnID]];
        for (int j = 0; j < twomol->natoms; j++) {
          int jj = equivalences[j][1][rxnID]-1;
          if (atom->map(update_mega_glove[jj+1][i]) >= 0 &&
              atom->map(update_mega_glove[jj+1][i]) < nlocal) {
            if (landlocked_atoms[j][rxnID] == 1)
              type[atom->map(update_mega_glove[jj+1][i])] = twomol->type[j];
            if (twomol->qflag && atom->q_flag && custom_charges[jj][rxnID] == 1) {
              double *q = atom->q;
              sim_total_charge += q[atom->map(update_mega_glove[jj+1][i])];
              mol_total_charge += twomol->q[j];
              n_custom_charge++;
            }
          }
        }
      }
      charge_rescale_addend = (sim_total_charge-mol_total_charge)/n_custom_charge;
    }

    // update charges and types of landlocked atoms
    // also keep track of 'stabilization' groups here
    for (int i = 0; i < update_num_mega; i++) {
      rxnID = update_mega_glove[0][i];
      twomol = atom->molecules[reacted_mol[rxnID]];
      for (int j = 0; j < twomol->natoms; j++) {
        int jj = equivalences[j][1][rxnID]-1;
        int ilocal = atom->map(update_mega_glove[jj+1][i]);
        if (ilocal >= 0 && ilocal < nlocal) {

          // update->ntimestep could be 0. so add 1 throughout
          i_limit_tags[ilocal] = update->ntimestep + 1;
          if (stabilization_flag == 1) i_statted_tags[ilocal] = 0;
          i_react_tags[ilocal] = rxnID;

          if (landlocked_atoms[j][rxnID] == 1)
            type[ilocal] = twomol->type[j];
          if (twomol->qflag && atom->q_flag && custom_charges[jj][rxnID] == 1) {
            double *q = atom->q;
            q[ilocal] = twomol->q[j]+charge_rescale_addend;
          }
        }
      }
    }

    int insert_num;
    // very nice and easy to completely overwrite special bond info for landlocked atoms
    for (int i = 0; i < update_num_mega; i++) {
      rxnID = update_mega_glove[0][i];
      onemol = atom->molecules[unreacted_mol[rxnID]];
      twomol = atom->molecules[reacted_mol[rxnID]];
      for (int j = 0; j < twomol->natoms; j++) {
        int jj = equivalences[j][1][rxnID]-1;
        int ilocal = atom->map(update_mega_glove[jj+1][i]);
        if (ilocal < nlocal && ilocal >= 0) {
          if (landlocked_atoms[j][rxnID] == 1) {
            for (int k = 0; k < 3; k++) {
              nspecial[ilocal][k] = twomol->nspecial[j][k];
            }
            for (int p = 0; p < twomol->nspecial[j][2]; p++) {
              special[ilocal][p] = update_mega_glove[equivalences[twomol->special[j][p]-1][1][rxnID]][i];
            }
          }
          // now delete and replace landlocked atoms from non-landlocked atoms' special info
          // delete 1-2, 1-3, 1-4 specials individually. only delete if special exists in pre-reaction template
          if (landlocked_atoms[j][rxnID] == 0) {
            int ispec, fspec, imolspec, fmolspec, nspecdel[3];
            for (int k = 0; k < 3; k++) nspecdel[k] = 0;
            for (int k = 0; k < atom->maxspecial; k++) delflag[k] = 0;
            for (int specn = 0; specn < 3; specn++) {
              if (specn == 0) {
                imolspec = 0;
                ispec = 0;
              } else {
                imolspec = onemol->nspecial[jj][specn-1];
                ispec = nspecial[ilocal][specn-1];
              }
              fmolspec = onemol->nspecial[jj][specn];
              fspec = nspecial[ilocal][specn];
              for (int k = ispec; k < fspec; k++) {
                for (int p = imolspec; p < fmolspec; p++) {
                  if (update_mega_glove[onemol->special[jj][p]][i] == special[ilocal][k]) {
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
            for (int k = 0; k < twomol->nspecial[j][2]; k++) {
              if (k > twomol->nspecial[j][1] - 1) {
                insert_num = nspecial[ilocal][2]++;
              } else if (k > twomol->nspecial[j][0] - 1) {
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
              special[ilocal][insert_num] = update_mega_glove[equivalences[twomol->special[j][k]-1][1][rxnID]][i];
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
      rxnID = update_mega_glove[0][i];
      twomol = atom->molecules[reacted_mol[rxnID]];
      // let's first delete all bond info about landlocked atoms
      for (int j = 0; j < twomol->natoms; j++) {
        int jj = equivalences[j][1][rxnID]-1;
        if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
          if (landlocked_atoms[j][rxnID] == 1) {
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
          if (landlocked_atoms[j][rxnID] == 0) {
            for (int p = num_bond[atom->map(update_mega_glove[jj+1][i])]-1; p > -1 ; p--) {
              for (int n = 0; n < twomol->natoms; n++) {
                int nn = equivalences[n][1][rxnID]-1;
                if (n!=j && bond_atom[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i] && landlocked_atoms[n][rxnID] == 1) {
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
      for (int j = 0; j < twomol->natoms; j++) {
        int jj = equivalences[j][1][rxnID]-1;
        if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
          if (landlocked_atoms[j][rxnID] == 1)  {
            num_bond[atom->map(update_mega_glove[jj+1][i])] = twomol->num_bond[j];
            delta_bonds += twomol->num_bond[j];
            for (int p = 0; p < twomol->num_bond[j]; p++) {
              bond_type[atom->map(update_mega_glove[jj+1][i])][p] = twomol->bond_type[j][p];
              bond_atom[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[equivalences[twomol->bond_atom[j][p]-1][1][rxnID]][i];
              // Check cached history data to see if bond regenerated
              if (n_histories > 0)
                for (auto &ihistory: histories)
                  dynamic_cast<FixBondHistory *>(ihistory)->check_cache(atom->map(update_mega_glove[jj+1][i]), p);
            }
          }
          if (landlocked_atoms[j][rxnID] == 0) {
            for (int p = 0; p < twomol->num_bond[j]; p++) {
              if (landlocked_atoms[twomol->bond_atom[j][p]-1][rxnID] == 1) {
                insert_num = num_bond[atom->map(update_mega_glove[jj+1][i])];
                bond_type[atom->map(update_mega_glove[jj+1][i])][insert_num] = twomol->bond_type[j][p];
                bond_atom[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->bond_atom[j][p]-1][1][rxnID]][i];
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
    if (force->angle && twomol->angleflag) {
      int *num_angle = atom->num_angle;
      int **angle_type = atom->angle_type;
      tagint **angle_atom1 = atom->angle_atom1;
      tagint **angle_atom2 = atom->angle_atom2;
      tagint **angle_atom3 = atom->angle_atom3;

      for (int i = 0; i < update_num_mega; i++) {
        rxnID = update_mega_glove[0][i];
        twomol = atom->molecules[reacted_mol[rxnID]];
        for (int j = 0; j < twomol->natoms; j++) {
          int jj = equivalences[j][1][rxnID]-1;
          if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
            if (landlocked_atoms[j][rxnID] == 1) {
              delta_angle -= num_angle[atom->map(update_mega_glove[jj+1][i])];
              num_angle[atom->map(update_mega_glove[jj+1][i])] = 0;
            }
            if (landlocked_atoms[j][rxnID] == 0) {
              for (int p = num_angle[atom->map(update_mega_glove[jj+1][i])]-1; p > -1; p--) {
                for (int n = 0; n < twomol->natoms; n++) {
                  int nn = equivalences[n][1][rxnID]-1;
                  if (n!=j && landlocked_atoms[n][rxnID] == 1 &&
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
        for (int j = 0; j < twomol->natoms; j++) {
          int jj = equivalences[j][1][rxnID]-1;
          if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
            if (landlocked_atoms[j][rxnID] == 1) {
              num_angle[atom->map(update_mega_glove[jj+1][i])] = twomol->num_angle[j];
              delta_angle += twomol->num_angle[j];
              for (int p = 0; p < twomol->num_angle[j]; p++) {
                angle_type[atom->map(update_mega_glove[jj+1][i])][p] = twomol->angle_type[j][p];
                angle_atom1[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[equivalences[twomol->angle_atom1[j][p]-1][1][rxnID]][i];
                angle_atom2[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[equivalences[twomol->angle_atom2[j][p]-1][1][rxnID]][i];
                angle_atom3[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[equivalences[twomol->angle_atom3[j][p]-1][1][rxnID]][i];
              }
            }
            if (landlocked_atoms[j][rxnID] == 0) {
              for (int p = 0; p < twomol->num_angle[j]; p++) {
                if (landlocked_atoms[twomol->angle_atom1[j][p]-1][rxnID] == 1 ||
                    landlocked_atoms[twomol->angle_atom2[j][p]-1][rxnID] == 1 ||
                    landlocked_atoms[twomol->angle_atom3[j][p]-1][rxnID] == 1) {
                  insert_num = num_angle[atom->map(update_mega_glove[jj+1][i])];
                  angle_type[atom->map(update_mega_glove[jj+1][i])][insert_num] = twomol->angle_type[j][p];
                  angle_atom1[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->angle_atom1[j][p]-1][1][rxnID]][i];
                  angle_atom2[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->angle_atom2[j][p]-1][1][rxnID]][i];
                  angle_atom3[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->angle_atom3[j][p]-1][1][rxnID]][i];
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

    // Dihedrals! first let's delete all dihedral info for landlocked atoms
    if (force->dihedral && twomol->dihedralflag) {
      int *num_dihedral = atom->num_dihedral;
      int **dihedral_type = atom->dihedral_type;
      tagint **dihedral_atom1 = atom->dihedral_atom1;
      tagint **dihedral_atom2 = atom->dihedral_atom2;
      tagint **dihedral_atom3 = atom->dihedral_atom3;
      tagint **dihedral_atom4 = atom->dihedral_atom4;

      for (int i = 0; i < update_num_mega; i++) {
        rxnID = update_mega_glove[0][i];
        twomol = atom->molecules[reacted_mol[rxnID]];
        for (int j = 0; j < twomol->natoms; j++) {
          int jj = equivalences[j][1][rxnID]-1;
          if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
            if (landlocked_atoms[j][rxnID] == 1) {
              delta_dihed -= num_dihedral[atom->map(update_mega_glove[jj+1][i])];
              num_dihedral[atom->map(update_mega_glove[jj+1][i])] = 0;
            }
            if (landlocked_atoms[j][rxnID] == 0) {
              for (int p = num_dihedral[atom->map(update_mega_glove[jj+1][i])]-1; p > -1; p--) {
                for (int n = 0; n < twomol->natoms; n++) {
                  int nn = equivalences[n][1][rxnID]-1;
                  if (n!=j && landlocked_atoms[n][rxnID] == 1 &&
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
        for (int j = 0; j < twomol->natoms; j++) {
          int jj = equivalences[j][1][rxnID]-1;
          if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
            if (landlocked_atoms[j][rxnID] == 1) {
              num_dihedral[atom->map(update_mega_glove[jj+1][i])] = twomol->num_dihedral[j];
              delta_dihed += twomol->num_dihedral[j];
              for (int p = 0; p < twomol->num_dihedral[j]; p++) {
                dihedral_type[atom->map(update_mega_glove[jj+1][i])][p] = twomol->dihedral_type[j][p];
                dihedral_atom1[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[equivalences[twomol->dihedral_atom1[j][p]-1][1][rxnID]][i];
                dihedral_atom2[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[equivalences[twomol->dihedral_atom2[j][p]-1][1][rxnID]][i];
                dihedral_atom3[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[equivalences[twomol->dihedral_atom3[j][p]-1][1][rxnID]][i];
                dihedral_atom4[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[equivalences[twomol->dihedral_atom4[j][p]-1][1][rxnID]][i];
              }
            }
            if (landlocked_atoms[j][rxnID] == 0) {
              for (int p = 0; p < twomol->num_dihedral[j]; p++) {
                if (landlocked_atoms[twomol->dihedral_atom1[j][p]-1][rxnID] == 1 ||
                    landlocked_atoms[twomol->dihedral_atom2[j][p]-1][rxnID] == 1 ||
                    landlocked_atoms[twomol->dihedral_atom3[j][p]-1][rxnID] == 1 ||
                    landlocked_atoms[twomol->dihedral_atom4[j][p]-1][rxnID] == 1) {
                  insert_num = num_dihedral[atom->map(update_mega_glove[jj+1][i])];
                  dihedral_type[atom->map(update_mega_glove[jj+1][i])][insert_num] = twomol->dihedral_type[j][p];
                  dihedral_atom1[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->dihedral_atom1[j][p]-1][1][rxnID]][i];
                  dihedral_atom2[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->dihedral_atom2[j][p]-1][1][rxnID]][i];
                  dihedral_atom3[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->dihedral_atom3[j][p]-1][1][rxnID]][i];
                  dihedral_atom4[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->dihedral_atom4[j][p]-1][1][rxnID]][i];
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

    // finally IMPROPERS!!!! first let's delete all improper info for landlocked atoms
    if (force->improper && twomol->improperflag) {
      int *num_improper = atom->num_improper;
      int **improper_type = atom->improper_type;
      tagint **improper_atom1 = atom->improper_atom1;
      tagint **improper_atom2 = atom->improper_atom2;
      tagint **improper_atom3 = atom->improper_atom3;
      tagint **improper_atom4 = atom->improper_atom4;

      for (int i = 0; i < update_num_mega; i++) {
        rxnID = update_mega_glove[0][i];
        twomol = atom->molecules[reacted_mol[rxnID]];
        for (int j = 0; j < twomol->natoms; j++) {
          int jj = equivalences[j][1][rxnID]-1;
          if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
            if (landlocked_atoms[j][rxnID] == 1) {
              delta_imprp -= num_improper[atom->map(update_mega_glove[jj+1][i])];
              num_improper[atom->map(update_mega_glove[jj+1][i])] = 0;
            }
            if (landlocked_atoms[j][rxnID] == 0) {
              for (int p = num_improper[atom->map(update_mega_glove[jj+1][i])]-1; p > -1; p--) {
                for (int n = 0; n < twomol->natoms; n++) {
                  int nn = equivalences[n][1][rxnID]-1;
                  if (n!=j && landlocked_atoms[n][rxnID] == 1 &&
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
        for (int j = 0; j < twomol->natoms; j++) {
          int jj = equivalences[j][1][rxnID]-1;
          if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
            if (landlocked_atoms[j][rxnID] == 1) {
              num_improper[atom->map(update_mega_glove[jj+1][i])] = twomol->num_improper[j];
              delta_imprp += twomol->num_improper[j];
              for (int p = 0; p < twomol->num_improper[j]; p++) {
                improper_type[atom->map(update_mega_glove[jj+1][i])][p] = twomol->improper_type[j][p];
                improper_atom1[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[equivalences[twomol->improper_atom1[j][p]-1][1][rxnID]][i];
                improper_atom2[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[equivalences[twomol->improper_atom2[j][p]-1][1][rxnID]][i];
                improper_atom3[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[equivalences[twomol->improper_atom3[j][p]-1][1][rxnID]][i];
                improper_atom4[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[equivalences[twomol->improper_atom4[j][p]-1][1][rxnID]][i];
              }
            }
            if (landlocked_atoms[j][rxnID] == 0) {
              for (int p = 0; p < twomol->num_improper[j]; p++) {
                if (landlocked_atoms[twomol->improper_atom1[j][p]-1][rxnID] == 1 ||
                    landlocked_atoms[twomol->improper_atom2[j][p]-1][rxnID] == 1 ||
                    landlocked_atoms[twomol->improper_atom3[j][p]-1][rxnID] == 1 ||
                    landlocked_atoms[twomol->improper_atom4[j][p]-1][rxnID] == 1) {
                  insert_num = num_improper[atom->map(update_mega_glove[jj+1][i])];
                  improper_type[atom->map(update_mega_glove[jj+1][i])][insert_num] = twomol->improper_type[j][p];
                  improper_atom1[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->improper_atom1[j][p]-1][1][rxnID]][i];
                  improper_atom2[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->improper_atom2[j][p]-1][1][rxnID]][i];
                  improper_atom3[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->improper_atom3[j][p]-1][1][rxnID]][i];
                  improper_atom4[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->improper_atom4[j][p]-1][1][rxnID]][i];
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

  memory->destroy(update_mega_glove);

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

  // reset mol ids
  if (reset_mol_ids_flag) reset_mol_ids->reset();

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
insert created atoms
------------------------------------------------------------------------- */

int FixBondReact::insert_atoms(tagint **my_mega_glove, int iupdate)
{
  // inserting atoms based off fix_deposit->pre_exchange
  int flag;
  imageint *imageflags;
  double **coords,lamda[3],rotmat[3][3];
  double *newcoord;
  double **v = atom->v;
  double t,delx,dely,delz,rsq;

  memory->create(coords,twomol->natoms,3,"bond/react:coords");
  memory->create(imageflags,twomol->natoms,"bond/react:imageflags");

  double *sublo,*subhi;
  if (domain->triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // find current max atom and molecule IDs
  tagint *tag = atom->tag;
  double **x = atom->x;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  tagint maxtag_all,maxmol_all;
  tagint max = 0;
  for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
  MPI_Allreduce(&max,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

  max = 0;
  for (int i = 0; i < nlocal; i++) max = MAX(max,molecule[i]);
  MPI_Allreduce(&max,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

  int dimension = domain->dimension;

  // only proc that owns reacting atom (use ibonding),
  // fits post-reaction template to reaction site, for creating atoms
  int n2superpose = 0;
  for (int j = 0; j < twomol->natoms; j++) {
    if (modify_create_fragid[rxnID] >= 0)
      if (!twomol->fragmentmask[modify_create_fragid[rxnID]][j]) continue;
    if (!create_atoms[j][rxnID] && !delete_atoms[equivalences[j][1][rxnID]][rxnID])
      n2superpose++;
  }

  int ifit = atom->map(my_mega_glove[ibonding[rxnID]+1][iupdate]); // use this local ID to find fitting proc
  Superpose3D<double, double **> superposer(n2superpose);
  int fitroot = 0;
  if (ifit >= 0 && ifit < atom->nlocal) {
    fitroot = comm->me;

    // get 'temperatere' averaged over site, used for created atoms' vels
    t = get_temperature(my_mega_glove,1,iupdate);

    double **xfrozen; // coordinates for the "frozen" target molecule
    double **xmobile; // coordinates for the "mobile" molecule
    memory->create(xfrozen,n2superpose,3,"bond/react:xfrozen");
    memory->create(xmobile,n2superpose,3,"bond/react:xmobile");
    tagint iatom;
    tagint iref = -1; // choose first atom as reference
    int fit_incr = 0;
    for (int j = 0; j < twomol->natoms; j++) {
      if (modify_create_fragid[rxnID] >= 0)
        if (!twomol->fragmentmask[modify_create_fragid[rxnID]][j]) continue;
      int ipre = equivalences[j][1][rxnID]-1; // equiv pre-reaction template index
      if (!create_atoms[j][rxnID] && !delete_atoms[ipre][rxnID]) {
        if (atom->map(my_mega_glove[ipre+1][iupdate]) < 0) {
          error->warning(FLERR," eligible atoms skipped for created-atoms fit on rank {}\n",
                         comm->me);
          continue;
        }
        iatom = atom->map(my_mega_glove[ipre+1][iupdate]);
        if (iref == -1) iref = iatom;
        iatom = domain->closest_image(iref,iatom);
        for (int k = 0; k < 3; k++) {
          xfrozen[fit_incr][k] = x[iatom][k];
          xmobile[fit_incr][k] = twomol->x[j][k];
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
  for (int m = 0; m < twomol->natoms; m++) {
    if (create_atoms[m][rxnID] == 1) {
      // apply optimal rotation/translation for created atom coords
      // also map coords back into simulation box
      if (fitroot == comm->me) {
        MathExtra::matvec(rotmat,twomol->x[m],coords[m]);
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
  if (overlapsq[rxnID] > 0) {
    int abortflag = 0;
    for (int m = 0; m < twomol->natoms; m++) {
      if (create_atoms[m][rxnID] == 1) {
        for (int i = 0; i < nlocal; i++) {
          delx = coords[m][0] - x[i][0];
          dely = coords[m][1] - x[i][1];
          delz = coords[m][2] - x[i][2];
          domain->minimum_image(delx,dely,delz);
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < overlapsq[rxnID]) {
            abortflag = 1;
            break;
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

  // clear ghost count and any ghost bonus data internal to AtomVec
  // same logic as beginning of Comm::exchange()
  // do it now b/c inserting atoms will overwrite ghost atoms
  atom->nghost = 0;
  atom->avec->clear_bonus();

  // check if new atoms are in my sub-box or above it if I am highest proc
  // if so, add atom to my list via create_atom()
  // initialize additional info about the atoms
  // set group mask to "all" plus fix group
  int preID; // new equivalences index
  int add_count = 0;
  for (int m = 0; m < twomol->natoms; m++) {
    if (create_atoms[m][rxnID] == 1) {
      // increase atom count
      add_count++;
      preID = onemol->natoms+add_count;

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
      if (flag) {
        root = comm->me;

        atom->avec->create_atom(twomol->type[m],coords[m]);
        int n = atom->nlocal - 1;
        atom->tag[n] = maxtag_all + add_count;

        // locally update mega_glove
        my_mega_glove[preID][iupdate] = atom->tag[n];

        if (atom->molecule_flag) {
          if (twomol->moleculeflag) {
            atom->molecule[n] = maxmol_all + twomol->molecule[m];
          } else {
            atom->molecule[n] = maxmol_all + 1;
          }
        }

        atom->mask[n] = 1 | groupbit;
        atom->image[n] = imageflags[m];

        // guess a somewhat reasonable initial velocity based on reaction site
        // further control is possible using bond_react_MASTER_group
        // compute |velocity| corresponding to a given temperature t, using specific atom's mass
        double vtnorm = sqrt(t / (force->mvv2e / (dimension * force->boltz)) / atom->mass[twomol->type[m]]);
        v[n][0] = random[rxnID]->uniform();
        v[n][1] = random[rxnID]->uniform();
        v[n][2] = random[rxnID]->uniform();
        double vnorm = sqrt(v[n][0]*v[n][0] + v[n][1]*v[n][1] + v[n][2]*v[n][2]);
        v[n][0] = v[n][0]/vnorm*vtnorm;
        v[n][1] = v[n][1]/vnorm*vtnorm;
        v[n][2] = v[n][2]/vnorm*vtnorm;
        modify->create_attribute(n);
      }
      // globally update mega_glove and equivalences
      MPI_Allreduce(MPI_IN_PLACE,&root,1,MPI_INT,MPI_SUM,world);
      MPI_Bcast(&my_mega_glove[preID][iupdate],1,MPI_LMP_TAGINT,root,world);
      equivalences[m][0][rxnID] = m+1;
      equivalences[m][1][rxnID] = preID;
      reverse_equiv[preID-1][0][rxnID] = preID;
      reverse_equiv[preID-1][1][rxnID] = m+1;
    }
  }

  // reset global natoms here
  // reset atom map elsewhere, after all calls to 'insert_atoms'
  atom->natoms += add_count;
  if (atom->natoms < 0)
    error->all(FLERR,"Too many total atoms");
  maxtag_all += add_count;
  if (maxtag_all >= MAXTAGINT)
    error->all(FLERR,"New atom IDs exceed maximum allowed ID");
  // atom creation successful
  memory->destroy(coords);
  memory->destroy(imageflags);
  return 1;
}

/* ----------------------------------------------------------------------
add equal-style variable to keyword argument list
------------------------------------------------------------------------- */

void FixBondReact::read_variable_keyword(const char *myarg, int keyword, int myrxn)
{
  var_id[keyword][myrxn] = input->variable->find(myarg);
  if (var_id[keyword][myrxn] < 0)
    error->all(FLERR,"Fix bond/react: Variable name {} does not exist",myarg);
  if (!input->variable->equalstyle(var_id[keyword][myrxn]))
    error->all(FLERR,"Fix bond/react: Variable {} is not equal-style",myarg);
  var_flag[keyword][myrxn] = 1;
}

/* ----------------------------------------------------------------------
read map file
------------------------------------------------------------------------- */

void FixBondReact::read_map_file(int myrxn)
{
  int rv;
  char line[MAXLINE],keyword[MAXLINE];
  char *eof,*ptr;

  // skip 1st line of file
  eof = fgets(line,MAXLINE,fp);
  if (eof == nullptr) error->one(FLERR,"Fix bond/react: Unexpected end of superimpose file");

  // read header lines
  // skip blank lines or lines that start with "#"
  // stop when read an unrecognized line

  ncreate = 0;
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
      if (nequivalent != onemol->natoms)
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
      rv = sscanf(line,"%d",&nconstraints[myrxn]);
      if (rv != 1) error->one(FLERR, "Map file header is incorrectly formatted");
      if (maxnconstraints < nconstraints[myrxn]) maxnconstraints = nconstraints[myrxn];
      constraints.resize(maxnconstraints, std::vector<Constraint>(nreacts));
    } else break;
  }

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
      rv = sscanf(line,"%d",&ibonding[myrxn]);
      if (rv != 1) error->one(FLERR, "InitiatorIDs section is incorrectly formatted");
      if (ibonding[myrxn] > onemol->natoms)
        error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
      readline(line);
      rv = sscanf(line,"%d",&jbonding[myrxn]);
      if (rv != 1) error->one(FLERR, "InitiatorIDs section is incorrectly formatted");
      if (jbonding[myrxn] > onemol->natoms)
        error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    } else if (strcmp(keyword,"EdgeIDs") == 0) {
      EdgeIDs(line, myrxn);
    } else if (strcmp(keyword,"Equivalences") == 0) {
      equivflag = 1;
      Equivalences(line, myrxn);
    } else if (strcmp(keyword,"DeleteIDs") == 0) {
      DeleteAtoms(line, myrxn);
    } else if (strcmp(keyword,"CreateIDs") == 0) {
      CreateAtoms(line, myrxn);
    } else if (strcmp(keyword,"ChiralIDs") == 0) {
      ChiralCenters(line, myrxn);
    } else if (strcmp(keyword,"Constraints") == 0) {
      ReadConstraints(line, myrxn);
    } else error->one(FLERR,"Fix bond/react: Unknown section in map file");

    parse_keyword(1,line,keyword);

  }

  // error check
  if (bondflag == 0 || equivflag == 0)
    error->all(FLERR,"Fix bond/react: Map file missing InitiatorIDs or Equivalences section\n");
}

void FixBondReact::EdgeIDs(char *line, int myrxn)
{
  // puts a 1 at edge(edgeID)

  int tmp,rv;
  for (int i = 0; i < nedge; i++) {
    readline(line);
    rv = sscanf(line,"%d",&tmp);
    if (rv != 1) error->one(FLERR, "EdgeIDs section is incorrectly formatted");
    if (tmp > onemol->natoms)
      error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    edge[tmp-1][myrxn] = 1;
  }
}

void FixBondReact::Equivalences(char *line, int myrxn)
{
  int tmp1,tmp2,rv;
  for (int i = 0; i < nequivalent; i++) {
    readline(line);
    rv = sscanf(line,"%d %d",&tmp1,&tmp2);
    if (rv != 2) error->one(FLERR, "Equivalences section is incorrectly formatted");
    if (tmp1 > onemol->natoms || tmp2 > twomol->natoms)
      error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    //equivalences is-> clmn 1: post-reacted, clmn 2: pre-reacted
    equivalences[tmp2-1][0][myrxn] = tmp2;
    equivalences[tmp2-1][1][myrxn] = tmp1;
    //reverse_equiv is-> clmn 1: pre-reacted, clmn 2: post-reacted
    reverse_equiv[tmp1-1][0][myrxn] = tmp1;
    reverse_equiv[tmp1-1][1][myrxn] = tmp2;
  }
}

void FixBondReact::DeleteAtoms(char *line, int myrxn)
{
  int tmp,rv;
  for (int i = 0; i < ndelete; i++) {
    readline(line);
    rv = sscanf(line,"%d",&tmp);
    if (rv != 1) error->one(FLERR, "DeleteIDs section is incorrectly formatted");
    if (tmp > onemol->natoms)
      error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    delete_atoms[tmp-1][myrxn] = 1;
  }
}

void FixBondReact::CreateAtoms(char *line, int myrxn)
{
  create_atoms_flag[myrxn] = 1;
  int tmp,rv;
  for (int i = 0; i < ncreate; i++) {
    readline(line);
    rv = sscanf(line,"%d",&tmp);
    if (rv != 1) error->one(FLERR, "CreateIDs section is incorrectly formatted");
    create_atoms[tmp-1][myrxn] = 1;
  }
  if (twomol->xflag == 0)
    error->one(FLERR,"Fix bond/react: 'Coords' section required in post-reaction template when creating new atoms");
}

void FixBondReact::CustomCharges(int ifragment, int myrxn)
{
  for (int i = 0; i < onemol->natoms; i++)
    if (onemol->fragmentmask[ifragment][i])
      custom_charges[i][myrxn] = 1;
    else
      custom_charges[i][myrxn] = 0;
}

void FixBondReact::ChiralCenters(char *line, int myrxn)
{
  int tmp,rv;
  for (int i = 0; i < nchiral; i++) {
    readline(line);
    rv = sscanf(line,"%d",&tmp);
    if (rv != 1) error->one(FLERR, "ChiralIDs section is incorrectly formatted");
    if (tmp > onemol->natoms)
      error->one(FLERR,"Fix bond/react: Invalid template atom ID in map file");
    chiral_atoms[tmp-1][0][myrxn] = 1;
    if (onemol->xflag == 0)
      error->one(FLERR,"Fix bond/react: Molecule template 'Coords' section required for chiralIDs keyword");
    if ((int) onemol_nxspecial[tmp-1][0] != 4)
      error->one(FLERR,"Fix bond/react: Chiral atoms must have exactly four first neighbors");
    for (int j = 0; j < 4; j++) {
      for (int k = j+1; k < 4; k++) {
        if (onemol->type[onemol_xspecial[tmp-1][j]-1] ==
            onemol->type[onemol_xspecial[tmp-1][k]-1])
          error->one(FLERR,"Fix bond/react: First neighbors of chiral atoms must be of mutually different types");
      }
    }
    // record order of atom types, and coords
    double my4coords[12];
    for (int j = 0; j < 4; j++) {
      chiral_atoms[tmp-1][j+2][myrxn] = onemol->type[onemol_xspecial[tmp-1][j]-1];
      for (int k = 0; k < 3; k++) {
        my4coords[3*j+k] = onemol->x[onemol_xspecial[tmp-1][j]-1][k];
      }
    }
    // get orientation
    chiral_atoms[tmp-1][1][myrxn] = get_chirality(my4coords);
  }
}

void FixBondReact::ReadConstraints(char *line, int myrxn)
{
  int rv;
  double tmp[MAXCONARGS];
  char **strargs,*ptr,*lptr;
  memory->create(strargs,MAXCONARGS,MAXLINE,"bond/react:strargs");
  auto constraint_type = new char[MAXLINE];
  strcpy(constraintstr[myrxn],"("); // string for boolean constraint logic
  for (int i = 0; i < nconstraints[myrxn]; i++) {
    readline(line);
    // find left parentheses, add to constraintstr, and update line
    for (int j = 0; j < (int)strlen(line); j++) {
      if (line[j] == '(') strcat(constraintstr[myrxn],"(");
      if (isalpha(line[j])) {
        line = line + j;
        break;
      }
    }
    // 'C' indicates where to sub in next constraint
    strcat(constraintstr[myrxn],"C");
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
      if (lptr[j] == ')') strcat(constraintstr[myrxn],")");
    // find logic symbols, and trim line via ptr
    if ((ptr = strstr(lptr,"&&"))) {
      strcat(constraintstr[myrxn],"&&");
      *ptr = '\0';
    } else if ((ptr = strstr(lptr,"||"))) {
      strcat(constraintstr[myrxn],"||");
      *ptr = '\0';
    } else if (i+1 < nconstraints[myrxn]) {
      strcat(constraintstr[myrxn],"&&");
    }
    if ((ptr = strchr(lptr,')')))
      *ptr = '\0';
    rv = sscanf(line,"%s",constraint_type);
    if (rv != 1) error->one(FLERR, "Constraints section is incorrectly formatted");
    if (strcmp(constraint_type,"distance") == 0) {
      constraints[i][myrxn].type = DISTANCE;
      rv = sscanf(line,"%*s %s %s %lg %lg",strargs[0],strargs[1],&tmp[0],&tmp[1]);
      if (rv != 4) error->one(FLERR, "Distance constraint is incorrectly formatted");
      readID(strargs[0], i, myrxn, 0);
      readID(strargs[1], i, myrxn, 1);
      // cutoffs
      constraints[i][myrxn].par[0] = tmp[0]*tmp[0]; // using square of distance
      constraints[i][myrxn].par[1] = tmp[1]*tmp[1];
    } else if (strcmp(constraint_type,"angle") == 0) {
      constraints[i][myrxn].type = ANGLE;
      rv = sscanf(line,"%*s %s %s %s %lg %lg",strargs[0],strargs[1],strargs[2],&tmp[0],&tmp[1]);
      if (rv != 5) error->one(FLERR, "Angle constraint is incorrectly formatted");
      readID(strargs[0], i, myrxn, 0);
      readID(strargs[1], i, myrxn, 1);
      readID(strargs[2], i, myrxn, 2);
      constraints[i][myrxn].par[0] = tmp[0]/180.0 * MY_PI;
      constraints[i][myrxn].par[1] = tmp[1]/180.0 * MY_PI;
    } else if (strcmp(constraint_type,"dihedral") == 0) {
      constraints[i][myrxn].type = DIHEDRAL;
      tmp[2] = 181.0; // impossible range
      tmp[3] = 182.0;
      rv = sscanf(line,"%*s %s %s %s %s %lg %lg %lg %lg",strargs[0],strargs[1],
             strargs[2],strargs[3],&tmp[0],&tmp[1],&tmp[2],&tmp[3]);
      if (!(rv == 6 || rv == 8)) error->one(FLERR, "Dihedral constraint is incorrectly formatted");
      readID(strargs[0], i, myrxn, 0);
      readID(strargs[1], i, myrxn, 1);
      readID(strargs[2], i, myrxn, 2);
      readID(strargs[3], i, myrxn, 3);
      constraints[i][myrxn].par[0] = tmp[0]/180.0 * MY_PI;
      constraints[i][myrxn].par[1] = tmp[1]/180.0 * MY_PI;
      constraints[i][myrxn].par[2] = tmp[2]/180.0 * MY_PI;
      constraints[i][myrxn].par[3] = tmp[3]/180.0 * MY_PI;
    } else if (strcmp(constraint_type,"arrhenius") == 0) {
      constraints[i][myrxn].type = ARRHENIUS;
      constraints[i][myrxn].par[0] = narrhenius++;
      rv = sscanf(line,"%*s %lg %lg %lg %lg",&tmp[0],&tmp[1],&tmp[2],&tmp[3]);
      if (rv != 4) error->one(FLERR, "Arrhenius constraint is incorrectly formatted");
      constraints[i][myrxn].par[1] = tmp[0];
      constraints[i][myrxn].par[2] = tmp[1];
      constraints[i][myrxn].par[3] = tmp[2];
      constraints[i][myrxn].par[4] = tmp[3];
    } else if (strcmp(constraint_type,"rmsd") == 0) {
      constraints[i][myrxn].type = RMSD;
      strcpy(strargs[0],"0");
      rv = sscanf(line,"%*s %lg %s",&tmp[0],strargs[0]);
      if (!(rv == 1 || rv == 2)) error->one(FLERR, "RMSD constraint is incorrectly formatted");
      constraints[i][myrxn].par[0] = tmp[0]; // RMSDmax
      constraints[i][myrxn].id[0] = -1; // optional molecule fragment
      if (isalpha(strargs[0][0])) {
        int ifragment = onemol->findfragment(strargs[0]);
        if (ifragment < 0) error->one(FLERR,"Fix bond/react: Molecule fragment does not exist");
        else constraints[i][myrxn].id[0] = ifragment;
      }
    } else if (strcmp(constraint_type,"custom") == 0) {
      constraints[i][myrxn].type = CUSTOM;
      std::vector<std::string> args = utils::split_words(line);
      constraints[i][myrxn].str = args[1];
    } else error->one(FLERR,"Fix bond/react: Illegal constraint type in 'Constraints' section of map file");
  }
  strcat(constraintstr[myrxn],")"); // close boolean constraint logic string
  delete [] constraint_type;
  memory->destroy(strargs);
}

/* ----------------------------------------------------------------------
if ID starts with character, assume it is a pre-reaction molecule fragment ID
otherwise, it is a pre-reaction atom ID
---------------------------------------------------------------------- */

void FixBondReact::readID(char *strarg, int iconstr, int myrxn, int i)
{
  if (isalpha(strarg[0])) {
    constraints[iconstr][myrxn].idtype[i] = FRAG; // fragment vs. atom ID flag
    int ifragment = onemol->findfragment(strarg);
    if (ifragment < 0)
      error->one(FLERR,"Fix bond/react: Molecule fragment {} does not exist", strarg);
    constraints[iconstr][myrxn].id[i] = ifragment;
  } else {
    constraints[iconstr][myrxn].idtype[i] = ATOM; // fragment vs. atom ID flag
    int iatom = utils::inumeric(FLERR, strarg, true, lmp);
    if (iatom > onemol->natoms)
      error->one(FLERR,"Fix bond/react: Invalid template atom ID {} in map file", strarg);
    constraints[iconstr][myrxn].id[i] = iatom;
  }
}

void FixBondReact::open(char *file)
{
  fp = fopen(file,"r");
  if (fp == nullptr) error->one(FLERR, "Fix bond/react: Cannot open map file {}", file);
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
  return (double) reaction_count_total[n];

}

/* ---------------------------------------------------------------------- */

void FixBondReact::post_integrate_respa(int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_integrate();
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
    if (closeneigh[rxnID] != 0)
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
    if (closeneigh[rxnID] != 0) {
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
  int revision = 1;
  set[0].nreacts = nreacts;
  set[0].max_rate_limit_steps = max_rate_limit_steps;

  for (int i = 0; i < nreacts; i++) {
    set[i].reaction_count_total = reaction_count_total[i];

    strncpy(set[i].rxn_name,rxn_name[i],MAXNAME-1);
    set[i].rxn_name[MAXNAME-1] = '\0';
  }

  int rbufcount = max_rate_limit_steps*nreacts;
  int *rbuf;
  if (rbufcount) {
    memory->create(rbuf,rbufcount,"bond/react:rbuf");
    memcpy(rbuf,&store_rxn_count[0][0],sizeof(int)*rbufcount);
  }

  if (comm->me == 0) {
    int size = nreacts*sizeof(Set)+(rbufcount+1)*sizeof(int);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(&revision,sizeof(int),1,fp);
    fwrite(set,sizeof(Set),nreacts,fp);
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
  int n,revision,r_nreacts,r_max_rate_limit_steps,ibufcount,n2cpy;
  int **ibuf;

  n = 0;
  if (lmp->restart_ver > utils::date2num("3 Nov 2022")) revision = buf[n++];
  else revision = 0;

  Set *set_restart = (Set *) &buf[n*sizeof(int)];
  r_nreacts = set_restart[0].nreacts;

  n2cpy = 0;
  if (revision > 0) {
    r_max_rate_limit_steps = set_restart[0].max_rate_limit_steps;
    if (r_max_rate_limit_steps > 0) {
      ibufcount = r_max_rate_limit_steps*r_nreacts;
      memory->create(ibuf,r_max_rate_limit_steps,r_nreacts,"bond/react:ibuf");
      memcpy(&ibuf[0][0],&buf[sizeof(int)+r_nreacts*sizeof(Set)],sizeof(int)*ibufcount);
      n2cpy = r_max_rate_limit_steps;
    }
  }

  if (max_rate_limit_steps < n2cpy) n2cpy = max_rate_limit_steps;
  for (int i = 0; i < r_nreacts; i++) {
    for (int j = 0; j < nreacts; j++) {
      if (strcmp(set_restart[i].rxn_name,rxn_name[j]) == 0) {
        reaction_count_total[j] = set_restart[i].reaction_count_total;
        // read rate_limit restart information
        for (int k = 0; k < n2cpy; k++)
          store_rxn_count[k][j] = ibuf[k][i];
      }
    }
  }
  if (revision > 0 && r_max_rate_limit_steps > 0) memory->destroy(ibuf);
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
