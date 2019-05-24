/* ----------------------------------------------------------------------
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov

Copyright (2003) Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software.  This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
Contributing Author: Jacob Gissinger (jacob.gissinger@colorado.edu)
------------------------------------------------------------------------- */

#include "fix_bond_react.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "update.h"
#include "modify.h"
#include "respa.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "random_mars.h"
#include "molecule.h"
#include "group.h"
#include "citeme.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#include <algorithm>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

static const char cite_fix_bond_react[] =
  "fix bond/react:\n\n"
  "@Article{Gissinger17,\n"
  " author = {J. R. Gissinger, B. D. Jensen, K. E. Wise},\n"
  " title = {Modeling chemical reactions in classical molecular dynamics simulations},\n"
  " journal = {Polymer},\n"
  " year =    2017,\n"
  " volume =  128,\n"
  " pages =   {211--217}\n"
  "}\n\n";

#define BIG 1.0e20
#define DELTA 16
#define MAXGUESS 20 // max # of guesses allowed by superimpose algorithm
#define MAXCONARGS 7 // max # of arguments for any type of constraint + rxnID

// various statuses of superimpose algorithm:
// ACCEPT: site successfully matched to pre-reacted template
// REJECT: site does not match pre-reacted template
// PROCEED: normal execution (non-guessing mode)
// CONTINUE: a neighbor has been assigned, skip to next neighbor
// GUESSFAIL: a guess has failed (if no more restore points, status = 'REJECT')
// RESTORE: restore mode, load most recent restore point
enum{ACCEPT,REJECT,PROCEED,CONTINUE,GUESSFAIL,RESTORE};

// types of available reaction constraints
enum{DISTANCE,ANGLE,ARRHENIUS};

/* ---------------------------------------------------------------------- */

FixBondReact::FixBondReact(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_bond_react);

  fix1 = NULL;
  fix2 = NULL;
  fix3 = NULL;

  if (narg < 8) error->all(FLERR,"Illegal fix bond/react command: "
                           "too few arguments");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  newton_bond = force->newton_bond;

  restart_global = 1;
  attempted_rxn = 0;
  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  global_freq = 1;
  extvector = 0;
  rxnID = 0;
  nconstraints = 0;
  narrhenius = 0;
  status = PROCEED;

  nxspecial = NULL;
  onemol_nxspecial = NULL;
  twomol_nxspecial = NULL;
  xspecial = NULL;
  onemol_xspecial = NULL;
  twomol_xspecial = NULL;

  // these group names are reserved for use exclusively by bond/react
  master_group = (char *) "bond_react_MASTER_group";

  // by using fixed group names, only one instance of fix bond/react is allowed.
  if (modify->find_fix_by_style("^bond/react") != -1)
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
  int num_common_keywords = 1;
  for (int m = 0; m < num_common_keywords; m++) {
    if (strcmp(arg[iarg],"stabilization") == 0) {
      if (strcmp(arg[iarg+1],"no") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'stabilization' keyword has too few arguments");
        iarg += 2;
      }
      if (strcmp(arg[iarg+1],"yes") == 0) {
        if (iarg+4 > narg) error->all(FLERR,"Illegal fix bond/react command:"
                                      "'stabilization' keyword has too few arguments");
        int n = strlen(arg[iarg+2]) + 1;
        exclude_group = new char[n];
        strcpy(exclude_group,arg[iarg+2]);
        stabilization_flag = 1;
        nve_limit_xmax = arg[iarg+3];
        iarg += 4;
      }
    } else if (strcmp(arg[iarg],"react") == 0) {
      break;
    } else error->all(FLERR,"Illegal fix bond/react command: unknown keyword");
  }

  // set up common variables as vectors of length 'nreacts'
  // nevery, cutoff, onemol, twomol, superimpose file

  // this looks excessive
  // the price of vectorization (all reactions in one command)?
  memory->create(rxn_name,nreacts,MAXLINE,"bond/react:rxn_name");
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
  memory->create(stabilize_steps_flag,nreacts,"bond/react:stabilize_steps_flag");
  memory->create(update_edges_flag,nreacts,"bond/react:update_edges_flag");
  memory->create(constraints,1,MAXCONARGS,"bond/react:constraints");
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
    stabilize_steps_flag[i] = 0;
    update_edges_flag[i] = 0;
    // set default limit duration to 60 timesteps
    limit_duration[i] = 60;
    reaction_count[i] = 0;
    local_rxn_count[i] = 0;
    ghostly_rxn_count[i] = 0;
    reaction_count_total[i] = 0;
  }

  char **files;
  files = new char*[nreacts];

  for (int rxn = 0; rxn < nreacts; rxn++) {

    if (strcmp(arg[iarg],"react") != 0) error->all(FLERR,"Illegal fix bond/react command: "
                                                   "'react' or 'stabilization' has incorrect arguments");

    iarg++;

    int n = strlen(arg[iarg]) + 1;
    if (n > MAXLINE) error->all(FLERR,"Reaction name (react-ID) is too long (limit: 256 characters)");
    strncpy(rxn_name[rxn],arg[iarg++],n);

    int igroup = group->find(arg[iarg++]);
    if (igroup == -1) error->all(FLERR,"Could not find fix group ID");
    groupbits[rxn] = group->bitmask[igroup];

    nevery[rxn] = force->inumeric(FLERR,arg[iarg++]);
    if (nevery[rxn] <= 0) error->all(FLERR,"Illegal fix bond/react command: "
                                     "'Nevery' must be a positive integer");

    double cutoff = force->numeric(FLERR,arg[iarg++]);
    if (cutoff < 0.0) error->all(FLERR,"Illegal fix bond/react command: "
                                 "'Rmin' cannot be negative");
    cutsq[rxn][0] = cutoff*cutoff;

    cutoff = force->numeric(FLERR,arg[iarg++]);
    if (cutoff < 0.0) error->all(FLERR,"Illegal fix bond/react command:"
                                 "'Rmax' cannot be negative");
    cutsq[rxn][1] = cutoff*cutoff;

    unreacted_mol[rxn] = atom->find_molecule(arg[iarg++]);
    if (unreacted_mol[rxn] == -1) error->all(FLERR,"Unreacted molecule template ID for "
                                             "fix bond/react does not exist");
    reacted_mol[rxn] = atom->find_molecule(arg[iarg++]);
    if (reacted_mol[rxn] == -1) error->all(FLERR,"Reacted molecule template ID for "
                                           "fix bond/react does not exist");

    //read superimpose file
    files[rxn] = new char[strlen(arg[iarg])+1];
    strcpy(files[rxn],arg[iarg]);
    iarg++;

    while (iarg < narg && strcmp(arg[iarg],"react") != 0 ) {
      if (strcmp(arg[iarg],"prob") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'prob' keyword has too few arguments");
        fraction[rxn] = force->numeric(FLERR,arg[iarg+1]);
        seed[rxn] = force->inumeric(FLERR,arg[iarg+2]);
        if (fraction[rxn] < 0.0 || fraction[rxn] > 1.0)
          error->all(FLERR,"Illegal fix bond/react command: "
                     "probability fraction must between 0 and 1, inclusive");
        if (seed[rxn] <= 0) error->all(FLERR,"Illegal fix bond/react command: "
                                       "probability seed must be positive");
        iarg += 3;
      } else if (strcmp(arg[iarg],"max_rxn") == 0) {
              if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                                        "'max_rxn' has too few arguments");
              max_rxn[rxn] = force->inumeric(FLERR,arg[iarg+1]);
              if (max_rxn[rxn] < 0) error->all(FLERR,"Illegal fix bond/react command: "
                                                                 "'max_rxn' cannot be negative");
              iarg += 2;
      } else if (strcmp(arg[iarg],"stabilize_steps") == 0) {
        if (stabilization_flag == 0) error->all(FLERR,"Stabilize_steps keyword "
                                                "used without stabilization keyword");
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'stabilize_steps' has too few arguments");
        limit_duration[rxn] = force->numeric(FLERR,arg[iarg+1]);
        stabilize_steps_flag[rxn] = 1;
        iarg += 2;
      } else if (strcmp(arg[iarg],"update_edges") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/react command: "
                                      "'update_edges' has too few arguments");
        if (strcmp(arg[iarg+1],"none") == 0) update_edges_flag[rxn] = 0;
        else if (strcmp(arg[iarg+1],"charges") == 0) update_edges_flag[rxn] = 1;
        else if (strcmp(arg[iarg+1],"custom") == 0) update_edges_flag[rxn] = 2;
        else error->all(FLERR,"Illegal value for 'update_edges' keyword'");
        iarg += 2;
      } else error->all(FLERR,"Illegal fix bond/react command: unknown keyword");
    }
  }

  max_natoms = 0; // the number of atoms in largest molecule template
  for (int myrxn = 0; myrxn < nreacts; myrxn++) {
    twomol = atom->molecules[reacted_mol[myrxn]];
    max_natoms = MAX(max_natoms,twomol->natoms);
  }

  memory->create(equivalences,max_natoms,2,nreacts,"bond/react:equivalences");
  memory->create(reverse_equiv,max_natoms,2,nreacts,"bond/react:reverse_equiv");
  memory->create(edge,max_natoms,nreacts,"bond/react:edge");
  memory->create(landlocked_atoms,max_natoms,nreacts,"bond/react:landlocked_atoms");
  memory->create(custom_edges,max_natoms,nreacts,"bond/react:custom_edges");
  memory->create(delete_atoms,max_natoms,nreacts,"bond/react:delete_atoms");

  for (int j = 0; j < nreacts; j++)
    for (int i = 0; i < max_natoms; i++) {
      edge[i][j] = 0;
      if (update_edges_flag[j] == 1) custom_edges[i][j] = 1;
      else custom_edges[i][j] = 0;
      delete_atoms[i][j] = 0;
    }

  // read all map files afterward
  for (int i = 0; i < nreacts; i++) {
    open(files[i]);
    onemol = atom->molecules[unreacted_mol[i]];
    twomol = atom->molecules[reacted_mol[i]];
    onemol->check_attributes(0);
    twomol->check_attributes(0);
    if (onemol->natoms != twomol->natoms)
      error->all(FLERR,"Bond/react: Reaction templates must contain the same number of atoms");
    get_molxspecials();
    read(i);
    fclose(fp);
    iatomtype[i] = onemol->type[ibonding[i]-1];
    jatomtype[i] = onemol->type[jbonding[i]-1];
    find_landlocked_atoms(i);
  }

  // initialize Marsaglia RNG with processor-unique seed (Arrhenius prob)

  rrhandom = new class RanMars*[narrhenius];
  int tmp = 0;
  for (int i = 0; i < nconstraints; i++) {
    if (constraints[i][1] == ARRHENIUS) {
      rrhandom[tmp++] = new RanMars(lmp,(int) constraints[i][6] + me);
    }
  }

  for (int i = 0; i < nreacts; i++) {
    delete [] files[i];
  }
  delete [] files;

  if (atom->molecular != 1)
    error->all(FLERR,"Bond/react: Cannot use fix bond/react with non-molecular systems");

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

  random = new class RanMars*[nreacts];
  for (int i = 0; i < nreacts; i++) {
    random[i] = new RanMars(lmp,seed[i] + me);
  }

  // set comm sizes needed by this fix
  // forward is big due to comm of broken bonds and 1-2 neighbors

  comm_forward = MAX(2,2+atom->maxspecial);
  comm_reverse = 2;

  // allocate arrays local to this fix
  nmax = 0;
  partner = finalpartner = NULL;
  distsq = NULL;
  probability = NULL;
  maxcreate = 0;
  created = NULL;
  ncreate = NULL;
  allncreate = 0;
  local_num_mega = 0;
  ghostly_num_mega = 0;
  restore =  NULL;

  // zero out stats
  global_megasize = 0;
  avail_guesses = 0;
  glove_counter = 0;
  guess_branch = new int[MAXGUESS]();
  pioneer_count = new int[max_natoms];
  local_mega_glove = NULL;
  ghostly_mega_glove = NULL;
  global_mega_glove = NULL;

  // these are merely loop indices that became important
  pion = neigh = trace = 0;

  id_fix1 = NULL;
  id_fix2 = NULL;
  id_fix3 = NULL;
  statted_id = NULL;
  custom_exclude_flag = 0;

  // used to store restart info
  set = new Set[nreacts];
  memset(set,0,nreacts*sizeof(Set));
}

/* ---------------------------------------------------------------------- */

FixBondReact::~FixBondReact()
{
  for (int i = 0; i < nreacts; i++) {
    delete random[i];
  }
  delete [] random;

  memory->destroy(partner);
  memory->destroy(finalpartner);
  memory->destroy(ncreate);
  memory->destroy(distsq);
  memory->destroy(probability);
  memory->destroy(created);
  memory->destroy(edge);
  memory->destroy(equivalences);
  memory->destroy(reverse_equiv);
  memory->destroy(landlocked_atoms);
  memory->destroy(custom_edges);
  memory->destroy(delete_atoms);

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
  memory->destroy(stabilize_steps_flag);
  memory->destroy(update_edges_flag);

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
    memory->destroy(local_mega_glove);
    memory->destroy(ghostly_mega_glove);
  }

  memory->destroy(global_mega_glove);

  if (stabilization_flag == 1) {
    // check nfix in case all fixes have already been deleted
    if (id_fix1 && modify->nfix) modify->delete_fix(id_fix1);
    delete [] id_fix1;

    if (id_fix3 && modify->nfix) modify->delete_fix(id_fix3);
    delete [] id_fix3;
  }

  if (id_fix2 && modify->nfix) modify->delete_fix(id_fix2);
  delete [] id_fix2;

  delete [] statted_id;
  delete [] guess_branch;
  delete [] pioneer_count;
  delete [] set;

  if (group) {
    char **newarg;
    newarg = new char*[2];
    newarg[0] = master_group;
    newarg[1] = (char *) "delete";
    group->assign(2,newarg);
    if (stabilization_flag == 1) {
      newarg[0] = exclude_group;
      group->assign(2,newarg);
      delete [] exclude_group;
    }
    delete [] newarg;
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
  int len = strlen("bond_react_props_internal") + 1;
  id_fix2 = new char[len];
  strcpy(id_fix2,"bond_react_props_internal");

  int ifix = modify->find_fix(id_fix2);
  if (ifix == -1) {
    char **newarg = new char*[7];
    newarg[0] = (char *) "bond_react_props_internal";
    newarg[1] = (char *) "all"; // group ID is ignored
    newarg[2] = (char *) "property/atom";
    newarg[3] = (char *) "i_limit_tags";
    newarg[4] = (char *) "i_react_tags";
    newarg[5] = (char *) "ghost";
    newarg[6] = (char *) "yes";
    modify->add_fix(7,newarg);
    delete [] newarg;
  }

  // create master_group if not already existing
  // NOTE: limit_tags and react_tags automaticaly intitialized to zero (unless read from restart)
  group->find_or_create(master_group);
  char **newarg;
  newarg = new char*[5];
  newarg[0] = master_group;
  newarg[1] = (char *) "dynamic";
  newarg[2] = (char *) "all";
  newarg[3] = (char *) "property";
  newarg[4] = (char *) "limit_tags";
  group->assign(5,newarg);
  delete [] newarg;

  if (stabilization_flag == 1) {
    int igroup = group->find(exclude_group);
    // create exclude_group if not already existing, or use as parent group if static
    if (igroup == -1 || group->dynamic[igroup] == 0) {
      // create stabilization per-atom property
      len = strlen("bond_react_stabilization_internal") + 1;
      id_fix3 = new char[len];
      strcpy(id_fix3,"bond_react_stabilization_internal");

      ifix = modify->find_fix(id_fix3);
      if (ifix == -1) {
        char **newarg = new char*[6];
        newarg[0] = (char *) id_fix3;
        newarg[1] = (char *) "all"; // group ID is ignored
        newarg[2] = (char *) "property/atom";
        newarg[3] = (char *) "i_statted_tags";
        newarg[4] = (char *) "ghost";
        newarg[5] = (char *) "yes";
        modify->add_fix(6,newarg);
        fix3 = modify->fix[modify->nfix-1];
        delete [] newarg;
      }

      len = strlen("statted_tags") + 1;
      statted_id = new char[len];
      strcpy(statted_id,"statted_tags");

      // if static group exists, use as parent group
      // also, rename dynamic exclude_group by appending '_REACT'
      char *exclude_PARENT_group;
      int n = strlen(exclude_group) + 1;
      exclude_PARENT_group = new char[n];
      strcpy(exclude_PARENT_group,exclude_group);
      n += strlen("_REACT");
      delete [] exclude_group;
      exclude_group = new char[n];
      strcpy(exclude_group,exclude_PARENT_group);
      strcat(exclude_group,"_REACT");

      group->find_or_create(exclude_group);
      char **newarg;
      newarg = new char*[5];
      newarg[0] = exclude_group;
      newarg[1] = (char *) "dynamic";
      if (igroup == -1) newarg[2] = (char *) "all";
      else newarg[2] = (char *) exclude_PARENT_group;
      newarg[3] = (char *) "property";
      newarg[4] = (char *) "statted_tags";
      group->assign(5,newarg);
      delete [] newarg;
      delete [] exclude_PARENT_group;

      // on to statted_tags (system-wide thermostat)
      // intialize per-atom statted_flags to 1
      // (only if not already initialized by restart)
      if (fix3->restart_reset != 1) {
        int flag;
        int index = atom->find_custom("statted_tags",flag);
        int *i_statted_tags = atom->ivector[index];

        for (int i = 0; i < atom->nlocal; i++)
          i_statted_tags[i] = 1;
      }
    } else {
        // sleeping code, for future capabilities
        custom_exclude_flag = 1;
        // first we have to find correct fix group reference
        int n = strlen("GROUP_") + strlen(exclude_group) + 1;
        char *fix_group = new char[n];
        strcpy(fix_group,"GROUP_");
        strcat(fix_group,exclude_group);
        int ifix = modify->find_fix(fix_group);
        Fix *fix = modify->fix[ifix];
        delete [] fix_group;

        // this returns names of corresponding property
        int unused;
        char * idprop;
        idprop = (char *) fix->extract("property",unused);
        if (idprop == NULL)
          error->all(FLERR,"Exclude group must be a per-atom property group");

        len = strlen(idprop) + 1;
        statted_id = new char[len];
        strcpy(statted_id,idprop);

        // intialize per-atom statted_tags to 1
        // need to correct for smooth restarts
        //int flag;
        //int index = atom->find_custom(statted_id,flag);
        //int *i_statted_tags = atom->ivector[index];
        //for (int i = 0; i < atom->nlocal; i++)
        //  i_statted_tags[i] = 1;
      }


    // let's create a new nve/limit fix to limit newly reacted atoms
    len = strlen("bond_react_MASTER_nve_limit") + 1;
    id_fix1 = new char[len];
    strcpy(id_fix1,"bond_react_MASTER_nve_limit");

    ifix = modify->find_fix(id_fix1);

    if (ifix == -1) {
      char **newarg = new char*[4];
      newarg[0] = id_fix1;
      newarg[1] = master_group;
      newarg[2] = (char *) "nve/limit";
      newarg[3] = nve_limit_xmax;
      modify->add_fix(4,newarg);
      fix1 = modify->fix[modify->nfix-1];
      delete [] newarg;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBondReact::init()
{

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // check cutoff for iatomtype,jatomtype
  for (int i = 0; i < nreacts; i++) {
    if (force->pair == NULL || cutsq[i][1] > force->pair->cutsq[iatomtype[i]][jatomtype[i]])
      error->all(FLERR,"Bond/react: Fix bond/react cutoff is longer than pairwise cutoff");
  }

  // need a half neighbor list, built every Nevery steps
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->occasional = 1;

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
  // check if any reactions could occur on this timestep
  int nevery_check = 1;
  for (int i = 0; i < nreacts; i++) {
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
    memory->destroy(ncreate);
    memory->destroy(probability);
    nmax = atom->nmax;
    memory->create(partner,nmax,"bond/react:partner");
    memory->create(finalpartner,nmax,"bond/react:finalpartner");
    memory->create(distsq,nmax,2,"bond/react:distsq");
    memory->create(ncreate,nreacts,"bond/react:ncreate");
    memory->create(probability,nmax,"bond/react:probability");
  }

  // reset create counts
  for (int i = 0; i < nreacts; i++) {
    ncreate[i] = 0;
  }

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
    if (max_rxn[rxnID] <= reaction_count_total[rxnID]) continue;
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
      if (force->newton_pair) comm->reverse_comm_fix(this);
    } else {
      close_partner();
      commflag = 2;
      comm->reverse_comm_fix(this);
    }

    // each atom now knows its winning partner
    // for prob check, generate random value for each atom with a bond partner
    // forward comm of partner and random value, so ghosts have it

    if (fraction[rxnID] < 1.0) {
      for (int i = 0; i < nlocal; i++)
      if (partner[i]) probability[i] = random[rxnID]->uniform();
    }

    commflag = 2;
    comm->forward_comm_fix(this,2);

    // consider for reaction:
    // only if both atoms list each other as winning bond partner
    //   and probability constraint is satisfied
    // if other atom is owned by another proc, it should do same thing

    int temp_ncreate = 0;
    for (int i = 0; i < nlocal; i++) {
      if (partner[i] == 0) {
        continue;
      }

      j = atom->map(partner[i]);
      if (partner[j] != tag[i]) {
        continue;
      }

      // apply probability constraint using RN for atom with smallest ID

      if (fraction[rxnID] < 1.0) {
        if (tag[i] < tag[j]) {
          if (probability[i] >= fraction[rxnID]) continue;
        } else {
          if (probability[j] >= fraction[rxnID]) continue;
        }
      }

      // store final created bond partners and count the rxn possibility once

      finalpartner[i] = tag[j];
      finalpartner[j] = tag[i];

      if (tag[i] < tag[j]) temp_ncreate++;
    }

    // cycle loop if no even eligible bonding atoms were found (on any proc)
    int some_chance;
    MPI_Allreduce(&temp_ncreate,&some_chance,1,MPI_INT,MPI_SUM,world);
    if (!some_chance) continue;

    // communicate final partner

    commflag = 3;
    comm->forward_comm_fix(this);

    // add instance to 'created' only if this processor
    // owns the atoms with smaller global ID
    // NOTE: we no longer care about ghost-ghost instances as bond/create did
    // this is because we take care of updating topology later (and differently)
    for (int i = 0; i < nlocal; i++) {

      if (finalpartner[i] == 0) continue;

      j = atom->map(finalpartner[i]);
      // if (j < 0 || tag[i] < tag[j]) {
      if (tag[i] < tag[j]) { //atom->map(std::min(tag[i],tag[j])) <= nlocal &&
        if (ncreate[rxnID] == maxcreate) {
          maxcreate += DELTA;
          // third column of 'created': bond/react integer ID
          memory->grow(created,maxcreate,2,nreacts,"bond/react:created");
        }
        // to ensure types remain in same order
        // unnecessary now taken from reaction map file
        if (iatomtype[rxnID] == type[i]) {
          created[ncreate[rxnID]][0][rxnID] = tag[i];
          created[ncreate[rxnID]][1][rxnID] = finalpartner[i];
        } else {
          created[ncreate[rxnID]][0][rxnID] = finalpartner[i];
          created[ncreate[rxnID]][1][rxnID] = tag[i];
        }
        ncreate[rxnID]++;
      }
    }
  }

  // break loop if no even eligible bonding atoms were found (on any proc)
  int some_chance;

  allncreate = 0;
  for (int i = 0; i < nreacts; i++)
    allncreate += ncreate[i];

  MPI_Allreduce(&allncreate,&some_chance,1,MPI_INT,MPI_SUM,world);
  if (!some_chance) {
    unlimit_bond();
    return;
  }

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
  int flag;
  int index1 = atom->find_custom("limit_tags",flag);
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
  int flag;
  int index1 = atom->find_custom("limit_tags",flag);
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

      delx = x[i1][0] - x[i2][0];
      dely = x[i1][1] - x[i2][1];
      delz = x[i1][2] - x[i2][2];
      domain->minimum_image(delx,dely,delz); // ghost location fix
      rsq = delx*delx + dely*dely + delz*delz;
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
    memory->destroy(local_mega_glove);
    memory->destroy(ghostly_mega_glove);
  }

  memory->create(glove,max_natoms,2,"bond/react:glove");
  memory->create(restore_pt,MAXGUESS,4,"bond/react:restore_pt");
  memory->create(pioneers,max_natoms,"bond/react:pioneers");
  memory->create(restore,max_natoms,MAXGUESS,"bond/react:restore");
  memory->create(local_mega_glove,max_natoms+1,allncreate,"bond/react:local_mega_glove");
  memory->create(ghostly_mega_glove,max_natoms+1,allncreate,"bond/react:ghostly_mega_glove");

  attempted_rxn = 1;

  for (int i = 0; i < max_natoms+1; i++) {
    for (int j = 0; j < allncreate; j++) {
      local_mega_glove[i][j] = 0;
      ghostly_mega_glove[i][j] = 0;
    }
  }

  // let's finally begin the superimpose loop
  for (rxnID = 0; rxnID < nreacts; rxnID++) {
    for (lcl_inst = 0; lcl_inst < ncreate[rxnID]; lcl_inst++) {

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
      glove[myibonding-1][1] = created[lcl_inst][0][rxnID];
      glove_counter++;
      glove[myjbonding-1][0] = myjbonding;
      glove[myjbonding-1][1] = created[lcl_inst][1][rxnID];
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
             check_constraints() ) {
          status = ACCEPT;
          glove_ghostcheck();
        } else
          status = REJECT;
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

        if (status == ACCEPT && check_constraints()) { // reaction site found successfully!
          glove_ghostcheck();
        }
        hang_catch++;
        // let's go ahead and catch the simplest of hangs
        //if (hang_catch > onemol->natoms*4)
        if (hang_catch > atom->nlocal*30) {
          error->one(FLERR,"Bond/react: Excessive iteration of superimpose algorithm");
        }
      }
    }
  }

  global_megasize = 0;

  ghost_glovecast(); // consolidate all mega_gloves to all processors
  dedup_mega_gloves(0); // make sure atoms aren't added to more than one reaction

  MPI_Allreduce(&local_rxn_count[0],&reaction_count[0],nreacts,MPI_INT,MPI_SUM,world);

  if (me == 0)
    for (int i = 0; i < nreacts; i++)
      reaction_count_total[i] += reaction_count[i] + ghostly_rxn_count[i];

  MPI_Bcast(&reaction_count_total[0], nreacts, MPI_INT, 0, world);

  // check if we overstepped our reaction limit
  for (int i = 0; i < nreacts; i++) {
    if (reaction_count_total[i] > max_rxn[i]) {
      // let's randomly choose rxns to skip, unbiasedly from local and ghostly
      int *local_rxncounts;
      int *all_localskips;
      memory->create(local_rxncounts,nprocs,"bond/react:local_rxncounts");
      memory->create(all_localskips,nprocs,"bond/react:all_localskips");
      MPI_Gather(&local_rxn_count[i],1,MPI_INT,local_rxncounts,1,MPI_INT,0,world);
      if (me == 0) {
        int overstep = reaction_count_total[i] - max_rxn[i];
        int delta_rxn = reaction_count[i] + ghostly_rxn_count[i];
        int *rxn_by_proc;
        memory->create(rxn_by_proc,delta_rxn,"bond/react:rxn_by_proc");
        for (int j = 0; j < delta_rxn; j++)
          rxn_by_proc[j] = -1; // corresponds to ghostly
        int itemp = 0;
        for (int j = 0; j < nprocs; j++)
          for (int k = 0; k < local_rxncounts[j]; k++)
            rxn_by_proc[itemp++] = j;
        std::random_shuffle(&rxn_by_proc[0],&rxn_by_proc[delta_rxn]);
        for (int j = 0; j < nprocs; j++)
          all_localskips[j] = 0;
        nghostlyskips[i] = 0;
        for (int j = 0; j < overstep; j++) {
          if (rxn_by_proc[j] == -1) nghostlyskips[i]++;
          else all_localskips[rxn_by_proc[j]]++;
        }
        memory->destroy(rxn_by_proc);
      }
      reaction_count_total[i] = max_rxn[i];
      MPI_Scatter(&all_localskips[0],1,MPI_INT,&nlocalskips[i],1,MPI_INT,0,world);
      MPI_Bcast(&nghostlyskips[i],1,MPI_INT,0,world);
      memory->destroy(local_rxncounts);
      memory->destroy(all_localskips);
    }
  }

  // this updates topology next step
  next_reneighbor = update->ntimestep;

  // call limit_bond in 'global_mega_glove mode.' oh, and local mode
  limit_bond(0); // add reacting atoms to nve/limit
  limit_bond(1);
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
  int flag;
  int index1 = atom->find_custom("limit_tags",flag);
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
  //  we accept this temporary and infrequent decrease in reaction occurences

  for (int i = 0; i < nxspecial[atom->map(glove[pion][1])][0]; i++) {
    if (atom->map(xspecial[atom->map(glove[pion][1])][i]) < 0) {
      error->one(FLERR,"Bond/react: Fix bond/react needs ghost atoms from further away"); // parallel issues.
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
              error->one(FLERR,"Bond/react: Fix bond/react needs ghost atoms from further away");
            }

            for (int j = 0; j < onemol_nxspecial[onemol_xspecial[pion][neigh]-1][0]; j++) {
              pioneer_count[onemol_xspecial[onemol_xspecial[pion][neigh]-1][j]-1]++;
            }

            glove_counter++;
            if (glove_counter == onemol->natoms) {
              status = ACCEPT;
              ring_check();
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
          error->one(FLERR,"Bond/react: Fix bond/react needs ghost atoms from further away");
        }

        for (int ii = 0; ii < onemol_nxspecial[onemol_xspecial[pion][neigh]-1][0]; ii++) {
          pioneer_count[onemol_xspecial[onemol_xspecial[pion][neigh]-1][ii]-1]++;
        }

        glove_counter++;
        if (glove_counter == onemol->natoms) {
          status = ACCEPT;
          ring_check();
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
        error->warning(FLERR,"Bond/react: Fix bond/react failed because MAXGUESS set too small. ask developer for info");
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

    int already_assigned = 0;
    for (int j = 0; j < onemol->natoms; j++) {
      if (glove[j][1] == xspecial[atom->map(glove[pion][1])][i]) {
        already_assigned = 1;
        break;
      }
    }

    if (already_assigned == 0 &&
        type[(int)atom->map(xspecial[atom->map(glove[pion][1])][i])] == onemol->type[(int)onemol_xspecial[pion][neigh]-1]) {
      if (num_choices > 5) { // here failed because too many identical first neighbors. but really no limit if situation arises
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

  //std::size_t size = sizeof(tag_choices) / sizeof(tag_choices[0]);
  std::sort(tag_choices, tag_choices + num_choices); //, std::greater<int>());
  glove[onemol_xspecial[pion][neigh]-1][0] = onemol_xspecial[pion][neigh];
  glove[onemol_xspecial[pion][neigh]-1][1] = tag_choices[guess_branch[avail_guesses-1]-1];
  guess_branch[avail_guesses-1]--;

  //another check for ghost atoms. perhaps remove the one in make_a_guess
  if (atom->map(glove[(int)onemol_xspecial[pion][neigh]-1][1]) < 0) {
    error->one(FLERR,"Bond/react: Fix bond/react needs ghost atoms from further away");
  }

  if (guess_branch[avail_guesses-1] == 0) avail_guesses--;

  for (int i = 0; i < onemol_nxspecial[onemol_xspecial[pion][neigh]-1][0]; i++) {
    pioneer_count[onemol_xspecial[onemol_xspecial[pion][neigh]-1][i]-1]++;
  }
  glove_counter++;
  if (glove_counter == onemol->natoms) {
    status = ACCEPT;
    ring_check();
    return;
  }
  status = CONTINUE;
}

/* ----------------------------------------------------------------------
  Check that newly assigned atoms have correct bonds
  Necessary for certain ringed structures
------------------------------------------------------------------------- */

void FixBondReact::ring_check()
{
  // ring_check can be made more efficient by re-introducing 'frozen' atoms
  // 'frozen' atoms have been assigned and also are no longer pioneers

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
      if (ring_fail == 1) {
        status = GUESSFAIL;
        return;
      }
    }
  }
}

/* ----------------------------------------------------------------------
evaluate constraints: return 0 if any aren't satisfied
------------------------------------------------------------------------- */

int FixBondReact::check_constraints()
{
  tagint atom1,atom2,atom3;
  double delx,dely,delz,rsq;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double rsq1,rsq2,r1,r2,c,t,prrhob;

  double **x = atom->x;

  for (int i = 0; i < nconstraints; i++) {
    if (constraints[i][0] == rxnID) {
      if (constraints[i][1] == DISTANCE) {
        atom1 = atom->map(glove[(int) constraints[i][2]-1][1]);
        atom2 = atom->map(glove[(int) constraints[i][3]-1][1]);
        delx = x[atom1][0] - x[atom2][0];
        dely = x[atom1][1] - x[atom2][1];
        delz = x[atom1][2] - x[atom2][2];
        domain->minimum_image(delx,dely,delz); // ghost location fix
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < constraints[i][4] || rsq > constraints[i][5]) return 0;
      } else if (constraints[i][1] == ANGLE) {
        atom1 = atom->map(glove[(int) constraints[i][2]-1][1]);
        atom2 = atom->map(glove[(int) constraints[i][3]-1][1]);
        atom3 = atom->map(glove[(int) constraints[i][4]-1][1]);

        // 1st bond
        delx1 = x[atom1][0] - x[atom2][0];
        dely1 = x[atom1][1] - x[atom2][1];
        delz1 = x[atom1][2] - x[atom2][2];
        rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
        r1 = sqrt(rsq1);

        // 2nd bond
        delx2 = x[atom3][0] - x[atom2][0];
        dely2 = x[atom3][1] - x[atom2][1];
        delz2 = x[atom3][2] - x[atom2][2];
        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
        r2 = sqrt(rsq2);

        // angle (cos and sin)
        c = delx1*delx2 + dely1*dely2 + delz1*delz2;
        c /= r1*r2;
        if (c > 1.0) c = 1.0;
        if (c < -1.0) c = -1.0;
        if (acos(c) < constraints[i][5] || acos(c) > constraints[i][6]) return 0;
      } else if (constraints[i][1] == ARRHENIUS) {
        t = get_temperature();
        prrhob = constraints[i][3]*pow(t,constraints[i][4])*
               exp(-constraints[i][5]/(force->boltz*t));
        if (prrhob < rrhandom[(int) constraints[i][2]]->uniform()) return 0;
      }
    }
  }
  return 1;
}

/* ----------------------------------------------------------------------
compute local temperature: average over all atoms in reaction template
------------------------------------------------------------------------- */

double FixBondReact::get_temperature()
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
      ilocal = atom->map(glove[i][1]);
      t += (v[ilocal][0]*v[ilocal][0] + v[ilocal][1]*v[ilocal][1] +
        v[ilocal][2]*v[ilocal][2]) * rmass[ilocal];
    }
  } else {
    for (i = 0; i < onemol->natoms; i++) {
      ilocal = atom->map(glove[i][1]);
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
  //   note: due to various usage/defintions of impropers, treated same as dihedrals
  // if angles exist: this means edge atoms not in their 1-2 neighbors list
  // if just bonds: this just means that edge atoms are not landlocked
  // Note: landlocked defined in terms of reacted template
  // if no edge atoms (small reacting molecule), all atoms are landlocked
  // we can delete all current topology of landlocked atoms and replace

  // always remove edge atoms from landlocked list
  for (int i = 0; i < twomol->natoms; i++) {
    if (edge[equivalences[i][1][myrxn]-1][myrxn] == 1) landlocked_atoms[i][myrxn] = 0;
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
    if (twomol->type[i] != onemol->type[equivalences[i][1][myrxn]-1] && landlocked_atoms[i][myrxn] == 0) {
      char str[128];
      snprintf(str,128,"Bond/react: Atom affected by reaction %s too close to template edge",rxn_name[myrxn]);
      error->all(FLERR,str);
    }
  }

  // additionally, if a bond changes type, but neither involved atom is landlocked, bad
  // would someone want to change an angle type but not bond or atom types? (etc.) ...hopefully not yet
  for (int i = 0; i < twomol->natoms; i++) {
    if (landlocked_atoms[i][myrxn] == 0) {
      for (int j = 0; j < twomol->num_bond[i]; j++) {
        int twomol_atomj = twomol->bond_atom[i][j];
        if (landlocked_atoms[twomol_atomj-1][myrxn] == 0) {
          int onemol_atomi = equivalences[i][1][myrxn];
          int onemol_batom;
          for (int m = 0; m < onemol->num_bond[onemol_atomi-1]; m++) {
            onemol_batom = onemol->bond_atom[onemol_atomi-1][m];
            if (onemol_batom == equivalences[twomol_atomj-1][1][myrxn]) {
              if (twomol->bond_type[i][j] != onemol->bond_type[onemol_atomi-1][m]) {
                char str[128];
                snprintf(str,128,"Bond/react: Atom affected by reaction %s too close to template edge",rxn_name[myrxn]);
                error->all(FLERR,str);
              }
            }
          }
          if (newton_bond) {
            int onemol_atomj = equivalences[twomol_atomj-1][1][myrxn];
            for (int m = 0; m < onemol->num_bond[onemol_atomj-1]; m++) {
              onemol_batom = onemol->bond_atom[onemol_atomj-1][m];
              if (onemol_batom == equivalences[i][1][myrxn]) {
                if (twomol->bond_type[i][j] != onemol->bond_type[onemol_atomj-1][m]) {
                  char str[128];
                  snprintf(str,128,"Bond/react: Atom affected by reaction %s too close to template edge",rxn_name[myrxn]);
                  error->all(FLERR,str);
                }
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
         error->all(FLERR,"Bond/react: A deleted atom cannot remain bonded to an atom that is not deleted");
        }
      }
    }
  }

  // also, if atoms change number of bonds, but aren't landlocked, that could be bad
  if (me == 0)
    for (int i = 0; i < twomol->natoms; i++) {
      if (twomol_nxspecial[i][0] != onemol_nxspecial[equivalences[i][1][myrxn]-1][0] && landlocked_atoms[i][myrxn] == 0) {
        char str[128];
        snprintf(str,128,"Bond/react: Atom affected by reaction %s too close to template edge",rxn_name[myrxn]);
        error->warning(FLERR,str);
        break;
      }
    }
}

/* ----------------------------------------------------------------------
let's dedup global_mega_glove
allows for same site undergoing different pathways, in parallel
------------------------------------------------------------------------- */

void FixBondReact::dedup_mega_gloves(int dedup_mode)
{
  // dedup_mode == 0 for local_dedup
  // dedup_mode == 1 for global_mega_glove
  for (int i = 0; i < nreacts; i++) {
    if (dedup_mode == 0) local_rxn_count[i] = 0;
    if (dedup_mode == 1) ghostly_rxn_count[i] = 0;
  }

  int dedup_size = 0;
  if (dedup_mode == 0) {
    dedup_size = local_num_mega;
  } else if (dedup_mode == 1) {
    dedup_size = global_megasize;
  }

  tagint **dedup_glove;
  memory->create(dedup_glove,max_natoms+1,dedup_size,"bond/react:dedup_glove");

  if (dedup_mode == 0) {
    for (int i = 0; i < dedup_size; i++) {
      for (int j = 0; j < max_natoms+1; j++) {
        dedup_glove[j][i] = local_mega_glove[j][i];
      }
    }
  } else if (dedup_mode == 1) {
    for (int i = 0; i < dedup_size; i++) {
      for (int j = 0; j < max_natoms+1; j++) {
        dedup_glove[j][i] = global_mega_glove[j][i];
      }
    }
  }

  // dedup_mask is size dedup_size and filters reactions that have been deleted
  // a value of 1 means this reaction instance has been deleted
  int *dedup_mask = new int[dedup_size];
  int *dup_list = new int[dedup_size];

  for (int i = 0; i < dedup_size; i++) {
    dedup_mask[i] = 0;
    dup_list[i] = 0;
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
      int num_dups = 0;
      int myrxnid1 = dedup_glove[0][i];
      onemol = atom->molecules[unreacted_mol[myrxnid1]];
      for (int j = 0; j < onemol->natoms; j++) {
        int check1 = dedup_glove[j+1][i];
        for (int ii = i + 1; ii < dedup_size; ii++) {
          int already_listed = 0;
          for (int jj = 0; jj < num_dups; jj++) {
            if (dup_list[jj] == ii) {
              already_listed = 1;
              break;
            }
          }
          if (dedup_mask[ii] == 0 && already_listed == 0) {
            int myrxnid2 = dedup_glove[0][ii];
            twomol = atom->molecules[unreacted_mol[myrxnid2]];
            for (int jj = 0; jj < twomol->natoms; jj++) {
              int check2 = dedup_glove[jj+1][ii];
              if (check2 == check1) {
                // add this rxn instance as well
                if (num_dups == 0) dup_list[num_dups++] = i;
                dup_list[num_dups++] = ii;
                break;
              }
            }
          }
        }
      }
      // here we choose random number and therefore reaction instance
      int myrand = 1;
      if (num_dups > 0) {
        myrand = floor(random[0]->uniform()*num_dups);
        for (int iii = 0; iii < num_dups; iii++) {
          if (iii != myrand) dedup_mask[dup_list[iii]] = 1;
        }
      }
    }
  }

  // we must update local_mega_glove and local_megasize
  // we can simply overwrite local_mega_glove column by column
  if (dedup_mode == 0) {
    int new_local_megasize = 0;
    for (int i = 0; i < local_num_mega; i++) {
      if (dedup_mask[i] == 0) {
        local_rxn_count[(int) dedup_glove[0][i]]++;
        for (int j = 0; j < max_natoms+1; j++) {
          local_mega_glove[j][new_local_megasize] = dedup_glove[j][i];
        }
        new_local_megasize++;
      }
    }

    local_num_mega = new_local_megasize;
  }

  // we must update global_mega_glove and global_megasize
  // we can simply overwrite global_mega_glove column by column
  if (dedup_mode == 1) {
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
  delete [] dup_list;
}

/* ----------------------------------------------------------------------
let's limit movement of newly bonded atoms
and exclude them from other thermostats via exclude_group
------------------------------------------------------------------------- */

void FixBondReact::limit_bond(int limit_bond_mode)
{
  //two types of passes: 1) while superimpose algorithm is iterating (only local atoms)
  //                     2) once more for global_mega_glove [after de-duplicating rxn instances]
  //in second case, only add local atoms to group
  //as with update_everything, we can pre-prepare these arrays, then run generic limit_bond code

  //create local, generic variables for onemol->natoms and glove
  //to be filled differently on respective passes

  int nlocal = atom->nlocal;
  int temp_limit_num = 0;
  tagint *temp_limit_glove;
  if (limit_bond_mode == 0) {
    int max_temp = local_num_mega * (max_natoms + 1);
    temp_limit_glove = new tagint[max_temp];
    for (int j = 0; j < local_num_mega; j++) {
      rxnID = local_mega_glove[0][j];
      onemol = atom->molecules[unreacted_mol[rxnID]];
      for (int i = 0; i < onemol->natoms; i++) {
        temp_limit_glove[temp_limit_num++] = local_mega_glove[i+1][j];
      }
    }

  } else if (limit_bond_mode == 1) {
    int max_temp = global_megasize * (max_natoms + 1);
    temp_limit_glove = new tagint[max_temp];
    for (int j = 0; j < global_megasize; j++) {
      rxnID = global_mega_glove[0][j];
      onemol = atom->molecules[unreacted_mol[rxnID]];
      for (int i = 0; i < onemol->natoms; i++) {
        if (atom->map(global_mega_glove[i+1][j]) >= 0 &&
            atom->map(global_mega_glove[i+1][j]) < nlocal)
          temp_limit_glove[temp_limit_num++] = global_mega_glove[i+1][j];
      }
    }
  }

  if (temp_limit_num == 0) {
    delete [] temp_limit_glove;
    return;
  }

  // we must keep our own list of limited atoms
  // this will be a new per-atom property!

  int flag;
  int index1 = atom->find_custom("limit_tags",flag);
  int *i_limit_tags = atom->ivector[index1];

  int *i_statted_tags;
  if (stabilization_flag == 1) {
    int index2 = atom->find_custom(statted_id,flag);
    i_statted_tags = atom->ivector[index2];
  }

  int index3 = atom->find_custom("react_tags",flag);
  int *i_react_tags = atom->ivector[index3];

  for (int i = 0; i < temp_limit_num; i++) {
    // update->ntimestep could be 0. so add 1 throughout
    i_limit_tags[atom->map(temp_limit_glove[i])] = update->ntimestep + 1;
    if (stabilization_flag == 1) i_statted_tags[atom->map(temp_limit_glove[i])] = 0;
    i_react_tags[atom->map(temp_limit_glove[i])] = rxnID;
  }

  delete [] temp_limit_glove;
}

/* ----------------------------------------------------------------------
let's unlimit movement of newly bonded atoms after n timesteps.
we give them back to the system thermostat
------------------------------------------------------------------------- */

void FixBondReact::unlimit_bond()
{
  //let's now unlimit in terms of i_limit_tags
  //we just run through all nlocal, looking for > limit_duration
  //then we return i_limit_tag to 0 (which removes from dynamic group)
  int flag;
  int index1 = atom->find_custom("limit_tags",flag);
  int *i_limit_tags = atom->ivector[index1];

  int *i_statted_tags;
  if (stabilization_flag == 1) {
    int index2 = atom->find_custom(statted_id,flag);
    i_statted_tags = atom->ivector[index2];
  }

  int index3 = atom->find_custom("react_tags",flag);
  int *i_react_tags = atom->ivector[index3];

  for (int i = 0; i < atom->nlocal; i++) {
    // unlimit atoms for next step! this resolves # of procs disparity, mostly
    // first '1': indexing offset, second '1': for next step
    if (i_limit_tags[i] != 0 && (update->ntimestep + 1 - i_limit_tags[i]) > limit_duration[i_react_tags[i]]) { // + 1
      i_limit_tags[i] = 0;
      if (stabilization_flag == 1) i_statted_tags[i] = 1;
      i_react_tags[i] = 0;
    }
  }

  //really should only communicate this per-atom property, not entire reneighboring
  next_reneighbor = update->ntimestep;
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

  int ghostly = 0;
  #if !defined(MPI_STUBS)
    if (comm->style == 0) {
      for (int i = 0; i < onemol->natoms; i++) {
        int ilocal = atom->map(glove[i][1]);
        if (ilocal >= atom->nlocal || localsendlist[ilocal] == 1) {
          ghostly = 1;
          break;
        }
      }
    } else {
      ghostly = 1;
    }
  #endif

  if (ghostly == 1) {
    ghostly_mega_glove[0][ghostly_num_mega] = rxnID;
    ghostly_rxn_count[rxnID]++; //for debuginng
    for (int i = 0; i < onemol->natoms; i++) {
      ghostly_mega_glove[i+1][ghostly_num_mega] = glove[i][1];
    }
    ghostly_num_mega++;
  } else {
    local_mega_glove[0][local_num_mega] = rxnID;
    local_rxn_count[rxnID]++; //for debuginng
    for (int i = 0; i < onemol->natoms; i++) {
      local_mega_glove[i+1][local_num_mega] = glove[i][1];
    }
    local_num_mega++;
  }
}

/* ----------------------------------------------------------------------
broadcast entries of mega_glove which contain nonlocal atoms for perusal by all processors
------------------------------------------------------------------------- */

void FixBondReact::ghost_glovecast()
{
#if !defined(MPI_STUBS)

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
  for (int i = 0; i < me; i++) {
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
  if (me == 0) {
    MPI_Gatherv(MPI_IN_PLACE, ghostly_num_mega, column, // Note: some values ignored for MPI_IN_PLACE
              &(global_mega_glove[0][0]), allncols, allstarts,
              column, 0, world);
  } else {
    MPI_Gatherv(&(global_mega_glove[0][start]), ghostly_num_mega, column,
              &(global_mega_glove[0][0]), allncols, allstarts,
              column, 0, world);
  }

  if (me == 0) dedup_mega_gloves(1); // global_mega_glove mode
  MPI_Bcast(&global_megasize,1,MPI_INT,0,world);
  MPI_Bcast(&(global_mega_glove[0][0]), global_megasize, column, 0, world);

  delete [] allstarts;
  delete [] allncols;

  MPI_Type_free(&column);
  MPI_Type_free(&columnunsized);
#endif
}

/* ----------------------------------------------------------------------
update charges, types, special lists and all topology
------------------------------------------------------------------------- */

void FixBondReact::update_everything()
{
  int *type = atom->type;

  int nlocal = atom->nlocal;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;

  // used when deleting atoms
  int ndel,ndelone;
  int *mark = new int[nlocal];
  for (int i = 0; i < nlocal; i++) mark[i] = 0;
  tagint *tag = atom->tag;
  AtomVec *avec = atom->avec;

  // update atom->nbonds, etc.
  // TODO: correctly tally with 'newton off'
  int delta_bonds = 0;
  int delta_angle = 0;
  int delta_dihed = 0;
  int delta_imprp = 0;

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
        for (int j = 0; j < max_natoms+1; j++)
          update_mega_glove[j][update_num_mega] = local_mega_glove[j][i];
        update_num_mega++;
      }
    } else if (pass == 1) {
      for (int i = 0; i < global_megasize; i++) {
        rxnID = global_mega_glove[0][i];
        // reactions already shuffled from dedup procedure, so can skip first N
        if (iskip[rxnID]++ < nghostlyskips[rxnID]) continue;
        for (int j = 0; j < max_natoms+1; j++)
          update_mega_glove[j][update_num_mega] = global_mega_glove[j][i];
        update_num_mega++;
      }
    }
    delete [] iskip;

    // mark to-delete atoms
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

    // update charges and types of landlocked atoms
    for (int i = 0; i < update_num_mega; i++) {
      rxnID = update_mega_glove[0][i];
      twomol = atom->molecules[reacted_mol[rxnID]];
      for (int j = 0; j < twomol->natoms; j++) {
        int jj = equivalences[j][1][rxnID]-1;
        if ((landlocked_atoms[j][rxnID] == 1 || custom_edges[jj][rxnID] == 1) &&
            atom->map(update_mega_glove[jj+1][i]) >= 0 &&
            atom->map(update_mega_glove[jj+1][i]) < nlocal) {
          type[atom->map(update_mega_glove[jj+1][i])] = twomol->type[j];
          if (twomol->qflag && atom->q_flag) {
            double *q = atom->q;
            q[atom->map(update_mega_glove[jj+1][i])] = twomol->q[j];
          }
        }
      }
    }

    //maybe add check that all 1-3 neighbors of a local atoms are at least ghosts -> unneeded --jg
    //okay, here goes:
    for (int i = 0; i < update_num_mega; i++) {
      rxnID = update_mega_glove[0][i];
      twomol = atom->molecules[reacted_mol[rxnID]];
      for (int j = 0; j < twomol->natoms; j++) {
        int jj = equivalences[j][1][rxnID]-1;
        if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
          if (landlocked_atoms[j][rxnID] == 1) {
            for (int k = 0; k < nspecial[atom->map(update_mega_glove[jj+1][i])][2]; k++) {
              if (atom->map(special[atom->map(update_mega_glove[jj+1][i])][k]) < 0) {
                error->all(FLERR,"Bond/react: Fix bond/react needs ghost atoms from further away - most likely too many processors");
              }
            }
          }
        }
      }
    }

    int insert_num;
    // very nice and easy to completely overwrite special bond info for landlocked atoms
    for (int i = 0; i < update_num_mega; i++) {
      rxnID = update_mega_glove[0][i];
      twomol = atom->molecules[reacted_mol[rxnID]];
      for (int j = 0; j < twomol->natoms; j++) {
        int jj = equivalences[j][1][rxnID]-1;
        if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
          if (landlocked_atoms[j][rxnID] == 1) {
            for (int k = 0; k < 3; k++) {
              nspecial[atom->map(update_mega_glove[jj+1][i])][k] = twomol->nspecial[j][k];
            }
            for (int p = 0; p < twomol->nspecial[j][2]; p++) {
              special[atom->map(update_mega_glove[jj+1][i])][p] = update_mega_glove[equivalences[twomol->special[j][p]-1][1][rxnID]][i];
            }
          }
          // now delete and replace landlocked atoms from non-landlocked atoms' special info
          if (landlocked_atoms[j][rxnID] == 0) {
            for (int k = nspecial[atom->map(update_mega_glove[jj+1][i])][2]-1; k > -1; k--) {
              for (int p = 0; p < twomol->natoms; p++) {
                int pp = equivalences[p][1][rxnID]-1;
                if (p!=j && special[atom->map(update_mega_glove[jj+1][i])][k] == update_mega_glove[pp+1][i]
                    && landlocked_atoms[p][rxnID] == 1 ) {
                  for (int n = k; n < nspecial[atom->map(update_mega_glove[jj+1][i])][2]-1; n++) {
                    special[atom->map(update_mega_glove[jj+1][i])][n] = special[atom->map(update_mega_glove[jj+1][i])][n+1];
                  }
                  if (k+1 > nspecial[atom->map(update_mega_glove[jj+1][i])][1]) {
                    nspecial[atom->map(update_mega_glove[jj+1][i])][2]--;
                  } else if (k+1 > nspecial[atom->map(update_mega_glove[jj+1][i])][0]) {
                    nspecial[atom->map(update_mega_glove[jj+1][i])][1]--;
                    nspecial[atom->map(update_mega_glove[jj+1][i])][2]--;
                  } else {
                    nspecial[atom->map(update_mega_glove[jj+1][i])][0]--;
                    nspecial[atom->map(update_mega_glove[jj+1][i])][1]--;
                    nspecial[atom->map(update_mega_glove[jj+1][i])][2]--;
                  }
                  break;
                }
              }
            }
            // now reassign from reacted template
            for (int k = 0; k < twomol->nspecial[j][2]; k++) {
              if (landlocked_atoms[twomol->special[j][k]-1][rxnID] == 1) {
                if (k > twomol->nspecial[j][1] - 1) {
                  insert_num = nspecial[atom->map(update_mega_glove[jj+1][i])][2]++;
                } else if (k > twomol->nspecial[j][0] - 1) {
                  insert_num = nspecial[atom->map(update_mega_glove[jj+1][i])][1]++;
                  nspecial[atom->map(update_mega_glove[jj+1][i])][2]++;
                } else {
                  insert_num = nspecial[atom->map(update_mega_glove[jj+1][i])][0]++;
                  nspecial[atom->map(update_mega_glove[jj+1][i])][1]++;
                  nspecial[atom->map(update_mega_glove[jj+1][i])][2]++;
                }
                if (nspecial[atom->map(update_mega_glove[jj+1][i])][2] > atom->maxspecial)
                  error->one(FLERR,"Bond/react special bond generation overflow");
                for (int n = nspecial[atom->map(update_mega_glove[jj+1][i])][2]-1; n > insert_num; n--) {
                  special[atom->map(update_mega_glove[jj+1][i])][n] = special[atom->map(update_mega_glove[jj+1][i])][n-1];
                }
                special[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->special[j][k]-1][1][rxnID]][i];
              }
            }
          }
        }
      }
    }

    // next let's update bond info
    // cool thing is, newton_bond issues are already taken care of in templates
    // same with class2 improper issues, which is why this fix started in the first place
    for (int i = 0; i < update_num_mega; i++) {
      rxnID = update_mega_glove[0][i];
      twomol = atom->molecules[reacted_mol[rxnID]];
      // let's first delete all bond info about landlocked atoms
      for (int j = 0; j < twomol->natoms; j++) {
        int jj = equivalences[j][1][rxnID]-1;
        if (atom->map(update_mega_glove[jj+1][i]) < nlocal && atom->map(update_mega_glove[jj+1][i]) >= 0) {
          if (landlocked_atoms[j][rxnID] == 1) {
            delta_bonds -= num_bond[atom->map(update_mega_glove[jj+1][i])];
            num_bond[atom->map(update_mega_glove[jj+1][i])] = 0;
          }
          if (landlocked_atoms[j][rxnID] == 0) {
            for (int p = num_bond[atom->map(update_mega_glove[jj+1][i])]-1; p > -1 ; p--) {
              for (int n = 0; n < twomol->natoms; n++) {
                int nn = equivalences[n][1][rxnID]-1;
                if (n!=j && bond_atom[atom->map(update_mega_glove[jj+1][i])][p] == update_mega_glove[nn+1][i] && landlocked_atoms[n][rxnID] == 1) {
                  for (int m = p; m < num_bond[atom->map(update_mega_glove[jj+1][i])]-1; m++) {
                    bond_type[atom->map(update_mega_glove[jj+1][i])][m] = bond_type[atom->map(update_mega_glove[jj+1][i])][m+1];
                    bond_atom[atom->map(update_mega_glove[jj+1][i])][m] = bond_atom[atom->map(update_mega_glove[jj+1][i])][m+1];
                  }
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
            }
          }
          if (landlocked_atoms[j][rxnID] == 0) {
            for (int p = 0; p < twomol->num_bond[j]; p++) {
              if (landlocked_atoms[twomol->bond_atom[j][p]-1][rxnID] == 1) {
                insert_num = num_bond[atom->map(update_mega_glove[jj+1][i])];
                bond_type[atom->map(update_mega_glove[jj+1][i])][insert_num] = twomol->bond_type[j][p];
                bond_atom[atom->map(update_mega_glove[jj+1][i])][insert_num] = update_mega_glove[equivalences[twomol->bond_atom[j][p]-1][1][rxnID]][i];
                num_bond[atom->map(update_mega_glove[jj+1][i])]++;
                if (num_bond[atom->map(update_mega_glove[jj+1][i])] > atom->bond_per_atom)
                  error->one(FLERR,"Bond/react topology/atom exceed system topology/atom");
                delta_bonds++;
              }
            }
          }
        }
      }
    }

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
                    error->one(FLERR,"Bond/react topology/atom exceed system topology/atom");
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
                    error->one(FLERR,"Bond/react topology/atom exceed system topology/atom");
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
                    error->one(FLERR,"Bond/react topology/atom exceed system topology/atom");
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
  delete [] mark;

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
}

/* ----------------------------------------------------------------------
read superimpose file
------------------------------------------------------------------------- */

void FixBondReact::read(int myrxn)
{
  char line[MAXLINE],keyword[MAXLINE];
  char *eof,*ptr;

  // skip 1st line of file
  eof = fgets(line,MAXLINE,fp);
  if (eof == NULL) error->one(FLERR,"Bond/react: Unexpected end of superimpose file");

  // read header lines
  // skip blank lines or lines that start with "#"
  // stop when read an unrecognized line

  while (1) {

    readline(line);

    // trim anything from '#' onward
    // if line is blank, continue

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    if (strstr(line,"edgeIDs")) sscanf(line,"%d",&nedge);
    else if (strstr(line,"equivalences")) {
      sscanf(line,"%d",&nequivalent);
      if (nequivalent != onemol->natoms)
        error->one(FLERR,"Bond/react: Number of equivalences in map file must "
                                  "equal number of atoms in reaction templates");
    }
    else if (strstr(line,"customIDs")) sscanf(line,"%d",&ncustom);
    else if (strstr(line,"deleteIDs")) sscanf(line,"%d",&ndelete);
    else if (strstr(line,"constraints")) {
      sscanf(line,"%d",&nconstr);
      memory->grow(constraints,nconstraints+nconstr,MAXCONARGS,"bond/react:constraints");
    } else break;
  }

  // grab keyword and skip next line

  parse_keyword(0,line,keyword);
  readline(line);

  // loop over sections of superimpose file

  int equivflag = 0, bondflag = 0, customedgesflag = 0;
  while (strlen(keyword)) {
    if (strcmp(keyword,"BondingIDs") == 0) {
      bondflag = 1;
      readline(line);
      sscanf(line,"%d",&ibonding[myrxn]);
      readline(line);
      sscanf(line,"%d",&jbonding[myrxn]);
    } else if (strcmp(keyword,"EdgeIDs") == 0) {
      EdgeIDs(line, myrxn);
    } else if (strcmp(keyword,"Equivalences") == 0) {
      equivflag = 1;
      Equivalences(line, myrxn);
    } else if (strcmp(keyword,"Custom Edges") == 0) {
      customedgesflag = 1;
      CustomEdges(line, myrxn);
    } else if (strcmp(keyword,"DeleteIDs") == 0) {
      DeleteAtoms(line, myrxn);
    } else if (strcmp(keyword,"Constraints") == 0) {
      Constraints(line, myrxn);
    } else error->one(FLERR,"Bond/react: Unknown section in map file");

    parse_keyword(1,line,keyword);

  }

  // error check
  if (bondflag == 0 || equivflag == 0)
    error->all(FLERR,"Bond/react: Map file missing BondingIDs or Equivalences section\n");

  if (update_edges_flag[myrxn] == 2 && customedgesflag == 0)
    error->all(FLERR,"Bond/react: Map file must have a Custom Edges section when using 'update_edges custom'\n");

  if (update_edges_flag[myrxn] != 2 && customedgesflag == 1)
    error->all(FLERR,"Bond/react: Specify 'update_edges custom' to include Custom Edges section in map file\n");
}

void FixBondReact::EdgeIDs(char *line, int myrxn)
{
  // puts a 1 at edge(edgeID)

  int tmp;
  for (int i = 0; i < nedge; i++) {
    readline(line);
    sscanf(line,"%d",&tmp);
    edge[tmp-1][myrxn] = 1;
  }
}

void FixBondReact::Equivalences(char *line, int myrxn)
{
  int tmp1;
  int tmp2;
  for (int i = 0; i < nequivalent; i++) {
    readline(line);
    sscanf(line,"%d %d",&tmp1,&tmp2);
    //equivalences is-> clmn 1: post-reacted, clmn 2: pre-reacted
    equivalences[tmp2-1][0][myrxn] = tmp2;
    equivalences[tmp2-1][1][myrxn] = tmp1;
    //reverse_equiv is-> clmn 1: pre-reacted, clmn 2: post-reacted
    reverse_equiv[tmp1-1][0][myrxn] = tmp1;
    reverse_equiv[tmp1-1][1][myrxn] = tmp2;
  }
}

void FixBondReact::CustomEdges(char *line, int myrxn)
{
  // 0 for 'none', 1 for 'charges'

  int tmp;
  int n = MAX(strlen("none"),strlen("charges")) + 1;
  char *edgemode = new char[n];
  for (int i = 0; i < ncustom; i++) {
    readline(line);
    sscanf(line,"%d %s",&tmp,edgemode);
    if (strcmp(edgemode,"none") == 0)
      custom_edges[tmp-1][myrxn] = 0;
    else if (strcmp(edgemode,"charges") == 0)
      custom_edges[tmp-1][myrxn] = 1;
    else
      error->one(FLERR,"Bond/react: Illegal value in 'Custom Edges' section of map file");
  }
  delete [] edgemode;
}

void FixBondReact::DeleteAtoms(char *line, int myrxn)
{
  int tmp;
  for (int i = 0; i < ndelete; i++) {
    readline(line);
    sscanf(line,"%d",&tmp);
    delete_atoms[tmp-1][myrxn] = 1;
  }
}

void FixBondReact::Constraints(char *line, int myrxn)
{
  double tmp[MAXCONARGS];
  int n = strlen("distance") + 1;
  char *constraint_type = new char[n];
  for (int i = 0; i < nconstr; i++) {
    readline(line);
    sscanf(line,"%s",constraint_type);
    constraints[nconstraints][0] = myrxn;
    if (strcmp(constraint_type,"distance") == 0) {
      constraints[nconstraints][1] = DISTANCE;
      sscanf(line,"%*s %lg %lg %lg %lg",&tmp[0],&tmp[1],&tmp[2],&tmp[3]);
      constraints[nconstraints][2] = tmp[0];
      constraints[nconstraints][3] = tmp[1];
      constraints[nconstraints][4] = tmp[2]*tmp[2]; // using square of distance
      constraints[nconstraints][5] = tmp[3]*tmp[3];
    } else if (strcmp(constraint_type,"angle") == 0) {
      constraints[nconstraints][1] = ANGLE;
      sscanf(line,"%*s %lg %lg %lg %lg %lg",&tmp[0],&tmp[1],&tmp[2],&tmp[3],&tmp[4]);
      constraints[nconstraints][2] = tmp[0];
      constraints[nconstraints][3] = tmp[1];
      constraints[nconstraints][4] = tmp[2];
      constraints[nconstraints][5] = tmp[3]/180.0 * MY_PI;
      constraints[nconstraints][6] = tmp[4]/180.0 * MY_PI;
    } else if (strcmp(constraint_type,"arrhenius") == 0) {
      constraints[nconstraints][1] = ARRHENIUS;
      constraints[nconstraints][2] = narrhenius++;
      sscanf(line,"%*s %lg %lg %lg %lg",&tmp[0],&tmp[1],&tmp[2],&tmp[3]);
      constraints[nconstraints][3] = tmp[0];
      constraints[nconstraints][4] = tmp[1];
      constraints[nconstraints][5] = tmp[2];
      constraints[nconstraints][6] = tmp[3];
    } else
      error->one(FLERR,"Bond/react: Illegal constraint type in 'Constraints' section of map file");
    nconstraints++;
  }
  delete [] constraint_type;
}

void FixBondReact::open(char *file)
{
  fp = fopen(file,"r");
  if (fp == NULL) {
    char str[128];
    snprintf(str,128,"Bond/react: Cannot open map file %s",file);
    error->one(FLERR,str);
  }
}

void FixBondReact::readline(char *line)
{
  int n;
  if (me == 0) {
    if (fgets(line,MAXLINE,fp) == NULL) n = 0;
    else n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  if (n == 0) error->all(FLERR,"Bond/react: Unexpected end of map file");
  MPI_Bcast(line,n,MPI_CHAR,0,world);
}

void FixBondReact::parse_keyword(int flag, char *line, char *keyword)
{
  if (flag) {

    // read upto non-blank line plus 1 following line
    // eof is set to 1 if any read hits end-of-file

    int eof = 0;
    if (me == 0) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
      while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
        if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
      }
      if (fgets(keyword,MAXLINE,fp) == NULL) eof = 1;
    }

    // if eof, set keyword empty and return

    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) {
      keyword[0] = '\0';
      return;
    }

    // bcast keyword line to all procs

    int n;
    if (me == 0) n = strlen(line) + 1;
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


void FixBondReact::skip_lines(int n, char *line)
{
  for (int i = 0; i < n; i++) readline(line);
}

int FixBondReact::parse(char *line, char **words, int max)
{
  char *ptr;
  char *r_token;

  r_token = line;
  int nwords = 0;
  words[nwords++] = strtok_r(r_token," \t\n\r\f",&r_token);

  while ((ptr = strtok_r(NULL," \t\n\r\f",&r_token))) {
    if (nwords < max) words[nwords] = ptr;
    nwords++;
  }

  return nwords;
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
      printf("hello you shouldn't be here\n");
      //buf[m++] = ubuf(bondcount[j]).d;
    }
    return m;
  }

  if (commflag == 2) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(partner[j]).d;
      buf[m++] = probability[j];
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
  int i,j,m,ns,last;

  m = 0;
  last = first + n;

  if (commflag == 1) {
    for (i = first; i < last; i++)
      printf("hello you shouldn't be here\n");
    // bondcount[i] = (int) ubuf(buf[m++]).i;
  } else if (commflag == 2) {
    for (i = first; i < last; i++) {
      partner[i] = (tagint) ubuf(buf[m++]).i;
      probability[i] = buf[m++];
    }
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

  if (commflag != 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      if (closeneigh[rxnID] != 0)
        if (buf[m+1] < distsq[j][1]) {
        partner[j] = (tagint) ubuf(buf[m++]).i;
          distsq[j][1] = buf[m++];
        } else m += 2;
      else
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
  set[0].nreacts = nreacts;
  for (int i = 0; i < nreacts; i++) {
    set[i].reaction_count_total = reaction_count_total[i];
    strncpy(set[i].rxn_name,rxn_name[i],MAXLINE);
    set[i].rxn_name[MAXLINE-1] = '\0';
  }

  if (me == 0) {
    int size = nreacts*sizeof(Set);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(set,sizeof(Set),nreacts,fp);
  }
}

/* ----------------------------------------------------------------------
   use selected state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixBondReact::restart(char *buf)
{
  Set *set_restart = (Set *) buf;
  for (int i = 0; i < set_restart[0].nreacts; i++) {
    for (int j = 0; j < nreacts; j++) {
      if (strcmp(set_restart[i].rxn_name,rxn_name[j]) == 0) {
        reaction_count_total[j] = set_restart[i].reaction_count_total;
      }
    }
  }
}

/* ----------------------------------------------------------------------
memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixBondReact::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes = 2*nmax * sizeof(tagint);
  bytes += nmax * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixBondReact::print_bb()
{

  //fix bond/create cargo code. eg nbonds needs to be added
  /*
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
  // printf(" " TAGINT_FORMAT,atom->special[i][j]);
  }
  // printf("\n");
}
*/
}
