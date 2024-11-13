// clang-format off
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

#include "fix_bond_dynamic.h"

#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "respa.h"
#include "update.h"
#include <iostream>

#include <cstring>
#include "math_const.h"
#include "random_mars.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBondDynamic::FixBondDynamic(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  distsq(nullptr), probabilities(nullptr), list(nullptr),
  random(nullptr), partners_possible_f(nullptr), partners_probs_f(nullptr),
  partners_possible(nullptr), partners_probs(nullptr), npos(nullptr), partners_success(nullptr)
{
  if (narg < 9) error->all(FLERR,"Illegal fix bond/dynamic command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix bond/dynamic command");

  dynamic_group_allow = 1;
  force_reneighbor = 1;
  next_reneighbor = -1;

  iatomtype = utils::inumeric(FLERR,arg[4],false,lmp);
  btype = utils::inumeric(FLERR,arg[5],false,lmp);
  ka = utils::numeric(FLERR,arg[6],false,lmp);
  kd = utils::numeric(FLERR,arg[7],false,lmp);
  double cutoff = utils::numeric(FLERR,arg[8],false,lmp);

  if (btype < 1 || btype > atom->nbondtypes)
    error->all(FLERR,"Invalid bond type in fix bond/dynamic command");
  if (cutoff < 0.0) error->all(FLERR,"Illegal fix bond/dynamic command");
  cutsq = cutoff*cutoff;

  // Default settings
  maxbond = atom->bond_per_atom;
  seed = 12345;
  jatomtype = iatomtype;

  // Flags for optional settings
  flag_prob = 0;
  flag_bell = 0;
  flag_catch = 0;
  flag_rouse = 0;
  flag_critical = 0;
  flag_mol = 0;

  // Parse remaining arguments
  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"prob") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix bond/dynamic command");
      prob_attach = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      prob_detach = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      flag_prob = 1;
      if (prob_attach < 0.0 || prob_attach > 1.0 || prob_detach < 0.0 || prob_detach > 1.0)
        error->all(FLERR,"Illegal fix bond/dynamic command");
      iarg += 3;
    } else if (strcmp(arg[iarg],"maxbond") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/dynamic command");
      maxbond = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (maxbond < 0) error->all(FLERR,"Illegal fix bond/dynamic command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"bell") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/dynamic command");
      f0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      flag_bell = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"catch") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix bond/dynamic command");
      fs0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      fc0 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      kc0_scale = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      flag_catch = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"rouse") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/dynamic command");
      double b0 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      b2 = b0*b0;
      flag_rouse = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"critical") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/dynamic command");
      double r_critical = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      r2_critical = r_critical*r_critical;
      flag_critical = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"jtype") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/dynamic command");
      jatomtype = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/dynamic command");
      flag_mol = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"seed") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/dynamic command");
      seed = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (seed <= 0) error->all(FLERR,"Illegal fix bond/dynamic command");
      iarg += 2; 
    } else error->all(FLERR,"Illegal fix bond/dynamic command");
  }

  // error checks
  if (flag_prob && flag_bell)
    error->all(FLERR,"Cannot use argument prob with argument bell");
  if (flag_prob && flag_catch)
    error->all(FLERR,"Cannot use argument prob with argument catch");
  if (flag_bell && flag_catch)
    error->all(FLERR,"Cannot use argument bell with argument catch");
  if (atom->molecular != Atom::MOLECULAR)
    error->all(FLERR,"Cannot use fix bond/dynamic with non-molecular systems");
  if (atom->bond_per_atom < maxbond)
    error->all(FLERR,"Maxbond too large in fix bond/dynamic - increase bonds/per/atom");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + me);

  // allocate values local to this fix
  nmax = 0;
  countflag = 0;
  distsq = probabilities = nullptr;
  partners_possible = partners_possible_f = nullptr;
  partners_probs = partners_probs_f = nullptr;
  npos = nullptr;
  partners_success = nullptr;

  comm_forward = 1+atom->maxspecial;
}

/* ---------------------------------------------------------------------- */

FixBondDynamic::~FixBondDynamic()
{
  delete random;

  // delete locally stored arrays

  memory->destroy(partners_possible);
  memory->destroy(partners_probs);
  memory->destroy(partners_possible_f);
  memory->destroy(partners_probs_f);
  memory->destroy(partners_success);
  memory->destroy(npos);

  if (new_fix_id && modify->nfix) modify->delete_fix(new_fix_id);
  delete [] new_fix_id;

}

/* ---------------------------------------------------------------------- */

int FixBondDynamic::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondDynamic::post_constructor()
{
  new_fix_id = utils::strdup(id + std::string("_FIX_PA"));
  modify->add_fix(fmt::format("{} {} property/atom i2_fbd_{} {} ghost yes",new_fix_id, group->names[igroup],id,std::to_string(maxbond)));

  int tmp1, tmp2;
  index = atom->find_custom(utils::strdup(std::string("fbd_")+id),tmp1,tmp2);
  tagint **fbd = atom->iarray[index];

  nmax = atom->nmax;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;

  for (int i = 0; i < nall; i++) {
    for (int m = 0; m < maxbond; m++) {
      if (mask[i] & groupbit) {
        fbd[i][m] = 0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBondDynamic::init()
{

  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // need a half neighbor list, built every Nevery steps
  neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);

}

/* ---------------------------------------------------------------------- */

void FixBondDynamic::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixBondDynamic::setup(int /*vflag*/)
{

  // compute initial bond neighbors if this is first run
  // can't do this earlier, in constructor or init, b/c need ghost info

  if (countflag) return;
  countflag = 1;

  // Initialize local atom array fbd
  // Custom list of dynamic bonds per atom
  //    0: open for bonding
  //   -1: to be broken this timestep
  //   -2: permanently broken
  //   >0: atom->tag of bonded atom

  tagint **fbd = atom->iarray[index];

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int nlocal = atom->nlocal;
  int **bond_type_raw = atom->bond_type;

  for (int i = 0; i < nlocal; i++) {
    if (num_bond[i] == 0) continue;
    for (int b = 0; b < num_bond[i]; b++) {
       //if (bond_type[i][b] == btype) {
      fbd[i][b] = bond_atom[i][b];
       //}
    }
  }

  // forward communication of fbd so ghost atoms have it
  commflag = 1;
  comm->forward_comm(this,maxbond);

  // Create initial memory allocations
  memory->create(partners_possible,nmax,maxbond,"bond/dynamic:possible_partners");
  memory->create(partners_probs,nmax,maxbond,"bond/dynamic:partners_probs");
  memory->create(partners_possible_f,nmax,maxbond,"bond/dynamic:possible_partners_f");
  memory->create(partners_probs_f,nmax,maxbond,"bond/dynamic:partners_probs_f");
  memory->create(partners_success,nmax,maxbond,"bond/dynamic:partners_success");
  memory->create(npos,nmax,"bond/dynamic:npos");
}

/* ---------------------------------------------------------------------- */

void FixBondDynamic::post_integrate()
{

  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;

  tagint **fbd = atom->iarray[index];

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  int  bondid = atom->bond_per_atom;
  tagint **bond_atom = atom->bond_atom;

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  int bondtype = 0;

  if (update->ntimestep % nevery) return;

  // acquire updated ghost atom positions
  // necessary b/c are calling this after integrate, but before Verlet comm
  comm->forward_comm();

  // basic atom information
  double **x = atom->x;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int *type = atom->type;
  Bond *bond = force->bond;
  double DT_EQ = (update->dt)*nevery;
  // JTC: Probably not worth worrying about, but this definition of DT_EQ won't be
  // compatible with a variable timestep like that used in fix dt/reset.
  // Not sure there's a great solution (maybe incrementing?) or a good error check

  /* BEGIN BREAKING PROCESS */
  // loop over local atoms
  // check for possible breaks

  //int size = sizeof(bond_type_raw);
  //printf("%i\n\n",size);

  //printf("%i\n\n",bondlist[0][2]);
  //printf("%i\n\n",bondlist[1][2]);
  //printf("%i\n\n",bondlist[2][2]);
  

  for (int i = 0; i < nlocal; i++) {

    // Skip atoms not in the desired group or of the wrong type
    if (!(mask[i] & groupbit)) continue;
    //if ((type[i] != iatomtype) && (type[i] != jatomtype)) continue;  //### TEMP (too restrictive)

    // Loop through each entry of fbd
    for (int b = 0; b < maxbond; b++) {
      
      // Tag of current bond pair
      tagint tagj = fbd[i][b];
      
      // tagj < 1 means bond is already detached or there is no bond
      if (tagj < 1) continue;

      // Skip bonds that don't belong to the right type
      //printf("nbonds %i\n",nbondlist);
      for (int n = 0; n < nbondlist; n++) {
        int iatom = bondlist[n][0];
        int jatom = bondlist[n][1];
        //printf("n:%i bondtype: %i \n",n,bondtype);

        //printf("iatom:%i \n",tag[iatom]);
        //printf("jatom:%i \n",tag[jatom]);
        //printf("tagi:%i \n",tag[i]);
        //printf("tagj:%i \n\n",tagj);

        if((tag[iatom]==tag[i] and tag[jatom]==tagj) || (tag[iatom]==tagj and tag[jatom]==tag[i])) {
          bondtype = bondlist[n][2];
          break;
        }

      }

      //printf("%i\n",fbd[i][b]);
      //printf("%i\n",tag[i]);
      //printf("bondtype %f\n",bondtype);
      //printf("conditional %i\n",(bondtype != btype));

      if (bondtype != btype) continue;
 
      // Local id of current bond pair
      int j = atom->map(tagj);

      //if (j < 0) continue;
      //  error->one(FLERR,"Fix bond/dynamic needs ghost atoms "    //### TEMP  
      //              "from further away 1");

      // Skip atoms not in the desired group or of the wrong type
      if (!(mask[j] & groupbit)) continue;
      //if ((type[j] != iatomtype) && (type[j] != jatomtype)) continue; //### TEMP (too restrictive)

      // Only consider each bond once - when my atom has the lower atom tag
      if (tag[i] > tagj) continue;

      // Random number for stochastic detachment
      double probability = random->uniform();

      // Detachment probability for constant kd
      double p_detach = 1 - exp(-kd*DT_EQ);

      // Flags that modify kd
      if (flag_bell) {

        // Find distance between two atoms
        double delx = x[i][0] - x[j][0];
        double dely = x[i][1] - x[j][1];
        double delz = x[i][2] - x[j][2];
        domain->minimum_image(delx, dely, delz);
        double rsq = delx*delx + dely*dely + delz*delz;

        // Find force in bond
        double fbond; // fbond is returned as f/r
        double engpot = bond->single(btype,rsq,i,j,fbond);
        double r = sqrt(rsq);
        double bondforce = fabs(fbond)*r;

        // Modify kd using Bell's law
        double kd_bell = kd*exp(fabs(bondforce)/f0);
        p_detach = 1 - exp(-kd_bell*DT_EQ);
      }
      if (flag_catch) {

         // Find distance between two atoms
        double delx = x[i][0] - x[j][0];
        double dely = x[i][1] - x[j][1];
        double delz = x[i][2] - x[j][2];
        domain->minimum_image(delx, dely, delz);
        double rsq = delx*delx + dely*dely + delz*delz;

        // Find force in bond
        double fbond; // fbond is returned as f/r
        double engpot = bond->single(btype,rsq,i,j,fbond);
        double r = sqrt(rsq);
        double bondforce = fabs(fbond)*r; 

        // Modify kd using two-path catch model
        // kd = slip + catch
        double kd_catch = kd*exp(fabs(bondforce)/fs0) + kd*kc0_scale*exp(-fabs(bondforce)/fc0);
        //printf("kd_catch %4.4f\n",kd_catch);
        //printf("fbond %4.4f\n",fbond);
        //printf("fs0 %4.4f\n",fs0);
        //printf("fc0 %4.4f\n",fc0);
        //printf("kco_scale %4.4f\n",kc0_scale);
        //printf("bond %4.4f\n",bondforce);
        
        p_detach = 1 - exp(-kd_catch*DT_EQ);
      }
      if (flag_critical) {

        // Find distance between two atoms
        double delx = x[i][0] - x[j][0];
        double dely = x[i][1] - x[j][1];
        double delz = x[i][2] - x[j][2];
        domain->minimum_image(delx, dely, delz);
        double rsq = delx*delx + dely*dely + delz*delz;

        // Compare to critical length for forced detachment
        if (rsq >= r2_critical) p_detach = 1.0;
      }
      if (flag_prob) {

        // Set detachment probability directly
        p_detach = prob_detach;
      }

      // Apply probability constraint
      if (probability > p_detach) continue;

      //double delx = x[i][0] - x[j][0];
      //double dely = x[i][1] - x[j][1];
      //double delz = x[i][2] - x[j][2];
      //domain->minimum_image(delx, dely, delz);
      //double rsq = delx*delx + dely*dely + delz*delz;
      //double fbond;
      //double engpot = bond->single(btype,rsq,i,j,fbond);
      //if (btype==1){
      //  printf("btype %i\n",btype);
      //}

      if (kd == 0) continue; 

      // if breaking was successful, update fbd to -tag
      fbd[i][b] *= -1;

      // find the entry of atom j and update its fbd as well
      // if j is a ghost atom, it will do this on its own processor
      // at the next step
      if (j < nlocal) {
        for (int bb = 0; bb < maxbond; bb++) {
          if (fbd[j][bb] == tag[i]) {
            fbd[j][bb] *= -1;
            break;
          }
        }
      }
    }
  }

  // forward communication of fbd so ghost atoms store their breaks
  commflag = 1;
  comm->forward_comm(this,maxbond);

  // Loop over ghost atoms, find corresponding entries of fdb and update
  // needed when a ghost atom bonded to an owned atom decides to break its bond
  for (int j = nlocal; j < nall; j++) {
    for (int b = 0; b < maxbond; b++) {
      tagint tagi = fbd[j][b];
      
      if (tagi > -1 || tagi == -INT_MAX) continue;

      int i = atom->map(-tagi);
      if (i < 0) continue;

      // find the entry of atom i and update its fbd as well
      for (int bb = 0; bb < maxbond; bb++) {
        if (fbd[i][bb] == tag[j]) {
          fbd[i][bb] *= -1;
          break;
        }
      }
    }
  }

  /* BEGIN CREATION PROCESS */

  // Possibly resize the possible_partners array
  if (atom->nmax > nmax) {
    memory->destroy(partners_possible);
    memory->destroy(partners_probs);
    memory->destroy(partners_possible_f);
    memory->destroy(partners_probs_f);
    memory->destroy(partners_success);
    memory->destroy(npos);
    nmax = atom->nmax;
    memory->create(partners_possible,nmax,maxbond,"bond/dynamic:possible_partners");
    memory->create(partners_probs,nmax,maxbond,"bond/dynamic:partners_probs");
    memory->create(partners_possible_f,nmax,maxbond,"bond/dynamic:possible_partners_f");
    memory->create(partners_probs_f,nmax,maxbond,"bond/dynamic:partners_probs_f");
    memory->create(partners_success,nmax,maxbond,"bond/dynamic:partners_success");
    memory->create(npos,nmax,"bond/dynamic:npos");
  }

  // Initialize arrays to zero
  for (int i = 0; i < nall; i++) {
    for (int j = 0; j < maxbond; j++) {
      partners_possible[i][j] = 0;
      partners_possible_f[i][j] = 0;
      partners_success[i][j] = 0;
      partners_probs[i][j] = 1.0;
      partners_probs_f[i][j] = 1.0;
    }
    npos[i] = 0;
  }

  // Determine how many open slots each atom has
  for (int i = 0; i < nlocal; i++) {
    for (int b = 0; b < maxbond; b++) {
      if (fbd[i][b] == 0) npos[i]++;
    }
  }

  // Forward communication of open slots
  commflag = 2;
  comm->forward_comm(this,1);

  // build temporary neighbor list to determine closest images
  // neighbor->build_one(list,1); OLD
  neighbor->build_one(list);
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  tagint *molecule = atom->molecule;

  // find potential bonding partners

  for (int i = 0; i < nlocal; i++) {

    // Skip if ka = 0;
    if (ka == 0) continue;

    // Skip irrelevant atoms
    if (!(mask[i] & groupbit)) continue;
    if ((type[i] != iatomtype) && (type[i] != jatomtype)) continue;
    if (npos[i] == 0) continue;

    // Neighbor list of atom i
    int *jlist = firstneigh[i];
    int jnum  = numneigh[i];

    // Looping through neighbor list of atom i
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      // Skip irrelevant atoms
      if (!(mask[j] & groupbit)) continue;
      if ((type[j] != iatomtype) && (type[j] != jatomtype)) continue;
      if (npos[j] == 0) continue;
      if (tag[i] == tag[j]) continue;

      // Skip bonds that don't belong to the right type
      //if (bond_type[i][j] != btype) continue;

      // flag_mol = 1 means only atoms on different molecules can bond
      if (flag_mol == 1) {
        if (molecule[i] == molecule[j]) continue;
      }
      // flag_mol = 2 means only atoms on the same molecule can bond
      if (flag_mol == 2) {
        if (molecule[i] != molecule[j]) continue;
      }

      // do not allow a duplicate bond to be created
      // check fbd matrix of atom i
      // abs() in case this bond was just broken
      int done = 0;
      for (int b = 0; b < maxbond; b++) {
        if (abs(fbd[i][b]) == tag[j]) {
          done = 1;
          break;
        }
      }
      if (done) continue;

      // check if this ghost atom was already seen
      // possible due to duplicate ghost atoms
      for (int n = 0; n < maxbond; n++) {
        if (partners_possible[i][n] == tag[j]) {
          done = 1;
          break;
        }
      }
      if (done) continue;

      // Find the distance between these two atoms
      double delx = x[i][0] - x[j][0];
      double dely = x[i][1] - x[j][1];
      double delz = x[i][2] - x[j][2];
      domain->minimum_image(delx, dely, delz);
      double rsq = delx*delx + dely*dely + delz*delz;

      // Skip if out of range
      if (rsq > cutsq) continue;

      // Determine probability of attachment
      double probability = random->uniform();

      // Attachment probability at constant ka
      double p_attach = 1.0 - exp(-ka*DT_EQ);

      // Flags that modify ka
      if (flag_rouse) {
        double ka_rouse = ka*pow(b2/rsq,2);
        p_attach = 1 - exp(-ka_rouse*DT_EQ);
      }
      if (flag_prob) {

        // Set attachment probability directly
        p_attach = prob_attach;
      }

      // No reason to consider this if it will be unsuccessful
      if (probability > p_attach) continue;

      // Next, check where to insert, if at all
      int bb = 0;
      for (int n = 0; n < maxbond; n++) {
        if (partners_probs[i][n] < probability) {
          bb++;
        } else {
          break;
        }
      }

      // continue if this is not possible
      if (bb > maxbond-1) continue; 

      // If bb is the last entry, no shifting required
      if (bb == maxbond-1) {
        partners_possible[i][bb] = tag[j];
        partners_probs[i][bb] = probability;
      } else {
        // Shift all elements to the right
        for (int n = maxbond-2; n >= bb; n--) {
          partners_possible[i][n+1] = partners_possible[i][n];
          partners_probs[i][n+1] = partners_probs[i][n];
          if (bb == maxbond-2) break;
        }

        partners_possible[i][bb] = tag[j];
        partners_probs[i][bb] = probability;
      }
    }
  }

  // forward communication of partners arrays
  commflag = 3;
  comm->forward_comm(this,2*maxbond);

  // compile list of atoms j that see my owned atoms i
  // could be ghost or local
  for (int j = 0; j < nall; j++) {

    // Skip irrelevant atoms
    if (!(mask[j] & groupbit)) continue;
    if ((type[j] != iatomtype) && (type[j] != jatomtype)) continue;
    if (npos[j] == 0) continue;

    // Loop through each entry of fbd
    for (int b = 0; b < maxbond; b++) {
      tagint tagi = partners_possible[j][b];
      
      if (tagi < 1) continue;

      int i = atom->map(tagi);
      if (i < 0 || i > nlocal) continue;

      // First case: atom i already has atom j as a final partner
      // check probabilities - update to more likely one if needed
      int done = 0;
      for (int bb = 0; bb < maxbond; bb++) {
        if (partners_possible_f[i][bb] == tag[j]) {
          done = 1;
          if (partners_probs[j][b] < partners_probs_f[i][bb]) partners_probs_f[i][bb] = partners_probs[j][b];
          break;
        }
      }
      if (done) continue;

      // Second case: atom i has open slots for final partners
      // just insert the new partner and probability in the first open slot
      for (int bb = 0; bb < maxbond; bb++) {
        if (partners_possible_f[i][bb] == 0) {
          partners_possible_f[i][bb] = tag[j];
          partners_probs_f[i][bb] = partners_probs[j][b];
          done = 1;
          break;
        }
      }
      if (done) continue;

      // Last case: atom i already has a full final partners matrix
      // find the least likely partner - check if this one is more likely and replace
      int imax = 0;
      int pmax = partners_probs_f[i][imax];
      for (int bb = 0; bb < maxbond; bb++) {
        if (partners_probs_f[i][bb] > pmax) {
          pmax = partners_probs_f[i][bb];
          imax = bb;
        }
      }
      if (partners_probs[j][b] < pmax) {
        partners_probs_f[i][imax] = partners_probs[j][b];
        partners_possible_f[i][imax] = tag[j];
      }
    }
  }

  // Add lists together - no need to communicate partners_probs_f
  for (int i = 0; i < nlocal; i++) {
    for (int b = 0; b < maxbond; b++) {

      double probability = partners_probs_f[i][b];
      tagint tagj = partners_possible_f[i][b];

      // Check where to insert
      int bb = 0;
      for (int n = 0; n < maxbond; n++) {
        if (partners_probs[i][n] < probability) {
          bb++;
        } else {
          break;
        }
      }

      // continue if this is not possible
      if (bb > maxbond-1) continue;

      // If bb is the last entry, no shifting required
      if (bb == maxbond-1) {
        partners_possible[i][bb] = tagj;
        partners_probs[i][bb] = probability;
      } else {
        // Shift all elements to the right
        for (int n = maxbond-2; n >= bb; n--) {
          partners_possible[i][n+1] = partners_possible[i][n];
          partners_probs[i][n+1] = partners_probs[i][n];
          if (bb == maxbond-2) break;
        }

        partners_possible[i][bb] = tagj;
        partners_probs[i][bb] = probability;
      }
    }
  }

  // forward communication of partners arrays
  commflag = 3;
  comm->forward_comm(this,2*maxbond);

  // At this point, the formation of bonds is completely determined
  // Just need to limit the formation to npos per atom and update arrays

  for (int i = 0; i < nlocal; i++) {

    // Skip irrelevant atoms
    if (!(mask[i] & groupbit)) continue;
    if ((type[i] != iatomtype) && (type[i] != jatomtype)) continue;
    if (npos[i] == 0) continue;

    // Loop through possibilites
    for (int n = 0; n < npos[i]; n++) {
      tagint tagj = partners_possible[i][n];

      if (tagj < 1) continue;

      // Check local id of potential partner
      int j = atom->map(tagj);
      if (j < 0)
        error->one(FLERR,"Fix bond/dynamic needs ghost atoms "
                    "from further away 2");

      // find where this bond is in atom j's list
      // if its location is past npos[j], then this was unsuccessful
      // this will be consistant across processors as we are looping through npos[i]
      int bb = 0;
      for (int b = 0; b < maxbond; b++) {
        if (partners_possible[j][b] == 0) {
          continue;
        } else if (partners_possible[j][b] == tag[i]) {
          break;
        } else {
          bb++;
        }
      }
      if (bb > npos[j]-1) continue;

      // Success! A bond was created. Mark as successful
      // atom j will also do this, whatever proc it's on
      partners_success[i][n] = 1;
    }
  }

  // forward communication of success array
  commflag = 4;
  comm->forward_comm(this,maxbond);

  for (int i = 0; i < nlocal; i++) {
    for (int b = 0; b < maxbond; b++) {

      // First, process broken bonds
      // check for negative fbd entry
      if (fbd[i][b] < 0) {
        tagint tagj = -fbd[i][b];
        if (tagj == INT_MAX) continue;
        int j = atom->map(tagj);

        if (j < 0)
          error->one(FLERR,"Fix bond/dynamic needs ghost atoms "
                      "from further away");

        // Check if this broke irreversibly
        int flag_remove = 0;
        if (flag_critical) {
          double delx = x[i][0] - x[j][0];
          double dely = x[i][1] - x[j][1];
          double delz = x[i][2] - x[j][2];
          domain->minimum_image(delx, dely, delz);
          double rsq = delx*delx + dely*dely + delz*delz;
          if (rsq >= r2_critical) flag_remove = 1;
        }

        // Update atom properties and fbd
        process_broken(i,j);
        fbd[i][b] = 0;

        // Only do this once if both local
        if (j < nlocal && i < j) {
          for (int bb = 0; bb < maxbond; bb++) {
            if (fbd[j][bb] == -tag[i]) {
              fbd[j][bb] = 0;

              if (flag_critical && flag_remove) {
                fbd[i][b] = -INT_MAX;
                fbd[j][bb] = -INT_MAX;
              }
              break;
            }
          }
        } else if (j >= nlocal) {
          if (flag_critical && flag_remove) {
            fbd[i][b] = -INT_MAX;
          }
        }
      }

      // Next, process created bonds
      // check for positive partners_success entry
      if (!partners_success[i][b]) continue;

      tagint tagj = partners_possible[i][b];
      int j = atom->map(tagj);

      // Only do this once if both local
      if (j < nlocal && i > j) continue;

      if (j < 0)
        error->one(FLERR,"Fix bond/dynamic needs ghost atoms "
                    "from further away 3");

      // do not allow a duplicate bond to be created
      // check fbd entry of atom i
      int done = 0;
      for (int k = 0; k < maxbond; k++) {
        if (fbd[i][k] == tag[j]) {
          done = 1;
          break;
        }
      }
      if (done) continue;

      // Update atom properties and fbd
      process_created(i,j);
      for (int bb = 0; bb < maxbond; bb++) {
        if (fbd[i][bb] == 0) {
          fbd[i][bb] = tagj;
          break;
        }
      }

      // Update atom j if owned
      if (j < nlocal) {
        for (int bb = 0; bb < maxbond; bb++) {
          if (fbd[j][bb] == 0) {
            fbd[j][bb] = tag[i];
            break; 
          }
        }
      }

    }
  }

  // forward communication of fbd
  commflag = 1;
  comm->forward_comm(this,maxbond);

  // forward communication of special lists
  commflag = 5;
  comm->forward_comm(this);

  update_topology();

  // trigger reneighboring
  next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixBondDynamic::process_broken(int i, int j)
{
  
  // First add the pair to new_broken_pairs
  auto tag_pair = std::make_pair(atom->tag[i], atom->tag[j]);
  new_broken_pairs.push_back(tag_pair);

  // Manually search and remove from atom arrays
  // need to remove in case special bonds arrays rebuilt
  int nlocal = atom->nlocal;

  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;
  

  if (i < nlocal) {
    int n = num_bond[i];

    int done = 0;
    for (int m = 0; m < n; m++) {
      if (bond_atom[i][m] == tag[j]) {
        for (int k = m; k < n - 1; k++) {
          bond_type[i][k] = bond_type[i][k + 1];
          bond_atom[i][k] = bond_atom[i][k + 1];
        }
        num_bond[i]--;
        break;
      }
      if (done) break;
    }
  }

  if (j < nlocal) {
    int n = num_bond[j];

    int done = 0;
    for (int m = 0; m < n; m++) {
      if (bond_atom[j][m] == tag[i]) {
        for (int k = m; k < n - 1; k++) {
          bond_type[j][k] = bond_type[j][k + 1];
          bond_atom[j][k] = bond_atom[j][k + 1];
        }
        num_bond[j]--;
        break;
      }
      if (done) break;
    }
  }

  // Update special neighbor list
  tagint *slist;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  // remove i from special bond list for atom j and vice versa
  // ignore n2, n3 since 1-3, 1-4 special factors required to be 1.0
  if (i < nlocal) {
    slist = special[i];
    int n1 = nspecial[i][0];
    int m;
    for (m = 0; m < n1; m++)
      if (slist[m] == tag[j]) break;
    for (; m < n1 - 1; m++) slist[m] = slist[m + 1];
    nspecial[i][0]--;
    nspecial[i][1] = nspecial[i][2] = nspecial[i][0];
  }

  if (j < nlocal) {
    slist = special[j];
    int n1 = nspecial[j][0];
    int m;
    for (int m = 0; m < n1; m++)
      if (slist[m] == tag[i]) break;
    for (; m < n1 - 1; m++) slist[m] = slist[m + 1];
    nspecial[j][0]--;
    nspecial[j][1] = nspecial[j][2] = nspecial[j][0];
  }

}

/* --------------------------------------------------------------------- */

void FixBondDynamic::process_created(int i, int j)
{

  // First add the pair to new_created_pairs
  auto tag_pair = std::make_pair(atom->tag[i], atom->tag[j]);
  new_created_pairs.push_back(tag_pair);

  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;

  int nlocal = atom->nlocal;

  // Add bonds to atom class for i and j
  if (i < nlocal) {
    if (num_bond[i] == atom->bond_per_atom)
      error->one(FLERR,"New bond exceeded bonds per atom in fix bond/dynamic");
    bond_type[i][num_bond[i]] = btype;
    bond_atom[i][num_bond[i]] = atom->tag[j];
    num_bond[i]++;
  }

  if (j < nlocal) {
    if (num_bond[j] == atom->bond_per_atom)
      error->one(FLERR,"New bond exceeded bonds per atom in fix bond/dynamic");
    bond_type[j][num_bond[j]] = btype;
    bond_atom[j][num_bond[j]] = atom->tag[i];
    num_bond[j]++;
  }

  // add i to special bond list for atom j and vice versa
  // ignore n2, n3 since 1-3, 1-4 special factors required to be 1.0

  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  if (i < nlocal) {
    int n1 = nspecial[i][0];
    if (n1 >= atom->maxspecial)
      error->one(FLERR, "Special list size exceeded in fix bond/dynamic");
    special[i][n1] = atom->tag[j];
    nspecial[i][0] += 1;
    nspecial[i][1] = nspecial[i][2] = nspecial[i][0];
  }

  if (j < nlocal) {
    int n1 = nspecial[j][0];
    if (n1 >= atom->maxspecial)
      error->one(FLERR, "Special list size exceeded in fix bond/dynamic");
    special[j][n1] = atom->tag[i];
    nspecial[j][0] += 1;
    nspecial[j][1] = nspecial[j][2] = nspecial[j][0];
  }

}

/* ----------------------------------------------------------------------
  Update special lists for recently broken/created bonds
  Assumes appropriate atom/bond arrays were updated, e.g. had called
      neighbor->add_temporary_bond(i1, i2, btype);
------------------------------------------------------------------------- */

void FixBondDynamic::update_topology()
{

  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;

  // In theory could communicate a list of broken bonds to neighboring processors here
  // to remove restriction that users use Newton bond off

  for (int ilist = 0; ilist < neighbor->nlist; ilist++) {
    NeighList *list = neighbor->lists[ilist];

    // Skip copied lists, will update original
    if (list->copy) continue;

    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;

    for (auto const &it : new_broken_pairs) {
      tagint tag1 = it.first;
      tagint tag2 = it.second;
      int i1 = atom->map(tag1);
      int i2 = atom->map(tag2);

      if (i1 < 0 || i2 < 0) {
        error->one(FLERR,"Fix bond/dynamic needs ghost atoms "
                    "from further away 4");
      }

      // Loop through atoms of owned atoms i j
      if (i1 < nlocal) {
        int *jlist = firstneigh[i1];
        int jnum = numneigh[i1];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= SPECIALMASK;    // Clear special bond bits
          if (tag[j] == tag2) jlist[jj] = j;
        }
      }

      if (i2 < nlocal) {
        int *jlist = firstneigh[i2];
        int jnum = numneigh[i2];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= SPECIALMASK;    // Clear special bond bits
          if (tag[j] == tag1) jlist[jj] = j;
        }
      }
    }
  }

  for (int ilist = 0; ilist < neighbor->nlist; ilist++) {
    NeighList *list = neighbor->lists[ilist];

    // Skip copied lists, will update original
    if (list->copy) continue;

    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;

    for (auto const &it : new_created_pairs) {
      tagint tag1 = it.first;
      tagint tag2 = it.second;
      int i1 = atom->map(tag1);
      int i2 = atom->map(tag2);

      if (i1 < 0 || i2 < 0) {
        error->one(FLERR,"Fix bond/dynamic needs ghost atoms "
                    "from further away 5");
      }

      // Loop through atoms of owned atoms i j
      if (i1 < nlocal) {
        int *jlist = firstneigh[i1];
        int jnum = numneigh[i1];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          if (((j >> SBBITS) & 3) != 0) continue;               // Skip bonded pairs
          if (tag[j] == tag2) jlist[jj] = j ^ (1 << SBBITS);    // Add 1-2 special bond bits
        }
      }

      if (i2 < nlocal) {
        int *jlist = firstneigh[i2];
        int jnum = numneigh[i2];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          if (((j >> SBBITS) & 3) != 0) continue;               // Skip bonded pairs
          if (tag[j] == tag1) jlist[jj] = j ^ (1 << SBBITS);    // Add 1-2 special bond bits
        }
      }
    }
  }

  new_broken_pairs.clear();
  new_created_pairs.clear();

}

/* ---------------------------------------------------------------------- */

void FixBondDynamic::post_integrate_respa(int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_integrate();
}

/* ---------------------------------------------------------------------- */

int FixBondDynamic::pack_forward_comm(int n, int *list, double *buf,
                                    int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;

  if (commflag == 1) {
      tagint **fbd = atom->iarray[index];

      for (int i = 0; i < n; i++) {
        int j = list[i];
        for (int k = 0; k < maxbond; k++) {
          buf[m++] = ubuf(fbd[j][k]).d;
        }
      }
      return m;
  }

  if (commflag == 2) {
      for (int i = 0; i < n; i++) {
        int j = list[i];
        buf[m++] = ubuf(npos[j]).d;
      }
      return m;
  }

  if (commflag == 3) {
      for (int i = 0; i < n; i++) {
        int j = list[i];
        for (int k = 0; k < maxbond; k++) {
          buf[m++] = ubuf(partners_possible[j][k]).d;
          buf[m++] = partners_probs[j][k];
        }
      }
      return m;
  }

  if (commflag == 4) {
      for (int i = 0; i < n; i++) {
        int j = list[i];
        for (int k = 0; k < maxbond; k++) {
          buf[m++] = ubuf(partners_success[j][k]).d;
        }
      }
      return m;
  }

  if (commflag == 5) {
    int **nspecial = atom->nspecial;
    tagint **special = atom->special;

    for (int i = 0; i < n; i++) {
        int j = list[i];
        int ns = nspecial[j][0];
        buf[m++] = ubuf(ns).d;
        for (int k = 0; k < ns; k++) {
            buf[m++] = ubuf(special[j][k]).d;
        }
      }
      return m;
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixBondDynamic::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;

  if (commflag == 1) {
    tagint **fbd = atom->iarray[index];

    for (int i = first; i < last; i++) {
        for (int j = 0; j < maxbond; j++) {
          fbd[i][j] = (tagint) ubuf(buf[m++]).i;
        }
    }

  } else if (commflag == 2) {
    for (int i = first; i < last; i++) {
      npos[i] = (int) ubuf(buf[m++]).i;
    }

  } else if (commflag == 3) {
    for (int i = first; i < last; i++) {
        for (int j = 0; j < maxbond; j++) {
          partners_possible[i][j] = (tagint) ubuf(buf[m++]).i;
          partners_probs[i][j] = buf[m++];
        }
    }

  } else if (commflag == 4) {
    for (int i = first; i < last; i++) {
        for (int j = 0; j < maxbond; j++) {
          partners_success[i][j] = (int) ubuf(buf[m++]).i;
        }
    }

  } else if (commflag == 5) {
    int **nspecial = atom->nspecial;
    tagint **special = atom->special;

    for (int i = first; i < last; i++) {
      int ns = (int) ubuf(buf[m++]).i;
      nspecial[i][0] = ns;
      for (int j = 0; j < ns; j++) {
          special[i][j] = (tagint) ubuf(buf[m++]).i;
      }
    }

  }

}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixBondDynamic::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 3*nmax * sizeof(tagint);
  bytes += (double)nmax*3 * sizeof(double);
  bytes += (double)nmax*4 * sizeof(int);
  return bytes;
}