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

#include <cstring>
#include <cstdlib>
#include "fix_drude.h"
#include "atom.h"
#include "comm.h"
#include "modify.h"
#include "error.h"
#include "memory.h"
#include "molecule.h"
#include "atom_vec.h"

#include <set>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDrude::FixDrude(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 3 + atom->ntypes) error->all(FLERR,"Illegal fix drude command");

  comm_border = 1; // drudeid
  special_alter_flag = 1;
  create_attribute = 1;
  is_reduced = false;

  memory->create(drudetype, atom->ntypes+1, "fix_drude::drudetype");
  for (int i=3; i<narg; i++) {
      if (arg[i][0] == 'n' || arg[i][0] == 'N' || arg[i][0] == '0')
          drudetype[i-2] = NOPOL_TYPE;
      else if (arg[i][0] == 'c' || arg[i][0] == 'C' || arg[i][0] == '1')
          drudetype[i-2] = CORE_TYPE;
      else if (arg[i][0] == 'd' || arg[i][0] == 'D' || arg[i][0] == '2')
          drudetype[i-2] = DRUDE_TYPE;
      else
          error->all(FLERR, "Illegal fix drude command");
  }

  drudeid = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);
  atom->add_callback(2);

  // one-time assignment of Drude partners

  build_drudeid();

  // set rebuildflag = 0 to indicate special lists have never been rebuilt

  rebuildflag = 0;
}

/* ---------------------------------------------------------------------- */

FixDrude::~FixDrude()
{
  atom->delete_callback(id,2);
  atom->delete_callback(id,1);
  atom->delete_callback(id,0);
  memory->destroy(drudetype);
  memory->destroy(drudeid);
}

/* ---------------------------------------------------------------------- */

void FixDrude::init()
{
  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"drude") == 0) count++;
  if (count > 1) error->all(FLERR,"More than one fix drude");

  if (!rebuildflag) rebuild_special();
}

/* ---------------------------------------------------------------------- */

int FixDrude::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   look in bond lists for Drude partner tags and fill drudeid
------------------------------------------------------------------------- */
void FixDrude::build_drudeid(){
  int nlocal = atom->nlocal;
  int *type = atom->type;

  std::vector<tagint> drude_vec; // list of my Drudes' tags
  std::vector<tagint> core_drude_vec;
  partner_set = new std::set<tagint>[nlocal]; // Temporary sets of bond partner tags

  if (atom->molecular == 1)
  {
    // Build list of my atoms' bond partners
    for (int i=0; i<nlocal; i++){
      if (drudetype[type[i]] == NOPOL_TYPE) continue;
      drudeid[i] = 0;
      for (int k=0; k<atom->num_bond[i]; k++){
        core_drude_vec.push_back(atom->tag[i]);
        core_drude_vec.push_back(atom->bond_atom[i][k]);
      }
    }
  }
  else
  {
    // Template case
    class Molecule **atommols;
    atommols = atom->avec->onemols;

    // Build list of my atoms' bond partners
    for (int i=0; i<nlocal; i++){
      int imol = atom->molindex[i];
      int iatom = atom->molatom[i];
      tagint *batom = atommols[imol]->bond_atom[iatom];
      tagint tagprev = atom->tag[i] - iatom - 1;
      int nbonds = atommols[imol]->num_bond[iatom];

      if (drudetype[type[i]] == NOPOL_TYPE) continue;
      drudeid[i] = 0;
      for (int k=0; k<nbonds; k++){
        core_drude_vec.push_back(atom->tag[i]);
        core_drude_vec.push_back(batom[k]+tagprev);
      }
    }
  }
  // Loop on procs to fill my atoms' sets of bond partners
  comm->ring(core_drude_vec.size(), sizeof(tagint),
             (char *) core_drude_vec.data(),
             4, ring_build_partner, NULL, (void *)this, 1);

  // Build the list of my Drudes' tags
  // The only bond partners of a Drude particle is its core,
  // so fill drudeid for my Drudes.
  for (int i=0; i<nlocal; i++){
    if (drudetype[type[i]] == DRUDE_TYPE){
      drude_vec.push_back(atom->tag[i]);
      drudeid[i] = *partner_set[i].begin(); // only one 1-2 neighbor, the core
    }
  }
  // At this point each of my Drudes knows its core.
  // Send my list of Drudes to other procs and myself
  // so that each core finds its Drude.
  comm->ring(drude_vec.size(), sizeof(tagint),
             (char *) drude_vec.data(),
             3, ring_search_drudeid, NULL, (void *)this, 1);
  delete [] partner_set;
}

/* ----------------------------------------------------------------------
 * when receive buffer, build the set of received Drude tags.
 * Look in my cores' bond partner tags if there is a Drude tag.
 * If so fill this core's dureid.
------------------------------------------------------------------------- */
void FixDrude::ring_search_drudeid(int size, char *cbuf, void *ptr){
  // Search for the drude partner of my cores
  FixDrude *fdptr = (FixDrude *) ptr;
  Atom *atom = fdptr->atom;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  std::set<tagint> *partner_set = fdptr->partner_set;
  tagint *drudeid = fdptr->drudeid;
  int *drudetype = fdptr->drudetype;

  tagint *first = (tagint *) cbuf;
  tagint *last = first + size;
  std::set<tagint> drude_set(first, last);
  std::set<tagint>::iterator it;

  for (int i=0; i<nlocal; i++) {
    if (drudetype[type[i]] != CORE_TYPE || drudeid[i] > 0) continue;
    for (it = partner_set[i].begin(); it != partner_set[i].end(); it++) { // Drude-core are 1-2 neighbors
      if (drude_set.count(*it) > 0){
        drudeid[i] = *it;
        break;
      }
    }
  }
}

/* ----------------------------------------------------------------------
 * buffer contains bond partners. Look for my atoms and add their partner's
 * tag in its set of bond partners.
------------------------------------------------------------------------- */
void FixDrude::ring_build_partner(int size, char *cbuf, void *ptr){
  // Add partners from incoming list
  FixDrude *fdptr = (FixDrude *) ptr;
  Atom *atom = fdptr->atom;
  int nlocal = atom->nlocal;
  std::set<tagint> *partner_set = fdptr->partner_set;
  tagint *it = (tagint *) cbuf;
  tagint *last = it + size;

  while (it < last) {
    int j = atom->map(*it);
    if (j >= 0 && j < nlocal)
      partner_set[j].insert(*(it+1));
    j = atom->map(*(it+1));
    if (j >= 0 && j < nlocal)
      partner_set[j].insert(*it);
    it += 2;
  }
}


/* ----------------------------------------------------------------------
   allocate atom-based array for drudeid
------------------------------------------------------------------------- */

void FixDrude::grow_arrays(int nmax)
{
  memory->grow(drudeid,nmax,"fix_drude:drudeid");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixDrude::copy_arrays(int i, int j, int /*delflag*/)
{
    drudeid[j] = drudeid[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixDrude::pack_exchange(int i, double *buf)
{
    int m = 0;
    buf[m++] = ubuf(drudeid[i]).d;
    return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixDrude::unpack_exchange(int nlocal, double *buf)
{
    int m = 0;
    drudeid[nlocal] = (tagint) ubuf(buf[m++]).i;
    return m;
}

/* ----------------------------------------------------------------------
   pack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixDrude::pack_border(int n, int *list, double *buf)
{
    int m = 0;
    for (int i=0; i<n; i++){
        int j = list[i];
        buf[m++] = ubuf(drudeid[j]).d;
    }
    return m;
}

/* ----------------------------------------------------------------------
   unpack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixDrude::unpack_border(int n, int first, double *buf)
{
    int m = 0;
    for (int i=first; i<first+n; i++){
        drudeid[i] = (tagint) ubuf(buf[m++]).i;
    }
    return m;
}

/* ----------------------------------------------------------------------
   Rebuild the list of special neighbors if atom_style is Drude
   so that each Drude particle is equivalent to its core atom.
------------------------------------------------------------------------- */

void FixDrude::rebuild_special(){
  rebuildflag = 1;

  int nlocal = atom->nlocal;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int *type = atom->type;

  if (atom->molecular != 1)
    return;

  // Make sure that drude partners know each other
  //build_drudeid();

  // Log info
  if (comm->me == 0) {
    if (screen) fprintf(screen, "Rebuild special list taking Drude particles into account\n");
    if (logfile) fprintf(logfile, "Rebuild special list taking Drude particles into account\n");
  }
  int nspecmax, nspecmax_old, nspecmax_loc;
  nspecmax_loc = 0;
  for (int i=0; i<nlocal; i++) {
    if (nspecmax_loc < nspecial[i][2]) nspecmax_loc = nspecial[i][2];
  }
  MPI_Allreduce(&nspecmax_loc, &nspecmax_old, 1, MPI_INT, MPI_MAX, world);
  if (comm->me == 0) {
    if (screen) fprintf(screen, "Old max number of 1-2 to 1-4 neighbors: %d\n", nspecmax_old);
    if (logfile) fprintf(logfile, "Old max number of 1-2 to 1-4 neighbors: %d\n", nspecmax_old);
  }

  // Build lists of drude and core-drude pairs
  std::vector<tagint> drude_vec, core_drude_vec, core_special_vec;
  for (int i=0; i<nlocal; i++) {
    if (drudetype[type[i]] == DRUDE_TYPE) {
      drude_vec.push_back(atom->tag[i]);
    } else if (drudetype[type[i]] == CORE_TYPE){
      core_drude_vec.push_back(atom->tag[i]);
      core_drude_vec.push_back(drudeid[i]);
    }
  }
  // Remove Drude particles from the special lists of each proc
  comm->ring(drude_vec.size(), sizeof(tagint),
             (char *) drude_vec.data(),
             9, ring_remove_drude, NULL, (void *)this, 1);
  // Add back Drude particles in the lists just after their core
  comm->ring(core_drude_vec.size(), sizeof(tagint),
             (char *) core_drude_vec.data(),
             10, ring_add_drude, NULL, (void *)this, 1);

  // Check size of special list
  nspecmax_loc = 0;
  for (int i=0; i<nlocal; i++) {
    if (nspecmax_loc < nspecial[i][2]) nspecmax_loc = nspecial[i][2];
  }
  MPI_Allreduce(&nspecmax_loc, &nspecmax, 1, MPI_INT, MPI_MAX, world);
  if (comm->me == 0) {
    if (screen) fprintf(screen, "New max number of 1-2 to 1-4 neighbors: %d (+%d)\n", nspecmax, nspecmax - nspecmax_old);
    if (logfile) fprintf(logfile, "New max number of 1-2 to 1-4 neighbors: %d (+%d)\n", nspecmax, nspecmax - nspecmax_old);
  }
  if (atom->maxspecial < nspecmax) {
    char str[1024];
    sprintf(str, "Not enough space in special: extra/special/per/atom should be at least %d", nspecmax - nspecmax_old);
    error->all(FLERR, str);
  }

  // Build list of cores' special lists to communicate to ghost drude particles
  for (int i=0; i<nlocal; i++) {
    if (drudetype[type[i]] != CORE_TYPE) continue;
    core_special_vec.push_back(atom->tag[i]);
    core_special_vec.push_back((tagint) nspecial[i][0]);
    core_special_vec.push_back((tagint) nspecial[i][1]);
    core_special_vec.push_back((tagint) nspecial[i][2]);
    for (int j=1; j<nspecial[i][2]; j++)
      core_special_vec.push_back(special[i][j]);
  }
  // Copy core's list into their drude list
  comm->ring(core_special_vec.size(), sizeof(tagint),
             (char *) core_special_vec.data(),
             11, ring_copy_drude, NULL, (void *)this, 1);
}

/* ----------------------------------------------------------------------
 * When receive buffer, build a set of drude tags, look into my atoms'
 * special list if some tags are drude particles. If so, remove it.
------------------------------------------------------------------------- */
void FixDrude::ring_remove_drude(int size, char *cbuf, void *ptr){
  // Remove all drude particles from special list
  FixDrude *fdptr = (FixDrude *) ptr;
  Atom *atom = fdptr->atom;
  int nlocal = atom->nlocal;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int *type = atom->type;
  tagint *first = (tagint *) cbuf;
  tagint *last = first + size;
  std::set<tagint> drude_set(first, last);
  int *drudetype = fdptr->drudetype;

  for (int i=0; i<nlocal; i++) {
    if (drudetype[type[i]] == DRUDE_TYPE) continue;
    for (int j=0; j<nspecial[i][2]; j++) {
      if (drude_set.count(special[i][j]) > 0) { // I identify a drude in my special list, remove it
        // left shift
        nspecial[i][2]--;
        for (int k=j; k<nspecial[i][2]; k++)
          special[i][k] = special[i][k+1];
        if (j < nspecial[i][1]) {
          nspecial[i][1]--;
          if (j < nspecial[i][0]) nspecial[i][0]--;
        }
        j--;
      }
    }
  }
}

/* ----------------------------------------------------------------------
 * When receive buffer, build a map core tag -> drude tag.
 * Loop on my atoms' special list to find core tags. Insert their Drude
 * particle if they have one.
------------------------------------------------------------------------- */
void FixDrude::ring_add_drude(int size, char *cbuf, void *ptr){
  // Assume special array size is big enough
  // Add all particle just after their core in the special list
  FixDrude *fdptr = (FixDrude *) ptr;
  Atom *atom = fdptr->atom;
  int nlocal = atom->nlocal;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int *type = atom->type;
  tagint *drudeid = fdptr->drudeid;
  int *drudetype = fdptr->drudetype;

  tagint *first = (tagint *) cbuf;
  tagint *last = first + size;
  std::map<tagint, tagint> core_drude_map;

  tagint *it = first;
  while (it < last) {
    tagint core_tag = *it;
    it++;
    core_drude_map[core_tag] = *it;
    it++;
  }

  for (int i=0; i<nlocal; i++) {
    if (drudetype[type[i]] == DRUDE_TYPE) continue;
    if (core_drude_map.count(atom->tag[i]) > 0) { // I identify myself as a core, add my own drude
      // right shift
      for (int k=nspecial[i][2]-1; k>=0; k--)
        special[i][k+1] = special[i][k];
      special[i][0] = drudeid[i];
      nspecial[i][0]++;
      nspecial[i][1]++;
      nspecial[i][2]++;
    }
    for (int j=0; j<nspecial[i][2]; j++) {
      if (core_drude_map.count(special[i][j]) > 0) { // I identify a core in my special list, add his drude
        // right shift
        for (int k=nspecial[i][2]-1; k>j; k--)
          special[i][k+1] = special[i][k];
        special[i][j+1] = core_drude_map[special[i][j]];
        nspecial[i][2]++;
        if (j < nspecial[i][1]) {
          nspecial[i][1]++;
          if (j < nspecial[i][0]) nspecial[i][0]++;
        }
        j++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
 * When receive buffer, build a map core tag -> pointer to special info
 * in the buffer. Loop on my Drude particles and copy their special
 * info from that of their core if the latter is found in the map.
------------------------------------------------------------------------- */
void FixDrude::ring_copy_drude(int size, char *cbuf, void *ptr){
  // Copy special list of drude from its core (except itself)
  FixDrude *fdptr = (FixDrude *) ptr;
  Atom *atom = fdptr->atom;
  int nlocal = atom->nlocal;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int *type = atom->type;
  tagint *drudeid = fdptr->drudeid;
  int *drudetype = fdptr->drudetype;

  tagint *first = (tagint *) cbuf;
  tagint *last = first + size;
  std::map<tagint, tagint*> core_special_map;

  tagint *it = first;
  while (it < last) {
    tagint core_tag = *it;
    it++;
    core_special_map[core_tag] = it;
    it += 2;
    it += (int) *it;
  }

  for (int i=0; i<nlocal; i++) {
    if (drudetype[type[i]] != DRUDE_TYPE) continue;
    if (core_special_map.count(drudeid[i]) > 0) { // My core is in this list, copy its special info
      it = core_special_map[drudeid[i]];
      nspecial[i][0] = (int) *it;
      it++;
      nspecial[i][1] = (int) *it;
      it++;
      nspecial[i][2] = (int) *it;
      it++;
      special[i][0] = drudeid[i];
      for (int k=1; k<nspecial[i][2]; k++) {
        special[i][k] = *it;
        it++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
 * Set drudeid when a new atom is created,
 * special list must be up-to-date
 * ----------------------------------------------------------------------*/
void FixDrude::set_arrays(int i){
    if (drudetype[atom->type[i]] != NOPOL_TYPE){
        if (atom->nspecial[i] ==0) error->all(FLERR, "Polarizable atoms cannot be inserted with special lists info from the molecule template");
        drudeid[i] = atom->special[i][0]; // Drude partner should be at first place in the special list
    } else {
        drudeid[i] = 0;
    }
}

