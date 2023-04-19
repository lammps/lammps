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

#include "fix_bond_relocate.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "citeme.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_fix_bond_swap[] =
    "fix bond/relocate command: doi:10.1063/1.1628670\n\n"
    "@Article{Auhl03,\n"
    " author = {R. Auhl and R. Everaers and G. S. Grest and K. Kremer and S. J. Plimpton},\n"
    " title = {Equilibration of Long Chain Polymer Melts in Computer Simulations},\n"
    " journal = {J.~Chem.\\ Phys.},\n"
    " year =    2003,\n"
    " volume =  119,\n"
    " number =  12,\n"
    " pages =   {12718--12728}\n"
    "}\n\n";

#define DELTA_PERMUTE 100

/* ---------------------------------------------------------------------- */

FixBondRelocate::FixBondRelocate(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), tflag(0), randomized_eligible_sourceatoms(nullptr),
    randomized_eligible_targetatoms(nullptr), id_temp(nullptr), type(nullptr), x(nullptr),
    list(nullptr), temperature(nullptr), random(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_bond_swap);

  if (narg != 11) error->all(FLERR, "Illegal fix bond/relocate command");

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, "Illegal fix bond/relocate command");

  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;

  relocate_btype = utils::inumeric(FLERR, arg[4], false, lmp);
  decide_btype = utils::inumeric(FLERR, arg[5], false, lmp);

  double low_cutoff = utils::numeric(FLERR, arg[6], false, lmp);
  double high_cutoff = utils::numeric(FLERR, arg[7], false, lmp);

  locutsq = low_cutoff * low_cutoff;
  hicutsq = high_cutoff * high_cutoff;

  fraction = utils::numeric(FLERR, arg[8], false, lmp);
  decide_functionality = utils::inumeric(FLERR, arg[10], false, lmp);

  // initialize Marsaglia RNG with processor-unique seed

  int seed = utils::inumeric(FLERR, arg[9], false, lmp);
  random = new RanMars(lmp, seed + comm->me);

  // error check

  if (atom->molecular != Atom::MOLECULAR)
    error->all(FLERR, "Cannot use fix bond/relocate with non-molecular systems");

  // create a new compute temp style
  // id = fix-ID + temp, compute group = fix group

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} all temp", id_temp));
  tflag = 1;

  // initialize two permutation lists

  nmax = 0;
  randomized_eligible_sourceatoms = nullptr;
  randomized_eligible_targetatoms = nullptr;

  maxpermute = 0;
  permute = nullptr;

  naccept = foursome = 0;
}

/* ---------------------------------------------------------------------- */

FixBondRelocate::~FixBondRelocate()
{
  delete random;

  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete[] id_temp;

  memory->destroy(randomized_eligible_sourceatoms);
  memory->destroy(randomized_eligible_targetatoms);
  delete[] permute;
}

/* ---------------------------------------------------------------------- */

int FixBondRelocate::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondRelocate::init()
{
  // require an atom style with molecule IDs

  if (atom->molecule == nullptr)
    error->all(FLERR, "Must use atom style with molecule IDs with fix bond/relocate");

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0) error->all(FLERR, "Temperature ID for fix bond/relocate does not exist");
  temperature = modify->compute[icompute];

  // pair and bonds must be defined
  // no dihedral or improper potentials allowed
  // special bonds must be 0 1 1

  if (force->pair == nullptr || force->bond == nullptr)
    error->all(FLERR, "Fix bond/relocate requires pair and bond styles");

  if (force->pair->single_enable == 0)
    error->all(FLERR, "Pair style does not support fix bond/relocate");

  if (force->angle == nullptr && atom->nangles > 0 && comm->me == 0)
    error->warning(FLERR,
                   "Fix bond/relocate will not preserve correct angle "
                   "topology because no angle_style is defined");

  if (force->dihedral || force->improper)
    error->all(FLERR, "Fix bond/relocate cannot use dihedral or improper styles");

  if (force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0)
    error->all(FLERR, "Fix bond/relocate requires special_bonds = {any},1,1");

  // need a half neighbor list, built every Nevery steps

  neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);

  // zero out stats

  naccept = foursome = 0;
  angleflag = 0;
  if (force->angle) angleflag = 1;
}

/* ---------------------------------------------------------------------- */

void FixBondRelocate::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   look for and perform swaps
   NOTE: used to do this every pre_neighbor(), but think that is a bug
         b/c was doing it after exchange() and before neighbor->build()
         which is when neigh lists are actually out-of-date or even bogus,
         now do it based on user-specified Nevery, and trigger reneigh
         if any swaps performed, like fix bond/create
------------------------------------------------------------------------- */

void FixBondRelocate::post_integrate()
{
  int ii, jj, m, inum, jnum;
  int inext, iprev, ilast, jnext, jprev, jlast, ibond, iangle, jbond, jangle;
  int ibondtype, jbondtype, iangletype, inextangletype, jangletype, jnextangletype;
  tagint itag, inexttag, iprevtag, ilasttag, jtag, jnexttag, jprevtag, jlasttag;
  tagint i1, i2, i3, j1, j2, j3;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double delta, factor;

  if (update->ntimestep % nevery) return;

  // compute current temp for Boltzmann factor test

  double t_current = temperature->compute_scalar();

  // local ptrs to atom arrays

  tagint *tag = atom->tag;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *num_angle = atom->num_angle;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int newton_bond = force->newton_bond;
  int nlocal = atom->nlocal;

  type = atom->type;
  x = atom->x;

  neighbor->build_one(list, 1);
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // randomize list of my owned atoms that are in fix group
  // grow atom list if necessary

  if (atom->nmax > nmax) {
    memory->destroy(randomized_eligible_sourceatoms);
    memory->destroy(randomized_eligible_targetatoms);
    nmax = atom->nmax;
    memory->create(randomized_eligible_sourceatoms, nmax, "bondswap:alist");
    memory->create(randomized_eligible_targetatoms, nmax, "bondswap:alist2");
  }

  // use randomized permutation of both I and J atoms in double loop below
  // this is to avoid any bias in accepted MC swaps based on
  //   ordering LAMMPS creates on a processor for atoms or their neighbors

  // create a random permutation of list of Neligible atoms
  // uses one-pass Fisher-Yates shuffle on an initial identity permutation
  // output: randomized randomized_eligible_sourceatoms[] vector, used in outer loop to select an I atom
  // similar randomized permutation is created for neighbors of each I atom

  int neligible_source = 0;
  int neligible_target = 0;
  for (ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (!(mask[i] & groupbit)) { continue; }
    bool is_eligible_source = false;
    bool is_eligible_target = false;

    int num_decide_bonds = 0;
    for (int bond_i = 0; bond_i < num_bond[i]; ++bond_i) {
      if (bond_type[i][bond_i] == decide_btype) { num_decide_bonds++; }
      if (bond_type[i][bond_i] == relocate_btype) { is_eligible_source = true; }
    }

    is_eligible_target = (num_decide_bonds == decide_functionality);

    if (is_eligible_source && is_eligible_target) {
      randomized_eligible_sourceatoms[neligible_source++] = i;
    }
    if (is_eligible_target) { randomized_eligible_targetatoms[neligible_target++] = i; }
  }

  int tmp;
  for (int i = 0; i < neligible_source - 1; i++) {
    int j = i + static_cast<int>(random->uniform() * (neligible_source - i));
    tmp = randomized_eligible_sourceatoms[i];
    randomized_eligible_sourceatoms[i] = randomized_eligible_sourceatoms[j];
    randomized_eligible_sourceatoms[j] = tmp;
  }
  for (int i = 0; i < neligible_target - 1; i++) {
    int j = i + static_cast<int>(random->uniform() * (neligible_target - i));
    tmp = randomized_eligible_targetatoms[i];
    randomized_eligible_targetatoms[i] = randomized_eligible_targetatoms[j];
    randomized_eligible_targetatoms[j] = tmp;
  }

  // examine ntest of my eligible atoms for potential swaps
  // atom I is randomly selected via atom list
  // atom I's bond of type relocate_btype is going to be possibly removed.
  // select another atom

  int ntest = static_cast<int>(fraction * neligible_source);
  int accept = 0;

  for (int itest_source = 0; itest_source < ntest; itest_source++) {
    int selected_source_a = randomized_eligible_sourceatoms[itest_source];

    for (int itest_target = 0; itest_target < neligible_target; ++itest_target) {
      int test_target_a1 = randomized_eligible_targetatoms[itest_target];
      jlist = firstneigh[test_target_a1];
      jnum = numneigh[test_target_a1];

      // look at all J neighbors of atom test_target_a1
      //   done in a randomized permutation, via neighbor_permutation()
      // J must be on-processor (J < nlocal)
      // I,J must be in fix group
      // I,J must have exactly decide_functionality

      neighbor_permutation(jnum);

      for (jj = 0; jj < jnum; jj++) {
        int test_target_a2 = jlist[permute[jj]];
        test_target_a2 &= NEIGHMASK;
        bool accept_move = false;
        if (test_target_a2 >= nlocal) continue;
        if ((mask[test_target_a2] & groupbit) == 0) continue;

        if (dist_rsq(test_target_a1, test_target_a2) >= hicutsq) continue;
        if (dist_rsq(test_target_a1, test_target_a2) < locutsq) continue;

        // potential new bond to form found: between test_target_a1 and test_target_a2

        // TODO: correct these
        delta = pair_eng(test_target_a1, inext) + pair_eng(test_target_a2, jnext) -
            pair_eng(test_target_a1, jnext) - pair_eng(inext, test_target_a2);
        delta += bond_eng(ibondtype, test_target_a1, jnext) +
            bond_eng(jbondtype, test_target_a2, inext) -
            bond_eng(ibondtype, test_target_a1, inext) - bond_eng(jbondtype, test_target_a2, jnext);
        if (angleflag)
          delta += angle_eng(iangletype, iprev, test_target_a1, jnext) +
              angle_eng(jnextangletype, test_target_a1, jnext, jlast) +
              angle_eng(jangletype, jprev, test_target_a2, inext) +
              angle_eng(inextangletype, test_target_a2, inext, ilast) -
              angle_eng(iangletype, iprev, test_target_a1, inext) -
              angle_eng(inextangletype, test_target_a1, inext, ilast) -
              angle_eng(jangletype, jprev, test_target_a2, jnext) -
              angle_eng(jnextangletype, test_target_a2, jnext, jlast);

        if (delta < 0.0)
          accept_move = true;
        else {
          factor = exp(-delta / force->boltz / t_current);
          if (random->uniform() < factor) accept_move = true;
        }
        if (!accept_move) { continue; }
        // assemble the bond to break
        bool atom_has_additional_bonds_to_relocate = false;
        int bond_i_to_remove = -1;
        for (int bond_i = 0; bond_i < num_bond[selected_source_a]; ++bond_i) {
          if (bond_type[selected_source_a][bond_i] == relocate_btype) {
            if (bond_i_to_remove >= 0) { atom_has_additional_bonds_to_relocate = true; }
            bond_i_to_remove = bond_i;
          }
        }
        // delete the bond if Metropolis criterion is fulfilled
        this->break_bond(selected_source_a, bond_i_to_remove);
        // create the new bond if Metropolis criterion is fulfilled
        this->create_bond(test_target_a1, test_target_a2);
        // check how to continue with this source atom – does it have additional candidates?
        if (atom_has_additional_bonds_to_relocate) {
          goto one_bond_relocated;
        } else {
          goto all_bonds_of_atom_relocated;
        }
      }
    one_bond_relocated:
    }
  all_bonds_of_atom_relocated:
  }

done:

  // trigger immediate reneighboring if swaps occurred on one or more procs

  int accept_any;
  MPI_Allreduce(&accept, &accept_any, 1, MPI_INT, MPI_SUM, world);
  if (accept_any) next_reneighbor = update->ntimestep;

  if (!accept) return;
  naccept++;
}

/* ---------------------------------------------------------------------- */

int FixBondRelocate::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "temp") == 0) {
    if (narg < 2) error->all(FLERR, "Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete[] id_temp;
    id_temp = utils::strdup(arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) error->all(FLERR, "Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,
                 "Fix_modify temperature ID does not "
                 "compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR, "Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- 
  Break the bond on source_atom with index bond_i
------------------------------------------------------------------------- */
void FixBondRelocate::break_bond(int source_atom, int bond_i)
{
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;
  int bond_target_atom_tag = bond_atom[source_atom][bond_i];
  int bond_source_atom_tag = atom->tag[source_atom];
  int bond_target_atom = atom->map(bond_target_atom_tag);
  auto histories = modify->get_fix_by_style("BOND_HISTORY");
  int n_histories = histories.size();

  if (n_histories > 0) {
    for (auto &ihistory : histories) {
      // TODO: check if need to use tags instead
      dynamic_cast<FixBondHistory *>(ihistory)->delete_history(source_atom, bond_target_atom);
    }
  }

  // first, remove it from the bonds of both contributors
  for (int i = num_bond[source_atom] - 2; i >= bond_i; i--) {
    bond_atom[source_atom][i] = bond_atom[source_atom][i + 1];
  }
  num_bond[source_atom]--;

  int bond_i_target = -1;
  int num_bonds_to_source = 0;
  for (int i = 0; i < num_bond[bond_target_atom]; ++i) {
    if (bond_atom[bond_target_atom][i] == bond_source_atom_tag) {
      if (bond_i_target >= 0) { num_bonds_to_source++; }
      if (bond_type[bond_target_atom][i] == relocate_btype) { bond_i_target = i; }
    }
  }
  for (int i = num_bond[bond_target_atom] - 2; i >= bond_i_target; i--) {
    bond_atom[bond_target_atom][i] = bond_atom[bond_target_atom][i + 1];
  }
  num_bond[bond_target_atom]--;

  // then, adjust the special list
  if (num_bonds_to_source == 1) {
    // TODO: adjust special list
  }
  // then, TODO: check for angles and stuff
}

/* ---------------------------------------------------------------------- 
  Create a bond between source_atom and target_atom
------------------------------------------------------------------------- */
void FixBondRelocate::create_bond(int source_atom, int target_atom) {}

/* ----------------------------------------------------------------------
   compute squared distance between atoms I,J
   must use minimum_image since J was found thru atom->map()
------------------------------------------------------------------------- */

double FixBondRelocate::dist_rsq(int i, int j)
{
  double delx = x[i][0] - x[j][0];
  double dely = x[i][1] - x[j][1];
  double delz = x[i][2] - x[j][2];
  domain->minimum_image(delx, dely, delz);
  return (delx * delx + dely * dely + delz * delz);
}

/* ----------------------------------------------------------------------
   return pairwise interaction energy between atoms I,J
   will always be full non-bond interaction, so factors = 1 in single() call
------------------------------------------------------------------------- */

double FixBondRelocate::pair_eng(int i, int j)
{
  double tmp;
  double rsq = dist_rsq(i, j);
  return force->pair->single(i, j, type[i], type[j], rsq, 1.0, 1.0, tmp);
}

/* ---------------------------------------------------------------------- */

double FixBondRelocate::bond_eng(int btype, int i, int j)
{
  double tmp;
  double rsq = dist_rsq(i, j);
  return force->bond->single(btype, rsq, i, j, tmp);
}

/* ---------------------------------------------------------------------- */

double FixBondRelocate::angle_eng(int atype, int i, int j, int k)
{
  // test for non-existent angle at end of chain

  if (i == -1 || k == -1) return 0.0;
  return force->angle->single(atype, i, j, k);
}

/* ----------------------------------------------------------------------
   create a random permutation of one atom's N neighbor list atoms
   uses one-pass Fisher-Yates shuffle on an initial identity permutation
   output: randomized permute[] vector, used to index neighbors
------------------------------------------------------------------------- */

void FixBondRelocate::neighbor_permutation(int n)
{
  int i, j, tmp;

  if (n > maxpermute) {
    delete[] permute;
    maxpermute = n + DELTA_PERMUTE;
    permute = new int[maxpermute];
  }

  // Fisher-Yates shuffle

  for (i = 0; i < n; i++) permute[i] = i;

  for (i = 0; i < n - 1; i++) {
    j = i + static_cast<int>(random->uniform() * (n - i));
    tmp = permute[i];
    permute[i] = permute[j];
    permute[j] = tmp;
  }
}

/* ----------------------------------------------------------------------
   return bond swapping stats
   n = 1 is # of swaps
   n = 2 is # of attempted swaps
------------------------------------------------------------------------- */

double FixBondRelocate::compute_vector(int n)
{
  double one, all;
  if (n == 0)
    one = naccept;
  else
    one = foursome;
  MPI_Allreduce(&one, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  return all;
}

/* ----------------------------------------------------------------------
   memory usage of randomized_eligible_sourceatoms
------------------------------------------------------------------------- */

double FixBondRelocate::memory_usage()
{
  double bytes = (double) nmax * sizeof(int);
  return bytes;
}
