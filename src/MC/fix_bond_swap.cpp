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

#include "fix_bond_swap.h"

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
  "fix bond/swap command: doi:10.1063/1.1628670\n\n"
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

FixBondSwap::FixBondSwap(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  tflag(0), alist(nullptr), id_temp(nullptr), type(nullptr), x(nullptr), list(nullptr),
  temperature(nullptr), random(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_bond_swap);

  if (narg != 7) error->all(FLERR,"Illegal fix bond/swap command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix bond/swap command");

  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;

  fraction = utils::numeric(FLERR,arg[4],false,lmp);
  double cutoff = utils::numeric(FLERR,arg[5],false,lmp);
  cutsq = cutoff*cutoff;

  // initialize Marsaglia RNG with processor-unique seed

  int seed = utils::inumeric(FLERR,arg[6],false,lmp);
  random = new RanMars(lmp,seed + comm->me);

  // error check

  if (atom->molecular != Atom::MOLECULAR)
    error->all(FLERR,"Cannot use fix bond/swap with non-molecular systems");

  // create a new compute temp style
  // id = fix-ID + temp, compute group = fix group

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} all temp",id_temp));
  tflag = 1;

  // initialize two permutation lists

  nmax = 0;
  alist = nullptr;

  maxpermute = 0;
  permute = nullptr;

  naccept = foursome = 0;
}

/* ---------------------------------------------------------------------- */

FixBondSwap::~FixBondSwap()
{
  delete random;

  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete[] id_temp;

  memory->destroy(alist);
  delete[] permute;
}

/* ---------------------------------------------------------------------- */

int FixBondSwap::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondSwap::init()
{
  // require an atom style with molecule IDs

  if (atom->molecule == nullptr)
    error->all(FLERR,
               "Must use atom style with molecule IDs with fix bond/swap");

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix bond/swap does not exist");
  temperature = modify->compute[icompute];

  // pair and bonds must be defined
  // no dihedral or improper potentials allowed
  // special bonds must be 0 1 1

  if (force->pair == nullptr || force->bond == nullptr)
    error->all(FLERR,"Fix bond/swap requires pair and bond styles");

  if (force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support fix bond/swap");

  if (force->angle == nullptr && atom->nangles > 0 && comm->me == 0)
    error->warning(FLERR,"Fix bond/swap will not preserve correct angle "
                   "topology because no angle_style is defined");

  if (force->dihedral || force->improper)
    error->all(FLERR,"Fix bond/swap cannot use dihedral or improper styles");

  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 ||
      force->special_lj[3] != 1.0)
    error->all(FLERR,"Fix bond/swap requires special_bonds = 0,1,1");

  // need a half neighbor list, built every Nevery steps

  neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);

  // zero out stats

  naccept = foursome = 0;
  angleflag = 0;
  if (force->angle) angleflag = 1;
}

/* ---------------------------------------------------------------------- */

void FixBondSwap::init_list(int /*id*/, NeighList *ptr)
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

void FixBondSwap::post_integrate()
{
  int i,j,ii,jj,m,inum,jnum;
  int inext,iprev,ilast,jnext,jprev,jlast,ibond,iangle,jbond,jangle;
  int ibondtype,jbondtype,iangletype,inextangletype,jangletype,jnextangletype;
  tagint itag,inexttag,iprevtag,ilasttag,jtag,jnexttag,jprevtag,jlasttag;
  tagint i1,i2,i3,j1,j2,j3;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double delta,factor;

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

  neighbor->build_one(list,1);
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // randomize list of my owned atoms that are in fix group
  // grow atom list if necessary

  if (atom->nmax > nmax) {
    memory->destroy(alist);
    nmax = atom->nmax;
    memory->create(alist,nmax,"bondswap:alist");
  }

  // use randomized permutation of both I and J atoms in double loop below
  // this is to avoid any bias in accepted MC swaps based on
  //   ordering LAMMPS creates on a processor for atoms or their neighbors

  // create a random permutation of list of Neligible atoms
  // uses one-pass Fisher-Yates shuffle on an initial identity permutation
  // output: randomized alist[] vector, used in outer loop to select an I atom
  // similar randomized permutation is created for neighbors of each I atom

  int neligible = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit)
      alist[neligible++] = i;
  }

  int tmp;
  for (i = 0; i < neligible-1; i++) {
    j = i + static_cast<int> (random->uniform() * (neligible-i));
    tmp = alist[i];
    alist[i] = alist[j];
    alist[j] = tmp;
  }

  // examine ntest of my eligible atoms for potential swaps
  // atom I is randomly selected via atom list
  // look at all J neighbors of atom I
  //   done in a randomized permutation, via neighbor_permutation()
  // J must be on-processor (J < nlocal)
  // I,J must be in fix group
  // I,J must have same molecule IDs
  //   use case 1 (see doc page):
  //     if user defines mol IDs appropriately for linear chains,
  //     this will mean they are same distance from (either) chain end
  //   use case 2 (see doc page):
  //     if user defines a unique mol ID for desired bond sites (on any chain)
  //     and defines the fix group as these sites,
  //     this will mean they are eligible bond sites

  int ntest = static_cast<int> (fraction * neligible);
  int accept = 0;

  for (int itest = 0; itest < ntest; itest++) {
    i = alist[itest];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    neighbor_permutation(jnum);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[permute[jj]];
      j &= NEIGHMASK;
      if (j >= nlocal) continue;
      if ((mask[j] & groupbit) == 0) continue;
      if (molecule[i] != molecule[j]) continue;

      // loop over all bond partners of atoms I and J
      // use num_bond for this, not special list, so also have bondtypes
      // inext,jnext = atoms bonded to I,J
      // inext,jnext must be on-processor (inext,jnext < nlocal)
      // inext,jnext must be in fix group
      // inext,jnext must have same molecule IDs
      //   in use cases above ...
      //   for case 1: this ensures chain length is preserved
      //   for case 2: always satisfied b/c fix group = bond-able atoms
      // 4 atoms must be unique (no duplicates): inext != jnext, inext != j
      //   already know i != inext, j != jnext
      // all 4 old and new bonds must have length < cutoff

      for (ibond = 0; ibond < num_bond[i]; ibond++) {
        inext = atom->map(bond_atom[i][ibond]);
        if (inext >= nlocal || inext < 0) continue;
        if ((mask[inext] & groupbit) == 0) continue;
        ibondtype = bond_type[i][ibond];

        for (jbond = 0; jbond < num_bond[j]; jbond++) {
          jnext = atom->map(bond_atom[j][jbond]);
          if (jnext >= nlocal || jnext < 0) continue;
          if ((mask[jnext] & groupbit) == 0) continue;
          jbondtype = bond_type[j][jbond];

          if (molecule[inext] != molecule[jnext]) continue;
          if (inext == jnext || inext == j) continue;
          if (dist_rsq(i,inext) >= cutsq) continue;
          if (dist_rsq(j,jnext) >= cutsq) continue;
          if (dist_rsq(i,jnext) >= cutsq) continue;
          if (dist_rsq(j,inext) >= cutsq) continue;

          // if angles are enabled:
          // find other atoms i,inext,j,jnext are in angles with
          //   and angletypes: i/j angletype, i/j nextangletype
          // use num_angle for this, not special list, so also have angletypes
          // 4 atoms consecutively along 1st chain: iprev,i,inext,ilast
          // 4 atoms consecutively along 2nd chain: jprev,j,jnext,jlast
          // prev or last atom can be non-existent at end of chain
          //   set prev/last = -1 in this case
          // if newton bond = 0, then angles are stored by all 4 atoms
          //   so require that iprev,ilast,jprev,jlast be owned by this proc
          //   so all copies of angles can be updated if a swap takes place

          if (angleflag) {
            itag = tag[i];
            inexttag = tag[inext];
            jtag = tag[j];
            jnexttag = tag[jnext];

            iprev = -1;
            for (iangle = 0; iangle < num_angle[i]; iangle++) {
              i1 = angle_atom1[i][iangle];
              i2 = angle_atom2[i][iangle];
              i3 = angle_atom3[i][iangle];
              if (i2 == itag && i3 == inexttag) iprev = atom->map(i1);
              else if (i1 == inexttag && i2 == itag) iprev = atom->map(i3);
              if (iprev >= 0) {
                iangletype = angle_type[i][iangle];
                break;
              }
            }
            if (!newton_bond && iprev >= nlocal) continue;

            ilast = -1;
            for (iangle = 0; iangle < num_angle[inext]; iangle++) {
              i1 = angle_atom1[inext][iangle];
              i2 = angle_atom2[inext][iangle];
              i3 = angle_atom3[inext][iangle];
              if (i1 == itag && i2 == inexttag) ilast = atom->map(i3);
              else if (i2 == inexttag && i3 == itag) ilast = atom->map(i1);
              if (ilast >= 0) {
                inextangletype = angle_type[inext][iangle];
                break;
              }
            }
            if (!newton_bond && ilast >= nlocal) continue;

            jprev = -1;
            for (jangle = 0; jangle < num_angle[j]; jangle++) {
              j1 = angle_atom1[j][jangle];
              j2 = angle_atom2[j][jangle];
              j3 = angle_atom3[j][jangle];
              if (j2 == jtag && j3 == jnexttag) jprev = atom->map(j1);
              else if (j1 == jnexttag && j2 == jtag) jprev = atom->map(j3);
              if (jprev >= 0) {
                jangletype = angle_type[j][jangle];
                break;
              }
            }
            if (!newton_bond && jprev >= nlocal) continue;

            jlast = -1;
            for (jangle = 0; jangle < num_angle[jnext]; jangle++) {
              j1 = angle_atom1[jnext][jangle];
              j2 = angle_atom2[jnext][jangle];
              j3 = angle_atom3[jnext][jangle];
              if (j1 == jtag && j2 == jnexttag) jlast = atom->map(j3);
              else if (j2 == jnexttag && j3 == jtag) jlast = atom->map(j1);
              if (jlast >= 0) {
                jnextangletype = angle_type[jnext][jangle];
                break;
              }
            }
            if (!newton_bond && jlast >= nlocal) continue;
          }

          // valid foursome found between 2 chains:
          //   chains = iprev-i-inext-ilast and jprev-j-jnext-jlast
          //   prev or last values are -1 if do not exist due to end of chain
          //   OK to call angle_eng with -1 atom, since just return 0.0
          // current energy of foursome =
          //   E_nb(i,j) + E_nb(i,jnext) + E_nb(inext,j) + E_nb(inext,jnext) +
          //   E_bond(i,inext) + E_bond(j,jnext) +
          //   E_angle(iprev,i,inext) + E_angle(i,inext,ilast) +
          //   E_angle(jprev,j,jnext) + E_angle(j,jnext,jlast)
          // new energy of foursome with swapped bonds =
          //   E_nb(i,j) + E_nb(i,inext) + E_nb(j,jnext) + E_nb(inext,jnext) +
          //   E_bond(i,jnext) + E_bond(j,inext) +
          //   E_angle(iprev,i,jnext) + E_angle(i,jnext,jlast) +
          //   E_angle(jprev,j,inext) + E_angle(j,inext,ilast)
          // energy delta = add/subtract differing terms between 2 formulas

          foursome++;

          delta = pair_eng(i,inext) + pair_eng(j,jnext) -
            pair_eng(i,jnext) - pair_eng(inext,j);
          delta += bond_eng(ibondtype,i,jnext) + bond_eng(jbondtype,j,inext) -
            bond_eng(ibondtype,i,inext) - bond_eng(jbondtype,j,jnext);
          if (angleflag)
            delta += angle_eng(iangletype,iprev,i,jnext) +
              angle_eng(jnextangletype,i,jnext,jlast) +
              angle_eng(jangletype,jprev,j,inext) +
              angle_eng(inextangletype,j,inext,ilast) -
              angle_eng(iangletype,iprev,i,inext) -
              angle_eng(inextangletype,i,inext,ilast) -
              angle_eng(jangletype,jprev,j,jnext) -
              angle_eng(jnextangletype,j,jnext,jlast);

          // if delta <= 0, accept swap
          // if delta > 0, compute Boltzmann factor with current temperature
          //   only accept if greater than random value
          // whether accept or not, exit test loop

          if (delta < 0.0) accept = 1;
          else {
            factor = exp(-delta/force->boltz/t_current);
            if (random->uniform() < factor) accept = 1;
          }

          goto done;
        }
      }
    }
  }

 done:

  // trigger immediate reneighboring if swaps occurred on one or more procs

  int accept_any;
  MPI_Allreduce(&accept,&accept_any,1,MPI_INT,MPI_SUM,world);
  if (accept_any) next_reneighbor = update->ntimestep;

  if (!accept) return;
  naccept++;

  // find instances of bond/history to reset history
  auto histories = modify->get_fix_by_style("BOND_HISTORY");
  int n_histories = histories.size();

  // change bond partners of affected atoms
  // on atom i: bond i-inext changes to i-jnext
  // on atom j: bond j-jnext changes to j-inext
  // on atom inext: bond inext-i changes to inext-j
  // on atom jnext: bond jnext-j changes to jnext-i

  for (ibond = 0; ibond < num_bond[i]; ibond++)
    if (bond_atom[i][ibond] == tag[inext]) {
      if (n_histories > 0)
        for (auto &ihistory: histories)
          dynamic_cast<FixBondHistory *>(ihistory)->delete_history(i,ibond);
      bond_atom[i][ibond] = tag[jnext];
    }
  for (jbond = 0; jbond < num_bond[j]; jbond++)
    if (bond_atom[j][jbond] == tag[jnext]) {
      if (n_histories > 0)
        for (auto &ihistory: histories)
          dynamic_cast<FixBondHistory *>(ihistory)->delete_history(j,jbond);
      bond_atom[j][jbond] = tag[inext];
    }
  for (ibond = 0; ibond < num_bond[inext]; ibond++)
    if (bond_atom[inext][ibond] == tag[i]) {
      if (n_histories > 0)
        for (auto &ihistory: histories)
          dynamic_cast<FixBondHistory *>(ihistory)->delete_history(inext,ibond);
      bond_atom[inext][ibond] = tag[j];
    }
  for (jbond = 0; jbond < num_bond[jnext]; jbond++)
    if (bond_atom[jnext][jbond] == tag[j]) {
      if (n_histories > 0)
        for (auto &ihistory: histories)
          dynamic_cast<FixBondHistory *>(ihistory)->delete_history(jnext,jbond);
      bond_atom[jnext][jbond] = tag[i];
    }

  // set global tags of 4 atoms in bonds

  itag = tag[i];
  inexttag = tag[inext];

  jtag = tag[j];
  jnexttag = tag[jnext];

  // change 1st special neighbors of affected atoms: i,j,inext,jnext
  // don't need to change 2nd/3rd special neighbors for any atom
  //   since special bonds = 0 1 1 means they are never used

  for (m = 0; m < nspecial[i][0]; m++)
    if (special[i][m] == inexttag) special[i][m] = jnexttag;
  for (m = 0; m < nspecial[j][0]; m++)
    if (special[j][m] == jnexttag) special[j][m] = inexttag;
  for (m = 0; m < nspecial[inext][0]; m++)
    if (special[inext][m] == itag) special[inext][m] = jtag;
  for (m = 0; m < nspecial[jnext][0]; m++)
    if (special[jnext][m] == jtag) special[jnext][m] = itag;

  // done if no angles

  if (!angleflag) return;

  // set global tags of 4 additional atoms in angles, 0 if no angle

  if (iprev >= 0) iprevtag = tag[iprev];
  else iprevtag = 0;
  if (ilast >= 0) ilasttag = tag[ilast];
  else ilasttag = 0;

  if (jprev >= 0) jprevtag = tag[jprev];
  else jprevtag = 0;
  if (jlast >= 0) jlasttag = tag[jlast];
  else jlasttag = 0;

  // change angle partners of affected atoms
  // must check if each angle is stored as a-b-c or c-b-a
  // on atom i:
  //   angle iprev-i-inext changes to iprev-i-jnext
  //   angle i-inext-ilast changes to i-jnext-jlast
  // on atom j:
  //   angle jprev-j-jnext changes to jprev-j-inext
  //   angle j-jnext-jlast changes to j-inext-ilast
  // on atom inext:
  //   angle iprev-i-inext changes to jprev-j-inext
  //   angle i-inext-ilast changes to j-inext-ilast
  // on atom jnext:
  //   angle jprev-j-jnext changes to iprev-i-jnext
  //   angle j-jnext-jlast changes to i-jnext-jlast

  for (iangle = 0; iangle < num_angle[i]; iangle++) {
    i1 = angle_atom1[i][iangle];
    i2 = angle_atom2[i][iangle];
    i3 = angle_atom3[i][iangle];

    if (i1 == iprevtag && i2 == itag && i3 == inexttag)
      angle_atom3[i][iangle] = jnexttag;
    else if (i1 == inexttag && i2 == itag && i3 == iprevtag)
      angle_atom1[i][iangle] = jnexttag;
    else if (i1 == itag && i2 == inexttag && i3 == ilasttag) {
      angle_atom2[i][iangle] = jnexttag;
      angle_atom3[i][iangle] = jlasttag;
    } else if (i1 == ilasttag && i2 == inexttag && i3 == itag) {
      angle_atom1[i][iangle] = jlasttag;
      angle_atom2[i][iangle] = jnexttag;
    }
  }

  for (jangle = 0; jangle < num_angle[j]; jangle++) {
    j1 = angle_atom1[j][jangle];
    j2 = angle_atom2[j][jangle];
    j3 = angle_atom3[j][jangle];

    if (j1 == jprevtag && j2 == jtag && j3 == jnexttag)
      angle_atom3[j][jangle] = inexttag;
    else if (j1 == jnexttag && j2 == jtag && j3 == jprevtag)
      angle_atom1[j][jangle] = inexttag;
    else if (j1 == jtag && j2 == jnexttag && j3 == jlasttag) {
      angle_atom2[j][jangle] = inexttag;
      angle_atom3[j][jangle] = ilasttag;
    } else if (j1 == jlasttag && j2 == jnexttag && j3 == jtag) {
      angle_atom1[j][jangle] = ilasttag;
      angle_atom2[j][jangle] = inexttag;
    }
  }

  for (iangle = 0; iangle < num_angle[inext]; iangle++) {
    i1 = angle_atom1[inext][iangle];
    i2 = angle_atom2[inext][iangle];
    i3 = angle_atom3[inext][iangle];

    if (i1 == iprevtag && i2 == itag && i3 == inexttag) {
      angle_atom1[inext][iangle] = jprevtag;
      angle_atom2[inext][iangle] = jtag;
    } else if (i1 == inexttag && i2 == itag && i3 == iprevtag) {
      angle_atom2[inext][iangle] = jtag;
      angle_atom3[inext][iangle] = jprevtag;
    } else if (i1 == itag && i2 == inexttag && i3 == ilasttag)
      angle_atom1[inext][iangle] = jtag;
    else if (i1 == ilasttag && i2 == inexttag && i3 == itag)
      angle_atom3[inext][iangle] = jtag;
  }

  for (jangle = 0; jangle < num_angle[jnext]; jangle++) {
    j1 = angle_atom1[jnext][jangle];
    j2 = angle_atom2[jnext][jangle];
    j3 = angle_atom3[jnext][jangle];

    if (j1 == jprevtag && j2 == jtag && j3 == jnexttag) {
      angle_atom1[jnext][jangle] = iprevtag;
      angle_atom2[jnext][jangle] = itag;
    } else if (j1 == jnexttag && j2 == jtag && j3 == jprevtag) {
      angle_atom2[jnext][jangle] = itag;
      angle_atom3[jnext][jangle] = iprevtag;
    } else if (j1 == jtag && j2 == jnexttag && j3 == jlasttag)
      angle_atom1[jnext][jangle] = itag;
    else if (j1 == jlasttag && j2 == jnexttag && j3 == jtag)
      angle_atom3[jnext][jangle] = itag;
  }

  // done if newton bond set

  if (newton_bond) return;

  // change angles stored by iprev,ilast,jprev,jlast
  // on atom iprev: angle iprev-i-inext changes to iprev-i-jnext
  // on atom jprev: angle jprev-j-jnext changes to jprev-j-inext
  // on atom ilast: angle i-inext-ilast changes to j-inext-ilast
  // on atom jlast: angle j-jnext-jlast changes to i-jnext-jlast

  for (iangle = 0; iangle < num_angle[iprev]; iangle++) {
    i1 = angle_atom1[iprev][iangle];
    i2 = angle_atom2[iprev][iangle];
    i3 = angle_atom3[iprev][iangle];

    if (i1 == iprevtag && i2 == itag && i3 == inexttag)
      angle_atom3[iprev][iangle] = jnexttag;
    else if (i1 == inexttag && i2 == itag && i3 == iprevtag)
      angle_atom1[iprev][iangle] = jnexttag;
  }

  for (jangle = 0; jangle < num_angle[jprev]; jangle++) {
    j1 = angle_atom1[jprev][jangle];
    j2 = angle_atom2[jprev][jangle];
    j3 = angle_atom3[jprev][jangle];

    if (j1 == jprevtag && j2 == jtag && j3 == jnexttag)
      angle_atom3[jprev][jangle] = inexttag;
    else if (j1 == jnexttag && j2 == jtag && j3 == jprevtag)
      angle_atom1[jprev][jangle] = inexttag;
  }

  for (iangle = 0; iangle < num_angle[ilast]; iangle++) {
    i1 = angle_atom1[ilast][iangle];
    i2 = angle_atom2[ilast][iangle];
    i3 = angle_atom3[ilast][iangle];

    if (i1 == itag && i2 == inexttag && i3 == ilasttag)
      angle_atom1[ilast][iangle] = jtag;
    else if (i1 == ilasttag && i2 == inexttag && i3 == itag)
      angle_atom3[ilast][iangle] = jtag;
  }

  for (jangle = 0; jangle < num_angle[jlast]; jangle++) {
    j1 = angle_atom1[jlast][jangle];
    j2 = angle_atom2[jlast][jangle];
    j3 = angle_atom3[jlast][jangle];

    if (j1 == jtag && j2 == jnexttag && j3 == jlasttag)
      angle_atom1[jlast][jangle] = itag;
    else if (j1 == jlasttag && j2 == jnexttag && j3 == jtag)
      angle_atom3[jlast][jangle] = itag;
  }
}

/* ---------------------------------------------------------------------- */

int FixBondSwap::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete[] id_temp;
    id_temp = utils::strdup(arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,"Fix_modify temperature ID does not "
                 "compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   compute squared distance between atoms I,J
   must use minimum_image since J was found thru atom->map()
------------------------------------------------------------------------- */

double FixBondSwap::dist_rsq(int i, int j)
{
  double delx = x[i][0] - x[j][0];
  double dely = x[i][1] - x[j][1];
  double delz = x[i][2] - x[j][2];
  domain->minimum_image(delx,dely,delz);
  return (delx*delx + dely*dely + delz*delz);
}

/* ----------------------------------------------------------------------
   return pairwise interaction energy between atoms I,J
   will always be full non-bond interaction, so factors = 1 in single() call
------------------------------------------------------------------------- */

double FixBondSwap::pair_eng(int i, int j)
{
  double tmp;
  double rsq = dist_rsq(i,j);
  return force->pair->single(i,j,type[i],type[j],rsq,1.0,1.0,tmp);
}

/* ---------------------------------------------------------------------- */

double FixBondSwap::bond_eng(int btype, int i, int j)
{
  double tmp;
  double rsq = dist_rsq(i,j);
  return force->bond->single(btype,rsq,i,j,tmp);
}

/* ---------------------------------------------------------------------- */

double FixBondSwap::angle_eng(int atype, int i, int j, int k)
{
  // test for non-existent angle at end of chain

  if (i == -1 || k == -1) return 0.0;
  return force->angle->single(atype,i,j,k);
}

/* ----------------------------------------------------------------------
   create a random permutation of one atom's N neighbor list atoms
   uses one-pass Fisher-Yates shuffle on an initial identity permutation
   output: randomized permute[] vector, used to index neighbors
------------------------------------------------------------------------- */

void FixBondSwap::neighbor_permutation(int n)
{
  int i,j,tmp;

  if (n > maxpermute) {
    delete[] permute;
    maxpermute = n + DELTA_PERMUTE;
    permute = new int[maxpermute];
  }

  // Fisher-Yates shuffle

  for (i = 0; i < n; i++) permute[i] = i;

  for (i = 0; i < n-1; i++) {
    j = i + static_cast<int> (random->uniform() * (n-i));
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

double FixBondSwap::compute_vector(int n)
{
  double one,all;
  if (n == 0) one = naccept;
  else one = foursome;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

/* ----------------------------------------------------------------------
   memory usage of alist
------------------------------------------------------------------------- */

double FixBondSwap::memory_usage()
{
  double bytes = (double)nmax * sizeof(int);
  return bytes;
}
