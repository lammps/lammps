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

#include "fix_bond_mcmove.h"

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

#define DELTA_PERMUTE 100

/* ---------------------------------------------------------------------- */

FixBondMcMove::FixBondMcMove(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    tflag(0), alist(nullptr), id_temp(nullptr), type(nullptr), x(nullptr), list(nullptr),
    temperature(nullptr), random(nullptr)
{
    if (narg != 7) error->all(FLERR,"Illegal fix bond/mcmove command");

    trackbondtype = utils::inumeric(FLERR,arg[3],false,lmp);
    if (trackbondtype <= 0) error->all(FLERR,"Illegal fix bond/mcmove command");

    movingbondtype = utils::inumeric(FLERR,arg[4],false,lmp);
    if (movingbondtype <= 0) error->all(FLERR,"Illegal fix bond/mcmove command");

    nevery = utils::inumeric(FLERR,arg[5],false,lmp);
    if (nevery <= 0) error->all(FLERR,"Illegal fix bond/mcmove command");

    force_reneighbor = 1;
    next_reneighbor = -1;
    vector_flag = 1;
    size_vector = 2;
    global_freq = 1;
    extvector = 0;

    fraction = utils::numeric(FLERR,arg[6],false,lmp);

    // initialize Marsaglia RNG with processor-unique seed

    int seed = utils::inumeric(FLERR,arg[7],false,lmp);
    random = new RanMars(lmp,seed + comm->me);

    // error check

    if (atom->molecular != Atom::MOLECULAR)
        error->all(FLERR,"Cannot use fix bond/mcmove with non-molecular systems");

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

    naccept = threesome = 0;
}

/* ---------------------------------------------------------------------- */

FixBondMcMove::~FixBondMcMove()
{
    delete random;

    // delete temperature if fix created it

    if (tflag) modify->delete_compute(id_temp);
    delete[] id_temp;

    memory->destroy(alist);
    delete[] permute;
}

/* ---------------------------------------------------------------------- */

int FixBondMcMove::setmask()
{
    int mask = 0;
    mask |= POST_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondMcMove::init()
{
    // require an atom style with molecule IDs

    if (atom->molecule == nullptr)
        error->all(FLERR,
                   "Must use atom style with molecule IDs with fix bond/mcmove");

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
        error->all(FLERR,"Temperature ID for fix bond/mcmove does not exist");
    temperature = modify->compute[icompute];

    // pair and bonds must be defined
    // no dihedral or improper potentials allowed
    // special bonds must be 0 1 1

    if (force->pair == nullptr || force->bond == nullptr)
        error->all(FLERR,"Fix bond/mcmove requires pair and bond styles");

    if (force->pair->single_enable == 0)
        error->all(FLERR,"Pair style does not support fix bond/mcmove");

    if (force->angle == nullptr && atom->nangles > 0 && comm->me == 0)
        error->warning(FLERR,"Fix bond/mcmove will not preserve correct angle "
                       "topology because no angle_style is defined");

    if (force->dihedral || force->improper)
        error->all(FLERR,"Fix bond/mcmove cannot use dihedral or improper styles");

    if (force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0)
        error->all(FLERR,"Fix bond/mcmove requires special_bonds = {any},1,1");

    // need a half neighbor list, built every Nevery steps

    neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);

    // zero out stats

    naccept = threesome = 0;
    angleflag = 0;
    if (force->angle) angleflag = 1;
}

/* ---------------------------------------------------------------------- */

void FixBondMcMove::init_list(int /*id*/, NeighList *ptr)
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

void FixBondMcMove::post_integrate()
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
        bool do_swap = false;
        int track_bond_i = -1;
        int move_bond_i = -1;

        // loop over all bonds of atom I
        int num_track_bonds = 0;
        int num_move_bonds = 0;
        for (int bond_i = 0; bond_i < num_bond[i]; ++bond_i) {
            if (bond_type[i][bond_i] == trackbondtype) {
                num_track_bonds++;
                track_bond_i = bond_i;
            }
            if (bond_type[i][bond_i] == movingbondtype) {
                num_move_bonds++;
                move_bond_i = bond_i;
            }
        }
        if (num_move_bonds < 1 || num_track_bonds < 1) {
            continue;
        }

        // default bonds for track & move of this atom is the last one
        // now we randomly select another one
        if (num_track_bonds > 1) {
            for (int bond_i = 0; bond_i < num_bond[i]; ++bond_i) {
                if (bond_type[i][bond_i] == trackbondtype) {
                    if (this->random->uniform()*static_cast<double>(num_track_bonds) < 1.){
                        track_bond_i = bond_i; // choose this randomly, not just the last one
                        break;
                    }
                }
            }
        }
        if (num_move_bonds > 1){
            for (int bond_i = 0; bond_i < num_bond[i]; ++bond_i) {
                if (bond_type[i][bond_i] == movingbondtype) {
                    if (this->random->uniform()*static_cast<double>(num_move_bonds) < 1.){
                        move_bond_i = bond_i; // choose this randomly, not just the last one
                        break;
                    }
                }
            }
        }

        // want to move the other end of the bond move_bond_i to
        // the other end of the track_bond_i
        // but accept only under conditions.
        j = atom->map(bond_atom[i][move_bond_i]);
        inext = atom->map(bond_atom[i][track_bond_i]);
        // now, it's: swap bond j-i to j-inext

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

        // valid threesome found between 2 chains:
        //   chains = iprev-i-inext-ilast and jprev-j-jnext-jlast
        //   prev or last values are -1 if do not exist due to end of chain
        //   OK to call angle_eng with -1 atom, since just return 0.0
        // current energy of threesome =
        //   E_nb(i,j) + E_nb(i,jnext) + E_nb(inext,j) + E_nb(inext,jnext) +
        //   E_bond(i,inext) + E_bond(j,jnext) +
        //   E_angle(iprev,i,inext) + E_angle(i,inext,ilast) +
        //   E_angle(jprev,j,jnext) + E_angle(j,jnext,jlast)
        // new energy of threesome with swapped bonds =
        //   E_nb(i,j) + E_nb(i,inext) + E_nb(j,jnext) + E_nb(inext,jnext) +
        //   E_bond(i,jnext) + E_bond(j,inext) +
        //   E_angle(iprev,i,jnext) + E_angle(i,jnext,jlast) +
        //   E_angle(jprev,j,inext) + E_angle(j,inext,ilast)
        // energy delta = add/subtract differing terms between 2 formulas

        threesome++;

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

        if (delta < 0.0) {
            accept += 1;
            do_swap = true;
        }
        else {
            factor = exp(-delta/force->boltz/t_current);
            if (random->uniform() < factor) {
                accept += 1;
                do_swap = true;
            }
        }

        // actually invoke swapping procedure.
        if (do_swap) {
          move_bond(j, jprev, jnext, i, iprev, inext, ilast);
        }
        // goto done;
    }

done:

    // trigger immediate reneighboring if swaps occurred on one or more procs

    int accept_any;
    MPI_Allreduce(&accept,&accept_any,1,MPI_INT,MPI_SUM,world);
    if (accept_any) next_reneighbor = update->ntimestep;

    if (!accept) return;
    naccept+= accept;
}

// move bond from
void FixBondMcMove::move_bond(int rotation_i, int rot_prev, int rot_next, int prev_partner, int preprev_partner, int next_partner, int postnext_partner) {
    int m, inum;
    int iangle,jangle,ibond,jbond;
    tagint i1,i2,i3,j1,j2,j3;
    int *ilist,*jlist,*numneigh,**firstneigh;

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

// find instances of bond/history to reset history
    auto histories = modify->get_fix_by_style("BOND_HISTORY");
    int n_histories = histories.size();

    // change bond partners of affected atoms
    // on atom rot: bond rot-prev_partner changes to rot-next_partner

    for (ibond = 0; ibond < num_bond[rotation_i]; ++ibond) {
        if (bond_atom[rotation_i][ibond] == tag[prev_partner]) {
            if (n_histories > 0) {
                for (auto &ihistory: histories) {
                    dynamic_cast<FixBondHistory *>(ihistory)->delete_history(rotation_i,prev_partner);
                }
            }
            bond_atom[rotation_i][ibond] = tag[next_partner];
        }
    }

    // on atom next_partner: new bond rot-next_partner
    if (atom->bond_per_atom < num_bond[next_partner]+1) {
        error->all(FLERR, "No more space to create bond. Use extra/bond/per/atom to increase.");
    }
    bond_atom[next_partner][num_bond[next_partner]] = tag[rotation_i];
    num_bond[next_partner]++;

    // on atom prev_partner: bond rot-prev_partner is removed
    for (ibond = 0; ibond < num_bond[prev_partner]; ++ibond) {
        if (bond_atom[prev_partner][ibond] == tag[rotation_i]) {
            for (int k = ibond; k < num_bond[prev_partner]-1; ++k) {

                bond_atom[prev_partner][k] = bond_atom[prev_partner][k+1];
                bond_type[prev_partner][k] = bond_type[prev_partner][k+1];
                if (n_histories > 0)
                    for (auto &ihistory: histories)
                        dynamic_cast<FixBondHistory *>(ihistory)->shift_history(prev_partner,k,k+1);
            }
            if (n_histories > 0)
                for (auto &ihistory: histories)
                    dynamic_cast<FixBondHistory *>(ihistory)->delete_history(prev_partner,num_bond[prev_partner]-1);
            num_bond[prev_partner]--;
            break;
        }
    }

    // set global tags of 3 atoms in bonds

    tagint rottag = tag[rotation_i];
    tagint prevtag = tag[prev_partner];
    tagint nexttag = tag[next_partner];

    // change 1st special neighbors of affected atoms: i,j,inext,jnext
    // don't need to change 2nd/3rd special neighbors for any atom
    //   since special bonds = 0 1 1 means they are never used
    // TODO: check / verify somehow
    if (force->special_lj[0] != 1.0){
    // change new 1st special neighbor on rotation atom
    for (m = 0; m < nspecial[rotation_i][0]; m++)
        if (special[rotation_i][m] == prevtag) special[rotation_i][m] = nexttag;
    
    // add the new primary partner to the new bond partner
    for (m = (nspecial[next_partner][0]+nspecial[next_partner][1]+nspecial[next_partner][2]); m > nspecial[next_partner][0];--m) {
      special[next_partner][m] = special[next_partner][m-1];
    }
    special[next_partner][nspecial[next_partner][0]] = rottag;
    nspecial[next_partner][0]++;

    // remove the old primary partner
    for (m = 0; m < nspecial[prev_partner][0]; m++)
        if (special[prev_partner][m] == rottag) {
            for (int k = m; k < nspecial[prev_partner][0]-1; ++k) {
                special[prev_partner][k] = special[prev_partner][k+1];
            }
            nspecial[prev_partner][0]--;
            break;
        }
}
    // done if no angles

    if (!angleflag) return;

    // set global tags of 4 additional atoms in angles, 0 if no angle
    int preprev_tag, postnext_tag, rotprev_tag, rotnext_tag;
    if (preprev_partner >= 0) preprev_tag = tag[preprev_partner];
    else preprev_tag = 0;
    if (postnext_partner >= 0) postnext_tag = tag[postnext_partner];
    else postnext_tag = 0;

    if (rot_prev >= 0) rotprev_tag = tag[rot_prev];
    else rotprev_tag = 0;
    if (rot_next >= 0) rotnext_tag = tag[rot_next];
    else rotnext_tag = 0;

    // change angle partners of affected atoms
    // must check if each angle is stored as a-b-c or c-b-a
    // TODO: verify/check somehow

    // on atom rot:
    //    angle rot_prev-rot-prev_partner changes to rot_prev-rot-next_partner
    //    angle rot_next-rot-prev_partner changes to rot_next-rot-next_partner

    for (iangle = 0; iangle < num_angle[rotation_i]; iangle++) {
        i1 = angle_atom1[rotation_i][iangle];
        i2 = angle_atom2[rotation_i][iangle];
        i3 = angle_atom3[rotation_i][iangle];

        if (i1 == rotprev_tag && i2 == rottag && i3 == prevtag)
            angle_atom3[rotation_i][iangle] = nexttag;
        else if (i1 == prevtag && i2 == rottag && i3 == rotprev_tag)
            angle_atom1[rotation_i][iangle] = nexttag;
        else if (i1 == rotnext_tag && i2 == rottag && i3 == prevtag)
            angle_atom3[rotation_i][iangle] = nexttag;
        else if (i1 == prevtag && i2 == rottag && i3 == rotnext_tag)
            angle_atom1[rotation_i][iangle] = nexttag; 
    }

    // on atom prev_partner:
    //    angle rot-prev_partner-next_partner:
    //        remove if not newton
    //        otherwise, change to rot-next_partner-prev_partner
    //    angle pre_prev-prev_partner-rot: remove
    //    angle next_partner-prev_partner-rot: remove
    // TODO: if there is another bond (type?) leading to these angles, are we destroying them?
    for (iangle = num_angle[prev_partner]-1; iangle >= 0; iangle--) {
        i1 = angle_atom1[prev_partner][iangle];
        i2 = angle_atom2[prev_partner][iangle];
        i3 = angle_atom3[prev_partner][iangle];
        bool remove_iangle = false;

        if (i1 == rottag && i2 == prevtag && i3 == nexttag) {
          remove_iangle = true;
        } else if (i1 == nexttag && i2 == prevtag && i3 == rottag) {
          remove_iangle = true;
        } else if (i1 == preprev_tag && i2 == prevtag && i3 == rottag) {
          remove_iangle = true;
        }else if (i1 == rottag && i2 == prevtag && i3 == preprev_tag) {
          remove_iangle = true;
        }

        if (remove_iangle) {
          for (m = iangle; m < num_angle[prev_partner]-1; ++m) {
            angle_atom1[prev_partner][m] = angle_atom1[prev_partner][m+1];
            angle_atom2[prev_partner][m] = angle_atom2[prev_partner][m+1];
            angle_atom3[prev_partner][m] = angle_atom3[prev_partner][m+1];
          }
          num_angle[prev_partner]--;
        }
      }

      // on atom next_partner:
      //    angle rot-next_partner-postnext: add
      // alternative: add rot-next_partner-prev_partner
      int n_next_partner_angles = num_angle[next_partner];
      angle_atom1[next_partner][n_next_partner_angles] = rottag;
      angle_atom2[next_partner][n_next_partner_angles] = nexttag;
      angle_atom3[next_partner][n_next_partner_angles] = postnext_tag;
      num_angle[next_partner]++;
      // TODO: do we need to add an angle type or something here?

    // done if newton bond set

    if (newton_bond) return;

    // change angles stored by other atoms
    // TODO: check/verfiy somehow, and correct
    // on atom rotprev: change potential rotprev-rot-prev_partner to rotprev-rot-next_partner
    // break out with the objective to keep the possibly destroyed angles small.
    // could happen e.g. if there is another bond leading to these angles
    for (iangle = 0; iangle < num_angle[rot_prev]; ++iangle) {
        i1 = angle_atom1[rot_prev][iangle];
        i2 = angle_atom2[rot_prev][iangle];
        i3 = angle_atom3[rot_prev][iangle];

        if (i1 == rotprev_tag && i2 == rottag && i3 == prevtag) {
          angle_atom3[rot_prev][iangle] = nexttag; break;
        } else if (i1 == prevtag && i2 == rottag && i3 == rotprev_tag) {
          angle_atom1[rot_prev][iangle] = nexttag; break;
        }
    }

    // on atom rotnext: change potential rotnext-rot-prev_partner to rotnext-rot-next_partner
    for (iangle = 0; iangle < num_angle[rot_next]; ++iangle) {
        i1 = angle_atom1[rot_next][iangle];
        i2 = angle_atom2[rot_next][iangle];
        i3 = angle_atom3[rot_next][iangle];

        if (i1 == rotnext_tag && i2 == rottag && i3 == prevtag) {
          angle_atom3[rot_next][iangle] = nexttag; break;
        } else if (i1 == prevtag && i2 == rottag && i3 == rotnext_tag) {
          angle_atom1[rot_next][iangle] = nexttag; break;
        }
    }

    // TODO: continue with other atoms
}

/* ---------------------------------------------------------------------- */

int FixBondMcMove::modify_param(int narg, char **arg)
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

double FixBondMcMove::dist_rsq(int i, int j)
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

double FixBondMcMove::pair_eng(int i, int j)
{
    double tmp;
    double rsq = dist_rsq(i,j);
    return force->pair->single(i,j,type[i],type[j],rsq,1.0,1.0,tmp);
}

/* ---------------------------------------------------------------------- */

double FixBondMcMove::bond_eng(int btype, int i, int j)
{
    double tmp;
    double rsq = dist_rsq(i,j);
    return force->bond->single(btype,rsq,i,j,tmp);
}

/* ---------------------------------------------------------------------- */

double FixBondMcMove::angle_eng(int atype, int i, int j, int k)
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

void FixBondMcMove::neighbor_permutation(int n)
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

double FixBondMcMove::compute_vector(int n)
{
    double one,all;
    if (n == 0) one = naccept;
    else one = threesome;
    MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
    return all;
}

/* ----------------------------------------------------------------------
   memory usage of alist
------------------------------------------------------------------------- */

double FixBondMcMove::memory_usage()
{
    double bytes = (double)nmax * sizeof(int);
    return bytes;
}
