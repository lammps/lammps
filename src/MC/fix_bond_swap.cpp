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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_bond_swap.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "group.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "compute.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBondSwap::FixBondSwap(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix bond/swap command");

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;

  fraction = atof(arg[3]);
  double cutoff = atof(arg[4]);
  cutsq = cutoff*cutoff;

  // initialize Marsaglia RNG with processor-unique seed

  int seed = atoi(arg[5]);
  random = new RanMars(lmp,seed + comm->me);

  // create a new compute temp style
  // id = fix-ID + temp, compute group = fix group

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "temp";
  modify->add_compute(3,newarg);
  delete [] newarg;
  tflag = 1;

  // initialize atom list

  nmax = 0;
  alist = NULL;

  naccept = foursome = 0;
}

/* ---------------------------------------------------------------------- */

FixBondSwap::~FixBondSwap()
{
  delete random;

  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete [] id_temp;

  memory->destroy(alist);
}

/* ---------------------------------------------------------------------- */

int FixBondSwap::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondSwap::init()
{
  // require an atom style with molecule IDs

  if (atom->molecule == NULL)
    error->all(FLERR,"Must use atom style with molecule IDs with fix bond/swap");

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0) 
    error->all(FLERR,"Temperature ID for fix bond/swap does not exist");
  temperature = modify->compute[icompute];

  // pair and bonds must be defined
  // no dihedral or improper potentials allowed
  // special bonds must be 0 1 1

  if (force->pair == NULL || force->bond == NULL)
    error->all(FLERR,"Fix bond/swap requires pair and bond styles");

  if (force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support fix bond/swap");

  if (force->angle == NULL && atom->nangles > 0 && comm->me == 0)
    error->warning(FLERR,"Fix bond/swap will ignore defined angles");

  if (force->dihedral || force->improper)
    error->all(FLERR,"Fix bond/swap cannot use dihedral or improper styles");

  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 ||
      force->special_lj[3] != 1.0)
    error->all(FLERR,"Fix bond/swap requires special_bonds = 0,1,1");

  // need a half neighbor list, built when ever re-neighboring occurs

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;

  // zero out stats

  naccept = foursome = 0;
  angleflag = 0;
  if (force->angle) angleflag = 1;
}

/* ---------------------------------------------------------------------- */

void FixBondSwap::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixBondSwap::pre_neighbor()
{
  int i,j,ii,jj,m,inum,jnum;
  int inext,iprev,ilast,jnext,jprev,jlast,ibond,iangle,jbond,jangle;
  int itag,inexttag,iprevtag,ilasttag,jtag,jnexttag,jprevtag,jlasttag;
  int ibondtype,jbondtype,iangletype,inextangletype,jangletype,jnextangletype;
  int i1,i2,i3,j1,j2,j3,tmp;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double delta,factor;

  // compute current temp for Boltzmann factor test

  double t_current = temperature->compute_scalar();

  // local ptrs to atom arrays

  int *tag = atom->tag;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *num_angle = atom->num_angle;
  int **angle_atom1 = atom->angle_atom1;
  int **angle_atom2 = atom->angle_atom2;
  int **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  int **nspecial = atom->nspecial;
  int **special = atom->special;
  int newton_bond = force->newton_bond;
  int nlocal = atom->nlocal;

  type = atom->type;
  x = atom->x;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // randomize list of my owned atoms that are in fix group
  // grow atom list if necessary

  if (nlocal > nmax) {
    memory->destroy(alist);
    nmax = atom->nmax;
    memory->create(alist,nmax,"bondswap:alist");
  }

  int neligible = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit)
      alist[neligible++] = i;
  }

  for (i = 0; i < neligible; i++) {
    j = static_cast<int> (random->uniform() * neligible);
    tmp = alist[i];
    alist[i] = alist[j];
    alist[j] = tmp;
  }

  // examine ntest of my eligible atoms for potential swaps
  // atom i is randomly selected via atom list
  // look at all j neighbors of atom i
  // atom j must be on-processor (j < nlocal)
  // atom j must be in fix group
  // i and j must be same distance from chain end (mol[i] = mol[j])
  // NOTE: must use extra parens in if test on mask[j] & groupbit

  int ntest = static_cast<int> (fraction * neligible);
  int accept = 0;

  for (int itest = 0; itest < ntest; itest++) {
    i = alist[itest];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      if (j >= nlocal) continue;
      if ((mask[j] & groupbit) == 0) continue;
      if (molecule[i] != molecule[j]) continue;

      // look at all bond partners of atoms i and j
      // use num_bond for this, not special list, so also find bondtypes
      // inext,jnext = bonded atoms
      // inext,jnext must be on-processor (inext,jnext < nlocal)
      // inext,jnext must be same dist from chain end (mol[inext] = mol[jnext])
      // since swaps may occur between two ends of a single chain, insure
      //   the 4 atoms are unique (no duplicates): inext != jnext, inext != j
      // all 4 old and new bonds must have length < cutoff

      for (ibond = 0; ibond < num_bond[i]; ibond++) {
	inext = atom->map(bond_atom[i][ibond]);
	if (inext >= nlocal || inext < 0) continue;
	ibondtype = bond_type[i][ibond];

	for (jbond = 0; jbond < num_bond[j]; jbond++) {
	  jnext = atom->map(bond_atom[j][jbond]);
	  if (jnext >= nlocal || jnext < 0) continue;
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
	  // use num_angle for this, not special list, so also find angletypes
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
  if (!accept) return;
  naccept++;
  
  // change bond partners of affected atoms
  // on atom i: bond i-inext changes to i-jnext
  // on atom j: bond j-jnext changes to j-inext
  // on atom inext: bond inext-i changes to inext-j
  // on atom jnext: bond jnext-j changes to jnext-i

  for (ibond = 0; ibond < num_bond[i]; ibond++)
    if (bond_atom[i][ibond] == tag[inext]) bond_atom[i][ibond] = tag[jnext];
  for (jbond = 0; jbond < num_bond[j]; jbond++)
    if (bond_atom[j][jbond] == tag[jnext]) bond_atom[j][jbond] = tag[inext];
  for (ibond = 0; ibond < num_bond[inext]; ibond++)
    if (bond_atom[inext][ibond] == tag[i]) bond_atom[inext][ibond] = tag[j];
  for (jbond = 0; jbond < num_bond[jnext]; jbond++)
    if (bond_atom[jnext][jbond] == tag[j]) bond_atom[jnext][jbond] = tag[i];

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
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,"Fix_modify temperature ID does not compute temperature");
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
  double rsq = dist_rsq(i,j);
  return force->bond->single(btype,rsq,i,j);
}

/* ---------------------------------------------------------------------- */

double FixBondSwap::angle_eng(int atype, int i, int j, int k)
{
  // test for non-existent angle at end of chain

  if (i == -1 || k == -1) return 0.0;
  return force->angle->single(atype,i,j,k);
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
  double bytes = nmax * sizeof(int);
  return bytes;
}
