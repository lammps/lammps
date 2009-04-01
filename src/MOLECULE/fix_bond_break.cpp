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
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_bond_break.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixBondBreak::FixBondBreak(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all("Illegal fix bond/break command");

  MPI_Comm_rank(world,&me);

  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix bond/break command");

  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  size_vector = 2;
  scalar_vector_freq = 1;
  extvector = 1;

  btype = atoi(arg[4]);
  double cutoff = atof(arg[5]);

  if (btype < 1 || btype > atom->nbondtypes)
    error->all("Invalid bond type in fix bond/break command");
  if (cutoff < 0.0) error->all("Illegal fix bond/break command");

  cutsq = cutoff*cutoff;

  // optional keywords

  fraction = 1.0;
  int seed = 12345;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"prob") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix bond/break command");
      fraction = atof(arg[iarg+1]);
      seed = atoi(arg[iarg+2]);
      if (fraction < 0.0 || fraction > 1.0)
	error->all("Illegal fix bond/break command");
      if (seed <= 0) error->all("Illegal fix bond/break command");
      iarg += 3;
    } else error->all("Illegal fix bond/break command");
  }

  // error check

  if (atom->molecular == 0)
    error->all("Cannot use fix bond/break with non-molecular systems");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + me);

  // set comm sizes needed by this fix

  comm_forward = 2;
  comm_reverse = 2;

  // allocate arrays local to this fix

  nmax = 0;
  partner = NULL;
  distsq = NULL;

  // zero out stats

  breakcount = 0;
  breakcounttotal = 0;
}

/* ---------------------------------------------------------------------- */

FixBondBreak::~FixBondBreak()
{
  delete random;

  // delete locally stored arrays

  memory->sfree(partner);
  memory->sfree(distsq);
}

/* ---------------------------------------------------------------------- */

int FixBondBreak::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondBreak::init()
{
  // require special bonds = 0,1,1

  int flag = 0;
  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 ||
      force->special_lj[3] != 1.0) flag = 1;
  if (force->special_coul[1] != 0.0 || force->special_coul[2] != 1.0 ||
      force->special_coul[3] != 1.0) flag = 1;
  if (flag) error->all("Fix bond/break requires special_bonds = 0,1,1");

  // warn if angles, dihedrals, impropers are being used

  if (force->angle || force->dihedral || force->improper) {
    if (me == 0) 
      error->warning("Broken bonds will not alter angles, "
		     "dihedrals, or impropers");
  }

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixBondBreak::post_integrate()
{
  int i,j,k,m,n,i1,i2,n1,n3,possible,type;
  double delx,dely,delz,rsq,min,max;
  int *slist;

  if (update->ntimestep % nevery) return;

  // need updated ghost atom positions

  comm->communicate();

  // resize bond partner list and initialize it
  // probability array overlays distsq array
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->sfree(partner);
    memory->sfree(distsq);
    nmax = atom->nmax;
    partner = (int *)
      memory->smalloc(nmax*sizeof(int),"bond/break:partner");
    distsq = (double *)
      memory->smalloc(nmax*sizeof(double),"bond/break:distsq");
    probability = distsq;
  }

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    partner[i] = 0;
    distsq[i] = 0.0;
  }

  // loop over bond list
  // setup possible partner list of bonds to break

  double **x = atom->x;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    if (!(mask[i1] & groupbit)) continue;
    if (!(mask[i2] & groupbit)) continue;
    if (type != btype) continue;

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    domain->minimum_image(delx,dely,delz);
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq <= cutsq) continue;

    if (rsq > distsq[i1]) {
      partner[i1] = tag[i2];
      distsq[i1] = rsq;
    }
    if (rsq > distsq[i2]) {
      partner[i2] = tag[i1];
      distsq[i2] = rsq;
    }
  }

  // reverse comm of partner info

  if (force->newton_bond) comm->reverse_comm_fix(this);

  // each atom now knows its winning partner
  // for prob check, generate random value for each atom with a bond partner
  // forward comm of partner and random value, so ghosts have it

  if (fraction < 1.0) {
    for (i = 0; i < nlocal; i++)
      if (partner[i]) probability[i] = random->uniform();
  }

  comm->comm_fix(this);

  // break bonds
  // if both atoms list each other as winning bond partner
  // and probability constraint is satisfied

  int **bond_type = atom->bond_type;
  int **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  int **nspecial = atom->nspecial;
  int **special = atom->special;

  int nbreak = 0;
  for (i = 0; i < nlocal; i++) {
    if (partner[i] == 0) continue;
    j = atom->map(partner[i]);
    if (partner[j] != tag[i]) continue;

    // apply probability constraint
    // MIN,MAX insures values are added in same order on different procs

    if (fraction < 1.0) {
      min = MIN(probability[i],probability[j]);
      max = MAX(probability[i],probability[j]);
      if (0.5*(min+max) >= fraction) continue;
    }

    // delete bond from atom I if I stores it
    // atom J will also do this

    for (m = 0; m < num_bond[i]; m++) {
      if (bond_atom[i][m] == partner[i]) {
	for (k = m; k < num_bond[i]-1; k++) {
	  bond_atom[i][k] = bond_atom[i][k+1];
	  bond_type[i][k] = bond_type[i][k+1];
	}
	num_bond[i]--;
	break;
      }
    }

    // remove J from special bond list for atom I
    // atom J will also do this

    slist = atom->special[i];
    n1 = nspecial[i][0];
    n3 = nspecial[i][2];
    for (m = 0; m < n1; m++)
      if (slist[m] == partner[i]) break;
    for (; m < n3-1; m++) slist[m] = slist[m+1];
    nspecial[i][0]--;
    nspecial[i][1]--;
    nspecial[i][2]--;
    
    // count the broken bond once

    if (tag[i] < tag[j]) nbreak++;
  }

  // tally stats

  MPI_Allreduce(&nbreak,&breakcount,1,MPI_INT,MPI_SUM,world);
  breakcounttotal += breakcount;
  atom->nbonds -= breakcount;

  // trigger reneighboring if any bonds were formed

  if (breakcount) next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixBondBreak::post_integrate_respa(int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_integrate();
}

/* ---------------------------------------------------------------------- */

int FixBondBreak::pack_comm(int n, int *list, double *buf,
			     int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = partner[j];
    buf[m++] = probability[j];
  }
  return 2;
}

/* ---------------------------------------------------------------------- */

void FixBondBreak::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    partner[i] = static_cast<int> (buf[m++]);
    probability[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int FixBondBreak::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = partner[i];
    buf[m++] = distsq[i];
  }
  return 2;
}

/* ---------------------------------------------------------------------- */

void FixBondBreak::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (buf[m+1] > distsq[j]) {
      partner[j] = static_cast<int> (buf[m++]);
      distsq[j] = buf[m++];
    } else m += 2;
  }
}

/* ---------------------------------------------------------------------- */

double FixBondBreak::compute_vector(int n)
{
  if (n == 1) return (double) breakcount;
  return (double) breakcounttotal;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double FixBondBreak::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax * sizeof(double);
  return bytes;
}
