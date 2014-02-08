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
#include "fix_bond_create.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixBondCreate::FixBondCreate(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8) error->all(FLERR,"Illegal fix bond/create command");

  MPI_Comm_rank(world,&me);

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix bond/create command");

  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;

  iatomtype = force->inumeric(FLERR,arg[4]);
  jatomtype = force->inumeric(FLERR,arg[5]);
  double cutoff = force->numeric(FLERR,arg[6]);
  btype = force->inumeric(FLERR,arg[7]);

  if (iatomtype < 1 || iatomtype > atom->ntypes ||
      jatomtype < 1 || jatomtype > atom->ntypes)
    error->all(FLERR,"Invalid atom type in fix bond/create command");
  if (cutoff < 0.0) error->all(FLERR,"Illegal fix bond/create command");
  if (btype < 1 || btype > atom->nbondtypes)
    error->all(FLERR,"Invalid bond type in fix bond/create command");

  cutsq = cutoff*cutoff;

  // optional keywords

  imaxbond = 0;
  inewtype = iatomtype;
  jmaxbond = 0;
  jnewtype = jatomtype;
  fraction = 1.0;
  int seed = 12345;

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"iparam") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix bond/create command");
      imaxbond = force->inumeric(FLERR,arg[iarg+1]);
      inewtype = force->inumeric(FLERR,arg[iarg+2]);
      if (imaxbond < 0) error->all(FLERR,"Illegal fix bond/create command");
      if (inewtype < 1 || inewtype > atom->ntypes)
        error->all(FLERR,"Invalid atom type in fix bond/create command");
      iarg += 3;
    } else if (strcmp(arg[iarg],"jparam") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix bond/create command");
      jmaxbond = force->inumeric(FLERR,arg[iarg+1]);
      jnewtype = force->inumeric(FLERR,arg[iarg+2]);
      if (jmaxbond < 0) error->all(FLERR,"Illegal fix bond/create command");
      if (jnewtype < 1 || jnewtype > atom->ntypes)
        error->all(FLERR,"Invalid atom type in fix bond/create command");
      iarg += 3;
    } else if (strcmp(arg[iarg],"prob") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix bond/create command");
      fraction = force->numeric(FLERR,arg[iarg+1]);
      seed = force->inumeric(FLERR,arg[iarg+2]);
      if (fraction < 0.0 || fraction > 1.0)
        error->all(FLERR,"Illegal fix bond/create command");
      if (seed <= 0) error->all(FLERR,"Illegal fix bond/create command");
      iarg += 3;
    } else error->all(FLERR,"Illegal fix bond/create command");
  }

  // error check

  if (atom->molecular != 1)
    error->all(FLERR,"Cannot use fix bond/create with non-molecular systems");
  if (iatomtype == jatomtype &&
      ((imaxbond != jmaxbond) || (inewtype != jnewtype)))
    error->all(FLERR,
               "Inconsistent iparam/jparam values in fix bond/create command");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + me);

  // perform initial allocation of atom-based arrays
  // register with Atom class
  // bondcount values will be initialized in setup()

  bondcount = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  countflag = 0;

  // set comm sizes needed by this fix

  comm_forward = 2;
  comm_reverse = 2;

  // allocate arrays local to this fix

  nmax = 0;
  partner = NULL;
  distsq = NULL;

  // zero out stats

  createcount = 0;
  createcounttotal = 0;
}

/* ---------------------------------------------------------------------- */

FixBondCreate::~FixBondCreate()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  delete random;

  // delete locally stored arrays

  memory->destroy(bondcount);
  memory->destroy(partner);
  memory->destroy(distsq);
}

/* ---------------------------------------------------------------------- */

int FixBondCreate::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::init()
{
  // check cutoff for iatomtype,jatomtype

  if (force->pair == NULL || cutsq > force->pair->cutsq[iatomtype][jatomtype])
    error->all(FLERR,"Fix bond/create cutoff is longer than pairwise cutoff");

  // require special bonds = *,1,1
  // [0] can be anything b/c duplicate bond is checked for
  // [1],[2] must be 1.0 b/c only special lists of I,J are updated when
  //   bond I-J is created, not special lists of neighbors of I,J,etc

  if (force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0)
    error->all(FLERR,"Fix bond/create requires special_bonds lj = *,1,1");

  if (atom->q_flag)
    if (force->special_coul[2] != 1.0 || force->special_coul[3] != 1.0)
      error->all(FLERR,"Fix bond/create requires special_bonds coul = *,1,1");

  // warn if angles, dihedrals, impropers are being used

  if (force->angle || force->dihedral || force->improper) {
    if (me == 0)
      error->warning(FLERR,"Created bonds will not create angles, "
                     "dihedrals, or impropers");
  }

  // need a half neighbor list, built when ever re-neighboring occurs

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::setup(int vflag)
{
  int i,j,m;

  // compute initial bondcount if this is first run
  // can't do this earlier, like in constructor or init, b/c need ghost info

  if (countflag) return;
  countflag = 1;

  // count bonds stored with each bond I own
  // if newton bond is not set, just increment count on atom I
  // if newton bond is set, also increment count on atom J even if ghost
  // bondcount is long enough to tally ghost atom counts

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  int newton_bond = force->newton_bond;

  for (i = 0; i < nall; i++) bondcount[i] = 0;

  for (i = 0; i < nlocal; i++)
    for (j = 0; j < num_bond[i]; j++) {
      if (bond_type[i][j] == btype) {
        bondcount[i]++;
        if (newton_bond) {
          m = atom->map(bond_atom[i][j]);
          if (m < 0)
            error->one(FLERR,
                       "Could not count initial bonds in fix bond/create");
          bondcount[m]++;
        }
      }
    }

  // if newton_bond is set, need to sum bondcount

  commflag = 0;
  if (newton_bond) comm->reverse_comm_fix(this);
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::post_integrate()
{
  int i,j,k,m,ii,jj,inum,jnum,itype,jtype,n1,n3,possible;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  tagint *slist;

  if (update->ntimestep % nevery) return;

  // need updated ghost atom positions

  comm->forward_comm();

  // forward comm of bondcount, so ghosts have it

  commflag = 0;
  comm->forward_comm_fix(this);

  // resize bond partner list and initialize it
  // probability array overlays distsq array
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(partner);
    memory->destroy(distsq);
    nmax = atom->nmax;
    memory->create(partner,nmax,"bond/create:partner");
    memory->create(distsq,nmax,"bond/create:distsq");
    probability = distsq;
  }

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    partner[i] = 0;
    distsq[i] = BIG;
  }

  // loop over neighbors of my atoms
  // each atom sets one closest eligible partner atom ID to bond with

  double **x = atom->x;
  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  int *mask = atom->mask;
  int *type = atom->type;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      if (!(mask[j] & groupbit)) continue;
      jtype = type[j];

      possible = 0;
      if (itype == iatomtype && jtype == jatomtype) {
        if ((imaxbond == 0 || bondcount[i] < imaxbond) &&
            (jmaxbond == 0 || bondcount[j] < jmaxbond))
          possible = 1;
      } else if (itype == jatomtype && jtype == iatomtype) {
        if ((jmaxbond == 0 || bondcount[i] < jmaxbond) &&
            (imaxbond == 0 || bondcount[j] < imaxbond))
          possible = 1;
      }
      if (!possible) continue;

      // do not allow a duplicate bond to be created
      // check existing bonds of both I and J

      for (k = 0; k < num_bond[i]; k++)
        if (bond_atom[i][k] == tag[j]) possible = 0;
      if (j < nlocal)
        for (k = 0; k < num_bond[j]; k++)
          if (bond_atom[j][k] == tag[i]) possible = 0;
      if (!possible) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq >= cutsq) continue;

      if (rsq < distsq[i]) {
        partner[i] = tag[j];
        distsq[i] = rsq;
      }
      if (rsq < distsq[j]) {
        partner[j] = tag[i];
        distsq[j] = rsq;
      }
    }
  }

  // reverse comm of distsq and partner
  // not needed if newton_pair off since I,J pair was seen by both procs

  commflag = 1;
  if (force->newton_pair) comm->reverse_comm_fix(this);

  // each atom now knows its winning partner
  // for prob check, generate random value for each atom with a bond partner
  // forward comm of partner and random value, so ghosts have it

  if (fraction < 1.0) {
    for (i = 0; i < nlocal; i++)
      if (partner[i]) probability[i] = random->uniform();
  }

  commflag = 1;
  comm->forward_comm_fix(this);

  // create bonds for atoms I own
  // if other atom is owned by another proc, it should create same bond
  // if both atoms list each other as winning bond partner
  // and probability constraint is satisfied

  int **bond_type = atom->bond_type;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int newton_bond = force->newton_bond;

  int ncreate = 0;
  for (i = 0; i < nlocal; i++) {
    if (partner[i] == 0) continue;
    j = atom->map(partner[i]);
    if (partner[j] != tag[i]) continue;

    // apply probability constraint using RN for atom with smallest ID

    if (fraction < 1.0) {
      if (tag[i] < tag[j]) {
        if (probability[i] >= fraction) continue;
      } else {
        if (probability[j] >= fraction) continue;
      }
    }

    // if newton_bond is set, only store with I or J
    // if not newton_bond, store bond with both I and J

    if (!newton_bond || tag[i] < tag[j]) {
      if (num_bond[i] == atom->bond_per_atom)
        error->one(FLERR,"New bond exceeded bonds per atom in fix bond/create");
      bond_type[i][num_bond[i]] = btype;
      bond_atom[i][num_bond[i]] = tag[j];
      num_bond[i]++;
    }

    // add a 1-2 neighbor to special bond list for atom I
    // atom J will also do this

    slist = special[i];
    n1 = nspecial[i][0];
    n3 = nspecial[i][2];
    if (n3 == atom->maxspecial)
      error->one(FLERR,
                 "New bond exceeded special list size in fix bond/create");
    for (m = n3; m > n1; m--) slist[m+1] = slist[m];
    slist[n1] = tag[j];
    nspecial[i][0]++;
    nspecial[i][1]++;
    nspecial[i][2]++;

    // increment bondcount, convert atom to new type if limit reached

    bondcount[i]++;
    if (type[i] == iatomtype)
      if (bondcount[i] == imaxbond) type[i] = inewtype;
    else
      if (bondcount[i] == jmaxbond) type[i] = jnewtype;

    // count the created bond once

    if (tag[i] < tag[j]) ncreate++;
  }

  // tally stats

  MPI_Allreduce(&ncreate,&createcount,1,MPI_INT,MPI_SUM,world);
  createcounttotal += createcount;
  atom->nbonds += createcount;

  // trigger reneighboring if any bonds were formed

  if (createcount) next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::post_integrate_respa(int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_integrate();
}

/* ---------------------------------------------------------------------- */

int FixBondCreate::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;

  if (commflag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = bondcount[j];
    }
    return 1;

  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = partner[j];
      buf[m++] = probability[j];
    }
    return 2;
  }
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (commflag == 0) {
    for (i = first; i < last; i++)
      bondcount[i] = static_cast<int> (buf[m++]);

  } else {
    for (i = first; i < last; i++) {
      partner[i] = static_cast<tagint> (buf[m++]);
      probability[i] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixBondCreate::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (commflag == 0) {
    for (i = first; i < last; i++)
      buf[m++] = bondcount[i];
    return 1;

  } else {
    for (i = first; i < last; i++) {
      buf[m++] = distsq[i];
      buf[m++] = partner[i];
    }
    return 2;
  }
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;

  if (commflag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      bondcount[j] += static_cast<int> (buf[m++]);
    }

  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      if (buf[m] < distsq[j]) {
        distsq[j] = buf[m++];
        partner[j] = static_cast<tagint> (buf[m++]);
      } else m += 2;
    }
  }
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixBondCreate::grow_arrays(int nmax)
{
  memory->grow(bondcount,nmax,"bond/create:bondcount");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixBondCreate::copy_arrays(int i, int j, int delflag)
{
  bondcount[j] = bondcount[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixBondCreate::pack_exchange(int i, double *buf)
{
  buf[0] = bondcount[i];
  return 1;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixBondCreate::unpack_exchange(int nlocal, double *buf)
{
  bondcount[nlocal] = static_cast<int> (buf[0]);
  return 1;
}

/* ---------------------------------------------------------------------- */

double FixBondCreate::compute_vector(int n)
{
  if (n == 1) return (double) createcount;
  return (double) createcounttotal;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixBondCreate::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax*2 * sizeof(int);
  bytes += nmax * sizeof(double);
  return bytes;
}
