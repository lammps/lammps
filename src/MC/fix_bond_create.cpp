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

#include "fix_bond_create.h"
#include <mpi.h>
#include <cstring>
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "force.h"
#include "modify.h"
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
#define DELTA 16

/* ---------------------------------------------------------------------- */

FixBondCreate::FixBondCreate(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  bondcount(NULL), partner(NULL), finalpartner(NULL), distsq(NULL),
  probability(NULL), created(NULL), copy(NULL), random(NULL), list(NULL)
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
  atype = dtype = itype = 0;

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
    } else if (strcmp(arg[iarg],"atype") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/create command");
      atype = force->inumeric(FLERR,arg[iarg+1]);
      if (atype < 0) error->all(FLERR,"Illegal fix bond/create command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"dtype") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/create command");
      dtype = force->inumeric(FLERR,arg[iarg+1]);
      if (dtype < 0) error->all(FLERR,"Illegal fix bond/create command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"itype") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/create command");
      itype = force->inumeric(FLERR,arg[iarg+1]);
      if (itype < 0) error->all(FLERR,"Illegal fix bond/create command");
      iarg += 2;
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
  // forward is big due to comm of broken bonds and 1-2 neighbors

  comm_forward = MAX(2,2+atom->maxspecial);
  comm_reverse = 2;

  // allocate arrays local to this fix

  nmax = 0;
  partner = finalpartner = NULL;
  distsq = NULL;

  maxcreate = 0;
  created = NULL;

  // copy = special list for one atom
  // size = ms^2 + ms is sufficient
  // b/c in rebuild_special_one() neighs of all 1-2s are added,
  //   then a dedup(), then neighs of all 1-3s are added, then final dedup()
  // this means intermediate size cannot exceed ms^2 + ms

  int maxspecial = atom->maxspecial;
  copy = new tagint[maxspecial*maxspecial + maxspecial];

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
  memory->destroy(finalpartner);
  memory->destroy(distsq);
  memory->destroy(created);
  delete [] copy;
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
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // check cutoff for iatomtype,jatomtype

  if (force->pair == NULL || cutsq > force->pair->cutsq[iatomtype][jatomtype])
    error->all(FLERR,"Fix bond/create cutoff is longer than pairwise cutoff");

  // warn if more than one fix bond/create or also a fix bond/break
  // because this fix stores per-atom state in bondcount
  //   if other fixes create/break bonds, this fix will not know about it

  int count = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"bond/create") == 0) count++;
    if (strcmp(modify->fix[i]->style,"bond/break") == 0) count++;
  }
  if (count > 1 && me == 0)
    error->warning(FLERR,"Fix bond/create is used multiple times "
                   " or with fix bond/break - may not work as expected");

  // enable angle/dihedral/improper creation if atype/dtype/itype
  //   option was used and a force field has been specified

  if (atype && force->angle) {
    angleflag = 1;
    if (atype > atom->nangletypes)
      error->all(FLERR,"Fix bond/create angle type is invalid");
  } else angleflag = 0;

  if (dtype && force->dihedral) {
    dihedralflag = 1;
    if (dtype > atom->ndihedraltypes)
      error->all(FLERR,"Fix bond/create dihedral type is invalid");
  } else dihedralflag = 0;

  if (itype && force->improper) {
    improperflag = 1;
    if (itype > atom->nimpropertypes)
      error->all(FLERR,"Fix bond/create improper type is invalid");
  } else improperflag = 0;

  if (force->improper) {
    if (force->improper_match("class2") || force->improper_match("ring"))
      error->all(FLERR,"Cannot yet use fix bond/create with this "
                 "improper style");
  }

  // need a half neighbor list, built every Nevery steps

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->occasional = 1;

  lastcheck = -1;
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::setup(int /*vflag*/)
{
  int i,j,m;

  // compute initial bondcount if this is first run
  // can't do this earlier, in constructor or init, b/c need ghost info

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
            error->one(FLERR,"Fix bond/create needs ghost atoms "
                       "from further away");
          bondcount[m]++;
        }
      }
    }

  // if newton_bond is set, need to sum bondcount

  commflag = 1;
  if (newton_bond) comm->reverse_comm_fix(this,1);
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::post_integrate()
{
  int i,j,k,m,n,ii,jj,inum,jnum,itype,jtype,n1,n2,n3,possible;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  tagint *slist;

  if (update->ntimestep % nevery) return;

  // check that all procs have needed ghost atoms within ghost cutoff
  // only if neighbor list has changed since last check
  // needs to be <= test b/c neighbor list could have been re-built in
  //   same timestep as last post_integrate() call, but afterwards
  // NOTE: no longer think is needed, due to error tests on atom->map()
  // NOTE: if delete, can also delete lastcheck and check_ghosts()

  //if (lastcheck <= neighbor->lastcall) check_ghosts();

  // acquire updated ghost atom positions
  // necessary b/c are calling this after integrate, but before Verlet comm

  comm->forward_comm();

  // forward comm of bondcount, so ghosts have it

  commflag = 1;
  comm->forward_comm_fix(this,1);

  // resize bond partner list and initialize it
  // probability array overlays distsq array
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(partner);
    memory->destroy(finalpartner);
    memory->destroy(distsq);
    nmax = atom->nmax;
    memory->create(partner,nmax,"bond/create:partner");
    memory->create(finalpartner,nmax,"bond/create:finalpartner");
    memory->create(distsq,nmax,"bond/create:distsq");
    probability = distsq;
  }

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    partner[i] = 0;
    finalpartner[i] = 0;
    distsq[i] = BIG;
  }

  // loop over neighbors of my atoms
  // each atom sets one closest eligible partner atom ID to bond with

  double **x = atom->x;
  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int *mask = atom->mask;
  int *type = atom->type;

  neighbor->build_one(list,1);
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
      // check 1-2 neighbors of atom I

      for (k = 0; k < nspecial[i][0]; k++)
        if (special[i][k] == tag[j]) possible = 0;
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

  commflag = 2;
  if (force->newton_pair) comm->reverse_comm_fix(this);

  // each atom now knows its winning partner
  // for prob check, generate random value for each atom with a bond partner
  // forward comm of partner and random value, so ghosts have it

  if (fraction < 1.0) {
    for (i = 0; i < nlocal; i++)
      if (partner[i]) probability[i] = random->uniform();
  }

  commflag = 2;
  comm->forward_comm_fix(this,2);

  // create bonds for atoms I own
  // only if both atoms list each other as winning bond partner
  //   and probability constraint is satisfied
  // if other atom is owned by another proc, it should do same thing

  int **bond_type = atom->bond_type;
  int newton_bond = force->newton_bond;

  ncreate = 0;
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
    // atom J will also do this consistently, whatever proc it is on

    if (!newton_bond || tag[i] < tag[j]) {
      if (num_bond[i] == atom->bond_per_atom)
        error->one(FLERR,"New bond exceeded bonds per atom in fix bond/create");
      bond_type[i][num_bond[i]] = btype;
      bond_atom[i][num_bond[i]] = tag[j];
      num_bond[i]++;
    }

    // add a 1-2 neighbor to special bond list for atom I
    // atom J will also do this, whatever proc it is on
    // need to first remove tag[j] from later in list if it appears
    // prevents list from overflowing, will be rebuilt in rebuild_special_one()

    slist = special[i];
    n1 = nspecial[i][0];
    n2 = nspecial[i][1];
    n3 = nspecial[i][2];
    for (m = n1; m < n3; m++)
      if (slist[m] == tag[j]) break;
    if (m < n3) {
      for (n = m; n < n3-1; n++) slist[n] = slist[n+1];
      n3--;
      if (m < n2) n2--;
    }
    if (n3 == atom->maxspecial)
      error->one(FLERR,
                 "New bond exceeded special list size in fix bond/create");
    for (m = n3; m > n1; m--) slist[m] = slist[m-1];
    slist[n1] = tag[j];
    nspecial[i][0] = n1+1;
    nspecial[i][1] = n2+1;
    nspecial[i][2] = n3+1;

    // increment bondcount, convert atom to new type if limit reached
    // atom J will also do this, whatever proc it is on

    bondcount[i]++;
    if (type[i] == iatomtype) {
      if (bondcount[i] == imaxbond) type[i] = inewtype;
    } else {
      if (bondcount[i] == jmaxbond) type[i] = jnewtype;
    }

    // store final created bond partners and count the created bond once

    finalpartner[i] = tag[j];
    finalpartner[j] = tag[i];
    if (tag[i] < tag[j]) ncreate++;
  }

  // tally stats

  MPI_Allreduce(&ncreate,&createcount,1,MPI_INT,MPI_SUM,world);
  createcounttotal += createcount;
  atom->nbonds += createcount;

  // trigger reneighboring if any bonds were formed
  // this insures neigh lists will immediately reflect the topology changes
  // done if any bonds created

  if (createcount) next_reneighbor = update->ntimestep;
  if (!createcount) return;

  // communicate final partner and 1-2 special neighbors
  // 1-2 neighs already reflect created bonds

  commflag = 3;
  comm->forward_comm_fix(this);

  // create list of broken bonds that influence my owned atoms
  //   even if between owned-ghost or ghost-ghost atoms
  // finalpartner is now set for owned and ghost atoms so loop over nall
  // OK if duplicates in broken list due to ghosts duplicating owned atoms
  // check J < 0 to insure a broken bond to unknown atom is included
  //   i.e. a bond partner outside of cutoff length

  ncreate = 0;
  for (i = 0; i < nall; i++) {
    if (finalpartner[i] == 0) continue;
    j = atom->map(finalpartner[i]);
    if (j < 0 || tag[i] < tag[j]) {
      if (ncreate == maxcreate) {
        maxcreate += DELTA;
        memory->grow(created,maxcreate,2,"bond/create:created");
      }
      created[ncreate][0] = tag[i];
      created[ncreate][1] = finalpartner[i];
      ncreate++;
    }
  }

  // update special neigh lists of all atoms affected by any created bond
  // also add angles/dihedrals/impropers induced by created bonds

  update_topology();

  // DEBUG
  //print_bb();
}

/* ----------------------------------------------------------------------
   insure all atoms 2 hops away from owned atoms are in ghost list
   this allows dihedral 1-2-3-4 to be properly created
     and special list of 1 to be properly updated
   if I own atom 1, but not 2,3,4, and bond 3-4 is added
     then 2,3 will be ghosts and 3 will store 4 as its finalpartner
------------------------------------------------------------------------- */

void FixBondCreate::check_ghosts()
{
  int i,j,n;
  tagint *slist;

  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int nlocal = atom->nlocal;

  int flag = 0;
  for (i = 0; i < nlocal; i++) {
    slist = special[i];
    n = nspecial[i][1];
    for (j = 0; j < n; j++)
      if (atom->map(slist[j]) < 0) flag = 1;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall)
    error->all(FLERR,"Fix bond/create needs ghost atoms from further away");
  lastcheck = update->ntimestep;
}

/* ----------------------------------------------------------------------
   double loop over my atoms and created bonds
   influenced = 1 if atom's topology is affected by any created bond
     yes if is one of 2 atoms in bond
     yes if either atom ID appears in as 1-2 or 1-3 in atom's special list
     else no
   if influenced by any created bond:
     rebuild the atom's special list of 1-2,1-3,1-4 neighs
     check for angles/dihedrals/impropers to create due modified special list
------------------------------------------------------------------------- */

void FixBondCreate::update_topology()
{
  int i,j,k,n,influence,influenced;
  tagint id1,id2;
  tagint *slist;

  tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int nlocal = atom->nlocal;

  nangles = 0;
  ndihedrals = 0;
  nimpropers = 0;
  overflow = 0;

  //printf("NCREATE %d: ",ncreate);
  //for (i = 0; i < ncreate; i++)
  //  printf(" %d %d,",created[i][0],created[i][1]);
  //printf("\n");

  for (i = 0; i < nlocal; i++) {
    influenced = 0;
    slist = special[i];

    for (j = 0; j < ncreate; j++) {
      id1 = created[j][0];
      id2 = created[j][1];

      influence = 0;
      if (tag[i] == id1 || tag[i] == id2) influence = 1;
      else {
        n = nspecial[i][1];
        for (k = 0; k < n; k++)
          if (slist[k] == id1 || slist[k] == id2) {
            influence = 1;
            break;
          }
      }
      if (!influence) continue;
      influenced = 1;
    }

    // rebuild_special_one() first, since used by create_angles, etc

    if (influenced) {
      rebuild_special_one(i);
      if (angleflag) create_angles(i);
      if (dihedralflag) create_dihedrals(i);
      if (improperflag) create_impropers(i);
    }
  }

  int overflowall;
  MPI_Allreduce(&overflow,&overflowall,1,MPI_INT,MPI_SUM,world);
  if (overflowall) error->all(FLERR,"Fix bond/create induced too many "
                              "angles/dihedrals/impropers per atom");

  int newton_bond = force->newton_bond;

  int all;
  if (angleflag) {
    MPI_Allreduce(&nangles,&all,1,MPI_INT,MPI_SUM,world);
    if (!newton_bond) all /= 3;
    atom->nangles += all;
  }
  if (dihedralflag) {
    MPI_Allreduce(&ndihedrals,&all,1,MPI_INT,MPI_SUM,world);
    if (!newton_bond) all /= 4;
    atom->ndihedrals += all;
  }
  if (improperflag) {
    MPI_Allreduce(&nimpropers,&all,1,MPI_INT,MPI_SUM,world);
    if (!newton_bond) all /= 4;
    atom->nimpropers += all;
  }
}

/* ----------------------------------------------------------------------
   re-build special list of atom M
   does not affect 1-2 neighs (already include effects of new bond)
   affects 1-3 and 1-4 neighs due to other atom's augmented 1-2 neighs
------------------------------------------------------------------------- */

void FixBondCreate::rebuild_special_one(int m)
{
  int i,j,n,n1,cn1,cn2,cn3;
  tagint *slist;

  tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  // existing 1-2 neighs of atom M

  slist = special[m];
  n1 = nspecial[m][0];
  cn1 = 0;
  for (i = 0; i < n1; i++)
    copy[cn1++] = slist[i];

  // new 1-3 neighs of atom M, based on 1-2 neighs of 1-2 neighs
  // exclude self
  // remove duplicates after adding all possible 1-3 neighs

  cn2 = cn1;
  for (i = 0; i < cn1; i++) {
    n = atom->map(copy[i]);
    if (n < 0)
      error->one(FLERR,"Fix bond/create needs ghost atoms from further away");
    slist = special[n];
    n1 = nspecial[n][0];
    for (j = 0; j < n1; j++)
      if (slist[j] != tag[m]) copy[cn2++] = slist[j];
  }

  cn2 = dedup(cn1,cn2,copy);
  if (cn2 > atom->maxspecial)
    error->one(FLERR,"Special list size exceeded in fix bond/create");

  // new 1-4 neighs of atom M, based on 1-2 neighs of 1-3 neighs
  // exclude self
  // remove duplicates after adding all possible 1-4 neighs

  cn3 = cn2;
  for (i = cn1; i < cn2; i++) {
    n = atom->map(copy[i]);
    if (n < 0)
      error->one(FLERR,"Fix bond/create needs ghost atoms from further away");
    slist = special[n];
    n1 = nspecial[n][0];
    for (j = 0; j < n1; j++)
      if (slist[j] != tag[m]) copy[cn3++] = slist[j];
  }

  cn3 = dedup(cn2,cn3,copy);
  if (cn3 > atom->maxspecial)
    error->one(FLERR,"Special list size exceeded in fix bond/create");

  // store new special list with atom M

  nspecial[m][0] = cn1;
  nspecial[m][1] = cn2;
  nspecial[m][2] = cn3;
  memcpy(special[m],copy,cn3*sizeof(int));
}

/* ----------------------------------------------------------------------
   create any angles owned by atom M induced by newly created bonds
   walk special list to find all possible angles to create
   only add an angle if a new bond is one of its 2 bonds (I-J,J-K)
   for newton_bond on, atom M is central atom
   for newton_bond off, atom M is any of 3 atoms in angle
------------------------------------------------------------------------- */

void FixBondCreate::create_angles(int m)
{
  int i,j,n,i2local,n1,n2;
  tagint i1,i2,i3;
  tagint *s1list,*s2list;

  tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  int num_angle = atom->num_angle[m];
  int *angle_type = atom->angle_type[m];
  tagint *angle_atom1 = atom->angle_atom1[m];
  tagint *angle_atom2 = atom->angle_atom2[m];
  tagint *angle_atom3 = atom->angle_atom3[m];

  // atom M is central atom in angle
  // double loop over 1-2 neighs
  // avoid double counting by 2nd loop as j = i+1,N not j = 1,N
  // consider all angles, only add if:
  //   a new bond is in the angle and atom types match

  i2 = tag[m];
  n2 = nspecial[m][0];
  s2list = special[m];

  for (i = 0; i < n2; i++) {
    i1 = s2list[i];
    for (j = i+1; j < n2; j++) {
      i3 = s2list[j];

      // angle = i1-i2-i3

      for (n = 0; n < ncreate; n++) {
        if (created[n][0] == i1 && created[n][1] == i2) break;
        if (created[n][0] == i2 && created[n][1] == i1) break;
        if (created[n][0] == i2 && created[n][1] == i3) break;
        if (created[n][0] == i3 && created[n][1] == i2) break;
      }
      if (n == ncreate) continue;

      // NOTE: this is place to check atom types of i1,i2,i3

      if (num_angle < atom->angle_per_atom) {
        angle_type[num_angle] = atype;
        angle_atom1[num_angle] = i1;
        angle_atom2[num_angle] = i2;
        angle_atom3[num_angle] = i3;
        num_angle++;
        nangles++;
      } else overflow = 1;
    }
  }

  atom->num_angle[m] = num_angle;
  if (force->newton_bond) return;

  // for newton_bond off, also consider atom M as atom 1 in angle

  i1 = tag[m];
  n1 = nspecial[m][0];
  s1list = special[m];

  for (i = 0; i < n1; i++) {
    i2 = s1list[i];
    i2local = atom->map(i2);
    if (i2local < 0)
      error->one(FLERR,"Fix bond/create needs ghost atoms from further away");
    s2list = special[i2local];
    n2 = nspecial[i2local][0];

    for (j = 0; j < n2; j++) {
      i3 = s2list[j];
      if (i3 == i1) continue;

      // angle = i1-i2-i3

      for (n = 0; n < ncreate; n++) {
        if (created[n][0] == i1 && created[n][1] == i2) break;
        if (created[n][0] == i2 && created[n][1] == i1) break;
        if (created[n][0] == i2 && created[n][1] == i3) break;
        if (created[n][0] == i3 && created[n][1] == i2) break;
      }
      if (n == ncreate) continue;

      // NOTE: this is place to check atom types of i1,i2,i3

      if (num_angle < atom->angle_per_atom) {
        angle_type[num_angle] = atype;
        angle_atom1[num_angle] = i1;
        angle_atom2[num_angle] = i2;
        angle_atom3[num_angle] = i3;
        num_angle++;
        nangles++;
      } else overflow = 1;
    }
  }

  atom->num_angle[m] = num_angle;
}

/* ----------------------------------------------------------------------
   create any dihedrals owned by atom M induced by newly created bonds
   walk special list to find all possible dihedrals to create
   only add a dihedral if a new bond is one of its 3 bonds (I-J,J-K,K-L)
   for newton_bond on, atom M is central atom
   for newton_bond off, atom M is any of 4 atoms in dihedral
------------------------------------------------------------------------- */

void FixBondCreate::create_dihedrals(int m)
{
  int i,j,k,n,i1local,i2local,i3local,n1,n2,n3;
  tagint i1,i2,i3,i4;
  tagint *s1list,*s2list,*s3list;

  tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  int num_dihedral = atom->num_dihedral[m];
  int *dihedral_type = atom->dihedral_type[m];
  tagint *dihedral_atom1 = atom->dihedral_atom1[m];
  tagint *dihedral_atom2 = atom->dihedral_atom2[m];
  tagint *dihedral_atom3 = atom->dihedral_atom3[m];
  tagint *dihedral_atom4 = atom->dihedral_atom4[m];

  // atom M is 2nd atom in dihedral
  // double loop over 1-2 neighs
  // two triple loops: one over neighs at each end of triplet
  // avoid double counting by 2nd loop as j = i+1,N not j = 1,N
  // avoid double counting due to another atom being 2nd atom in same dihedral
  //   by requiring ID of 2nd atom < ID of 3rd atom
  //   don't do this if newton bond off since want to double count
  // consider all dihedrals, only add if:
  //   a new bond is in the dihedral and atom types match

  i2 = tag[m];
  n2 = nspecial[m][0];
  s2list = special[m];

  for (i = 0; i < n2; i++) {
    i1 = s2list[i];

    for (j = i+1; j < n2; j++) {
      i3 = s2list[j];
      if (force->newton_bond && i2 > i3) continue;
      i3local = atom->map(i3);
      if (i3local < 0)
        error->one(FLERR,"Fix bond/create needs ghost atoms from further away");
      s3list = special[i3local];
      n3 = nspecial[i3local][0];

      for (k = 0; k < n3; k++) {
        i4 = s3list[k];
        if (i4 == i1 || i4 == i2 || i4 == i3) continue;

        // dihedral = i1-i2-i3-i4

        for (n = 0; n < ncreate; n++) {
          if (created[n][0] == i1 && created[n][1] == i2) break;
          if (created[n][0] == i2 && created[n][1] == i1) break;
          if (created[n][0] == i2 && created[n][1] == i3) break;
          if (created[n][0] == i3 && created[n][1] == i2) break;
          if (created[n][0] == i3 && created[n][1] == i4) break;
          if (created[n][0] == i4 && created[n][1] == i3) break;
        }
        if (n < ncreate) {
          // NOTE: this is place to check atom types of i3,i2,i1,i4
          if (num_dihedral < atom->dihedral_per_atom) {
            dihedral_type[num_dihedral] = dtype;
            dihedral_atom1[num_dihedral] = i1;
            dihedral_atom2[num_dihedral] = i2;
            dihedral_atom3[num_dihedral] = i3;
            dihedral_atom4[num_dihedral] = i4;
            num_dihedral++;
            ndihedrals++;
          } else overflow = 1;
        }
      }
    }
  }

  for (i = 0; i < n2; i++) {
    i1 = s2list[i];
    if (force->newton_bond && i2 > i1) continue;
    i1local = atom->map(i1);
    if (i1local < 0)
      error->one(FLERR,"Fix bond/create needs ghost atoms from further away");
    s3list = special[i1local];
    n3 = nspecial[i1local][0];

    for (j = i+1; j < n2; j++) {
      i3 = s2list[j];

      for (k = 0; k < n3; k++) {
        i4 = s3list[k];
        if (i4 == i1 || i4 == i2 || i4 == i3) continue;

        // dihedral = i3-i2-i1-i4

        for (n = 0; n < ncreate; n++) {
          if (created[n][0] == i3 && created[n][1] == i2) break;
          if (created[n][0] == i2 && created[n][1] == i3) break;
          if (created[n][0] == i2 && created[n][1] == i1) break;
          if (created[n][0] == i1 && created[n][1] == i2) break;
          if (created[n][0] == i1 && created[n][1] == i4) break;
          if (created[n][0] == i4 && created[n][1] == i1) break;
        }
        if (n < ncreate) {
          // NOTE: this is place to check atom types of i3,i2,i1,i4
          if (num_dihedral < atom->dihedral_per_atom) {
            dihedral_type[num_dihedral] = dtype;
            dihedral_atom1[num_dihedral] = i3;
            dihedral_atom2[num_dihedral] = i2;
            dihedral_atom3[num_dihedral] = i1;
            dihedral_atom4[num_dihedral] = i4;
            num_dihedral++;
            ndihedrals++;
          } else overflow = 1;
        }
      }
    }
  }

  atom->num_dihedral[m] = num_dihedral;
  if (force->newton_bond) return;

  // for newton_bond off, also consider atom M as atom 1 in dihedral

  i1 = tag[m];
  n1 = nspecial[m][0];
  s1list = special[m];

  for (i = 0; i < n1; i++) {
    i2 = s1list[i];
    i2local = atom->map(i2);
    if (i2local < 0)
      error->one(FLERR,"Fix bond/create needs ghost atoms from further away");
    s2list = special[i2local];
    n2 = nspecial[i2local][0];

    for (j = 0; j < n2; j++) {
      i3 = s2list[j];
      if (i3 == i1) continue;
      i3local = atom->map(i3);
      if (i3local < 0)
        error->one(FLERR,"Fix bond/create needs ghost atoms from further away");
      s3list = special[i3local];
      n3 = nspecial[i3local][0];

      for (k = 0; k < n3; k++) {
        i4 = s3list[k];
        if (i4 == i1 || i4 == i2 || i4 == i3) continue;

        // dihedral = i1-i2-i3-i4

        for (n = 0; n < ncreate; n++) {
          if (created[n][0] == i1 && created[n][1] == i2) break;
          if (created[n][0] == i2 && created[n][1] == i1) break;
          if (created[n][0] == i2 && created[n][1] == i3) break;
          if (created[n][0] == i3 && created[n][1] == i2) break;
          if (created[n][0] == i3 && created[n][1] == i4) break;
          if (created[n][0] == i4 && created[n][1] == i3) break;
        }
        if (n < ncreate) {
          // NOTE: this is place to check atom types of i3,i2,i1,i4
          if (num_dihedral < atom->dihedral_per_atom) {
            dihedral_type[num_dihedral] = dtype;
            dihedral_atom1[num_dihedral] = i1;
            dihedral_atom2[num_dihedral] = i2;
            dihedral_atom3[num_dihedral] = i3;
            dihedral_atom4[num_dihedral] = i4;
            num_dihedral++;
            ndihedrals++;
          } else overflow = 1;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   create any impropers owned by atom M induced by newly created bonds
   walk special list to find all possible impropers to create
   only add an improper if a new bond is one of its 3 bonds (I-J,I-K,I-L)
   for newton_bond on, atom M is central atom
   for newton_bond off, atom M is any of 4 atoms in improper
------------------------------------------------------------------------- */

void FixBondCreate::create_impropers(int m)
{
  int i,j,k,n,i1local,n1,n2;
  tagint i1,i2,i3,i4;
  tagint *s1list,*s2list;

  tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  int num_improper = atom->num_improper[m];
  int *improper_type = atom->improper_type[m];
  tagint *improper_atom1 = atom->improper_atom1[m];
  tagint *improper_atom2 = atom->improper_atom2[m];
  tagint *improper_atom3 = atom->improper_atom3[m];
  tagint *improper_atom4 = atom->improper_atom4[m];

  // atom M is central atom in improper
  // triple loop over 1-2 neighs
  // avoid double counting by 2nd loop as j = i+1,N not j = 1,N
  // consider all impropers, only add if:
  //   a new bond is in the improper and atom types match

  i1 = tag[m];
  n1 = nspecial[m][0];
  s1list = special[m];

  for (i = 0; i < n1; i++) {
    i2 = s1list[i];
    for (j = i+1; j < n1; j++) {
      i3 = s1list[j];
      for (k = j+1; k < n1; k++) {
        i4 = s1list[k];

        // improper = i1-i2-i3-i4

        for (n = 0; n < ncreate; n++) {
          if (created[n][0] == i1 && created[n][1] == i2) break;
          if (created[n][0] == i2 && created[n][1] == i1) break;
          if (created[n][0] == i1 && created[n][1] == i3) break;
          if (created[n][0] == i3 && created[n][1] == i1) break;
          if (created[n][0] == i1 && created[n][1] == i4) break;
          if (created[n][0] == i4 && created[n][1] == i1) break;
        }
        if (n == ncreate) continue;

        // NOTE: this is place to check atom types of i1,i2,i3,i4

        if (num_improper < atom->improper_per_atom) {
          improper_type[num_improper] = itype;
          improper_atom1[num_improper] = i1;
          improper_atom2[num_improper] = i2;
          improper_atom3[num_improper] = i3;
          improper_atom4[num_improper] = i4;
          num_improper++;
          nimpropers++;
        } else overflow = 1;
      }
    }
  }

  atom->num_improper[m] = num_improper;
  if (force->newton_bond) return;

  // for newton_bond off, also consider atom M as atom 2 in improper

  i2 = tag[m];
  n2 = nspecial[m][0];
  s2list = special[m];

  for (i = 0; i < n2; i++) {
    i1 = s2list[i];
    i1local = atom->map(i1);
    if (i1local < 0)
      error->one(FLERR,"Fix bond/create needs ghost atoms from further away");
    s1list = special[i1local];
    n1 = nspecial[i1local][0];

    for (j = 0; j < n1; j++) {
      i3 = s1list[j];
      if (i3 == i1 || i3 == i2) continue;

      for (k = j+1; k < n1; k++) {
        i4 = s1list[k];
        if (i4 == i1 || i4 == i2) continue;

        // improper = i1-i2-i3-i4

        for (n = 0; n < ncreate; n++) {
          if (created[n][0] == i1 && created[n][1] == i2) break;
          if (created[n][0] == i2 && created[n][1] == i1) break;
          if (created[n][0] == i1 && created[n][1] == i3) break;
          if (created[n][0] == i3 && created[n][1] == i1) break;
          if (created[n][0] == i1 && created[n][1] == i4) break;
          if (created[n][0] == i4 && created[n][1] == i1) break;
        }
        if (n < ncreate) {
          // NOTE: this is place to check atom types of i3,i2,i1,i4
          if (num_improper < atom->improper_per_atom) {
            improper_type[num_improper] = itype;
            improper_atom1[num_improper] = i1;
            improper_atom2[num_improper] = i2;
            improper_atom3[num_improper] = i3;
            improper_atom4[num_improper] = i4;
            num_improper++;
            nimpropers++;
          } else overflow = 1;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   remove all ID duplicates in copy from Nstart:Nstop-1
   compare to all previous values in copy
   return N decremented by any discarded duplicates
------------------------------------------------------------------------- */

int FixBondCreate::dedup(int nstart, int nstop, tagint *copy)
{
  int i;

  int m = nstart;
  while (m < nstop) {
    for (i = 0; i < m; i++)
      if (copy[i] == copy[m]) {
        copy[m] = copy[nstop-1];
        nstop--;
        break;
      }
    if (i == m) m++;
  }

  return nstop;
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::post_integrate_respa(int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_integrate();
}

/* ---------------------------------------------------------------------- */

int FixBondCreate::pack_forward_comm(int n, int *list, double *buf,
                                     int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,k,m,ns;

  m = 0;

  if (commflag == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(bondcount[j]).d;
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

  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = ubuf(finalpartner[j]).d;
    ns = nspecial[j][0];
    buf[m++] = ubuf(ns).d;
    for (k = 0; k < ns; k++)
      buf[m++] = ubuf(special[j][k]).d;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::unpack_forward_comm(int n, int first, double *buf)
{
  int i,j,m,ns,last;

  m = 0;
  last = first + n;

  if (commflag == 1) {
    for (i = first; i < last; i++)
      bondcount[i] = (int) ubuf(buf[m++]).i;

  } else if (commflag == 2) {
    for (i = first; i < last; i++) {
      partner[i] = (tagint) ubuf(buf[m++]).i;
      probability[i] = buf[m++];
    }

  } else {
    int **nspecial = atom->nspecial;
    tagint **special = atom->special;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      finalpartner[i] = (tagint) ubuf(buf[m++]).i;
      ns = (int) ubuf(buf[m++]).i;
      nspecial[i][0] = ns;
      for (j = 0; j < ns; j++)
        special[i][j] = (tagint) ubuf(buf[m++]).i;
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixBondCreate::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (commflag == 1) {
    for (i = first; i < last; i++)
      buf[m++] = ubuf(bondcount[i]).d;
    return m;
  }

  for (i = first; i < last; i++) {
    buf[m++] = ubuf(partner[i]).d;
    buf[m++] = distsq[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;

  if (commflag == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      bondcount[j] += (int) ubuf(buf[m++]).i;
    }

  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      if (buf[m+1] < distsq[j]) {
        partner[j] = (tagint) ubuf(buf[m++]).i;
        distsq[j] = buf[m++];
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

void FixBondCreate::copy_arrays(int i, int j, int /*delflag*/)
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
  if (n == 0) return (double) createcount;
  return (double) createcounttotal;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixBondCreate::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes = 2*nmax * sizeof(tagint);
  bytes += nmax * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::print_bb()
{
  for (int i = 0; i < atom->nlocal; i++) {
    printf("TAG " TAGINT_FORMAT ": %d nbonds: ",atom->tag[i],atom->num_bond[i]);
    for (int j = 0; j < atom->num_bond[i]; j++) {
      printf(" " TAGINT_FORMAT,atom->bond_atom[i][j]);
    }
    printf("\n");
    printf("TAG " TAGINT_FORMAT ": %d nangles: ",atom->tag[i],atom->num_angle[i]);
    for (int j = 0; j < atom->num_angle[i]; j++) {
      printf(" " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT ",",
             atom->angle_atom1[i][j], atom->angle_atom2[i][j],
             atom->angle_atom3[i][j]);
    }
    printf("\n");
    printf("TAG " TAGINT_FORMAT ": %d ndihedrals: ",atom->tag[i],atom->num_dihedral[i]);
    for (int j = 0; j < atom->num_dihedral[i]; j++) {
      printf(" " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " "
             TAGINT_FORMAT ",", atom->dihedral_atom1[i][j],
             atom->dihedral_atom2[i][j],atom->dihedral_atom3[i][j],
             atom->dihedral_atom4[i][j]);
    }
    printf("\n");
    printf("TAG " TAGINT_FORMAT ": %d nimpropers: ",atom->tag[i],atom->num_improper[i]);
    for (int j = 0; j < atom->num_improper[i]; j++) {
      printf(" " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " "
             TAGINT_FORMAT ",",atom->improper_atom1[i][j],
             atom->improper_atom2[i][j],atom->improper_atom3[i][j],
             atom->improper_atom4[i][j]);
    }
    printf("\n");
    printf("TAG " TAGINT_FORMAT ": %d %d %d nspecial: ",atom->tag[i],
           atom->nspecial[i][0],atom->nspecial[i][1],atom->nspecial[i][2]);
    for (int j = 0; j < atom->nspecial[i][2]; j++) {
      printf(" " TAGINT_FORMAT,atom->special[i][j]);
    }
    printf("\n");
  }
}

/* ---------------------------------------------------------------------- */

void FixBondCreate::print_copy(const char *str, tagint m,
                              int n1, int n2, int n3, int *v)
{
  printf("%s " TAGINT_FORMAT ": %d %d %d nspecial: ",str,m,n1,n2,n3);
  for (int j = 0; j < n3; j++) printf(" %d",v[j]);
  printf("\n");
}
