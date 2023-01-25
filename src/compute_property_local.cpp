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

#include "compute_property_local.h"

#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

enum { NONE, NEIGH, PAIR, BOND, ANGLE, DIHEDRAL, IMPROPER };
enum { TYPE, RADIUS };

#define DELTA 10000

/* ---------------------------------------------------------------------- */

ComputePropertyLocal::ComputePropertyLocal(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), vlocal(nullptr), alocal(nullptr), indices(nullptr),
    pack_choice(nullptr)
{
  if (narg < 4) error->all(FLERR, "Illegal compute property/local command");

  local_flag = 1;
  nvalues = narg - 3;
  pack_choice = new FnPtrPack[nvalues];

  kindflag = NONE;

  int i;
  nvalues = 0;
  int iarg = 3;
  while (iarg < narg) {
    i = iarg - 3;

    if (strcmp(arg[iarg], "natom1") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_patom1;
      if (kindflag != NONE && kindflag != NEIGH)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = NEIGH;
    } else if (strcmp(arg[iarg], "natom2") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_patom2;
      if (kindflag != NONE && kindflag != NEIGH)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = NEIGH;
    } else if (strcmp(arg[iarg], "ntype1") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_ptype1;
      if (kindflag != NONE && kindflag != NEIGH)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = NEIGH;
    } else if (strcmp(arg[iarg], "ntype2") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_ptype2;
      if (kindflag != NONE && kindflag != NEIGH)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = NEIGH;

    } else if (strcmp(arg[iarg], "patom1") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_patom1;
      if (kindflag != NONE && kindflag != PAIR)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = PAIR;
    } else if (strcmp(arg[iarg], "patom2") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_patom2;
      if (kindflag != NONE && kindflag != PAIR)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = PAIR;
    } else if (strcmp(arg[iarg], "ptype1") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_ptype1;
      if (kindflag != NONE && kindflag != PAIR)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = PAIR;
    } else if (strcmp(arg[iarg], "ptype2") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_ptype2;
      if (kindflag != NONE && kindflag != PAIR)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = PAIR;

    } else if (strcmp(arg[iarg], "batom1") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_batom1;
      if (kindflag != NONE && kindflag != BOND)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = BOND;
    } else if (strcmp(arg[iarg], "batom2") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_batom2;
      if (kindflag != NONE && kindflag != BOND)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = BOND;
    } else if (strcmp(arg[iarg], "btype") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_btype;
      if (kindflag != NONE && kindflag != BOND)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = BOND;

    } else if (strcmp(arg[iarg], "aatom1") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_aatom1;
      if (kindflag != NONE && kindflag != ANGLE)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = ANGLE;
    } else if (strcmp(arg[iarg], "aatom2") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_aatom2;
      if (kindflag != NONE && kindflag != ANGLE)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = ANGLE;
    } else if (strcmp(arg[iarg], "aatom3") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_aatom3;
      if (kindflag != NONE && kindflag != ANGLE)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = ANGLE;
    } else if (strcmp(arg[iarg], "atype") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_atype;
      if (kindflag != NONE && kindflag != ANGLE)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = ANGLE;

    } else if (strcmp(arg[iarg], "datom1") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_datom1;
      if (kindflag != NONE && kindflag != DIHEDRAL)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = DIHEDRAL;
    } else if (strcmp(arg[iarg], "datom2") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_datom2;
      if (kindflag != NONE && kindflag != DIHEDRAL)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = DIHEDRAL;
    } else if (strcmp(arg[iarg], "datom3") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_datom3;
      if (kindflag != NONE && kindflag != DIHEDRAL)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = DIHEDRAL;
    } else if (strcmp(arg[iarg], "datom4") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_datom4;
      if (kindflag != NONE && kindflag != DIHEDRAL)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = DIHEDRAL;
    } else if (strcmp(arg[iarg], "dtype") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_dtype;
      if (kindflag != NONE && kindflag != DIHEDRAL)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = DIHEDRAL;

    } else if (strcmp(arg[iarg], "iatom1") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_iatom1;
      if (kindflag != NONE && kindflag != IMPROPER)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = IMPROPER;
    } else if (strcmp(arg[iarg], "iatom2") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_iatom2;
      if (kindflag != NONE && kindflag != IMPROPER)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = IMPROPER;
    } else if (strcmp(arg[iarg], "iatom3") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_iatom3;
      if (kindflag != NONE && kindflag != IMPROPER)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = IMPROPER;
    } else if (strcmp(arg[iarg], "iatom4") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_iatom4;
      if (kindflag != NONE && kindflag != IMPROPER)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = IMPROPER;
    } else if (strcmp(arg[iarg], "itype") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_itype;
      if (kindflag != NONE && kindflag != IMPROPER)
        error->all(FLERR, "Compute property/local cannot use these inputs together");
      kindflag = IMPROPER;

    } else
      break;

    nvalues++;
    iarg++;
  }

  if (nvalues == 1)
    size_local_cols = 0;
  else
    size_local_cols = nvalues;

  // optional args

  cutstyle = TYPE;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute property/local command");
      if (strcmp(arg[iarg + 1], "type") == 0)
        cutstyle = TYPE;
      else if (strcmp(arg[iarg + 1], "radius") == 0)
        cutstyle = RADIUS;
      else
        error->all(FLERR, "Illegal compute property/local command");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal compute property/local command");
  }

  // error check

  if (atom->molecular == 2 &&
      (kindflag == BOND || kindflag == ANGLE || kindflag == DIHEDRAL || kindflag == IMPROPER))
    error->all(FLERR,
               "Compute property/local does not (yet) work "
               "with atom_style template");

  if (kindflag == BOND && atom->avec->bonds_allow == 0)
    error->all(FLERR, "Compute property/local for property that isn't allocated");
  if (kindflag == ANGLE && atom->avec->angles_allow == 0)
    error->all(FLERR, "Compute property/local for property that isn't allocated");
  if (kindflag == DIHEDRAL && atom->avec->dihedrals_allow == 0)
    error->all(FLERR, "Compute property/local for property that isn't allocated");
  if (kindflag == IMPROPER && atom->avec->impropers_allow == 0)
    error->all(FLERR, "Compute property/local for property that isn't allocated");
  if (cutstyle == RADIUS && !atom->radius_flag)
    error->all(FLERR, "Compute property/local requires atom attribute radius");

  nmax = 0;
  vlocal = nullptr;
  alocal = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputePropertyLocal::~ComputePropertyLocal()
{
  delete[] pack_choice;
  memory->destroy(vlocal);
  memory->destroy(alocal);
  memory->destroy(indices);
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::init()
{
  if (kindflag == NEIGH || kindflag == PAIR) {
    if (force->pair == nullptr)
      error->all(FLERR, "No pair style is defined for compute property/local");
    if (force->pair->single_enable == 0)
      error->all(FLERR, "Pair style does not support compute property/local");
  }

  // for NEIGH/PAIR need an occasional half neighbor list
  // set size to same value as request made by force->pair
  // this should enable it to always be a copy list  (e.g. for granular pstyle)

  if (kindflag == NEIGH || kindflag == PAIR) {
    int neighflags = NeighConst::REQ_OCCASIONAL;
    auto pairrequest = neighbor->find_request(force->pair);
    if (pairrequest && pairrequest->get_size()) neighflags |= NeighConst::REQ_SIZE;
    neighbor->add_request(this, neighflags);
  }

  // do initial memory allocation so that memory_usage() is correct
  // cannot be done yet for NEIGH/PAIR, since neigh list does not exist

  if (kindflag == NEIGH)
    ncount = 0;
  else if (kindflag == PAIR)
    ncount = 0;
  else if (kindflag == BOND)
    ncount = count_bonds(0);
  else if (kindflag == ANGLE)
    ncount = count_angles(0);
  else if (kindflag == DIHEDRAL)
    ncount = count_dihedrals(0);
  else if (kindflag == IMPROPER)
    ncount = count_impropers(0);

  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and generate list of indices

  if (kindflag == NEIGH)
    ncount = count_pairs(0, 0);
  else if (kindflag == PAIR)
    ncount = count_pairs(0, 1);
  else if (kindflag == BOND)
    ncount = count_bonds(0);
  else if (kindflag == ANGLE)
    ncount = count_angles(0);
  else if (kindflag == DIHEDRAL)
    ncount = count_dihedrals(0);
  else if (kindflag == IMPROPER)
    ncount = count_impropers(0);

  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;

  if (kindflag == NEIGH)
    ncount = count_pairs(1, 0);
  else if (kindflag == PAIR)
    ncount = count_pairs(1, 1);
  else if (kindflag == BOND)
    ncount = count_bonds(1);
  else if (kindflag == ANGLE)
    ncount = count_angles(1);
  else if (kindflag == DIHEDRAL)
    ncount = count_dihedrals(1);
  else if (kindflag == IMPROPER)
    ncount = count_impropers(1);

  // fill vector or array with local values

  if (nvalues == 1) {
    buf = vlocal;
    (this->*pack_choice[0])(0);
  } else {
    if (alocal) buf = &alocal[0][0];
    for (int n = 0; n < nvalues; n++) (this->*pack_choice[n])(n);
  }
}

/* ----------------------------------------------------------------------
   count pairs and compute pair info on this proc
   only count pair once if newton_pair is off
   both atom I,J must be in group
   if allflag is set, compute requested info about pair
   if forceflag = 1, pair must be within force cutoff, else neighbor cutoff
------------------------------------------------------------------------- */

int ComputePropertyLocal::count_pairs(int allflag, int forceflag)
{
  int i, j, m, ii, jj, inum, jnum, itype, jtype;
  tagint itag, jtag;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, radsum;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double **x = atom->x;
  double *radius = atom->radius;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)

  if (allflag == 0) neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // skip if I or J are not in group
  // for newton = 0 and J = ghost atom,
  //   need to ensure I,J pair is only output by one proc
  //   use same itag,jtag logic as in Neighbor::neigh_half_nsq()

  double **cutsq = force->pair->cutsq;

  m = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itag = tag[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;

      // itag = jtag is possible for long cutoffs that include images of self

      if (newton_pair == 0 && j >= nlocal) {
        jtag = tag[j];
        if (itag > jtag) {
          if ((itag + jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag + jtag) % 2 == 1) continue;
        } else {
          if (x[j][2] < ztmp) continue;
          if (x[j][2] == ztmp) {
            if (x[j][1] < ytmp) continue;
            if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
          }
        }
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      if (forceflag) {
        if (cutstyle == TYPE) {
          if (rsq >= cutsq[itype][jtype]) continue;
        } else {
          radsum = radius[i] + radius[j];
          if (rsq >= radsum * radsum) continue;
        }
      }

      if (allflag) {
        indices[m][0] = i;
        indices[m][1] = j;
      }
      m++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   count bonds on this proc
   only count bond once if newton_bond is off
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if bond is deleted (type = 0), do not count
   if bond is turned off (type < 0), still count
------------------------------------------------------------------------- */

int ComputePropertyLocal::count_bonds(int flag)
{
  int i, atom1, atom2;

  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int m = 0;
  for (atom1 = 0; atom1 < nlocal; atom1++) {
    if (!(mask[atom1] & groupbit)) continue;
    for (i = 0; i < num_bond[atom1]; i++) {
      atom2 = atom->map(bond_atom[atom1][i]);
      if (atom2 < 0 || !(mask[atom2] & groupbit)) continue;
      if (newton_bond == 0 && tag[atom1] > tag[atom2]) continue;
      if (bond_type[atom1][i] == 0) continue;

      if (flag) {
        indices[m][0] = atom1;
        indices[m][1] = i;
      }
      m++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   count angles on this proc
   only count if 2nd atom is the one storing the angle
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if angle is deleted (type = 0), do not count
   if angle is turned off (type < 0), still count
------------------------------------------------------------------------- */

int ComputePropertyLocal::count_angles(int flag)
{
  int i, atom1, atom2, atom3;

  int *num_angle = atom->num_angle;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (atom2 = 0; atom2 < nlocal; atom2++) {
    if (!(mask[atom2] & groupbit)) continue;
    for (i = 0; i < num_angle[atom2]; i++) {
      if (tag[atom2] != angle_atom2[atom2][i]) continue;
      atom1 = atom->map(angle_atom1[atom2][i]);
      if (atom1 < 0 || !(mask[atom1] & groupbit)) continue;
      atom3 = atom->map(angle_atom3[atom2][i]);
      if (atom3 < 0 || !(mask[atom3] & groupbit)) continue;
      if (angle_type[atom2][i] == 0) continue;

      if (flag) {
        indices[m][0] = atom2;
        indices[m][1] = i;
      }
      m++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   count dihedrals on this proc
   only count if 2nd atom is the one storing the dihedral
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
------------------------------------------------------------------------- */

int ComputePropertyLocal::count_dihedrals(int flag)
{
  int i, atom1, atom2, atom3, atom4;

  int *num_dihedral = atom->num_dihedral;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom3 = atom->dihedral_atom3;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (atom2 = 0; atom2 < nlocal; atom2++) {
    if (!(mask[atom2] & groupbit)) continue;
    for (i = 0; i < num_dihedral[atom2]; i++) {
      if (tag[atom2] != dihedral_atom2[atom2][i]) continue;
      atom1 = atom->map(dihedral_atom1[atom2][i]);
      if (atom1 < 0 || !(mask[atom1] & groupbit)) continue;
      atom3 = atom->map(dihedral_atom3[atom2][i]);
      if (atom3 < 0 || !(mask[atom3] & groupbit)) continue;
      atom4 = atom->map(dihedral_atom4[atom2][i]);
      if (atom4 < 0 || !(mask[atom4] & groupbit)) continue;

      if (flag) {
        indices[m][0] = atom2;
        indices[m][1] = i;
      }
      m++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   count impropers on this proc
   only count if 2nd atom is the one storing the improper
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
------------------------------------------------------------------------- */

int ComputePropertyLocal::count_impropers(int flag)
{
  int i, atom1, atom2, atom3, atom4;

  int *num_improper = atom->num_improper;
  tagint **improper_atom1 = atom->improper_atom1;
  tagint **improper_atom2 = atom->improper_atom2;
  tagint **improper_atom3 = atom->improper_atom3;
  tagint **improper_atom4 = atom->improper_atom4;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (atom2 = 0; atom2 < nlocal; atom2++) {
    if (!(mask[atom2] & groupbit)) continue;
    for (i = 0; i < num_improper[atom2]; i++) {
      if (tag[atom2] != improper_atom2[atom2][i]) continue;
      atom1 = atom->map(improper_atom1[atom2][i]);
      if (atom1 < 0 || !(mask[atom1] & groupbit)) continue;
      atom3 = atom->map(improper_atom3[atom2][i]);
      if (atom3 < 0 || !(mask[atom3] & groupbit)) continue;
      atom4 = atom->map(improper_atom4[atom2][i]);
      if (atom4 < 0 || !(mask[atom4] & groupbit)) continue;

      if (flag) {
        indices[m][0] = atom2;
        indices[m][1] = i;
      }
      m++;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::reallocate(int n)
{
  // grow vector_local or array_local, also indices

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vlocal);
    memory->create(vlocal, nmax, "property/local:vector_local");
    vector_local = vlocal;
  } else {
    memory->destroy(alocal);
    memory->create(alocal, nmax, nvalues, "property/local:array_local");
    array_local = alocal;
  }

  memory->destroy(indices);
  memory->create(indices, nmax, 2, "property/local:indices");
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputePropertyLocal::memory_usage()
{
  double bytes = (double) nmax * nvalues * sizeof(double);
  bytes += (double) nmax * 2 * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword compute property/local can output
   the atom property is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_patom1(int n)
{
  int i;
  tagint *tag = atom->tag;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    buf[n] = tag[i];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_patom2(int n)
{
  int i;
  tagint *tag = atom->tag;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][1];
    buf[n] = tag[i];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_ptype1(int n)
{
  int i;
  int *type = atom->type;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    buf[n] = type[i];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_ptype2(int n)
{
  int i;
  int *type = atom->type;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][1];
    buf[n] = type[i];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_batom1(int n)
{
  int i;
  tagint *tag = atom->tag;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    buf[n] = tag[i];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_batom2(int n)
{
  int i, j;
  tagint **bond_atom = atom->bond_atom;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = bond_atom[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_btype(int n)
{
  int i, j;
  int **bond_type = atom->bond_type;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = bond_type[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_aatom1(int n)
{
  int i, j;
  tagint **angle_atom1 = atom->angle_atom1;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = angle_atom1[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_aatom2(int n)
{
  int i, j;
  tagint **angle_atom2 = atom->angle_atom2;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = angle_atom2[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_aatom3(int n)
{
  int i, j;
  tagint **angle_atom3 = atom->angle_atom3;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = angle_atom3[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_atype(int n)
{
  int i, j;
  int **angle_type = atom->angle_type;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = angle_type[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_datom1(int n)
{
  int i, j;
  tagint **dihedral_atom1 = atom->dihedral_atom1;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = dihedral_atom1[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_datom2(int n)
{
  int i, j;
  tagint **dihedral_atom2 = atom->dihedral_atom2;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = dihedral_atom2[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_datom3(int n)
{
  int i, j;
  tagint **dihedral_atom3 = atom->dihedral_atom3;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = dihedral_atom3[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_datom4(int n)
{
  int i, j;
  tagint **dihedral_atom4 = atom->dihedral_atom4;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = dihedral_atom4[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_dtype(int n)
{
  int i, j;
  int **dihedral_type = atom->dihedral_type;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = dihedral_type[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_iatom1(int n)
{
  int i, j;
  tagint **improper_atom1 = atom->improper_atom1;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = improper_atom1[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_iatom2(int n)
{
  int i, j;
  tagint **improper_atom2 = atom->improper_atom2;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = improper_atom2[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_iatom3(int n)
{
  int i, j;
  tagint **improper_atom3 = atom->improper_atom3;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = improper_atom3[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_iatom4(int n)
{
  int i, j;
  tagint **improper_atom4 = atom->improper_atom4;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = improper_atom4[i][j];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_itype(int n)
{
  int i, j;
  int **improper_type = atom->improper_type;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = improper_type[i][j];
    n += nvalues;
  }
}
