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

#include "string.h"
#include "compute_property_local.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

// customize by adding keyword
// same list as in dump_custom.cpp

enum{ID,MOL,TYPE,MASS,
       X,Y,Z,XS,YS,ZS,XSTRI,YSTRI,ZSTRI,XU,YU,ZU,XUTRI,YUTRI,ZUTRI,IX,IY,IZ,
       VX,VY,VZ,FX,FY,FZ,
       Q,MUX,MUY,MUZ,RADIUS,OMEGAX,OMEGAY,OMEGAZ,ANGMOMX,ANGMOMY,ANGMOMZ,
       QUATW,QUATI,QUATJ,QUATK,TQX,TQY,TQZ};
enum{NONE,BOND,ANGLE,DIHEDRAL,IMPROPER};

#define DELTA 10000

/* ---------------------------------------------------------------------- */

ComputePropertyLocal::ComputePropertyLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal compute property/atom command");

  local_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  pack_choice = new FnPtrPack[nvalues];

  kindflag = NONE;

  int i;
  for (int iarg = 3; iarg < narg; iarg++) {
    i = iarg-3;

    if (strcmp(arg[iarg],"batom1") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_batom1;
      if (kindflag != NONE && kindflag != BOND) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = BOND;
    } else if (strcmp(arg[iarg],"batom2") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_batom2;
      if (kindflag != NONE && kindflag != BOND) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = BOND;
    } else if (strcmp(arg[iarg],"btype") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_btype;
      if (kindflag != NONE && kindflag != BOND) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = BOND;

    } else if (strcmp(arg[iarg],"aatom1") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_aatom1;
      if (kindflag != NONE && kindflag != ANGLE) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = ANGLE;
    } else if (strcmp(arg[iarg],"aatom2") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_aatom2;
      if (kindflag != NONE && kindflag != ANGLE) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = ANGLE;
    } else if (strcmp(arg[iarg],"aatom3") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_aatom3;
      if (kindflag != NONE && kindflag != ANGLE) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = ANGLE;
    } else if (strcmp(arg[iarg],"atype") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_atype;
      if (kindflag != NONE && kindflag != ANGLE) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = ANGLE;

    } else if (strcmp(arg[iarg],"datom1") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_datom1;
      if (kindflag != NONE && kindflag != DIHEDRAL) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = DIHEDRAL;
    } else if (strcmp(arg[iarg],"datom2") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_datom2;
      if (kindflag != NONE && kindflag != DIHEDRAL) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = DIHEDRAL;
    } else if (strcmp(arg[iarg],"datom3") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_datom3;
      if (kindflag != NONE && kindflag != DIHEDRAL) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = DIHEDRAL;
    } else if (strcmp(arg[iarg],"datom4") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_datom4;
      if (kindflag != NONE && kindflag != DIHEDRAL) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = DIHEDRAL;
    } else if (strcmp(arg[iarg],"dtype") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_dtype;
      if (kindflag != NONE && kindflag != DIHEDRAL) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = DIHEDRAL;

    } else if (strcmp(arg[iarg],"iatom1") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_iatom1;
      if (kindflag != NONE && kindflag != IMPROPER) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = IMPROPER;
    } else if (strcmp(arg[iarg],"iatom2") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_iatom2;
      if (kindflag != NONE && kindflag != IMPROPER) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = IMPROPER;
    } else if (strcmp(arg[iarg],"iatom3") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_iatom3;
      if (kindflag != NONE && kindflag != IMPROPER) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = IMPROPER;
    } else if (strcmp(arg[iarg],"iatom4") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_iatom4;
      if (kindflag != NONE && kindflag != IMPROPER) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = IMPROPER;
    } else if (strcmp(arg[iarg],"itype") == 0) {
      pack_choice[i] = &ComputePropertyLocal::pack_itype;
      if (kindflag != NONE && kindflag != IMPROPER) 
	error->all("Compute property/local cannot use these inputs together");
      kindflag = IMPROPER;

    } else error->all("Invalid keyword in compute property/local command");
  }

  // error check

  if (kindflag == BOND && atom->avec->bonds_allow == 0)
    error->all("Compute property/local for property that isn't allocated");
  if (kindflag == ANGLE && atom->avec->angles_allow == 0)
    error->all("Compute property/local for property that isn't allocated");
  if (kindflag == DIHEDRAL && atom->avec->dihedrals_allow == 0)
    error->all("Compute property/local for property that isn't allocated");
  if (kindflag == IMPROPER && atom->avec->impropers_allow == 0)
    error->all("Compute property/local for property that isn't allocated");

  nmax = 0;
  vector = NULL;
  array = NULL;
  indices = NULL;
}

/* ---------------------------------------------------------------------- */

ComputePropertyLocal::~ComputePropertyLocal()
{
  delete [] pack_choice;
  memory->sfree(vector);
  memory->destroy_2d_double_array(array);
  memory->destroy_2d_int_array(indices);
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::init()
{
  // do initial memory allocation so that memory_usage() is correct

  if (kindflag == BOND) ncount = count_bonds(0);
  else if (kindflag == ANGLE) ncount = count_angles(0);
  else if (kindflag == DIHEDRAL) ncount = count_dihedrals(0);
  else if (kindflag == IMPROPER) ncount = count_impropers(0);

  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and generate list of indices

  if (kindflag == BOND) ncount = count_bonds(0);
  else if (kindflag == ANGLE) ncount = count_angles(0);
  else if (kindflag == DIHEDRAL) ncount = count_dihedrals(0);
  else if (kindflag == IMPROPER) ncount = count_impropers(0);

  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;

  if (kindflag == BOND) ncount = count_bonds(1);
  else if (kindflag == ANGLE) ncount = count_angles(1);
  else if (kindflag == DIHEDRAL) ncount = count_dihedrals(1);
  else if (kindflag == IMPROPER) ncount = count_impropers(1);

  // fill vector or array with local values

  if (nvalues == 1) {
    buf = vector;
    (this->*pack_choice[0])(0);
  } else {
    if (array) buf = array[0];
    for (int n = 0; n < nvalues; n++)
      (this->*pack_choice[n])(n);
  }
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
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int i,atom1,atom2;

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
------------------------------------------------------------------------- */

int ComputePropertyLocal::count_angles(int flag)
{
  int *num_angle = atom->num_angle;
  int **angle_atom1 = atom->angle_atom1;
  int **angle_atom2 = atom->angle_atom2;
  int **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int i,atom1,atom2,atom3;

  int m = 0;
  for (atom2 = 0; atom2 < nlocal; atom2++) {
    if (!(mask[atom2] & groupbit)) continue;
    for (i = 0; i < num_angle[atom2]; i++) {
      if (tag[atom2] != angle_atom2[atom2][i]) continue;
      atom1 = atom->map(angle_atom1[atom2][i]);
      if (atom1 < 0 || !(mask[atom1] & groupbit)) continue;
      atom3 = atom->map(angle_atom3[atom2][i]);
      if (atom3 < 0 || !(mask[atom3] & groupbit)) continue;

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
  int *num_dihedral = atom->num_dihedral;
  int **dihedral_atom1 = atom->dihedral_atom1;
  int **dihedral_atom2 = atom->dihedral_atom2;
  int **dihedral_atom3 = atom->dihedral_atom3;
  int **dihedral_atom4 = atom->dihedral_atom4;
  int **dihedral_type = atom->dihedral_type;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int i,atom1,atom2,atom3,atom4;

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
  int *num_improper = atom->num_improper;
  int **improper_atom1 = atom->improper_atom1;
  int **improper_atom2 = atom->improper_atom2;
  int **improper_atom3 = atom->improper_atom3;
  int **improper_atom4 = atom->improper_atom4;
  int **improper_type = atom->improper_type;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int i,atom1,atom2,atom3,atom4;

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
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;
  if (nvalues == 1) {
    memory->sfree(vector);
    vector = (double *) memory->smalloc(nmax*sizeof(double),
					"property/local:vector");
    vector_local = vector;
  } else {
    memory->destroy_2d_double_array(array);
    array = memory->create_2d_double_array(nmax,nvalues,
					   "property/local:array");
    array_local = array;
  }

  memory->destroy_2d_int_array(indices);
  indices = memory->create_2d_int_array(nmax,2,"property/local:indices");
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputePropertyLocal::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  bytes += nmax*2 * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword compute property/local can output
   the atom property is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_batom1(int n)
{
  int i;
  int *tag = atom->tag;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    buf[n] = tag[i];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyLocal::pack_batom2(int n)
{
  int i,j;
  int **bond_atom = atom->bond_atom;

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
  int i,j;
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
  int i,j;
  int **angle_atom1 = atom->angle_atom1;

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
  int i,j;
  int **angle_atom2 = atom->angle_atom2;

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
  int i,j;
  int **angle_atom3 = atom->angle_atom3;

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
  int i,j;
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
  int i,j;
  int **dihedral_atom1 = atom->dihedral_atom1;

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
  int i,j;
  int **dihedral_atom2 = atom->dihedral_atom2;

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
  int i,j;
  int **dihedral_atom3 = atom->dihedral_atom3;

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
  int i,j;
  int **dihedral_atom4 = atom->dihedral_atom4;

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
  int i,j;
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
  int i,j;
  int **improper_atom1 = atom->improper_atom1;

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
  int i,j;
  int **improper_atom2 = atom->improper_atom2;

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
  int i,j;
  int **improper_atom3 = atom->improper_atom3;

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
  int i,j;
  int **improper_atom4 = atom->improper_atom4;

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
  int i,j;
  int **improper_type = atom->improper_type;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    buf[n] = improper_type[i][j];
    n += nvalues;
  }
}
