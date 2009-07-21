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

#include "stdlib.h"
#include "atom_vec_molecular.h"
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecMolecular::AtomVecMolecular(LAMMPS *lmp, int narg, char **arg) : 
  AtomVec(lmp, narg, arg)
{
  molecular = 1;
  bonds_allow = 1;
  angles_allow = 1;
  dihedrals_allow = 1;
  impropers_allow = 1;
  mass_type = 1;
  size_comm = 3;
  size_reverse = 3;
  size_border = 7;
  size_data_atom = 6;
  size_data_vel = 4;
  xcol_data = 4;

  atom->molecule_flag = 1;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n 
------------------------------------------------------------------------- */

void AtomVecMolecular::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  atom->nmax = nmax;

  tag = atom->tag = (int *) 
    memory->srealloc(atom->tag,nmax*sizeof(int),"atom:tag");
  type = atom->type = (int *)
    memory->srealloc(atom->type,nmax*sizeof(int),"atom:type");
  mask = atom->mask = (int *) 
    memory->srealloc(atom->mask,nmax*sizeof(int),"atom:mask");
  image = atom->image = (int *) 
    memory->srealloc(atom->image,nmax*sizeof(int),"atom:image");
  x = atom->x = memory->grow_2d_double_array(atom->x,nmax,3,"atom:x");
  v = atom->v = memory->grow_2d_double_array(atom->v,nmax,3,"atom:v");
  f = atom->f = memory->grow_2d_double_array(atom->f,nmax,3,"atom:f");

  molecule = atom->molecule = (int *) 
    memory->srealloc(atom->molecule,nmax*sizeof(int),"atom:molecule");

  nspecial = atom->nspecial =
    memory->grow_2d_int_array(atom->nspecial,nmax,3,"atom:nspecial");
  special = atom->special =
    memory->grow_2d_int_array(atom->special,nmax,atom->maxspecial,
			      "atom:special");

  num_bond = atom->num_bond = (int *) 
    memory->srealloc(atom->num_bond,nmax*sizeof(int),"atom:num_bond");
  bond_type = atom->bond_type = 
    memory->grow_2d_int_array(atom->bond_type,nmax,atom->bond_per_atom,
			      "atom:bond_type");
  bond_atom = atom->bond_atom = 
    memory->grow_2d_int_array(atom->bond_atom,nmax,atom->bond_per_atom,
			      "atom:bond_atom");

  num_angle = atom->num_angle = (int *) 
    memory->srealloc(atom->num_angle,nmax*sizeof(int),"atom:num_angle");
  angle_type = atom->angle_type = 
    memory->grow_2d_int_array(atom->angle_type,nmax,atom->angle_per_atom,
			      "atom:angle_type");
  angle_atom1 = atom->angle_atom1 = 
    memory->grow_2d_int_array(atom->angle_atom1,nmax,atom->angle_per_atom,
			      "atom:angle_atom1");
  angle_atom2 = atom->angle_atom2 = 
    memory->grow_2d_int_array(atom->angle_atom2,nmax,atom->angle_per_atom,
			      "atom:angle_atom2");
  angle_atom3 = atom->angle_atom3 = 
    memory->grow_2d_int_array(atom->angle_atom3,nmax,atom->angle_per_atom,
			      "atom:angle_atom3");

  num_dihedral = atom->num_dihedral = (int *) 
    memory->srealloc(atom->num_dihedral,nmax*sizeof(int),"atom:num_dihedral");
  dihedral_type = atom->dihedral_type = 
    memory->grow_2d_int_array(atom->dihedral_type,nmax,atom->dihedral_per_atom,
			      "atom:dihedral_type");
  dihedral_atom1 = atom->dihedral_atom1 = 
    memory->grow_2d_int_array(atom->dihedral_atom1,nmax,
			      atom->dihedral_per_atom,"atom:dihedral_atom1");
  dihedral_atom2 = atom->dihedral_atom2 = 
    memory->grow_2d_int_array(atom->dihedral_atom2,nmax,
			      atom->dihedral_per_atom,"atom:dihedral_atom2");
  dihedral_atom3 = atom->dihedral_atom3 = 
    memory->grow_2d_int_array(atom->dihedral_atom3,nmax,
			      atom->dihedral_per_atom,"atom:dihedral_atom3");
  dihedral_atom4 = atom->dihedral_atom4 = 
    memory->grow_2d_int_array(atom->dihedral_atom4,nmax,
			      atom->dihedral_per_atom,"atom:dihedral_atom4");

  num_improper = atom->num_improper = (int *) 
    memory->srealloc(atom->num_improper,nmax*sizeof(int),"atom:num_improper");
  improper_type = atom->improper_type = 
    memory->grow_2d_int_array(atom->improper_type,nmax,atom->improper_per_atom,
			      "atom:improper_type");
  improper_atom1 = atom->improper_atom1 = 
    memory->grow_2d_int_array(atom->improper_atom1,nmax,
			      atom->improper_per_atom,"atom:improper_atom1");
  improper_atom2 = atom->improper_atom2 = 
    memory->grow_2d_int_array(atom->improper_atom2,nmax,
			      atom->improper_per_atom,"atom:improper_atom2");
  improper_atom3 = atom->improper_atom3 = 
    memory->grow_2d_int_array(atom->improper_atom3,nmax,
			      atom->improper_per_atom,"atom:improper_atom3");
  improper_atom4 = atom->improper_atom4 = 
    memory->grow_2d_int_array(atom->improper_atom4,nmax,
			      atom->improper_per_atom,"atom:improper_atom4");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ---------------------------------------------------------------------- */

void AtomVecMolecular::reset_special()
{
  special = atom->special;
}

/* ---------------------------------------------------------------------- */

void AtomVecMolecular::copy(int i, int j)
{
  int k;

  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  molecule[j] = molecule[i];

  num_bond[j] = num_bond[i];
  for (k = 0; k < num_bond[j]; k++) {
    bond_type[j][k] = bond_type[i][k];
    bond_atom[j][k] = bond_atom[i][k];
  }

  num_angle[j] = num_angle[i];
  for (k = 0; k < num_angle[j]; k++) {
    angle_type[j][k] = angle_type[i][k];
    angle_atom1[j][k] = angle_atom1[i][k];
    angle_atom2[j][k] = angle_atom2[i][k];
    angle_atom3[j][k] = angle_atom3[i][k];
  }

  num_dihedral[j] = num_dihedral[i];
  for (k = 0; k < num_dihedral[j]; k++) {
    dihedral_type[j][k] = dihedral_type[i][k];
    dihedral_atom1[j][k] = dihedral_atom1[i][k];
    dihedral_atom2[j][k] = dihedral_atom2[i][k];
    dihedral_atom3[j][k] = dihedral_atom3[i][k];
    dihedral_atom4[j][k] = dihedral_atom4[i][k];
  }

  num_improper[j] = num_improper[i];
  for (k = 0; k < num_improper[j]; k++) {
    improper_type[j][k] = improper_type[i][k];
    improper_atom1[j][k] = improper_atom1[i][k];
    improper_atom2[j][k] = improper_atom2[i][k];
    improper_atom3[j][k] = improper_atom3[i][k];
    improper_atom4[j][k] = improper_atom4[i][k];
  }

  nspecial[j][0] = nspecial[i][0];
  nspecial[j][1] = nspecial[i][1];
  nspecial[j][2] = nspecial[i][2];
  for (k = 0; k < nspecial[j][2]; k++) special[j][k] = special[i][k];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j);
}

/* ---------------------------------------------------------------------- */

int AtomVecMolecular::pack_comm(int n, int *list, double *buf, 
				int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecMolecular::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecMolecular::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecMolecular::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecMolecular::pack_border(int n, int *list, double *buf,
				  int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = molecule[j];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = molecule[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMolecular::pack_border_one(int i, double *buf)
{
  buf[0] = molecule[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void AtomVecMolecular::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    molecule[i] = static_cast<int> (buf[m++]);
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecMolecular::unpack_border_one(int i, double *buf)
{
  molecule[i] = static_cast<int> (buf[0]);
  return 1;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them 
------------------------------------------------------------------------- */

int AtomVecMolecular::pack_exchange(int i, double *buf)
{
  int k;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  buf[m++] = image[i];

  buf[m++] = molecule[i];

  buf[m++] = num_bond[i];
  for (k = 0; k < num_bond[i]; k++) {
    buf[m++] = bond_type[i][k];
    buf[m++] = bond_atom[i][k];
  }

  buf[m++] = num_angle[i];
  for (k = 0; k < num_angle[i]; k++) {
    buf[m++] = angle_type[i][k];
    buf[m++] = angle_atom1[i][k];
    buf[m++] = angle_atom2[i][k];
    buf[m++] = angle_atom3[i][k];
  }

  buf[m++] = num_dihedral[i];
  for (k = 0; k < num_dihedral[i]; k++) {
    buf[m++] = dihedral_type[i][k];
    buf[m++] = dihedral_atom1[i][k];
    buf[m++] = dihedral_atom2[i][k];
    buf[m++] = dihedral_atom3[i][k];
    buf[m++] = dihedral_atom4[i][k];
  }

  buf[m++] = num_improper[i];
  for (k = 0; k < num_improper[i]; k++) {
    buf[m++] = improper_type[i][k];
    buf[m++] = improper_atom1[i][k];
    buf[m++] = improper_atom2[i][k];
    buf[m++] = improper_atom3[i][k];
    buf[m++] = improper_atom4[i][k];
  }

  buf[m++] = nspecial[i][0];
  buf[m++] = nspecial[i][1];
  buf[m++] = nspecial[i][2];
  for (k = 0; k < nspecial[i][2]; k++) buf[m++] = special[i][k];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMolecular::unpack_exchange(double *buf)
{
  int k;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = static_cast<int> (buf[m++]);

  molecule[nlocal] = static_cast<int> (buf[m++]);

  num_bond[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = static_cast<int> (buf[m++]);
    bond_atom[nlocal][k] = static_cast<int> (buf[m++]);
  }

  num_angle[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_angle[nlocal]; k++) {
    angle_type[nlocal][k] = static_cast<int> (buf[m++]);
    angle_atom1[nlocal][k] = static_cast<int> (buf[m++]);
    angle_atom2[nlocal][k] = static_cast<int> (buf[m++]);
    angle_atom3[nlocal][k] = static_cast<int> (buf[m++]);
  }

  num_dihedral[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_dihedral[nlocal]; k++) {
    dihedral_type[nlocal][k] = static_cast<int> (buf[m++]);
    dihedral_atom1[nlocal][k] = static_cast<int> (buf[m++]);
    dihedral_atom2[nlocal][k] = static_cast<int> (buf[m++]);
    dihedral_atom3[nlocal][k] = static_cast<int> (buf[m++]);
    dihedral_atom4[nlocal][k] = static_cast<int> (buf[m++]);
  }

  num_improper[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_improper[nlocal]; k++) {
    improper_type[nlocal][k] = static_cast<int> (buf[m++]);
    improper_atom1[nlocal][k] = static_cast<int> (buf[m++]);
    improper_atom2[nlocal][k] = static_cast<int> (buf[m++]);
    improper_atom3[nlocal][k] = static_cast<int> (buf[m++]);
    improper_atom4[nlocal][k] = static_cast<int> (buf[m++]);
  }

  nspecial[nlocal][0] = static_cast<int> (buf[m++]);
  nspecial[nlocal][1] = static_cast<int> (buf[m++]);
  nspecial[nlocal][2] = static_cast<int> (buf[m++]);
  for (k = 0; k < nspecial[nlocal][2]; k++)
    special[nlocal][k] = static_cast<int> (buf[m++]);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) 
      m += modify->fix[atom->extra_grow[iextra]]->
	unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecMolecular::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 0;
  for (i = 0; i < nlocal; i++)
    n += 16 + 2*num_bond[i] + 4*num_angle[i] +
      5*num_dihedral[i] + 5*num_improper[i];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++) 
      for (i = 0; i < nlocal; i++)
	n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive   
------------------------------------------------------------------------- */

int AtomVecMolecular::pack_restart(int i, double *buf)
{
  int k;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  buf[m++] = image[i];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = molecule[i];
 
  buf[m++] = num_bond[i];
  for (k = 0; k < num_bond[i]; k++) {
    buf[m++] = MAX(bond_type[i][k],-bond_type[i][k]);
    buf[m++] = bond_atom[i][k];
  }
  
  buf[m++] = num_angle[i];
  for (k = 0; k < num_angle[i]; k++) {
    buf[m++] = MAX(angle_type[i][k],-angle_type[i][k]);
    buf[m++] = angle_atom1[i][k];
    buf[m++] = angle_atom2[i][k];
    buf[m++] = angle_atom3[i][k];
  }

  buf[m++] = num_dihedral[i];
  for (k = 0; k < num_dihedral[i]; k++) {
    buf[m++] = MAX(dihedral_type[i][k],-dihedral_type[i][k]);
    buf[m++] = dihedral_atom1[i][k];
    buf[m++] = dihedral_atom2[i][k];
    buf[m++] = dihedral_atom3[i][k];
    buf[m++] = dihedral_atom4[i][k];
  }

  buf[m++] = num_improper[i];
  for (k = 0; k < num_improper[i]; k++) {
    buf[m++] = MAX(improper_type[i][k],-improper_type[i][k]);
    buf[m++] = improper_atom1[i][k];
    buf[m++] = improper_atom2[i][k];
    buf[m++] = improper_atom3[i][k];
    buf[m++] = improper_atom4[i][k];
  }

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++) 
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecMolecular::unpack_restart(double *buf)
{
  int k;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      atom->extra = memory->grow_2d_double_array(atom->extra,nmax,
						 atom->nextra_store,
						 "atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = static_cast<int> (buf[m++]);
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  molecule[nlocal] = static_cast<int> (buf[m++]);
    
  num_bond[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = static_cast<int> (buf[m++]);
    bond_atom[nlocal][k] = static_cast<int> (buf[m++]);
  }

  num_angle[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_angle[nlocal]; k++) {
    angle_type[nlocal][k] = static_cast<int> (buf[m++]);
    angle_atom1[nlocal][k] = static_cast<int> (buf[m++]);
    angle_atom2[nlocal][k] = static_cast<int> (buf[m++]);
    angle_atom3[nlocal][k] = static_cast<int> (buf[m++]);
  }

  num_dihedral[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_dihedral[nlocal]; k++) {
    dihedral_type[nlocal][k] = static_cast<int> (buf[m++]);
    dihedral_atom1[nlocal][k] = static_cast<int> (buf[m++]);
    dihedral_atom2[nlocal][k] = static_cast<int> (buf[m++]);
    dihedral_atom3[nlocal][k] = static_cast<int> (buf[m++]);
    dihedral_atom4[nlocal][k] = static_cast<int> (buf[m++]);
  }

  num_improper[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_improper[nlocal]; k++) {
    improper_type[nlocal][k] = static_cast<int> (buf[m++]);
    improper_atom1[nlocal][k] = static_cast<int> (buf[m++]);
    improper_atom2[nlocal][k] = static_cast<int> (buf[m++]);
    improper_atom3[nlocal][k] = static_cast<int> (buf[m++]);
    improper_atom4[nlocal][k] = static_cast<int> (buf[m++]);
  }

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecMolecular::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = (512 << 20) | (512 << 10) | 512;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  molecule[nlocal] = 0;
  num_bond[nlocal] = 0;
  num_angle[nlocal] = 0;
  num_dihedral[nlocal] = 0;
  num_improper[nlocal] = 0;
  nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecMolecular::data_atom(double *coord, int imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = atoi(values[0]);
  if (tag[nlocal] <= 0)
    error->one("Invalid atom ID in Atoms section of data file");

  molecule[nlocal] = atoi(values[1]);

  type[nlocal] = atoi(values[2]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one("Invalid atom type in Atoms section of data file");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  num_bond[nlocal] = 0;
  num_angle[nlocal] = 0;
  num_dihedral[nlocal] = 0;
  num_improper[nlocal] = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecMolecular::data_atom_hybrid(int nlocal, char **values)
{
  molecule[nlocal] = atoi(values[0]);

  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  num_bond[nlocal] = 0;
  num_angle[nlocal] = 0;
  num_dihedral[nlocal] = 0;
  num_improper[nlocal] = 0;

  return 1;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory 
------------------------------------------------------------------------- */

double AtomVecMolecular::memory_usage()
{
  double bytes = 0.0;

  if (atom->memcheck("tag")) bytes += nmax * sizeof(int);
  if (atom->memcheck("type")) bytes += nmax * sizeof(int);
  if (atom->memcheck("mask")) bytes += nmax * sizeof(int);
  if (atom->memcheck("image")) bytes += nmax * sizeof(int);
  if (atom->memcheck("x")) bytes += nmax*3 * sizeof(double);
  if (atom->memcheck("v")) bytes += nmax*3 * sizeof(double);
  if (atom->memcheck("f")) bytes += nmax*3 * sizeof(double);

  if (atom->memcheck("molecule")) bytes += nmax * sizeof(int);
  if (atom->memcheck("nspecial")) bytes += nmax*3 * sizeof(int);
  if (atom->memcheck("special")) bytes += nmax*atom->maxspecial * sizeof(int);

  if (atom->memcheck("num_bond")) bytes += nmax * sizeof(int);
  if (atom->memcheck("bond_type"))
    bytes += nmax*atom->bond_per_atom * sizeof(int);
  if (atom->memcheck("bond_atom"))
    bytes += nmax*atom->bond_per_atom * sizeof(int);

  if (atom->memcheck("num_angle")) bytes += nmax * sizeof(int);
  if (atom->memcheck("angle_type"))
    bytes += nmax*atom->angle_per_atom * sizeof(int);
  if (atom->memcheck("angle_atom1"))
    bytes += nmax*atom->angle_per_atom * sizeof(int);
  if (atom->memcheck("angle_atom2"))
    bytes += nmax*atom->angle_per_atom * sizeof(int);
  if (atom->memcheck("angle_atom3"))
    bytes += nmax*atom->angle_per_atom * sizeof(int);

  if (atom->memcheck("num_dihedral")) bytes += nmax * sizeof(int);
  if (atom->memcheck("dihedral_type"))
    bytes += nmax*atom->dihedral_per_atom * sizeof(int);
  if (atom->memcheck("dihedral_atom1"))
    bytes += nmax*atom->dihedral_per_atom * sizeof(int);
  if (atom->memcheck("dihedral_atom2"))
    bytes += nmax*atom->dihedral_per_atom * sizeof(int);
  if (atom->memcheck("dihedral_atom3"))
    bytes += nmax*atom->dihedral_per_atom * sizeof(int);
  if (atom->memcheck("dihedral_atom4"))
    bytes += nmax*atom->dihedral_per_atom * sizeof(int);

  if (atom->memcheck("num_improper")) bytes += nmax * sizeof(int);
  if (atom->memcheck("improper_type"))
    bytes += nmax*atom->improper_per_atom * sizeof(int);
  if (atom->memcheck("improper_atom1"))
    bytes += nmax*atom->improper_per_atom * sizeof(int);
  if (atom->memcheck("improper_atom2"))
    bytes += nmax*atom->improper_per_atom * sizeof(int);
  if (atom->memcheck("improper_atom3"))
    bytes += nmax*atom->improper_per_atom * sizeof(int);
  if (atom->memcheck("improper_atom4"))
    bytes += nmax*atom->improper_per_atom * sizeof(int);

  return bytes;
}
