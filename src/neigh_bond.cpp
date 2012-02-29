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

#include "neighbor.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BONDDELTA 10000

/* ---------------------------------------------------------------------- */

void Neighbor::bond_all()
{
  int i,m,atom1;

  int nlocal = atom->nlocal;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *tag = atom->tag;
  int newton_bond = force->newton_bond;

  nbondlist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_bond[i]; m++) {
      atom1 = atom->map(bond_atom[i][m]);
      if (atom1 == -1) {
	char str[128];
	sprintf(str,
		"Bond atoms %d %d missing on proc %d at step " BIGINT_FORMAT,
		tag[i],bond_atom[i][m],me,update->ntimestep);
	error->one(FLERR,str);
      }
      if (newton_bond || i < atom1) {
	if (nbondlist == maxbond) {
	  maxbond += BONDDELTA;
	  memory->grow(bondlist,maxbond,3,"neighbor:bondlist");
	}
	bondlist[nbondlist][0] = i;
	bondlist[nbondlist][1] = atom1;
	bondlist[nbondlist][2] = bond_type[i][m];
	nbondlist++;
      }
    }
}

/* ---------------------------------------------------------------------- */

void Neighbor::bond_partial()
{
  int i,m,atom1;

  int nlocal = atom->nlocal;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *tag = atom->tag;
  int newton_bond = force->newton_bond;

  nbondlist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_bond[i]; m++) {
      if (bond_type[i][m] <= 0) continue;
      atom1 = atom->map(bond_atom[i][m]);
      if (atom1 == -1) {
	char str[128];
	sprintf(str,
		"Bond atoms %d %d missing on proc %d at step " BIGINT_FORMAT,
		tag[i],bond_atom[i][m],me,update->ntimestep);
	error->one(FLERR,str);
      }
      if (newton_bond || i < atom1) {
	if (nbondlist == maxbond) {
	  maxbond += BONDDELTA;
	  memory->grow(bondlist,maxbond,3,"neighbor:bondlist");
	}
	bondlist[nbondlist][0] = i;
	bondlist[nbondlist][1] = atom1;
	bondlist[nbondlist][2] = bond_type[i][m];
	nbondlist++;
      }
    }
}

/* ---------------------------------------------------------------------- */

void Neighbor::angle_all()
{
  int i,m,atom1,atom2,atom3;

  int nlocal = atom->nlocal;
  int *num_angle = atom->num_angle;
  int **angle_atom1 = atom->angle_atom1;
  int **angle_atom2 = atom->angle_atom2;
  int **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  int newton_bond = force->newton_bond;

  nanglelist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_angle[i]; m++) {
      atom1 = atom->map(angle_atom1[i][m]);
      atom2 = atom->map(angle_atom2[i][m]);
      atom3 = atom->map(angle_atom3[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1) {
	char str[128];
	sprintf(str,
		"Angle atoms %d %d %d missing on proc %d at step " 
		BIGINT_FORMAT,
		angle_atom1[i][m],angle_atom2[i][m],angle_atom3[i][m],
		me,update->ntimestep);
	error->one(FLERR,str);
      }
      if (newton_bond || (i <= atom1 && i <= atom2 && i <= atom3)) {
	if (nanglelist == maxangle) {
	  maxangle += BONDDELTA;
	  memory->grow(anglelist,maxangle,4,"neighbor:anglelist");
	}
	anglelist[nanglelist][0] = atom1;
	anglelist[nanglelist][1] = atom2;
	anglelist[nanglelist][2] = atom3;
	anglelist[nanglelist][3] = angle_type[i][m];
	nanglelist++;
      }
    }
}

/* ---------------------------------------------------------------------- */

void Neighbor::angle_partial()
{
  int i,m,atom1,atom2,atom3;

  int nlocal = atom->nlocal;
  int *num_angle = atom->num_angle;
  int **angle_atom1 = atom->angle_atom1;
  int **angle_atom2 = atom->angle_atom2;
  int **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  int newton_bond = force->newton_bond;

  nanglelist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_angle[i]; m++) {
      if (angle_type[i][m] <= 0) continue;
      atom1 = atom->map(angle_atom1[i][m]);
      atom2 = atom->map(angle_atom2[i][m]);
      atom3 = atom->map(angle_atom3[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1) {
	char str[128];
	sprintf(str,
		"Angle atoms %d %d %d missing on proc %d at step " 
		BIGINT_FORMAT,
		angle_atom1[i][m],angle_atom2[i][m],angle_atom3[i][m],
		me,update->ntimestep);
	error->one(FLERR,str);
      }
      if (newton_bond || (i <= atom1 && i <= atom2 && i <= atom3)) {
	if (nanglelist == maxangle) {
	  maxangle += BONDDELTA;
	  memory->grow(anglelist,maxangle,4,"neighbor:anglelist");
	}
	anglelist[nanglelist][0] = atom1;
	anglelist[nanglelist][1] = atom2;
	anglelist[nanglelist][2] = atom3;
	anglelist[nanglelist][3] = angle_type[i][m];
	nanglelist++;
      }
    }
}

/* ---------------------------------------------------------------------- */

void Neighbor::dihedral_all()
{
  int i,m,atom1,atom2,atom3,atom4;

  int nlocal = atom->nlocal;
  int *num_dihedral = atom->num_dihedral;
  int **dihedral_atom1 = atom->dihedral_atom1;
  int **dihedral_atom2 = atom->dihedral_atom2;
  int **dihedral_atom3 = atom->dihedral_atom3;
  int **dihedral_atom4 = atom->dihedral_atom4;
  int **dihedral_type = atom->dihedral_type;
  int newton_bond = force->newton_bond;

  ndihedrallist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_dihedral[i]; m++) {
      atom1 = atom->map(dihedral_atom1[i][m]);
      atom2 = atom->map(dihedral_atom2[i][m]);
      atom3 = atom->map(dihedral_atom3[i][m]);
      atom4 = atom->map(dihedral_atom4[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
	char str[128];
	sprintf(str,
		"Dihedral atoms %d %d %d %d missing on proc %d at step " 
		BIGINT_FORMAT,
		dihedral_atom1[i][m],dihedral_atom2[i][m],
		dihedral_atom3[i][m],dihedral_atom4[i][m],
		me,update->ntimestep);
	error->one(FLERR,str);
      }
      if (newton_bond || 
	  (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
	if (ndihedrallist == maxdihedral) {
	  maxdihedral += BONDDELTA;
	  memory->grow(dihedrallist,maxdihedral,5,"neighbor:dihedrallist");
	}
	dihedrallist[ndihedrallist][0] = atom1;
	dihedrallist[ndihedrallist][1] = atom2;
	dihedrallist[ndihedrallist][2] = atom3;
	dihedrallist[ndihedrallist][3] = atom4;
	dihedrallist[ndihedrallist][4] = dihedral_type[i][m];
	ndihedrallist++;
      }
    }
}

/* ---------------------------------------------------------------------- */

void Neighbor::dihedral_partial()
{
  int i,m,atom1,atom2,atom3,atom4;

  int nlocal = atom->nlocal;
  int *num_dihedral = atom->num_dihedral;
  int **dihedral_atom1 = atom->dihedral_atom1;
  int **dihedral_atom2 = atom->dihedral_atom2;
  int **dihedral_atom3 = atom->dihedral_atom3;
  int **dihedral_atom4 = atom->dihedral_atom4;
  int **dihedral_type = atom->dihedral_type;
  int newton_bond = force->newton_bond;

  ndihedrallist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_dihedral[i]; m++) {
      if (dihedral_type[i][m] <= 0) continue;
      atom1 = atom->map(dihedral_atom1[i][m]);
      atom2 = atom->map(dihedral_atom2[i][m]);
      atom3 = atom->map(dihedral_atom3[i][m]);
      atom4 = atom->map(dihedral_atom4[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
	char str[128];
	sprintf(str,
		"Dihedral atoms %d %d %d %d missing on proc %d at step " 
		BIGINT_FORMAT,
		dihedral_atom1[i][m],dihedral_atom2[i][m],
		dihedral_atom3[i][m],dihedral_atom4[i][m],
		me,update->ntimestep);
	error->one(FLERR,str);
      }
      if (newton_bond || 
	  (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
	if (ndihedrallist == maxdihedral) {
	  maxdihedral += BONDDELTA;
	  memory->grow(dihedrallist,maxdihedral,5,"neighbor:dihedrallist");
	}
	dihedrallist[ndihedrallist][0] = atom1;
	dihedrallist[ndihedrallist][1] = atom2;
	dihedrallist[ndihedrallist][2] = atom3;
	dihedrallist[ndihedrallist][3] = atom4;
	dihedrallist[ndihedrallist][4] = dihedral_type[i][m];
	ndihedrallist++;
      }
    }
}

/* ---------------------------------------------------------------------- */

void Neighbor::improper_all()
{
  int i,m,atom1,atom2,atom3,atom4;

  int nlocal = atom->nlocal;
  int *num_improper = atom->num_improper;
  int **improper_atom1 = atom->improper_atom1;
  int **improper_atom2 = atom->improper_atom2;
  int **improper_atom3 = atom->improper_atom3;
  int **improper_atom4 = atom->improper_atom4;
  int **improper_type = atom->improper_type;
  int newton_bond = force->newton_bond;

  nimproperlist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_improper[i]; m++) {
      atom1 = atom->map(improper_atom1[i][m]);
      atom2 = atom->map(improper_atom2[i][m]);
      atom3 = atom->map(improper_atom3[i][m]);
      atom4 = atom->map(improper_atom4[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
	char str[128];
	sprintf(str,
		"Improper atoms %d %d %d %d missing on proc %d at step " 
		BIGINT_FORMAT,
		improper_atom1[i][m],improper_atom2[i][m],
		improper_atom3[i][m],improper_atom4[i][m],
		me,update->ntimestep);
	error->one(FLERR,str);
      }
      if (newton_bond || 
	  (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
	if (nimproperlist == maximproper) {
	  maximproper += BONDDELTA;
	  memory->grow(improperlist,maximproper,5,"neighbor:improperlist");
	}
	improperlist[nimproperlist][0] = atom1;
	improperlist[nimproperlist][1] = atom2;
	improperlist[nimproperlist][2] = atom3;
	improperlist[nimproperlist][3] = atom4;
	improperlist[nimproperlist][4] = improper_type[i][m];
	nimproperlist++;
      }
    }
}

/* ---------------------------------------------------------------------- */

void Neighbor::improper_partial()
{
  int i,m,atom1,atom2,atom3,atom4;

  int nlocal = atom->nlocal;
  int *num_improper = atom->num_improper;
  int **improper_atom1 = atom->improper_atom1;
  int **improper_atom2 = atom->improper_atom2;
  int **improper_atom3 = atom->improper_atom3;
  int **improper_atom4 = atom->improper_atom4;
  int **improper_type = atom->improper_type;
  int newton_bond = force->newton_bond;

  nimproperlist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_improper[i]; m++) {
      if (improper_type[i][m] <= 0) continue;
      atom1 = atom->map(improper_atom1[i][m]);
      atom2 = atom->map(improper_atom2[i][m]);
      atom3 = atom->map(improper_atom3[i][m]);
      atom4 = atom->map(improper_atom4[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
	char str[128];
	sprintf(str,
		"Improper atoms %d %d %d %d missing on proc %d at step " 
		BIGINT_FORMAT,
		improper_atom1[i][m],improper_atom2[i][m],
		improper_atom3[i][m],improper_atom4[i][m],
		me,update->ntimestep);
	error->one(FLERR,str);
      }
      if (newton_bond || 
	  (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)) {
	if (nimproperlist == maximproper) {
	  maximproper += BONDDELTA;
	  memory->grow(improperlist,maximproper,5,"neighbor:improperlist");
	}
	improperlist[nimproperlist][0] = atom1;
	improperlist[nimproperlist][1] = atom2;
	improperlist[nimproperlist][2] = atom3;
	improperlist[nimproperlist][3] = atom4;
	improperlist[nimproperlist][4] = improper_type[i][m];
	nimproperlist++;
      }
    }
}
