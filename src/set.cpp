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

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "set.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "group.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "random_park.h"
#include "error.h"

using namespace LAMMPS_NS;

#define ONEATOM  0
#define ATOM     1
#define BOND     2
#define ANGLE    3
#define DIHEDRAL 4
#define IMPROPER 5
#define CHARGE   6
#define DIPOLE   7

#define WARMUP 100

/* ---------------------------------------------------------------------- */

Set::Set(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Set::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all("Set command before simulation box is defined");
  if (atom->natoms == 0)
    error->all("Set command with no atoms existing");
  if (narg != 3) error->all("Illegal set command");

  // set style

  int style;
  if (strcmp(arg[1],"atom") == 0) style = ATOM;
  else if (strcmp(arg[1],"bond") == 0) style = BOND;
  else if (strcmp(arg[1],"angle") == 0) style = ANGLE;
  else if (strcmp(arg[1],"dihedral") == 0) style = DIHEDRAL;
  else if (strcmp(arg[1],"improper") == 0) style = IMPROPER;
  else if (strcmp(arg[1],"charge") == 0) style = CHARGE;
  else if (strcmp(arg[1],"dipole") == 0) style = DIPOLE;
  else style = ONEATOM;

  // set atom or group

  int atomid,igroup,groupbit;
  if (style == ONEATOM) atomid = atoi(arg[0]);
  else {
    igroup = group->find(arg[0]);
    if (igroup == -1) error->all("Cannot find set command group ID");
    groupbit = group->bitmask[igroup];
  }
  
  // consistency checks

  if ((style == BOND && atom->avec->bonds_allow == 0) ||
      (style == ANGLE && atom->avec->angles_allow == 0) ||
      (style == DIHEDRAL && atom->avec->dihedrals_allow == 0) ||
      (style == IMPROPER && atom->avec->impropers_allow == 0) ||
      (style == CHARGE && atom->q == NULL) ||
      (style == DIPOLE && atom->dipole == NULL))
    error->all("Cannot set this attribute for this atom style");

  if (style == ONEATOM && strcmp(arg[1],"q") == 0 && atom->q == NULL)
    error->all("Cannot set this attribute for this atom style");

  if (style == ONEATOM && strcmp(arg[1],"mol") == 0 && atom->molecular == 0)
    error->all("Cannot set this attribute for this atom style");

  // border swap to insure type and mask is current for off-proc atoms
  // only needed for BOND, ANGLE, etc types
  // enforce PBC before in case atoms are outside box
  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  
  if (style == BOND || style == ANGLE || 
      style == DIHEDRAL || style == IMPROPER) {
    if (comm->me == 0 && screen) fprintf(screen,"System init for set ...\n");
    lmp->init();

    domain->pbc();
    domain->reset_box();
    comm->setup();
    comm->exchange();
    comm->borders();
  }

  if (comm->me == 0 && screen) fprintf(screen,"Setting atom values ...\n");
  
  // set new values

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int m,atom1,atom2,atom3,atom4;

  int count = 0;

  // for style ATOM, atom must be in group

  if (style == ATOM) {
    int ivalue = atoi(arg[2]);
    if (ivalue < 1 || ivalue > atom->ntypes)
      error->all("Invalid type in set command");
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	atom->type[i] = ivalue;
	count++;
      }
  }

  // for style BOND, each of 2 atoms must be in group

  if (style == BOND) {
    int ivalue = atoi(arg[2]);
    if (ivalue < 1 || ivalue > atom->nbondtypes)
      error->all("Invalid type in set command");
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_bond[i]; m++) {
	atom1 = atom->map(atom->bond_atom[i][m]);
	if (atom1 == -1) error->one("Bond atom missing in set command");
	if (mask[i] & groupbit && mask[atom1] & groupbit) {
	  atom->bond_type[i][m] = ivalue;
	  count++;
	}
      }
  }

  // for style ANGLE, each of 3 atoms must be in group

  if (style == ANGLE) {
    int ivalue = atoi(arg[2]);
    if (ivalue < 1 || ivalue > atom->nangletypes)
      error->all("Invalid type in set command");
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_angle[i]; m++) {
	atom1 = atom->map(atom->angle_atom1[i][m]);
	atom2 = atom->map(atom->angle_atom2[i][m]);
	atom3 = atom->map(atom->angle_atom3[i][m]);
	if (atom1 == -1 || atom2 == -1 || atom3 == -1)
	  error->one("Angle atom missing in set command");
	if (mask[atom1] & groupbit && mask[atom2] & groupbit &&
	    mask[atom3] & groupbit) {
	  atom->angle_type[i][m] = ivalue;
	  count++;
	}
      }
  }

  // for style DIHEDRAL, each of 4 atoms must be in group

  if (style == DIHEDRAL) {
    int ivalue = atoi(arg[2]);
    if (ivalue < 1 || ivalue > atom->ndihedraltypes)
      error->all("Invalid type in set command");
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_dihedral[i]; m++) {
	atom1 = atom->map(atom->dihedral_atom1[i][m]);
	atom2 = atom->map(atom->dihedral_atom2[i][m]);
	atom3 = atom->map(atom->dihedral_atom3[i][m]);
	atom4 = atom->map(atom->dihedral_atom4[i][m]);
	if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
	  error->one("Dihedral atom missing in set command");
	if (mask[atom1] & groupbit && mask[atom2] & groupbit &&
	    mask[atom3] & groupbit && mask[atom4] & groupbit) {
	  atom->dihedral_type[i][m] = ivalue;
	  count++;
	}
      }
  }

  // for style IMPROPER, each of 4 atoms must be in group

  if (style == IMPROPER) {
    int ivalue = atoi(arg[2]);
    if (ivalue < 1 || ivalue > atom->nimpropertypes)
      error->all("Invalid type in set command");
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_improper[i]; m++) {
	atom1 = atom->map(atom->improper_atom1[i][m]);
	atom2 = atom->map(atom->improper_atom2[i][m]);
	atom3 = atom->map(atom->improper_atom3[i][m]);
	atom4 = atom->map(atom->improper_atom4[i][m]);
	if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
	  error->one("Improper atom missing in set command");
	if (mask[atom1] & groupbit && mask[atom2] & groupbit &&
	    mask[atom3] & groupbit && mask[atom4] & groupbit) {
	  atom->improper_type[i][m] = ivalue;
	  count++;
	}
      }
  }

  // for style CHARGE, set charge of all atoms in group

  if (style == CHARGE) {
    double value = atof(arg[2]);
    double *q = atom->q;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	q[i] = value;
	count++;
      }
  }

  // for style DIPOLE, create unique random number stream for each proc
  // set dipole moment of each atom in group to random orientation
  // dipole length is determined by dipole type array

  if (style == DIPOLE) {
    int ivalue = atoi(arg[2]);
    if (ivalue <= 0) error->all("Invalid random number seed in set command");
    double msq,scale;
    RanPark *random = new RanPark(lmp,ivalue + comm->me);
    for (int i = 0; i < WARMUP; i++) random->uniform();
    int *type = atom->type;
    double *dipole = atom->dipole;
    double **mu = atom->mu;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	mu[i][0] = random->uniform() - 0.5;
	mu[i][1] = random->uniform() - 0.5;
	mu[i][2] = random->uniform() - 0.5;
	msq = mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2];
	scale = dipole[type[i]]/sqrt(msq);
	mu[i][0] *= scale;
	mu[i][1] *= scale;
	mu[i][2] *= scale;
	count++;
      }
    delete random;
  }

  // for style ONEATOM, set attribute of single atom ID

  if (style == ONEATOM) {
    int *tag = atom->tag;
    for (int i = 0; i < nlocal; i++)
      if (atomid == tag[i]) {
	if (strcmp(arg[1],"type") == 0) {
	  atom->type[i] = atoi(arg[2]);
	  if (atom->type[i] < 1 || atom->type[i] > atom->ntypes)
	    error->one("Invalid type in set command");
	} else if (strcmp(arg[1],"mol") == 0) atom->molecule[i] = atoi(arg[2]);
	else if (strcmp(arg[1],"x") == 0) atom->x[i][0] = atof(arg[2]);
	else if (strcmp(arg[1],"y") == 0) atom->x[i][1] = atof(arg[2]);
	else if (strcmp(arg[1],"z") == 0) atom->x[i][2] = atof(arg[2]);
	else if (strcmp(arg[1],"vx") == 0) atom->v[i][0] = atof(arg[2]);
	else if (strcmp(arg[1],"vy") == 0) atom->v[i][1] = atof(arg[2]);
	else if (strcmp(arg[1],"vz") == 0) atom->v[i][2] = atof(arg[2]);
	else if (strcmp(arg[1],"q") == 0) atom->q[i] = atof(arg[2]);
	else error->one("Illegal set command");
	count++;
      }
  }

  // statistics

  int allcount;
  MPI_Allreduce(&count,&allcount,1,MPI_INT,MPI_SUM,world);
  
  if (comm->me == 0) {
    if (screen) fprintf(screen,"  %d settings made\n",allcount);
    if (logfile) fprintf(logfile,"  %d settings made\n",allcount);
  }
}
