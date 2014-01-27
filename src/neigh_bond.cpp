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
#include "atom_vec.h"
#include "molecule.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "output.h"
#include "thermo.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BONDDELTA 10000

enum{IGNORE,WARN,ERROR};           // same as thermo.cpp

// bondlist, anglelist, dihedrallist, improperlist
//   no longer store atom->map() of the bond partners
// instead store domain->closest_image() of the bond partners of atom I
// this enables distances between list atoms to be calculated
//   w/out invoking domain->minimium_image(), e.g. in bond->compute()

/* ---------------------------------------------------------------------- */

void Neighbor::bond_all()
{
  int i,m,atom1;

  int nlocal = atom->nlocal;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  tagint *tag = atom->tag;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nbondlist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_bond[i]; m++) {
      atom1 = atom->map(bond_atom[i][m]);
      if (atom1 == -1) {
        nmissing++;
        if (lostbond == ERROR) {
          char str[128];
          sprintf(str,"Bond atoms " TAGINT_FORMAT " " TAGINT_FORMAT 
                  " missing on proc %d at step " BIGINT_FORMAT,
                  tag[i],bond_atom[i][m],me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
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

  if (cluster_check) bond_check();
  if (lostbond == IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Bond atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::bond_template()
{
  int i,m,atom1;
  int imol,iatom;
  tagint tagprev;
  int *num_bond;
  int **bond_atom,**bond_type;

  Molecule **onemols = atom->avec->onemols;

  tagint *tag = atom->tag;
  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nbondlist = 0;

  for (i = 0; i < nlocal; i++) {
    if (molindex[i] < 0) continue;
    imol = molindex[i];
    iatom = molatom[i];
    tagprev = tag[i] - iatom - 1;
    num_bond = onemols[imol]->num_bond;
    bond_atom = onemols[imol]->bond_atom;
    bond_type = onemols[imol]->bond_type;

    for (m = 0; m < num_bond[iatom]; m++) {
      if (bond_type[iatom][m] <= 0) continue;
      atom1 = atom->map(bond_atom[iatom][m]+tagprev);
      if (atom1 == -1) {
        nmissing++;
        if (lostbond == ERROR) {
          char str[128];
          sprintf(str,"Bond atoms " TAGINT_FORMAT " " TAGINT_FORMAT 
                  " missing on proc %d at step " BIGINT_FORMAT,
                  tag[i],bond_atom[iatom][m]+tagprev,me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
      if (newton_bond || i < atom1) {
        if (nbondlist == maxbond) {
          maxbond += BONDDELTA;
          memory->grow(bondlist,maxbond,3,"neighbor:bondlist");
        }
        bondlist[nbondlist][0] = i;
        bondlist[nbondlist][1] = atom1;
        bondlist[nbondlist][2] = bond_type[iatom][m];
        nbondlist++;
      }
    }
  }

  if (cluster_check) bond_check();
  if (lostbond == IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Bond atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::bond_partial()
{
  int i,m,atom1;

  int nlocal = atom->nlocal;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  tagint *tag = atom->tag;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nbondlist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_bond[i]; m++) {
      if (bond_type[i][m] <= 0) continue;
      atom1 = atom->map(bond_atom[i][m]);
      if (atom1 == -1) {
        nmissing++;
        if (lostbond == ERROR) {
          char str[128];
          sprintf(str,"Bond atoms " TAGINT_FORMAT " " TAGINT_FORMAT 
                  " missing on proc %d at step " BIGINT_FORMAT,
                  tag[i],bond_atom[i][m],me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
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

  if (cluster_check) bond_check();
  if (lostbond == IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Bond atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::bond_check()
{
  int i,j;
  double dx,dy,dz,dxstart,dystart,dzstart;
  
  double **x = atom->x;
  int flag = 0;

  for (int m = 0; m < nbondlist; m++) {
    i = bondlist[m][0];
    j = bondlist[m][1];
    dxstart = dx = x[i][0] - x[j][0];
    dystart = dy = x[i][1] - x[j][1];
    dzstart = dz = x[i][2] - x[j][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Bond extent > half of periodic box length");
}

/* ---------------------------------------------------------------------- */

void Neighbor::angle_all()
{
  int i,m,atom1,atom2,atom3;

  int nlocal = atom->nlocal;
  int *num_angle = atom->num_angle;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nanglelist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_angle[i]; m++) {
      atom1 = atom->map(angle_atom1[i][m]);
      atom2 = atom->map(angle_atom2[i][m]);
      atom3 = atom->map(angle_atom3[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1) {
        nmissing++;
        if (lostbond == ERROR) {
          char str[128];
          sprintf(str,"Angle atoms "
                  TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT 
                  " missing on proc %d at step " BIGINT_FORMAT,
                  angle_atom1[i][m],angle_atom2[i][m],angle_atom3[i][m],
                  me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
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

  if (cluster_check) angle_check();
  if (lostbond == IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Angle atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::angle_template()
{
  int i,m,atom1,atom2,atom3;
  int imol,iatom;
  tagint tagprev;
  int *num_angle;
  int **angle_atom1,**angle_atom2,**angle_atom3,**angle_type;

  Molecule **onemols = atom->avec->onemols;

  tagint *tag = atom->tag;
  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nanglelist = 0;

  for (i = 0; i < nlocal; i++) {
    if (molindex[i] < 0) continue;
    imol = molindex[i];
    iatom = molatom[i];
    tagprev = tag[i] - iatom - 1;
    num_angle = onemols[imol]->num_angle;
    angle_atom1 = onemols[imol]->angle_atom1;
    angle_atom2 = onemols[imol]->angle_atom2;
    angle_atom3 = onemols[imol]->angle_atom3;
    angle_type = onemols[imol]->angle_type;

    for (m = 0; m < num_angle[iatom]; m++) {
      if (angle_type[iatom][m] <= 0) continue;
      atom1 = atom->map(angle_atom1[iatom][m]+tagprev);
      atom2 = atom->map(angle_atom2[iatom][m]+tagprev);
      atom3 = atom->map(angle_atom3[iatom][m]+tagprev);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1) {
        nmissing++;
        if (lostbond == ERROR) {
          char str[128];
          sprintf(str,"Angle atoms "
                  TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT 
                  " missing on proc %d at step " BIGINT_FORMAT,
                  angle_atom1[iatom][m]+tagprev,angle_atom2[iatom][m]+tagprev,
                  angle_atom3[iatom][m]+tagprev,
                  me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      if (newton_bond || (i <= atom1 && i <= atom2 && i <= atom3)) {
        if (nanglelist == maxangle) {
          maxangle += BONDDELTA;
          memory->grow(anglelist,maxangle,4,"neighbor:anglelist");
        }
        anglelist[nanglelist][0] = atom1;
        anglelist[nanglelist][1] = atom2;
        anglelist[nanglelist][2] = atom3;
        anglelist[nanglelist][3] = angle_type[iatom][m];
        nanglelist++;
      }
    }
  }

  if (cluster_check) angle_check();
  if (lostbond == IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Angle atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::angle_partial()
{
  int i,m,atom1,atom2,atom3;

  int nlocal = atom->nlocal;
  int *num_angle = atom->num_angle;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  int **angle_type = atom->angle_type;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nanglelist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_angle[i]; m++) {
      if (angle_type[i][m] <= 0) continue;
      atom1 = atom->map(angle_atom1[i][m]);
      atom2 = atom->map(angle_atom2[i][m]);
      atom3 = atom->map(angle_atom3[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1) {
        nmissing++;
        if (lostbond == ERROR) {
          char str[128];
          sprintf(str,"Angle atoms "
                  TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT 
                  " missing on proc %d at step " BIGINT_FORMAT,
                  angle_atom1[i][m],angle_atom2[i][m],angle_atom3[i][m],
                  me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
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

  if (cluster_check) angle_check();
  if (lostbond == IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Angle atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::angle_check()
{
  int i,j,k;
  double dx,dy,dz,dxstart,dystart,dzstart;
  
  double **x = atom->x;
  int flag = 0;

  // check all 3 distances
  // in case angle potential computes any of them

  for (int m = 0; m < nanglelist; m++) {
    i = anglelist[m][0];
    j = anglelist[m][1];
    k = anglelist[m][2];
    dxstart = dx = x[i][0] - x[j][0];
    dystart = dy = x[i][1] - x[j][1];
    dzstart = dz = x[i][2] - x[j][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[i][0] - x[k][0];
    dystart = dy = x[i][1] - x[k][1];
    dzstart = dz = x[i][2] - x[k][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[j][0] - x[k][0];
    dystart = dy = x[j][1] - x[k][1];
    dzstart = dz = x[j][2] - x[k][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Angle extent > half of periodic box length");
}

/* ---------------------------------------------------------------------- */

void Neighbor::dihedral_all()
{
  int i,m,atom1,atom2,atom3,atom4;

  int nlocal = atom->nlocal;
  int *num_dihedral = atom->num_dihedral;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom3 = atom->dihedral_atom3;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  int **dihedral_type = atom->dihedral_type;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  ndihedrallist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_dihedral[i]; m++) {
      atom1 = atom->map(dihedral_atom1[i][m]);
      atom2 = atom->map(dihedral_atom2[i][m]);
      atom3 = atom->map(dihedral_atom3[i][m]);
      atom4 = atom->map(dihedral_atom4[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
        nmissing++;
        if (lostbond == ERROR) {
          char str[128];
          sprintf(str,"Dihedral atoms "
                  TAGINT_FORMAT " " TAGINT_FORMAT " " 
                  TAGINT_FORMAT " " TAGINT_FORMAT 
                  " missing on proc %d at step " BIGINT_FORMAT,
                  dihedral_atom1[i][m],dihedral_atom2[i][m],
                  dihedral_atom3[i][m],dihedral_atom4[i][m],
                  me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain->closest_image(i,atom4);
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

  if (cluster_check) dihedral_check(ndihedrallist,dihedrallist);
  if (lostbond == IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Dihedral atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::dihedral_template()
{
  int i,m,atom1,atom2,atom3,atom4;
  int imol,iatom;
  tagint tagprev;
  int *num_dihedral;
  int **dihedral_atom1,**dihedral_atom2,**dihedral_atom3,**dihedral_atom4;
  int **dihedral_type;

  Molecule **onemols = atom->avec->onemols;

  tagint *tag = atom->tag;
  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  ndihedrallist = 0;

  for (i = 0; i < nlocal; i++) {
    if (molindex[i] < 0) continue;
    imol = molindex[i];
    iatom = molatom[i];
    tagprev = tag[i] - iatom - 1;
    num_dihedral = onemols[imol]->num_dihedral;
    dihedral_atom1 = onemols[imol]->dihedral_atom1;
    dihedral_atom2 = onemols[imol]->dihedral_atom2;
    dihedral_atom3 = onemols[imol]->dihedral_atom3;
    dihedral_atom4 = onemols[imol]->dihedral_atom4;
    dihedral_type = onemols[imol]->dihedral_type;

    for (m = 0; m < num_dihedral[iatom]; m++) {
      atom1 = atom->map(dihedral_atom1[iatom][m]+tagprev);
      atom2 = atom->map(dihedral_atom2[iatom][m]+tagprev);
      atom3 = atom->map(dihedral_atom3[iatom][m]+tagprev);
      atom4 = atom->map(dihedral_atom4[iatom][m]+tagprev);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
        nmissing++;
        if (lostbond == ERROR) {
          char str[128];
          sprintf(str,"Dihedral atoms "
                  TAGINT_FORMAT " " TAGINT_FORMAT " " 
                  TAGINT_FORMAT " " TAGINT_FORMAT 
                  " missing on proc %d at step " BIGINT_FORMAT,
                  dihedral_atom1[iatom][m]+tagprev,
                  dihedral_atom2[iatom][m]+tagprev,
                  dihedral_atom3[iatom][m]+tagprev,
                  dihedral_atom4[iatom][m]+tagprev,
                  me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain->closest_image(i,atom4);
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
        dihedrallist[ndihedrallist][4] = dihedral_type[iatom][m];
        ndihedrallist++;
      }
    }
  }

  if (cluster_check) dihedral_check(ndihedrallist,dihedrallist);
  if (lostbond == IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Dihedral atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::dihedral_partial()
{
  int i,m,atom1,atom2,atom3,atom4;

  int nlocal = atom->nlocal;
  int *num_dihedral = atom->num_dihedral;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom3 = atom->dihedral_atom3;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  int **dihedral_type = atom->dihedral_type;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  ndihedrallist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_dihedral[i]; m++) {
      if (dihedral_type[i][m] <= 0) continue;
      atom1 = atom->map(dihedral_atom1[i][m]);
      atom2 = atom->map(dihedral_atom2[i][m]);
      atom3 = atom->map(dihedral_atom3[i][m]);
      atom4 = atom->map(dihedral_atom4[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
        nmissing++;
        if (lostbond == ERROR) {
          char str[128];
          sprintf(str,"Dihedral atoms "
                  TAGINT_FORMAT " " TAGINT_FORMAT " " 
                  TAGINT_FORMAT " " TAGINT_FORMAT 
                  " missing on proc %d at step " BIGINT_FORMAT,
                  dihedral_atom1[i][m],dihedral_atom2[i][m],
                  dihedral_atom3[i][m],dihedral_atom4[i][m],
                  me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain->closest_image(i,atom4);
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

  if (cluster_check) dihedral_check(ndihedrallist,dihedrallist);
  if (lostbond == IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Dihedral atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::dihedral_check(int nlist, int **list)
{
  int i,j,k,l;
  double dx,dy,dz,dxstart,dystart,dzstart;
  
  double **x = atom->x;
  int flag = 0;

  // check all 6 distances
  // in case dihedral/improper potential computes any of them

  for (int m = 0; m < nlist; m++) {
    i = list[m][0];
    j = list[m][1];
    k = list[m][2];
    l = list[m][3];
    dxstart = dx = x[i][0] - x[j][0];
    dystart = dy = x[i][1] - x[j][1];
    dzstart = dz = x[i][2] - x[j][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[i][0] - x[k][0];
    dystart = dy = x[i][1] - x[k][1];
    dzstart = dz = x[i][2] - x[k][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[i][0] - x[l][0];
    dystart = dy = x[i][1] - x[l][1];
    dzstart = dz = x[i][2] - x[l][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[j][0] - x[k][0];
    dystart = dy = x[j][1] - x[k][1];
    dzstart = dz = x[j][2] - x[k][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[j][0] - x[l][0];
    dystart = dy = x[j][1] - x[l][1];
    dzstart = dz = x[j][2] - x[l][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
    dxstart = dx = x[k][0] - x[l][0];
    dystart = dy = x[k][1] - x[l][1];
    dzstart = dz = x[k][2] - x[l][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dxstart || dy != dystart || dz != dzstart) flag = 1;
  }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) 
    error->all(FLERR,"Dihedral/improper extent > half of periodic box length");
}

/* ---------------------------------------------------------------------- */

void Neighbor::improper_all()
{
  int i,m,atom1,atom2,atom3,atom4;

  int nlocal = atom->nlocal;
  int *num_improper = atom->num_improper;
  tagint **improper_atom1 = atom->improper_atom1;
  tagint **improper_atom2 = atom->improper_atom2;
  tagint **improper_atom3 = atom->improper_atom3;
  tagint **improper_atom4 = atom->improper_atom4;
  int **improper_type = atom->improper_type;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nimproperlist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_improper[i]; m++) {
      atom1 = atom->map(improper_atom1[i][m]);
      atom2 = atom->map(improper_atom2[i][m]);
      atom3 = atom->map(improper_atom3[i][m]);
      atom4 = atom->map(improper_atom4[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
        nmissing++;
        if (lostbond == ERROR) {
          char str[128];
          sprintf(str,"Improper atoms "
                  TAGINT_FORMAT " " TAGINT_FORMAT " " 
                  TAGINT_FORMAT " " TAGINT_FORMAT 
                  " missing on proc %d at step " BIGINT_FORMAT,
                  improper_atom1[i][m],improper_atom2[i][m],
                  improper_atom3[i][m],improper_atom4[i][m],
                  me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain-> closest_image(i,atom4);
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

  if (cluster_check) dihedral_check(nimproperlist,improperlist);
  if (lostbond == IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Improper atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::improper_template()
{
  int i,m,atom1,atom2,atom3,atom4;
  int imol,iatom;
  tagint tagprev;
  int *num_improper;
  int **improper_atom1,**improper_atom2,**improper_atom3,**improper_atom4;
  int **improper_type;

  Molecule **onemols = atom->avec->onemols;

  tagint *tag = atom->tag;
  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nimproperlist = 0;

  for (i = 0; i < nlocal; i++) {
    if (molindex[i] < 0) continue;
    imol = molindex[i];
    iatom = molatom[i];
    tagprev = tag[i] - iatom - 1;
    num_improper = onemols[imol]->num_improper;
    improper_atom1 = onemols[imol]->improper_atom1;
    improper_atom2 = onemols[imol]->improper_atom2;
    improper_atom3 = onemols[imol]->improper_atom3;
    improper_atom4 = onemols[imol]->improper_atom4;
    improper_type = onemols[imol]->improper_type;

    for (m = 0; m < num_improper[iatom]; m++) {
      atom1 = atom->map(improper_atom1[iatom][m]+tagprev);
      atom2 = atom->map(improper_atom2[iatom][m]+tagprev);
      atom3 = atom->map(improper_atom3[iatom][m]+tagprev);
      atom4 = atom->map(improper_atom4[iatom][m]+tagprev);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
        nmissing++;
        if (lostbond == ERROR) {
          char str[128];
          sprintf(str,"Improper atoms "
                  TAGINT_FORMAT " " TAGINT_FORMAT " " 
                  TAGINT_FORMAT " " TAGINT_FORMAT 
                  " missing on proc %d at step " BIGINT_FORMAT,
                  improper_atom1[iatom][m]+tagprev,
                  improper_atom2[iatom][m]+tagprev,
                  improper_atom3[iatom][m]+tagprev,
                  improper_atom4[iatom][m]+tagprev,
                  me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain-> closest_image(i,atom4);
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
        improperlist[nimproperlist][4] = improper_type[iatom][m];
        nimproperlist++;
      }
    }
  }

  if (cluster_check) dihedral_check(nimproperlist,improperlist);
  if (lostbond == IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Improper atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void Neighbor::improper_partial()
{
  int i,m,atom1,atom2,atom3,atom4;

  int nlocal = atom->nlocal;
  int *num_improper = atom->num_improper;
  tagint **improper_atom1 = atom->improper_atom1;
  tagint **improper_atom2 = atom->improper_atom2;
  tagint **improper_atom3 = atom->improper_atom3;
  tagint **improper_atom4 = atom->improper_atom4;
  int **improper_type = atom->improper_type;
  int newton_bond = force->newton_bond;

  int lostbond = output->thermo->lostbond;
  int nmissing = 0;
  nimproperlist = 0;

  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_improper[i]; m++) {
      if (improper_type[i][m] <= 0) continue;
      atom1 = atom->map(improper_atom1[i][m]);
      atom2 = atom->map(improper_atom2[i][m]);
      atom3 = atom->map(improper_atom3[i][m]);
      atom4 = atom->map(improper_atom4[i][m]);
      if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1) {
        nmissing++;
        if (lostbond == ERROR) {
          char str[128];
          sprintf(str,"Improper atoms "
                  TAGINT_FORMAT " " TAGINT_FORMAT " " 
                  TAGINT_FORMAT " " TAGINT_FORMAT 
                  " missing on proc %d at step " BIGINT_FORMAT,
                  improper_atom1[i][m],improper_atom2[i][m],
                  improper_atom3[i][m],improper_atom4[i][m],
                  me,update->ntimestep);
          error->one(FLERR,str);
        }
        continue;
      }
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain->closest_image(i,atom4);
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

  if (cluster_check) dihedral_check(nimproperlist,improperlist);
  if (lostbond == IGNORE) return;

  int all;
  MPI_Allreduce(&nmissing,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,
            "Improper atoms missing at step " BIGINT_FORMAT,update->ntimestep);
    if (me == 0) error->warning(FLERR,str);
  }
}
