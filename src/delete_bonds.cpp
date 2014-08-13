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
#include "stdlib.h"
#include "string.h"
#include "delete_bonds.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "neighbor.h"
#include "comm.h"
#include "force.h"
#include "group.h"
#include "special.h"
#include "error.h"

#include <stdlib.h>

using namespace LAMMPS_NS;

enum{MULTI,ATOM,BOND,ANGLE,DIHEDRAL,IMPROPER,STATS};

/* ---------------------------------------------------------------------- */

DeleteBonds::DeleteBonds(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void DeleteBonds::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Delete_bonds command before simulation box is defined");
  if (atom->natoms == 0)
    error->all(FLERR,"Delete_bonds command with no atoms existing");
  if (atom->molecular != 1)
    error->all(FLERR,"Cannot use delete_bonds with non-molecular system");

  if (narg < 2) error->all(FLERR,"Illegal delete_bonds command");

  // init entire system since comm->borders is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  if (comm->me == 0 && screen)
    fprintf(screen,"System init for delete_bonds ...\n");
  lmp->init();

  if (comm->me == 0 && screen) fprintf(screen,"Deleting bonds ...\n");

  // identify group

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Cannot find delete_bonds group ID");
  int groupbit = group->bitmask[igroup];

  // set style and which = type value

  int style = -1;
  if (strcmp(arg[1],"multi") == 0) style = MULTI;
  else if (strcmp(arg[1],"atom") == 0) style = ATOM;
  else if (strcmp(arg[1],"bond") == 0) style = BOND;
  else if (strcmp(arg[1],"angle") == 0) style = ANGLE;
  else if (strcmp(arg[1],"dihedral") == 0) style = DIHEDRAL;
  else if (strcmp(arg[1],"improper") == 0) style = IMPROPER;
  else if (strcmp(arg[1],"stats") == 0) style = STATS;
  else error->all(FLERR,"Illegal delete_bonds command");

  // setup list of types (atom,bond,etc) to consider
  // use force->bounds() to allow setting of range of types
  // range can be 0 to ntypes inclusive

  int *tlist = NULL;

  int iarg = 2;
  if (style != MULTI && style != STATS) {
    if (narg < 3) error->all(FLERR,"Illegal delete_bonds command");

    int n = -1;
    if (style == ATOM) n = atom->ntypes;
    if (style == BOND) n = atom->nbondtypes;
    if (style == ANGLE) n = atom->nangletypes;
    if (style == DIHEDRAL) n = atom->ndihedraltypes;
    if (style == IMPROPER) n = atom->nimpropertypes;

    tlist = new int[n+1];
    for (int i = 0; i <= n; i++) tlist[i] = 0;
    int nlo,nhi;
    force->bounds(arg[2],n,nlo,nhi,0);
    for (int i = nlo; i <= nhi; i++) tlist[i] = 1;

    iarg++;
  }

  // grab optional keywords

  int any_flag = 0;
  int undo_flag = 0;
  int remove_flag = 0;
  int special_flag = 0;
  int induce_flag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"any") == 0) any_flag = 1;
    else if (strcmp(arg[iarg],"undo") == 0) undo_flag = 1;
    else if (strcmp(arg[iarg],"remove") == 0) remove_flag = 1;
    else if (strcmp(arg[iarg],"special") == 0) special_flag = 1;
    else if (strcmp(arg[iarg],"induce") == 0) induce_flag = 1;
    else error->all(FLERR,"Illegal delete_bonds command");
    iarg++;
  }

  // border swap to insure type and mask is current for off-proc atoms
  // enforce PBC before in case atoms are outside box

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  // set topology interactions either off or on
  // criteria for an interaction to potentially be changed (set flag = 1)
  //   all atoms or any atom in interaction must be in group, based on any_flag
  //   for style = MULTI, all bond/angle/dihedral/improper, no other criteria
  //   for style = ATOM, same as MULTI, plus at least one atom is specified type
  //   for style = BOND/ANGLE/DIHEDRAL/IMPROPER, interaction is specified type
  //   for style = STATS only compute stats, flag is always 0
  // if flag = 1
  //   set interaction type negative if undo_flag = 0
  //   set interaction type positive if undo_flag = 1

  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  int i,m,n,consider,flag,itype;
  int atom1,atom2,atom3,atom4;

  if (atom->avec->bonds_allow && 
      (style == BOND || style == MULTI || style == ATOM)) {
    int *num_bond = atom->num_bond;
    int **bond_type = atom->bond_type;

    for (i = 0; i < nlocal; i++) {
      for (m = 0; m < num_bond[i]; m++) {
        atom1 = atom->map(atom->bond_atom[i][m]);
        if (atom1 == -1) error->one(FLERR,"Bond atom missing in delete_bonds");
        consider = 0;
        if (!any_flag && mask[i] & groupbit && mask[atom1] & groupbit)
          consider = 1;
        if (any_flag && (mask[i] & groupbit || mask[atom1] & groupbit))
          consider = 1;
        if (consider) {
          flag = 0;
          if (style == MULTI) flag = 1;
          else if (style == ATOM) {
            if (tlist[type[i]] || tlist[type[atom1]]) flag = 1;
          } else if (style == BOND) {
            itype = abs(bond_type[i][m]);
            if (tlist[itype]) flag = 1;
          }
          if (flag) {
            if (undo_flag == 0 && bond_type[i][m] > 0)
              bond_type[i][m] = -bond_type[i][m];
            if (undo_flag == 1 && bond_type[i][m] < 0)
              bond_type[i][m] = -bond_type[i][m];
          }
        }
      }
    }
  }

  if (atom->avec->angles_allow &&
      (style == ANGLE || style == MULTI || style == ATOM)) {
    int *num_angle = atom->num_angle;
    int **angle_type = atom->angle_type;

    for (i = 0; i < nlocal; i++) {
      for (m = 0; m < num_angle[i]; m++) {
        atom1 = atom->map(atom->angle_atom1[i][m]);
        atom2 = atom->map(atom->angle_atom2[i][m]);
        atom3 = atom->map(atom->angle_atom3[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1)
          error->one(FLERR,"Angle atom missing in delete_bonds");
        consider = 0;
        if (!any_flag && mask[atom1] & groupbit && mask[atom2] & groupbit &&
            mask[atom3] & groupbit) consider = 1;
        if (any_flag && (mask[atom1] & groupbit || mask[atom2] & groupbit ||
                          mask[atom3] & groupbit)) consider = 1;
        if (consider) {
          flag = 0;
          if (style == MULTI) flag = 1;
          else if (style == ATOM) {
            if (tlist[type[atom1]] || tlist[type[atom2]] ||
                tlist[type[atom3]]) flag = 1;
          } else if (style == ANGLE) {
            itype = abs(angle_type[i][m]);
            if (tlist[itype]) flag = 1;
          }
          if (flag) {
            if (undo_flag == 0 && angle_type[i][m] > 0)
              angle_type[i][m] = -angle_type[i][m];
            if (undo_flag == 1 && angle_type[i][m] < 0)
              angle_type[i][m] = -angle_type[i][m];
          }
        }
      }
    }
  }

  if (atom->avec->dihedrals_allow &&
      (style == DIHEDRAL || style == MULTI || style == ATOM)) {
    int *num_dihedral = atom->num_dihedral;
    int **dihedral_type = atom->dihedral_type;

    for (i = 0; i < nlocal; i++) {
      for (m = 0; m < num_dihedral[i]; m++) {
        atom1 = atom->map(atom->dihedral_atom1[i][m]);
        atom2 = atom->map(atom->dihedral_atom2[i][m]);
        atom3 = atom->map(atom->dihedral_atom3[i][m]);
        atom4 = atom->map(atom->dihedral_atom4[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
          error->one(FLERR,"Dihedral atom missing in delete_bonds");
        consider = 0;
        if (!any_flag && mask[atom1] & groupbit && mask[atom2] & groupbit &&
            mask[atom3] & groupbit && mask[atom4] & groupbit) consider = 1;
        if (any_flag && (mask[atom1] & groupbit || mask[atom2] & groupbit ||
                         mask[atom3] & groupbit || mask[atom4] & groupbit))
          consider = 1;
        if (consider) {
          flag = 0;
          if (style == MULTI) flag = 1;
          else if (style == ATOM) {
              if (tlist[type[atom1]] || tlist[type[atom2]] ||
                  tlist[type[atom3]] || tlist[type[atom4]]) flag = 1;
          } else if (style == DIHEDRAL) {
            itype = abs(dihedral_type[i][m]);
            if (tlist[itype]) flag = 1;
          }
          if (flag) {
            if (undo_flag == 0 && dihedral_type[i][m] > 0)
              dihedral_type[i][m] = -dihedral_type[i][m];
            if (undo_flag == 1 && dihedral_type[i][m] < 0)
              dihedral_type[i][m] = -dihedral_type[i][m];
          }
        }
      }
    }
  }

  if (atom->avec->impropers_allow &&
      (style == IMPROPER || style == MULTI || style == ATOM)) {
    int *num_improper = atom->num_improper;
    int **improper_type = atom->improper_type;

    for (i = 0; i < nlocal; i++) {
      for (m = 0; m < num_improper[i]; m++) {
        atom1 = atom->map(atom->improper_atom1[i][m]);
        atom2 = atom->map(atom->improper_atom2[i][m]);
        atom3 = atom->map(atom->improper_atom3[i][m]);
        atom4 = atom->map(atom->improper_atom4[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
          error->one(FLERR,"Improper atom missing in delete_bonds");
        consider = 0;
        if (!any_flag && mask[atom1] & groupbit && mask[atom2] & groupbit &&
            mask[atom3] & groupbit && mask[atom4] & groupbit) consider = 1;
        if (any_flag && (mask[atom1] & groupbit || mask[atom2] & groupbit ||
                         mask[atom3] & groupbit || mask[atom4] & groupbit))
          consider = 1;
        if (consider) {
          flag = 0;
          if (style == MULTI) flag = 1;
          else if (style == ATOM) {
              if (tlist[type[atom1]] || tlist[type[atom2]] ||
                  tlist[type[atom3]] || tlist[type[atom4]]) flag = 1;
          } else if (style == IMPROPER) {
            itype = abs(improper_type[i][m]);
            if (tlist[itype]) flag = 1;
          }
          if (flag) {
            if (undo_flag == 0 && improper_type[i][m] > 0)
              improper_type[i][m] = -improper_type[i][m];
            if (undo_flag == 1 && improper_type[i][m] < 0)
              improper_type[i][m] = -improper_type[i][m];
          }
        }
      }
    }
  }

  delete [] tlist;
    
  // induce turn off of angles, dihedral, impropers due to turned off bonds
  // induce turn off of dihedrals due to turned off angles
  // all atoms or any atom in interaction must be in group, based on any_flag

  if (induce_flag) {

    // circulate list of turned off bonds around ring of procs

    // circulate list of turned off angles around ring of procs

  }

  // remove interactions if requested
  // all atoms or any atom in interaction must be in group, based on any_flag

  if (remove_flag) {

    if (atom->avec->bonds_allow) {
      for (i = 0; i < nlocal; i++) {
        m = 0;
        while (m < atom->num_bond[i]) {
          if (atom->bond_type[i][m] <= 0) {
            atom1 = atom->map(atom->bond_atom[i][m]);
            flag = 0;
            if (!any_flag && mask[i] & groupbit && mask[atom1] & groupbit)
              flag = 1;
            if (any_flag && (mask[i] & groupbit || mask[atom1] & groupbit))
              flag = 1;
            if (flag) {
              n = atom->num_bond[i];
              atom->bond_type[i][m] = atom->bond_type[i][n-1];
              atom->bond_atom[i][m] = atom->bond_atom[i][n-1];
              atom->num_bond[i]--;
            } else m++;
          } else m++;
        }
      }
    }

    if (atom->avec->angles_allow) {
      for (i = 0; i < nlocal; i++) {
        m = 0;
        while (m < atom->num_angle[i]) {
          if (atom->angle_type[i][m] <= 0) {
            atom1 = atom->map(atom->angle_atom1[i][m]);
            atom2 = atom->map(atom->angle_atom2[i][m]);
            atom3 = atom->map(atom->angle_atom3[i][m]);
            flag = 0;
            if (!any_flag && mask[atom1] & groupbit && mask[atom2] & groupbit &&
                mask[atom3] & groupbit) flag = 1;
            if (any_flag && (mask[atom1] & groupbit || mask[atom2] & groupbit ||
                             mask[atom3] & groupbit)) flag = 1;
            if (flag) {
              n = atom->num_angle[i];
              atom->angle_type[i][m] = atom->angle_type[i][n-1];
              atom->angle_atom1[i][m] = atom->angle_atom1[i][n-1];
              atom->angle_atom2[i][m] = atom->angle_atom2[i][n-1];
              atom->angle_atom3[i][m] = atom->angle_atom3[i][n-1];
              atom->num_angle[i]--;
            } else m++;
          } else m++;
        }
      }
    }

    if (atom->avec->dihedrals_allow) {
      for (i = 0; i < nlocal; i++) {
        m = 0;
        while (m < atom->num_dihedral[i]) {
          if (atom->dihedral_type[i][m] <= 0) {
            atom1 = atom->map(atom->dihedral_atom1[i][m]);
            atom2 = atom->map(atom->dihedral_atom2[i][m]);
            atom3 = atom->map(atom->dihedral_atom3[i][m]);
            atom4 = atom->map(atom->dihedral_atom4[i][m]);
            flag = 0;
            if (!any_flag && mask[atom1] & groupbit && mask[atom2] & groupbit &&
                mask[atom3] & groupbit && mask[atom4] & groupbit) flag = 1;
            if (any_flag && (mask[atom1] & groupbit || mask[atom2] & groupbit ||
                             mask[atom3] & groupbit || mask[atom4] & groupbit))
              flag = 1;
            if (flag) {
              n = atom->num_dihedral[i];
              atom->dihedral_type[i][m] = atom->dihedral_type[i][n-1];
              atom->dihedral_atom1[i][m] = atom->dihedral_atom1[i][n-1];
              atom->dihedral_atom2[i][m] = atom->dihedral_atom2[i][n-1];
              atom->dihedral_atom3[i][m] = atom->dihedral_atom3[i][n-1];
              atom->dihedral_atom4[i][m] = atom->dihedral_atom4[i][n-1];
              atom->num_dihedral[i]--;
            } else m++;
          } else m++;
        }
      }
    }

    if (atom->avec->impropers_allow) {
      for (i = 0; i < nlocal; i++) {
        m = 0;
        while (m < atom->num_improper[i]) {
          if (atom->improper_type[i][m] <= 0) {
            atom1 = atom->map(atom->improper_atom1[i][m]);
            atom2 = atom->map(atom->improper_atom2[i][m]);
            atom3 = atom->map(atom->improper_atom3[i][m]);
            atom4 = atom->map(atom->improper_atom4[i][m]);
            flag = 0;
            if (!any_flag && mask[atom1] & groupbit && mask[atom2] & groupbit &&
                mask[atom3] & groupbit && mask[atom4] & groupbit) flag = 1;
            if (any_flag && (mask[atom1] & groupbit || mask[atom2] & groupbit ||
                             mask[atom3] & groupbit || mask[atom4] & groupbit))
              flag = 1;
            if (flag) {
              n = atom->num_improper[i];
              atom->improper_type[i][m] = atom->improper_type[i][n-1];
              atom->improper_atom1[i][m] = atom->improper_atom1[i][n-1];
              atom->improper_atom2[i][m] = atom->improper_atom2[i][n-1];
              atom->improper_atom3[i][m] = atom->improper_atom3[i][n-1];
              atom->improper_atom4[i][m] = atom->improper_atom4[i][n-1];
              atom->num_improper[i]--;
            } else m++;
          } else m++;
        }
      }
    }

  }

  // if interactions were removed, recompute global counts

  if (remove_flag) {

    if (atom->avec->bonds_allow) {
      bigint nbonds = 0;
      for (i = 0; i < nlocal; i++) nbonds += atom->num_bond[i];
      MPI_Allreduce(&nbonds,&atom->nbonds,1,MPI_LMP_BIGINT,
                    MPI_SUM,world);
      if (force->newton_bond == 0) atom->nbonds /= 2;
    }

    if (atom->avec->angles_allow) {
      bigint nangles = 0;
      for (i = 0; i < nlocal; i++) nangles += atom->num_angle[i];
      MPI_Allreduce(&nangles,&atom->nangles,1,MPI_LMP_BIGINT,
                    MPI_SUM,world);
      if (force->newton_bond == 0) atom->nangles /= 3;
    }

    if (atom->avec->dihedrals_allow) {
      bigint ndihedrals = 0;
      for (i = 0; i < nlocal; i++) ndihedrals += atom->num_dihedral[i];
      MPI_Allreduce(&ndihedrals,&atom->ndihedrals,
                    1,MPI_LMP_BIGINT,MPI_SUM,world);
      if (force->newton_bond == 0) atom->ndihedrals /= 4;
    }

    if (atom->avec->impropers_allow) {
      bigint nimpropers = 0;
      for (i = 0; i < nlocal; i++) nimpropers += atom->num_improper[i];
      MPI_Allreduce(&nimpropers,&atom->nimpropers,
                    1,MPI_LMP_BIGINT,MPI_SUM,world);
      if (force->newton_bond == 0) atom->nimpropers /= 4;
    }

  }

  // compute and print stats

  bigint tmp;
  bigint bond_on,bond_off;
  bigint angle_on,angle_off;
  bigint dihedral_on,dihedral_off;
  bigint improper_on,improper_off;

  if (atom->avec->bonds_allow) {
    bond_on = bond_off = 0;
    for (i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_bond[i]; m++)
        if (atom->bond_type[i][m] > 0) bond_on++;
        else bond_off++;
    MPI_Allreduce(&bond_on,&tmp,1,MPI_LMP_BIGINT,MPI_SUM,world);
    bond_on = tmp;
    MPI_Allreduce(&bond_off,&tmp,1,MPI_LMP_BIGINT,MPI_SUM,world);
    bond_off = tmp;
    if (force->newton_bond == 0) {
      bond_on /= 2;
      bond_off /= 2;
    }
  }

  if (atom->avec->angles_allow) {
    angle_on = angle_off = 0;
    for (i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_angle[i]; m++)
        if (atom->angle_type[i][m] > 0) angle_on++;
        else angle_off++;
    MPI_Allreduce(&angle_on,&tmp,1,MPI_LMP_BIGINT,MPI_SUM,world);
    angle_on = tmp;
    MPI_Allreduce(&angle_off,&tmp,1,MPI_LMP_BIGINT,MPI_SUM,world);
    angle_off = tmp;
    if (force->newton_bond == 0) {
      angle_on /= 3;
      angle_off /= 3;
    }
  }

  if (atom->avec->dihedrals_allow) {
    dihedral_on = dihedral_off = 0;
    for (i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_dihedral[i]; m++)
        if (atom->dihedral_type[i][m] > 0) dihedral_on++;
        else dihedral_off++;
    MPI_Allreduce(&dihedral_on,&tmp,1,MPI_LMP_BIGINT,MPI_SUM,world);
    dihedral_on = tmp;
    MPI_Allreduce(&dihedral_off,&tmp,1,MPI_LMP_BIGINT,MPI_SUM,world);
    dihedral_off = tmp;
    if (force->newton_bond == 0) {
      dihedral_on /= 4;
      dihedral_off /= 4;
    }
  }

  if (atom->avec->impropers_allow) {
    improper_on = improper_off = 0;
    for (i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_improper[i]; m++)
        if (atom->improper_type[i][m] > 0) improper_on++;
        else improper_off++;
    MPI_Allreduce(&improper_on,&tmp,1,MPI_LMP_BIGINT,MPI_SUM,world);
    improper_on = tmp;
    MPI_Allreduce(&improper_off,&tmp,1,MPI_LMP_BIGINT,MPI_SUM,world);
    improper_off = tmp;
    if (force->newton_bond == 0) {
      improper_on /= 4;
      improper_off /= 4;
    }
  }

  if (comm->me == 0) {
    if (atom->avec->bonds_allow) {
      if (screen) fprintf(screen,
                          "  " BIGINT_FORMAT " total bonds, " BIGINT_FORMAT
                          " turned on, " BIGINT_FORMAT " turned off\n",
                          atom->nbonds,bond_on,bond_off);
      if (logfile) fprintf(logfile,
                           "  " BIGINT_FORMAT " total bonds, " BIGINT_FORMAT
                           " turned on, " BIGINT_FORMAT " turned off\n",
                           atom->nbonds,bond_on,bond_off);
    }
    if (atom->avec->angles_allow) {
      if (screen) fprintf(screen,
                          "  " BIGINT_FORMAT " total angles, " BIGINT_FORMAT
                          " turned on, " BIGINT_FORMAT " turned off\n",
                          atom->nangles,angle_on,angle_off);
      if (logfile) fprintf(logfile,
                          "  " BIGINT_FORMAT " total angles, " BIGINT_FORMAT
                           " turned on, " BIGINT_FORMAT " turned off\n",
                           atom->nangles,angle_on,angle_off);
    }
    if (atom->avec->dihedrals_allow) {
      if (screen) fprintf(screen,
                          "  " BIGINT_FORMAT " total dihedrals, "
                          BIGINT_FORMAT " turned on, " BIGINT_FORMAT
                          " turned off\n",
                          atom->ndihedrals,dihedral_on,dihedral_off);
      if (logfile) fprintf(logfile,
                          "  " BIGINT_FORMAT " total dihedrals, "
                          BIGINT_FORMAT " turned on, " BIGINT_FORMAT
                          " turned off\n",
                          atom->ndihedrals,dihedral_on,dihedral_off);
    }
    if (atom->avec->impropers_allow) {
      if (screen) fprintf(screen,
                          "  " BIGINT_FORMAT " total impropers, "
                          BIGINT_FORMAT " turned on, " BIGINT_FORMAT
                          " turned off\n",
                          atom->nimpropers,improper_on,improper_off);
      if (logfile) fprintf(logfile,
                          "  " BIGINT_FORMAT " total impropers, "
                          BIGINT_FORMAT " turned on, " BIGINT_FORMAT
                          " turned off\n",
                          atom->nimpropers,improper_on,improper_off);
    }
  }

  // re-compute special list if requested

  if (special_flag) {
    Special special(lmp);
    special.build();
  }
}
