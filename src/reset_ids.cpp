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

#include "reset_ids.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "comm.h"
#include "special.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ResetIDs::ResetIDs(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void ResetIDs::command(int narg, char **/*arg*/)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Reset_ids command before simulation box is defined");
  if (narg != 0) error->all(FLERR,"Illegal reset_ids command");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use reset_ids unless atoms have IDs");

  // NOTE: check if any fixes exist which store atom IDs?
  // if so, this operation will mess up the fix

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Resetting atom IDs ...\n");
    if (logfile) fprintf(logfile,"Resetting atom IDs ...\n");
  }

  // create an atom map if one doesn't exist already

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // initialize system since comm->borders() will be invoked

  lmp->init();

  // setup domain, communication
  // acquire ghosts - is that really necessary?
  // exchange will clear map, borders will reset
  // this is the map needed to lookup current global IDs for bond topology

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  // oldIDs = copy of current owned IDs

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  tagint *oldIDs;
  memory->create(oldIDs,nlocal,"reset_ids:oldIDs");

  for (int i = 0; i < nlocal; i++) {
    oldIDs[i] = tag[i];
    tag[i] = 0;
  }

  // assign new contigous IDs to owned atoms via tag_extend()

  atom->tag_extend();

  // newIDs = copy of new IDs
  // restore old IDs, consistent with existing atom map
  // forward_comm_array acquires new IDs for ghost atoms

  double **newIDs;
  memory->create(newIDs,nall,1,"reset_ids:newIDs");

  for (int i = 0; i < nlocal; i++) {
    newIDs[i][0] = tag[i];
    tag[i] = oldIDs[i];
  }

  comm->forward_comm_array(1,newIDs);

  // loop over bonds, angles, etc and reset IDs in stored topology arrays
  // only necessary for molecular = 1, not molecular = 2
  // badcount = atom IDs that could not be found

  int badcount = 0;

  if (atom->molecular == 1) {
    int j,m;
    tagint oldID;

    if (atom->avec->bonds_allow) {
      int *num_bond = atom->num_bond;
      tagint **bond_atom = atom->bond_atom;
      for (int i = 0; i < nlocal; i++) {
        for (j = 0; j < num_bond[i]; j++) {
          oldID = bond_atom[i][j];
          m = atom->map(oldID);
          if (m >= 0) bond_atom[i][j] = static_cast<tagint> (newIDs[m][0]);
          else badcount++;
        }
      }
    }

    if (atom->avec->angles_allow) {
      int *num_angle = atom->num_angle;
      tagint **angle_atom1 = atom->angle_atom1;
      tagint **angle_atom2 = atom->angle_atom2;
      tagint **angle_atom3 = atom->angle_atom3;
      for (int i = 0; i < nlocal; i++) {
        for (j = 0; j < num_angle[i]; j++) {
          oldID = angle_atom1[i][j];
          m = atom->map(oldID);
          if (m >= 0) angle_atom1[i][j] = static_cast<tagint> (newIDs[m][0]);
          else badcount++;

          oldID = angle_atom2[i][j];
          m = atom->map(oldID);
          if (m >= 0) angle_atom2[i][j] = static_cast<tagint> (newIDs[m][0]);
          else badcount++;

          oldID = angle_atom3[i][j];
          m = atom->map(oldID);
          if (m >= 0) angle_atom3[i][j] = static_cast<tagint> (newIDs[m][0]);
          else badcount++;
        }
      }
    }

    if (atom->avec->dihedrals_allow) {
      int *num_dihedral = atom->num_dihedral;
      tagint **dihedral_atom1 = atom->dihedral_atom1;
      tagint **dihedral_atom2 = atom->dihedral_atom2;
      tagint **dihedral_atom3 = atom->dihedral_atom3;
      tagint **dihedral_atom4 = atom->dihedral_atom4;
      for (int i = 0; i < nlocal; i++) {
        for (j = 0; j < num_dihedral[i]; j++) {
          oldID = dihedral_atom1[i][j];
          m = atom->map(oldID);
          if (m >= 0) dihedral_atom1[i][j] = static_cast<tagint> (newIDs[m][0]);
          else badcount++;

          oldID = dihedral_atom2[i][j];
          m = atom->map(oldID);
          if (m >= 0) dihedral_atom2[i][j] = static_cast<tagint> (newIDs[m][0]);
          else badcount++;

          oldID = dihedral_atom3[i][j];
          m = atom->map(oldID);
          if (m >= 0) dihedral_atom3[i][j] = static_cast<tagint> (newIDs[m][0]);
          else badcount++;

          oldID = dihedral_atom4[i][j];
          m = atom->map(oldID);
          if (m >= 0) dihedral_atom4[i][j] = static_cast<tagint> (newIDs[m][0]);
          else badcount++;
        }
      }
    }

    if (atom->avec->impropers_allow) {
      int *num_improper = atom->num_improper;
      tagint **improper_atom1 = atom->improper_atom1;
      tagint **improper_atom2 = atom->improper_atom2;
      tagint **improper_atom3 = atom->improper_atom3;
      tagint **improper_atom4 = atom->improper_atom4;
      for (int i = 0; i < nlocal; i++) {
        for (j = 0; j < num_improper[i]; j++) {
          oldID = improper_atom1[i][j];
          m = atom->map(oldID);
          if (m >= 0) improper_atom1[i][j] = static_cast<tagint> (newIDs[m][0]);
          else badcount++;

          oldID = improper_atom2[i][j];
          m = atom->map(oldID);
          if (m >= 0) improper_atom2[i][j] = static_cast<tagint> (newIDs[m][0]);
          else badcount++;

          oldID = improper_atom3[i][j];
          m = atom->map(oldID);
          if (m >= 0) improper_atom3[i][j] = static_cast<tagint> (newIDs[m][0]);
          else badcount++;

          oldID = improper_atom4[i][j];
          m = atom->map(oldID);
          if (m >= 0) improper_atom4[i][j] = static_cast<tagint> (newIDs[m][0]);
          else badcount++;
        }
      }
    }
  }

  // error check

  int all;
  MPI_Allreduce(&badcount,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,"Reset_ids missing %d bond topology atom IDs - "
            "use comm_modify cutoff",all);
    error->all(FLERR,str);
  }

  // reset IDs and atom map for owned atoms

  atom->map_clear();
  atom->nghost = 0;
  for (int i = 0; i < nlocal; i++) tag[i] = static_cast<tagint> (newIDs[i][0]);
  atom->map_init();
  atom->map_set();

  // need to update exclusions with new atom IDs

  if (atom->molecular == 1) {
    Special special(lmp);
    special.build();
  }

  // delete temporary atom map

  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  }

  // clean up

  memory->destroy(oldIDs);
  memory->destroy(newIDs);
}
