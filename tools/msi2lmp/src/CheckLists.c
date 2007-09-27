#include "Msi2LMP2.h"

void CheckLists() {
  int i,j,k,ic_type;

  for (i=0; i < total_no_bonds; i++) {
    if ((atoms[bonds[i].members[0]].type != bondtypes[bonds[i].type].types[0])
	|| (atoms[bonds[i].members[1]].type != bondtypes[bonds[i].type].types[1])) {
      fprintf(stderr,"Warning atom types in bond %d are inconsistent with bond type %d\n",i,bonds[i].type);
    }
  }

  for (i=0; i < total_no_angles;i++) {
    if (   (atoms[angles[i].members[0]].type != angletypes[angles[i].type].types[0])
	|| (atoms[angles[i].members[1]].type != angletypes[angles[i].type].types[1])
	|| (atoms[angles[i].members[2]].type != angletypes[angles[i].type].types[2]))
      {
	fprintf(stderr,"Warning atom types in angle %d are inconsistent with angle type %d\n", i,angles[i].type);
      }
  }

  for (i=0; i < total_no_dihedrals; i++) {
    if (   (atoms[dihedrals[i].members[0]].type != dihedraltypes[dihedrals[i].type].types[0])
	|| (atoms[dihedrals[i].members[1]].type != dihedraltypes[dihedrals[i].type].types[1])
	|| (atoms[dihedrals[i].members[2]].type != dihedraltypes[dihedrals[i].type].types[2])
	|| (atoms[dihedrals[i].members[3]].type != dihedraltypes[dihedrals[i].type].types[3])) {
      fprintf(stderr,"Warning atom types in dihedral %d are inconsistent with dihedral type %d\n",i,dihedrals[i].type);
    }
  }

  for (i=0; i < total_no_oops; i++) {

    if (   (atoms[oops[i].members[0]].type != ooptypes[oops[i].type].types[0])
	|| (atoms[oops[i].members[1]].type != ooptypes[oops[i].type].types[1])
	|| (atoms[oops[i].members[2]].type != ooptypes[oops[i].type].types[2])
	|| (atoms[oops[i].members[3]].type != ooptypes[oops[i].type].types[3]))
      {
	fprintf(stderr,"Warning atom types in oop %d are inconsistent with oop type %d\n",i,oops[i].type);
      }
  }
}
