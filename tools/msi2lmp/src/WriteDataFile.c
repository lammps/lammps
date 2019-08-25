/*
*  This function creates and writes the data file to be used with LAMMPS
*/

#include "msi2lmp.h"
#include "Forcefield.h"

#include <stdlib.h>

void WriteDataFile(char *nameroot)
{
  int i,j,k,m;
  char line[MAX_LINE_LENGTH];
  FILE *DatF;

  /* Open data file */

  sprintf(line,"%s.data",rootname);
  if (pflag > 0) {
    printf(" Writing LAMMPS data file %s.data",rootname);
    if (forcefield & FF_TYPE_CLASS1) puts(" for Class I force field");
    if (forcefield & FF_TYPE_CLASS2) puts(" for Class II force field");
    if (forcefield & FF_TYPE_OPLSAA) puts(" for OPLS-AA force field");
  }

  if ((DatF = fopen(line,"w")) == NULL ) {
    printf("Cannot open %s\n",line);
    exit(62);
  }

  if (forcefield & (FF_TYPE_CLASS1|FF_TYPE_OPLSAA)) total_no_angle_angles = 0;

  if (hintflag) fprintf(DatF, "LAMMPS data file. msi2lmp " MSI2LMP_VERSION
                        " / CGCMM for %s\n\n", nameroot);
  else fprintf(DatF, "LAMMPS data file. msi2lmp " MSI2LMP_VERSION
               " for %s\n\n", nameroot);
  fprintf(DatF, " %6d atoms\n", total_no_atoms);
  fprintf(DatF, " %6d bonds\n", total_no_bonds);
  fprintf(DatF, " %6d angles\n",total_no_angles);
  fprintf(DatF, " %6d dihedrals\n", total_no_dihedrals);
  fprintf(DatF, " %6d impropers\n", total_no_oops+total_no_angle_angles);
  fputs("\n",DatF);


  fprintf(DatF, " %3d atom types\n", no_atom_types);
  if (no_bond_types > 0)
    fprintf(DatF, " %3d bond types\n", no_bond_types);
  if (no_angle_types> 0)
    fprintf(DatF, " %3d angle types\n", no_angle_types);
  if (no_dihedral_types > 0) fprintf  (DatF," %3d dihedral types\n",
                                       no_dihedral_types);
  if (forcefield & FF_TYPE_CLASS1) {
    if (no_oop_types > 0)
      fprintf  (DatF, " %3d improper types\n", no_oop_types);
  }

  if (forcefield & FF_TYPE_CLASS2) {
    if ((no_oop_types + no_angleangle_types) > 0)
      fprintf  (DatF, " %3d improper types\n",
                no_oop_types + no_angleangle_types);
  }

  /* Modified by SLTM to print out triclinic box types 10/05/10 - lines 56-68 */

  if (TriclinicFlag == 0) {
    fputs("\n",DatF);
    fprintf(DatF, " %15.9f %15.9f xlo xhi\n", box[0][0], box[1][0]);
    fprintf(DatF, " %15.9f %15.9f ylo yhi\n", box[0][1], box[1][1]);
    fprintf(DatF, " %15.9f %15.9f zlo zhi\n", box[0][2], box[1][2]);
  } else {
    fputs("\n",DatF);
    fprintf(DatF, " %15.9f %15.9f xlo xhi\n", box[0][0], box[1][0]);
    fprintf(DatF, " %15.9f %15.9f ylo yhi\n", box[0][1], box[1][1]);
    fprintf(DatF, " %15.9f %15.9f zlo zhi\n", box[0][2], box[1][2]);
    fprintf(DatF, " %15.9f %15.9f %15.9f xy xz yz\n",box[2][0], box[2][1], box[2][2]);
  }

  /* MASSES */

  fprintf(DatF, "\nMasses\n\n");
  for(k=0; k < no_atom_types; k++) {
    if (hintflag) fprintf(DatF, " %3d %10.6f # %s\n",k+1,atomtypes[k].mass,atomtypes[k].potential);
    else fprintf(DatF, " %3d %10.6f\n",k+1,atomtypes[k].mass);
  }
  fputs("\n",DatF);


  /* COEFFICIENTS */

  fputs("Pair Coeffs",DatF);
  if (hintflag) {
    if (forcefield & (FF_TYPE_CLASS1|FF_TYPE_OPLSAA))
      fputs(" # lj/cut/coul/long\n\n",DatF);
    else if (forcefield & FF_TYPE_CLASS2)
      fputs(" # lj/class2/coul/long\n\n",DatF);
  } else fputs("\n\n",DatF);

  for (i=0; i < no_atom_types; i++) {
    fprintf(DatF, " %3i ", i+1);
    for ( j = 0; j < 2; j++)
      fprintf(DatF, "%14.10f ",atomtypes[i].params[j]);

    if (hintflag) fprintf(DatF, "# %s\n",atomtypes[i].potential);
    else fputs("\n",DatF);
  }
  fputs("\n",DatF);

  if (no_bond_types > 0) {
    m = 0;
    if (forcefield & FF_TYPE_CLASS1) m = 2;
    if (forcefield & FF_TYPE_OPLSAA) m = 2;
    if (forcefield & FF_TYPE_CLASS2) m = 4;

    fputs("Bond Coeffs",DatF);
    if (hintflag) {
      if (forcefield & (FF_TYPE_CLASS1|FF_TYPE_OPLSAA))
        fputs(" # harmonic\n\n",DatF);
      else if (forcefield & FF_TYPE_CLASS2)
        fputs(" # class2\n\n",DatF);
    } else fputs("\n\n",DatF);

    for (i=0; i < no_bond_types; i++) {
      fprintf(DatF, " %3i", i+1);
      for ( j = 0; j < m; j++)
        fprintf(DatF, " %10.4f", bondtypes[i].params[j]);

      if (hintflag) fprintf(DatF," # %s-%s\n",atomtypes[bondtypes[i].types[0]].potential,
                            atomtypes[bondtypes[i].types[1]].potential);
      else fputs("\n",DatF);
    }
    fputs("\n",DatF);
  }

  if (no_angle_types > 0) {
    m = 0;
    if (forcefield & FF_TYPE_CLASS1) m = 2;
    if (forcefield & FF_TYPE_OPLSAA) m = 2;
    if (forcefield & FF_TYPE_CLASS2) m = 4;

    fputs("Angle Coeffs",DatF);
    if (hintflag) {
      if (forcefield & (FF_TYPE_CLASS1|FF_TYPE_OPLSAA))
        fputs(" # harmonic\n\n",DatF);
      else if (forcefield & FF_TYPE_CLASS2)
        fputs(" # class2\n\n",DatF);
    } else fputs("\n\n",DatF);

    for (i=0; i < no_angle_types; i++) {
      fprintf(DatF, " %3i", i+1);
      for ( j = 0; j < m; j++)
        fprintf(DatF, " %10.4f", angletypes[i].params[j]);

      if (hintflag) fprintf(DatF," # %s-%s-%s\n",
                            atomtypes[angletypes[i].types[0]].potential,
                            atomtypes[angletypes[i].types[1]].potential,
                            atomtypes[angletypes[i].types[2]].potential);
      else fputs("\n",DatF);
    }
    fputs("\n",DatF);
  }

  if (no_dihedral_types > 0) {

    fputs("Dihedral Coeffs",DatF);
    if (forcefield & FF_TYPE_CLASS1) {

      if (hintflag) fputs(" # harmonic\n\n",DatF);
      else fputs("\n\n",DatF);

      for (i=0; i < no_dihedral_types; i++) {
        fprintf(DatF, "%3i %10.4f %3i %3i", i+1,
                dihedraltypes[i].params[0],
                (int) dihedraltypes[i].params[1],
                (int) dihedraltypes[i].params[2]);
        if (hintflag) fprintf(DatF," # %s-%s-%s-%s\n",
                              atomtypes[dihedraltypes[i].types[0]].potential,
                              atomtypes[dihedraltypes[i].types[1]].potential,
                              atomtypes[dihedraltypes[i].types[2]].potential,
                              atomtypes[dihedraltypes[i].types[3]].potential);
        else fputs("\n",DatF);
      }
      fputs("\n",DatF);
    } else if (forcefield & FF_TYPE_OPLSAA) {

      if (hintflag) fputs(" # opls\n\n",DatF);
      else fputs("\n\n",DatF);

      for (i=0; i < no_dihedral_types; i++) {
        fprintf(DatF, " %3i",i+1);
        for ( j = 0; j < 4; j++)
          fprintf(DatF, " %10.4f",dihedraltypes[i].params[j]);

        if (hintflag) fprintf(DatF," # %s-%s-%s-%s\n",
                              atomtypes[dihedraltypes[i].types[0]].potential,
                              atomtypes[dihedraltypes[i].types[1]].potential,
                              atomtypes[dihedraltypes[i].types[2]].potential,
                              atomtypes[dihedraltypes[i].types[3]].potential);
        else fputs("\n",DatF);
      }
      fputs("\n",DatF);
    } else if (forcefield & FF_TYPE_CLASS2) {

      if (hintflag) fputs(" # class2\n\n",DatF);
      else fputs("\n\n",DatF);

      for (i=0; i < no_dihedral_types; i++) {
        fprintf(DatF, " %3i",i+1);
        for ( j = 0; j < 6; j++)
          fprintf(DatF, " %10.4f",dihedraltypes[i].params[j]);

        if (hintflag) fprintf(DatF,"# %s-%s-%s-%s\n",
                              atomtypes[dihedraltypes[i].types[0]].potential,
                              atomtypes[dihedraltypes[i].types[1]].potential,
                              atomtypes[dihedraltypes[i].types[2]].potential,
                              atomtypes[dihedraltypes[i].types[3]].potential);
        else fputs("\n",DatF);
      }
      fputs("\n",DatF);
    }
  }

  if (forcefield & FF_TYPE_CLASS1) {
    if (no_oop_types > 0) {
      /* cvff improper coeffs are: type K0 d n */
      if (hintflag) fputs("Improper Coeffs # cvff\n\n",DatF);
      else fputs("Improper Coeffs\n\n",DatF);

      for (i=0; i < no_oop_types; i++) {
        fprintf(DatF,"%5i %10.4f %3i %3i ",i+1,
                ooptypes[i].params[0], (int) ooptypes[i].params[1],
                (int) ooptypes[i].params[2]);

        if (hintflag) fprintf(DatF,"# %s-%s-%s-%s\n",
                              atomtypes[ooptypes[i].types[0]].potential,
                              atomtypes[ooptypes[i].types[1]].potential,
                              atomtypes[ooptypes[i].types[2]].potential,
                              atomtypes[ooptypes[i].types[3]].potential);
        else fputs("\n",DatF);
      }
      fputs("\n",DatF);
    }
  } else if (forcefield & FF_TYPE_OPLSAA) {
    if (no_oop_types > 0) {
      /* opls improper coeffs are like cvff: type K0 d(=-1) n(=2) */
      if (hintflag) fputs("Improper Coeffs # cvff\n\n",DatF);
      else fputs("Improper Coeffs\n\n",DatF);

      for (i=0; i < no_oop_types; i++) {
        fprintf(DatF,"%5i %10.4f %3i %3i ",i+1,
                ooptypes[i].params[0], (int) ooptypes[i].params[1],
                (int) ooptypes[i].params[2]);

        if (hintflag) fprintf(DatF,"# %s-%s-%s-%s\n",
                              atomtypes[ooptypes[i].types[0]].potential,
                              atomtypes[ooptypes[i].types[1]].potential,
                              atomtypes[ooptypes[i].types[2]].potential,
                              atomtypes[ooptypes[i].types[3]].potential);
        else fputs("\n",DatF);
      }
      fputs("\n",DatF);
    }
  } else if (forcefield & FF_TYPE_CLASS2) {
    if ((no_oop_types + no_angleangle_types) > 0) {
      if (hintflag) fputs("Improper Coeffs # class2\n\n",DatF);
      else fputs("Improper Coeffs\n\n",DatF);

      for (i=0; i < no_oop_types; i++) {
        fprintf(DatF, "%3i ", i+1);
        for ( j = 0; j < 2; j++)
          fprintf(DatF, "%10.4f ", ooptypes[i].params[j]);

        if (hintflag) fprintf(DatF,"# %s-%s-%s-%s\n",
                              atomtypes[ooptypes[i].types[0]].potential,
                              atomtypes[ooptypes[i].types[1]].potential,
                              atomtypes[ooptypes[i].types[2]].potential,
                              atomtypes[ooptypes[i].types[3]].potential);
        else fputs("\n",DatF);
      }
      for (i=0; i < no_angleangle_types; i++) {
        fprintf(DatF, "%3i ", i+no_oop_types+1);
        for ( j = 0; j < 2; j++)
          fprintf(DatF, "%10.4f ", 0.0);
        fputs("\n",DatF);
      }
      fputs("\n",DatF);
    }
  }

  if (forcefield & FF_TYPE_CLASS2) {

    if (no_angle_types > 0) {
      fprintf(DatF,"BondBond Coeffs\n\n");
      for (i=0; i < no_angle_types; i++) {
        fprintf(DatF, "%3i ", i+1);
        for ( j = 0; j < 3; j++)
          fprintf(DatF, "%10.4f ", angletypes[i].bondbond_cross_term[j]);
        fputs("\n",DatF);
      }
      fputs("\n",DatF);

      fprintf(DatF,"BondAngle Coeffs\n\n");

      for (i=0; i < no_angle_types; i++) {
        fprintf(DatF, "%3i ", i+1);
        for ( j = 0; j < 4; j++)
          fprintf(DatF, "%10.4f ",angletypes[i].bondangle_cross_term[j]);
        fputs("\n",DatF);
      }
      fputs("\n",DatF);
    }

    if ((no_oop_types+no_angleangle_types) > 0) {
      fprintf(DatF,"AngleAngle Coeffs\n\n");
      for (i=0; i < no_oop_types; i++) {
        fprintf(DatF, "%3i ", i+1);
        for ( j = 0; j < 6; j++)
          fprintf(DatF, "%10.4f ", ooptypes[i].angleangle_params[j]);
        fputs("\n",DatF);
      }
      for (i=0; i < no_angleangle_types; i++) {
        fprintf(DatF, "%3i ", i+no_oop_types+1);
        for ( j = 0; j < 6; j++)
          fprintf(DatF, "%10.4f ", angleangletypes[i].params[j]);
        fputs("\n",DatF);
      }
      fputs("\n",DatF);
    }

    if (no_dihedral_types > 0) {
      fprintf(DatF,"AngleAngleTorsion Coeffs\n\n");
      for (i=0; i < no_dihedral_types; i++) {
        fprintf(DatF, "%3i ", i+1);
        for ( j = 0; j < 3; j++)
          fprintf(DatF,"%10.4f ",
                  dihedraltypes[i].angleangledihedral_cross_term[j]);
        fputs("\n",DatF);
      }
      fputs("\n",DatF);

      fprintf(DatF,"EndBondTorsion Coeffs\n\n");
      for (i=0; i < no_dihedral_types; i++) {
        fprintf(DatF, "%i ", i+1);
        for ( j = 0; j < 8; j++)
          fprintf(DatF, "%10.4f ",
                  dihedraltypes[i].endbonddihedral_cross_term[j]);
        fputs("\n",DatF);
      }
      fputs("\n",DatF);

      fprintf(DatF,"MiddleBondTorsion Coeffs\n\n");
      for (i=0; i < no_dihedral_types; i++) {
        fprintf(DatF, "%3i ", i+1);
        for ( j = 0; j < 4; j++)
          fprintf(DatF,"%10.4f ",
                  dihedraltypes[i].midbonddihedral_cross_term[j]);
        fputs("\n",DatF);
      }
      fputs("\n",DatF);


      fprintf(DatF,"BondBond13 Coeffs\n\n");
      for (i=0; i < no_dihedral_types; i++) {
        fprintf(DatF, "%3i ", i+1);
        for ( j = 0; j < 3; j++)
          fprintf(DatF, "%10.4f ",
                  dihedraltypes[i].bond13_cross_term[j]);
        fputs("\n",DatF);
      }
      fputs("\n",DatF);

      fprintf(DatF,"AngleTorsion Coeffs\n\n");
      for (i=0; i < no_dihedral_types; i++) {
        fprintf(DatF, "%3i ", i+1);
        for ( j = 0; j < 8; j++)
          fprintf(DatF, "%10.4f ",
                  dihedraltypes[i].angledihedral_cross_term[j]);
        fputs("\n",DatF);
      }
      fputs("\n",DatF);
    }
  }

  /*--------------------------------------------------------------------*/

  /* ATOMS */

  if (hintflag) fputs("Atoms # full\n\n",DatF);
  else fputs("Atoms\n\n",DatF);

  for(k=0; k < total_no_atoms; k++) {
    int typ = atoms[k].type;
    fprintf(DatF," %6i %6i %3i %9.6f %15.9f %15.9f %15.9f %3i %3i %3i",
            k+1,
            atoms[k].molecule,
            typ+1,
            atoms[k].q,
            atoms[k].x[0],
            atoms[k].x[1],
            atoms[k].x[2],
            atoms[k].image[0],
            atoms[k].image[1],
            atoms[k].image[2]);
    if (hintflag) fprintf(DatF," # %s\n",atomtypes[typ].potential);
    else fputs("\n",DatF);
  }
  fputs("\n",DatF);

  /***** BONDS *****/

  if (total_no_bonds > 0) {
    fprintf(DatF, "Bonds\n\n");
    for(k=0; k < total_no_bonds; k++)
      fprintf(DatF, "%6i %3i %6i %6i\n",k+1,
              bonds[k].type+1,
              bonds[k].members[0]+1,
              bonds[k].members[1]+1);
    fputs("\n",DatF);
  }

  /***** ANGLES *****/

  if (total_no_angles > 0) {
    fprintf(DatF, "Angles\n\n");
    for(k=0; k < total_no_angles; k++)
      fprintf(DatF, "%6i %3i %6i %6i %6i\n",k+1,
              angles[k].type+1,
              angles[k].members[0]+1,
              angles[k].members[1]+1,
              angles[k].members[2]+1);
    fputs("\n",DatF);
  }


  /***** TORSIONS *****/

  if (total_no_dihedrals > 0)   {
    fprintf(DatF,"Dihedrals\n\n");
    for(k=0; k < total_no_dihedrals; k++)
      fprintf(DatF, "%6i %3i %6i %6i %6i %6i\n",k+1,
              dihedrals[k].type+1,
              dihedrals[k].members[0]+1,
              dihedrals[k].members[1]+1,
              dihedrals[k].members[2]+1,
              dihedrals[k].members[3]+1);
    fputs("\n",DatF);
  }

  /***** OUT-OF-PLANES *****/

  if (total_no_oops+total_no_angle_angles > 0) {
    fprintf(DatF,"Impropers\n\n");
    for (k=0; k < total_no_oops; k++)
      fprintf(DatF, "%6i %3i %6i %6i %6i %6i \n", k+1,
              oops[k].type+1,
              oops[k].members[0]+1,
              oops[k].members[1]+1,
              oops[k].members[2]+1,
              oops[k].members[3]+1);
    if (forcefield & FF_TYPE_CLASS2) {
      for (k=0; k < total_no_angle_angles; k++)
        fprintf(DatF, "%6i %3i %6i %6i %6i %6i \n",k+total_no_oops+1,
                angleangles[k].type+no_oop_types+1,
                angleangles[k].members[0]+1,
                angleangles[k].members[1]+1,
                angleangles[k].members[2]+1,
                angleangles[k].members[3]+1);
    }
    fputs("\n",DatF);
  }

  /* Close data file */

  if (fclose(DatF) !=0) {
    printf("Error closing %s.lammps05\n", rootname);
    exit(61);
  }
}

