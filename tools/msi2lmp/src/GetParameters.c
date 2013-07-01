
#include "msi2lmp.h"
#include "Forcefield.h"

#include <string.h>
#include <stdlib.h>
#include <math.h>

static int find_improper_body_data(char [][5],struct FrcFieldItem,int *);
static void rearrange_improper(int,int);
static int find_trigonal_body_data(char [][5],struct FrcFieldItem);
static int find_angleangle_data(char [][5],struct FrcFieldItem,int[]);
static int find_match(int, char [][5],struct FrcFieldItem,int *);
static int match_types(int,int,char [][5],char [][5],int *);
static double get_r0(int,int);
static double get_t0(int,int,int);
static int quo_cp();
static void get_equivs(int,char [][5],char[][5]);
static int find_equiv_type(char[]);

static void condexit(int val)
{
    if (iflag == 0) exit(val);
}

/**********************************************************************/
/*                                                                    */
/*  GetParameters is a long routine for searching the forcefield      */
/*  parameters (read in by ReadFrcFile) for parameters corresponding  */
/*  to the different internal coordinate types derived by MakeLists   */
/*                                                                    */
/**********************************************************************/

void GetParameters(int Forcefield)
{
  int i,j,k,backwards,cp_type,rearrange;
  int kloc[3],multiplicity;
  char potential_types[4][5];
  char equiv_types[4][5];
  double rab,rbc,rcd,tabc,tbcd,tabd,tcbd;

  if (pflag > 1) fprintf(stderr," Trying Atom Equivalences if needed\n");

  /**********************************************************************/
  /*                                                                    */
  /*   Find masses  of atom types                                       */
  /*                                                                    */
  /**********************************************************************/

  for (i=0; i < no_atom_types; i++) {
    backwards = -1;
    strncpy(potential_types[0],atomtypes[i].potential,5);
    k = find_match(1,potential_types,ff_atomtypes,&backwards);
    if (k < 0) {
      printf(" Unable to find mass for %s\n",atomtypes[i].potential);
      condexit(10);
    } else {
      atomtypes[i].mass = ff_atomtypes.data[k].ff_param[0];
    }
  }

  /**********************************************************************/
  /*                                                                    */
  /*   Find VDW parameters for atom types                               */
  /*                                                                    */
  /**********************************************************************/

  for (i=0; i < no_atom_types; i++) {
    backwards = 0;
    for (j=0; j < 2; j++) atomtypes[i].params[j] = 0.0;
    strncpy(potential_types[0],atomtypes[i].potential,5);
    k = find_match(1,potential_types,ff_vdw,&backwards);
    if (k < 0) {
      get_equivs(1,potential_types,equiv_types);

      if (pflag > 2) printf("Using equivalences for VDW %s -> %s\n",
                            potential_types[0],equiv_types[0]);

      k = find_match(1,equiv_types,ff_vdw,&backwards);
    }
    if (k < 0) {
      printf(" Unable to find vdw data for %s\n",atomtypes[i].potential);
      condexit(11);
    } else {
      if (Forcefield == 1) {
        if((ff_vdw.data[k].ff_param[0] != 0.0 ) &&
           (ff_vdw.data[k].ff_param[1] != 0.0)) {
          atomtypes[i].params[0] =
            (ff_vdw.data[k].ff_param[1]*
             ff_vdw.data[k].ff_param[1])/(4.0*ff_vdw.data[k].ff_param[0]);
          atomtypes[i].params[1] = pow((ff_vdw.data[k].ff_param[0]/
                                        ff_vdw.data[k].ff_param[1]),
                                       (1.0/6.0));
        }
      } else {
        atomtypes[i].params[0] = ff_vdw.data[k].ff_param[1];
        atomtypes[i].params[1] = ff_vdw.data[k].ff_param[0];
      }
    }
  }

  if (pflag > 2) {
    printf("\n Atom Types, Masses and VDW Parameters\n");
    for (i=0; i < no_atom_types; i++) {
      printf(" %3s %8.4f %8.4f %8.4f\n",
             atomtypes[i].potential,atomtypes[i].mass, atomtypes[i].params[0],atomtypes[i].params[1]);
    }
  }

  /**********************************************************************/
  /*                                                                    */
  /*   Find parameters for bond types                                   */
  /*                                                                    */
  /**********************************************************************/

  for (i=0; i < no_bond_types; i++) {
    backwards = 0;
    for (j=0; j < 4; j++) bondtypes[i].params[j] = 0.0;
    for (j=0; j < 2; j++)
      strncpy(potential_types[j],
              atomtypes[bondtypes[i].types[j]].potential,5);
    k = find_match(2,potential_types,ff_bond,&backwards);
    if (k < 0) {
      get_equivs(2,potential_types,equiv_types);

      if (pflag > 2) {
        printf("Using equivalences for bond %s %s -> %s %s\n",
               potential_types[0],potential_types[1],
               equiv_types[0],equiv_types[1]);
      }
      k = find_match(2,equiv_types,ff_bond,&backwards);
    }
    if (k < 0) {
      printf(" Unable to find bond data for %s %s\n",
             potential_types[0],potential_types[1]);
      condexit(12);
    } else {
      if (Forcefield == 1) {
        bondtypes[i].params[0] = ff_bond.data[k].ff_param[1];
        bondtypes[i].params[1] = ff_bond.data[k].ff_param[0];
      } else {
        for (j=0; j < 4; j++)
          bondtypes[i].params[j] = ff_bond.data[k].ff_param[j];
      }
    }
  }

  if (pflag > 2) {
    printf("\n Bond Types and  Parameters\n");
    for (i=0; i < no_bond_types; i++) {
      for (j=0; j < 2; j++)
        printf("%-3s",atomtypes[bondtypes[i].types[j]].potential);
      for (j=0; j < 4; j++)
        printf(" %8.4f",bondtypes[i].params[j]);
      printf("\n");
    }
  }


  /**********************************************************************/
  /*                                                                    */
  /*   Find parameters for angle types including bondbond,              */
  /*   and bondangle parameters if Class II                             */
  /*                                                                    */
  /*   Each of the cross terms are searched separately even though      */
  /*   they share a given angle type. This allows parameters to be      */
  /*   in different order in the forcefield for each cross term or      */
  /*   maybe not even there.                                            */
  /*                                                                    */
  /**********************************************************************/
  for (i=0; i < no_angle_types; i++) {
    backwards = 0;
    for (j=0; j < 4; j++) angletypes[i].params[j] = 0.0;
    for (j=0; j < 3; j++)
      strncpy(potential_types[j],atomtypes[angletypes[i].types[j]].potential,5);
    k = find_match(3,potential_types,ff_ang,&backwards);
    if (k < 0) {
      get_equivs(3,potential_types,equiv_types);
      if (pflag > 2) {
        printf("Using equivalences for angle %s %s %s -> %s %s %s\n",
               potential_types[0],potential_types[1],
               potential_types[2],
               equiv_types[0],equiv_types[1],
               equiv_types[2]);
      }
      k = find_match(3,equiv_types,ff_ang,&backwards);
    }
    if (k < 0) {
      printf(" Unable to find angle data for %s %s %s\n",
             potential_types[0],potential_types[1],potential_types[2]);
      condexit(13);
    } else {
      if (Forcefield == 1) {
        angletypes[i].params[0] = ff_ang.data[k].ff_param[1];
        angletypes[i].params[1] = ff_ang.data[k].ff_param[0];
      } else {
        for (j=0; j < 4; j++)
          angletypes[i].params[j] = ff_ang.data[k].ff_param[j];
      }
    }
    if (Forcefield > 1) {
      get_equivs(3,potential_types,equiv_types);
      if (pflag > 2) {
        printf("Using equivalences for 3 body cross terms %s %s %s -> %s %s %s\n",
               potential_types[0],potential_types[1],potential_types[2],
               equiv_types[0],equiv_types[1],equiv_types[2]);
      }
      for (j=0; j < 3; j++) angletypes[i].bondbond_cross_term[j] = 0.0;
      for (j=0; j < 4; j++) angletypes[i].bondangle_cross_term[j] = 0.0;

      rab = get_r0(angletypes[i].types[0],angletypes[i].types[1]);
      rbc = get_r0(angletypes[i].types[1],angletypes[i].types[2]);

      angletypes[i].bondbond_cross_term[1] = rab;
      angletypes[i].bondbond_cross_term[2] = rbc;
      angletypes[i].bondangle_cross_term[2] = rab;
      angletypes[i].bondangle_cross_term[3] = rbc;

      k = find_match(3,potential_types,ff_bonbon,&backwards);
      if (k < 0) {
        k = find_match(3,equiv_types,ff_bonbon,&backwards);
      }
      if (k < 0) {
        printf(" Unable to find bondbond data for %s %s %s\n",
               potential_types[0],potential_types[1],potential_types[2]);
        condexit(14);
      } else {
        angletypes[i].bondbond_cross_term[0] = ff_bonbon.data[k].ff_param[0];
      }
      k = find_match(3,potential_types,ff_bonang,&backwards);
      if (k < 0) {
        k = find_match(3,equiv_types,ff_bonang,&backwards);
      }
      if (k < 0) {
        printf(" Unable to find bondangle data for %s %s %s\n",
               potential_types[0],potential_types[1],potential_types[2]);
        condexit(15);
      } else {
        if (backwards) {
          angletypes[i].bondangle_cross_term[0] = ff_bonang.data[k].ff_param[1];
          angletypes[i].bondangle_cross_term[1] = ff_bonang.data[k].ff_param[0];
        } else {
          angletypes[i].bondangle_cross_term[0] = ff_bonang.data[k].ff_param[0];
          angletypes[i].bondangle_cross_term[1] = ff_bonang.data[k].ff_param[1];
        }
      }
    }
  }

  if (pflag > 2) {
    printf("\n Angle Types and Parameters\n");
    for (i=0; i < no_angle_types; i++) {
      for (j=0; j < 3; j++)
        printf(" %-3s", atomtypes[angletypes[i].types[j]].potential);
      for (j=0; j < 4; j++) printf(" %8.4f",angletypes[i].params[j]);
      printf("\n");
    }

    if (forcefield > 1) {
      printf("\n BondBond Types and  Parameters\n");
      for (i=0; i < no_angle_types; i++) {
        for (j=0; j < 3; j++)
          printf("%-3s",atomtypes[angletypes[i].types[j]].potential);
        for (j=0; j < 3; j++)
          printf(" %8.4f",angletypes[i].bondbond_cross_term[j]);
        printf("\n");
      }
      printf("\n BondAngle Types and  Parameters\n");
      for (i=0; i < no_angle_types; i++) {
        for (j=0; j < 3; j++)
          printf("%-3s",atomtypes[angletypes[i].types[j]].potential);
        for (j=0; j < 4; j++)
          printf(" %8.4f",angletypes[i].bondangle_cross_term[j]);
        printf("\n");
      }
    }
  }

  /**********************************************************************/
  /*                                                                    */
  /*   Find parameters for dihedral types including endbonddihedral,    */
  /*   midbonddihedral, angledihedral, angleangledihedral and           */
  /*   bondbond13 parameters if Class II                                */
  /*                                                                    */
  /*   Each of the cross terms are searched separately even though      */
  /*   they share a given dihedral type. This allows parameters to be   */
  /*   in different order in the forcefield for each cross term or      */
  /*   maybe not even there.                                            */
  /*                                                                    */
  /**********************************************************************/

  for (i=0; i < no_dihedral_types; i++) {
    for (j=0; j < 6; j++)
      dihedraltypes[i].params[j] = 0.0;
    for (j=0; j < 4; j++)
      strncpy(potential_types[j],
              atomtypes[dihedraltypes[i].types[j]].potential,5);
    backwards = 0;
    k = find_match(4,potential_types,ff_tor,&backwards);

    if (k < 0) {
      get_equivs(4,potential_types,equiv_types);

      if (pflag > 2) {
        printf("Using equivalences for dihedral %s %s %s %s -> %s %s %s %s\n",
               potential_types[0],potential_types[1],
               potential_types[2],potential_types[3],
               equiv_types[0],equiv_types[1],
               equiv_types[2],equiv_types[3]);
      }
      k = find_match(4,equiv_types,ff_tor,&backwards);
    }
    if (k < 0) {
      printf(" Unable to find torsion data for %s %s %s %s\n",
             potential_types[0],
             potential_types[1],
             potential_types[2],
             potential_types[3]);
      condexit(16);
    } else {
      if (Forcefield == 1) {
        multiplicity = 1;
        if (ff_tor.data[k].ff_types[0][0] == '*')
          multiplicity =
            atomtypes[dihedraltypes[i].types[1]].no_connect-1;
        if (ff_tor.data[k].ff_types[3][0] == '*')
          multiplicity *=
            atomtypes[dihedraltypes[i].types[2]].no_connect-1;

        dihedraltypes[i].params[0] = ff_tor.data[k].ff_param[0]/(double) multiplicity;
        if (ff_tor.data[k].ff_param[2] == 0.0)
          dihedraltypes[i].params[1] = 1.0;
        else if (ff_tor.data[k].ff_param[2] == 180.0)
          dihedraltypes[i].params[1] = -1.0;
        else {
          printf("Non planar phi0 for %s %s %s %s\n",
                 potential_types[0],potential_types[1],
                 potential_types[2],potential_types[3]);
          dihedraltypes[i].params[1] = 0.0;
        }
        dihedraltypes[i].params[2] = ff_tor.data[k].ff_param[1];
      }
      else {
        for (j=0; j < 6; j++)
          dihedraltypes[i].params[j] = ff_tor.data[k].ff_param[j];
      }
    }

    if (Forcefield > 1) {
      get_equivs(4,potential_types,equiv_types);
      if (pflag > 2) {
        printf("Using equivalences for linear 4 body cross terms  %s %s %s %s -> %s %s %s %s\n",
               potential_types[0],potential_types[1],
               potential_types[2],potential_types[3],
               equiv_types[0],equiv_types[1],
               equiv_types[2],equiv_types[3]);
      }

      for (j=0; j < 8; j++)
        dihedraltypes[i].endbonddihedral_cross_term[j] = 0.0;
      for (j=0; j < 4; j++)
        dihedraltypes[i].midbonddihedral_cross_term[j] = 0.0;
      for (j=0; j < 8; j++)
        dihedraltypes[i].angledihedral_cross_term[j] = 0.0;
      for (j=0; j < 3; j++)
        dihedraltypes[i].angleangledihedral_cross_term[j] = 0.0;
      for (j=0; j < 3; j++)
        dihedraltypes[i].bond13_cross_term[j] = 0.0;

      rab = get_r0(dihedraltypes[i].types[0],dihedraltypes[i].types[1]);
      rbc = get_r0(dihedraltypes[i].types[1],dihedraltypes[i].types[2]);
      rcd = get_r0(dihedraltypes[i].types[2],dihedraltypes[i].types[3]);
      tabc = get_t0(dihedraltypes[i].types[0],
                    dihedraltypes[i].types[1],
                    dihedraltypes[i].types[2]);

      tbcd = get_t0(dihedraltypes[i].types[1],
                    dihedraltypes[i].types[2],
                    dihedraltypes[i].types[3]);

      dihedraltypes[i].endbonddihedral_cross_term[6] = rab;
      dihedraltypes[i].endbonddihedral_cross_term[7] = rcd;
      dihedraltypes[i].midbonddihedral_cross_term[3] = rbc;
      dihedraltypes[i].angledihedral_cross_term[6] = tabc;
      dihedraltypes[i].angledihedral_cross_term[7] = tbcd;
      dihedraltypes[i].angleangledihedral_cross_term[1] = tabc;
      dihedraltypes[i].angleangledihedral_cross_term[2] = tbcd;
      dihedraltypes[i].bond13_cross_term[1] = rab;
      dihedraltypes[i].bond13_cross_term[2] = rcd;

      backwards = 0;
      k = find_match(4,potential_types,ff_endbontor,&backwards);
      if (k < 0) {
        k = find_match(4,equiv_types,ff_endbontor,&backwards);
      }
      if (k < 0) {
        printf(" Unable to find endbonddihedral data for %s %s %s %s\n",
               potential_types[0],potential_types[1],
               potential_types[2],potential_types[3]);
        condexit(17);
      } else {
        if (backwards) {
          dihedraltypes[i].endbonddihedral_cross_term[0] =
            ff_endbontor.data[k].ff_param[3];
          dihedraltypes[i].endbonddihedral_cross_term[1] =
            ff_endbontor.data[k].ff_param[4];
          dihedraltypes[i].endbonddihedral_cross_term[2] =
            ff_endbontor.data[k].ff_param[5];
          dihedraltypes[i].endbonddihedral_cross_term[3] =
            ff_endbontor.data[k].ff_param[0];
          dihedraltypes[i].endbonddihedral_cross_term[4] =
            ff_endbontor.data[k].ff_param[1];
          dihedraltypes[i].endbonddihedral_cross_term[5] =
            ff_endbontor.data[k].ff_param[2];
        }
        else {
          dihedraltypes[i].endbonddihedral_cross_term[0] =
            ff_endbontor.data[k].ff_param[0];
          dihedraltypes[i].endbonddihedral_cross_term[1] =
            ff_endbontor.data[k].ff_param[1];
          dihedraltypes[i].endbonddihedral_cross_term[2] =
            ff_endbontor.data[k].ff_param[2];
          dihedraltypes[i].endbonddihedral_cross_term[3] =
            ff_endbontor.data[k].ff_param[3];
          dihedraltypes[i].endbonddihedral_cross_term[4] =
            ff_endbontor.data[k].ff_param[4];
          dihedraltypes[i].endbonddihedral_cross_term[5] =
            ff_endbontor.data[k].ff_param[5];
        }
      }
      backwards = 0;
      k = find_match(4,potential_types,ff_midbontor,&backwards);
      if (k < 0) {
        k = find_match(4,equiv_types,ff_midbontor,&backwards);
      }
      if (k < 0) {
        printf(" Unable to find midbonddihedral data for %s %s %s %s\n",
               potential_types[0],potential_types[1],
               potential_types[2],potential_types[3]);
        condexit(18);
      } else {
        dihedraltypes[i].midbonddihedral_cross_term[0] =
          ff_midbontor.data[k].ff_param[0];
        dihedraltypes[i].midbonddihedral_cross_term[1] =
          ff_midbontor.data[k].ff_param[1];
        dihedraltypes[i].midbonddihedral_cross_term[2] =
          ff_midbontor.data[k].ff_param[2];
      }

      backwards = 0;
      k = find_match(4,potential_types,ff_angtor,&backwards);
      if (k < 0) {
        k = find_match(4,equiv_types,ff_angtor,&backwards);
      }
      if (k < 0) {
        printf(" Unable to find angledihedral data for %s %s %s %s\n",
               potential_types[0],potential_types[1],
               potential_types[2],potential_types[3]);
        condexit(19);
      } else {
        if (backwards) {
          dihedraltypes[i].angledihedral_cross_term[0] =
            ff_angtor.data[k].ff_param[3];
          dihedraltypes[i].angledihedral_cross_term[1] =
            ff_angtor.data[k].ff_param[4];
          dihedraltypes[i].angledihedral_cross_term[2] =
            ff_angtor.data[k].ff_param[5];
          dihedraltypes[i].angledihedral_cross_term[3] =
            ff_angtor.data[k].ff_param[0];
          dihedraltypes[i].angledihedral_cross_term[4] =
            ff_angtor.data[k].ff_param[1];
          dihedraltypes[i].angledihedral_cross_term[5] =
            ff_angtor.data[k].ff_param[2];
        }
        else {
          dihedraltypes[i].angledihedral_cross_term[0] =
            ff_angtor.data[k].ff_param[0];
          dihedraltypes[i].angledihedral_cross_term[1] =
            ff_angtor.data[k].ff_param[1];
          dihedraltypes[i].angledihedral_cross_term[2] =
            ff_angtor.data[k].ff_param[2];
          dihedraltypes[i].angledihedral_cross_term[3] =
            ff_angtor.data[k].ff_param[3];
          dihedraltypes[i].angledihedral_cross_term[4] =
            ff_angtor.data[k].ff_param[4];
          dihedraltypes[i].angledihedral_cross_term[5] =
            ff_angtor.data[k].ff_param[5];
        }
      }
      backwards = 0;
      k = find_match(4,potential_types,ff_angangtor,&backwards);
      if (k < 0) {
        k = find_match(4,equiv_types,ff_angangtor,&backwards);
      }
      if (k < 0) {
        printf(" Unable to find angleangledihedral data for %s %s %s %s\n",
               potential_types[0],potential_types[1],
               potential_types[2],potential_types[3]);
        condexit(20);
      } else {
        dihedraltypes[i].angleangledihedral_cross_term[0] =
          ff_angangtor.data[k].ff_param[0];
      }
      cp_type = quo_cp();
      if ((cp_type >= 0) &&
          ((dihedraltypes[i].types[0] == cp_type) ||
           (dihedraltypes[i].types[1] == cp_type) ||
           (dihedraltypes[i].types[2] == cp_type) ||
           (dihedraltypes[i].types[3] == cp_type)   )) {
        backwards = 0;
        k = find_match(4,potential_types,ff_bonbon13,&backwards);
        if (k < 0) {
          k = find_match(4,equiv_types,ff_bonbon13,&backwards);
        }
        if (k < 0) {
          printf(" Unable to find bond13 data for %s %s %s %s\n",
                 potential_types[0],potential_types[1],
                 potential_types[2],potential_types[3]);
          condexit(21);
        } else {
          dihedraltypes[i].bond13_cross_term[0] =
            ff_bonbon13.data[k].ff_param[0];
        }
      }
    }
  }

  if (pflag > 2) {
    printf("\n Dihedral Types and  Parameters\n");
    for (i=0; i < no_dihedral_types; i++) {
      for (j=0; j < 4; j++)
        printf("%-3s",atomtypes[dihedraltypes[i].types[j]].potential);
      for (j=0; j < 6; j++)
        printf(" %8.4f",dihedraltypes[i].params[j]);
      printf("\n");
    }

    if (forcefield > 1) {

      printf("\n EndBondDihedral Types and Parameters\n");
      for (i=0; i < no_dihedral_types; i++) {
        for (j=0; j < 4; j++)
          printf("%-3s",atomtypes[dihedraltypes[i].types[j]].potential);
        for (j=0; j < 8; j++)
          printf(" %8.4f",dihedraltypes[i].endbonddihedral_cross_term[j]);
        printf("\n");
      }
      printf("\n MidBondDihedral Types and Parameters\n");
      for (i=0; i < no_dihedral_types; i++) {
        for (j=0; j < 4; j++)
          printf(" %-3s",atomtypes[dihedraltypes[i].types[j]].potential);
        for (j=0; j < 4; j++)
          printf(" %8.4f",dihedraltypes[i].midbonddihedral_cross_term[j]);
        printf("\n");
      }

      printf("\n AngleDihedral Types and Parameters\n");
      for (i=0; i < no_dihedral_types; i++) {
        for (j=0; j < 4; j++)
          printf("%-3s",atomtypes[dihedraltypes[i].types[j]].potential);
        for (j=0; j < 8; j++)
          printf(" %8.4f",dihedraltypes[i].angledihedral_cross_term[j]);
        printf("\n");
      }

      printf("\n AngleAngleDihedral Types and Parameters\n");
      for (i=0; i < no_dihedral_types; i++) {
        for (j=0; j < 4; j++)
          printf(" %-3s",atomtypes[dihedraltypes[i].types[j]].potential);
        for (j=0; j < 3; j++)
          printf("%8.4f",dihedraltypes[i].angleangledihedral_cross_term[j]);
        printf("\n");
      }

      printf("\n Bond13 Types and  Parameters\n");

      for (i=0; i < no_dihedral_types; i++) {
        for (j=0; j < 4; j++)
          printf(" %-3s",atomtypes[dihedraltypes[i].types[j]].potential);
        for (j=0; j < 3; j++)
          printf(" %8.4f",dihedraltypes[i].bond13_cross_term[j]);
        printf("\n");
      }
    }
  }


  /**********************************************************************/
  /*                                                                    */
  /*   Find parameters for oop types                                    */
  /*                                                                    */
  /*   This is the most complicated of all the types because the        */
  /*   the class I oop is actually an improper torsion and does         */
  /*   not have the permutation symmetry of a well defined oop          */
  /*   The net result is that if one does not find the current          */
  /*   atom type ordering in the Forcefield file then one must try each */
  /*   of the next permutations (6 in total) and when a match is found  */
  /*   the program must go back and rearrange the oop type AND the atom */
  /*   ordering in the oop lists for those with the current type        */
  /*                                                                    */
  /*   The Class II oop types are easier but also tedious since the     */
  /*   program has to try all permutations of the a c and d atom        */
  /*   types to find a match. A special routine is used to do this.     */
  /*                                                                    */
  /*   Fortunately, there are typically few oop types                   */
  /*                                                                    */
  /**********************************************************************/

  if (forcefield == 1) {
    for (i=0; i < no_oop_types; i++) {
      for (j=0; j < 3; j++) ooptypes[i].params[j] = 0.0;
      for (j=0; j < 4; j++)
        strncpy(potential_types[j],
                atomtypes[ooptypes[i].types[j]].potential,5);

      k = find_improper_body_data(potential_types,ff_oop,&rearrange);
      if (k < 0) {
        get_equivs(5,potential_types,equiv_types);

        if (pflag > 2) {
          printf("Using equivalences for oop %s %s %s %s -> %s %s %s %s\n",
                 potential_types[0],potential_types[1],
                 potential_types[2],potential_types[3],
                 equiv_types[0],equiv_types[1],
                 equiv_types[2],equiv_types[3]);
        }
        k = find_improper_body_data(equiv_types,ff_oop,&rearrange);
      }
      if (k < 0) {
        printf(" Unable to find oop data for %s %s %s %s\n",
               potential_types[0],
               potential_types[1],potential_types[2],potential_types[3]);
        condexit(22);
      } else {
        ooptypes[i].params[0] = ff_oop.data[k].ff_param[0];
        if (ff_oop.data[k].ff_param[2] == 0.0)
          ooptypes[i].params[1] = 1.0;
        else if (ff_oop.data[k].ff_param[2] == 180.0)
          ooptypes[i].params[1] = -1.0;
        else {
          printf("Non planar phi0 for %s %s %s %s\n",
                 potential_types[0],potential_types[1],
                 potential_types[2],potential_types[3]);
          ooptypes[i].params[1] = 0.0;
        }
        ooptypes[i].params[2] = ff_oop.data[k].ff_param[1];
        if (rearrange > 0) rearrange_improper(i,rearrange);
      }
    }
  } else {
    for (i=0; i < no_oop_types; i++) {
      for (j=0; j < 3; j++)
        ooptypes[i].params[j] = 0.0;
      for (j=0; j < 4; j++)
        strncpy(potential_types[j],
                atomtypes[ooptypes[i].types[j]].potential,5);
      k = find_trigonal_body_data(potential_types,ff_oop);
      if (k < 0) {
        get_equivs(5,potential_types,equiv_types);
        if (pflag > 2) {
          printf("Using equivalences for oop %s %s %s %s -> %s %s %s %s\n",
                 potential_types[0],potential_types[1],
                 potential_types[2],potential_types[3],
                 equiv_types[0],equiv_types[1],
                 equiv_types[2],equiv_types[3]);
        }
        k = find_trigonal_body_data(equiv_types,ff_oop);
      }
      if (k < 0) {
        printf(" Unable to find oop data for %s %s %s %s\n",
               potential_types[0],
               potential_types[1],potential_types[2],potential_types[3]);
        condexit(23);
      } else {
        for (j=0; j < 2; j++)
          ooptypes[i].params[j] = ff_oop.data[k].ff_param[j];
      }
    }
  }

  if (pflag > 2) {
    printf("\n OOP Types and  Parameters\n");
    for (i=0; i < no_oop_types; i++) {
      for (j=0; j < 4; j++)
        printf("%-3s",atomtypes[ooptypes[i].types[j]].potential);
      for (j=0; j < 3; j++)
        printf(" %8.4f",ooptypes[i].params[j]);
      printf("\n");
    }
  }


  /**********************************************************************/
  /*                                                                    */
  /*   Find parameters for angleangle types (Class II only)             */
  /*                                                                    */
  /*   This is somewhat complicated in that one set of four types       */
  /*   a b c d has three angleangle combinations so for each type       */
  /*   the program needs to find three sets of parameters by            */
  /*   progressively looking for data for different permutations of     */
  /*   a c and d                                                        */
  /*                                                                    */
  /**********************************************************************/

  if (forcefield > 1) {

    for (i=0; i < no_oop_types; i++) {

      for (j=0; j < 6; j++) ooptypes[i].angleangle_params[j] = 0.0;

      for (j=0; j < 4; j++)
        strncpy(potential_types[j],
                atomtypes[ooptypes[i].types[j]].potential,5);


      tabc = get_t0(ooptypes[i].types[0],
                    ooptypes[i].types[1],
                    ooptypes[i].types[2]);

      tabd = get_t0(ooptypes[i].types[0],
                    ooptypes[i].types[1],
                    ooptypes[i].types[3]);
      tcbd = get_t0(ooptypes[i].types[2],
                    ooptypes[i].types[1],

                    ooptypes[i].types[3]);

      ooptypes[i].angleangle_params[3] = tabc;
      ooptypes[i].angleangle_params[4] = tcbd;
      ooptypes[i].angleangle_params[5] = tabd;

      k = find_angleangle_data(potential_types,ff_angang,kloc);
      if (k < 0) {
        get_equivs(5,potential_types,equiv_types);
        if (pflag > 2) {
          printf("Using equivalences for angleangle %s %s %s %s -> %s %s %s %s\n",
                 potential_types[0],potential_types[1],
                 potential_types[2],potential_types[3],
                 equiv_types[0],equiv_types[1],
                 equiv_types[2],equiv_types[3]);
          k = find_angleangle_data(equiv_types,ff_angang,kloc);
        }
      }
      if (k < 0) {
        printf(" Unable to find angleangle data for %s %s %s %s\n",
               potential_types[0],
               potential_types[1],potential_types[2],potential_types[3]);
        condexit(24);
      } else {
        for (j=0; j < 3; j++) {
          if (kloc[j] > -1)
            ooptypes[i].angleangle_params[j] = ff_angang.data[kloc[j]].ff_param[0];
        }
      }
    }

    for (i=0; i < no_angleangle_types; i++) {
      for (j=0; j < 6; j++) angleangletypes[i].params[j] = 0.0;
      for (j=0; j < 4; j++)
        strncpy(potential_types[j],
                atomtypes[angleangletypes[i].types[j]].potential,5);

      tabc = get_t0(angleangletypes[i].types[0],
                    angleangletypes[i].types[1],
                    angleangletypes[i].types[2]);
      tabd = get_t0(angleangletypes[i].types[0],
                    angleangletypes[i].types[1],
                    angleangletypes[i].types[3]);
      tcbd = get_t0(angleangletypes[i].types[2],
                    angleangletypes[i].types[1],
                    angleangletypes[i].types[3]);

      angleangletypes[i].params[3] = tabc;
      angleangletypes[i].params[4] = tcbd;
      angleangletypes[i].params[5] = tabd;

      k = find_angleangle_data(potential_types,ff_angang,kloc);
      if (k < 0) {
        get_equivs(5,potential_types,equiv_types);
        if (pflag > 2) {
          printf("Using equivalences for angleangle %s %s %s %s -> %s %s %s %s\n",
                 potential_types[0],potential_types[1],
                 potential_types[2],potential_types[3],
                 equiv_types[0],equiv_types[1],
                 equiv_types[2],equiv_types[3]);
        }
        k = find_angleangle_data(equiv_types,ff_angang,kloc);
      }
      if (k < 0) {
        printf(" Unable to find angleangle data for %s %s %s %s\n",
               potential_types[0],
               potential_types[1],potential_types[2],potential_types[3]);
        condexit(25);
      } else {
        for (j=0; j < 3; j++) {
          if (kloc[j] > -1)
            angleangletypes[i].params[j] =
              ff_angang.data[kloc[j]].ff_param[0];
        }
      }
    }
    if (pflag > 2) {
      printf("\n AngleAngle Types and  Parameters\n");
      for (i=0; i < no_oop_types; i++) {
        for (j=0; j < 4; j++)
          printf("%-3s",atomtypes[ooptypes[i].types[j]].potential);
        for (j=0; j < 6; j++)
          printf(" %8.4f",ooptypes[i].angleangle_params[j]);
        printf("\n");
      }
      for (i=0; i < no_angleangle_types; i++) {
        for (j=0; j < 4; j++)
          printf(" %-3s",atomtypes[angleangletypes[i].types[j]].potential);
        for (j=0; j < 6; j++) printf(" %8.4f",angleangletypes[i].params[j]);
        printf("\n");
      }
    }
  }
}

int find_improper_body_data(char types1[][5],struct FrcFieldItem item,
                            int *rearrange_ptr)
{
  int k,backwards;
  char mirror_types[4][5];

  backwards = 0;

  /* a b c d */

  *rearrange_ptr = 0;
  k = find_match(4,types1,item,&backwards);
  if (k >= 0) return k;

  /* a b d c */

  *rearrange_ptr = 1;
  strncpy(mirror_types[0],types1[0],5);
  strncpy(mirror_types[1],types1[1],5);
  strncpy(mirror_types[2],types1[3],5);
  strncpy(mirror_types[3],types1[2],5);
  k = find_match(4,mirror_types,item,&backwards);
  if (k >= 0) return k;

  /* d b a c */

  *rearrange_ptr = 2;
  strncpy(mirror_types[0],types1[3],5);
  strncpy(mirror_types[2],types1[0],5);
  strncpy(mirror_types[3],types1[2],5);
  k = find_match(4,mirror_types,item,&backwards);
  if (k >= 0) return k;

  /* d b c a */

  *rearrange_ptr = 3;
  strncpy(mirror_types[2],types1[2],5);
  strncpy(mirror_types[3],types1[0],5);
  k = find_match(4,mirror_types,item,&backwards);
  if (k >= 0) return k;

  /* c b a d */

  *rearrange_ptr = 4;
  strncpy(mirror_types[0],types1[2],5);
  strncpy(mirror_types[2],types1[0],5);
  strncpy(mirror_types[3],types1[3],5);
  k = find_match(4,mirror_types,item,&backwards);
  if (k >= 0) return k;

  /* c b d a */

  *rearrange_ptr = 5;
  strncpy(mirror_types[2],types1[3],5);
  strncpy(mirror_types[3],types1[0],5);
  k = find_match(4,mirror_types,item,&backwards);
  return k;
}

void rearrange_improper(int ooptype,int rearrange)
{
  int i,j,temp[4];

  for (i=0; i < 4; i++) temp[i] = ooptypes[ooptype].types[i];

  switch (rearrange) {
  case 1:
    ooptypes[ooptype].types[0] = temp[0];
    ooptypes[ooptype].types[2] = temp[3];
    ooptypes[ooptype].types[3] = temp[2];
    for (i=0; i < total_no_oops; i++) {
      if (oops[i].type == ooptype) {
        for (j=0; j < 4; j++) temp[j] = oops[i].members[j];
        oops[i].members[2] = temp[3];
        oops[i].members[3] = temp[2];
      }
    }
    break;
  case 2:
    ooptypes[ooptype].types[0] = temp[3];
    ooptypes[ooptype].types[2] = temp[0];
    ooptypes[ooptype].types[3] = temp[2];
    for (i=0; i < total_no_oops; i++) {
      if (oops[i].type == ooptype) {
        for (j=0; j < 4; j++) temp[j] = oops[i].members[j];
        oops[i].members[0] = temp[3];
        oops[i].members[2] = temp[0];
        oops[i].members[3] = temp[2];
      }
    }
    break;
  case 3:
    ooptypes[ooptype].types[0] = temp[3];
    ooptypes[ooptype].types[2] = temp[2];
    ooptypes[ooptype].types[3] = temp[0];
    for (i=0; i < total_no_oops; i++) {
      if (oops[i].type == ooptype) {
        for (j=0; j < 4; j++) temp[j] = oops[i].members[j];
        oops[i].members[0] = temp[3];
        oops[i].members[2] = temp[2];
        oops[i].members[3] = temp[0];
      }
    }
    break;
  case 4:
    ooptypes[ooptype].types[0] = temp[2];
    ooptypes[ooptype].types[2] = temp[0];
    ooptypes[ooptype].types[3] = temp[3];
    for (i=0; i < total_no_oops; i++) {
      if (oops[i].type == ooptype) {
        for (j=0; j < 4; j++) temp[j] = oops[i].members[j];
        oops[i].members[0] = temp[2];
        oops[i].members[2] = temp[0];
        oops[i].members[3] = temp[3];
      }
    }
    break;
  case 5:
    ooptypes[ooptype].types[0] = temp[2];
    ooptypes[ooptype].types[2] = temp[3];
    ooptypes[ooptype].types[3] = temp[0];
    for (i=0; i < total_no_oops; i++) {
      if (oops[i].type == ooptype) {
        for (j=0; j < 4; j++) temp[j] = oops[i].members[j];
        oops[i].members[0] = temp[2];
        oops[i].members[2] = temp[3];
        oops[i].members[3] = temp[0];
      }
    }
    break;
  default:
    break;
  }
}

int find_trigonal_body_data(char types1[][5],struct FrcFieldItem item)
{
  int k,backwards;
  char mirror_types[4][5];

  backwards = -1;

  /* a b c d */

  k = find_match(4,types1,item,&backwards);
  if (k >= 0) return k;

  /* a b d c */

  strncpy(mirror_types[0],types1[0],5);
  strncpy(mirror_types[1],types1[1],5);
  strncpy(mirror_types[2],types1[3],5);
  strncpy(mirror_types[3],types1[2],5);
  k = find_match(4,mirror_types,item,&backwards);
  if (k >= 0) return k;

  /* d b a c */

  strncpy(mirror_types[0],types1[3],5);
  strncpy(mirror_types[2],types1[0],5);
  strncpy(mirror_types[3],types1[2],5);
  k = find_match(4,mirror_types,item,&backwards);
  if (k >= 0) return k;

  /* d b c a */

  strncpy(mirror_types[2],types1[2],5);
  strncpy(mirror_types[3],types1[0],5);
  k = find_match(4,mirror_types,item,&backwards);
  if (k >= 0) return k;
  /* c b a d */

  strncpy(mirror_types[0],types1[2],5);
  strncpy(mirror_types[2],types1[0],5);
  strncpy(mirror_types[3],types1[3],5);
  k = find_match(4,mirror_types,item,&backwards);
  if (k >= 0) return k;

  /* c b d a */

  strncpy(mirror_types[2],types1[3],5);
  strncpy(mirror_types[3],types1[0],5);
  k = find_match(4,mirror_types,item,&backwards);
  return k;
}

int find_angleangle_data(char types1[][5],struct FrcFieldItem item,int kloc[3])
{
  int k,backwards = -1;
  char mirror_types[4][5];

  strncpy(mirror_types[1],types1[1],5);

  /* go for first parameter a b c d or d b c a */

  k = find_match(4,types1,item,&backwards);
  if (k < 0) {
    strncpy(mirror_types[0],types1[3],5);
    strncpy(mirror_types[2],types1[2],5);
    strncpy(mirror_types[3],types1[0],5);
    k = find_match(4,mirror_types,item,&backwards);
  }
  kloc[0] = k;

  /* go for second parameter d b a c or c b a d */

  strncpy(mirror_types[0],types1[3],5);
  strncpy(mirror_types[2],types1[0],5);
  strncpy(mirror_types[3],types1[2],5);
  k = find_match(4,mirror_types,item,&backwards);
  if (k < 0) {
    strncpy(mirror_types[0],types1[2],5);
    strncpy(mirror_types[3],types1[3],5);
    k = find_match(4,mirror_types,item,&backwards);
  }
  kloc[1] = k;

  /* go for third parameter a b d c or c b d a */

  strncpy(mirror_types[0],types1[0],5);
  strncpy(mirror_types[2],types1[3],5);
  strncpy(mirror_types[3],types1[2],5);
  k = find_match(4,mirror_types,item,&backwards);
  if (k < 0) {
    strncpy(mirror_types[0],types1[2],5);
    strncpy(mirror_types[3],types1[0],5);
    k = find_match(4,mirror_types,item,&backwards);
  }
  kloc[2] = k;
  k = 0;
  if ((kloc[0] < 0) && (kloc[1] < 0) && (kloc[2] < 0)) k = -1;
  return k;
}

int find_match(int n, char types1[][5],struct FrcFieldItem item,int
               *backwards_ptr)
{
  int k,match;

  match = 0;
  k=0;

  /* Try for an exact match (no wildcards) first */

  while (!match && (k < item.entries)) {
    if (match_types(n, 0,types1,item.data[k].ff_types,backwards_ptr) == 1)
      match = 1;
    else
      k++;
  }

  /* Try again - allow wildcard matching  */

  if (!match) {
    k=0;
    while (!match && (k < item.entries)) {
      if (match_types(n,1,types1,item.data[k].ff_types,backwards_ptr) == 1)
        match = 1;
      else
        k++;
    }
  }
  if (match) return k;
  else return -1;
}

int match_types(int n,int wildcard,char types1[][5],char types2[][5],
                int *backwards_ptr)
{
  int k,match;

  /* Routine to match short arrays of characters strings which contain
     atom potential types. The arrays range from 1 to 4 (VDW or equivalences,
     bond, angle, dihedrals or oops). There are potentially four ways the
     arrays can match: exact match (forwards), exact match when one array is
     run backwards (backwards), forwards with wildcard character match allowed
     (forwards *) and finally backwards with wildcard character match
     (backwards *). If the variable, backwards (pointed by backwards_ptr)
     is -1, then the backwards options are not to be used (such when
     matching oop types)
  */


  if (wildcard == 0) {

  /* forwards */

    k=0;
    match = 1;
    while (match && (k < n)) {
      if (strncmp(types1[k],types2[k],5) == 0)
        k++;
      else
        match = 0;
    }
  }
  else {

  /* forwards * */

    k=0;

    match = 1;
    while (match && (k < n)) {
      if ((strncmp(types1[k],types2[k],5) == 0) ||
          (types2[k][0] == '*'))
        k++;
      else
        match = 0;
    }
  }

  if (match) {
    *backwards_ptr = 0;
    return 1;
  }
  if ((n < 2) || (*backwards_ptr == -1)) return 0;

  if (wildcard == 0) {

  /* backwards */

    k=0;
    match = 1;
    while (match && (k < n)) {
      if (strncmp(types1[n-k-1],types2[k],5) == 0)
        k++;
      else
        match = 0;
    }
  }
  else {

  /* backwards * */

    k=0;
    match = 1;
    while (match && (k < n)) {
      if ((strncmp(types1[n-k-1],types2[k],5) == 0) ||
          (types2[k][0] == '*')                   )
        k++;
      else
        match = 0;
    }
  }

  if (match) {
    *backwards_ptr = 1;
    return 1;
  }
  else return 0;
}

double get_r0(int typei,int typej)
{
  int k,match;
  double r;

  k=0;
  match=0;
  r = 0.0;

  while (!match && (k < no_bond_types)) {
    if (((typei == bondtypes[k].types[0]) &&
         (typej == bondtypes[k].types[1])) ||
        ((typej == bondtypes[k].types[0]) &&
         (typei == bondtypes[k].types[1]))   ) {
      r = bondtypes[k].params[0];
      match = 1;
    }
    else
      k++;
  }

  if (match == 0)
    printf("Unable to find r0 for types %d %d\n",typei,typej);
  return r;
}

double get_t0(int typei,int typej,int typek)
{
  int k,match;
  double theta;

  k=0;
  match=0;
  theta = 0.0;

  while (!match && (k < no_angle_types)) {
    if (((typei == angletypes[k].types[0]) &&
         (typej == angletypes[k].types[1]) &&
         (typek == angletypes[k].types[2])) ||
        ((typek == angletypes[k].types[0]) &&
         (typej == angletypes[k].types[1]) &&
         (typei == angletypes[k].types[2]))   ) {
      theta = angletypes[k].params[0];
      match = 1;
    }
    else
      k++;
  }

  if (match == 0)
    printf(" Unable to find t0 for types %d %d %d\n",
           typei,typej,typek);
  return theta;
}

int quo_cp()
{
  char cp[] = "cp  ";
  int i,type,found;

  i = 0;
  type = -1;
  found = 0;

  while (!found && (i < no_atom_types)) {
    if (strncmp(atomtypes[i].potential,cp,2) == 0) {
      found = 1;
      type = i;
    }
    else
      i++;
  }

  return type;
}

void get_equivs(int ic,char potential_types[][5],char equiv_types[][5])
{
  int i,k;
  switch (ic) {
  case 1:
    k = find_equiv_type(potential_types[0]);
    if (k > -1) strncpy(equiv_types[0],equivalence.data[k].ff_types[1],5);
    break;

  case 2:
    for (i=0; i < 2; i++) {
      k = find_equiv_type(potential_types[i]);
      if (k > -1) strncpy(equiv_types[i],equivalence.data[k].ff_types[2],5);
    }
    break;
  case 3:
    for (i=0; i < 3; i++) {
      k = find_equiv_type(potential_types[i]);
      if (k > -1) strncpy(equiv_types[i],equivalence.data[k].ff_types[3],5);
    }
    break;
  case 4:
    for (i=0; i < 4; i++) {
      k = find_equiv_type(potential_types[i]);
      if (k > -1) strncpy(equiv_types[i],equivalence.data[k].ff_types[4],5);
    }
    break;

  case 5:
    for (i=0; i < 4; i++) {
      k = find_equiv_type(potential_types[i]);
      if (k > -1)
        /* XXX: this leads to an out-of-bounds access.
           the change below does not seem to make a difference.
           strncpy(equiv_types[i],equivalence.data[k].ff_types[5],5); */
        strncpy(equiv_types[i],equivalence.data[k].ff_types[4],5);
    }
    break;
  default:
    break;
  }
  return;
}

int find_equiv_type(char potential_type[5])
{
  int j,k,match;

  j = -1;
  k = 0;
  match = 0;

  while (!match && (k < equivalence.entries)) {
    if (strncmp(potential_type,
                equivalence.data[k].ff_types[0],5) == 0) {
      match = 1;
      j = k;
    } else {
      k++;
    }
  }
  if (j < 0)
    printf(" Unable to find equivalent type for %s\n",potential_type);
  return j;
}
