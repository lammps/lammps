/*
 *   This function fills in the keyword field, the number of members for each
 *   item and the number of parameters for each item
 *
 */

#include "msi2lmp.h"
#include "Forcefield.h"

#include <string.h>

void InitializeItems(void)
{
  /* ATOM TYPES */
  strcpy(ff_atomtypes.keyword,"#atom_types");
  ff_atomtypes.number_of_members = 1;
  ff_atomtypes.number_of_parameters = 1;

  /* EQUIVALENCE */

  strcpy(equivalence.keyword,"#equivalence");
  equivalence.number_of_members = 6;
  equivalence.number_of_parameters = 0;

  /* NON-BOND */

  strcpy(ff_vdw.keyword,"#nonbond");
  ff_vdw.number_of_members = 1;
  ff_vdw.number_of_parameters = 2;

  /* BOND */

  ff_bond.number_of_members = 2;
  if (forcefield & (FF_TYPE_CLASS1|FF_TYPE_OPLSAA)) {
    strcpy(ff_bond.keyword,"#quadratic_bond");
    ff_bond.number_of_parameters = 2;
  }

  if (forcefield & FF_TYPE_CLASS2) {
    strcpy(ff_bond.keyword,"#quartic_bond");
    ff_bond.number_of_parameters = 4;
  }

  /* MORSE */

  if (forcefield & FF_TYPE_CLASS1) {
    ff_morse.number_of_members = 2;
    strcpy(ff_morse.keyword,"#morse_bond");
    ff_morse.number_of_parameters = 3;
  }

  /* ANGLE */

  ff_ang.number_of_members = 3;
  if (forcefield & (FF_TYPE_CLASS1|FF_TYPE_OPLSAA)) {
    strcpy(ff_ang.keyword,"#quadratic_angle");
    ff_ang.number_of_parameters = 2;
  }

  if (forcefield & FF_TYPE_CLASS2) {
    strcpy(ff_ang.keyword,"#quartic_angle");
    ff_ang.number_of_parameters = 4;
  }

  /* TORSION */

  ff_tor.number_of_members = 4;
  if (forcefield & (FF_TYPE_CLASS1|FF_TYPE_OPLSAA)) {
    strcpy(ff_tor.keyword,"#torsion_1");
    ff_tor.number_of_parameters = 3;
  }

  if (forcefield & FF_TYPE_CLASS2) {
    strcpy(ff_tor.keyword,"#torsion_3");
    ff_tor.number_of_parameters = 6;
  }

  /* OOP */

  ff_oop.number_of_members = 4;
  if (forcefield & (FF_TYPE_CLASS1|FF_TYPE_OPLSAA)) {
    strcpy(ff_oop.keyword,"#out_of_plane");
    ff_oop.number_of_parameters = 3;
  }

  if (forcefield & FF_TYPE_CLASS2) {
    strcpy(ff_oop.keyword,"#wilson_out_of_plane");
    ff_oop.number_of_parameters = 2;
  }

  if (forcefield & FF_TYPE_CLASS2) {
    /* BOND-BOND */

    strcpy(ff_bonbon.keyword,"#bond-bond");
    ff_bonbon.number_of_members = 3;
    ff_bonbon.number_of_parameters = 1;

    /* BOND-ANGLE */

    strcpy(ff_bonang.keyword,"#bond-angle");
    ff_bonang.number_of_members = 3;
    ff_bonang.number_of_parameters = 2;

    /* ANGLE-TORSION */

    strcpy(ff_angtor.keyword,"#angle-torsion_3");
    ff_angtor.number_of_members = 4;
    ff_angtor.number_of_parameters = 6;

    /* ANGLE-ANGLE-TORSION */

    strcpy(ff_angangtor.keyword,"#angle-angle-torsion_1");
    ff_angangtor.number_of_members = 4;
    ff_angangtor.number_of_parameters = 1;

    /* END-BOND-TORSION */

    strcpy(ff_endbontor.keyword,"#end_bond-torsion_3");
    ff_endbontor.number_of_members = 4;
    ff_endbontor.number_of_parameters = 6;

    /* MID-BOND-TORSION */

    strcpy(ff_midbontor.keyword,"#middle_bond-torsion_3");
    ff_midbontor.number_of_members = 4;
    ff_midbontor.number_of_parameters = 3;

    /* ANGLE-ANGLE */

    strcpy(ff_angang.keyword,"#angle-angle");
    ff_angang.number_of_members = 4;
    ff_angang.number_of_parameters = 1;

    /* BOND-BOND-1-3 */

    strcpy(ff_bonbon13.keyword,"#bond-bond_1_3");
    ff_bonbon13.number_of_members = 4;
    ff_bonbon13.number_of_parameters = 1;
  }
}
