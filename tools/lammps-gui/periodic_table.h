// clang-format off
/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: periodic_table.h,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.12 $       $Date: 2009/01/21 17:45:41 $
 *
 ***************************************************************************/

/*
 * periodic table of elements and helper functions to convert
 * ordinal numbers to labels and back.
 * all tables and functions are declared static, so that it
 * can be safely included by all plugins that may need it.
 *
 * 2002-2009 akohlmey@cmm.chem.upenn.edu, vmd@ks.uiuc.edu
 */

#include <string.h>
#include <ctype.h>

/* periodic table of elements for translation of ordinal to atom type */
static const char *pte_label[] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P" , "S",  "Cl", "Ar", "K",  "Ca", "Sc",
    "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc",
    "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
    "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg"
};
static const int nr_pte_entries = sizeof(pte_label) / sizeof(char *);

/* corresponding table of masses. */
static const double pte_mass[] = {
    /* X  */ 0.00000, 1.00794, 4.00260, 6.941, 9.012182, 10.811,
    /* C  */ 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
    /* Na */ 22.989770, 24.3050, 26.981538, 28.0855, 30.973761,
    /* S  */ 32.065, 35.453, 39.948, 39.0983, 40.078, 44.955910,
    /* Ti */ 47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.9332,
    /* Ni */ 58.6934, 63.546, 65.409, 69.723, 72.64, 74.92160,
    /* Se */ 78.96, 79.904, 83.798, 85.4678, 87.62, 88.90585,
    /* Zr */ 91.224, 92.90638, 95.94, 98.0, 101.07, 102.90550,
    /* Pd */ 106.42, 107.8682, 112.411, 114.818, 118.710, 121.760,
    /* Te */ 127.60, 126.90447, 131.293, 132.90545, 137.327,
    /* La */ 138.9055, 140.116, 140.90765, 144.24, 145.0, 150.36,
    /* Eu */ 151.964, 157.25, 158.92534, 162.500, 164.93032,
    /* Er */ 167.259, 168.93421, 173.04, 174.967, 178.49, 180.9479,
    /* W  */ 183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,
    /* Hg */ 200.59, 204.3833, 207.2, 208.98038, 209.0, 210.0, 222.0,
    /* Fr */ 223.0, 226.0, 227.0, 232.0381, 231.03588, 238.02891,
    /* Np */ 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0,
    /* Md */ 258.0, 259.0, 262.0, 261.0, 262.0, 266.0, 264.0, 269.0,
    /* Mt */ 268.0, 271.0, 272.0
};

/*
 * corresponding table of VDW radii.
 * van der Waals radii are taken from A. Bondi,
 * J. Phys. Chem., 68, 441 - 452, 1964,
 * except the value for H, which is taken from R.S. Rowland & R. Taylor,
 * J.Phys.Chem., 100, 7384 - 7391, 1996. Radii that are not available in
 * either of these publications have RvdW = 2.00 \AA
 * The radii for Ions (Na, K, Cl, Ca, Mg, and Cs are based on the CHARMM27
 * Rmin/2 parameters for (SOD, POT, CLA, CAL, MG, CES) by default.
 */
static const double pte_vdw_radius[] = {
    /* X  */ 1.5, 1.2, 1.4, 1.82, 2.0, 2.0,
    /* C  */ 1.7, 1.55, 1.52, 1.47, 1.54,
    /* Na */ 1.36, 1.18, 2.0, 2.1, 1.8,
    /* S  */ 1.8, 2.27, 1.88, 1.76, 1.37, 2.0,
    /* Ti */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Ni */ 1.63, 1.4, 1.39, 1.07, 2.0, 1.85,
    /* Se */ 1.9, 1.85, 2.02, 2.0, 2.0, 2.0,
    /* Zr */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Pd */ 1.63, 1.72, 1.58, 1.93, 2.17, 2.0,
    /* Te */ 2.06, 1.98, 2.16, 2.1, 2.0,
    /* La */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Eu */ 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Er */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* W  */ 2.0, 2.0, 2.0, 2.0, 1.72, 1.66,
    /* Hg */ 1.55, 1.96, 2.02, 2.0, 2.0, 2.0, 2.0,
    /* Fr */ 2.0, 2.0, 2.0, 2.0, 2.0, 1.86,
    /* Np */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Md */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Mt */ 2.0, 2.0, 2.0
};

/* lookup functions */

static const char *get_pte_label(const int idx)
{
    if ((idx < 1) || (idx >= nr_pte_entries)) return pte_label[0];

    return pte_label[idx];
}

static double get_pte_mass(const int idx)
{
    if ((idx < 1) || (idx >= nr_pte_entries)) return pte_mass[0];

    return pte_mass[idx];
}

static double get_pte_vdw_radius(const int idx)
{
    if ((idx < 1) || (idx >= nr_pte_entries)) return pte_vdw_radius[0];

#if 1
    /* Replace with Hydrogen radius with an "all-atom" radius */
    if (idx == 1)
      return 1.0;    /* H  */
#else
    /* Replace with old VMD atom radii values */
    switch (idx) {
      case  1: return 1.0;    /* H  */
      case  6: return 1.5;    /* C  */
      case  7: return 1.4;    /* N  */
      case  8: return 1.3;    /* O  */
      case  9: return 1.2;    /* F  */
      case 15: return 1.5;    /* P  */
      case 16: return 1.9;    /* S  */
    }
#endif

    return pte_vdw_radius[idx];
}

static int get_pte_idx(const char *label)
{
    int i;
    char atom[3];

    /* zap string */
    atom[0] = (char) 0;
    atom[1] = (char) 0;
    atom[2] = (char) 0;
    /* if we don't have a null-pointer, there must be at least two
     * chars, which is all we need. we convert to the capitalization
     * convention of the table above during assignment. */
    if (label != NULL) {
        atom[0] = (char) toupper((int) label[0]);
        atom[1] = (char) tolower((int) label[1]);
    }
    /* discard numbers in atom label */
    if (isdigit(atom[1])) atom[1] = (char) 0;

    for (i=0; i < nr_pte_entries; ++i) {
        if ( (pte_label[i][0] == atom[0])
             && (pte_label[i][1] == atom[1]) ) return i;
    }

    return 0;
}

static int get_pte_idx_from_string(const char *label) {
  int i, ind;
  char atom[3];

  if (label != NULL) {
    /* zap string */
    atom[0] = atom[1] = atom[2] = '\0';

    for (ind=0,i=0; (ind<2) && (label[i]!='\0'); i++) {
      if (label[i] != ' ') {
        atom[ind] = toupper(label[i]);
        ind++;
      }
    }

    if (ind < 1)
      return 0; /* no non-whitespace characters */

    for (i=0; i < nr_pte_entries; ++i) {
      if ((toupper(pte_label[i][0]) == atom[0]) && (toupper(pte_label[i][1]) == atom[1]))
        return i;
    }
  }

  return 0;
}

#if 0
#include <stdio.h>

int main() {
  int i;

  printf("Periodic table check/dump\n");
  printf("  Table contains data for %d elements\n", nr_pte_entries);
  printf("   Mass table size check: %d\n", sizeof(pte_mass) / sizeof(double));
  printf("    VDW table size check: %d\n", sizeof(pte_vdw_radius) / sizeof(double));
  printf("\n");
  printf("Symbol Num    Mass   rVDW\n");
  for (i=0; i<nr_pte_entries; i++) {
    printf("   %-2s  %3d  %6.2f  %4.2f\n",
      get_pte_label(i), i, get_pte_mass(i), get_pte_vdw_radius(i));
  }
  return 0;
}
#endif
