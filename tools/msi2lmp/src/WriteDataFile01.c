/*
*  This function creates and writes the data file to be used with LAMMPS
*/

#include "Msi2LMP2.h"
#include "Forcefield.h"

void WriteDataFile01(char *nameroot,int forcefield)
{
   int i,j,k,m;
   char line[MAX_LINE_LENGTH];
   FILE *DatF;

/* Open data file */

   sprintf(line,"%s.lammps01",rootname);
   if (pflag > 0) fprintf(stderr," Writing LAMMPS 2001 data file: %s\n",line);
   if( (DatF = fopen(line,"w")) == NULL ) {
      fprintf(stderr,"Cannot open %s\n",line);
      exit(2);
   }



   if (forcefield == 1) total_no_angle_angles = 0;

   fprintf(DatF, "LAMMPS data file for %s\n\n", nameroot);
   fprintf(DatF, " %6d atoms\n", total_no_atoms);
   fprintf(DatF, " %6d bonds\n", total_no_bonds);
   fprintf(DatF, " %6d angles\n",total_no_angles);
   fprintf(DatF, " %6d dihedrals\n", total_no_dihedrals);
   fprintf(DatF, " %6d impropers\n", total_no_oops+total_no_angle_angles);
   fprintf(DatF, "\n");


   fprintf(DatF, " %3d atom types\n", no_atom_types);
   if (no_bond_types > 0)
     fprintf(DatF, " %3d bond types\n", no_bond_types);
   if (no_angle_types> 0)
     fprintf(DatF, " %3d angle types\n", no_angle_types);
   if (no_dihedral_types > 0) fprintf  (DatF," %3d dihedral types\n",
					no_dihedral_types);
   if (forcefield > 1) {
     if ((no_oop_types + no_angleangle_types) > 0)
       fprintf  (DatF, " %3d improper types\n", 
		 no_oop_types + no_angleangle_types);
   }
   else {
     if (no_oop_types > 0)
       fprintf  (DatF, " %3d improper types\n", no_oop_types);
   }

   fprintf(DatF, "\n");
   fprintf(DatF, " %15.9f %15.9f xlo xhi\n", pbc[0], pbc[3]);
   fprintf(DatF, " %15.9f %15.9f ylo yhi\n", pbc[1], pbc[4]);
   fprintf(DatF, " %15.9f %15.9f zlo zhi\n", pbc[2], pbc[5]);

   /* MASSES */


   fprintf(DatF, "\nMasses\n\n");
   for(k=0; k < no_atom_types; k++)
     fprintf(DatF, " %3d %10.6f\n",k+1,atomtypes[k].mass);
   fprintf(DatF, "\n");


   /* COEFFICIENTS */

   fprintf(DatF,"Nonbond Coeffs\n\n");
   for (i=0; i < no_atom_types; i++) {
     fprintf(DatF, " %3i ", i+1);
     for ( j = 0; j < 2; j++)
       fprintf(DatF, "%14.10f ", atomtypes[i].params[j]);
     fprintf(DatF, "\n");
   }
   fprintf(DatF, "\n");

   if (no_bond_types > 0) {
     if (forcefield == 1) m = 2;
     if (forcefield == 2) m = 4;

     fprintf(DatF,"Bond Coeffs\n\n");
     for (i=0; i < no_bond_types; i++) {
       fprintf(DatF, "%3i ", i+1);
       for ( j = 0; j < m; j++)
	 fprintf(DatF, "%10.4f ", bondtypes[i].params[j]);
       fprintf(DatF,"\n");
     }
     fprintf(DatF, "\n");
   }

   if (no_angle_types > 0) {

     if (forcefield == 1) m = 2;
     if (forcefield == 2) m = 4;

     fprintf(DatF,"Angle Coeffs\n\n");
     for (i=0; i < no_angle_types; i++) {
       fprintf(DatF, "%3i ", i+1);
       for ( j = 0; j < m; j++)
	
	 fprintf(DatF, "%10.4f ", angletypes[i].params[j]);
       fprintf(DatF,"\n");
     }
     fprintf(DatF, "\n");
   }

   if (no_dihedral_types > 0) {
     if (forcefield == 1) m = 3;
     if (forcefield == 2) m = 6;

     fprintf(DatF,"Dihedral Coeffs\n\n");
     for (i=0; i < no_dihedral_types; i++) {
       fprintf(DatF, "%3i ", i+1);
       for ( j = 0; j < m; j++)
	 fprintf(DatF, "%10.4f ", dihedraltypes[i].params[j]);
       fprintf(DatF,"\n");
     }
     fprintf(DatF, "\n");
   }
   if (forcefield == 1) {
     if (no_oop_types > 0) {
       fprintf(DatF,"Improper Coeffs\n\n");
       for (i=0; i < no_oop_types; i++) {
	 fprintf(DatF, "%3i ", i+1);
	 for ( j = 0; j < 3; j++)
	   fprintf(DatF, "%10.4f ", ooptypes[i].params[j]);
	 fprintf(DatF, "\n");
       }
       fprintf(DatF, "\n");
     }
   }
   else if (forcefield == 2) {
     if ((no_oop_types + no_angleangle_types) > 0) {
       fprintf(DatF,"Improper Coeffs\n\n");
       for (i=0; i < no_oop_types; i++) {
	 fprintf(DatF, "%3i ", i+1);
	 for ( j = 0; j < 2; j++)
	   fprintf(DatF, "%10.4f ", ooptypes[i].params[j]);
	 fprintf(DatF, "\n");
       }
       for (i=0; i < no_angleangle_types; i++) {
	 fprintf(DatF, "%3i ", i+no_oop_types+1);
	 for ( j = 0; j < 2; j++)
	   fprintf(DatF, "%10.4f ", 0.0);
	 fprintf(DatF, "\n");
       }
       fprintf(DatF, "\n");
     }
   }
   if (forcefield == 2) {

     if (no_angle_types > 0) {
       fprintf(DatF,"BondBond Coeffs\n\n");
       for (i=0; i < no_angle_types; i++) {
	 fprintf(DatF, "%3i ", i+1);
	 for ( j = 0; j < 3; j++)
	   fprintf(DatF, "%10.4f ", angletypes[i].bondbond_cross_term[j]);
	 fprintf(DatF, "\n");
       }
       fprintf(DatF, "\n");

       fprintf(DatF,"BondAngle Coeffs\n\n");

       for (i=0; i < no_angle_types; i++) {
	 fprintf(DatF, "%3i ", i+1);
	 for ( j = 0; j < 4; j++)
	   fprintf(DatF, "%10.4f ",angletypes[i].bondangle_cross_term[j]);
	 fprintf(DatF, "\n");
       }
       fprintf(DatF, "\n");
     }

     if ((no_oop_types+no_angleangle_types) > 0) {
       fprintf(DatF,"AngleAngle Coeffs\n\n");
       for (i=0; i < no_oop_types; i++) {
	 fprintf(DatF, "%3i ", i+1);
	 for ( j = 0; j < 6; j++)
	   fprintf(DatF, "%10.4f ", ooptypes[i].angleangle_params[j]);
	 fprintf(DatF, "\n");
       }
       for (i=0; i < no_angleangle_types; i++) {
	 fprintf(DatF, "%3i ", i+no_oop_types+1);
	 for ( j = 0; j < 6; j++)
	   fprintf(DatF, "%10.4f ", angleangletypes[i].params[j]);
	 fprintf(DatF, "\n");
       }
       fprintf(DatF, "\n");
     }

     if (no_dihedral_types > 0) {
       fprintf(DatF,"AngleAngleTorsion Coeffs\n\n");
       for (i=0; i < no_dihedral_types; i++) {
	 fprintf(DatF, "%3i ", i+1);
	 for ( j = 0; j < 3; j++)
	   fprintf(DatF,"%10.4f ", 
		   dihedraltypes[i].angleangledihedral_cross_term[j]);
	 fprintf(DatF, "\n");
       }
       fprintf(DatF, "\n");

       fprintf(DatF,"EndBondTorsion Coeffs\n\n");
       for (i=0; i < no_dihedral_types; i++) {
	 fprintf(DatF, "%i ", i+1);
	 for ( j = 0; j < 8; j++)
	   fprintf(DatF, "%10.4f ", 
		   dihedraltypes[i].endbonddihedral_cross_term[j]);
	 fprintf(DatF, "\n");
       }
       fprintf(DatF, "\n");

       fprintf(DatF,"MiddleBondTorsion Coeffs\n\n");
       for (i=0; i < no_dihedral_types; i++) {
	 fprintf(DatF, "%3i ", i+1);
	 for ( j = 0; j < 4; j++)
	   fprintf(DatF,"%10.4f ", 
		   dihedraltypes[i].midbonddihedral_cross_term[j]);
	 fprintf(DatF, "\n");
       }
       fprintf(DatF, "\n");


       fprintf(DatF,"BondBond13 Coeffs\n\n");
       for (i=0; i < no_dihedral_types; i++) {
	 fprintf(DatF, "%3i ", i+1);
	 for ( j = 0; j < 3; j++)
	   fprintf(DatF, "%10.4f ",
		   dihedraltypes[i].bond13_cross_term[j]);
	 fprintf(DatF, "\n");
       }
       fprintf(DatF, "\n");

       fprintf(DatF,"AngleTorsion Coeffs\n\n");
       for (i=0; i < no_dihedral_types; i++) {
	 fprintf(DatF, "%3i ", i+1);
	 for ( j = 0; j < 8; j++)
	   fprintf(DatF, "%10.4f ",
		   dihedraltypes[i].angledihedral_cross_term[j]);
	 fprintf(DatF, "\n");
       }
       fprintf(DatF, "\n");
     }
   }

   /*--------------------------------------------------------------------*/

   /* ATOMS */

   fprintf(DatF, "Atoms\n\n");
   for(k=0; k < total_no_atoms; k++) {
      fprintf(DatF, " %6i %6i %3i %9.6f %15.9f %15.9f %15.9f\n",
              k+1,
              atoms[k].molecule,
	      atoms[k].type+1,
              atoms[k].q,
              atoms[k].x[0],
	      atoms[k].x[1],
              atoms[k].x[2]);
   }
   fprintf(DatF, "\n");

   /***** BONDS *****/

   if (total_no_bonds > 0) {
     fprintf(DatF, "Bonds\n\n");
     for(k=0; k < total_no_bonds; k++)
       fprintf(DatF, "%6i %3i %6i %6i\n",k+1,
	       bonds[k].type+1,
	       bonds[k].members[0]+1,
	       bonds[k].members[1]+1);
     fprintf(DatF,"\n");
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
     fprintf(DatF, "\n");
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
      fprintf(DatF, "\n");
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
     if (forcefield == 2) {
       for (k=0; k < total_no_angle_angles; k++)
	 fprintf(DatF, "%6i %3i %6i %6i %6i %6i \n",k+total_no_oops+1,
		 angleangles[k].type+no_oop_types+1,
		 angleangles[k].members[0]+1,
		 angleangles[k].members[1]+1,
		 angleangles[k].members[2]+1,
		 angleangles[k].members[3]+1);
     }
     fprintf(DatF, "\n");
   }

/* Close data file */

   if (fclose(DatF) !=0) {
      fprintf(stderr,"Error closing %s.lammps01\n", rootname);
      exit(1);
   }

}
