/******************************
*
*  This function opens the .mdf file and extracts connectivity information
*  into the atoms Atom structure.  It also updates the charge from the .car
*  file because the charge in the .mdf file has more significant figures.
*
*/

#include "msi2lmp.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>


/* Prototype for function to process a single atom
   Returns int that flags end of data file */

static int get_molecule(char line[], int connect_col_no,
                        int q_col_no, int *counter);

/* Prototype for function that takes connectivty record as stated in
   .mdf file and fills in any default values */
static void MakeConnectFullForm(int *counter);

/* prototype for function to clean strange characters out of strings */

static void clean_string(char *);

static int blank_line(char *line)
{
  while (*line != '\0') {
    if (isalnum((int) *line)) return 0;
    ++line;
  }
  return 1;
}

void ReadMdfFile(void)
{
  char line[MAX_LINE_LENGTH];        /* Temporary storage for reading lines */
  char *col_no;                /* Pointer to column number stored as char */
  char *col_name;                /* Pointer to column name */
  int connect_col_no = 0;        /* Column number where connection info begins */
  int q_col_no = 0;                /* Column number containg charge information */
  int atom_counter=0;          /* Keeps track of current atom number */

  int i,j,k,kk,l,n,match,match2,status;
  char *temp_string;
  char *temp_residue;
  char *temp_atom_name;
  char *sptr;
  unsigned char at_end = 0;

  /* Open .mdf file for reading */


  sprintf(line,"%s.mdf",rootname);
  if (pflag > 0) printf(" Reading mdf file: %s\n",line);
  if ((MdfF = fopen(line,"r")) == NULL ) {
    printf("Cannot open %s\n",line);
    exit(41);
  }

  while (!at_end) {

    sptr = fgets(line,MAX_LINE_LENGTH,MdfF);

    if (sptr != NULL) {

      clean_string(line);

      if (strncmp(line,"#end",4) == 0) {
        at_end = 1;
      } else if (strncmp(line,"@column",7) == 0) {

        temp_string = strtok(line,WHITESPACE);
        col_no = strtok(NULL,WHITESPACE);
        col_name = strtok(NULL,WHITESPACE);
        if (strncmp(col_name,"charge",6) == 0) {
          if (strlen(col_name) < 8) {
            q_col_no = atoi(col_no);
          }
        } else if (strncmp(col_name,"connect",7) == 0) {
          connect_col_no = atoi(col_no);
        }
      } else if (strncmp(line,"@molecule",9) == 0) {

        if ((q_col_no == 0) | (connect_col_no == 0)) {
          printf("Unable to process molecule without knowing charge\n");
          printf("and connections columns\n");
          exit(42);
        }
        sptr = fgets(line,MAX_LINE_LENGTH,MdfF);
        status = get_molecule(line,connect_col_no,q_col_no,&atom_counter);
        if (status == 0) {
          printf("Trouble reading molecule - exiting\n");
          exit(43);
        }
      }
    } else {
      printf("End of File found or error reading line\n");
      at_end = 1;
    }
  }

  /* Next build list of residues for each molecule This will
     facilitate assigning connections numbers as well as figuring
     out bonds, angles, etc.  This first loop just figures out the
     number of residues in each molecule and allocates memory to
     store information for each residue. The second loop fills
     in starting and ending atom positions for each residue
  */

  temp_string = (char *)calloc(MAX_STRING,sizeof(char));

  for (n=0; n < no_molecules; n++) {
    molecule[n].no_residues = 1;

    strncpy(temp_string,atoms[molecule[n].start].residue_string,MAX_NAME);
    for (i=molecule[n].start+1; i < molecule[n].end; i++) {
      if (strncmp(temp_string,atoms[i].residue_string,MAX_NAME) != 0) {
        molecule[n].no_residues++;
        strncpy(temp_string,atoms[i].residue_string,MAX_NAME);
      }
    }

    molecule[n].residue = (struct ResidueList *)
      calloc(molecule[n].no_residues, sizeof(struct ResidueList));

    if (molecule[n].residue == NULL) {
      printf("Unable to allocate memory for residue list - molecule %d\n",n);
      exit(44);
    }
  }
  for (n=0; n < no_molecules; n++) {
    j = 0;
    strncpy(molecule[n].residue[j].name,
            atoms[molecule[n].start].residue_string,MAX_NAME);

    molecule[n].residue[j].start = molecule[n].start;
    for (i=molecule[n].start+1; i < molecule[n].end; i++) {
      if (strncmp(molecule[n].residue[j].name,
                  atoms[i].residue_string,MAX_NAME) != 0) {

        molecule[n].residue[j].end = i;
        molecule[n].residue[++j].start = i;
        strncpy(molecule[n].residue[j].name,atoms[i].residue_string,MAX_NAME);
      }
    }
    molecule[n].residue[j].end = molecule[n].end;
    /*
      printf("Molecule %d has %d residues",n,molecule[n].no_residues);
      for (i=0; i < molecule[n].no_residues; i++) {
      printf(" %s",molecule[n].residue[i].name);
      }
      printf("\n");
      for (i=molecule[n].start; i < molecule[n].end; i++) {
      printf(" atom %d residue %s\n",i,atoms[i].residue_string);
      }
      printf(" residue %s start %d end %d\n",molecule[n].residue[i].name,
      molecule[n].residue[i].start,molecule[n].residue[i].end);
      }
    */
  }

  /* Assign atom names in connections[] to corresponding atom numbers */

  for (n=0; n < no_molecules; n++) {
    for (j=0; j < molecule[n].no_residues; j++) {
      for (i=molecule[n].residue[j].start; i < molecule[n].residue[j].end;
           i++) {
        for (l=0; l < atoms[i].no_connect; l++) {
          strncpy(temp_string,atoms[i].connections[l],MAX_STRING);
          temp_residue = strtok(temp_string,":");
          temp_atom_name = strtok(NULL,"%");

          if (strcmp(temp_residue,molecule[n].residue[j].name) == 0) {

            /* atom and connection are part of same residue
               Search names on just that residue            */

            k = molecule[n].residue[j].start;
            match = 0;
            while (!match && (k < molecule[n].residue[j].end)) {
              if (strcmp(atoms[k].name,temp_atom_name) == 0) {
                atoms[i].conn_no[l] = k;
                match = 1;
              } else
                k++;
            }
            if (match == 0) {
              printf("Unable to resolve atom number of atom %d conn %d string %s:%s\n"
                     " Something is wrong in the MDF file\n",
                     i,l,temp_residue,temp_atom_name);
              exit(45);
            }
          } else {

            /* atom and connection are on different residues
               First find the residue that the connection is
               on then loop over its atoms
            */

            k=0;
            match = 0;
            while (!match && (k < molecule[n].no_residues)) {
              if (strcmp(temp_residue,molecule[n].residue[k].name) == 0) {
                kk = molecule[n].residue[k].start;
                match2 = 0;
                while (!match2 && (kk < molecule[n].residue[k].end)) {
                  if (strcmp(atoms[kk].name,temp_atom_name) == 0) {
                    atoms[i].conn_no[l] = kk;
                    match2 = 1;
                  } else
                    kk++;
                }
                if (match2 == 0) {
                  printf("Unable to resolve atom number of atom %d conn %d string %s\n"
                         " Something is wrong in the MDF file\n",
                         i,l,atoms[i].connections[l]);
                  exit(46);
                }
                match = 1;
              } else
                k++;
            }
            if (match == 0) {
              printf("Unable to find residue associated with conn %d %s on atom %d\n"
                     " Something is wrong in the MDF file\n", l,atoms[i].connections[l],i);
              exit(47);
            }
          } /* end if */
        } /* l - loop over connections on atom i */
      } /* i - loop on atoms in residue j molecule n */
    } /* j - loop on residues in molecule n */
  } /* n - loop over molecules */

  free(temp_string);
  /*
    for (n=0; n < no_molecules; n++) {

    printf("Molecule %d has %d residues\n",n,molecule[n].no_residues);
    for (j=0; j < molecule[n].no_residues; j++) {
    printf(" Residue %d named %s\n",j,molecule[n].residue[j].name);
    for (i=molecule[n].residue[j].start; i < molecule[n].residue[j].end;
    i++) {
    printf("  Atom %d type %s connected to ",i,atoms[i].potential);
    for (l=0; l < atoms[i].no_connect; l++) printf(" %d ",
    atoms[i].conn_no[l]);
    printf("\n");
    }
    }
    }

  */

  /* Close .mdf file */

  if (fclose(MdfF) !=0) {
    printf("Error closing %s.car\n", rootname);
    exit(1);
  }

} /* End ReadMdfFile function */

/*--------------------- get_molecule Function-----------------------*/

int get_molecule(char *line, int connect_col_no, int q_col_no,
                 int *counter)
{
  char *cur_field;        /* For storing current string token */
  int i;                 /* Used in loop counters */
  int connect_no;      /* Connection number within atom */
  int r_val = 1;        /* Return value.  1 = successful
                           0 = EOF encountered */
  /* Loop over atoms */

  /* blank line signals end of molecule*/
  while(!blank_line(fgets(line,MAX_LINE_LENGTH,MdfF))) {
    /*  while(strlen(fgets(line,MAX_LINE_LENGTH,MdfF)) > 2) { */

    clean_string(line);

    /* Get atom name */
    cur_field = strtok(line,":");
    sscanf(cur_field, "%s", atoms[*counter].residue_string);
    cur_field = strtok(NULL,WHITESPACE);
    /* Compare atom name with that in .car file */
    if (strcmp(atoms[*counter].name, cur_field)) {
      printf("Names %s from .car file and %s from .mdf file do not match\n",
             atoms[*counter].name, cur_field);
      printf("counter = %d\n",*counter);
      printf("Program Terminating\n");
      exit(4);
    }

    /* Skip unwanted fields until charge column, then update charge */

    for (i=1; i < q_col_no; i++) strtok(NULL,WHITESPACE);
    cur_field = strtok(NULL, WHITESPACE);
    atoms[*counter].q = atof(cur_field);

    /* Continue skipping unwanted fields until connectivity records begin */

    for ( i = (q_col_no + 1); i < connect_col_no; i++) strtok(NULL,WHITESPACE);

    /* Process connections */

    connect_no = 0; /* reset connections counter */
    while ((cur_field = strtok(NULL,WHITESPACE)) && (connect_no < MAX_CONNECTIONS)) {
      sscanf(cur_field, "%s", atoms[*counter].connections[connect_no++]);
    }
    atoms[*counter].no_connect = connect_no;
    MakeConnectFullForm(counter);
    (*counter)++;

  } /* End atom processing loop */

  return r_val;

} /* End get_molecule function */



 /*------------------------MakeConnectFullForm Function--------------------*/

void MakeConnectFullForm(int *counter) {

  /* This function processes the connection names after all connections
     for an atom have been read in.
     It replaces any short forms that use implied default values
     with the full form connectivity record */

  int i;                /* Counter for character array */
  int j;                /* loop counter */
  char tempname[2*MAX_STRING];   /* name of connection */
  char tempcell[10];   /* Values from connectivity record */
  char tempsym[5];              /* " " */
  char tempbo[6];               /* " " */
  char *charptr;

  for ( j = 0; j < atoms[*counter].no_connect; j++) {
    /* If not full name, make name full */
    if (strchr(atoms[*counter].connections[j],':') == NULL) {
      strcpy(tempname,atoms[*counter].residue_string);
      strcat(tempname,":");
      strcat(tempname,
             atoms[*counter].connections[j]);
      /* truncate at capacity of target storage */
      tempname[MAX_STRING-1] = '\0';
      sscanf(tempname, "%s",
             atoms[*counter].connections[j]);
    } else sscanf(atoms[*counter].connections[j], "%s", tempname);
    /* Set cell variables */

    i=0;
    charptr = (strchr(tempname,'%'));
    if (charptr != NULL) {
      while ( *charptr!='#' && *charptr!='/' && *charptr!='\000')
        tempcell[i++] = *(charptr++);
      tempcell[i] = '\000';
    } else strcpy(tempcell, "%000");

    /* Set symmetry variables
       -- If not 1, cannot handle at this time */

    i = 0;
    charptr = (strchr(tempname,'#'));
    if (charptr != NULL)  {
      while (*charptr != '/' && *charptr !='\000') {
        tempsym[i++] = *(charptr++);
        if ((i==2) && (tempsym[1] != '1')) {
          printf("Msi2LMP is not equipped to handle symmetry operations\n");
          exit(5);
        }
      }
      tempsym[i] = '\000';
    } else strcpy(tempsym, "#1");

    /* Set bond order and record in data structure */

    i = 0;
    charptr = strchr(tempname,'/');
    if (charptr != NULL) {
      charptr++;
      while (*charptr != '\000')
        tempbo[i++] = *(charptr++);
      tempbo[i] = '\000';
    } else strcpy(tempbo, "1.0");

    atoms[*counter].bond_order[j] = atof(tempbo);

    /* Build connection name and store in atoms data structure */

    strtok( tempname, "%#/");
    strcat( tempname, tempcell);
    strcat( tempname, tempsym);
    strcat( tempname, "/");
    strcat( tempname, tempbo);
    if (strlen(tempname) >= MAX_STRING) {
      printf("tempname overrun %s. truncating to ",tempname);
      /* truncate at capacity of target storage */
      tempname[MAX_STRING-1] = '\0';
      puts(tempname);
    }
    sscanf(tempname, "%s", atoms[*counter].connections[j]);
  }/*End for loop*/
}/* End function MakeNameLong
  */

void clean_string(char *string) {
  int i,n;
  short k;

  n = strlen(string);
  for (i=0; i < n; i++) {
    k = (short)string[i];
    if ((k<32) | (k>127)) string[i] = '\0';
  }
}

