/****************************** Lmp2arc.c ******************************

   This program converts a contiguous set of LAMMPS dump (position) files 
   to a MSI .arc file using a .car file as the template. This enables one to
   visualize a LAMMPS trajectory using Insight II.

   by Michael Peachey & John Carpenter  (SGI/Cray)
   currently maintained by Kelly L. Anderson (AFRL/ML)
   
   Version 0.1:  8/14/97  MJP
   Version 1.0:  7/27/98  JEC
   Version 1.1: 10/22/98  JEC
   Version 1.2: 11/04/98  JEC
   Version 1.3:  4/10/05  KLA


   1. To build the program

     % cd src
     % make
     % cp lmp2arc.exe ../bin

     The program is only a few files so should compile
     quickly. 

   2. Testing the program

     There are three test jobs in the test directory
     that one can run

     % ../bin/lmp2arc.exe -car water8.car <<EOF > water8.arc
     > water8.pos
     > EOF
     % diff water8.arc correct/water8.arc

     % ../bin/lmp2arc.exe -move_mol -car water8.car <<EOF > water8.arc
     > water8.pos
     > EOF
     % diff water8.arc correct/water8.arc.unwrap

     % ../bin/lmp2arc.exe -trueflags -move_mol -car water8.car \
       <<EOF > water8.arc
     > water8.pos
     > EOF
     % diff water8.arc correct/water8.arc.trueflags

     The only differences you should see are date differences.

  3. Usage
    
     Usage: lmp2arc.exe [options] -car carfile < infile > outfile

      options

	      -trueflags  (optional)

	         trueflags are present in the LAMMPS position file.

	      -move_mol   (optional)

                 Unwraps molecules. All atoms in the position file are in
                 the simulation cell. This means that molecules which
                 are sticking out of the cell will have some of its
                 atoms "wrapped" to the other side of the periodic cell.
                 This leads to very messy visualizations. Specification
                 of -move_mol will attempt to "unwrap" or straighten 
                 molecules and yield a reasonable visualization. However,
                 molecules which oscillate around half in or out of the
                 box may appear to jump from one side of the box to
                 the other. 

                 Two algorithms are used depending upon the presence of
                 trueflags. The algorithm that uses trueflags is the
                 most robust, but the other geometric based algorithm 
                 should be adequate.

	      -car filename

	         the name of the .car file that corresponds to the
		 position information. This is required.

	      -skip n   (default is 0)

	         Skip every n timesteps in the position file.

	      -npico n  (default is 2000)

	         Number of timesteps in 1 picosecond of simulation.

              -2001   (default is 2005)
 
                 Version of LAMMPS to convert from.
                 

	      stdin

	         file containing one or more names of LAMMPS position files

	      stdout

	         the name of the new .arc file

   Examples:

    % lmp2arc.exe -2001 -trueflags -move_mol -skip 4 -car water.car <<EOF > water.arc
      water1.pos
      water2.pos
      EOF

>>>> Program output <<<<<

 lmp2arc v1.2 - LAMMPS MD trajectory to MSI .arc file conversion


 Car file name is water.car

 Number of Atoms       = 24
 Number of Molecules   = 8

 Position file names:
 water1.pos
 water2.pos

 Processing Timesteps

 20 40 60 80 100

 102 frames were written to the ArcFile


 Program Exiting Normally


 4. Questions or comments?

    Send them to Kelly L. Anderson (kelly.anderson@cantab.net)


*/
#define MAIN
#include "lmp2.h" 

 /* External function prototypes */
   
   extern void ReadCarFile(FILE *,struct Sys *);
   extern void ProcessPosFile01(int,char **,struct Sys *,FILE *);
   extern void ProcessPosFile05(int,char **,struct Sys *,FILE *);

int main(int argc, char *argv[])
{
   int got_car,ierr,i;
   char *carfilename;
   struct Sys *sysinfo;
   int num_posfiles;
   int outv;
   char *posnames[MAX_POS_FILES];
   char line[MAX_LINE_LENGTH];

   FILE *CarFile,*ArcFile;

/* default input values */

  trueflag = 0;
  move_molecules = 0;
  nskip = 0;
  got_car = 0;
  npico = 2000;
  outv = 2005;
  
/* read input options from command line
   should do more error checking for missing args */

  for (i = 1; i < argc; i++) {

    if (!strcmp(argv[i],"-car")) {
      carfilename = (char *) malloc(strlen(argv[i+1]));
      strcpy(carfilename,argv[i+1]);
      i++;
      got_car = 1; 
      }
    else if (!strcmp(argv[i],"-trueflags"))
      trueflag = 1;
    else if (!strcmp(argv[i],"-move_mol"))
      move_molecules = 1;
    else if (!strcmp(argv[i],"-2001")) {
      outv = 2001;
      fprintf(stderr,"\n LAMMPS 2001 file used.\n"); 
     }
    else if (!strcmp(argv[i],"-skip")) {
      sscanf(argv[i+1],"%d",&nskip);
      i++;
      }
    else if (!strcmp(argv[i],"-npico")) {
      sscanf(argv[i+1],"%d",&npico);
      i++;
      }
      else {
      fprintf(stderr,"\nSyntax error: lmp2arc.exe [-2001] [-trueflags] [-move_mol] [-skip n]\n");
      fprintf(stderr,"                      [-npico m] -car carfile\n");
      fprintf(stderr,"                      < posfile_list > ArcFile\n");
      exit(1);
      }
  }

/* Check for option inconsistancies */

  if (!got_car) {
     fprintf(stderr,"Must specify car file\n");
     exit(1);
  }

   fprintf(stderr,"\n lmp2arc v1.3 - LAMMPS MD trajectory to ACCELRYS .arc file conversion\n\n");

   sysinfo = (struct Sys *) calloc(1,sizeof(struct Sys));

   /* Open CarFile */

   fprintf(stderr,"\n Car file name is %s\n",carfilename);
   if( (CarFile = fopen(carfilename,"r")) == NULL )
   {
      fprintf(stderr,"Cannot open %s\n",carfilename);
      exit(2);
   }

   ReadCarFile(CarFile, sysinfo);

   fclose(CarFile);
   free(carfilename);

   fprintf(stderr,"\n Number of Atoms       = %d\n",sysinfo->natoms);
   fprintf(stderr," Number of Molecules   = %d\n",sysinfo->no_molecules);

   /* Get position file names */

   num_posfiles = 0;
   while (gets(line) != NULL) {
     if (ferror(stdin) == 0) {
       posnames[num_posfiles] = (char *)calloc(strlen(line),sizeof(char));
       posnames[num_posfiles] = strcpy(posnames[num_posfiles],line);
       num_posfiles++;
     }
     else
       perror(strerror(ferror(stdin)));
   }

   fprintf(stderr,"\n Position file names:\n");
   for (i=0;i < num_posfiles; i++) fprintf(stderr," %s\n",posnames[i]);

   /* Assign ArcFile to stdout */

   ArcFile = stdout;

   if (outv == 2001) {  
       ProcessPosFile01(num_posfiles,posnames,sysinfo,ArcFile); 
   }
   else if (outv == 2005) {  
       ProcessPosFile05(num_posfiles,posnames,sysinfo,ArcFile); 
   }

   fclose(ArcFile);

   free(sysinfo);

   fprintf(stderr,"\n\n Program Exiting Normally\n\n");
}
