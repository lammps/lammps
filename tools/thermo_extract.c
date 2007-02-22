/* Written by Vikas Varshney
  Dated Feb 22, 2007
  Contact Info: vikas.varshney@gmail.com

The main purpose of this simple code is to extract thermodynamic
information such as Temperature, Pressure, Energy, etc from LAMMPS log
files, to plot them as a function of timestep in a simple plotting
programs such as xmgrace.

Such plots are useful at initial stages of a simulation to make sure
the system is equilibrated.

To run this program please compile it as
 gcc -O3 thermo_extract.c -lm -o thermo_extract

To extract temperature as a function of timestep, run the code as

thermo_extract -p Temp -s/-m logfile1 logfile2 ... (all the log files)

Time and temperature is output as 2 columns of data.  It can
be re-directed to a file for plotting.

 -p for thermo parameter to extract
 -s if the output format is single line, or
 -m is the output format is multi line

Few points

1. The log file must contain a "thermo N" command where N is the
  frequency of the thermo output.

2. The name of the property "Temp" (in previous example) has to match
  the thermo property in the logfile.

3. The log file(s) can contain thermo info for multiple runs.

4. Currently, you can only use one parameter as input, so for extracting
  multiple properties, you have to run the code again.

5. The Program uses the following keywords to keep track of the data in
  the file. In general, these keywords are not used in input script (which
  in turn get printed in log file)  but are worth mentioning. 

  PLEASE avoid using  
  "Dangerous", "Memory", "Step" and "Loop" in the input script

  "thermo" and "run" except when they are used to specify frequency of 
   thermo output and total number of steps to be run, respectively.

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>
#include <stddef.h>

FILE *input;
char filename_input[100];
char string1[100],*string2;
char temp[80];
int j,k,check,timestep;
double datapoint;
int thermostylecounter,thermocounter,runcounter,totalsteps;
int narg1;
char *token;
char templl[80];
int thermo[100];
int run[100];

void scan_first_part()
{
 while(strcmp(temp,"Memory"))
  fscanf(input,"%s ",&temp);
}

void scan_data_multi(int temp1,char* chartemp)
{
 totalsteps=run[temp1]/thermo[temp1]+1;
 for(j=0;j<totalsteps;j++)
  {
    while(strcmp(temp,"Step"))
      fscanf(input,"%s ",&temp);

    fscanf(input,"%s ",&temp);
    timestep=atoi(temp);

    while(strcmp(temp,chartemp))
      fscanf(input,"%s ",&temp);

    fscanf(input,"%s ",&temp);
    fscanf(input,"%s ",&temp);
    datapoint=atof(temp);
    if ((temp1!=0) && (j==0))
      {}
    else
      printf("%d %f\n",timestep,datapoint);
  }
}

void scan_data_single(int temp1,char* chartemp)
{
 fgets(string1,100,input);
 fgets(string1,100,input);
 narg1=0;
 for(j=0;j<strlen(string1);j++)
  if(string1[j]==' ')
    narg1++;

 token=strtok(string1," ");
 check=0;
 for(j=1;j<narg1;j++)
  {
    token=strtok(NULL," ");
    if (!strcmp(token,chartemp))
      {
        check=j;
        break;
      }
  }

 if (check>0)
  {}
 else
  {
    printf("Couldn't recognize parameter\n");
    exit(1);
  }

 totalsteps=run[temp1]/thermo[temp1]+1;
 for(j=0;j<totalsteps;j++)
  {
    fscanf(input,"%s ",&temp);
    timestep=atoi(temp);
    for(k=0;k<check;k++)
      fscanf(input,"%s ",&temp);
    datapoint=atof(temp);
    for(k=0;k<narg1-1-check;k++)
      fscanf(input,"%s ",&temp);
    if ((temp1!=0) && (j==0))
      {}
    else
      printf("%d %f\n",timestep,datapoint);
  }
}

int main(int arg, char **argv)
{
 int singleoutput,next,eofpointer,i,file;
 char firsttemp[256];
 int iterations;
 if (arg<5)
  {
    printf("Insufficient number of arguements\n");
    printf("Please make sure the systax has the following format\n");
    printf("thermo_extract -p parameter -m/-s logfile1 logfile2 ...\n");
    exit(1);
  }
 for(file=4;file<arg;file++)
  {
    sprintf(filename_input,argv[file]);
    if ((input=fopen(filename_input,"r"))==NULL)
      {
        printf("Error reading %dth file.\n",file);
        exit(1);
      }
    iterations=0;
    while(!feof(input))
      {
        fscanf(input,"%s ",&firsttemp);
        if (!strcmp(firsttemp,"Dangerous"))
          iterations++;
      }
    rewind(input);
    for(i=0;i<100;i++)
      {
        thermo[i]=0;
        run[i]=0;
      }

    for(i=0;i<iterations;i++)
      {
        while(strcmp(firsttemp,"run"))
          {
            fscanf(input,"%s ",&firsttemp);
            if (!strcmp(firsttemp,"thermo"))
              {
                fscanf(input,"%s ",&firsttemp);
                thermo[i]=atoi(firsttemp);
              }
          }
        fscanf(input,"%s ",&firsttemp);
        run[i]=atoi(firsttemp);
        if ((i==0) && (thermo[i]==0))
          {
            printf("thermo info not found in the file.. Exiting\n");
            exit(1);
          }
        if (thermo[i]==0)
          thermo[i]=thermo[i-1];
      }

    rewind(input);

    if (!strcmp(argv[3],"-m"))
      thermostylecounter=1;
    else if (!strcmp(argv[3],"-s"))
      thermostylecounter=0;
    else
      {
        printf("Couldn't recognize file format. Please use either -m or -s\n");
        exit(1);
      }

    if (thermostylecounter==1)
      for(i=0;i<iterations;i++)
        {
	  scan_first_part();
          scan_data_multi(i,argv[2]);
          while(strcmp(temp,"Loop"))
            fscanf(input,"%s ",&temp);
          for(j=0;j<19;j++)
            fgets(temp,80,input);
        }
    else
      for(i=0;i<iterations;i++)
        {
	  scan_first_part();
          scan_data_single(i,argv[2]);
          while(strcmp(temp,"Loop"))
            fscanf(input,"%s ",&temp);
          for(j=0;j<19;j++)
            fgets(temp,80,input);
        }
    fclose(input);
  }
 return 0;
}
