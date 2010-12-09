/*** modified on 12/09/2010 by Aidan Thompson               ***/
/*** modified on 10/21/2010 by avv/cfdrc                    ***/
/*** this code was given to LAMMPS by Sergey Zybin, CalTech ***/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Cutoff "Cutoff.dic"
#define Snapshots "bonds.reax"
#define Output1 "fra.out"
#define Output2 "mol.out"
#define Output3 "num.out"
#define Output4 "dipol.out"
#define Output5 "vibtemp.out"
#define Output6 "conc.out"

#define MaxNumBonds 10 //Max Number of bonds per atom
#define MaxMolTypes 80000 //Max Number of molecule types 
#define Natomtypes 4 //C,H,O,N

#define Carbon 0
#define Hydrogen 1
#define Nitrogen 3
#define Oxygen 2

// Prototypes for functions

void FileOpenError(FILE *fp,char *name);
void Allocation();
int AtomIndicator(char *symbol);
void SkipLine(FILE *fp);

int Natom,Nmolecule;//total numbers of atoms and molecules
int Nmoltype=0;//number of species
int *Atomdic;//this array stores atom type accessible by atom id-1.
int *MolName;//an array that stores number of atoms of each type for each molecule by using counters of Natomtypes*imoleculetypes+atomtype so if H2 is the second molecule entry the array will have 4=0, 5=2, 6=0, 7=0.
int *MolType;//array with the same info as MolName but accessed by Nmoltype, not molecule id.
int *Mol;//an array that stores the molecule id for each atom
int *NMol;//an array that stores the number of molecules of this species for this frame
double *Netc;//an array that stores charge
double Mass[Natomtypes];
double *rx,*ry,*rz,*Charge;
double *vx,*vy,*vz;
double COdic[Natomtypes][Natomtypes];
double latc_x[3],latc_y[3],latc_z[3];//lattice parameters of the box read from config.out
double Ilatc_x[3],Ilatc_y[3],Ilatc_z[3];//normalized box lengths
double Origin[3];//origin coordinates of the box
FILE *trj,*dpol,*vtemp,*conc;

main() {

  void GetAtomlist();
  void GetAtomlist_mod();
  void GetCOdict();
  void AnalyzeTrajectory();
  void PostProcess();
  
  GetAtomlist_mod(Snapshots);
  GetCOdict(Cutoff);
  AnalyzeTrajectory(Snapshots,Output3);
  PostProcess(Output1,Output2,Output3);

}

void GetAtomlist_mod(char *name) {
/*This subroutine reads the bonds.petn file to set the number of atoms.  It also checks that a full frame exists because it loops over the first frame to populate the Atomdic with the proper type for each atom id.*/
  FILE *fp;
  char buffer[1000];
  int  i,j, Nlink, iatom, itemp, id;
  double BO,dtemp;//BO is bond order
  
  fp=fopen(name,"r"); 
  FileOpenError(fp,name);  
  
 while(fscanf(fp,"%s",buffer)!=EOF) {

   if(strcmp(buffer,"particles")==0) {
      fscanf(fp,"%d",&Natom);
      printf("# of atoms : %d\n",Natom);
      Allocation();
    }
    if(strcmp(buffer,"q")==0) {    
     
     for(i=0;i<Natom;i++) {
        fscanf(fp,"%d",&iatom);
        iatom--;
        fscanf(fp,"%d",&itemp);
	itemp--;
	if (itemp == Carbon){
	   sprintf(buffer, "C");
	}else if(itemp == Hydrogen) {
	  sprintf(buffer, "H");
	}else if(itemp == Oxygen) {
	  sprintf(buffer, "O");  
        }else if(itemp == Nitrogen) {
	  sprintf(buffer, "N");
	} 
	 
        fscanf(fp,"%d",&Nlink);
        for(j=0;j<Nlink;j++) {
          fscanf(fp,"%d",&itemp);
          itemp--;    
        }
        fscanf(fp,"%d",&itemp);

        for(j=0;j<Nlink;j++) {
          fscanf(fp,"%lf",&BO);   
        }
        fscanf(fp,"%lf",&dtemp);
        fscanf(fp,"%lf",&dtemp);
        fscanf(fp,"%lf",&dtemp);
	
	id=AtomIndicator(buffer);
        Atomdic[iatom]=id;
      }
      if (i == Natom) goto Ending;  
  }  
}
Ending:;
}
  
  

void GetAtomlist(char *name) {
/*This subroutine reads the head of the configuration file to set the lattice parameters and the number of atoms.  It also checks that a full frame exists because it loops over the first frame to populate the Atomdic with the proper type for each atom id.*/
  
  FILE *fp;
  int i,iatom,count=0,itemp,id,flag;
  double dtemp;
  double a,b,c,d,e,f,g,h,ii,det;
  char buffer[1000];

  FILE *temp;
  fp=fopen(name,"r"); 
  FileOpenError(fp,name);

  while(fscanf(fp,"%s",buffer)!=EOF) {
    if(strcmp(buffer,"vectors:")==0) {
      fscanf(fp,"%s",buffer); 
      for(i=0;i<3;i++) fscanf(fp,"%lf",&latc_x[i]); 
      fscanf(fp,"%s",buffer); 
      for(i=0;i<3;i++) fscanf(fp,"%lf",&latc_y[i]); 
      fscanf(fp,"%s",buffer); 
      for(i=0;i<3;i++) fscanf(fp,"%lf",&latc_z[i]); 
    }
    else if(strcmp(buffer,"origin:")==0) {
      fscanf(fp,"%s",buffer); 
      for(i=0;i<3;i++) fscanf(fp,"%lf",&Origin[i]); 
    }
    else if(strcmp(buffer,"particles")==0) {
      fscanf(fp,"%d",&Natom);
      printf("# of atoms : %d\n",Natom);
      Allocation();
    }
    else if(strcmp(buffer,"positions")==0) {
      for(i=0,flag=0;i<Natom;i++) {
        fscanf(fp,"%d",&id);
        id--;
        fscanf(fp,"%s",buffer);
        iatom=AtomIndicator(buffer);
        Atomdic[id]=iatom;
        count++;
//        fscanf(fp,"%d",&itemp);
//        fprintf(temp,"\t%d",itemp);
//        fscanf(fp,"%lf",&dtemp);
//        fprintf(temp,"\t%lf",dtemp);
//        fscanf(fp,"%lf",&dtemp);
//        fprintf(temp,"\t%lf",dtemp);
//        fscanf(fp,"%lf",&dtemp);
//        fprintf(temp,"\t%lf",dtemp);
//        fscanf(fp,"%lf",&dtemp);
//        fprintf(temp,"\t%lf\n",dtemp);
	fgets(buffer,1000,fp);
        flag=1;
/*         fscanf(fp,"%lf",&dtemp); */
/*         fscanf(fp,"%lf",&dtemp); */
/*         fscanf(fp,"%lf",&dtemp); */
/*         fscanf(fp,"%lf",&dtemp); */
      }
      if(flag==1) break;
    }
  }
  if(Natom!=count) {
    printf("Can't read whole # of atoms.\n");
    exit(0);
  }


  a=latc_x[0]; b=latc_y[0]; c=latc_z[0];
  d=latc_x[1]; e=latc_y[1]; f=latc_z[1];
  g=latc_x[2]; h=latc_y[2]; ii=latc_z[2];
  det=a*e*ii-a*f*h-d*b*ii+d*c*h+g*b*f-g*c*e;
  if (det != 0) {
  Ilatc_x[0]=(e*ii-f*h)/det;
  Ilatc_x[1]=-(b*ii-c*h)/det;
  Ilatc_x[2]=(b*f-c*e)/det;

  Ilatc_y[0]=-(d*ii-f*g)/det;
  Ilatc_y[1]=(a*ii-c*g)/det;
  Ilatc_y[2]=-(a*f-c*d)/det;

  Ilatc_z[0]=(e*g-d*h)/det;
  Ilatc_z[1]=-(b*g-a*h)/det;
  Ilatc_z[2]=(a*e-b*d)/det;
  }

  close(fp);
}

void GetCOdict(char *name) {
/* This subroutine populates the cut-off dictionary for each bond and throws an error if any possible combination of atoms do not have a specified bond*/
  FILE *fp;
  int i,j,iatom,jatom;
  double cut;
#define BUFLEN 1000
  char buffer[BUFLEN];
  fp=fopen(name,"r"); 
  FileOpenError(fp,name);

  for(i=0;i<Natomtypes;i++) {
    for(j=0;j<Natomtypes;j++)  COdic[i][j]=-1;
  }

  while(fscanf(fp,"%s",buffer)!=EOF) {
    if (strcmp(buffer,"#") == 0) {
      SkipLine(fp);
      continue;
    }
    iatom=AtomIndicator(buffer);
    fscanf(fp,"%s",buffer);
    jatom=AtomIndicator(buffer);
    fscanf(fp,"%lf",&cut);
    COdic[iatom][jatom]=cut;
    COdic[jatom][iatom]=cut;
  }
  for(i=0;i<Natomtypes;i++) {
    for(j=0;j<Natomtypes;j++) {
      if(COdic[i][j]==-1) {
        printf("Can't fine a certain pair for cut off dictionary!\n");
        exit(0);
      }
    }
  }


  close(fp);
}
void AnalyzeTrajectory(char *name,char *name2) {
/***************************************************************************/
/*  This subroutine loops through the data both with the bond lists and the position files.  
    Step 1: It checks the timestep for the bond list and then reads in the position data ensuring that the timesteps match.
    Step 2: It loops over the bonds to identify molecules 
    Step 3: It loops over the molecules to determine name and species.
    Step 4: It calculates and writes dipole.out and vibtemp.out.
    Step 5: It writes num.out */
/**************************************************************************/

  FILE *fp,*out3;//fp is the file with the bond lists.  out3 refers to a scratch file to keep track of how many of each type of molecules exist.
  int i,j,k,l,iatom,jatom,Nlink,count,iatomtype,jatomtype;
/* jatom is the neighbor, iatom is the self, Nlink is the number of neighbors, count is a running index to identify order of atoms (e.g. atom 5 is on line 38)*/
  int frame;
  int itemp,Molid;
/*itemp is a generic variable as needed, Molid is the line of the atom (e.g. atom 5 is on line 38 means Molid[4]=37 because of the zero).*/

  double BO,dtemp;//BO is bond order
  char buffer[1000];
  int Cnct[MaxNumBonds];//connectivity table for a given atom

  int Nspec,flag_mol,flag_spec,flag_identity,nd,timestep;
/*  Nspec is the number of species during this frame.  flag_mol is a binary flag with 1 indicating that this atom belongs to this molecule.  flag_spec is a binary flag with 1 indicating that a new species has been found.  I have no idea what purpose flag_identity serves, but a 1 indicates that flag_spec is 1.
nd is a return value from CheckExistence.*/

  int Name[Natomtypes],MolTemp[Natom];
/*Name is a counter for how many atoms of a given type are in this particular molecule. MolTemp is an array that keeps track of molecule id by atom id.  For example, if atom 5 belongs to molecule 3 MolTemp[4]=3*/

  int CheckExistence(); 
  void GetAtomdata();
  void GetCOMDipole();

  fp=fopen(name,"r"); 
  FileOpenError(fp,name);
  out3=fopen(name2,"w"); 
  FileOpenError(out3,name2);

  frame=0;

  while(fscanf(fp,"%s",buffer)!=EOF) {
/* Step 1 */
    if(strcmp(buffer,"Timestep")==0) {
      fscanf(fp,"%d",&timestep);
      //      GetAtomdata(timestep);
    }

/* Step 2 */
    if(strcmp(buffer,"q")==0) {
      frame++;
      count=0; 
/*             zeroing out the counters        */
      for(i=0;i<Natom;i++) {
        MolTemp[i]=-1;
        NMol[i]=1;
      }
/*          Step 2A reads the atomid as iatom, the first itemp is a discard of the atom type, Nlink is the number of neighbors, and then the loop over all the bonded neighbors is stored in Cnct, and the last itemp is the ignored 0 to mark the end.*/
      for(i=0;i<Natom;i++) {
        fscanf(fp,"%d",&iatom);
        iatom--;
        if(MolTemp[iatom]==-1) {
          MolTemp[iatom]=count;
          count++;
        }
        fscanf(fp,"%d",&itemp);
        fscanf(fp,"%d",&Nlink);
        for(j=0;j<Nlink;j++) {
          fscanf(fp,"%d",&itemp);
          itemp--;
          Cnct[j]=itemp;
        }
        fscanf(fp,"%d",&itemp);

/*         Step 2B reads the bond orders from the bond file and assigns atoms to molecules.  The first check is if a bond exists.  The second check is if i and j already belong to the same molecule.  If they don't and j has no molecule, assign it to i's molecule.  If they do, then j was assigned to a different number so loop through the molecule list and set all the atoms on that molecule to the new i molecule.  Then read and discard the remainder of the line.*/

        for(j=0;j<Nlink;j++) {
          fscanf(fp,"%lf",&BO);
          jatom=Cnct[j];
          Molid=MolTemp[jatom];
          iatomtype=Atomdic[iatom];
          jatomtype=Atomdic[jatom];
          if(BO>COdic[iatomtype][jatomtype]) {
            if(Molid!=MolTemp[iatom]) {
              if(MolTemp[jatom]==-1) MolTemp[jatom]=MolTemp[iatom];
              else {
                for(k=0;k<Natom;k++) {
                  if(MolTemp[k]==Molid) MolTemp[k]=MolTemp[iatom];
                }
              }
            }
          } 
        }
        fscanf(fp,"%lf",&dtemp);
        fscanf(fp,"%lf",&dtemp);
        fscanf(fp,"%lf",&dtemp);
      }

/*       Step 3 builds the chemical formula for each molecule.  For each molecule (the for loop over i with count as the maximum), loop over all the atoms.*/

      for(i=0;i<MaxMolTypes;i++) Netc[i]=0.0;
      for(i=0,Nspec=0,Nmolecule=0;i<count;i++) {

 /*     Step 3A: if the atom belongs to this molecule, increment the appropriate type counter (e.g. we have another carbon), flip flag_mol to 1, set the Mol value to this molecule, and add the atomic charge to the running charge total. */

        for(j=0;j<Natomtypes;j++) Name[j]=0;
        for(j=0,flag_mol=0;j<Natom;j++) {
          if(MolTemp[j]==i) {
            jatom=Atomdic[j]; 
            Name[jatom]++;
            flag_mol=1; 
            Mol[j]=Nmolecule; 
            Netc[Nmolecule]+=Charge[j];
          }
        }

/*     Step 3B: If this molecule contains any atoms (flag_mol==1), then determine if it is a new species.  MolName is populated as new species are found.  In one of the least efficient matching mechanisms ever, we walk through all the known species names.  If no match is found, the flag_spec is set to 1.  A zero flag_spec means that we have a match, increment the appropriate species counter.  Otherwise, if we have a zero number of species or a flipped flag_spec, we have a new species and need to populate its name.  Then we increment the number of species as appropriate and always add a new molecule.*/

        if(flag_mol==1) {
          flag_identity=1;
          for(k=0;k<Nspec;k++) {
            flag_spec=0;
            for(l=0;l<Natomtypes;l++) {
              if(MolName[Natomtypes*k+l]!=Name[l]) flag_spec=1;
            }
            if(flag_spec==0) NMol[k]++;  
            flag_identity*=flag_spec; 
          }
          if(Nspec==0 || flag_identity==1) {
            for(l=0;l<Natomtypes;l++) {
              MolName[Natomtypes*Nspec+l]=Name[l];
            } 
            Nspec++; 
          }
          Nmolecule++;
        }//ends the if(flag_mol==1) loop

      }//ends the for i up to count loop

/*  Step 4 does the end of frame cleanup including the writes for COM */
      fprintf(out3,"%d\t",timestep);
      printf("Timestep = %d\n",timestep);
      fprintf(dpol,"%d\n",Nmolecule);
      fprintf(dpol,"TSTEP %d\n",timestep);
      fprintf(vtemp,"%d\n",Nmolecule);
      fprintf(vtemp,"TSTEP %d\n",timestep);
      fprintf(conc,"TSTEP %d\n",timestep);
      fprintf(conc,"%d\n",Nmolecule);
      GetCOMDipole(Nmolecule,timestep,frame);

/*  Step 5 prints the list of species and their quantities   */
      fprintf(out3,"%d\t",Nspec);
      for(i=0;i<Nspec;i++) {
        nd=CheckExistence(i);
        fprintf(out3,"%d ",nd);
        fprintf(out3,"%d\t",NMol[i]);
      }
      fprintf(out3,"\n");
    }
  }
  fflush(out3);
  fclose(out3);
  fclose(fp);
}

int CheckExistence(int id) {
/* The purpose of this function is to have a list of only the species that occur in this frame. */

  int i,j,num,molid,itemp;
  int flag,flag_identity=1;

  for(i=0;i<Nmoltype;i++) {
    flag=0;
    for(j=0;j<Natomtypes;j++) {
      molid=MolType[Natomtypes*i+j];
      if(molid!=MolName[Natomtypes*id+j]) flag=1;
    }
    if(flag==0) return i; 
    flag_identity*=flag;
  }

  for(i=0;i<Natomtypes;i++) {
    MolType[Natomtypes*Nmoltype+i]=MolName[Natomtypes*id+i];
  }
  Nmoltype++;

  return Nmoltype-1;
}

void FileOpenError(FILE *fp,char *name) {
  if(fp==NULL) {
    printf("file open error(%s).\n",name);
    exit(0);
  }
}

void Allocation() {
/* This subroutine sets up most of the important adjustable arrays based on the number of atoms.  It also opens the dipole and vibtemp files for writing.*/
  int flag=0;
  char name[100];
  Atomdic=(int *)malloc(Natom*sizeof(int));
  Mol=(int *)malloc(Natom*sizeof(int));
  MolName=(int *)malloc(Natom*Natomtypes*sizeof(int));
//20Jun07 JLB: I have made the absurdly huge array much smaller
//  MolType=(int *)malloc(Natom*MaxMolTypes*sizeof(int));
  MolType=(int *)malloc((Natomtypes+20)*MaxMolTypes*sizeof(int));
  NMol=(int *)malloc(Natom*sizeof(int));
  rx=(double *)malloc(Natom*sizeof(double));
  ry=(double *)malloc(Natom*sizeof(double));
  rz=(double *)malloc(Natom*sizeof(double));
  vx=(double *)malloc(Natom*sizeof(double));
  vy=(double *)malloc(Natom*sizeof(double));
  vz=(double *)malloc(Natom*sizeof(double));
  Charge=(double *)malloc(Natom*sizeof(double));
  Netc=(double *)malloc(MaxMolTypes*sizeof(double));
  
  if(Atomdic==NULL || NMol==NULL || Mol==NULL || MolType==NULL) flag=1;
  if(rx==NULL || ry==NULL || rz==NULL || Charge==NULL) flag=1;
  if(vx==NULL || vy==NULL || vz==NULL) flag=1;
  if(Netc==NULL) flag=1;

  if(flag==1) {
    printf("malloc error.\n");
    exit(0);
  }
  //  trj=fopen(Atomlist,"r");
  //  FileOpenError(trj,Atomlist);
  dpol=fopen(Output4,"w");
  FileOpenError(dpol,Output4);
  vtemp=fopen(Output5,"w");
  FileOpenError(vtemp,Output5);
  conc=fopen(Output6,"w");
  FileOpenError(conc,Output6);
}

int AtomIndicator(char *symbol) {
/* This subroutine identifies the atom type from the config file and sets the masses.  Why set the masses here, I don't know.*/

  int value;
  Mass[Hydrogen]=1.0079;
  Mass[Carbon]=12.0107;
  Mass[Oxygen]=15.9994;
  Mass[Nitrogen]=14.0067;
  if(strcmp(symbol,"H")==0) value=Hydrogen;
  else if(strcmp(symbol,"C")==0) value=Carbon;
  else if(strcmp(symbol,"O")==0) value=Oxygen;
  else if(strcmp(symbol,"N")==0) value=Nitrogen;
  else {
    printf("Can't recognize the atom type : %s.\n",symbol);
    exit(0);
  }

  return value;
}
void PostProcess(char *name,char *name2,char *name3) {

/*  This subroutine does the global lists after all the frames have been read.  out is the fra.out with the list of species and number of each per frame.  out2 is mol.out with the chemical identities of each species. out3 is num.out which is the scratch file with number of species per frame in a non-plottable format.*/

  int i,j;
  FILE *out,*out2,*out3;
  int step,molnum,id,num,itemp;
  int numcount[Nmoltype];

  out=fopen(name,"w");
  FileOpenError(out,name);
  out2=fopen(name2,"w");
  FileOpenError(out2,name2);
  out3=fopen(name3,"r");
  FileOpenError(out3,name3);

  fprintf(out,"Timestep\t");
  fprintf(out2,"%d\n\n",Nmoltype);

/*  To get the chemical identity of a species, this process loops over each atom type and writes the number of atoms of that type in the formula.    */
  for(i=0;i<Nmoltype;i++) {
    fprintf(out2,"%d\t",i);
    for(j=0;j<Natomtypes;j++) {
      itemp=MolType[Natomtypes*i+j];
      if(itemp!=0) {
        if(j==Carbon) {
          fprintf(out,"C");
          fprintf(out2,"C");
        }
        else if(j==Hydrogen) {
          fprintf(out,"H");
          fprintf(out2,"H");
        }
        else if(j==Oxygen) {
          fprintf(out,"O");
          fprintf(out2,"O");
        }
        else if(j==Nitrogen) {
          fprintf(out,"N");
          fprintf(out2,"N");
        }
        if(itemp!=1) {
          fprintf(out,"%d",itemp);
          fprintf(out2,"%d",itemp);
        }
      }
    }
    fprintf(out,"\t");
    fprintf(out2,"\n");
  }
  fprintf(out,"\n");

/*    for each molecule type, read the num.out scratch file and find the number of those molecules.  Then write it out in a nice column format suitable for plotting. */
  while(fscanf(out3,"%d",&step)!=EOF) {
    fscanf(out3,"%d",&num);
    for(i=0;i<Nmoltype;i++) numcount[i]=0;
    for(i=0;i<num;i++) {
      fscanf(out3,"%d",&id);
      fscanf(out3,"%d",&molnum);
      numcount[id]=molnum; 
    }
    fprintf(out,"%d\t",step);
    for(i=0;i<Nmoltype;i++) {
      fprintf(out,"%d\t",numcount[i]);
    }
    fprintf(out,"\n");
  }
}
void GetAtomdata(int tstep) {
/* This subroutine reads in the position data for a frame*/

  int i,iatom,count=0,itemp,id,flag;
  double dtemp;
  int step;
  char buffer[1000];

  while(fscanf(trj,"%s",buffer)!=EOF) {
    if(strcmp(buffer,"Timestep")==0) {
      fscanf(trj,"%d",&step);
      if(step!=tstep) {
        printf("The sequences are not matched between bond data file and coordinate file\n");
        printf("step in config file : %d\tStep in bond file : %d\n",tstep,step);
        exit(0);
      }
    }
    else if(strcmp(buffer,"positions")==0) {
      for(i=0,flag=0;i<Natom;i++) {
        fscanf(trj,"%d",&id);
        id--;
        fscanf(trj,"%s",buffer);
        iatom=AtomIndicator(buffer);
        Atomdic[id]=iatom;
        fscanf(trj,"%d",&itemp);
        fscanf(trj,"%lf",&dtemp);
        Charge[id]=dtemp;
        fscanf(trj,"%lf",&dtemp);
        rx[id]=dtemp;
        fscanf(trj,"%lf",&dtemp);
        ry[id]=dtemp;
        fscanf(trj,"%lf",&dtemp);
        rz[id]=dtemp;
        fscanf(trj,"%lf",&dtemp);
        vx[id]=dtemp;
        fscanf(trj,"%lf",&dtemp);
        vy[id]=dtemp;
        fscanf(trj,"%lf",&dtemp);
        vz[id]=dtemp;
        fscanf(trj,"%lf",&dtemp);
        flag=1;
      }
      if(flag==1) break;
    }
  }
}
void GetCOMDipole(int max,int tstep,int frame) {
/* max is the number of molecules */
  int i;
  void AnalCOMDipole();

  for(i=0;i<max;i++) AnalCOMDipole(i);

}
void AnalCOMDipole(int molid) {

  int i,j,ix,iy,iz,count,itemp;
  int ntot;
  double mass,chg;
  int atomid;
  double totmass;
  double x0,y0,z0,xtemp,ytemp,ztemp,x,y,z,r2,minr2,dx,dy,dz;
  double COMx,COMy,COMz;
  double vCOMx,vCOMy,vCOMz;
  double Dpolx,Dpoly,Dpolz,Dpol;
  double VibTempx,VibTempy,VibTempz,VibTemp;
  int Name[Natomtypes];

  double tempx,tempy,tempz;

  totmass=0.0;
  vCOMx=vCOMy=vCOMz=0.0;

/* Again in an inefficient way, we loop over all the atoms to find the ones that are members of this molecule (molid). Then calculate the desired quantities, including the fact that the first atom sets our standard for measurement.  Then, we write the quantities to the COM files.  

JB 26Sep07: to facilitate my write out of the chemical formula, redo the molecule determination here. */
  for(j=0;j<Natomtypes;j++) Name[j]=0;
  for(i=0,count=0;i<Natom;i++) {
    if(Mol[i]==molid) {
      atomid=Atomdic[i];
      Name[atomid]++;
      mass=Mass[atomid];
      chg=Charge[i];
      totmass+=mass;
      vCOMx+=mass*vx[i];
      vCOMy+=mass*vy[i];
      vCOMz+=mass*vz[i];
      if(count==0) {  	
//        printf("i,  rx[i], ry[i], rz[i] , mass %d %e %e %e %e\n", i,  rx[i], ry[i], rz[i] , mass);
        x0=rx[i];
        y0=ry[i];
        z0=rz[i];
        COMx=mass*x0;
        COMy=mass*y0;
        COMz=mass*z0;
      }
      if(count!=0) {
        minr2=1.0e8;
        for(ix=-1;ix<2;ix++) {
          for(iy=-1;iy<2;iy++) {
            for(iz=-1;iz<2;iz++) {
              xtemp=rx[i]+ix*latc_x[0]+iy*latc_y[0]+iz*latc_z[0];
              ytemp=ry[i]+ix*latc_x[1]+iy*latc_y[1]+iz*latc_z[1];
              ztemp=rz[i]+ix*latc_x[2]+iy*latc_y[2]+iz*latc_z[2];
              dx=xtemp-x0;
              dy=ytemp-y0;
              dz=ztemp-z0;
              r2=dx*dx+dy*dy+dz*dz;
              if(r2<minr2) {
                minr2=r2;
                x=xtemp;
                y=ytemp;
                z=ztemp;
              }
            }
          }
        }
        rx[i]=x;
        ry[i]=y;
        rz[i]=z;
//        printf("AAA i,  rx[i], ry[i], rz[i] , mass %d %e %e %e %e\n", i,  rx[i], ry[i], rz[i] , mass);
        COMx+=mass*x;
        COMy+=mass*y;
        COMz+=mass*z;
      }
      count++;
    }
  }
//  printf ("AAA totmass %e\n",totmass);
  COMx/=totmass;
  COMy/=totmass;
  COMz/=totmass;

  vCOMx/=totmass;
  vCOMy/=totmass;
  vCOMz/=totmass;

  VibTempx=VibTempy=VibTempz=0.0;
  for(i=0,ntot=0,Dpolx=0.0,Dpoly=0.0,Dpolz=0.0;i<Natom;i++) {
    if(Mol[i]==molid) {
      atomid=Atomdic[i];
      mass=Mass[atomid];
      chg=Charge[i];
      totmass+=mass;
      x=rx[i]-COMx;
      y=ry[i]-COMy;
      z=rz[i]-COMz;
      Dpolx+=chg*x;
      Dpoly+=chg*y;
      Dpolz+=chg*z;
      VibTempx+=mass*(vx[i]-vCOMx)*(vx[i]-vCOMx);
      VibTempy+=mass*(vy[i]-vCOMy)*(vy[i]-vCOMy);
      VibTempz+=mass*(vz[i]-vCOMz)*(vz[i]-vCOMz);
      ntot++;
    }
  }
  VibTempx/=(double)ntot;
  VibTempy/=(double)ntot;
  VibTempz/=(double)ntot;

  VibTemp=(VibTempx+VibTempy+VibTempz)/3.0;

  VibTemp*=1.2027E+06;


  COMx-=Origin[0]; COMy-=Origin[1]; COMz-=Origin[2];
//  printf("VVVV %e %e %e\n",Ilatc_x[0],Ilatc_x[1],Ilatc_x[2]);
  tempx=Ilatc_x[0]*COMx;
  tempx+=Ilatc_x[1]*COMy;
  tempx+=Ilatc_x[2]*COMz;
 
  tempy=Ilatc_y[0]*COMx;
  tempy+=Ilatc_y[1]*COMy;
  tempy+=Ilatc_y[2]*COMz;

  tempz=Ilatc_z[0]*COMx;
  tempz+=Ilatc_z[1]*COMy;
  tempz+=Ilatc_z[2]*COMz;

  tempx-=floor(tempx);
  tempy-=floor(tempy);
  tempz-=floor(tempz);

  COMx=latc_x[0]*tempx+latc_y[0]*tempy+latc_z[0]*tempz;
  COMy=latc_x[1]*tempx+latc_y[1]*tempy+latc_z[1]*tempz;
  COMz=latc_x[2]*tempx+latc_y[2]*tempy+latc_z[2]*tempz;
//  printf ("LATICE %e %e %e %e %e %e \n",latc_x[0], tempx, latc_y[0], tempy,latc_z[0],tempz );
  COMx+=Origin[0]; COMy+=Origin[1]; COMz+=Origin[2];

  Dpolx/=0.2082267;
  Dpoly/=0.2082267;
  Dpolz/=0.2082267;

  Dpol=Dpolx*Dpolx+Dpoly*Dpoly+Dpolz*Dpolz;
  Dpol=sqrt(Dpol);
//  printf("BBB %e %e %e\n",COMx,COMy,COMz);
  fprintf(dpol,"He %e %e %e ",COMx,COMy,COMz);
  fprintf(dpol,"%e %e %e %e\n",Dpol,Dpolx,Dpoly,Dpolz);
  fprintf(vtemp,"He %e %e %e",COMx,COMy,COMz);
  fprintf(vtemp,"%e %e %e %e\n",VibTemp,vCOMx,vCOMy,vCOMz);

/* JB 26Sep07: I am adding the routines that will write the chemical formula and the COM to a separate file for use with the composition profile.  Step 1: determine the chemical formula for this molecule.*/
    for(j=0;j<Natomtypes;j++) {
      itemp=Name[j];
      if(itemp!=0) {
        if(j==Carbon) {
          fprintf(conc,"C");
        }
        else if(j==Hydrogen) {
          fprintf(conc,"H");
        }
        else if(j==Oxygen) {
          fprintf(conc,"O");
        }
        else if(j==Nitrogen) {
          fprintf(conc,"N");
        }
        if(itemp!=1) {
          fprintf(conc,"%d",itemp);
        }
      }
    }
  fprintf(conc,"\t");
  fprintf(conc,"%e %e %e\n",COMx,COMy,COMz);

}

void SkipLine(FILE *fp) {
  int c;
  while ((c = getc(fp)) != EOF)
    if (c == '\n') break;
}

