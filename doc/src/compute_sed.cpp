/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Tianli Feng, Divya Chalise, Xiulin Ruan
------------------------------------------------------------------------- */

#include <mpi.h>
#include "stdlib.h"
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include "memory.h"
#include "compute_sed.h"
#include "atom.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;


// global variables
const double kB = 1.3806E-23; // Boltzmann constant
const double J_to_eV = 1.60219E-19;
const double m_to_A = 1E-10;	

/* ---------------------------------------------------------------------- */



ComputeSed::ComputeSed(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  kpt(NULL), eigenv(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal compute sed command");	
	
  std::ifstream setup_in("info.txt");
  std::ifstream kpnts("kpoints.txt");
  std::ifstream recpvec("recipvec.txt");
  std::ifstream cellmap("cellbasismap.txt");
  std::ifstream coord("lammps.txt");
  std::ifstream eigen_txt("eigenvectors.txt");
	
	
	
  if (!setup_in.is_open()) {
    error->all(FLERR,"Error opening the file info.txt");
  }
  if (!kpnts.is_open()) {
    error->all(FLERR,"Error opening the file kpoints.txt");
  }
  if (!recpvec.is_open()) {
    error->all(FLERR,"Error opening the file recipvec.txt");
  }
  if (!cellmap.is_open()) {
    error->all(FLERR,"Error opening the file cellbasismap.txt");
  }
  if (!coord.is_open()) {
    error->all(FLERR,"Error opening the coordinate file lammps.txt");
  }
	
  double xtmp[3], frac[3], rec[3][3], chg;
  int i,j,l,s;
  char tempstr[255];
  int nothing,id,nk, prepared;
  struct timeval beg;
	
  for(i=0;i<3;i++) {
		recpvec>>rec[0][i]>>rec[1][i]>>rec[2][i];
  }
	
  kpnts>>nk;// read in number of kpoints

  memory->create(kpt,nk,3,"sed:kpt");
	
  for(i=0;i<nk;i++) {
    for(j=0;j<3;j++) {
      kpnts>>frac[j];
    }
    for(l=0;l<3;l++) {
      kpt[i][l]=0;
      for(s=0;s<3;s++) {
        //put atoms coord back to Ang for convenience
		kpt[i][l]+=rec[l][s]*frac[s];
      }
    }
  }
  setup_in>>natoms;
  setup_in>>nbasis;
	
  eigenv = new EIGEN [nk*(nbasis*3)*nbasis];
  at2bs = new int [natoms];	//map atom id to basis atoms
  indat2cell = new int [natoms];
  mass = new double [nbasis];
	
  for(i=0;i<nbasis;i++)
    setup_in>>mass[i];
  setup_in>>nlines;
  setup_in>>nlinemin;
  nlinemax = nlinemin + nlines;
  setup_in>>delt;
  setup_in>>gap;
	
  nkst = 0;
  nked = nk;
  nkpp = nked - nkst; // #kpoints per processor
  memory->create(ratom,natoms,3,"sed:ratom");
	
  for(i=0;i<nk*3*nbasis*nbasis;i++) {
    for(l=0;l<3;l++)
      eigen_txt>>eigenv[i].re[l]>>eigenv[i].im[l];
  }
  int ncell=natoms/nbasis;
  int ind[nbasis][ncell];
  int cellid, basisid;
  int indexcell[ncell];
  int cellindex[natoms],basisindex[natoms];
  for (i=0;i<natoms;i++) {
    cellmap>>id>>cellid>>basisid;
    basisid=basisid-1;
    ind[basisid][cellid-1]=id-1;
    cellindex[i]=cellid-1;
    basisindex[i]=basisid;
    at2bs[id-1]=basisid;
	//find the id of the first basis atom in this cell C
    if (basisid==0) indexcell[cellid-1]=id-1;
  }
  for (i=0;i<natoms;i++)
    indat2cell[i]=i;
  for (i=0;i<natoms;i++) {
    if (basisindex[i]!=0) {
	  //find the id of the first basis atom sharing the same cell with A	
      indat2cell[i]=indexcell[cellindex[i]];
    }
  }
	/*read the initial positions from lammps.txt ; 
	should use averaged atomic positions instead if possible */
  i=0;// check if lammps.txt is in the right format
  char *pch;
  while(1) {
    coord.getline(tempstr,256);
    if(!(strstr(tempstr,"Atoms"))) continue;
    coord.getline(tempstr,256);
    coord.getline(tempstr,256);
    break;
  }
  pch= strtok (tempstr," ,\t");
  while(pch!=NULL) {
    pch= strtok (NULL," ,\t");
    i++;
  }
  if(i!=7) {
    error->all(FLERR,"Wrong lammps.txt format");
    exit(0);
  }
  else {
    coord.clear();
    coord.seekg(0, std::ios::beg);
  }
  prepared=0;
  while(!prepared) {
    coord.getline(tempstr,256);
    if(!(strstr(tempstr,"Atoms"))) continue;
    coord.getline(tempstr,256);
    for(i=0;i<natoms;i++) {
      coord>>id>>nothing>>j>>chg>>xtmp[0]>>xtmp[1]>>xtmp[2];
      for(l=0;l<3;l++) {
        ratom[id-1][l]=xtmp[l];// id starts from 0
      }
    }
    prepared=1;
  }
 
  setup_in.close();
  kpnts.close();
  recpvec.close();
  cellmap.close();
  coord.close();
  eigen_txt.close();
	
  int veclength=(nbasis*3*2)*nked;
  vector_flag = 1;
  size_vector = veclength+1;
  extvector = 1;
  vector = new double[veclength+1];
}
  
ComputeSed::~ComputeSed()
{
  delete [] eigenv;
  delete [] mass;
  delete [] at2bs;
  delete [] indat2cell;
  memory->destroy(kpt);
  memory->destroy(ratom);
  delete [] vector;
}
  
/* ---------------------------------------------------------------------- */  
void ComputeSed::init()
{
	
}
/* ---------------------------------------------------------------------- */

void ComputeSed::compute_vector()
{
  int k,m; //index
  int i,l,n,p,q,count; //looping variables
  int id, atid,brn,ine; //brn=branch, ine=index of eigenvector
  double eiqr, cs, sn,sqm; //sqm=sqrt(mass)
  double sre[nkpp][nbasis*3], sim[nkpp][nbasis*3];
  double vntemp[3*natoms], vnf[3*natoms];
  double	t0 = 1.0E-12;
	
  //Information from MD 
  double **v=atom->v;
  int tstep=update->ntimestep;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  m =	tstep/gap+1;

  //time calculator
  vector[0]=(tstep+1-nlinemin)*t0*delt*1e12; 
	 
  //initialization
  for (count=0;count<3*natoms;count++) {
    vntemp[count] =0.0;	 
    vnf[count] =0.0;
  }
  for(p=0;p<nkpp;p++) {
    for(brn=0;brn<nbasis*3;brn++) {
      sre[p][brn]=0.0;
      sim[p][brn]=0.0;
    }
  }
	
  if (m>=nlinemin) {
	//Assigning velocities according to atom id
    for (count=0;count<nlocal;count++) {
      atid=tag[count];
	  //because tag starts from 1 and vntemp has to start from 0 
      for(l=0;l<3;l++) vntemp[3*(atid-1)+l] =v[count][l] ;	 	
    }
		
    MPI_Allreduce(vntemp,vnf,3*natoms,MPI_DOUBLE,MPI_SUM,world);
    //vnf=vntemp;

    for(p=nkst;p<nked;p++) {
      for(q=0;q<natoms;q++) {
        //the only knowledge of tag of basis atom is needed, not cell info.
		k=at2bs[q];
        sqm=sqrt(mass[k]/natoms*nbasis);//ncell=natoms/nbasis
        id=q;//id of atom, starts from 0
        eiqr=0.;
        for(l=0;l<3;l++) eiqr+=kpt[p][l]*ratom[indat2cell[id]][l];
        cs=cos(eiqr);
        sn=sin(eiqr);
        for(brn=0;brn<nbasis*3;brn++) {
          ine=p*nbasis*3*nbasis+brn*nbasis+k;
          for(n=0;n<3;n++) {
            sre[p-nkst][brn]+=vnf[3*id+n]*(eigenv[ine].re[n]*cs-eigenv[ine].im[n]*sn)*sqm;
            sim[p-nkst][brn]+=vnf[3*id+n]*(-eigenv[ine].re[n]*sn-eigenv[ine].im[n]*cs)*sqm;		
          }
        }
      }//q
						
      for(brn=0;brn<nbasis*3;brn++) {
        vector[(nbasis*3*2)*p+2*brn+1]=sre[p-nkst][brn];
        vector[(nbasis*3*2)*p+2*brn+2]=sim[p-nkst][brn];
        sre[p-nkst][brn]=0.0;
        sim[p-nkst][brn]=0.0;			
      }
							
    }// for p	
		
  } //if 
}
		

		
		
	 	
		