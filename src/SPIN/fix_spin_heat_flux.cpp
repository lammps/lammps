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
   Contributing authors:Yanguang Zhou (UCLA, USA). v1.0
   Data: 05.08.2017
   Email_address: yanguangzhou@ucla.edu
   
   Mentions:
   1, This programm can be used to calculate the thermal conductivity of magnons
   2, If you wan to to use this code, please cite our papers.
   3, If you have any problems, please send me email.
------------------------------------------------------------------------- */
#ifdef FIX_CLASS
FixStyle(spin/heat/flux,FixSpinheatFlux)
#else
#include <math.h>
#include <string.h>
#include "fix_spin_heat_flux.h"
#include "atom.h"
#include "update.h"
#include "atom.h"
#include "memory.h"
#include "compute.h"
#include "modify.h"
#include "force.h"
#include "error.h"
#include "lattice.h"
#include "domain.h"
#include "update.h"
#include "group.h"
#include "stdlib.h"
#include "output.h"
#include "string.h"
#include "input.h"
#include <iostream> //The oupt function

const double PI    = 3.14159265358979; // define the global variables
#define INVOKED_PERATOM 8
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSpinheatFlux::FixSpinheatFlux(LAMMPS *lmp, int narg, char **arg) : 
Fix(lmp, narg, arg) //The definition has been done in lammps
{
 scalar_flag = 1;
 global_freq = 1;
 extscalar = 1;
 /* Parameter description
 arg[0       ] ID
 arg[1       ] group-ID-1
 arg[2       ] fix name
 arg[3       ] group-ID-2  
 arg[4       ] call every N time steps
               Here,I choose to output every N time steps
 arg[5       ] alpha
 arg[6       ] gamma
 arg[7       ] sigma
 arg[8       ] cut_off
 arg[last    ] The output file name*/
 
 int n = strlen(arg[3]) + 1;
 group2 = new char[n];
 strcpy(group2,arg[3]);  
 jgroup = group->find(group2);
 if (jgroup == -1) error->all(FLERR,"Compute group/group group ID does not exist"); 
 if (narg < 9) error->all(FLERR,"Illegal fix heat current spectrum command: number of arguments < 9");

nevery  = utils::inumeric(FLERR,arg[4],false,lmp);
if ((nevery<1)) error->all(FLERR,"Illegal fix heat current spectrum command");

 MPI_Comm_rank(world,&me);   // Making the current processor as the head node
 MPI_Comm_size(world,&nump);
 
 atompgroupi = group->count(igroup);
 atompgroupj = group->count(jgroup);  
 jgroupbit = group->bitmask[jgroup];
 
 alpha  = atof(arg[5]);
 gamma  = atof(arg[6]);
 sigma  = atof(arg[7]);
 r_cut  = atof(arg[8]);
  
 Jd  = new int [nump];
 Jd0 = new int [nump];
 mpi_Jd = new int [nump];
 for(j=0; j<nump; j++) {
   Jd[j] = 0;
   Jd0[j] = 0;
   mpi_Jd[j] = 0;}
   
 Jin  = new double [6*atompgroupj]; // changed 4*atompgroupj to 6*atompgroupj 
 Jout = new double [6*atompgroupj]; // changed 4*atompgroupj to 6*atompgroupj
 for(j=0; j<6*atompgroupj; j++) {
   Jin[j] = 0.0;
   Jout[j] = 0.0;}    

//========================================================================//
//==== Output file pointers, The head node is used as the output node ====//

 if( me==0 ){
   char *Modename1;
   n = strlen(arg[narg-1]) + 1;
   Modename1 = new char[n]; //The name of the eachmode VS. files
   strcpy(Modename1,arg[narg-1]);
   fp_W  = fopen(Modename1,"wt");  
   delete [] Modename1;}

// nvalid = next step on which end_of_step does something
 lastindex = -1;
 nsample = 0;
 nvalid = nextvalid();}

/* ---------------------------------------------------------------------- */

FixSpinheatFlux::~FixSpinheatFlux(){
 delete [] Jd;
 delete [] Jd0;
 delete [] mpi_Jd;
 delete [] Jin;
 delete [] Jout;
 // Close files
 if(me==0){
   fclose(fp_W);}  
}

/* ---------------------------------------------------------------------- */

int FixSpinheatFlux::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */
void FixSpinheatFlux::init()
{
 if (nvalid < update->ntimestep) {
   lastindex = -1;
   nvalid = nextvalid(); }}

/*---------------------------------------------------------------------- */
void FixSpinheatFlux::setup(int vflag)
{
 end_of_step();}

/* --------------------------------------------------------------------- */
void FixSpinheatFlux::end_of_step()
{    	
 ntimestep = update->ntimestep;
 if(ntimestep != nvalid) return;    
 x   = atom->x; 
 sp  = atom->sp;  
 fm  = atom->fm;
 mask  = atom->mask;	    
 
// Initialization: make sure that all the arrary is 0 
 J_mag = 0.0;
 for(j=0; j<nump; j++) {
   Jd[j] = 0;
   Jd0[j] = 0;
   mpi_Jd[j] = 0;}  

 for(j=0; j<6*atompgroupj; j++) {
   Jin[j] = 0.0;
   Jout[j] = 0.0;}    
      
// Make all the atoms in the J group locate in Process  
 j = 0;
 for(n=0; n<atom->nlocal; n++){	   
   if((mask[n] & jgroupbit)) {
     j++; } }
 mpi_Jd[me] = j; 	

 MPI_Allreduce(mpi_Jd,Jd0,nump,MPI_INT,MPI_SUM,world);	
 	
 if(me==0){
   for(i=1;i<nump;i++) {
     temp = 0;
     for(j=0;j<i;j++){temp += Jd0[j];}
     mpi_Jd[i] = temp;}	  
   mpi_Jd[0] = 0;}
 if(me!=0){
     for(i=0;i<nump;i++) mpi_Jd[i] = 0; }
	   
 MPI_Allreduce(mpi_Jd,Jd,nump,MPI_INT,MPI_SUM,world);
   
// Spin force and spin in the I and J group
 i = 0;
 for (n = 0; n < atom->nlocal; n++) {
   if (mask[n] & jgroupbit) {
    Jin[i+Jd[me]] = x[n][0];
    Jin[i+Jd[me]+atompgroupj] = x[n][1];
    Jin[i+Jd[me]+2*atompgroupj] = x[n][2];
    Jin[i+Jd[me]+3*atompgroupj] = sp[n][0];
    Jin[i+Jd[me]+4*atompgroupj] = sp[n][1]; // add spy component
    Jin[i+Jd[me]+5*atompgroupj] = sp[n][2]; // add spz component
	i++; }}
	
 MPI_Allreduce(Jin,Jout,6*atompgroupj,MPI_DOUBLE,MPI_SUM,world); // changed 4*atomgroupcountj to 6*atomgroupcountj
 
// Heat flux from magnon: Jij*(omega_i cross product spin_i)*spin_j 
 jv = 0.0;
 for (i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      for (int j = 0; j < atompgroupj; j++){
        dx = x[i][0] - Jout[j];    
        dy = x[i][1] - Jout[j+atompgroupj];
        dz = x[i][2] - Jout[j+2*atompgroupj];
      
        if(dx > domain->xprd_half) {dx = dx - domain->xprd;}
	    else if (dx < -domain->xprd_half) {dx = dx + domain->xprd;}
    
        if(dy > domain->yprd_half) {dy = dy - domain->yprd;}
	    else if (dy < -domain->yprd_half) {dy = dy + domain->yprd;}
 
        if(dz > domain->zprd_half) {dz = dz - domain->zprd;}
      else if (dz < -domain->zprd_half) {dz = dz + domain->zprd;}
      
        r_delta = sqrt(dx*dx + dy*dy + dz*dz)/sigma;
        if (r_delta*sigma < r_cut) hevi_step = 1;
        else hevi_step = 0;
        Jij = 4*alpha*r_delta*r_delta*(1-gamma*r_delta*r_delta)*exp(-r_delta*r_delta)*hevi_step;
      
        jv += Jij * (
      (fm[i][1]*sp[i][2]-fm[i][2]*sp[i][1]) * Jout[j+3*atompgroupj]      // x
    + (fm[i][2]*sp[i][0]-fm[i][0]*sp[i][2]) * Jout[j+4*atompgroupj]      // Added y component
    + (fm[i][0]*sp[i][1]-fm[i][1]*sp[i][0]) * Jout[j+5*atompgroupj]);}}} // Added z component

 MPI_Allreduce(&jv,&J_mag,1,MPI_DOUBLE,MPI_SUM,world);
 if(me==0){
    fprintf(fp_W,"%d %8.6f\n",nvalid,J_mag);
	if (ntimestep%(10*nevery)==0) fflush(fp_W); }
 
 nvalid += nevery; 
}

/* ---------------------------------------------------------------------- */
bigint FixSpinheatFlux::nextvalid(){
 nvalid = update->ntimestep;
 //if ((nvalid % nevery) || (nvalid ==0)) 
 nvalid = (nvalid/nevery + 1)*nevery;
 return nvalid;  }

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

// Added functionality to make the heat flux a usuable variable that can be passed 
// to other functions. 
double FixSpinheatFlux::compute_scalar()
{

	return J_mag;
}

#endif


















