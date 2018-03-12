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
   Contributing authors:
   Roberto Guerra, e-mail: roberto.guerra at unimi dot it
   Silvia Bonfanti, e-mail: silviafisica at gmail dot com

   This is the potential described in:
     J.Chem. Phys. 140, 104106 (2014) ; DOI: 10.1063/1.4867272
   and reparametrized in:
     J. Phys. Chem. C, 2017, 121 (41), pp 22826–22835 ; DOI: 10.1021/acs.jpcc.7b07091

   Assumptions:
   - The structure is assumed to be stacked graphitic layers.
   - Each layer must have different atom types
   - Pair interaction is inter-layer, i.e. between different layers:
       pair_coeff  1*2 3*4  ilp  ILP.par.2017  B N B N
   - Each atom must have one to three INTRA-layer neighbors within a cutoff of sqrt(2.5).
   - The intra-layer neighbors of each atom must be static: no melting!
   - The normal vector of atom i is defined by the normalized cross product (a-b)^(a-c),
     where a-b and a-c are the bond vectors among the three nearest neighbors of i.
     If i has only two intra-neighbors, the cross product is (i-a)^(i-b).
     if i it has just one intra-neighbor j, then j is expected to have three neighbors and
     the normal vector of j will be used instead.
   - Coulomb interaction among partial charges is currently not implemented.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_ilp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include <iostream>
#include "atom_vec.h"
#include "domain.h"

using namespace LAMMPS_NS;
// to print
using namespace std;

#define MAXLINE 1024
#define TOL 1.0e-9
#define DELTA 4


/* ---------------------------------------------------------------------- */
PairILP::PairILP(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;

  // initialize element to parameter maps
  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;
  map = NULL;

  // never compute energy offset
  offset_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairILP::~PairILP()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(offset);
  }

  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);
  if (allocated) delete [] map;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairILP::init_style()
{

  if(force->newton_pair == 0)
    error->all(FLERR,"Pair style ilp requires newton pair on");

  // map array is required to properly map from globalid to localid
  if(atom->map_style==0)
    error->all(FLERR,"Pair style ilp requires atom_modify map array");

  int irequest = neighbor->request(this,instance_me);

}

/* ----------------------------------------------------------------------
   workhorse routine that computes pairwise interactions
---------------------------------------------------------------------------*/

void PairILP::compute(int eflag, int vflag)
{
  #define ILP_d 15.0
  int i,j,ii,jj,ll,kk,a,b,c,iparam_ij;
  double xtmp,ytmp,ztmp,tmp,delx,dely,delz,evdwl;
  double r,rsq,r3,r4,r6,aasimp,asimp,bsimp,csimp,dsimp,esimp;
  double ndotr_ij,ndotr_ji,rho2_ij,rho2_ji;
  double temp_b,temp_c,sumCfrho,dsimp_ij,dsimp_ji;
  static int once=1,*intra_numneigh,**intra_neigh;
  static double *nx=NULL,*ny=NULL,*nz=NULL,*nnorm=NULL,**bx=NULL,**by=NULL,**bz=NULL; //https://stackoverflow.com/questions/12134315/realloc-on-null-valued-or-undefined-pointer
  // consider using: static double[2] *bx=NULL,*by=NULL,*bz=NULL; https://stackoverflow.com/questions/49014310
  static double TAP_A,TAP_B,TAP_C,TAP_D,dTAP_A,dTAP_B,dTAP_C,dTAP_D;
  double PBCx,PBCy,PBCz;

  PBCx = domain->boxhi[0] - domain->boxlo[0];
  PBCy = domain->boxhi[1] - domain->boxlo[1];
  PBCz = domain->boxhi[2] - domain->boxlo[2];

  evdwl = 0.;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type  = atom->type;
  int nlocal = atom->nlocal; // number of atoms in this thread
  int nghost = atom->nghost; // number of ghost atoms in this thread
  int ntotal = nlocal+nghost;// max # of owned+ghost in arrays on this proc
  static int ntotal0 = 0;    // tracks changes of ntotal
  int inum,jnum,itype,jtype,*ilist,*jlist,*numneigh,**firstneigh;

  inum = list->inum;         //number of atoms in the inter-neighbors list (should be equal to nlocal)
  ilist = list->ilist;       //inter-neighbors list
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  //Here we alloc and build the static intra-neighbors list and other arrays for internal use
  //BEGIN ONCE
  if(once){

    int natoms; // total number of atoms of the whole system
    int nprocs=comm->nprocs;   // number of threads
    double *xx,*yy,*zz,Rcut;
    int imin,imax,*tt,*mm,global_id,max_global_id;

    //finding the maximum global index
    imax=0;
    for(i=0;i<ntotal;i++){
      global_id = atom->tag[i];
      if(global_id>imax) imax=global_id;
    }

    MPI_Reduce(&imax, &max_global_id, 1, MPI_INT, MPI_MAX, 0, world);
    MPI_Bcast(&max_global_id, 1, MPI_INT, 0, world);
    natoms=max_global_id;

    Rcut = cut_global;
    if(comm->me == 0) cout << "ILP: Rcut = " << Rcut << endl;

    TAP_A  = 20./pow(Rcut,7);
    TAP_B  = 70./pow(Rcut,6);
    TAP_C  = 84./pow(Rcut,5);
    TAP_D  = 35./pow(Rcut,4);

    dTAP_A = 7.*TAP_A;
    dTAP_B = 6.*TAP_B;
    dTAP_C = 5.*TAP_C;
    dTAP_D = 4.*TAP_D;

    if(comm->me == 0){
      cout << "ILP: nlocal0 = " << nlocal << " , nghost0 = " << nghost << " , ntotal0 = " << ntotal << " , natoms = " << atom->natoms << " , max global id = " << max_global_id << endl;
      cout << "ILP: box = " << PBCx << " x " << PBCy << " x " << PBCz << endl;
      cout << "ILP: building the inter-layer neighbors list..." << endl;
    }
    MPI_Barrier(world);

    if( (xx=(double *)calloc(natoms,sizeof(double))) == NULL ) { error->all(FLERR,"allocation error."); }
    if( (yy=(double *)calloc(natoms,sizeof(double))) == NULL ) { error->all(FLERR,"allocation error."); }
    if( (zz=(double *)calloc(natoms,sizeof(double))) == NULL ) { error->all(FLERR,"allocation error."); }
    if( (tt=(int    *)calloc(natoms,sizeof(int)   )) == NULL ) { error->all(FLERR,"allocation error."); }
    if( (mm=(int    *)calloc(natoms,sizeof(int)   )) == NULL ) { error->all(FLERR,"allocation error."); }

    if (comm->me == 0) cout << "  1) memory allocated" << endl;

    /**************NOTE: mostly copied from WriteData::atoms() in write_data.cpp **************/

    // communication buffer for all my Atom info
    // max_size = largest buffer needed by any proc

    int ncol = atom->avec->size_data_atom + 3;
    int sendrow = atom->nlocal;
    int maxrow;
    MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

    double **buf;
    if (comm->me ==0) memory->create(buf,MAX(1,maxrow ),ncol,"write_data:buf");
    else              memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

    // pack my atom data into buf
    atom->avec->pack_data(buf);

    // proc 0 pings each proc, receives its chunk, collects to buffin
    // all other procs wait for ping, send their chunk to proc 0
    int tmp,recvrow;

    if (comm->me == 0) {
      MPI_Status status;
      MPI_Request request;

      for (int iproc = 0; iproc < nprocs; iproc++) {
        if (iproc){
          MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
          MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
          MPI_Wait(&request,&status);
          MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
          recvrow /= ncol;
        }
        else recvrow = sendrow;

        //atom->avec->write_data(fp,recvrow,buf);
        for(ii=0;ii<recvrow;ii++){
             j =(tagint) ubuf(buf[ii][0]).i - 1; // atom id starts from 1, arrays starts from 0
          mm[j]=(tagint) ubuf(buf[ii][1]).i;     // NOTE: depends on atom_style: absent for atomic/charge, 1 for full
          tt[j]=(tagint) ubuf(buf[ii][2]).i;     // NOTE: depends on atom_style: 1      for atomic/charge, 2 for full
          xx[j]=buf[ii][4];                      // NOTE: depends on atom_style: 3      for atomic/charge, 4 for full
          yy[j]=buf[ii][5];                      //
          zz[j]=buf[ii][6];                      //
        }

      }

    }
    else {
      MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
      MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
    }

    memory->destroy(buf);
    /*******************************************************************/
    MPI_Barrier(world);
    if (comm->me == 0) cout << "  2) coordinates collected" << endl;

    // we send xx,yy,zz from proc0 to all the other procs
    MPI_Bcast(xx, natoms, MPI_DOUBLE, 0, world);
    MPI_Bcast(yy, natoms, MPI_DOUBLE, 0, world);
    MPI_Bcast(zz, natoms, MPI_DOUBLE, 0, world);
    MPI_Bcast(tt, natoms, MPI_INT   , 0, world);
    MPI_Bcast(mm, natoms, MPI_INT   , 0, world);
    MPI_Barrier(world);
    if (comm->me == 0) cout << "  3) coordinates shared" << endl;

    // now that we have ALL the atom coordinates, we
    // allocate the intra-neighbor list...
    intra_numneigh=(int *)calloc(natoms+1,sizeof(int));
    intra_neigh=(int **)malloc((natoms+1)*sizeof(int*));
    for (ii = 0; ii < (natoms+1); ii++) {
      intra_neigh[ii]=(int *)malloc(3*sizeof(int));
    }

    // ...and we fill it
    for (i = 0; i < natoms; i++) {
      xtmp  = xx[i];
      ytmp  = yy[i];
      ztmp  = zz[i];
      itype = tt[i];
      if( itype==0 || map[itype] == -1 ) continue; // atoms not belonging to this pair_potential

      for (j = 0; j < natoms; j++) {
        if( i==j || itype==0 || map[itype] == -1 || mm[i]!=mm[j] ) continue; // atom j not in the same layer as i

        delx = xtmp - xx[j]; delx -= round(delx/PBCx)*PBCx;
        dely = ytmp - yy[j]; dely -= round(dely/PBCy)*PBCy;
        delz = ztmp - zz[j]; delz -= round(delz/PBCz)*PBCz;

        rsq = delx*delx + dely*dely + delz*delz;

        // sqrt(2.5)=1.58 should be ok for both bulk graphene and h-BN
        if (rsq < 2.5) {

          if(comm->me==0 && intra_numneigh[i+1]==3){
            cout << "atom " << i+1 << " has more than three intra-layer neighbors" << endl;
            error->all(FLERR,"all atoms must have 1..3 intra-layer neighbors (cutoff^2=2.5)");
          }
          intra_neigh[i+1][intra_numneigh[i+1]]=j+1;
          intra_numneigh[i+1]++;

        }
      } //end of jj loop
    }//end of ii loop
     MPI_Barrier(world);
    if (comm->me==0) cout << "  4) neighbors list built" << endl;

    // we check that the number of intra-neighbors is 1..3
    imin=3;imax=0;
    for (i = 0; i < natoms; i++) {
      itype = tt[i];
      if( itype==0 || map[itype]==-1 ) continue; // atoms not belonging to this pair_potential
      if(comm->me==0 && intra_numneigh[i+1]<3) cout << "atom " << i+1 << " has " << intra_numneigh[i+1] << " intra-neighbors" << endl;
      if(intra_numneigh[i+1]<imin) imin=intra_numneigh[i+1];
      if(intra_numneigh[i+1]>imax) imax=intra_numneigh[i+1];
    }
    if(comm->me == 0) cout << "ILP: min-max " << imin << "-" << imax << " intra-neighbors found." << endl;
    if(imin<1 || imax>3)
      error->all(FLERR,"all atoms must have 1..3 intra-layer neighbors [cutoff=sqrt(2.5)]");
    if(comm->me==0) cout << "ILP: init done." << endl;

    // freeing temporary arrays
    free(xx);free(yy);free(zz);free(tt);free(mm);

    once=0;
    MPI_Barrier(world);

  }
  //END ONCE


  //BEGIN update normal vectors
  // TODO: allocate with memory->create(bx,ntotal,2,"pair:ilp");
  if(ntotal>ntotal0){ //realloc if ntotal grows upon reneighboring
    if( (  nx=(double *)realloc(nx,ntotal*sizeof(double ))) == NULL ) { error->all(FLERR,"allocation error."); }
    if( (  ny=(double *)realloc(ny,ntotal*sizeof(double ))) == NULL ) { error->all(FLERR,"allocation error."); }
    if( (  nz=(double *)realloc(nz,ntotal*sizeof(double ))) == NULL ) { error->all(FLERR,"allocation error."); }
    if( (nnorm=(double*)realloc(nnorm,ntotal*sizeof(double ))) == NULL ) { error->all(FLERR,"allocation error."); }
    if( ( bx=(double **)realloc(bx,ntotal*sizeof(double*))) == NULL ) { error->all(FLERR,"allocation error."); }
    if( ( by=(double **)realloc(by,ntotal*sizeof(double*))) == NULL ) { error->all(FLERR,"allocation error."); }
    if( ( bz=(double **)realloc(bz,ntotal*sizeof(double*))) == NULL ) { error->all(FLERR,"allocation error."); }
    for(i=ntotal0;i<ntotal;i++){ // https://stackoverflow.com/questions/49014310
      if( (   bx[i]=(double *)malloc(2*sizeof(double))) == NULL ) { error->all(FLERR,"allocation error."); }
      if( (   by[i]=(double *)malloc(2*sizeof(double))) == NULL ) { error->all(FLERR,"allocation error."); }
      if( (   bz[i]=(double *)malloc(2*sizeof(double))) == NULL ) { error->all(FLERR,"allocation error."); }
    }
    ntotal0=ntotal;
  }

  for(i=0;i<ntotal;i++){ //i=0..ntotal since we need normals also for ghost atoms

    itype=type[i];
    if( itype==0 || map[itype] == -1 ) continue; // atoms not belonging to this pair_potential

    ll=atom->tag[i];
    if(intra_numneigh[ll]==3){
      // a,b,c are the three local neighbors of i, mapped from the global (static) neighbors
      a=atom->map( intra_neigh[ll][0] );
      b=atom->map( intra_neigh[ll][1] );
      c=atom->map( intra_neigh[ll][2] );
    }
    else if(intra_numneigh[ll]==2){
      a=atom->map( intra_neigh[ll][0] );
      b=atom->map( intra_neigh[ll][1] );
      c=i;
    }
    else{ // only one intra-neighbor: it must have three intra-neighbors
      kk=intra_neigh[ll][0];
      if(intra_numneigh[kk]!=3) error->all(FLERR,"cannot determine normal vector of atom chain.");
      a=atom->map( intra_neigh[kk][0] );
      b=atom->map( intra_neigh[kk][1] );
      c=atom->map( intra_neigh[kk][2] );
    }

    if(a==-1 || b==-1 || c==-1) continue; //edge neighbors may not belong to this proc

    // the two bonds are a-b and a-c
    bx[i][0]=x[a][0]-x[b][0]; bx[i][0] -= round(bx[i][0]/PBCx)*PBCx; //TODO TEST wrapping needed?
    by[i][0]=x[a][1]-x[b][1]; by[i][0] -= round(by[i][0]/PBCy)*PBCy;
    bz[i][0]=x[a][2]-x[b][2]; bz[i][0] -= round(bz[i][0]/PBCz)*PBCz;
    bx[i][1]=x[a][0]-x[c][0]; bx[i][1] -= round(bx[i][1]/PBCx)*PBCx;
    by[i][1]=x[a][1]-x[c][1]; by[i][1] -= round(by[i][1]/PBCy)*PBCy;
    bz[i][1]=x[a][2]-x[c][2]; bz[i][1] -= round(bz[i][1]/PBCz)*PBCz;

    //cross product between the bonds
    nx[i]=(by[i][0]*bz[i][1]-bz[i][0]*by[i][1]);
    ny[i]=(bz[i][0]*bx[i][1]-bx[i][0]*bz[i][1]);
    nz[i]=(bx[i][0]*by[i][1]-by[i][0]*bx[i][1]);
    nnorm[i]=sqrt( nx[i]*nx[i] + ny[i]*ny[i] + nz[i]*nz[i] );
    nx[i]/=nnorm[i]; ny[i]/=nnorm[i]; nz[i]/=nnorm[i];

  }
  //END update normal vectors


  //////////////////////////////////////////////
  // loop over full neighbor list of my atoms
  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j][0]; //delx -= round(delx/PBCx)*PBCx;
      dely = ytmp - x[j][1]; //dely -= round(dely/PBCy)*PBCy;
      delz = ztmp - x[j][2]; //delz -= round(delz/PBCz)*PBCz;

      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq > cutsq[itype][jtype]) continue;

      r  = sqrt(rsq);
      r3 = r*rsq;
      r4 = rsq*rsq;
      r6 = r3*r3;

      iparam_ij = elem2param[map[itype]][map[jtype]];
      Param& p = params[iparam_ij];

      aasimp= exp(-ILP_d*( (r/(p.sR*p.reff)) - 1.));
      bsimp= r4*(   TAP_A*r3 -  TAP_B*rsq +  TAP_C*r -  TAP_D ) + 1.; //Tapering function
      csimp= rsq*( dTAP_A*r3 - dTAP_B*rsq + dTAP_C*r - dTAP_D );      //Tapering derivative (already divided by r)
      dsimp= p.C6/r6;

      esimp= csimp/(1.+aasimp)*dsimp;                                      //dispersive term 1
      esimp+= bsimp*dsimp*( 1./pow(1.+aasimp,2)*aasimp*ILP_d/(p.sR*p.reff)  //dispersive term 2
                           -1./(1.+aasimp)*6./r
                          )/r;
      tmp=esimp*delx;
      f[i][0] += tmp;
      f[j][0] -= tmp;
      tmp=esimp*dely;
      f[i][1] += tmp;
      f[j][1] -= tmp;
      tmp=esimp*delz;
      f[i][2] += tmp;
      f[j][2] -= tmp;

      ndotr_ij=(nx[i]*delx + ny[i]*dely + nz[i]*delz);
      rho2_ij = rsq - ndotr_ij*ndotr_ij;
      asimp= exp(p.alpha*(1.-r/p.beta) );
      dsimp_ij= exp( -rho2_ij/p.gamma2 );
      esimp= bsimp*asimp*p.C*dsimp_ij*2./p.gamma2;                         //repulsive term 3  (ij)

      tmp=esimp*(delx-ndotr_ij*nx[i]);
      f[i][0] += tmp;
      f[j][0] -= tmp;
      tmp=esimp*(dely-ndotr_ij*ny[i]);
      f[i][1] += tmp;
      f[j][1] -= tmp;
      tmp=esimp*(delz-ndotr_ij*nz[i]);
      f[i][2] += tmp;
      f[j][2] -= tmp;


      // evaluate and apply the forces on i neighbors due to its normal vector variation
      ll=atom->tag[i];
      if(intra_numneigh[ll]==3){
        // a,b,c are the three local neighbors of i, mapped from the global (static) neighbors
        a=atom->map( intra_neigh[ll][0] );
        b=atom->map( intra_neigh[ll][1] );
        c=atom->map( intra_neigh[ll][2] );
      }
      else if(intra_numneigh[ll]==2){
        a=atom->map( intra_neigh[ll][0] );
        b=atom->map( intra_neigh[ll][1] );
        c=i;
      }
      else{ // only one intra-neighbor: it must have three intra-neighbors
        kk=intra_neigh[ll][0];
        a=atom->map( intra_neigh[kk][0] );
        b=atom->map( intra_neigh[kk][1] );
        c=atom->map( intra_neigh[kk][2] );
      }

      if(a==-1 || b==-1 || c==-1) error->all(FLERR,"negative id mapping, try to increase the skin");

      temp_b = esimp*(-ndotr_ij)*( dely*bz[i][1] - delz*by[i][1] - ndotr_ij*( ny[i]*bz[i][1] - nz[i]*by[i][1]) )/nnorm[i];
      temp_c = esimp*(-ndotr_ij)*(-dely*bz[i][0] + delz*by[i][0] - ndotr_ij*(-ny[i]*bz[i][0] + nz[i]*by[i][0]) )/nnorm[i];
      f[a][0]-=( temp_b + temp_c );
      f[b][0]+=temp_b;
      f[c][0]+=temp_c;

      temp_b = esimp*(-ndotr_ij)*(-delx*bz[i][1] + delz*bx[i][1] - ndotr_ij*(-nx[i]*bz[i][1] + nz[i]*bx[i][1]) )/nnorm[i];
      temp_c = esimp*(-ndotr_ij)*( delx*bz[i][0] - delz*bx[i][0] - ndotr_ij*( nx[i]*bz[i][0] - nz[i]*bx[i][0]) )/nnorm[i];
      f[a][1]-=( temp_b + temp_c );
      f[b][1]+=temp_b;
      f[c][1]+=temp_c;

      temp_b = esimp*(-ndotr_ij)*( delx*by[i][1] - dely*bx[i][1] - ndotr_ij*( nx[i]*by[i][1] - ny[i]*bx[i][1]) )/nnorm[i];
      temp_c = esimp*(-ndotr_ij)*(-delx*by[i][0] + dely*bx[i][0] - ndotr_ij*(-nx[i]*by[i][0] + ny[i]*bx[i][0]) )/nnorm[i];
      f[a][2]-=( temp_b + temp_c );
      f[b][2]+=temp_b;
      f[c][2]+=temp_c;


      //JI
      ndotr_ji=(nx[j]*delx + ny[j]*dely + nz[j]*delz);
      rho2_ji = rsq - ndotr_ji*ndotr_ji;
      dsimp_ji= exp( -rho2_ji/p.gamma2 );
      esimp= bsimp*asimp*p.C*dsimp_ji*2./p.gamma2;                         //repulsive term 3  (ji)

      tmp=esimp*(delx-ndotr_ji*nx[i]);
      f[i][0] += tmp;
      f[j][0] -= tmp;
      tmp=esimp*(dely-ndotr_ji*ny[i]);
      f[i][1] += tmp;
      f[j][1] -= tmp;
      tmp=esimp*(delz-ndotr_ji*nz[i]);
      f[i][2] += tmp;
      f[j][2] -= tmp;

      // evaluate and apply the forces on j neighbors due to its normal vector variation
      ll=atom->tag[j];
      if(intra_numneigh[ll]==3){
        // a,b,c are the three local neighbors of i, mapped from the global (static) neighbors
        a=atom->map( intra_neigh[ll][0] );
        b=atom->map( intra_neigh[ll][1] );
        c=atom->map( intra_neigh[ll][2] );
      }
      else if(intra_numneigh[ll]==2){
        a=atom->map( intra_neigh[ll][0] );
        b=atom->map( intra_neigh[ll][1] );
        c=i;
      }
      else{ // only one intra-neighbor: it must have three intra-neighbors
        kk=intra_neigh[ll][0];
        a=atom->map( intra_neigh[kk][0] );
        b=atom->map( intra_neigh[kk][1] );
        c=atom->map( intra_neigh[kk][2] );
      }
      if(a==-1 || b==-1 || c==-1) error->all(FLERR,"negative id mapping, try to increase the skin");

      temp_b = esimp*(-ndotr_ji)*( dely*bz[j][1] - delz*by[j][1] - ndotr_ji*( ny[j]*bz[j][1] - nz[j]*by[j][1]) )/nnorm[j];
      temp_c = esimp*(-ndotr_ji)*(-dely*bz[j][0] + delz*by[j][0] - ndotr_ji*(-ny[j]*bz[j][0] + nz[j]*by[j][0]) )/nnorm[j];
      f[a][0]-=( temp_b + temp_c );
      f[b][0]+=temp_b;
      f[c][0]+=temp_c;

      temp_b = esimp*(-ndotr_ji)*(-delx*bz[j][1] + delz*bx[j][1] - ndotr_ji*(-nx[j]*bz[j][1] + nz[j]*bx[j][1]) )/nnorm[j];
      temp_c = esimp*(-ndotr_ji)*( delx*bz[j][0] - delz*bx[j][0] - ndotr_ji*( nx[j]*bz[j][0] - nz[j]*bx[j][0]) )/nnorm[j];
      f[a][1]-=( temp_b + temp_c );
      f[b][1]+=temp_b;
      f[c][1]+=temp_c;

      temp_b = esimp*(-ndotr_ji)*( delx*by[j][1] - dely*bx[j][1] - ndotr_ji*( nx[j]*by[j][1] - ny[j]*bx[j][1]) )/nnorm[j];
      temp_c = esimp*(-ndotr_ji)*(-delx*by[j][0] + dely*bx[j][0] - ndotr_ji*(-nx[j]*by[j][0] + ny[j]*bx[j][0]) )/nnorm[j];
      f[a][2]-=( temp_b + temp_c );
      f[b][2]+=temp_b;
      f[c][2]+=temp_c;


      // RADIAL FORCES
      dsimp= p.epsilon + p.C*( dsimp_ij + dsimp_ji );
      esimp= -csimp*asimp*dsimp;                       //repulsive term 1
      esimp+= bsimp*asimp*dsimp*p.alpha/p.beta/r;      //repulsive term 2

      f[i][0] += delx*esimp;
      f[i][1] += dely*esimp;
      f[i][2] += delz*esimp;
      f[j][0] -= delx*esimp;
      f[j][1] -= dely*esimp;
      f[j][2] -= delz*esimp;


      if (eflag) {
        //DISPERSIVE TERM
        evdwl = bsimp*(-p.C6/(1.+aasimp)/r6);
        //REPULSIVE TERM
        evdwl+= bsimp*exp(p.alpha*(1.-r/p.beta))
		*( p.epsilon + p.C*(
				    exp(-rho2_ij/p.gamma2)
				  + exp(-rho2_ji/p.gamma2)
				    )
		);
      }
      if (evflag){
        ev_tally(i,j,nlocal,force->newton_pair,evdwl,0.,0.,delx,dely,delz);
      }


    }//END LOOP ON jj
  }//END LOOP ON ii

  if (vflag_fdotr) virial_fdotr_compute();      //NOTE: virial not implemented!

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairILP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(offset,n+1,n+1,"pair:offset");
  map = new int[atom->ntypes+1];

}


/* ----------------------------------------------------------------------
   global settings
   reads the input script line with arguments you define
------------------------------------------------------------------------- */

void PairILP::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  if (strcmp(force->pair_style,"hybrid/overlay")!=0)
    error->all(FLERR,"ERROR: requires hybrid/overlay pair_style");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairILP::coeff(int narg, char **arg)
{
  int i,j,n;

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }


  read_file(arg[2]);

  double cut_one = cut_global;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   perform initialization for one i,j type pair (cutoff)
------------------------------------------------------------------------- */

double PairILP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  offset[j][i] = offset[i][j] = 0.0;

  return cut[i][j];
}


/* ----------------------------------------------------------------------
   read potential parameters file
------------------------------------------------------------------------- */

void PairILP::read_file(char *filename)
{
  int params_per_line = 13;
  char **words = new char*[params_per_line+1];
  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(filename);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open ILP potential file %s",filename);
      error->one(FLERR,str);
    }
  }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  int i,j,n,m,nwords,ielement,jelement;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Insufficient format in ILP potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement = 1st args
    // if these 2 args are in element list, then parse this line
    // else skip to next line (continue)

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),"pair:params");
    }

    //ILP.par lines contains: spec1, spec2, Reff(kcal/mol), C6(meV*Ang^6), alpha, beta(Ang), gamma(Ang), epsilon(kcal/mol), lambda(1/Ang), csi(meV), eta(meV), C(kcal/mol), sR
    params[nparams].ielement  = ielement;
    params[nparams].jelement  = jelement;
    params[nparams].reff      = atof(words[2]);
    params[nparams].C6        = atof(words[3]);
    params[nparams].alpha     = atof(words[4]);
    params[nparams].beta      = atof(words[5]);
    params[nparams].gamma     = atof(words[6]);
    params[nparams].epsilon   = atof(words[7]);
    params[nparams].lambda    = atof(words[8]);
    params[nparams].csi       = atof(words[9]);
    params[nparams].eta       = atof(words[10]);
    params[nparams].C         = atof(words[11]);
    params[nparams].sR        = atof(words[12]);

    // here rescale parameters to eV and Angstrom
    double kcalmol2eV = 43.3634420165735e-3;

    params[nparams].C6       *= kcalmol2eV;
    params[nparams].epsilon  *= kcalmol2eV;
    params[nparams].C        *= kcalmol2eV;

    // precompute some quantities
    params[nparams].gamma2    = pow(params[nparams].gamma,2);
    params[nparams].oolambda3 = 1./pow(params[nparams].lambda,3);

    nparams++;
    //if(nparams >= pow(atom->ntypes,3)) break;
  }

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,"pair:elem2param");
  for (i = 0; i < nelements; i++) {
    for (j = 0; j < nelements; j++) {
      n = -1;
      for (m = 0; m < nparams; m++) {
        if (i == params[m].ielement && j == params[m].jelement) {
          if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
          n = m;
        }
      }
      if (n < 0) error->all(FLERR,"Potential file is missing an entry");
      elem2param[i][j] = n;
    }
  }

  delete [] words;

}

/* ---------------------------------------------------------------------- */
