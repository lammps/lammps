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
   Contributing authors: Valeriu Smirichinski, Ryan S. Elliott,
                         Ellad B. Tadmor (U Minn)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_kim.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h" 
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "memory.h"
#include "error.h"

#include <iostream>
#include "KIMservice.h"

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */

PairKIM::PairKIM(LAMMPS *lmp) : Pair(lmp)
{
  nelements = 0;
  elements = NULL;

  virialGlobal_ind = -1;
  virialPerAtom_ind = -1;
  process_d1Edr_ind = -1;

  modelname = NULL;
  testname = NULL;

  // reverse comm even if newton off

  comm_reverse_off = 9;
}

/* ---------------------------------------------------------------------- */

PairKIM::~PairKIM()
{
  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;

  delete [] modelname;
  delete [] testname;

  my_free();

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
  }
}

/* ---------------------------------------------------------------------- */

void PairKIM::compute(int eflag , int vflag)
{
  if (eflag || vflag) {
    if (vflag) ev_setup(eflag,5);
    else ev_setup(eflag,vflag);
  }
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 
	 vflag_global = vflag_atom = 0;
  
  set_volatiles();

  long long nall = (long long) (atom->nghost + atom->nlocal);

  int energyPerAtom_flag =0, virialPerAtom_flag=0, 
    process_d1Edr_flag=0, virialGlobal_flag=0;
  int kimerr;

  if (eflag_atom==1)  energyPerAtom_flag=1;
  else                energyPerAtom_flag=0;

  if (vflag_atom==1)  virialPerAtom_flag=1;
  else                virialPerAtom_flag=0;

  if (vflag_global ==1) virialGlobal_flag =1;
  else                  virialGlobal_flag =0;

 if (process_d1Edr_ind >=0) if (vflag_atom==1 || vflag_global==1) {
     process_d1Edr_flag =1;
  }else {
      process_d1Edr_flag=0;
  }

  pkim->set_compute_byI_multiple(&kimerr,12, energyPerAtom_ind,
				 energyPerAtom_flag,1,
     virialPerAtom_ind,virialPerAtom_flag,1,
				 virialGlobal_ind,virialGlobal_flag,1,
				 process_d1Edr_ind,process_d1Edr_flag,1);
  err_treat(__LINE__,"set_compute_byI_multiple",kimerr);

  init2zero(pkim,&kimerr);
  err_treat(__LINE__,"PairKIM::init2zero(..)",kimerr);

  pkim->model_compute(&kimerr);

  err_treat(__LINE__,"PairKIM::pkim->model_compute() error",kimerr);

  static int runs=0;
  if (runs < 1) {
    runs++;
    cout<<"nghost = "<<atom->nghost<<", nlocals = "<<atom->nlocal<<endl;
    cout<<"NBC method: "<<(char *) pkim->get_NBC_method(&kimerr)<<endl;
  }
 
  comm->reverse_comm_pair(this);
  
  if(virialGlobal_flag==1 && !pkim->virialGlobal_need2add ){
      for(int i=0; i<6;i++) virial[i]=-1.0*virial[i];
  }
  if(virialPerAtom_flag==1 && !pkim->virialPerAtom_need2add ){
      for(int i=0; i<atom->nlocal;i++) 
         for(int j=0; j<6;j++) vatom[i][j]=-1.0*vatom[i][j];
  }
}
 
/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */

void PairKIM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairKIM::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  delete [] modelname;
  int n = strlen(arg[0]) + 1;
  modelname = new char[n];
  strcpy(modelname,arg[0]);

  delete [] testname;
  char *test = "test_LAMMPS";
  n = strlen(test) + 1;
  testname = new char[n];
  strcpy(testname,test);

  ev_setup(3,6);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairKIM::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (narg != 2 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to KIM elements
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
  for (i = 2; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-1] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-1] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // KIM initialization

  my_init();
  if (!pkim->model_init())
    err_treat(__LINE__,"KIM API:model_init() failed",-1);

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
	setflag[i][j] = 1;
	count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairKIM::init_style()
{
  if (force->newton_pair != 0)
    error->all(FLERR,"Pair style pair_KIM requires newton pair Off");

  // request half or full neighbor list depending on KIM model requirement

  int irequest = neighbor->request(this);
  if (pkim->requiresFullNeighbors()) {
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairKIM::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cut_global;
}

/* ---------------------------------------------------------------------- */

int PairKIM::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  double *fp = &(atom->f[0][0]);
  double *va=&(vatom[0][0]);

  m = 0;
  last = first + n;
  if (vflag_atom == 1) {
    for (i = first; i < last; i++) {
        buf[m++] = fp[3*i+0];
        buf[m++] = fp[3*i+1];
        buf[m++] = fp[3*i+2];

        buf[m++] = va[6*i+0];
        buf[m++] = va[6*i+1];
        buf[m++] = va[6*i+2];
        buf[m++] = va[6*i+3];
        buf[m++] = va[6*i+4];
        buf[m++] = va[6*i+5];
    }
    return 9;

  } else {
    for (i = first; i < last; i++){
      buf[m++] = fp[3*i+0];
      buf[m++] = fp[3*i+1];
      buf[m++] = fp[3*i+2];
    }
    return 3;
  }
}

/* ---------------------------------------------------------------------- */

void PairKIM::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  double *fp = &(atom->f[0][0]);
 
  double *va=&(vatom[0][0]);
  m = 0;
  if (vflag_atom == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      fp[3*j+0]+= buf[m++];
      fp[3*j+1]+= buf[m++];
      fp[3*j+2]+= buf[m++];

      va[j*6+0]+=buf[m++];
      va[j*6+1]+=buf[m++];
      va[j*6+2]+=buf[m++];
      va[j*6+3]+=buf[m++];
      va[j*6+4]+=buf[m++];
      va[j*6+5]+=buf[m++];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      fp[3*j+0]+= buf[m++];
      fp[3*j+1]+= buf[m++];
      fp[3*j+2]+= buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   KIM-specific interface
------------------------------------------------------------------------- */

int PairKIM::get_neigh(void **kimmdl,int *mode,int *request,
        int *atom, int *numnei, int **nei1atom, double **pRij)
{
  static int neighObject_ind=-1, numberOfAtoms_ind=-1, coordinates_ind=-1;
  static bool key1strun=true;
  static double Rij[3*KIM_API_MAX_NEIGHBORS];
   
  double * x;
  *pRij=&Rij[0];
   
  KIM_API_model *pkim = (KIM_API_model *) *kimmdl;
  //get neighObj from KIM API obj
  if (key1strun) {  //this block will be changed in case of hybride/composite KIM potential
    int kimerr;
    
    neighObject_ind = pkim->get_index("neighObject",&kimerr);
    err_treat(__LINE__,"get_neigh:get_index of : neighObject",kimerr);
    
    
    NeighList * neiobj = (NeighList * ) (*pkim)[neighObject_ind].data;
    cout<<"inum="<<neiobj->inum;
    cout<<", ghostflag ="<<neiobj->ghostflag<<" requires full neighbor? :"<<pkim->requiresFullNeighbors()<<endl;
    
    numberOfAtoms_ind= pkim->get_index("numberOfAtoms",&kimerr);
    err_treat(__LINE__,"get_neigh:get_index of : neighObject",kimerr);
    key1strun=false;
    
    coordinates_ind= pkim->get_index("coordinates",&kimerr);
    err_treat(__LINE__,"get_neigh:get_index of : coordinates",kimerr);
    key1strun=false;
  }
    
  NeighList * neiobj = (NeighList * ) (*pkim)[neighObject_ind].data;
  
  int * pnAtoms = (int *)(*pkim)[numberOfAtoms_ind].data;
  if(pkim->support_Rij) x = (double *)(*pkim)[coordinates_ind].data;
  int nAtoms = *pnAtoms;
  int j,jj,inum,*ilist,*numneigh,**firstneigh;
  inum = neiobj->inum;             //# of I atoms neighbors are stored for
  ilist =neiobj->ilist ;          //local indices of I atoms
  numneigh =neiobj-> numneigh ;   // # of J neighbors for each I atom
  firstneigh =neiobj->firstneigh ;// ptr to 1st J int value of each I atom

  static int pointsto=0;
  if (*mode==0){ //iterator mode
    if (*request== 0){ //restart iterator
      pointsto =0;
      *numnei=0;
      return KIM_STATUS_NEIGH_ITER_INIT_OK; //succsesful restart
    } else if(*request==1){//increment iterator
      if (pointsto > inum || inum <0){
	printf("pair_KIM: iterator is >= then numnei or < 0\n");
	cout<<" pointsto ="<<pointsto<<" inum="<<inum<<endl;
	exit (237);
      }else if(pointsto == inum) {
	pointsto ==0;
	*numnei=0;
	return KIM_STATUS_NEIGH_ITER_PAST_END; //reached end by iterator
      }else{
	*numnei = numneigh[pointsto] ;
	*atom   =ilist[pointsto];
	*nei1atom = firstneigh[pointsto];
	pointsto++;
	if (*numnei > KIM_API_MAX_NEIGHBORS) 
	  return KIM_STATUS_NEIGH_TOO_MANY_NEIGHBORS;
	if (pkim->support_Rij){
	  for( jj=0; jj < *numnei; jj++){
	    int i = *atom;
	    j = (*nei1atom)[jj];
	    Rij[jj*3 +0] = -x[i*3+0] + x[j*3+0];
	    Rij[jj*3 +1] = -x[i*3+1] + x[j*3+1];
	    Rij[jj*3 +2] = -x[i*3+2] + x[j*3+2];
	    
	  }
	}
	return KIM_STATUS_OK;//successful increment
      }
    }
  }else if (*mode ==1){//locator mode
    //...
    if (*request >= nAtoms || inum <0) return -1;
    if (*request >= inum) {
      *atom=*request;
      *numnei=0;
      return KIM_STATUS_OK;//successfull but no neighbors in the list
    }
    *atom=ilist[*request];
    *numnei=numneigh[*request];
    *nei1atom = firstneigh[*request];
    if(*numnei > KIM_API_MAX_NEIGHBORS) return KIM_STATUS_NEIGH_TOO_MANY_NEIGHBORS;
    if (pkim->support_Rij){
      for(int jj=0; jj < *numnei; jj++){
	int i = *atom;
	int j = (*nei1atom)[jj];
	Rij[jj*3 +0]=-x[i*3+0] + x[j*3+0];
	Rij[jj*3 +1]=-x[i*3+1] + x[j*3+1];
	Rij[jj*3 +2]=-x[i*3+2] + x[j*3+2];
        
      }
    }
    return KIM_STATUS_OK;//successful end
  }else{
    return KIM_STATUS_NEIGH_INVALID_MODE;//invalid mode
  }
  return -16;//should not get here: unspecified error
}

/* ---------------------------------------------------------------------- */

void PairKIM::my_free(){
    int kimerr;
    pkim->model_destroy(&kimerr);
  
    pkim->free();
    // cout<<(*pkim);
  
    delete  pkim;
    if (atypeMapKIM != NULL) delete [] atypeMapKIM;
    if (atomTypes!=NULL) delete [] atomTypes;
}

/* ---------------------------------------------------------------------- */

void PairKIM::my_init()
{
  support_atypes=true;
  //if (narg == 0) support_atypes=false;

  // the section has to be done once?*********************
  // get paths

  strcpy(modelfile,"");
  strcpy(testfile,"");
  char * kim_dir =getenv("KIM_DIR");
  char * kim_models_dir = getenv("KIM_MODELS_DIR");
  char * current_path = getenv("PWD");
  if (kim_dir == NULL)
    cout<<"pair_KIM: error in initialization: env   KIM_DIR is not set"<<endl;

    if (current_path == NULL)
      cout<<"pair_KIM: error in initialization: env   PWD is not supported on OS system"<<endl;

    if (kim_models_dir == NULL) {
      strcat(modelfile,kim_dir);strcat(modelfile,"MODELs/");
      strcat(modelfile,modelname); strcat(modelfile,"/");
      strcat(modelfile,modelname); strcat(modelfile,".kim");
    }else{
      strcat(modelfile,kim_models_dir);strcat(modelfile,modelname);
      strcat(modelfile,"/");
      strcat(modelfile,modelname);
      strcat(modelfile,".kim");
    }
    strcat(testfile,current_path); strcat(testfile,"/");
    strcat(testfile,testname);strcat(testfile,".kim");

    // now testfile and modelfile are set

    int kimerr;

    if (!(strcmp(update->unit_style,"metal")==0))
      err_treat(__LINE__,"LAMMPS unit_style must be metal to work with KIM models",-12);
  
    // create temporary kim file
    test_descriptor_string = new char[80*60];
    for(int i=50;i<80*60;i++) test_descriptor_string[i]='\0';

    // 1. write atomic header and atomic type spec

    int i_s;
    sprintf(test_descriptor_string,
	    "# This file is automatically generated from LAMMPS "
	    "pair_style PairKIM command\n");
    i_s=strlen(test_descriptor_string);
    
    sprintf(&test_descriptor_string[i_s],"TEST_NAME :=%s\n" ,testname);
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],"#Unit_Handling := flexible\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],"Unit_length := A\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],"Unit_energy := eV \n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],"Unit_charge := e\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],"Unit_temperature := K\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],"Unit_time := fs\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],"SUPPORTED_ATOM/PARTICLES_TYPES:\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s], 
	    "# Simbol/name           Type            code\n");
    int code=1;
    for (int i=0; i < nelements; i++){
      i_s=strlen(test_descriptor_string);
      sprintf(&test_descriptor_string[i_s],"%s\t\t\tspec\t\t%d\n",
	      elements[i],code++);
    }
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],"\n");
    
    // 2. conventions

    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s], \
	    "CONVENTIONS:\n# Name                  Type\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s], "ZeroBasedLists           flag\n\n");
    i_s=strlen(test_descriptor_string);

    // can use iterator or locator neighbor mode, unless hybrid

    int iterator_only = 0;
    if (force->pair_match("hybrid",0)) iterator_only = 1;
    if (iterator_only)
      sprintf(&test_descriptor_string[i_s],"Neigh_IterAccess         flag\n\n");
    else
      sprintf(&test_descriptor_string[i_s],"Neigh_BothAccess         flag\n\n");

    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],"NEIGH-PURE-H             flag\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],"NEIGH-PURE-F             flag\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],"NEIGH-RVEC-F             flag\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],"CLUSTER                  flag\n\n");
    
    // 3. input-output

    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],"MODEL_INPUT:\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "# Name                  Type         Unit       Shape\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "numberOfAtoms           integer      none       []\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "numberContributingAtoms    integer      none    []\n\n");
    
    if (support_atypes){
      i_s=strlen(test_descriptor_string);
      sprintf(&test_descriptor_string[i_s],
	      "numberAtomTypes         integer      none []\n\n");
      i_s=strlen(test_descriptor_string);
      sprintf(&test_descriptor_string[i_s],
	      "atomTypes               integer      none [numberOfAtoms]\n\n");
    }
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "coordinates             real*8       length     "
	    "[numberOfAtoms,3]\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "neighObject             pointer      none       []\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "get_half_neigh          method       none       []\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "get_full_neigh          method       none       []\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s], "MODEL_OUPUT:\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "# Name                  Type         Unit       Shape\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "compute                 method       none       []\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "destroy                 method       none       []\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "cutoff                  real*8       length     []\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "energy                  real*8       energy     []\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "forces                  real*8       force      "
	    "[numberOfAtoms,3]\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "energyPerAtom           real*8       energy    "
	    "[numberOfAtoms]\n\n");
    //virial and virial per atom will be added here
    //  i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s], 
	    "virialGlobal          real*8       pressure    [6] \n\n\0");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "process_d1Edr        method       none   [] \n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "virialPerAtom         real*8       pressure   "
	    "[numberOfAtoms,6] \n\n\0");
    
    //.kim file created now init and maptypes
    
    my_init2();
    
    create_atypeMapKIM();

    //check if Rij needed for get_neigh
    char * NBC_method =(char *) pkim->get_NBC_method(&kimerr);
    err_treat(__LINE__,"NBC method not set",kimerr);
    support_Rij=false;
    if (strcmp(NBC_method,"NEIGH-RVEC-F")==0) support_Rij=true;
    delete [] NBC_method;
}

/* ---------------------------------------------------------------------- */

void PairKIM::my_init2()
{
  pkim = new KIM_API_model();
  if(!pkim->init_str_testname(test_descriptor_string,modelname)){
    cout<<"pkim: error in initialization "<<endl;
    exit (357);
  }
  
  delete [] test_descriptor_string;
  
  int kimerr=0;
    
  // get correct index of each variable in kim_api object
  
  pkim->get_index_multiple(&kimerr,36,
			   "coordinates", &coordinates_ind,1,
			   "cutoff",  &cutoff_ind,1,
			   "energyPerAtom",&energyPerAtom_ind,1,
			   "energy", &energy_ind,1,
			   "forces",  &forces_ind,1,
			   "neighObject", &neighObject_ind,1,
			   "numberOfAtoms",&numberOfAtoms_ind,1,
			   "get_half_neigh", &get_half_neigh_ind,1,
			   "get_full_neigh",&get_full_neigh_ind,1,
			   "virialPerAtom", &virialPerAtom_ind,1,
			   "virialGlobal", &virialGlobal_ind,1,
			   "process_d1Edr",&process_d1Edr_ind,1 );
  err_treat(__LINE__,"get_index_multiple",kimerr);
  
  if(support_atypes){
    numberAtomTypes_ind = pkim->get_index("numberAtomTypes",&kimerr);
    err_treat(__LINE__,"numberAtomTypes",kimerr);
    atomTypes_ind = pkim->get_index("atomTypes",&kimerr);
    err_treat(__LINE__,"atomTypes",kimerr);
  }
  
  set_statics();
}

/* ---------------------------------------------------------------------- */

void PairKIM:: create_atypeMapKIM()
{
  atypeMapKIM = new int[2*nelements];
  for(int i = 0; i < nelements; i++){
    atypeMapKIM[i*2 + 0] = i;
    int kimerr;
    atypeMapKIM[i*2 + 1] = pkim->get_aTypeCode(elements[i],&kimerr);
    err_treat(__LINE__,"create_atypeMapKIM: symbol not found ",kimerr);
  }
  atomTypes=NULL;
}

/* ---------------------------------------------------------------------- */

void PairKIM::atypeMap2KIM(){
    static long long n=0;
    if (n!= (long long)localnall) {
        n=(long long)localnall;
        if (atomTypes != NULL) delete [] atomTypes;
        atomTypes = new int[n];
    }

  int itype;
  int *type = atom->type;

  for (int i= 0; i < (int) n; i++) {
    itype = map[type[i]];
    itype = MAX(itype,0);
    atomTypes[i] = atypeMapKIM[2*itype+1];
  }
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_statics(){ //set pointers that are fixed (never reallocated)

  if(support_atypes)
    pkim->set_data_byi(numberAtomTypes_ind,1,(void *) &(atom->ntypes )) ;
  if ( process_d1Edr_ind >= 0 )
    pkim->set_data("process_d1Edr",1, (void*) &process_d1Edr);
  localnall = (int) (atom->nghost + atom->nlocal);
  int kimerr;
  pkim->set_data_byI_multiple(&kimerr,24,
			      energy_ind,1,  (void *)&(this->eng_vdwl),1,
			      cutoff_ind,1,  (void *)&(this->cut_global),1,
			      get_half_neigh_ind,1,(void *) &get_neigh,1,
			      get_full_neigh_ind,1,(void *) &get_neigh,1,
			      numberOfAtoms_ind,1,  (void *) &localnall,1,
			      virialGlobal_ind,1,(void *)&(virial[0]),1);
  err_treat(__LINE__,"set_data_byI_multiple",kimerr);
  
  pkim->set_data("numberContributingAtoms",1,(void *)&(atom->nlocal));
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_volatiles()
{ //set pointers that might be reallocated
  localnall = (int) (atom->nghost + atom->nlocal);
  long long nall= (long long) localnall;
    
  pkim->set_data_byi(coordinates_ind,nall*3,(void *) &(atom->x[0][0]) );
  pkim->set_data_byi(forces_ind,nall*3,(void *) &(atom->f[0][0]) ) ;
  pkim->set_data_byi(energyPerAtom_ind,nall,  (void *)this->eatom) ;
  pkim->set_data_byi(neighObject_ind,1,  (void *)this->list) ;
  
  if(support_atypes) {
    atypeMap2KIM();
    pkim->set_data_byi(atomTypes_ind,nall,  (void *) atomTypes) ;
  }
  
  if(vflag_atom != 1) {
    (*pkim)[virialPerAtom_ind].flag->calculate = 0;
  } else {
    (*pkim)[virialPerAtom_ind].flag->calculate = 1;
    pkim->set_data_byi(virialPerAtom_ind,(intptr_t)nall,
		       (void *)&vatom[0][0]) ;
  }
  if (vflag_global != 1) {
    (*pkim)[virialGlobal_ind].flag->calculate = 0;
  } else {
    (*pkim)[virialGlobal_ind].flag->calculate = 1;
  }
}

/* ---------------------------------------------------------------------- */

namespace process_fij_4_pair_KIM
{
  double *virialGlobal;
  double *virialPerAtom;
  int virialGlobal_flag;
  int virialPerAtom_flag;
  int *numberOfAtoms;
  bool halfNeighbors;
  double *de, *r;
  double ** pdx;
  int *i, *j, *ier;
  KIM_API_model **ppkim;
}

/* ---------------------------------------------------------------------- */

void PairKIM:: init2zero(KIM_API_model *pkim,int *kimerr)
{
  using namespace process_fij_4_pair_KIM;
  int p1_ind=-1;bool process_d1=false;
  virialGlobal_flag = 0;
  virialPerAtom_flag =0;
  int ierGlobal,ierPerAtom;
  virialGlobal = (double *) pkim->get_data("virialGlobal",&ierGlobal);
  virialPerAtom = (double *) pkim->get_data("virialPerAtom",&ierPerAtom);
  numberOfAtoms = (int *) pkim->get_data("numberOfAtoms",kimerr);
  //halfNeighbors = !pkim->requiresFullNeighbors();
  p1_ind = pkim->get_index("process_d1Edr");
  if (*kimerr !=KIM_STATUS_OK) return;
  if (ierGlobal == KIM_STATUS_OK && virialGlobal != NULL) {
    virialGlobal_flag = pkim->isit_compute("virialGlobal");
    if (virialGlobal_flag==1 && pkim->virialGlobal_need2add){
      virialGlobal[0] =0.0;  virialGlobal[1] =0.0;  virialGlobal[2] =0.0;
      virialGlobal[3] =0.0;  virialGlobal[4] =0.0;  virialGlobal[5] =0.0;
      process_d1=true;
    }
  }
  
  if (ierPerAtom == KIM_STATUS_OK && virialPerAtom != NULL) {
    virialPerAtom_flag = pkim->isit_compute("virialPerAtom");
    if (virialPerAtom_flag==1 && pkim->virialPerAtom_need2add) {
      for (int i =0;i<(*numberOfAtoms)*6 ;i++) virialPerAtom[i]=0.0;
      process_d1=true;
    }
  }
  if (p1_ind >=0){
    if (process_d1) {
      pkim->set2_compute("process_d1Edr");
    } else {
      pkim->set2_donotcompute("process_d1Edr");
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairKIM::process_d1Edr(KIM_API_model **ppkim,
			    double *de,double *r,double ** pdx,
			    int *i,int *j,int *ier)
{
  using namespace process_fij_4_pair_KIM;
  
  KIM_API_model *pkim = *ppkim;
  *ier=KIM_STATUS_FAIL;
  double vir[6],v;
  double *dx = *pdx;
  //v=(*de)/((*r)*(*r))/3.0;
  //v=(*de)/((*r)*(*r));
  
  v=-(*de)/(*r);
  
  if (virialGlobal_flag ==1 && pkim->virialGlobal_need2add) {
    static int numruns1 =0;
    if (numruns1 <1) cout<<"*********** global virial computing from pair_KIM!!!*********"<<endl;numruns1++;
    vir[0] = v * dx[0] * dx[0];
    vir[1] = v * dx[1] * dx[1];
    vir[2] = v * dx[2] * dx[2];
    vir[3] = v * dx[1] * dx[2];
    vir[4] = v * dx[0] * dx[2];
    vir[5] = v * dx[0] * dx[1];
    virialGlobal[0] += vir[0];
    virialGlobal[1] += vir[1];
    virialGlobal[2] += vir[2];
    virialGlobal[3] += vir[3];
    virialGlobal[4] += vir[4];
    virialGlobal[5] += vir[5];
  }

  if (virialPerAtom_flag==1 && pkim->virialPerAtom_need2add ){
    static int numruns =0; 
    if (numruns <1) cout<<"*********** virial computing from pair_KIM!!!*********"<<endl;numruns++;
    vir[0] =0.5 * v * dx[0] * dx[0];
    vir[1] =0.5 * v * dx[1] * dx[1];
    vir[2] =0.5 * v * dx[2] * dx[2];
    vir[3] =0.5 * v * dx[1] * dx[2];
    vir[4] =0.5 * v * dx[0] * dx[2];
    vir[5] =0.5 * v * dx[0] * dx[1];
    virialPerAtom[(*i)*6 + 0] += vir[0];
    virialPerAtom[(*i)*6 + 1] += vir[1];
    virialPerAtom[(*i)*6 + 2] += vir[2];
    virialPerAtom[(*i)*6 + 3] += vir[3];
    virialPerAtom[(*i)*6 + 4] += vir[4];
    virialPerAtom[(*i)*6 + 5] += vir[5];
    
    virialPerAtom[(*j)*6 + 0] += vir[0];
    virialPerAtom[(*j)*6 + 1] += vir[1];
    virialPerAtom[(*j)*6 + 2] += vir[2];
    virialPerAtom[(*j)*6 + 3] += vir[3];
    virialPerAtom[(*j)*6 + 4] += vir[4];
    virialPerAtom[(*j)*6 + 5] += vir[5];
  }
  
  *ier = KIM_STATUS_OK;
}

/* ---------------------------------------------------------------------- */

 void PairKIM::err_treat(char* msg, int errcode)
{
  if (KIM_API_model::report_error(__LINE__,__FILE__,msg,errcode) < 1) exit(327);
}

/* ---------------------------------------------------------------------- */

void PairKIM::err_treat(int ln,char* msg, int errcode)
{
  if (KIM_API_model::report_error(ln,__FILE__,msg,errcode) < 1) exit(327);
}
