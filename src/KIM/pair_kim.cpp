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
   Contributing authors: Valeriu Smirichinski, Ryan Elliott, 
                         Ellad Tadmor (U Minn)
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

#include "KIM_API.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairKIM::PairKIM(LAMMPS *lmp) : Pair(lmp)
{
  nelements = 0;
  elements = NULL;

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  no_virial_fdotr_compute = 1;

  virialGlobal_ind = -1;
  particleVirial_ind = -1;
  process_dEdr_ind = -1;

  modelname = NULL;
  testname = NULL;

  maxall = 0;
  kimtype = NULL;
  fcopy = NULL;
  onebuf = NULL;

  pointsto = 0;
  memory->create(Rij,3*KIM_API_MAX_NEIGHBORS,"pair:Rij");

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

  memory->destroy(kimtype);
  memory->destroy(fcopy);
  memory->destroy(onebuf);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
  }

  memory->destroy(Rij);

  kim_free();
}

/* ---------------------------------------------------------------------- */

void PairKIM::compute(int eflag , int vflag)
{
  int kimerr;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = eflag_global = eflag_atom = vflag_global = vflag_atom = 0;

  // grow kimtype array if necessary
  // needs to be atom->nmax in length

  if (atom->nmax > maxall) {
    memory->destroy(kimtype);
    memory->destroy(fcopy);
    maxall = atom->nmax;
    memory->create(kimtype,maxall,"pair:kimtype");
    memory->create(fcopy,maxall,3,"pair:fcopy");
  }

  // kimtype = KIM atom type for each LAMMPS atom
  // set ielement to valid 0 if map[] stores an un-used -1

  int *type = atom->type;
  int nall = atom->nlocal + atom->nghost;
  int ielement;

  for (int i = 0; i < nall; i++) {
    ielement = map[type[i]];
    ielement = MAX(ielement,0);
    kimtype[i] = atypeMapKIM[2*ielement+1];
  }

  // pass current atom pointers to KIM

  set_volatiles();

  // set callback for virial tallying if necessary

  int process_dEdr_flag = 0; 
  if (process_dEdr_ind >= 0 && (vflag_atom || vflag_global))
    process_dEdr_flag = 1;

  pkim->setm_compute_by_index(&kimerr,12, particleEnergy_ind,
				 eflag_atom,1,
				 particleVirial_ind,vflag_atom,1,
				 virialGlobal_ind,vflag_global,1,
				 process_dEdr_ind,process_dEdr_flag,1);
  kim_error(__LINE__,"setm_compute_by_index",kimerr);

  // KIM initialization

  init2zero(pkim,&kimerr);
  kim_error(__LINE__,"PairKIM::init2zero(..)",kimerr);

  // compute energy, forces, virial via KIM model

  pkim->model_compute(&kimerr);
  kim_error(__LINE__,"PairKIM::pkim->model_compute() error",kimerr);

  // reverse comm for ghost atoms
  // required even when newton flag is off

  comm->reverse_comm_pair(this);

  // sum fcopy to f if running in hybrid mode

  if (hybrid) {
    double **f = atom->f;
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) {
      f[i][0] += fcopy[i][0];
      f[i][1] += fcopy[i][1];
      f[i][2] += fcopy[i][2];
    }
  }
  
  // flip sign of virial

  if (vflag_global && !pkim->virial_need2add)
    for (int i = 0; i < 6; i++) virial[i] = -1.0*virial[i];
  if (vflag_atom && !pkim->particleVirial_need2add)
    for (int i = 0; i < atom->nlocal; i++) 
      for (int j = 0; j < 6; j++) vatom[i][j] = -1.0*vatom[i][j];
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
  const char *test = "test_LAMMPS";
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

  kim_init();
  if (!pkim->model_init())
    kim_error(__LINE__,"KIM API:model_init() failed",-1);

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
    error->all(FLERR,"Pair style kim requires newton pair off");

  // setup onebuf for neighbors of one atom if needed

  molecular = atom->molecular;
  if (molecular) {
    memory->destroy(onebuf);
    memory->create(onebuf,neighbor->oneatom,"pair:onebuf");
  }

  // hybrid = 1 if running with pair hybrid

  hybrid = 0;
  if (force->pair_match("hybrid",0)) hybrid = 1;

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
  double *fp;
  if (hybrid) fp = &(fcopy[0][0]);
  else fp = &(atom->f[0][0]);
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
  double *fp;
  if (hybrid) fp = &(fcopy[0][0]);
  else fp = &(atom->f[0][0]);
 
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
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double PairKIM::memory_usage()
{
  double bytes = maxall * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   KIM-specific interface
------------------------------------------------------------------------- */

void PairKIM::kim_error(int ln, const char* msg, int errcode)
{
  if (errcode == 1) return;
  KIM_API_model::report_error(ln,(char *) __FILE__,
			      (char *) msg,errcode);
  error->all(__FILE__,ln,"Internal KIM error");
}

/* ---------------------------------------------------------------------- */

int PairKIM::get_neigh(void **kimmdl,int *mode,int *request,
		       int *atom, int *numnei, int **nei1atom, double **pRij)
{
  double * x;
   
  KIM_API_model *pkim = (KIM_API_model *) *kimmdl;

  // get neighObj from KIM API obj
  // this block will be changed in case of hybride/composite KIM potential

  int kimerr;
  PairKIM *self = (PairKIM *) pkim->get_test_buffer(&kimerr);
  
  *pRij = &(self->Rij[0]);
  
  NeighList * neiobj = (NeighList * ) (*pkim)[self->neighObject_ind].data;
  
  int * pnAtoms = (int *)(*pkim)[self->numberOfParticles_ind].data;
  if(pkim->support_Rij) x = (double *)(*pkim)[self->coordinates_ind].data;
  int nAtoms = *pnAtoms;
  int j,jj,inum,*ilist,*numneigh,**firstneigh;
  inum = neiobj->inum;             //# of I atoms neighbors are stored for
  ilist =neiobj->ilist ;          //local indices of I atoms
  numneigh =neiobj-> numneigh ;   // # of J neighbors for each I atom
  firstneigh =neiobj->firstneigh ;// ptr to 1st J int value of each I atom

  if (*mode==0){ //iterator mode
    if (*request== 0){ //restart iterator
      self->pointsto =0;
      *numnei=0;
      return KIM_STATUS_NEIGH_ITER_INIT_OK; //succsesful restart
    } else if (*request==1) { //increment iterator
      if (self->pointsto > inum || inum <0){
	self->error->one(FLERR,"KIM neighbor iterator exceeded range");
      } else if (self->pointsto == inum) {
	self->pointsto = 0;
	*numnei = 0;
	return KIM_STATUS_NEIGH_ITER_PAST_END; //reached end by iterator
      } else {
	*atom = ilist[self->pointsto];
	*numnei = numneigh[*atom];

	// strip off neighbor mask for molecular systems

	if (self->molecular) {
	  int n = *numnei;
	  int *ptr = firstneigh[*atom];
	  int *onebuf = self->onebuf;
	  for (int i = 0; i < n; i++)
	    onebuf[i] = *(ptr++) & NEIGHMASK;
	  *nei1atom = onebuf;
	} else *nei1atom = firstneigh[*atom];

	self->pointsto++;
	if (*numnei > KIM_API_MAX_NEIGHBORS) 
	  return KIM_STATUS_NEIGH_TOO_MANY_NEIGHBORS;
	if (pkim->support_Rij) {
	  for (jj=0; jj < *numnei; jj++) {
	    int i = *atom;
	    j = (*nei1atom)[jj];
	    self->Rij[jj*3 +0] = -x[i*3+0] + x[j*3+0];
	    self->Rij[jj*3 +1] = -x[i*3+1] + x[j*3+1];
	    self->Rij[jj*3 +2] = -x[i*3+2] + x[j*3+2];
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
    *atom = *request;
    *numnei = numneigh[*atom];
    
    // strip off neighbor mask for molecular systems

    if (self->molecular) {
      int n = *numnei;
      int *ptr = firstneigh[*atom];
      int *onebuf = self->onebuf;
      for (int i = 0; i < n; i++)
	onebuf[i] = *(ptr++) & NEIGHMASK;
      *nei1atom = onebuf;
    } else *nei1atom = firstneigh[*atom];

    if (*numnei > KIM_API_MAX_NEIGHBORS) 
      return KIM_STATUS_NEIGH_TOO_MANY_NEIGHBORS;
    if (pkim->support_Rij){
      for(int jj=0; jj < *numnei; jj++){
	int i = *atom;
	int j = (*nei1atom)[jj];
	self->Rij[jj*3 +0]=-x[i*3+0] + x[j*3+0];
	self->Rij[jj*3 +1]=-x[i*3+1] + x[j*3+1];
	self->Rij[jj*3 +2]=-x[i*3+2] + x[j*3+2];   
      }
    }
    return KIM_STATUS_OK;//successful end
  } else return KIM_STATUS_NEIGH_INVALID_MODE;//invalid mode

  return -16;//should not get here: unspecified error
}

/* ---------------------------------------------------------------------- */

void PairKIM::kim_free()
{
  int kimerr;
  pkim->model_destroy(&kimerr);
  pkim->free();
  delete pkim;
  delete [] atypeMapKIM;
}

/* ---------------------------------------------------------------------- */

void PairKIM::kim_init()
{
  support_atypes = true;
  strcpy(modelfile,"");
  strcpy(testfile,"");
  char *kim_dir = getenv("KIM_DIR");
  char *kim_models_dir = getenv("KIM_MODELS_DIR");
  char *current_path = getenv("PWD");
  if (kim_dir == NULL)
    error->all(FLERR,"KIM_DIR environment variable is unset");
  if (current_path == NULL)
    error->all(FLERR,"PWD environment variable is unset");

  if (kim_models_dir == NULL) {
    strcat(modelfile,kim_dir);strcat(modelfile,"MODELs/");
    strcat(modelfile,modelname); strcat(modelfile,"/");
    strcat(modelfile,modelname); strcat(modelfile,".kim");
  } else {
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
    this->kim_error(__LINE__,
		    "LAMMPS unit_style must be metal to "
		    "work with KIM models",-12);
  
  // create temporary kim file

  test_descriptor_string = new char[80*60];
  for (int i=50;i<80*60;i++) test_descriptor_string[i]='\0';
  
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
  sprintf(&test_descriptor_string[i_s],				\
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
  sprintf(&test_descriptor_string[i_s],"NEIGH_PURE_H             flag\n\n");
  i_s=strlen(test_descriptor_string);
  sprintf(&test_descriptor_string[i_s],"NEIGH_PURE_F             flag\n\n");
  i_s=strlen(test_descriptor_string);
  sprintf(&test_descriptor_string[i_s],"NEIGH_RVEC_F             flag\n\n");
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
	  "numberOfParticles           integer      none       []\n\n");
  i_s=strlen(test_descriptor_string);
  sprintf(&test_descriptor_string[i_s],
	  "numberContributingParticles    integer      none    []\n\n");
  
  if (support_atypes){
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "numberParticleTypes         integer      none []\n\n");
    i_s=strlen(test_descriptor_string);
    sprintf(&test_descriptor_string[i_s],
	    "particleTypes               integer      none [numberOfParticles]\n\n");
  }
  i_s=strlen(test_descriptor_string);
  sprintf(&test_descriptor_string[i_s],
	  "coordinates             real*8       length     "
	  "[numberOfParticles,3]\n\n");
  i_s=strlen(test_descriptor_string);
  sprintf(&test_descriptor_string[i_s],
	  "neighObject             pointer      none       []\n\n");
  i_s=strlen(test_descriptor_string);
  sprintf(&test_descriptor_string[i_s],
	  "get_neigh          method       none       []\n\n");
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
	  "[numberOfParticles,3]\n\n");
  i_s=strlen(test_descriptor_string);
  sprintf(&test_descriptor_string[i_s],
	  "particleEnergy           real*8       energy    "
	  "[numberOfParticles]\n\n");
  //virial and virial per atom will be added here
  //  i_s=strlen(test_descriptor_string);
  sprintf(&test_descriptor_string[i_s], 
	  "virial          real*8       energy    [6] \n\n\0");
  i_s=strlen(test_descriptor_string);
  sprintf(&test_descriptor_string[i_s],
	  "process_dEdr        method       none   [] \n\n");
  i_s=strlen(test_descriptor_string);
  sprintf(&test_descriptor_string[i_s],
	  "particleVirial         real*8       energy   "
	  "[numberOfParticles,6] \n\n\0");
  
  // kim file created now init and maptypes

  pkim = new KIM_API_model();
  if (!pkim->init_str_testname(test_descriptor_string,modelname))
    error->all(FLERR,"KIM initialization failed");
  
  delete [] test_descriptor_string;
    
  // get correct index of each variable in kim_api object
  
  pkim->getm_index(&kimerr,33,
			   "coordinates", &coordinates_ind,1,
			   "cutoff",  &cutoff_ind,1,
			   "particleEnergy",&particleEnergy_ind,1,
			   "energy", &energy_ind,1,
			   "forces",  &forces_ind,1,
			   "neighObject", &neighObject_ind,1,
			   "numberOfParticles",&numberOfParticles_ind,1,
			   "get_neigh", &get_neigh_ind,1,
			   "particleVirial", &particleVirial_ind,1,
			   "virial", &virialGlobal_ind,1,
			   "process_dEdr",&process_dEdr_ind,1 );
  this->kim_error(__LINE__,"getm_index",kimerr);
  
  if(support_atypes){
    numberParticleTypes_ind = pkim->get_index((char *) "numberParticleTypes",&kimerr);
    this->kim_error(__LINE__,"numberParticleTypes",kimerr);
    particleTypes_ind = pkim->get_index((char *) "particleTypes",&kimerr);
    this->kim_error(__LINE__,"particleTypes",kimerr);
  }
  
  set_statics();

  // setup mapping between LAMMPS and KIM atom types

  atypeMapKIM = new int[2*nelements];
  for(int i = 0; i < nelements; i++){
    atypeMapKIM[i*2 + 0] = i;
    int kimerr;
    atypeMapKIM[i*2 + 1] = pkim->get_partcl_type_code(elements[i],&kimerr);
    this->kim_error(__LINE__,"create_atypeMapKIM: symbol not found ",kimerr);
  }

  // check if Rij needed for get_neigh
  
  char * NBC_method =(char *) pkim->get_NBC_method(&kimerr);
  this->kim_error(__LINE__,"NBC method not set",kimerr);
  support_Rij=false;
  if (strcmp(NBC_method,"NEIGH_RVEC_F")==0) support_Rij=true;
  free((void*)NBC_method);
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_statics()
{ 
  if(support_atypes)
    pkim->set_data_by_index(numberParticleTypes_ind,1,(void *) &(atom->ntypes )) ;
  if ( process_dEdr_ind >= 0 )
    pkim->set_data((char *) "process_dEdr",1, (void*) &process_dEdr);
  localnall = (int) (atom->nghost + atom->nlocal);
  int kimerr;
  pkim->setm_data_by_index(&kimerr,20,
			      energy_ind,1,  (void *)&(this->eng_vdwl),1,
			      cutoff_ind,1,  (void *)&(this->cut_global),1,
			      get_neigh_ind,1,(void *) &get_neigh,1,
			      numberOfParticles_ind,1,  (void *) &localnall,1,
			      virialGlobal_ind,1,(void *)&(virial[0]),1);
  this->kim_error(__LINE__,"setm_data_by_index",kimerr);
  
  pkim->set_data((char *) "numberContributingParticles",1,(void *)&(atom->nlocal));

  pkim->set_test_buffer( (void *)this, &kimerr);

  this->kim_error(__LINE__,"set_test_buffer",kimerr);
  
}

/* ---------------------------------------------------------------------- */

void PairKIM::set_volatiles()
{
  localnall = (int) (atom->nghost + atom->nlocal);
  bigint nall = (bigint) localnall;
    
  pkim->set_data_by_index(coordinates_ind,nall*3,(void *) &(atom->x[0][0]) );
  if (hybrid)
    pkim->set_data_by_index(forces_ind,nall*3,(void *) &(fcopy[0][0]) ) ;
  else
    pkim->set_data_by_index(forces_ind,nall*3,(void *) &(atom->f[0][0]) ) ;
  pkim->set_data_by_index(particleEnergy_ind,nall,  (void *)this->eatom) ;
  pkim->set_data_by_index(neighObject_ind,1,  (void *)this->list) ;
  
  if (support_atypes)
    pkim->set_data_by_index(particleTypes_ind,nall,  (void *) kimtype) ;
  
  if(vflag_atom != 1) {
    (*pkim)[particleVirial_ind].flag->calculate = 0;
  } else {
    (*pkim)[particleVirial_ind].flag->calculate = 1;
    pkim->set_data_by_index(particleVirial_ind,(intptr_t)nall,
		       (void *)&vatom[0][0]) ;
  }
  if (vflag_global != 1) {
    (*pkim)[virialGlobal_ind].flag->calculate = 0;
  } else {
    (*pkim)[virialGlobal_ind].flag->calculate = 1;
  }
}

/* ---------------------------------------------------------------------- */

void PairKIM::init2zero(KIM_API_model *pkim, int *kimerr)
{
  // fetch pointer to this = PairKIM, so that callbacks can access class members
  // if error, have no ptr to PairKIM, so let KIM generate error message

  PairKIM *self = (PairKIM *) pkim->get_test_buffer(kimerr);
  if (*kimerr < 1)
    KIM_API_model::report_error(__LINE__,(char *) __FILE__,
				(char *) "Self ptr to PairKIM is invalid",
				*kimerr);

  int p1_ind=-1;bool process_d1=false;
  self->process_dE.virialGlobal_flag = 0;
  self->process_dE.particleVirial_flag =0;
  int ierGlobal,ierPerAtom;
  self->process_dE.virialGlobal = (double *) pkim->get_data((char *) "virial",&ierGlobal);
  self->process_dE.particleVirial = (double *) pkim->get_data((char *) "particleVirial",
					    &ierPerAtom);
  self->process_dE.numberOfParticles = (int *) pkim->get_data((char *) "numberOfParticles",kimerr);
  //halfNeighbors = !pkim->requiresFullNeighbors();
  p1_ind = pkim->get_index((char *) "process_dEdr");
  if (*kimerr !=KIM_STATUS_OK) return;
  if (ierGlobal == KIM_STATUS_OK && self->process_dE.virialGlobal != NULL) {
    self->process_dE.virialGlobal_flag = pkim->get_compute((char *) "virial");
    if (self->process_dE.virialGlobal_flag==1 && pkim->virial_need2add){
      self->process_dE.virialGlobal[0] =0.0;  self->process_dE.virialGlobal[1] =0.0;
      self->process_dE.virialGlobal[2] =0.0;  self->process_dE.virialGlobal[3] =0.0;
      self->process_dE.virialGlobal[4] =0.0;  self->process_dE.virialGlobal[5] =0.0;
      process_d1=true;
    }
  }
  
  if (ierPerAtom == KIM_STATUS_OK && self->process_dE.particleVirial != NULL) {
    self->process_dE.particleVirial_flag = pkim->get_compute((char *) "particleVirial");
    if (self->process_dE.particleVirial_flag==1 && pkim->particleVirial_need2add) {
      for (int i =0;i<(*self->process_dE.numberOfParticles)*6 ;i++) self->process_dE.particleVirial[i]=0.0;
      process_d1=true;
    }
  }
  if (p1_ind >= 0) {
    if (process_d1) pkim->set_compute((char *) "process_dEdr",1,kimerr);
    else pkim->set_compute((char *) "process_dEdr",0,kimerr);
  }
}

/* ---------------------------------------------------------------------- */

void PairKIM::process_dEdr(KIM_API_model **ppkim,
			    double *de,double *r,double ** pdx,
			    int *i,int *j,int *ier)
{
  KIM_API_model *pkim = *ppkim;
  PairKIM *self = (PairKIM *) pkim->get_test_buffer(ier);

  *ier=KIM_STATUS_FAIL;
  double vir[6],v;
  double *dx = *pdx;
  //v=(*de)/((*r)*(*r))/3.0;
  //v=(*de)/((*r)*(*r));
  
  v=-(*de)/(*r);
  
  if (self->process_dE.virialGlobal_flag ==1 && pkim->virial_need2add) {
    vir[0] = v * dx[0] * dx[0];
    vir[1] = v * dx[1] * dx[1];
    vir[2] = v * dx[2] * dx[2];
    vir[3] = v * dx[1] * dx[2];
    vir[4] = v * dx[0] * dx[2];
    vir[5] = v * dx[0] * dx[1];
    self->process_dE.virialGlobal[0] += vir[0];
    self->process_dE.virialGlobal[1] += vir[1];
    self->process_dE.virialGlobal[2] += vir[2];
    self->process_dE.virialGlobal[3] += vir[3];
    self->process_dE.virialGlobal[4] += vir[4];
    self->process_dE.virialGlobal[5] += vir[5];
  }

  if (self->process_dE.particleVirial_flag==1 && pkim->particleVirial_need2add ){
    vir[0] =0.5 * v * dx[0] * dx[0];
    vir[1] =0.5 * v * dx[1] * dx[1];
    vir[2] =0.5 * v * dx[2] * dx[2];
    vir[3] =0.5 * v * dx[1] * dx[2];
    vir[4] =0.5 * v * dx[0] * dx[2];
    vir[5] =0.5 * v * dx[0] * dx[1];
    self->process_dE.particleVirial[(*i)*6 + 0] += vir[0];
    self->process_dE.particleVirial[(*i)*6 + 1] += vir[1];
    self->process_dE.particleVirial[(*i)*6 + 2] += vir[2];
    self->process_dE.particleVirial[(*i)*6 + 3] += vir[3];
    self->process_dE.particleVirial[(*i)*6 + 4] += vir[4];
    self->process_dE.particleVirial[(*i)*6 + 5] += vir[5];
    
    self->process_dE.particleVirial[(*j)*6 + 0] += vir[0];
    self->process_dE.particleVirial[(*j)*6 + 1] += vir[1];
    self->process_dE.particleVirial[(*j)*6 + 2] += vir[2];
    self->process_dE.particleVirial[(*j)*6 + 3] += vir[3];
    self->process_dE.particleVirial[(*j)*6 + 4] += vir[4];
    self->process_dE.particleVirial[(*j)*6 + 5] += vir[5];
  }
  
  *ier = KIM_STATUS_OK;
}
