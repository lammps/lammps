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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fix_rx.h"
#include "atom.h"
#include "error.h"
#include "group.h"
#include "modify.h"
#include "force.h"
#include "memory.h"
#include "comm.h"
#include "update.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair_dpd_fdt_energy.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,HARMONIC};
enum{LUCY};

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

FixRX::FixRX(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6 || narg > 7) error->all(FLERR,"Illegal fix rx command");
  nevery = 1;

  nreactions = maxparam = 0;
  params = NULL;
  mol2param = NULL;
  pairDPDE = NULL;

  // Keep track of the argument list.
  int iarg = 3;

  // Read the kinetic file in arg[3].
  kineticsFile = arg[iarg++];

  // Determine the local temperature averaging method in arg[4].
  wtFlag = 0;
  localTempFlag = NONE;

  {
    char *word = arg[iarg++];
    if (strcmp(word,"none") == 0)
    {
      wtFlag = 0;
      localTempFlag = NONE;
    }
    else if (strcmp(word,"lucy") == 0)
    {
      wtFlag = LUCY;
      localTempFlag = HARMONIC;
    }
    else
      error->all(FLERR,"Illegal fix rx local temperature weighting technique");
  }

  // Determine the ODE solver/stepper strategy in arg[5].
  odeIntegrationFlag = ODE_LAMMPS_RK4;

  {
    char *word = arg[iarg++];
    if (strcmp(word,"lammps_rk4") == 0 || strcmp(word,"rk4") == 0)
      odeIntegrationFlag = ODE_LAMMPS_RK4;
    else
      error->all(FLERR,"Illegal ODE integration type");
  }

  /// Set the default ODE parameters here. Modify with arg[6].
  /// RK4:   This is the # of steps that will be taken with h = dt_dpd / minSteps;
  minSteps = 1;

  if (odeIntegrationFlag == ODE_LAMMPS_RK4 && narg==7)
  {
    // Only one option: the # of steps to take.
    char *word = arg[iarg++];
    minSteps = atoi( word );
  }
}

/* ---------------------------------------------------------------------- */

FixRX::~FixRX()
{
  // De-Allocate memory to prevent memory leak
  for (int ii = 0; ii < nreactions; ii++){
    delete [] stoich[ii];
    delete [] stoichReactants[ii];
    delete [] stoichProducts[ii];
  }
  delete [] stoich;
  delete [] stoichReactants;
  delete [] stoichProducts;
  delete [] kR;
}  

/* ---------------------------------------------------------------------- */

void FixRX::post_constructor()
{
  int maxspecies = 1000;
  int nUniqueSpecies = 0;
  bool match;

  for (int i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"property/atom",13) == 0) 
      error->all(FLERR,"fix rx cannot be combined with fix property/atom");  

  char **tmpspecies = new char*[maxspecies];  
  for(int jj=0; jj < maxspecies; jj++)
    tmpspecies[jj] = NULL;

  // open file on proc 0

  FILE *fp;
  fp = NULL;
  if (comm->me == 0) {
    fp = force->open_potential(kineticsFile);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open rx file %s",kineticsFile);
      error->one(FLERR,str);
    }
  }

  // Assign species names to tmpspecies array and determine the number of unique species

  int n,nwords;
  char line[MAXLINE],*ptr;
  int eof = 0;
  char * word;

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

    // words = ptrs to all words in line

    nwords = 0;
    word = strtok(line," \t\n\r\f");
    while (word != NULL){
      word = strtok(NULL, " \t\n\r\f");
      match=false;
      for(int jj=0;jj<nUniqueSpecies;jj++){
        if(strcmp(word,tmpspecies[jj])==0){
          match=true;
          break;
        }
      }
      if(!match){ 
        if(nUniqueSpecies+1>=maxspecies)
          error->all(FLERR,"Exceeded the maximum number of species permitted in fix rx.");
        tmpspecies[nUniqueSpecies] = new char[strlen(word)];
        strcpy(tmpspecies[nUniqueSpecies],word);
        nUniqueSpecies++;
      }
      word = strtok(NULL, " \t\n\r\f");
      if(strcmp(word,"+") != 0 && strcmp(word,"=") != 0) break;
      word = strtok(NULL, " \t\n\r\f");
    }
  }
  atom->nspecies_dpd = nUniqueSpecies;
  nspecies = atom->nspecies_dpd;

  // new id = fix-ID + FIX_STORE_ATTRIBUTE
  // new fix group = group for this fix

  id_fix_species = NULL;
  id_fix_species_old = NULL;

  n = strlen(id) + strlen("_SPECIES") + 1;
  id_fix_species = new char[n];
  n = strlen(id) + strlen("_SPECIES_OLD") + 1;
  id_fix_species_old = new char[n];

  strcpy(id_fix_species,id);
  strcat(id_fix_species,"_SPECIES");
  strcpy(id_fix_species_old,id);
  strcat(id_fix_species_old,"_SPECIES_OLD");

  char **newarg = new char*[nspecies+5];
  char **newarg2 = new char*[nspecies+5];
  newarg[0] = id_fix_species;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "property/atom";
  newarg2[0] = id_fix_species_old;
  newarg2[1] = group->names[igroup];
  newarg2[2] = (char *) "property/atom";
  for(int ii=0; ii<nspecies; ii++){
    char str1[2+strlen(tmpspecies[ii])+1];
    char str2[2+strlen(tmpspecies[ii])+4];
    strcpy(str1,"d_");
    strcpy(str2,"d_");
    strncat(str1,tmpspecies[ii],strlen(tmpspecies[ii]));
    strncat(str2,tmpspecies[ii],strlen(tmpspecies[ii]));
    strncat(str2,"Old",3);
    newarg[ii+3] = new char[strlen(str1)];
    newarg2[ii+3] = new char[strlen(str2)];
    strcpy(newarg[ii+3],str1);
    strcpy(newarg2[ii+3],str2);
  }
  newarg[nspecies+3] = (char *) "ghost";
  newarg[nspecies+4] = (char *) "yes";
  newarg2[nspecies+3] = (char *) "ghost";
  newarg2[nspecies+4] = (char *) "yes";

  modify->add_fix(nspecies+5,newarg);
  fix_species = (FixPropertyAtom *) modify->fix[modify->nfix-1];

  modify->add_fix(nspecies+5,newarg2);
  fix_species_old = (FixPropertyAtom *) modify->fix[modify->nfix-1];

  if(nspecies==0) error->all(FLERR,"There are no rx species specified.");

  for(int jj=0;jj<maxspecies;jj++) delete tmpspecies[jj];
  delete tmpspecies;

  delete newarg;
  delete newarg2;

  read_file( kineticsFile );

  // set comm size needed by this Pair
  comm_forward = nspecies*2;
  comm_reverse = 2;
}

/* ---------------------------------------------------------------------- */

int FixRX::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRX::init()
{
  pairDPDE = (PairDPDfdtEnergy *) force->pair_match("dpd/fdt/energy",1);
  if (pairDPDE == NULL)
    error->all(FLERR,"Must use pair_style dpd/fdt/energy with fix rx");

  bool eos_flag = false;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"eos/table/rx") == 0) eos_flag = true;
  if(!eos_flag) error->all(FLERR,"fix rx requires fix eos/table/rx to be specified");  
}

/* ---------------------------------------------------------------------- */

void FixRX::setup_pre_force(int vflag)
{
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int *mask = atom->mask;
  int newton_pair = force->newton_pair;
  double tmp;
  int ii;

  if(localTempFlag){
    if (newton_pair) {
      dpdThetaLocal = new double[nlocal+nghost];
      for (ii = 0; ii < nlocal+nghost; ii++)
    dpdThetaLocal[ii] = double(0.0);
    } else {
      dpdThetaLocal = new double[nlocal];
      for (ii = 0; ii < nlocal; ii++)
    dpdThetaLocal[ii] = double(0.0);
    }
    computeLocalTemperature();
  }

  for (int id = 0; id < nlocal; id++)
    for (int ispecies=0; ispecies<nspecies; ispecies++){
      tmp = atom->dvector[ispecies][id];
      atom->dvector[ispecies+nspecies][id] = tmp;
    }
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){

      // Set the reaction rate constants to zero:  no reactions occur at step 0
      for(int irxn=0;irxn<nreactions;irxn++)
        kR[irxn] = double(0.0);
      if(odeIntegrationFlag==ODE_LAMMPS_RK4) rk4(i);
    }

  // Communicate the updated momenta and velocities to all nodes
  comm->forward_comm_fix(this);
  if(localTempFlag) delete [] dpdThetaLocal;
}

/* ---------------------------------------------------------------------- */

void FixRX::pre_force(int vflag)
{
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int *mask = atom->mask;
  double *dpdTheta = atom->dpdTheta;
  int newton_pair = force->newton_pair;
  int ii;
  double theta;

  if(localTempFlag){
    if (newton_pair) {
      dpdThetaLocal = new double[nlocal+nghost];
      for (ii = 0; ii < nlocal+nghost; ii++)
        dpdThetaLocal[ii] = double(0.0);
    } else {
      dpdThetaLocal = new double[nlocal];
      for (ii = 0; ii < nlocal; ii++)
        dpdThetaLocal[ii] = double(0.0);
    }
    computeLocalTemperature();
  }

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      if(localTempFlag) theta = dpdThetaLocal[i];
      else theta=dpdTheta[i];

      //Compute the reaction rate constants
      for(int irxn=0;irxn<nreactions;irxn++)
        kR[irxn] = Arr[irxn]*pow(theta,nArr[irxn])*exp(-Ea[irxn]/force->boltz/theta);
      if(odeIntegrationFlag==ODE_LAMMPS_RK4) rk4(i);
    }

  // Communicate the updated momenta and velocities to all nodes
  comm->forward_comm_fix(this);
  if(localTempFlag) delete [] dpdThetaLocal;
}

/* ---------------------------------------------------------------------- */

void FixRX::read_file(char *file)
{
  nreactions = 0;

  // open file on proc 0

  FILE *fp;
  fp = NULL;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open rx file %s",file);
      error->one(FLERR,str);
    }
  }

  // Count the number of reactions from kinetics file

  int n,nwords,ispecies;
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

    nreactions++;
  }

  // open file on proc 0
  if (comm->me == 0) fp = force->open_potential(file);

  // read each reaction from kinetics file
  eof=0;
  char * word;
  double tmpStoich;
  double sign;

  Arr = new double[nreactions];
  nArr = new double[nreactions];
  Ea = new double[nreactions];
  tempExp = new double[nreactions];
  stoich = new double*[nreactions];
  stoichReactants = new double*[nreactions];
  stoichProducts = new double*[nreactions];
  for (int ii=0;ii<nreactions;ii++){
    stoich[ii] = new double[nspecies];
    stoichReactants[ii] = new double[nspecies];
    stoichProducts[ii] = new double[nspecies];
  }
  kR = new double[nreactions];
  for (int ii=0;ii<nreactions;ii++){
    for (int jj=0;jj<nspecies;jj++){
      stoich[ii][jj] = double(0.0);
      stoichReactants[ii][jj] = double(0.0);
      stoichProducts[ii][jj] = double(0.0);
    }
  }

  nreactions=0;
  sign = double(-1.0);
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

    // words = ptrs to all words in line

    nwords = 0;
    word = strtok(line," \t\n\r\f");
    while (word != NULL){
      tmpStoich = atof(word);
      word = strtok(NULL, " \t\n\r\f");
      for (ispecies = 0; ispecies < nspecies; ispecies++){
        if (strcmp(word,&atom->dname[ispecies][0]) == 0){
          stoich[nreactions][ispecies] += sign*tmpStoich;
          if(sign<double(0.0))
            stoichReactants[nreactions][ispecies] += tmpStoich;
          else stoichProducts[nreactions][ispecies] += tmpStoich;
          break;
        }
      }
      if(ispecies==nspecies){
        if (comm->me) {
          fprintf(stderr,"%s mol fraction is not found in data file\n",word);
          fprintf(stderr,"nspecies=%d ispecies=%d\n",nspecies,ispecies);
        }
        error->all(FLERR,"Illegal fix rx command");
      }
      word = strtok(NULL, " \t\n\r\f");
      if(strcmp(word,"=") == 0) sign = double(1.0);
      if(strcmp(word,"+") != 0 && strcmp(word,"=") != 0){
        if(word==NULL)
          error->all(FLERR,"Missing parameters in reaction kinetic equation");
        Arr[nreactions] = atof(word);
        word = strtok(NULL, " \t\n\r\f");
        if(word==NULL)
          error->all(FLERR,"Missing parameters in reaction kinetic equation");
        nArr[nreactions]  = atof(word);
        word = strtok(NULL, " \t\n\r\f");
        if(word==NULL)
          error->all(FLERR,"Missing parameters in reaction kinetic equation");
        Ea[nreactions]  = atof(word);
        sign = double(-1.0);
        break;
      }
      word = strtok(NULL, " \t\n\r\f");
    }
    nreactions++;
  }
}

/* ---------------------------------------------------------------------- */

void FixRX::setupParams()
{
  int i,j,n;

  // set mol2param for all combinations
  // must be a single exact match to lines read from file

  memory->destroy(mol2param);
  memory->create(mol2param,nspecies,"pair:mol2param");

  for (i = 0; i < nspecies; i++) {
    n = -1;
    for (j = 0; j < nreactions; j++) {
      if (i == params[j].ispecies) {
    if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
    n = j;
      }
    }
    mol2param[i] = n;
  }
}

/* ---------------------------------------------------------------------- */

void FixRX::rk4(int id)
{
  double *k1 = new double[6*nspecies + nreactions];
  double *k2 = k1 + nspecies;
  double *k3 = k2 + nspecies;
  double *k4 = k3 + nspecies;
  double *y  = k4 + nspecies;
  double *yp = y  + nspecies;

  double *dummyArray = yp + nspecies; // Passed to the rhs function.

  const int numSteps = minSteps;

  const double h = update->dt / double(numSteps);

  // Update ConcOld
  for (int ispecies = 0; ispecies < nspecies; ispecies++)
  {
    const double tmp = atom->dvector[ispecies][id];
    atom->dvector[ispecies+nspecies][id] = tmp;
    y[ispecies] = tmp;
  }

  // Run the requested steps with h.
  for (int step = 0; step < numSteps; step++)
  {
    // k1
    rhs(0.0,y,k1,dummyArray);

    // k2
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      yp[ispecies] = y[ispecies] + double(0.5)*h*k1[ispecies];

    rhs(0.0,yp,k2,dummyArray);

    // k3
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      yp[ispecies] = y[ispecies] + double(0.5)*h*k2[ispecies];

    rhs(0.0,yp,k3,dummyArray);

    // k4
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      yp[ispecies] = y[ispecies] + h*k3[ispecies];

    rhs(0.0,yp,k4,dummyArray);

    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      y[ispecies] += h*(k1[ispecies]/6.0 + k2[ispecies]/3.0 + 
                    k3[ispecies]/3.0 + k4[ispecies]/6.0);

  } // end for (int step...

  // Store the solution back in atom->dvector.
  for (int ispecies = 0; ispecies < nspecies; ispecies++){
    if(y[ispecies] < double(-1.0e-10))
      error->one(FLERR,"Computed concentration in RK4 solver is < -1.0e-10");
    else if(y[ispecies] < double(0.0))
      y[ispecies] = double(0.0);
    atom->dvector[ispecies][id] = y[ispecies];
  }
  delete [] k1;
}

/* ---------------------------------------------------------------------- */

int FixRX::rhs(double t, const double *y, double *dydt, void *params)
{
  double rxnRateLawForward;
  double *rxnRateLaw = (double *) params;
  double VDPD = domain->xprd * domain->yprd * domain->zprd / atom->natoms;
  double concentration;
  int nspecies = atom->nspecies_dpd;

  for(int ispecies=0; ispecies<nspecies; ispecies++)
    dydt[ispecies] = double(0.0);

  // Construct the reaction rate laws
  for(int jrxn=0; jrxn<nreactions; jrxn++){
    rxnRateLawForward = kR[jrxn];

    for(int ispecies=0; ispecies<nspecies; ispecies++){
      concentration = y[ispecies]/VDPD;
      rxnRateLawForward *= pow(concentration,stoichReactants[jrxn][ispecies]);
    }
    rxnRateLaw[jrxn] = rxnRateLawForward;
  }
  
  // Construct the reaction rates for each species
  for(int ispecies=0; ispecies<nspecies; ispecies++)
    for(int jrxn=0; jrxn<nreactions; jrxn++)
      dydt[ispecies] += stoich[jrxn][ispecies]*VDPD*rxnRateLaw[jrxn];

  return 0;
}

/* ---------------------------------------------------------------------- */

void FixRX::computeLocalTemperature()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;

  // local temperature variables
  double wij;
  double *dpdTheta = atom->dpdTheta;

  // Initialize the local density and local temperature arrays
  if (newton_pair) {
    sumWeights = new double[nlocal+nghost];
    for (ii = 0; ii < nlocal+nghost; ii++)
      sumWeights[ii] = double(0.0);
  } else {
    sumWeights = new double[nlocal];
    for (ii = 0; ii < nlocal; ii++)
      dpdThetaLocal[ii] = double(0.0);
  }

  inum = pairDPDE->list->inum;
  ilist = pairDPDE->list->ilist;
  numneigh = pairDPDE->list->numneigh;
  firstneigh = pairDPDE->list->firstneigh;

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < pairDPDE->cutsq[itype][jtype]) {
        double rcut = sqrt(pairDPDE->cutsq[itype][jtype]);
    double rij = sqrt(rsq);
    double ratio = rij/rcut;

    // Lucy's Weight Function
    if(wtFlag==LUCY){
      wij = (double(1.0)+double(3.0)*ratio) * (double(1.0)-ratio)*(double(1.0)-ratio)*(double(1.0)-ratio);
      dpdThetaLocal[i] += wij/dpdTheta[j];
      if (newton_pair || j < nlocal) 
        dpdThetaLocal[j] += wij/dpdTheta[i];
    }

    sumWeights[i] += wij;
        if (newton_pair || j < nlocal) {
      sumWeights[j] += wij;
    }
      
      }
    }
  }
  if (newton_pair) comm->reverse_comm_fix(this);

  // self-interaction for local temperature
  for (i = 0; i < nlocal; i++){

    // Lucy Weight Function
    if(wtFlag==LUCY){
      wij = double(1.0);
      dpdThetaLocal[i] += wij / dpdTheta[i];
    }
    sumWeights[i] += wij;
    
    // Normalized local temperature
    dpdThetaLocal[i] = dpdThetaLocal[i] / sumWeights[i];
    
    if(localTempFlag == HARMONIC)
      dpdThetaLocal[i] = double(1.0) / dpdThetaLocal[i];
    
  }

  delete [] sumWeights;
}

/* ---------------------------------------------------------------------- */

int FixRX::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int ii,jj,m;
  double tmp;

  m = 0;
  for (ii = 0; ii < n; ii++) {
    jj = list[ii];
    for(int ispecies=0;ispecies<nspecies;ispecies++){
      tmp = atom->dvector[ispecies][jj];
      buf[m++] = tmp;
      tmp = atom->dvector[ispecies+nspecies][jj];
      buf[m++] = tmp;
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRX::unpack_forward_comm(int n, int first, double *buf)
{
  int ii,m,last;
  double tmp;

  m = 0;
  last = first + n ;
  for (ii = first; ii < last; ii++){
    for(int ispecies=0;ispecies<nspecies;ispecies++){
      tmp = buf[m++];
      atom->dvector[ispecies][ii] = tmp;
      tmp = buf[m++];
      atom->dvector[ispecies+nspecies][ii] = tmp;
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixRX::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = dpdThetaLocal[i];
    buf[m++] = sumWeights[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRX::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    dpdThetaLocal[j] += buf[m++];
    sumWeights[j] += buf[m++];
  }
}
