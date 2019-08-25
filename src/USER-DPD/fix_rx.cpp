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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat> // DBL_EPSILON
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
#include "neigh_request.h"
#include "math_special.h"
#include "pair_dpd_fdt_energy.h"

#include <vector> // std::vector<>
#include <algorithm> // std::max

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathSpecial;

enum{NONE,HARMONIC};
enum{LUCY};

#define MAXLINE 1024
#define DELTA 4

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

#define SparseKinetics_enableIntegralReactions (true)
#define SparseKinetics_invalidIndex (-1)

namespace /* anonymous */
{

typedef double TimerType;
TimerType getTimeStamp(void) { return MPI_Wtime(); }
double getElapsedTime( const TimerType &t0, const TimerType &t1) { return t1-t0; }

} // end namespace

/* ---------------------------------------------------------------------- */

FixRX::FixRX(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), mol2param(NULL), nreactions(0),
  params(NULL), Arr(NULL), nArr(NULL), Ea(NULL), tempExp(NULL),
  stoich(NULL), stoichReactants(NULL), stoichProducts(NULL), kR(NULL),
  pairDPDE(NULL), dpdThetaLocal(NULL), sumWeights(NULL), sparseKinetics_nu(NULL),
  sparseKinetics_nuk(NULL), sparseKinetics_inu(NULL), sparseKinetics_isIntegralReaction(NULL),
  kineticsFile(NULL), id_fix_species(NULL),
  id_fix_species_old(NULL), fix_species(NULL), fix_species_old(NULL)
{
  if (narg < 7 || narg > 12) error->all(FLERR,"Illegal fix rx command");
  nevery = 1;

  nreactions = maxparam = 0;
  params = NULL;
  mol2param = NULL;
  pairDPDE = NULL;
  id_fix_species = NULL;
  id_fix_species_old = NULL;

  const int Verbosity = 1;

  // Keep track of the argument list.
  int iarg = 3;

  // Read the kinetic file in arg[3].
  kineticsFile = arg[iarg++];

  // Determine the local temperature averaging method in arg[4].
  wtFlag = 0;
  localTempFlag = NONE;

  {
    char *word = arg[iarg++];
    if (strcmp(word,"none") == 0){
      wtFlag = 0;
      localTempFlag = NONE;
    }
    else if (strcmp(word,"lucy") == 0){
      wtFlag = LUCY;
      localTempFlag = HARMONIC;
    }
    else
      error->all(FLERR,"Illegal fix rx local temperature weighting technique");
  }

  // Select either sparse and dense matrix
  // representations of the stoichiometric matrix.
  useSparseKinetics = true;
  {
    char *word = arg[iarg++];

    if (strcmp(word,"sparse") == 0)
      useSparseKinetics = true;
    else if (strcmp(word,"dense") == 0)
      useSparseKinetics = false;
    else {
      std::string errmsg = "Illegal command " + std::string(word)
                             + " expected \"sparse\" or \"dense\"\n";
      error->all(FLERR, errmsg.c_str());
    }

    if (comm->me == 0 and Verbosity > 1){
      std::string msg = "FixRX: matrix format is ";
      if (useSparseKinetics)
         msg += std::string("sparse");
      else
         msg += std::string("dense");

      error->message(FLERR, msg.c_str());
    }
  }

  // Determine the ODE solver/stepper strategy in arg[6].
  odeIntegrationFlag = ODE_LAMMPS_RK4;

  {
    char *word = arg[iarg++];
    if (strcmp(word,"lammps_rk4") == 0 || strcmp(word,"rk4") == 0)
      odeIntegrationFlag = ODE_LAMMPS_RK4;
    else if (strcmp(word,"lammps_rkf45") == 0 || strcmp(word,"rkf45") == 0)
      odeIntegrationFlag = ODE_LAMMPS_RKF45;
    else {
      std::string errmsg = "Illegal ODE integration type: " + std::string(word);
      error->all(FLERR, errmsg.c_str());
    }
  }

  /// Set the default ODE parameters here. Modify with arg[].
  /// 'minSteps' has a different meaning for RK4 and RKF45.
  /// RK4:   This is the # of steps that will be taken with h = dt_dpd / minSteps;
  /// RKF45: This sets as h0 = dt_dpd / minSteps. If minSteps == 0, RKF45 will
  ///        estimate h0 internally. h will be adjusted as needed on subsequent steps.
  minSteps = 1;
  maxIters = 100;
  relTol   = 1.0e-6;
  absTol   = 1.0e-8;

  diagnosticFrequency = 0;
  for (int i = 0; i < numDiagnosticCounters; ++i){
    diagnosticCounter[i] = 0;
    diagnosticCounterPerODE[i] = NULL;
  }

  if (odeIntegrationFlag == ODE_LAMMPS_RK4 && narg==8){
    char *word = arg[iarg++];
    minSteps = atoi( word );

    if (comm->me == 0 and Verbosity > 1){
      char msg[128];
      sprintf(msg, "FixRX: RK4 numSteps= %d", minSteps);
      error->message(FLERR, msg);
    }
  }
  else if (odeIntegrationFlag == ODE_LAMMPS_RK4 && narg>8){
    error->all(FLERR,"Illegal fix rx command.  Too many arguments for RK4 solver.");
  }
  else if (odeIntegrationFlag == ODE_LAMMPS_RKF45){
    // Must have four options.
    if (narg < 11)
      error->all(FLERR,"Illegal fix rx command.  Too few arguments for RKF45 solver.");

    minSteps = atoi( arg[iarg++] );
    maxIters = atoi( arg[iarg++] );
    relTol   = strtod( arg[iarg++], NULL);
    absTol   = strtod( arg[iarg++], NULL);

    if (iarg < narg)
      diagnosticFrequency = atoi( arg[iarg++] );

    // maxIters must be at least minSteps.
    maxIters = std::max( minSteps, maxIters );

    if (comm->me == 0 and Verbosity > 1){
      //printf("FixRX: RKF45 minSteps= %d maxIters= %d absTol= %e relTol= %e\n", minSteps, maxIters, absTol, relTol);
      char msg[128];
      sprintf(msg, "FixRX: RKF45 minSteps= %d maxIters= %d relTol= %.1e absTol= %.1e diagnosticFrequency= %d", minSteps, maxIters, relTol, absTol, diagnosticFrequency);
      error->message(FLERR, msg);
    }
  }

  // Initialize/Create the sparse matrix database.
  sparseKinetics_nu = NULL;
  sparseKinetics_nuk = NULL;
  sparseKinetics_inu = NULL;
  sparseKinetics_isIntegralReaction = NULL;
  sparseKinetics_maxReactants = 0;
  sparseKinetics_maxProducts = 0;
  sparseKinetics_maxSpecies = 0;
}

/* ---------------------------------------------------------------------- */

FixRX::~FixRX()
{
  //printf("Inside FixRX::~FixRX copymode= %d\n", copymode);
  if (copymode) return;

  // De-Allocate memory to prevent memory leak
  for (int ii = 0; ii < nreactions; ii++){
    delete [] stoich[ii];
    delete [] stoichReactants[ii];
    delete [] stoichProducts[ii];
  }
  delete [] Arr;
  delete [] nArr;
  delete [] Ea;
  delete [] tempExp;
  delete [] stoich;
  delete [] stoichReactants;
  delete [] stoichProducts;
  delete [] kR;
  delete [] id_fix_species;
  delete [] id_fix_species_old;

  if (useSparseKinetics){
     memory->destroy( sparseKinetics_nu );
     memory->destroy( sparseKinetics_nuk );
     memory->destroy( sparseKinetics_inu );
     memory->destroy( sparseKinetics_isIntegralReaction );
  }
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
  int tmpmaxstrlen = 0;
  for(int jj=0; jj < maxspecies; jj++)
    tmpspecies[jj] = NULL;

  // open file on proc 0

  FILE *fp;
  fp = NULL;
  if (comm->me == 0) {
    fp = force->open_potential(kineticsFile);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open rx file %s",kineticsFile);
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
        tmpspecies[nUniqueSpecies] = new char[strlen(word)+1];
        strcpy(tmpspecies[nUniqueSpecies],word);
        tmpmaxstrlen = MAX(tmpmaxstrlen,strlen(word));
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
  char *str1 = new char[tmpmaxstrlen+3];
  char *str2 = new char[tmpmaxstrlen+6];
  for(int ii=0; ii<nspecies; ii++){
    strcpy(str1,"d_");
    strcpy(str2,"d_");
    strcat(str1,tmpspecies[ii]);
    strcat(str2,tmpspecies[ii]);
    strcat(str2,"Old");
    newarg[ii+3] = new char[strlen(str1)+1];
    newarg2[ii+3] = new char[strlen(str2)+1];
    strcpy(newarg[ii+3],str1);
    strcpy(newarg2[ii+3],str2);
  }
  delete[] str1;
  delete[] str2;
  newarg[nspecies+3] = (char *) "ghost";
  newarg[nspecies+4] = (char *) "yes";
  newarg2[nspecies+3] = (char *) "ghost";
  newarg2[nspecies+4] = (char *) "yes";

  modify->add_fix(nspecies+5,newarg,1);
  fix_species = (FixPropertyAtom *) modify->fix[modify->nfix-1];
  restartFlag = modify->fix[modify->nfix-1]->restart_reset;

  modify->add_fix(nspecies+5,newarg2,1);
  fix_species_old = (FixPropertyAtom *) modify->fix[modify->nfix-1];

  if(nspecies==0) error->all(FLERR,"There are no rx species specified.");

  for(int jj=0;jj<nspecies;jj++) {
    delete[] tmpspecies[jj];
    delete[] newarg[jj+3];
    delete[] newarg2[jj+3];
  }

  delete[] newarg;
  delete[] newarg2;
  delete[] tmpspecies;

  read_file( kineticsFile );

  if (useSparseKinetics)
    this->initSparse();

  // set comm size needed by this Pair
  comm_forward = nspecies*2;
  comm_reverse = 2;
}

/* ---------------------------------------------------------------------- */

void FixRX::initSparse()
{
  const int Verbosity = 1;

  if (comm->me == 0 and Verbosity > 1){
    for (int k = 0; k < nspecies; ++k)
      printf("atom->dname[%d]= %s\n", k, atom->dname[k]);

    printf("stoich[][]\n");
    for (int i = 0; i < nreactions; ++i){
      int nreac_i = 0, nprod_i = 0;
      printf("%d: ", i);
      for (int k = 0; k < nspecies; ++k){
         printf(" %g", stoich[i][k]);
         if (stoich[i][k] < 0.0) nreac_i++;
         else if (stoich[i][k] > 0.0) nprod_i++;
      }
      printf(" : %d %d\n", nreac_i, nprod_i);
    }

    printf("stoichReactants[][]\n");
    for (int i = 0; i < nreactions; ++i){
      int nreac_i = 0;
      printf("%d: ", i);
      for (int k = 0; k < nspecies; ++k){
         printf(" %g", stoichReactants[i][k]);
         if (stoichReactants[i][k] > 0.0) nreac_i++;
      }
      printf(" : %d\n", nreac_i);
    }

    printf("stoichProducts[][]\n");
    for (int i = 0; i < nreactions; ++i){
      int nprod_i = 0;
      printf("%d: ", i);
      for (int k = 0; k < nspecies; ++k){
         printf(" %g", stoichProducts[i][k]);
         if (stoichProducts[i][k] > 0.0) nprod_i++;
      }
      printf(" : %d\n", nprod_i);
    }
  } // if (Verbose)

  // 1) Measure the sparsity of stoich[][]
  int nzeros = 0;
  int mxprod = 0;
  int mxreac = 0;
  int mxspec = 0;
  int nIntegral = 0;
  for (int i = 0; i < nreactions; ++i){
    int nreac_i = 0, nprod_i = 0;
    std::string pstr, rstr;
    bool allAreIntegral = true;
    for (int k = 0; k < nspecies; ++k){
      if (stoichReactants[i][k] == 0 and stoichProducts[i][k] == 0)
        nzeros++;

      if (stoichReactants[i][k] > 0.0){
        allAreIntegral &= (std::fmod( stoichReactants[i][k], 1.0 ) == 0.0);

        nreac_i++;
        if (rstr.length() > 0)
          rstr += " + ";

        char digit[6];
        sprintf(digit, "%4.1f ", stoichReactants[i][k]); rstr += digit;
        rstr += atom->dname[k];
      }
      if (stoichProducts[i][k] > 0.0){
        allAreIntegral &= (std::fmod( stoichProducts[i][k], 1.0 ) == 0.0);

        nprod_i++;
        if (pstr.length() > 0)
          pstr += " + ";

        char digit[6];
        sprintf(digit, "%4.1f ", stoichProducts[i][k]); pstr += digit;

        pstr += atom->dname[k];
      }
    }
    if (comm->me == 0 and Verbosity > 1)
      printf("rx%3d: %d %d %d ... %s %s %s\n", i, nreac_i, nprod_i, allAreIntegral, rstr.c_str(), /*reversible[i]*/ (false) ? "<=>" : "=", pstr.c_str());

    mxreac = std::max( mxreac, nreac_i );
    mxprod = std::max( mxprod, nprod_i );
    mxspec = std::max( mxspec, nreac_i + nprod_i );
    if (allAreIntegral) nIntegral++;
  }

  if (comm->me == 0 and Verbosity > 1){
    char msg[256];
    sprintf(msg, "FixRX: Sparsity of Stoichiometric Matrix= %.1f%% non-zeros= %d nspecies= %d nreactions= %d maxReactants= %d maxProducts= %d maxSpecies= %d integralReactions= %d", 100*(double(nzeros) / (nspecies * nreactions)), nzeros, nspecies, nreactions, mxreac, mxprod, (mxreac + mxprod), SparseKinetics_enableIntegralReactions);
    error->message(FLERR, msg);
  }

  // Allocate the sparse matrix data.
  {
     sparseKinetics_maxSpecies   = (mxreac + mxprod);
     sparseKinetics_maxReactants = mxreac;
     sparseKinetics_maxProducts  = mxprod;

     memory->create( sparseKinetics_nu , nreactions, sparseKinetics_maxSpecies, "sparseKinetics_nu");
     memory->create( sparseKinetics_nuk, nreactions, sparseKinetics_maxSpecies, "sparseKinetics_nuk");

     for (int i = 0; i < nreactions; ++i)
        for (int k = 0; k < sparseKinetics_maxSpecies; ++k){
           sparseKinetics_nu [i][k] = 0.0;
           sparseKinetics_nuk[i][k] = SparseKinetics_invalidIndex; // Initialize with an invalid index.
        }

     if (SparseKinetics_enableIntegralReactions){
        memory->create( sparseKinetics_inu, nreactions, sparseKinetics_maxSpecies, "sparseKinetics_inu");
        memory->create( sparseKinetics_isIntegralReaction, nreactions, "sparseKinetics_isIntegralReaction");

        for (int i = 0; i < nreactions; ++i){
           sparseKinetics_isIntegralReaction[i] = false;
           for (int k = 0; k < sparseKinetics_maxSpecies; ++k)
              sparseKinetics_inu[i][k] = 0;
        }
     }
  }

  // Measure the distribution of the # of moles for the ::fastpowi function.
  std::vector<int> nu_bin(10);

  for (int i = 0; i < nreactions; ++i){
    int nreac_i = 0, nprod_i = 0;
    bool isIntegral_i = true;
    for (int k = 0; k < nspecies; ++k){
      if (stoichReactants[i][k] > 0.0){
        const int idx = nreac_i;
        sparseKinetics_nu [i][idx] = stoichReactants[i][k];
        sparseKinetics_nuk[i][idx] = k;

        isIntegral_i &= (std::fmod( stoichReactants[i][k], 1.0 ) == 0.0);
        if (SparseKinetics_enableIntegralReactions){
          sparseKinetics_inu[i][idx] = (int)sparseKinetics_nu[i][idx];
          if (isIntegral_i){
            if (sparseKinetics_inu[i][idx] >= nu_bin.size())
               nu_bin.resize( sparseKinetics_inu[i][idx] );

            nu_bin[ sparseKinetics_inu[i][idx] ] ++;
          }
        }

        nreac_i++;
      }
      if (stoichProducts[i][k] > 0.0){
        const int idx = sparseKinetics_maxReactants + nprod_i;
        sparseKinetics_nu [i][idx] = stoichProducts[i][k];
        sparseKinetics_nuk[i][idx] = k;

        isIntegral_i &= (std::fmod( sparseKinetics_nu[i][idx], 1.0 ) == 0.0);
        if (SparseKinetics_enableIntegralReactions){
          sparseKinetics_inu[i][idx] = (int) sparseKinetics_nu[i][idx];
          if (isIntegral_i){
            if (sparseKinetics_inu[i][idx] >= nu_bin.size())
               nu_bin.resize( sparseKinetics_inu[i][idx] );

            nu_bin[ sparseKinetics_inu[i][idx] ] ++;
          }
        }

        nprod_i++;
      }
    }

    if (SparseKinetics_enableIntegralReactions)
       sparseKinetics_isIntegralReaction[i] = isIntegral_i;
  }

  if (comm->me == 0 and Verbosity > 1){
    for (int i = 1; i < nu_bin.size(); ++i)
      if (nu_bin[i] > 0)
        printf("nu_bin[%d] = %d\n", i, nu_bin[i]);

    for (int i = 0; i < nreactions; ++i){
      std::string pstr, rstr;

      for (int kk = 0; kk < sparseKinetics_maxReactants; kk++){
        const int k = sparseKinetics_nuk[i][kk];
        if (k != SparseKinetics_invalidIndex){
          if (rstr.length() > 0)
            rstr += " + ";

          char digit[6];
          if (SparseKinetics_enableIntegralReactions and sparseKinetics_isIntegralReaction[i])
            sprintf(digit,"%d ", sparseKinetics_inu[i][kk]);
          else
            sprintf(digit,"%4.1f ", sparseKinetics_nu[i][kk]);
          rstr += digit;
          rstr += atom->dname[k];
        }
      }

      for (int kk = sparseKinetics_maxReactants; kk < sparseKinetics_maxSpecies; kk++){
        const int k = sparseKinetics_nuk[i][kk];
        if (k != SparseKinetics_invalidIndex){
          if (pstr.length() > 0)
            pstr += " + ";

          char digit[6];
          if (SparseKinetics_enableIntegralReactions and sparseKinetics_isIntegralReaction[i])
            sprintf(digit,"%d ", sparseKinetics_inu[i][kk]);
          else
            sprintf(digit,"%4.1f ", sparseKinetics_nu[i][kk]);
          pstr += digit;
          pstr += atom->dname[k];
        }
      }
      if (comm->me == 0 and Verbosity > 1)
        printf("rx%3d: %s %s %s\n", i, rstr.c_str(), /*reversible[i]*/ (false) ? "<=>" : "=", pstr.c_str());
    }
    // end for nreactions
  }
  // end if Verbose
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
    pairDPDE = (PairDPDfdtEnergy *) force->pair_match("dpd/fdt/energy/kk",1);

  if (pairDPDE == NULL)
    error->all(FLERR,"Must use pair_style dpd/fdt/energy with fix rx");

  bool eos_flag = false;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"eos/table/rx") == 0) eos_flag = true;
  if(!eos_flag) error->all(FLERR,"fix rx requires fix eos/table/rx to be specified");

  // need a half neighbor list
  // built whenever re-neighboring occurs

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
}

/* ---------------------------------------------------------------------- */

void FixRX::init_list(int, class NeighList* ptr)
{
  this->list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixRX::setup_pre_force(int /*vflag*/)
{
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int *mask = atom->mask;
  int newton_pair = force->newton_pair;
  double tmp;

  if(restartFlag){
    restartFlag = 0;
  }
  else
  {
    int ode_counter[4] = {0};

    UserRHSData userData;
    userData.kFor = new double[nreactions];
    userData.rxnRateLaw = new double[nreactions];

    double *rwork = new double[8*nspecies];

    if(localTempFlag){
      int count = nlocal + (newton_pair ? nghost : 0);
      dpdThetaLocal = new double[count];
      memset(dpdThetaLocal, 0, sizeof(double)*count);
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
          userData.kFor[irxn] = 0.0;

        if (odeIntegrationFlag == ODE_LAMMPS_RK4)
          rk4(i, rwork, &userData);
        else if (odeIntegrationFlag == ODE_LAMMPS_RKF45)
          rkf45(i, rwork, &userData, ode_counter);
      }

    // Communicate the updated momenta and velocities to all nodes
    comm->forward_comm_fix(this);
    if(localTempFlag) delete [] dpdThetaLocal;

    delete [] userData.kFor;
    delete [] userData.rxnRateLaw;
    delete [] rwork;
  }
}

/* ---------------------------------------------------------------------- */

void FixRX::pre_force(int /*vflag*/)
{
  //TimerType timer_start = getTimeStamp();

  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int *mask = atom->mask;
  double *dpdTheta = atom->dpdTheta;
  int newton_pair = force->newton_pair;

  if(localTempFlag){
    int count = nlocal + (newton_pair ? nghost : 0);
    dpdThetaLocal = new double[count];
    memset(dpdThetaLocal, 0, sizeof(double)*count);
    computeLocalTemperature();
  }

  TimerType timer_localTemperature = getTimeStamp();

  // Zero the counters for the ODE solvers.
  int nSteps = 0;
  int nIters = 0;
  int nFuncs = 0;
  int nFails = 0;

  if (odeIntegrationFlag == ODE_LAMMPS_RKF45 && diagnosticFrequency == 1)
  {
    memory->create( diagnosticCounterPerODE[StepSum], nlocal, "FixRX::diagnosticCounterPerODE");
    memory->create( diagnosticCounterPerODE[FuncSum], nlocal, "FixRX::diagnosticCounterPerODE");
  }

#if 0
  #pragma omp parallel \
     reduction(+: nSteps, nIters, nFuncs, nFails )
#endif
  {
    double *rwork = new double[8*nspecies];

    UserRHSData userData;
    userData.kFor = new double[nreactions];
    userData.rxnRateLaw = new double[nreactions];

    int ode_counter[4] = { 0 };

    //#pragma omp for schedule(runtime)
    for (int i = 0; i < nlocal; i++)
    {
      if (mask[i] & groupbit)
      {
        double theta;
        if (localTempFlag)
          theta = dpdThetaLocal[i];
        else
          theta = dpdTheta[i];

        //Compute the reaction rate constants
        for (int irxn = 0; irxn < nreactions; irxn++)
          userData.kFor[irxn] = Arr[irxn]*pow(theta,nArr[irxn])*exp(-Ea[irxn]/force->boltz/theta);

        if (odeIntegrationFlag == ODE_LAMMPS_RK4)
          rk4(i, rwork, &userData);
        else if (odeIntegrationFlag == ODE_LAMMPS_RKF45)
          rkf45(i, rwork, &userData, ode_counter);
      }
    }

    nSteps += ode_counter[0];
    nIters += ode_counter[1];
    nFuncs += ode_counter[2];
    nFails += ode_counter[3];

    delete [] rwork;
    delete [] userData.kFor;
    delete [] userData.rxnRateLaw;

  } // end parallel region

  TimerType timer_ODE = getTimeStamp();

  // Communicate the updated momenta and velocities to all nodes
  comm->forward_comm_fix(this);
  if(localTempFlag) delete [] dpdThetaLocal;

  //TimerType timer_stop = getTimeStamp();

  double time_ODE = getElapsedTime(timer_localTemperature, timer_ODE);

  //printf("me= %d total= %g temp= %g ode= %g comm= %g nlocal= %d nfc= %d %d\n", comm->me,
  //                       getElapsedTime(timer_start, timer_stop),
  //                       getElapsedTime(timer_start, timer_localTemperature),
  //                       getElapsedTime(timer_localTemperature, timer_ODE),
  //                       getElapsedTime(timer_ODE, timer_stop), nlocal, nFuncs, nSteps);

  // Warn the user if a failure was detected in the ODE solver.
  if (nFails > 0){
    char sbuf[128];
    sprintf(sbuf,"in FixRX::pre_force, ODE solver failed for %d atoms.", nFails);
    error->warning(FLERR, sbuf);
  }

  // Compute and report ODE diagnostics, if requested.
  if (odeIntegrationFlag == ODE_LAMMPS_RKF45 && diagnosticFrequency != 0){
    // Update the counters.
    diagnosticCounter[StepSum] += nSteps;
    diagnosticCounter[FuncSum] += nFuncs;
    diagnosticCounter[TimeSum] += time_ODE;
    diagnosticCounter[AtomSum] += nlocal;
    diagnosticCounter[numDiagnosticCounters-1] ++;

    if ( (diagnosticFrequency > 0 &&
               ((update->ntimestep - update->firststep) % diagnosticFrequency) == 0) ||
         (diagnosticFrequency < 0 && update->ntimestep == update->laststep) )
      this->odeDiagnostics();

    for (int i = 0; i < numDiagnosticCounters; ++i)
      if (diagnosticCounterPerODE[i])
        memory->destroy( diagnosticCounterPerODE[i] );
  }
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
      snprintf(str,128,"Cannot open rx file %s",file);
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
      stoich[ii][jj] = 0.0;
      stoichReactants[ii][jj] = 0.0;
      stoichProducts[ii][jj] = 0.0;
    }
  }

  nreactions=0;
  sign = -1.0;
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
          if(sign<0.0)
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
      if(word==NULL) error->all(FLERR,"Missing parameters in reaction kinetic equation");
      if(strcmp(word,"=") == 0) sign = 1.0;
      if(strcmp(word,"+") != 0 && strcmp(word,"=") != 0){
        if(word==NULL) error->all(FLERR,"Missing parameters in reaction kinetic equation");
        Arr[nreactions] = atof(word);
        word = strtok(NULL, " \t\n\r\f");
        if(word==NULL) error->all(FLERR,"Missing parameters in reaction kinetic equation");
        nArr[nreactions]  = atof(word);
        word = strtok(NULL, " \t\n\r\f");
        if(word==NULL) error->all(FLERR,"Missing parameters in reaction kinetic equation");
        Ea[nreactions]  = atof(word);
        sign = -1.0;
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

void FixRX::rk4(int id, double *rwork, void* v_params)
{
  double *k1 = rwork;
  double *k2 = k1 + nspecies;
  double *k3 = k2 + nspecies;
  double *k4 = k3 + nspecies;
  double *y  = k4 + nspecies;
  double *yp = y  + nspecies;

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
    rhs(0.0,y,k1,v_params);

    // k2
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      yp[ispecies] = y[ispecies] + 0.5*h*k1[ispecies];

    rhs(0.0,yp,k2,v_params);

    // k3
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      yp[ispecies] = y[ispecies] + 0.5*h*k2[ispecies];

    rhs(0.0,yp,k3,v_params);

    // k4
    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      yp[ispecies] = y[ispecies] + h*k3[ispecies];

    rhs(0.0,yp,k4,v_params);

    for (int ispecies = 0; ispecies < nspecies; ispecies++)
      y[ispecies] += h*(k1[ispecies]/6.0 + k2[ispecies]/3.0 + k3[ispecies]/3.0 + k4[ispecies]/6.0);

  } // end for (int step...

  // Store the solution back in atom->dvector.
  for (int ispecies = 0; ispecies < nspecies; ispecies++){
    if(y[ispecies] < -MY_EPSILON)
      error->one(FLERR,"Computed concentration in RK4 solver is < -10*DBL_EPSILON");
    else if(y[ispecies] < MY_EPSILON)
      y[ispecies] = 0.0;
    atom->dvector[ispecies][id] = y[ispecies];
  }
}

/* ---------------------------------------------------------------------- */

//     f1 = dt*f(t,x)
//     f2 = dt*f(t+ c20*dt,x + c21*f1)
//     f3 = dt*f(t+ c30*dt,x + c31*f1 + c32*f2)
//     f4 = dt*f(t+ c40*dt,x + c41*f1 + c42*f2 + c43*f3)
//     f5 = dt*f(t+dt,x + c51*f1 + c52*f2 + c53*f3 + c54*f4)
//     f6 = dt*f(t+ c60*dt,x + c61*f1 + c62*f2 + c63*f3 + c64*f4 + c65*f5)
//
//     fifth-order runge-kutta integration
//        x5 = x + b1*f1 + b3*f3 + b4*f4 + b5*f5 + b6*f6
//     fourth-order runge-kutta integration
//        x  = x + a1*f1 + a3*f3 + a4*f4 + a5*f5

void FixRX::rkf45_step (const int neq, const double h, double y[], double y_out[], double rwk[], void* v_param)
{
   const double c21=0.25;
   const double c31=0.09375;
   const double c32=0.28125;
   const double c41=0.87938097405553;
   const double c42=-3.2771961766045;
   const double c43=3.3208921256258;
   const double c51=2.0324074074074;
   const double c52=-8.0;
   const double c53=7.1734892787524;
   const double c54=-0.20589668615984;
   const double c61=-0.2962962962963;
   const double c62=2.0;
   const double c63=-1.3816764132554;
   const double c64=0.45297270955166;
   const double c65=-0.275;
   const double a1=0.11574074074074;
   const double a3=0.54892787524366;
   const double a4=0.5353313840156;
   const double a5=-0.2;
   const double b1=0.11851851851852;
   const double b3=0.51898635477583;
   const double b4=0.50613149034201;
   const double b5=-0.18;
   const double b6=0.036363636363636;

   // local dependent variables (5 total)
   double* f1 = &rwk[    0];
   double* f2 = &rwk[  neq];
   double* f3 = &rwk[2*neq];
   double* f4 = &rwk[3*neq];
   double* f5 = &rwk[4*neq];
   double* f6 = &rwk[5*neq];

   // scratch for the intermediate solution.
   //double* ytmp = &rwk[6*neq];
   double* ytmp = y_out;

   // 1)
   rhs (0.0, y, f1, v_param);

   for (int k = 0; k < neq; k++){
      f1[k] *= h;
      ytmp[k] = y[k] + c21 * f1[k];
   }

   // 2)
   rhs(0.0, ytmp, f2, v_param);

   for (int k = 0; k < neq; k++){
      f2[k] *= h;
      ytmp[k] = y[k] + c31 * f1[k] + c32 * f2[k];
   }

   // 3)
   rhs(0.0, ytmp, f3, v_param);

   for (int k = 0; k < neq; k++) {
      f3[k] *= h;
      ytmp[k] = y[k] + c41 * f1[k] + c42 * f2[k] + c43 * f3[k];
   }

   // 4)
   rhs(0.0, ytmp, f4, v_param);

   for (int k = 0; k < neq; k++) {
      f4[k] *= h;
      ytmp[k] = y[k] + c51 * f1[k] + c52 * f2[k] + c53 * f3[k] + c54 * f4[k];
   }

   // 5)
   rhs(0.0, ytmp, f5, v_param);

   for (int k = 0; k < neq; k++) {
      f5[k] *= h;
      ytmp[k] = y[k] + c61*f1[k] + c62*f2[k] + c63*f3[k] + c64*f4[k] + c65*f5[k];
   }

   // 6)
   rhs(0.0, ytmp, f6, v_param);

   for (int k = 0; k < neq; k++)
   {
      //const double f6 = h * ydot[k];
      f6[k] *= h;

      // 5th-order solution.
      const double r5 = b1*f1[k] + b3*f3[k] + b4*f4[k] + b5*f5[k] + b6*f6[k];

      // 4th-order solution.
      const double r4 = a1*f1[k] + a3*f3[k] + a4*f4[k] + a5*f5[k];

      // Truncation error: difference between 4th and 5th-order solutions.
      rwk[k] = fabs(r5 - r4);

      // Update solution.
    //y_out[k] = y[k] + r5; // Local extrapolation
      y_out[k] = y[k] + r4;
   }

   return;
}

int FixRX::rkf45_h0 (const int neq, const double t, const double /*t_stop*/,
                     const double hmin, const double hmax,
                     double& h0, double y[], double rwk[], void* v_params)
{
   // Set lower and upper bounds on h0, and take geometric mean as first trial value.
   // Exit with this value if the bounds cross each other.

   // Adjust upper bound based on ydot ...
   double hg = sqrt(hmin*hmax);

   //if (hmax < hmin)
   //{
   //   h0 = hg;
   //   return;
   //}

   // Start iteration to find solution to ... {WRMS norm of (h0^2 y'' / 2)} = 1

   double *ydot  = rwk;
   double *y1    = ydot + neq;
   double *ydot1 = y1 + neq;

   const int max_iters = 10;
   bool hnew_is_ok = false;
   double hnew = hg;
   int iter = 0;

   // compute ydot at t=t0
   rhs (t, y, ydot, v_params);

   while(1)
   {
      // Estimate y'' with finite-difference ...

      for (int k = 0; k < neq; k++)
         y1[k] = y[k] + hg * ydot[k];

      // compute y' at t1
      rhs (t + hg, y1, ydot1, v_params);

      // Compute WRMS norm of y''
      double yddnrm = 0.0;
      for (int k = 0; k < neq; k++){
         double ydd = (ydot1[k] - ydot[k]) / hg;
         double wterr = ydd / (relTol * fabs( y[k] ) + absTol);
         yddnrm += wterr * wterr;
      }

      yddnrm = sqrt( yddnrm / double(neq) );

      //std::cout << "iter " << _iter << " hg " << hg << " y'' " << yddnrm << std::endl;
      //std::cout << "ydot " << ydot[neq-1] << std::endl;

      // should we accept this?
      if (hnew_is_ok || iter == max_iters){
         hnew = hg;
         if (iter == max_iters)
            fprintf(stderr, "ERROR_HIN_MAX_ITERS\n");
         break;
      }

      // Get the new value of h ...
      hnew = (yddnrm*hmax*hmax > 2.0) ? sqrt(2.0 / yddnrm) : sqrt(hg * hmax);

      // test the stopping conditions.
      double hrat = hnew / hg;

      // Accept this value ... the bias factor should bring it within range.
      if ( (hrat > 0.5) && (hrat < 2.0) )
         hnew_is_ok = true;

      // If y'' is still bad after a few iterations, just accept h and give up.
      if ( (iter > 1) && hrat > 2.0 ) {
         hnew = hg;
         hnew_is_ok = true;
      }

      //printf("iter=%d, yddnrw=%e, hnew=%e, hmin=%e, hmax=%e\n", iter, yddnrm, hnew, hmin, hmax);

      hg = hnew;
      iter ++;
   }

   // bound and bias estimate
   h0 = hnew * 0.5;
   h0 = fmax(h0, hmin);
   h0 = fmin(h0, hmax);
   //printf("h0=%e, hmin=%e, hmax=%e\n", h0, hmin, hmax);

   return (iter + 1);
}

void FixRX::odeDiagnostics(void)
{
  TimerType timer_start = getTimeStamp();

  // Compute:
  // 1) Average # of ODE integrator steps and RHS evaluations per atom globally.
  // 2) RMS     # of  ...
  // 3) Average # of ODE steps and RHS evaluations per MPI task.
  // 4) RMS     # of ODE steps and RHS evaluations per MPI task.
  // 5) MAX     # of ODE steps and RHS evaluations per MPI task.
  //
  // ... 1,2 are for ODE control diagnostics.
  // ... 3-5 are for load balancing diagnostics.
  //
  // To do this, we'll need to
  // a) Allreduce (sum) the sum of nSteps / nFuncs. Dividing by atom->natoms
  //    gives the avg # of steps/funcs per atom globally.
  // b) Reduce (sum) to root the sum of squares of the differences.
  //    i) Sum_i (steps_i - avg_steps_global)^2
  //   ii) Sum_i (funcs_i - avg_funcs_global)^2
  //  iii) (avg_steps_local - avg_steps_global)^2
  //   iv) (avg_funcs_local - avg_funcs_global)^2

  const int numCounters = numDiagnosticCounters-1;

  // # of time-steps for averaging.
  const int nTimes = this->diagnosticCounter[numDiagnosticCounters-1];

  // # of ODE's per time-step (on average).
  //const int nODEs  = this->diagnosticCounter[AtomSum] / nTimes;

  // Sum up the sums from each task.
  double sums[numCounters];
  double my_vals[numCounters];
  double max_per_proc[numCounters];
  double min_per_proc[numCounters];

  if(1)
  {
     static bool firstStep = true;

     static TimerType oldTimeStamp (-1);

     TimerType now = getTimeStamp();

     // Query the fix database and look for rx_weight for the balance fix.
     int type_flag = -1;
     int rx_weight_index = atom->find_custom( "rx_weight", /*0:int, 1:float*/ type_flag );

     // Compute the average # of neighbors.
     double averageNumNeighbors = 0;
     {
        const int inum = pairDPDE->list->inum;
        const int* ilist = pairDPDE->list->ilist;
        const int* numneigh = pairDPDE->list->numneigh;

        for (int ii = 0; ii < inum; ++ii)
        {
           const int i = ilist[ii];
           averageNumNeighbors += numneigh[i];
        }

        averageNumNeighbors /= inum;
     }

     printf("me= %d nst= %g nfc= %g time= %g nlocal= %g lmpnst= %g weight_idx= %d 1st= %d aveNeigh= %g\n", comm->me, this->diagnosticCounter[0], this->diagnosticCounter[1], this->diagnosticCounter[2], this->diagnosticCounter[3], this->diagnosticCounter[4], rx_weight_index, firstStep, averageNumNeighbors);

     if (rx_weight_index != -1 && !firstStep && 0)
     {
        double *rx_weight = atom->dvector[rx_weight_index];

        const int nlocal = atom->nlocal;
        const int *mask = atom->mask;

        if (odeIntegrationFlag == ODE_LAMMPS_RKF45 && diagnosticFrequency == 1)
        {
          const double total_time = getElapsedTime( oldTimeStamp, now );
          const double fixrx_time = this->diagnosticCounter[TimeSum];
          const double time_ratio = fixrx_time / total_time;

          double tsum = 0.0;
          double tmin = 100000, tmax = 0;
          for (int i = 0; i < nlocal; ++i)
            if (mask[i] & groupbit)
            {
              double nfunc_ratio = double( diagnosticCounterPerODE[FuncSum][i] ) / diagnosticCounter[FuncSum];
              rx_weight[i] = nfunc_ratio * fixrx_time + (total_time - fixrx_time) / nlocal;
              tmin = fmin( tmin, rx_weight[i] );
              tmax = fmax( tmax, rx_weight[i] );
              tsum += rx_weight[i];
              //rx_weight[i] = (double) diagnosticCounterPerODE[FuncSum][i];
            }

          printf("me= %d total= %g fixrx= %g ratio= %g tsum= %g %g %g %g\n", comm->me, total_time, fixrx_time, time_ratio, tsum, (total_time - fixrx_time) / nlocal, tmin, tmax);
        }
        else
        {
          error->warning(FLERR, "Dynamic load balancing enabled but per-atom weights not available.");

          for (int i = 0; i < nlocal; ++i)
            if (mask[i] & groupbit)
              rx_weight[i] = 1.0;
        }
     }

     firstStep = false;
     oldTimeStamp = now;
  }

  // Compute counters per dpd time-step.
  for (int i = 0; i < numCounters; ++i){
    my_vals[i] = this->diagnosticCounter[i] / nTimes;
    //printf("my sum[%d] = %f %d\n", i, my_vals[i], comm->me);
  }

  MPI_Allreduce (my_vals, sums, numCounters, MPI_DOUBLE, MPI_SUM, world);

  MPI_Reduce (my_vals, max_per_proc, numCounters, MPI_DOUBLE, MPI_MAX, 0, world);
  MPI_Reduce (my_vals, min_per_proc, numCounters, MPI_DOUBLE, MPI_MIN, 0, world);

  const double nODEs = sums[numCounters-1];

  double avg_per_atom[numCounters], avg_per_proc[numCounters];

  // Averages per-ODE and per-proc per time-step.
  for (int i = 0; i < numCounters; ++i){
    avg_per_atom[i] = sums[i] / nODEs;
    avg_per_proc[i] = sums[i] / comm->nprocs;
  }

  // Sum up the differences from each task.
  double sum_sq[2*numCounters];
  double my_sum_sq[2*numCounters];
  for (int i = 0; i < numCounters; ++i){
    double diff_i = my_vals[i] - avg_per_proc[i];
    my_sum_sq[i] = diff_i * diff_i;
  }

  double max_per_ODE[numCounters], min_per_ODE[numCounters];

  // Process the per-ODE RMS of the # of steps/funcs
  if (diagnosticFrequency == 1){
    double my_max[numCounters], my_min[numCounters];

    const int nlocal = atom->nlocal;
    const int *mask = atom->mask;

    for (int i = 0; i < numCounters; ++i){
      my_sum_sq[i+numCounters] = 0;
      my_max[i] = 0;
      my_min[i] = DBL_MAX;

      if (diagnosticCounterPerODE[i] != NULL){
        for (int j = 0; j < nlocal; ++j)
          if (mask[j] & groupbit){
            double diff = double(diagnosticCounterPerODE[i][j]) - avg_per_atom[i];
            my_sum_sq[i+numCounters] += diff*diff;

            my_max[i] = std::max( my_max[i], (double)diagnosticCounterPerODE[i][j] );
            my_min[i] = std::min( my_min[i], (double)diagnosticCounterPerODE[i][j] );
          }
      }
    }

    MPI_Reduce (my_sum_sq, sum_sq, 2*numCounters, MPI_DOUBLE, MPI_SUM, 0, world);

    MPI_Reduce (my_max, max_per_ODE, numCounters, MPI_DOUBLE, MPI_MAX, 0, world);
    MPI_Reduce (my_min, min_per_ODE, numCounters, MPI_DOUBLE, MPI_MIN, 0, world);
  }
  else
    MPI_Reduce (my_sum_sq, sum_sq, numCounters, MPI_DOUBLE, MPI_SUM, 0, world);

  TimerType timer_stop = getTimeStamp();
  double time_local = getElapsedTime( timer_start, timer_stop );

  if (comm->me == 0){
    char smesg[128];

#define print_mesg(smesg) {\
    if (screen)  fprintf(screen,"%s\n", smesg); \
    if (logfile) fprintf(logfile,"%s\n", smesg); }

    sprintf(smesg, "FixRX::ODE Diagnostics:  # of iters  |# of rhs evals| run-time (sec) | # atoms");
    print_mesg(smesg);

    sprintf(smesg, "         AVG per ODE  : %-12.5g | %-12.5g | %-12.5g", avg_per_atom[0], avg_per_atom[1], avg_per_atom[2]);
    print_mesg(smesg);

    // only valid for single time-step!
    if (diagnosticFrequency == 1){
      double rms_per_ODE[numCounters];
      for (int i = 0; i < numCounters; ++i)
        rms_per_ODE[i] = sqrt( sum_sq[i+numCounters] / nODEs );

      sprintf(smesg, "         RMS per ODE  : %-12.5g | %-12.5g ", rms_per_ODE[0], rms_per_ODE[1]);
      print_mesg(smesg);

      sprintf(smesg, "         MAX per ODE  : %-12.5g | %-12.5g ", max_per_ODE[0], max_per_ODE[1]);
      print_mesg(smesg);

      sprintf(smesg, "         MIN per ODE  : %-12.5g | %-12.5g ", min_per_ODE[0], min_per_ODE[1]);
      print_mesg(smesg);
    }

    sprintf(smesg, "         AVG per Proc : %-12.5g | %-12.5g | %-12.5g | %-12.5g", avg_per_proc[StepSum], avg_per_proc[FuncSum], avg_per_proc[TimeSum], avg_per_proc[AtomSum]);
    print_mesg(smesg);

    if (comm->nprocs > 1){
      double rms_per_proc[numCounters];
      for (int i = 0; i < numCounters; ++i)
        rms_per_proc[i] = sqrt( sum_sq[i] / comm->nprocs );

      sprintf(smesg, "         RMS per Proc : %-12.5g | %-12.5g | %-12.5g | %-12.5g", rms_per_proc[0], rms_per_proc[1], rms_per_proc[2], rms_per_proc[AtomSum]);
      print_mesg(smesg);

      sprintf(smesg, "         MAX per Proc : %-12.5g | %-12.5g | %-12.5g | %-12.5g", max_per_proc[0], max_per_proc[1], max_per_proc[2], max_per_proc[AtomSum]);
      print_mesg(smesg);

      sprintf(smesg, "         MIN per Proc : %-12.5g | %-12.5g | %-12.5g | %-12.5g", min_per_proc[0], min_per_proc[1], min_per_proc[2], min_per_proc[AtomSum]);
      print_mesg(smesg);
    }

    sprintf(smesg, "  AVG'd over %d time-steps", nTimes);
    print_mesg(smesg);
    sprintf(smesg, "  AVG'ing took %g sec", time_local);
    print_mesg(smesg);

#undef print_mesg

  }

  // Reset the counters.
  for (int i = 0; i < numDiagnosticCounters; ++i)
    diagnosticCounter[i] = 0;

  return;
}

void FixRX::rkf45(int id, double *rwork, void *v_param, int ode_counter[])
{
  // Rounding coefficient.
  const double uround = DBL_EPSILON;

  // Adaption limit (shrink or grow)
  const double adaption_limit = 4.0;

  //double *y = new double[8*nspecies + nreactions];
  double *y = rwork;

  const int neq = nspecies;

  // Update ConcOld and initialize the ODE solution vector y[].
  for (int ispecies = 0; ispecies < nspecies; ispecies++){
    const double tmp = atom->dvector[ispecies][id];
    atom->dvector[ispecies+nspecies][id] = tmp;
    y[ispecies] = tmp;
  }

  // Integration length.
  const double t_stop = update->dt; // DPD time-step.

  // Safety factor on the adaption. very specific but not necessary .. 0.9 is common.
  const double hsafe = 0.840896415;

  // Time rounding factor.
  const double tround = t_stop * uround;

  // Counters for diagnostics.
  int nst = 0; // # of steps (accepted)
  int nit = 0; // # of iterations total
  int nfe = 0; // # of RHS evaluations

  // Min/Max step-size limits.
  const double h_min = 100.0 * tround;
  const double h_max = (minSteps > 0) ? t_stop / double(minSteps) : t_stop;

  // Set the initial step-size. 0 forces an internal estimate ... stable Euler step size.
  double h = (minSteps > 0) ? t_stop / double(minSteps) : 0.0;

  double t = 0.0;

  if (h < h_min){
    //fprintf(stderr,"hin not implemented yet\n");
    //exit(-1);
    nfe = rkf45_h0 (neq, t, t_stop, h_min, h_max, h, y, y + neq, v_param);
  }

  //printf("t= %e t_stop= %e h= %e\n", t, t_stop, h);

  // Integrate until we reach the end time.
  while (fabs(t - t_stop) > tround){
    double *yout = y + neq;
    double *eout = yout + neq;

    // Take a trial step.
    rkf45_step (neq, h, y, yout, eout, v_param);

    // Estimate the solution error.
      // ... weighted 2-norm of the error.
      double err2 = 0.0;
      for (int k = 0; k < neq; k++){
        const double wterr = eout[k] / (relTol * fabs( y[k] ) + absTol);
        err2 += wterr * wterr;
      }

    double err = fmax( uround, sqrt( err2 / double(nspecies) ));

    // Accept the solution?
    if (err <= 1.0 || h <= h_min){
      t += h;
      nst++;

      for (int k = 0; k < neq; k++)
        y[k] = yout[k];
    }

    // Adjust h for the next step.
    double hfac = hsafe * sqrt( sqrt( 1.0 / err ) );

    // Limit the adaption.
    hfac = fmax( hfac, 1.0 / adaption_limit );
    hfac = fmin( hfac,       adaption_limit );

    // Apply the adaption factor...
    h *= hfac;

    // Limit h.
    h = fmin( h, h_max );
    h = fmax( h, h_min );

    // Stretch h if we're within 5% ... and we didn't just fail.
    if (err <= 1.0 && (t + 1.05*h) > t_stop)
      h = t_stop - t;

    // And don't overshoot the end.
    if (t + h > t_stop)
      h = t_stop - t;

    nit++;
    nfe += 6;

    if (maxIters && nit > maxIters){
      //fprintf(stderr,"atom[%d] took too many iterations in rkf45 %d %e %e\n", id, nit, t, t_stop);
      //nFails ++;
      ode_counter[3] ++;
      break;
      // We should set an error here so that the solution is not used!
    }

  } // end while

  ode_counter[0] += nst;
  ode_counter[1] += nit;
  ode_counter[2] += nfe;

  //if (diagnosticFrequency == 1 && diagnosticCounterPerODE[StepSum] != NULL)
  if (diagnosticCounterPerODE[StepSum] != NULL){
    diagnosticCounterPerODE[StepSum][id] = nst;
    diagnosticCounterPerODE[FuncSum][id] = nfe;
  }
  //printf("id= %d nst= %d nit= %d\n", id, nst, nit);

  // Store the solution back in atom->dvector.
  for (int ispecies = 0; ispecies < nspecies; ispecies++){
    if(y[ispecies] < -1.0e-10)
      error->one(FLERR,"Computed concentration in RKF45 solver is < -1.0e-10");
    else if(y[ispecies] < MY_EPSILON)
      y[ispecies] = 0.0;
    atom->dvector[ispecies][id] = y[ispecies];
  }
}

/* ---------------------------------------------------------------------- */

int FixRX::rhs(double t, const double *y, double *dydt, void *params)
{
  // Use the sparse format instead.
  if (useSparseKinetics)
    return this->rhs_sparse( t, y, dydt, params);
  else
    return this->rhs_dense ( t, y, dydt, params);
}

/* ---------------------------------------------------------------------- */

int FixRX::rhs_dense(double /*t*/, const double *y, double *dydt, void *params)
{
  UserRHSData *userData = (UserRHSData *) params;

  double *rxnRateLaw = userData->rxnRateLaw;
  double *kFor       = userData->kFor;

  const double VDPD = domain->xprd * domain->yprd * domain->zprd / atom->natoms;
  const int nspecies = atom->nspecies_dpd;

  for(int ispecies=0; ispecies<nspecies; ispecies++)
    dydt[ispecies] = 0.0;

  // Construct the reaction rate laws
  for(int jrxn=0; jrxn<nreactions; jrxn++){
    double rxnRateLawForward = kFor[jrxn];

    for(int ispecies=0; ispecies<nspecies; ispecies++){
      const double concentration = y[ispecies]/VDPD;
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

int FixRX::rhs_sparse(double /*t*/, const double *y, double *dydt, void *v_params) const
{
   UserRHSData *userData = (UserRHSData *) v_params;

   const double VDPD = domain->xprd * domain->yprd * domain->zprd / atom->natoms;

   #define kFor         (userData->kFor)
   #define kRev         (NULL)
   #define rxnRateLaw   (userData->rxnRateLaw)
   #define conc         (dydt)
   #define maxReactants (this->sparseKinetics_maxReactants)
   #define maxSpecies   (this->sparseKinetics_maxSpecies)
   #define nuk          (this->sparseKinetics_nuk)
   #define nu           (this->sparseKinetics_nu)
   #define inu          (this->sparseKinetics_inu)
   #define isIntegral(idx) (SparseKinetics_enableIntegralReactions \
                             && this->sparseKinetics_isIntegralReaction[idx])

   for (int k = 0; k < nspecies; ++k)
      conc[k] = y[k] / VDPD;

   // Construct the reaction rate laws
   for (int i = 0; i < nreactions; ++i)
   {
      double rxnRateLawForward;
      if (isIntegral(i)){
         rxnRateLawForward = kFor[i] * powint( conc[ nuk[i][0] ], inu[i][0]);
         for (int kk = 1; kk < maxReactants; ++kk){
            const int k = nuk[i][kk];
            if (k == SparseKinetics_invalidIndex) break;
            //if (k != SparseKinetics_invalidIndex)
               rxnRateLawForward *= powint( conc[k], inu[i][kk] );
         }
      } else {
         rxnRateLawForward = kFor[i] * pow( conc[ nuk[i][0] ], nu[i][0]);
         for (int kk = 1; kk < maxReactants; ++kk){
            const int k = nuk[i][kk];
            if (k == SparseKinetics_invalidIndex) break;
            //if (k != SparseKinetics_invalidIndex)
               rxnRateLawForward *= pow( conc[k], nu[i][kk] );
         }
      }

      rxnRateLaw[i] = rxnRateLawForward;
   }

   // Construct the reaction rates for each species from the
   // Stoichiometric matrix and ROP vector.
   for (int k = 0; k < nspecies; ++k)
      dydt[k] = 0.0;

   for (int i = 0; i < nreactions; ++i){
      // Reactants ...
      dydt[ nuk[i][0] ] -= nu[i][0] * rxnRateLaw[i];
      for (int kk = 1; kk < maxReactants; ++kk){
         const int k = nuk[i][kk];
         if (k == SparseKinetics_invalidIndex) break;
         //if (k != SparseKinetics_invalidIndex)
            dydt[k] -= nu[i][kk] * rxnRateLaw[i];
      }

      // Products ...
      dydt[ nuk[i][maxReactants] ] += nu[i][maxReactants] * rxnRateLaw[i];
      for (int kk = maxReactants+1; kk < maxSpecies; ++kk){
         const int k = nuk[i][kk];
         if (k == SparseKinetics_invalidIndex) break;
         //if (k != SparseKinetics_invalidIndex)
            dydt[k] += nu[i][kk] * rxnRateLaw[i];
      }
   }

   // Add in the volume factor to convert to the proper units.
   for (int k = 0; k < nspecies; ++k)
      dydt[k] *= VDPD;

   #undef kFor
   #undef kRev
   #undef rxnRateLaw
   #undef conc
   #undef maxReactants
   #undef maxSpecies
   #undef nuk
   #undef nu
   #undef inu
   #undef isIntegral
   //#undef invalidIndex

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
  double wij=0.0;
  double *dpdTheta = atom->dpdTheta;

  // Initialize the local temperature weight array
  int sumWeightsCt = nlocal + (newton_pair ? nghost : 0);
  sumWeights = new double[sumWeightsCt];
  memset(sumWeights, 0, sizeof(double)*sumWeightsCt);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

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
          wij = (1.0+3.0*ratio) * (1.0-ratio)*(1.0-ratio)*(1.0-ratio);
          dpdThetaLocal[i] += wij/dpdTheta[j];
          if (newton_pair || j < nlocal)
            dpdThetaLocal[j] += wij/dpdTheta[i];
        }

        sumWeights[i] += wij;
        if (newton_pair || j < nlocal)
          sumWeights[j] += wij;
      }
    }
  }
  if (newton_pair) comm->reverse_comm_fix(this);

  // self-interaction for local temperature
  for (i = 0; i < nlocal; i++){

    // Lucy Weight Function
    if(wtFlag==LUCY){
      wij = 1.0;
      dpdThetaLocal[i] += wij / dpdTheta[i];
    }
    sumWeights[i] += wij;

    // Normalized local temperature
    dpdThetaLocal[i] = dpdThetaLocal[i] / sumWeights[i];

    if(localTempFlag == HARMONIC)
      dpdThetaLocal[i] = 1.0 / dpdThetaLocal[i];

  }

  delete [] sumWeights;
}

/* ---------------------------------------------------------------------- */

int FixRX::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
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
