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
   The SMTBQ code has been developed with the financial support of  CNRS and
   of the Regional Council of Burgundy (Convention n¡ 2010-9201AAO037S03129)

   Copyright (2015)
   Universite de Bourgogne : Nicolas SALLES, Olivier POLITANO
   Universite de Paris-Sud Orsay : R. Tetot
   Aalto University (Finland) : E. Maras

   Please cite the related publication:
   N. Salles, O. Politano, E. Amzallag and R. Tetot,
   Comput. Mater. Sci., 111 (2016) 181-189

   Contact : lammps@u-bourgogne.fr

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of
   the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   See the GNU General Public License for more details:
   <http://www.gnu.org/licenses/>.
   ------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_smtbq.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "group.h"
#include "update.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "error.h"
#include "domain.h"

#include <fstream>
#include <iomanip>

using namespace std;

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define MAXLINE 2048
#define MAXTOKENS 2048
#define DELTA 4
#define PGDELTA 1
#define MAXNEIGH 24

/* ------------------------------------------------------------------------------------

   Calculates the factorial of an integer n via recursion.

   ------------------------------------------------------------------------------------ */
static double factorial(int n)
{
  if (n <= 1) return 1.0;
  else return static_cast<double>(n)*factorial(n-1);
}

/* ---------------------------------------------------------------------- */

PairSMTBQ::PairSMTBQ(LAMMPS *lmp) : Pair(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nproc);

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;

  nmax = 0;
  rmin = 0.0;
  dr = 0.0;
  ds = 0.0;
  kmax = 0;

  nelements = 0;
  elements = NULL;
  nparams = 0;
  maxparam = 0;
  params = NULL;
  intparams = NULL;

  intype = NULL;
  coultype = NULL;
  fafb = NULL;
  dfafb = NULL;
  potqn = NULL;
  dpotqn = NULL;
  Vself = 0.0;
  tabsmb = NULL;
  tabsmr = NULL;
  dtabsmb = NULL;
  dtabsmr = NULL;

  sbcov = NULL;
  coord = NULL;
  sbmet = NULL;
  chimet = NULL;
  ecov = NULL;

  potmad = NULL;
  potself = NULL;
  potcov = NULL;
  qf = NULL;
  q1 = NULL;
  q2 = NULL;
  tab_comm = NULL;

  nvsm = NULL;
  vsm = NULL;
  flag_QEq = NULL;
  nQEqaall = NULL;
  nQEqcall = NULL;
  nQEqall = NULL;
  nteam = 0;
  cluster = 0;

  Nevery = 0.0;
  Neverypot = 0.0;

  fct = NULL;

  maxpage = 0;

  // set comm size needed by this Pair

  comm_forward = 1;
  comm_reverse = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
   ------------------------------------------------------------------------- */

PairSMTBQ::~PairSMTBQ()
{
  int i;
  if (elements) {
    for ( i = 0; i < nelements; i++) delete [] elements[i];

    for( i = 0; i < atom->ntypes ; i++ ) free( params[i].nom );
    for( i = 1; i <= maxintparam ; i++ ) free( intparams[i].typepot );
    for( i = 1; i <= maxintparam ; i++ ) free( intparams[i].mode );
  }

  free(QEqMode);
  free(QInitMode);
  free(writepot);
  free(writeenerg);
  free(Bavard);

  delete [] elements;
  memory->sfree(params);
  memory->sfree(intparams);

  memory->destroy(intype);
  memory->destroy(coultype);
  memory->destroy(fafb);
  memory->destroy(dfafb);
  memory->destroy(potqn);
  memory->destroy(dpotqn);

  memory->destroy(fafbOxOxSurf);
  memory->destroy(dfafbOxOxSurf);
  memory->destroy(fafbTiOxSurf);
  memory->destroy(dfafbTiOxSurf);

  memory->destroy(fafbOxOxBB);
  memory->destroy(dfafbOxOxBB);
  memory->destroy(fafbTiOxBB);
  memory->destroy(dfafbTiOxBB);

  memory->destroy(ecov);
  memory->destroy(sbcov);
  memory->destroy(coord);
  memory->destroy(sbmet);
  memory->destroy(tabsmb);
  memory->destroy(tabsmr);
  memory->destroy(dtabsmb);
  memory->destroy(dtabsmr);

  memory->destroy(potmad);
  memory->destroy(potself);
  memory->destroy(potcov);
  memory->destroy(chimet);

  memory->destroy(nvsm);
  memory->destroy(vsm);;
  memory->destroy(flag_QEq);

  memory->destroy(nQEqall);
  memory->destroy(nQEqcall);
  memory->destroy(nQEqaall);

  memory->destroy(qf);
  memory->destroy(q1);
  memory->destroy(q2);
  memory->destroy(tab_comm);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
    delete [] esm;
  }

  memory->destroy(fct);
}

/* ---------------------------------------------------------------------- */

void PairSMTBQ::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
  esm = new double[n];
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairSMTBQ::settings(int narg, char **/*arg*/)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

void PairSMTBQ::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (strstr(force->pair_style,"hybrid"))
    error->all(FLERR,"Pair style SMTBQ is not compatible with hybrid styles");

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

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

  // read potential file and initialize potential parameters

  read_file(arg[2]);

  n = atom->ntypes;

  // generate Coulomb 1/r energy look-up table

  if (comm->me == 0 && screen) fprintf(screen,"Pair SMTBQ:\n");
  if (comm->me == 0 && screen)
    fprintf(screen,"  generating Coulomb integral lookup table ...\n");

  tabqeq();

  // ------------


  if (comm->me == 0 && screen)
    fprintf(screen,"  generating Second Moment integral lookup table ...\n");

  tabsm();

  // ------------

  // clear setflag since coeff() called once with I,J = * *

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

void PairSMTBQ::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style SMTBQ requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style SMTBQ requires newton pair on");
  if (!atom->q_flag)
    error->all(FLERR,"Pair style SMTBQ requires atom attribute q");


  // need a full neighbor list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  pgsize = neighbor->pgsize;
  oneatom = neighbor->oneatom;
  //  if (maxpage == 0) add_pages();

}

/* ---------------------------------------------------------------------- */

double PairSMTBQ::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  return cutmax;
}

/* ----------------------------------------------------------------------
   ---------------------------------------------------------------------- */

void PairSMTBQ::read_file(char *file)
{
  int num_atom_types,i,k,m,test,j,verbose;
  char **words;

  memory->sfree(params);
  params = NULL;
  memory->sfree(intparams);
  intparams = NULL;
  nparams = 0;
  maxparam = 0;
  maxintparam = 0;

  verbose = 1;
  verbose = 0;

  // open file on all processors
  FILE *fp;
  fp = force->open_potential(file);
  if ( fp  == NULL ) {
    char str[128];
    snprintf(str,128,"Cannot open SMTBQ potential file %s",file);
    error->one(FLERR,str);
  }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  char *ptr;

  ptr = (char*) malloc(sizeof(char)*MAXLINE);
  words = (char**) malloc(sizeof(char*)*MAXTOKENS);
  for (i=0; i < MAXTOKENS; i++)
    words[i] = (char*) malloc(sizeof(char)*MAXTOKENS);


  /* strip comment, skip line if blank */

  if (verbose) printf ("\n");
  fgets(ptr,MAXLINE,fp);
  while (strchr(ptr,'#')) {
    if (verbose) printf ("%s",ptr);
    fgets(ptr,MAXLINE,fp);
  }


  // Nombre d'atome different dans la structure
  //  ===============================================
  Tokenize( ptr, &words );
  num_atom_types = atoi(words[1]);
  if (verbose) printf (" %s %d\n", words[0], num_atom_types);

  memory->create(intype,num_atom_types,num_atom_types,"pair:intype");

  m = 0;
  for (i = 0; i < num_atom_types; i++) {
    for (j = 0; j < num_atom_types; j++) {
      if (j < i) { intype[i][j] = intype[j][i];}
      else       { intype[i][j] = 0;
        m = m + 1;   }
      if (verbose) printf ("i %d, j %d, intype %d - nb pair %d\n",i,j,intype[i][j],m);
    }
  }

  // load up parameter settings and error check their values

  if (nparams == maxparam) {
    maxparam += DELTA;
    params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                        "pair:params");
    maxintparam += m;
    intparams = (Intparam *) memory->srealloc(intparams,(maxintparam+1)*sizeof(Intparam),
                                              "pair:intparams");
  }

  for (i=0; i < num_atom_types; i++)
    params[i].nom = (char*) malloc(sizeof(char)*3);

  for (i=1; i <= maxintparam; i++)
    intparams[i].typepot = (char*) malloc(sizeof(char)*15);

  for (i=1; i <= maxintparam; i++)
    intparams[i].mode = (char*) malloc(sizeof(char)*6);

  QEqMode = (char*) malloc(sizeof(char)*19);
  Bavard = (char*) malloc(sizeof(char)*6);
  QInitMode = (char*) malloc(sizeof(char)*19);
  writepot = (char*) malloc(sizeof(char)*6);
  writeenerg = (char*) malloc(sizeof(char)*6);


  //  Little loop for ion's parameters
  // ================================================
  for (i=0; i<num_atom_types; i++) {

    fgets(ptr,MAXLINE,fp);  if (verbose) printf ("%s",ptr);

    // Line 2 - Al

    fgets( ptr, MAXLINE, fp);
    Tokenize( ptr, &words );
    strcpy(params[i].nom , words[1]);
    params[i].sto = atof(words[2]);
    if (verbose) printf (" %s %s %f\n", words[0],params[i].nom,params[i].sto);

    //Line 3 - Charges

    fgets( ptr, MAXLINE, fp);
    Tokenize( ptr, &words );

    params[i].qform = atof(words[1]);
    params[i].masse = atof(words[2]);
    if (verbose) printf (" %s %f %f \n", words[0],params[i].qform, params[i].masse);

    // Line 4 - Parametres QEq

    fgets( ptr, MAXLINE, fp);
    Tokenize ( ptr, &words );
    params[i].ne = atof(words[1]) ;
    params[i].chi = atof(words[2])  ;
    params[i].dj = atof(words[3]) ;

    if(strcmp(params[i].nom,"O")!=0){
      params[i].R = atof(words[4]) ;
      if (verbose) printf(" %s %f %f %f %f\n",words[0],params[i].ne,params[i].chi,
                          params[i].dj,params[i].R);
    } else {
      if (verbose) printf(" %s %f %f %f\n",words[0],params[i].ne,params[i].chi,params[i].dj);
    }


    // Line 4bis - Coordinance et rayon pour Ox
    if(strcmp(params[i].nom,"O")==0){

      fgets( ptr, MAXLINE, fp);
      Tokenize ( ptr, &words );

      coordOxBB=   atof(words[1]) ;
      coordOxBulk=  atof(words[2]) ;
      coordOxSurf=  atof(words[3]) ;
      ROxBB = atof(words[4]) ;
      params[i].R = atof(words[5]) ;
      ROxSurf = atof(words[6]) ;
      if (verbose) printf(" %s %f %f %f %f %f %f\n",words[0],coordOxBB,coordOxBulk,coordOxSurf,ROxBB,params[i].R,ROxSurf);
    }

    // Ligne 5 - Nombre d'etats partages

    fgets( ptr, MAXLINE, fp);
    Tokenize ( ptr, &words );
    params[i].n0 = atof(words[1]);
    if (verbose) printf(" %s %f\n",words[0],params[i].n0);

    // Parametres de Slater
    params[i].dzeta = (2.0*params[i].ne + 1.0)/(4.0*params[i].R);
    if (verbose) printf (" Parametre dzeta (Slater) : %f\n",params[i].dzeta);

  } // Fin elements i

  /* =====================================================================
     reading the interaction's parameters
     ===================================================================== */

  m = 0; maxintsm = 0;        //
  for (k=0 ; k<=maxintparam ; k++){intparams[k].intsm = 0;}
  //  ---------------------------------
  for (k = 0; k < maxintparam; k++) {
    //  ---------------------------------
    m += 1;

    // Ligne 5 - parametre des potentials
    fgets(ptr,MAXLINE,fp);  if (verbose) printf ("%s",ptr);

    // Lecture des protagonistes
    fgets( ptr, MAXLINE, fp);
    Tokenize( ptr, &words );

    test = 0;
    for (i = 0; i <num_atom_types; i++)
      {
        if (strcmp(params[i].nom,words[1])==0) break;
        if (i == num_atom_types - 1) test = 1;
      }
    //   if (test == 0) printf (" on a %s -> %d = %s\n",words[1],i,params[i].nom);

    for (j = 0; j <num_atom_types; j++)
      {
        if (strcmp(params[j].nom,words[2])==0) break;
        if (j == num_atom_types - 1) test = 1;
      }
    //    if (test == 0) printf (" on a %s -> %d = %s\n",words[2],j,params[j].nom);


    if ( test == 1 ) {
      if (verbose) printf ("========== fin des interaction ==========\n");
      break ; }


    intype[i][j] = m;
    intype[j][i] = intype[i][j];
    strcpy( intparams[m].typepot , words[3] );
    intparams[m].intsm = 0;
    if (verbose) printf (" itype %d jtype %d - intype %d\n",i,j,intype[i][j]);

    if (strcmp(intparams[m].typepot,"second_moment") !=0 &&
        strcmp(intparams[m].typepot,"buck") != 0 &&
        strcmp(intparams[m].typepot,"buckPlusAttr") != 0) {
      error->all(FLERR,"the potential other than second_moment or buckingham have not treated here\n");}


    // On detemrine le type d'interaction
    // -----------------------------------
    if (strcmp(intparams[m].typepot,"second_moment") == 0) {
      maxintsm += 1;
      strcpy( intparams[m].mode , words[4] );
      intparams[m].intsm = maxintsm;

      if (strcmp(intparams[m].mode,"oxide") != 0 &&
          strcmp(intparams[m].mode,"metal") != 0){
        error->all(FLERR,"needs mode to second moment interaction : oxide or metal"); }

      //      if (strcmp(intparams[m].mode,"oxide") == 0)
      //             intparams[m].ncov = min((params[i].sto)*(params[i].n0),(params[j].sto)*(params[j].n0));

      if (verbose) printf(" %s %s %s %s %s \n",words[0],words[1],words[2],
                          intparams[m].typepot,intparams[m].mode);

      fgets( ptr, MAXLINE, fp);
      Tokenize( ptr, &words );

      intparams[m].a = atof(words[1])   ;
      intparams[m].p = atof(words[2])   ;
      intparams[m].ksi = atof(words[3]) ;
      intparams[m].q = atof(words[4])   ;
      if (verbose) printf (" %s %f %f %f %f\n",words[0],
                           intparams[m].a,intparams[m].p,intparams[m].ksi,intparams[m].q);

      // Ligne 6 - rayon de coupure potential SM

      fgets( ptr, MAXLINE, fp);
      Tokenize( ptr, &words );

      intparams[m].dc1 = atof(words[1]) ;
      intparams[m].dc2 = atof(words[2]) ;
      intparams[m].r0 = atof(words[3]) ;


      if (strcmp(intparams[m].mode,"metal") == 0) {
        if (verbose) printf (" %s %f %f %f\n",words[0],
                             intparams[m].dc1,intparams[m].dc2,intparams[m].r0);
      } else {
        if (verbose) printf (" %s %f %f %f\n",words[0],
                             intparams[m].dc1,intparams[m].dc2,intparams[m].r0);
      }


    } else if (strcmp(intparams[m].typepot,"buck") == 0) {

      if (verbose) printf(" %s %s %s %s\n",words[0],words[1],words[2],
                          intparams[m].typepot);

      fgets( ptr, MAXLINE, fp);
      Tokenize( ptr, &words );

      intparams[m].abuck = atof(words[1]) ; intparams[m].rhobuck = atof(words[2]) ;
      if (verbose) printf (" %s %f %f\n",words[0],intparams[m].abuck,intparams[m].rhobuck);

    }

    else if (strcmp(intparams[m].typepot,"buckPlusAttr") == 0) {

      if (verbose) printf(" %s %s %s %s\n",words[0],words[1],words[2],
                          intparams[m].typepot);

      fgets( ptr, MAXLINE, fp);
      Tokenize( ptr, &words );

      intparams[m].abuck = atof(words[1]) ; intparams[m].rhobuck = atof(words[2]) ;
      if (verbose) printf (" %s %f %f\n",words[0],intparams[m].abuck,intparams[m].rhobuck);


      fgets( ptr, MAXLINE, fp);
      Tokenize( ptr, &words );

      intparams[m].aOO = atof(words[1]) ; intparams[m].bOO = atof(words[2]) ;
      intparams[m].r1OO = atof(words[3]) ;intparams[m].r2OO = atof(words[4]) ;
      if (verbose) printf (" %s %f %f %f %f \n",words[0],intparams[m].aOO,
                           intparams[m].bOO,intparams[m].r1OO,intparams[m].r2OO);

    }
    if (verbose) printf (" intsm %d \n",intparams[m].intsm);

  } // for maxintparam


  /* ====================================================================
     tables Parameters
     ==================================================================== */

  // Ligne 9 - rayon de coupure Electrostatique
  if (test == 0) {
    fgets(ptr,MAXLINE,fp);
    if (verbose) printf ("%s\n",ptr);

    fgets( ptr, MAXLINE, fp);
  }
  Tokenize( ptr, &words );

  for (i=0 ; i<num_atom_types; i++) { params[i].cutsq = atof(words[1]); }
  cutmax = atof(words[1]);
  if (verbose) printf (" %s %f\n",words[0],params[0].cutsq);

  // Ligne 9 - parametre pour les tableaux

  fgets( ptr, MAXLINE, fp);
  Tokenize( ptr, &words );

  rmin = atof(words[1]) ; dr = atof(words[2]);
  if (verbose) printf (" %s %f %f\n",words[0],rmin,dr);

  kmax = int(cutmax*cutmax/(2.0*dr*rmin));
  ds = cutmax*cutmax/static_cast<double>(kmax) ;
  if (verbose) printf (" kmax %d et ds %f\n",kmax,ds);

  /* ======================================================== */
  fgets( ptr, MAXLINE, fp);
  if (verbose) printf ("%s",ptr);

  fgets( ptr, MAXLINE, fp);
  Tokenize( ptr, &words );
  Qstep = atoi(words[1]);
  if (verbose) printf (" %s " BIGINT_FORMAT "\n",words[0],Qstep);

  fgets( ptr, MAXLINE, fp);
  Tokenize( ptr, &words );
  loopmax = atoi(words[1]);
  precision = atof(words[2]);
  if (verbose) printf (" %s %d %f\n",words[0],loopmax,precision);

  /* Param de coordination ============================================= */

  fgets( ptr, MAXLINE, fp);
  if (verbose) printf ("%s",ptr);

  fgets( ptr, MAXLINE, fp);
  Tokenize( ptr, &words );
  r1Coord = atof(words[1]);
  r2Coord = atof(words[2]);
  if (verbose) printf (" %s %f %f\n",words[0],r1Coord,r2Coord);

  /* Mode for QInit============================================= */
  fgets( ptr, MAXLINE, fp);
  if (verbose) printf ("%s",ptr);

  fgets( ptr, MAXLINE, fp);
  Tokenize( ptr, &words );
  strcpy( QInitMode , words[1] );
  if (strcmp(QInitMode,"true") == 0) QOxInit= atof(words[2]);
  else QOxInit = 0.0;
  if (verbose) printf (" %s %s %f\n",words[0],QInitMode,QOxInit);


  /* Mode for QEq============================================= */

  fgets( ptr, MAXLINE, fp);
  if (verbose) printf ("%s",ptr);

  fgets( ptr, MAXLINE, fp);
  Tokenize( ptr, &words );
  strcpy( QEqMode , words[1] );
  if (verbose) printf (" %s %s\n",words[0],QEqMode);

  fgets( ptr, MAXLINE, fp);

  if (strcmp(QEqMode,"BulkFromSlab") == 0) {
    Tokenize( ptr, &words );
    zlim1QEq = atof(words[1]);
    zlim2QEq = atof(words[2]);
    if (verbose) printf (" %s %f %f\n",words[0],zlim1QEq,zlim2QEq);

  } else if (strcmp(QEqMode,"Surface") == 0) {
    Tokenize( ptr, &words );
    zlim1QEq = atof(words[1]);
    if (verbose) printf (" %s %f \n",words[0],zlim1QEq);

  } else if (strcmp(QEqMode,"QEqAll") != 0         &&
             strcmp(QEqMode,"QEqAllParallel") != 0 &&
             strcmp(QEqMode,"Surface") != 0 ) {
    error->all(FLERR,"The QEq Mode is not known. QEq mode should be :\n"
               "  Possible QEq  modes    |   parameters\n"
               "  QEqAll                      |   no parameters\n"
               "  QEqAllParallel        |   no parameters\n"
               "  Surface                |   zlim   (QEq only for z>zlim)\n"
               "  BulkFromSlab                |   zlim1  zlim2  (QEq only for zlim1<z<zlim2)\n");
  }

  /* Bavard============================================= */

  fgets( ptr, MAXLINE, fp);
  if (verbose) printf ("%s",ptr);

  fgets( ptr, MAXLINE, fp);
  Tokenize( ptr, &words );
  strcpy( Bavard , words[1] );
  if (verbose) printf (" %s %s\n",words[0],Bavard);

  // ---------------------------------------
  //  Writing the energy component.

  fgets( ptr, MAXLINE, fp);
  Tokenize( ptr, &words );
  strcpy( writeenerg, words[1] );
  if (strcmp (writeenerg,"true") == 0) { Nevery = atof(words[2]); }
  else { Nevery = 0.0; }
  if (verbose) printf (" %s %s %f\n",words[0],writeenerg,Nevery);

  // ---------------------------------------
  //  Writing the chimical electronic potential.

  fgets( ptr, MAXLINE, fp);
  Tokenize( ptr, &words );
  strcpy( writepot, words[1] );
  if (strcmp (writepot,"true") == 0) { Neverypot = atof(words[2]); }
  else { Neverypot = 0.0; }
  if (verbose) printf (" %s %s %f\n",words[0],writepot,Neverypot);


  /* ======================================================== */

  /* deallocate helper storage */
  for( i = 0; i < MAXTOKENS ; i++ ) free( words[i] );
  free( words );
  free( ptr );
  fclose(fp);

  // === Rayon de coupure premier voisins : 1,2*r0
  for (i=0 ; i<num_atom_types ; i++) {
    for (j=0 ; j<=i ; j++) {
      m = intype[i][j];
      if (m == 0) continue;
      if (intparams[m].intsm == 0) continue;

      intparams[m].neig_cut = 1.2*intparams[m].r0;
      if (strcmp(intparams[m].typepot,"second_moment") == 0 )
        if (verbose) printf (" Rc 1er voisin, typepot %s -> %f Ang\n",
                             intparams[m].typepot,intparams[m].neig_cut);
    }
  }

  //A adapter au STO
  ncov = min((params[0].sto)*(params[0].n0),(params[1].sto)*(params[1].n0));

  if (verbose) printf (" Parametre ncov = %f\n",ncov);
  if (verbose) printf (" ********************************************* \n");
}

/* ----------------------------------------------------------------------
 *                           COMPUTE
 ---------------------------------------------------------------------- */

void PairSMTBQ::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,m,gp;
  tagint itag,jtag;
  int itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,iq,jq,Eself;
  double ecovtot,ErepOO,ErepMO,Eion,Ecoh;
  double **tmp,**tmpAll,*nmol;
  double dq,dqcov;

  if (atom->nmax > nmax) {
    memory->destroy(ecov);
    memory->destroy(potmad);
    memory->destroy(potself);
    memory->destroy(potcov);
    memory->destroy(sbcov);
    memory->destroy(coord);
    memory->destroy(sbmet);
    memory->destroy(chimet);
    memory->destroy(flag_QEq);
    memory->destroy(qf);
    memory->destroy(q1);
    memory->destroy(q2);
    memory->destroy(tab_comm);

    nmax = atom->nmax;

    memory->create(ecov,nmax,"pair:ecov");
    memory->create(potmad,nmax,"pair:potmad");
    memory->create(potself,nmax,"pair:potself");
    memory->create(potcov,nmax,"pair:potcov");
    memory->create(sbcov,nmax,"pair:sbcov");
    memory->create(coord,nmax,"pair:coord");
    memory->create(sbmet,nmax,"pair:sbmet");
    memory->create(chimet,nmax,"pair:chimet");
    memory->create(flag_QEq,nmax,"pair:flag_QEq");
    memory->create(qf,nmax,"pair:qf");
    memory->create(q1,nmax,"pair:q1");
    memory->create(q2,nmax,"pair:q2");
    memory->create(tab_comm,nmax,"pair:tab_comm");
  }


  evdwl = ecoul = ecovtot = ErepOO = ErepMO = Eion = 0.0;
  Eself = 0.0;

  if (eflag || vflag) { ev_setup(eflag,vflag); }
  else { evflag = vflag_fdotr = vflag_atom = 0; }

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int newton_pair = force->newton_pair;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  const bigint step = update->ntimestep;
  if (step == 0 || ((Qstep > 0) && (step % Qstep == 0))) Charge();
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // this is necessary to get sbcov or sbmet table in order to caclulate the covalent or metal bonding
  if (Qstep == 0 || (step % Qstep != 0)) QForce_charge(0);


  //   Charges Communication
  //   ----------------------
  forward(q) ; // reverse(q);

  memory->create(nmol,nteam+1,"pair:nmol");
  memory->create(tmp,nteam+1,7,"pair:tmp");
  memory->create(tmpAll,nteam+1,7,"pair:tmpAll");


  for (i=0; i<nteam+1; i++) {
    nmol[i] = static_cast<double>(nQEqall[i]);
    for (j=0; j<7; j++) { tmp[i][j] = 0.0; tmpAll[i][j] = 0.0; }
  }


  /* ------------------------------------------------------------------------
     Energy component store in tmp[gp][:] with gp is # QEq group
     0 -> ionic energy
     1 -> coulombian energy
     2 -> Electrosatic energy (ionic + Coulombian)
     3 -> Short int. Ox-Ox
     4 -> Short int. SMTB (repulsion)
     5 -> Covalent energy SMTB
     6 -> Somme des Q(i)²
     ------------------------------------------------------------------------- */

  /* -------------- N-body forces Calcul --------------- */

  for (ii = 0; ii < inum; ii++) {
    //  ===============================
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    iq = q[i];
    gp = flag_QEq[i];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // --- For atom i
    tmp[gp][6] += iq*iq;


    // self energy, only on i atom
    // ---------------------------
    Eself = self(&params[itype],iq);
    tmp[gp][0] += Eself;
    tmp[gp][2] += Eself;

    if (evflag) ev_tally_full (i,0.0,2.0*Eself,0.0,0.0,0.0,0.0);


    //     N-body energy of i
    //    ---------------------
    dq = fabs(params[itype].qform) - fabs(iq);
    dqcov = dq*(2.0*ncov/params[itype].sto - dq);

    ecov[i] = - sqrt(sbcov[i]*dqcov + sbmet[i]);
    ecovtot += ecov[i];
    tmp[gp][5] += ecov[i];

    if (evflag) ev_tally_full(i,0.0,2.0*ecov[i],0.0,0.0,0.0,0.0);


    //   Coulombian Interaction
    //   -----------------------
    evdwl = 2.0*Vself*iq*iq ;
    tmp[gp][1] += Vself*iq*iq;
    tmp[gp][2] += Vself*iq*iq;

    if (evflag) ev_tally_full (i,0.0,evdwl,0.0,0.0,0.0,0.0);
    evdwl = 0.0 ;


    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      //  ===============================
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      jtag = tag[j]; jq = q[j];


      //   .......................................................................
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < x[i][2]) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }
      //   .......................................................................


      // # of interaction
      // ----------------
      m = intype[itype][jtype];


      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;


      //    ---------------------------------
      if (sqrt(rsq) > cutmax) continue;
      //    ---------------------------------


      // Coulombian Energy
      // ------------------
      evdwl = 0.0 ; fpair = 0.0;
      potqeq(i,j,iq,jq,rsq,fpair,eflag,evdwl);

      tmp[gp][1] += evdwl;
      tmp[gp][2] += evdwl;


      // Coulombian Force
      // -----------------
      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;


      if (evflag)
        ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
      evdwl = 0.0; fpair = 0.0 ;



      //    ---------------------
      if (m == 0) continue;
      //    ---------------------

      //    ----------------------------------------------
      if ( strcmp(intparams[m].typepot,"buck") == 0 ||
           strcmp(intparams[m].typepot,"buckPlusAttr") ==0 ) {
        //    ----------------------------------------------

        evdwl = 0.0; fpair =0.0;
        rep_OO (&intparams[m],rsq,fpair,eflag,evdwl);
        ErepOO += evdwl ;
        tmp[gp][3] += evdwl;


        // repulsion is pure two-body, sums up pair repulsive forces
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        f[j][0] -= delx*fpair;
        f[j][1] -= dely*fpair;
        f[j][2] -= delz*fpair;


        if (evflag)
          ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
        evdwl = 0.0; fpair = 0.0 ;

      }  // ----------------------------------- Rep O-O

      if (strcmp(intparams[m].typepot,"buckPlusAttr") == 0 ) {
        //    ----------------------------------------------

        evdwl = 0.0; fpair =0.0;
        Attr_OO (&intparams[m],rsq,fpair,eflag,evdwl);
        ErepOO += evdwl ;
        tmp[gp][3] += evdwl;


        // repulsion is pure two-body, sums up pair repulsive forces
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        f[j][0] -= delx*fpair;
        f[j][1] -= dely*fpair;
        f[j][2] -= delz*fpair;


        if (evflag)
          ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
        evdwl = 0.0; fpair = 0.0 ;

      }  // ----------------------------------- Attr O-O


      //    -----------------------------------------------------------------
      if (strcmp(intparams[m].typepot,"second_moment") != 0 ) continue;
      //    -----------------------------------------------------------------


      if (sqrt(rsq) > intparams[m].dc2) continue;
      //    -------------------------------------------

      //   Repulsion : Energy + force
      //   ----------------------------
      evdwl = 0.0; fpair = 0.0 ;
      repulsive(&intparams[m],rsq,i,j,fpair,eflag,evdwl);
      ErepMO += evdwl;
      tmp[gp][4] += 2.0*evdwl;

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;


      if (evflag)
        ev_tally(i,j,nlocal,newton_pair,2.0*evdwl,0.0,fpair,delx,dely,delz);

      evdwl = 0.0 ; fpair = 0.0;
      //     ----- ----- ----- ----- ----- -----

      //    Attraction : force
      //    ------------------
      fpair = 0.0;
      f_att(&intparams[m], i, j, rsq, fpair) ;

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag)
        ev_tally(i,j,nlocal,newton_pair,0.0,0.0,fpair,delx,dely,delz);


    } // --------------------------------- End j

  } // ---------------------------------- End i


  if (vflag_fdotr) virial_fdotr_compute();


  for (i = 0; i < nteam+1; i++) {
    MPI_Allreduce(tmp[i],tmpAll[i],7,MPI_DOUBLE,MPI_SUM,world);
  }

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  if (me == 0 && fmod(static_cast<double>(step), Nevery) == 0.0 && strcmp(writeenerg,"true") == 0) {

    ofstream fichierE;

    if (step == 0) { fichierE.open ("Energy_component.txt", ios::out | ios::trunc) ;}
    else {           fichierE.open ("Energy_component.txt", ios::out | ios::app) ;}

    if (fichierE) fichierE<< setprecision(9) <<step;

    for (gp = 0; gp < nteam+1; gp++) {
      if (nmol[gp] == 0) continue;
      if (fichierE) fichierE<< setprecision(9) <<" "<<gp<<" "<<nmol[gp]
                            <<" "<<tmpAll[gp][2]<<" "<<tmpAll[gp][3]<<" "<<tmpAll[gp][4]+tmpAll[gp][5];
    }
    if (fichierE) fichierE<<endl;
    if (fichierE) fichierE.close();
  }
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  if (me == 0&& strcmp(Bavard,"false") != 0) {
    printf ("A la fin de Compute\n");

    printf ("Nemton_pair : %d, evflag %d, tail_flag %d,vflag_fdotr %d\n",
            newton_pair,evflag,force->pair->tail_flag,vflag_fdotr);
    printf ("neighbor->includegroup %d\n",neighbor->includegroup);



    for (gp=0; gp<nteam+1; gp++) { // -------------------------- Boucle sur les gp

      printf ("Energie de la team %d -- %d atome -------(nteam = %d)  \n ",gp,int(nmol[gp]),nteam);

      if (nmol[gp] == 0) {
        printf (" ============================================ \n \n");
        continue;
      }
      printf ("Vself %f, Som q2 : %f, nmol %f\n",Vself,tmpAll[gp][6],nmol[gp]);
      //   printf ("nmol %f\n",nmol[gp]);
      printf ("Energie coul tot   : %f | %f par mol\n",tmpAll[gp][1],tmpAll[gp][1]/nmol[gp]);
      printf ("Energie ionique    : %f | %f par mol\n",tmpAll[gp][0],tmpAll[gp][0]/nmol[gp]);
      printf ("Energie elect tot  : %f | %f par mol\n",tmpAll[gp][2],tmpAll[gp][2]/nmol[gp]);
      printf ("Energie cp pair ox : %f | %f par mol\n",tmpAll[gp][3],tmpAll[gp][3]/nmol[gp]);
      printf ("Energie cp pair sm : %f | %f par mol\n",tmpAll[gp][4],tmpAll[gp][4]/nmol[gp]);
      printf ("Energie cov sm     : %f | %f par mol\n",tmpAll[gp][5],tmpAll[gp][5]/nmol[gp]);

      Ecoh = tmpAll[gp][2] + tmpAll[gp][3] + tmpAll[gp][4] + tmpAll[gp][5];
      printf ("Energie totale     : %f | %f par mol\n",Ecoh,Ecoh/nmol[gp]);
      printf ("================================================= \n");
      printf ("    \n");

    }  // ----------------------------------------------------- Boucle sur les gp




  } // ------------ Call me == 0

  memory->destroy(nmol);
  memory->destroy(tmp);
  memory->destroy(tmpAll);
}

/* ----------------------------------------------------------------------
   Partie Electrostatique
   ----------------------------------------------------------------------*/

double PairSMTBQ::self(Param *param, double qi)
{
  double self_tmp;
  double s1=param->chi, s2=param->dj;

  self_tmp = qi*(s1+0.5*qi*s2);

  return self_tmp;
}

/* ---------------------------------------------------------------------- */

double PairSMTBQ::qfo_self(Param *param, double qi)
{
  double self_d;
  double s1 = param->chi;
  double s2 = param->dj;

  self_d = 0.0 ;
  self_d = s1+qi*s2;

  return self_d;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

void PairSMTBQ::tabqeq()
{
  int i,j,k,m,verbose;
  int nntype;
  double rc,s,r;
  double alf;

  int ii;
  double za,zb,ra,rb,gam,dgam,dza,dzb,
    d2zaa,d2zab,d2zbb,d2zra,d2zrb,d2gamr2,na,nb;
  double aCoeff,bCoeff,rcoupe,nang;

  int n = atom->ntypes;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  nmax = atom->nmax;

  verbose = 1;
  verbose = 0;

  nntype = int((n+1)*n/2);

  rc = cutmax ;
  alf = 0.3 ;
  //  alf = 0.2 ;


  if (verbose) printf ("kmax %d, ds %f, nmax %d\n",kmax,ds,nmax);
  if (verbose) printf ("nlocal = %d, nghost = %d\n",nlocal,nghost);
  if (verbose) printf ("nntypes %d, kmax %d, rc %f, n %d\n",nntype,kmax,rc,n);

  // allocate arrays

  memory->create(coultype,n,n,"pair:intype");
  memory->create(potqn,kmax+5,"pair:potqn");
  memory->create(dpotqn,kmax+5,"pair:dpotqn");
  memory->create(fafb,kmax+5,nntype,"pair:fafb");
  memory->create(dfafb,kmax+5,nntype,"pair:dfafb");
  memory->create(fafbOxOxSurf,kmax+5,"pair:fafbOxOxSurf");
  memory->create(dfafbOxOxSurf,kmax+5,"pair:dfafbOxOxSurf");
  memory->create(fafbTiOxSurf,kmax+5,"pair:fafbTiOxSurf");
  memory->create(dfafbTiOxSurf,kmax+5,"pair:dfafbTiOxSurf");

  memory->create(fafbOxOxBB,kmax+5,"pair:fafbOxOxBB");
  memory->create(dfafbOxOxBB,kmax+5,"pair:dfafbOxOxBB");
  memory->create(fafbTiOxBB,kmax+5,"pair:fafbTiOxB");
  memory->create(dfafbTiOxBB,kmax+5,"pair:dfafbTiOxBB");


  memory->create(ecov,nmax,"pair:ecov");
  memory->create(potmad,nmax,"pair:potmad");
  memory->create(potself,nmax,"pair:potself");
  memory->create(potcov,nmax,"pair:potcov");
  memory->create(sbcov,nmax,"pair:sbcov");
  memory->create(coord,nmax,"pair:coord");
  memory->create(sbmet,nmax,"pair:sbmet");
  memory->create(chimet,nmax,"pair:chimet");

  //  memory->create(nvsm,nmax,"pair:nvsm");
  //  memory->create(vsm,nmax,nmax,"pair:vsm");
  memory->create(flag_QEq,nmax,"pair:flag_QEq");

  memory->create(qf,nmax,"pair:qf");
  memory->create(q1,nmax,"pair:q1");
  memory->create(q2,nmax,"pair:q2");
  memory->create(tab_comm,nmax,"pair:tab_comm");

  memory->create(fct,31,"pair:fct");

  // set interaction number: 0-0=0, 1-1=1, 0-1=1-0=2

  m = 0; k = n;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j == i) {
        coultype[i][j] = m;
        m += 1;
      } else if (j != i && j > i) {
        coultype[i][j] = k;
        k += 1;
      } else if (j != i && j < i) {
        coultype[i][j] = coultype[j][i];
      }
      if (verbose) printf ("i %d, j %d, coultype %d\n",i,j,coultype[i][j]);
    }
  }

  // -------- Tabqn --------

  // -------------------
  //   Ouverture du fichier
  //   ofstream fichier("tabqeq.txt", ios::out | ios::trunc) ;
  // -------------------

  double mu;

  mu = erfc(alf*rc)/rc ;

  //if (fichier) fichier <<" r  -  potqn " <<endl ;

  //-------------------------
  for (k=0; k < kmax+5; k++)
    //-------------------------
    {
      s = static_cast<double>(k)*ds ; r = sqrt(s);
      if (k==0) r=10e-30;
      potqn[k] = 14.4*(erfc(alf*r)/r - mu) ;

      // $$$ Here is (1/r)*dE/dr
      dpotqn[k] = -14.4*( (erfc(alf*r)/(r*r) + 2.0*alf/MY_PIS/r*exp(-alf*alf*r*r))/r  ) ;
    }

  Vself = -14.4*(alf/MY_PIS + mu*0.5) ;

  // --------------------
  // default arrays to zero

  for (i = 0; i < kmax+5; i ++) {
    for (j = 0; j < nntype; j ++) {
      fafb[i][j] = 0.0;
      dfafb[i][j] = 0.0;
    }
    fafbOxOxSurf[i] = 0.0;
    fafbTiOxSurf[i] = 0.0;
    dfafbOxOxSurf[i] = 0.0;
    dfafbTiOxSurf[i] = 0.0;

    fafbOxOxBB[i] = 0.0;
    fafbTiOxBB[i] = 0.0;
    dfafbOxOxBB[i] = 0.0;
    dfafbTiOxBB[i] = 0.0;
  }


  // Set Tabqeq
  double dij,ddij;

  // -------------------
  //   Ouverture du fichier
  //ofstream fichier("dtabqeq.txt", ios::out | ios::trunc) ;
  // -------------------

  //if (fichier) fichier <<" k , r , fafb , dfafb , dfafb2 , dgam , d(1/r) , dpotqn" <<endl ;


  rcoupe = cutmax ;
  double cang ;

  for (i = 0; i < n ; i++){
    for (j = i; j < n ; j++){

      rc = cutmax; if (verbose) printf ("cutmax %f\n",cutmax);
      m = coultype[i][j] ;
      na = params[i].ne ;
      nb = params[j].ne ;
      za = params[i].dzeta ;
      zb = params[j].dzeta ;
      ra = params[i].R;
      rb = params[j].R;

      ii = 0 ; nang =cang= 5.0 ;
      // --------------------------
      for (k = 0; k < kmax+5; k++)
        // --------------------------
        {
          gam = dgam = dza = dzb = d2zaa = d2zab =
            d2zbb = d2zra = d2zrb = d2gamr2 = 0.0 ;
          dij = 0.0 ;

          s = static_cast<double>(k)*ds ; r = sqrt(s) ;
          if (k==0) r=10e-30;

          gammas(na,nb,za,zb,r,gam,dgam,dza,dzb,
                 d2zaa,d2zab,d2zbb,d2zra,d2zrb,d2gamr2) ;

          // --- Jij

          dij = 14.4 * (1.0/r - static_cast<double>(gam));
          ddij = 14.4 * (-1.0/(r*r) - static_cast<double>(dgam)) ;

          // Cutting Fonction

          if (dij < 0.01 && ii==0)
            {
              ii=2;
              if (ii==2) if (verbose) printf ("rc : %f\n",r);
              rc = r ; ii=1 ;
              if ((rc+nang)>rcoupe) nang = rcoupe - rc ;
              bCoeff =  (2*dij+ddij*nang)/(dij*nang);
              aCoeff = dij*exp(-bCoeff*rc) /square(nang);
            }
          if (r > rc) {dij = aCoeff *square(r- rc-nang) *exp(bCoeff*r);
            ddij = aCoeff*(r- rc-nang) *(2+bCoeff*(r-rc-nang))*exp(bCoeff*r);
          }

          if (r > (rc+nang)) {dij = 0.0 ; ddij = 0.0;}

          fafb[k][m] = potqn[k] - dij ;
          if (k == 1) fafb[0][m] = fafb[k][m] ;

          dfafb[k][m] = dpotqn[k] - ddij/r ;
        }

      // Make the table fafbOxOxSurf
      rc = cutmax;
      if(strcmp(params[i].nom,"O")==0 || strcmp(params[j].nom,"O")==0){
        if(strcmp(params[i].nom,"O")==0) {
          ra = ROxSurf;
          za = (2.0*params[i].ne + 1.0)/(4.0*ra);}
        if(strcmp(params[j].nom,"O")==0) {
          rb = ROxSurf;
          zb = (2.0*params[j].ne + 1.0)/(4.0*rb); }

        ii = 0 ; nang =cang= 5.0 ;
        // --------------------------
        for (k = 0; k < kmax+5; k++)
          // --------------------------
          {
            gam = dgam = dza = dzb = d2zaa = d2zab =
              d2zbb = d2zra = d2zrb = d2gamr2 = 0.0 ;
            dij = 0.0 ;

            s = static_cast<double>(k)*ds ; r = sqrt(s) ;
            if (k==0) r=10e-30;

            gammas(na,nb,za,zb,r,gam,dgam,dza,dzb,
                   d2zaa,d2zab,d2zbb,d2zra,d2zrb,d2gamr2) ;

            // --- Jij

            dij = 14.4 * (1.0/r - static_cast<double>(gam));
            ddij = 14.4 * (-1.0/(r*r) - static_cast<double>(dgam)) ;

            if (dij < 0.01 && ii==0)
              {
                ii=2;
                if (ii==2) if (verbose) printf ("rc : %f\n",r);
                rc = r ; ii=1 ;
                if ((rc+nang)>rcoupe) nang = rcoupe - rc ;
                bCoeff =  (2*dij+ddij*nang)/(dij*nang);
                aCoeff = dij*exp(-bCoeff*rc) /square(nang);
              }
            if (r > rc) {dij = aCoeff *square(r- rc-nang) *exp(bCoeff*r);
              ddij = aCoeff*(r- rc-nang) *(2+bCoeff*(r-rc-nang))*exp(bCoeff*r);
            }


            if (r > (rc+nang)) {dij = 0.0 ; ddij = 0.0;}

            if(strcmp(params[i].nom,"O")==0 && strcmp(params[j].nom,"O")==0){
              fafbOxOxSurf[k] = potqn[k] - dij ;
              if (k == 1) fafbOxOxSurf[0] = fafbOxOxSurf[k] ;

              dfafbOxOxSurf[k] = dpotqn[k] - ddij/r ;
            }
            else { fafbTiOxSurf[k] = potqn[k] - dij ;
              if (k == 1) fafbTiOxSurf[0] = fafbTiOxSurf[k] ;

              dfafbTiOxSurf[k] = dpotqn[k] - ddij/r ;}

          }


      } //for k
      //end of make the table fafbOxOxSurf

      // Makes the table fafbOxOxBB
      rc = cutmax;
      if(strcmp(params[i].nom,"O")==0 || strcmp(params[j].nom,"O")==0){
        if(strcmp(params[i].nom,"O")==0) {
          ra = ROxBB;
          za = (2.0*params[i].ne + 1.0)/(4.0*ra);}
        if(strcmp(params[j].nom,"O")==0) {
          rb = ROxBB;
          zb = (2.0*params[j].ne + 1.0)/(4.0*rb); }


        ii = 0 ; nang =cang= 5.0 ;
        // --------------------------
        for (k = 0; k < kmax+5; k++)
          // --------------------------
          {
            gam = dgam = dza = dzb = d2zaa = d2zab =
              d2zbb = d2zra = d2zrb = d2gamr2 = 0.0 ;
            dij = 0.0 ;

            s = static_cast<double>(k)*ds ; r = sqrt(s) ;
            if (k==0) r=10e-30;

            gammas(na,nb,za,zb,r,gam,dgam,dza,dzb,
                   d2zaa,d2zab,d2zbb,d2zra,d2zrb,d2gamr2) ;

            // --- Jij

            dij = 14.4 * (1.0/r - static_cast<double>(gam));
            ddij = 14.4 * (-1.0/(r*r) - static_cast<double>(dgam)) ;

            if (dij < 0.01 && ii==0)  {
              ii=2;
              if (ii==2) if (verbose) printf ("rc : %f\n",r);
              rc = r ; ii=1 ;
              if ((rc+nang)>rcoupe) nang = rcoupe - rc ;
              bCoeff =  (2*dij+ddij*nang)/(dij*nang);
              aCoeff = dij*exp(-bCoeff*rc) /square(nang);
            }
            if (r > rc) {
              dij = aCoeff *square(r- rc-nang) *exp(bCoeff*r);
              ddij = aCoeff*(r- rc-nang) *(2+bCoeff*(r-rc-nang))*exp(bCoeff*r);
            }


            if (r > (rc+nang)) {dij = 0.0 ; ddij = 0.0;}

            if(strcmp(params[i].nom,"O")==0 && strcmp(params[j].nom,"O")==0){
              fafbOxOxBB[k] = potqn[k] - dij ;
              if (k == 1) fafbOxOxBB[0] = fafbOxOxBB[k] ;
              dfafbOxOxBB[k] = dpotqn[k] - ddij/r ; }
            else { fafbTiOxBB[k] = potqn[k] - dij ;
              if (k == 1) fafbTiOxBB[0] = fafbTiOxBB[k] ;
              dfafbTiOxBB[k] = dpotqn[k] - ddij/r ;
            }
          }

      } //for k
      //end of make the table fafbOxOxBB

    }
  } //for i,j

  //if (fichier) fichier.close() ;

}

/* ---------------------------------------------------------------------*/

void PairSMTBQ::potqeq(int i, int j, double qi, double qj, double rsq,
                       double &fforce, int /*eflag*/, double &eng)
{

  /* ===================================================================
     Coulombian energy calcul between i and j atoms
     with fafb table make in sm_table().
     fafb[i][j] : i is the table's step (r)
     j is the interaction's # (in intype[itype][jtype])
     dfafb is derivate energy (force)
     ==================================================================== */

  int itype,jtype,l,m;
  double r,t1,t2,sds,xi,engBulk,engSurf,fforceBulk,fforceSurf,dcoordloc,dcoupureloc;
  double engBB,fforceBB, dIntfcoup2loc,iCoord,jCoord,iIntfCoup2,jIntfCoup2;

  int *type = atom->type;
  //  int n = atom->ntypes;

  itype = map[type[i]];
  jtype = map[type[j]];
  m = coultype[itype][jtype];

  r = rsq;
  sds = r/ds ;  l = int(sds) ;
  xi = sds - static_cast<double>(l) ;


  iCoord=coord[i];
  jCoord=coord[j];
  iIntfCoup2= Intfcoup2(iCoord,coordOxBulk,0.15);
  jIntfCoup2= Intfcoup2(jCoord,coordOxBulk,0.15);

  // ---- Energies Interpolation

  t1 = fafb[l][m] + (fafb[l+1][m] - fafb[l][m])*xi;
  t2 = fafb[l+1][m] + (fafb[l+2][m] - fafb[l+1][m])*(xi-1.0);

  engBulk = qi*qj*(t1 + (t2 - t1)*xi/2.0);
  eng=engBulk;


  // ---- Forces Interpolation

  t1 = dfafb[l][m] + (dfafb[l+1][m] - dfafb[l][m])*xi;
  t2 = dfafb[l+1][m] + (dfafb[l+2][m] - dfafb[l+1][m])*(xi-1);


  fforce = - qi*qj*(t1 + (t2 - t1)*xi/2.0) ;


  if(strcmp(params[itype].nom,"O")==0 || strcmp(params[jtype].nom,"O")==0){

    if(strcmp(params[itype].nom,"O")==0 && strcmp(params[jtype].nom,"O")==0){
      // between two oxygens

      t1 = fafbOxOxSurf[l] + (fafbOxOxSurf[l+1] - fafbOxOxSurf[l])*xi;
      t2 = fafbOxOxSurf[l+1] + (fafbOxOxSurf[l+2] - fafbOxOxSurf[l+1])*(xi-1.0);
      engSurf = qi*qj*(t1 + (t2 - t1)*xi/2.0);

      t1 = fafbOxOxBB[l] + (fafbOxOxBB[l+1] - fafbOxOxBB[l])*xi;
      t2 = fafbOxOxBB[l+1] + (fafbOxOxBB[l+2] - fafbOxOxBB[l+1])*(xi-1.0);
      engBB = qi*qj*(t1 + (t2 - t1)*xi/2.0);

      eng= engBulk + (iCoord+jCoord-2*coordOxBulk)/(2*(coordOxBB-coordOxBulk)) *(engBB-engBulk)
        + (iIntfCoup2+jIntfCoup2)*((engBulk-engSurf)/(2*(coordOxBulk-coordOxSurf))
                                   - (engBB-engBulk)/(2*(coordOxBB-coordOxBulk))) ;


      // ---- Interpolation des forces

      fforceBulk=fforce;
      t1 = dfafbOxOxSurf[l] + (dfafbOxOxSurf[l+1] - dfafbOxOxSurf[l])*xi;
      t2 = dfafbOxOxSurf[l+1] + (dfafbOxOxSurf[l+2] - dfafbOxOxSurf[l+1])*(xi-1);
      fforceSurf = - qi*qj*(t1 + (t2 - t1)*xi/2.0) ;

      t1 = dfafbOxOxBB[l] + (dfafbOxOxBB[l+1] - dfafbOxOxBB[l])*xi;
      t2 = dfafbOxOxBB[l+1] + (dfafbOxOxBB[l+2] - dfafbOxOxBB[l+1])*(xi-1);
      fforceBB = - qi*qj*(t1 + (t2 - t1)*xi/2.0) ;

      fforce= fforceBulk + (iCoord+jCoord-2*coordOxBulk)/(2*(coordOxBB-coordOxBulk))*(fforceBB-fforceBulk)
        + (iIntfCoup2+jIntfCoup2)*((fforceBulk-fforceSurf)/(2*(coordOxBulk-coordOxSurf))
                                   - (fforceBB-fforceBulk)/(2*(coordOxBB-coordOxBulk))) ;

    } else {    // between metal and oxygen

      t1 = fafbTiOxSurf[l] + (fafbTiOxSurf[l+1] - fafbTiOxSurf[l])*xi;
      t2 = fafbTiOxSurf[l+1] + (fafbTiOxSurf[l+2] - fafbTiOxSurf[l+1])*(xi-1.0);
      engSurf = qi*qj*(t1 + (t2 - t1)*xi/2.0);
      t1 = fafbTiOxBB[l] + (fafbTiOxBB[l+1] - fafbTiOxBB[l])*xi;
      t2 = fafbTiOxBB[l+1] + (fafbTiOxBB[l+2] - fafbTiOxBB[l+1])*(xi-1.0);
      engBB = qi*qj*(t1 + (t2 - t1)*xi/2.0);

      if(strcmp(params[jtype].nom,"O")==0) //the atom j is an oxygen
        {         iIntfCoup2=jIntfCoup2;
          iCoord=jCoord;        }

      eng = engBulk + (engBulk-engSurf)/(coordOxBulk-coordOxSurf) * iIntfCoup2
        + (engBB-engBulk)/(coordOxBB-coordOxBulk) * (iCoord-coordOxBulk-iIntfCoup2);


      // ---- Forces Interpolation

      fforceBulk=fforce;
      t1 = dfafbTiOxSurf[l] + (dfafbTiOxSurf[l+1] - dfafbTiOxSurf[l])*xi;
      t2 = dfafbTiOxSurf[l+1] + (dfafbTiOxSurf[l+2] - dfafbTiOxSurf[l+1])*(xi-1);
      fforceSurf = - qi*qj*(t1 + (t2 - t1)*xi/2.0) ;

      t1 = dfafbTiOxBB[l] + (dfafbTiOxBB[l+1] - dfafbTiOxBB[l])*xi;
      t2 = dfafbTiOxBB[l+1] + (dfafbTiOxBB[l+2] - dfafbTiOxBB[l+1])*(xi-1);
      fforceBB = - qi*qj*(t1 + (t2 - t1)*xi/2.0) ;

      dcoordloc =  fcoupured(sqrt(r),r1Coord,r2Coord) ;


      dcoupureloc =  fcoupured(iCoord,coordOxSurf,coordOxBulk) ;
      dIntfcoup2loc=   fcoup2(iCoord,coordOxBulk,0.15)*dcoupureloc ;
      fforce = fforceBulk + 1/(coordOxBulk-coordOxSurf) * ((fforceBulk-fforceSurf)* iIntfCoup2
                                                           - (engBulk-engSurf) *dIntfcoup2loc)
        + 1/(coordOxBB-coordOxBulk) * ((fforceBB-fforceBulk)*(iCoord-coordOxBulk- iIntfCoup2)
                                       - (engBB-engBulk) *(dcoordloc-dIntfcoup2loc));


    }
  }
}

/* -------------------------------------------------------------------- */

void PairSMTBQ::pot_ES (int i, int j, double rsq, double &eng)
{

  /* ===================================================================
     Coulombian potential energy calcul between i and j atoms
     with fafb table make in sm_table().
     fafb[i][j] : i is the table's step (r)
     j is the interaction's # (in intype[itype][jtype])
     dfafb is derivate energy (force)
     ==================================================================== */

  int itype,jtype,l,m;
  double r,t1,t2,sds,xi,engBulk,engSurf;
  double engBB,iCoord,jCoord,iIntfCoup2,jIntfCoup2;

  int *type = atom->type;
  //  int n = atom->ntypes;

  itype = map[type[i]];
  jtype = map[type[j]];
  m = coultype[itype][jtype];

  r = rsq;
  sds = r/ds ;  l = int(sds) ;
  xi = sds - static_cast<double>(l) ;


  iCoord=coord[i];
  jCoord=coord[j];
  iIntfCoup2= Intfcoup2(iCoord,coordOxBulk,0.15);
  jIntfCoup2= Intfcoup2(jCoord,coordOxBulk,0.15);

  // ---- Energies Interpolation

  t1 = fafb[l][m] + (fafb[l+1][m] - fafb[l][m])*xi;
  t2 = fafb[l+1][m] + (fafb[l+2][m] - fafb[l+1][m])*(xi-1.0);


  eng = (t1 + (t2 - t1)*xi/2.0);
  engBulk=eng;


  if(itype==0 || jtype==0){

    if(itype==0 && jtype==0){   // between two oxygens

      t1 = fafbOxOxSurf[l] + (fafbOxOxSurf[l+1] - fafbOxOxSurf[l])*xi;
      t2 = fafbOxOxSurf[l+1] + (fafbOxOxSurf[l+2] - fafbOxOxSurf[l+1])*(xi-1.0);
      engSurf = (t1 + (t2 - t1)*xi/2.0);

      t1 = fafbOxOxBB[l] + (fafbOxOxBB[l+1] - fafbOxOxBB[l])*xi;
      t2 = fafbOxOxBB[l+1] + (fafbOxOxBB[l+2] - fafbOxOxBB[l+1])*(xi-1.0);
      engBB = (t1 + (t2 - t1)*xi/2.0);

      eng= engBulk + (iCoord+jCoord-2*coordOxBulk)/(2*(coordOxBB-coordOxBulk))*(engBB-engBulk)
        + (iIntfCoup2+jIntfCoup2)*((engBulk-engSurf)/(2*(coordOxBulk-coordOxSurf))
                                   - (engBB-engBulk)/(2*(coordOxBB-coordOxBulk))) ;

    } else {    // between metal and oxygen

      t1 = fafbTiOxSurf[l] + (fafbTiOxSurf[l+1] - fafbTiOxSurf[l])*xi;
      t2 = fafbTiOxSurf[l+1] + (fafbTiOxSurf[l+2] - fafbTiOxSurf[l+1])*(xi-1.0);
      engSurf = (t1 + (t2 - t1)*xi/2.0);

      t1 = fafbTiOxBB[l] + (fafbTiOxBB[l+1] - fafbTiOxBB[l])*xi;
      t2 = fafbTiOxBB[l+1] + (fafbTiOxBB[l+2] - fafbTiOxBB[l+1])*(xi-1.0);
      engBB = (t1 + (t2 - t1)*xi/2.0);

      if (jtype==0) {   //the atom j is an oxygen
        iIntfCoup2=jIntfCoup2;
        iCoord=jCoord;
      }

      eng = engBulk + (engBulk-engSurf)/(coordOxBulk-coordOxSurf)*iIntfCoup2
        + (engBB-engBulk)/(coordOxBB-coordOxBulk) * (iCoord-coordOxBulk-iIntfCoup2);


    }



  }


}

/* -------------------------------------------------------------------- */

void PairSMTBQ::pot_ES2 (int i, int j, double rsq, double &pot)
{
  int l,m,itype,jtype;
  double sds,xi,t1,t2,r;

  int *type = atom->type;


  if (sqrt(rsq) > cutmax) return ;

  itype = map[type[i]];
  jtype = map[type[j]];
  m = coultype[itype][jtype];

  r = rsq ;
  sds = r/ds ;  l = int(sds) ;
  xi = sds - static_cast<double>(l) ;

  // ---- Energies Interpolation

  t1 = fafb[l][m] + (fafb[l+1][m] - fafb[l][m])*xi;
  t2 = fafb[l+1][m] + (fafb[l+2][m] - fafb[l+1][m])*(xi-1.0);

  pot = (t1 + (t2 - t1)*xi/2.0) ;

}

/* --------------------------------------------------------------------
   Oxygen-Oxygen Interaction
   -------------------------------------------------------------------- */

void PairSMTBQ::rep_OO(Intparam *intparam, double rsq, double &fforce,
                       int /*eflag*/, double &eng)
{
  double r,tmp_exp,tmp;
  double A = intparam->abuck ;
  double rho = intparam->rhobuck ;

  r = sqrt(rsq);
  tmp = - r/rho ;
  tmp_exp = exp( tmp );

  eng = A * tmp_exp ;

  fforce = A/rho * tmp_exp / r ; //( - )

}


void PairSMTBQ::Attr_OO(Intparam *intparam, double rsq, double &fforce,
                        int /*eflag*/, double &eng)
{
  double r,tmp_exp;
  double aOO = intparam->aOO ;
  double bOO = intparam->bOO ;
  double r1OO = intparam->r1OO ;
  double r2OO = intparam->r2OO ;

  r = sqrt(rsq);
  tmp_exp= exp( bOO* r);
  eng = aOO * tmp_exp* fcoupure(r,r1OO,r2OO);

  fforce = - (aOO*bOO * tmp_exp * fcoupure(r,r1OO,r2OO)+ aOO*tmp_exp *fcoupured(r,r1OO,r2OO))/ r ; //( - )

}


/* ----------------------------------------------------------------------
   covalente Interaction
   ----------------------------------------------------------------------*/


void PairSMTBQ::tabsm()
{
  int k,m;
  double s,r,tmpb,tmpr,fcv,fcdv;

  memory->create(tabsmb,kmax,maxintsm+1,"pair:tabsmb");
  memory->create(tabsmr,kmax,maxintsm+1,"pair:tabsmr");
  memory->create(dtabsmb,kmax,maxintsm+1,"pair:dtabsmb");
  memory->create(dtabsmr,kmax,maxintsm+1,"pair:dtabsmr");


  for (m = 0; m <= maxintparam; m++) {

    if (intparams[m].intsm == 0) continue;

    double rc1 = intparams[m].dc1;
    double rc2 = intparams[m].dc2;
    double A = intparams[m].a;
    double p = intparams[m].p;
    double Ksi = intparams[m].ksi;
    double q = intparams[m].q;
    double rzero = intparams[m].r0;
    int sm = intparams[m].intsm;


    for (k=0; k < kmax; k++)
      {
        s = static_cast<double>(k)*ds ; r = sqrt(s);
        if (k==0) r=10e-30;
        tmpb = exp( -2.0*q*(r/rzero - 1.0));
        tmpr = exp( -p*(r/rzero - 1.0));

        if (r <= rc1)
          {

            // -- Energy
            tabsmb[k][sm] = Ksi*Ksi * tmpb ;
            tabsmr[k][sm] = A * tmpr ;

            // -- Force
            /*  dtabsmb ne correspond pas vraiment a une force puisqu'il y a le /r
                (on a donc une unite force/distance). Le programme multiplie ensuite
                (dans le  PairSMTBQ::compute ) dtabsmb par la projection du vecteur r
                sur un axe x (ou y ou z) pour determiner la composante de la force selon
                cette direction. Donc tout est ok au final.  */

            dtabsmb[k][sm] = - 2.0 *Ksi*Ksi* q/rzero * tmpb /r;
            dtabsmr[k][sm] = - A * p/rzero * tmpr/r ;

          } // if

        else if (r > rc1 && r <= rc2)
          {

            // -- Energie
            fcv = fcoupure(r,intparams[sm].dc1,intparams[sm].dc2);
            tabsmb[k][sm] = fcv* Ksi*Ksi * tmpb ;
            tabsmr[k][sm] = fcv* A * tmpr ;

            // -- Force
            /*   dtabsmb ne correspond pas vraiment a une force puisqu'il y a le /r
                 (on a donc une unite force/distance). Le programme multiplie ensuite
                 (dans le  PairSMTBQ::compute ) d tabsmb par la projection du vecteur
                 r sur un axe x (ou y ou z) pour determiner la composante de la force
                 selon cette direction. Donc tout est ok au final. */

            fcdv = fcoupured(r,intparams[sm].dc1,intparams[sm].dc2);
            dtabsmb[k][sm] = (fcv*( - 2.0 *Ksi*Ksi* q/rzero * tmpb )+fcdv* Ksi*Ksi * tmpb )/r ;
            dtabsmr[k][sm] = (fcv*( - A * p/rzero * tmpr )+fcdv*A * tmpr  )/r ;

          }

        else
          {

            // -- Energie
            tabsmb[k][sm] = 0.0;
            tabsmr[k][sm] = 0.0;

            // -- Force
            dtabsmb[k][sm] = 0.0;
            dtabsmr[k][sm] = 0.0;

          }



      }   // for kmax


  } // for maxintparam

}





/* -------------------------------------------------------------- */

void PairSMTBQ::repulsive(Intparam *intparam, double rsq, int /*i*/, int /*j*/,
                          double &fforce, int /*eflag*/, double &eng)
{

  /* ================================================
     rsq    : square of ij distance
     fforce : repulsion force
     eng    : repulsion energy
     eflag  : Si oui ou non on calcule l'energie
     =================================================*/

  int l;
  double r,sds,xi,t1,t2,dt1,dt2,sweet;

  double rrcs = intparam->dc2;
  int sm = intparam->intsm;

  //  printf ("On rentre dans repulsive\n");


  r = rsq;
  if (sqrt(r) > rrcs) return ;

  sds = r/ds ;  l = int(sds) ;
  xi = sds - static_cast<double>(l) ;

  t1 = tabsmr[l][sm] + (tabsmr[l+1][sm] - tabsmr[l][sm])*xi ;
  t2 = tabsmr[l+1][sm] + (tabsmr[l+2][sm] - tabsmr[l+1][sm])*(xi-1.0) ;

  dt1 = dtabsmr[l][sm] + (dtabsmr[l+1][sm] - dtabsmr[l][sm])*xi ;
  dt2 = dtabsmr[l+1][sm] + (dtabsmr[l+2][sm] - dtabsmr[l+1][sm])*(xi-1.0) ;

  if (strcmp(intparam->mode,"oxide") == 0)
    {
      fforce = - 2.0*(dt1 + (dt2 - dt1)*xi/2.0);
      eng = (t1 + (t2 - t1)*xi/2.0) ;
    }
  else if (strcmp(intparam->mode,"metal") == 0)
    {
      sweet = 1.0;
      fforce = - 2.0*(dt1 + (dt2 - dt1)*xi/2.0) * sweet ;
      eng = (t1 + (t2 - t1)*xi/2.0) * sweet ;
    }

}


/* --------------------------------------------------------------------------------- */


void PairSMTBQ::attractive(Intparam *intparam, double rsq,
                           int /*eflag*/, int i, double /*iq*/, int /*j*/, double /*jq*/)
{
  int itype,l;
  double r,t1,t2,xi,sds;
  double sweet,mu;

  double rrcs = intparam->dc2;
  int *type = atom->type;
  int sm = intparam->intsm;

  itype = map[type[i]];

  r = rsq;
  if (sqrt(r) > rrcs) return ;


  sds = r/ds ;  l = int(sds) ;
  xi = sds - static_cast<double>(l) ;

  t1 = tabsmb[l][sm] + (tabsmb[l+1][sm] - tabsmb[l][sm])*xi ;
  t2 = tabsmb[l+1][sm] + (tabsmb[l+2][sm] - tabsmb[l+1][sm])*(xi-1.0) ;



  if (strcmp(intparam->mode,"oxide") == 0) {
    mu = 0.5*(sqrt(params[1].sto) + sqrt(params[0].sto));

    //      dq = fabs(params[itype].qform) - fabs(iq);
    //      dqcov = dq*(2.0*ncov/params[itype].sto - dq);

    sbcov[i] += (t1 + (t2 - t1)*xi/2.0) *params[itype].sto*mu*mu;

    //      if (i < 10) printf ("i %d, iq %f sbcov %f \n",i,iq,sbcov[i]);

    if (sqrt(r)<r1Coord) { coord[i] +=  1 ; }
    else if (sqrt(r)<r2Coord){ coord[i] +=  fcoupure(sqrt(r),r1Coord,r2Coord) ;}


  }
  else if (strcmp(intparam->mode,"metal") == 0) {
    sweet = 1.0;
    sbmet[i] += (t1 + (t2 - t1)*xi/2.0) * sweet ;
  }

}

/* ---------------------------------------------------------------------- */

void PairSMTBQ::f_att(Intparam *intparam, int i, int j,double rsq, double &fforce)
{
  int itype,jtype,l;
  int *type = atom->type;

  double r,sds,xi,dt1,dt2,dq,dqcovi,dqcovj;
  double fcov_ij,fcov_ji,sweet,iq,jq,mu;

  int sm = intparam->intsm;
  double *q = atom->q;

  itype = map[type[i]];
  jtype = map[type[j]];
  iq = q[i] ; jq = q[j];

  r = rsq;

  sds = r/ds ;  l = int(sds) ;
  xi = sds - static_cast<double>(l) ;

  dt1 = dtabsmb[l][sm] + (dtabsmb[l+1][sm] - dtabsmb[l][sm])*xi ;
  dt2 = dtabsmb[l+1][sm] + (dtabsmb[l+2][sm] - dtabsmb[l+1][sm])*(xi-1.0) ;

  dq = fabs(params[itype].qform) - fabs(iq);
  dqcovi = dq*(2.0*ncov/params[itype].sto - dq);

  dq = fabs(params[jtype].qform) - fabs(jq);
  dqcovj = dq*(2.0*ncov/params[jtype].sto - dq);

  if (strcmp(intparam->mode,"oxide") == 0) {
    //------------------------------------------
    mu = 0.5*(sqrt(params[1].sto) + sqrt(params[0].sto));
    fcov_ij = (dt1 + (dt2 - dt1)*xi/2.0) * dqcovi *params[itype].sto*mu*mu;
    fcov_ji = (dt1 + (dt2 - dt1)*xi/2.0) * dqcovj *params[jtype].sto*mu*mu;

    fforce = 0.5 * ( fcov_ij/sqrt(sbcov[i]*dqcovi + sbmet[i])
                     + fcov_ji/sqrt(sbcov[j]*dqcovj + sbmet[j]) ) ;
  }

  else if (strcmp(intparam->mode,"metal") == 0) {
    //-----------------------------------------------
    sweet = 1.0;
    fcov_ij = (dt1 + (dt2 - dt1)*xi/2.0) * sweet ;

    fforce = 0.5 * fcov_ij*( 1.0/sqrt(sbcov[i]*dqcovi + sbmet[i])
                             + 1.0/sqrt(sbcov[j]*dqcovj + sbmet[j]) ) ;
  }

}

/* ---------------------------------------------------------------------- */

void PairSMTBQ::pot_COV(Param *param, int i, double &qforce)
{
  double iq,dq,DQ,sign;

  double *q = atom->q;
  double qform = param->qform;
  double sto = param->sto;

  sign = qform / fabs(qform);
  iq = q[i];

  dq = fabs(qform) - fabs(iq);
  DQ = dq*(2.0*ncov/sto - dq);

  if (fabs(iq) < 1.0e-7 || fabs(sbcov[i]) < 1.0e-7) {
    qforce = 0.0; }
  else {
    qforce = sign*sbcov[i]/sqrt(sbcov[i]*DQ + sbmet[i])*(ncov/sto - dq) ;
  }

}

/* ---------------------------------------------------------------------- */

double PairSMTBQ::potmet(Intparam *intparam, double rsq,
                         int i, double iq, int j, double jq)
{
  int l,itype,jtype;
  int *type = atom->type;
  double chi,sds,xi,t1,t2,r,dsweet,dq,dqcovi,dqcovj;

  int sm = intparam->intsm;
  itype = map[type[i]];
  jtype = map[type[j]];

  r = rsq;
  sds = r/ds ;  l = int(sds) ;
  xi = sds - static_cast<double>(l) ;

  t1 = tabsmb[l][sm] + (tabsmb[l+1][sm] - tabsmb[l][sm])*xi ;
  t2 = tabsmb[l+1][sm] + (tabsmb[l+2][sm] - tabsmb[l+1][sm])*(xi-1.0) ;

  dq = fabs(params[itype].qform) - fabs(iq);
  dqcovi = dq*(2.0*ncov/params[itype].sto - dq);

  dq = fabs(params[jtype].qform) - fabs(jq);
  dqcovj = dq*(2.0*ncov/params[jtype].sto - dq);

  dsweet = 0.0;
  chi = (t1 + (t2 - t1)*xi/2.0) * dsweet *( 1.0/(2.0*sqrt(sbcov[i]*dqcovi+sbmet[i]))
                                            + 1.0/(2.0*sqrt(sbcov[j]*dqcovj+sbmet[j])) );

  return chi;
}


/* ----------------------------------------------------------------------
   Cutting Function
   ----------------------------------------------------------------------*/


/* -------------------------------------------------------------------- */

double PairSMTBQ::fcoupure(double r, double rep_dc1, double rep_dc2)
{
  double ddc,a3,a4,a5,x;

  if (r<rep_dc1)
    {return 1;}
  else if (r> rep_dc2)
    {return 0;}
  else
    {

      ddc = rep_dc2 - rep_dc1;
      x = r - rep_dc2;

      a3 = -10/(ddc*ddc*ddc);
      a4 = -15/(ddc*ddc*ddc*ddc);
      a5 = -6/(ddc*ddc*ddc*ddc*ddc);

      return x*x*x*(a3 + x*(a4 + x*a5));}
}

/* ----------------------------------------------------------------------
   Derivate of cutting function
   ----------------------------------------------------------------------*/


/* ----------------------------------------------------------------------- */

double PairSMTBQ::fcoupured(double r,  double rep_dc1, double rep_dc2)
{

  double ddc,a3,a4,a5,x;

  if ( r>rep_dc1 && r<rep_dc2) {
    ddc = rep_dc2 - rep_dc1;
    x = r - rep_dc2;

    a3 = -10/(ddc*ddc*ddc);
    a4 = -15/(ddc*ddc*ddc*ddc);
    a5 = -6/(ddc*ddc*ddc*ddc*ddc);

    return  x*x*(3*a3 + x*(4*a4 + 5*x*a5));}
  else {return 0;}
}


/* ----------------------------------------------------------------------
   cutting function for derive (force)
   ----------------------------------------------------------------------*/


/* -------------------------------------------------------------------- */



double PairSMTBQ::fcoup2(double c, double x, double delta)
{
  double dc;

  if (c<x-delta)
    {return 1;}
  else if (c> x+delta)
    {return 0;}
  else
    {
      dc = c - x-delta;
      return dc*dc*(3*delta+dc)/(4*delta*delta*delta);
    }
}

/* ----------------------------------------------------------------------
   Primitive of cutting function for derive (force)
   ----------------------------------------------------------------------*/


/* -------------------------------------------------------------------- */

double PairSMTBQ::Primfcoup2(double c, double x, double delta)
{

  return (c*(c*c*c - 4* c*c* x - 4* (x - 2 *delta) * (x+delta)*(x+delta) +
             6* c *(x*x - delta*delta)))/(16* delta*delta*delta);

}


/* ----------------------------------------------------------------------
   Integral of cutting function for derive (force) between x and c
   ----------------------------------------------------------------------*/


/* -------------------------------------------------------------------- */

double PairSMTBQ::Intfcoup2(double c, double x, double delta)
{

  if (c<x-delta)
    {return c - x + delta + Primfcoup2(x-delta,x,delta) - Primfcoup2(x,x,delta) ;}
  else if (c> x+delta)
    {return  Primfcoup2(x+delta,x,delta) - Primfcoup2(x,x,delta) ;}
  else
    {
      return  Primfcoup2(c,x,delta) - Primfcoup2(x,x,delta) ;}
}


/* ---------------------------------------------------------------------
   Energy derivation respect charge Q
   --------------------------------------------------------------------- */

void PairSMTBQ::QForce_charge(int loop)
{
  int i,j,ii,jj,jnum;
  int itype,jtype,m,gp;
  double xtmp,ytmp,ztmp;
  double rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double iq,jq,fqi,fqj,fqij,fqij2,fqjj;
  int eflag = 0;


  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int step = update->ntimestep;

  int inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over full neighbor list of my atoms

  fqi = fqj = fqij = fqij2 = fqjj = 0.0;


  //   ==================
  if (loop == 0) {
    //   ==================

    memset(sbcov,0,sizeof(double)*atom->nmax);
    memset(coord,0,sizeof(double)*atom->nmax);
    memset(sbmet,0,sizeof(double)*atom->nmax);

    for (ii = 0; ii < inum; ii ++) {
      //--------------------------------
      i = ilist[ii];
      itype = map[type[i]];

      gp = flag_QEq[i];

      itype = map[type[i]];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      iq = q[i];




      // two-body interactions

      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        //  -------------------------------
        j = jlist[jj];
        j &= NEIGHMASK;

        jtype = map[type[j]];
        jq = q[j];

        m = intype[itype][jtype];
        if (intparams[m].intsm == 0) continue ;


        const double delx = x[j][0] - xtmp;
        const double dely = x[j][1] - ytmp;
        const double delz = x[j][2] - ztmp;
        rsq = delx*delx + dely*dely + delz*delz;

        //    Covalente charge forces - sbcov initialization
        //     ------------------------------------------------
        if (sqrt(rsq) > intparams[m].dc2) continue;

        attractive (&intparams[m],rsq,eflag,i,iq,j,jq);


      } // ---------------------------------------------- for jj


    } // -------------------------------------------- ii

    // ===============================================
    //    Communicates the tables *sbcov and *sbmet
    //    to calculate the N-body forces
    // ================================================

    forward (sbcov) ; // reverse (sbcov);
    forward (coord) ; // reverse (coord);
    forward (sbmet) ; // reverse (sbmet);


    if (nteam == 0) return; //no oxide
    if (Qstep == 0 || (step % Qstep != 0)) return;

    // =======================
  } // End of If(loop)
    // =======================


  // ===============================================

  for (ii=0; ii<inum; ii++)
    {
      i = ilist[ii]; itype = map[type[i]];
      gp = flag_QEq[i]; iq = q[i];

      potmad[i] = potself[i] = potcov[i] =  0.0 ;

      if (gp == 0) continue;

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      fqi = 0.0 ;

      //  Madelung potential
      // --------------------
      potmad[i] += 2.0*Vself*iq ;

      // charge force from self energy
      // -----------------------------
      fqi = qfo_self(&params[itype],iq);
      potself[i] = fqi ;



      // Loop on Second moment neighbor
      // -------------------------------

      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++)
        {
          j = jlist[jj];
          j &= NEIGHMASK;
          jtype = map[type[j]];
          m = intype[itype][jtype];
          jq = q[j];


          const double delx = x[j][0] - xtmp;
          const double dely = x[j][1] - ytmp;
          const double delz = x[j][2] - ztmp;
          rsq = delx*delx + dely*dely + delz*delz;

          // long range q-dependent
          if (sqrt(rsq) > cutmax) continue;

          //   1/r charge forces
          //  --------------------
          fqij = 0.0;
          //                  pot_ES2 (i,j,rsq,fqij2);
          pot_ES (i,j,rsq,fqij);

          potmad[i] += jq*fqij ;


        } // ------ jj

      fqi = 0.0;
      pot_COV (&params[itype],i,fqi) ;
      potcov[i] = fqi ;


    } // ------- ii


}

/* ---------------------------------------------------------------------- */

void PairSMTBQ::Charge()
{
  int i,ii,iloop,itype,gp,m;
  int *ilist;
  double heatpq,qmass,dtq,dtq2,qtot,qtotll;
  double t_init,t_end,dt;

  double *Transf,*TransfAll;

  double *q = atom->q;
  double **x = atom->x;
  int *type = atom->type;
  int step = update->ntimestep;

  int inum = list->inum;
  ilist = list->ilist;


  if (me == 0) t_init = MPI_Wtime();
  if (step == 0) cluster = 0;

  // ---------------------------
  //  Mise en place des groupes
  // ---------------------------

  if (strcmp(QEqMode,"BulkFromSlab") == 0)
    groupBulkFromSlab_QEq();
  else if (strcmp(QEqMode,"QEqAll") == 0)
    groupQEqAll_QEq();
  else if (strcmp(QEqMode,"QEqAllParallel") == 0)
    groupQEqAllParallel_QEq();
  else if (strcmp(QEqMode,"Surface") == 0)
    groupSurface_QEq();



  if (nteam+1 != cluster) {
    memory->destroy(nQEqall);
    memory->destroy(nQEqaall);
    memory->destroy(nQEqcall);

    cluster = nteam+1;
    memory->create(nQEqall,nteam+1,"pair:nQEq");
    memory->create(nQEqaall,nteam+1,"pair:nQEqa");
    memory->create(nQEqcall,nteam+1,"pair:nQEqc");
  }
  // ---------------------------


  double enegtotall[nteam+1],enegchkall[nteam+1],enegmaxall[nteam+1],qtota[nteam+1],qtotc[nteam+1];
  double qtotcll[nteam+1],qtotall[nteam+1];
  double sigmaa[nteam+1],sigmac[nteam+1],sigmaall[nteam+1],sigmacll[nteam+1];
  int end[nteam+1], nQEq[nteam+1],nQEqc[nteam+1],nQEqa[nteam+1];


  iloop = 0;


  heatpq = 0.07;
  qmass  = 0.000548580;
  dtq    = 0.0006; // 0.0006
  dtq2   = 0.5*dtq*dtq/qmass;

  double enegchk[nteam+1];
  double enegtot[nteam+1];
  double enegmax[nteam+1];



  for (i=0; i<nteam+1; i++) {
    nQEq[i] = nQEqa[i] = nQEqc[i] = 0;
    nQEqall[i] = nQEqcall[i] = nQEqaall[i] = end[i]= 0;
    enegchk[i] = enegtot[i] = enegmax[i] = 0.0;
    qtota[i] = qtotc[i] = qtotall[i] = qtotcll[i] = 0.0;
    sigmaall[i] = sigmacll[i] = 0.0;
  }
  qtot = qtotll = 0.0 ;


  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii]; itype = map[type[i]];
    gp = flag_QEq[i];
    q1[i] = q2[i] = qf[i] = 0.0;

    qtot += q[i] ;
    nQEq[gp] += 1;
    if (itype != 0) { qtotc[gp] += q[i]; nQEqc[gp] += 1; }
    if (itype == 0) { qtota[gp] += q[i]; nQEqa[gp] += 1; }
  }

  MPI_Allreduce(nQEq,nQEqall,nteam+1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(nQEqc,nQEqcall,nteam+1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(nQEqa,nQEqaall,nteam+1,MPI_INT,MPI_SUM,world);

  MPI_Allreduce(qtotc,qtotcll,nteam+1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(qtota,qtotall,nteam+1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&qtot,&qtotll,1,MPI_DOUBLE,MPI_SUM,world);


  //  ---------------------------------
  //     Init_charge(nQEq,nQEqa,nQEqc);
  //  ---------------------------------

  if (update->ntimestep == 0 && (strcmp(QInitMode,"true") == 0)   ) {
    //Carefull here it won't be fine if there are more than 2 species!!!
    QOxInit=max(QOxInit, -0.98* params[1].qform *nQEqcall[gp]/nQEqaall[gp])   ;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii]; itype = map[type[i]];
      gp = flag_QEq[i];

      //    if (gp == 0) continue;

      if (itype == 0) q[i] = QOxInit ;
      if (itype > 0) q[i] = -QOxInit * nQEqaall[gp]/nQEqcall[gp];
    }
  }

  if (nteam == 0 || Qstep == 0) return;
  if (step % Qstep != 0) return;
  //  --------------------------------------

  // ----

  // ----


  if (me == 0 && strcmp(Bavard,"false") != 0) {
    for (gp = 0; gp < nteam+1; gp++) {
      printf (" ++++ Groupe %d - Nox %d Ncat %d\n",gp,nQEqaall[gp],nQEqcall[gp]);
      if (nQEqcall[gp] !=0 && nQEqaall[gp] !=0 )
        printf (" neutralite des charges %f\n qtotc %f qtota %f\n",
                qtotll,qtotcll[gp]/nQEqcall[gp],qtotall[gp]/nQEqaall[gp]);
      printf (" ---------------------------- \n");}
  }

  // ======================= Tab transfert ==================
  //  Transf[gp] = enegtot[gp]
  //  Transf[gp+cluster] = Qtotc[gp]
  //  Transf[gp+2cluster] = Qtota[gp]
  //  Transf[3cluster] = Qtot
  //  -------------------------------------------------------
  memory->create(Transf,3*cluster+1,"pair:Tranf");
  memory->create(TransfAll,3*cluster+1,"pair:TranfAll");
  // ========================================================


  // --------------------------------------------
  for (iloop = 0; iloop < loopmax; iloop ++ ) {
    // --------------------------------------------

    qtot = qtotll = Transf[3*cluster] = 0.0 ;
    for (gp=0; gp<nteam+1; gp++) {
      Transf[gp] = Transf[gp+cluster] = Transf[gp+2*cluster] = 0.0;
      TransfAll[gp] = TransfAll[gp+cluster] = TransfAll[gp+2*cluster] = 0.0;
      enegtot[gp] = enegtotall[gp] = enegchkall[gp] = enegmaxall[gp] = 0.0 ;
    }

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      q1[i] += qf[i]*dtq2 - heatpq*q1[i];
      q[i]  += q1[i];

      Transf[3*cluster] += q[i];
      itype = map[type[i]];
      gp = flag_QEq[i];
      if (itype == 0) Transf[gp+2*cluster] += q[i];
      if (itype != 0) Transf[gp+cluster] += q[i];
    }

    //   Communication des charges
    //  ---------------------------
    forward(q) ; // reverse(q);


    //   Calcul des potential
    //  ----------------------
    QForce_charge(iloop);

    //     exit(1);

    for (ii = 0; ii < inum; ii++)
      {
        i = ilist[ii]; itype = map[type[i]];
        gp = flag_QEq[i];

        qf[i] = 0.0;
        qf[i] = potself[i]+potmad[i]+potcov[i]+chimet[i] ;
        Transf[gp] += qf[i];
      }

    MPI_Allreduce(Transf,TransfAll,3*cluster+1,MPI_DOUBLE,MPI_SUM,world);

    for (i = 0; i < nteam+1; i++) {

      if (nQEqall[i] !=0) TransfAll[i] /= static_cast<double>(nQEqall[i]);
      enegchk[i] = enegmax[i] = 0.0;
    }

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      itype = map[type[i]];
      gp = flag_QEq[i];

      if (gp == 0) continue;

      q2[i] = TransfAll[gp] - qf[i];

      enegmax[gp] = MAX(enegmax[gp],fabs(q2[i]));
      enegchk[gp] += fabs(q2[i]);
      qf[i] = q2[i];

    } // Boucle local

    MPI_Allreduce(enegchk,enegchkall,nteam+1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(enegmax,enegmaxall,nteam+1,MPI_DOUBLE,MPI_MAX,world);


    for (gp = 0; gp < nteam+1; gp++) {
      if(nQEqall[gp] !=0) {
        enegchk[gp] = enegchkall[gp]/static_cast<double>(nQEqall[gp]);
        enegmax[gp] = enegmaxall[gp];
      }
    }

    // -----------------------------------------------------
    //                Convergence Test
    // -----------------------------------------------------

    m = 0;
    for (gp = 1; gp < nteam+1; gp++) {
      if (enegchk[gp] <= precision && enegmax[gp] <= 100.0*precision) end[gp] =  1;
    }
    for (gp = 1; gp < nteam+1; gp++) { m += end[gp] ; }

    if (m == nteam) {        break;  }
    // -----------------------------------------------------
    // -----------------------------------------------------

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      q1[i] += qf[i]*dtq2 - heatpq*q1[i];
    }

    // --------------------------------------------
  } // -------------------------------- iloop
    // --------------------------------------------


  // =======================================
  //    Charge Communication.
  // =======================================
  forward(q); // reverse(q);

  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ==========================================
  //   Ecriture des potentials dans un fichier
  // ==========================================

  if (strcmp(writepot,"true") == 0 && fmod(static_cast<double>(step), Neverypot) == 0.0) {

    ofstream fichierpot("Electroneg_component.txt", ios::out | ios::trunc) ;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      itype = map[type[i]];
      gp = flag_QEq[i];


      if (fichierpot) fichierpot<< setprecision(9) <<i <<" "<<itype<<" "<<x[i][0]<<" "<<x[i][1]
                                <<" "<<x[i][2]<<" "<<q[i]<<" "<<potself[i] + potmad[i]<<" "<<potcov[i]
                                <<" "<<sbcov[i]<<" "<<TransfAll[gp]<<endl;

    }
    if (fichierpot) fichierpot.close() ;
  }

  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


  //   Statistique (ecart type)
  //   ------------------------
  for (i=0; i<nteam+1; i++) {
    if(nQEqcall[i] !=0)
      { TransfAll[i+cluster] /= static_cast<double>(nQEqcall[i]) ;
        TransfAll[i+2*cluster] /= static_cast<double>(nQEqaall[i]) ;}
    sigmaa[i] = sigmac[i] = 0.0;
  }

  qtot = 0.0 ;
  for (ii = 0; ii < inum; ii++)
    {
      i = ilist[ii];
      itype = map[type[i]];
      gp = flag_QEq[i];
      //      qtot += q[i];

      if (gp == 0) continue;

      if (itype == 0) sigmaa[gp] += (q[i]-TransfAll[gp+2*cluster])*(q[i]-TransfAll[gp+2*cluster]);
      if (itype == 1) sigmac[gp] += (q[i]-TransfAll[gp+cluster])*(q[i]-TransfAll[gp+cluster]);
    }

  MPI_Allreduce(sigmaa,sigmaall,nteam+1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(sigmac,sigmacll,nteam+1,MPI_DOUBLE,MPI_SUM,world);

  for (gp = 1; gp < nteam+1; gp++) {
    sigmaall[gp] = sqrt(sigmaall[gp]/static_cast<double>(nQEqaall[gp])) ;
    sigmacll[gp] = sqrt(sigmacll[gp]/static_cast<double>(nQEqcall[gp])) ;
  }



  if (me == 0 && strcmp(Bavard,"false") != 0){
    for (gp = 0; gp < nteam+1; gp++) {
      printf (" -------------- Groupe %d -----------------\n",gp);
      printf (" qtotc %f(+- %f) qtota %f(+- %f)\n",
              TransfAll[gp+cluster],sigmacll[gp],TransfAll[gp+2*cluster],sigmaall[gp]);
      printf (" Potentiel elec total : %f\n iloop %d, qtot %f\n",TransfAll[gp],iloop,TransfAll[3*cluster]);
      printf (" convergence : %f - %f\n",enegchk[gp],enegmax[gp]);
    }

    t_end = MPI_Wtime();
    dt = t_end - t_init;
    printf (" temps dans charges : %f seconde. \n",dt);
    printf (" ======================================================== \n");
  }

  // ============== Destroy Tab
  memory->destroy(Transf);
  memory->destroy(TransfAll);

}

/* ---------------------------------------------------------------------- */

void PairSMTBQ::groupBulkFromSlab_QEq()
{ int ii,i;
  double **x=atom->x;
  int *ilist;
  double ztmp;
  int inum = list->inum;
  ilist = list->ilist;

  for (ii = 0; ii < inum; ii++)
    {
      i = ilist[ii];
      ztmp = x[i][2];
      if (ztmp>zlim1QEq && ztmp< zlim2QEq)
        flag_QEq[i]=1;
      else
        flag_QEq[i]=0;

      nteam=1;

    }

}

// ----------------------------------------------

void PairSMTBQ::groupQEqAll_QEq()
{ int ii,i;
  int *ilist;
  int inum = list->inum;
  ilist = list->ilist;

  nteam=1;

  for (ii = 0; ii < inum; ii++)
    {
      i= ilist[ii];
      flag_QEq[i]=1;
    }

}

// ----------------------------------------------

void PairSMTBQ::groupSurface_QEq()
{ int ii,i;
  double **x=atom->x;
  int *ilist;
  double ztmp;
  int inum = list->inum;
  ilist = list->ilist;

  for (ii = 0; ii < inum; ii++)
    {
      i = ilist[ii];
      ztmp = x[i][2];
      if (ztmp>zlim1QEq)
        flag_QEq[i]=1;
      else
        flag_QEq[i]=0;

      nteam=1;

    }

}


void PairSMTBQ::groupQEqAllParallel_QEq()
{
  int ii,i,jj,j,kk,k,itype,jtype,ktype,jnum,m,gp,zz,z,kgp;
  int iproc,team_elt[10][nproc],team_QEq[10][nproc][5];
  int *ilist,*jlist,*numneigh,**firstneigh,ngp,igp;
  double delr[3],xtmp,ytmp,ztmp,rsq;
  int **flag_gp, *nelt, **tab_gp;
  int QEq,QEqall[nproc];

  double **x = atom->x;
  int *type = atom->type;
  const int nlocal = atom->nlocal;
  const int nghost = atom->nghost;
  const int nall = nlocal + nghost;
  int inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;


  // +++++++++++++++++++++++++++++++++++++++++++++++++
  //  On declare et initialise nos p'tits tableaux
  // +++++++++++++++++++++++++++++++++++++++++++++++++

  int **tabtemp,**Alltabtemp, *gptmp, *Allgptmp;

  memory->create(tabtemp,10*nproc+10,nproc,"pair:tabtemp");
  memory->create(Alltabtemp,10*nproc+10,nproc,"pair:Alltabtemp");
  memory->create(gptmp,10*nproc+10,"pair:gptmp");
  memory->create(Allgptmp,10*nproc+10,"pair:Allgptmp");

  memory->create(flag_gp,nproc,nall,"pair:flag_gp");
  memory->create(nelt,nall,"pair:nelt");
  memory->create(tab_gp,10,nall,"pair:flag_gp");


  for (i = 0; i < nall ; i++) { flag_QEq[i] = 0; }
  for (i = 0; i < 10*nproc; i++) {
    gptmp[i] = 0; Allgptmp[i] = 0;
    for (j=0;j<nproc;j++) { tabtemp[i][j] = 0;
      Alltabtemp[i][j] = 0;}
  }
  for (i = 0; i < 10; i++) {
    for (k = 0; k < nall; k++) { tab_gp[i][k] = 0;
      if (i == 0) nelt[k] = 0;
    }
    for (j = 0; j < nproc; j++) {
      team_elt[i][j] = 0;
      for (k = 0; k < 5; k++) { team_QEq[i][j][k] = 0; }
    }
  }

  QEq = 0;


  //   printf ("groupeQEq me %d - nloc %d nghost %d boite %d\n",
  //             me,nlocal,nghost,nall);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  On identifie les atomes rentrant dans le schema QEq +
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for (ii = 0; ii < inum; ii++)
    {
      i = ilist[ii] ; itype = map[type[i]] ;

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];


      jlist = firstneigh[i];
      jnum = numneigh[i];
      for (jj = 0; jj < jnum; jj++ )
        {
          j = jlist[jj] ;
          j &= NEIGHMASK;
          jtype = map[type[j]];
          if (jtype == itype) continue;
          m = intype[itype][jtype];

          delr[0] = x[j][0] - xtmp;
          delr[1] = x[j][1] - ytmp;
          delr[2] = x[j][2] - ztmp;
          rsq = vec3_dot(delr,delr);

          if (sqrt(rsq) <= intparams[m].dc2) {
            flag_QEq[i] = 1; flag_QEq[j] = 1;
          }
        }
      if (flag_QEq[i] == 1) {
        QEq = 1;
      }
    }

  // ::::::::::::::: Testing the presence of oxide :::::::::::
  m = 0;
  MPI_Allgather(&QEq,1,MPI_INT,QEqall,1,MPI_INT,world);
  for (iproc = 0; iproc < nproc; iproc++) {
    if (QEqall[iproc] == 1) m = 1;
  }
  if (m == 0) {
    memory->destroy(tabtemp);
    memory->destroy(Alltabtemp);
    memory->destroy(gptmp);
    memory->destroy(Allgptmp);
    memory->destroy(flag_gp);
    memory->destroy(tab_gp);
    memory->destroy(nelt);

    return;
  }
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::


  for (m = 0; m < nproc; m++) {
    for (i = 0; i < nall; i++) { flag_gp[m][i] = 0; }
  }

  // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  //  It includes oxygens entering the QEq scheme          O
  // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO


  ngp = igp = 0; nelt[ngp] = 0;

  // On prend un oxygène
  //   printf ("[me %d] On prend un oxygene\n",me);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii] ; itype = map[type[i]];
    if (itype != 0 || flag_QEq[i] == 0) continue;

    m = 0;

    if (ngp != 0 && flag_gp[me][i] == ngp) continue;


    //   Grouping Initialisation
    //  -----------------------------
    if (flag_gp[me][i] == 0) {
      ngp += 1; nelt[ngp] = 0;
      tab_gp[ngp][nelt[ngp]] = i;
      flag_gp[me][i] = ngp;
      nelt[ngp] += 1;
    }
    //  -------------------------------


    //       Loop on the groups
    //      ----------------------
    for (kk = 0; kk < nelt[ngp]; kk++)
      {
        k = tab_gp[ngp][kk];
        ktype = map[type[k]];
        //              printf ("[me %d] kk - gp %d elemt %d : atom %d(%d)\n",me,ngp,kk,k,ktype);
        if (k >= nlocal) continue;

        xtmp = x[k][0];
        ytmp = x[k][1];
        ztmp = x[k][2];

        //  Loop on the oxygen's neighbor of the group
        //  ---------------------------------------------
        jlist = firstneigh[k];
        jnum = numneigh[k];
        for (j = 0; j < nall; j++ )
          {
            jtype = map[type[j]];
            if (jtype == ktype) continue;
            m = intype[itype][jtype];

            if (jtype == 0 && flag_QEq[j] == 0) continue;


            if (flag_gp[me][j] == ngp)  continue;

            delr[0] = x[j][0] - xtmp;
            delr[1] = x[j][1] - ytmp;
            delr[2] = x[j][2] - ztmp;
            rsq = vec3_dot(delr,delr);

            //     -------------------------------------
            if (sqrt(rsq) <= cutmax) {

              flag_QEq[j] = 1; //Entre dans le schema QEq

              //   :::::::::::::::::::: Meeting of two group in the same proc :::::::::::::::::::::

              if (flag_gp[me][j] != 0 && flag_gp[me][j] != ngp && nelt[flag_gp[me][j]] != 0) {
                printf("[me %d] (atom %d) %d [elt %d] rencontre un nouveau groupe %d [elt %d] (atom %d)\n",
                       me,k,ngp,nelt[ngp],flag_gp[me][j],nelt[flag_gp[me][j]],j);

                //         On met a jours les tableaux
                //        -----------------------------
                igp = flag_gp[me][j];
                z = min(igp,ngp);

                if (z == igp) { igp = z; }
                else if (z == ngp) {
                  ngp = igp ; igp = z;
                  flag_gp[me][j] = ngp;
                }

                for (zz = 0; zz < nelt[ngp]; zz++) {
                  z = tab_gp[ngp][zz];
                  tab_gp[igp][nelt[igp]] = z;
                  nelt[igp] += 1;
                  flag_gp[me][z] = igp;
                  tab_gp[ngp][zz] = 0;
                }

                nelt[ngp] = 0;
                for (z = nlocal; z < nall; z++) {
                  if (flag_gp[me][z] == ngp) flag_gp[me][z] = igp;
                }

                m = 1; kk = 0;
                ngp = igp;
                break;
              }
              //   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

              flag_gp[me][j] = ngp;
              if (j < nlocal)
                {
                  tab_gp[ngp][nelt[ngp]] = j;
                  nelt[ngp] += 1;
                }
            }
          } // for j
      } // for k
  } // for ii

  // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  //   Groups communication
  // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  for (i = 0; i < nproc; i++) {
    forward_int(flag_gp[i]); //reverse_int(flag_gp[i]);
  }
  // ---


  // =======================================================
  //  Loop on the cation to make them joined in the oxygen's
  //  group which it interacts
  // =======================================================
  igp = 0;
  for (ii = 0; ii < inum; ii++)
    {
      i = ilist[ii] ; itype = map[type[i]];
      if (itype == 0) continue;

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      jlist = firstneigh[i];
      jnum = numneigh[i];
      for (jj = 0; jj < jnum; jj++ )
        {
          j = jlist[jj] ;
          j &= NEIGHMASK;
          jtype = map[type[j]];
          if (jtype != 0) continue;

          m = 0;
          for (iproc = 0; iproc < nproc; iproc++) {
            if (flag_gp[iproc][j] != 0) m = flag_gp[iproc][j];
          }
          if (m == 0) continue;

          delr[0] = x[j][0] - xtmp;
          delr[1] = x[j][1] - ytmp;
          delr[2] = x[j][2] - ztmp;
          rsq = vec3_dot(delr,delr);

          //    ----------------------------------------
          if (sqrt(rsq) <= cutmax) {
            //       if (sqrt(rsq) <= intparams[m].dc2) {
            //    ----------------------------------------

            flag_QEq[i] = 1; igp = flag_gp[me][j];

            if (flag_gp[me][i] == 0) flag_gp[me][i] = igp;

            if (flag_gp[me][i] != igp && igp != 0) {
              printf ("[me %d] Cation i %d gp %d [nelt %d] rencontre j %d(%d)du groupe %d [nelt %d]\n",
                      me,i,flag_gp[me][i],nelt[flag_gp[me][i]],j,jtype,igp,nelt[igp]);

              igp = min(flag_gp[me][i],flag_gp[me][j]);
              if (igp == flag_gp[me][i]) { kgp = flag_gp[me][j]; }
              else { kgp = flag_gp[me][i]; }

              for (k = 0; k < nelt[kgp]; k++) {
                z = tab_gp[kgp][k];
                tab_gp[igp][nelt[igp]] = z;
                nelt[igp] += 1;
                flag_gp[me][z] = igp;
                tab_gp[kgp][k] = 0;
              }
              nelt[kgp] = 0;

              for (k = 0; k < nall; k++) {
                if (flag_gp[me][k] == kgp) flag_gp[me][k] = igp;
              }

            }
            m = 0;
            for (k = 0; k < nelt[igp]; k++) {
              if (tab_gp[igp][k] == i) m = 1;
            }

            if (i >= nlocal || m == 1 ) continue;
            //          printf ("[me %d] igp %d - nelt %d atom %d\n",me,igp,nelt[igp],i);
            tab_gp[igp][nelt[igp]] = i;
            nelt[igp] += 1;
            break;
          }

        } // voisin j

    } // atom i

  /* ==================================================
     Group Communication between proc for unification
     ================================================== */
  for (i = 0; i < nproc; i++) {
    forward_int(flag_gp[i]);// reverse_int(flag_gp[i]);
  }

  //  =============== End of COMM =================


  for (i = 0; i < nall; i++) {

    m = 10*me + flag_gp[me][i];
    if (m == 10*me) continue; // Pas de groupe zero
    gptmp[m] = 1;
    for (k = 0; k < nproc; k++) {

      if (k == me) continue;
      if (tabtemp[m][k] != 0) continue;

      if (flag_gp[k][i] != 0) {
        tabtemp[m][k] = 10*k + flag_gp[k][i];
      }
    }

  }

  for (k = 0; k < 10*nproc; k++) {
    MPI_Allreduce(tabtemp[k],Alltabtemp[k],nproc,MPI_INT,MPI_SUM,world); }
  MPI_Allreduce(gptmp,Allgptmp,10*nproc,MPI_INT,MPI_SUM,world);

  nteam = 0; iproc = 0;
  for (igp = 0; igp < 10*nproc; igp++) {
    if (Allgptmp[igp] == 0) continue;
    iproc = int(static_cast<double>(igp)/10.0);
    ngp = igp - 10*iproc;
    if (nteam == 0) {

      nteam += 1;
      team_elt[nteam][iproc] = 0;
      team_QEq[nteam][iproc][team_elt[nteam][iproc]] = ngp;
      team_elt[nteam][iproc] += 1;
    } else {
      m = 0;
      for (i = 1; i < nteam+1; i++) {
        for (k = 0; k < team_elt[i][iproc]; k++) {
          if (ngp == team_QEq[i][iproc][k]) m = 1;
        } }
      if (m == 1) continue;
      //         create a new team!!
      //      ---------------------------
      if (m == 0) {
        nteam += 1;
        team_elt[nteam][iproc] = 0;
        team_QEq[nteam][iproc][team_elt[nteam][iproc]] = ngp;
        team_elt[nteam][iproc] += 1;
      }
    }
    //  -------
    //   On a mis une graine dans le groupe et nous allons
    //  remplir se groupe en questionnant le "tabtemp[igp][iproc]=jgp"
    //
    for (kk = 0; kk < nproc; kk++) {
      for (k = 0; k < team_elt[nteam][kk]; k++) {

        // On prend le gp le k ieme element de la team nteam sur le proc iproc
        //         ngp = 0;
        ngp = team_QEq[nteam][kk][k];
        kgp = 10*kk + ngp;

        // On regarde sur les autre proc si ce gp ne pointe pas vers un autre.
        for (i = 0; i < nproc; i++) {
          if (i == kk) continue;
          if (Alltabtemp[kgp][i] == 0) continue;

          if (Alltabtemp[kgp][i] != 0) ngp = Alltabtemp[kgp][i];
          ngp = ngp - 10*i;

          //  Est ce que ce groupe est nouveau?
          m = 0;
          for (j = 0; j < team_elt[nteam][i]; j++) {
            if (team_QEq[nteam][i][j] == ngp) m = 1;
          }

          if (m == 0) {
            iproc = i; k = 0;
            team_QEq[nteam][i][team_elt[nteam][i]] = ngp ;
            team_elt[nteam][i] += 1;
          }
        } // regard sur les autre proc

      } // On rempli de proche en proche
    } // boucle kk sur les proc
  }

  //  Finalement on met le numero de la team en indice du flag_QEq, c mieu!
  //  ---------------------------------------------------------------------

  for (ii = 0; ii < inum; ii++)
    {
      i = ilist[ii]; m = 0; itype = map[type[i]];
      if (flag_QEq[i] == 0) continue;

      gp = flag_gp[me][i];
      for (j = 1; j <= nteam; j++) {
        for (k = 0; k < team_elt[j][me]; k++) {
          if (gp == team_QEq[j][me][k]) {
            flag_QEq[i] = j; m = 1;
            break;
          }
        }
        if (m == 1) break;
      }
    }

  memory->destroy(tabtemp);
  memory->destroy(Alltabtemp);
  memory->destroy(gptmp);
  memory->destroy(Allgptmp);
  memory->destroy(flag_gp);
  memory->destroy(tab_gp);
  memory->destroy(nelt);

}

/* ---------------------------------------------------------------------- */

void PairSMTBQ::Init_charge(int * /*nQEq*/, int * /*nQEqa*/, int * /*nQEqc*/)
{
  int ii,i,gp,itype;
  int *ilist,test[nteam],init[nteam];
  double bound,tot,totll;

  int inum = list->inum;
  int *type = atom->type;
  double *q = atom->q;
  ilist = list->ilist;

  if (nteam == 0) return;

  if (me == 0) printf (" ======== Init_charge ======== \n");
  for (gp = 0; gp < cluster; gp++) {
    test[gp] = 0; init[gp] = 0;
  }


  //  On fait un test sur les charges pour voir sont
  //  elles sont dans le domaine delimiter par DeltaQ
  // -------------------------------------------------
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii]; itype = map[type[i]];
    gp = flag_QEq[i];

    if (gp != 0 && itype == 0) {
      bound = fabs(2.0*ncov/params[itype].sto - fabs(params[itype].qform)) ;
      if (bound == fabs(params[itype].qform)) continue;
      if (fabs(q[i]) < bound) test[gp] = 1;
    }
  }

  MPI_Allreduce(test,init,nteam+1,MPI_INT,MPI_SUM,world);

  //  On fait que sur les atomes hybrides!!!
  // ----------------------------------------

  tot = totll = 0.0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii]; itype = map[type[i]];
    gp = flag_QEq[i];

    if (gp != 0 && init[gp] != 0) {
      if (itype == 0) q[i] = -1.96;
      if (itype != 0) q[i] = 1.96 * static_cast<double>(nQEqaall[gp]) / static_cast<double>(nQEqcall[gp]);
    }
    tot += q[i];
  }
  MPI_Allreduce(&tot,&totll,1,MPI_INT,MPI_SUM,world);
  if (me == 0) printf (" === Fin de init_charge qtot %20.15f ====\n",totll);

}
/* ----------------------------------------------------------------------
 *                        COMMUNICATION
 * ---------------------------------------------------------------------- */

int PairSMTBQ::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i ++) {
    j = list[i];
    buf[m++] = tab_comm[j];
    //    if (j < 3) printf ("[%d] %d pfc %d %d buf_send = %f \n",me,n,i,m-1,buf[m-1]);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSMTBQ::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n ;
  for (i = first; i < last; i++) {
    tab_comm[i] = buf[m++];
    //    if (i<first+3) printf ("[%d] ufc %d %d buf_recv = %f \n",me,i,m-1,buf[m-1]);
  }
}

/* ---------------------------------------------------------------------- */

int PairSMTBQ::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = tab_comm[i];
    //    if (i<first+3) printf ("[%d] prc %d %d buf_send = %f \n",me,i,m-1,buf[m-1]);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSMTBQ::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    //  tab_comm[j] += buf[m++];
    tab_comm[j] = buf[m++];
    //    if (j<3) printf ("[%d] %d urc %d %d buf_recv = %f \n",me,n,i,m-1,buf[m-1]);
  }
}

/* ---------------------------------------------------------------------- */

void PairSMTBQ::forward(double *tab)
{
  int i;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  for (i=0; i<nlocal+nghost; i++) tab_comm[i] = tab[i];

  comm->forward_comm_pair(this);

  for (i=0; i<nlocal+nghost; i++) tab[i] = tab_comm[i];
}

/* ---------------------------------------------------------------------- */

void PairSMTBQ::reverse(double *tab)
{
  int i;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  for (i=0; i<nlocal+nghost; i++)  tab_comm[i] = tab[i];

  comm->reverse_comm_pair(this);

  for (i=0; i<nlocal+nghost; i++)  tab[i] = tab_comm[i];
}

/* ---------------------------------------------------------------------- */

void PairSMTBQ::forward_int(int *tab)
{
  int i;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  for (i=0; i<nlocal+nghost; i++) { tab_comm[i] = static_cast<double>(tab[i]);}

  comm->forward_comm_pair(this);

  for (i=0; i<nlocal+nghost; i++) {
    if (fabs(tab_comm[i]) > 0.1) tab[i] = int(tab_comm[i]) ; }
}

/* ---------------------------------------------------------------------- */

void PairSMTBQ::reverse_int(int *tab)
{
  int i;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  for (i=0; i<nlocal+nghost; i++) { tab_comm[i] = static_cast<double>(tab[i]);}

  comm->reverse_comm_pair(this);

  for (i=0; i<nlocal+nghost; i++) {
    if (fabs(tab_comm[i]) > 0.1) tab[i] = int(tab_comm[i]); }
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

double PairSMTBQ::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += nmax * sizeof(int);
  bytes += MAXNEIGH * nmax * sizeof(int);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void PairSMTBQ::add_pages(int howmany)
{
  int toppage = maxpage;
  maxpage += howmany*PGDELTA;

  pages = (int **)
    memory->srealloc(pages,maxpage*sizeof(int *),"pair:pages");
  for (int i = toppage; i < maxpage; i++)
    memory->create(pages[i],pgsize,"pair:pages[i]");
}

/* ---------------------------------------------------------------------- */


int PairSMTBQ::Tokenize( char* s, char*** tok )
{
  char test[MAXLINE];
  const char *sep = "' ";
  char *mot;
  int count=0;
  mot = NULL;


  strncpy( test, s, MAXLINE-1 );

  for( mot = strtok(test, sep); mot; mot = strtok(NULL, sep) ) {
    strncpy( (*tok)[count], mot, MAXLINE );
    count++;
  }

  return count;
}



void PairSMTBQ::CheckEnergyVSForce()
{
  double drL,iq,jq,rsq,evdwlCoul,fpairCoul,eflag=0,ErepR,frepR,fpair,evdwl;
  int i,j,iiiMax,iii,iCoord;
  int itype,jtype,l,m;
  double r,t1,t2,sds,xi,engSurf,fforceSurf;
  double eng,fforce,engBB,fforceBB;

  double za,zb,gam,dgam,dza,dzb,
    d2zaa,d2zab,d2zbb,d2zra,d2zrb,d2gamr2,na,nb;
  int *type = atom->type;
  const char *NameFile;


  i=0;
  j=1;
  map[type[i]]=0;  //ox
  itype=map[type[i]];
  iq=-1;


  map[type[j]]=0;  //ox
  jtype=map[type[j]];
  jq=-1;
  coord[i]=coordOxBulk;
  coord[j]=coordOxBulk;
  m = intype[itype][jtype];


  na = params[itype].ne ;
  nb = params[jtype].ne ;
  za = params[itype].dzeta ;
  zb = params[jtype].dzeta ;


  //   Ouverture du fichier

  for (iCoord=1;iCoord<5; iCoord++)
    {

      if (iCoord==1)
        {coord[i]=2.2;
          coord[j]=2.1;
          NameFile=(const char *)"energyandForceOxOxUnderCoord.txt";
        }
      if (iCoord==2)
        {coord[i]=coordOxBulk;
          coord[j]=coordOxBulk;
          NameFile=(const char *)"energyandForceOxOxCoordBulk.txt";
        }
      if (iCoord==3)
        {coord[i]=3.8;
          coord[j]=4;
          NameFile=(const char *)"energyandForceOxOxOverCoord.txt";
        }
      if (iCoord==4)
        {coord[i]=2.2;
          coord[j]=3.5;
          NameFile=(const char *)"energyandForceOxOxUnderOverCoord.txt";
        }


      ofstream fichierOxOx(NameFile, ios::out | ios::trunc) ;

      drL=0.0001;
      iiiMax=int((cutmax-1.2)/drL);
      for (iii=1; iii< iiiMax ; iii++){
        r=1.2+drL*iii;
        rsq=r*r;
        evdwlCoul = 0.0 ; fpairCoul = 0.0;
        potqeq(i,j,iq,jq,rsq,fpairCoul,eflag,evdwlCoul);
        fpairCoul=fpairCoul*r;

        rep_OO (&intparams[m],rsq,fpair,eflag,evdwl);
        ErepR = evdwl;
        frepR= fpair*r;

        gam = dgam = dza = dzb = d2zaa = d2zab =
          d2zbb = d2zra = d2zrb = d2gamr2 = 0.0 ;


        //        gammas_(na,nb,za,zb,r,gam,dgam,dza,dzb,
        //                d2zaa,d2zab,d2zbb,d2zra,d2zrb,d2gamr2) ;
        gammas(na,nb,za,zb,r,gam,dgam,dza,dzb,
               d2zaa,d2zab,d2zbb,d2zra,d2zrb,d2gamr2) ;


        sds = rsq/ds ;  l = int(sds) ;
        xi = sds - static_cast<double>(l) ;


        t1 = fafb[l][m] + (fafb[l+1][m] - fafb[l][m])*xi;
        t2 = fafb[l+1][m] + (fafb[l+2][m] - fafb[l+1][m])*(xi-1.0);
        eng = iq*jq*(t1 + (t2 - t1)*xi/2.0);

        t1 = dfafb[l][m] + (dfafb[l+1][m] - dfafb[l][m])*xi;
        t2 = dfafb[l+1][m] + (dfafb[l+2][m] - dfafb[l+1][m])*(xi-1);
        fforce = - iq*jq*(t1 + (t2 - t1)*xi/2.0)*r ;


        t1 = fafbOxOxSurf[l] + (fafbOxOxSurf[l+1] - fafbOxOxSurf[l])*xi;
        t2 = fafbOxOxSurf[l+1] + (fafbOxOxSurf[l+2] - fafbOxOxSurf[l+1])*(xi-1.0);
        engSurf = iq*jq*(t1 + (t2 - t1)*xi/2.0);

        t1 = fafbOxOxBB[l] + (fafbOxOxBB[l+1] - fafbOxOxBB[l])*xi;
        t2 = fafbOxOxBB[l+1] + (fafbOxOxBB[l+2] - fafbOxOxBB[l+1])*(xi-1.0);
        engBB = iq*jq*(t1 + (t2 - t1)*xi/2.0);

        t1 = dfafbOxOxSurf[l] + (dfafbOxOxSurf[l+1] - dfafbOxOxSurf[l])*xi;
        t2 = dfafbOxOxSurf[l+1] + (dfafbOxOxSurf[l+2] - dfafbOxOxSurf[l+1])*(xi-1);
        fforceSurf = - iq*jq*(t1 + (t2 - t1)*xi/2.0)*r ;

        t1 = dfafbOxOxBB[l] + (dfafbOxOxBB[l+1] - dfafbOxOxBB[l])*xi;
        t2 = dfafbOxOxBB[l+1] + (dfafbOxOxBB[l+2] - dfafbOxOxBB[l+1])*(xi-1);
        fforceBB = - iq*jq*(t1 + (t2 - t1)*xi/2.0)*r ;

        if (fichierOxOx) { fichierOxOx<< setprecision (9)  <<r <<"  "<<evdwlCoul <<"  " <<fpairCoul <<"  "<<eng <<"  " <<fforce <<"  "<<engSurf <<"  " <<fforceSurf <<"  "<<engBB <<"  " <<fforceBB <<"  "<<ErepR<<"  "<<frepR<<"  "<<gam<<"  "<<dgam<<endl ;}

      }


      if (fichierOxOx) fichierOxOx.close() ;
    }



  map[type[j]]=1;  //met
  jtype=map[type[j]];
  jq=1;
  coord[i]=coordOxBulk;
  coord[j]=6;
  m = intype[itype][jtype];


  na = params[itype].ne ;
  nb = params[jtype].ne ;
  za = params[itype].dzeta ;
  zb = params[jtype].dzeta ;


  //   Ouverture du fichier

  for (iCoord=1;iCoord<4; iCoord++)
    {

      if (iCoord==1)
        {coord[i]=2.2;
          coord[j]=2.1;
          NameFile="energyandForceOxTiUnderCoord.txt";
        }
      if (iCoord==2)
        {coord[i]=coordOxBulk;
          coord[j]=coordOxBulk;
          NameFile="energyandForceOxTiCoordBulk.txt";
        }
      if (iCoord==3)
        {coord[i]=3.8;
          coord[j]=4;
          NameFile="energyandForceOxTiOverCoord.txt";
        }


      ofstream fichierOxTi(NameFile, ios::out | ios::trunc) ;

      drL=0.0001;
      iiiMax=int((cutmax-1.2)/drL);
      for (iii=1; iii< iiiMax ; iii++){
        r=1.2+drL*iii;
        rsq=r*r;
        evdwlCoul = 0.0 ; fpairCoul = 0.0;
        potqeq(i,j,iq,jq,rsq,fpairCoul,eflag,evdwlCoul);
        fpairCoul=fpairCoul*r;

        rep_OO (&intparams[m],rsq,fpair,eflag,evdwl);
        ErepR = evdwl;
        frepR= fpair*r;

        gam = dgam = dza = dzb = d2zaa = d2zab =
          d2zbb = d2zra = d2zrb = d2gamr2 = 0.0 ;


        //        gammas_(na,nb,za,zb,r,gam,dgam,dza,dzb,
        //                d2zaa,d2zab,d2zbb,d2zra,d2zrb,d2gamr2) ;
        gammas(na,nb,za,zb,r,gam,dgam,dza,dzb,
               d2zaa,d2zab,d2zbb,d2zra,d2zrb,d2gamr2) ;


        if (fichierOxTi) { fichierOxTi<< setprecision (9)  <<r <<"  "<<evdwlCoul <<"  " <<fpairCoul <<"  "<<ErepR<<"  "<<frepR<<"  "<<gam<<"  "<<dgam<<endl ;}

      }


      if (fichierOxTi) fichierOxTi.close() ;
    }

  exit(0);
}

/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   GAMMAS FUNCTION (GALE)
   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

void PairSMTBQ::gammas(double &na, double &nb, double &za, double &zb, double &r, double &gam,
                       double &dgam, double &dza, double &dzb, double &d2zaa, double &d2zab, double &d2zbb,
                       double &d2zra, double &d2zrb, double &d2gamr2)
{
  /*  ---------------------------------------------------------------
      Subroutine calculates the integral over two s orbtials
      for Slater functions

      On input :

      na = principle quantum number of A
      nb = principle quantum number of B
      za = orbital exponent of A
      zb = orbital exponent of B
      r  = distance between A and B

      On exit :

      gam     = integral
      dgam    = first derivative of gam with respect to r
      dza     = first derivative of gam with respect to za
      dzb     = first derivative of gam with respect to zb
      d2zaa   = second derivative of gam with respect to za/za
      d2zab   = second derivative of gam with respect to za/zb
      d2zbb   = second derivative of gam with respect to zb/zb
      d2zra   = second derivative of gam with respect to r/za
      d2zrb   = second derivative of gam with respect to r/zb
      d2gamr2 = second derivative of gam with respect to r

      Julian Gale, Imperial College, December 1997
      ---------------------------------------------------------------- */

  int i;
  double z2ra,z2rb,na2,nb2,halfr,drtrm,rtrm,ss,deriv,
    dzeta1,dzeta2,d2zeta11,d2zeta12,d2zeta22,d2zeta1r,
    d2zeta2r,deriv2,trm1,trm2,trm3,d2rtrm,ztrm,ctrm,
    rfct1,rgam1,rdza1,rgam2;



  gam=0.0;
  dgam=0.0;
  dza=0.0;
  dzb=0.0;
  d2zaa=0.0;
  d2zab=0.0;
  d2zbb=0.0;
  d2zra=0.0;
  d2zrb=0.0;
  d2gamr2=0.0;

  //  This routine only handles two centre integrals


  if (r < 1.0e-10) return;


  //  Create local variables

  z2ra=2.0*za*r;
  z2rb=2.0*zb*r;
  na2=2*na;
  nb2=2*nb;
  halfr=0.5*r;
  d2rtrm=powint(halfr,na2-2);
  drtrm=d2rtrm*halfr;
  rtrm=drtrm*halfr;

  //  First term


  css(ss,na2-1,0,z2ra,0.0,r,deriv,dzeta1,dzeta2,
      d2zeta11,d2zeta12,d2zeta22,d2zeta1r,d2zeta2r,
      deriv2)        ;

  gam=rtrm*ss;
  dgam=rtrm*deriv+static_cast<double>(na)*drtrm*ss;
  dza=rtrm*dzeta1;
  d2zaa=rtrm*d2zeta11;
  d2zra=rtrm*d2zeta1r+static_cast<double>(na)*drtrm*dzeta1;
  d2gamr2=d2gamr2+0.5*static_cast<double>(na*(na2-1))*d2rtrm*ss + 2.0*static_cast<double>(na)*drtrm*deriv+rtrm*deriv2;

  //  Sum over 2*nb

  rtrm=drtrm;
  drtrm=d2rtrm;
  ztrm=0.5/(zb*static_cast<double>(nb2));

  for (i = nb2; i >= 1; i--) {
    rtrm=rtrm*halfr;
    drtrm=drtrm*halfr;
    ztrm=ztrm*2.0*zb;
    ctrm=ztrm/factorial(static_cast<int>(nb2-i));

    css(ss,na2-1,nb2-i,z2ra,z2rb,r,deriv,dzeta1,dzeta2,
        d2zeta11,d2zeta12,d2zeta22,d2zeta1r,d2zeta2r,deriv2);

    trm1=static_cast<double>(i)*ctrm;
    trm2=trm1*rtrm;
    gam=gam-trm2*ss;
    trm3=trm1*static_cast<double>(na2+nb2-i)*drtrm;
    dgam=dgam-trm2*deriv-0.5*trm3*ss;
    d2gamr2=d2gamr2-trm2*deriv2-trm3*deriv-0.5*trm3*static_cast<double>(na2+nb2-i-1)*ss/r;
    dza=dza-trm2*dzeta1;
    dzb=dzb-(trm2/zb)*((static_cast<double>(nb2-i))*ss+zb*dzeta2);
    d2zaa=d2zaa-trm2*d2zeta11;
    d2zab=d2zab-(trm2/zb)*((static_cast<double>(nb2-i))*dzeta1+zb*d2zeta12);
    d2zbb=d2zbb-(trm2/zb)*(2.0*(static_cast<double>(nb2-i))*dzeta2+zb*d2zeta22 +
                           (static_cast<double>((nb2-i-1)*(nb2-i))*ss/zb));
    d2zra=d2zra-trm2*d2zeta1r-0.5*trm3*dzeta1;
    d2zrb=d2zrb-(trm2/zb)*((static_cast<double>(nb2-i))*deriv+zb*d2zeta2r) -
      0.5*(trm3/zb)*((static_cast<double>(nb2-i))*ss+zb*dzeta2);
  }

  //  Multiply by coefficients

  trm3=powint(2.0*za,na2+1)/factorial(static_cast<int>(na2));
  gam=gam*trm3;
  dgam=dgam*trm3;
  rfct1=((static_cast<double>(na2+1))/za);
  rgam1=rfct1*gam;
  dza=dza*trm3;
  rdza1=2.0*rfct1*dza;
  dza=dza+rgam1;
  dzb=dzb*trm3;
  rgam2=rgam1*static_cast<double>(na2)/za;
  d2zaa=d2zaa*trm3+rgam2+rdza1;
  d2zab=d2zab*trm3+rfct1*dzb;
  d2zbb=d2zbb*trm3;
  d2zra=d2zra*trm3+rfct1*dgam;
  d2zrb=d2zrb*trm3;
  d2gamr2=d2gamr2*trm3;
  return;

}
/* --------------------------------------------------------------------------------
   Css
   -------------------------------------------------------------------------------- */
void PairSMTBQ::css(double &s, double nn1, double nn2, double alpha, double beta, double r,
                    double &deriv, double &dzeta1, double &dzeta2, double &d2zeta11, double &d2zeta12,
                    double &d2zeta22, double &d2zeta1r, double &d2zeta2r, double &deriv2)
{
  //      implicit real (a-h,o-z)
  //      common /fctrl/ fct(30) // A RAJOUTER DANS Pair_SMTBQ.h
  /* ------------------------------------------------------------------
     Modified integral calculation routine for Slater orbitals
     including derivatives. This version is for S orbitals only.

     dzeta1 and dzeta2 are the first derivatives with respect to zetas
     and d2zeta11/d2zeta12/d2zeta22 are the second.
     d2zeta1r and d2zeta2r are the mixed zeta/r second derivatives
     deriv2 is the second derivative with respect to r

     Julian Gale, Imperial College, December 1997
     ------------------------------------------------------------------- */
  int i,i1,nni1;  //ulim
  double ulim,n1,n2,p,pt,x,k,dpdz1,dpdz2,dptdz1,dptdz2,dpdr,dptdr,d2pdz1r,
    d2pdz2r,d2ptdz1r,d2ptdz2r,zeta1,zeta2,sumzeta,difzeta,coff;

  double da1[30],da2[30],db1[30],db2[30];
  double d2a11[30],d2a12[30],d2a22[30],dar[30];
  double d2b11[30],d2b12[30],d2b22[30],dbr[30];
  double d2a1r[30],d2a2r[30],d2b1r[30],d2b2r[30];
  double d2ar2[30],d2br2[30];

  double *a,*b;
  memory->create(a,31,"pair:a");
  memory->create(b,31,"pair:a");


  //  Set up factorials - stored as factorial(n) in location(n+1)

  for (i = 1; i <= 30; i++) {
    fct[i]=factorial(i-1);
  }
  dzeta1=0.0;
  dzeta2=0.0;
  d2zeta11=0.0;
  d2zeta12=0.0;
  d2zeta22=0.0;
  d2zeta1r=0.0;
  d2zeta2r=0.0;
  deriv=0.0;
  deriv2=0.0;
  n1=nn1;
  n2=nn2;
  p =(alpha + beta)*0.5;
  pt=(alpha - beta)*0.5;
  x = 0.0;
  zeta1=alpha/r;
  zeta2=beta/r;
  sumzeta=zeta1+zeta2;
  difzeta=zeta1-zeta2;

  //  Partial derivative terms for zeta derivatives

  dpdz1=r;
  dpdz2=r;
  dptdz1=r;
  dptdz2=-r;
  dpdr=0.5*sumzeta;
  dptdr=0.5*difzeta;
  d2pdz1r=1.0;
  d2pdz2r=1.0;
  d2ptdz1r=1.0;
  d2ptdz2r=-1.0;

  //  Reverse quantum numbers if necessary -
  //  also change the sign of difzeta to match
  //  change in sign of pt

  if (n2 < n1) {
    k = n1;
    n1= n2;
    n2= k;
    pt=-pt;
    difzeta=-difzeta;
    dptdr=-dptdr;
    dptdz1=-dptdz1;
    dptdz2=-dptdz2;
    d2ptdz1r=-d2ptdz1r;
    d2ptdz2r=-d2ptdz2r;
  }

  //  Trap for enormously long distances which would cause
  //  caintgs or cbintgs to crash with an overflow

  if (p > 86.0 || pt > 86.0) {
    s=0.0;
    return;
  }
  //***************************
  //  Find a and b integrals  *
  //***************************
  caintgs(p,n1+n2+3,a);
  cbintgs(pt,n1+n2+3,b);

  //  Convert derivatives with respect to p and pt
  //  into derivatives with respect to zeta1 and
  //  zeta2

  ulim=n1+n2+1;
  for (i = 1; i <= int(ulim); i++) {
    da1[i]=-a[i+1]*dpdz1;
    da2[i]=-a[i+1]*dpdz2;
    db1[i]=-b[i+1]*dptdz1;
    db2[i]=-b[i+1]*dptdz2;
    d2a11[i]=a[i+2]*dpdz1*dpdz1;
    d2a12[i]=a[i+2]*dpdz1*dpdz2;
    d2a22[i]=a[i+2]*dpdz2*dpdz2;
    d2b11[i]=b[i+2]*dptdz1*dptdz1;
    d2b12[i]=b[i+2]*dptdz1*dptdz2;
    d2b22[i]=b[i+2]*dptdz2*dptdz2;
    dar[i]=-a[i+1]*dpdr;
    dbr[i]=-b[i+1]*dptdr;
    d2a1r[i]=a[i+2]*dpdz1*dpdr-a[i+1]*d2pdz1r;
    d2a2r[i]=a[i+2]*dpdz2*dpdr-a[i+1]*d2pdz2r;
    d2b1r[i]=b[i+2]*dptdz1*dptdr-b[i+1]*d2ptdz1r;
    d2b2r[i]=b[i+2]*dptdz2*dptdr-b[i+1]*d2ptdz2r;
    d2ar2[i]=a[i+2]*dpdr*dpdr;
    d2br2[i]=b[i+2]*dptdr*dptdr;
  }

  //  Begin section used for overlap integrals involving s functions

  for (i1 = 1; i1 <= int(ulim); i1++) {
    nni1=n1+n2-i1+2;
    coff=coeffs(n1,n2,i1-1);
    deriv=deriv+coff*(dar[i1]*b[nni1]+a[i1]*dbr[nni1]);
    x=x+coff*a[i1]*b[nni1];
    dzeta1=dzeta1+coff*(da1[i1]*b[nni1]+a[i1]*db1[nni1]);
    dzeta2=dzeta2+coff*(da2[i1]*b[nni1]+a[i1]*db2[nni1]);
    d2zeta11=d2zeta11+coff*(d2a11[i1]*b[nni1]+a[i1]*d2b11[nni1]+
                            2.0*da1[i1]*db1[nni1]);
    d2zeta12=d2zeta12+coff*(d2a12[i1]*b[nni1]+a[i1]*d2b12[nni1]+
                            da1[i1]*db2[nni1]+da2[i1]*db1[nni1]);
    d2zeta22=d2zeta22+coff*(d2a22[i1]*b[nni1]+a[i1]*d2b22[nni1]+
                            2.0*da2[i1]*db2[nni1]);
    d2zeta1r=d2zeta1r+coff*(d2a1r[i1]*b[nni1]+dar[i1]*db1[nni1]+
                            da1[i1]*dbr[nni1]+a[i1]*d2b1r[nni1]);
    d2zeta2r=d2zeta2r+coff*(d2a2r[i1]*b[nni1]+dar[i1]*db2[nni1]+
                            da2[i1]*dbr[nni1]+a[i1]*d2b2r[nni1]);
    deriv2=deriv2+coff*(d2ar2[i1]*b[nni1]+a[i1]*d2br2[nni1]+
                        2.0*dar[i1]*dbr[nni1]);
  }
  s=x*0.5;
  deriv=0.5*deriv;
  deriv2=0.5*deriv2;
  dzeta1=0.5*dzeta1;
  dzeta2=0.5*dzeta2;
  d2zeta11=0.5*d2zeta11;
  d2zeta12=0.5*d2zeta12;
  d2zeta22=0.5*d2zeta22;
  d2zeta1r=0.5*d2zeta1r;
  d2zeta2r=0.5*d2zeta2r;

  memory->destroy(a);
  memory->destroy(b);

  return;
}
/* -------------------------------------------------------------------------------
   coeffs
   ------------------------------------------------------------------------------- */
double PairSMTBQ::coeffs(int na, int nb, int k)
{
  //     implicit real (a-h,o-z)
  //     common /fctrl/ fct(30)

  int il,je,ia,i,j,ie,l;
  double coeffs;

  // Statement function
  //      binm(n,i)=fct(n+1)/(fct(n-i+1)*fct(i+1));

  coeffs=0.0;
  l=na+nb-k;
  ie=min(l,na)+1;
  je=min(l,nb);
  ia=l-je+1;
  for (il = ia; il <= ie; il++) {
    i=il-1;
    j=l-i;  // D'ou vient le i
    coeffs=coeffs + binm(na,i)*binm(nb,j)*powint(-1.,j);
  }
  return coeffs;
}

// ============================================

double PairSMTBQ::binm(int n, int i)
{
  return fct[n+1]/(fct[n-i+1]*fct[i+1]);
}

/* ---------------------------------------------------------------------------------
   Caintgs
   --------------------------------------------------------------------------------  */
void PairSMTBQ::caintgs (double x, int k, double *a)
{
  //      implicit real (a-h,o-z)
  //      dimension a(30)
  int i;
  double cste,rx;

  cste=exp(-x);
  rx=1.0/x;
  a[1]=cste*rx;
  for (i = 1; i <= k; i++) {
    a[i+1]=(a[i]*static_cast<double>(i)+cste)*rx;
  }
  return;
}
/* -----------------------------------------------------------------------------------
   Cbintgs
   ----------------------------------------------------------------------------------- */
void PairSMTBQ::cbintgs( double x, int k, double *b)
{
  //      implicit real (a-h,o-z)
  /* *******************************************************************
     ! Fills array of b-integrals. note that b(i) is b(i-1) in the
     ! usual notation
     ! for x.gt.3                          exponential formula is used
     ! for 2.lt.x.le.3 and k.le.10   exponential formula is used
     ! for 2.lt.x.le.3 and k.gt.10   15 term series is used
     ! for 1.lt.x .e.2 and k.le.7    exponential formula is used
     ! for 1.lt.x.le.2 and k.gt.7    12 term series is used
     ! for .5.lt.x.le.1 and k.le.5   exponential formula is used
     ! for .5.lt.x.le.1 and k.gt.5    7 term series is used
     ! for x.le..5                    6 term series is used
     !******************************************************************* */
  //      dimension b(30)
  //      common /fctrl/ fct(30)

  int i0,m,last,i;
  double absx,expx,expmx,ytrm,y,rx;

  i0=0;
  absx=fabs(x);

  if (absx > 3.0) goto g120;
  if (absx > 2.0) goto g20;
  if (absx > 1.0) goto g50;
  if (absx > 0.5) goto g80;
  if (absx > 1.0e-8) goto g110;
  goto g170;
 g110: last=6;
  goto g140;
 g80: if (k <= 5) goto g120;
  last=7;
  goto g140;
 g50: if (k <= 7) goto g120;
  last=12;
  goto g140;
 g20: if (k <= 10) goto g120;
  last=15;
  goto g140;

 g120: expx=exp(x);
  expmx=1./expx;
  rx=1.0/x;
  b[1]=(expx-expmx)*rx;
  for (i = 1; i <= k ; i++) {
    b[i+1]=(static_cast<double>(i)*b[i]+ powint(-1.0,i)*expx-expmx)*rx;
  }
  goto g190;
  //
  //  Series to calculate b(i)
  //
 g140: for (i = i0; i <= k ; i++) {
    y=0.;
    for (m = i0; m <= last; m++) {
      ytrm = powint(-x,m-1)*(1. - powint(-1.,m+i+1))/(fct[m+1]*static_cast<double>(m+i+1));
      y = y + ytrm*(-x);
    }
    b[i+1] = y;
  }
  goto g190;
  //
  //  x extremely small
  //
 g170: for (i = i0; i <= k; i++) {
    b[i+1] = (1.-powint(-1.,i+1))/static_cast<double>(i+1);
  }
 g190:
  return;
}

/* ============================== This is the END... ================================== */
