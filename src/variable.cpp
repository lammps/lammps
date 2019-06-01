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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <unistd.h>
#include "variable.h"
#include "universe.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "comm.h"
#include "region.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "fix_store.h"
#include "force.h"
#include "output.h"
#include "thermo.h"
#include "random_mars.h"
#include "math_const.h"
#include "atom_masks.h"
#include "lmppython.h"
#include "memory.h"
#include "info.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define VARDELTA 4
#define MAXLEVEL 4
#define MAXLINE 256
#define CHUNK 1024
#define VALUELENGTH 64               // also in python.cpp
#define MAXFUNCARG 6

#define MYROUND(a) (( a-floor(a) ) >= .5) ? ceil(a) : floor(a)

enum{INDEX,LOOP,WORLD,UNIVERSE,ULOOP,STRING,GETENV,
     SCALARFILE,ATOMFILE,FORMAT,EQUAL,ATOM,VECTOR,PYTHON,INTERNAL};
enum{ARG,OP};

// customize by adding a function
// if add before XOR:
// also set precedence level in constructor and precedence length in *.h

enum{DONE,ADD,SUBTRACT,MULTIPLY,DIVIDE,CARAT,MODULO,UNARY,
     NOT,EQ,NE,LT,LE,GT,GE,AND,OR,XOR,
     SQRT,EXP,LN,LOG,ABS,SIN,COS,TAN,ASIN,ACOS,ATAN,ATAN2,
     RANDOM,NORMAL,CEIL,FLOOR,ROUND,RAMP,STAGGER,LOGFREQ,LOGFREQ2,
     LOGFREQ3,STRIDE,STRIDE2,VDISPLACE,SWIGGLE,CWIGGLE,GMASK,RMASK,
     GRMASK,IS_ACTIVE,IS_DEFINED,IS_AVAILABLE,
     VALUE,ATOMARRAY,TYPEARRAY,INTARRAY,BIGINTARRAY,VECTORARRAY};

// customize by adding a special function

enum{SUM,XMIN,XMAX,AVE,TRAP,SLOPE};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PERATOM 8

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

Variable::Variable(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);

  nvar = maxvar = 0;
  names = NULL;
  style = NULL;
  num = NULL;
  which = NULL;
  pad = NULL;
  reader = NULL;
  data = NULL;
  dvalue = NULL;
  vecs = NULL;

  eval_in_progress = NULL;

  randomequal = NULL;
  randomatom = NULL;

  // customize by assigning a precedence level

  precedence[DONE] = 0;
  precedence[OR] = precedence[XOR] = 1;
  precedence[AND] = 2;
  precedence[EQ] = precedence[NE] = 3;
  precedence[LT] = precedence[LE] = precedence[GT] = precedence[GE] = 4;
  precedence[ADD] = precedence[SUBTRACT] = 5;
  precedence[MULTIPLY] = precedence[DIVIDE] = precedence[MODULO] = 6;
  precedence[CARAT] = 7;
  precedence[UNARY] = precedence[NOT] = 8;
}

/* ---------------------------------------------------------------------- */

Variable::~Variable()
{
  for (int i = 0; i < nvar; i++) {
    delete [] names[i];
    delete reader[i];
    if (style[i] == LOOP || style[i] == ULOOP) delete [] data[i][0];
    else for (int j = 0; j < num[i]; j++) delete [] data[i][j];
    delete [] data[i];
    if (style[i] == VECTOR) memory->destroy(vecs[i].values);
  }
  memory->sfree(names);
  memory->destroy(style);
  memory->destroy(num);
  memory->destroy(which);
  memory->destroy(pad);
  memory->sfree(reader);
  memory->sfree(data);
  memory->sfree(dvalue);
  memory->sfree(vecs);

  memory->destroy(eval_in_progress);

  delete randomequal;
  delete randomatom;

}

/* ----------------------------------------------------------------------
   called by variable command in input script
------------------------------------------------------------------------- */

void Variable::set(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal variable command");

  int replaceflag = 0;

  // DELETE
  // doesn't matter if variable no longer exists

  if (strcmp(arg[1],"delete") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal variable command");
    if (find(arg[0]) >= 0) remove(find(arg[0]));
    return;

  // INDEX
  // num = listed args, which = 1st value, data = copied args

  } else if (strcmp(arg[1],"index") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal variable command");
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) grow();
    style[nvar] = INDEX;
    num[nvar] = narg - 2;
    which[nvar] = 0;
    pad[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(num[nvar],&arg[2],data[nvar]);

  // LOOP
  // 1 arg + pad: num = N, which = 1st value, data = single string
  // 2 args + pad: num = N2, which = N1, data = single string

  } else if (strcmp(arg[1],"loop") == 0) {
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) grow();
    style[nvar] = LOOP;
    int nfirst = 0,nlast = 0;
    if (narg == 3 || (narg == 4 && strcmp(arg[3],"pad") == 0)) {
      nfirst = 1;
      nlast = force->inumeric(FLERR,arg[2]);
      if (nlast <= 0) error->all(FLERR,"Illegal variable command");
      if (narg == 4 && strcmp(arg[3],"pad") == 0) {
        char digits[12];
        sprintf(digits,"%d",nlast);
        pad[nvar] = strlen(digits);
      } else pad[nvar] = 0;
    } else if (narg == 4 || (narg == 5 && strcmp(arg[4],"pad") == 0)) {
      nfirst = force->inumeric(FLERR,arg[2]);
      nlast = force->inumeric(FLERR,arg[3]);
      if (nfirst > nlast || nlast < 0)
        error->all(FLERR,"Illegal variable command");
      if (narg == 5 && strcmp(arg[4],"pad") == 0) {
        char digits[12];
        sprintf(digits,"%d",nlast);
        pad[nvar] = strlen(digits);
      } else pad[nvar] = 0;
    } else error->all(FLERR,"Illegal variable command");
    num[nvar] = nlast;
    which[nvar] = nfirst-1;
    data[nvar] = new char*[1];
    data[nvar][0] = NULL;

  // WORLD
  // num = listed args, which = partition this proc is in, data = copied args
  // error check that num = # of worlds in universe

  } else if (strcmp(arg[1],"world") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal variable command");
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) grow();
    style[nvar] = WORLD;
    num[nvar] = narg - 2;
    if (num[nvar] != universe->nworlds)
      error->all(FLERR,"World variable count doesn't match # of partitions");
    which[nvar] = universe->iworld;
    pad[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(num[nvar],&arg[2],data[nvar]);

  // UNIVERSE and ULOOP
  // for UNIVERSE: num = listed args, data = copied args
  // for ULOOP: num = N, data = single string
  // which = partition this proc is in
  // universe proc 0 creates lock file
  // error check that all other universe/uloop variables are same length

  } else if (strcmp(arg[1],"universe") == 0 || strcmp(arg[1],"uloop") == 0) {
    if (strcmp(arg[1],"universe") == 0) {
      if (narg < 3) error->all(FLERR,"Illegal variable command");
      if (find(arg[0]) >= 0) return;
      if (nvar == maxvar) grow();
      style[nvar] = UNIVERSE;
      num[nvar] = narg - 2;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      copy(num[nvar],&arg[2],data[nvar]);
    } else if (strcmp(arg[1],"uloop") == 0) {
      if (narg < 3 || narg > 4 || (narg == 4 && strcmp(arg[3],"pad") != 0))
        error->all(FLERR,"Illegal variable command");
      if (find(arg[0]) >= 0) return;
      if (nvar == maxvar) grow();
      style[nvar] = ULOOP;
      num[nvar] = force->inumeric(FLERR,arg[2]);
      data[nvar] = new char*[1];
      data[nvar][0] = NULL;
      if (narg == 4) {
        char digits[12];
        sprintf(digits,"%d",num[nvar]);
        pad[nvar] = strlen(digits);
      } else pad[nvar] = 0;
    }

    if (num[nvar] < universe->nworlds)
      error->all(FLERR,"Universe/uloop variable count < # of partitions");
    which[nvar] = universe->iworld;

    if (universe->me == 0) {
      FILE *fp = fopen("tmp.lammps.variable","w");
      if (fp == NULL)
        error->one(FLERR,"Cannot open temporary file for world counter.");
      fprintf(fp,"%d\n",universe->nworlds);
      fclose(fp);
      fp = NULL;
    }

    for (int jvar = 0; jvar < nvar; jvar++)
      if (num[jvar] && (style[jvar] == UNIVERSE || style[jvar] == ULOOP) &&
          num[nvar] != num[jvar])
        error->all(FLERR,
                   "All universe/uloop variables must have same # of values");

  // STRING
  // replace pre-existing var if also style STRING (allows it to be reset)
  // num = 1, which = 1st value
  // data = 1 value, string to eval

  } else if (strcmp(arg[1],"string") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");

    int maxcopy = strlen(arg[2]) + 1;
    int maxwork = maxcopy;
    char *scopy = new char[maxcopy];
    char *work = new char[maxwork];
    strcpy(scopy,arg[2]);
    input->substitute(scopy,work,maxcopy,maxwork,1);
    delete [] work;

    int ivar = find(arg[0]);
    if (ivar >= 0) {
      if (style[ivar] != STRING)
        error->all(FLERR,"Cannot redefine variable as a different style");
      delete [] data[ivar][0];
      copy(1,&scopy,data[ivar]);
      replaceflag = 1;
    } else {
      if (nvar == maxvar) grow();
      style[nvar] = STRING;
      num[nvar] = 1;
      which[nvar] = 0;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      copy(1,&scopy,data[nvar]);
    }
    delete [] scopy;

  // GETENV
  // remove pre-existing var if also style GETENV (allows it to be reset)
  // num = 1, which = 1st value
  // data = 1 value, string to eval

  } else if (strcmp(arg[1],"getenv") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    if (find(arg[0]) >= 0) {
      if (style[find(arg[0])] != GETENV)
        error->all(FLERR,"Cannot redefine variable as a different style");
      remove(find(arg[0]));
    }
    if (nvar == maxvar) grow();
    style[nvar] = GETENV;
    num[nvar] = 1;
    which[nvar] = 0;
    pad[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(1,&arg[2],data[nvar]);
    data[nvar][1] = new char[VALUELENGTH];
    strcpy(data[nvar][1],"(undefined)");

  // SCALARFILE for strings or numbers
  // which = 1st value
  // data = 1 value, string to eval

  } else if (strcmp(arg[1],"file") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) grow();
    style[nvar] = SCALARFILE;
    num[nvar] = 1;
    which[nvar] = 0;
    pad[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    data[nvar][0] = new char[MAXLINE];
    reader[nvar] = new VarReader(lmp,arg[0],arg[2],SCALARFILE);
    int flag = reader[nvar]->read_scalar(data[nvar][0]);
    if (flag) error->all(FLERR,"File variable could not read value");

  // ATOMFILE for numbers
  // which = 1st value
  // data = NULL

  } else if (strcmp(arg[1],"atomfile") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) grow();
    style[nvar] = ATOMFILE;
    num[nvar] = 1;
    which[nvar] = 0;
    pad[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    data[nvar][0] = NULL;
    reader[nvar] = new VarReader(lmp,arg[0],arg[2],ATOMFILE);
    int flag = reader[nvar]->read_peratom();
    if (flag) error->all(FLERR,"Atomfile variable could not read values");

  // FORMAT
  // num = 3, which = 1st value
  // data = 3 values
  //   1st is name of variable to eval, 2nd is format string,
  //   3rd is filled on retrieval

  } else if (strcmp(arg[1],"format") == 0) {
    if (narg != 4) error->all(FLERR,"Illegal variable command");
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) grow();
    style[nvar] = FORMAT;
    num[nvar] = 3;
    which[nvar] = 0;
    pad[nvar] = 0;
    if (!utils::strmatch(arg[3],"%[0-9 ]*\\.[0-9]+[efgEFG]"))
      error->all(FLERR,"Incorrect conversion in format string");
    data[nvar] = new char*[num[nvar]];
    copy(2,&arg[2],data[nvar]);
    data[nvar][2] = new char[VALUELENGTH];
    strcpy(data[nvar][2],"(undefined)");

  // EQUAL
  // replace pre-existing var if also style EQUAL (allows it to be reset)
  // num = 2, which = 1st value
  // data = 2 values, 1st is string to eval, 2nd is filled on retrieval

  } else if (strcmp(arg[1],"equal") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    int ivar = find(arg[0]);
    if (ivar >= 0) {
      if (style[ivar] != EQUAL)
        error->all(FLERR,"Cannot redefine variable as a different style");
      delete [] data[ivar][0];
      copy(1,&arg[2],data[ivar]);
      replaceflag = 1;
    } else {
      if (nvar == maxvar) grow();
      style[nvar] = EQUAL;
      num[nvar] = 2;
      which[nvar] = 0;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      copy(1,&arg[2],data[nvar]);
      data[nvar][1] = new char[VALUELENGTH];
      strcpy(data[nvar][1],"(undefined)");
    }

  // ATOM
  // replace pre-existing var if also style ATOM (allows it to be reset)
  // num = 1, which = 1st value
  // data = 1 value, string to eval

  } else if (strcmp(arg[1],"atom") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    int ivar = find(arg[0]);
    if (ivar >= 0) {
      if (style[ivar] != ATOM)
        error->all(FLERR,"Cannot redefine variable as a different style");
      delete [] data[ivar][0];
      copy(1,&arg[2],data[ivar]);
      replaceflag = 1;
    } else {
      if (nvar == maxvar) grow();
      style[nvar] = ATOM;
      num[nvar] = 1;
      which[nvar] = 0;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      copy(1,&arg[2],data[nvar]);
    }

  // VECTOR
  // replace pre-existing var if also style VECTOR (allows it to be reset)
  // num = 1, which = 1st value
  // data = 1 value, string to eval

  } else if (strcmp(arg[1],"vector") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    int ivar = find(arg[0]);
    if (ivar >= 0) {
      if (style[ivar] != VECTOR)
        error->all(FLERR,"Cannot redefine variable as a different style");
      delete [] data[ivar][0];
      copy(1,&arg[2],data[ivar]);
      replaceflag = 1;
    } else {
      if (nvar == maxvar) grow();
      style[nvar] = VECTOR;
      num[nvar] = 1;
      which[nvar] = 0;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      copy(1,&arg[2],data[nvar]);
    }

  // PYTHON
  // replace pre-existing var if also style PYTHON (allows it to be reset)
  // num = 2, which = 1st value
  // data = 2 values, 1st is Python func to invoke, 2nd is filled by invoke

  } else if (strcmp(arg[1],"python") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    if (!python->is_enabled())
      error->all(FLERR,"LAMMPS is not built with Python embedded");
    int ivar = find(arg[0]);
    if (ivar >= 0) {
      if (style[ivar] != PYTHON)
        error->all(FLERR,"Cannot redefine variable as a different style");
      delete [] data[ivar][0];
      copy(1,&arg[2],data[ivar]);
      replaceflag = 1;
    } else {
      if (nvar == maxvar) grow();
      style[nvar] = PYTHON;
      num[nvar] = 2;
      which[nvar] = 1;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      copy(1,&arg[2],data[nvar]);
      data[nvar][1] = new char[VALUELENGTH];
      strcpy(data[nvar][1],"(undefined)");
    }

  // INTERNAL
  // replace pre-existing var if also style INTERNAL (allows it to be reset)
  // num = 1, for string representation of dvalue, set by retrieve()
  // dvalue = numeric initialization from 2nd arg, reset by internal_set()

  } else if (strcmp(arg[1],"internal") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    int ivar = find(arg[0]);
    if (ivar >= 0) {
      if (style[ivar] != INTERNAL)
        error->all(FLERR,"Cannot redefine variable as a different style");
      dvalue[nvar] = force->numeric(FLERR,arg[2]);
      replaceflag = 1;
    } else {
      if (nvar == maxvar) grow();
      style[nvar] = INTERNAL;
      num[nvar] = 1;
      which[nvar] = 0;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      data[nvar][0] = new char[VALUELENGTH];
      dvalue[nvar] = force->numeric(FLERR,arg[2]);
    }

  } else error->all(FLERR,"Illegal variable command");

  // set name of variable, if not replacing one flagged with replaceflag
  // name must be all alphanumeric chars or underscores

  if (replaceflag) return;

  int n = strlen(arg[0]) + 1;
  names[nvar] = new char[n];
  strcpy(names[nvar],arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(names[nvar][i]) && names[nvar][i] != '_') {
      char errmsg[128];
      snprintf(errmsg,128,"Variable name '%s' must have only alphanumeric "
              "characters or underscore",names[nvar]);
      error->all(FLERR,errmsg);
    }
  nvar++;
}

/* ----------------------------------------------------------------------
   INDEX variable created by command-line argument
   make it INDEX rather than STRING so cannot be re-defined in input script
------------------------------------------------------------------------- */

void Variable::set(char *name, int narg, char **arg)
{
  char **newarg = new char*[2+narg];
  newarg[0] = name;
  newarg[1] = (char *) "index";
  for (int i = 0; i < narg; i++) newarg[2+i] = arg[i];
  set(2+narg,newarg);
  delete [] newarg;
}

/* ----------------------------------------------------------------------
   set existing STRING variable to str
   return 0 if successful
   return -1 if variable doesn't exist or isn't a STRING variable
   called via library interface, so external programs can set variables
------------------------------------------------------------------------- */

int Variable::set_string(char *name, char *str)
{
  int ivar = find(name);
  if (ivar < 0) return -1;
  if (style[ivar] != STRING) return -1;
  delete [] data[ivar][0];
  copy(1,&str,data[ivar]);
  return 0;
}

/* ----------------------------------------------------------------------
   increment variable(s)
   return 0 if OK if successfully incremented
   return 1 if any variable is exhausted, free the variable to allow re-use
------------------------------------------------------------------------- */

int Variable::next(int narg, char **arg)
{
  int ivar;

  if (narg == 0) error->all(FLERR,"Illegal next command");

  // check that variables exist and are all the same style
  // exception: UNIVERSE and ULOOP variables can be mixed in same next command

  for (int iarg = 0; iarg < narg; iarg++) {
    ivar = find(arg[iarg]);
    if (ivar < 0) {
      char errmsg[128];
      snprintf(errmsg,128,"Invalid variable '%s' in next command",arg[iarg]);
      error->all(FLERR,errmsg);
    }
    if (style[ivar] == ULOOP && style[find(arg[0])] == UNIVERSE) continue;
    else if (style[ivar] == UNIVERSE && style[find(arg[0])] == ULOOP) continue;
    else if (style[ivar] != style[find(arg[0])])
      error->all(FLERR,"All variables in next command must have same style");
  }

  // invalid styles: STRING, EQUAL, WORLD, ATOM, VECTOR, GETENV,
  //                 FORMAT, PYTHON, INTERNAL

  int istyle = style[find(arg[0])];
  if (istyle == STRING || istyle == EQUAL || istyle == WORLD ||
      istyle == GETENV || istyle == ATOM || istyle == VECTOR ||
      istyle == FORMAT || istyle == PYTHON || istyle == INTERNAL)
    error->all(FLERR,"Invalid variable style with next command");

  // if istyle = UNIVERSE or ULOOP, insure all such variables are incremented

  if (istyle == UNIVERSE || istyle == ULOOP)
    for (int i = 0; i < nvar; i++) {
      if (style[i] != UNIVERSE && style[i] != ULOOP) continue;
      int iarg = 0;
      for (iarg = 0; iarg < narg; iarg++)
        if (strcmp(arg[iarg],names[i]) == 0) break;
      if (iarg == narg)
        error->universe_one(FLERR,"Next command must list all "
                            "universe and uloop variables");
    }

  // increment all variables in list
  // if any variable is exhausted, set flag = 1 and remove var to allow re-use

  int flag = 0;

  if (istyle == INDEX || istyle == LOOP) {
    for (int iarg = 0; iarg < narg; iarg++) {
      ivar = find(arg[iarg]);
      which[ivar]++;
      if (which[ivar] >= num[ivar]) {
        flag = 1;
        remove(ivar);
      }
    }

  } else if (istyle == SCALARFILE) {

    for (int iarg = 0; iarg < narg; iarg++) {
      ivar = find(arg[iarg]);
      int done = reader[ivar]->read_scalar(data[ivar][0]);
      if (done) {
        flag = 1;
        remove(ivar);
      }
    }

  } else if (istyle == ATOMFILE) {

    for (int iarg = 0; iarg < narg; iarg++) {
      ivar = find(arg[iarg]);
      int done = reader[ivar]->read_peratom();
      if (done) {
        flag = 1;
        remove(ivar);
      }
    }

  } else if (istyle == UNIVERSE || istyle == ULOOP) {

    // wait until lock file can be created and owned by proc 0 of this world
    // rename() is not atomic in practice, but no known simple fix
    //   means multiple procs can read/write file at the same time (bad!)
    // random delays help
    // delay for random fraction of 1 second before first rename() call
    // delay for random fraction of 1 second before subsequent tries
    // when successful, read next available index and Bcast it within my world

    int nextindex;
    if (me == 0) {
      int seed = 12345 + universe->me + which[find(arg[0])];
      RanMars *random = new RanMars(lmp,seed);
      int delay = (int) (1000000*random->uniform());
      usleep(delay);
      while (1) {
        if (!rename("tmp.lammps.variable","tmp.lammps.variable.lock")) break;
        delay = (int) (1000000*random->uniform());
        usleep(delay);
      }
      delete random;

      FILE *fp = fopen("tmp.lammps.variable.lock","r");
      fscanf(fp,"%d",&nextindex);
      //printf("READ %d %d\n",universe->me,nextindex);
      fclose(fp);
      fp = fopen("tmp.lammps.variable.lock","w");
      fprintf(fp,"%d\n",nextindex+1);
      //printf("WRITE %d %d\n",universe->me,nextindex+1);
      fclose(fp);
      fp = NULL;
      rename("tmp.lammps.variable.lock","tmp.lammps.variable");
      if (universe->uscreen)
        fprintf(universe->uscreen,
                "Increment via next: value %d on partition %d\n",
                nextindex+1,universe->iworld);
      if (universe->ulogfile)
        fprintf(universe->ulogfile,
                "Increment via next: value %d on partition %d\n",
                nextindex+1,universe->iworld);
    }
    MPI_Bcast(&nextindex,1,MPI_INT,0,world);

    // set all variables in list to nextindex
    // must increment all UNIVERSE and ULOOP variables here
    // error check above tested for this

    for (int iarg = 0; iarg < narg; iarg++) {
      ivar = find(arg[iarg]);
      which[ivar] = nextindex;
      if (which[ivar] >= num[ivar]) {
        flag = 1;
        remove(ivar);
      }
    }
  }

  return flag;
}

/* ----------------------------------------------------------------------
   search for name in list of variables names
   return index or -1 if not found
------------------------------------------------------------------------- */

int Variable::find(char *name)
{
  if(name==NULL) return -1;
  for (int i = 0; i < nvar; i++)
    if (strcmp(name,names[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values in all VarReaders via fix STORE
   called when atom is created
------------------------------------------------------------------------- */

void Variable::set_arrays(int /*i*/)
{
  for (int i = 0; i < nvar; i++)
    if (reader[i] && style[i] == ATOMFILE)
      reader[i]->fixstore->vstore[i] = 0.0;
}

/* ----------------------------------------------------------------------
   called by python command in input script
   simply pass input script line args to Python class
------------------------------------------------------------------------- */

void Variable::python_command(int narg, char **arg)
{
  if (!python->is_enabled())
    error->all(FLERR,"LAMMPS is not built with Python embedded");
  python->command(narg,arg);
}

/* ----------------------------------------------------------------------
   return 1 if variable is EQUAL or INTERNAL or PYTHON numeric style, 0 if not
   this is checked before call to compute_equal() to return a double
------------------------------------------------------------------------- */

int Variable::equalstyle(int ivar)
{
  if (style[ivar] == EQUAL || style[ivar] == INTERNAL) return 1;
  if (style[ivar] == PYTHON) {
    int ifunc = python->variable_match(data[ivar][0],names[ivar],1);
    if (ifunc < 0) return 0;
    else return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   return 1 if variable is ATOM or ATOMFILE style, 0 if not
   this is checked before call to compute_atom() to return a vector of doubles
------------------------------------------------------------------------- */

int Variable::atomstyle(int ivar)
{
  if (style[ivar] == ATOM || style[ivar] == ATOMFILE) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return 1 if variable is VECTOR style, 0 if not
   this is checked before call to compute_vector() to return a vector of doubles
------------------------------------------------------------------------- */

int Variable::vectorstyle(int ivar)
{
  if (style[ivar] == VECTOR) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check if variable with name is PYTHON and matches funcname
   called by Python class before it invokes a Python function
   return data storage so Python function can return a value for this variable
   return NULL if not a match
------------------------------------------------------------------------- */

char *Variable::pythonstyle(char *name, char *funcname)
{
  int ivar = find(name);
  if (ivar < 0) return NULL;
  if (style[ivar] != PYTHON) return NULL;
  if (strcmp(data[ivar][0],funcname) != 0) return NULL;
  return data[ivar][1];
}

/* ----------------------------------------------------------------------
   return 1 if variable is INTERNAL style, 0 if not
   this is checked before call to set_internal() to assure it can be set
------------------------------------------------------------------------- */

int Variable::internalstyle(int ivar)
{
  if (style[ivar] == INTERNAL) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return ptr to the data text associated with a variable
   if INDEX or WORLD or UNIVERSE or STRING or SCALARFILE,
     return ptr to stored string
   if LOOP or ULOOP, write int to data[0] and return ptr to string
   if EQUAL, evaluate variable and put result in str
   if FORMAT, evaluate its variable and put formatted result in str
   if GETENV, query environment and put result in str
   if PYTHON, evaluate Python function, it will put result in str
   if INTERNAL, convert dvalue and put result in str
   if ATOM or ATOMFILE or VECTOR, return NULL
   return NULL if no variable with name, or which value is bad,
     caller must respond
------------------------------------------------------------------------- */

char *Variable::retrieve(char *name)
{
  int ivar = find(name);
  if (ivar < 0) return NULL;
  if (which[ivar] >= num[ivar]) return NULL;

  if (eval_in_progress[ivar])
    print_var_error(FLERR,"Variable has a circular dependency",ivar);

  eval_in_progress[ivar] = 1;

  char *str = NULL;
  if (style[ivar] == INDEX || style[ivar] == WORLD ||
      style[ivar] == UNIVERSE || style[ivar] == STRING ||
      style[ivar] == SCALARFILE) {
    str = data[ivar][which[ivar]];
  } else if (style[ivar] == LOOP || style[ivar] == ULOOP) {
    char result[16];
    if (pad[ivar] == 0) sprintf(result,"%d",which[ivar]+1);
    else {
      char padstr[16];
      sprintf(padstr,"%%0%dd",pad[ivar]);
      sprintf(result,padstr,which[ivar]+1);
    }
    int n = strlen(result) + 1;
    delete [] data[ivar][0];
    data[ivar][0] = new char[n];
    strcpy(data[ivar][0],result);
    str = data[ivar][0];
  } else if (style[ivar] == EQUAL) {
    double answer = evaluate(data[ivar][0],NULL,ivar);
    sprintf(data[ivar][1],"%.15g",answer);
    str = data[ivar][1];
  } else if (style[ivar] == FORMAT) {
    int jvar = find(data[ivar][0]);
    if (jvar == -1) return NULL;
    if (!equalstyle(jvar)) return NULL;
    double answer = compute_equal(jvar);
    sprintf(data[ivar][2],data[ivar][1],answer);
    str = data[ivar][2];
  } else if (style[ivar] == GETENV) {
    const char *result = getenv(data[ivar][0]);
    if (result == NULL) result = (const char *) "";
    int n = strlen(result) + 1;
    if (n > VALUELENGTH) {
      delete [] data[ivar][1];
      data[ivar][1] = new char[n];
    }
    strcpy(data[ivar][1],result);
    str = data[ivar][1];
  } else if (style[ivar] == PYTHON) {
    int ifunc = python->variable_match(data[ivar][0],name,0);
    if (ifunc < 0) {
      char errmsg[128];
      snprintf(errmsg,128,"Python variable '%s' does not match Python function",name);
      error->all(FLERR,errmsg);
    }
    python->invoke_function(ifunc,data[ivar][1]);
    str = data[ivar][1];
    // if Python func returns a string longer than VALUELENGTH
    // then the Python class stores the result, query it via long_string()
    char *strlong = python->long_string(ifunc);
    if (strlong) str = strlong;
  } else if (style[ivar] == INTERNAL) {
    sprintf(data[ivar][0],"%.15g",dvalue[ivar]);
    str = data[ivar][0];
  } else if (style[ivar] == ATOM || style[ivar] == ATOMFILE ||
             style[ivar] == VECTOR) return NULL;

  eval_in_progress[ivar] = 0;

  return str;
}

/* ----------------------------------------------------------------------
   return result of equal-style variable evaluation
   can be EQUAL or INTERNAL style or PYTHON numeric style
   for PYTHON, don't need to check python->variable_match() error return,
     since caller will have already checked via equalstyle()
------------------------------------------------------------------------- */

double Variable::compute_equal(int ivar)
{
  if (eval_in_progress[ivar])
    print_var_error(FLERR,"Variable has a circular dependency",ivar);

  eval_in_progress[ivar] = 1;

  double value = 0.0;
  if (style[ivar] == EQUAL) value = evaluate(data[ivar][0],NULL,ivar);
  else if (style[ivar] == INTERNAL) value = dvalue[ivar];
  else if (style[ivar] == PYTHON) {
    int ifunc = python->find(data[ivar][0]);
    if (ifunc < 0)
      print_var_error(FLERR,"Python variable has no function",ivar);
    python->invoke_function(ifunc,data[ivar][1]);
    value = atof(data[ivar][1]);
  }

  eval_in_progress[ivar] = 0;
  return value;
}

/* ----------------------------------------------------------------------
   return result of immediate equal-style variable evaluation
   called from Input::substitute()
   don't need to flag eval_in_progress since is an immediate variable
------------------------------------------------------------------------- */

double Variable::compute_equal(char *str)
{
  return evaluate(str,NULL,-1);
}

/* ----------------------------------------------------------------------
   compute result of atom-style and atomfile-style variable evaluation
   only computed for atoms in igroup, else result is 0.0
   answers are placed every stride locations into result
   if sumflag, add variable values to existing result
------------------------------------------------------------------------- */

void Variable::compute_atom(int ivar, int igroup,
                            double *result, int stride, int sumflag)
{
  Tree *tree = NULL;
  double *vstore;

  if (eval_in_progress[ivar])
    print_var_error(FLERR,"Variable has a circular dependency",ivar);

  eval_in_progress[ivar] = 1;

  if (style[ivar] == ATOM) {
    treetype = ATOM;
    evaluate(data[ivar][0],&tree,ivar);
    collapse_tree(tree);
  } else vstore = reader[ivar]->fixstore->vstore;

  if (result == NULL) {
    if (style[ivar] == ATOM) free_tree(tree);
    eval_in_progress[ivar] = 0;
    return;
  }

  int groupbit = group->bitmask[igroup];
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (style[ivar] == ATOM) {
    if (sumflag == 0) {
      int m = 0;
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) result[m] = eval_tree(tree,i);
        else result[m] = 0.0;
        m += stride;
      }

    } else {
      int m = 0;
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) result[m] += eval_tree(tree,i);
        m += stride;
      }
    }

  } else {
    if (sumflag == 0) {
      int m = 0;
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) result[m] = vstore[i];
        else result[m] = 0.0;
        m += stride;
      }

    } else {
      int m = 0;
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) result[m] += vstore[i];
        m += stride;
      }
    }
  }

  if (style[ivar] == ATOM) free_tree(tree);
  eval_in_progress[ivar] = 0;
}

/* ----------------------------------------------------------------------
   compute result of vector-style variable evaluation
   return length of vector and result pointer to vector values
     if length == 0 or -1 (mismatch), generate an error
   if variable already computed on this timestep, just return
   else evaluate the formula and its length, store results in VecVar entry
------------------------------------------------------------------------- */

int Variable::compute_vector(int ivar, double **result)
{
  Tree *tree = NULL;
  if (vecs[ivar].currentstep == update->ntimestep) {
    *result = vecs[ivar].values;
    return vecs[ivar].n;
  }

  if (eval_in_progress[ivar])
    print_var_error(FLERR,"Variable has a circular dependency",ivar);

  eval_in_progress[ivar] = 1;

  treetype = VECTOR;
  evaluate(data[ivar][0],&tree,ivar);
  collapse_tree(tree);
  int nlen = size_tree_vector(tree);
  if (nlen == 0)
    print_var_error(FLERR,"Vector-style variable has zero length",ivar);

  if (nlen < 0)
    print_var_error(FLERR,"Inconsistent lengths in vector-style variable",ivar);

  // (re)allocate space for results if necessary

  if (nlen > vecs[ivar].nmax) {
    memory->destroy(vecs[ivar].values);
    vecs[ivar].nmax = nlen;
    memory->create(vecs[ivar].values,vecs[ivar].nmax,"variable:values");
  }

  vecs[ivar].n = nlen;
  vecs[ivar].currentstep = update->ntimestep;
  double *vec = vecs[ivar].values;
  for (int i = 0; i < nlen; i++)
    vec[i] = eval_tree(tree,i);

  free_tree(tree);
  eval_in_progress[ivar] = 0;

  *result = vec;
  return nlen;
}

/* ----------------------------------------------------------------------
   set value stored by INTERNAL style ivar
------------------------------------------------------------------------- */

void Variable::internal_set(int ivar, double value)
{
  dvalue[ivar] = value;
}

/* ----------------------------------------------------------------------
   remove Nth variable from list and compact list
   delete reader explicitly if it exists
------------------------------------------------------------------------- */

void Variable::remove(int n)
{
  delete [] names[n];
  if (style[n] == LOOP || style[n] == ULOOP) delete [] data[n][0];
  else for (int i = 0; i < num[n]; i++) delete [] data[n][i];
  delete [] data[n];
  delete reader[n];

  for (int i = n+1; i < nvar; i++) {
    names[i-1] = names[i];
    style[i-1] = style[i];
    num[i-1] = num[i];
    which[i-1] = which[i];
    pad[i-1] = pad[i];
    reader[i-1] = reader[i];
    data[i-1] = data[i];
    dvalue[i-1] = dvalue[i];
  }
  nvar--;
}

/* ----------------------------------------------------------------------
  make space in arrays for new variable
------------------------------------------------------------------------- */

void Variable::grow()
{
  int old = maxvar;
  maxvar += VARDELTA;
  names = (char **) memory->srealloc(names,maxvar*sizeof(char *),"var:names");
  memory->grow(style,maxvar,"var:style");
  memory->grow(num,maxvar,"var:num");
  memory->grow(which,maxvar,"var:which");
  memory->grow(pad,maxvar,"var:pad");

  reader = (VarReader **)
    memory->srealloc(reader,maxvar*sizeof(VarReader *),"var:reader");
  for (int i = old; i < maxvar; i++) reader[i] = NULL;

  data = (char ***) memory->srealloc(data,maxvar*sizeof(char **),"var:data");
  memory->grow(dvalue,maxvar,"var:dvalue");

  vecs = (VecVar *) memory->srealloc(vecs,maxvar*sizeof(VecVar),"var:vecvar");
  for (int i = old; i < maxvar; i++) {
    vecs[i].nmax = 0;
    vecs[i].currentstep = -1;
    vecs[i].values = NULL;
  }

  memory->grow(eval_in_progress,maxvar,"var:eval_in_progress");
  for (int i = 0; i < maxvar; i++) eval_in_progress[i] = 0;
}

/* ----------------------------------------------------------------------
   copy narg strings from **from to **to, and allocate space for them
------------------------------------------------------------------------- */

void Variable::copy(int narg, char **from, char **to)
{
  int n;
  for (int i = 0; i < narg; i++) {
    n = strlen(from[i]) + 1;
    to[i] = new char[n];
    strcpy(to[i],from[i]);
  }
}

/* ----------------------------------------------------------------------
   recursive evaluation of a string str
   str is an equal-style or atom-style or vector-style formula
     containing one or more items:
     number = 0.0, -5.45, 2.8e-4, ...
     constant = PI, version, yes, no, on, off
     thermo keyword = ke, vol, atoms, ...
     math operation = (),-x,x+y,x-y,x*y,x/y,x^y,
                      x==y,x!=y,x<y,x<=y,x>y,x>=y,x&&y,x||y,
                      sqrt(x),exp(x),ln(x),log(x),abs(x),
                      sin(x),cos(x),tan(x),asin(x),atan2(y,x),...
     group function = count(group), mass(group), xcm(group,x), ...
     special function = sum(x),min(x), ...
     atom value = x[i], y[i], vx[i], ...
     atom vector = x, y, vx, ...
     compute = c_ID, c_ID[i], c_ID[i][j]
     fix = f_ID, f_ID[i], f_ID[i][j]
     variable = v_name, v_name[i]
   equal-style variables passes in tree = NULL:
     evaluate the formula, return result as a double
   atom-style and vector-style variables pass in tree = non-NULL:
     parse the formula but do not evaluate it
     create a parse tree and return it
------------------------------------------------------------------------- */

double Variable::evaluate(char *str, Tree **tree, int ivar)
{
  int op,opprevious;
  double value1,value2;
  char onechar;
  char *ptr;

  double argstack[MAXLEVEL];
  Tree *treestack[MAXLEVEL];
  int opstack[MAXLEVEL];
  int nargstack = 0;
  int ntreestack = 0;
  int nopstack = 0;

  int i = 0;
  int expect = ARG;

  while (1) {
    onechar = str[i];

    // whitespace: just skip

    if (isspace(onechar)) i++;

    // ----------------
    // parentheses: recursively evaluate contents of parens
    // ----------------

    else if (onechar == '(') {
      if (expect == OP)
        print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
      expect = OP;

      char *contents;
      i = find_matching_paren(str,i,contents,ivar);
      i++;

      // evaluate contents and push on stack

      if (tree) {
        Tree *newtree = NULL;
        evaluate(contents,&newtree,ivar);
        treestack[ntreestack++] = newtree;
      } else argstack[nargstack++] = evaluate(contents,NULL,ivar);

      delete [] contents;

    // ----------------
    // number: push value onto stack
    // ----------------

    } else if (isdigit(onechar) || onechar == '.') {
      if (expect == OP)
        print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
      expect = OP;

      // istop = end of number, including scientific notation

      int istart = i;
      while (isdigit(str[i]) || str[i] == '.') i++;
      if (str[i] == 'e' || str[i] == 'E') {
        i++;
        if (str[i] == '+' || str[i] == '-') i++;
        while (isdigit(str[i])) i++;
      }
      int istop = i - 1;

      int n = istop - istart + 1;
      char *number = new char[n+1];
      strncpy(number,&str[istart],n);
      number[n] = '\0';

      if (tree) {
        Tree *newtree = new Tree();
        newtree->type = VALUE;
        newtree->value = atof(number);
        newtree->first = newtree->second = NULL;
        newtree->extra = 0;
        treestack[ntreestack++] = newtree;
      } else argstack[nargstack++] = atof(number);

      delete [] number;

    // ----------------
    // letter: c_ID, c_ID[], c_ID[][], f_ID, f_ID[], f_ID[][],
    //         v_name, v_name[], exp(), xcm(,), x, x[], PI, vol
    // ----------------

    } else if (isalpha(onechar)) {
      if (expect == OP)
        print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
      expect = OP;

      // istop = end of word
      // word = all alphanumeric or underscore

      int istart = i;
      while (isalnum(str[i]) || str[i] == '_') i++;
      int istop = i-1;

      int n = istop - istart + 1;
      char *word = new char[n+1];
      strncpy(word,&str[istart],n);
      word[n] = '\0';

      // ----------------
      // compute
      // ----------------

      if (strncmp(word,"c_",2) == 0 || strncmp(word,"C_",2) == 0) {
        if (domain->box_exist == 0)
          print_var_error(FLERR,"Variable evaluation before "
                          "simulation box is defined",ivar);

        // uppercase used to force access of
        // global vector vs global scalar, and global array vs global vector

        int lowercase = 1;
        if (word[0] == 'C') lowercase = 0;

        int icompute = modify->find_compute(word+2);
        if (icompute < 0)
          print_var_error(FLERR,"Invalid compute ID in variable formula",ivar);
        Compute *compute = modify->compute[icompute];

        // parse zero or one or two trailing brackets
        // point i beyond last bracket
        // nbracket = # of bracket pairs
        // index1,index2 = int inside each bracket pair, possibly an atom ID

        int nbracket;
        tagint index1,index2;
        if (str[i] != '[') nbracket = 0;
        else {
          nbracket = 1;
          ptr = &str[i];
          index1 = int_between_brackets(ptr,1);
          i = ptr-str+1;
          if (str[i] == '[') {
            nbracket = 2;
            ptr = &str[i];
            index2 = int_between_brackets(ptr,1);
            i = ptr-str+1;
          }
        }

        // c_ID = scalar from global scalar, must be lowercase

        if (nbracket == 0 && compute->scalar_flag && lowercase) {

          if (update->whichflag == 0) {
            if (compute->invoked_scalar != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between "
                              "runs is not current",ivar);
          } else if (!(compute->invoked_flag & INVOKED_SCALAR)) {
            compute->compute_scalar();
            compute->invoked_flag |= INVOKED_SCALAR;
          }

          value1 = compute->scalar;
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->first = newtree->second = NULL;
            newtree->nextra = 0;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // c_ID[i] = scalar from global vector, must be lowercase

        } else if (nbracket == 1 && compute->vector_flag && lowercase) {

          if (index1 > compute->size_vector &&
              compute->size_vector_variable == 0)
            print_var_error(FLERR,"Variable formula compute vector "
                            "is accessed out-of-range",ivar,0);
          if (update->whichflag == 0) {
            if (compute->invoked_vector != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between runs "
                              "is not current",ivar);
          } else if (!(compute->invoked_flag & INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= INVOKED_VECTOR;
          }

          if (compute->size_vector_variable &&
              index1 > compute->size_vector) value1 = 0.0;
          else value1 = compute->vector[index1-1];
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->first = newtree->second = NULL;
            newtree->nextra = 0;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // c_ID[i][j] = scalar from global array, must be lowercase

        } else if (nbracket == 2 && compute->array_flag && lowercase) {

          if (index1 > compute->size_array_rows &&
              compute->size_array_rows_variable == 0)
            print_var_error(FLERR,"Variable formula compute array "
                            "is accessed out-of-range",ivar,0);
          if (index2 > compute->size_array_cols)
            print_var_error(FLERR,"Variable formula compute array "
                            "is accessed out-of-range",ivar,0);
          if (update->whichflag == 0) {
            if (compute->invoked_array != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between runs "
                              "is not current",ivar);
          } else if (!(compute->invoked_flag & INVOKED_ARRAY)) {
            compute->compute_array();
            compute->invoked_flag |= INVOKED_ARRAY;
          }

          if (compute->size_array_rows_variable &&
              index1 > compute->size_array_rows) value1 = 0.0;
          else value1 = compute->array[index1-1][index2-1];
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->first = newtree->second = NULL;
            newtree->nextra = 0;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // c_ID = vector from global vector, lowercase or uppercase

        } else if (nbracket == 0 && compute->vector_flag) {

          if (tree == NULL)
            print_var_error(FLERR,"Compute global vector in "
                            "equal-style variable formula",ivar);
          if (treetype == ATOM)
            print_var_error(FLERR,"Compute global vector in "
                            "atom-style variable formula",ivar);
          if (compute->size_vector == 0)
            print_var_error(FLERR,"Variable formula compute "
                            "vector is zero length",ivar);
          if (update->whichflag == 0) {
            if (compute->invoked_vector != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between "
                              "runs is not current",ivar);
          } else if (!(compute->invoked_flag & INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= INVOKED_VECTOR;
          }

          Tree *newtree = new Tree();
          newtree->type = VECTORARRAY;
          newtree->array = compute->vector;
          newtree->nvector = compute->size_vector;
          newtree->nstride = 1;
          newtree->selfalloc = 0;
          newtree->first = newtree->second = NULL;
          newtree->nextra = 0;
          treestack[ntreestack++] = newtree;

        // c_ID[i] = vector from global array, lowercase or uppercase

        } else if (nbracket == 1 && compute->array_flag) {

          if (tree == NULL)
            print_var_error(FLERR,"Compute global vector in "
                            "equal-style variable formula",ivar);
          if (treetype == ATOM)
            print_var_error(FLERR,"Compute global vector in "
                            "atom-style variable formula",ivar);
          if (compute->size_array_rows == 0)
            print_var_error(FLERR,"Variable formula compute array "
                            "is zero length",ivar);
          if (update->whichflag == 0) {
            if (compute->invoked_array != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between "
                              "runs is not current",ivar);
          } else if (!(compute->invoked_flag & INVOKED_ARRAY)) {
            compute->compute_array();
            compute->invoked_flag |= INVOKED_ARRAY;
          }

          Tree *newtree = new Tree();
          newtree->type = VECTORARRAY;
          newtree->array = &compute->array[0][index1-1];
          newtree->nvector = compute->size_array_rows;
          newtree->nstride = compute->size_array_cols;
          newtree->selfalloc = 0;
          newtree->first = newtree->second = NULL;
          newtree->nextra = 0;
          treestack[ntreestack++] = newtree;

        // c_ID[i] = scalar from per-atom vector

        } else if (nbracket == 1 && compute->peratom_flag &&
                   compute->size_peratom_cols == 0) {

          if (update->whichflag == 0) {
            if (compute->invoked_peratom != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable "
                              "between runs is not current",ivar);
          } else if (!(compute->invoked_flag & INVOKED_PERATOM)) {
            compute->compute_peratom();
            compute->invoked_flag |= INVOKED_PERATOM;
          }

          peratom2global(1,NULL,compute->vector_atom,1,index1,
                         tree,treestack,ntreestack,argstack,nargstack);

        // c_ID[i][j] = scalar from per-atom array

        } else if (nbracket == 2 && compute->peratom_flag &&
                   compute->size_peratom_cols > 0) {

          if (index2 > compute->size_peratom_cols)
            print_var_error(FLERR,"Variable formula compute "
                            "array is accessed out-of-range",ivar,0);
          if (update->whichflag == 0) {
            if (compute->invoked_peratom != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable "
                              "between runs is not current",ivar);
          } else if (!(compute->invoked_flag & INVOKED_PERATOM)) {
            compute->compute_peratom();
            compute->invoked_flag |= INVOKED_PERATOM;
          }

          if (compute->array_atom)
            peratom2global(1,NULL,&compute->array_atom[0][index2-1],
                           compute->size_peratom_cols,index1,
                           tree,treestack,ntreestack,argstack,nargstack);
          else
            peratom2global(1,NULL,NULL,
                           compute->size_peratom_cols,index1,
                           tree,treestack,ntreestack,argstack,nargstack);

        // c_ID = vector from per-atom vector

        } else if (nbracket == 0 && compute->peratom_flag &&
                   compute->size_peratom_cols == 0) {

          if (tree == NULL)
            print_var_error(FLERR,"Per-atom compute in "
                            "equal-style variable formula",ivar);
          if (treetype == VECTOR)
            print_var_error(FLERR,"Per-atom compute in "
                            "vector-style variable formula",ivar);
          if (update->whichflag == 0) {
            if (compute->invoked_peratom != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable "
                              "between runs is not current",ivar);
          } else if (!(compute->invoked_flag & INVOKED_PERATOM)) {
            compute->compute_peratom();
            compute->invoked_flag |= INVOKED_PERATOM;
          }

          Tree *newtree = new Tree();
          newtree->type = ATOMARRAY;
          newtree->array = compute->vector_atom;
          newtree->nstride = 1;
          newtree->selfalloc = 0;
          newtree->first = newtree->second = NULL;
          newtree->nextra = 0;
          treestack[ntreestack++] = newtree;

        // c_ID[i] = vector from per-atom array

        } else if (nbracket == 1 && compute->peratom_flag &&
                   compute->size_peratom_cols > 0) {

          if (tree == NULL)
            print_var_error(FLERR,"Per-atom compute in "
                            "equal-style variable formula",ivar);
          if (treetype == VECTOR)
            print_var_error(FLERR,"Per-atom compute in "
                            "vector-style variable formula",ivar);
          if (index1 > compute->size_peratom_cols)
            print_var_error(FLERR,"Variable formula compute array "
                            "is accessed out-of-range",ivar,0);
          if (update->whichflag == 0) {
            if (compute->invoked_peratom != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable "
                              "between runs is not current",ivar);
          } else if (!(compute->invoked_flag & INVOKED_PERATOM)) {
            compute->compute_peratom();
            compute->invoked_flag |= INVOKED_PERATOM;
          }

          Tree *newtree = new Tree();
          newtree->type = ATOMARRAY;
          if (compute->array_atom)
            newtree->array = &compute->array_atom[0][index1-1];
          else
            newtree->array = NULL;
          newtree->nstride = compute->size_peratom_cols;
          newtree->selfalloc = 0;
          newtree->first = newtree->second = NULL;
          newtree->nextra = 0;
          treestack[ntreestack++] = newtree;

        } else if (nbracket == 1 && compute->local_flag) {
          print_var_error(FLERR,"Cannot access local data via indexing",ivar);
        } else print_var_error(FLERR,
                               "Mismatched compute in variable formula",ivar);

      // ----------------
      // fix
      // ----------------

      } else if (strncmp(word,"f_",2) == 0 || strncmp(word,"F_",2) == 0) {
        if (domain->box_exist == 0)
          print_var_error(FLERR,"Variable evaluation before "
                          "simulation box is defined",ivar);

        // uppercase used to force access of
        // global vector vs global scalar, and global array vs global vector

        int lowercase = 1;
        if (word[0] == 'F') lowercase = 0;

        int ifix = modify->find_fix(word+2);
        if (ifix < 0) {
          char msg[128];
          snprintf(msg,128,"Invalid fix ID '%s' in variable formula",word+2);
          print_var_error(FLERR,msg,ivar);
        }
        Fix *fix = modify->fix[ifix];

        // parse zero or one or two trailing brackets
        // point i beyond last bracket
        // nbracket = # of bracket pairs
        // index1,index2 = int inside each bracket pair, possibly an atom ID

        int nbracket;
        tagint index1,index2;
        if (str[i] != '[') nbracket = 0;
        else {
          nbracket = 1;
          ptr = &str[i];
          index1 = int_between_brackets(ptr,1);
          i = ptr-str+1;
          if (str[i] == '[') {
            nbracket = 2;
            ptr = &str[i];
            index2 = int_between_brackets(ptr,1);
            i = ptr-str+1;
          }
        }

        // f_ID = scalar from global scalar, must be lowercase

        if (nbracket == 0 && fix->scalar_flag && lowercase) {

          if (update->whichflag > 0 && update->ntimestep % fix->global_freq)
            print_var_error(FLERR,"Fix in variable not computed "
                            "at a compatible time",ivar);

          value1 = fix->compute_scalar();
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->first = newtree->second = NULL;
            newtree->nextra = 0;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // f_ID[i] = scalar from global vector, must be lowercase

        } else if (nbracket == 1 && fix->vector_flag && lowercase) {

          if (index1 > fix->size_vector &&
              fix->size_vector_variable == 0)
            print_var_error(FLERR,"Variable formula fix vector is "
                            "accessed out-of-range",ivar,0);
          if (update->whichflag > 0 && update->ntimestep % fix->global_freq)
            print_var_error(FLERR,"Fix in variable not computed "
                            "at a compatible time",ivar);

          value1 = fix->compute_vector(index1-1);
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->first = newtree->second = NULL;
            newtree->nextra = 0;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // f_ID[i][j] = scalar from global array, must be lowercase

        } else if (nbracket == 2 && fix->array_flag && lowercase) {

          if (index1 > fix->size_array_rows &&
              fix->size_array_rows_variable == 0)
            print_var_error(FLERR,"Variable formula fix array is "
                            "accessed out-of-range",ivar,0);
          if (index2 > fix->size_array_cols)
            print_var_error(FLERR,"Variable formula fix array is "
                            "accessed out-of-range",ivar,0);
          if (update->whichflag > 0 && update->ntimestep % fix->global_freq)
            print_var_error(FLERR,"Fix in variable not computed at a "
                            "compatible time",ivar);

          value1 = fix->compute_array(index1-1,index2-1);
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->first = newtree->second = NULL;
            newtree->nextra = 0;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // f_ID = vector from global vector, lowercase or uppercase

        } else if (nbracket == 0 && fix->vector_flag) {

          if (update->whichflag > 0 && update->ntimestep % fix->global_freq)
            print_var_error(FLERR,"Fix in variable not computed at "
                            "compatible time",ivar);
          if (tree == NULL)
            print_var_error(FLERR,"Fix global vector in "
                            "equal-style variable formula",ivar);
          if (treetype == ATOM)
            print_var_error(FLERR,"Fix global vector in "
                            "atom-style variable formula",ivar);
          if (fix->size_vector == 0)
            print_var_error(FLERR,"Variable formula "
                            "fix vector is zero length",ivar);

          int nvec = fix->size_vector;
          double *vec;
          memory->create(vec,nvec,"variable:values");
          for (int m = 0; m < nvec; m++)
            vec[m] = fix->compute_vector(m);

          Tree *newtree = new Tree();
          newtree->type = VECTORARRAY;
          newtree->array = vec;
          newtree->nvector = nvec;
          newtree->nstride = 1;
          newtree->selfalloc = 1;
          newtree->first = newtree->second = NULL;
          newtree->nextra = 0;
          treestack[ntreestack++] = newtree;

        // f_ID[i] = vector from global array, lowercase or uppercase

        } else if (nbracket == 1 && fix->array_flag) {

          if (update->whichflag > 0 && update->ntimestep % fix->global_freq)
            print_var_error(FLERR,"Fix in variable not computed "
                            "at a compatible time",ivar);
          if (tree == NULL)
            print_var_error(FLERR,"Fix global vector in "
                            "equal-style variable formula",ivar);
          if (treetype == ATOM)
            print_var_error(FLERR,"Fix global vector in "
                            "atom-style variable formula",ivar);
          if (fix->size_array_rows == 0)
            print_var_error(FLERR,"Variable formula fix array is "
                            "zero length",ivar);

          int nvec = fix->size_array_rows;
          double *vec;
          memory->create(vec,nvec,"variable:values");
          for (int m = 0; m < nvec; m++)
            vec[m] = fix->compute_array(m,index1-1);

          Tree *newtree = new Tree();
          newtree->type = VECTORARRAY;
          newtree->array = vec;
          newtree->nvector = nvec;
          newtree->nstride = 1;
          newtree->selfalloc = 1;
          newtree->first = newtree->second = NULL;
          newtree->nextra = 0;
          treestack[ntreestack++] = newtree;

        // f_ID[i] = scalar from per-atom vector

        } else if (nbracket == 1 && fix->peratom_flag &&
                   fix->size_peratom_cols == 0) {

          if (update->whichflag > 0 &&
              update->ntimestep % fix->peratom_freq)
            print_var_error(FLERR,"Fix in variable not computed "
                            "at a compatible time",ivar);

          peratom2global(1,NULL,fix->vector_atom,1,index1,
                         tree,treestack,ntreestack,argstack,nargstack);

        // f_ID[i][j] = scalar from per-atom array

        } else if (nbracket == 2 && fix->peratom_flag &&
                   fix->size_peratom_cols > 0) {

          if (index2 > fix->size_peratom_cols)
            print_var_error(FLERR,"Variable formula fix array is "
                            "accessed out-of-range",ivar,0);
          if (update->whichflag > 0 &&
              update->ntimestep % fix->peratom_freq)
            print_var_error(FLERR,"Fix in variable not computed "
                            "at a compatible time",ivar);

          if (fix->array_atom)
            peratom2global(1,NULL,&fix->array_atom[0][index2-1],
                           fix->size_peratom_cols,index1,
                           tree,treestack,ntreestack,argstack,nargstack);
          else
            peratom2global(1,NULL,NULL,
                           fix->size_peratom_cols,index1,
                           tree,treestack,ntreestack,argstack,nargstack);

        // f_ID = vector from per-atom vector

        } else if (nbracket == 0 && fix->peratom_flag &&
                   fix->size_peratom_cols == 0) {

          if (tree == NULL)
            print_var_error(FLERR,"Per-atom fix in "
                            "equal-style variable formula",ivar);
          if (update->whichflag > 0 &&
              update->ntimestep % fix->peratom_freq)
            print_var_error(FLERR,"Fix in variable not computed at "
                            "compatible time",ivar);

          Tree *newtree = new Tree();
          newtree->type = ATOMARRAY;
          newtree->array = fix->vector_atom;
          newtree->nstride = 1;
          newtree->selfalloc = 0;
          newtree->first = newtree->second = NULL;
          newtree->nextra = 0;
          treestack[ntreestack++] = newtree;

        // f_ID[i] = vector from per-atom array

        } else if (nbracket == 1 && fix->peratom_flag &&
                   fix->size_peratom_cols > 0) {

          if (tree == NULL)
            print_var_error(FLERR,"Per-atom fix in "
                            "equal-style variable formula",ivar);
          if (index1 > fix->size_peratom_cols)
            print_var_error(FLERR,"Variable formula fix array "
                            "is accessed out-of-range",ivar,0);
          if (update->whichflag > 0 &&
              update->ntimestep % fix->peratom_freq)
            print_var_error(FLERR,"Fix in variable not computed at "
                            "compatible time",ivar);

          Tree *newtree = new Tree();
          newtree->type = ATOMARRAY;
          if (fix->array_atom)
            newtree->array = &fix->array_atom[0][index1-1];
          else
            newtree->array = NULL;
          newtree->nstride = fix->size_peratom_cols;
          newtree->selfalloc = 0;
          newtree->first = newtree->second = NULL;
          newtree->nextra = 0;
          treestack[ntreestack++] = newtree;

        } else print_var_error(FLERR,"Mismatched fix in variable formula",ivar);

      // ----------------
      // variable
      // ----------------

      } else if (strncmp(word,"v_",2) == 0) {

        int ivar = find(word+2);
        if (ivar < 0)
          print_var_error(FLERR,"Invalid variable reference "
                          "in variable formula",ivar);
        if (eval_in_progress[ivar])
          print_var_error(FLERR,"Variable has circular dependency",ivar);

        // parse zero or one trailing brackets
        // point i beyond last bracket
        // nbracket = # of bracket pairs
        // index = int inside bracket, possibly an atom ID

        int nbracket;
        tagint index;
        if (str[i] != '[') nbracket = 0;
        else {
          nbracket = 1;
          ptr = &str[i];
          index = int_between_brackets(ptr,1);
          i = ptr-str+1;
        }

        // v_name = scalar from internal-style variable
        // access value directly

        if (nbracket == 0 && style[ivar] == INTERNAL) {

          value1 = dvalue[ivar];
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->first = newtree->second = NULL;
            newtree->nextra = 0;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // v_name = scalar from non atom/atomfile & non vector-style variable
        // access value via retrieve()

        } else if (nbracket == 0 && style[ivar] != ATOM &&
                   style[ivar] != ATOMFILE && style[ivar] != VECTOR) {

          char *var = retrieve(word+2);
          if (var == NULL)
            print_var_error(FLERR,"Invalid variable evaluation in "
                            "variable formula",ivar);
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = atof(var);
            newtree->first = newtree->second = NULL;
            newtree->nextra = 0;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = atof(var);

        // v_name = per-atom vector from atom-style variable
        // evaluate the atom-style variable as newtree

        } else if (nbracket == 0 && style[ivar] == ATOM) {

          if (tree == NULL)
            print_var_error(FLERR,"Atom-style variable in "
                            "equal-style variable formula",ivar);
          if (treetype == VECTOR)
            print_var_error(FLERR,"Atom-style variable in "
                            "vector-style variable formula",ivar);

          Tree *newtree = NULL;
          evaluate(data[ivar][0],&newtree,ivar);
          treestack[ntreestack++] = newtree;

        // v_name = per-atom vector from atomfile-style variable

        } else if (nbracket == 0 && style[ivar] == ATOMFILE) {

          if (tree == NULL)
            print_var_error(FLERR,"Atomfile-style variable in "
                            "equal-style variable formula",ivar);
          if (treetype == VECTOR)
            print_var_error(FLERR,"Atomfile-style variable in "
                            "vector-style variable formula",ivar);

          Tree *newtree = new Tree();
          newtree->type = ATOMARRAY;
          newtree->array = reader[ivar]->fixstore->vstore;
          newtree->nstride = 1;
          newtree->selfalloc = 0;
          newtree->first = newtree->second = NULL;
          newtree->nextra = 0;
          treestack[ntreestack++] = newtree;

        // v_name = vector from vector-style variable
        // evaluate the vector-style variable, put result in newtree

        } else if (nbracket == 0 && style[ivar] == VECTOR) {

          if (tree == NULL)
            print_var_error(FLERR,"Vector-style variable in "
                            "equal-style variable formula",ivar);
          if (treetype == ATOM)
            print_var_error(FLERR,"Vector-style variable in "
                            "atom-style variable formula",ivar);

          double *vec;
          int nvec = compute_vector(ivar,&vec);

          Tree *newtree = new Tree();
          newtree->type = VECTORARRAY;
          newtree->array = vec;
          newtree->nvector = nvec;
          newtree->nstride = 1;
          newtree->selfalloc = 0;
          newtree->first = newtree->second = NULL;
          newtree->nextra = 0;
          treestack[ntreestack++] = newtree;

        // v_name[N] = scalar from atom-style variable
        // compute the per-atom variable in result
        // use peratom2global to extract single value from result

        } else if (nbracket && style[ivar] == ATOM) {

          double *result;
          memory->create(result,atom->nlocal,"variable:result");
          compute_atom(ivar,0,result,1,0);
          peratom2global(1,NULL,result,1,index,
                         tree,treestack,ntreestack,argstack,nargstack);
          memory->destroy(result);

        // v_name[N] = scalar from atomfile-style variable

        } else if (nbracket && style[ivar] == ATOMFILE) {

          peratom2global(1,NULL,reader[ivar]->fixstore->vstore,1,index,
                         tree,treestack,ntreestack,argstack,nargstack);

        // v_name[N] = scalar from vector-style variable
        // compute the vector-style variable, extract single value

        } else if (nbracket && style[ivar] == VECTOR) {

          double *vec;
          int nvec = compute_vector(ivar,&vec);
          if (index <= 0 || index > nvec)
            print_var_error(FLERR,"Invalid index into "
                            "vector-style variable",ivar);
          int m = index;   // convert from tagint to int

          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = vec[m-1];
            newtree->first = newtree->second = NULL;
            newtree->nextra = 0;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = vec[m-1];

        } else print_var_error(FLERR,"Mismatched variable in "
                               "variable formula",ivar);

      // ----------------
      // math/group/special function or atom value/vector or
      // constant or thermo keyword
      // ----------------

      } else {

        // ----------------
        // math or group or special function
        // ----------------

        if (str[i] == '(') {
          char *contents;
          i = find_matching_paren(str,i,contents,ivar);
          i++;

          if (math_function(word,contents,tree,treestack,ntreestack,
                            argstack,nargstack,ivar));
          else if (group_function(word,contents,tree,treestack,ntreestack,
                                  argstack,nargstack,ivar));
          else if (special_function(word,contents,tree,treestack,ntreestack,
                                    argstack,nargstack,ivar));
          else {
            char msg[128];
            snprintf(msg,128,"Invalid math/group/special function '%s()'"
                     "in variable formula", word);
            print_var_error(FLERR,msg,ivar);
          }
          delete [] contents;

        // ----------------
        // atom value
        // ----------------

        } else if (str[i] == '[') {
          if (domain->box_exist == 0)
            print_var_error(FLERR,"Variable evaluation before "
                            "simulation box is defined",ivar);

          ptr = &str[i];
          tagint id = int_between_brackets(ptr,1);
          i = ptr-str+1;

          peratom2global(0,word,NULL,0,id,
                         tree,treestack,ntreestack,argstack,nargstack);

        // ----------------
        // atom vector
        // ----------------

        } else if (is_atom_vector(word)) {
          if (domain->box_exist == 0)
            print_var_error(FLERR,"Variable evaluation before "
                            "simulation box is defined",ivar);

          atom_vector(word,tree,treestack,ntreestack);

        // ----------------
        // constant
        // ----------------

        } else if (is_constant(word)) {
          value1 = constant(word);
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->first = newtree->second = NULL;
            newtree->nextra = 0;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // ----------------
        // thermo keyword
        // ----------------

        } else {
          if (domain->box_exist == 0)
            print_var_error(FLERR,"Variable evaluation before "
                            "simulation box is defined",ivar);

          int flag = output->thermo->evaluate_keyword(word,&value1);
          if (flag) {
            char msg[128];
            snprintf(msg,128,"Invalid thermo keyword '%s' in variable formula",word);
            print_var_error(FLERR,msg,ivar);
          }
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->first = newtree->second = NULL;
            newtree->nextra = 0;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;
        }
      }

      delete [] word;

    // ----------------
    // math operator, including end-of-string
    // ----------------

    } else if (strchr("+-*/^<>=!&|%\0",onechar)) {
      if (onechar == '+') op = ADD;
      else if (onechar == '-') op = SUBTRACT;
      else if (onechar == '*') op = MULTIPLY;
      else if (onechar == '/') op = DIVIDE;
      else if (onechar == '%') op = MODULO;
      else if (onechar == '^') op = CARAT;
      else if (onechar == '=') {
        if (str[i+1] != '=')
          print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
        op = EQ;
        i++;
      } else if (onechar == '!') {
        if (str[i+1] == '=') {
          op = NE;
          i++;
        } else op = NOT;
      } else if (onechar == '<') {
        if (str[i+1] != '=') op = LT;
        else {
          op = LE;
          i++;
        }
      } else if (onechar == '>') {
        if (str[i+1] != '=') op = GT;
        else {
          op = GE;
          i++;
        }
      } else if (onechar == '&') {
        if (str[i+1] != '&')
          print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
        op = AND;
        i++;
      } else if (onechar == '|') {
        if (str[i+1] == '|') op = OR;
        else if (str[i+1] == '^') op = XOR;
        else print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
        i++;
      } else op = DONE;

      i++;

      if (op == SUBTRACT && expect == ARG) {
        opstack[nopstack++] = UNARY;
        continue;
      }
      if (op == NOT && expect == ARG) {
        opstack[nopstack++] = op;
        continue;
      }

      if (expect == ARG)
        print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
      expect = ARG;

      // evaluate stack as deep as possible while respecting precedence
      // before pushing current op onto stack

      while (nopstack && precedence[opstack[nopstack-1]] >= precedence[op]) {
        opprevious = opstack[--nopstack];

        if (tree) {
          Tree *newtree = new Tree();
          newtree->type = opprevious;
          if ((opprevious == UNARY) || (opprevious == NOT)) {
            newtree->first = treestack[--ntreestack];
            newtree->second = NULL;
            newtree->nextra = 0;
          } else {
            newtree->second = treestack[--ntreestack];
            newtree->first = treestack[--ntreestack];
            newtree->nextra = 0;
          }
          treestack[ntreestack++] = newtree;

        } else {
          value2 = argstack[--nargstack];
          if (opprevious != UNARY && opprevious != NOT)
            value1 = argstack[--nargstack];

          if (opprevious == ADD)
            argstack[nargstack++] = value1 + value2;
          else if (opprevious == SUBTRACT)
            argstack[nargstack++] = value1 - value2;
          else if (opprevious == MULTIPLY)
            argstack[nargstack++] = value1 * value2;
          else if (opprevious == DIVIDE) {
            if (value2 == 0.0)
              print_var_error(FLERR,"Divide by 0 in variable formula",ivar,0);
            argstack[nargstack++] = value1 / value2;
          } else if (opprevious == MODULO) {
            if (value2 == 0.0)
              print_var_error(FLERR,"Modulo 0 in variable formula",ivar,0);
            argstack[nargstack++] = fmod(value1,value2);
          } else if (opprevious == CARAT) {
            if (value2 == 0.0)
              argstack[nargstack++] = 1.0;
            else if ((value1 == 0.0) && (value2 < 0.0))
              print_var_error(FLERR,"Invalid power expression in "
                              "variable formula",ivar,0);
            else argstack[nargstack++] = pow(value1,value2);
          } else if (opprevious == UNARY) {
            argstack[nargstack++] = -value2;
          } else if (opprevious == NOT) {
            if (value2 == 0.0) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == EQ) {
            if (value1 == value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == NE) {
            if (value1 != value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == LT) {
            if (value1 < value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == LE) {
            if (value1 <= value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == GT) {
            if (value1 > value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == GE) {
            if (value1 >= value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == AND) {
            if (value1 != 0.0 && value2 != 0.0) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == OR) {
            if (value1 != 0.0 || value2 != 0.0) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == XOR) {
            if ((value1 == 0.0 && value2 != 0.0) ||
                (value1 != 0.0 && value2 == 0.0)) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          }
        }
      }

      // if end-of-string, break out of entire formula evaluation loop

      if (op == DONE) break;

      // push current operation onto stack

      opstack[nopstack++] = op;

    } else print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
  }

  if (nopstack) print_var_error(FLERR,"Invalid syntax in variable formula",ivar);

  // for atom-style variable, return remaining tree
  // for equal-style variable, return remaining arg

  if (tree) {
    if (ntreestack != 1)
      print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
    *tree = treestack[0];
    return 0.0;
  } else {
    if (nargstack != 1)
      print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
    return argstack[0];
  }
}

/* ----------------------------------------------------------------------
   one-time collapse of an atom-style variable parse tree
   tree was created by one-time parsing of formula string via evaluate()
   only keep tree nodes that depend on
     ATOMARRAY, TYPEARRAY, INTARRAY, BIGINTARRAY, VECTOR
   remainder is converted to single VALUE
   this enables optimal eval_tree loop over atoms
   customize by adding a function:
     sqrt(),exp(),ln(),log(),abs(),sin(),cos(),tan(),asin(),acos(),atan(),
     atan2(y,x),random(x,y,z),normal(x,y,z),ceil(),floor(),round(),
     ramp(x,y),stagger(x,y),logfreq(x,y,z),logfreq2(x,y,z),
     logfreq3(x,y,z),stride(x,y,z),vdisplace(x,y),swiggle(x,y,z),
     cwiggle(x,y,z),gmask(x),rmask(x),grmask(x,y)
---------------------------------------------------------------------- */

double Variable::collapse_tree(Tree *tree)
{
  double arg1,arg2;

  if (tree->type == VALUE) return tree->value;
  if (tree->type == ATOMARRAY) return 0.0;
  if (tree->type == TYPEARRAY) return 0.0;
  if (tree->type == INTARRAY) return 0.0;
  if (tree->type == BIGINTARRAY) return 0.0;
  if (tree->type == VECTORARRAY) return 0.0;

  if (tree->type == ADD) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = arg1 + arg2;
    return tree->value;
  }

  if (tree->type == SUBTRACT) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = arg1 - arg2;
    return tree->value;
  }

  if (tree->type == MULTIPLY) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = arg1 * arg2;
    return tree->value;
  }

  if (tree->type == DIVIDE) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg2 == 0.0) error->one(FLERR,"Divide by 0 in variable formula");
    tree->value = arg1 / arg2;
    return tree->value;
  }

  if (tree->type == MODULO) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg2 == 0.0) error->one(FLERR,"Modulo 0 in variable formula");
    tree->value = fmod(arg1,arg2);
    return tree->value;
  }

  if (tree->type == CARAT) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg2 == 0.0) error->one(FLERR,"Power by 0 in variable formula");
    tree->value = pow(arg1,arg2);
    return tree->value;
  }

  if (tree->type == UNARY) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = -arg1;
    return tree->value;
  }

  if (tree->type == NOT) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 == 0.0) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == EQ) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 == arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == NE) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 != arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == LT) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 < arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == LE) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 <= arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == GT) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 > arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == GE) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 >= arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == AND) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 != 0.0 && arg2 != 0.0) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == OR) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 != 0.0 || arg2 != 0.0) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == XOR) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if ((arg1 == 0.0 && arg2 != 0.0) || (arg1 != 0.0 && arg2 == 0.0))
      tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == SQRT) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 < 0.0)
      error->one(FLERR,"Sqrt of negative value in variable formula");
    tree->value = sqrt(arg1);
    return tree->value;
  }

  if (tree->type == EXP) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = exp(arg1);
    return tree->value;
  }

  if (tree->type == LN) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 <= 0.0)
      error->one(FLERR,"Log of zero/negative value in variable formula");
    tree->value = log(arg1);
    return tree->value;
  }

  if (tree->type == LOG) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 <= 0.0)
      error->one(FLERR,"Log of zero/negative value in variable formula");
    tree->value = log10(arg1);
    return tree->value;
  }

  if (tree->type == ABS) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = fabs(arg1);
    return tree->value;
  }

  if (tree->type == SIN) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = sin(arg1);
    return tree->value;
  }

  if (tree->type == COS) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = cos(arg1);
    return tree->value;
  }

  if (tree->type == TAN) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = tan(arg1);
    return tree->value;
  }

  if (tree->type == ASIN) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 < -1.0 || arg1 > 1.0)
      error->one(FLERR,"Arcsin of invalid value in variable formula");
    tree->value = asin(arg1);
    return tree->value;
  }

  if (tree->type == ACOS) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 < -1.0 || arg1 > 1.0)
      error->one(FLERR,"Arccos of invalid value in variable formula");
    tree->value = acos(arg1);
    return tree->value;
  }

  if (tree->type == ATAN) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = atan(arg1);
    return tree->value;
  }

  if (tree->type == ATAN2) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = atan2(arg1,arg2);
    return tree->value;
  }

  // random() or normal() do not become a single collapsed value

  if (tree->type == RANDOM) {
    collapse_tree(tree->first);
    collapse_tree(tree->second);
    if (randomatom == NULL) {
      int seed = static_cast<int> (collapse_tree(tree->extra[0]));
      if (seed <= 0)
        error->one(FLERR,"Invalid math function in variable formula");
      randomatom = new RanMars(lmp,seed+me);
    }
    return 0.0;
  }

  if (tree->type == NORMAL) {
    collapse_tree(tree->first);
    double sigma = collapse_tree(tree->second);
    if (sigma < 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    if (randomatom == NULL) {
      int seed = static_cast<int> (collapse_tree(tree->extra[0]));
      if (seed <= 0)
        error->one(FLERR,"Invalid math function in variable formula");
      randomatom = new RanMars(lmp,seed+me);
    }
    return 0.0;
  }

  if (tree->type == CEIL) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = ceil(arg1);
    return tree->value;
  }

  if (tree->type == FLOOR) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = floor(arg1);
    return tree->value;
  }

  if (tree->type == ROUND) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = MYROUND(arg1);
    return tree->value;
  }

  if (tree->type == RAMP) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    double delta = update->ntimestep - update->beginstep;
    if (delta != 0.0) delta /= update->endstep - update->beginstep;
    tree->value = arg1 + delta*(arg2-arg1);
    return tree->value;
  }

  if (tree->type == STAGGER) {
    int ivalue1 = static_cast<int> (collapse_tree(tree->first));
    int ivalue2 = static_cast<int> (collapse_tree(tree->second));
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue1 <= ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    int lower = update->ntimestep/ivalue1 * ivalue1;
    int delta = update->ntimestep - lower;
    if (delta < ivalue2) tree->value = lower+ivalue2;
    else tree->value = lower+ivalue1;
    return tree->value;
  }

  if (tree->type == LOGFREQ) {
    int ivalue1 = static_cast<int> (collapse_tree(tree->first));
    int ivalue2 = static_cast<int> (collapse_tree(tree->second));
    int ivalue3 = static_cast<int> (collapse_tree(tree->extra[0]));
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 || ivalue2 >= ivalue3)
      error->one(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) tree->value = ivalue1;
    else {
      int lower = ivalue1;
      while (update->ntimestep >= ivalue3*lower) lower *= ivalue3;
      int multiple = update->ntimestep/lower;
      if (multiple < ivalue2) tree->value = (multiple+1)*lower;
      else tree->value = lower*ivalue3;
    }
    return tree->value;
  }

  if (tree->type == LOGFREQ2) {
    int ivalue1 = static_cast<int> (collapse_tree(tree->first));
    int ivalue2 = static_cast<int> (collapse_tree(tree->second));
    int ivalue3 = static_cast<int> (collapse_tree(tree->extra[0]));
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 )
      error->all(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) tree->value = ivalue1;
    else {
      tree->value = ivalue1;
      double delta = ivalue1*(ivalue3-1.0)/ivalue2;
      int count = 0;
      while (update->ntimestep >= tree->value) {
        tree->value += delta;
        count++;
        if (count % ivalue2 == 0) delta *= ivalue3;
      }
    }
    tree->value = ceil(tree->value);
    return tree->value;
  }

  if (tree->type == LOGFREQ3) {
    int ivalue1 = static_cast<int> (collapse_tree(tree->first));
    int ivalue2 = static_cast<int> (collapse_tree(tree->second));
    int ivalue3 = static_cast<int> (collapse_tree(tree->extra[0]));
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (ivalue1 <= 0 || ivalue2 <= 1 || ivalue3 <= 0 ||
        ivalue3-ivalue1+1 < ivalue2 )
      error->all(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) tree->value = ivalue1;
    //else if (update->ntimestep <= ivalue3){
    else {
      tree->value = ivalue1;
      double logsp = ivalue1;
      double factor = pow(((double)ivalue3)/ivalue1, 1.0/(ivalue2-1));
      int linsp = ivalue1;
      while (update->ntimestep >= (tree->value)) {
        logsp *= factor;
        linsp++;
        if (linsp > logsp) tree->value = linsp;
        else tree->value = ceil(logsp)-(((int)ceil(logsp)-1)/ivalue3);
      }
    }
    if (update->ntimestep > ivalue3)
      error->all(FLERR,"Calls to variable exceeded limit");
    return tree->value;
  }

  if (tree->type == STRIDE) {
    int ivalue1 = static_cast<int> (collapse_tree(tree->first));
    int ivalue2 = static_cast<int> (collapse_tree(tree->second));
    int ivalue3 = static_cast<int> (collapse_tree(tree->extra[0]));
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) tree->value = ivalue1;
    else if (update->ntimestep < ivalue2) {
      int offset = update->ntimestep - ivalue1;
      tree->value = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
      if (tree->value > ivalue2) tree->value = MAXBIGINT;
    } else tree->value = MAXBIGINT;
    return tree->value;
  }

  if (tree->type == STRIDE2) {
    int ivalue1 = static_cast<int> (collapse_tree(tree->first));
    int ivalue2 = static_cast<int> (collapse_tree(tree->second));
    int ivalue3 = static_cast<int> (collapse_tree(tree->extra[0]));
    int ivalue4 = static_cast<int> (collapse_tree(tree->extra[1]));
    int ivalue5 = static_cast<int> (collapse_tree(tree->extra[2]));
    int ivalue6 = static_cast<int> (collapse_tree(tree->extra[3]));
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE || tree->extra[1]->type != VALUE ||
        tree->extra[2]->type != VALUE || tree->extra[3]->type != VALUE)
      return 0.0;
    tree->type = VALUE;
    if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    if (ivalue4 < 0 || ivalue5 < 0 || ivalue6 <= 0 || ivalue4 > ivalue5)
      error->one(FLERR,"Invalid math function in variable formula");
    if (ivalue4 < ivalue1 || ivalue5 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    bigint istep;
    if (update->ntimestep < ivalue1) istep = ivalue1;
    else if (update->ntimestep < ivalue2) {
      if (update->ntimestep < ivalue4 || update->ntimestep > ivalue5) {
        int offset = update->ntimestep - ivalue1;
        istep = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
        if (update->ntimestep < ivalue2 && istep > ivalue4)
          tree->value = ivalue4;
      } else {
        int offset = update->ntimestep - ivalue4;
        istep = ivalue4 + (offset/ivalue6)*ivalue6 + ivalue6;
        if (istep > ivalue5) {
          int offset = ivalue5 - ivalue1;
          istep = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
          if (istep > ivalue2) istep = MAXBIGINT;
        }
      }
    } else istep = MAXBIGINT;
    tree->value = istep;
    return tree->value;
  }

  if (tree->type == VDISPLACE) {
    double arg1 = collapse_tree(tree->first);
    double arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    double delta = update->ntimestep - update->beginstep;
    tree->value = arg1 + arg2*delta*update->dt;
    return tree->value;
  }

  if (tree->type == SWIGGLE) {
    double arg1 = collapse_tree(tree->first);
    double arg2 = collapse_tree(tree->second);
    double arg3 = collapse_tree(tree->extra[0]);
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg3 == 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    double delta = update->ntimestep - update->beginstep;
    double omega = 2.0*MY_PI/arg3;
    tree->value = arg1 + arg2*sin(omega*delta*update->dt);
    return tree->value;
  }

  if (tree->type == CWIGGLE) {
    double arg1 = collapse_tree(tree->first);
    double arg2 = collapse_tree(tree->second);
    double arg3 = collapse_tree(tree->extra[0]);
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg3 == 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    double delta = update->ntimestep - update->beginstep;
    double omega = 2.0*MY_PI/arg3;
    tree->value = arg1 + arg2*(1.0-cos(omega*delta*update->dt));
    return tree->value;
  }

  // mask functions do not become a single collapsed value

  if (tree->type == GMASK) return 0.0;
  if (tree->type == RMASK) return 0.0;
  if (tree->type == GRMASK) return 0.0;

  return 0.0;
}

/* ----------------------------------------------------------------------
   evaluate an atom-style or vector-style variable parse tree
   index I = atom I or vector index I
   tree was created by one-time parsing of formula string via evaluate()
   customize by adding a function:
     sqrt(),exp(),ln(),log(),sin(),cos(),tan(),asin(),acos(),atan(),
     atan2(y,x),random(x,y,z),normal(x,y,z),ceil(),floor(),round(),
     ramp(x,y),stagger(x,y),logfreq(x,y,z),logfreq2(x,y,z),
     logfreq3(x,y,z),stride(x,y,z),stride2(x,y,z),vdisplace(x,y),
     swiggle(x,y,z),cwiggle(x,y,z),gmask(x),rmask(x),grmask(x,y)
---------------------------------------------------------------------- */

double Variable::eval_tree(Tree *tree, int i)
{
  double arg,arg1,arg2,arg3;

  if (tree->type == VALUE) return tree->value;
  if (tree->type == ATOMARRAY) return tree->array[i*tree->nstride];
  if (tree->type == TYPEARRAY) return tree->array[atom->type[i]];
  if (tree->type == INTARRAY) return (double) tree->iarray[i*tree->nstride];
  if (tree->type == BIGINTARRAY) return (double) tree->barray[i*tree->nstride];
  if (tree->type == VECTORARRAY) return tree->array[i*tree->nstride];

  if (tree->type == ADD)
    return eval_tree(tree->first,i) + eval_tree(tree->second,i);
  if (tree->type == SUBTRACT)
    return eval_tree(tree->first,i) - eval_tree(tree->second,i);
  if (tree->type == MULTIPLY)
    return eval_tree(tree->first,i) * eval_tree(tree->second,i);
  if (tree->type == DIVIDE) {
    double denom = eval_tree(tree->second,i);
    if (denom == 0.0) error->one(FLERR,"Divide by 0 in variable formula");
    return eval_tree(tree->first,i) / denom;
  }
  if (tree->type == MODULO) {
    double denom = eval_tree(tree->second,i);
    if (denom == 0.0) error->one(FLERR,"Modulo 0 in variable formula");
    return fmod(eval_tree(tree->first,i),denom);
  }
  if (tree->type == CARAT) {
    double exponent = eval_tree(tree->second,i);
    if (exponent == 0.0) error->one(FLERR,"Power by 0 in variable formula");
    return pow(eval_tree(tree->first,i),exponent);
  }
  if (tree->type == UNARY) return -eval_tree(tree->first,i);

  if (tree->type == NOT) {
    if (eval_tree(tree->first,i) == 0.0) return 1.0;
    else return 0.0;
  }
  if (tree->type == EQ) {
    if (eval_tree(tree->first,i) == eval_tree(tree->second,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == NE) {
    if (eval_tree(tree->first,i) != eval_tree(tree->second,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == LT) {
    if (eval_tree(tree->first,i) < eval_tree(tree->second,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == LE) {
    if (eval_tree(tree->first,i) <= eval_tree(tree->second,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == GT) {
    if (eval_tree(tree->first,i) > eval_tree(tree->second,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == GE) {
    if (eval_tree(tree->first,i) >= eval_tree(tree->second,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == AND) {
    if (eval_tree(tree->first,i) != 0.0 && eval_tree(tree->second,i) != 0.0)
      return 1.0;
    else return 0.0;
  }
  if (tree->type == OR) {
    if (eval_tree(tree->first,i) != 0.0 || eval_tree(tree->second,i) != 0.0)
      return 1.0;
    else return 0.0;
  }
  if (tree->type == XOR) {
    if ((eval_tree(tree->first,i) == 0.0 && eval_tree(tree->second,i) != 0.0)
        ||
        (eval_tree(tree->first,i) != 0.0 && eval_tree(tree->second,i) == 0.0))
      return 1.0;
    else return 0.0;
  }

  if (tree->type == SQRT) {
    arg1 = eval_tree(tree->first,i);
    if (arg1 < 0.0)
      error->one(FLERR,"Sqrt of negative value in variable formula");
    return sqrt(arg1);
  }
  if (tree->type == EXP)
    return exp(eval_tree(tree->first,i));
  if (tree->type == LN) {
    arg1 = eval_tree(tree->first,i);
    if (arg1 <= 0.0)
      error->one(FLERR,"Log of zero/negative value in variable formula");
    return log(arg1);
  }
  if (tree->type == LOG) {
    arg1 = eval_tree(tree->first,i);
    if (arg1 <= 0.0)
      error->one(FLERR,"Log of zero/negative value in variable formula");
    return log10(arg1);
  }
  if (tree->type == ABS)
    return fabs(eval_tree(tree->first,i));

  if (tree->type == SIN)
    return sin(eval_tree(tree->first,i));
  if (tree->type == COS)
    return cos(eval_tree(tree->first,i));
  if (tree->type == TAN)
    return tan(eval_tree(tree->first,i));

  if (tree->type == ASIN) {
    arg1 = eval_tree(tree->first,i);
    if (arg1 < -1.0 || arg1 > 1.0)
      error->one(FLERR,"Arcsin of invalid value in variable formula");
    return asin(arg1);
  }
  if (tree->type == ACOS) {
    arg1 = eval_tree(tree->first,i);
    if (arg1 < -1.0 || arg1 > 1.0)
      error->one(FLERR,"Arccos of invalid value in variable formula");
    return acos(arg1);
  }
  if (tree->type == ATAN)
    return atan(eval_tree(tree->first,i));
  if (tree->type == ATAN2)
    return atan2(eval_tree(tree->first,i),eval_tree(tree->second,i));

  if (tree->type == RANDOM) {
    double lower = eval_tree(tree->first,i);
    double upper = eval_tree(tree->second,i);
    if (randomatom == NULL) {
      int seed = static_cast<int> (eval_tree(tree->extra[0],i));
      if (seed <= 0)
        error->one(FLERR,"Invalid math function in variable formula");
      randomatom = new RanMars(lmp,seed+me);
    }
    return randomatom->uniform()*(upper-lower)+lower;
  }
  if (tree->type == NORMAL) {
    double mu = eval_tree(tree->first,i);
    double sigma = eval_tree(tree->second,i);
    if (sigma < 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    if (randomatom == NULL) {
      int seed = static_cast<int> (eval_tree(tree->extra[0],i));
      if (seed <= 0)
        error->one(FLERR,"Invalid math function in variable formula");
      randomatom = new RanMars(lmp,seed+me);
    }
    return mu + sigma*randomatom->gaussian();
  }

  if (tree->type == CEIL)
    return ceil(eval_tree(tree->first,i));
  if (tree->type == FLOOR)
    return floor(eval_tree(tree->first,i));
  if (tree->type == ROUND)
    return MYROUND(eval_tree(tree->first,i));

  if (tree->type == RAMP) {
    arg1 = eval_tree(tree->first,i);
    arg2 = eval_tree(tree->second,i);
    double delta = update->ntimestep - update->beginstep;
    if (delta != 0.0) delta /= update->endstep - update->beginstep;
    arg = arg1 + delta*(arg2-arg1);
    return arg;
  }

  if (tree->type == STAGGER) {
    int ivalue1 = static_cast<int> (eval_tree(tree->first,i));
    int ivalue2 = static_cast<int> (eval_tree(tree->second,i));
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue1 <= ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    int lower = update->ntimestep/ivalue1 * ivalue1;
    int delta = update->ntimestep - lower;
    if (delta < ivalue2) arg = lower+ivalue2;
    else arg = lower+ivalue1;
    return arg;
  }

  if (tree->type == LOGFREQ) {
    int ivalue1 = static_cast<int> (eval_tree(tree->first,i));
    int ivalue2 = static_cast<int> (eval_tree(tree->second,i));
    int ivalue3 = static_cast<int> (eval_tree(tree->extra[0],i));
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 || ivalue2 >= ivalue3)
      error->one(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) arg = ivalue1;
    else {
      int lower = ivalue1;
      while (update->ntimestep >= ivalue3*lower) lower *= ivalue3;
      int multiple = update->ntimestep/lower;
      if (multiple < ivalue2) arg = (multiple+1)*lower;
      else arg = lower*ivalue3;
    }
    return arg;
  }

  if (tree->type == LOGFREQ2) {
    int ivalue1 = static_cast<int> (eval_tree(tree->first,i));
    int ivalue2 = static_cast<int> (eval_tree(tree->second,i));
    int ivalue3 = static_cast<int> (eval_tree(tree->extra[0],i));
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 )
      error->all(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) arg = ivalue1;
    else {
      arg = ivalue1;
      double delta = ivalue1*(ivalue3-1.0)/ivalue2;
      int count = 0;
      while (update->ntimestep >= arg) {
        arg += delta;
        count++;
        if (count % ivalue2 == 0) delta *= ivalue3;
      }
    }
    arg = ceil(arg);
    return arg;
  }

  if (tree->type == STRIDE) {
    int ivalue1 = static_cast<int> (eval_tree(tree->first,i));
    int ivalue2 = static_cast<int> (eval_tree(tree->second,i));
    int ivalue3 = static_cast<int> (eval_tree(tree->extra[0],i));
    if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) arg = ivalue1;
    else if (update->ntimestep < ivalue2) {
      int offset = update->ntimestep - ivalue1;
      arg = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
      if (arg > ivalue2) arg = MAXBIGINT;
    } else arg = MAXBIGINT;
    return arg;
  }

  if (tree->type == STRIDE2) {
    int ivalue1 = static_cast<int> (eval_tree(tree->first,i));
    int ivalue2 = static_cast<int> (eval_tree(tree->second,i));
    int ivalue3 = static_cast<int> (eval_tree(tree->extra[0],i));
    int ivalue4 = static_cast<int> (eval_tree(tree->extra[1],i));
    int ivalue5 = static_cast<int> (eval_tree(tree->extra[2],i));
    int ivalue6 = static_cast<int> (eval_tree(tree->extra[3],i));
    if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    if (ivalue4 < 0 || ivalue5 < 0 || ivalue6 <= 0 || ivalue4 > ivalue5)
      error->one(FLERR,"Invalid math function in variable formula");
    if (ivalue4 < ivalue1 || ivalue5 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    bigint istep;
    if (update->ntimestep < ivalue1) istep = ivalue1;
    else if (update->ntimestep < ivalue2) {
      if (update->ntimestep < ivalue4 || update->ntimestep > ivalue5) {
        int offset = update->ntimestep - ivalue1;
        istep = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
        if (update->ntimestep < ivalue2 && istep > ivalue4)
          tree->value = ivalue4;
      } else {
        int offset = update->ntimestep - ivalue4;
        istep = ivalue4 + (offset/ivalue6)*ivalue6 + ivalue6;
        if (istep > ivalue5) {
          int offset = ivalue5 - ivalue1;
          istep = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
          if (istep > ivalue2) istep = MAXBIGINT;
        }
      }
    } else istep = MAXBIGINT;
    arg = istep;
    return arg;
  }

  if (tree->type == VDISPLACE) {
    arg1 = eval_tree(tree->first,i);
    arg2 = eval_tree(tree->second,i);
    double delta = update->ntimestep - update->beginstep;
    arg = arg1 + arg2*delta*update->dt;
    return arg;
  }

  if (tree->type == SWIGGLE) {
    arg1 = eval_tree(tree->first,i);
    arg2 = eval_tree(tree->second,i);
    arg3 = eval_tree(tree->extra[0],i);
    if (arg3 == 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    double delta = update->ntimestep - update->beginstep;
    double omega = 2.0*MY_PI/arg3;
    arg = arg1 + arg2*sin(omega*delta*update->dt);
    return arg;
  }

  if (tree->type == CWIGGLE) {
    arg1 = eval_tree(tree->first,i);
    arg2 = eval_tree(tree->second,i);
    arg3 = eval_tree(tree->extra[0],i);
    if (arg3 == 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    double delta = update->ntimestep - update->beginstep;
    double omega = 2.0*MY_PI/arg3;
    arg = arg1 + arg2*(1.0-cos(omega*delta*update->dt));
    return arg;
  }

  if (tree->type == GMASK) {
    if (atom->mask[i] & tree->ivalue1) return 1.0;
    else return 0.0;
  }

  if (tree->type == RMASK) {
    if (domain->regions[tree->ivalue1]->match(atom->x[i][0],
                                              atom->x[i][1],
                                              atom->x[i][2])) return 1.0;
    else return 0.0;
  }

  if (tree->type == GRMASK) {
    if ((atom->mask[i] & tree->ivalue1) &&
        (domain->regions[tree->ivalue2]->match(atom->x[i][0],
                                               atom->x[i][1],
                                               atom->x[i][2]))) return 1.0;
    else return 0.0;
  }

  return 0.0;
}

/* ----------------------------------------------------------------------
   scan entire tree, find size of vectors for vector-style variable
   return N for consistent vector size
   return 0 for no vector size, caller flags as error
   return -1 for inconsistent vector size, caller flags as error
------------------------------------------------------------------------- */

int Variable::size_tree_vector(Tree *tree)
{
  int nsize = 0;
  if (tree->type == VECTORARRAY) nsize = tree->nvector;
  if (tree->first) nsize = compare_tree_vector(nsize,
                                               size_tree_vector(tree->first));
  if (tree->second) nsize = compare_tree_vector(nsize,
                                                size_tree_vector(tree->second));
  if (tree->nextra) {
    for (int i = 0; i < tree->nextra; i++)
      nsize = compare_tree_vector(nsize,size_tree_vector(tree->extra[i]));
  }
  return nsize;
}

/* ----------------------------------------------------------------------
   compare size of two vectors for vector-style variable
   return positive size if same or one has no size 0
   return -1 error if one is already error or not same positive size
------------------------------------------------------------------------- */

int Variable::compare_tree_vector(int i, int j)
{
  if (i < 0 || j < 0) return -1;
  if (i == 0 || j == 0) return MAX(i,j);
  if (i != j) return -1;
  return i;
}

/* ---------------------------------------------------------------------- */

void Variable::free_tree(Tree *tree)
{
  if (tree->first) free_tree(tree->first);
  if (tree->second) free_tree(tree->second);
  if (tree->nextra) {
    for (int i = 0; i < tree->nextra; i++) free_tree(tree->extra[i]);
    delete [] tree->extra;
  }

  if (tree->selfalloc) memory->destroy(tree->array);
  delete tree;
}

/* ----------------------------------------------------------------------
   find matching parenthesis in str, allocate contents = str between parens
   i = left paren
   return loc or right paren
------------------------------------------------------------------------- */

int Variable::find_matching_paren(char *str, int i, char *&contents, int ivar)
{
  // istop = matching ')' at same level, allowing for nested parens

  int istart = i;
  int ilevel = 0;
  while (1) {
    i++;
    if (!str[i]) break;
    if (str[i] == '(') ilevel++;
    else if (str[i] == ')' && ilevel) ilevel--;
    else if (str[i] == ')') break;
  }
  if (!str[i]) print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
  int istop = i;

  int n = istop - istart - 1;
  contents = new char[n+1];
  strncpy(contents,&str[istart+1],n);
  contents[n] = '\0';

  return istop;
}

/* ----------------------------------------------------------------------
   find int between brackets and return it
   return a tagint, since value can be an atom ID
   ptr initially points to left bracket
   return it pointing to right bracket
   error if no right bracket or brackets are empty or index = 0
   if varallow = 0: error if any between-bracket chars are non-digits
   if varallow = 1: also allow for v_name, where name is variable name
------------------------------------------------------------------------- */

tagint Variable::int_between_brackets(char *&ptr, int varallow)
{
  int varflag;
  tagint index;

  char *start = ++ptr;

  if (varallow && strstr(ptr,"v_") == ptr) {
    varflag = 1;
    while (*ptr && *ptr != ']') {
      if (!isalnum(*ptr) && *ptr != '_')
        error->all(FLERR,"Variable name between brackets must be "
                   "alphanumeric or underscore characters");
      ptr++;
    }

  } else {
    varflag = 0;
    while (*ptr && *ptr != ']') {
      if (!isdigit(*ptr))
        error->all(FLERR,"Non digit character between brackets in variable");
      ptr++;
    }
  }

  if (*ptr != ']') error->all(FLERR,"Mismatched brackets in variable");
  if (ptr == start) error->all(FLERR,"Empty brackets in variable");

  *ptr = '\0';

  // evaluate index as floating point variable or as tagint via ATOTAGINT()

  if (varflag) {
    char *id = start+2;
    int ivar = find(id);
    if (ivar < 0)
      error->all(FLERR,"Invalid variable name in variable formula");

    char *var = retrieve(id);
    if (var == NULL)
      error->all(FLERR,"Invalid variable evaluation in variable formula");
    index = static_cast<tagint> (atof(var));

  } else index = ATOTAGINT(start);

  *ptr = ']';

  if (index == 0)
    error->all(FLERR,"Index between variable brackets must be positive");
  return index;
}

/* ----------------------------------------------------------------------
   process a math function in formula
   push result onto tree or arg stack
   word = math function
   contents = str between parentheses with comma-separated args
   return 0 if not a match, 1 if successfully processed
   customize by adding a math function:
     sqrt(),exp(),ln(),log(),abs(),sin(),cos(),tan(),asin(),acos(),atan(),
     atan2(y,x),random(x,y,z),normal(x,y,z),ceil(),floor(),round(),
     ramp(x,y),stagger(x,y),logfreq(x,y,z),logfreq2(x,y,z),
     logfreq3(x,y,z),stride(x,y,z),stride2(x,y,z,a,b,c),vdisplace(x,y),
     swiggle(x,y,z),cwiggle(x,y,z)
------------------------------------------------------------------------- */

int Variable::math_function(char *word, char *contents, Tree **tree,
                            Tree **treestack, int &ntreestack,
                            double *argstack, int &nargstack, int ivar)
{
  // word not a match to any math function

  if (strcmp(word,"sqrt") && strcmp(word,"exp") &&
      strcmp(word,"ln") && strcmp(word,"log") &&
      strcmp(word,"abs") &&
      strcmp(word,"sin") && strcmp(word,"cos") &&
      strcmp(word,"tan") && strcmp(word,"asin") &&
      strcmp(word,"acos") && strcmp(word,"atan") &&
      strcmp(word,"atan2") && strcmp(word,"random") &&
      strcmp(word,"normal") && strcmp(word,"ceil") &&
      strcmp(word,"floor") && strcmp(word,"round") &&
      strcmp(word,"ramp") && strcmp(word,"stagger") &&
      strcmp(word,"logfreq") && strcmp(word,"logfreq2") &&
      strcmp(word,"logfreq3") && strcmp(word,"stride") &&
      strcmp(word,"stride2") && strcmp(word,"vdisplace") &&
      strcmp(word,"swiggle") && strcmp(word,"cwiggle"))
    return 0;

  // parse contents for comma-separated args
  // narg = number of args, args = strings between commas

  char *args[MAXFUNCARG];
  int narg = parse_args(contents,args);

  Tree *newtree = NULL;
  double value1,value2;
  double values[MAXFUNCARG-2];

  if (tree) {
    newtree = new Tree();
    newtree->first = newtree->second = NULL;
    newtree->nextra = 0;
    Tree *argtree = NULL;
    evaluate(args[0],&argtree,ivar);
    newtree->first = argtree;
    if (narg > 1) {
      evaluate(args[1],&argtree,ivar);
      newtree->second = argtree;
      if (narg > 2) {
        newtree->nextra = narg-2;
        newtree->extra = new Tree*[narg-2];
        for (int i = 2; i < narg; i++) {
          evaluate(args[i],&argtree,ivar);
          newtree->extra[i-2] = argtree;
        }
      }
    }
    treestack[ntreestack++] = newtree;

  } else {
    value1 = evaluate(args[0],NULL,ivar);
    if (narg > 1) {
      value2 = evaluate(args[1],NULL,ivar);
      if (narg > 2) {
        for (int i = 2; i < narg; i++)
          values[i-2] = evaluate(args[i],NULL,ivar);
      }
    }
  }

  // individual math functions
  // customize by adding a function

  if (strcmp(word,"sqrt") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = SQRT;
    else {
      if (value1 < 0.0)
        print_var_error(FLERR,"Sqrt of negative value in "
                        "variable formula",ivar,0);
      argstack[nargstack++] = sqrt(value1);
    }

  } else if (strcmp(word,"exp") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = EXP;
    else argstack[nargstack++] = exp(value1);
  } else if (strcmp(word,"ln") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = LN;
    else {
      if (value1 <= 0.0)
        print_var_error(FLERR,"Log of zero/negative value in "
                        "variable formula",ivar,0);
      argstack[nargstack++] = log(value1);
    }
  } else if (strcmp(word,"log") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = LOG;
    else {
      if (value1 <= 0.0)
        print_var_error(FLERR,"Log of zero/negative value in "
                        "variable formula",ivar,0);
      argstack[nargstack++] = log10(value1);
    }
  } else if (strcmp(word,"abs") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = ABS;
    else argstack[nargstack++] = fabs(value1);

  } else if (strcmp(word,"sin") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = SIN;
    else argstack[nargstack++] = sin(value1);
  } else if (strcmp(word,"cos") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = COS;
    else argstack[nargstack++] = cos(value1);
  } else if (strcmp(word,"tan") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = TAN;
    else argstack[nargstack++] = tan(value1);

  } else if (strcmp(word,"asin") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = ASIN;
    else {
      if (value1 < -1.0 || value1 > 1.0)
        print_var_error(FLERR,"Arcsin of invalid value in variable formula",ivar,0);
      argstack[nargstack++] = asin(value1);
    }
  } else if (strcmp(word,"acos") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = ACOS;
    else {
      if (value1 < -1.0 || value1 > 1.0)
        print_var_error(FLERR,"Arccos of invalid value in variable formula",ivar,0);
      argstack[nargstack++] = acos(value1);
    }
  } else if (strcmp(word,"atan") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = ATAN;
    else argstack[nargstack++] = atan(value1);
  } else if (strcmp(word,"atan2") == 0) {
    if (narg != 2)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = ATAN2;
    else argstack[nargstack++] = atan2(value1,value2);

  } else if (strcmp(word,"random") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = RANDOM;
    else {
      if (randomequal == NULL) {
        int seed = static_cast<int> (values[0]);
        if (seed <= 0)
          print_var_error(FLERR,"Invalid math function in variable formula",ivar);
        randomequal = new RanMars(lmp,seed);
      }
      argstack[nargstack++] = randomequal->uniform()*(value2-value1) + value1;
    }
  } else if (strcmp(word,"normal") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = NORMAL;
    else {
      if (value2 < 0.0)
        print_var_error(FLERR,"Invalid math function in variable formula",ivar);
      if (randomequal == NULL) {
        int seed = static_cast<int> (values[0]);
        if (seed <= 0)
          print_var_error(FLERR,"Invalid math function in variable formula",ivar);
        randomequal = new RanMars(lmp,seed);
      }
      argstack[nargstack++] = value1 + value2*randomequal->gaussian();
    }

  } else if (strcmp(word,"ceil") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = CEIL;
    else argstack[nargstack++] = ceil(value1);

  } else if (strcmp(word,"floor") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = FLOOR;
    else argstack[nargstack++] = floor(value1);

  } else if (strcmp(word,"round") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = ROUND;
    else argstack[nargstack++] = MYROUND(value1);

  } else if (strcmp(word,"ramp") == 0) {
    if (narg != 2)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (update->whichflag == 0)
      print_var_error(FLERR,"Cannot use ramp in "
                      "variable formula between runs",ivar);
    if (tree) newtree->type = RAMP;
    else {
      double delta = update->ntimestep - update->beginstep;
      if (delta != 0.0) delta /= update->endstep - update->beginstep;
      double value = value1 + delta*(value2-value1);
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"stagger") == 0) {
    if (narg != 2)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = STAGGER;
    else {
      int ivalue1 = static_cast<int> (value1);
      int ivalue2 = static_cast<int> (value2);
      if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue1 <= ivalue2)
        print_var_error(FLERR,"Invalid math function in variable formula",ivar);
      int lower = update->ntimestep/ivalue1 * ivalue1;
      int delta = update->ntimestep - lower;
      double value;
      if (delta < ivalue2) value = lower+ivalue2;
      else value = lower+ivalue1;
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"logfreq") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = LOGFREQ;
    else {
      int ivalue1 = static_cast<int> (value1);
      int ivalue2 = static_cast<int> (value2);
      int ivalue3 = static_cast<int> (values[0]);
      if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 || ivalue2 >= ivalue3)
        print_var_error(FLERR,"Invalid math function in variable formula",ivar);
      double value;
      if (update->ntimestep < ivalue1) value = ivalue1;
      else {
        int lower = ivalue1;
        while (update->ntimestep >= ivalue3*lower) lower *= ivalue3;
        int multiple = update->ntimestep/lower;
        if (multiple < ivalue2) value = (multiple+1)*lower;
        else value = lower*ivalue3;
      }
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"logfreq2") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = LOGFREQ2;
    else {
      int ivalue1 = static_cast<int> (value1);
      int ivalue2 = static_cast<int> (value2);
      int ivalue3 = static_cast<int> (values[0]);
      if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 )
        print_var_error(FLERR,"Invalid math function in variable formula",ivar);
      double value;
      if (update->ntimestep < ivalue1) value = ivalue1;
      else {
        value = ivalue1;
        double delta = ivalue1*(ivalue3-1.0)/ivalue2;
        int count = 0;
        while (update->ntimestep >= value) {
          value += delta;
          count++;
          if (count % ivalue2 == 0) delta *= ivalue3;
        }
      }
      argstack[nargstack++] = ceil(value);
    }

  } else if (strcmp(word,"logfreq3") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = LOGFREQ3;
    else {
      int ivalue1 = static_cast<int> (value1);
      int ivalue2 = static_cast<int> (value2);
      int ivalue3 = static_cast<int> (values[0]);
      if (ivalue1 <= 0 || ivalue2 <= 1 || ivalue3 <= 0 ||
          ivalue3-ivalue1+1 < ivalue2 )
        print_var_error(FLERR,"Invalid math function in variable formula",ivar);
      double value;
      if (update->ntimestep < ivalue1) value = ivalue1;
      //else if (update->ntimestep <= ivalue3){
      else {
        value = ivalue1;
        double logsp = ivalue1;
        double factor = pow(((double)ivalue3)/ivalue1, 1.0/(ivalue2-1));
        int linsp = ivalue1;
        while (update->ntimestep >= value) {
          logsp *= factor;
          linsp++;
          if (linsp > logsp) value = linsp;
          else value = ceil(logsp)-(((int)ceil(logsp)-1)/ivalue3);
        }
      }
      if (update->ntimestep > ivalue3)
        error->all(FLERR,"Calls to variable exceeded limit");
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"stride") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = STRIDE;
    else {
      int ivalue1 = static_cast<int> (value1);
      int ivalue2 = static_cast<int> (value2);
      int ivalue3 = static_cast<int> (values[0]);
      if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
        error->one(FLERR,"Invalid math function in variable formula");
      double value;
      if (update->ntimestep < ivalue1) value = ivalue1;
      else if (update->ntimestep < ivalue2) {
        int offset = update->ntimestep - ivalue1;
        value = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
        if (value > ivalue2) value = MAXBIGINT;
      } else value = MAXBIGINT;
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"stride2") == 0) {
    if (narg != 6)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = STRIDE2;
    else {
      int ivalue1 = static_cast<int> (value1);
      int ivalue2 = static_cast<int> (value2);
      int ivalue3 = static_cast<int> (values[0]);
      int ivalue4 = static_cast<int> (values[1]);
      int ivalue5 = static_cast<int> (values[2]);
      int ivalue6 = static_cast<int> (values[3]);
      if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
        error->one(FLERR,"Invalid math function in variable formula");
      if (ivalue4 < 0 || ivalue5 < 0 || ivalue6 <= 0 || ivalue4 > ivalue5)
        error->one(FLERR,"Invalid math function in variable formula");
      if (ivalue4 < ivalue1 || ivalue5 > ivalue2)
        error->one(FLERR,"Invalid math function in variable formula");
      bigint istep;
      if (update->ntimestep < ivalue1) istep = ivalue1;
      else if (update->ntimestep < ivalue2) {
        if (update->ntimestep < ivalue4 || update->ntimestep > ivalue5) {
          int offset = update->ntimestep - ivalue1;
          istep = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
          if (update->ntimestep < ivalue4 && istep > ivalue4) istep = ivalue4;
        } else {
          int offset = update->ntimestep - ivalue4;
          istep = ivalue4 + (offset/ivalue6)*ivalue6 + ivalue6;
          if (istep > ivalue5) {
            int offset = ivalue5 - ivalue1;
            istep = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
            if (istep > ivalue2) istep = MAXBIGINT;
          }
        }
      } else istep = MAXBIGINT;
      double value = istep;
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"vdisplace") == 0) {
    if (narg != 2)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (update->whichflag == 0)
      print_var_error(FLERR,"Cannot use vdisplace in "
                      "variable formula between runs",ivar);
    if (tree) newtree->type = VDISPLACE;
    else {
      double delta = update->ntimestep - update->beginstep;
      double value = value1 + value2*delta*update->dt;
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"swiggle") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (update->whichflag == 0)
      print_var_error(FLERR,"Cannot use swiggle in "
                      "variable formula between runs",ivar);
    if (tree) newtree->type = CWIGGLE;
    else {
      if (values[0] == 0.0)
        print_var_error(FLERR,"Invalid math function in variable formula",ivar);
      double delta = update->ntimestep - update->beginstep;
      double omega = 2.0*MY_PI/values[0];
      double value = value1 + value2*sin(omega*delta*update->dt);
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"cwiggle") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (update->whichflag == 0)
      print_var_error(FLERR,"Cannot use cwiggle in "
                      "variable formula between runs",ivar);
    if (tree) newtree->type = CWIGGLE;
    else {
      if (values[0] == 0.0)
        print_var_error(FLERR,"Invalid math function in variable formula",ivar);
      double delta = update->ntimestep - update->beginstep;
      double omega = 2.0*MY_PI/values[0];
      double value = value1 + value2*(1.0-cos(omega*delta*update->dt));
      argstack[nargstack++] = value;
    }
  }

  // delete stored args

  for (int i = 0; i < narg; i++) delete [] args[i];

  return 1;
}

/* ----------------------------------------------------------------------
   process a group function in formula with optional region arg
   push result onto tree or arg stack
   word = group function
   contents = str between parentheses with one,two,three args
   return 0 if not a match, 1 if successfully processed
   customize by adding a group function with optional region arg:
     count(group),mass(group),charge(group),
     xcm(group,dim),vcm(group,dim),fcm(group,dim),
     bound(group,xmin),gyration(group),ke(group),angmom(group,dim),
     torque(group,dim),inertia(group,dim),omega(group,dim)
------------------------------------------------------------------------- */

int Variable::group_function(char *word, char *contents, Tree **tree,
                             Tree **treestack, int &ntreestack,
                             double *argstack, int &nargstack, int ivar)
{
  // word not a match to any group function

  if (strcmp(word,"count") && strcmp(word,"mass") &&
      strcmp(word,"charge") && strcmp(word,"xcm") &&
      strcmp(word,"vcm") && strcmp(word,"fcm") &&
      strcmp(word,"bound") && strcmp(word,"gyration") &&
      strcmp(word,"ke") && strcmp(word,"angmom") &&
      strcmp(word,"torque") && strcmp(word,"inertia") &&
      strcmp(word,"omega"))
    return 0;

  // parse contents for comma-separated args
  // narg = number of args, args = strings between commas

  char *args[MAXFUNCARG];
  int narg = parse_args(contents,args);

  // group to operate on

  int igroup = group->find(args[0]);
  if (igroup == -1)
    print_var_error(FLERR,"Group ID in variable formula does not exist",ivar);

  // match word to group function

  double value;

  if (strcmp(word,"count") == 0) {
    if (narg == 1) value = group->count(igroup);
    else if (narg == 2)
      value = group->count(igroup,region_function(args[1],ivar));
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);

  } else if (strcmp(word,"mass") == 0) {
    if (narg == 1) value = group->mass(igroup);
    else if (narg == 2) value = group->mass(igroup,region_function(args[1],ivar));
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);

  } else if (strcmp(word,"charge") == 0) {
    if (narg == 1) value = group->charge(igroup);
    else if (narg == 2)
      value = group->charge(igroup,region_function(args[1],ivar));
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);

  } else if (strcmp(word,"xcm") == 0) {
    atom->check_mass(FLERR);
    double xcm[3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
    } else if (narg == 3) {
      int iregion = region_function(args[2],ivar);
      double masstotal = group->mass(igroup,iregion);
      group->xcm(igroup,masstotal,xcm,iregion);
    } else print_var_error(FLERR,"Invalid group function in "
                           "variable formula",ivar);
    if (strcmp(args[1],"x") == 0) value = xcm[0];
    else if (strcmp(args[1],"y") == 0) value = xcm[1];
    else if (strcmp(args[1],"z") == 0) value = xcm[2];
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);

  } else if (strcmp(word,"vcm") == 0) {
    atom->check_mass(FLERR);
    double vcm[3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->vcm(igroup,masstotal,vcm);
    } else if (narg == 3) {
      int iregion = region_function(args[2],ivar);
      double masstotal = group->mass(igroup,iregion);
      group->vcm(igroup,masstotal,vcm,iregion);
    } else print_var_error(FLERR,"Invalid group function in "
                           "variable formula",ivar);
    if (strcmp(args[1],"x") == 0) value = vcm[0];
    else if (strcmp(args[1],"y") == 0) value = vcm[1];
    else if (strcmp(args[1],"z") == 0) value = vcm[2];
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);

  } else if (strcmp(word,"fcm") == 0) {
    double fcm[3];
    if (narg == 2) group->fcm(igroup,fcm);
    else if (narg == 3) group->fcm(igroup,fcm,region_function(args[2],ivar));
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);
    if (strcmp(args[1],"x") == 0) value = fcm[0];
    else if (strcmp(args[1],"y") == 0) value = fcm[1];
    else if (strcmp(args[1],"z") == 0) value = fcm[2];
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);

  } else if (strcmp(word,"bound") == 0) {
    double minmax[6];
    if (narg == 2) group->bounds(igroup,minmax);
    else if (narg == 3)
      group->bounds(igroup,minmax,region_function(args[2],ivar));
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);
    if (strcmp(args[1],"xmin") == 0) value = minmax[0];
    else if (strcmp(args[1],"xmax") == 0) value = minmax[1];
    else if (strcmp(args[1],"ymin") == 0) value = minmax[2];
    else if (strcmp(args[1],"ymax") == 0) value = minmax[3];
    else if (strcmp(args[1],"zmin") == 0) value = minmax[4];
    else if (strcmp(args[1],"zmax") == 0) value = minmax[5];
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);

  } else if (strcmp(word,"gyration") == 0) {
    atom->check_mass(FLERR);
    double xcm[3];
    if (narg == 1) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
      value = group->gyration(igroup,masstotal,xcm);
    } else if (narg == 2) {
      int iregion = region_function(args[1],ivar);
      double masstotal = group->mass(igroup,iregion);
      group->xcm(igroup,masstotal,xcm,iregion);
      value = group->gyration(igroup,masstotal,xcm,iregion);
    } else print_var_error(FLERR,"Invalid group function in "
                           "variable formula",ivar);

  } else if (strcmp(word,"ke") == 0) {
    if (narg == 1) value = group->ke(igroup);
    else if (narg == 2) value = group->ke(igroup,region_function(args[1],ivar));
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);

  } else if (strcmp(word,"angmom") == 0) {
    atom->check_mass(FLERR);
    double xcm[3],lmom[3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
      group->angmom(igroup,xcm,lmom);
    } else if (narg == 3) {
      int iregion = region_function(args[2],ivar);
      double masstotal = group->mass(igroup,iregion);
      group->xcm(igroup,masstotal,xcm,iregion);
      group->angmom(igroup,xcm,lmom,iregion);
    } else print_var_error(FLERR,"Invalid group function in "
                           "variable formula",ivar);
    if (strcmp(args[1],"x") == 0) value = lmom[0];
    else if (strcmp(args[1],"y") == 0) value = lmom[1];
    else if (strcmp(args[1],"z") == 0) value = lmom[2];
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);

  } else if (strcmp(word,"torque") == 0) {
    atom->check_mass(FLERR);
    double xcm[3],tq[3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
      group->torque(igroup,xcm,tq);
    } else if (narg == 3) {
      int iregion = region_function(args[2],ivar);
      double masstotal = group->mass(igroup,iregion);
      group->xcm(igroup,masstotal,xcm,iregion);
      group->torque(igroup,xcm,tq,iregion);
    } else print_var_error(FLERR,"Invalid group function in "
                           "variable formula",ivar);
    if (strcmp(args[1],"x") == 0) value = tq[0];
    else if (strcmp(args[1],"y") == 0) value = tq[1];
    else if (strcmp(args[1],"z") == 0) value = tq[2];
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);

  } else if (strcmp(word,"inertia") == 0) {
    atom->check_mass(FLERR);
    double xcm[3],inertia[3][3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
      group->inertia(igroup,xcm,inertia);
    } else if (narg == 3) {
      int iregion = region_function(args[2],ivar);
      double masstotal = group->mass(igroup,iregion);
      group->xcm(igroup,masstotal,xcm,iregion);
      group->inertia(igroup,xcm,inertia,iregion);
    } else print_var_error(FLERR,"Invalid group function in "
                           "variable formula",ivar);
    if (strcmp(args[1],"xx") == 0) value = inertia[0][0];
    else if (strcmp(args[1],"yy") == 0) value = inertia[1][1];
    else if (strcmp(args[1],"zz") == 0) value = inertia[2][2];
    else if (strcmp(args[1],"xy") == 0) value = inertia[0][1];
    else if (strcmp(args[1],"yz") == 0) value = inertia[1][2];
    else if (strcmp(args[1],"xz") == 0) value = inertia[0][2];
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);

  } else if (strcmp(word,"omega") == 0) {
    atom->check_mass(FLERR);
    double xcm[3],angmom[3],inertia[3][3],omega[3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
      group->angmom(igroup,xcm,angmom);
      group->inertia(igroup,xcm,inertia);
      group->omega(angmom,inertia,omega);
    } else if (narg == 3) {
      int iregion = region_function(args[2],ivar);
      double masstotal = group->mass(igroup,iregion);
      group->xcm(igroup,masstotal,xcm,iregion);
      group->angmom(igroup,xcm,angmom,iregion);
      group->inertia(igroup,xcm,inertia,iregion);
      group->omega(angmom,inertia,omega);
    } else print_var_error(FLERR,"Invalid group function in "
                           "variable formula",ivar);
    if (strcmp(args[1],"x") == 0) value = omega[0];
    else if (strcmp(args[1],"y") == 0) value = omega[1];
    else if (strcmp(args[1],"z") == 0) value = omega[2];
    else print_var_error(FLERR,"Invalid group function in variable formula",ivar);
  }

  // delete stored args

  for (int i = 0; i < narg; i++) delete [] args[i];

  // save value in tree or on argstack

  if (tree) {
    Tree *newtree = new Tree();
    newtree->type = VALUE;
    newtree->value = value;
    newtree->first = newtree->second = NULL;
    newtree->nextra = 0;
    treestack[ntreestack++] = newtree;
  } else argstack[nargstack++] = value;

  return 1;
}

/* ---------------------------------------------------------------------- */

int Variable::region_function(char *id, int ivar)
{
  int iregion = domain->find_region(id);
  if (iregion == -1)
    print_var_error(FLERR,"Region ID in variable formula does not exist",ivar);

  // init region in case sub-regions have been deleted

  domain->regions[iregion]->init();

  return iregion;
}

/* ----------------------------------------------------------------------
   process a special function in formula
   push result onto tree or arg stack
   word = special function
   contents = str between parentheses with one,two,three args
   return 0 if not a match, 1 if successfully processed
   customize by adding a special function:
     sum(x),min(x),max(x),ave(x),trap(x),slope(x),
     gmask(x),rmask(x),grmask(x,y),next(x)
------------------------------------------------------------------------- */

int Variable::special_function(char *word, char *contents, Tree **tree,
                               Tree **treestack, int &ntreestack,
                               double *argstack, int &nargstack, int ivar)
{
  bigint sx,sxx;
  double value,sy,sxy;

  // word not a match to any special function

  if (strcmp(word,"sum") && strcmp(word,"min") && strcmp(word,"max") &&
      strcmp(word,"ave") && strcmp(word,"trap") && strcmp(word,"slope") &&
      strcmp(word,"gmask") && strcmp(word,"rmask") &&
      strcmp(word,"grmask") && strcmp(word,"next") &&
      strcmp(word,"is_active") && strcmp(word,"is_defined") &&
      strcmp(word,"is_available"))
    return 0;

  // parse contents for comma-separated args
  // narg = number of args, args = strings between commas

  char *args[MAXFUNCARG];
  int narg = parse_args(contents,args);

  // special functions that operate on global vectors

  if (strcmp(word,"sum") == 0 || strcmp(word,"min") == 0 ||
      strcmp(word,"max") == 0 || strcmp(word,"ave") == 0 ||
      strcmp(word,"trap") == 0 || strcmp(word,"slope") == 0) {

    int method = 0;
    if (strcmp(word,"sum") == 0) method = SUM;
    else if (strcmp(word,"min") == 0) method = XMIN;
    else if (strcmp(word,"max") == 0) method = XMAX;
    else if (strcmp(word,"ave") == 0) method = AVE;
    else if (strcmp(word,"trap") == 0) method = TRAP;
    else if (strcmp(word,"slope") == 0) method = SLOPE;

    if (narg != 1)
      print_var_error(FLERR,"Invalid special function in variable formula",ivar);

    Compute *compute = NULL;
    Fix *fix = NULL;
    int ivar = -1;
    int index,nvec,nstride;
    char *ptr1,*ptr2;

    // argument is compute

    if (strstr(args[0],"c_") == args[0]) {
      ptr1 = strchr(args[0],'[');
      if (ptr1) {
        ptr2 = ptr1;
        index = (int) int_between_brackets(ptr2,0);
        *ptr1 = '\0';
      } else index = 0;

      int icompute = modify->find_compute(&args[0][2]);
      if (icompute < 0) {
        char msg[128];
        snprintf(msg,128,"Invalid compute ID '%s' in variable formula",word+2);
        print_var_error(FLERR,msg,ivar);
      }
      compute = modify->compute[icompute];
      if (index == 0 && compute->vector_flag) {
        if (update->whichflag == 0) {
          if (compute->invoked_vector != update->ntimestep)
            print_var_error(FLERR,"Compute used in variable between runs "
                            "is not current",ivar);
        } else if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        nvec = compute->size_vector;
        nstride = 1;
      } else if (index && compute->array_flag) {
        if (index > compute->size_array_cols)
          print_var_error(FLERR,"Variable formula compute array "
                          "is accessed out-of-range",ivar,0);
        if (update->whichflag == 0) {
          if (compute->invoked_array != update->ntimestep)
            print_var_error(FLERR,"Compute used in variable between runs "
                            "is not current",ivar);
        } else if (!(compute->invoked_flag & INVOKED_ARRAY)) {
          compute->compute_array();
          compute->invoked_flag |= INVOKED_ARRAY;
        }
        nvec = compute->size_array_rows;
        nstride = compute->size_array_cols;
      } else print_var_error(FLERR,"Mismatched compute in variable formula",ivar);

    // argument is fix

    } else if (strstr(args[0],"f_") == args[0]) {
      ptr1 = strchr(args[0],'[');
      if (ptr1) {
        ptr2 = ptr1;
        index = (int) int_between_brackets(ptr2,0);
        *ptr1 = '\0';
      } else index = 0;

      int ifix = modify->find_fix(&args[0][2]);
      if (ifix < 0)
        print_var_error(FLERR,"Invalid fix ID in variable formula",ivar);
      fix = modify->fix[ifix];
      if (index == 0 && fix->vector_flag) {
        if (update->whichflag > 0 && update->ntimestep % fix->global_freq)
          print_var_error(FLERR,"Fix in variable not computed at "
                          "compatible time",ivar);
        nvec = fix->size_vector;
        nstride = 1;
      } else if (index && fix->array_flag) {
        if (index > fix->size_array_cols)
          print_var_error(FLERR,"Variable formula fix array "
                          "is accessed out-of-range",ivar);
        if (update->whichflag > 0 && update->ntimestep % fix->global_freq)
          print_var_error(FLERR,"Fix in variable not computed at "
                          "compatible time",ivar);
        nvec = fix->size_array_rows;
        nstride = fix->size_array_cols;
      } else print_var_error(FLERR,"Mismatched fix in variable formula",ivar);

    // argument is vector-style variable

    } else if (strstr(args[0],"v_") == args[0]) {
      ptr1 = strchr(args[0],'[');
      if (ptr1) {
        ptr2 = ptr1;
        index = (int) int_between_brackets(ptr2,0);
        *ptr1 = '\0';
      } else index = 0;

      if (index)
        print_var_error(FLERR,"Invalid special function in "
                        "variable formula",ivar);
      ivar = find(&args[0][2]);
      if (ivar < 0)
        print_var_error(FLERR,"Invalid special function in "
                        "variable formula",ivar);
      if (style[ivar] != VECTOR)
        print_var_error(FLERR,"Mis-matched special function variable "
                        "in variable formula",ivar);
      if (eval_in_progress[ivar])
        print_var_error(FLERR,"Variable has circular dependency",ivar);

      double *vec;
      nvec = compute_vector(ivar,&vec);
      nstride = 1;

    } else print_var_error(FLERR,"Invalid special function in "
                           "variable formula",ivar);

    value = 0.0;
    if (method == SLOPE) {
      sx = sxx = 0;
      sy = sxy = 0.0;
    }
    else if (method == XMIN) value = BIG;
    else if (method == XMAX) value = -BIG;

    if (compute) {
      double *vec;
      if (index) {
        if (compute->array) vec = &compute->array[0][index-1];
        else vec = NULL;
      } else vec = compute->vector;

      int j = 0;
      for (int i = 0; i < nvec; i++) {
        if (method == SUM) value += vec[j];
        else if (method == XMIN) value = MIN(value,vec[j]);
        else if (method == XMAX) value = MAX(value,vec[j]);
        else if (method == AVE) value += vec[j];
        else if (method == TRAP) value += vec[j];
        else if (method == SLOPE) {
          sx += i;
          sy += vec[j];
          sxx += i*i;
          sxy += i*vec[j];
        }
        j += nstride;
      }
      if (method == TRAP) value -= 0.5*vec[0] + 0.5*vec[nvec-1];
    }

    if (fix) {
      double one;
      for (int i = 0; i < nvec; i++) {
        if (index) one = fix->compute_array(i,index-1);
        else one = fix->compute_vector(i);
        if (method == SUM) value += one;
        else if (method == XMIN) value = MIN(value,one);
        else if (method == XMAX) value = MAX(value,one);
        else if (method == AVE) value += one;
        else if (method == TRAP) value += one;
        else if (method == SLOPE) {
          sx += i;
          sy += one;
          sxx += i*i;
          sxy += i*one;
        }
      }
      if (method == TRAP) {
        if (index) value -= 0.5*fix->compute_array(0,index-1) +
                     0.5*fix->compute_array(nvec-1,index-1);
        else value -= 0.5*fix->compute_vector(0) +
               0.5*fix->compute_vector(nvec-1);
      }
    }

    if (ivar >= 0) {
      double one;
      double *vec = vecs[ivar].values;
      for (int i = 0; i < nvec; i++) {
        one = vec[i];
        if (method == SUM) value += one;
        else if (method == XMIN) value = MIN(value,one);
        else if (method == XMAX) value = MAX(value,one);
        else if (method == AVE) value += one;
        else if (method == TRAP) value += one;
        else if (method == SLOPE) {
          sx += i;
          sy += one;
          sxx += i*i;
          sxy += i*one;
        }
      }
      if (method == TRAP) value -= 0.5*vec[0] + 0.5*vec[nvec-1];
    }

    if (method == AVE) value /= nvec;

    if (method == SLOPE) {
      double numerator = nvec*sxy - sx*sy;
      double denominator = nvec*sxx - sx*sx;
      if (denominator != 0.0) value = numerator/denominator;
      else value = BIG;
    }

    // save value in tree or on argstack

    if (tree) {
      Tree *newtree = new Tree();
      newtree->type = VALUE;
      newtree->value = value;
      newtree->first = newtree->second = NULL;
      newtree->nextra = 0;
      treestack[ntreestack++] = newtree;
    } else argstack[nargstack++] = value;

  // mask special functions

  } else if (strcmp(word,"gmask") == 0) {
    if (tree == NULL)
      print_var_error(FLERR,"Gmask function in equal-style "
                      "variable formula",ivar);
    if (narg != 1)
      print_var_error(FLERR,"Invalid special function in variable formula",ivar);

    int igroup = group->find(args[0]);
    if (igroup == -1)
      print_var_error(FLERR,"Group ID in variable formula does not exist",ivar);

    Tree *newtree = new Tree();
    newtree->type = GMASK;
    newtree->ivalue1 = group->bitmask[igroup];
    newtree->first = newtree->second = NULL;
    newtree->nextra = 0;
    treestack[ntreestack++] = newtree;

  } else if (strcmp(word,"rmask") == 0) {
    if (tree == NULL)
      print_var_error(FLERR,"Rmask function in equal-style "
                      "variable formula",ivar);
    if (narg != 1)
      print_var_error(FLERR,"Invalid special function in variable formula",ivar);

    int iregion = region_function(args[0],ivar);
    domain->regions[iregion]->prematch();

    Tree *newtree = new Tree();
    newtree->type = RMASK;
    newtree->ivalue1 = iregion;
    newtree->first = newtree->second = NULL;
    newtree->nextra = 0;
    treestack[ntreestack++] = newtree;

  } else if (strcmp(word,"grmask") == 0) {
    if (tree == NULL)
      print_var_error(FLERR,"Grmask function in equal-style "
                      "variable formula",ivar);
    if (narg != 2)
      print_var_error(FLERR,"Invalid special function in variable formula",ivar);

    int igroup = group->find(args[0]);
    if (igroup == -1)
      print_var_error(FLERR,"Group ID in variable formula does not exist",ivar);
    int iregion = region_function(args[1],ivar);
    domain->regions[iregion]->prematch();

    Tree *newtree = new Tree();
    newtree->type = GRMASK;
    newtree->ivalue1 = group->bitmask[igroup];
    newtree->ivalue2 = iregion;
    newtree->first = newtree->second = NULL;
    newtree->nextra = 0;
    treestack[ntreestack++] = newtree;

  // special function for file-style or atomfile-style variables

  } else if (strcmp(word,"next") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid special function in variable formula",ivar);

    int ivar = find(args[0]);
    if (ivar < 0)
      print_var_error(FLERR,"Variable ID in "
                      "variable formula does not exist",ivar);

    // SCALARFILE has single current value, read next one
    // save value in tree or on argstack

    if (style[ivar] == SCALARFILE) {
      double value = atof(data[ivar][0]);
      int done = reader[ivar]->read_scalar(data[ivar][0]);
      if (done) remove(ivar);

      if (tree) {
        Tree *newtree = new Tree();
        newtree->type = VALUE;
        newtree->value = value;
        newtree->first = newtree->second = NULL;
        newtree->nextra = 0;
        treestack[ntreestack++] = newtree;
      } else argstack[nargstack++] = value;

    // ATOMFILE has per-atom values, save values in tree
    // copy current per-atom values into result so can read next ones
    // set selfalloc = 1 so result will be deleted by free_tree() after eval

    } else if (style[ivar] == ATOMFILE) {
      if (tree == NULL)
        print_var_error(FLERR,"Atomfile variable in "
                        "equal-style variable formula",ivar);

      double *result;
      memory->create(result,atom->nlocal,"variable:result");
      memcpy(result,reader[ivar]->fixstore->vstore,atom->nlocal*sizeof(double));

      int done = reader[ivar]->read_peratom();
      if (done) remove(ivar);

      Tree *newtree = new Tree();
      newtree->type = ATOMARRAY;
      newtree->array = result;
      newtree->nstride = 1;
      newtree->selfalloc = 1;
      newtree->first = newtree->second = NULL;
      newtree->nextra = 0;
      treestack[ntreestack++] = newtree;

    } else print_var_error(FLERR,"Invalid variable style in "
                           "special function next",ivar);

  } else if (strcmp(word,"is_active") == 0) {
    if (narg != 2)
      print_var_error(FLERR,"Invalid is_active() function in "
                      "variable formula",ivar);

    Info info(lmp);
    value = (info.is_active(args[0],args[1])) ? 1.0 : 0.0;

    // save value in tree or on argstack

    if (tree) {
      Tree *newtree = new Tree();
      newtree->type = VALUE;
      newtree->value = value;
      newtree->first = newtree->second = NULL;
      newtree->nextra = 0;
      treestack[ntreestack++] = newtree;
    } else argstack[nargstack++] = value;

  } else if (strcmp(word,"is_available") == 0) {
    if (narg != 2)
      print_var_error(FLERR,"Invalid is_available() function in "
                      "variable formula",ivar);

    Info info(lmp);
    value = (info.is_available(args[0],args[1])) ? 1.0 : 0.0;

    // save value in tree or on argstack

    if (tree) {
      Tree *newtree = new Tree();
      newtree->type = VALUE;
      newtree->value = value;
      newtree->first = newtree->second = NULL;
      newtree->nextra = 0;
      treestack[ntreestack++] = newtree;
    } else argstack[nargstack++] = value;

  } else if (strcmp(word,"is_defined") == 0) {
    if (narg != 2)
      print_var_error(FLERR,"Invalid is_defined() function in "
                      "variable formula",ivar);

    Info info(lmp);
    value = (info.is_defined(args[0],args[1])) ? 1.0 : 0.0;

    // save value in tree or on argstack

    if (tree) {
      Tree *newtree = new Tree();
      newtree->type = VALUE;
      newtree->value = value;
      newtree->first = newtree->second = NULL;
      newtree->nextra = 0;
      treestack[ntreestack++] = newtree;
    } else argstack[nargstack++] = value;
  }

  // delete stored args

  for (int i = 0; i < narg; i++) delete [] args[i];

  return 1;
}

/* ----------------------------------------------------------------------
   extract a global value from a per-atom quantity in a formula
   flag = 0 -> word is an atom vector
   flag = 1 -> vector is a per-atom compute or fix quantity with nstride
   id = global ID of atom, converted to local index
   push result onto tree or arg stack
   customize by adding an atom vector:
     id,mass,type,mol,x,y,z,vx,vy,vz,fx,fy,fz,q
------------------------------------------------------------------------- */

void Variable::peratom2global(int flag, char *word,
                              double *vector, int nstride, tagint id,
                              Tree **tree, Tree **treestack, int &ntreestack,
                              double *argstack, int &nargstack)
{
  // error check for ID larger than any atom
  // int_between_brackets() already checked for ID <= 0

  if (atom->map_style == 0)
    error->all(FLERR,
               "Indexed per-atom vector in variable formula without atom map");

  if (id > atom->map_tag_max)
    error->all(FLERR,"Variable atom ID is too large");

  // if ID does not exist, index will be -1 for all procs,
  // and mine will be set to 0.0

  int index = atom->map(id);

  double mine;
  if (index >= 0 && index < atom->nlocal) {

    if (flag == 0) {
      if (strcmp(word,"id") == 0) mine = atom->tag[index];
      else if (strcmp(word,"mass") == 0) {
        if (atom->rmass) mine = atom->rmass[index];
        else mine = atom->mass[atom->type[index]];
      }
      else if (strcmp(word,"type") == 0) mine = atom->type[index];
      else if (strcmp(word,"mol") == 0) {
        if (!atom->molecule_flag)
          error->one(FLERR,"Variable uses atom property that isn't allocated");
        mine = atom->molecule[index];
      }
      else if (strcmp(word,"x") == 0) mine = atom->x[index][0];
      else if (strcmp(word,"y") == 0) mine = atom->x[index][1];
      else if (strcmp(word,"z") == 0) mine = atom->x[index][2];
      else if (strcmp(word,"vx") == 0) mine = atom->v[index][0];
      else if (strcmp(word,"vy") == 0) mine = atom->v[index][1];
      else if (strcmp(word,"vz") == 0) mine = atom->v[index][2];
      else if (strcmp(word,"fx") == 0) mine = atom->f[index][0];
      else if (strcmp(word,"fy") == 0) mine = atom->f[index][1];
      else if (strcmp(word,"fz") == 0) mine = atom->f[index][2];
      else if (strcmp(word,"q") == 0) {
        if (!atom->q_flag)
          error->one(FLERR,"Variable uses atom property that isn't allocated");
        mine = atom->q[index];
      }
      else error->one(FLERR,"Invalid atom vector in variable formula");

    } else mine = vector[index*nstride];

  } else mine = 0.0;

  double value;
  MPI_Allreduce(&mine,&value,1,MPI_DOUBLE,MPI_SUM,world);

  if (tree) {
    Tree *newtree = new Tree();
    newtree->type = VALUE;
    newtree->value = value;
    newtree->first = newtree->second = NULL;
    newtree->nextra = 0;
    treestack[ntreestack++] = newtree;
  } else argstack[nargstack++] = value;
}

/* ----------------------------------------------------------------------
   check if word matches an atom vector
   return 1 if yes, else 0
   customize by adding an atom vector:
     id,mass,type,mol,x,y,z,vx,vy,vz,fx,fy,fz,q
------------------------------------------------------------------------- */

int Variable::is_atom_vector(char *word)
{
  if (strcmp(word,"id") == 0) return 1;
  if (strcmp(word,"mass") == 0) return 1;
  if (strcmp(word,"type") == 0) return 1;
  if (strcmp(word,"mol") == 0) return 1;
  if (strcmp(word,"x") == 0) return 1;
  if (strcmp(word,"y") == 0) return 1;
  if (strcmp(word,"z") == 0) return 1;
  if (strcmp(word,"vx") == 0) return 1;
  if (strcmp(word,"vy") == 0) return 1;
  if (strcmp(word,"vz") == 0) return 1;
  if (strcmp(word,"fx") == 0) return 1;
  if (strcmp(word,"fy") == 0) return 1;
  if (strcmp(word,"fz") == 0) return 1;
  if (strcmp(word,"q") == 0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   process an atom vector in formula
   push result onto tree
   word = atom vector
   customize by adding an atom vector:
     id,mass,type,mol,x,y,z,vx,vy,vz,fx,fy,fz,q
------------------------------------------------------------------------- */

void Variable::atom_vector(char *word, Tree **tree,
                           Tree **treestack, int &ntreestack)
{
  if (tree == NULL)
    error->all(FLERR,"Atom vector in equal-style variable formula");

  Tree *newtree = new Tree();
  newtree->type = ATOMARRAY;
  newtree->nstride = 3;
  newtree->selfalloc = 0;
  newtree->first = newtree->second = NULL;
  newtree->nextra = 0;
  treestack[ntreestack++] = newtree;

  if (strcmp(word,"id") == 0) {
    if (sizeof(tagint) == sizeof(smallint)) {
      newtree->type = INTARRAY;
      newtree->iarray = (int *) atom->tag;
    } else {
      newtree->type = BIGINTARRAY;
      newtree->barray = (bigint *) atom->tag;
    }
    newtree->nstride = 1;

  } else if (strcmp(word,"mass") == 0) {
    if (atom->rmass) {
      newtree->nstride = 1;
      newtree->array = atom->rmass;
    } else {
      newtree->type = TYPEARRAY;
      newtree->array = atom->mass;
    }

  } else if (strcmp(word,"type") == 0) {
    newtree->type = INTARRAY;
    newtree->nstride = 1;
    newtree->iarray = atom->type;

  } else if (strcmp(word,"mol") == 0) {
    if (!atom->molecule_flag)
      error->one(FLERR,"Variable uses atom property that isn't allocated");
    if (sizeof(tagint) == sizeof(smallint)) {
      newtree->type = INTARRAY;
      newtree->iarray = (int *) atom->molecule;
    } else {
      newtree->type = BIGINTARRAY;
      newtree->barray = (bigint *) atom->molecule;
    }
    newtree->nstride = 1;
  }

  else if (strcmp(word,"x") == 0) newtree->array = &atom->x[0][0];
  else if (strcmp(word,"y") == 0) newtree->array = &atom->x[0][1];
  else if (strcmp(word,"z") == 0) newtree->array = &atom->x[0][2];
  else if (strcmp(word,"vx") == 0) newtree->array = &atom->v[0][0];
  else if (strcmp(word,"vy") == 0) newtree->array = &atom->v[0][1];
  else if (strcmp(word,"vz") == 0) newtree->array = &atom->v[0][2];
  else if (strcmp(word,"fx") == 0) newtree->array = &atom->f[0][0];
  else if (strcmp(word,"fy") == 0) newtree->array = &atom->f[0][1];
  else if (strcmp(word,"fz") == 0) newtree->array = &atom->f[0][2];

  else if (strcmp(word,"q") == 0) {
    newtree->nstride = 1;
    newtree->array = atom->q;
  }
}

/* ----------------------------------------------------------------------
   check if word matches a constant
   return 1 if yes, else 0
   customize by adding a constant: PI, version
------------------------------------------------------------------------- */

int Variable::is_constant(char *word)
{
  if (strcmp(word,"PI") == 0) return 1;
  if (strcmp(word,"version") == 0) return 1;
  if (strcmp(word,"yes") == 0) return 1;
  if (strcmp(word,"no") == 0) return 1;
  if (strcmp(word,"on") == 0) return 1;
  if (strcmp(word,"off") == 0) return 1;
  if (strcmp(word,"true") == 0) return 1;
  if (strcmp(word,"false") == 0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   process a constant in formula
   customize by adding a constant: PI, version
------------------------------------------------------------------------- */

double Variable::constant(char *word)
{
  if (strcmp(word,"PI") == 0) return MY_PI;
  if (strcmp(word,"version") == 0) return atof(universe->num_ver);
  if (strcmp(word,"yes") == 0) return 1.0;
  if (strcmp(word,"no") == 0) return 0.0;
  if (strcmp(word,"on") == 0) return 1.0;
  if (strcmp(word,"off") == 0) return 0.0;
  if (strcmp(word,"true") == 0) return 1.0;
  if (strcmp(word,"false") == 0) return 0.0;
  return 0.0;
}

/* ----------------------------------------------------------------------
   parse string for comma-separated args
   store copy of each arg in args array
   max allowed # of args = MAXFUNCARG
------------------------------------------------------------------------- */

int Variable::parse_args(char *str, char **args)
{
  int n;
  char *ptrnext;

  int narg = 0;
  char *ptr = str;

  while (ptr && narg < MAXFUNCARG) {
    ptrnext = find_next_comma(ptr);
    if (ptrnext) *ptrnext = '\0';
    n = strlen(ptr) + 1;
    args[narg] = new char[n];
    strcpy(args[narg],ptr);
    narg++;
    ptr = ptrnext;
    if (ptr) ptr++;
  }

  if (ptr) error->all(FLERR,"Too many args in variable function");
  return narg;
}


/* ----------------------------------------------------------------------
   find next comma in str
   skip commas inside one or more nested parenthesis
   only return ptr to comma at level 0, else NULL if not found
------------------------------------------------------------------------- */

char *Variable::find_next_comma(char *str)
{
  int level = 0;
  for (char *p = str; *p; ++p) {
    if ('(' == *p) level++;
    else if (')' == *p) level--;
    else if (',' == *p && !level) return p;
  }
  return NULL;
}


/* ----------------------------------------------------------------------
   helper routine for printing variable name with error message
------------------------------------------------------------------------- */

void Variable::print_var_error(const char *srcfile, int lineno,
                               const char *errmsg, int ivar, int global)
{
  if ((ivar >= 0) && (ivar < nvar)) {
    char msg[128];

    snprintf(msg,128,"Variable %s: %s",names[ivar],errmsg);
    if (global)
      error->all(srcfile,lineno,msg);
    else
      error->one(srcfile,lineno,msg);
  } else {
    if (global)
      error->all(srcfile,lineno,errmsg);
    else
      error->one(srcfile,lineno,errmsg);
  }
}

/* ----------------------------------------------------------------------
   debug routine for printing formula tree recursively
------------------------------------------------------------------------- */

void Variable::print_tree(Tree *tree, int level)
{
  printf("TREE %d: %d %g\n",level,tree->type,tree->value);
  if (tree->first) print_tree(tree->first,level+1);
  if (tree->second) print_tree(tree->second,level+1);
  if (tree->nextra)
    for (int i = 0; i < tree->nextra; i++) print_tree(tree->extra[i],level+1);
  return;
}

/* ----------------------------------------------------------------------
   recursive evaluation of string str
   called from "if" command in input script
   str is a boolean expression containing one or more items:
     number = 0.0, -5.45, 2.8e-4, ...
     math operation = (),x==y,x!=y,x<y,x<=y,x>y,x>=y,x&&y,x||y
------------------------------------------------------------------------- */

double Variable::evaluate_boolean(char *str)
{
  int op,opprevious,flag1,flag2;
  double value1,value2;
  char onechar;
  char *str1,*str2;

  struct Arg {
    int flag;          // 0 for numeric value, 1 for string
    double value;      // stored numeric value
    char *str;         // stored string
  };

  Arg argstack[MAXLEVEL];
  int opstack[MAXLEVEL];
  int nargstack = 0;
  int nopstack = 0;

  int i = 0;
  int expect = ARG;

  while (1) {
    onechar = str[i];

    // whitespace: just skip

    if (isspace(onechar)) i++;

    // ----------------
    // parentheses: recursively evaluate contents of parens
    // ----------------

    else if (onechar == '(') {
      if (expect == OP)
        error->all(FLERR,"Invalid Boolean syntax in if command");
      expect = OP;

      char *contents;
      i = find_matching_paren(str,i,contents,-1);
      i++;

      // evaluate contents and push on stack

      argstack[nargstack].value = evaluate_boolean(contents);
      argstack[nargstack].flag = 0;
      nargstack++;

      delete [] contents;

    // ----------------
    // number: push value onto stack
    // ----------------

    } else if (isdigit(onechar) || onechar == '.' || onechar == '-') {
      if (expect == OP)
        error->all(FLERR,"Invalid Boolean syntax in if command");
      expect = OP;

      // set I to end of number, including scientific notation

      int istart = i++;
      while (isdigit(str[i]) || str[i] == '.') i++;
      if (str[i] == 'e' || str[i] == 'E') {
        i++;
        if (str[i] == '+' || str[i] == '-') i++;
        while (isdigit(str[i])) i++;
      }

      onechar = str[i];
      str[i] = '\0';
      argstack[nargstack].value = atof(&str[istart]);
      str[i] = onechar;

      argstack[nargstack++].flag = 0;

    // ----------------
    // string: push string onto stack
    // ----------------

    } else if (isalpha(onechar)) {
      if (expect == OP)
        error->all(FLERR,"Invalid Boolean syntax in if command");
      expect = OP;

      // set I to end of string

      int istart = i++;
      while (isalnum(str[i]) || str[i] == '_') i++;

      int n = i - istart + 1;
      argstack[nargstack].str = new char[n];
      onechar = str[i];
      str[i] = '\0';
      strcpy(argstack[nargstack].str,&str[istart]);
      str[i] = onechar;

      argstack[nargstack++].flag = 1;

    // ----------------
    // Boolean operator, including end-of-string
    // ----------------

    } else if (strchr("<>=!&|\0",onechar)) {
      if (onechar == '=') {
        if (str[i+1] != '=')
          error->all(FLERR,"Invalid Boolean syntax in if command");
        op = EQ;
        i++;
      } else if (onechar == '!') {
        if (str[i+1] == '=') {
          op = NE;
          i++;
        } else op = NOT;
      } else if (onechar == '<') {
        if (str[i+1] != '=') op = LT;
        else {
          op = LE;
          i++;
        }
      } else if (onechar == '>') {
        if (str[i+1] != '=') op = GT;
        else {
          op = GE;
          i++;
        }
      } else if (onechar == '&') {
        if (str[i+1] != '&')
          error->all(FLERR,"Invalid Boolean syntax in if command");
        op = AND;
        i++;
      } else if (onechar == '|') {
        if (str[i+1] == '|') op = OR;
        else if (str[i+1] == '^') op = XOR;
        else error->all(FLERR,"Invalid Boolean syntax in if command");
        i++;
      } else op = DONE;

      i++;

      if (op == NOT && expect == ARG) {
        opstack[nopstack++] = op;
        continue;
      }

      if (expect == ARG)
        error->all(FLERR,"Invalid Boolean syntax in if command");
      expect = ARG;

      // evaluate stack as deep as possible while respecting precedence
      // before pushing current op onto stack

      while (nopstack && precedence[opstack[nopstack-1]] >= precedence[op]) {
        opprevious = opstack[--nopstack];

        nargstack--;
        flag2 = argstack[nargstack].flag;
        value2 = argstack[nargstack].value;
        str2 = argstack[nargstack].str;
        if (opprevious != NOT) {
          nargstack--;
          flag1 = argstack[nargstack].flag;
          value1 = argstack[nargstack].value;
          str1 = argstack[nargstack].str;
        }

        if (opprevious == NOT) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value2 == 0.0) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        } else if (opprevious == EQ) {
          if (flag1 != flag2)
            error->all(FLERR,"Invalid Boolean syntax in if command");
          if (flag2 == 0) {
            if (value1 == value2) argstack[nargstack].value = 1.0;
            else argstack[nargstack].value = 0.0;
          } else {
            if (strcmp(str1,str2) == 0) argstack[nargstack].value = 1.0;
            else argstack[nargstack].value = 0.0;
            delete [] str1;
            delete [] str2;
          }
        } else if (opprevious == NE) {
          if (flag1 != flag2)
            error->all(FLERR,"Invalid Boolean syntax in if command");
          if (flag2 == 0) {
            if (value1 != value2) argstack[nargstack].value = 1.0;
            else argstack[nargstack].value = 0.0;
          } else {
            if (strcmp(str1,str2) != 0) argstack[nargstack].value = 1.0;
            else argstack[nargstack].value = 0.0;
            delete [] str1;
            delete [] str2;
          }
        } else if (opprevious == LT) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value1 < value2) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        } else if (opprevious == LE) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value1 <= value2) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        } else if (opprevious == GT) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value1 > value2) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        } else if (opprevious == GE) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value1 >= value2) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        } else if (opprevious == AND) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value1 != 0.0 && value2 != 0.0) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        } else if (opprevious == OR) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value1 != 0.0 || value2 != 0.0) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        } else if (opprevious == XOR) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if ((value1 == 0.0 && value2 != 0.0) ||
              (value1 != 0.0 && value2 == 0.0))
            argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        }

        argstack[nargstack++].flag = 0;
      }

      // if end-of-string, break out of entire formula evaluation loop

      if (op == DONE) break;

      // push current operation onto stack

      opstack[nopstack++] = op;

    } else error->all(FLERR,"Invalid Boolean syntax in if command");
  }

  if (nopstack) error->all(FLERR,"Invalid Boolean syntax in if command");
  if (nargstack != 1) error->all(FLERR,"Invalid Boolean syntax in if command");
  return argstack[0].value;
}

/* ----------------------------------------------------------------------
   class to read variable values from a file
   for flag = SCALARFILE, reads one value per line
   for flag = ATOMFILE, reads set of one value per atom
------------------------------------------------------------------------- */

VarReader::VarReader(LAMMPS *lmp, char *name, char *file, int flag) :
  Pointers(lmp)
{
  me = comm->me;
  style = flag;
  fp = NULL;

  if (me == 0) {
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open file variable file %s",file);
      error->one(FLERR,str);
    }
  }

  // if atomfile-style variable, must store per-atom values read from file
  // allocate a new fix STORE, so they persist
  // id = variable-ID + VARIABLE_STORE, fix group = all

  fixstore = NULL;
  id_fix = NULL;
  buffer = NULL;

  if (style == ATOMFILE) {
    if (atom->map_style == 0)
      error->all(FLERR,
                 "Cannot use atomfile-style variable unless atom map exists");

    int n = strlen(name) + strlen("_VARIABLE_STORE") + 1;
    id_fix = new char[n];
    strcpy(id_fix,name);
    strcat(id_fix,"_VARIABLE_STORE");

    char **newarg = new char*[6];
    newarg[0] = id_fix;
    newarg[1] = (char *) "all";
    newarg[2] = (char *) "STORE";
    newarg[3] = (char *) "peratom";
    newarg[4] = (char *) "0";
    newarg[5] = (char *) "1";
    modify->add_fix(6,newarg);
    fixstore = (FixStore *) modify->fix[modify->nfix-1];
    delete [] newarg;

    buffer = new char[CHUNK*MAXLINE];
  }
}

/* ---------------------------------------------------------------------- */

VarReader::~VarReader()
{
  if (me == 0) {
    fclose(fp);
    fp = NULL;
  }

  // check modify in case all fixes have already been deleted

  if (fixstore) {
    if (modify) modify->delete_fix(id_fix);
    delete [] id_fix;
    delete [] buffer;
  }
}

/* ----------------------------------------------------------------------
   read for SCALARFILE style
   read next value from file into str for file-style variable
   strip comments, skip blank lines
   return 0 if successful, 1 if end-of-file
------------------------------------------------------------------------- */

int VarReader::read_scalar(char *str)
{
  int n;
  char *ptr;

  // read one string from file

  if (me == 0) {
    while (1) {
      if (fgets(str,MAXLINE,fp) == NULL) n = 0;
      else n = strlen(str);
      if (n == 0) break;                                 // end of file
      str[n-1] = '\0';                                   // strip newline
      if ((ptr = strchr(str,'#'))) *ptr = '\0';          // strip comment
      if (strtok(str," \t\n\r\f") == NULL) continue;     // skip if blank
      n = strlen(str) + 1;
      break;
    }
  }

  MPI_Bcast(&n,1,MPI_INT,0,world);
  if (n == 0) return 1;
  MPI_Bcast(str,n,MPI_CHAR,0,world);
  return 0;
}

/* ----------------------------------------------------------------------
   read snapshot of per-atom values from file
   into str for atomfile-style variable
   return 0 if successful, 1 if end-of-file
------------------------------------------------------------------------- */

int VarReader::read_peratom()
{
  int i,m,n,nchunk,eof;
  tagint tag;
  char *ptr,*next;
  double value;

  // set all per-atom values to 0.0
  // values that appear in file will overwrite this

  double *vstore = fixstore->vstore;

  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) vstore[i] = 0.0;

  // read one string from file, convert to Nlines

  char str[MAXLINE];
  if (me == 0) {
    while (1) {
      if (fgets(str,MAXLINE,fp) == NULL) n = 0;
      else n = strlen(str);
      if (n == 0) break;                                 // end of file
      str[n-1] = '\0';                                   // strip newline
      if ((ptr = strchr(str,'#'))) *ptr = '\0';          // strip comment
      if (strtok(str," \t\n\r\f") == NULL) continue;     // skip if blank
      n = strlen(str) + 1;
      break;
    }
  }

  MPI_Bcast(&n,1,MPI_INT,0,world);
  if (n == 0) return 1;

  MPI_Bcast(str,n,MPI_CHAR,0,world);
  bigint nlines = force->bnumeric(FLERR,str);
  tagint map_tag_max = atom->map_tag_max;

  bigint nread = 0;
  while (nread < nlines) {
    nchunk = MIN(nlines-nread,CHUNK);
    eof = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eof) return 1;

    char *buf = buffer;
    for (i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');
      *next = '\0';
      sscanf(buf,TAGINT_FORMAT " %lg",&tag,&value);
      if (tag <= 0 || tag > map_tag_max)
        error->one(FLERR,"Invalid atom ID in variable file");
      if ((m = atom->map(tag)) >= 0) vstore[m] = value;
      buf = next + 1;
    }

    nread += nchunk;
  }

  return 0;
}
