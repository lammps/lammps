/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "unistd.h"
#include "variable.h"
#include "universe.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "output.h"
#include "thermo.h"
#include "memory.h"
#include "error.h"

#define VARDELTA 4

enum{INDEX,LOOP,EQUAL,WORLD,UNIVERSE};

/* ---------------------------------------------------------------------- */

Variable::Variable()
{
  MPI_Comm_rank(world,&me);

  nvar = maxvar = 0;
  names = NULL;
  style = NULL;
  num = NULL;
  index = NULL;
  data = NULL;
}

/* ---------------------------------------------------------------------- */

Variable::~Variable()
{
  for (int i = 0; i < nvar; i++) {
    delete [] names[i];
    for (int j = 0; j < num[i]; j++) delete [] data[i][j];
    delete [] data[i];
  }
  memory->sfree(names);
  memory->sfree(style);
  memory->sfree(num);
  memory->sfree(index);
  memory->sfree(data);
}

/* ----------------------------------------------------------------------
   called by variable command in input script
------------------------------------------------------------------------- */

void Variable::set(int narg, char **arg)
{
  if (narg < 3) error->all("Illegal variable command");

  // if var already exists, just skip (except EQUAL vars)

  if (find(arg[0]) >= 0 && strcmp(arg[1],"equal") != 0) return;
      
  // make space for new variable

  if (nvar == maxvar) {
    maxvar += VARDELTA;
    names = (char **)
      memory->srealloc(names,maxvar*sizeof(char *),"var:names");
    style = (int *) memory->srealloc(style,maxvar*sizeof(int),"var:style");
    num = (int *) memory->srealloc(num,maxvar*sizeof(int),"var:num");
    index = (int *) memory->srealloc(index,maxvar*sizeof(int),"var:index");
    data = (char ***) 
      memory->srealloc(data,maxvar*sizeof(char **),"var:data");
  }

  // INDEX
  // num = listed args, index = 1st value, data = copied args

  if (strcmp(arg[1],"index") == 0) {
    style[nvar] = INDEX;
    num[nvar] = narg - 2;
    index[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(num[nvar],&arg[2],data[nvar]);

  // LOOP
  // num = N, index = 1st value, data = list of NULLS since never used

  } else if (strcmp(arg[1],"loop") == 0) {
    style[nvar] = LOOP;
    num[nvar] = atoi(arg[2]);
    index[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    for (int i = 0; i < num[nvar]; i++) data[nvar][i] = NULL;
    
  // EQUAL
  // remove pre-existing var if also style EQUAL (allows it to be reset)
  // num = 2, index = 1st value
  // data = 2 values, 1st is string to eval, 2nd is filled on retrieval

  } else if (strcmp(arg[1],"equal") == 0) {
    if (find(arg[0]) >= 0) {
      if (style[find(arg[0])] != EQUAL)
	error->all("Cannot redefine variable as a different style");
      remove(find(arg[0]));
    }
    style[nvar] = EQUAL;
    num[nvar] = 2;
    index[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(1,&arg[2],data[nvar]);
    data[nvar][1] = NULL;
    
  // WORLD
  // num = listed args, index = partition this proc is in, data = copied args
  // error check that num = # of worlds in universe

  } else if (strcmp(arg[1],"world") == 0) {
    style[nvar] = WORLD;
    num[nvar] = narg - 2;
    if (num[nvar] != universe->nworlds)
      error->all("World variable count doesn't match # of partitions");
    index[nvar] = universe->iworld;
    data[nvar] = new char*[num[nvar]];
    copy(num[nvar],&arg[2],data[nvar]);

  // UNIVERSE
  // num = listed args, index = partition this proc is in, data = copied args
  // universe proc 0 creates lock file
  // error check that all other universe variables are same length

  } else if (strcmp(arg[1],"universe") == 0) {
    style[nvar] = UNIVERSE;
    num[nvar] = narg - 2;
    if (num[nvar] < universe->nworlds)
      error->all("Universe variable count < # of partitions");
    index[nvar] = universe->iworld;
    data[nvar] = new char*[num[nvar]];
    copy(num[nvar],&arg[2],data[nvar]);

    if (universe->me == 0) {
      FILE *fp = fopen("tmp.lammps.variable","w");
      fprintf(fp,"%d\n",universe->nworlds);
      fclose(fp);
    }

    for (int jvar = 0; jvar < nvar; jvar++)
      if (num[jvar] && style[jvar] == UNIVERSE && num[nvar] != num[jvar])
	error->all("All universe variables must have same # of values");

    if (me == 0) {
      if (universe->uscreen)
	fprintf(universe->uscreen,
		"Initial ${%s} setting: value %d on partition %d\n",
		arg[0],index[nvar]+1,universe->iworld);
      if (universe->ulogfile)
	fprintf(universe->ulogfile,
		"Initial ${%s} setting: value %d on partition %d\n",
		arg[0],index[nvar]+1,universe->iworld);
    }

  } else error->all("Illegal variable command");

  // set variable name (after possible EQUAL remove)

  names[nvar] = new char[strlen(arg[0]) + 1];
  strcpy(names[nvar],arg[0]);
  nvar++;
}

/* ----------------------------------------------------------------------
   single-value INDEX variable created by command-line argument
------------------------------------------------------------------------- */

void Variable::set(char *name, char *value)
{
  int ivar = find(name);
  if (ivar >= 0) error->all("Command-line variable already exists");

  if (nvar == maxvar) {
    maxvar += VARDELTA;
    names = (char **)
      memory->srealloc(names,maxvar*sizeof(char *),"var:names");
    style = (int *) memory->srealloc(style,maxvar*sizeof(int),"var:style");
    num = (int *) memory->srealloc(num,maxvar*sizeof(int),"var:num");
    index = (int *) memory->srealloc(index,maxvar*sizeof(int),"var:index");
    data = (char ***) 
      memory->srealloc(data,maxvar*sizeof(char **),"var:data");
  }

  names[nvar] = new char[strlen(name) + 1];
  strcpy(names[nvar],name);
  style[nvar] = INDEX;
  num[nvar] = 1;
  index[nvar] = 0;
  data[nvar] = new char*[num[nvar]];
  copy(1,&value,data[nvar]);
  nvar++;
}

/* ----------------------------------------------------------------------
   increment variable(s)
   return 0 if OK if successfully incremented
   return 1 if any variable is exhausted, free the variable to allow re-use
------------------------------------------------------------------------- */

int Variable::next(int narg, char **arg)
{
  int ivar;

  if (narg == 0) error->all("Illegal next command");

  // check that variables exist and are all the same style

  for (int iarg = 0; iarg < narg; iarg++) {
    ivar = find(arg[iarg]);
    if (ivar == -1) error->all("Invalid variable in next command");
    if (style[ivar] != style[find(arg[0])])
      error->all("All variables in next command must be same style");
  }

  // check for invalid styles EQUAL or WORLD

  int istyle = style[find(arg[0])];
  if (istyle == EQUAL || istyle == WORLD)
    error->all("Invalid variable style with next command");

  // increment all variables in list
  // if any variable is exhausted, set flag = 1 and remove var to allow re-use

  int flag = 0;

  if (istyle == INDEX || istyle == LOOP) {
    for (int iarg = 0; iarg < narg; iarg++) {
      ivar = find(arg[iarg]);
      index[ivar]++;
      if (index[ivar] >= num[ivar]) {
	flag = 1;
	remove(ivar);
      }
    }

  } else if (istyle == UNIVERSE) {

    // wait until lock file can be created and owned by proc 0 of this world
    // read next available index and Bcast it within my world
    // set all variables in list to nextindex

    int nextindex;
    if (me == 0) {
      while (1) {
	if (!rename("tmp.lammps.variable","tmp.lammps.variable.lock")) break;
	usleep(100000);
      }
      FILE *fp = fopen("tmp.lammps.variable.lock","r");
      fscanf(fp,"%d",&nextindex);
      fclose(fp);
      fp = fopen("tmp.lammps.variable.lock","w");
      fprintf(fp,"%d\n",nextindex+1);
      fclose(fp);
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

    for (int iarg = 0; iarg < narg; iarg++) {
      ivar = find(arg[iarg]);
      index[ivar] = nextindex;
      if (index[ivar] >= num[ivar]) {
	flag = 1;
	remove(ivar);
      }
    }
  }

  return flag;
}

/* ----------------------------------------------------------------------
   return ptr to the data text associated with a variable
   return NULL if no variable or index is bad, caller must respond
------------------------------------------------------------------------- */

char *Variable::retrieve(char *name)
{
  int ivar = find(name);
  if (ivar == -1) return NULL;
  if (index[ivar] >= num[ivar]) return NULL;

  char *str;
  if (style[ivar] == INDEX) {
    str = data[ivar][index[ivar]];
  } else if (style[ivar] == LOOP) {
    char *value = new char[16];
    sprintf(value,"%d",index[ivar]+1);
    int n = strlen(value) + 1;
    if (data[ivar][0]) delete [] data[ivar][0];
    data[ivar][0] = new char[n];
    strcpy(data[ivar][0],value);
    delete [] value;
    str = data[ivar][0];
  } else if (style[ivar] == EQUAL) {
    char *value = evaluate(data[ivar][0]);
    int n = strlen(value) + 1;
    if (data[ivar][1]) delete [] data[ivar][1];
    data[ivar][1] = new char[n];
    strcpy(data[ivar][1],value);
    delete [] value;
    str = data[ivar][1];
  } else if (style[ivar] == WORLD) {
    str = data[ivar][index[ivar]];
  } else if (style[ivar] == UNIVERSE) {
    str = data[ivar][index[ivar]];
  }
  return str;
}

/* ----------------------------------------------------------------------
   search for name in list of variables names
   return index or -1 if not found
------------------------------------------------------------------------- */
  
int Variable::find(char *name)
{
  for (int i = 0; i < nvar; i++)
    if (strcmp(name,names[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   remove Nth variable from list and compact list
------------------------------------------------------------------------- */
  
void Variable::remove(int n)
{
  delete [] names[n];
  for (int i = 0; i < num[n]; i++) delete [] data[n][i];
  delete [] data[n];

  for (int i = n+1; i < nvar; i++) {
    names[i-1] = names[i];
    style[i-1] = style[i];
    num[i-1] = num[i];
    index[i-1] = index[i];
    data[i-1] = data[i];
  }
  nvar--;
}

/* ----------------------------------------------------------------------
   copy narg strings from **from to **to 
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
   recursive method to evaluate a string
   string can be a number: 0.0, -5.45, etc
   string can be a keyword: lx, ly, lz, vol, etc
   string can be a vector: x[123], y[3], vx[34], etc
     value inside brackets must be global ID from 1 to Natoms
   string can be a function: div(x,y), mult(x,y), add(x,y), xcm(group,x), etc
     function args can be combo of  number, keyword, vector, fn(), group, etc
   see lists of valid keywords, vectors, and functions below
   when string is evaluated, put result in a newly allocated string
   return address of result string (will be freed by caller)
------------------------------------------------------------------------- */

char *Variable::evaluate(char *str)
{
  // allocate a new string for the eventual result

  char *result = new char[32];
  double answer;

  // if string is "function", grab one or two args, evaulate args recursively

  if (strchr(str,'(')) {
    if (str[strlen(str)-1] != ')')
      error->all("Cannot evaluate variable equal command");

    char *ptr = strchr(str,'(');
    int n = ptr - str;
    char *func = new char[n+1];
    strncpy(func,str,n);
    func[n] = '\0';

    char *comma = ++ptr;
    int level = 0;
    while (1) {
      if (*comma == '\0')
	error->all("Cannot evaluate variable equal command");
      else if (*comma == ',' && level == 0) break;
      else if (*comma == ')' && level == 0) break;
      else if (*comma == '(') level++;
      else if (*comma == ')') level--;
      comma++;
    }

    char *arg1,*arg2;
    n = comma - ptr;
    arg1 = new char[n+1];
    strncpy(arg1,ptr,n);
    arg1[n] = '\0';

    if (*comma == ',') {
      ptr = comma + 1;
      comma = &str[strlen(str)-1];
      n = comma - ptr;
      arg2 = new char[n+1];
      strncpy(arg2,ptr,n);
      arg2[n] = '\0';
    } else arg2 = NULL;
    
    double value1,value2;
    char *strarg1 = NULL;
    char *strarg2 = NULL;

    // customize by adding function to this list and to if statement
    // math functions: add(x,y),sub(x,y),mult(x,y),div(x,y),
    //                 neg(x),pow(x,y),exp(x),ln(x),sqrt(x)
    // group functions: mass(group),charge(group),xcm(group,dim),
    //                  vcm(group,dim),bound(group,xmin),gyration(group)

    if (strcmp(func,"add") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      strarg2 = evaluate(arg2);
      value1 = atof(strarg1);
      value2 = atof(strarg2);
      answer = value1 + value2;

    } else if (strcmp(func,"sub") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      strarg2 = evaluate(arg2);
      value1 = atof(strarg1);
      value2 = atof(strarg2);
      answer = value1 - value2;

    } else if (strcmp(func,"mult") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      strarg2 = evaluate(arg2);
      value1 = atof(strarg1);
      value2 = atof(strarg2);
      answer = value1 * value2;

    } else if (strcmp(func,"div") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      strarg2 = evaluate(arg2);
      value1 = atof(strarg1);
      value2 = atof(strarg2);
      if (value2 == 0.0)
	error->all("Cannot evaluate variable equal command");
      answer = value1 / value2;

    } else if (strcmp(func,"neg") == 0) {
      if (arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      value1 = atof(strarg1);
      answer = -value1;

    } else if (strcmp(func,"pow") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      strarg2 = evaluate(arg2);
      value1 = atof(strarg1);
      value2 = atof(strarg2);
      if (value2 == 0.0)
	error->all("Cannot evaluate variable equal command");
      answer = pow(value1,value2);

    } else if (strcmp(func,"exp") == 0) {
      if (arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      value1 = atof(strarg1);
      answer = exp(value1);

    } else if (strcmp(func,"ln") == 0) {
      if (arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      value1 = atof(strarg1);
      if (value1 == 0.0)
	error->all("Cannot evaluate variable equal command");
      answer = log(value1);

    } else if (strcmp(func,"sqrt") == 0) {
      if (arg2) error->all("Cannot evaluate variable equal command");
      strarg1 = evaluate(arg1);
      value1 = atof(strarg1);
      if (value1 == 0.0)
	error->all("Cannot evaluate variable equal command");
      answer = sqrt(value1);

    } else if (strcmp(func,"mass") == 0) {
      if (arg2) error->all("Cannot evaluate variable equal command");
      int igroup = group->find(arg1);
      if (igroup == -1) error->all("Variable equal group ID does not exist");
      atom->check_mass();
      answer = group->mass(igroup);

    } else if (strcmp(func,"charge") == 0) {
      if (arg2) error->all("Cannot evaluate variable equal command");
      int igroup = group->find(arg1);
      if (igroup == -1) error->all("Variable equal group ID does not exist");
      answer = group->charge(igroup);

    } else if (strcmp(func,"xcm") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      int igroup = group->find(arg1);
      if (igroup == -1) error->all("Variable equal group ID does not exist");
      atom->check_mass();
      double masstotal = group->mass(igroup);
      double xcm[3];
      group->xcm(igroup,masstotal,xcm);
      if (strcmp(arg2,"x") == 0) answer = xcm[0];
      else if (strcmp(arg2,"y") == 0) answer = xcm[1];
      else if (strcmp(arg2,"z") == 0) answer = xcm[2];
      else error->all("Cannot evaluate variable equal command");

    } else if (strcmp(func,"vcm") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      int igroup = group->find(arg1);
      if (igroup == -1) error->all("Variable equal group ID does not exist");
      atom->check_mass();
      double masstotal = group->mass(igroup);
      double vcm[3];
      group->vcm(igroup,masstotal,vcm);
      if (strcmp(arg2,"x") == 0) answer = vcm[0];
      else if (strcmp(arg2,"y") == 0) answer = vcm[1];
      else if (strcmp(arg2,"z") == 0) answer = vcm[2];
      else error->all("Cannot evaluate variable equal command");

    } else if (strcmp(func,"bound") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      int igroup = group->find(arg1);
      if (igroup == -1) error->all("Variable equal group ID does not exist");
      double minmax[6];
      group->bounds(igroup,minmax);
      if (strcmp(arg2,"xmin") == 0) answer = minmax[0];
      else if (strcmp(arg2,"xmax") == 0) answer = minmax[1];
      else if (strcmp(arg2,"ymin") == 0) answer = minmax[2];
      else if (strcmp(arg2,"ymax") == 0) answer = minmax[3];
      else if (strcmp(arg2,"zmin") == 0) answer = minmax[4];
      else if (strcmp(arg2,"zmax") == 0) answer = minmax[5];
      else error->all("Cannot evaluate variable equal command");

    } else if (strcmp(func,"gyration") == 0) {
      if (!arg2) error->all("Cannot evaluate variable equal command");
      int igroup = group->find(arg1);
      if (igroup == -1) error->all("Variable equal group ID does not exist");
      atom->check_mass();
      double masstotal = group->mass(igroup);
      double xcm[3];
      group->xcm(igroup,masstotal,xcm);
      answer = group->gyration(igroup,masstotal,xcm);

    } else error->all("Cannot evaluate variable equal command");

    delete [] func;
    delete [] arg1;
    delete [] strarg1;
    if (arg2) {
      delete [] arg2;
      delete [] strarg2;
    }
    sprintf(result,"%g",answer);

  // if string is "vector", find which proc owns atom, grab vector value

  } else if (strchr(str,'[')) {
    if (str[strlen(str)-1] != ']')
      error->all("Cannot evaluate variable equal command");
    if (atom->map_style == 0)
      error->all("Cannot use vectors in variables unless atom map exists");

    char *ptr = strchr(str,'[');
    int n = ptr - str;
    char *vector = new char[n+1];
    strncpy(vector,str,n);
    vector[n] = '\0';

    char *arg;
    ptr++;
    char *ptr2 = &str[strlen(str)-1];
    n = ptr2 - ptr;
    arg = new char[n+1];
    strncpy(arg,ptr,n);
    arg[n] = '\0';
    int i = atom->map(atoi(arg));

    double mine;
    if (i >= 0 && i < atom->nlocal) {

      // customize by adding vector to this list and to if statement
      // x,y,z,vx,vy,vz,fx,fy,fz

      if (strcmp(vector,"x") == 0) mine = atom->x[i][0];
      else if (strcmp(vector,"y") == 0) mine = atom->x[i][1];
      else if (strcmp(vector,"z") == 0) mine = atom->x[i][2];
      else if (strcmp(vector,"vx") == 0) mine = atom->v[i][0];
      else if (strcmp(vector,"vy") == 0) mine = atom->v[i][1];
      else if (strcmp(vector,"vz") == 0) mine = atom->v[i][2];
      else if (strcmp(vector,"fx") == 0) mine = atom->f[i][0];
      else if (strcmp(vector,"fy") == 0) mine = atom->f[i][1];
      else if (strcmp(vector,"fz") == 0) mine = atom->f[i][2];

      else error->one("Invalid vector in variable equal command");

    } else mine = 0.0;
    MPI_Allreduce(&mine,&answer,1,MPI_DOUBLE,MPI_SUM,world);

    delete [] vector;
    delete [] arg;
    sprintf(result,"%g",answer);

  // if string is "keyword", compute appropriate value via thermo

  } else if (str[0] - 'a' >= 0 && str[0] - 'a' < 26) {

    // customize by adding keyword to this list and to if statement
    // keywords defined by thermo_style custom are passed to compute_value()

    if (domain->box_exist == 0)
      error->all("Using variable equal keyword before simulation box is defined");

    // process these keywords here, not in thermo
    // they are allowed anytime after simulation box is defined

    if (strcmp(str,"step") == 0)
      answer = update->ntimestep;
    else if (strcmp(str,"atoms") == 0)
      answer = atom->natoms;
    else if (strcmp(str,"lx") == 0)
      answer = domain->xprd;
    else if (strcmp(str,"ly") == 0)
      answer = domain->yprd;
    else if (strcmp(str,"lz") == 0)
      answer = domain->zprd;
    else if (strcmp(str,"vol") == 0)
      answer = domain->xprd * domain->yprd * domain->zprd;

    // process other keywords via thermo
    // error if 1st run has not yet started
    // warning if out-of-sync with thermo, since values can be bogus

    else {
      if (update->first_update == 0)
	error->all("Using variable equal keyword before initial run");
      if (update->ntimestep != output->next_thermo && me == 0)
	error->warning("Using variable equal keyword with non-current thermo");
      int flag = output->thermo->compute_value(str,&answer);
      if (flag) error->all("Invalid keyword in variable equal command");
    }
    
    sprintf(result,"%g",answer);

  // string is a number, just copy to result

  } else strcpy(result,str);

  // return newly allocated string

  return result;
}
