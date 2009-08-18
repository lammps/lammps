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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "unistd.h"
#include "variable.h"
#include "universe.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "output.h"
#include "thermo.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define VARDELTA 4
#define MAXLEVEL 4

#define MYROUND(a) (( a-floor(a) ) >= .5) ? ceil(a) : floor(a)

enum{DUMMY0,INVOKED_SCALAR,INVOKED_VECTOR,DUMMMY3,INVOKED_PERATOM};
enum{INDEX,LOOP,EQUAL,WORLD,UNIVERSE,ULOOP,ATOM};
enum{ARG,OP};
enum{DONE,ADD,SUBTRACT,MULTIPLY,DIVIDE,CARAT,UNARY,
       SQRT,EXP,LN,LOG,SIN,COS,TAN,ASIN,ACOS,ATAN,
       CEIL,FLOOR,ROUND,VALUE,ATOMARRAY,TYPEARRAY,INTARRAY};

/* ---------------------------------------------------------------------- */

Variable::Variable(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);

  nvar = maxvar = 0;
  names = NULL;
  style = NULL;
  num = NULL;
  index = NULL;
  data = NULL;

  precedence[DONE] = 0;
  precedence[ADD] = precedence[SUBTRACT] = 1;
  precedence[MULTIPLY] = precedence[DIVIDE] = 2;
  precedence[CARAT] = 3;
  precedence[UNARY] = 4;
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
  if (narg < 2) error->all("Illegal variable command");

  // DELETE
  // doesn't matter if variable no longer exists

  if (strcmp(arg[1],"delete") == 0) {
    if (narg != 2) error->all("Illegal variable command");
    if (find(arg[0]) >= 0) remove(find(arg[0]));
    return;

  // INDEX
  // num = listed args, index = 1st value, data = copied args

  } else if (strcmp(arg[1],"index") == 0) {
    if (narg < 3) error->all("Illegal variable command");
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) extend();
    style[nvar] = INDEX;
    num[nvar] = narg - 2;
    index[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(num[nvar],&arg[2],data[nvar]);

  // LOOP
  // num = N, index = 1st value, data = list of NULLS since never used

  } else if (strcmp(arg[1],"loop") == 0) {
    if (narg != 3) error->all("Illegal variable command");
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) extend();
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
    if (narg != 3) error->all("Illegal variable command");
    if (find(arg[0]) >= 0) {
      if (style[find(arg[0])] != EQUAL)
	error->all("Cannot redefine variable as a different style");
      remove(find(arg[0]));
    }
    if (nvar == maxvar) extend();
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
    if (narg < 3) error->all("Illegal variable command");
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) extend();
    style[nvar] = WORLD;
    num[nvar] = narg - 2;
    if (num[nvar] != universe->nworlds)
      error->all("World variable count doesn't match # of partitions");
    index[nvar] = universe->iworld;
    data[nvar] = new char*[num[nvar]];
    copy(num[nvar],&arg[2],data[nvar]);

  // UNIVERSE and ULOOP
  // for UNIVERSE: num = listed args, data = copied args
  // for ULOOP: num = N, data = list of NULLS since never used
  // index = partition this proc is in
  // universe proc 0 creates lock file
  // error check that all other universe/uloop variables are same length

  } else if (strcmp(arg[1],"universe") == 0 || strcmp(arg[1],"uloop") == 0) {
    if (strcmp(arg[1],"universe") == 0) {
      if (narg < 3) error->all("Illegal variable command");
      if (find(arg[0]) >= 0) return;
      if (nvar == maxvar) extend();
      style[nvar] = UNIVERSE;
      num[nvar] = narg - 2;
      data[nvar] = new char*[num[nvar]];
      copy(num[nvar],&arg[2],data[nvar]);
    } else {
      if (narg != 3) error->all("Illegal variable command");
      if (find(arg[0]) >= 0) return;
      if (nvar == maxvar) extend();
      style[nvar] = ULOOP;
      num[nvar] = atoi(arg[2]);
      data[nvar] = new char*[num[nvar]];
      for (int i = 0; i < num[nvar]; i++) data[nvar][i] = NULL;
    }

    if (num[nvar] < universe->nworlds)
      error->all("Universe/uloop variable count < # of partitions");
    index[nvar] = universe->iworld;

    if (universe->me == 0) {
      FILE *fp = fopen("tmp.lammps.variable","w");
      fprintf(fp,"%d\n",universe->nworlds);
      fclose(fp);
    }

    for (int jvar = 0; jvar < nvar; jvar++)
      if (num[jvar] && (style[jvar] == UNIVERSE || style[jvar] == ULOOP) && 
	  num[nvar] != num[jvar])
	error->all("All universe/uloop variables must have same # of values");

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

  // ATOM
  // remove pre-existing var if also style ATOM (allows it to be reset)
  // num = 1, index = 1st value
  // data = 1 value, string to eval

  } else if (strcmp(arg[1],"atom") == 0) {
    if (narg != 3) error->all("Illegal variable command");
    if (find(arg[0]) >= 0) {
      if (style[find(arg[0])] != ATOM)
	error->all("Cannot redefine variable as a different style");
      remove(find(arg[0]));
    }
    if (nvar == maxvar) extend();
    style[nvar] = ATOM;
    num[nvar] = 1;
    index[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(1,&arg[2],data[nvar]);
    
  } else error->all("Illegal variable command");

  // set name of variable
  // must come at end, since EQUAL/ATOM reset may have removed name
  // name must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  names[nvar] = new char[n];
  strcpy(names[nvar],arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(names[nvar][i]) && names[nvar][i] != '_')
      error->all("Variable name must be alphanumeric or underscore characters");

  nvar++;
}

/* ----------------------------------------------------------------------
   single-value INDEX variable created by command-line argument
------------------------------------------------------------------------- */

void Variable::set(char *name, char *value)
{
  char **newarg = new char*[3];
  newarg[0] = name;
  newarg[1] = (char *) "index";
  newarg[2] = value;
  set(3,newarg);
  delete [] newarg;
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
  // exception: UNIVERSE and ULOOP variables can be mixed in same next command

  for (int iarg = 0; iarg < narg; iarg++) {
    ivar = find(arg[iarg]);
    if (ivar == -1) error->all("Invalid variable in next command");
    if (style[ivar] == ULOOP && style[find(arg[0])] == UNIVERSE) continue;
    else if (style[ivar] == UNIVERSE && style[find(arg[0])] == ULOOP) continue;
    else if (style[ivar] != style[find(arg[0])])
      error->all("All variables in next command must be same style");
  }

  // invalid styles EQUAL or WORLD or ATOM

  int istyle = style[find(arg[0])];
  if (istyle == EQUAL || istyle == WORLD || istyle == ATOM)
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

  } else if (istyle == UNIVERSE || istyle == ULOOP) {

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
   if EQUAL var, evaluates variable and puts result in str
   return NULL if no variable or index is bad, caller must respond
------------------------------------------------------------------------- */

char *Variable::retrieve(char *name)
{
  int ivar = find(name);
  if (ivar == -1) return NULL;
  if (index[ivar] >= num[ivar]) return NULL;

  char *str;
  if (style[ivar] == INDEX || style[ivar] == WORLD || 
      style[ivar] == UNIVERSE) {
    str = data[ivar][index[ivar]];
  } else if (style[ivar] == LOOP || style[ivar] == ULOOP) {
    char *value = new char[16];
    sprintf(value,"%d",index[ivar]+1);
    int n = strlen(value) + 1;
    if (data[ivar][0]) delete [] data[ivar][0];
    data[ivar][0] = new char[n];
    strcpy(data[ivar][0],value);
    delete [] value;
    str = data[ivar][0];
  } else if (style[ivar] == EQUAL) {
    char result[32];
    double answer = evaluate(data[ivar][0],NULL);
    sprintf(result,"%.10g",answer);
    int n = strlen(result) + 1;
    if (data[ivar][1]) delete [] data[ivar][1];
    data[ivar][1] = new char[n];
    strcpy(data[ivar][1],result);
    str = data[ivar][1];
  } else if (style[ivar] == ATOM) return NULL;

  return str;
}

/* ----------------------------------------------------------------------
   return result of equal-style variable evaluation
------------------------------------------------------------------------- */

double Variable::compute_equal(int ivar)
{
  return evaluate(data[ivar][0],NULL);
}

/* ----------------------------------------------------------------------
   compute result of atom-style variable evaluation
   stride used since result may not be contiguous memory locs
   if sumflag, add variable values to existing result
------------------------------------------------------------------------- */

void Variable::compute_atom(int ivar, int igroup,
			    double *result, int stride, int sumflag)
{
  Tree *tree;
  double tmp = evaluate(data[ivar][0],&tree);

  int groupbit = group->bitmask[igroup];
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (sumflag == 0) {
    int m = 0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] && groupbit) result[m] = eval_tree(tree,i);
      else result[m] = 0.0;
      m += stride;
    }

  } else {
    int m = 0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] && groupbit) result[m] += eval_tree(tree,i);
      m += stride;
    }
  }

  free_tree(tree);
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
   return 1 if variable is EQUAL style, 0 if not
------------------------------------------------------------------------- */
  
int Variable::equalstyle(int ivar)
{
  if (style[ivar] == EQUAL) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return 1 if variable is ATOM style, 0 if not
------------------------------------------------------------------------- */
  
int Variable::atomstyle(int ivar)
{
  if (style[ivar] == ATOM) return 1;
  return 0;
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
  make space in arrays for new variable
------------------------------------------------------------------------- */

void Variable::extend()
{
  maxvar += VARDELTA;
  names = (char **)
    memory->srealloc(names,maxvar*sizeof(char *),"var:names");
  style = (int *) memory->srealloc(style,maxvar*sizeof(int),"var:style");
  num = (int *) memory->srealloc(num,maxvar*sizeof(int),"var:num");
  index = (int *) memory->srealloc(index,maxvar*sizeof(int),"var:index");
  data = (char ***) 
    memory->srealloc(data,maxvar*sizeof(char **),"var:data");
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
   recursive evaluation of a string str
   string is a equal-style or atom-style formula containing one or more items:
     number = 0.0, -5.45, 2.8e-4, ...
     thermo keyword = ke, vol, atoms, ...
     math operation = (),-x,x+y,x-y,x*y,x/y,x^y,
                      sqrt(x),exp(x),ln(x),log(x),
		      sin(x),cos(x),tan(x),asin(x),acos(x),atan(x)
     group function = count(group), mass(group), xcm(group,x), ...
     atom value = x[123], y[3], vx[34], ...
     atom vector = x[], y[], vx[], ...
     compute = c_ID, c_ID[2], c_ID[], c_ID[][2], c_ID[N], c_ID[N][2]
     fix = f_ID, f_ID[2], f_ID[], f_ID[][2], f_ID[N], f_ID[N][2]
     variable = v_name, v_name[], v_name[N]
   if tree ptr passed in, return ptr to parse tree created for formula
   if no tree ptr passed in, return value of string as double
------------------------------------------------------------------------- */

double Variable::evaluate(char *str, Tree **tree)
{
  int op,opprevious;
  double value1,value2;
  char onechar;

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
      if (expect == OP) error->all("Invalid syntax in variable formula");
      expect = OP;

      char *contents;
      i = find_matching_paren(str,i,contents);
      i++;

      // evaluate contents and push on stack

      if (tree) {
	Tree *newtree;
	double tmp = evaluate(contents,&newtree);
	treestack[ntreestack++] = newtree;
      } else argstack[nargstack++] = evaluate(contents,NULL);

      delete [] contents;

    // ----------------
    // number: push value onto stack
    // ----------------

    } else if (isdigit(onechar) || onechar == '.') {
      if (expect == OP) error->all("Invalid syntax in variable formula");
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
	newtree->left = newtree->right = NULL;
	treestack[ntreestack++] = newtree;
      } else argstack[nargstack++] = atof(number);

      delete [] number;

    // ----------------
    // letter: c_ID, f_ID, v_name, exp(), xcm(,), x[234], x[], vol
    // ----------------

    } else if (islower(onechar)) {
      if (expect == OP) error->all("Invalid syntax in variable formula");
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

      if (strncmp(word,"c_",2) == 0) {
	n = strlen(word) - 2 + 1;
	char *id = new char[n];
	strcpy(id,&word[2]);

	int icompute = modify->find_compute(id);
	if (icompute < 0) error->all("Invalid compute ID in variable formula");
	Compute *compute = modify->compute[icompute];
	delete [] id;

	if (domain->box_exist == 0)
	  error->all("Variable evaluation before simulation box is defined");
 
	// parse zero or one or two trailing brackets
	// point i beyond last bracket
	// nbracket = # of bracket pairs
	// index1,index2 = int inside each bracket pair, 0 if empty

	int nbracket,index1,index2;
	if (str[i] != '[') nbracket = 0;
	else {
	  nbracket = 1;
	  i = int_between_brackets(str,i,index1,1);
	  i++;
	  if (str[i] == '[') {
	    nbracket = 2;
	    i = int_between_brackets(str,i,index2,0);
	    i++;
	  }
	}

        // c_ID = global scalar

	if (nbracket == 0 && compute->scalar_flag) {

	  if (update->whichflag == 0) {
	    if (compute->invoked_scalar != update->ntimestep)
	      error->all("Compute used in variable between runs "
			 "is not current");
	  } else if (!(compute->invoked_flag & INVOKED_SCALAR)) {
	    compute->compute_scalar();
	    compute->invoked_flag |= INVOKED_SCALAR;
	  }

	  value1 = compute->scalar;
	  if (tree) {
	    Tree *newtree = new Tree();
	    newtree->type = VALUE;
	    newtree->value = value1;
	    newtree->left = newtree->right = NULL;
	    treestack[ntreestack++] = newtree;
	  } else argstack[nargstack++] = value1;

        // c_ID[2] = global vector

	} else if (nbracket == 1 && index1 > 0 && compute->vector_flag) {

	  if (index1 > compute->size_vector)
	      error->all("Compute vector in variable formula is too small");
	  if (update->whichflag == 0) {
	    if (compute->invoked_vector != update->ntimestep)
	      error->all("Compute used in variable between runs "
			 "is not current");
	  } else if (!(compute->invoked_flag & INVOKED_VECTOR)) {
	    compute->compute_vector();
	    compute->invoked_flag |= INVOKED_VECTOR;
	  }

	  value1 = compute->vector[index1-1];
	  if (tree) {
	    Tree *newtree = new Tree();
	    newtree->type = VALUE;
	    newtree->value = value1;
	    newtree->left = newtree->right = NULL;
	    treestack[ntreestack++] = newtree;
	  } else argstack[nargstack++] = value1;

        // c_ID[] = per-atom scalar

	} else if (nbracket == 1 && index1 == 0 && 
		   compute->peratom_flag && compute->size_peratom == 0) {

	  if (tree == NULL)
	    error->all("Per-atom compute in equal-style variable formula");
	  if (update->whichflag == 0) {
	    if (compute->invoked_peratom != update->ntimestep)
	      error->all("Compute used in variable between runs "
			 "is not current");
	  } else if (!(compute->invoked_flag & INVOKED_PERATOM)) {
	    compute->compute_peratom();
	    compute->invoked_flag |= INVOKED_PERATOM;
	  }

	  Tree *newtree = new Tree();
	  newtree->type = ATOMARRAY;
	  newtree->array = compute->scalar_atom;
	  newtree->nstride = 1;
	  newtree->left = newtree->right = NULL;
	  treestack[ntreestack++] = newtree;

        // c_ID[N] = global value from per-atom scalar

	} else if (nbracket == 1 && index1 > 0 && 
		   compute->peratom_flag && compute->size_peratom == 0) {

	  if (update->whichflag == 0) {
	    if (compute->invoked_peratom != update->ntimestep)
	      error->all("Compute used in variable between runs "
			 "is not current");
	  } else if (!(compute->invoked_flag & INVOKED_PERATOM)) {
	    compute->compute_peratom();
	    compute->invoked_flag |= INVOKED_PERATOM;
	  }

	  peratom2global(1,NULL,compute->scalar_atom,1,index1,
			 tree,treestack,ntreestack,argstack,nargstack);

        // c_ID[][2] = per-atom vector

	} else if (nbracket == 2 && index1 == 0 && index2 > 0 &&
		   compute->peratom_flag) {

	  if (tree == NULL)
	    error->all("Per-atom compute in equal-style variable formula");
	  if (index2 > compute->size_peratom)
	    error->all("Compute vector in variable formula is too small");
	  if (update->whichflag == 0) {
	    if (compute->invoked_peratom != update->ntimestep)
	      error->all("Compute used in variable between runs "
			 "is not current");
	  } else if (!(compute->invoked_flag & INVOKED_PERATOM)) {
	    compute->compute_peratom();
	    compute->invoked_flag |= INVOKED_PERATOM;
	  }

	  Tree *newtree = new Tree();
	  newtree->type = ATOMARRAY;
	  newtree->array = &compute->vector_atom[0][index2-1];
	  newtree->nstride = compute->size_peratom;
	  newtree->left = newtree->right = NULL;
	  treestack[ntreestack++] = newtree;

        // c_ID[N][2] = global value from per-atom vector

	} else if (nbracket == 2 && index1 > 0 && index2 > 0 &&
		   compute->peratom_flag) {

	  if (index2 > compute->size_peratom)
	    error->all("Compute vector in variable formula is too small");
	  if (update->whichflag == 0) {
	    if (compute->invoked_peratom != update->ntimestep)
	      error->all("Compute used in variable between runs "
			 "is not current");
	  } else if (!(compute->invoked_flag & INVOKED_PERATOM)) {
	    compute->compute_peratom();
	    compute->invoked_flag |= INVOKED_PERATOM;
	  }

	  peratom2global(1,NULL,&compute->vector_atom[0][index2-1],
			 compute->size_peratom,index1,
			 tree,treestack,ntreestack,argstack,nargstack);

	} else error->all("Mismatched compute in variable formula");

      // ----------------
      // fix
      // ----------------

      } else if (strncmp(word,"f_",2) == 0) {
	n = strlen(word) - 2 + 1;
	char *id = new char[n];
	strcpy(id,&word[2]);

	int ifix = modify->find_fix(id);
	if (ifix < 0) error->all("Invalid fix ID in variable formula");
	Fix *fix = modify->fix[ifix];
	delete [] id;

	if (domain->box_exist == 0)
	  error->all("Variable evaluation before simulation box is defined");
 
	// parse zero or one or two trailing brackets
	// point i beyond last bracket
	// nbracket = # of bracket pairs
	// index1,index2 = int inside each bracket pair, 0 if empty

	int nbracket,index1,index2;
	if (str[i] != '[') nbracket = 0;
	else {
	  nbracket = 1;
	  i = int_between_brackets(str,i,index1,1);
	  i++;
	  if (str[i] == '[') {
	    nbracket = 2;
	    i = int_between_brackets(str,i,index2,0);
	    i++;
	  }
	}

        // f_ID = global scalar

	if (nbracket == 0 && fix->scalar_flag) {

	  if (update->whichflag > 0 &&
	      update->ntimestep % fix->scalar_vector_freq)
	    error->all("Fix in variable not computed at compatible time");
	  value1 = fix->compute_scalar();
	  if (tree) {
	    Tree *newtree = new Tree();
	    newtree->type = VALUE;
	    newtree->value = value1;
	    newtree->left = newtree->right = NULL;
	    treestack[ntreestack++] = newtree;
	  } else argstack[nargstack++] = value1;

        // f_ID[2] = global vector

	} else if (nbracket == 1 && index1 > 0 && fix->vector_flag) {

	  if (index1 > fix->size_vector)
	      error->all("Fix vector in variable formula is too small");
	  if (update->whichflag > 0 && 
	      update->ntimestep % fix->scalar_vector_freq)
	    error->all("Fix in variable not computed at compatible time");
	  value1 = fix->compute_vector(index1-1);
	  if (tree) {
	    Tree *newtree = new Tree();
	    newtree->type = VALUE;
	    newtree->value = value1;
	    newtree->left = newtree->right = NULL;
	    treestack[ntreestack++] = newtree;
	  } else argstack[nargstack++] = value1;

        // f_ID[] = per-atom scalar

	} else if (nbracket == 1 && index1 == 0 && 
		   fix->peratom_flag && fix->size_peratom == 0) {

	  if (tree == NULL)
	    error->all("Per-atom fix in equal-style variable formula");
	  if (update->whichflag > 0 && 
	      update->ntimestep % fix->peratom_freq)
	    error->all("Fix in variable not computed at compatible time");
	  Tree *newtree = new Tree();
	  newtree->type = ATOMARRAY;
	  newtree->array = fix->scalar_atom;
	  newtree->nstride = 1;
	  newtree->left = newtree->right = NULL;
	  treestack[ntreestack++] = newtree;

        // f_ID[N] = global value from per-atom scalar

	} else if (nbracket == 1 && index1 > 0 && 
		   fix->peratom_flag && fix->size_peratom == 0) {

	  if (update->whichflag > 0 && 
	      update->ntimestep % fix->peratom_freq)
	    error->all("Fix in variable not computed at compatible time");
	  peratom2global(1,NULL,fix->scalar_atom,1,index1,
			 tree,treestack,ntreestack,argstack,nargstack);

        // f_ID[][2] = per-atom vector

	} else if (nbracket == 2 && index1 == 0 && index2 > 0 &&
		   fix->peratom_flag) {

	  if (tree == NULL)
	    error->all("Per-atom fix in equal-style variable formula");
	  if (index2 > fix->size_peratom)
	    error->all("Fix vector in variable formula is too small");
	  if (update->whichflag > 0 && 
	      update->ntimestep % fix->peratom_freq)
	    error->all("Fix in variable not computed at compatible time");
	  Tree *newtree = new Tree();
	  newtree->type = ATOMARRAY;
	  newtree->array = &fix->vector_atom[0][index2-1];
	  newtree->nstride = fix->size_peratom;
	  newtree->left = newtree->right = NULL;
	  treestack[ntreestack++] = newtree;

        // f_ID[N][2] = global value from per-atom vector

	} else if (nbracket == 2 && index1 > 0 && index2 > 0 &&
		   fix->peratom_flag) {

	  if (index2 > fix->size_peratom)
	    error->all("Fix vector in variable formula is too small");
	  if (update->whichflag > 0 && 
	      update->ntimestep % fix->peratom_freq)
	    error->all("Fix in variable not computed at compatible time");
	  peratom2global(1,NULL,&fix->vector_atom[0][index2-1],
			 fix->size_peratom,index1,
			 tree,treestack,ntreestack,argstack,nargstack);

	} else error->all("Mismatched fix in variable formula");

      // ----------------
      // variable
      // ----------------

      } else if (strncmp(word,"v_",2) == 0) {
	n = strlen(word) - 2 + 1;
	char *id = new char[n];
	strcpy(id,&word[2]);

	int ivar = find(id);
	if (ivar < 0) error->all("Invalid variable name in variable formula");

	// parse zero or one trailing brackets
	// point i beyond last bracket
	// nbracket = # of bracket pairs
	// index = int inside bracket, 0 if empty

	int nbracket,index;
	if (str[i] != '[') nbracket = 0;
	else {
	  nbracket = 1;
	  i = int_between_brackets(str,i,index,1);
	  i++;
	}

        // v_name = non atom-style variable = global value

	if (nbracket == 0 && style[ivar] != ATOM) {

	  char *var = retrieve(id);
	  if (var == NULL)
	    error->all("Invalid variable evaluation in variable formula");
	  if (tree) {
	    Tree *newtree = new Tree();
	    newtree->type = VALUE;
	    newtree->value = atof(var);
	    newtree->left = newtree->right = NULL;
	    treestack[ntreestack++] = newtree;
	  } else argstack[nargstack++] = atof(var);

        // v_name[] = atom-style variable

	} else if (nbracket && index == 0 && style[ivar] == ATOM) {

	  if (tree == NULL)
	    error->all("Atom-style variable in equal-style variable formula");
	  Tree *newtree;
	  double tmp = evaluate(data[ivar][0],&newtree);
	  treestack[ntreestack++] = newtree;

        // v_name[N] = global value from atom-style variable
	// compute the per-atom variable in result
	// use peratom2global to extract single value from result

	} else if (nbracket && index > 0 && style[ivar] == ATOM) {

	  double *result = (double *)
	    memory->smalloc(atom->nlocal*sizeof(double),"variable:result");
	  compute_atom(ivar,0,result,1,0);
	  peratom2global(1,NULL,result,1,index,
			 tree,treestack,ntreestack,argstack,nargstack);
	  memory->sfree(result);

	} else error->all("Mismatched variable in variable formula");

	delete [] id;

      // ----------------
      // math/group function or atom value/vector or thermo keyword
      // ----------------

      } else {

	// math or group function

	if (str[i] == '(') {
	  char *contents;
	  i = find_matching_paren(str,i,contents);
	  i++;

	  if (math_function(word,contents,tree,
			    treestack,ntreestack,argstack,nargstack));
	  else if (group_function(word,contents,tree,
				  treestack,ntreestack,argstack,nargstack));
	  else error->all("Invalid math or group function in variable formula");

	  delete [] contents;

	// ----------------
	// atom value or vector
	// ----------------

	} else if (str[i] == '[') {
	  int id;
	  i = int_between_brackets(str,i,id,1);
	  i++;

	  if (domain->box_exist == 0)
	    error->all("Variable evaluation before simulation box is defined");
 
	  // ID between brackets exists: atom value
	  // empty brackets: atom vector

	  if (id > 0)
	    peratom2global(0,word,NULL,0,id,
			   tree,treestack,ntreestack,argstack,nargstack);
	  else
	    atom_vector(word,tree,treestack,ntreestack);

	// ----------------
	// thermo keyword
	// ----------------

	} else {
	  if (domain->box_exist == 0)
	    error->all("Variable evaluation before simulation box is defined");
 
	  int flag = output->thermo->evaluate_keyword(word,&value1);
	  if (flag) error->all("Invalid thermo keyword in variable formula");
	  if (tree) {
	    Tree *newtree = new Tree();
	    newtree->type = VALUE;
	    newtree->value = value1;
	    newtree->left = newtree->right = NULL;
	    treestack[ntreestack++] = newtree;
	  } else argstack[nargstack++] = value1;
	}
      }

      delete [] word;

    // ----------------
    // math operator, including end-of-string
    // ----------------

    } else if (strchr("+-*/^\0",onechar)) {
      if (onechar == '+') op = ADD;
      else if (onechar == '-') op = SUBTRACT;
      else if (onechar == '*') op = MULTIPLY;
      else if (onechar == '/') op = DIVIDE;
      else if (onechar == '^') op = CARAT;
      else op = DONE;
      i++;

      if (op == SUBTRACT && expect == ARG) {
	opstack[nopstack++] = UNARY;
	continue;
      }

      if (expect == ARG) error->all("Invalid syntax in variable formula");
      expect = ARG;

      // evaluate stack as deep as possible while respecting precedence
      // before pushing current op onto stack

      while (nopstack && precedence[opstack[nopstack-1]] >= precedence[op]) {
	opprevious = opstack[--nopstack];

	if (tree) {
	  Tree *newtree = new Tree();
	  newtree->type = opprevious;
	  if (opprevious == UNARY) {
	    newtree->left = treestack[--ntreestack];
	    newtree->right = NULL;
	  } else {
	    newtree->right = treestack[--ntreestack];
	    newtree->left = treestack[--ntreestack];
	  }
	  treestack[ntreestack++] = newtree;

	} else {
	  value2 = argstack[--nargstack];
	  if (opprevious != UNARY) value1 = argstack[--nargstack];

	  if (opprevious == ADD)
	    argstack[nargstack++] = value1 + value2;
	  else if (opprevious == SUBTRACT)
	    argstack[nargstack++] = value1 - value2;
	  else if (opprevious == MULTIPLY)
	    argstack[nargstack++] = value1 * value2;
	  else if (opprevious == DIVIDE) {
	    if (value2 == 0.0) error->all("Divide by 0 in variable formula");
	    argstack[nargstack++] = value1 / value2;
	  } else if (opprevious == CARAT) {
	    if (value2 == 0.0) error->all("Power by 0 in variable formula");
	    argstack[nargstack++] = pow(value1,value2);
	  } else if (opprevious == UNARY)
	    argstack[nargstack++] = -value2;
	}
      }

      // if end-of-string, break out of entire formula evaluation loop

      if (op == DONE) break;

      // push current operation onto stack

      opstack[nopstack++] = op;

    } else error->all("Invalid syntax in variable formula");
  }

  if (nopstack) error->all("Invalid syntax in variable formula");

  // for atom-style variable, return remaining tree
  // for equal-style variable, return remaining arg

  if (tree) {
    if (ntreestack != 1) error->all("Invalid syntax in variable formula");
    *tree = treestack[0];
    return 0.0;
  } else {
    if (nargstack != 1) error->all("Invalid syntax in variable formula");
    return argstack[0];
  }
}

/* ----------------------------------------------------------------------
   process an evaulation tree
   customize by adding a math function:
     sqrt(),exp(),ln(),log(),sin(),cos(),tan(),asin(),acos(),atan()
     ceil(),floor(),round()
---------------------------------------------------------------------- */

double Variable::eval_tree(Tree *tree, int i)
{
  if (tree->type == VALUE) return tree->value;
  if (tree->type == ATOMARRAY) return tree->array[i*tree->nstride];
  if (tree->type == TYPEARRAY) return tree->array[atom->type[i]];
  if (tree->type == INTARRAY) return tree->iarray[i*tree->nstride];

  if (tree->type == ADD)
    return eval_tree(tree->left,i) + eval_tree(tree->right,i);
  if (tree->type == SUBTRACT)
    return eval_tree(tree->left,i) - eval_tree(tree->right,i);
  if (tree->type == MULTIPLY)
    return eval_tree(tree->left,i) * eval_tree(tree->right,i);
  if (tree->type == DIVIDE) {
    double denom = eval_tree(tree->right,i);
    if (denom == 0.0) error->all("Divide by 0 in variable formula");
    return eval_tree(tree->left,i) / denom;
  }
  if (tree->type == CARAT) {
    double exponent = eval_tree(tree->right,i);
    if (exponent == 0.0) error->all("Power by 0 in variable formula");
    return pow(eval_tree(tree->left,i),exponent);
  }
  if (tree->type == UNARY)
    return -eval_tree(tree->left,i);

  if (tree->type == SQRT) {
    double arg = eval_tree(tree->left,i);
    if (arg < 0.0) error->all("Sqrt of negative in variable formula");
    return sqrt(arg);
  }
  if (tree->type == EXP)
    return exp(eval_tree(tree->left,i));
  if (tree->type == LN) {
    double arg = eval_tree(tree->left,i);
    if (arg <= 0.0) error->all("Log of zero/negative in variable formula");
    return log(arg);
  }
  if (tree->type == LOG) {
    double arg = eval_tree(tree->left,i);
    if (arg <= 0.0) error->all("Log of zero/negative in variable formula");
    return log10(arg);
  }

  if (tree->type == SIN)
    return sin(eval_tree(tree->left,i));
  if (tree->type == COS)
    return cos(eval_tree(tree->left,i));
  if (tree->type == TAN)
    return tan(eval_tree(tree->left,i));

  if (tree->type == ASIN) {
    double arg = eval_tree(tree->left,i);
    if (arg < -1.0 || arg > 1.0)
      error->all("Arcsin of invalid value in variable formula");
    return asin(arg);
  }
  if (tree->type == ACOS) {
    double arg = eval_tree(tree->left,i);
    if (arg < -1.0 || arg > 1.0)
      error->all("Arccos of invalid value in variable formula");
    return acos(arg);
  }
  if (tree->type == ATAN)
    return atan(eval_tree(tree->left,i));

  if (tree->type == CEIL)
    return ceil(eval_tree(tree->left,i));
  if (tree->type == FLOOR)
    return floor(eval_tree(tree->left,i));
  if (tree->type == ROUND)
    return MYROUND(eval_tree(tree->left,i));

  return 0.0;
}

/* ---------------------------------------------------------------------- */

void Variable::free_tree(Tree *tree)
{
  if (tree->left) free_tree(tree->left);
  if (tree->right) free_tree(tree->right);
  delete tree;
}

/* ----------------------------------------------------------------------
   find matching parenthesis in str, allocate contents = str between parens
   i = left paren
   return loc or right paren
------------------------------------------------------------------------- */

int Variable::find_matching_paren(char *str, int i,char *&contents)
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
  if (!str[i]) error->all("Invalid syntax in variable formula");
  int istop = i;

  int n = istop - istart - 1;
  contents = new char[n+1];
  strncpy(contents,&str[istart+1],n);
  contents[n] = '\0';

  return istop;
}

/* ----------------------------------------------------------------------
   find int between brackets and set index to its value
   if emptyflag, then brackets can be empty (index = 0)
   else if not emptyflag and brackets empty, is an error
   else contents of brackets must be positive int
   i = left bracket
   return loc of right bracket
------------------------------------------------------------------------- */

int Variable::int_between_brackets(char *str, int i, int &index, int emptyflag)
{
  // istop = matching ']'

  int istart = i;
  while (!str[i] || str[i] != ']') i++;
  if (!str[i]) error->all("Invalid syntax in variable formula");
  int istop = i;

  int n = istop - istart - 1;
  char *arg = new char[n+1];
  strncpy(arg,&str[istart+1],n);
  arg[n] = '\0';

  if (n == 0 && emptyflag) index = 0;
  else if (n == 0) error->all("Empty brackets in variable formula");
  else {
    index = atoi(arg);
    if (index <= 0) error->all("Invalid index in variable formula");
  }
  delete [] arg;
 
  return istop;
}

/* ----------------------------------------------------------------------
   process a math function in formula
   push result onto tree or arg stack
   word = math function
   contents = str bewteen parentheses
   return 0 if not a match, 1 if successfully processed
   customize by adding a math function in 2 places:
     sqrt(),exp(),ln(),log(),sin(),cos(),tan(),asin(),acos(),atan()
     ceil(),floor(),round()
------------------------------------------------------------------------- */

int Variable::math_function(char *word, char *contents, Tree **tree,
			    Tree **treestack, int &ntreestack,
			    double *argstack, int &nargstack)
{
  // word not a match to any math function

  if (strcmp(word,"sqrt") && strcmp(word,"exp") && 
      strcmp(word,"ln") && strcmp(word,"log") &&
      strcmp(word,"sin") && strcmp(word,"cos") &&
      strcmp(word,"tan") && strcmp(word,"asin") &&
      strcmp(word,"acos") && strcmp(word,"atan") &&
      strcmp(word,"ceil") && strcmp(word,"floor") && strcmp(word,"round"))
    return 0;
    
  Tree *newtree;
  double value;

  if (tree) {
    newtree = new Tree();
    Tree *argtree;
    double tmp = evaluate(contents,&argtree);
    newtree->left = argtree;
    newtree->right = NULL;
    treestack[ntreestack++] = newtree;
  } else value = evaluate(contents,NULL);
    
  if (strcmp(word,"sqrt") == 0) {
    if (tree) newtree->type = SQRT;
    else {
      if (value < 0.0) error->all("Sqrt of negative in variable formula");
      argstack[nargstack++] = sqrt(value);
    }

  } else if (strcmp(word,"exp") == 0) {
    if (tree) newtree->type = EXP;
    else argstack[nargstack++] = exp(value);
  } else if (strcmp(word,"ln") == 0) {
    if (tree) newtree->type = LN;
    else {
      if (value <= 0.0) error->all("Log of zero/negative in variable formula");
      argstack[nargstack++] = log(value);
    }
  } else if (strcmp(word,"log") == 0) {
    if (tree) newtree->type = LOG;
    else {
      if (value <= 0.0) error->all("Log of zero/negative in variable formula");
      argstack[nargstack++] = log10(value);
    }

  } else if (strcmp(word,"sin") == 0) {
    if (tree) newtree->type = SIN;
    else argstack[nargstack++] = sin(value);
  } else if (strcmp(word,"cos") == 0) {
    if (tree) newtree->type = COS;
    else argstack[nargstack++] = cos(value);
  } else if (strcmp(word,"tan") == 0) {
    if (tree) newtree->type = TAN;
    else argstack[nargstack++] = tan(value);

  } else if (strcmp(word,"asin") == 0) {
    if (tree) newtree->type = ASIN;
    else {
      if (value < -1.0 || value > 1.0) 
	error->all("Arcsin of invalid value in variable formula");
      argstack[nargstack++] = asin(value);
    }
  } else if (strcmp(word,"acos") == 0) {
    if (tree) newtree->type = ACOS;
    else {
      if (value < -1.0 || value > 1.0) 
	error->all("Arccos of invalid value in variable formula");
      argstack[nargstack++] = acos(value);
    }
  } else if (strcmp(word,"atan") == 0) {
    if (tree) newtree->type = ATAN;
    else argstack[nargstack++] = atan(value);

  } else if (strcmp(word,"ceil") == 0) {
    if (tree) newtree->type = CEIL;
    else argstack[nargstack++] = ceil(value);

  } else if (strcmp(word,"floor") == 0) {
    if (tree) newtree->type = FLOOR;
    else argstack[nargstack++] = floor(value);

  } else if (strcmp(word,"round") == 0) {
    if (tree) newtree->type = ROUND;
    else argstack[nargstack++] = MYROUND(value);
  }

  return 1;
}

/* ----------------------------------------------------------------------
   process a group function in formula with optional region arg
   push result onto tree or arg stack
   word = group function
   contents = str bewteen parentheses with one,two,three args
   return 0 if not a match, 1 if successfully processed
   customize by adding a group function:
     count(group),mass(group),charge(group),
     xcm(group,dim),vcm(group,dim),fcm(group,dim),
     bound(group,xmin),gyration(group),ke(group)
------------------------------------------------------------------------- */

int Variable::group_function(char *word, char *contents, Tree **tree,
			     Tree **treestack, int &ntreestack,
			     double *argstack, int &nargstack)
{
  // parse contents for arg1,arg2,arg3 separated by commas
  // ptr1,ptr2 = location of 1st and 2nd comma, NULL if none

  char *arg1,*arg2,*arg3;
  char *ptr1,*ptr2;

  ptr1 = strchr(contents,',');
  if (ptr1) {
    *ptr1 = '\0';
    ptr2 = strchr(ptr1+1,',');
    if (ptr2) *ptr2 = '\0';
  } else ptr2 = NULL;

  int n = strlen(contents) + 1;
  arg1 = new char[n];
  strcpy(arg1,contents);
  int narg = 1;
  if (ptr1) {
    n = strlen(ptr1+1) + 1;
    arg2 = new char[n];
    strcpy(arg2,ptr1+1);
    narg = 2;
  } else arg2 = NULL;
  if (ptr2) {
    n = strlen(ptr2+1) + 1;
    arg3 = new char[n];
    strcpy(arg3,ptr2+1);
    narg = 3;
  } else arg3 = NULL;

  int igroup = group->find(arg1);
  if (igroup == -1)
    error->all("Group ID in variable formula does not exist");

  Tree *newtree;
  double value;
  
  if (tree) {
    newtree = new Tree();
    newtree->type = VALUE;
    newtree->left = newtree->right = NULL;
    treestack[ntreestack++] = newtree;
  }

  // match word to group function

  if (strcmp(word,"count") == 0) {
    if (narg == 1) value = group->count(igroup);
    else if (narg == 2) value = group->count(igroup,region_function(arg2));
    else error->all("Invalid group function in variable formula");

  } else if (strcmp(word,"mass") == 0) {
    if (narg == 1) value = group->mass(igroup);
    else if (narg == 2) value = group->mass(igroup,region_function(arg2));
    else error->all("Invalid group function in variable formula");

  } else if (strcmp(word,"charge") == 0) {
    if (narg == 1) value = group->charge(igroup);
    else if (narg == 2) value = group->charge(igroup,region_function(arg2));
    else error->all("Invalid group function in variable formula");

  } else if (strcmp(word,"xcm") == 0) {
    atom->check_mass();
    double xcm[3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
    } else if (narg == 3) {
      int iregion = region_function(arg3);
      double masstotal = group->mass(igroup,iregion);
      group->xcm(igroup,masstotal,xcm,iregion);
    } else error->all("Invalid group function in variable formula");
    if (strcmp(arg2,"x") == 0) value = xcm[0];
    else if (strcmp(arg2,"y") == 0) value = xcm[1];
    else if (strcmp(arg2,"z") == 0) value = xcm[2];
    else error->all("Invalid group function in variable formula");

  } else if (strcmp(word,"vcm") == 0) {
    atom->check_mass();
    double vcm[3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->vcm(igroup,masstotal,vcm);
    } else if (narg == 3) {
      int iregion = region_function(arg3);
      double masstotal = group->mass(igroup,iregion);
      group->vcm(igroup,masstotal,vcm,iregion);
    } else error->all("Invalid group function in variable formula");
    if (strcmp(arg2,"x") == 0) value = vcm[0];
    else if (strcmp(arg2,"y") == 0) value = vcm[1];
    else if (strcmp(arg2,"z") == 0) value = vcm[2];
    else error->all("Invalid group function in variable formula");

  } else if (strcmp(word,"fcm") == 0) {
    double fcm[3];
    if (narg == 2) group->fcm(igroup,fcm);
    else if (narg == 3) group->fcm(igroup,fcm,region_function(arg3));
    else error->all("Invalid group function in variable formula");
    if (strcmp(arg2,"x") == 0) value = fcm[0];
    else if (strcmp(arg2,"y") == 0) value = fcm[1];
    else if (strcmp(arg2,"z") == 0) value = fcm[2];
    else error->all("Invalid group function in variable formula");

  } else if (strcmp(word,"bound") == 0) {
    double minmax[6];
    if (narg == 2) group->bounds(igroup,minmax);
    else if (narg == 3) group->bounds(igroup,minmax,region_function(arg3));
    else error->all("Invalid group function in variable formula");
    if (strcmp(arg2,"xmin") == 0) value = minmax[0];
    else if (strcmp(arg2,"xmax") == 0) value = minmax[1];
    else if (strcmp(arg2,"ymin") == 0) value = minmax[2];
    else if (strcmp(arg2,"ymax") == 0) value = minmax[3];
    else if (strcmp(arg2,"zmin") == 0) value = minmax[4];
    else if (strcmp(arg2,"zmax") == 0) value = minmax[5];
    else error->all("Invalid group function in variable formula");

  } else if (strcmp(word,"gyration") == 0) {
    atom->check_mass();
    double xcm[3];
    if (narg == 1) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
      value = group->gyration(igroup,masstotal,xcm);
    } else if (narg == 2) {
      int iregion = region_function(arg2);
      double masstotal = group->mass(igroup,iregion);
      group->xcm(igroup,masstotal,xcm,iregion);
      value = group->gyration(igroup,masstotal,xcm,iregion);
    } else error->all("Invalid group function in variable formula");

  } else if (strcmp(word,"ke") == 0) {
    if (narg == 1) value = group->ke(igroup);
    else if (narg == 2) value = group->ke(igroup,region_function(arg2));
    else error->all("Invalid group function in variable formula");

  } else return 0;
    
  delete [] arg1;
  delete [] arg2;
    
  if (tree) newtree->value= value;
  else argstack[nargstack++] = value;

  return 1;
}

/* ----------------------------------------------------------------------
   process a group function in formula with optional region arg
   push result onto tree or arg stack
   word = group function
   contents = str bewteen parentheses with one,two,three args
   return 0 if not a match, 1 if successfully processed
   customize by adding a group function:
     count(group),mass(group),charge(group),
     xcm(group,dim),vcm(group,dim),fcm(group,dim),
     bound(group,xmin),gyration(group)
------------------------------------------------------------------------- */

int Variable::region_function(char *id)
{
  int iregion = domain->find_region(id);
  if (iregion == -1)
    error->all("Invalid region in group function in variable formula");
  return iregion;
}

/* ----------------------------------------------------------------------
   extract a global value from a per-atom quantity in a formula
   flag = 0 -> word is an atom vector
   flag = 1 -> vector is a per-atom compute or fix quantity with nstride
   id = positive global ID of atom, converted to local index
   push result onto tree or arg stack
   customize by adding an atom vector: mass,type,x,y,z,vx,vy,vz,fx,fy,fz
------------------------------------------------------------------------- */

void Variable::peratom2global(int flag, char *word,
			      double *vector, int nstride, int id,
			      Tree **tree, Tree **treestack, int &ntreestack,
			      double *argstack, int &nargstack)
{
  if (atom->map_style == 0)
    error->all("Indexed per-atom vector in variable formula without atom map");

  int index = atom->map(id);

  double mine;
  if (index >= 0 && index < atom->nlocal) {

    if (flag == 0) {
      if (strcmp(word,"mass") == 0) {
	if (atom->rmass) mine = atom->rmass[index];
	else mine = atom->mass[atom->type[index]];
      }
      else if (strcmp(word,"type") == 0) mine = atom->type[index];
      else if (strcmp(word,"x") == 0) mine = atom->x[index][0];
      else if (strcmp(word,"y") == 0) mine = atom->x[index][1];
      else if (strcmp(word,"z") == 0) mine = atom->x[index][2];
      else if (strcmp(word,"vx") == 0) mine = atom->v[index][0];
      else if (strcmp(word,"vy") == 0) mine = atom->v[index][1];
      else if (strcmp(word,"vz") == 0) mine = atom->v[index][2];
      else if (strcmp(word,"fx") == 0) mine = atom->f[index][0];
      else if (strcmp(word,"fy") == 0) mine = atom->f[index][1];
      else if (strcmp(word,"fz") == 0) mine = atom->f[index][2];
      
      else error->one("Invalid atom vector in variable formula");

    } else mine = vector[index*nstride];
    
  } else mine = 0.0;

  double value;
  MPI_Allreduce(&mine,&value,1,MPI_DOUBLE,MPI_SUM,world);
  
  if (tree) {
    Tree *newtree = new Tree();
    newtree->type = VALUE;
    newtree->value = value;
    newtree->left = newtree->right = NULL;
    treestack[ntreestack++] = newtree;
  } else argstack[nargstack++] = value;
}

/* ----------------------------------------------------------------------
   process an atom vector in formula
   push result onto tree
   word = atom vector
   customize by adding an atom vector: mass,type,x,y,z,vx,vy,vz,fx,fy,fz
------------------------------------------------------------------------- */

void Variable::atom_vector(char *word, Tree **tree,
			   Tree **treestack, int &ntreestack)
{
  if (tree == NULL)
    error->all("Atom vector in equal-style variable formula");

  Tree *newtree = new Tree();
  newtree->type = ATOMARRAY;
  newtree->nstride = 3;
  newtree->left = newtree->right = NULL;
  treestack[ntreestack++] = newtree;
	    
  if (strcmp(word,"mass") == 0) {
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
  
  else error->all("Invalid atom vector in variable formula");
}
