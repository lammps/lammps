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
#include "output.h"
#include "thermo.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define VARDELTA 4

enum{INDEX,LOOP,EQUAL,WORLD,UNIVERSE,ULOOP,ATOM};
enum{VALUE,ATOMARRAY,TYPEARRAY,ADD,SUB,MULT,DIV,NEG,POW,EXP,LN,SQRT};

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

  // if var already exists, just skip, except EQUAL vars

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

  // UNIVERSE and ULOOP
  // for UNIVERSE: num = listed args, data = copied args
  // for ULOOP: num = N, data = list of NULLS since never used
  // index = partition this proc is in
  // universe proc 0 creates lock file
  // error check that all other universe/uloop variables are same length

  } else if (strcmp(arg[1],"universe") == 0 || strcmp(arg[1],"uloop") == 0) {
    if (strcmp(arg[1],"universe") == 0) {
      style[nvar] = UNIVERSE;
      num[nvar] = narg - 2;
      data[nvar] = new char*[num[nvar]];
      copy(num[nvar],&arg[2],data[nvar]);
    } else {
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
  // num = 1, index = 1st value
  // data = 1 value, string to eval

  } else if (strcmp(arg[1],"atom") == 0) {
    style[nvar] = ATOM;
    num[nvar] = 1;
    index[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(1,&arg[2],data[nvar]);
    
  } else error->all("Illegal variable command");

  // set name of variable
  // must come at end, since EQUAL reset may have removed name

  names[nvar] = new char[strlen(arg[0]) + 1];
  strcpy(names[nvar],arg[0]);
  nvar++;
}

/* ----------------------------------------------------------------------
   single-value INDEX variable created by command-line argument
------------------------------------------------------------------------- */

void Variable::set(char *name, char *value)
{
  char **newarg = new char*[3];
  newarg[0] = name;
  newarg[1] = "index";
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
   recursive evaluation of a "string"
   string is defined as one of several items:
     number =  0.0, -5.45, 2.8e-4, ...
     thermo keyword = ke, vol, atoms, ...
     variable keyword = v_a, v_myvar, ...
     math function = div(x,y), mult(x,y), add(x,y), ...
     group function = mass(group), xcm(group,x), ...
     atom vector = x[123], y[3], vx[34], ...
     compute vector = c_mytemp[0], c_thermo_press[3], ...
   numbers start with a digit or "." or "-" (no parens or brackets)
   keywords start with a lowercase letter (no parens or brackets)
   functions contain ()
     can have 1 or 2 args, each of which can be a "string"
   vectors contain []
     single arg must be integer
     for atom vectors, arg is global ID of atom
     for compute vectors, 0 is the scalar value, 1-N are vector values
   see lists of valid functions & vectors below
   return answer = value of string
------------------------------------------------------------------------- */

double Variable::evaluate(char *str, Tree *tree)
{
  double answer = 0.0;
  if (tree) {
    tree->type = VALUE;
    tree->left = tree->right = NULL;
  }

  // string is a "function" since contains ()
  // grab one or two args, separated by ","
  // evaulate args recursively
  // if tree = NULL, evaluate math or group function
  // else store as leaf in tree

  if (strchr(str,'(')) {
    if (str[strlen(str)-1] != ')') error->all("Cannot evaluate variable");

    char *ptr = strchr(str,'(');
    int n = ptr - str;
    char *func = new char[n+1];
    strncpy(func,str,n);
    func[n] = '\0';

    char *comma = ++ptr;
    int level = 0;
    while (1) {
      if (*comma == '\0') error->all("Cannot evaluate variable");
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

    // customize by adding math function to this list and to if statement
    //   add(x,y),sub(x,y),mult(x,y),div(x,y),neg(x),
    //   pow(x,y),exp(x),ln(x),sqrt(x)

    if (strcmp(func,"add") == 0) {
      if (!arg2) error->all("Cannot evaluate variable");
      if (tree) {
	tree->type = ADD;
	tree->left = new Tree();
	tree->right = new Tree();
	value1 = evaluate(arg1,tree->left);
	value2 = evaluate(arg2,tree->right);
      } else answer = evaluate(arg1,NULL) + evaluate(arg2,NULL);

    } else if (strcmp(func,"sub") == 0) {
      if (!arg2) error->all("Cannot evaluate variable");
      if (tree) {
	tree->type = SUB;
	tree->left = new Tree();
	tree->right = new Tree();
	value1 = evaluate(arg1,tree->left);
	value2 = evaluate(arg2,tree->right);
      } else answer = evaluate(arg1,NULL) - evaluate(arg2,NULL);

    } else if (strcmp(func,"mult") == 0) {
      if (!arg2) error->all("Cannot evaluate variable");
      if (tree) {
	tree->type = MULT;
	tree->left = new Tree();
	tree->right = new Tree();
	value1 = evaluate(arg1,tree->left);
	value2 = evaluate(arg2,tree->right);
      } else answer = evaluate(arg1,NULL) * evaluate(arg2,NULL);

    } else if (strcmp(func,"div") == 0) {
      if (!arg2) error->all("Cannot evaluate variable");
      if (tree) {
	tree->type = DIV;
	tree->left = new Tree();
	tree->right = new Tree();
	value1 = evaluate(arg1,tree->left);
	value2 = evaluate(arg2,tree->right);
      } else {
	value1 = evaluate(arg1,NULL);
	value2 = evaluate(arg2,NULL);
	if (value2 == 0.0) error->all("Cannot evaluate variable");
	answer = value1 / value2;
      }

    } else if (strcmp(func,"neg") == 0) {
      if (arg2) error->all("Cannot evaluate variable");
      if (tree) {
	tree->type = NEG;
	tree->left = new Tree();
	value1 = evaluate(arg1,tree->left);
      } else answer = -evaluate(arg1,NULL);

    } else if (strcmp(func,"pow") == 0) {
      if (!arg2) error->all("Cannot evaluate variable");
      if (tree) {
	tree->type = POW;
	tree->left = new Tree();
	tree->right = new Tree();
	value1 = evaluate(arg1,tree->left);
	value2 = evaluate(arg2,tree->right);
      } else {
	value1 = evaluate(arg1,NULL);
	value2 = evaluate(arg2,NULL);
	if (value2 == 0.0) error->all("Cannot evaluate variable");
	answer = pow(value1,value2);
      }

    } else if (strcmp(func,"exp") == 0) {
      if (arg2) error->all("Cannot evaluate variable");
      if (tree) {
	tree->type = EXP;
	tree->left = new Tree();
	value1 = evaluate(arg1,tree->left);
      } else answer = exp(evaluate(arg1,NULL));

    } else if (strcmp(func,"ln") == 0) {
      if (arg2) error->all("Cannot evaluate variable");
      if (tree) {
	tree->type = LN;
	tree->left = new Tree();
	value1 = evaluate(arg1,tree->left);
      } else {
	value1 = evaluate(arg1,NULL);
	if (value1 == 0.0) error->all("Cannot evaluate variable");
	answer = log(value1);
      }

    } else if (strcmp(func,"sqrt") == 0) {
      if (arg2) error->all("Cannot evaluate variable");
      if (tree) {
	tree->type = SQRT;
	tree->left = new Tree();
	value1 = evaluate(arg1,tree->left);
      } else {
	value1 = evaluate(arg1,NULL);
	if (value1 == 0.0) error->all("Cannot evaluate variable");
	answer = sqrt(value1);
      }

    // customize by adding group function to this list and to if statement
    //   mass(group),charge(group),xcm(group,dim),vcm(group,dim),
    //   fcm(group,dim),bound(group,xmin),gyration(group)

    } else if (strcmp(func,"mass") == 0) {
      if (arg2) error->all("Cannot evaluate variable");
      int igroup = group->find(arg1);
      if (igroup == -1) error->all("Variable group ID does not exist");
      atom->check_mass();
      answer = group->mass(igroup);

    } else if (strcmp(func,"charge") == 0) {
      if (arg2) error->all("Cannot evaluate variable");
      int igroup = group->find(arg1);
      if (igroup == -1) error->all("Variable group ID does not exist");
      answer = group->charge(igroup);

    } else if (strcmp(func,"xcm") == 0) {
      if (!arg2) error->all("Cannot evaluate variable");
      int igroup = group->find(arg1);
      if (igroup == -1) error->all("Variable group ID does not exist");
      atom->check_mass();
      double masstotal = group->mass(igroup);
      double xcm[3];
      group->xcm(igroup,masstotal,xcm);
      if (strcmp(arg2,"x") == 0) answer = xcm[0];
      else if (strcmp(arg2,"y") == 0) answer = xcm[1];
      else if (strcmp(arg2,"z") == 0) answer = xcm[2];
      else error->all("Cannot evaluate variable");

    } else if (strcmp(func,"vcm") == 0) {
      if (!arg2) error->all("Cannot evaluate variable");
      int igroup = group->find(arg1);
      if (igroup == -1) error->all("Variable group ID does not exist");
      atom->check_mass();
      double masstotal = group->mass(igroup);
      double vcm[3];
      group->vcm(igroup,masstotal,vcm);
      if (strcmp(arg2,"x") == 0) answer = vcm[0];
      else if (strcmp(arg2,"y") == 0) answer = vcm[1];
      else if (strcmp(arg2,"z") == 0) answer = vcm[2];
      else error->all("Cannot evaluate variable");

    } else if (strcmp(func,"fcm") == 0) {
      if (!arg2) error->all("Cannot evaluate variable");
      int igroup = group->find(arg1);
      if (igroup == -1) error->all("Variable group ID does not exist");
      double fcm[3];
      group->fcm(igroup,fcm);
      if (strcmp(arg2,"x") == 0) answer = fcm[0];
      else if (strcmp(arg2,"y") == 0) answer = fcm[1];
      else if (strcmp(arg2,"z") == 0) answer = fcm[2];
      else error->all("Cannot evaluate variable");

    } else if (strcmp(func,"bound") == 0) {
      if (!arg2) error->all("Cannot evaluate variable");
      int igroup = group->find(arg1);
      if (igroup == -1) error->all("Variable group ID does not exist");
      double minmax[6];
      group->bounds(igroup,minmax);
      if (strcmp(arg2,"xmin") == 0) answer = minmax[0];
      else if (strcmp(arg2,"xmax") == 0) answer = minmax[1];
      else if (strcmp(arg2,"ymin") == 0) answer = minmax[2];
      else if (strcmp(arg2,"ymax") == 0) answer = minmax[3];
      else if (strcmp(arg2,"zmin") == 0) answer = minmax[4];
      else if (strcmp(arg2,"zmax") == 0) answer = minmax[5];
      else error->all("Cannot evaluate variable");

    } else if (strcmp(func,"gyration") == 0) {
      if (!arg2) error->all("Cannot evaluate variable");
      int igroup = group->find(arg1);
      if (igroup == -1) error->all("Variable group ID does not exist");
      atom->check_mass();
      double masstotal = group->mass(igroup);
      double xcm[3];
      group->xcm(igroup,masstotal,xcm);
      answer = group->gyration(igroup,masstotal,xcm);

    } else error->all("Invalid math/group function in variable");

    delete [] func;
    delete [] arg1;
    delete [] arg2;

  // string is a "vector" since contains []
  // index = everything between [] evaluated as integer
  // if vector name starts with "c_", trailing chars are compute ID
  //   check if compute ID exists, invoke it with index as arg
  // else is atom vector
  //   find which proc owns atom via atom->map()
  //   grab atom-based value with index as global atom ID

  } else if (strchr(str,'[')) {
    if (str[strlen(str)-1] != ']') error->all("Cannot evaluate variable");

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

    if (strncmp(vector,"c_",2) == 0) {
      n = strlen(vector) - 2 + 1;
      char *id = new char[n];
      strcpy(id,&vector[2]);

      int icompute;
      for (icompute = 0; icompute < modify->ncompute; icompute++)
	if (strcmp(id,modify->compute[icompute]->id) == 0) break;
      if (icompute == modify->ncompute)
	error->all("Invalid compute ID in variable");
      delete [] id;
      modify->compute[icompute]->init();

      // call compute() if index = 0, else compute_vector()
      // make pre-call to Compute object's id_pre if it is defined

      int index = atoi(arg);
      if (index == 0) {
	if (modify->compute[icompute]->scalar_flag == 0)
	  error->all("Variable compute ID does not compute scalar info");
	if (modify->compute[icompute]->id_pre) {
	  int ipre = modify->find_compute(modify->compute[icompute]->id_pre);
	  if (ipre < 0) error->all("Could not pre-compute in variable");
	  answer = modify->compute[ipre]->compute_scalar();
	}
	answer = modify->compute[icompute]->compute_scalar();
      } else if (index > 0) {
	if (modify->compute[icompute]->vector_flag == 0)
	  error->all("Variable compute ID does not compute scalar info");
	if (index > modify->compute[icompute]->size_vector)
	  error->all("Variable compute ID vector is not large enough");
	if (modify->compute[icompute]->id_pre) {
	  int ipre = modify->find_compute(modify->compute[icompute]->id_pre);
	  if (ipre < 0) error->all("Could not pre-compute in variable");
	  modify->compute[ipre]->compute_vector();
	}
	modify->compute[icompute]->compute_vector();
	answer = modify->compute[icompute]->vector[index-1];
      } else error->all("Invalid compute ID index in variable");

    } else if (tree) {

      if (strlen(arg)) error->all("Invalid atom vector in variable");

      // customize by adding atom vector to this list and to if statement
      // mass,x,y,z,vx,vy,vz,fx,fy,fz

      tree->type = ATOMARRAY;
      tree->nstride = 3;

      if (strcmp(vector,"mass") == 0) {
	if (atom->mass) {
	  tree->type = TYPEARRAY;
	  tree->array = atom->mass;
	} else {
	  tree->nstride = 1;
	  tree->array = atom->rmass;
	}
      }
      else if (strcmp(vector,"x") == 0) tree->array = &atom->x[0][0];
      else if (strcmp(vector,"y") == 0) tree->array = &atom->x[0][1];
      else if (strcmp(vector,"z") == 0) tree->array = &atom->x[0][2];
      else if (strcmp(vector,"vx") == 0) tree->array = &atom->v[0][0];
      else if (strcmp(vector,"vy") == 0) tree->array = &atom->v[0][1];
      else if (strcmp(vector,"vz") == 0) tree->array = &atom->v[0][2];
      else if (strcmp(vector,"fx") == 0) tree->array = &atom->f[0][0];
      else if (strcmp(vector,"fy") == 0) tree->array = &atom->f[0][1];
      else if (strcmp(vector,"fz") == 0) tree->array = &atom->f[0][2];
      
      else error->all("Invalid atom vector in variable");

    } else {

      if (atom->map_style == 0)
	error->all("Cannot use atom vector in variable unless atom map exists");

      int index = atom->map(atoi(arg));

      // customize by adding atom vector to this list and to if statement
      // mass,x,y,z,vx,vy,vz,fx,fy,fz

      double mine;
      if (index >= 0 && index < atom->nlocal) {

	if (strcmp(vector,"mass") == 0) {
	  if (atom->mass) mine = atom->mass[atom->type[index]];
	  else mine = atom->rmass[index];
	}
        else if (strcmp(vector,"x") == 0) mine = atom->x[index][0];
	else if (strcmp(vector,"y") == 0) mine = atom->x[index][1];
	else if (strcmp(vector,"z") == 0) mine = atom->x[index][2];
	else if (strcmp(vector,"vx") == 0) mine = atom->v[index][0];
	else if (strcmp(vector,"vy") == 0) mine = atom->v[index][1];
	else if (strcmp(vector,"vz") == 0) mine = atom->v[index][2];
	else if (strcmp(vector,"fx") == 0) mine = atom->f[index][0];
	else if (strcmp(vector,"fy") == 0) mine = atom->f[index][1];
	else if (strcmp(vector,"fz") == 0) mine = atom->f[index][2];
	
	else error->one("Invalid atom vector in variable");
	
      } else mine = 0.0;

      MPI_Allreduce(&mine,&answer,1,MPI_DOUBLE,MPI_SUM,world);
    }

    delete [] vector;
    delete [] arg;

  // string is "keyword" since starts with lowercase letter
  // if keyword starts with "v_", trailing chars are variable name
  //   evaluate it via retrieve(), convert it to double
  // else is thermo keyword
  //   evaluate it via evaluate_keyword()

  } else if (islower(str[0])) {

    if (strncmp(str,"v_",2) == 0) {
      int n = strlen(str) - 2 + 1;
      char *id = new char[n];
      strcpy(id,&str[2]);

      char *v = retrieve(id);
      if (v == NULL) error->all("Invalid variable name in variable");
      delete [] id;

      answer = atof(v);

    } else {
      int flag = output->thermo->evaluate_keyword(str,&answer);
      if (flag) error->all("Invalid thermo keyword in variable");
    }

  // string is a number since starts with digit or "." or "-"
  // just copy to result

  } else if (isdigit(str[0]) || str[0] == '.' || str[0] == '-') {

    answer = atof(str);

  // string is an error

  } else error->all("Cannot evaluate variable");

  // store answer in tree and return it

  if (tree) tree->value = answer;
  return answer;
}

/* ---------------------------------------------------------------------- */

void Variable::build_parse_tree(int ivar)
{
  if (style[ivar] != ATOM)
    error->all("Cannot build parse tree for variable that is not atom style");
  ptree = new Tree();
  double tmp = evaluate(data[ivar][0],ptree);
}

/* ---------------------------------------------------------------------- */

void Variable::evaluate_parse_tree(int igroup, double *result)
{
  int groupbit = group->bitmask[igroup];
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] && groupbit) result[i] = eval_tree(ptree,i);
    else result[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void Variable::free_parse_tree()
{
  free_tree(ptree);
}

/* ---------------------------------------------------------------------- */

double Variable::eval_tree(Tree *tree, int i)
{
  if (tree->type == VALUE) return tree->value;
  if (tree->type == ATOMARRAY) return tree->array[i*tree->nstride];
  if (tree->type == TYPEARRAY) return tree->array[atom->type[i]];

  if (tree->type == ADD)
    return eval_tree(tree->left,i) + eval_tree(tree->right,i);
  if (tree->type == SUB)
    return eval_tree(tree->left,i) - eval_tree(tree->right,i);
  if (tree->type == MULT)
    return eval_tree(tree->left,i) * eval_tree(tree->right,i);
  if (tree->type == DIV)
    return eval_tree(tree->left,i) / eval_tree(tree->right,i);
  if (tree->type == NEG)
    return -eval_tree(tree->left,i);
  if (tree->type == POW)
    return pow(eval_tree(tree->left,i),eval_tree(tree->right,i));
  if (tree->type == EXP)
    return exp(eval_tree(tree->left,i));
  if (tree->type == LN)
    return log(eval_tree(tree->left,i));
  if (tree->type == SQRT)
    return sqrt(eval_tree(tree->left,i));

  return 0.0;
}

/* ---------------------------------------------------------------------- */

void Variable::free_tree(Tree *tree)
{
  if (tree->left) free_tree(tree->left);
  if (tree->right) free_tree(tree->right);
  delete tree;
}
