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

#include "stdlib.h"
#include "string.h"
#include "dump_custom.h"
#include "atom.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_centro.h"
#include "fix_energy.h"
#include "fix_stress.h"
#include "memory.h"
#include "error.h"

// customize by adding to 1st enum

enum{TAG,MOL,TYPE,X,Y,Z,XS,YS,ZS,XU,YU,ZU,IX,IY,IZ,VX,VY,VZ,FX,FY,FZ,
       Q,MUX,MUY,MUZ,TQX,TQY,TQZ,CENTRO,ENG,SXX,SYY,SZZ,SXY,SXZ,SYZ};
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE};

/* ---------------------------------------------------------------------- */

DumpCustom::DumpCustom(int narg, char **arg) : Dump(narg, arg)
{
  if (narg == 5) error->all("No dump custom arguments specified");

  nevery = atoi(arg[3]);

  size_one = narg-5;
  pack_choice = new FnPtrPack[size_one];
  vtype = new int[size_one];

  iregion = -1;
  nthresh = 0;
  thresh_array = NULL;
  thresh_op = NULL;
  thresh_value = NULL;

  fixflag_centro = 0;
  fixflag_energy = 0;
  fixflag_stress = 0;

  // customize by adding to if statement

  int i;
  for (int iarg = 5; iarg < narg; iarg++) {
    i = iarg-5;
    if (strcmp(arg[iarg],"tag") == 0) {
      pack_choice[i] = &DumpCustom::pack_tag;
      vtype[i] = INT;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (atom->molecule == NULL)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_molecule;
      vtype[i] = INT;
    } else if (strcmp(arg[iarg],"type") == 0) {
      pack_choice[i] = &DumpCustom::pack_type;
      vtype[i] = INT;
    } else if (strcmp(arg[iarg],"x") == 0) {
      pack_choice[i] = &DumpCustom::pack_x;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"y") == 0) {
      pack_choice[i] = &DumpCustom::pack_y;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"z") == 0) {
      pack_choice[i] = &DumpCustom::pack_z;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xs") == 0) {
      pack_choice[i] = &DumpCustom::pack_xs;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"ys") == 0) {
      pack_choice[i] = &DumpCustom::pack_ys;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"zs") == 0) {
      pack_choice[i] = &DumpCustom::pack_zs;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xu") == 0) {
      pack_choice[i] = &DumpCustom::pack_xu;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"yu") == 0) {
      pack_choice[i] = &DumpCustom::pack_yu;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"zu") == 0) {
      pack_choice[i] = &DumpCustom::pack_zu;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"ix") == 0) {
      pack_choice[i] = &DumpCustom::pack_ix;
      vtype[i] = INT;
    } else if (strcmp(arg[iarg],"iy") == 0) {
      pack_choice[i] = &DumpCustom::pack_iy;
      vtype[i] = INT;
    } else if (strcmp(arg[iarg],"iz") == 0) {
      pack_choice[i] = &DumpCustom::pack_iz;
      vtype[i] = INT;
    } else if (strcmp(arg[iarg],"vx") == 0) {
      pack_choice[i] = &DumpCustom::pack_vx;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      pack_choice[i] = &DumpCustom::pack_vy;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      pack_choice[i] = &DumpCustom::pack_vz;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"fx") == 0) {
      pack_choice[i] = &DumpCustom::pack_fx;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"fy") == 0) {
      pack_choice[i] = &DumpCustom::pack_fy;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"fz") == 0) {
      pack_choice[i] = &DumpCustom::pack_fz;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"q") == 0) {
      if (atom->q == NULL)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_q;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"mux") == 0) {
      if (atom->mu == NULL)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_mux;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"muy") == 0) {
      if (atom->mu == NULL)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_muy;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"muz") == 0) {
      if (atom->mu == NULL)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_muz;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"tqx") == 0) {
      if (atom->torque == NULL)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_tqx;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"tqy") == 0) {
      if (atom->torque == NULL)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_tqy;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"tqz") == 0) {
      if (atom->torque == NULL)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_tqz;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"centro") == 0) {
      pack_choice[i] = &DumpCustom::pack_centro;
      vtype[i] = DOUBLE;
      fixflag_centro = 1;
    } else if (strcmp(arg[iarg],"eng") == 0) {
      pack_choice[i] = &DumpCustom::pack_energy;
      vtype[i] = DOUBLE;
      fixflag_energy = 1;
    } else if (strcmp(arg[iarg],"sxx") == 0) {
      pack_choice[i] = &DumpCustom::pack_sxx;
      vtype[i] = DOUBLE;
      fixflag_stress = 1;
    } else if (strcmp(arg[iarg],"syy") == 0) {
      pack_choice[i] = &DumpCustom::pack_syy;
      vtype[i] = DOUBLE;
      fixflag_stress = 1;
    } else if (strcmp(arg[iarg],"szz") == 0) {
      pack_choice[i] = &DumpCustom::pack_szz;
      vtype[i] = DOUBLE;
      fixflag_stress = 1;
    } else if (strcmp(arg[iarg],"sxy") == 0) {
      pack_choice[i] = &DumpCustom::pack_sxy;
      vtype[i] = DOUBLE;
      fixflag_stress = 1;
    } else if (strcmp(arg[iarg],"sxz") == 0) {
      pack_choice[i] = &DumpCustom::pack_sxz;
      vtype[i] = DOUBLE;
      fixflag_stress = 1;
    } else if (strcmp(arg[iarg],"syz") == 0) {
      pack_choice[i] = &DumpCustom::pack_syz;
      vtype[i] = DOUBLE;
      fixflag_stress = 1;
    } else error->all("Invalid keyword in dump custom command");
  }

  if (fixflag_centro) fix_create("CENTRO");
  if (fixflag_energy) fix_create("ENERGY");
  if (fixflag_stress) fix_create("STRESS");

  vformat = new char*[size_one];

  format_default = new char[3*size_one+1];
  format_default[0] = '\0';

  for (i = 0; i < size_one; i++) {
    if (vtype[i] == INT) format_default = strcat(format_default,"%d ");
    else format_default = strcat(format_default,"%g ");
    vformat[i] = NULL;
  }

  maxchoose = 0;
  choose = NULL;
  dchoose = NULL;

  // one-time file open

  if (multifile == 0) openfile();
}

/* ----------------------------------------------------------------------
   free memory
------------------------------------------------------------------------- */

DumpCustom::~DumpCustom()
{
  delete [] vtype;
  for (int i = 0; i < size_one; i++) delete [] vformat[i];
  delete [] vformat;
  delete [] pack_choice;

  if (nthresh) {
    memory->sfree(thresh_array);
    memory->sfree(thresh_op);
    memory->sfree(thresh_value);
  }

  // delete fixes for computing per-atom custom quantities if they exist

  if (fixflag_centro) fix_delete("CENTRO");
  if (fixflag_energy) fix_delete("ENERGY");
  if (fixflag_stress) fix_delete("STRESS");

  if (maxchoose) {
    delete [] choose;
    delete [] dchoose;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::init()
{
  delete [] format;
  char *str;
  if (format_user) str = format_user;
  else str = format_default;

  int n = strlen(str) + 1;
  format = new char[n];
  strcpy(format,str);

  // tokenize the format string and add space at end of each format element

  char *ptr;
  for (int i = 0; i < size_one; i++) {
    if (i == 0) ptr = strtok(format," \0");
    else ptr = strtok(NULL," \0");
    delete [] vformat[i];
    vformat[i] = new char[strlen(ptr) + 2];
    strcpy(vformat[i],ptr);
    vformat[i] = strcat(vformat[i]," ");
  }

  // setup function ptrs

  if (binary) header_choice = &DumpCustom::header_binary;
  else header_choice = &DumpCustom::header_item;

  if (binary) write_choice = &DumpCustom::write_binary;
  else write_choice = &DumpCustom::write_text;

  // set fixflag if dump invokes any fixes

  fixflag = 0;
  if (fixflag_centro || fixflag_energy || fixflag_stress) fixflag = 1;

  // find associated fixes that must exist
  // could have changed locations in fix list since created

  if (fixflag_centro) ifix_centro = fix_match("CENTRO");
  if (fixflag_energy) ifix_energy = fix_match("ENERGY");
  if (fixflag_stress) ifix_stress = fix_match("STRESS");

  // allocate choose arrays if first time

  if (maxchoose == 0) {
    maxchoose = atom->nmax;
    choose = new int[maxchoose];
    dchoose = new double[maxchoose];
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf and choose arrays
------------------------------------------------------------------------- */

int DumpCustom::memory_usage()
{
  int bytes = maxbuf * sizeof(double);
  bytes += maxchoose * sizeof(int);
  bytes += maxchoose * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

int DumpCustom::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"region") == 0) {
    if (narg < 2) error->all("Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) iregion = -1;
    else {
      for (iregion = 0; iregion < domain->nregion; iregion++)
	if (strcmp(arg[1],domain->regions[iregion]->id) == 0) break;
      if (iregion == domain->nregion)
	error->all("Dump_modify region ID does not exist");
    }
    return 2;

  } else if (strcmp(arg[0],"thresh") == 0) {
    if (narg < 2) error->all("Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) {
      if (nthresh) {
	memory->sfree(thresh_array);
	memory->sfree(thresh_op);
	memory->sfree(thresh_value);
	thresh_array = NULL;
	thresh_op = NULL;
	thresh_value = NULL;
      }
      nthresh = 0;
      return 2;
    }
    
    if (narg < 4) error->all("Illegal dump_modify command");
    thresh_array = (int *)
      memory->srealloc(thresh_array,(nthresh+1)*sizeof(int),
		       "dump:thresh_array");
    thresh_op = (int *)
      memory->srealloc(thresh_op,(nthresh+1)*sizeof(int),
		       "dump:thresh_op");
    thresh_value = (double *)
      memory->srealloc(thresh_value,(nthresh+1)*sizeof(double),
		       "dump:thresh_value");

    // customize by adding to if statement

    if (strcmp(arg[1],"tag") == 0) thresh_array[nthresh] = TAG;
    else if (strcmp(arg[1],"mol") == 0) thresh_array[nthresh] = MOL;
    else if (strcmp(arg[1],"type") == 0) thresh_array[nthresh] = TYPE;
    else if (strcmp(arg[1],"x") == 0) thresh_array[nthresh] = X;
    else if (strcmp(arg[1],"y") == 0) thresh_array[nthresh] = Y;
    else if (strcmp(arg[1],"z") == 0) thresh_array[nthresh] = Z;
    else if (strcmp(arg[1],"xs") == 0) thresh_array[nthresh] = XS;
    else if (strcmp(arg[1],"ys") == 0) thresh_array[nthresh] = YS;
    else if (strcmp(arg[1],"zs") == 0) thresh_array[nthresh] = ZS;
    else if (strcmp(arg[1],"xu") == 0) thresh_array[nthresh] = XU;
    else if (strcmp(arg[1],"yu") == 0) thresh_array[nthresh] = YU;
    else if (strcmp(arg[1],"zu") == 0) thresh_array[nthresh] = ZU;
    else if (strcmp(arg[1],"ix") == 0) thresh_array[nthresh] = IX;
    else if (strcmp(arg[1],"iy") == 0) thresh_array[nthresh] = IY;
    else if (strcmp(arg[1],"iz") == 0) thresh_array[nthresh] = IZ;
    else if (strcmp(arg[1],"vx") == 0) thresh_array[nthresh] = VX;
    else if (strcmp(arg[1],"vy") == 0) thresh_array[nthresh] = VY;
    else if (strcmp(arg[1],"vz") == 0) thresh_array[nthresh] = VZ;
    else if (strcmp(arg[1],"fx") == 0) thresh_array[nthresh] = FX;
    else if (strcmp(arg[1],"fy") == 0) thresh_array[nthresh] = FY;
    else if (strcmp(arg[1],"fz") == 0) thresh_array[nthresh] = FZ;
    else if (strcmp(arg[1],"q") == 0) thresh_array[nthresh] = Q;
    else if (strcmp(arg[1],"mux") == 0) thresh_array[nthresh] = MUX;
    else if (strcmp(arg[1],"muy") == 0) thresh_array[nthresh] = MUY;
    else if (strcmp(arg[1],"muz") == 0) thresh_array[nthresh] = MUZ;
    else if (strcmp(arg[1],"tqx") == 0) thresh_array[nthresh] = TQX;
    else if (strcmp(arg[1],"tqy") == 0) thresh_array[nthresh] = TQY;
    else if (strcmp(arg[1],"tqz") == 0) thresh_array[nthresh] = TQZ;
    else if (strcmp(arg[1],"centro") == 0) thresh_array[nthresh] = CENTRO;
    else if (strcmp(arg[1],"eng") == 0) thresh_array[nthresh] = ENG;
    else if (strcmp(arg[1],"sxx") == 0) thresh_array[nthresh] = SXX;
    else if (strcmp(arg[1],"syy") == 0) thresh_array[nthresh] = SYY;
    else if (strcmp(arg[1],"szz") == 0) thresh_array[nthresh] = SZZ;
    else if (strcmp(arg[1],"sxy") == 0) thresh_array[nthresh] = SXY;
    else if (strcmp(arg[1],"sxz") == 0) thresh_array[nthresh] = SXZ;
    else if (strcmp(arg[1],"syz") == 0) thresh_array[nthresh] = SYZ;
    else error->all("Invalid dump_modify threshhold operator");

    if (strcmp(arg[2],">=") == 0) thresh_op[nthresh] = LT;
    else if (strcmp(arg[2],">") == 0) thresh_op[nthresh] = LE;
    else if (strcmp(arg[2],"<=") == 0) thresh_op[nthresh] = GT;
    else if (strcmp(arg[2],"<") == 0) thresh_op[nthresh] = GE;
    else if (strcmp(arg[2],"=") == 0) thresh_op[nthresh] = NEQ;
    else if (strcmp(arg[2],"<>") == 0) thresh_op[nthresh] = EQ;
    else error->all("Invalid dump_modify threshhold operator");

    thresh_value[nthresh] = atof(arg[3]);

    // add new fix if necessary

    if (thresh_array[nthresh] == CENTRO && fixflag_centro == 0) {
      fix_create("CENTRO");
      fixflag_centro = 1;
    }
    if (thresh_array[nthresh] == ENG && fixflag_energy == 0) {
      fix_create("ENERGY");
      fixflag_energy = 1;
    }
    if ((thresh_array[nthresh] == SXX || thresh_array[nthresh] == SYY ||
	thresh_array[nthresh] == SZZ || thresh_array[nthresh] == SXY ||
	thresh_array[nthresh] == SXZ || thresh_array[nthresh] == SYZ) &&
	fixflag_stress == 0) {
      fix_create("STRESS");
      fixflag_stress = 1;
    }

    nthresh++;
    return 4;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_header(int ndump)
{
  (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

int DumpCustom::count()
{
  int i;

  // invoke any fixes that compute dump quantities

  if (fixflag) {
    if (fixflag_centro) modify->fix[ifix_centro]->dump();
    if (fixflag_energy) modify->fix[ifix_energy]->dump();
    if (fixflag_stress) modify->fix[ifix_stress]->dump();
  }

  // grow choose and dchoose arrays if needed

  int nlocal = atom->nlocal;
  if (nlocal > maxchoose) {
    delete [] choose;
    delete [] dchoose;
    maxchoose = atom->nmax;
    choose = new int[maxchoose];
    dchoose = new double[maxchoose];
  }

  // choose all local atoms for output

  for (i = 0; i < nlocal; i++) choose[i] = 1;
  nmine = nlocal;

  // un-choose if not in group

  if (igroup) {
    int *mask = atom->mask;
    for (i = 0; i < nlocal; i++)
      if (!(mask[i] & groupbit)) {
	choose[i] = 0;
	nmine--;
      }
  }

  // un-choose if not in region

  if (iregion >= 0) {
    Region *region = domain->regions[iregion];
    double **x = atom->x;
    for (i = 0; i < nlocal; i++)
      if (choose[i] && region->match(x[i][0],x[i][1],x[i][2]) == 0) {
	choose[i] = 0;
	nmine--;
      }
  }

  // un-choose if any threshhold criterion isn't met

  if (nthresh) {
    double *ptr;
    double value;
    int nstride;
    int nlocal = atom->nlocal;
    
    for (int ithresh = 0; ithresh < nthresh; ithresh++) {

      // customize by adding to if statement

      if (thresh_array[ithresh] == TAG) {
	int *tag = atom->tag;
	for (i = 0; i < nlocal; i++) dchoose[i] = tag[i];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == MOL) {
	int *molecule = atom->molecule;
	for (i = 0; i < nlocal; i++) dchoose[i] = molecule[i];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == TYPE) {
	int *type = atom->type;
	for (i = 0; i < nlocal; i++) dchoose[i] = type[i];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == X) {
	ptr = &atom->x[0][0];
	nstride = 3;
      } else if (thresh_array[ithresh] == Y) {
	ptr = &atom->x[0][1];
	nstride = 3;
      } else if (thresh_array[ithresh] == Z) {
	ptr = &atom->x[0][2];
	nstride = 3;
      } else if (thresh_array[ithresh] == XS) {
        double **x = atom->x;
	double boxxlo = domain->boxxlo;
	double invxprd = 1.0/domain->xprd;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = (x[i][0] - boxxlo) * invxprd;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == YS) {
        double **x = atom->x;
	double boxylo = domain->boxylo;
	double invyprd = 1.0/domain->yprd;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = (x[i][1] - boxylo) * invyprd;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == ZS) {
        double **x = atom->x;
	double boxzlo = domain->boxzlo;
	double invzprd = 1.0/domain->zprd;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = (x[i][2] - boxzlo) * invzprd;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == XU) {
        double **x = atom->x;
	int *image = atom->image;
	double xprd = domain->xprd;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = x[i][0] + ((image[i] & 1023) - 512) * xprd;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == YU) {
        double **x = atom->x;
	int *image = atom->image;
	double yprd = domain->yprd;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = x[i][1] + ((image[i] >> 10 & 1023) - 512) * yprd;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == ZU) {
        double **x = atom->x;
	int *image = atom->image;
	double zprd = domain->zprd;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = x[i][2] + ((image[i] >> 20) - 512) * zprd;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == IX) {
	int *image = atom->image;
	for (i = 0; i < nlocal; i++)
	  dchoose[i] = (image[i] & 1023) - 512;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == IY) {
	int *image = atom->image;
	for (i = 0; i < nlocal; i++)
	  dchoose[i] = (image[i] >> 10 & 1023) - 512;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == IX) {
	int *image = atom->image;
	for (i = 0; i < nlocal; i++)
	  dchoose[i] = (image[i] >> 20) - 512;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == VX) {
	ptr = &atom->v[0][0];
	nstride = 3;
      } else if (thresh_array[ithresh] == VY) {
	ptr = &atom->v[0][1];
	nstride = 3;
      } else if (thresh_array[ithresh] == VZ) {
	ptr = &atom->v[0][2];
	nstride = 3;
      } else if (thresh_array[ithresh] == FX) {
	ptr = &atom->f[0][0];
	nstride = 3;
      } else if (thresh_array[ithresh] == FY) {
	ptr = &atom->f[0][1];
	nstride = 3;
      } else if (thresh_array[ithresh] == FZ) {
	ptr = &atom->f[0][2];
	nstride = 3;
      } else if (thresh_array[ithresh] == Q) {
	ptr = atom->q;
	nstride = 1;
      } else if (thresh_array[ithresh] == MUX) {
	ptr = &atom->mu[0][0];
	nstride = 3;
      } else if (thresh_array[ithresh] == MUY) {
	ptr = &atom->mu[0][1];
	nstride = 3;
      } else if (thresh_array[ithresh] == MUZ) {
	ptr = &atom->mu[0][2];
	nstride = 3;
      } else if (thresh_array[ithresh] == TQX) {
	ptr = &atom->torque[0][0];
	nstride = 3;
      } else if (thresh_array[ithresh] == TQY) {
	ptr = &atom->torque[0][1];
	nstride = 3;
      } else if (thresh_array[ithresh] == TQZ) {
	ptr = &atom->torque[0][2];
	nstride = 3;
      } else if (thresh_array[ithresh] == CENTRO) {
	ptr = ((FixCentro *) modify->fix[ifix_centro])->centro;
	nstride = 3;
      } else if (thresh_array[ithresh] == ENG) {
	ptr = ((FixEnergy *) modify->fix[ifix_energy])->energy;
	nstride = 3;
      } else if (thresh_array[ithresh] == SXX) {
	ptr = &((FixStress *) modify->fix[ifix_stress])->stress[0][0];
	nstride = 6;
      } else if (thresh_array[ithresh] == SYY) {
	ptr = &((FixStress *) modify->fix[ifix_stress])->stress[0][1];
	nstride = 6;
      } else if (thresh_array[ithresh] == SZZ) {
	ptr = &((FixStress *) modify->fix[ifix_stress])->stress[0][2];
	nstride = 6;
      } else if (thresh_array[ithresh] == SXY) {
	ptr = &((FixStress *) modify->fix[ifix_stress])->stress[0][3];
	nstride = 6;
      } else if (thresh_array[ithresh] == SXZ) {
	ptr = &((FixStress *) modify->fix[ifix_stress])->stress[0][4];
	nstride = 6;
      } else if (thresh_array[ithresh] == SYZ) {
	ptr = &((FixStress *) modify->fix[ifix_stress])->stress[0][5];
	nstride = 6;
      }

      value = thresh_value[ithresh];

      if (thresh_op[ithresh] == GE) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr >= value) {
	    choose[i] = 0;
	    nmine--;
	  }
      } else if (thresh_op[ithresh] == GT) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr > value) {
	    choose[i] = 0;
	    nmine--;
	  }
      } else if (thresh_op[ithresh] == LE) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr <= value) {
	    choose[i] = 0;
	    nmine--;
	  }
      } else if (thresh_op[ithresh] == LT) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr < value) {
	    choose[i] = 0;
	    nmine--;
	  }
      } else if (thresh_op[ithresh] == EQ) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr == value) {
	    choose[i] = 0;
	    nmine--;
	  }
      } else if (thresh_op[ithresh] == NEQ) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr != value) {
	    choose[i] = 0;
	    nmine--;
	  }
      }
    }
  }

  return nmine;
}

/* ---------------------------------------------------------------------- */

int DumpCustom::pack()
{
  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
  return nmine*size_one;
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_data(int n, double *buf)
{
  (this->*write_choice)(n,buf);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_binary(int ndump)
{
  fwrite(&update->ntimestep,sizeof(int),1,fp);
  fwrite(&ndump,sizeof(int),1,fp);
  fwrite(&domain->boxxlo,sizeof(double),1,fp);
  fwrite(&domain->boxxhi,sizeof(double),1,fp);
  fwrite(&domain->boxylo,sizeof(double),1,fp);
  fwrite(&domain->boxyhi,sizeof(double),1,fp);
  fwrite(&domain->boxzlo,sizeof(double),1,fp);
  fwrite(&domain->boxzhi,sizeof(double),1,fp);
  fwrite(&size_one,sizeof(int),1,fp);
  if (multiproc) {
    int one = 1;
    fwrite(&one,sizeof(int),1,fp);
  } else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_item(int ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,"%d\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS\n");
  fprintf(fp,"%g %g\n",domain->boxxlo,domain->boxxhi);
  fprintf(fp,"%g %g\n",domain->boxylo,domain->boxyhi);
  fprintf(fp,"%g %g\n",domain->boxzlo,domain->boxzhi);
  fprintf(fp,"ITEM: ATOMS\n");
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_binary(int n, double *buf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(buf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_text(int n, double *buf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < size_one; j++) {
      if (vtype[j] == INT) fprintf(fp,vformat[j],static_cast<int> (buf[m]));
      else fprintf(fp,vformat[j],buf[m]);
      m++;
    }
    fprintf(fp,"\n");
  }
}

/* ----------------------------------------------------------------------
   create a fix of style keyword for computing a per-atom custom quantity
   fix-ID = dump-ID + keyword
------------------------------------------------------------------------- */

void DumpCustom::fix_create(char *keyword)
{
  char **fixarg = new char*[3];

  int n = strlen(id) + strlen(keyword) + 1;
  fixarg[0] = new char[n];
  strcpy(fixarg[0],id);
  strcat(fixarg[0],keyword);

  n = strlen(group->names[igroup]) + 1;
  fixarg[1] = new char[n];
  strcpy(fixarg[1],group->names[igroup]);

  n = strlen(keyword) + 1;
  fixarg[2] = new char[n];
  strcpy(fixarg[2],keyword);

  modify->add_fix(3,fixarg);
  for (int i = 0; i < 3; i++) delete [] fixarg[i];
  delete [] fixarg;
}

/* ----------------------------------------------------------------------
   delete a fix of style keyword for computing a per-atom custom quantity 
------------------------------------------------------------------------- */

void DumpCustom::fix_delete(char *keyword)
{
  int n = strlen(id) + strlen(keyword) + 1;
  char *name = new char[n];
  strcpy(name,id);
  strcat(name,keyword);
  modify->delete_fix(name);
  delete [] name;
}

/* ----------------------------------------------------------------------
   find which fix matches keyword
   fix was created by dump custom, so must exist
------------------------------------------------------------------------- */

int DumpCustom::fix_match(char *keyword)
{
  int n = strlen(id) + strlen(keyword) + 1;
  char *name = new char[n];
  strcpy(name,id);
  strcat(name,keyword);

  int i;
  for (i = 0; i < modify->nfix; i++)
    if (strcmp(name,modify->fix[i]->id) == 0) break;

  delete [] name;
  return i;
}

// ----------------------------------------------------------------------
// one routine for every kind of quantity dump custom can output
// the atom quantity is packed into buf starting at n with stride size_one
// customize by adding a method
// ----------------------------------------------------------------------

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_tag(int n)
{
  int *tag = atom->tag;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = tag[i];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_molecule(int n)
{
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = molecule[i];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_type(int n)
{
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = type[i];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_x(int n)
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = x[i][0];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_y(int n)
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = x[i][1];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_z(int n)
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = x[i][2];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_xs(int n)
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double boxxlo = domain->boxxlo;
  double invxprd = 1.0/domain->xprd;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = (x[i][0] - boxxlo) * invxprd;
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_ys(int n)
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double boxylo = domain->boxylo;
  double invyprd = 1.0/domain->yprd;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = (x[i][1] - boxylo) * invyprd;
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_zs(int n)
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double boxzlo = domain->boxzlo;
  double invzprd = 1.0/domain->zprd;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = (x[i][2] - boxzlo) * invzprd;
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_xu(int n)
{
  double **x = atom->x;
  int *image = atom->image;
  double xprd = domain->xprd;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = x[i][0] + ((image[i] & 1023) - 512) * xprd;
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_yu(int n)
{
  double **x = atom->x;
  int *image = atom->image;
  double yprd = domain->yprd;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = x[i][1] + ((image[i] >> 10 & 1023) - 512) * yprd;
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_zu(int n)
{
  double **x = atom->x;
  int *image = atom->image;
  double zprd = domain->zprd;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = x[i][2] + ((image[i] >> 20) - 512) * zprd;
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_ix(int n)
{
  int *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = (image[i] & 1023) - 512;
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_iy(int n)
{
  int *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = (image[i] >> 10 & 1023) - 512;
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_iz(int n)
{
  int *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = (image[i] >> 20) - 512;
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_vx(int n)
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = v[i][0];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_vy(int n)
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = v[i][1];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_vz(int n)
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = v[i][2];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fx(int n)
{
  double **f = atom->f;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = f[i][0];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fy(int n)
{
  double **f = atom->f;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = f[i][1];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fz(int n)
{
  double **f = atom->f;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = f[i][2];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_q(int n)
{
  double *q = atom->q;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = q[i];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_mux(int n)
{
  double **mu = atom->mu;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = mu[i][0];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_muy(int n)
{
  double **mu = atom->mu;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = mu[i][1];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_muz(int n)
{
  double **mu = atom->mu;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = mu[i][2];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_tqx(int n)
{
  double **torque = atom->torque;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = torque[i][0];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_tqy(int n)
{
  double **torque = atom->torque;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = torque[i][1];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_tqz(int n)
{
  double **torque = atom->torque;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = torque[i][2];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_centro(int n)
{
  double *centro = ((FixCentro *) modify->fix[ifix_centro])->centro;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = centro[i];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_energy(int n)
{
  double *energy = ((FixEnergy *) modify->fix[ifix_energy])->energy;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = energy[i];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_sxx(int n)
{
  double **stress = ((FixStress *) modify->fix[ifix_stress])->stress;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = stress[i][0];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_syy(int n)
{
  double **stress = ((FixStress *) modify->fix[ifix_stress])->stress;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = stress[i][1];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_szz(int n)
{
  double **stress = ((FixStress *) modify->fix[ifix_stress])->stress;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = stress[i][2];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_sxy(int n)
{
  double **stress = ((FixStress *) modify->fix[ifix_stress])->stress;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = stress[i][3];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_sxz(int n)
{
  double **stress = ((FixStress *) modify->fix[ifix_stress])->stress;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = stress[i][4];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_syz(int n)
{
  double **stress = ((FixStress *) modify->fix[ifix_stress])->stress;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = stress[i][5];
      n += size_one;
    }
}
