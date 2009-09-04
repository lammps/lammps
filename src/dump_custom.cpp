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

#include "stdlib.h"
#include "string.h"
#include "dump_custom.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

// customize by adding keyword to 1st enum

enum{ID,MOL,TYPE,MASS,
       X,Y,Z,XS,YS,ZS,XSTRI,YSTRI,ZSTRI,XU,YU,ZU,XUTRI,YUTRI,ZUTRI,IX,IY,IZ,
       VX,VY,VZ,FX,FY,FZ,
       Q,MUX,MUY,MUZ,RADIUS,OMEGAX,OMEGAY,OMEGAZ,ANGMOMX,ANGMOMY,ANGMOMZ,
       QUATW,QUATI,QUATJ,QUATK,TQX,TQY,TQZ,
       COMPUTE,FIX,VARIABLE};
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE};
enum{DUMMY0,INVOKED_SCALAR,INVOKED_VECTOR,DUMMMY3,INVOKED_PERATOM};

/* ---------------------------------------------------------------------- */

DumpCustom::DumpCustom(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg)
{
  if (narg == 5) error->all("No dump custom arguments specified");

  nevery = atoi(arg[3]);

  size_one = nfield = narg-5;
  pack_choice = new FnPtrPack[nfield];
  vtype = new int[nfield];

  iregion = -1;
  nthresh = 0;
  thresh_array = NULL;
  thresh_op = NULL;
  thresh_value = NULL;

  // computes, fixes, variables which the dump accesses

  field2index = (int *) memory->smalloc(nfield*sizeof(int),"dump:field2index");
  argindex = (int *) memory->smalloc(nfield*sizeof(int),"dump:argindex");

  ncompute = 0;
  id_compute = NULL;
  compute = NULL;

  nfix = 0;
  id_fix = NULL;
  fix = NULL;

  nvariable = 0;
  id_variable = NULL;
  variable = NULL;
  vbuf = NULL;

  // process keywords

  parse_fields(narg,arg);

  // atom selection arrays

  maxlocal = 0;
  choose = NULL;
  dchoose = NULL;

  // setup format strings

  vformat = new char*[size_one];

  format_default = new char[3*size_one+1];
  format_default[0] = '\0';

  for (int i = 0; i < size_one; i++) {
    if (vtype[i] == INT) format_default = strcat(format_default,"%d ");
    else format_default = strcat(format_default,"%g ");
    vformat[i] = NULL;
  }

  // setup column string

  int n = 0;
  for (int iarg = 5; iarg < narg; iarg++) n += strlen(arg[iarg]) + 2;
  columns = new char[n];
  columns[0] = '\0';
  for (int iarg = 5; iarg < narg; iarg++) {
    strcat(columns,arg[iarg]);
    strcat(columns," ");
  }
}

/* ---------------------------------------------------------------------- */

DumpCustom::~DumpCustom()
{
  delete [] pack_choice;
  delete [] vtype;

  memory->sfree(thresh_array);
  memory->sfree(thresh_op);
  memory->sfree(thresh_value);

  memory->sfree(field2index);
  memory->sfree(argindex);

  for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
  memory->sfree(id_compute);
  delete [] compute;

  for (int i = 0; i < nfix; i++) delete [] id_fix[i];
  memory->sfree(id_fix);
  delete [] fix;

  for (int i = 0; i < nvariable; i++) delete [] id_variable[i];
  memory->sfree(id_variable);
  delete [] variable;
  for (int i = 0; i < nvariable; i++) memory->sfree(vbuf[i]);
  delete [] vbuf;

  memory->sfree(choose);
  memory->sfree(dchoose);

  for (int i = 0; i < size_one; i++) delete [] vformat[i];
  delete [] vformat;

  delete [] columns;
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

  if (binary && domain->triclinic == 0)
    header_choice = &DumpCustom::header_binary;
  else if (binary && domain->triclinic == 1)
    header_choice = &DumpCustom::header_binary_triclinic;
  else if (!binary && domain->triclinic == 0)
    header_choice = &DumpCustom::header_item;
  else if (!binary && domain->triclinic == 1)
    header_choice = &DumpCustom::header_item_triclinic;

  if (binary) write_choice = &DumpCustom::write_binary;
  else write_choice = &DumpCustom::write_text;

  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) error->all("Could not find dump custom compute ID");
    compute[i] = modify->compute[icompute];
  }

  int ifix;
  for (int i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all("Could not find dump custom fix ID");
    fix[i] = modify->fix[ifix];
    if (nevery % modify->fix[ifix]->peratom_freq)
      error->all("Dump custom and fix not computed at compatible times");
  }

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0) error->all("Could not find dump custom variable name");
    variable[i] = ivariable;
  }

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_header(int ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_binary(int ndump)
{
  fwrite(&update->ntimestep,sizeof(int),1,fp);
  fwrite(&ndump,sizeof(int),1,fp);
  fwrite(&domain->triclinic,sizeof(int),1,fp);
  fwrite(&boxxlo,sizeof(double),1,fp);
  fwrite(&boxxhi,sizeof(double),1,fp);
  fwrite(&boxylo,sizeof(double),1,fp);
  fwrite(&boxyhi,sizeof(double),1,fp);
  fwrite(&boxzlo,sizeof(double),1,fp);
  fwrite(&boxzhi,sizeof(double),1,fp);
  fwrite(&size_one,sizeof(int),1,fp);
  if (multiproc) {
    int one = 1;
    fwrite(&one,sizeof(int),1,fp);
  } else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_binary_triclinic(int ndump)
{
  fwrite(&update->ntimestep,sizeof(int),1,fp);
  fwrite(&ndump,sizeof(int),1,fp);
  fwrite(&domain->triclinic,sizeof(int),1,fp);
  fwrite(&boxxlo,sizeof(double),1,fp);
  fwrite(&boxxhi,sizeof(double),1,fp);
  fwrite(&boxylo,sizeof(double),1,fp);
  fwrite(&boxyhi,sizeof(double),1,fp);
  fwrite(&boxzlo,sizeof(double),1,fp);
  fwrite(&boxzhi,sizeof(double),1,fp);
  fwrite(&boxxy,sizeof(double),1,fp);
  fwrite(&boxxz,sizeof(double),1,fp);
  fwrite(&boxyz,sizeof(double),1,fp);
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
  fprintf(fp,"%g %g\n",boxxlo,boxxhi);
  fprintf(fp,"%g %g\n",boxylo,boxyhi);
  fprintf(fp,"%g %g\n",boxzlo,boxzhi);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_item_triclinic(int ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%d\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,"%d\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS xy xz yz\n");
  fprintf(fp,"%g %g %g\n",boxxlo,boxxhi,boxxy);
  fprintf(fp,"%g %g %g\n",boxylo,boxyhi,boxxz);
  fprintf(fp,"%g %g %g\n",boxzlo,boxzhi,boxyz);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

int DumpCustom::count()
{
  int i;

  // grow choose and variable vbuf arrays if needed

  int nlocal = atom->nlocal;
  if (nlocal > maxlocal) {
    maxlocal = atom->nmax;

    memory->sfree(choose);
    memory->sfree(dchoose);
    choose = (int *) memory->smalloc(maxlocal*sizeof(int),"dump:choose");
    dchoose = (double *) 
      memory->smalloc(maxlocal*sizeof(double),"dump:dchoose");

    for (i = 0; i < nvariable; i++) {
      memory->sfree(vbuf[i]);
      vbuf[i] = (double *) 
	memory->smalloc(maxlocal*sizeof(double),"dump:vbuf");
    }
  }

  // invoke Computes for per-atom dump quantities

  if (ncompute) {
    int ntimestep = update->ntimestep;
    for (i = 0; i < ncompute; i++)
      if (!(compute[i]->invoked_flag & INVOKED_PERATOM)) {
	compute[i]->compute_peratom();
	compute[i]->invoked_flag |= INVOKED_PERATOM;
      }
  }

  // evaluate atom-style Variables for per-atom dump quantities

  if (nvariable)
    for (i = 0; i < nvariable; i++)
      input->variable->compute_atom(variable[i],igroup,vbuf[i],1,0);
  
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

      if (thresh_array[ithresh] == ID) {
	int *tag = atom->tag;
	for (i = 0; i < nlocal; i++) dchoose[i] = tag[i];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == MOL) {
	if (!atom->molecule_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	int *molecule = atom->molecule;
	for (i = 0; i < nlocal; i++) dchoose[i] = molecule[i];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == TYPE) {
	int *type = atom->type;
	for (i = 0; i < nlocal; i++) dchoose[i] = type[i];
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == MASS) {
	if (atom->rmass) {
	  ptr = atom->rmass;
	  nstride = 1;
	} else {
	  double *mass = atom->mass;
	  int *type = atom->type;
	  for (i = 0; i < nlocal; i++) dchoose[i] = mass[type[i]];
	  ptr = dchoose;
	  nstride = 1;
	}

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
	double boxxlo = domain->boxlo[0];
	double invxprd = 1.0/domain->xprd;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = (x[i][0] - boxxlo) * invxprd;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == YS) {
        double **x = atom->x;
	double boxylo = domain->boxlo[1];
	double invyprd = 1.0/domain->yprd;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = (x[i][1] - boxylo) * invyprd;
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == ZS) {
        double **x = atom->x;
	double boxzlo = domain->boxlo[2];
	double invzprd = 1.0/domain->zprd;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = (x[i][2] - boxzlo) * invzprd;
	ptr = dchoose;
	nstride = 1;

      } else if (thresh_array[ithresh] == XSTRI) {
        double **x = atom->x;
	double *boxlo = domain->boxlo;
	double *h_inv = domain->h_inv;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = h_inv[0]*(x[i][0]-boxlo[0]) + 
	    h_inv[5]*(x[i][1]-boxlo[1]) + h_inv[4]*(x[i][2]-boxlo[2]);
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == YSTRI) {
        double **x = atom->x;
	double *boxlo = domain->boxlo;
	double *h_inv = domain->h_inv;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = h_inv[1]*(x[i][1]-boxlo[1]) + 
	    h_inv[3]*(x[i][2]-boxlo[2]);
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == ZSTRI) {
        double **x = atom->x;
	double *boxlo = domain->boxlo;
	double *h_inv = domain->h_inv;
	for (i = 0; i < nlocal; i++) 
	  dchoose[i] = h_inv[2]*(x[i][2]-boxlo[2]);
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

      } else if (thresh_array[ithresh] == XUTRI) {
        double **x = atom->x;
	int *image = atom->image;
	double *h = domain->h;
	int xbox,ybox,zbox;
	for (i = 0; i < nlocal; i++) {
	  xbox = (image[i] & 1023) - 512;
	  ybox = (image[i] >> 10 & 1023) - 512;
	  zbox = (image[i] >> 20) - 512;
	  dchoose[i] = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
	}
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == YUTRI) {
        double **x = atom->x;
	int *image = atom->image;
	double *h = domain->h;
	int ybox,zbox;
	for (i = 0; i < nlocal; i++) {
	  ybox = (image[i] >> 10 & 1023) - 512;
	  zbox = (image[i] >> 20) - 512;
	  dchoose[i] = x[i][1] + h[1]*ybox + h[3]*zbox;
	}
	ptr = dchoose;
	nstride = 1;
      } else if (thresh_array[ithresh] == ZUTRI) {
        double **x = atom->x;
	int *image = atom->image;
	double *h = domain->h;
	int zbox;
	for (i = 0; i < nlocal; i++) {
	  zbox = (image[i] >> 20) - 512;
	  dchoose[i] = x[i][2] + h[2]*zbox;
	}
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
      } else if (thresh_array[ithresh] == IZ) {
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
	if (!atom->q_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = atom->q;
	nstride = 1;
      } else if (thresh_array[ithresh] == MUX) {
	if (!atom->mu_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->mu[0][0];
	nstride = 3;
      } else if (thresh_array[ithresh] == MUY) {
	if (!atom->mu_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->mu[0][1];
	nstride = 3;
      } else if (thresh_array[ithresh] == MUZ) {
	if (!atom->mu_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->mu[0][2];
	nstride = 3;
      } else if (thresh_array[ithresh] == RADIUS) {
	if (!atom->radius_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = atom->radius;
	nstride = 1;
      } else if (thresh_array[ithresh] == OMEGAX) {
	if (!atom->omega_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->omega[0][0];
	nstride = 3;
      } else if (thresh_array[ithresh] == OMEGAY) {
	if (!atom->omega_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->omega[0][1];
	nstride = 3;
      } else if (thresh_array[ithresh] == OMEGAZ) {
	if (!atom->omega_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->omega[0][2];
	nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMX) {
	if (!atom->angmom_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->angmom[0][0];
	nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMY) {
	if (!atom->angmom_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->angmom[0][1];
	nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMZ) {
	if (!atom->angmom_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->angmom[0][2];
	nstride = 3;
      } else if (thresh_array[ithresh] == QUATW) {
	if (!atom->quat_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->quat[0][0];
	nstride = 4;
      } else if (thresh_array[ithresh] == QUATI) {
	if (!atom->quat_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->quat[0][1];
	nstride = 4;
      } else if (thresh_array[ithresh] == QUATJ) {
	if (!atom->quat_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->quat[0][2];
	nstride = 4;
      } else if (thresh_array[ithresh] == QUATK) {
	if (!atom->quat_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->quat[0][3];
	nstride = 4;
      } else if (thresh_array[ithresh] == TQX) {
	if (!atom->torque_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->torque[0][0];
	nstride = 3;
      } else if (thresh_array[ithresh] == TQY) {
	if (!atom->torque_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->torque[0][1];
	nstride = 3;
      } else if (thresh_array[ithresh] == TQZ) {
	if (!atom->torque_flag)
	  error->all("Threshhold for an atom quantity that isn't allocated");
	ptr = &atom->torque[0][2];
	nstride = 3;

      } else if (thresh_array[ithresh] == COMPUTE) {
	i = nfield + ithresh;
	if (argindex[i] == 0) {
	  ptr = compute[field2index[i]]->scalar_atom;
	  nstride = 1;
	} else {
	  ptr = &compute[field2index[i]]->vector_atom[0][argindex[i]-1];
	  nstride = compute[field2index[i]]->size_peratom;
	}

      } else if (thresh_array[ithresh] == FIX) {
	i = nfield + ithresh;
	if (argindex[i] == 0) {
	  ptr = fix[field2index[i]]->scalar_atom;
	  nstride = 1;
	} else {
	  ptr = &fix[field2index[i]]->vector_atom[0][argindex[i]-1];
	  nstride = fix[field2index[i]]->size_peratom;
	}

      } else if (thresh_array[ithresh] == VARIABLE) {
	i = nfield + ithresh;
	ptr = vbuf[field2index[i]];
	nstride = 1;
      }

      // unselect atoms that don't meet threshhold criterion

      value = thresh_value[ithresh];

      if (thresh_op[ithresh] == LT) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr >= value) {
	    choose[i] = 0;
	    nmine--;
	  }
      } else if (thresh_op[ithresh] == LE) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr > value) {
	    choose[i] = 0;
	    nmine--;
	  }
      } else if (thresh_op[ithresh] == GT) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr <= value) {
	    choose[i] = 0;
	    nmine--;
	  }
      } else if (thresh_op[ithresh] == GE) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr < value) {
	    choose[i] = 0;
	    nmine--;
	  }
      } else if (thresh_op[ithresh] == EQ) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr != value) {
	    choose[i] = 0;
	    nmine--;
	  }
      } else if (thresh_op[ithresh] == NEQ) {
	for (i = 0; i < nlocal; i++, ptr += nstride)
	  if (choose[i] && *ptr == value) {
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

/* ---------------------------------------------------------------------- */

void DumpCustom::parse_fields(int narg, char **arg)
{
  // customize by adding to if statement

  int i;
  for (int iarg = 5; iarg < narg; iarg++) {
    i = iarg-5;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[i] = &DumpCustom::pack_id;
      vtype[i] = INT;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (!atom->molecule_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_molecule;
      vtype[i] = INT;
    } else if (strcmp(arg[iarg],"type") == 0) {
      pack_choice[i] = &DumpCustom::pack_type;
      vtype[i] = INT;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      pack_choice[i] = &DumpCustom::pack_mass;
      vtype[i] = DOUBLE;

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
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_xs_triclinic;
      else pack_choice[i] = &DumpCustom::pack_xs;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"ys") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_ys_triclinic;
      else pack_choice[i] = &DumpCustom::pack_ys;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"zs") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_zs_triclinic;
      else pack_choice[i] = &DumpCustom::pack_zs;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"xu") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_xu_triclinic;
      else pack_choice[i] = &DumpCustom::pack_xu;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"yu") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_yu_triclinic;
      else pack_choice[i] = &DumpCustom::pack_yu;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"zu") == 0) {
      if (domain->triclinic) pack_choice[i] = &DumpCustom::pack_zu_triclinic;
      else pack_choice[i] = &DumpCustom::pack_zu;
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
      if (!atom->q_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_q;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"mux") == 0) {
      if (!atom->mu_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_mux;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"muy") == 0) {
      if (!atom->mu_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_muy;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"muz") == 0) {
      if (!atom->mu_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_muz;
      vtype[i] = DOUBLE;

    } else if (strcmp(arg[iarg],"radius") == 0) {
      if (!atom->radius_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_radius;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"omegax") == 0) {
      if (!atom->omega_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_omegax;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"omegay") == 0) {
      if (!atom->omega_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_omegay;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"omegaz") == 0) {
      if (!atom->omega_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_omegaz;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"angmomx") == 0) {
      if (!atom->angmom_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_angmomx;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"angmomy") == 0) {
      if (!atom->angmom_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_angmomy;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"angmomz") == 0) {
      if (!atom->angmom_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_angmomz;
      vtype[i] = DOUBLE;

    } else if (strcmp(arg[iarg],"quatw") == 0) {
      if (!atom->quat_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_quatw;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"quati") == 0) {
      if (!atom->quat_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_quati;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"quatj") == 0) {
      if (!atom->quat_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_quatj;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"quatk") == 0) {
      if (!atom->quat_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_quatk;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"tqx") == 0) {
      if (!atom->torque_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_tqx;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"tqy") == 0) {
      if (!atom->torque_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_tqy;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"tqz") == 0) {
      if (!atom->torque_flag)
	error->all("Dumping an atom quantity that isn't allocated");
      pack_choice[i] = &DumpCustom::pack_tqz;
      vtype[i] = DOUBLE;

    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is between []

    } else if (strncmp(arg[iarg],"c_",2) == 0) {
      pack_choice[i] = &DumpCustom::pack_compute;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all("Invalid keyword in dump custom command");
	argindex[i] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[i] = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all("Could not find dump custom compute ID");
      if (modify->compute[n]->peratom_flag == 0)
	error->all("Dump custom compute ID does not compute peratom info");
      if (argindex[i] == 0 && modify->compute[n]->size_peratom > 0)
	error->all("Dump custom compute ID does not compute scalar per atom");
      if (argindex[i] > 0 && modify->compute[n]->size_peratom == 0)
	error->all("Dump custom compute ID does not compute vector per atom");
      if (argindex[i] > 0 && 
	  argindex[i] > modify->compute[n]->size_peratom)
	error->all("Dump custom compute ID vector is not large enough");

      field2index[i] = add_compute(suffix);
      delete [] suffix;
      
    // fix value = f_ID
    // if no trailing [], then arg is set to 0, else arg is between []

    } else if (strncmp(arg[iarg],"f_",2) == 0) {
      pack_choice[i] = &DumpCustom::pack_fix;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all("Invalid keyword in dump custom command");
	argindex[i] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[i] = 0;

      n = modify->find_fix(suffix);
      if (n < 0) error->all("Could not find dump custom fix ID");
      if (modify->fix[n]->peratom_flag == 0)
	error->all("Dump custom fix ID does not compute peratom info");
      if (argindex[i] == 0 && modify->fix[n]->size_peratom > 0)
	error->all("Dump custom fix ID does not compute scalar per atom");
      if (argindex[i] > 0 && modify->fix[n]->size_peratom == 0)
	error->all("Dump custom fix ID does not compute vector per atom");
      if (argindex[i] > 0 && 
	  argindex[i] > modify->fix[n]->size_peratom)
	error->all("Dump custom fix ID vector is not large enough");

      field2index[i] = add_fix(suffix);
      delete [] suffix;

    // variable value = v_name

    } else if (strncmp(arg[iarg],"v_",2) == 0) {
      pack_choice[i] = &DumpCustom::pack_variable;
      vtype[i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      argindex[i] = 0;

      n = input->variable->find(suffix);
      if (n < 0) error->all("Could not find dump custom variable name");
      if (input->variable->atomstyle(n) == 0)
	error->all("Dump custom variable is not atom-style variable");

      field2index[i] = add_variable(suffix);
      delete [] suffix;

    } else error->all("Invalid keyword in dump custom command");
  }
}

/* ----------------------------------------------------------------------
   add Compute to list of Compute objects used by dump
   return index of where this Compute is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustom::add_compute(char *id)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,id_compute[icompute]) == 0) break;
  if (icompute < ncompute) return icompute;
  
  id_compute = (char **)
    memory->srealloc(id_compute,(ncompute+1)*sizeof(char *),"dump:id_compute");
  delete [] compute;
  compute = new Compute*[ncompute+1];

  int n = strlen(id) + 1;
  id_compute[ncompute] = new char[n];
  strcpy(id_compute[ncompute],id);
  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   add Fix to list of Fix objects used by dump
   return index of where this Fix is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustom::add_fix(char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,id_fix[ifix]) == 0) break;
  if (ifix < nfix) return ifix;
  
  id_fix = (char **)
    memory->srealloc(id_fix,(nfix+1)*sizeof(char *),"dump:id_fix");
  delete [] fix;
  fix = new Fix*[nfix+1];

  int n = strlen(id) + 1;
  id_fix[nfix] = new char[n];
  strcpy(id_fix[nfix],id);
  nfix++;
  return nfix-1;
}

/* ----------------------------------------------------------------------
   add Variable to list of Variables used by dump
   return index of where this Variable is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustom::add_variable(char *id)
{
  int ivariable;
  for (ivariable = 0; ivariable < nvariable; ivariable++)
    if (strcmp(id,id_variable[ivariable]) == 0) break;
  if (ivariable < nvariable) return ivariable;
  
  id_variable = (char **)
    memory->srealloc(id_variable,(nvariable+1)*sizeof(char *),
		     "dump:id_variable");
  delete [] variable;
  variable = new int[nvariable+1];
  delete [] vbuf;
  vbuf = new double*[nvariable+1];
  for (int i = 0; i <= nvariable; i++) vbuf[i] = NULL;

  int n = strlen(id) + 1;
  id_variable[nvariable] = new char[n];
  strcpy(id_variable[nvariable],id);
  nvariable++;
  return nvariable-1;
}

/* ---------------------------------------------------------------------- */

int DumpCustom::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"region") == 0) {
    if (narg < 2) error->all("Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) iregion = -1;
    else {
      iregion = domain->find_region(arg[1]);
      if (iregion == -1) error->all("Dump_modify region ID does not exist");
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
    
    // grow threshhold arrays
    
    thresh_array = (int *)
      memory->srealloc(thresh_array,(nthresh+1)*sizeof(int),
		       "dump:thresh_array");
    thresh_op = (int *)
      memory->srealloc(thresh_op,(nthresh+1)*sizeof(int),
		       "dump:thresh_op");
    thresh_value = (double *)
      memory->srealloc(thresh_value,(nthresh+1)*sizeof(double),
		       "dump:thresh_value");

    // set keyword type of threshhold
    // customize by adding to if statement
    
    if (strcmp(arg[1],"id") == 0) thresh_array[nthresh] = ID;
    else if (strcmp(arg[1],"mol") == 0) thresh_array[nthresh] = MOL;
    else if (strcmp(arg[1],"type") == 0) thresh_array[nthresh] = TYPE;
    else if (strcmp(arg[1],"mass") == 0) thresh_array[nthresh] = MASS;

    else if (strcmp(arg[1],"x") == 0) thresh_array[nthresh] = X;
    else if (strcmp(arg[1],"y") == 0) thresh_array[nthresh] = Y;
    else if (strcmp(arg[1],"z") == 0) thresh_array[nthresh] = Z;

    else if (strcmp(arg[1],"xs") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = XS;
    else if (strcmp(arg[1],"xs") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = XSTRI;
    else if (strcmp(arg[1],"ys") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = YS;
    else if (strcmp(arg[1],"ys") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = YSTRI;
    else if (strcmp(arg[1],"zs") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = ZS;
    else if (strcmp(arg[1],"zs") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = ZSTRI;

    else if (strcmp(arg[1],"xu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = XU;
    else if (strcmp(arg[1],"xu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = XUTRI;
    else if (strcmp(arg[1],"yu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = YU;
    else if (strcmp(arg[1],"yu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = YUTRI;
    else if (strcmp(arg[1],"zu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = ZU;
    else if (strcmp(arg[1],"zu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = ZUTRI;

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
    else if (strcmp(arg[1],"radius") == 0) thresh_array[nthresh] = RADIUS;
    else if (strcmp(arg[1],"omegax") == 0) thresh_array[nthresh] = OMEGAX;
    else if (strcmp(arg[1],"omegay") == 0) thresh_array[nthresh] = OMEGAY;
    else if (strcmp(arg[1],"omegaz") == 0) thresh_array[nthresh] = OMEGAZ;
    else if (strcmp(arg[1],"angmomx") == 0) thresh_array[nthresh] = ANGMOMX;
    else if (strcmp(arg[1],"angmomy") == 0) thresh_array[nthresh] = ANGMOMY;
    else if (strcmp(arg[1],"angmomz") == 0) thresh_array[nthresh] = ANGMOMZ;
    else if (strcmp(arg[1],"quatw") == 0) thresh_array[nthresh] = QUATW;
    else if (strcmp(arg[1],"quati") == 0) thresh_array[nthresh] = QUATI;
    else if (strcmp(arg[1],"quatj") == 0) thresh_array[nthresh] = QUATJ;
    else if (strcmp(arg[1],"quatk") == 0) thresh_array[nthresh] = QUATK;
    else if (strcmp(arg[1],"tqx") == 0) thresh_array[nthresh] = TQX;
    else if (strcmp(arg[1],"tqy") == 0) thresh_array[nthresh] = TQY;
    else if (strcmp(arg[1],"tqz") == 0) thresh_array[nthresh] = TQZ;
    
    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is between []
    // must grow field2index and argindex arrays, since access is beyond nfield

    else if (strncmp(arg[1],"c_",2) == 0) {
      thresh_array[nthresh] = COMPUTE;
      field2index = (int *) memory->srealloc(field2index,
					     (nfield+nthresh+1)*sizeof(int),
					     "dump:field2index");
      argindex = (int *) memory->srealloc(argindex,
					  (nfield+nthresh+1)*sizeof(int),
					  "dump:argindex");
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);
    
      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all("Invalid keyword in dump custom command");
	argindex[nfield+nthresh] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[nfield+nthresh] = 0;
      
      n = modify->find_compute(suffix);
      if (n < 0) error->all("Could not find dump custom compute ID");

      if (modify->compute[n]->peratom_flag == 0)
	error->all("Dump custom compute ID does not compute peratom info");
      if (argindex[nfield+nthresh] == 0 && 
	  modify->compute[n]->size_peratom > 0)
	error->all("Dump custom compute ID does not compute scalar per atom");
      if (argindex[nfield+nthresh] > 0 && 
	  modify->compute[n]->size_peratom == 0)
	error->all("Dump custom compute ID does not compute vector per atom");
      if (argindex[nfield+nthresh] > 0 && 
	  argindex[nfield+nthresh] > modify->compute[n]->size_peratom)
	error->all("Dump custom compute ID vector is not large enough");

      field2index[nfield+nthresh] = add_compute(suffix);
      delete [] suffix;

    // fix value = f_ID
    // if no trailing [], then arg is set to 0, else arg is between []
    // must grow field2index and argindex arrays, since access is beyond nfield

    } else if (strncmp(arg[1],"f_",2) == 0) {
      thresh_array[nthresh] = FIX;
      field2index = (int *) memory->srealloc(field2index,
					     (nfield+nthresh+1)*sizeof(int),
					     "dump:field2index");
      argindex = (int *) memory->srealloc(argindex,
					  (nfield+nthresh+1)*sizeof(int),
					  "dump:argindex");
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);
    
      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all("Invalid keyword in dump custom command");
	argindex[nfield+nthresh] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[nfield+nthresh] = 0;
      
      n = modify->find_fix(suffix);
      if (n < 0) error->all("Could not find dump custom fix ID");

      if (modify->fix[n]->peratom_flag == 0)
	error->all("Dump custom fix ID does not compute peratom info");
      if (argindex[nfield+nthresh] == 0 && 
	  modify->fix[n]->size_peratom > 0)
	error->all("Dump custom fix ID does not compute scalar per atom");
      if (argindex[nfield+nthresh] > 0 && 
	  modify->fix[n]->size_peratom == 0)
	error->all("Dump custom fix ID does not compute vector per atom");
      if (argindex[nfield+nthresh] > 0 && 
	  argindex[nfield+nthresh] > modify->fix[n]->size_peratom)
	error->all("Dump custom fix ID vector is not large enough");

      field2index[nfield+nthresh] = add_fix(suffix);
      delete [] suffix;

    // variable value = v_ID
    // must grow field2index and argindex arrays, since access is beyond nfield

    } else if (strncmp(arg[1],"v_",2) == 0) {
      thresh_array[nthresh] = VARIABLE;
      field2index = (int *) memory->srealloc(field2index,
					     (nfield+nthresh+1)*sizeof(int),
					     "dump:field2index");
      argindex = (int *) memory->srealloc(argindex,
					  (nfield+nthresh+1)*sizeof(int),
					  "dump:argindex");
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);
    
      argindex[nfield+nthresh] = 0;
      
      n = input->variable->find(suffix);
      if (n < 0) error->all("Could not find dump custom variable name");
      if (input->variable->atomstyle(n) == 0)
	error->all("Dump custom variable is not atom-style variable");

      field2index[nfield+nthresh] = add_variable(suffix);
      delete [] suffix;

    } else error->all("Invalid dump_modify threshhold operator");

    // set operation type of threshhold

    if (strcmp(arg[2],"<") == 0) thresh_op[nthresh] = LT;
    else if (strcmp(arg[2],"<=") == 0) thresh_op[nthresh] = LE;
    else if (strcmp(arg[2],">") == 0) thresh_op[nthresh] = GT;
    else if (strcmp(arg[2],">=") == 0) thresh_op[nthresh] = GE;
    else if (strcmp(arg[2],"==") == 0) thresh_op[nthresh] = EQ;
    else if (strcmp(arg[2],"!=") == 0) thresh_op[nthresh] = NEQ;
    else error->all("Invalid dump_modify threshhold operator");

    // set threshhold value

    thresh_value[nthresh] = atof(arg[3]);

    nthresh++;
    return 4;

  // pass along params to child class

  } else {
    int n = modify_param2(narg,arg);
    return n;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf, choose, variable arrays
------------------------------------------------------------------------- */

double DumpCustom::memory_usage()
{
  double bytes = maxbuf * sizeof(double);
  bytes += maxlocal * sizeof(int);
  bytes += maxlocal * sizeof(double);
  bytes += maxlocal * nvariable * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
------------------------------------------------------------------------- */

void DumpCustom::pack_compute(int n)
{
  double *vector = compute[field2index[n]]->scalar_atom;
  double **array = compute[field2index[n]]->vector_atom;
  int index = argindex[n];

  int nlocal = atom->nlocal;

  if (index == 0) {
    for (int i = 0; i < nlocal; i++)
      if (choose[i]) {
	buf[n] = vector[i];
	n += size_one;
      }
  } else {
    index--;
    for (int i = 0; i < nlocal; i++)
      if (choose[i]) {
	buf[n] = array[i][index];
	n += size_one;
      }
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fix(int n)
{
  double *vector = fix[field2index[n]]->scalar_atom;
  double **array = fix[field2index[n]]->vector_atom;
  int index = argindex[n];

  int nlocal = atom->nlocal;

  if (index == 0) {
    for (int i = 0; i < nlocal; i++)
      if (choose[i]) {
	buf[n] = vector[i];
	n += size_one;
      }
  } else {
    index--;
    for (int i = 0; i < nlocal; i++)
      if (choose[i]) {
	buf[n] = array[i][index];
	n += size_one;
      }
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_variable(int n)
{
  double *vector = vbuf[field2index[n]];

  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = vector[i];
      n += size_one;
    }
}

/* ----------------------------------------------------------------------
   one method for every keyword dump custom can output
   the atom quantity is packed into buf starting at n with stride size_one
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_id(int n)
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

void DumpCustom::pack_mass(int n)
{
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (choose[i]) {
	buf[n] = rmass[i];
	n += size_one;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (choose[i]) {
	buf[n] = mass[type[i]];
	n += size_one;
      }
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

  double boxxlo = domain->boxlo[0];
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

  double boxylo = domain->boxlo[1];
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

  double boxzlo = domain->boxlo[2];
  double invzprd = 1.0/domain->zprd;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = (x[i][2] - boxzlo) * invzprd;
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_xs_triclinic(int n)
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = h_inv[0]*(x[i][0]-boxlo[0]) + 
	h_inv[5]*(x[i][1]-boxlo[1]) + h_inv[4]*(x[i][2]-boxlo[2]);
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_ys_triclinic(int n)
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = h_inv[1]*(x[i][1]-boxlo[1]) + h_inv[3]*(x[i][2]-boxlo[2]);
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_zs_triclinic(int n)
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = h_inv[2]*(x[i][2]-boxlo[2]);
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_xu(int n)
{
  double **x = atom->x;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;

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
  int nlocal = atom->nlocal;

  double yprd = domain->yprd;

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
  int nlocal = atom->nlocal;

  double zprd = domain->zprd;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = x[i][2] + ((image[i] >> 20) - 512) * zprd;
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_xu_triclinic(int n)
{
  double **x = atom->x;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  int xbox,ybox,zbox;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      buf[n] = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_yu_triclinic(int n)
{
  double **x = atom->x;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  int ybox,zbox;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      buf[n] = x[i][1] + h[1]*ybox + h[3]*zbox;
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_zu_triclinic(int n)
{
  double **x = atom->x;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  int zbox;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      zbox = (image[i] >> 20) - 512;
      buf[n] = x[i][2] + h[2]*zbox;
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

void DumpCustom::pack_radius(int n)
{
  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = radius[i];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_omegax(int n)
{
  double **omega = atom->omega;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = omega[i][0];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_omegay(int n)
{
  double **omega = atom->omega;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = omega[i][1];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_omegaz(int n)
{
  double **omega = atom->omega;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = omega[i][2];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_angmomx(int n)
{
  double **angmom = atom->angmom;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = angmom[i][0];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_angmomy(int n)
{
  double **angmom = atom->angmom;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = angmom[i][1];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_angmomz(int n)
{
  double **angmom = atom->angmom;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = angmom[i][2];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_quatw(int n)
{
  double **quat = atom->quat;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = quat[i][0];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_quati(int n)
{
  double **quat = atom->quat;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = quat[i][1];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_quatj(int n)
{
  double **quat = atom->quat;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = quat[i][2];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_quatk(int n)
{
  double **quat = atom->quat;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = quat[i][3];
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
