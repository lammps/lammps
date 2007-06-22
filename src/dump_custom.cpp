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
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

// customize by adding keyword to 1st enum

enum{TAG,MOL,TYPE,X,Y,Z,XS,YS,ZS,XU,YU,ZU,IX,IY,IZ,
     VX,VY,VZ,FX,FY,FZ,
     Q,MUX,MUY,MUZ,
     QUATW,QUATI,QUATJ,QUATK,TQX,TQY,TQZ,
     EPAIR,KE,ETOTAL,CENTRO,SXX,SYY,SZZ,SXY,SXZ,SYZ,
     COMPUTE};
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE};

/* ---------------------------------------------------------------------- */

DumpCustom::DumpCustom(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg)
{
  if (narg == 5) error->all("No dump custom arguments specified");

  size_one = nfield = narg-5;
  pack_choice = new FnPtrPack[nfield];
  vtype = new int[nfield];

  iregion = -1;
  nthresh = 0;
  thresh_array = NULL;
  thresh_op = NULL;
  thresh_value = NULL;

  // flags, IDs, and memory for compute objects dump may create

  index_epair = index_ke = index_etotal = index_centro = index_stress = -1;

  style_epair = "epair/atom";
  style_ke = "ke/atom";
  style_etotal = "etotal/atom";
  style_centro = "centro/atom";
  style_stress = "stress/atom";

  ncompute = 0;
  id_compute = NULL;
  compute = NULL;
  field2compute = (int *) memory->smalloc(nfield*sizeof(int),
					  "dump:field2compute");
  arg_compute = (int *) memory->smalloc(nfield*sizeof(int),"dump:arg_compute");

  // process keywords

  parse_fields(narg,arg);

  // create the requested Computes

  if (index_epair >= 0) create_compute(style_epair,NULL);
  if (index_ke >= 0) create_compute(style_ke,NULL);
  if (index_etotal >= 0) create_compute(style_etotal,style_epair);
  if (index_centro >= 0) create_compute(style_centro,NULL);
  if (index_stress >= 0) create_compute(style_stress,NULL);

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

  // one-time file open

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

DumpCustom::~DumpCustom()
{
  delete [] pack_choice;
  delete [] vtype;

  memory->sfree(thresh_array);
  memory->sfree(thresh_op);
  memory->sfree(thresh_value);

  // delete Compute classes if dump custom created them

  if (index_epair >= 0)  modify->delete_compute(id_compute[index_epair]);
  if (index_ke >= 0)  modify->delete_compute(id_compute[index_ke]);
  if (index_etotal >= 0)  modify->delete_compute(id_compute[index_etotal]);
  if (index_centro >= 0)  modify->delete_compute(id_compute[index_centro]);
  if (index_stress >= 0)  modify->delete_compute(id_compute[index_stress]);

  for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
  memory->sfree(id_compute);
  memory->sfree(compute);
  memory->sfree(field2compute);
  memory->sfree(arg_compute);

  delete [] choose;
  delete [] dchoose;

  for (int i = 0; i < size_one; i++) delete [] vformat[i];
  delete [] vformat;
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

  // find current ptr for each Compute ID

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) error->all("Could not find dump custom compute ID");
    compute[i] = modify->compute[icompute];
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_header(int ndump)
{
  (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_binary(int ndump)
{
  fwrite(&update->ntimestep,sizeof(int),1,fp);
  fwrite(&ndump,sizeof(int),1,fp);
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
  fprintf(fp,"ITEM: ATOMS\n");
}

/* ---------------------------------------------------------------------- */

int DumpCustom::count()
{
  int i;

  // grow choose arrays if needed

  int nlocal = atom->nlocal;
  if (nlocal > maxlocal) {
    delete [] choose;
    delete [] dchoose;
    maxlocal = atom->nmax;
    choose = new int[maxlocal];
    dchoose = new double[maxlocal];
  }

  // invoke Computes for per-atom dump quantities

  if (ncompute)
    for (i = 0; i < ncompute; i++) compute[i]->compute_peratom();

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
      } else if (thresh_array[ithresh] == QUATW) {
	ptr = &atom->quat[0][0];
	nstride = 4;
      } else if (thresh_array[ithresh] == QUATI) {
	ptr = &atom->quat[0][1];
	nstride = 4;
      } else if (thresh_array[ithresh] == QUATJ) {
	ptr = &atom->quat[0][2];
	nstride = 4;
      } else if (thresh_array[ithresh] == QUATK) {
	ptr = &atom->quat[0][3];
	nstride = 4;
      } else if (thresh_array[ithresh] == TQX) {
	ptr = &atom->torque[0][0];
	nstride = 3;
      } else if (thresh_array[ithresh] == TQY) {
	ptr = &atom->torque[0][1];
	nstride = 3;
      } else if (thresh_array[ithresh] == TQZ) {
	ptr = &atom->torque[0][2];
	nstride = 3;
      } else if (thresh_array[ithresh] == EPAIR) {
	ptr = compute[index_epair]->scalar_atom;
	nstride = 1;
      } else if (thresh_array[ithresh] == KE) {
	ptr = compute[index_ke]->scalar_atom;
	nstride = 1;
      } else if (thresh_array[ithresh] == ETOTAL) {
	ptr = compute[index_etotal]->scalar_atom;
	nstride = 1;
      } else if (thresh_array[ithresh] == CENTRO) {
	ptr = compute[index_centro]->scalar_atom;
	nstride = 1;
      } else if (thresh_array[ithresh] == SXX) {
	ptr = &compute[index_stress]->vector_atom[0][0];
	nstride = 6;
      } else if (thresh_array[ithresh] == SYY) {
	ptr = &compute[index_stress]->vector_atom[0][1];
	nstride = 6;
      } else if (thresh_array[ithresh] == SZZ) {
	ptr = &compute[index_stress]->vector_atom[0][2];
	nstride = 6;
      } else if (thresh_array[ithresh] == SXY) {
	ptr = &compute[index_stress]->vector_atom[0][3];
	nstride = 6;
      } else if (thresh_array[ithresh] == SXZ) {
	ptr = &compute[index_stress]->vector_atom[0][4];
	nstride = 6;
      } else if (thresh_array[ithresh] == SYZ) {
	ptr = &compute[index_stress]->vector_atom[0][5];
	nstride = 6;
      } else if (thresh_array[ithresh] == COMPUTE) {
	i = nfield + ithresh;
	if (arg_compute[i] == 0) {
	  ptr = compute[field2compute[i]]->scalar_atom;
	  nstride = 1;
	} else {
	  ptr = &compute[field2compute[i]]->vector_atom[0][arg_compute[i]-1];
	  nstride = compute[field2compute[i]]->size_peratom;
	}
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

    if (strcmp(arg[iarg],"tag") == 0) {
      pack_choice[i] = &DumpCustom::pack_tag;
      vtype[i] = INT;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (!atom->molecule_flag)
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
      if (domain->triclinic)
	error->all("Cannot dump scaled coords with triclinic box");
      pack_choice[i] = &DumpCustom::pack_xs;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"ys") == 0) {
      if (domain->triclinic)
	error->all("Cannot dump scaled coords with triclinic box");
      pack_choice[i] = &DumpCustom::pack_ys;
      vtype[i] = DOUBLE;
    } else if (strcmp(arg[iarg],"zs") == 0) {
      if (domain->triclinic)
	error->all("Cannot dump scaled coords with triclinic box");
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

    } else if (strcmp(arg[iarg],"epair") == 0) {
      pack_choice[i] = &DumpCustom::pack_epair;
      vtype[i] = DOUBLE;
      index_epair = add_compute(style_epair,1);
    } else if (strcmp(arg[iarg],"ke") == 0) {
      pack_choice[i] = &DumpCustom::pack_ke;
      vtype[i] = DOUBLE;
      index_ke = add_compute(style_ke,1);
    } else if (strcmp(arg[iarg],"etotal") == 0) {
      pack_choice[i] = &DumpCustom::pack_etotal;
      vtype[i] = DOUBLE;
      index_epair = add_compute(style_epair,1);
      index_etotal = add_compute(style_etotal,1);
    } else if (strcmp(arg[iarg],"centro") == 0) {
      pack_choice[i] = &DumpCustom::pack_centro;
      vtype[i] = DOUBLE;
      index_centro = add_compute(style_centro,1);

    } else if (strcmp(arg[iarg],"sxx") == 0) {
      pack_choice[i] = &DumpCustom::pack_sxx;
      vtype[i] = DOUBLE;
      index_stress = add_compute(style_stress,1);
    } else if (strcmp(arg[iarg],"syy") == 0) {
      pack_choice[i] = &DumpCustom::pack_syy;
      vtype[i] = DOUBLE;
      index_stress = add_compute(style_stress,1);
    } else if (strcmp(arg[iarg],"szz") == 0) {
      pack_choice[i] = &DumpCustom::pack_szz;
      vtype[i] = DOUBLE;
      index_stress = add_compute(style_stress,1);
    } else if (strcmp(arg[iarg],"sxy") == 0) {
      pack_choice[i] = &DumpCustom::pack_sxy;
      vtype[i] = DOUBLE;
      index_stress = add_compute(style_stress,1);
    } else if (strcmp(arg[iarg],"sxz") == 0) {
      pack_choice[i] = &DumpCustom::pack_sxz;
      vtype[i] = DOUBLE;
      index_stress = add_compute(style_stress,1);
    } else if (strcmp(arg[iarg],"syz") == 0) {
      pack_choice[i] = &DumpCustom::pack_syz;
      vtype[i] = DOUBLE;
      index_stress = add_compute(style_stress,1);

    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is between []
    // if Compute has pre-compute, first add it to list

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
	arg_compute[i] = atoi(ptr+1);
	*ptr = '\0';
      } else arg_compute[i] = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all("Could not find dump custom compute ID");
      if (modify->compute[n]->peratom_flag == 0)
	error->all("Dump custom compute ID does not compute peratom info");
      if (arg_compute[i] == 0 && modify->compute[n]->size_peratom > 0)
	error->all("Dump custom compute ID does not compute scalar per atom");
      if (arg_compute[i] > 0 && modify->compute[n]->size_peratom == 0)
	error->all("Dump custom compute ID does not compute vector per atom");
      if (arg_compute[i] > 0 && 
	  arg_compute[i] > modify->compute[n]->size_peratom)
	error->all("Dump custom compute ID vector is not large enough");
      if (modify->compute[n]->id_pre)
	int tmp = add_compute(modify->compute[n]->id_pre,0);
      field2compute[i] = add_compute(suffix,0);
      delete [] suffix;

    } else error->all("Invalid keyword in dump custom command");
  }
}

/* ----------------------------------------------------------------------
   add Compute to list of Compute objects to call
   return index of where this Compute is in call list
   compute ID = dump-ID + "_" + keyword if appendflag is set, else just keyword
   if already in call list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustom::add_compute(char *keyword, int appendflag)
{
  int n = strlen(id) + strlen(keyword) + 2;
  char *name = new char[n];
  if (appendflag) {
    strcpy(name,id);
    strcat(name,"_");
    strcat(name,keyword);
  } else strcpy(name,keyword);

  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(name,id_compute[icompute]) == 0) break;
  if (icompute < ncompute) {
    delete [] name;
    return icompute;
  }
  
  id_compute = (char **)
    memory->srealloc(id_compute,(ncompute+1)*sizeof(char *),"dump:id_compute");
  compute = (Compute **) 
    memory->srealloc(compute,(ncompute+1)*sizeof(Compute *),"dump:compute");

  n = strlen(name) + 1;
  id_compute[ncompute] = new char[n];
  strcpy(id_compute[ncompute],name);
  delete [] name;
  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   create a compute
   compute ID = dump-ID + "_" + keyword, compute style = keyword
   pass additional extra arg to Modify::add_compute() if defined
------------------------------------------------------------------------- */

void DumpCustom::create_compute(char *keyword, char *extra)
{
  int n = strlen(id) + strlen(keyword) + 2;
  char *name = new char[n];
  strcpy(name,id);
  strcat(name,"_");
  strcat(name,keyword);

  char **newarg = new char*[4];
  newarg[0] = name;
  newarg[1] = group->names[igroup];
  newarg[2] = keyword;

  if (extra) {
    n = strlen(id) + strlen(extra) + 2;
    newarg[3] = new char[n];
    strcpy(newarg[3],id);
    strcat(newarg[3],"_");
    strcat(newarg[3],extra);
  } else newarg[3] = NULL;

  if (extra) modify->add_compute(4,newarg);
  else modify->add_compute(3,newarg);

  delete [] name;
  delete [] newarg[3];
  delete [] newarg;
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

    // error check for triclinic box

    if (domain->triclinic &&
	(strcmp(arg[1],"xs") == 0 || strcmp(arg[1],"ys") == 0 || 
	 strcmp(arg[1],"zs") == 0))
	error->all("Cannot dump scaled coords with triclinic box");

    // set keyword type of threshhold
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
    else if (strcmp(arg[1],"quatw") == 0) thresh_array[nthresh] = QUATW;
    else if (strcmp(arg[1],"quati") == 0) thresh_array[nthresh] = QUATI;
    else if (strcmp(arg[1],"quatj") == 0) thresh_array[nthresh] = QUATJ;
    else if (strcmp(arg[1],"quatk") == 0) thresh_array[nthresh] = QUATK;
    else if (strcmp(arg[1],"tqx") == 0) thresh_array[nthresh] = TQX;
    else if (strcmp(arg[1],"tqy") == 0) thresh_array[nthresh] = TQY;
    else if (strcmp(arg[1],"tqz") == 0) thresh_array[nthresh] = TQZ;
    else if (strcmp(arg[1],"epair") == 0) {
      thresh_array[nthresh] = EPAIR;
      if (index_epair < 0) {
	index_epair = add_compute(style_epair,1);
	create_compute(style_epair,NULL);
      }
    } else if (strcmp(arg[1],"ke") == 0) {
      thresh_array[nthresh] = KE;
      if (index_ke < 0) {
	index_ke = add_compute(style_ke,1);
	create_compute(style_ke,NULL);
      }
    } else if (strcmp(arg[1],"etotal") == 0) {
      thresh_array[nthresh] = ETOTAL;
      if (index_etotal < 0) {
	if (index_epair < 0) {
	  index_epair = add_compute(style_epair,1);
	  create_compute(style_epair,NULL);
	}
	index_etotal = add_compute(style_etotal,1);
	create_compute(style_etotal,style_epair);
      }
    } else if (strcmp(arg[1],"centro") == 0) {
      thresh_array[nthresh] = CENTRO;
      if (index_centro < 0) {
	index_centro = add_compute(style_centro,1);
	create_compute(style_centro,NULL);
      }
    } else if (strcmp(arg[1],"sxx") == 0) {
      thresh_array[nthresh] = SXX;
      if (index_stress < 0) {
	index_stress = add_compute(style_stress,1);
	create_compute(style_stress,NULL);
      }
    } else if (strcmp(arg[1],"syy") == 0) {
      thresh_array[nthresh] = SYY;
      if (index_stress < 0) {
	index_stress = add_compute(style_stress,1);
	create_compute(style_stress,NULL);
      }
    } else if (strcmp(arg[1],"szz") == 0) {
      thresh_array[nthresh] = SZZ;
      if (index_stress < 0) {
	index_stress = add_compute(style_stress,1);
	create_compute(style_stress,NULL);
      }
    } else if (strcmp(arg[1],"sxy") == 0) {
      thresh_array[nthresh] = SXY;
      if (index_stress < 0) {
	index_stress = add_compute(style_stress,1);
	create_compute(style_stress,NULL);
      }
    } else if (strcmp(arg[1],"sxz") == 0) {
      thresh_array[nthresh] = SXZ;
      if (index_stress < 0) {
	index_stress = add_compute(style_stress,1);
	create_compute(style_stress,NULL);
      }
    } else if (strcmp(arg[1],"syz") == 0) {
      thresh_array[nthresh] = SYZ;
      if (index_stress < 0) {
	index_stress = add_compute(style_stress,1);
	create_compute(style_stress,NULL);
      }
    
    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is between []
    // must grow field2compute and arg_compute arrays,
    //   since access is beyond nfield
    // if Compute has pre-compute, first add it to list

    } else if (strncmp(arg[1],"c_",2) == 0) {
      thresh_array[nthresh] = COMPUTE;
      field2compute = (int *) memory->srealloc(field2compute,
					       (nfield+nthresh+1)*sizeof(int),
					       "dump:field2compute");
      arg_compute = (int *) memory->srealloc(arg_compute,
					     (nfield+nthresh+1)*sizeof(int),
					     "dump:arg_compute");
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);
    
      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all("Invalid keyword in dump custom command");
	arg_compute[nfield+nthresh] = atoi(ptr+1);
	*ptr = '\0';
      } else arg_compute[nfield+nthresh] = 0;
      
      n = modify->find_compute(suffix);
      if (n < 0) error->all("Could not find dump custom compute ID");

      if (modify->compute[n]->peratom_flag == 0)
	error->all("Dump custom compute ID does not compute peratom info");
      if (arg_compute[nfield+nthresh] == 0 && 
	  modify->compute[n]->size_peratom > 0)
	error->all("Dump custom compute ID does not compute scalar per atom");
      if (arg_compute[nfield+nthresh] > 0 && 
	  modify->compute[n]->size_peratom == 0)
	error->all("Dump custom compute ID does not compute vector per atom");
      if (arg_compute[nfield+nthresh] > 0 && 
	  arg_compute[nfield+nthresh] > modify->compute[n]->size_peratom)
	error->all("Dump custom compute ID vector is not large enough");
      if (modify->compute[n]->id_pre)
	int tmp = add_compute(modify->compute[n]->id_pre,0);
      field2compute[nfield+nthresh] = add_compute(suffix,0);
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
  }
  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf and choose and local arrays
------------------------------------------------------------------------- */

int DumpCustom::memory_usage()
{
  int bytes = maxbuf * sizeof(double);
  bytes += maxlocal * sizeof(int);
  bytes += maxlocal * sizeof(double);
  return bytes;
}

// ----------------------------------------------------------------------
// one method for every keyword dump custom can output
// the atom quantity is packed into buf starting at n with stride size_one
// customize by adding a method
// ----------------------------------------------------------------------

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_compute(int n)
{
  double *vector = compute[field2compute[n]]->scalar_atom;
  double **array = compute[field2compute[n]]->vector_atom;
  int index = arg_compute[n];
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

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_epair(int n)
{
  double *epair = compute[index_epair]->scalar_atom;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = epair[i];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_ke(int n)
{
  double *ke = compute[index_ke]->scalar_atom;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = ke[i];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_etotal(int n)
{
  double *etotal = compute[index_etotal]->scalar_atom;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = etotal[i];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_centro(int n)
{
  double *centro = compute[index_centro]->scalar_atom;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = centro[i];
      n += size_one;
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_sxx(int n)
{
  double **stress = compute[index_stress]->vector_atom;
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
  double **stress = compute[index_stress]->vector_atom;
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
  double **stress = compute[index_stress]->vector_atom;
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
  double **stress = compute[index_stress]->vector_atom;
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
  double **stress = compute[index_stress]->vector_atom;
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
  double **stress = compute[index_stress]->vector_atom;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (choose[i]) {
      buf[n] = stress[i][5];
      n += size_one;
    }
}
