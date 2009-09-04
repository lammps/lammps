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

#include "string.h"
#include "dump_atom.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpAtom::DumpAtom(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5) error->all("Illegal dump atom command");

  scale_flag = 1;
  image_flag = 0;
  format_default = NULL;
}

/* ---------------------------------------------------------------------- */

void DumpAtom::init()
{
  if (image_flag == 0) size_one = 5;
  else size_one = 8;

  // default format depends on image flags

  delete [] format;
  if (format_user) {
    int n = strlen(format_user) + 2;
    format = new char[n];
    strcpy(format,format_user);
    strcat(format,"\n");
  } else {
    char *str;
    if (image_flag == 0) str = (char *) "%d %d %g %g %g";
    else str = (char *) "%d %d %g %g %g %d %d %d";
    int n = strlen(str) + 2;
    format = new char[n];
    strcpy(format,str);
    strcat(format,"\n");
  }

  // setup column string

  if (scale_flag == 0 && image_flag == 0)
    columns = "id type x y z";
  else if (scale_flag == 0 && image_flag == 1)
    columns = "id type x y z ix iy iz";
  else if (scale_flag == 1 && image_flag == 0)
    columns = "id type xs ys zs";
  else if (scale_flag == 1 && image_flag == 1)
    columns = "id type xs ys zs ix iy iz";

  // setup function ptrs

  if (binary && domain->triclinic == 0)
    header_choice = &DumpAtom::header_binary;
  else if (binary && domain->triclinic == 1)
    header_choice = &DumpAtom::header_binary_triclinic;
  else if (!binary && domain->triclinic == 0)
    header_choice = &DumpAtom::header_item;
  else if (!binary && domain->triclinic == 1)
    header_choice = &DumpAtom::header_item_triclinic;

  if (scale_flag == 1 && image_flag == 0 && domain->triclinic == 0)
    pack_choice = &DumpAtom::pack_scale_noimage;
  else if (scale_flag == 1 && image_flag == 1 && domain->triclinic == 0)
    pack_choice = &DumpAtom::pack_scale_image;
  else if (scale_flag == 1 && image_flag == 0 && domain->triclinic == 1)
    pack_choice = &DumpAtom::pack_scale_noimage_triclinic;
  else if (scale_flag == 1 && image_flag == 1 && domain->triclinic == 1)
    pack_choice = &DumpAtom::pack_scale_image_triclinic;
  else if (scale_flag == 0 && image_flag == 0)
    pack_choice = &DumpAtom::pack_noscale_noimage;
  else if (scale_flag == 0 && image_flag == 1)
    pack_choice = &DumpAtom::pack_noscale_image;

  if (binary) write_choice = &DumpAtom::write_binary;
  else if (image_flag == 0) write_choice = &DumpAtom::write_noimage;
  else if (image_flag == 1) write_choice = &DumpAtom::write_image;

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

int DumpAtom::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"scale") == 0) {
    if (narg < 2) error->all("Illegal dump_modify command");
    if (strcmp(arg[1],"yes") == 0) scale_flag = 1;
    else if (strcmp(arg[1],"no") == 0) scale_flag = 0;
    else error->all("Illegal dump_modify command");
    return 2;
  } else if (strcmp(arg[0],"image") == 0) {
    if (narg < 2) error->all("Illegal dump_modify command");
    if (strcmp(arg[1],"yes") == 0) image_flag = 1;
    else if (strcmp(arg[1],"no") == 0) image_flag = 0;
    else error->all("Illegal dump_modify command");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_header(int ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

int DumpAtom::count()
{
  if (igroup == 0) return atom->nlocal;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) m++;
  return m;
}

/* ---------------------------------------------------------------------- */

int DumpAtom::pack()
{
  return (this->*pack_choice)();
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_data(int n, double *buf)
{
  (this->*write_choice)(n,buf);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::header_binary(int ndump)
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

void DumpAtom::header_binary_triclinic(int ndump)
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

void DumpAtom::header_item(int ndump)
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

void DumpAtom::header_item_triclinic(int ndump)
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

int DumpAtom::pack_scale_image()
{
  int *tag = atom->tag;
  int *type = atom->type;
  int *image = atom->image;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double invxprd = 1.0/domain->xprd;
  double invyprd = 1.0/domain->yprd;
  double invzprd = 1.0/domain->zprd;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = (x[i][0] - boxxlo) * invxprd;
      buf[m++] = (x[i][1] - boxylo) * invyprd;
      buf[m++] = (x[i][2] - boxzlo) * invzprd;
      buf[m++] = (image[i] & 1023) - 512;
      buf[m++] = (image[i] >> 10 & 1023) - 512;
      buf[m++] = (image[i] >> 20) - 512;
    }

  return m;
}

/* ---------------------------------------------------------------------- */

int DumpAtom::pack_scale_noimage()
{
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double invxprd = 1.0/domain->xprd;
  double invyprd = 1.0/domain->yprd;
  double invzprd = 1.0/domain->zprd;
  
  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = (x[i][0] - boxxlo) * invxprd;
      buf[m++] = (x[i][1] - boxylo) * invyprd;
      buf[m++] = (x[i][2] - boxzlo) * invzprd;
    }

  return m;
}

/* ---------------------------------------------------------------------- */

int DumpAtom::pack_scale_image_triclinic()
{
  int *tag = atom->tag;
  int *type = atom->type;
  int *image = atom->image;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double lamda[3];

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      domain->x2lamda(x[i],lamda);
      buf[m++] = lamda[0];
      buf[m++] = lamda[1];
      buf[m++] = lamda[2];
      buf[m++] = (image[i] & 1023) - 512;
      buf[m++] = (image[i] >> 10 & 1023) - 512;
      buf[m++] = (image[i] >> 20) - 512;
    }

  return m;
}

/* ---------------------------------------------------------------------- */

int DumpAtom::pack_scale_noimage_triclinic()
{
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double lamda[3];
  
  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      domain->x2lamda(x[i],lamda);
      buf[m++] = lamda[0];
      buf[m++] = lamda[1];
      buf[m++] = lamda[2];
    }

  return m;
}

/* ---------------------------------------------------------------------- */

int DumpAtom::pack_noscale_image()
{
  int *tag = atom->tag;
  int *type = atom->type;
  int *image = atom->image;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      buf[m++] = (image[i] & 1023) - 512;
      buf[m++] = (image[i] >> 10 & 1023) - 512;
      buf[m++] = (image[i] >> 20) - 512;
    }

  return m;
}

/* ---------------------------------------------------------------------- */

int DumpAtom::pack_noscale_noimage()
{
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
    }

  return m;
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_binary(int n, double *buf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(buf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_image(int n, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
	    static_cast<int> (buf[m]), static_cast<int> (buf[m+1]),
	    buf[m+2],buf[m+3],buf[m+4], static_cast<int> (buf[m+5]),
	    static_cast<int> (buf[m+6]), static_cast<int> (buf[m+7]));
    m += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_noimage(int n, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
	    static_cast<int> (buf[m]), static_cast<int> (buf[m+1]),
	    buf[m+2],buf[m+3],buf[m+4]);
    m += size_one;
  }
}
