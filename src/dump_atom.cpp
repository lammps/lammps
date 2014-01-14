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
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define ONELINE 256
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpAtom::DumpAtom(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal dump atom command");

  scale_flag = 1;
  image_flag = 0;
  buffer_allow = 1;
  buffer_flag = 1;
  format_default = NULL;
}

/* ---------------------------------------------------------------------- */

void DumpAtom::init_style()
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

  // setup boundary string

  domain->boundary_string(boundstr);

  // setup column string

  if (scale_flag == 0 && image_flag == 0)
    columns = (char *) "id type x y z";
  else if (scale_flag == 0 && image_flag == 1)
    columns = (char *) "id type x y z ix iy iz";
  else if (scale_flag == 1 && image_flag == 0)
    columns = (char *) "id type xs ys zs";
  else if (scale_flag == 1 && image_flag == 1)
    columns = (char *) "id type xs ys zs ix iy iz";

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

  if (image_flag == 0) convert_choice = &DumpAtom::convert_noimage;
  else convert_choice = &DumpAtom::convert_image;

  if (binary) write_choice = &DumpAtom::write_binary;
  else if (buffer_flag == 1) write_choice = &DumpAtom::write_string;
  else if (image_flag == 0) write_choice = &DumpAtom::write_lines_noimage;
  else if (image_flag == 1) write_choice = &DumpAtom::write_lines_image;

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

int DumpAtom::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"scale") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"yes") == 0) scale_flag = 1;
    else if (strcmp(arg[1],"no") == 0) scale_flag = 0;
    else error->all(FLERR,"Illegal dump_modify command");
    return 2;
  } else if (strcmp(arg[0],"image") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"yes") == 0) image_flag = 1;
    else if (strcmp(arg[1],"no") == 0) image_flag = 0;
    else error->all(FLERR,"Illegal dump_modify command");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::pack(int *ids)
{
  (this->*pack_choice)(ids);
}

/* ---------------------------------------------------------------------- */

int DumpAtom::convert_string(int n, double *mybuf)
{
  return (this->*convert_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::header_binary(bigint ndump)
{
  fwrite(&update->ntimestep,sizeof(bigint),1,fp);
  fwrite(&ndump,sizeof(bigint),1,fp);
  fwrite(&domain->triclinic,sizeof(int),1,fp);
  fwrite(&domain->boundary[0][0],6*sizeof(int),1,fp);
  fwrite(&boxxlo,sizeof(double),1,fp);
  fwrite(&boxxhi,sizeof(double),1,fp);
  fwrite(&boxylo,sizeof(double),1,fp);
  fwrite(&boxyhi,sizeof(double),1,fp);
  fwrite(&boxzlo,sizeof(double),1,fp);
  fwrite(&boxzhi,sizeof(double),1,fp);
  fwrite(&size_one,sizeof(int),1,fp);
  if (multiproc) fwrite(&nclusterprocs,sizeof(int),1,fp);
  else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::header_binary_triclinic(bigint ndump)
{
  fwrite(&update->ntimestep,sizeof(bigint),1,fp);
  fwrite(&ndump,sizeof(bigint),1,fp);
  fwrite(&domain->triclinic,sizeof(int),1,fp);
  fwrite(&domain->boundary[0][0],6*sizeof(int),1,fp);
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
  if (multiproc) fwrite(&nclusterprocs,sizeof(int),1,fp);
  else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::header_item(bigint ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
  fprintf(fp,"%g %g\n",boxxlo,boxxhi);
  fprintf(fp,"%g %g\n",boxylo,boxyhi);
  fprintf(fp,"%g %g\n",boxzlo,boxzhi);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::header_item_triclinic(bigint ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
  fprintf(fp,"%g %g %g\n",boxxlo,boxxhi,boxxy);
  fprintf(fp,"%g %g %g\n",boxylo,boxyhi,boxxz);
  fprintf(fp,"%g %g %g\n",boxzlo,boxzhi,boxyz);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::pack_scale_image(int *ids)
{
  int m,n;

  int *tag = atom->tag;
  int *type = atom->type;
  imageint *image = atom->image;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double invxprd = 1.0/domain->xprd;
  double invyprd = 1.0/domain->yprd;
  double invzprd = 1.0/domain->zprd;

  m = n = 00;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = (x[i][0] - boxxlo) * invxprd;
      buf[m++] = (x[i][1] - boxylo) * invyprd;
      buf[m++] = (x[i][2] - boxzlo) * invzprd;
      buf[m++] = (image[i] & IMGMASK) - IMGMAX;
      buf[m++] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      buf[m++] = (image[i] >> IMG2BITS) - IMGMAX;
      if (ids) ids[n++] = tag[i];
    }
}

/* ---------------------------------------------------------------------- */

void DumpAtom::pack_scale_noimage(int *ids)
{
  int m,n;

  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double invxprd = 1.0/domain->xprd;
  double invyprd = 1.0/domain->yprd;
  double invzprd = 1.0/domain->zprd;

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = (x[i][0] - boxxlo) * invxprd;
      buf[m++] = (x[i][1] - boxylo) * invyprd;
      buf[m++] = (x[i][2] - boxzlo) * invzprd;
      if (ids) ids[n++] = tag[i];
    }
}

/* ---------------------------------------------------------------------- */

void DumpAtom::pack_scale_image_triclinic(int *ids)
{
  int m,n;

  int *tag = atom->tag;
  int *type = atom->type;
  imageint *image = atom->image;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double lamda[3];

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      domain->x2lamda(x[i],lamda);
      buf[m++] = lamda[0];
      buf[m++] = lamda[1];
      buf[m++] = lamda[2];
      buf[m++] = (image[i] & IMGMASK) - IMGMAX;
      buf[m++] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      buf[m++] = (image[i] >> IMG2BITS) - IMGMAX;
      if (ids) ids[n++] = tag[i];
    }
}

/* ---------------------------------------------------------------------- */

void DumpAtom::pack_scale_noimage_triclinic(int *ids)
{
  int m,n;

  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double lamda[3];

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      domain->x2lamda(x[i],lamda);
      buf[m++] = lamda[0];
      buf[m++] = lamda[1];
      buf[m++] = lamda[2];
      if (ids) ids[n++] = tag[i];
    }
}

/* ---------------------------------------------------------------------- */

void DumpAtom::pack_noscale_image(int *ids)
{
  int m,n;

  int *tag = atom->tag;
  int *type = atom->type;
  imageint *image = atom->image;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      buf[m++] = (image[i] & IMGMASK) - IMGMAX;
      buf[m++] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      buf[m++] = (image[i] >> IMG2BITS) - IMGMAX;
      if (ids) ids[n++] = tag[i];
    }
}

/* ---------------------------------------------------------------------- */

void DumpAtom::pack_noscale_noimage(int *ids)
{
  int m,n;

  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      if (ids) ids[n++] = tag[i];
    }
}

/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpAtom::convert_image(int n, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }

    offset += sprintf(&sbuf[offset],format,
                      static_cast<int> (mybuf[m]), 
                      static_cast<int> (mybuf[m+1]),
                      mybuf[m+2],mybuf[m+3],mybuf[m+4], 
                      static_cast<int> (mybuf[m+5]),
                      static_cast<int> (mybuf[m+6]),
                      static_cast<int> (mybuf[m+7]));
    m += size_one;
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

int DumpAtom::convert_noimage(int n, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }

    offset += sprintf(&sbuf[offset],format,
                      static_cast<int> (mybuf[m]),
                      static_cast<int> (mybuf[m+1]),
                      mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += size_one;
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_binary(int n, double *mybuf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(mybuf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_string(int n, double *mybuf)
{
  fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_lines_image(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
            static_cast<int> (mybuf[m]), static_cast<int> (mybuf[m+1]),
            mybuf[m+2],mybuf[m+3],mybuf[m+4], static_cast<int> (mybuf[m+5]),
            static_cast<int> (mybuf[m+6]), static_cast<int> (mybuf[m+7]));
    m += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
            static_cast<int> (mybuf[m]), static_cast<int> (mybuf[m+1]),
            mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += size_one;
  }
}
