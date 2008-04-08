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
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
------------------------------------------------------------------------- */

#include "math.h"
#include "inttypes.h"
#include "stdio.h"
#include "time.h"
#include "string.h"
#include "dump_dcd.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "output.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define NFILE_POS 8L
#define NSTEP_POS 20L

// necessary to set SEEK params b/c MPI-2 messes with these settings

#ifndef SEEK_SET
#define SEEK_SET        0
#define SEEK_CUR        1
#define SEEK_END        2
#endif

/* ---------------------------------------------------------------------- */

static inline void fwrite_int32(FILE* fd, uint32_t i)
{
  fwrite(&i,sizeof(uint32_t),1,fd);
}

/* ---------------------------------------------------------------------- */

DumpDCD::DumpDCD(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5) error->all("Illegal dump dcd command");
  if (igroup != group->find("all")) error->all("Dump dcd must use group all");
  if (binary || compressed || multifile || multiproc)
    error->all("Invalid dump dcd filename");

  size_one = 4;
  unwrap_flag = 0;
  format_default = NULL;
    
  // allocate global array for atom coords

  natoms = static_cast<int> (atom->natoms);
  if (natoms <= 0) error->all("Invalid natoms for dump dcd");
  if (atom->tag_consecutive() == 0)
    error->all("Atom IDs must be consecutive for dump dcd");
  coords = (float *) memory->smalloc(3*natoms*sizeof(float),"dump:coords");
  xf = &coords[0*natoms];
  yf = &coords[1*natoms];
  zf = &coords[2*natoms];

  openfile();
  headerflag = 0;
  nevery_save = 0;
  ntotal = 0;
}

/* ---------------------------------------------------------------------- */

DumpDCD::~DumpDCD()
{
  memory->sfree(coords);
}

/* ---------------------------------------------------------------------- */

void DumpDCD::init()
{
  if (unwrap_flag == 1 && domain->triclinic)
    error->all("Dump dcd cannot dump unwrapped coords with triclinic box");
  
  // check that dump frequency has not changed

  if (nevery_save == 0) {
    int idump;
    for (idump = 0; idump < output->ndump; idump++)
      if (strcmp(id,output->dump[idump]->id) == 0) break;
    nevery_save = output->dump_every[idump];
  } else {
    int idump;
    for (idump = 0; idump < output->ndump; idump++)
      if (strcmp(id,output->dump[idump]->id) == 0) break;
    if (nevery_save != output->dump_every[idump])
      error->all("Cannot change dump_modify every for dump dcd");
  }
}

/* ---------------------------------------------------------------------- */

int DumpDCD::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"unwrap") == 0) {
    if (narg < 2) error->all("Illegal dump_modify command");
    if (strcmp(arg[1],"yes") == 0) unwrap_flag = 1;
    else if (strcmp(arg[1],"no") == 0) unwrap_flag = 0;
    else error->all("Illegal dump_modify command");
    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf and global coords array
------------------------------------------------------------------------- */

double DumpDCD::memory_usage()
{
  double bytes = maxbuf * sizeof(double);
  bytes += 3*natoms * sizeof(float);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void DumpDCD::openfile()
{
  if (me == 0) {
    fp = fopen(filename,"wb");
    if (fp == NULL) error->one("Cannot open dump file");
  }
}

/* ---------------------------------------------------------------------- */

void DumpDCD::write_header(int n)
{
  if (n != natoms) error->all("Dump dcd of non-matching # of atoms");

  // first time, write header for entire file

  if (headerflag == 0) {
    if (me == 0) write_dcd_header("Written by LAMMPS");
    headerflag = 1;
    nframes = 0;
  }

  // dim[] = size and angle cosines of orthogonal or triclinic box
  // dim[1] = alpha = cosine of angle between b and c
  // dim[3] = beta = cosine of angle between a and c
  // dim[4] = gamma = cosine of angle between a and b
  // 48 = 6 doubles

  double dim[6];
  dim[0] = domain->xprd;
  dim[2] = domain->yprd;
  dim[5] = domain->zprd;
  if (domain->triclinic == 0) dim[1] = dim[3] = dim[4] = 0.0;
  else {
    double *h = domain->h;
    double alen = h[0];
    double blen = sqrt(h[5]*h[5] + h[1]*h[1]);
    double clen = sqrt(h[4]*h[4] + h[3]*h[3] + h[2]*h[2]);
    dim[1] = (h[5]*h[4] + h[1]*h[3]) / blen/clen;
    dim[3] = (h[0]*h[4]) / alen/clen;
    dim[4] = (h[0]*h[5]) / alen/blen;
  }
  
  if (me == 0) {
    uint32_t out_integer = 48;
    fwrite_int32(fp,out_integer);
    fwrite(dim,out_integer,1,fp); 
    fwrite_int32(fp,out_integer);
    if (flush_flag) fflush(fp);
  }
}

/* ---------------------------------------------------------------------- */

int DumpDCD::count()
{
  return atom->nlocal;
}

/* ---------------------------------------------------------------------- */

int DumpDCD::pack()
{
  int *tag = atom->tag;
  double **x = atom->x;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  // assume group all, so no need to perform mask check

  int m = 0;
  if (unwrap_flag) {
    double xprd = domain->xprd;
    double yprd = domain->yprd;
    double zprd = domain->zprd;

    for (int i = 0; i < nlocal; i++) {
      buf[m++] = tag[i];
      buf[m++] = x[i][0] + ((image[i] & 1023) - 512) * xprd;
      buf[m++] = x[i][1] + ((image[i] >> 10 & 1023) - 512) * yprd;
      buf[m++] = x[i][2] + ((image[i] >> 20) - 512) * zprd;
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      buf[m++] = tag[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void DumpDCD::write_data(int n, double *buf)
{
  // spread buf atom coords into global arrays

  int tag;
  int m = 0;
  for (int i = 0; i < n; i++) {
    tag = static_cast<int> (buf[m]) - 1;
    xf[tag] = buf[m+1];
    yf[tag] = buf[m+2];
    zf[tag] = buf[m+3];
    m += size_one;
  }

  // if last chunk of atoms in this snapshot, write global arrays to file

  ntotal += n;
  if (ntotal == natoms) {
    write_frame();
    ntotal = 0;
  }
}

/* ---------------------------------------------------------------------- */

void DumpDCD::write_frame()
{
  // write coords

  uint32_t out_integer = natoms*sizeof(float);
  fwrite_int32(fp,out_integer);
  fwrite(xf,out_integer,1,fp);
  fwrite_int32(fp,out_integer);
  fwrite_int32(fp,out_integer);
  fwrite(yf,out_integer,1,fp);
  fwrite_int32(fp,out_integer);
  fwrite_int32(fp,out_integer);
  fwrite(zf,out_integer,1,fp);
  fwrite_int32(fp,out_integer);
  
  // update NFILE and NSTEP fields in DCD header

  nframes++;
  out_integer = nframes;
  fseek(fp,NFILE_POS,SEEK_SET);
  fwrite_int32(fp,out_integer);
  out_integer = update->ntimestep;
  fseek(fp,NSTEP_POS,SEEK_SET);
  fwrite_int32(fp,out_integer);
  fseek(fp,0,SEEK_END);
}

/* ---------------------------------------------------------------------- */

void DumpDCD::write_dcd_header(const char *remarks)
{
  uint32_t out_integer;
  float out_float;
  char title_string[200];
  time_t cur_time;
  struct tm *tmbuf;

  out_integer = 84;
  fwrite_int32(fp,out_integer);
  strcpy(title_string,"CORD");
  fwrite(title_string,4,1,fp);
  fwrite_int32(fp,0);                    // NFILE = # of snapshots in file
  fwrite_int32(fp,update->ntimestep);    // START = timestep of first snapshot
  fwrite_int32(fp,nevery_save);          // SKIP = interval between snapshots
  fwrite_int32(fp,update->ntimestep);    // NSTEP = timestep of last snapshot
  fwrite_int32(fp,0);			 // NAMD writes NSTEP or ISTART
  fwrite_int32(fp,0);
  fwrite_int32(fp,0);
  fwrite_int32(fp,0);
  fwrite_int32(fp,0);
  out_float = update->dt;
  fwrite(&out_float,sizeof(float),1,fp);
  fwrite_int32(fp,1);
  fwrite_int32(fp,0);
  fwrite_int32(fp,0);
  fwrite_int32(fp,0);
  fwrite_int32(fp,0);
  fwrite_int32(fp,0);
  fwrite_int32(fp,0);
  fwrite_int32(fp,0);
  fwrite_int32(fp,0);
  fwrite_int32(fp,24);                   // pretend to be Charmm version 24
  fwrite_int32(fp,84);
  fwrite_int32(fp,164);
  fwrite_int32(fp,2);
  strncpy(title_string,remarks,80);
  title_string[79] = '\0';
  fwrite(title_string,80,1,fp);
  cur_time=time(NULL);
  tmbuf=localtime(&cur_time);
  memset(title_string,' ',81);
  strftime(title_string,80,"REMARKS Created %d %B,%Y at %R",tmbuf);
  fwrite(title_string,80,1,fp);
  fwrite_int32(fp,164);
  fwrite_int32(fp,4);
  fwrite_int32(fp,natoms);                // number of atoms in each snapshot
  fwrite_int32(fp,4);
  if (flush_flag) fflush(fp);
}
