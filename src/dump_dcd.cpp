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
                        Axel Kohlmeyer (Temple U), support for groups
------------------------------------------------------------------------- */

#include <cmath>
#include <inttypes.h> // <cinttypes> requires C++-11
#include <cstdio>
#include <ctime>
#include <cstring>
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

DumpDCD::DumpDCD(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg),
  coords(NULL)
{
  if (narg != 5) error->all(FLERR,"Illegal dump dcd command");
  if (binary || compressed || multifile || multiproc)
    error->all(FLERR,"Invalid dump dcd filename");

  size_one = 3;
  sort_flag = 1;
  sortcol = 0;

  unwrap_flag = 0;
  format_default = NULL;

  // allocate global array for atom coords

  bigint n = group->count(igroup);
  if (n > static_cast<bigint>(MAXSMALLINT/3/sizeof(float)))
    error->all(FLERR,"Too many atoms for dump dcd");
  natoms = static_cast<int> (n);

  memory->create(coords,3*natoms,"dump:coords");
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
  memory->destroy(coords);
}

/* ---------------------------------------------------------------------- */

void DumpDCD::init_style()
{
  if (sort_flag == 0 || sortcol != 0)
    error->all(FLERR,"Dump dcd requires sorting by atom ID");

  // check that dump frequency has not changed and is not a variable

  int idump;
  for (idump = 0; idump < output->ndump; idump++)
    if (strcmp(id,output->dump[idump]->id) == 0) break;
  if (output->every_dump[idump] == 0)
    error->all(FLERR,"Cannot use variable every setting for dump dcd");

  if (nevery_save == 0) nevery_save = output->every_dump[idump];
  else if (nevery_save != output->every_dump[idump])
    error->all(FLERR,"Cannot change dump_modify every for dump dcd");
}

/* ---------------------------------------------------------------------- */

void DumpDCD::openfile()
{
  if (me == 0) {
    fp = fopen(filename,"wb");
    if (fp == NULL) error->one(FLERR,"Cannot open dump file");
  }
}

/* ---------------------------------------------------------------------- */

void DumpDCD::write_header(bigint n)
{
  if (n != natoms) error->all(FLERR,"Dump dcd of non-matching # of atoms");
  if (update->ntimestep > MAXSMALLINT)
    error->one(FLERR,"Too big a timestep for dump dcd");

  // first time, write header for entire file

  if (headerflag == 0) {
    if (me == 0) write_dcd_header("Written by LAMMPS");
    headerflag = 1;
    nframes = 0;
  }

  // dim[] = size and angle cosines of orthogonal or triclinic box
  // dim[0] = a = length of unit cell vector along x-axis
  // dim[1] = gamma = cosine of angle between a and b
  // dim[2] = b = length of unit cell vector in xy-plane
  // dim[3] = beta = cosine of angle between a and c
  // dim[4] = alpha = cosine of angle between b and c
  // dim[5] = c = length of final unit cell vector
  // 48 = 6 doubles

  double dim[6];
  if (domain->triclinic) {
    double *h = domain->h;
    double alen = h[0];
    double blen = sqrt(h[5]*h[5] + h[1]*h[1]);
    double clen = sqrt(h[4]*h[4] + h[3]*h[3] + h[2]*h[2]);
    dim[0] = alen;
    dim[2] = blen;
    dim[5] = clen;
    dim[4] = (h[5]*h[4] + h[1]*h[3]) / blen/clen; // alpha
    dim[3] = (h[0]*h[4]) / alen/clen;             // beta
    dim[1] = (h[0]*h[5]) / alen/blen;             // gamma
  } else {
    dim[0] = domain->xprd;
    dim[2] = domain->yprd;
    dim[5] = domain->zprd;
    dim[1] = dim[3] = dim[4] = 0.0;
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

void DumpDCD::pack(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  m = n = 0;
  if (unwrap_flag) {
    double xprd = domain->xprd;
    double yprd = domain->yprd;
    double zprd = domain->zprd;
    double xy = domain->xy;
    double xz = domain->xz;
    double yz = domain->yz;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        int ix = (image[i] & IMGMASK) - IMGMAX;
        int iy = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        int iz = (image[i] >> IMG2BITS) - IMGMAX;

        if (domain->triclinic) {
          buf[m++] = x[i][0] + ix * xprd + iy * xy + iz * xz;
          buf[m++] = x[i][1] + iy * yprd + iz * yz;
          buf[m++] = x[i][2] + iz * zprd;
        } else {
          buf[m++] = x[i][0] + ix * xprd;
          buf[m++] = x[i][1] + iy * yprd;
          buf[m++] = x[i][2] + iz * zprd;
        }
        ids[n++] = tag[i];
      }
    }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        buf[m++] = x[i][0];
        buf[m++] = x[i][1];
        buf[m++] = x[i][2];
        ids[n++] = tag[i];
      }
  }
}

/* ---------------------------------------------------------------------- */

void DumpDCD::write_data(int n, double *mybuf)
{
  // copy buf atom coords into 3 global arrays

  int m = 0;
  for (int i = 0; i < n; i++) {
    xf[ntotal] = mybuf[m++];
    yf[ntotal] = mybuf[m++];
    zf[ntotal] = mybuf[m++];
    ntotal++;
  }

  // if last chunk of atoms in this snapshot, write global arrays to file

  if (ntotal == natoms) {
    write_frame();
    ntotal = 0;
  }
}

/* ---------------------------------------------------------------------- */

int DumpDCD::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"unwrap") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"yes") == 0) unwrap_flag = 1;
    else if (strcmp(arg[1],"no") == 0) unwrap_flag = 0;
    else error->all(FLERR,"Illegal dump_modify command");
    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf and global coords array
------------------------------------------------------------------------- */

bigint DumpDCD::memory_usage()
{
  bigint bytes = Dump::memory_usage();
  bytes += memory->usage(coords,natoms*3);
  return bytes;
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

  int ntimestep = update->ntimestep;

  out_integer = 84;
  fwrite_int32(fp,out_integer);
  strcpy(title_string,"CORD");
  fwrite(title_string,4,1,fp);
  fwrite_int32(fp,0);                    // NFILE = # of snapshots in file
  fwrite_int32(fp,ntimestep);            // START = timestep of first snapshot
  fwrite_int32(fp,nevery_save);          // SKIP = interval between snapshots
  fwrite_int32(fp,ntimestep);            // NSTEP = timestep of last snapshot
  fwrite_int32(fp,0);                         // NAMD writes NSTEP or ISTART
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
  strftime(title_string,80,"REMARKS Created %d %B,%Y at %H:%M",tmbuf);
  fwrite(title_string,80,1,fp);
  fwrite_int32(fp,164);
  fwrite_int32(fp,4);
  fwrite_int32(fp,natoms);                // number of atoms in each snapshot
  fwrite_int32(fp,4);
  if (flush_flag) fflush(fp);
}
