// -*- c++ -*-

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
   Contributing author:  Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "ndx_group.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "group.h"
#include "memory.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;
#define BUFLEN 4096
#define DELTA 16384

static char *find_section(FILE *fp, const char *name)
{
  char linebuf[BUFLEN];
  char *n,*p,*t,*r;

  while ((p = fgets(linebuf,BUFLEN,fp))) {
    t = strtok(p," \t\n\r\f");
    if ((t != NULL) && *t == '[') {
      t = strtok(NULL," \t\n\r\f");
      if (t != NULL) {
        n = t;
        t = strtok(NULL," \t\n\r\f");
        if ((t != NULL) && *t == ']') {
          if ((name == NULL) || strcmp(name,n) == 0) {
            int l = strlen(n);
            r = new char[l+1];
            strncpy(r,n,l+1);
            return r;
          }
        }
      }
    }
  }
  return NULL;
}

static tagint *read_section(FILE *fp, bigint &num)
{
  char linebuf[BUFLEN];
  char *p,*t;
  tagint *tagbuf;
  bigint nmax;

  num = 0;
  nmax = DELTA;
  tagbuf = (tagint *)malloc(sizeof(tagint)*nmax);

  while ((p = fgets(linebuf,BUFLEN,fp))) {
    t = strtok(p," \t\n\r\f");
    while (t != NULL) {
      // start of a new section. we are done here.
      if (*t == '[') return tagbuf;

      tagbuf[num++] = ATOTAGINT(t);
      if (num == nmax) {
        nmax += DELTA;
        tagbuf = (tagint *)realloc(tagbuf,sizeof(tagint)*nmax);
      }
      t = strtok(NULL," \t\n\r\f");
    }
  }
  return tagbuf;
}

/* ---------------------------------------------------------------------- */

void Ndx2Group::command(int narg, char **arg)
{
  int len;
  bigint num;
  FILE *fp;
  char *name = NULL;
  tagint *tags;

  if (narg < 1) error->all(FLERR,"Illegal ndx2group command");

  if (atom->tag_enable == 0)
      error->all(FLERR,"Must have atom IDs for ndx2group command");

  if (comm->me == 0) {
    fp = fopen(arg[0], "r");
    if (fp == NULL)
      error->one(FLERR,"Cannot open index file for reading");

    if (screen)
      fprintf(screen, "Reading groups from index file %s:\n",arg[0]);
    if (logfile)
      fprintf(logfile,"Reading groups from index file %s:\n",arg[0]);
  }

  if (narg == 1) {    // restore all groups

    do {
      if (comm->me == 0) {
        len = 0;

        // find the next section.
        // if we had processed a section, before we need to step back
        if (name != NULL) {
          rewind(fp);
          char *tmp = find_section(fp,name);
          delete[] tmp;
          delete[] name;
          name = NULL;
        }
        name = find_section(fp,NULL);
        if (name != NULL) {
          len=strlen(name);

          // skip over group "all", which is called "System" in gromacs
          if (strcmp(name,"System") == 0) continue;

          if (screen)
            fprintf(screen," Processing group '%s'\n",name);
          if (logfile)
            fprintf(logfile," Processing group '%s'\n",name);
        }
        MPI_Bcast(&len,1,MPI_INT,0,world);
        if (len > 0) {
          MPI_Bcast(name,len,MPI_CHAR,0,world);

          // read tags for atoms in group and broadcast
          num = 0;
          tags = read_section(fp,num);
          MPI_Bcast(&num,1,MPI_LMP_BIGINT,0,world);
          MPI_Bcast(tags,num,MPI_LMP_TAGINT,0,world);
          create(name,num,tags);
          free(tags);
        }
      } else {
        MPI_Bcast(&len,1,MPI_INT,0,world);
        if (len > 0) {
          delete[] name;
          name = new char[len+1];
          MPI_Bcast(name,len+1,MPI_CHAR,0,world);

          MPI_Bcast(&num,1,MPI_LMP_BIGINT,0,world);
          tags = (tagint *)malloc(sizeof(tagint)*(num ? num : 1));
          MPI_Bcast(tags,num,MPI_LMP_TAGINT,0,world);
          create(name,num,tags);
          free(tags);
        }
      }
    } while (len);

  } else {            // restore selected groups
    for (int idx=1; idx < narg; ++idx) {

      if (comm->me == 0) {
        len = 0;

        // find named section, search from beginning of file
        if (name != NULL) delete[] name;
        rewind(fp);
        name = find_section(fp,arg[idx]);
        if (name != NULL) len=strlen(name);

        if (screen)
          fprintf(screen," %s group '%s'\n",
                  len ? "Processing" : "Skipping",arg[idx]);
        if (logfile)
          fprintf(logfile,"%s group '%s'\n",
                  len ? "Processing" : "Skipping",arg[idx]);

        MPI_Bcast(&len,1,MPI_INT,0,world);
        if (len > 0) {
          MPI_Bcast(name,len+1,MPI_CHAR,0,world);
          // read tags for atoms in group and broadcast
          num = 0;
          tags = read_section(fp,num);
          MPI_Bcast(&num,1,MPI_LMP_BIGINT,0,world);
          MPI_Bcast(tags,num,MPI_LMP_TAGINT,0,world);
          create(name,num,tags);
          free(tags);
        }
      } else {

        MPI_Bcast(&len,1,MPI_INT,0,world);
        if (len > 0) {
          delete[] name;
          name = new char[len+1];
          MPI_Bcast(name,len+1,MPI_CHAR,0,world);

          MPI_Bcast(&num,1,MPI_LMP_BIGINT,0,world);
          tags = (tagint *)malloc(sizeof(tagint)*(num ? num : 1));
          MPI_Bcast(tags,num,MPI_LMP_TAGINT,0,world);
          create(name,num,tags);
          free(tags);
        }
      }
    }
  }

  delete[] name;
  if (comm->me == 0) {
    if (screen) fputs("\n",screen);
    if (logfile) fputs("\n",logfile);
    fclose(fp);
  }
}

/* ---------------------------------------------------------------------- */

void Ndx2Group::create(char *name, bigint num, tagint *tags)
{
  // wipe out all members if the group exists. gid==0 is group "all"
  int gid = group->find(name);
  if (gid > 0) {
    char *cmd[2];
    cmd[0] = name;
    cmd[1] = (char *)"clear";
    group->assign(2,cmd);
  }

  // map from global to local
  const int nlocal = atom->nlocal;
  int *flags = (int *)calloc(nlocal,sizeof(int));
  for (bigint i=0; i < num; ++i) {
    const int id = atom->map(tags[i]);
    if (id < nlocal && id >= 0)
      flags[id] = 1;
  }
  group->create(name,flags);
  free(flags);
}

