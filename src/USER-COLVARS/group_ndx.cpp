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

#include "group_ndx.h"
#include "atom.h"
#include "comm.h"
#include "group.h"
#include "memory.h"
#include "error.h"

#include <stdio.h>
#include <stdlib.h>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   helper function. integer comparison for qsort()
   ---------------------------------------------------------------------- */

static int cmptagint(const void *p1, const void *p2)
{
  const tagint i1 = * static_cast<const tagint *>(p1);
  const tagint i2 = * static_cast<const tagint *>(p2);
  if (i1 == i2) return 0;
  else {
    if (i1 < i2) return -1;
    else return 1;
  }
}

/* ----------------------------------------------------------------------
   helper function. writes out one group to a gromacs style index file
   ---------------------------------------------------------------------- */

static void write_group(FILE *fp, int gid, Atom *atom, Group *group, int me,
                        int np, MPI_Comm world, FILE *screen, FILE *logfile)
{
  char fmt[8];
  tagint *sendlist, *recvlist;
  bigint num = group->count(gid);
  int lnum, cols;

  if (me == 0) {
    if (screen) fprintf(screen, " writing group %s... ", group->names[gid]);
    if (logfile) fprintf(logfile, " writing group %s... ", group->names[gid]);

    // the "all" group in LAMMPS is called "System" in gromacs
    if (gid == 0) {
      fputs("[ System ]\n", fp);
    } else {
      fprintf(fp,"[ %s ]\n", group->names[gid]);
    }

    // derive format string for index lists
    bigint j = atom->natoms;
    int i = 0;
    while (j > 0) {
      ++i;
      j /= 10;
    }
    sprintf(fmt,"%%%dd ", i);
    cols = 80 / (i+1);
  }

  if (num > 0) {
    const int * const mask = atom->mask;
    const tagint * const tag = atom->tag;
    const int groupbit = group->bitmask[gid];
    const int nlocal = atom->nlocal;
    int i,j;

    sendlist = new tagint[nlocal];
    recvlist = new tagint[num];
    lnum = 0;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) sendlist[lnum++] = tag[i];

    MPI_Status status;
    MPI_Request request;
    int nrecv,allrecv;
    if (me == 0) {
      for (i = 0; i < lnum; i++)
        recvlist[i] = sendlist[i];

      allrecv = lnum;
      for (int i=1; i < np; ++i) {
        MPI_Irecv(recvlist+allrecv,num-allrecv,MPI_LMP_TAGINT,i,0, world,&request);
        MPI_Send(&nrecv,0,MPI_INT,i,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&nrecv);
        allrecv += nrecv;
      }

      // sort received list
      qsort((void *)recvlist, num, sizeof(tagint), cmptagint);
    } else {
      MPI_Recv(&nrecv,0,MPI_INT,0,0,world,&status);
      MPI_Rsend(sendlist,lnum,MPI_LMP_TAGINT,0,0,world);
    }
    delete [] sendlist;
  }

  if (me == 0) {
    int i, j;
    for(i = 0, j = 0; i < num; ++i) {
      fprintf(fp,fmt,recvlist[i]);
      ++j;
      if (j == cols) {
        fputs("\n",fp);
        j = 0;
      }
    }
    if (j > 0) fputs("\n",fp);
    if (screen) fputs("done\n",screen);
    if (logfile) fputs("done\n",logfile);
  }
  if (num > 0) delete[] recvlist;
}

/* ---------------------------------------------------------------------- */

void Group2Ndx::command(int narg, char **arg)
{
  FILE *fp;

  if (narg < 1) error->all(FLERR,"Illegal group2ndx command");

  if (atom->tag_enable == 0)
      error->all(FLERR,"Must have atom IDs for group2ndx command");

  if (comm->me == 0) {
    fp = fopen(arg[0], "w");
    if (fp == NULL)
      error->one(FLERR,"Cannot open index file for writing");

    if (screen)
      fprintf(screen, "Writing groups to index file %s:\n",arg[0]);
    if (logfile)
      fprintf(logfile,"Writing groups to index file %s:\n",arg[0]);
  }

  if (narg == 1) { // write out all groups
    for (int i=0; i < group->ngroup; ++i) {
      write_group(fp,i,atom,group,comm->me,comm->nprocs,world,screen,logfile);
    }

  } else { // write only selected groups
    for (int i=1; i < narg; ++i) {
      int gid = group->find(arg[i]);
      if (gid < 0) error->all(FLERR, "Non-existing group requested");
      write_group(fp,gid,atom,group,comm->me,comm->nprocs,world,screen,logfile);
    }
  }

  if (comm->me == 0) {
    if (screen) fputs("\n",screen);
    if (logfile) fputs("\n",logfile);
    fclose(fp);
  }
}

