// clang-format off
// -*- c++ -*-

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
#include "error.h"
#include "group.h"

#include <cmath>

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

/* ---------------------------------------------------------------------- */

void Group2Ndx::command(int narg, char **arg)
{
  FILE *fp = nullptr;

  if (narg < 1) error->all(FLERR,"Illegal group2ndx command");

  if (atom->tag_enable == 0)
      error->all(FLERR,"Must have atom IDs for group2ndx command");

  if (comm->me == 0) {
    fp = fopen(arg[0], "w");
    if (fp == nullptr)
      error->one(FLERR,"Cannot open index file for writing: {}", utils::getsyserror());
    utils::logmesg(lmp,"Writing groups to index file {}:\n",arg[0]);
  }

  if (narg == 1) { // write out all groups
    for (int i=0; i < group->ngroup; ++i) {
      write_group(fp,i);
    }
  } else { // write only selected groups
    for (int i=1; i < narg; ++i) {
      int gid = group->find(arg[i]);
      if (gid < 0) error->all(FLERR, "Non-existing group requested");
      write_group(fp,gid);
    }
  }

  if (comm->me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
  write out one group to a Gromacs style index file, fp is non-null only on MPI rank 0
   ---------------------------------------------------------------------- */
void Group2Ndx::write_group(FILE *fp, int gid)
{
  tagint *sendlist, *recvlist;
  bigint gcount = group->count(gid);
  int lnum, width, cols;

  if (fp) {
    utils::logmesg(lmp," writing group {}...",group->names[gid]);

    // the "all" group in LAMMPS is called "System" in Gromacs
    if (gid == 0) {
      fputs("[ System ]\n", fp);
    } else {
      fmt::print(fp,"[ {} ]\n", group->names[gid]);
    }
    width = log10((double) atom->natoms)+2;
    cols = 80 / width;
  }

  if (gcount > 0) {
    const int * const mask = atom->mask;
    const tagint * const tag = atom->tag;
    const int groupbit = group->bitmask[gid];
    const int nlocal = atom->nlocal;
    int i;

    sendlist = new tagint[nlocal];
    recvlist = new tagint[gcount];
    lnum = 0;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) sendlist[lnum++] = tag[i];

    int nrecv=0;
    bigint allrecv;
    if (comm->me == 0) {
      MPI_Status status;
      MPI_Request request;

      for (i=0; i < lnum; i++)
        recvlist[i] = sendlist[i];

      allrecv = lnum;
      for (i=1; i < comm->nprocs; ++i) {
        MPI_Irecv(recvlist+allrecv,gcount-allrecv,MPI_LMP_TAGINT,i,0, world,&request);
        MPI_Send(&nrecv,0,MPI_INT,i,0,world); // block rank "i" until we are ready to receive
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&nrecv);
        allrecv += nrecv;
      }

      // sort received list
      qsort((void *)recvlist, allrecv, sizeof(tagint), cmptagint);
    } else {
      MPI_Recv(&nrecv,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
      MPI_Rsend(sendlist,lnum,MPI_LMP_TAGINT,0,0,world);
    }
    delete [] sendlist;
  }

  if (fp) {
    int i, j;
    for (i=0, j=0; i < gcount; ++i) {
      fmt::print(fp,"{:>{}}",recvlist[i],width);
      ++j;
      if (j == cols) {
        fputs("\n",fp);
        j = 0;
      }
    }
    if (j > 0) fputs("\n",fp);
    utils::logmesg(lmp,"done\n");
  }
  if (gcount > 0) delete[] recvlist;
}
