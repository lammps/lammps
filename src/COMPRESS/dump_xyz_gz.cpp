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

#include "dump_xyz_gz.h"
#include "error.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

DumpXYZGZ::DumpXYZGZ(LAMMPS *lmp, int narg, char **arg) :
  DumpXYZ(lmp, narg, arg)
{
  gzFp = NULL;

  if (!compressed)
    error->all(FLERR,"Dump xyz/gz only writes compressed files");
}


/* ---------------------------------------------------------------------- */

DumpXYZGZ::~DumpXYZGZ()
{
  if (gzFp) gzclose(gzFp);
  gzFp = NULL;
  fp = NULL;
}


/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or gzipped
   some derived classes override this function
------------------------------------------------------------------------- */

void DumpXYZGZ::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  // if one file per timestep, replace '*' with current timestep

  char *filecurrent = filename;
  if (multiproc) filecurrent = multiname;

  if (multifile) {
    char *filestar = filecurrent;
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    if (padflag == 0)
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
              filestar,update->ntimestep,ptr+1);
    else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
    }
    *ptr = '*';
    if (maxfiles > 0) {
      if (numfiles < maxfiles) {
        nameslist[numfiles] = new char[strlen(filecurrent)+1];
        strcpy(nameslist[numfiles],filecurrent);
        ++numfiles;
      } else {
        remove(nameslist[fileidx]);
        delete[] nameslist[fileidx];
        nameslist[fileidx] = new char[strlen(filecurrent)+1];
        strcpy(nameslist[fileidx],filecurrent);
        fileidx = (fileidx + 1) % maxfiles;
      }
    }
  }

  // each proc with filewriter = 1 opens a file

  if (filewriter) {
    if (append_flag) {
      gzFp = gzopen(filecurrent,"ab9");
    } else {
      gzFp = gzopen(filecurrent,"wb9");
    }

    if (gzFp == NULL) error->one(FLERR,"Cannot open dump file");
  } else gzFp = NULL;

  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}

void DumpXYZGZ::write_header(bigint ndump)
{
  if (me == 0) {
    gzprintf(gzFp,BIGINT_FORMAT "\n",ndump);
    gzprintf(gzFp,"Atoms. Timestep: " BIGINT_FORMAT "\n",update->ntimestep);
  }
}

/* ---------------------------------------------------------------------- */

void DumpXYZGZ::write_data(int n, double *mybuf)
{
  gzwrite(gzFp,mybuf,sizeof(char)*n);
}

/* ---------------------------------------------------------------------- */

void DumpXYZGZ::write()
{
  DumpXYZ::write();
  if (filewriter) {
    if (multifile) {
      gzclose(gzFp);
      gzFp = NULL;
    } else {
      if (flush_flag)
        gzflush(gzFp,Z_SYNC_FLUSH);
    }
  }
}
