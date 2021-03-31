/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "dump_cfg_gz.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "update.h"


#include <cstring>


using namespace LAMMPS_NS;
#define UNWRAPEXPAND 10.0

DumpCFGGZ::DumpCFGGZ(LAMMPS *lmp, int narg, char **arg) :
  DumpCFG(lmp, narg, arg)
{
  gzFp = nullptr;

  compression_level = Z_BEST_COMPRESSION;

  if (!compressed)
    error->all(FLERR,"Dump cfg/gz only writes compressed files");
}


/* ---------------------------------------------------------------------- */

DumpCFGGZ::~DumpCFGGZ()
{
  if (gzFp) gzclose(gzFp);
  gzFp = nullptr;
  fp = nullptr;
}


/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or gzipped
   some derived classes override this function
------------------------------------------------------------------------- */

void DumpCFGGZ::openfile()
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
        nameslist[numfiles] = utils::strdup(filecurrent);
        ++numfiles;
      } else {
        remove(nameslist[fileidx]);
        delete[] nameslist[fileidx];
        nameslist[fileidx] = utils::strdup(filecurrent);
        fileidx = (fileidx + 1) % maxfiles;
      }
    }
  }

  // each proc with filewriter = 1 opens a file

  if (filewriter) {
    std::string mode;
    if (append_flag) {
      mode = fmt::format("ab{}", compression_level);
    } else {
      mode = fmt::format("wb{}", compression_level);
    }

    gzFp = gzopen(filecurrent, mode.c_str());

    if (gzFp == nullptr) error->one(FLERR,"Cannot open dump file");
  } else gzFp = nullptr;

  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */

void DumpCFGGZ::write_header(bigint n)
{
  // set scale factor used by AtomEye for CFG viz
  // default = 1.0
  // for peridynamics, set to pre-computed PD scale factor
  //   so PD particles mimic C atoms
  // for unwrapped coords, set to UNWRAPEXPAND (10.0)
  //   so molecules are not split across periodic box boundaries

  double scale = 1.0;
  if (atom->peri_flag) scale = atom->pdscale;
  else if (unwrapflag == 1) scale = UNWRAPEXPAND;

  char str[64];
  sprintf(str,"Number of particles = %s\n",BIGINT_FORMAT);
  gzprintf(gzFp,str,n);
  gzprintf(gzFp,"A = %g Angstrom (basic length-scale)\n",scale);
  gzprintf(gzFp,"H0(1,1) = %g A\n",domain->xprd);
  gzprintf(gzFp,"H0(1,2) = 0 A \n");
  gzprintf(gzFp,"H0(1,3) = 0 A \n");
  gzprintf(gzFp,"H0(2,1) = %g A \n",domain->xy);
  gzprintf(gzFp,"H0(2,2) = %g A\n",domain->yprd);
  gzprintf(gzFp,"H0(2,3) = 0 A \n");
  gzprintf(gzFp,"H0(3,1) = %g A \n",domain->xz);
  gzprintf(gzFp,"H0(3,2) = %g A \n",domain->yz);
  gzprintf(gzFp,"H0(3,3) = %g A\n",domain->zprd);
  gzprintf(gzFp,".NO_VELOCITY.\n");
  gzprintf(gzFp,"entry_count = %d\n",nfield-2);
  for (int i = 0; i < nfield-5; i++)
    gzprintf(gzFp,"auxiliary[%d] = %s\n",i,auxname[i]);
}

/* ---------------------------------------------------------------------- */

void DumpCFGGZ::write_data(int n, double *mybuf)
{
  gzwrite(gzFp,mybuf,sizeof(char)*n);
}

/* ---------------------------------------------------------------------- */

void DumpCFGGZ::write()
{
  DumpCFG::write();
  if (filewriter) {
    if (multifile) {
      gzclose(gzFp);
      gzFp = nullptr;
    } else {
      if (flush_flag)
        gzflush(gzFp,Z_SYNC_FLUSH);
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpCFGGZ::modify_param(int narg, char **arg)
{
  int consumed = DumpCFG::modify_param(narg, arg);
  if (consumed == 0) {
    if (strcmp(arg[0],"compression_level") == 0) {
      if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
      int min_level = Z_DEFAULT_COMPRESSION;
      int max_level = Z_BEST_COMPRESSION;
      compression_level = utils::inumeric(FLERR, arg[1], false, lmp);
      if (compression_level < min_level || compression_level > max_level)
        error->all(FLERR, fmt::format("Illegal dump_modify command: compression level must in the range of [{}, {}]", min_level, max_level));
      return 2;
    }
  }
  return consumed;
}

