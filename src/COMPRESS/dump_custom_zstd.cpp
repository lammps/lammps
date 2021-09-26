// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#ifdef LAMMPS_ZSTD

#include "dump_custom_zstd.h"

#include "file_writer.h"
#include "domain.h"
#include "error.h"
#include "update.h"

#include <cstring>


using namespace LAMMPS_NS;

DumpCustomZstd::DumpCustomZstd(LAMMPS *lmp, int narg, char **arg) :
  DumpCustom(lmp, narg, arg)
{
  if (!compressed)
    error->all(FLERR,"Dump custom/zstd only writes compressed files");
}

/* ---------------------------------------------------------------------- */

DumpCustomZstd::~DumpCustomZstd()
{
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or gzipped
   some derived classes override this function
------------------------------------------------------------------------- */

void DumpCustomZstd::openfile()
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
        if (remove(nameslist[fileidx]) != 0) {
          error->warning(FLERR, "Could not delete {}", nameslist[fileidx]);
        }
        delete[] nameslist[fileidx];
        nameslist[fileidx] = utils::strdup(filecurrent);
        fileidx = (fileidx + 1) % maxfiles;
      }
    }
  }

  // each proc with filewriter = 1 opens a file

  if (filewriter) {
    if (append_flag) {
      error->one(FLERR, "dump/zstd currently doesn't support append");
    }

    try {
      writer.open(filecurrent);
    } catch (FileWriterException &e) {
      error->one(FLERR, e.what());
    }
  }

  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */

void DumpCustomZstd::write_header(bigint ndump)
{
  std::string header;

  if ((multiproc) || (!multiproc && me == 0)) {
    if (unit_flag && !unit_count) {
      ++unit_count;
      header = fmt::format("ITEM: UNITS\n{}\n",update->unit_style);
    }

    if (time_flag) {
      header += fmt::format("ITEM: TIME\n{0:.16g}\n", compute_time());
    }

    header += fmt::format("ITEM: TIMESTEP\n{}\n", update->ntimestep);
    header += fmt::format("ITEM: NUMBER OF ATOMS\n{}\n", ndump);
    if (domain->triclinic == 0) {
      header += fmt::format("ITEM: BOX BOUNDS {}\n", boundstr);
      header += fmt::format("{0:-1.16e} {1:-1.16e}\n", boxxlo, boxxhi);
      header += fmt::format("{0:-1.16e} {1:-1.16e}\n", boxylo, boxyhi);
      header += fmt::format("{0:-1.16e} {1:-1.16e}\n", boxzlo, boxzhi);
    } else {
      header += fmt::format("ITEM: BOX BOUNDS xy xz yz {}\n", boundstr);
      header += fmt::format("{0:-1.16e} {1:-1.16e} {2:-1.16e}\n", boxxlo, boxxhi, boxxy);
      header += fmt::format("{0:-1.16e} {1:-1.16e} {2:-1.16e}\n", boxylo, boxyhi, boxxz);
      header += fmt::format("{0:-1.16e} {1:-1.16e} {2:-1.16e}\n", boxzlo, boxzhi, boxyz);
    }
    header += fmt::format("ITEM: ATOMS {}\n", columns);

    writer.write(header.c_str(), header.length());
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomZstd::write_data(int n, double *mybuf)
{
  if (buffer_flag == 1) {
    writer.write(mybuf, n);
  } else {
    constexpr size_t VBUFFER_SIZE = 256;
    char vbuffer[VBUFFER_SIZE];
    int m = 0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < nfield; j++) {
        int written = 0;
        if (vtype[j] == Dump::INT) {
          written = snprintf(vbuffer, VBUFFER_SIZE, vformat[j], static_cast<int> (mybuf[m]));
        } else if (vtype[j] == Dump::DOUBLE) {
          written = snprintf(vbuffer, VBUFFER_SIZE, vformat[j], mybuf[m]);
        } else if (vtype[j] == Dump::STRING) {
          written = snprintf(vbuffer, VBUFFER_SIZE, vformat[j], typenames[(int) mybuf[m]]);
        } else if (vtype[j] == Dump::BIGINT) {
          written = snprintf(vbuffer, VBUFFER_SIZE, vformat[j], static_cast<bigint> (mybuf[m]));
        }

        if (written > 0) {
          writer.write(vbuffer, written);
        } else if (written < 0) {
          error->one(FLERR, "Error while writing dump custom/gz output");
        }
        m++;
      }
      writer.write("\n", 1);
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomZstd::write()
{
  DumpCustom::write();
  if (filewriter) {
    if (multifile) {
      writer.close();
    } else {
      if (flush_flag && writer.isopen()) {
        writer.flush();
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpCustomZstd::modify_param(int narg, char **arg)
{
  int consumed = DumpCustom::modify_param(narg, arg);
  if (consumed == 0) {
    try {
      if (strcmp(arg[0],"checksum") == 0) {
        if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
        if (strcmp(arg[1],"yes") == 0) writer.setChecksum(true);
        else if (strcmp(arg[1],"no") == 0) writer.setChecksum(false);
        else error->all(FLERR,"Illegal dump_modify command");
        return 2;
      } else if (strcmp(arg[0],"compression_level") == 0) {
        if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
        int compression_level = utils::inumeric(FLERR, arg[1], false, lmp);
        writer.setCompressionLevel(compression_level);
        return 2;
      }
    } catch (FileWriterException &e) {
      error->one(FLERR,"Illegal dump_modify command: {}", e.what());
    }
  }
  return consumed;
}

#endif
