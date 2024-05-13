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
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#ifdef LAMMPS_ZSTD

#include "dump_xyz_zstd.h"

#include "error.h"
#include "file_writer.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

DumpXYZZstd::DumpXYZZstd(LAMMPS *lmp, int narg, char **arg) : DumpXYZ(lmp, narg, arg)
{
  if (!compressed) error->all(FLERR, "Dump xyz/zstd only writes compressed files");
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or compressed
   some derived classes override this function
------------------------------------------------------------------------- */

void DumpXYZZstd::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  // if one file per timestep, replace '*' with current timestep

  char *filecurrent = filename;
  if (multiproc) filecurrent = multiname;

  if (multifile) {
    filecurrent = utils::strdup(utils::star_subst(filecurrent, update->ntimestep, padflag));
    if (maxfiles > 0) {
      if (numfiles < maxfiles) {
        nameslist[numfiles] = utils::strdup(filecurrent);
        ++numfiles;
      } else {
        if (remove(nameslist[fileidx]) != 0) {
          error->warning(FLERR, fmt::format("Could not delete {}", nameslist[fileidx]));
        }
        delete[] nameslist[fileidx];
        nameslist[fileidx] = utils::strdup(filecurrent);
        fileidx = (fileidx + 1) % maxfiles;
      }
    }
  }

  // each proc with filewriter = 1 opens a file

  if (filewriter) {
    if (append_flag) { error->one(FLERR, "dump cfg/zstd currently doesn't support append"); }

    try {
      writer.open(filecurrent);
    } catch (FileWriterException &e) {
      error->one(FLERR, e.what());
    }
  }

  // delete string with timestep replaced

  if (multifile) delete[] filecurrent;
}

/* ---------------------------------------------------------------------- */

void DumpXYZZstd::write_header(bigint ndump)
{
  if (me == 0) {
    auto header = fmt::format("{}\n Atoms. Timestep: {}", ndump, update->ntimestep);
    if (time_flag) header += fmt::format(" Time: {:.6f}", compute_time());
    header += "\n";
    writer.write(header.c_str(), header.length());
  }
}

/* ---------------------------------------------------------------------- */

void DumpXYZZstd::write_data(int n, double *mybuf)
{
  if (buffer_flag) {
    writer.write(mybuf, n);
  } else {
    constexpr size_t VBUFFER_SIZE = 256;
    char vbuffer[VBUFFER_SIZE];
    int m = 0;
    for (int i = 0; i < n; i++) {
      int written =
          snprintf(vbuffer, VBUFFER_SIZE, format, typenames[static_cast<int>(mybuf[m + 1])],
                   mybuf[m + 2], mybuf[m + 3], mybuf[m + 4]);
      if (written > 0) {
        writer.write(vbuffer, written);
      } else if (written < 0) {
        error->one(FLERR, "Error while writing dump xyz/gz output");
      }
      m += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpXYZZstd::write()
{
  DumpXYZ::write();
  if (filewriter) {
    if (multifile) {
      writer.close();
    } else {
      if (flush_flag && writer.isopen()) { writer.flush(); }
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpXYZZstd::modify_param(int narg, char **arg)
{
  int consumed = DumpXYZ::modify_param(narg, arg);
  if (consumed == 0) {
    try {
      if (strcmp(arg[0], "checksum") == 0) {
        if (narg < 2) error->all(FLERR, "Illegal dump_modify command");
        writer.setChecksum(utils::logical(FLERR, arg[1], false, lmp) == 1);
        return 2;
      } else if (strcmp(arg[0], "compression_level") == 0) {
        if (narg < 2) error->all(FLERR, "Illegal dump_modify command");
        writer.setCompressionLevel(utils::inumeric(FLERR, arg[1], false, lmp));
        return 2;
      }
    } catch (FileWriterException &e) {
      error->one(FLERR, "Illegal dump_modify command: {}", e.what());
    }
  }
  return consumed;
}

#endif
