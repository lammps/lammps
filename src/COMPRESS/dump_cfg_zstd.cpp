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

#include "dump_cfg_zstd.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "file_writer.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
#define UNWRAPEXPAND 10.0

DumpCFGZstd::DumpCFGZstd(LAMMPS *lmp, int narg, char **arg) : DumpCFG(lmp, narg, arg)
{
  if (!compressed) error->all(FLERR, "Dump cfg/zstd only writes compressed files");
}

/* ---------------------------------------------------------------------- */

DumpCFGZstd::~DumpCFGZstd() {}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or compressed
   some derived classes override this function
------------------------------------------------------------------------- */

void DumpCFGZstd::openfile()
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
    char *ptr = strchr(filestar, '*');
    *ptr = '\0';
    if (padflag == 0)
      sprintf(filecurrent, "%s" BIGINT_FORMAT "%s", filestar, update->ntimestep, ptr + 1);
    else {
      char bif[8], pad[16];
      strcpy(bif, BIGINT_FORMAT);
      sprintf(pad, "%%s%%0%d%s%%s", padflag, &bif[1]);
      sprintf(filecurrent, pad, filestar, update->ntimestep, ptr + 1);
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

void DumpCFGZstd::write_header(bigint n)
{
  // set scale factor used by AtomEye for CFG viz
  // default = 1.0
  // for peridynamics, set to pre-computed PD scale factor
  //   so PD particles mimic C atoms
  // for unwrapped coords, set to UNWRAPEXPAND (10.0)
  //   so molecules are not split across periodic box boundaries

  double scale = 1.0;
  if (atom->peri_flag)
    scale = atom->pdscale;
  else if (unwrapflag == 1)
    scale = UNWRAPEXPAND;

  std::string header = fmt::format("Number of particles = {}\n", n);
  header += fmt::format("A = {0:g} Angstrom (basic length-scale)\n", scale);
  header += fmt::format("H0(1,1) = {0:g} A\n", domain->xprd);
  header += fmt::format("H0(1,2) = 0 A \n");
  header += fmt::format("H0(1,3) = 0 A \n");
  header += fmt::format("H0(2,1) = {0:g} A \n", domain->xy);
  header += fmt::format("H0(2,2) = {0:g} A\n", domain->yprd);
  header += fmt::format("H0(2,3) = 0 A \n");
  header += fmt::format("H0(3,1) = {0:g} A \n", domain->xz);
  header += fmt::format("H0(3,2) = {0:g} A \n", domain->yz);
  header += fmt::format("H0(3,3) = {0:g} A\n", domain->zprd);
  header += fmt::format(".NO_VELOCITY.\n");
  header += fmt::format("entry_count = {}\n", nfield - 2);
  for (int i = 0; i < nfield - 5; i++) header += fmt::format("auxiliary[{}] = {}\n", i, auxname[i]);

  writer.write(header.c_str(), header.length());
}

/* ---------------------------------------------------------------------- */

void DumpCFGZstd::write_data(int n, double *mybuf)
{
  if (buffer_flag) {
    writer.write(mybuf, n);
  } else {
    constexpr size_t VBUFFER_SIZE = 256;
    char vbuffer[VBUFFER_SIZE];
    if (unwrapflag == 0) {
      int m = 0;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < size_one; j++) {
          int written = 0;
          if (j == 0) {
            written = snprintf(vbuffer, VBUFFER_SIZE, "%f \n", mybuf[m]);
          } else if (j == 1) {
            written = snprintf(vbuffer, VBUFFER_SIZE, "%s \n", typenames[(int) mybuf[m]]);
          } else if (j >= 2) {
            if (vtype[j] == Dump::INT)
              written = snprintf(vbuffer, VBUFFER_SIZE, vformat[j], static_cast<int>(mybuf[m]));
            else if (vtype[j] == Dump::DOUBLE)
              written = snprintf(vbuffer, VBUFFER_SIZE, vformat[j], mybuf[m]);
            else if (vtype[j] == Dump::STRING)
              written = snprintf(vbuffer, VBUFFER_SIZE, vformat[j], typenames[(int) mybuf[m]]);
            else if (vtype[j] == Dump::BIGINT)
              written = snprintf(vbuffer, VBUFFER_SIZE, vformat[j], static_cast<bigint>(mybuf[m]));
          }
          if (written > 0) {
            writer.write(vbuffer, written);
          } else if (written < 0) {
            error->one(FLERR, "Error while writing dump cfg/gz output");
          }
          m++;
        }
        writer.write("\n", 1);
      }
    } else if (unwrapflag == 1) {
      int m = 0;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < size_one; j++) {
          int written = 0;
          if (j == 0) {
            written = snprintf(vbuffer, VBUFFER_SIZE, "%f \n", mybuf[m]);
          } else if (j == 1) {
            written = snprintf(vbuffer, VBUFFER_SIZE, "%s \n", typenames[(int) mybuf[m]]);
          } else if (j >= 2 && j <= 4) {
            double unwrap_coord = (mybuf[m] - 0.5) / UNWRAPEXPAND + 0.5;
            written = snprintf(vbuffer, VBUFFER_SIZE, vformat[j], unwrap_coord);
          } else if (j >= 5) {
            if (vtype[j] == Dump::INT)
              written = snprintf(vbuffer, VBUFFER_SIZE, vformat[j], static_cast<int>(mybuf[m]));
            else if (vtype[j] == Dump::DOUBLE)
              written = snprintf(vbuffer, VBUFFER_SIZE, vformat[j], mybuf[m]);
            else if (vtype[j] == Dump::STRING)
              written = snprintf(vbuffer, VBUFFER_SIZE, vformat[j], typenames[(int) mybuf[m]]);
            else if (vtype[j] == Dump::BIGINT)
              written = snprintf(vbuffer, VBUFFER_SIZE, vformat[j], static_cast<bigint>(mybuf[m]));
          }
          if (written > 0) {
            writer.write(vbuffer, written);
          } else if (written < 0) {
            error->one(FLERR, "Error while writing dump cfg/gz output");
          }
          m++;
        }
        writer.write("\n", 1);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpCFGZstd::write()
{
  DumpCFG::write();
  if (filewriter) {
    if (multifile) {
      writer.close();
    } else {
      if (flush_flag && writer.isopen()) { writer.flush(); }
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpCFGZstd::modify_param(int narg, char **arg)
{
  int consumed = DumpCFG::modify_param(narg, arg);
  if (consumed == 0) {
    try {
      if (strcmp(arg[0], "checksum") == 0) {
        if (narg < 2) error->all(FLERR, "Illegal dump_modify command");
        if (strcmp(arg[1], "yes") == 0)
          writer.setChecksum(true);
        else if (strcmp(arg[1], "no") == 0)
          writer.setChecksum(false);
        else
          error->all(FLERR, "Illegal dump_modify command");
        return 2;
      } else if (strcmp(arg[0], "compression_level") == 0) {
        if (narg < 2) error->all(FLERR, "Illegal dump_modify command");
        int compression_level = utils::inumeric(FLERR, arg[1], false, lmp);
        writer.setCompressionLevel(compression_level);
        return 2;
      }
    } catch (FileWriterException &e) {
      error->one(FLERR, e.what());
    }
  }
  return consumed;
}
#endif
