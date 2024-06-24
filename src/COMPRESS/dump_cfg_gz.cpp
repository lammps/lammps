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

#include "dump_cfg_gz.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "file_writer.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
static constexpr double UNWRAPEXPAND = 10.0;

DumpCFGGZ::DumpCFGGZ(LAMMPS *lmp, int narg, char **arg) : DumpCFG(lmp, narg, arg)
{
  if (!compressed) error->all(FLERR, "Dump cfg/gz only writes compressed files");
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or compressed
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
    filecurrent = utils::strdup(utils::star_subst(filecurrent, update->ntimestep, padflag));
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
    try {
      writer.open(filecurrent, append_flag);
    } catch (FileWriterException &e) {
      error->one(FLERR, e.what());
    }
  }

  // delete string with timestep replaced

  if (multifile) delete[] filecurrent;
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
  if (atom->peri_flag)
    scale = atom->pdscale;
  else if (unwrapflag == 1)
    scale = UNWRAPEXPAND;

  std::string header = fmt::format("Number of particles = {}\n", n);
  header += fmt::format("A = {:g} Angstrom (basic length-scale)\n", scale);
  header += fmt::format("H0(1,1) = {:g} A\n", domain->xprd);
  header += fmt::format("H0(1,2) = 0 A\n");
  header += fmt::format("H0(1,3) = 0 A\n");
  header += fmt::format("H0(2,1) = {:g} A\n", domain->xy);
  header += fmt::format("H0(2,2) = {:g} A\n", domain->yprd);
  header += fmt::format("H0(2,3) = 0 A\n");
  header += fmt::format("H0(3,1) = {:g} A\n", domain->xz);
  header += fmt::format("H0(3,2) = {:g} A\n", domain->yz);
  header += fmt::format("H0(3,3) = {:g} A\n", domain->zprd);
  header += fmt::format(".NO_VELOCITY.\n");
  header += fmt::format("entry_count = {}\n", nfield - 2);
  for (int i = 0; i < nfield - 5; i++) header += fmt::format("auxiliary[{}] = {}\n", i, auxname[i]);

  writer.write(header.c_str(), header.length());
}

/* ---------------------------------------------------------------------- */

void DumpCFGGZ::write_data(int n, double *mybuf)
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

void DumpCFGGZ::write()
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

int DumpCFGGZ::modify_param(int narg, char **arg)
{
  int consumed = DumpCFG::modify_param(narg, arg);
  if (consumed == 0) {
    try {
      if (strcmp(arg[0], "compression_level") == 0) {
        if (narg < 2) error->all(FLERR, "Illegal dump_modify command");
        int compression_level = utils::inumeric(FLERR, arg[1], false, lmp);
        writer.setCompressionLevel(compression_level);
        return 2;
      }
    } catch (FileWriterException &e) {
      error->one(FLERR, "Illegal dump_modify command: {}", e.what());
    }
  }
  return consumed;
}
