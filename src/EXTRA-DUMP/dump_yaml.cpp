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

#include "dump_yaml.h"

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "output.h"
#include "thermo.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
DumpYAML::DumpYAML(class LAMMPS *_lmp, int narg, char **args) :
    DumpCustom(_lmp, narg, args), thermo(false)
{
  buffer_allow = 0;
  buffer_flag = 0;
}

/* ---------------------------------------------------------------------- */

void DumpYAML::init_style()
{
  if (binary) error->all(FLERR, "Dump style yaml does not support binary output");
  if (multiproc) error->all(FLERR, "Dump style yaml does not support multi-processor output");

  DumpCustom::init_style();
}

/* ---------------------------------------------------------------------- */

void DumpYAML::write()
{
  // temporarily enable so write_header() is called
  // by all MPI ranks to compute thermo data
  if (thermo) filewriter = 1;

  Dump::write();
}

/* ---------------------------------------------------------------------- */

void DumpYAML::write_header(bigint ndump)
{
  std::string thermo_data;
  if (thermo) {
    Thermo *th = output->thermo;
    // output thermo data only on timesteps where it was computed
    if (update->ntimestep == *th->get_timestep()) {
      int nfield = *th->get_nfield();
      const auto &keywords = th->get_keywords();
      const auto &fields = th->get_fields();

      thermo_data += "thermo:\n  - keywords: [ ";
      for (int i = 0; i < nfield; ++i) thermo_data += fmt::format("{}, ", keywords[i]);
      thermo_data += "]\n  - data: [ ";

      for (int i = 0; i < nfield; ++i) {
        if (fields[i].type == multitype::DOUBLE)
          thermo_data += fmt::format("{}, ", fields[i].data.d);
        else if (fields[i].type == multitype::INT)
          thermo_data += fmt::format("{}, ", fields[i].data.i);
        else if (fields[i].type == multitype::BIGINT)
          thermo_data += fmt::format("{}, ", fields[i].data.b);
        else
          thermo_data += ", ";
      }
      thermo_data += "]\n";
      MPI_Barrier(world);
    }
  }

  if (comm->me == 0) {
    const std::string boundary(boundstr);
    fmt::print(fp, "---\ncreator: LAMMPS\ntimestep: {}\n", update->ntimestep);
    if (unit_flag) fmt::print(fp, "units: {}\n", update->unit_style);
    if (time_flag) fmt::print(fp, "time: {:.16g}\n", compute_time());

    fmt::print(fp, "natoms: {}\n", ndump);
    fputs("boundary: [ ", fp);
    for (const auto &bflag : boundary) {
      if (bflag == ' ') continue;
      fmt::print(fp, "{}, ", bflag);
    }
    fputs("]\n", fp);

    if (thermo) fmt::print(fp, thermo_data);

    fmt::print(fp, "box:\n  - [ {}, {} ]\n", boxxlo, boxxhi);
    fmt::print(fp, "  - [ {}, {} ]\n", boxylo, boxyhi);
    fmt::print(fp, "  - [ {}, {} ]\n", boxzlo, boxzhi);
    if (domain->triclinic) fmt::print(fp, "  - [ {}, {}, {} ]\n", boxxy, boxxz, boxyz);

    fmt::print(fp, "keywords: [ ");
    for (const auto &item : utils::split_words(columns)) fmt::print(fp, "{}, ", item);
    fputs(" ]\ndata:\n", fp);
  } else    // reset so that the remainder of the output is not multi-proc
    filewriter = 0;
}

/* ---------------------------------------------------------------------- */

void DumpYAML::write_data(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fputs("  - [ ", fp);
    for (int j = 0; j < nfield; j++) {
      if (vtype[j] == Dump::INT)
        fprintf(fp, vformat[j], static_cast<int>(mybuf[m]));
      else if (vtype[j] == Dump::DOUBLE)
        fprintf(fp, vformat[j], mybuf[m]);
      else if (vtype[j] == Dump::STRING)
        fprintf(fp, vformat[j], typenames[(int) mybuf[m]]);
      else if (vtype[j] == Dump::BIGINT)
        fprintf(fp, vformat[j], static_cast<bigint>(mybuf[m]));
      m++;
      fputs(", ", fp);
    }
    fputs("]\n", fp);
  }
}

/* ---------------------------------------------------------------------- */

void DumpYAML::write_footer()
{
  fputs("...\n", fp);
}

/* ---------------------------------------------------------------------- */

int DumpYAML::modify_param(int narg, char **arg)
{
  int n = DumpCustom::modify_param(narg, arg);
  if (n > 0) return n;

  if (strcmp(arg[0], "thermo") == 0) {
    if (narg < 2) error->all(FLERR, "expected 'yes' or 'no' after 'thermo' keyword.");
    thermo = utils::logical(FLERR, arg[1], false, lmp) == 1;
    return 2;
  } else
    return 0;
}
