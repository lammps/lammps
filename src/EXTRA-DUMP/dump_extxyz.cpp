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

#include "dump_extxyz.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "label_map.h"
#include "memory.h"
#include "modify.h"
#include "output.h"
#include "thermo.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

static constexpr int ONELINE = 512;
static constexpr int DELTA = 1048576;

/* ---------------------------------------------------------------------- */

DumpExtXYZ::DumpExtXYZ(LAMMPS *lmp, int narg, char **arg) :
    DumpXYZ(lmp, narg, arg), properties_string(nullptr)
{
  // style specific customizable settings
  with_vel = 1;
  with_forces = 1;
  with_mass = 0;
  with_pe = 1;
  with_temp = 1;
  with_press = 0;

  update_properties();

  // We want simulation time by default
  time_flag = 1;

  // dump may invoke computes
  clearstep = 1;

  // use type labels by default if present
  if (atom->labelmapflag) {
    typenames = new char *[ntypes + 1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = utils::strdup(atom->lmap->typelabel[itype - 1]);
    }
  }
}

/* ---------------------------------------------------------------------- */

DumpExtXYZ::~DumpExtXYZ()
{
  delete[] properties_string;
}

/* ---------------------------------------------------------------------- */

void DumpExtXYZ::update_properties()
{
  // How many per-atom elements we buffer
  size_one = 5 + (with_vel ? 3 : 0) + (with_forces ? 3 : 0) + (with_mass ? 1 : 0);

  // The properties string
  delete[] properties_string;
  properties_string = utils::strdup(
      fmt::format("species:S:1:pos:R:3{}{}{}", (with_vel ? ":vel:R:3" : ""),
                  (with_forces ? ":forces:R:3" : ""), (with_mass ? ":mass:R:1" : "")));

  // The output printf-style format
  delete[] format;
  if (format_line_user)
    format = utils::strdup(fmt::format("{}\n", format_line_user));
  else {
    format = utils::strdup(fmt::format("%s %g %g %g{}{}{}\n", (with_vel ? " %g %g %g" : ""),
                                       (with_forces ? " %g %g %g" : ""), (with_mass ? " %g" : "")));
  }
}

/* ---------------------------------------------------------------------- */

void DumpExtXYZ::init_style()
{
  if (!typenames)
    error->all(FLERR, Error::NOLASTLINE,
               "Must use either type lables or dump_modify element with dump style extxyz");

  DumpXYZ::init_style();
  update_properties();
}

/* ---------------------------------------------------------------------- */

int DumpExtXYZ::modify_param(int narg, char **arg)
{
  int rv = DumpXYZ::modify_param(narg, arg);
  if (rv > 0) return rv;

  if (strcmp(arg[0], "vel") == 0) {
    if (narg < 2) error->all(FLERR, "Illegal dump_modify command");
    with_vel = utils::logical(FLERR, arg[1], false, lmp);
    update_properties();
    return 2;
  }

  if (strcmp(arg[0], "forces") == 0) {
    if (narg < 2) error->all(FLERR, "Illegal dump_modify command");
    with_forces = utils::logical(FLERR, arg[1], false, lmp);
    update_properties();
    return 2;
  }

  if (strcmp(arg[0], "mass") == 0) {
    if (narg < 2) error->all(FLERR, "Illegal dump_modify command");
    with_mass = utils::logical(FLERR, arg[1], false, lmp);
    update_properties();
    return 2;
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpExtXYZ::write_header(bigint n)
{
  if (me == 0) {
    if (!fp)
      error->one(FLERR, Error::NOLASTLINE, "Must not use 'run pre no' after creating a new dump");

    std::string header = fmt::format("{}\nTimestep={}", n, update->ntimestep);
    if (time_flag) header += fmt::format(" Time={:.6f}", compute_time());
    header += fmt::format(" pbc=\"{} {} {}\"", domain->xperiodic ? "T" : "F",
                          domain->yperiodic ? "T" : "F", domain->zperiodic ? "T" : "F");
    header +=
        fmt::format(" Lattice=\"{:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g}\"", domain->xprd, 0.,
                    0., domain->xy, domain->yprd, 0., domain->xz, domain->yz, domain->zprd);

    if (output && output->thermo) {
      auto *pe = output->thermo->pe;
      if (pe) header += fmt::format(" Potential_energy={}", pe->compute_scalar());

      auto *temp = output->thermo->temperature;
      if (temp) header += fmt::format(" Temperature={}", temp->compute_scalar());

      auto *press = output->thermo->pressure;
      if (press) {
        press->compute_vector();
        header +=
            fmt::format(" Stress=\"{} {} {} {} {} {} {} {} {}\"", press->vector[0],
                        press->vector[3], press->vector[4], press->vector[3], press->vector[1],
                        press->vector[5], press->vector[4], press->vector[5], press->vector[2]);
      }
    }

    header += fmt::format(" Properties={}", properties_string);
    utils::print(fp, header + "\n");
  }
}

/* ---------------------------------------------------------------------- */

void DumpExtXYZ::pack(tagint *ids)
{
  int m, n;

  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      if (with_vel) {
        buf[m++] = v[i][0];
        buf[m++] = v[i][1];
        buf[m++] = v[i][2];
      }
      if (with_forces) {
        buf[m++] = f[i][0];
        buf[m++] = f[i][1];
        buf[m++] = f[i][2];
      }
      if (with_mass) {
        if (rmass) {
          buf[m++] = rmass[i];
        } else {
          buf[m++] = mass[type[i]];
        }
      }

      if (ids) ids[n++] = tag[i];
    }
}

/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpExtXYZ::convert_string(int n, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf, maxsbuf, "dump:sbuf");
    }

    if (size_one == 5) {
      offset += snprintf(&sbuf[offset], maxsbuf - offset, format,
                         typenames[static_cast<int>(mybuf[m + 1])], mybuf[m + 2], mybuf[m + 3],
                         mybuf[m + 4]);
    } else if (size_one == 6) {
      offset += snprintf(&sbuf[offset], maxsbuf - offset, format,
                         typenames[static_cast<int>(mybuf[m + 1])], mybuf[m + 2], mybuf[m + 3],
                         mybuf[m + 4], mybuf[m + 5]);
    } else if (size_one == 8) {
      offset += snprintf(&sbuf[offset], maxsbuf - offset, format,
                         typenames[static_cast<int>(mybuf[m + 1])], mybuf[m + 2], mybuf[m + 3],
                         mybuf[m + 4], mybuf[m + 5], mybuf[m + 6], mybuf[m + 7]);
    } else if (size_one == 9) {
      offset += snprintf(&sbuf[offset], maxsbuf - offset, format,
                         typenames[static_cast<int>(mybuf[m + 1])], mybuf[m + 2], mybuf[m + 3],
                         mybuf[m + 4], mybuf[m + 5], mybuf[m + 6], mybuf[m + 7], mybuf[m + 8]);
    } else if (size_one == 11) {
      offset += snprintf(&sbuf[offset], maxsbuf - offset, format,
                         typenames[static_cast<int>(mybuf[m + 1])], mybuf[m + 2], mybuf[m + 3],
                         mybuf[m + 4], mybuf[m + 5], mybuf[m + 6], mybuf[m + 7], mybuf[m + 8],
                         mybuf[m + 9], mybuf[m + 10]);
    } else if (size_one == 12) {
      offset += snprintf(&sbuf[offset], maxsbuf - offset, format,
                         typenames[static_cast<int>(mybuf[m + 1])], mybuf[m + 2], mybuf[m + 3],
                         mybuf[m + 4], mybuf[m + 5], mybuf[m + 6], mybuf[m + 7], mybuf[m + 8],
                         mybuf[m + 9], mybuf[m + 10], mybuf[m + 11]);
    } else {
      error->all(FLERR, "Invalid value of size_one for dump extxyz format.");
    }
    m += size_one;
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

void DumpExtXYZ::write_lines(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (size_one == 5) {
      fprintf(fp, format, typenames[static_cast<int>(mybuf[m + 1])], mybuf[m + 2], mybuf[m + 3],
              mybuf[m + 4]);
    } else if (size_one == 6) {
      fprintf(fp, format, typenames[static_cast<int>(mybuf[m + 1])], mybuf[m + 2], mybuf[m + 3],
              mybuf[m + 4], mybuf[m + 5]);
    } else if (size_one == 8) {
      fprintf(fp, format, typenames[static_cast<int>(mybuf[m + 1])], mybuf[m + 2], mybuf[m + 3],
              mybuf[m + 4], mybuf[m + 5], mybuf[m + 6], mybuf[m + 7]);
    } else if (size_one == 9) {
      fprintf(fp, format, typenames[static_cast<int>(mybuf[m + 1])], mybuf[m + 2], mybuf[m + 3],
              mybuf[m + 4], mybuf[m + 5], mybuf[m + 6], mybuf[m + 7], mybuf[m + 8]);
    } else if (size_one == 11) {
      fprintf(fp, format, typenames[static_cast<int>(mybuf[m + 1])], mybuf[m + 2], mybuf[m + 3],
              mybuf[m + 4], mybuf[m + 5], mybuf[m + 6], mybuf[m + 7], mybuf[m + 8], mybuf[m + 9],
              mybuf[m + 10]);
    } else if (size_one == 12) {
      fprintf(fp, format, typenames[static_cast<int>(mybuf[m + 1])], mybuf[m + 2], mybuf[m + 3],
              mybuf[m + 4], mybuf[m + 5], mybuf[m + 6], mybuf[m + 7], mybuf[m + 8], mybuf[m + 9],
              mybuf[m + 10], mybuf[m + 11]);
    } else {
      error->all(FLERR, "Invalid value of size_one for dump extxyz format.");
    }
    m += size_one;
  }
}
