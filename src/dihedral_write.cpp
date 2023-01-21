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
   Contributing authors:  Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "dihedral_write.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "math_const.h"
#include "update.h"

#include <cmath>
using namespace LAMMPS_NS;
using MathConst::DEG2RAD;
using MathConst::RAD2DEG;

static constexpr double epsilon = 6.5e-6;
#define MAXLINE 1024
/* ---------------------------------------------------------------------- */

void DihedralWrite::command(int narg, char **arg)
{
  // sanity checks

  if (domain->box_exist == 0)
    error->all(FLERR, "Dihedral_write command before simulation box is defined");
  if (atom->avec->dihedrals_allow == 0)
    error->all(FLERR, "Dihedral_write command when no dihedrals allowed");
  auto dihedral = force->dihedral;
  if (dihedral == nullptr)
    error->all(FLERR, "Dihedral_write command before an dihedral_style is defined");
  if (dihedral && (force->dihedral->writedata == 0))
    error->all(FLERR, "Dihedral style must support writing coeffs to data file for dihedral_write");

  if (dihedral &&
      (utils::strmatch(force->dihedral_style, "^charmm") ||
       utils::strmatch(force->dihedral_style, "^class2")))
    error->all(FLERR, "Dihedral_write command is not compatible with dihedral style {}",
               force->dihedral_style);

  // parse arguments

  if (narg != 4) error->all(FLERR, "Illegal dihedral_write command");

  int dtype = utils::inumeric(FLERR, arg[0], false, lmp);
  if ((dtype <= 0) || (dtype > atom->ndihedraltypes))
    error->all(FLERR, "Invalid dihedral type {} in dihedral_write command", dtype);

  int n = utils::inumeric(FLERR, arg[1], false, lmp);
  std::string table_file = arg[2];
  std::string keyword = arg[3];
  if (n < 2) error->all(FLERR, "Must have at least 2 table values");

  // make sure system is initialized before calling any functions

  lmp->init();

  // write out all dihedral_coeff settings to file. use function from write_data.
  // open table file in append mode if it already exists
  // add line with DATE: and UNITS: tag when creating new file
  // otherwise make certain that units are consistent
  // print header in format used by dihedral_style table

  FILE *fp = nullptr;
  std::string coeffs_file = table_file + ".tmp.coeffs";
  if (comm->me == 0) {

    fp = fopen(coeffs_file.c_str(), "w");
    force->dihedral->write_data(fp);
    fclose(fp);

    // units sanity check:
    // - if this is the first time we write to this potential file,
    //   write out a line with "DATE:" and "UNITS:" tags
    // - if the file already exists, print a message about appending
    //   while printing the date and check that units are consistent.
    if (platform::file_is_readable(table_file)) {
      std::string units = utils::get_potential_units(table_file, "table");
      if (!units.empty() && (units != update->unit_style)) {
        error->one(FLERR, "Trying to append to a table file with UNITS: {} while units are {}",
                   units, update->unit_style);
      }
      std::string date = utils::get_potential_date(table_file, "table");
      utils::logmesg(lmp, "Appending to table file {} with DATE: {}\n", table_file, date);
      fp = fopen(table_file.c_str(), "a");
    } else {
      utils::logmesg(lmp, "Creating table file {} with DATE: {}\n", table_file,
                     utils::current_date());
      fp = fopen(table_file.c_str(), "w");
      if (fp)
        fmt::print(fp, "# DATE: {} UNITS: {} Created by dihedral_write\n", utils::current_date(),
                   update->unit_style);
    }
    if (fp == nullptr)
      error->one(FLERR, "Cannot open dihedral_write file {}: {}", table_file, utils::getsyserror());
  }

  // split communicator so that we can run a new LAMMPS class instance only on comm->me == 0

  MPI_Comm singlecomm;
  int color = (comm->me == 0) ? 1 : MPI_UNDEFINED;
  int key = comm->me;
  MPI_Comm_split(world, color, key, &singlecomm);

  if (comm->me == 0) {
    // set up new LAMMPS instance with dummy system to evaluate dihedral potential
    //    const char *args[] = {"DihedralWrite", "-nocite", "-echo",   "none",
    //                          "-log",          "none",    "-screen", "none"};
    const char *args[] = {"DihedralWrite", "-nocite", "-echo", "screen", "-log", "none"};
    char **argv = (char **) args;
    int argc = sizeof(args) / sizeof(char *);
    LAMMPS *writer = new LAMMPS(argc, argv, singlecomm);

    // create dummy system replicating dihedral style settings
    writer->input->one(fmt::format("units {}", update->unit_style));
    writer->input->one("atom_style molecular");
    writer->input->one("atom_modify map array");
    writer->input->one("boundary f f f");
    writer->input->one("region box block -2 2 -2 2 -2 2");
    writer->input->one(fmt::format("create_box {} box dihedral/types {} "
                                   "extra/dihedral/per/atom 1 "
                                   "extra/special/per/atom 4",
                                   atom->ntypes, atom->ndihedraltypes));
    writer->input->one("create_atoms 1 single  1.0  0.0 -1.0");
    writer->input->one("create_atoms 1 single  0.0  0.0 -1.0");
    writer->input->one("create_atoms 1 single  0.0  0.0  1.0");
    writer->input->one("create_atoms 1 single  1.0  0.0  1.0");
    writer->input->one(fmt::format("create_bonds single/dihedral {} 1 2 3 4", dtype));

    writer->input->one("pair_style zero 10.0");
    writer->input->one("pair_coeff * *");
    writer->input->one("mass * 1.0");
    writer->input->one(fmt::format("dihedral_style {}", force->dihedral_style));
    FILE *coeffs;
    char line[MAXLINE];
    coeffs = fopen(coeffs_file.c_str(), "r");
    for (int i = 0; i < atom->ndihedraltypes; ++i) {
      fgets(line, MAXLINE, coeffs);
      writer->input->one(fmt::format("dihedral_coeff {}", line));
    }
    fclose(coeffs);
    platform::unlink(coeffs_file);

    // must complete a full setup() to initialize system including neighbor and dihedral lists.

    writer->input->one("run 0 post no");

    // move third atom to reproduce dihedrals

    double theta, phi, phi1, phi2, f;
    dihedral = writer->force->dihedral;
    auto atom4 = writer->atom->x[writer->atom->map(4)];

    // evaluate energy and force at each of N distances

    fmt::print(fp, "# Dihedral potential {} for dihedral type {}: i,theta,energy,force\n",
               force->dihedral_style, dtype);
    fmt::print(fp, "\n{}\nN {} DEGREES\n\n", keyword, n);

#define GET_ENERGY(myphi, mytheta)     \
  theta = mytheta;                     \
  atom4[0] = cos(theta * DEG2RAD);     \
  atom4[1] = sin(theta * DEG2RAD);     \
  dihedral->energy = 0.0;              \
  dihedral->compute(ENERGY_GLOBAL, 0); \
  myphi = dihedral->energy

    const double dtheta = 360.0 / static_cast<double>(n);
    for (int i = 0; i < n; i++) {
      GET_ENERGY(phi1, dtheta * static_cast<double>(i) - epsilon);
      GET_ENERGY(phi2, dtheta * static_cast<double>(i) + epsilon);
      GET_ENERGY(phi, dtheta * static_cast<double>(i));

      if (!std::isfinite(phi)) phi = 0.0;

      // get force from numerical differentiation
      f = -0.5 * (phi2 - phi1) / epsilon;
      if (!std::isfinite(f)) f = 0.0;
      fprintf(fp, "%8d %- 22.15g %- 22.15g %- 22.15g\n", i + 1, theta, phi, f);
    }

    // clean up
    delete writer;
    fclose(fp);
  }
  MPI_Comm_free(&singlecomm);
}
