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

#include "angle_write.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
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
static constexpr int MAXLINE = 1024;
/* ---------------------------------------------------------------------- */

void AngleWrite::command(int narg, char **arg)
{
  // sanity checks

  if (domain->box_exist == 0)
    error->all(FLERR, "Angle_write command before simulation box is defined");
  if (atom->avec->angles_allow == 0)
    error->all(FLERR, "Angle_write command when no angles allowed");
  auto angle = force->angle;
  if (angle == nullptr) error->all(FLERR, "Angle_write command before an angle_style is defined");
  if (angle && (force->angle->writedata == 0))
    error->all(FLERR, "Angle style must support writing coeffs to data file for angle_write");

  if (angle && (utils::strmatch(force->angle_style, "^class2")))
    error->all(FLERR, "Angle_write command is not compatible with angle style {}",
               force->angle_style);

  // parse arguments

  if (narg != 4) error->all(FLERR, "Illegal angle_write command");

  int atype = utils::inumeric(FLERR, arg[0], false, lmp);
  if ((atype <= 0) || (atype > atom->nangletypes))
    error->all(FLERR, "Invalid angle type {} in angle_write command", atype);

  int n = utils::inumeric(FLERR, arg[1], false, lmp);
  std::string table_file = arg[2];
  std::string keyword = arg[3];
  if (n < 2) error->all(FLERR, "Must have at least 2 table values");

  // make sure system is initialized before calling any functions

  lmp->init();

  double theta0 = angle->equilibrium_angle(atype) * RAD2DEG;

  // write out all angle_coeff settings to file. use function from write_data.
  // open table file in append mode if it already exists
  // add line with DATE: and UNITS: tag when creating new file
  // otherwise make certain that units are consistent
  // print header in format used by angle_style table

  FILE *fp = nullptr;
  std::string coeffs_file = table_file + ".tmp.coeffs";
  if (comm->me == 0) {

    fp = fopen(coeffs_file.c_str(), "w");
    force->angle->write_data(fp);
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
        fmt::print(fp, "# DATE: {} UNITS: {} Created by angle_write\n", utils::current_date(),
                   update->unit_style);
    }
    if (fp == nullptr)
      error->one(FLERR, "Cannot open angle_write file {}: {}", table_file, utils::getsyserror());
  }

  // split communicator so that we can run a new LAMMPS class instance only on comm->me == 0

  MPI_Comm singlecomm;
  int color = (comm->me == 0) ? 1 : MPI_UNDEFINED;
  int key = comm->me;
  MPI_Comm_split(world, color, key, &singlecomm);

  if (comm->me == 0) {
    // set up new LAMMPS instance with dummy system to evaluate angle potential
    LAMMPS::argv args = {"AngleWrite", "-nocite", "-echo",   "none",
                         "-log",       "none",    "-screen", "none"};
    LAMMPS *writer = new LAMMPS(args, singlecomm);

    // create dummy system replicating angle style settings
    writer->input->one(fmt::format("units {}", update->unit_style));
    writer->input->one("atom_style angle");
    writer->input->one("atom_modify map array");
    writer->input->one("boundary f f f");
    writer->input->one("region box block -2 2 -2 2 -1 1");
    writer->input->one(fmt::format("create_box {} box angle/types {} "
                                   "extra/angle/per/atom 1 "
                                   "extra/special/per/atom 4",
                                   atom->ntypes, atom->nangletypes));
    writer->input->one("create_atoms 1 single  0.0  0.0  0.0");
    writer->input->one("create_atoms 1 single  1.0  0.0  0.0");
    writer->input->one("create_atoms 1 single -1.0  0.0  0.0");
    writer->input->one(fmt::format("create_bonds single/angle {} 2 1 3", atype));

    writer->input->one("pair_style zero 10.0");
    writer->input->one("pair_coeff * *");
    writer->input->one("mass * 1.0");
    writer->input->one(fmt::format("angle_style {}", force->angle_style));
    FILE *coeffs;
    char line[MAXLINE] = {'\0'};
    coeffs = fopen(coeffs_file.c_str(), "r");
    for (int i = 0; i < atom->nangletypes; ++i) {
      fgets(line, MAXLINE, coeffs);
      writer->input->one(fmt::format("angle_coeff {}", line));
    }
    fclose(coeffs);
    platform::unlink(coeffs_file);

    // initialize system

    writer->init();

    // move third atom to reproduce angles

    double theta, phi, phi1, phi2, f;
    angle = writer->force->angle;
    int i1, i2, i3;
    i1 = writer->atom->map(1);
    i2 = writer->atom->map(2);
    i3 = writer->atom->map(3);
    auto atom3 = writer->atom->x[i3];

    // evaluate energy and force at each of N distances

    fmt::print(fp, "# Angle potential {} for angle type {}: i,theta,energy,force\n",
               force->angle_style, atype);
    fmt::print(fp, "\n{}\nN {} EQ {:.15g}\n\n", keyword, n, theta0);

#define GET_ENERGY(myphi, mytheta) \
  theta = mytheta;                 \
  atom3[0] = cos(theta * DEG2RAD); \
  atom3[1] = sin(theta * DEG2RAD); \
  myphi = angle->single(atype, i2, i1, i3)

    const double dtheta = 180.0 / static_cast<double>(n - 1);

    // get force for divergent 0 degree angle from interpolation to the right

    GET_ENERGY(phi, 0.0);
    GET_ENERGY(phi1, epsilon);
    GET_ENERGY(phi2, 2.0 * epsilon);

    f = (1.5 * phi - 2.0 * phi1 + 0.5 * phi2) / epsilon;
    if (!std::isfinite(f)) f = 0.0;
    if (!std::isfinite(phi)) phi = 0.0;
    fprintf(fp, "%8d %- 22.15g %- 22.15g %- 22.15g\n", 1, 0.0, phi, f);

    for (int i = 1; i < n - 1; i++) {
      GET_ENERGY(phi1, dtheta * static_cast<double>(i) - epsilon);
      GET_ENERGY(phi2, dtheta * static_cast<double>(i) + epsilon);
      GET_ENERGY(phi, dtheta * static_cast<double>(i));

      if (!std::isfinite(phi)) phi = 0.0;

      // get force from numerical differentiation
      f = -0.5 * (phi2 - phi1) / epsilon;
      if (!std::isfinite(f)) f = 0.0;
      fprintf(fp, "%8d %- 22.15g %- 22.15g %- 22.15g\n", i + 1, theta, phi, f);
    }

    // get force for divergent 180 degree angle from interpolation to the left
    GET_ENERGY(phi, 180.0);
    GET_ENERGY(phi1, 180.0 - epsilon);
    GET_ENERGY(phi2, 180.0 - 2.0 * epsilon);

    f = (2.0 * phi1 - 1.5 * phi - 0.5 * phi2) / epsilon;
    if (!std::isfinite(f)) f = 0.0;
    if (!std::isfinite(phi)) phi = 0.0;
    fprintf(fp, "%8d %- 22.15g %- 22.15g %- 22.15g\n", 1, 180.0, phi, f);

    // clean up
    delete writer;
    fclose(fp);
  }
  MPI_Comm_free(&singlecomm);
}
