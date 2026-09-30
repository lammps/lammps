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
   Contributing author:  Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "dump_movie.h"

#include "comm.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpMovie::DumpMovie(LAMMPS *lmp, int narg, char **arg) : DumpImage(lmp, narg, arg)
{
  if (multiproc || compressed || multifile) error->all(FLERR, "Invalid dump movie filename");

  filetype = PPM;
  bitrate = 2000;
  framerate = 24;
}

/* ---------------------------------------------------------------------- */

void DumpMovie::openfile()
{
  if ((comm->me == 0) && (fp == nullptr)) {
    auto ffmpeg = platform::find_exe_path("ffmpeg");
    if (ffmpeg.empty()) {
      error->one(FLERR, Error::NOLASTLINE, "Dump movie requires 'ffmpeg' but it could not be found");
    } else {
      auto moviecmd = fmt::format("{} -v error -y -r {:.2f} -f image2pipe -c:v ppm -i - -r 24.0 -b:v {}k {}",
                                  ffmpeg, framerate, bitrate, filename);
      fp.set_pclose();
      fp = platform::popen(moviecmd, "w");
    }
  }
}
/* ---------------------------------------------------------------------- */

void DumpMovie::init_style()
{
  // initialize image style circumventing multifile check

  multifile_override = 1;
  DumpImage::init_style();
}

/* ---------------------------------------------------------------------- */

int DumpMovie::modify_param(int narg, char **arg)
{
  int n = DumpImage::modify_param(narg, arg);
  if (n) return n;

  if (strcmp(arg[0], "bitrate") == 0) {
    if (narg < 2) error->all(FLERR, "Illegal dump_modify command");
    bitrate = utils::inumeric(FLERR, arg[1], false, lmp);
    if (bitrate <= 0.0) error->all(FLERR, "Illegal dump_modify command");
    return 2;
  }

  if (strcmp(arg[0], "framerate") == 0) {
    if (narg < 2) error->all(FLERR, "Illegal dump_modify command");
    framerate = utils::numeric(FLERR, arg[1], false, lmp);
    if ((framerate <= 0.1) || (framerate > 24.0))
      error->all(FLERR, "Illegal dump_modify framerate command");
    return 2;
  }

  return 0;
}
