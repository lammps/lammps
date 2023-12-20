/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(movie,DumpMovie);
// clang-format on
#else

#ifndef LMP_DUMP_MOVIE_H
#define LMP_DUMP_MOVIE_H

#include "dump_image.h"

namespace LAMMPS_NS {

class DumpMovie : public DumpImage {
 public:
  DumpMovie(LAMMPS *, int, char **);
  ~DumpMovie() override;

  void openfile() override;
  void init_style() override;
  int modify_param(int, char **) override;

 protected:
  double framerate;    // frame rate of animation
  int bitrate;         // bitrate of video file in kbps
};

}    // namespace LAMMPS_NS

#endif
#endif
