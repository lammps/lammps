/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

  virtual void openfile();
  virtual void init_style();
  virtual int modify_param(int, char **);

 protected:
  double framerate;    // frame rate of animation
  int bitrate;         // bitrate of video file in kbps
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid dump movie filename

The file produced by dump movie cannot be binary or compressed
and must be a single file for a single processor.

E: Support for writing movies not included

LAMMPS was not built with the -DLAMMPS_FFMPEG switch in the Makefile

E: Failed to open FFmpeg pipeline to file %s

The specified file cannot be opened.  Check that the path and name are
correct and writable and that the FFmpeg executable can be found and run.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
