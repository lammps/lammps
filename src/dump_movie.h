/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(movie,DumpMovie)

#else

#ifndef LMP_DUMP_MOVIE_H
#define LMP_DUMP_MOVIE_H

#include "dump_image.h"

namespace LAMMPS_NS {

class DumpMovie : public DumpImage {
 public:
  DumpMovie(LAMMPS *, int, char**);

  virtual void openfile();
  virtual void init_style();

 protected:
  int bitrate;                  // bitrate of video file in kbps
  int framerate;                // frame rate of animation
};

}

#endif
#endif

/* ERROR/WARNING messages:


*/
