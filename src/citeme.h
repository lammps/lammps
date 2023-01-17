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

#ifndef LMP_CITEME_H
#define LMP_CITEME_H

#include "pointers.h"
#include <set>

namespace LAMMPS_NS {

class CiteMe : protected Pointers {
 public:
  CiteMe(class LAMMPS *, int, int, const char *);
  ~CiteMe() override;
  void add(const std::string &);    // register publication for output
  void flush();                     // flush buffers to screen and logfile
  enum { VERBOSE, TERSE };

 private:
  FILE *fp;                 // explicit citation file pointer or NULL
  std::string citefile;     // name of the explicit citation file.
  int screen_flag;          // determine whether verbose or terse output
  int logfile_flag;         // determine whether verbose or terse output
  std::string scrbuffer;    // output buffer for screen
  std::string logbuffer;    // output buffer for logfile
  typedef std::set<std::size_t> citeset;
  citeset *cs;    // registered set of publications
};
}    // namespace LAMMPS_NS

#endif
