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

  /** Register a publication for citation output */

  void add(const std::string &);

  /** Flush accumulated citation text to screen and log file */

  void flush();

  enum {
    VERBOSE, /**< Display full BibTeX entries */
    TERSE    /**< Display only first line (DOI/brief description) */
  };

 private:
  FILE *fp;                 /**< File pointer for optional BibTeX citation file or NULL */
  std::string citefile;     /**< Name of the explicit citation file */
  int screen_flag;          /**< Output mode for screen (VERBOSE or TERSE) */
  int logfile_flag;         /**< Output mode for log file (VERBOSE or TERSE) */
  std::string scrbuffer;    /**< Output buffer for screen */
  std::string logbuffer;    /**< Output buffer for log file */
  std::set<std::size_t> cs; /**< Set of citation hashes to track uniqueness */
};
}    // namespace LAMMPS_NS

#endif
