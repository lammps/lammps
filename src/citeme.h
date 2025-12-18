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

/** \class LAMMPS_NS::CiteMe
 *
 * The CiteMe class provides a mechanism for LAMMPS to remind users to cite
 * relevant publications when they use specific contributed features. This
 * ensures proper attribution for the scientific work underlying LAMMPS
 * implementations. */
class CiteMe : protected Pointers {
 public:
  /** Constructor for CiteMe class
   *
   * \param lmp        Pointer to the main LAMMPS object
   * \param _screen    Output mode for screen (VERBOSE or TERSE)
   * \param _logfile   Output mode for log file (VERBOSE or TERSE)
   * \param _file      Optional filename for BibTeX citation output file (can be NULL) */
  CiteMe(class LAMMPS *, int, int, const char *);

  /** Destructor flushes any remaining citations */
  ~CiteMe() override;

  /** Register a publication for citation output

\verbatim embed:rst

Adds a citation to the set of publications to be cited.  Each citation
should contain a BibTeX entry and is output only once, even if add()
is called multiple times with the same citation.  The citation string
should start with a brief description including a DOI,  followed by a
complete BibTeX entry.

This method should typically be called in the constructor of a style
that implements a published method or algorithm.

\endverbatim

   * \param reference  String containing the citation in BibTeX format with DOI header */
  void add(const std::string &reference);

  /** Flush accumulated citation buffers to screen and log file

\verbatim embed:rst

Outputs all pending citations to the screen and log file with appropriate
formatting.  Called automatically by the destructor when LAMMPS terminates
or is reset by the :doc:`clear <clear>` command and at the end of a
:doc:`run <run>` or :doc:`minimize <minimize>` command.

\endverbatim
  */
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
