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
 *  \brief Manages citation reminders for contributed features in LAMMPS
 *
 * The CiteMe class provides a mechanism for LAMMPS to remind users to cite
 * relevant publications when they use specific contributed features. This
 * ensures proper attribution for the scientific work underlying LAMMPS
 * implementations.
 *
 * ## Overview
 *
 * When users enable and use contributed packages or special features in
 * LAMMPS, the CiteMe class tracks which publications should be cited and
 * displays appropriate reminders during the simulation run.
 *
 * ## Usage in Contributed Code
 *
 * To add a citation reminder to a contributed feature, developers should:
 *
 * 1. Define a static string containing the citation in BibTeX format at
 *    the top of the implementation file:
 * \code{.cpp}
 *    static const char cite_my_feature[] =
 *      "my_feature command: doi:10.1234/example.doi\n\n"
 *      "@Article{AuthorYear,\n"
 *      " author = {First Author and Second Author},\n"
 *      " title = {Title of Paper},\n"
 *      " journal = {Journal Name},\n"
 *      " year = 2024,\n"
 *      " volume = 100,\n"
 *      " pages = {1-10}\n"
 *      "}\n\n";
 * \endcode
 *
 * 2. Call the add() method in the style's constructor to register the citation:
 * \code{.cpp}
 *    MyStyle::MyStyle(LAMMPS *lmp) : Parent(lmp)
 *    {
 *      if (lmp->citeme) lmp->citeme->add(cite_my_feature);
 *      // ... rest of constructor
 *    }
 * \endcode
 *
 * ## Citation Output
 *
 * Citations are output in three ways:
 * - To the screen (controlled by the VERBOSE/TERSE flag)
 * - To the log file (controlled by a separate VERBOSE/TERSE flag)
 * - To an optional citation file in BibTeX format (if specified)
 *
 * The verbosity can be controlled via command-line arguments when starting
 * LAMMPS. In VERBOSE mode, the full BibTeX entry is displayed. In TERSE mode,
 * only the first line (typically containing the DOI) is shown.
 *
 * ## Example Output
 *
 * \code{.unparsed}
 * CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE
 *
 * Your simulation uses code contributions which should be cited:
 *
 * - pair gayberne command: doi:10.1063/1.3058435
 *
 * The file log.cite lists these citations in BibTeX format.
 *
 * CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE
 * \endcode
 *
 * \see add()
 * \see flush()
 */
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
   *
   * Adds a citation to the set of publications to be cited. Each citation
   * should contain a BibTeX entry and is output only once, even if add()
   * is called multiple times with the same citation. The citation string
   * should start with a brief description including a DOI, followed by a
   * complete BibTeX entry.
   *
   * This method should typically be called in the constructor of a style
   * that implements a published method or algorithm.
   *
   * \param reference  String containing the citation in BibTeX format with DOI header */
  void add(const std::string &);

  /** Flush accumulated citation buffers to screen and log file
   *
   * Outputs all pending citations to the screen and log file with
   * appropriate formatting. Called automatically by the destructor,
   * but can also be called explicitly to output citations at specific
   * points during execution. */
  void flush();

  /** Output verbosity modes */
  enum { VERBOSE, /**< Display full BibTeX entries */
         TERSE   /**< Display only first line (DOI/brief description) */ };

 private:
  FILE *fp;                 /**< File pointer for optional BibTeX citation file or NULL */
  std::string citefile;     /**< Name of the explicit citation file */
  int screen_flag;          /**< Output mode for screen (VERBOSE or TERSE) */
  int logfile_flag;         /**< Output mode for log file (VERBOSE or TERSE) */
  std::string scrbuffer;    /**< Output buffer for screen */
  std::string logbuffer;    /**< Output buffer for log file */
  using citeset = std::set<std::size_t>;
  citeset *cs;              /**< Set of citation hashes to track uniqueness */
};
}    // namespace LAMMPS_NS

#endif
