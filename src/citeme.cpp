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

#include "citeme.h"
#include "comm.h"
#include "universe.h"

#include <functional>    // IWYU pragma: keep

using namespace LAMMPS_NS;

namespace {
/** Separator line for citation output blocks */
constexpr char cite_separator[] =
    "CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE\n\n";

/** Header message for citation output */
constexpr char cite_nagline[] =
    "Your simulation uses code contributions which should be cited:\n";

/** Format string for citation file reference */
constexpr char cite_file[] = "The {} {} lists these citations in BibTeX format.\n\n";

/** Hash function for deduplicating citations */
std::hash<std::string> get_hash;
}    // namespace
/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Constructor initializes citation tracking and optionally opens citation file

   Only MPI rank 0 in the universe opens and writes to the citation file.
   All other ranks simply initialize the data structures but do not perform I/O.

   \param lmp      Pointer to main LAMMPS object
   \param _screen  Output mode for screen (VERBOSE or TERSE)
   \param _logfile Output mode for log file (VERBOSE or TERSE)
   \param _file    Optional filename for BibTeX output (NULL for no file)
------------------------------------------------------------------------- */

CiteMe::CiteMe(LAMMPS *lmp, int _screen, int _logfile, const char *_file) : Pointers(lmp)
{
  fp = nullptr;

  screen_flag = _screen;
  scrbuffer.clear();
  logfile_flag = _logfile;
  logbuffer.clear();

  if (_file && universe->me == 0) {
    citefile = _file;
    fp = fopen(_file, "w");
    if (fp) {
      fputs(cite_nagline, fp);
      fflush(fp);
    } else {
      utils::logmesg(
          lmp, "Unable to open citation file '" + citefile + "': " + utils::getsyserror() + "\n");
    }
  }
}

/* ----------------------------------------------------------------------
   Destructor flushes any remaining citations and closes citation file

   Ensures all pending citation output is written before cleanup.
------------------------------------------------------------------------- */

CiteMe::~CiteMe()
{
  flush();
  if (fp) fclose(fp);
}

/* ----------------------------------------------------------------------
   Add a citation to the output queue

   This method:
   - Uses hash-based deduplication to show each citation only once
   - Writes immediately to the BibTeX file (if enabled)
   - Buffers citations for screen and log file output
   - Handles both VERBOSE (full BibTeX) and TERSE (first line only) modes

   Only MPI rank 0 in the communicator performs the actual work to avoid
   duplicate output in parallel runs.

   \param reference String containing citation in BibTeX format with DOI header
------------------------------------------------------------------------- */

void CiteMe::add(const std::string &reference)
{
  if (comm->me != 0) return;

  std::size_t crc = get_hash(reference);
  if (cs.find(crc) != cs.end()) return;
  cs.insert(crc);

  if (fp) {
    fputs(reference.c_str(), fp);
    fflush(fp);
  }

  if (scrbuffer.empty()) {
    scrbuffer += "\n";
    scrbuffer += cite_separator;
    scrbuffer += cite_nagline;
    if (screen_flag == VERBOSE) scrbuffer += "\n";
  }

  if (logbuffer.empty()) {
    logbuffer += "\n";
    logbuffer += cite_separator;
    logbuffer += cite_nagline;
    if (logfile_flag == VERBOSE) logbuffer += "\n";
  }

  std::size_t found = reference.find_first_of('\n');
  std::string header = reference.substr(0, found + 1);
  if (screen_flag == VERBOSE) scrbuffer += "- " + reference;
  if (screen_flag == TERSE) scrbuffer += "- " + header;
  if (logfile_flag == VERBOSE) logbuffer += "- " + reference;
  if (logfile_flag == TERSE) logbuffer += "- " + header;
}

/* ----------------------------------------------------------------------
   Flush accumulated citation buffers to output streams

   Writes the buffered citations to screen and log file with appropriate
   formatting and separator lines. Clears buffers after output.

   Only MPI rank 0 performs output to avoid duplication in parallel runs.
------------------------------------------------------------------------- */

void CiteMe::flush()
{
  if (comm->me == 0) {
    if (!scrbuffer.empty()) {
      if (!citefile.empty()) scrbuffer += fmt::format(cite_file, "file", citefile);
      if (logfile_flag == VERBOSE) scrbuffer += fmt::format(cite_file, "log", "file");
      scrbuffer += cite_separator;
      if (screen) fputs(scrbuffer.c_str(), screen);
      scrbuffer.clear();
    }
    if (!logbuffer.empty()) {
      if (!citefile.empty()) logbuffer += fmt::format(cite_file, "file", citefile);
      if (screen_flag == VERBOSE) logbuffer += fmt::format(cite_file, "screen", "output");
      logbuffer += cite_separator;
      if (logfile) fputs(logbuffer.c_str(), logfile);
      logbuffer.clear();
    }
  }
}
