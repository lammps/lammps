// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "citeme.h"
#include "comm.h"
#include "universe.h"

#include <functional>           // IWYU pragma: keep

using namespace LAMMPS_NS;

static const char cite_separator[] =
  "CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE\n\n";

static const char cite_nagline[] =
  "Your simulation uses code contributions which should be cited:\n";

static const char cite_file[] = "The {} {} lists these citations in "
                                "BibTeX format.\n\n";

// define hash function
static std::hash<std::string> get_hash;

/* ---------------------------------------------------------------------- */

CiteMe::CiteMe(LAMMPS *lmp, int _screen, int _logfile, const char *_file)
  : Pointers(lmp)
{
  fp = nullptr;
  cs = new citeset();

  screen_flag = _screen;
  scrbuffer.clear();
  logfile_flag = _logfile;
  logbuffer.clear();

  if (_file && universe->me == 0) {
    citefile = _file;
    fp = fopen(_file,"w");
    if (fp) {
      fputs(cite_nagline,fp);
      fflush(fp);
    } else {
      utils::logmesg(lmp, "Unable to open citation file '" + citefile
                     + "': " + utils::getsyserror() + "\n");
    }
  }
}

/* ----------------------------------------------------------------------
   write out remaining citations at end of the regular output and clean up
------------------------------------------------------------------------- */

CiteMe::~CiteMe()
{
  flush();
  delete cs;

  if (fp) fclose(fp);
}

/* ----------------------------------------------------------------------
   process an added citation so it will be shown only once and as requested
------------------------------------------------------------------------- */

void CiteMe::add(const std::string &reference)
{
  if (comm->me != 0) return;

  std::size_t crc = get_hash(reference);
  if (cs->find(crc) != cs->end()) return;
  cs->insert(crc);

  if (fp) {
    fputs(reference.c_str(),fp);
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
  std::string header = reference.substr(0,found+1);
  if (screen_flag == VERBOSE) scrbuffer += "- " + reference;
  if (screen_flag == TERSE) scrbuffer += "- " + header;
  if (logfile_flag == VERBOSE) logbuffer += "- " + reference;
  if (logfile_flag == TERSE) logbuffer += "- " + header;
}

void CiteMe::flush()
{
  if (comm->me == 0) {
    if (!scrbuffer.empty()) {
      if (!citefile.empty())
        scrbuffer += fmt::format(cite_file,"file",citefile);
      if (logfile_flag == VERBOSE)
        scrbuffer += fmt::format(cite_file,"log","file");
      scrbuffer += cite_separator;
      if (screen) fputs(scrbuffer.c_str(),screen);
      scrbuffer.clear();
    }
    if (!logbuffer.empty()) {
      if (!citefile.empty())
        logbuffer += fmt::format(cite_file,"file",citefile);
      if (screen_flag == VERBOSE)
        logbuffer += fmt::format(cite_file,"screen","output");
      logbuffer += cite_separator;
      if (logfile) fputs(logbuffer.c_str(),logfile);
      logbuffer.clear();
    }
  }
  return;
}

