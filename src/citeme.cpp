/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
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

using namespace LAMMPS_NS;

static const char cite_header[] =
  "This LAMMPS simulation made specific use of work described in the\n"
  "following references.  See https://lammps.sandia.gov/cite.html\n"
  "for details.\n\n";

static const char cite_nagline[] = "Please see the log.cite file "
  "for references relevant to this simulation\n\n";

static const char cite_seefile[] = "Please see the citation file "
  "for references relevant to this simulation\n\n";

static const char cite_separator[] =
  "\nCITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE-CITE\n\n";

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
    fp = fopen(_file,"w");
    if (fp) {
      fputs(cite_header,fp);
      fflush(fp);
    } else {
      utils::logmesg(lmp, "Unable to open citation file '" + std::string(_file)
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

void CiteMe::add(const char *ref)
{
  if (comm->me != 0) return;
  if (cs->find(ref) != cs->end()) return;
  cs->insert(ref);

  if (fp) {
    fputs(ref,fp);
    fflush(fp);
  }

  if (scrbuffer.empty()) {
    scrbuffer += cite_separator;
    if (screen_flag == VERBOSE) scrbuffer += cite_header;
    if (screen_flag == TERSE) scrbuffer += cite_nagline;
  }
  if (logbuffer.empty()) {
    logbuffer += cite_separator;
    if (logfile_flag == VERBOSE) logbuffer += cite_header;
    if (logfile_flag == TERSE) logbuffer += cite_nagline;
  }

  std::string reference = ref;
  std::size_t found = reference.find_first_of("\n");
  std::string header = reference.substr(0,found+1);
  if (screen_flag == VERBOSE) scrbuffer += reference;
  if (screen_flag == TERSE) scrbuffer += header;
  if (logfile_flag == VERBOSE) logbuffer += reference;
  if (logfile_flag == TERSE) logbuffer += header;
}

void CiteMe::flush()
{
  if (comm->me == 0) {
    if (!scrbuffer.empty()) {
      scrbuffer += cite_separator;
      if (screen) fputs(scrbuffer.c_str(),screen);
      scrbuffer.clear();
    }
    if (!scrbuffer.empty()) {
      logbuffer += cite_separator;
      if (logfile) fputs(logbuffer.c_str(),logfile);
      logbuffer.clear();
    }
  }
  return;
}

