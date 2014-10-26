/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "citeme.h"
#include "version.h"
#include "universe.h"
#include "error.h"

using namespace LAMMPS_NS;

static const char cite_header[] = 
  "This LAMMPS simulation made specific use of work described in the\n"
  "following references.  See http://lammps.sandia.gov/cite.html\n"
  "for details.\n\n";

static const char cite_nagline[] = "\nPlease see the log.cite file "
  "for references relevant to this simulation\n\n";

/* ---------------------------------------------------------------------- */

CiteMe::CiteMe(LAMMPS *lmp) : Pointers(lmp)
{
  fp = NULL;
  cs = new citeset();
}

/* ---------------------------------------------------------------------- 
   write out nag-line at the end of the regular output and clean up
------------------------------------------------------------------------- */

CiteMe::~CiteMe()
{
  if (universe->me || cs->size() == 0) {
    delete cs;
    return;
  }

  delete cs;

  if (fp) {
    if (screen) fprintf(screen,cite_nagline);
    if (logfile) fprintf(logfile,cite_nagline);

    fclose(fp);
  }
}

/* ----------------------------------------------------------------------
   write out and register a citation so it will be written only once
------------------------------------------------------------------------- */

void CiteMe::add(const char *ref)
{
  if (universe->me) return;
  if (cs->find(ref) != cs->end()) return;
  cs->insert(ref);

  if (!fp) {
    fp = fopen("log.cite","w");
    if (!fp) return;
    fputs(cite_header,fp);
    fflush(fp);
  }

  fputs(ref,fp);
  fflush(fp);
}
