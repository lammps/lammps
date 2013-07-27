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
#include "comm.h"

#include <stdio.h>
#include <set>

using namespace LAMMPS_NS;

// to simplify the code below
typedef std::set<const char *> citeset;

static const char nagline[] = "\n"
  "------------------------------------------------------------------------\n"
  "This simulation made use of algorithms and methodologies described\n"
  "in the references listed in the file 'log.cite'\n"
  "------------------------------------------------------------------------\n";

static const char cite_header[] = "\n"
  "------------------------------------------------------------------------\n"
  "This simulation made use of algorithms and methodologies described\n"
  "in the following references:\n\n";

static const char lammps_version[] =
  "The LAMMPS Molecular Dynamics Simulator, Version " LAMMPS_VERSION "\n"
  "    http://lammps.sandia.gov\n\n";

/* ---------------------------------------------------------------------- */

CiteMe::CiteMe(LAMMPS *lmp) : Pointers(lmp) {

  FILE *fp;
  citeset *c = new citeset;

  _pubs = (void *)c;

  if ((universe->me == 0) && ((fp = fopen("log.cite","w")))) {
    fputs(cite_header,fp);
    fputs(lammps_version,fp);
    fflush(fp);
  } else {
    fp = NULL;
  }
  _fp = (void *) fp;

}

/* ----------------------------------------------------------------------
   Write out and register a citation so it will be written only once
------------------------------------------------------------------------- */

void CiteMe::add(const char *ref)
{
  citeset *c = (citeset *) _pubs;

  if (c->find(ref) == c->end()) {

    FILE *fp = (FILE *)_fp;
    fputs(ref,fp);
    fflush(fp);

    c->insert(ref);
  }
}

/* ---------------------------------------------------------------------- 
   Write out nag-line at the end of the regular output and clean up.
------------------------------------------------------------------------- */

CiteMe::~CiteMe()
{

  if (comm->me == 0) {
    if (screen)
      fputs(nagline,screen);

    if (logfile)
      fputs(nagline,logfile);
  } 

  FILE *fp = (FILE *)_fp;
  if (fp) fclose(fp)

  citeset *c = (citeset *) _pubs;
  delete c;
}
