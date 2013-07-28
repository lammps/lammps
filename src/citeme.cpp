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

typedef std::set<const char *> citeset;

static const char dashline[] =   "----------------------------------"
  "----------------------------------------------\n";

static const char nagline[] = "\nPlease see the file 'log.cite'"
  " for references relevant to this simulation\n\n";

static const char cite_header[] = "\nThis simulation made use of "
  "methodologies described in the following references:\n\n";

static const char lammps_version[] =
  "The LAMMPS Molecular Dynamics Simulator, Version " LAMMPS_VERSION
  ", lammps.sandia.gov\n";

/* ---------------------------------------------------------------------- */

CiteMe::CiteMe(LAMMPS *lmp) : Pointers(lmp) {

  FILE *fp;
  citeset *c;

  _pubs = (void *)c;

  if (lmp->cite_enable && (universe->me == 0) && ((fp = fopen("log.cite","w")))) {
    fputs(dashline,fp);
    fputs(lammps_version,fp);
    fputs(dashline,fp);
    fputs(cite_header,fp);
    fflush(fp);
    c = new citeset;
  } else {
    fp = NULL;
    c  = NULL;
  }
  _fp = (void *) fp;
  _pubs = (void *) c;
}

/* ----------------------------------------------------------------------
   Write out and register a citation so it will be written only once
------------------------------------------------------------------------- */

void CiteMe::add(const char *ref)
{
  if (_fp == NULL) return;

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

  if (lmp->cite_enable && (comm->me == 0)) {
    if (screen)
      fputs(nagline,screen);

    if (logfile)
      fputs(nagline,logfile);
  } 

  if (_fp != NULL) {
    
    FILE *fp = (FILE *)_fp;
    fputs(dashline,fp);
    fclose(fp);

    citeset *c = (citeset *) _pubs;
    delete c;
  }
}
