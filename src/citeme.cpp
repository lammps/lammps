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

#include <stdio.h>
#include <set>

using namespace LAMMPS_NS;

// the list of publications is below
static const char * const publication[] = {
  /* PLIMPTON_1995 */  "S. Plimpton,"
  " Fast Parallel Algorithms for Short-Range Molecular Dynamics,\n"
  "  J Comp Phys, 117, 1-19 (1995)\n\n",
  NULL 
};

/* ---------------------------------------------------------------------- */

CiteMe::CiteMe(LAMMPS *lmp) : Pointers(lmp) {
  list = static_cast<std::set<int> *>(new std::set<int>);
}

/* ---------------------------------------------------------------------- */

void CiteMe::add(int ref)
{
  std::set<int> *c = static_cast< std::set<int> *>(list);
  c->insert(ref);
}

/* ---------------------------------------------------------------------- */

static const char cite_header[] = "\n"
  "------------------------------------------------------------------------\n"
  "This simulation made use of algorithms and methodologies described\n"
  "in the following references:\n\n"
  "The LAMMPS Molecular Dynamics Simulator, Version " LAMMPS_VERSION "\n"
  "    http://lammps.sandia.gov\n\n";


CiteMe::~CiteMe(){
  std::set<int> *c = static_cast< std::set<int> *>(list);

  if (screen)
    fputs(cite_header,screen);

  if (logfile)
    fputs(cite_header,logfile);

  for (std::set<int>::const_iterator i = c->begin(); i != c->end(); ++i) {
    if (screen)
      fputs(publication[*i],screen);
    if (logfile)
      fputs(publication[*i],logfile);
  }

  delete c;
}


