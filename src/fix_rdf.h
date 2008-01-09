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

#ifndef FIX_RDF_H
#define FIX_RDF_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixRDF : public Fix {
 public:
  FixRDF(class LAMMPS *, int, char **);
  ~FixRDF();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void end_of_step();

 private:
  int me,first;
  FILE *fp;
  int maxbin;			 // # of rdf bins
  int n_rdf_pairs;            	 // # of rdf pairs
  int nframes;                   // # of rdf snapshots taken
  double delr,delrinv;		 // bin width and its inverse
  int **rdfpair;              	 // mapping from 2 types to rdf pair
  int **hist,**hist_all;	 // histogram bins
  int *nrdfatoms;             	 // # of atoms of each type in the group
  double **gr_ave,**ncoord_ave;  // accumulators for average rdf statistics
  class NeighList *list;         // half neighbor list
};

}

#endif
