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

#ifndef COMPUTE_RDF_H
#define COMPUTE_RDF_H

#include "stdio.h"
#include "compute.h"

namespace LAMMPS_NS {

class ComputeRDF : public Compute {
 public:
  ComputeRDF(class LAMMPS *, int, char **);
  ~ComputeRDF();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();

 private:
  int first;
  int nbin;			 // # of rdf bins
  int npairs;            	 // # of rdf pairs
  double delr,delrinv;		 // bin width and its inverse
  int ***rdfpair;              	 // map 2 type pair to rdf pair for each histo
  int **nrdfpair;                // # of histograms for each type pair
  int *ilo,*ihi,*jlo,*jhi;
  double **hist;	         // histogram bins
  double **histall;	         // summed histogram bins across all procs

  int *typecount;
  int *icount,*jcount;

  class NeighList *list;         // half neighbor list
};

}

#endif
