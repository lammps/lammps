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

#ifdef PAIR_CLASS                    

PairStyle(peri/lps,PairPeriLPS)      

#else                                 

#ifndef LMP_PAIR_PERI_LPS_H          
#define LMP_PAIR_PERI_LPS_H          

#include "pair.h"

namespace LAMMPS_NS {

class PairPeriLPS : public Pair {    
 public:
  PairPeriLPS(class LAMMPS *);       
  virtual ~PairPeriLPS();                                   
  int pack_comm(int, int *, double *, int, int *);  
  void unpack_comm(int, int, double *);             

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *) {}
  void read_restart_settings(FILE *) {}
  double single(int, int, int, int, double, double, double, double &);
  double memory_usage();
  double influence_function(double, double, double);
  void compute_dilatation();   

 protected:
  int ifix_peri;
  double **bulkmodulus;                
  double **shearmodulus;               
  double **s00, **alpha;                  
  double **cut;
 
  double *s0_new;                  
  double *theta;                      
  int nmax;

  void allocate();
};

}

#endif
#endif
