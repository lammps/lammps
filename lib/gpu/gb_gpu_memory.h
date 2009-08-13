/***************************************************************************
                               gb_gpu_memory.h
                             -------------------
                               W. Michael Brown

  Global variables for GPU Gayberne Library

 __________________________________________________________________________
    This file is part of the LAMMPS GPU Library
 __________________________________________________________________________

    begin                : Thu Jun 25 2009
    copyright            : (C) 2009 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
   ----------------------------------------------------------------------- */

#ifndef GB_GPU_MEMORY_H
#define GB_GPU_MEMORY_H

#define MAX_GPU_THREADS 4
#include "lj_gpu_memory.h"

enum{SPHERE_SPHERE,SPHERE_ELLIPSE,ELLIPSE_SPHERE,ELLIPSE_ELLIPSE};

template <class numtyp, class acctyp>
class GB_GPU_Memory : public LJ_GPU_Memory<numtyp,acctyp> {
 public:
  GB_GPU_Memory();
  ~GB_GPU_Memory(); 
 
  int* init(const int ij_size, const int ntypes, const double gamma,
            const double upsilon, const double mu, double **host_shape,
            double **host_well, double **host_cutsq, double **host_sigma, 
            double **host_epsilon, double *host_lshape, int **h_form,
            double **host_lj1, double **host_lj2, double **host_lj3, 
            double **host_lj4, double **host_offset, double *host_special_lj,
            const int max_nbors, const bool force_d, const int me);

  void clear();
 
  double host_memory_usage();
  
  // ----------------------------  DATA  ----------------------------

  // ilist with particles sorted by type
  NVC_HostI host_olist;

  // --------------- Const Data for Atoms
  NVC_ConstMatT shape, well;
  NVC_ConstMatI form;
  NVC_VecT lshape, gamma_upsilon_mu;
  

  // --------------- Timing Stuff
  NVCTimer time_kernel, time_gayberne, time_kernel2, time_gayberne2;
  
  // True if we want to use fast GB-sphere or sphere-sphere calculations 
  bool multiple_forms;
  int **host_form;
  int last_ellipse;
   
 private:
};

#endif

