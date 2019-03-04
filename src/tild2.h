/* ------------------------------------------------------------------
    TILD - Theoretically Informed Langevan Dynamics
    This replicates the TILD coe done by the Riggleman group, 
    previously known as Dynamical Mean Field Theory. 
    
    Copyright (2019) Christian Tabedzki and Zachariah Vicars.
    tabedzki@seas.upenn.edu zvicars@seas.upenn.edu
-------------------------------------------------------------------- */

#ifdef KSPACE_CLASS

KSpaceStyle(tild,TILD)

#else

#ifndef LMP_TILD_H
#define LMP_TILD_H

#include "pppm.h"

namespace LAMMPS_NS {

class TILD : public PPPM {
 public:
  TILD (class LAMMPS *);
  virtual ~TILD();
  virtual void init();
  virtual void setup();
  virtual void settings(int, char **);
  virtual void compute(int, int);
  double memory_usage();

  void compute_group_group(int, int, int);

 protected:
  static double *uG;
  static double a_squared;
  double **grad_uG, **grad_uG_hat;
  int kxmax,kymax,kzmax;
  int kcount,kmax,kmax3d,kmax_created;
  double gsqmx,volume;
  int nmax;

  double unitk[3];
  int *kxvecs,*kyvecs,*kzvecs;
  int kxmax_orig,kymax_orig,kzmax_orig;
  double *ug;
  double **eg,**vg;
  double **ek;
  double *sfacrl,*sfacim,*sfacrl_all,*sfacim_all;
  double ***cs,***sn;

  // group-group interactions

  int group_allocate_flag;
  double *sfacrl_A,*sfacim_A,*sfacrl_A_all,*sfacim_A_all;
  double *sfacrl_B,*sfacim_B,*sfacrl_B_all,*sfacim_B_all;

  double rms(int, double, bigint, double);
  // virtual void eik_dot_r();
  void coeffs();
  virtual void allocate();
  void deallocate();
  void slabcorr();
  void init_gauss();

  // group-group interactions

//   void slabcorr_groups(int,int,int);
//   void allocate_groups();
//   void deallocate_groups();
};

}

#endif
#endif

