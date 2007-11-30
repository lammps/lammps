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

#ifndef PPPM_H
#define PPPM_H

#include "kspace.h"

namespace LAMMPS_NS {

class PPPM : public KSpace {
 public:
  PPPM(class LAMMPS *, int, char **);
  ~PPPM();
  void init();
  void setup();
  void compute(int, int);
  void timing(int, double &, double &);
  double memory_usage();

 protected:
  int me,nprocs;
  double PI;
  double precision;
  int nfactors;
  int *factors;
  double qsum,qsqsum;
  double qqrd2e;
  double cutoff;
  double volume;
  double delxinv,delyinv,delzinv,delvolinv;
  double shift,shiftone;

  int nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in;
  int nxlo_out,nylo_out,nzlo_out,nxhi_out,nyhi_out,nzhi_out;
  int nxlo_ghost,nxhi_ghost,nylo_ghost,nyhi_ghost,nzlo_ghost,nzhi_ghost;
  int nxlo_fft,nylo_fft,nzlo_fft,nxhi_fft,nyhi_fft,nzhi_fft;
  int nlower,nupper;
  int ngrid,nfft,nbuf,nfft_both;

  double ***density_brick;
  double ***vdx_brick,***vdy_brick,***vdz_brick;
  double *greensfn;
  double **vg;
  double *fkx,*fky,*fkz;
  double *density_fft;
  double *work1,*work2;
  double *buf1,*buf2;

  double *gf_b;
  double **rho1d,**rho_coeff;

  class FFT3d *fft1,*fft2;
  class Remap *remap;

  int **part2grid;             // storage for particle -> grid mapping
  int nmax;

  int triclinic;               // domain settings, orthog or triclinic
  double *boxlo;
                               // TIP4P settings
  int typeH,typeO;             // atom types of TIP4P water H and O atoms
  double qdist;                // distance from O site to negative charge
  double alpha;                // geometric factor

  void set_grid();
  void allocate();
  void deallocate();
  int factorable(int);
  double rms(double, double, double, double, double **);
  double diffpr(double, double, double, double, double **);
  void compute_gf_denom();
  double gf_denom(double, double, double);
  virtual void particle_map();
  virtual void make_rho();
  void brick2fft();
  void fillbrick();
  void poisson(int, int);
  virtual void fieldforce();
  void procs2grid2d(int,int,int,int *, int*);
  void compute_rho1d(double, double, double);
  void compute_rho_coeff();
  void slabcorr(int);
};

}

#endif
