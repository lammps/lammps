/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS

KSpaceStyle(pppm,PPPM)

#else

#ifndef LMP_PPPM_H
#define LMP_PPPM_H

#include "lmptype.h"
#include "mpi.h"

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_FLOAT
#else
typedef double FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

#include "kspace.h"

namespace LAMMPS_NS {

class PPPM : public KSpace {
 public:
  PPPM(class LAMMPS *, int, char **);
  virtual ~PPPM();
  virtual void init();
  virtual void setup();
  virtual void compute(int, int);
  virtual void timing(int, double &, double &);
  virtual double memory_usage();

 protected:
  int me,nprocs;
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

  FFT_SCALAR ***density_brick;
  FFT_SCALAR ***vdx_brick,***vdy_brick,***vdz_brick;
  double *greensfn;
  double **vg;
  double *fkx,*fky,*fkz;
  FFT_SCALAR *density_fft;
  FFT_SCALAR *work1,*work2;
  FFT_SCALAR *buf1,*buf2;

  double *gf_b;
  FFT_SCALAR **rho1d,**rho_coeff;

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
  virtual void allocate();
  virtual void deallocate();
  int factorable(int);
  double rms(double, double, bigint, double, double **);
  double diffpr(double, double, double, double, double **);
  void compute_gf_denom();
  double gf_denom(double, double, double);
  virtual void particle_map();
  virtual void make_rho();
  virtual void brick2fft();
  virtual void fillbrick();
  virtual void poisson(int, int);
  virtual void fieldforce();
  void procs2grid2d(int,int,int,int *, int*);
  void compute_rho1d(const FFT_SCALAR &, const FFT_SCALAR &, 
		     const FFT_SCALAR &);
  void compute_rho_coeff();
  void slabcorr(int);
};

}

#endif
#endif
