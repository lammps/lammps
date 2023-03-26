/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(rheo/kernel,ComputeRHEOKernel)
// clang-format on
#else

#ifndef LMP_COMPUTE_RHEO_KERNEL_H
#define LMP_COMPUTE_RHEO_KERNEL_H

#include "compute.h"
#include <unordered_set>

namespace LAMMPS_NS {

class ComputeRHEOKernel : public Compute {
 public:
  ComputeRHEOKernel(class LAMMPS *, int, char **);
  ~ComputeRHEOKernel() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double memory_usage() override;
  double calc_w(int,int,double,double,double,double);
  double calc_dw(int,int,double,double,double,double);
  double calc_w_quintic(int,int,double,double,double,double);
  double calc_dw_quintic(int,int,double,double,double,double,double *,double *);

  double dWij[3], dWji[3], Wij, Wji;
  int correction_order;

 private:
  int solid_flag;
  int gsl_error_flag;
  std::unordered_set<tagint> gsl_error_tags;

  int kernel_style, zmin, dim, Mdim, ncor;
  int nmax_old, index_coord;
  double h, hsq, hinv, hsqinv, pre_w, pre_wp;
  double ***C;
  double *C0;

  class NeighList *list;
  class ComputeRHEOInterface *compute_interface;
  class FixRHEO *fix_rheo;

  int check_corrections(int);

  double calc_w_crk0(int,int,double,double,double,double);
  double calc_w_crk1(int,int,double,double,double,double);
  double calc_w_crk2(int,int,double,double,double,double);
  void calc_dw_crk1(int,int,double,double,double,double,double *);
  void calc_dw_crk2(int,int,double,double,double,double,double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
