/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef THERMO_H
#define THERMO_H

#include "lammps.h"

class Temperature;
class Pressure;

class Thermo : public LAMMPS {
  friend class WriteRestart;           // accesses lostflag
  friend class Pressure;               // accesses temperature,tensorflag
  friend class Temper;                 // accesses peflag,potential_energy
  friend class MinCG;                  // accesses compute_pe,potential_energy

 public:
  char *style;

  Thermo(int, char **);
  ~Thermo();
  void init();
  double lost_check();
  void modify_params(int, char **);
  void header();
  void compute(int);
  int compute_value(char *, double *);
  void fix_compute_pe();

 private:
  int me;
  char *line;
  int nfield,nfield_initial;
  char **keyword;
  int *vtype;
  int freeze_group_bit;
  int tensorflag,granflag,dipoleflag;
  double *inertia;

  char *format_int_def,*format_g_def,*format_f_def,*format_multi;
  int format_iflag,format_fflag;
  char *format_int_user,*format_float_user;
  int nformat;
  char **format_user;
  int *format_index;
  char **format;

  int normflag;          // 0 if do not normalize by atoms, 1 if normalize
  int normvalue;         // use this for normflag unless natoms = 0
  int normuserflag;      // 0 if user has not set, 1 if has
  int normuser;

  int firststep;
  int lostflag,lostbefore;
  int flushflag,lineflag;

  Pressure *pressure;
  Temperature *temperature;

  int ivalue,tempflag,pressflag,peflag;
  double dvalue,natoms,tempvalue,potential_energy;

  int ntemp,itemp_print;
  int *tempindices;
  int nfix,ifix_print,nfix_print,nfix_energy;
  double *fixvalues;
  int *fixflags,*fixenergy,*fixprint;
  int nextra;

  int nwindow;
  int npartial_t,ncurrent_t,npartial_p,ncurrent_p;
  int npartial_e,ncurrent_e,npartial_pe,ncurrent_pe;
  double tsum,psum,esum,pesum;
  double *tpast,*ppast,*epast,*pepast;

  typedef void (Thermo::*FnPtr)();
  void parse_fields(char *);
  void addfield(char *, FnPtr, int);
  FnPtr *vfunc;                      // list of ptrs to functions

  // customize by adding a method prototype

  void compute_step();      // functions that each compute one value to print
  void compute_atoms();
  void compute_cpu();
  void compute_temp();
  void compute_press();
  void compute_pe();
  void compute_ke();
  void compute_etotal();
  void compute_evdwl();
  void compute_ecoul();
  void compute_epair();
  void compute_ebond();
  void compute_eangle();
  void compute_edihed();
  void compute_eimp();
  void compute_emol();
  void compute_elong();
  void compute_etail();
  void compute_erot();
  void compute_vol();
  void compute_lx();
  void compute_ly();
  void compute_lz();
  void compute_pxx();
  void compute_pyy();
  void compute_pzz();
  void compute_pxy();
  void compute_pyz();
  void compute_pxz();
  void compute_gke();
  void compute_grot();
  void compute_fix();
  void compute_tave();
  void compute_pave();
  void compute_eave();
  void compute_peave();
  void compute_t_id();
};

#endif
