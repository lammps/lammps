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

#ifndef THERMO_H
#define THERMO_H

#include "pointers.h"

namespace LAMMPS_NS {

class Thermo : protected Pointers {
  friend class WriteRestart;           // accesses lostflag
  friend class MinCG;                  // accesses compute_pe

 public:
  char *style;
  int normflag;          // 0 if do not normalize by atoms, 1 if normalize
  int modified;          // 1 if thermo_modify has been used, else 0

  Thermo(class LAMMPS *, int, char **);
  ~Thermo();
  void init();
  double lost_check();
  void modify_params(int, char **);
  void header();
  void compute(int);
  int evaluate_keyword(char *, double *);

 private:
  int me;

  char *line;
  int nfield,nfield_initial;
  char **keyword;
  int *vtype;

  char *format_multi,*format_int_one_def,*format_int_multi_def;
  char *format_g_def,*format_f_def;
  char *format_int_user,*format_float_user;
  char **format,**format_user;

  int normvalue;         // use this for normflag unless natoms = 0
  int normuserflag;      // 0 if user has not set, 1 if has
  int normuser;

  int firststep;
  int lostflag,lostbefore;
  int flushflag,lineflag;

                         // data used by routines that compute single values
  int ivalue;            // integer value to print
  double dvalue,natoms;  // dvalue = double value to print
  int ifield;            // which field in thermo output is being computed
  int *field2index;      // which compute,fix,variable calcs this field
  int *argindex;         // index into compute,fix scalar,vector

                         // data for keyword-specific Compute objects
                         // index = where they are in computes list
                         // id = ID of Compute objects
                         // Compute * = ptrs to the Compute objects
  int index_temp,index_press_scalar,index_press_vector,index_pe;
  char *id_temp,*id_press,*id_pe;
  class Compute *temperature,*pressure,*pe;

  int ncompute;                // # of Compute objects called by thermo
  char **id_compute;           // their IDs
  int *compute_which;          // 0/1 if should call scalar() or vector()
  class Compute **computes;    // list of ptrs to the Compute objects

  int nfix;                    // # of Fix objects called by thermo
  char **id_fix;               // their IDs
  class Fix **fixes;           // list of ptrs to the Fix objects

  int nvariable;               // # of variables evaulated by thermo
  char **id_variable;          // list of variable names
  int *variables;              // list of Variable indices

  // private methods

  void allocate();
  void deallocate();

  void parse_fields(char *);
  int add_compute(const char *, int);
  int add_fix(const char *);
  int add_variable(const char *);

  typedef void (Thermo::*FnPtr)();
  void addfield(const char *, FnPtr, int);
  FnPtr *vfunc;                // list of ptrs to functions

  void compute_compute();      // functions that compute a single value
  void compute_fix();          // via calls to  Compute,Fix,Variable classes
  void compute_variable();

  // functions that compute a single value
  // customize a new keyword by adding a method prototype

  void compute_step();    
  void compute_atoms();
  void compute_cpu();

  void compute_temp();
  void compute_press();
  void compute_pe();
  void compute_ke();
  void compute_etotal();
  void compute_enthalpy();

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

  void compute_vol();
  void compute_lx();
  void compute_ly();
  void compute_lz();

  void compute_xlo();
  void compute_xhi();
  void compute_ylo();
  void compute_yhi();
  void compute_zlo();
  void compute_zhi();

  void compute_xy();
  void compute_xz();
  void compute_yz();

  void compute_pxx();
  void compute_pyy();
  void compute_pzz();
  void compute_pxy();
  void compute_pyz();
  void compute_pxz();
};

}

#endif
