/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_THERMO_H
#define LMP_THERMO_H

#include "pointers.h"
#include <map>

namespace LAMMPS_NS {

class Thermo : protected Pointers {
  friend class MinCG;              // accesses compute_pe
  friend class DumpNetCDF;         // accesses thermo properties
  friend class DumpNetCDFMPIIO;    // accesses thermo properties
  friend class DumpYAML;           // accesses thermo properties

 public:
  char *style;
  int normflag;    // 0 if do not normalize by atoms, 1 if normalize
  int modified;    // 1 if thermo_modify has been used, else 0
  int lostflag;    // IGNORE,WARN,ERROR
  int lostbond;    // ditto for atoms in bonds

  enum { IGNORE, WARN, ERROR };
  enum { INT, FLOAT, BIGINT };

  Thermo(class LAMMPS *, int, char **);
  ~Thermo() override;
  void init();
  bigint lost_check();
  void modify_params(int, char **);
  void header();
  void footer();
  void compute(int);
  int evaluate_keyword(const std::string &, double *);

 private:
  int nfield, nfield_initial;
  int *vtype;
  std::string line;
  std::vector<std::string> keyword, format, format_column_user, keyword_user;
  std::string format_line_user, format_float_user, format_int_user, format_bigint_user;
  std::map<std::string, int> key2col;

  int normvalue;       // use this for normflag unless natoms = 0
  int normuserflag;    // 0 if user has not set, 1 if has
  int normuser;

  int firststep;
  int lostbefore, warnbefore;
  int flushflag, lineflag;

  double last_tpcpu, last_spcpu;
  double last_time;
  bigint last_step;

  bigint natoms;

  // data used by routines that compute single values
  int ivalue;          // integer value to print
  double dvalue;       // double value to print
  bigint bivalue;      // big integer value to print
  int ifield;          // which field in thermo output is being computed
  int *field2index;    // which compute,fix,variable calcs this field
  int *argindex1;      // indices into compute,fix scalar,vector
  int *argindex2;

  // data for keyword-specific Compute objects
  // index = where they are in computes list
  // id = ID of Compute objects
  // Compute * = ptrs to the Compute objects
  int index_temp, index_press_scalar, index_press_vector, index_pe;
  class Compute *temperature, *pressure, *pe;

  int ncompute;                // # of Compute objects called by thermo
  char **id_compute;           // their IDs
  int *compute_which;          // 0/1/2 if should call scalar,vector,array
  class Compute **computes;    // list of ptrs to the Compute objects

  int nfix;             // # of Fix objects called by thermo
  char **id_fix;        // their IDs
  class Fix **fixes;    // list of ptrs to the Fix objects

  int nvariable;         // # of variables evaluated by thermo
  char **id_variable;    // list of variable names
  int *variables;        // list of Variable indices

  // private methods

  void allocate();
  void deallocate();

  void parse_fields(const std::string &);
  int add_compute(const char *, int);
  int add_fix(const char *);
  int add_variable(const char *);

  void check_temp(const std::string &);
  void check_pe(const std::string &);
  void check_press_scalar(const std::string &);
  void check_press_vector(const std::string &);

  typedef void (Thermo::*FnPtr)();
  void addfield(const char *, FnPtr, int);
  FnPtr *vfunc;    // list of ptrs to functions
  void call_vfunc(int ifield);

  void compute_compute();    // functions that compute a single value
  void compute_fix();        // via calls to  Compute,Fix,Variable classes
  void compute_variable();

  // functions that compute a single value
  // customize a new keyword by adding a method prototype

  void compute_step();
  void compute_elapsed();
  void compute_elapsed_long();
  void compute_dt();
  void compute_time();
  void compute_cpu();
  void compute_tpcpu();
  void compute_spcpu();
  void compute_cpuremain();
  void compute_part();
  void compute_timeremain();

  void compute_atoms();
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

  void compute_enthalpy();
  void compute_ecouple();
  void compute_econserve();

  void compute_vol();
  void compute_density();
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

  void compute_xlat();
  void compute_ylat();
  void compute_zlat();

  void compute_bonds();
  void compute_angles();
  void compute_dihedrals();
  void compute_impropers();

  void compute_pxx();
  void compute_pyy();
  void compute_pzz();
  void compute_pxy();
  void compute_pyz();
  void compute_pxz();

  void compute_fmax();
  void compute_fnorm();

  void compute_nbuild();
  void compute_ndanger();

  void compute_cella();
  void compute_cellb();
  void compute_cellc();
  void compute_cellalpha();
  void compute_cellbeta();
  void compute_cellgamma();
};

}    // namespace LAMMPS_NS

#endif
