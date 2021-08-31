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

#ifdef FIX_CLASS
// clang-format off
FixStyle(ttm/grid,FixTTMGrid);
// clang-format on
#else

#ifndef LMP_FIX_TTM_GRID_H
#define LMP_FIX_TTM_GRID_H

#include "fix_ttm.h"

namespace LAMMPS_NS {

class FixTTMGrid : public FixTTM {
 public:
  FixTTMGrid(class LAMMPS *, int, char **);
  ~FixTTMGrid();
  void post_constructor();
  void init();
  void post_force(int);
  void end_of_step();

  // grid communication

  void pack_forward_grid(int, void *, int, int *);
  void unpack_forward_grid(int, void *, int, int *);
  void pack_reverse_grid(int, void *, int, int *);
  void unpack_reverse_grid(int, void *, int, int *);
  void pack_gather_grid(int, void *);
  void unpack_gather_grid(int, void *, void *, int, int, int, int, int, int);

  void write_restart(FILE *);
  void restart(char *);
  double compute_vector(int);
  double memory_usage();

 private:
  int ngridmine,ngridout;
  int nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in;
  int nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out;
  double delxinv,delyinv,delzinv;
  double skin_original;
  FILE *FPout;

  class GridComm *gc;
  int ngc_buf1,ngc_buf2;
  double *gc_buf1,*gc_buf2;

  void allocate_grid();
  void deallocate_grid();
  void read_electron_temperatures(const std::string &);
  void write_electron_temperatures(const std::string &);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

E: Cannot open fix ttm file %s

The output file for the fix ttm command cannot be opened.  Check that
the path and name are correct.

E: Invalid random number seed in fix ttm command

Random number seed must be > 0.

E: Fix ttm electronic_specific_heat must be > 0.0

Self-explanatory.

E: Fix ttm electronic_density must be > 0.0

Self-explanatory.

E: Fix ttm electronic_thermal_conductivity must be >= 0.0

Self-explanatory.

E: Fix ttm gamma_p must be > 0.0

Self-explanatory.

E: Fix ttm gamma_s must be >= 0.0

Self-explanatory.

E: Fix ttm v_0 must be >= 0.0

Self-explanatory.

E: Fix ttm number of nodes must be > 0

Self-explanatory.

E: Cannot use fix ttm with 2d simulation

This is a current restriction of this fix due to the grid it creates.

E: Cannot use non-periodic boundares with fix ttm

This fix requires a fully periodic simulation box.

E: Cannot use fix ttm with triclinic box

This is a current restriction of this fix due to the grid it creates.

E: Electronic temperature dropped below zero

Something has gone wrong with the fix ttm electron temperature model.

E: Fix ttm electron temperatures must be > 0.0

Self-explanatory.

E: Initial temperatures not all set in fix ttm

Self-explanatory.

W: Too many inner timesteps in fix ttm

Self-explanatory.

*/
