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

#ifndef LMP_UPDATE_H
#define LMP_UPDATE_H

#include "pointers.h"

#include <map>

namespace LAMMPS_NS {

class Update : protected Pointers {
 public:
  double dt;                     // timestep
  double etol, ftol;             // minimizer tolerances on energy/force
  bigint ntimestep;              // current step (dynamics or min iterations)
  int nsteps;                    // # of steps to run (dynamics or min iter)
  int whichflag;                 // 0 for unset, 1 for dynamics, 2 for min
  double atime;                  // simulation time at atimestep
  bigint atimestep;              // last timestep atime was updated
  bigint firststep, laststep;    // 1st & last step of this run
  bigint beginstep, endstep;     // 1st and last step of multiple runs
  int first_update;              // 0 before initial update, 1 after
  int max_eval;                  // max force evaluations for minimizer
  int restrict_output;           // 1 if output should not write dump/restart
  int setupflag;                 // set when setup() is computing forces
  int multireplica;              // 1 if min across replicas, else 0
  int dt_default;                // 1 if dt is at default value, else 0

  bigint eflag_global, eflag_atom;    // timestep global/peratom eng is tallied on
  bigint vflag_global, vflag_atom;    // ditto for virial

  char *unit_style;

  class Integrate *integrate;
  char *integrate_style;

  class Min *minimize;
  char *minimize_style;

  typedef Integrate *(*IntegrateCreator)(LAMMPS *, int, char **);
  typedef Min *(*MinimizeCreator)(LAMMPS *);

  typedef std::map<std::string, IntegrateCreator> IntegrateCreatorMap;
  typedef std::map<std::string, MinimizeCreator> MinimizeCreatorMap;

  IntegrateCreatorMap *integrate_map;
  MinimizeCreatorMap *minimize_map;

  Update(class LAMMPS *);
  ~Update() override;
  void init();
  void set_units(const char *);
  void create_integrate(int, char **, int);
  void create_minimize(int, char **, int);
  void reset_timestep(int, char **);
  void reset_timestep(bigint, bool);
  void update_time();
  double memory_usage();

 private:
  void new_integrate(char *, int, char **, int, int &);
  void new_minimize(char *, int, char **, int, int &);
};

}    // namespace LAMMPS_NS

#endif
