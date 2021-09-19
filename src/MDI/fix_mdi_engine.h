/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(mdi/engine, FixMDIEngine);
// clang-format on
#else

#ifndef LMP_FIX_MDI_ENGINE_H
#define LMP_FIX_MDI_ENGINE_H

#include "fix.h"
#include <mdi.h>                // IWYU pragma: export

namespace LAMMPS_NS {

class FixMDIEngine : public Fix {
 public:
  FixMDIEngine(class LAMMPS *, int, char **);
  ~FixMDIEngine();
  int setmask();
  void init();

  int execute_command(const char *command, MDI_Comm driver_socket);
  char *engine_mode(const char *node);

  // receive and update forces

  void min_setup(int);
  void post_integrate();
  void post_force(int);
  void min_pre_force(int);     //@COORDS
  void min_post_force(int);    //@FORCES

  double *add_force;          // stores forces added using +FORCE command
  double potential_energy;    // stores potential energy
  double kinetic_energy;      // stores kinetic energy

  // current command

  char *command;

 protected:
  void exchange_forces();

 private:
  int lmpunits;    // REAL or METAL
  int master, ierr;
  int driver_socket;
  int most_recent_init;    // which MDI init command was most recently received?
                           // 0 - none
                           // 1 - MD
                           // 2 - OPTG
  bool exit_flag;
  bool local_exit_flag;
  char *current_node;
  char *target_node;    // is the code supposed to advance to a particular node?
                        // 0 - none
                        // 1 - @COORDS (before pre-force calculation)
                        // 2 - @PRE-FORCES (before final force calculation)
                        // 3 - @FORCES (before time integration)
                        // -1 - after MD_INIT command
                        // -2 - after MD_INIT command followed by @PRE-FORCES (actually @INIT_OPTG?)

  // command to be executed at the target node

  char *target_command;

  char *id_pe;
  char *id_ke;
  class Irregular *irregular;
  class Minimize *minimizer;
  class Compute *pe;
  class Compute *ke;

  void send_types(Error *);
  void send_labels(Error *);
  void send_masses(Error *);
  void receive_coordinates(Error *);
  void send_coordinates(Error *);
  void send_charges(Error *);
  void send_energy(Error *);
  void send_forces(Error *);
  void send_pe(Error *);
  void send_ke(Error *);
  void receive_forces(Error *, int);
  void send_cell(Error *);
  void receive_cell(Error *);
  void send_celldispl(Error *);
  void receive_celldispl(Error *);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.

E: Potential energy ID for fix mdi does not exist

Self-explanatory.

E: Cannot use MDI command without atom IDs

Self-explanatory.

E: MDI command requires consecutive atom IDs

Self-explanatory.

E: Unable to connect to driver

Self-explanatory.

E: Unable to ... driver

Self-explanatory.

E: Unknown command from driver

The driver sent a command that is not supported by the LAMMPS
interface.  In some cases this might be because a nonsensical
command was sent (i.e. "SCF").  In other cases, the LAMMPS
interface might benefit from being expanded.

*/
