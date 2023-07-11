/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MDI_ENGINE_H
#define LMP_MDI_ENGINE_H

#include "mdi.h"
#include "pointers.h"

namespace LAMMPS_NS {

class MDIEngine : protected Pointers {
 public:
  MDIEngine(class LAMMPS *, int, char **);
  void engine_node(const char *node);

 private:
  int lmpunits;    // REAL or METAL or NATIVE
  int root;        // 1 for proc 0, otherwise 0

  MDI_Comm mdicomm;    // MDI communicator

  // state of MDI engine

  int mode;             // which mode engine is in (DEFAULT,MD,OPTG)
  char *mdicmd;         // current MDI command being processed
  char *node_engine;    // which node engine is at
  char *node_driver;    // which node driver has requested
  bool node_match;      // true if driver and engine node currently match
  bool exit_command;    // true if EXIT command received from driver

  // unit conversion factors

  double lmp2mdi_length, mdi2lmp_length;
  double lmp2mdi_energy, mdi2lmp_energy;
  double lmp2mdi_velocity, mdi2lmp_velocity;
  double lmp2mdi_force, mdi2lmp_force;
  double lmp2mdi_pressure, mdi2lmp_pressure;
  double lmp2mdi_virial, mdi2lmp_virial;

  // flags for data received by engine
  // not acted on until a request to send <ENERGY,<FORCES,<PE,<STRESS

  int flag_natoms, flag_types;
  int flag_cell, flag_cell_displ;
  int flag_charges, flag_coords, flag_velocities;

  int sys_natoms;
  int *sys_types;
  double *sys_charges, *sys_coords, *sys_velocities;
  double sys_cell[9], sys_cell_displ[3];

  int nsteps;           // timesteps for MD
  double etol, ftol;    // 4 minimization params for OPTG
  int niterate, max_eval;

  int nbytes;    // NBYTES command value used for length by other commands

  int actionflag;    // 1 if MD or OPTG just completed, else 0

  int *elements;

  // buffers for MDI comm

  int maxatom;
  double *buf1, *buf1all;
  double *buf3, *buf3all;
  int *ibuf1, *ibuf1all;

  // other classes used by MDI

  char *id_ke, *id_pe, *id_press;    // computes invoked by MDI
  class Compute *ke, *pe, *press;

  class Irregular *irregular;    // irregular comm if new COORDS
                                 // are highly displaced

  // static method for MDI to callback to, when LAMMPS used as plugin engine

  static int execute_command_plugin_wrapper(const char *, MDI_Comm, void *);

  // class methods

  int execute_command(const char *, MDI_Comm);
  void mdi_commands();

  void mdi_md();
  void md();
  void mdi_optg();
  void optg();

  void evaluate();
  void create_system();
  void adjust_box();
  void adjust_charges();
  void adjust_coords();
  void adjust_velocities();

  void receive_cell();
  void receive_cell_displ();
  void receive_charges();
  void receive_coords();
  void receive_elements();
  void receive_natoms();
  void receive_nsteps();
  void receive_tolerance();
  void receive_types();
  void receive_velocities();

  void receive_double3(int);

  void send_cell();
  void send_cell_displ();
  void send_total_energy();
  void send_labels();
  void send_natoms();
  void send_pe();
  void send_ke_elec();
  void send_stress();

  void send_double1(int);
  void send_int1(int);
  void send_double3(int);

  void nbytes_command();
  void single_command();
  void many_commands();
  void infile();
  void send_ke();

  void unit_conversions();
  void reallocate();
  void deallocate();
  void allocate();
};

}    // namespace LAMMPS_NS

#endif
