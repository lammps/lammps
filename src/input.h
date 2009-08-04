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

#ifndef INPUT_H
#define INPUT_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class Input : protected Pointers {
 public:
  int narg;                    // # of command args
  char **arg;                  // parsed args for command
  class Variable *variable;    // defined variables

  Input(class LAMMPS *, int, char **);
  ~Input();
  void file();                   // process all input
  void file(const char *);       // process an input script
  char *one(const char *);       // process a single command
  void substitute(char *, int);  // substitute for variables in a string

 private:
  int me;                      // proc ID
  char *command;               // ptr to current command
  int maxarg;                  // max # of args in arg
  char *line,*copy,*work;      // input line & copy of it
  int echo_screen;             // 0 = no, 1 = yes
  int echo_log;                // 0 = no, 1 = yes
  int nfile,maxfile;           // current # and max # of open input files
  int label_active;            // 0 = no label, 1 = looking for label
  char *labelstr;              // label string being looked for
  int jump_skip;               // 1 if skipping next jump, 0 otherwise

  FILE **infiles;              // list of open input files

  void parse();                // parse an input text line
  int execute_command();       // execute a single command

  void clear();                // input script commands
  void echo();
  void ifthenelse();
  void include();
  void jump();
  void label();
  void log();
  void next_command();
  void print();
  void variable_command();

  void angle_coeff();          // LAMMPS commands
  void angle_style();
  void atom_modify();
  void atom_style();
  void bond_coeff();
  void bond_style();
  void boundary();
  void communicate();
  void compute();
  void compute_modify();
  void dielectric();
  void dihedral_coeff();
  void dihedral_style();
  void dimension();
  void dipole();
  void dump();
  void dump_modify();
  void fix();
  void fix_modify();
  void group_command();
  void improper_coeff();
  void improper_style();
  void kspace_modify();
  void kspace_style();
  void lattice();
  void mass();
  void min_modify();
  void min_style();
  void neigh_modify();
  void neighbor_command();
  void newton();
  void pair_coeff();
  void pair_modify();
  void pair_style();
  void pair_write();
  void processors();
  void region();
  void reset_timestep();
  void restart();
  void run_style();
  void shape();
  void special_bonds();
  void thermo();
  void thermo_modify();
  void thermo_style();
  void timestep();
  void uncompute();
  void undump();
  void unfix();
  void units();
};

}

#endif
