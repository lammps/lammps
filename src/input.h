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

#ifndef INPUT_H
#define INPUT_H

#include "stdio.h"
#include "lammps.h"

class Variable;

class Input : public LAMMPS {
 public:
  int narg;                    // # of command args
  char **arg;                  // parsed args for command
  Variable *variable;          // defined variables

  Input(int, char **);
  ~Input();
  void file();                   // process all input
  void file(char *);             // process an input script
  char *one(char *);             // process a single command
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

  void angle_coeff();          // individual commands
  void angle_style();
  void atom_modify();
  void atom_style();
  void bond_coeff();
  void bond_style();
  void boundary();
  void clear();
  void dielectric();
  void dihedral_coeff();
  void dihedral_style();
  void dimension();
  void dipole();
  void dump();
  void dump_modify();
  void echo();
  void fix();
  void fix_modify();
  void group_command();
  void improper_coeff();
  void improper_style();
  void include();
  void jump();
  void kspace_modify();
  void kspace_style();
  void label();
  void lattice();
  void log();
  void mass();
  void min_modify();
  void min_style();
  void neigh_modify();
  void neighbor_command();
  void newton();
  void next_command();
  void pair_coeff();
  void pair_modify();
  void pair_style();
  void pair_write();
  void print();
  void processors();
  void region();
  void reset_timestep();
  void restart();
  void run_style();
  void special_bonds();
  void temp_modify();
  void temperature();
  void thermo();
  void thermo_modify();
  void thermo_style();
  void timestep();
  void undump();
  void unfix();
  void units();
  void variable_command();
};

#endif
