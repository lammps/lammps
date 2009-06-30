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

#ifndef READ_DATA_H
#define READ_DATA_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class ReadData : protected Pointers {
 public:
  ReadData(class LAMMPS *);
  ~ReadData();
  void command(int, char **);

 private:
  int me;
  char *line,*keyword,*buffer;
  FILE *fp;
  int narg,maxarg,compressed;
  char **arg;

  void open(char *);
  void scan(int &, int &, int &, int &);
  int reallocate(int **, int, int);
  void header(int);
  void parse_keyword(int, int);
  void skip_lines(int);
  void parse_coeffs(int, char *);

  void atoms();
  void velocities();
  void bonds();
  void angles();
  void dihedrals();
  void impropers();

  void mass();
  void shape();
  void dipole();

  void paircoeffs();
  void bondcoeffs();
  void anglecoeffs(int);
  void dihedralcoeffs(int);
  void impropercoeffs(int);
};

}

#endif
