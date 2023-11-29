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

#ifndef HELPERS_H
#define HELPERS_H

#include <QString>
#include <string>

// duplicate string
extern char *mystrdup(const std::string &text);
extern char *mystrdup(const char *text);
extern char *mystrdup(const QString &text);

// find if executable is in path
extern bool has_exe(const QString &exe);

#endif
// Local Variables:
// c-basic-offset: 4
// End:
