/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple)
------------------------------------------------------------------------- */

#include <cstdio>

#include "molfile_interface.h"

using namespace LAMMPS_NS;

int main(int, char **)
{
  const char pdir[] = ".";
  const char ptype[] = "xyz";

  MolfileInterface mif;
  
  int rv = mif.find_plugin(pdir,ptype,MolfileInterface::M_WRITE);
  printf("plugin finder returns: %d\n", rv);
  printf("plugin now: %s\n",mif.get_plugin_name());
  
  rv = mif.find_plugin(pdir,ptype,MolfileInterface::M_WRITE);
  printf("plugin finder returns: %d\n", rv);
  if (rv == MolfileInterface::E_MATCH)
    printf("plugin now: %s\n",mif.get_plugin_name());
  else
    printf("plugin still: %s\n",mif.get_plugin_name());

  mif.forget_plugin();
  rv = mif.find_plugin(pdir,ptype,MolfileInterface::M_WRITE);
  printf("plugin finder returns: %d\n", rv);
  printf("plugin now: %s\n",mif.get_plugin_name());
  
  rv = mif.find_plugin(pdir,ptype,MolfileInterface::M_READ);
  printf("plugin finder returns: %d\n", rv);
  printf("plugin now: %s\n",mif.get_plugin_name());
  
  return 0;
  
}

