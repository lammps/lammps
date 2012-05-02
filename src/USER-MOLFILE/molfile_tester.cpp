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

#include <mpi.h>
#include "lammps.h"
#include "input.h"
#include "dump_molfile.h"
#include "molfile_interface.h"

int main(int, char **)
{
  const char pdir[] = ".";
  const char ptype[] = "xyz";

  LAMMPS_NS::MolfileInterface mif;
  
  int rv = mif.find_plugin(pdir,"trr",LAMMPS_NS::MolfileInterface::M_WRITE);
  printf("plugin finder returns: %d\n", rv);
  printf("plugin now: %s\n",mif.get_plugin_name());
  
  rv = mif.find_plugin(pdir,"g96",LAMMPS_NS::MolfileInterface::M_WRITE);
  printf("plugin finder returns: %d\n", rv);
  if (rv == LAMMPS_NS::MolfileInterface::E_MATCH)
    printf("plugin now: %s\n",mif.get_plugin_name());
  else
    printf("plugin still: %s\n",mif.get_plugin_name());

  mif.forget_plugin();
  rv = mif.find_plugin(pdir,ptype,LAMMPS_NS::MolfileInterface::M_WRITE);
  printf("plugin finder returns: %d\n", rv);
  printf("plugin now: %s\n",mif.get_plugin_name());
  
  rv = mif.find_plugin(pdir,ptype,LAMMPS_NS::MolfileInterface::M_READ);
  printf("plugin finder returns: %d\n", rv);
  printf("plugin now: %s\n",mif.get_plugin_name());

  LAMMPS_NS::LAMMPS *lmp;
  char *argv[] = {"lammps", "-log", "none", "-echo", "screen", NULL };
  
  lmp = new LAMMPS_NS::LAMMPS(5,argv,MPI_COMM_WORLD);
  lmp->input->file("in.test");
  lmp->input->one("dump xy1 all xyz 50 melt-native1.xyz");
  lmp->input->one("dump xy2 all xyz 50 melt-native2.xyz");
  lmp->input->one("dump_modify xy2 element C");
  
  char *molf1[] = {"mf1", "all", "molfile", "10", "melt1.dcd", "dcd", ".:.:.", NULL };
  LAMMPS_NS::DumpMolfile *dump1 = new LAMMPS_NS::DumpMolfile(lmp,8,molf1);
  dump1->init();
  lmp->input->one("run 10 pre no post no");
  dump1->write();
  lmp->input->one("run 10 pre no post no");
  dump1->write();
  lmp->input->one("run 10 pre no post no");
  dump1->write();
  lmp->input->one("run 10 pre no post no");
  dump1->write();
  lmp->input->one("run 10 pre no post no");
    
  char *molf2[] = {"mf2", "all", "molfile", "10", "melt2-*.dcd", "dcd", ".", NULL };
  LAMMPS_NS::DumpMolfile *dump2 = new LAMMPS_NS::DumpMolfile(lmp,8,molf2);
  dump2->write();
  lmp->input->one("run 10 pre no post no");
  dump2->write();
  lmp->input->one("run 10 pre no post no");
  dump2->write();
  lmp->input->one("run 10 pre no post no");
  dump2->write();
  lmp->input->one("run 10 pre no post no");
  dump2->write();
  lmp->input->one("run 10 pre no post no");
  dump2->write();
  
  char *molf3[] = {"mf3", "all", "molfile", "50", "melt3.xyz", "xyz", NULL };
  LAMMPS_NS::DumpMolfile *dump3 = new LAMMPS_NS::DumpMolfile(lmp,6,molf3);
  

  delete dump1;
  delete dump2;
  delete dump3;
  delete lmp;
  
  return 0;
  
}

