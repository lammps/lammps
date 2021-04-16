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
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

   Heavily modified and adapted for LAMMPS by the LAMMPS developers.
------------------------------------------------------------------------- */

#ifndef LMP_REAXFF_API_H
#define LMP_REAXFF_API_H

#include "reaxc_types.h"

namespace ReaxFF 
{
  struct API {
    control_params *control;
    reax_system *system;
    output_controls *out_control;
    simulation_data *data;
    storage *workspace;
    reax_list *lists;
  };
}

#endif
