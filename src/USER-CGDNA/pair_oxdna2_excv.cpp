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
/* ----------------------------------------------------------------------
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_oxdna2_excv.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairOxdna2Excv::PairOxdna2Excv(LAMMPS *lmp) : PairOxdnaExcv(lmp)
{

}

/* ---------------------------------------------------------------------- */

PairOxdna2Excv::~PairOxdna2Excv()
{

}

/* ----------------------------------------------------------------------
    compute vector COM-excluded volume interaction sites in oxDNA2
------------------------------------------------------------------------- */
void PairOxdna2Excv::compute_interaction_sites(double e1[3], double e2[3], 
    double /*e3*/[3], double rs[3], double rb[3])
{
  double d_cs_x=-0.34, d_cs_y=+0.3408, d_cb=+0.4;

  rs[0] = d_cs_x*e1[0] + d_cs_y*e2[0];
  rs[1] = d_cs_x*e1[1] + d_cs_y*e2[1];
  rs[2] = d_cs_x*e1[2] + d_cs_y*e2[2];

  rb[0] = d_cb*e1[0];
  rb[1] = d_cb*e1[1];
  rb[2] = d_cb*e1[2];

}
