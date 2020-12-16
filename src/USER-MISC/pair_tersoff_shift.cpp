/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Wengen Ouyang (Tel Aviv University)
   and Davide Mandelli (Istituto Italiano di Tecnologia)
   e-mail: w.g.ouyang at gmail dot com
------------------------------------------------------------------------- */

#include "pair_tersoff_shift.h"
#include "error.h"
#include <cstring>
#include "citeme.h"

using namespace LAMMPS_NS;

static const char cite_tersoff_shift[] =
  "@Article{Mandelli2019\n"
  " author = {D. Mandelli, W. Ouyang, M. Urbakh, and O. Hod},\n"
  " title = {he Princess and the Nanoscale Pea: Long-Range Penetration of Surface Distortions into Layered Materials Stacks},\n"
  " journal = {ACS Nano},\n"
  " volume =  13,\n"
  " pages =   {7603}\n"
  " year =    2019,\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairTersoffShift::PairTersoffShift(LAMMPS *lmp) : PairTersoff(lmp)
{
  shift_flag = 1;
  if (lmp->citeme) lmp->citeme->add(cite_tersoff_shift);
}
