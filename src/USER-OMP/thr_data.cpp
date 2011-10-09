/* -------------------------------------------------------------------------
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
   per-thread data management for LAMMPS
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "thr_data.h"

#include <string.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ThrData::ThrData(int tid) : _tid(tid)
{
  eatom_pair = vatom_pair = NULL;
  eatom_bond = vatom_bond = NULL;
  eatom_angle = vatom_angle = NULL;
  eatom_dihed = vatom_dihed = NULL;
  eatom_imprp = vatom_imprp = NULL;
  eatom_kspce = vatom_kspce = NULL;
  _maxeatom = _maxvatom = 0;
}

/* ---------------------------------------------------------------------- */

ThrData::~ThrData()
{
  if (eatom_pair) delete[] eatom_pair;
  if (vatom_pair) delete[] vatom_pair;
}

/* ---------------------------------------------------------------------- */

void ThrData::clear(int nall)
{
  eng_vdwl=eng_coul=eng_bond=eng_angle=eng_dihed=eng_imprp=eng_kspce=0.0;
  memset(virial_pair,0,6*sizeof(double));
  memset(virial_bond,0,6*sizeof(double));
  memset(virial_angle,0,6*sizeof(double));
  memset(virial_dihed,0,6*sizeof(double));
  memset(virial_imprp,0,6*sizeof(double));
  memset(virial_kspce,0,6*sizeof(double));

  if (_accflags & THR_EATOM) {
    if (_accflags & THR_PAIR)
      memset(eatom_pair,0,nall*sizeof(double));
    if (_accflags & THR_BOND)
      memset(eatom_bond,0,nall*sizeof(double));
    if (_accflags & THR_ANGLE)
      memset(eatom_angle,0,nall*sizeof(double));
    if (_accflags & THR_DIHEDRAL)
      memset(eatom_dihed,0,nall*sizeof(double));
    if (_accflags & THR_IMPROPER)
      memset(eatom_imprp,0,nall*sizeof(double));
    if (_accflags & THR_KSPACE)
      memset(eatom_kspce,0,nall*sizeof(double));
  }

  if (_accflags & THR_VATOM) {
    if (_accflags & THR_PAIR)
      memset(vatom_pair,0,6*nall*sizeof(double));
    if (_accflags & THR_BOND)
      memset(vatom_bond,0,6*nall*sizeof(double));
    if (_accflags & THR_ANGLE)
      memset(vatom_angle,0,6*nall*sizeof(double));
    if (_accflags & THR_DIHEDRAL)
      memset(vatom_dihed,0,6*nall*sizeof(double));
    if (_accflags & THR_IMPROPER)
      memset(vatom_imprp,0,6*nall*sizeof(double));
    if (_accflags & THR_KSPACE)
      memset(vatom_kspce,0,6*nall*sizeof(double));
  }

  _redflags = THR_NONE;
}

/* ---------------------------------------------------------------------- */
#define GrowMe(array,type,flag)				 \
  if (_redflags & flag) {				 \
    if (array ## _ ## type) delete[] array ## _ ## type; \
    array ## _ ## type = new double[_max ## array];	 \
  }

void ThrData::grow_arrays(int nmax)
{
  if (_accflags & THR_EATOM) {
    if (_maxeatom < nmax) {
      _maxeatom = nmax;
      GrowMe(eatom,pair,THR_PAIR);
      GrowMe(eatom,bond,THR_BOND);
      GrowMe(eatom,angle,THR_ANGLE);
      GrowMe(eatom,dihed,THR_DIHEDRAL);
      GrowMe(eatom,imprp,THR_IMPROPER);
      GrowMe(eatom,kspce,THR_KSPACE);
    }
  }
  if (_accflags & THR_VATOM) {
    if (_maxvatom < nmax) {
      _maxvatom = nmax;
      GrowMe(vatom,pair,THR_PAIR);
      GrowMe(vatom,bond,THR_BOND);
      GrowMe(vatom,angle,THR_ANGLE);
      GrowMe(vatom,dihed,THR_DIHEDRAL);
      GrowMe(vatom,imprp,THR_IMPROPER);
      GrowMe(vatom,kspce,THR_KSPACE);
    }
  }
}
#undef GrowMe
/* ---------------------------------------------------------------------- */

double ThrData::memory_usage() 
{
  double bytes = (7 + 6*6) * sizeof(double);
  bytes += 2 * sizeof(double*);
  bytes += 4 * sizeof(int);

  int count = 0;
  if (_redflags & THR_PAIR) ++count;
  if (_redflags & THR_BOND) ++count;
  if (_redflags & THR_ANGLE) ++count;
  if (_redflags & THR_DIHEDRAL) ++count;
  if (_redflags & THR_IMPROPER) ++count;
  if (_redflags & THR_KSPACE) ++count;
  bytes += count * _maxeatom * sizeof(double);
  bytes += count * 6 * _maxvatom * sizeof(double);

  return bytes;
}

