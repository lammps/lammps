/* ----------------------------------------------------------------------
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov

Copyright (2003) Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software.  This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
Contributing authors: Leo Silbert (SNL), Gary Grest (SNL),
		      Jeremy Lechman (SNL), Dan Bolintineanu (SNL), Ishan Srivastava (SNL)
----------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pair_granular.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_neigh_history.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define PI27SQ 266.47931882941264802866    // 27*PI**2
#define THREEROOT3 5.19615242270663202362  // 3*sqrt(3)
#define SIXROOT6 14.69693845669906728801   // 6*sqrt(6)
#define INVROOT6 0.40824829046386307274    // 1/sqrt(6)
#define FOURTHIRDS 1.333333333333333       // 4/3
#define TWOPI 6.28318530717959             // 2*PI

#define EPSILON 1e-10

enum {HOOKE, HERTZ, HERTZ_MATERIAL, DMT, JKR};
enum {VELOCITY, VISCOELASTIC, TSUJI};
enum {TANGENTIAL_NOHISTORY, TANGENTIAL_MINDLIN};
enum {TWIST_NONE, TWIST_NOHISTORY, TWIST_SDS, TWIST_MARSHALL};
enum {ROLL_NONE, ROLL_NOHISTORY, ROLL_SDS};

/* ---------------------------------------------------------------------- */

PairGranular::PairGranular(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  no_virial_fdotr_compute = 1;
  fix_history = NULL;

  single_extra = 9;
  svector = new double[single_extra];

  neighprev = 0;

  nmax = 0;
  mass_rigid = NULL;

  onerad_dynamic = NULL;
  onerad_frozen = NULL;
  maxrad_dynamic = NULL;
  maxrad_frozen = NULL;

  dt = update->dt;

  // set comm size needed by this Pair if used with fix rigid

  comm_forward = 1;

  use_history = 0;
  beyond_contact = 0;
  nondefault_history_transfer = 0;
  tangential_history_index = 0;
  roll_history_index = twist_history_index = 0;

}

/* ---------------------------------------------------------------------- */
PairGranular::~PairGranular()
{
  delete [] svector;
  if (fix_history) modify->delete_fix("NEIGH_HISTORY");

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);

    memory->destroy(normal_coeffs);
    memory->destroy(tangential_coeffs);
    memory->destroy(roll_coeffs);
    memory->destroy(twist_coeffs);

    delete [] onerad_dynamic;
    delete [] onerad_frozen;
    delete [] maxrad_dynamic;
    delete [] maxrad_frozen;
  }
  memory->destroy(mass_rigid);
}

void PairGranular::compute(int eflag, int vflag){
#ifdef TEMPLATED_PAIR_GRANULAR
  if (normal == HOOKE){
    if (damping == VELOCITY){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<0,0,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,0,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,0,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<0,0,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,0,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,0,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<0,0,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,0,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,0,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<0,0,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,0,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,0,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<0,0,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,0,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,0,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<0,0,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,0,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,0,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<0,0,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,0,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,0,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<0,0,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,0,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,0,1,3,2>(eflag, vflag);
        }
      }
    }
    else if (damping == VISCOELASTIC){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<0,1,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,1,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,1,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<0,1,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,1,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,1,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<0,1,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,1,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,1,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<0,1,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,1,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,1,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<0,1,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,1,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,1,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<0,1,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,1,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,1,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<0,1,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,1,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,1,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<0,1,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,1,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,1,1,3,2>(eflag, vflag);
        }
      }
    }
    else if (damping == TSUJI){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<0,2,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,2,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,2,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<0,2,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,2,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,2,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<0,2,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,2,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,2,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<0,2,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,2,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,2,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<0,2,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,2,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,2,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<0,2,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,2,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,2,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<0,2,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,2,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,2,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<0,2,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<0,2,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<0,2,1,3,2>(eflag, vflag);
        }
      }
    }
  }
  else if (normal == HERTZ){
    if (damping == VELOCITY){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<1,0,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,0,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,0,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<1,0,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,0,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,0,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<1,0,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,0,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,0,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<1,0,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,0,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,0,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<1,0,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,0,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,0,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<1,0,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,0,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,0,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<1,0,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,0,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,0,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<1,0,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,0,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,0,1,3,2>(eflag, vflag);
        }
      }
    }
    else if (damping == VISCOELASTIC){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<1,1,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,1,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,1,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<1,1,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,1,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,1,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<1,1,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,1,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,1,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<1,1,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,1,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,1,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<1,1,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,1,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,1,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<1,1,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,1,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,1,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<1,1,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,1,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,1,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<1,1,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,1,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,1,1,3,2>(eflag, vflag);
        }
      }
    }
    else if (damping == TSUJI){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<1,2,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,2,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,2,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<1,2,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,2,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,2,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<1,2,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,2,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,2,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<1,2,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,2,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,2,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<1,2,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,2,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,2,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<1,2,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,2,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,2,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<1,2,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,2,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,2,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<1,2,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<1,2,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<1,2,1,3,2>(eflag, vflag);
        }
      }
    }
  }
  else if (normal == HERTZ_MATERIAL){
    if (damping == VELOCITY){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<2,0,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,0,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,0,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<2,0,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,0,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,0,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<2,0,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,0,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,0,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<2,0,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,0,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,0,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<2,0,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,0,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,0,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<2,0,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,0,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,0,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<2,0,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,0,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,0,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<2,0,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,0,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,0,1,3,2>(eflag, vflag);
        }
      }
    }
    else if (damping == VISCOELASTIC){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<2,1,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,1,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,1,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<2,1,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,1,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,1,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<2,1,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,1,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,1,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<2,1,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,1,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,1,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<2,1,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,1,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,1,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<2,1,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,1,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,1,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<2,1,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,1,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,1,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<2,1,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,1,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,1,1,3,2>(eflag, vflag);
        }
      }
    }
    else if (damping == TSUJI){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<2,2,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,2,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,2,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<2,2,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,2,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,2,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<2,2,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,2,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,2,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<2,2,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,2,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,2,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<2,2,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,2,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,2,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<2,2,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,2,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,2,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<2,2,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,2,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,2,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<2,2,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<2,2,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<2,2,1,3,2>(eflag, vflag);
        }
      }
    }
  }
  else if (normal == DMT){
    if (damping == VELOCITY){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<3,0,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,0,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,0,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<3,0,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,0,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,0,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<3,0,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,0,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,0,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<3,0,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,0,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,0,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<3,0,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,0,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,0,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<3,0,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,0,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,0,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<3,0,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,0,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,0,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<3,0,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,0,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,0,1,3,2>(eflag, vflag);
        }
      }
    }
    else if (damping == VISCOELASTIC){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<3,1,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,1,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,1,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<3,1,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,1,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,1,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<3,1,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,1,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,1,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<3,1,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,1,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,1,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<3,1,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,1,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,1,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<3,1,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,1,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,1,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<3,1,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,1,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,1,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<3,1,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,1,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,1,1,3,2>(eflag, vflag);
        }
      }
    }
    else if (damping == TSUJI){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<3,2,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,2,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,2,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<3,2,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,2,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,2,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<3,2,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,2,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,2,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<3,2,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,2,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,2,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<3,2,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,2,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,2,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<3,2,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,2,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,2,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<3,2,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,2,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,2,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<3,2,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<3,2,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<3,2,1,3,2>(eflag, vflag);
        }
      }
    }
  }
  else if (normal == JKR){
    if (damping == VELOCITY){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<4,0,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,0,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,0,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<4,0,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,0,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,0,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<4,0,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,0,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,0,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<4,0,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,0,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,0,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<4,0,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,0,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,0,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<4,0,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,0,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,0,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<4,0,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,0,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,0,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<4,0,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,0,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,0,1,3,2>(eflag, vflag);
        }
      }
    }
    else if (damping == VISCOELASTIC){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<4,1,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,1,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,1,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<4,1,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,1,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,1,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<4,1,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,1,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,1,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<4,1,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,1,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,1,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<4,1,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,1,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,1,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<4,1,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,1,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,1,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<4,1,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,1,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,1,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<4,1,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,1,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,1,1,3,2>(eflag, vflag);
        }
      }
    }
    else if (damping == TSUJI){
      if (tangential == TANGENTIAL_NOHISTORY){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<4,2,0,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,2,0,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,2,0,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<4,2,0,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,2,0,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,2,0,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<4,2,0,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,2,0,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,2,0,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<4,2,0,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,2,0,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,2,0,3,2>(eflag, vflag);
        }
      }
      else if (tangential == TANGENTIAL_MINDLIN){
        if (twist == TWIST_NONE){
          if (roll == ROLL_NONE)              compute_templated<4,2,1,0,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,2,1,0,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,2,1,0,2>(eflag, vflag);
        }
        else if (twist == TWIST_NOHISTORY){
          if (roll == ROLL_NONE)              compute_templated<4,2,1,1,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,2,1,1,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,2,1,1,2>(eflag, vflag);
        }
        else if (twist == TWIST_SDS){
          if (roll == ROLL_NONE)              compute_templated<4,2,1,2,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,2,1,2,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,2,1,2,2>(eflag, vflag);
        }
        else if (twist == TWIST_MARSHALL){
          if (roll == ROLL_NONE)              compute_templated<4,2,1,3,0>(eflag, vflag);
          else if (roll == ROLL_NOHISTORY)    compute_templated<4,2,1,3,1>(eflag, vflag);
          else if (roll == ROLL_SDS)          compute_templated<4,2,1,3,2>(eflag, vflag);
        }
      }
    }
  }

#else
  compute_untemplated(Tp_normal, Tp_damping, Tp_tangential,
      Tp_roll, Tp_twist,
      eflag, vflag);
#endif
}

#ifdef TEMPLATED_PAIR_GRANULAR
template < int Tp_normal, int Tp_damping, int Tp_tangential,
           int Tp_twist, int Tp_roll >
void PairGranular::compute_templated(int eflag, int vflag)
#else
void PairGranular::compute_untemplated
  (int Tp_normal, int Tp_damping, int Tp_tangential,
   int Tp_twist, int Tp_roll, int eflag, int vflag)
#endif
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,nx,ny,nz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv;
  double Reff, delta, dR, dR2;

  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;

  double knfac, damp_normal;
  double k_tangential, damp_tangential;
  double Fne, Ft, Fdamp, Fntot, Fcrit, Fscrit, Frcrit;
  double fs, fs1, fs2, fs3;

  //For JKR
  double R2, coh, F_pulloff, delta_pulloff, dist_pulloff, a, a2, E;
  double t0, t1, t2, t3, t4, t5, t6;
  double sqrt1, sqrt2, sqrt3, sqrt4;

  double mi,mj,meff,damp,ccel,tor1,tor2,tor3;
  double relrot1,relrot2,relrot3,vrl1,vrl2,vrl3,vrlmag,vrlmaginv;

  //Rolling
  double k_roll, damp_roll;
  double roll1, roll2, roll3, torroll1, torroll2, torroll3;
  double rollmag, rolldotn, scalefac;
  double fr, fr1, fr2, fr3;

  //Twisting
  double k_twist, damp_twist, mu_twist;
  double signtwist, magtwist, magtortwist, Mtcrit;
  double tortwist1, tortwist2, tortwist3;

  double shrmag,rsht;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *history,*allhistory,**firsthistory;

  bool touchflag;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int historyupdate = 1;
  if (update->setupflag) historyupdate = 0;

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0){
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm_pair(this);
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firsthistory = fix_history->firstvalue;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    touch = firsttouch[i];
    allhistory = firsthistory[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++){
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      jtype = type[j];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      E = normal_coeffs[itype][jtype][0];
      Reff = radi*radj/(radi+radj);
      touchflag = false;

      if (Tp_normal == JKR){
        if (touch[jj]){
          R2 = Reff*Reff;
          coh = normal_coeffs[itype][jtype][3];
          a = cbrt(9.0*M_PI*coh*R2/(4*E));
          delta_pulloff = a*a/Reff - 2*sqrt(M_PI*coh*a/E);
          dist_pulloff = radsum-delta_pulloff;
          touchflag = (rsq < dist_pulloff*dist_pulloff);
        }
        else{
          touchflag = (rsq < radsum*radsum);
        }
      }
      else{
        touchflag = (rsq < radsum*radsum);
      }

      if (!touchflag){
        // unset non-touching neighbors
        touch[jj] = 0;
        history = &allhistory[size_history*jj];
        for (int k = 0; k < size_history; k++) history[k] = 0.0;
      }
      else{
        r = sqrt(rsq);
        rinv = 1.0/r;

        nx = delx*rinv;
        ny = dely*rinv;
        nz = delz*rinv;

        // relative translational velocity

        vr1 = v[i][0] - v[j][0];
        vr2 = v[i][1] - v[j][1];
        vr3 = v[i][2] - v[j][2];

        // normal component

        vnnr = vr1*nx + vr2*ny + vr3*nz; //v_R . n
        vn1 = nx*vnnr;
        vn2 = ny*vnnr;
        vn3 = nz*vnnr;

        // meff = effective mass of pair of particles
        // if I or J part of rigid body, use body mass
        // if I or J is frozen, meff is other particle

        mi = rmass[i];
        mj = rmass[j];
        if (fix_rigid) {
          if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
          if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
        }

        meff = mi*mj / (mi+mj);
        if (mask[i] & freeze_group_bit) meff = mj;
        if (mask[j] & freeze_group_bit) meff = mi;

        delta = radsum - r;
        dR = delta*Reff;
        if (Tp_normal == JKR){
          touch[jj] = 1;
          R2=Reff*Reff;
          coh = normal_coeffs[itype][jtype][3];
          dR2 = dR*dR;
          t0 = coh*coh*R2*R2*E;
          t1 = PI27SQ*t0;
          t2 = 8*dR*dR2*E*E*E;
          t3 = 4*dR2*E;
          sqrt1 = MAX(0, t0*(t1+2*t2)); //In case of sqrt(0) < 0 due to precision issues
          t4 = cbrt(t1+t2+THREEROOT3*M_PI*sqrt(sqrt1));
          t5 = t3/t4 + t4/E;
          sqrt2 = MAX(0, 2*dR + t5);
          t6 = sqrt(sqrt2);
          sqrt3 = MAX(0, 4*dR - t5 + SIXROOT6*coh*M_PI*R2/(E*t6));
          a = INVROOT6*(t6 + sqrt(sqrt3));
          a2 = a*a;
          knfac = FOURTHIRDS*E*a;
          Fne = knfac*a2/Reff - TWOPI*a2*sqrt(4*coh*E/(M_PI*a));
        }
        else{
          knfac = E; //Hooke
          Fne = knfac*delta;
          if (Tp_normal != HOOKE)
            a = sqrt(dR);
            Fne *= a;
          if (Tp_normal == DMT)
            Fne -= 4*MY_PI*normal_coeffs[itype][jtype][3]*Reff;
        }

        //Consider restricting Hooke to only have 'velocity' as an option for damping?
        if (Tp_damping == VELOCITY){
          damp_normal = 1;
        }
        else if (Tp_damping == VISCOELASTIC){
          if (Tp_normal == HOOKE) a = sqrt(dR);
          damp_normal = a*meff;
        }
        else if (Tp_damping == TSUJI){
          damp_normal = sqrt(meff*knfac);
        }

        Fdamp = -normal_coeffs[itype][jtype][1]*damp_normal*vnnr;

        Fntot = Fne + Fdamp;

        //****************************************
        //Tangential force, including history effects
        //****************************************

        // tangential component
        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity
        wr1 = (radi*omega[i][0] + radj*omega[j][0]);
        wr2 = (radi*omega[i][1] + radj*omega[j][1]);
        wr3 = (radi*omega[i][2] + radj*omega[j][2]);

        // relative tangential velocities
        vtr1 = vt1 - (nz*wr2-ny*wr3);
        vtr2 = vt2 - (nx*wr3-nz*wr1);
        vtr3 = vt3 - (ny*wr1-nx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // If any history is needed:
        if (use_history){
          touch[jj] = 1;
          history = &allhistory[size_history*jj];
        }


        if (Tp_normal == JKR){
          F_pulloff = 3*M_PI*coh*Reff;
          Fcrit = fabs(Fne + 2*F_pulloff);
        }
        else{
          Fcrit = fabs(Fne);
        }

        //------------------------------
        //Tangential forces
        //------------------------------
        k_tangential = tangential_coeffs[itype][jtype][0];
        damp_tangential = tangential_coeffs[itype][jtype][1]*damp_normal;

        if (Tp_tangential > 0){
          shrmag = sqrt(history[0]*history[0] + history[1]*history[1] +
              history[2]*history[2]);

          // Rotate and update displacements.
          // See e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
          if (historyupdate) {
            rsht = history[0]*nx + history[1]*ny + history[2]*nz;
            if (fabs(rsht) < EPSILON) rsht = 0;
            if (rsht > 0){
              scalefac = shrmag/(shrmag - rsht); //if rhst == shrmag, contacting pair has rotated 90 deg. in one step, in which case you deserve a crash!
              history[0] -= rsht*nx;
              history[1] -= rsht*ny;
              history[2] -= rsht*nz;
              //Also rescale to preserve magnitude
              history[0] *= scalefac;
              history[1] *= scalefac;
              history[2] *= scalefac;
            }
            //Update history
            history[0] += vtr1*dt;
            history[1] += vtr2*dt;
            history[2] += vtr3*dt;
          }

          // tangential forces = history + tangential velocity damping
          fs1 = -k_tangential*history[0] - damp_tangential*vtr1;
          fs2 = -k_tangential*history[1] - damp_tangential*vtr2;
          fs3 = -k_tangential*history[2] - damp_tangential*vtr3;

          // rescale frictional displacements and forces if needed
          Fscrit = tangential_coeffs[itype][jtype][2] * Fcrit;
          fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
          if (fs > Fscrit) {
            if (shrmag != 0.0) {
              history[0] = -1.0/k_tangential*(Fscrit*fs1/fs + damp_tangential*vtr1);
              history[1] = -1.0/k_tangential*(Fscrit*fs2/fs + damp_tangential*vtr2);
              history[2] = -1.0/k_tangential*(Fscrit*fs3/fs + damp_tangential*vtr3);
              fs1 *= Fscrit/fs;
              fs2 *= Fscrit/fs;
              fs3 *= Fscrit/fs;
            } else fs1 = fs2 = fs3 = 0.0;
          }
        }
        else{ //Classic pair gran/hooke (no history)
          fs = meff*damp_tangential*vrel;
          if (vrel != 0.0) Ft = MIN(Fne,fs) / vrel;
          else Ft = 0.0;
          fs1 = -Ft*vtr1;
          fs2 = -Ft*vtr2;
          fs3 = -Ft*vtr3;
        }

        //****************************************
        // Rolling resistance
        //****************************************

        if (Tp_roll != ROLL_NONE){
          relrot1 = omega[i][0] - omega[j][0];
          relrot2 = omega[i][1] - omega[j][1];
          relrot3 = omega[i][2] - omega[j][2];

          // rolling velocity, see eq. 31 of Wang et al, Particuology v 23, p 49 (2015)
          // This is different from the Marshall papers, which use the Bagi/Kuhn formulation
          // for rolling velocity (see Wang et al for why the latter is wrong)
          vrl1 = Reff*(relrot2*nz - relrot3*ny); //- 0.5*((radj-radi)/radsum)*vtr1;
          vrl2 = Reff*(relrot3*nx - relrot1*nz); //- 0.5*((radj-radi)/radsum)*vtr2;
          vrl3 = Reff*(relrot1*ny - relrot2*nx); //- 0.5*((radj-radi)/radsum)*vtr3;
          vrlmag = sqrt(vrl1*vrl1+vrl2*vrl2+vrl3*vrl3);
          if (vrlmag != 0.0) vrlmaginv = 1.0/vrlmag;
          else vrlmaginv = 0.0;

          if (Tp_roll > 1){
            int rhist0 = roll_history_index;
            int rhist1 = rhist0 + 1;
            int rhist2 = rhist1 + 1;

            // Rolling displacement
            rollmag = sqrt(history[rhist0]*history[rhist0] +
                history[rhist1]*history[rhist1] +
                history[rhist2]*history[rhist2]);

            rolldotn = history[rhist0]*nx + history[rhist1]*ny + history[rhist2]*nz;

            if (historyupdate){
              if (fabs(rolldotn) < EPSILON) rolldotn = 0;
              if (rolldotn > 0){ //Rotate into tangential plane
                scalefac = rollmag/(rollmag - rolldotn);
                history[rhist0] -= rolldotn*nx;
                history[rhist1] -= rolldotn*ny;
                history[rhist2] -= rolldotn*nz;
                //Also rescale to preserve magnitude
                history[rhist0] *= scalefac;
                history[rhist1] *= scalefac;
                history[rhist2] *= scalefac;
              }
              history[rhist0] += vrl1*dt;
              history[rhist1] += vrl2*dt;
              history[rhist2] += vrl3*dt;
            }

            k_roll = roll_coeffs[itype][jtype][0];
            damp_roll = roll_coeffs[itype][jtype][1];
            fr1 = -k_roll*history[rhist0] - damp_roll*vrl1;
            fr2 = -k_roll*history[rhist1] - damp_roll*vrl2;
            fr3 = -k_roll*history[rhist2] - damp_roll*vrl3;

            // rescale frictional displacements and forces if needed
            Frcrit = roll_coeffs[itype][jtype][2] * Fcrit;

            fr = sqrt(fr1*fr1 + fr2*fr2 + fr3*fr3);
            if (fr > Frcrit) {
              if (rollmag != 0.0) {
                history[rhist0] = -1.0/k_roll*(Frcrit*fr1/fr + damp_roll*vrl1);
                history[rhist1] = -1.0/k_roll*(Frcrit*fr2/fr + damp_roll*vrl2);
                history[rhist2] = -1.0/k_roll*(Frcrit*fr3/fr + damp_roll*vrl3);
                fr1 *= Frcrit/fr;
                fr2 *= Frcrit/fr;
                fr3 *= Frcrit/fr;
              } else fr1 = fr2 = fr3 = 0.0;
            }
          }
          else{ //
            fr = meff*roll_coeffs[itype][jtype][1]*vrlmag;
            if (vrlmag != 0.0) fr = MIN(Fne, fr) / vrlmag;
            else fr = 0.0;
            fr1 = -fr*vrl1;
            fr2 = -fr*vrl2;
            fr3 = -fr*vrl3;
          }
        }

        //****************************************
        // Twisting torque, including history effects
        //****************************************
        if (Tp_twist != TWIST_NONE){
          magtwist = relrot1*nx + relrot2*ny + relrot3*nz; //Omega_T (eq 29 of Marshall)
          if (Tp_twist == TWIST_MARSHALL){
            k_twist = 0.5*k_tangential*a*a;; //eq 32
            damp_twist = 0.5*damp_tangential*a*a;
            mu_twist = TWOTHIRDS*a;
          }
          else{
            k_twist = twist_coeffs[itype][jtype][0];
            damp_twist = twist_coeffs[itype][jtype][1];
            mu_twist = twist_coeffs[itype][jtype][2];
          }
          if (Tp_twist > 1){
            if (historyupdate){
              history[twist_history_index] += magtwist*dt;
            }
            magtortwist = -k_twist*history[twist_history_index] - damp_twist*magtwist;//M_t torque (eq 30)
            signtwist = (magtwist > 0) - (magtwist < 0);
            Mtcrit = TWOTHIRDS*a*Fscrit;//critical torque (eq 44)
            if (fabs(magtortwist) > Mtcrit) {
              history[twist_history_index] = 1.0/k_twist*(Mtcrit*signtwist - damp_twist*magtwist);
              magtortwist = -Mtcrit * signtwist; //eq 34
            }
          }
          else{
            if (magtwist > 0) magtortwist = -damp_twist*magtwist;
            else magtortwist = 0;
          }
        }
        // Apply forces & torques

        fx = nx*Fntot + fs1;
        fy = ny*Fntot + fs2;
        fz = nz*Fntot + fs3;

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        tor1 = ny*fs3 - nz*fs2;
        tor2 = nz*fs1 - nx*fs3;
        tor3 = nx*fs2 - ny*fs1;

        torque[i][0] -= radi*tor1;
        torque[i][1] -= radi*tor2;
        torque[i][2] -= radi*tor3;

        if (Tp_twist != TWIST_NONE){
          tortwist1 = magtortwist * nx;
          tortwist2 = magtortwist * ny;
          tortwist3 = magtortwist * nz;

          torque[i][0] += tortwist1;
          torque[i][1] += tortwist2;
          torque[i][2] += tortwist3;
        }

        if (Tp_roll != ROLL_NONE){
          torroll1 = Reff*(ny*fr3 - nz*fr2); //n cross fr
          torroll2 = Reff*(nz*fr1 - nx*fr3);
          torroll3 = Reff*(nx*fr2 - ny*fr1);

          torque[i][0] += torroll1;
          torque[i][1] += torroll2;
          torque[i][2] += torroll3;
        }

        if (force->newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;

          torque[j][0] -= radj*tor1;
          torque[j][1] -= radj*tor2;
          torque[j][2] -= radj*tor3;

          if (Tp_twist != TWIST_NONE){
            torque[j][0] -= tortwist1;
            torque[j][1] -= tortwist2;
            torque[j][2] -= tortwist3;
          }
          if (Tp_roll != ROLL_NONE){
            torque[j][0] -= torroll1;
            torque[j][1] -= torroll2;
            torque[j][2] -= torroll3;
          }
        }
        if (evflag) ev_tally_xyz(i,j,nlocal,0,
            0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }
}


/* ----------------------------------------------------------------------
allocate all arrays
------------------------------------------------------------------------- */

void PairGranular::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(normal_coeffs,n+1,n+1,4,"pair:normal_coeffs");
  memory->create(tangential_coeffs,n+1,n+1,3,"pair:tangential_coeffs");
  memory->create(roll_coeffs,n+1,n+1,3,"pair:roll_coeffs");
  memory->create(twist_coeffs,n+1,n+1,3,"pair:twist_coeffs");

  onerad_dynamic = new double[n+1];
  onerad_frozen = new double[n+1];
  maxrad_dynamic = new double[n+1];
  maxrad_frozen = new double[n+1];
}

/* ----------------------------------------------------------------------
  global settings
------------------------------------------------------------------------- */

void PairGranular::settings(int narg, char **arg)
{
  if (narg == 1){
    cutoff_global = force->numeric(FLERR,arg[0]);
  }
  else{
    cutoff_global = -1; //Will be set based on particle sizes, model choice
  }
  tangential_history = 0;
  roll_history = twist_history = 0;
  normal_set = tangential_set = damping_set = roll_set = twist_set = 0;
}

/* ----------------------------------------------------------------------
  set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGranular::coeff(int narg, char **arg)
{
  double normal_coeffs_local[4];
  double tangential_coeffs_local[4];
  double roll_coeffs_local[4];
  double twist_coeffs_local[4];

  if (narg < 2)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  int iarg = 2;
  while (iarg < narg){
    if (strcmp(arg[iarg], "hooke") == 0){
      if (iarg + 2 >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for Hooke option");
      if (!normal_set) normal = HOOKE;
      else if (normal != HOOKE) error->all(FLERR, "Illegal pair_coeff command, choice of normal contact model must be the same for all types");
      normal_coeffs_local[0] = force->numeric(FLERR,arg[iarg+1]); //kn
      normal_coeffs_local[1] = force->numeric(FLERR,arg[iarg+2]); //damping
      normal_set = 1;
      iarg += 3;
    }
    else if (strcmp(arg[iarg], "hertz") == 0){
      int num_coeffs = 2;
      if (iarg + num_coeffs >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for Hertz option");
      if (!normal_set) normal = HERTZ;
      else if (normal_set && normal != HERTZ) if (normal != HOOKE) error->all(FLERR, "Illegal pair_coeff command, choice of normal contact model must be the same for all types");
      normal_coeffs_local[0] = force->numeric(FLERR,arg[iarg+1]); //kn
      normal_coeffs_local[1] = force->numeric(FLERR,arg[iarg+2]); //damping
      normal_set = 1;
      iarg += num_coeffs+1;
    }
    else if (strcmp(arg[iarg], "hertz/material") == 0){
      int num_coeffs = 3;
      if (iarg + num_coeffs >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for Hertz option");
      if (!normal_set) normal = HERTZ_MATERIAL;
      else if (normal != HERTZ) if (normal != HOOKE) error->all(FLERR, "Illegal pair_coeff command, choice of normal contact model must be the same for all types");
      normal_coeffs_local[0] = force->numeric(FLERR,arg[iarg+1])*FOURTHIRDS; //E
      normal_coeffs_local[1] = force->numeric(FLERR,arg[iarg+2]); //damping
      normal_coeffs_local[2] = force->numeric(FLERR,arg[iarg+3]); //G
      normal_set = 1;
      iarg += num_coeffs+1;
    }
    else if (strcmp(arg[iarg], "dmt") == 0){
      if (iarg + 4 >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for Hertz option");
      if (!normal_set) normal = DMT;
      else if (normal != DMT) error->all(FLERR, "Illegal pair_coeff command, choice of normal contact model must be the same for all types");
      normal_coeffs_local[0] = force->numeric(FLERR,arg[iarg+1])*FOURTHIRDS; //E
      normal_coeffs_local[1] = force->numeric(FLERR,arg[iarg+2]); //damping
      normal_coeffs_local[2] = force->numeric(FLERR,arg[iarg+3]); //G
      normal_coeffs_local[3] = force->numeric(FLERR,arg[iarg+3]); //cohesion
      normal_set = 1;
      iarg += 5;
    }
    else if (strcmp(arg[iarg], "jkr") == 0){
      if (iarg + 4 >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for JKR option");
      beyond_contact = 1;
      if (!normal_set) normal = JKR;
      else if (normal != JKR) error->all(FLERR, "Illegal pair_coeff command, choice of normal contact model must be the same for all types");
      normal_coeffs_local[0] = force->numeric(FLERR,arg[iarg+1]); //E
      normal_coeffs_local[1] = force->numeric(FLERR,arg[iarg+2]); //damping
      normal_coeffs_local[2] = force->numeric(FLERR,arg[iarg+3]); //G
      normal_coeffs_local[3] = force->numeric(FLERR,arg[iarg+4]); //cohesion
      normal_set = 1;
      iarg += 5;
    }
    else if (strcmp(arg[iarg], "damp") == 0){
      if (iarg+1 >= narg) error->all(FLERR, "Illegal pair_coeff command, not enough parameters provided for damping model");
      if (strcmp(arg[iarg+1], "velocity") == 0){
        if (!damping_set) damping = VELOCITY;
        else if (damping != VELOCITY) error->all(FLERR, "Illegal pair_coeff command, choice of damping contact model must be the same for all types");
      }
      else if (strcmp(arg[iarg+1], "viscoelastic") == 0){
        if (!damping_set) damping = VISCOELASTIC;
        else if (damping != VISCOELASTIC) error->all(FLERR, "Illegal pair_coeff command, choice of damping contact model must be the same for all types");
      }
      else if (strcmp(arg[iarg+1], "tsuji") == 0){
        if (!damping_set) damping = TSUJI;
        if (damping != TSUJI) error->all(FLERR, "Illegal pair_coeff command, choice of damping contact model must be the same for all types");
      }
      damping_set = 1;
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "tangential") == 0){
      if (iarg + 4 >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for tangential model");
      if (strcmp(arg[iarg+1], "nohistory") == 0){
        if (!tangential_set) tangential = TANGENTIAL_NOHISTORY;
        else if (tangential != TANGENTIAL_NOHISTORY) error->all(FLERR, "Illegal pair_coeff command, choice of tangential contact model must be the same for all types");
      }
      else if (strcmp(arg[iarg+1], "mindlin") == 0){
        if (!tangential_set) tangential = TANGENTIAL_MINDLIN;
        else if (tangential != TANGENTIAL_MINDLIN) error->all(FLERR, "Illegal pair_coeff command, choice of tangential contact model must be the same for all types");;
        tangential_history = 1;
      }
      else{
        error->all(FLERR, "Illegal pair_coeff command, tangential model not recognized");
      }
      tangential_set = 1;
      tangential_coeffs_local[0] = force->numeric(FLERR,arg[iarg+2]); //kt
      tangential_coeffs_local[1] = force->numeric(FLERR,arg[iarg+3]); //gammat
      tangential_coeffs_local[2] = force->numeric(FLERR,arg[iarg+4]); //friction coeff.
      iarg += 5;
    }
    else if (strcmp(arg[iarg], "rolling") == 0){
      if (iarg + 1 >= narg) error->all(FLERR, "Illegal pair_coeff command, not enough parameters");
      if (strcmp(arg[iarg+1], "none") == 0){
        if (!roll_set) roll = ROLL_NONE;
        else if (roll != ROLL_NONE) error->all(FLERR, "Illegal pair_coeff command, choice of rolling friction model must be the same for all types");
        iarg += 2;
      }
      else{
        if (iarg + 4 >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for rolling model");
        if (strcmp(arg[iarg+1], "nohistory") == 0){
          if (!roll_set) roll = ROLL_NOHISTORY;
          else if (roll != ROLL_NOHISTORY) error->all(FLERR, "Illegal pair_coeff command, choice of rolling friction model must be the same for all types");
        }
        else if (strcmp(arg[iarg+1], "sds") == 0){
          if (!roll_set) roll = ROLL_SDS;
          else if (roll != ROLL_SDS) error->all(FLERR, "Illegal pair_coeff command, choice of rolling friction model must be the same for all types");
          roll_history = 1;
        }
        else{
          error->all(FLERR, "Illegal pair_coeff command, rolling friction model not recognized");
        }
        roll_set =1 ;
        roll_coeffs_local[0] = force->numeric(FLERR,arg[iarg+2]); //kR
        roll_coeffs_local[1] = force->numeric(FLERR,arg[iarg+3]); //gammaR
        roll_coeffs_local[2] = force->numeric(FLERR,arg[iarg+4]); //rolling friction coeff.
        iarg += 5;
      }
    }
    else if (strcmp(arg[iarg], "twisting") == 0){
      if (iarg + 1 >= narg) error->all(FLERR, "Illegal pair_coeff command, not enough parameters");
      if (strcmp(arg[iarg+1], "none") == 0){
        if (!twist_set) twist = TWIST_NONE;
        else if (twist != TWIST_NONE) error->all(FLERR, "Illegal pair_coeff command, choice of twisting friction model must be the same for all types");
        iarg += 2;
      }
      else if (strcmp(arg[iarg+1], "marshall") == 0){
        if (!twist_set) twist = TWIST_MARSHALL;
        else if (twist != TWIST_MARSHALL) error->all(FLERR, "Illegal pair_coeff command, choice of twisting friction model must be the same for all types");
        twist_history = 1;
        twist_set = 1;
        iarg += 2;
      }
      else{
        if (iarg + 4 >= narg) error->all(FLERR,"Illegal pair_coeff command, not enough parameters provided for twist model");
        else if (strcmp(arg[iarg+1], "nohistory") == 0){
          if (!twist_set) twist = TWIST_NOHISTORY;
          if (twist != TWIST_NOHISTORY) error->all(FLERR, "Illegal pair_coeff command, choice of twisting friction model must be the same for all types");
        }
        else if (strcmp(arg[iarg+1], "sds") == 0){
          if (!twist_set) twist = TWIST_SDS;
          else if (twist != TWIST_SDS) error->all(FLERR, "Illegal pair_coeff command, choice of twisting friction model must be the same for all types");
          twist_history = 1;
        }
        else{
          error->all(FLERR, "Illegal pair_coeff command, twisting friction model not recognized");
        }
        twist_set = 1;
        twist_coeffs_local[0] = force->numeric(FLERR,arg[iarg+2]); //kt
        twist_coeffs_local[1] = force->numeric(FLERR,arg[iarg+3]); //gammat
        twist_coeffs_local[2] = force->numeric(FLERR,arg[iarg+4]); //friction coeff.
        iarg += 5;
      }
    }
    else error->all(FLERR, "Illegal pair coeff command");
  }

  //It is an error not to specify normal or tangential model
  if (!normal_set) error->all(FLERR, "Illegal pair coeff command, must specify normal contact model");
  if (!tangential_set) error->all(FLERR, "Illegal pair coeff command, must specify tangential contact model");

  //If unspecified, set damping to VISCOELASTIC, twist/roll to NONE (cannot be changed by subsequent pair_coeff commands)
  if (!damping_set) damping = VISCOELASTIC;
  if (!roll_set) roll = ROLL_NONE;
  if (!twist_set) twist = TWIST_NONE;
  damping_set = roll_set = twist_set = 1;

  int count = 0;
  double damp;
  if (damping == TSUJI){
    double cor;
    cor = normal_coeffs_local[1];
    damp = 1.2728-4.2783*cor+11.087*pow(cor,2)-22.348*pow(cor,3)+
        27.467*pow(cor,4)-18.022*pow(cor,5)+
        4.8218*pow(cor,6);
  }
  else damp = normal_coeffs_local[1];

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      normal_coeffs[i][j][0] = normal_coeffs_local[0];
      normal_coeffs[i][j][1] = damp;
      if (normal != HERTZ && normal != HOOKE) normal_coeffs[i][j][2] = normal_coeffs_local[2];
      if ((normal == JKR) || (normal == DMT))
        normal_coeffs[i][j][3] = normal_coeffs_local[3];

      for (int k = 0; k < 3; k++)
        tangential_coeffs[i][j][k] = tangential_coeffs_local[k];

      if (roll != ROLL_NONE)
        for (int k = 0; k < 3; k++)
          roll_coeffs[i][j][k] = roll_coeffs_local[k];

      if (twist != TWIST_NONE && twist != TWIST_MARSHALL)
        for (int k = 0; k < 3; k++)
          twist_coeffs[i][j][k] = twist_coeffs_local[k];

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
	 init specific to this pair style
------------------------------------------------------------------------- */

void PairGranular::init_style()
{
  int i;

  // error and warning checks

  if (!atom->radius_flag || !atom->rmass_flag)
    error->all(FLERR,"Pair granular requires atom attributes radius, rmass");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair granular requires ghost atoms store velocity");

  // Determine whether we need a granular neigh list, how large it needs to be
  use_history = tangential_history || roll_history || twist_history;

  //For JKR, will need fix/neigh/history to keep track of touch arrays
  if (normal == JKR) use_history = 1;

  size_history = 3*tangential_history + 3*roll_history + twist_history;

  //Determine location of tangential/roll/twist histories in array
  if (roll_history){
    if (tangential_history) roll_history_index = 3;
    else roll_history_index = 0;
  }
  if (twist_history){
    if (tangential_history){
      if (roll_history) twist_history_index = 6;
      else twist_history_index = 3;
    }
    else{
      if (roll_history) twist_history_index = 3;
      else twist_history_index = 0;
    }
  }

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->size = 1;
  if (use_history) neighbor->requests[irequest]->history = 1;

  dt = update->dt;

  // if history is stored:
  // if first init, create Fix needed for storing history

  if (use_history && fix_history == NULL) {
    char dnumstr[16];
    sprintf(dnumstr,"%d",size_history);
    char **fixarg = new char*[4];
    fixarg[0] = (char *) "NEIGH_HISTORY";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "NEIGH_HISTORY";
    fixarg[3] = dnumstr;
    modify->add_fix(4,fixarg,1);
    delete [] fixarg;
    fix_history = (FixNeighHistory *) modify->fix[modify->nfix-1];
    fix_history->pair = this;
  }

  // check for FixFreeze and set freeze_group_bit

  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"freeze") == 0) break;
  if (i < modify->nfix) freeze_group_bit = modify->fix[i]->groupbit;
  else freeze_group_bit = 0;

  // check for FixRigid so can extract rigid body masses

  fix_rigid = NULL;
  for (i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) break;
  if (i < modify->nfix) fix_rigid = modify->fix[i];

  // check for FixPour and FixDeposit so can extract particle radii

  int ipour;
  for (ipour = 0; ipour < modify->nfix; ipour++)
    if (strcmp(modify->fix[ipour]->style,"pour") == 0) break;
  if (ipour == modify->nfix) ipour = -1;

  int idep;
  for (idep = 0; idep < modify->nfix; idep++)
    if (strcmp(modify->fix[idep]->style,"deposit") == 0) break;
  if (idep == modify->nfix) idep = -1;

  // set maxrad_dynamic and maxrad_frozen for each type
  // include future FixPour and FixDeposit particles as dynamic

  int itype;
  for (i = 1; i <= atom->ntypes; i++) {
    onerad_dynamic[i] = onerad_frozen[i] = 0.0;
    if (ipour >= 0) {
      itype = i;
      double radmax = *((double *) modify->fix[ipour]->extract("radius",itype));
      if (normal == JKR) radmax = radmax - 0.5*pulloff_distance(radmax, itype);
      onerad_dynamic[i] = radmax;
    }
    if (idep >= 0) {
      itype = i;
      double radmax = *((double *) modify->fix[idep]->extract("radius",itype));
      if (normal == JKR) radmax = radmax - 0.5*pulloff_distance(radmax, itype);
      onerad_dynamic[i] = radmax;
    }
  }

  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++){
    double radius_cut = radius[i];
    if (normal == JKR){
      radius_cut = radius[i] - 0.5*pulloff_distance(radius[i], type[i]);
    }
    if (mask[i] & freeze_group_bit){
      onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]],radius_cut);
    }
    else{
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]],radius_cut);
    }
  }

  MPI_Allreduce(&onerad_dynamic[1],&maxrad_dynamic[1],atom->ntypes,
      MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&onerad_frozen[1],&maxrad_frozen[1],atom->ntypes,
      MPI_DOUBLE,MPI_MAX,world);

  // set fix which stores history info

  if (size_history > 0){
    int ifix = modify->find_fix("NEIGH_HISTORY");
    if (ifix < 0) error->all(FLERR,"Could not find pair fix neigh history ID");
    fix_history = (FixNeighHistory *) modify->fix[ifix];
  }
}

/* ----------------------------------------------------------------------
	 init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGranular::init_one(int i, int j)
{
  double cutoff;
  if (setflag[i][j] == 0) {

    if (normal != HOOKE && normal != HERTZ){
      normal_coeffs[i][j][0] = mix_stiffnessE(normal_coeffs[i][i][0], normal_coeffs[j][j][0],
          normal_coeffs[i][i][2], normal_coeffs[j][j][2]);
      normal_coeffs[i][j][2] = mix_stiffnessG(normal_coeffs[i][i][0], normal_coeffs[j][j][0],
          normal_coeffs[i][i][2], normal_coeffs[j][j][2]);
    }
    else{
      normal_coeffs[i][j][0] = mix_geom(normal_coeffs[i][i][0], normal_coeffs[j][j][0]);
    }

    normal_coeffs[i][j][1] = mix_geom(normal_coeffs[i][i][1], normal_coeffs[j][j][1]);
    if ((normal == JKR) || (normal == DMT))
      normal_coeffs[i][j][3] = mix_geom(normal_coeffs[i][i][3], normal_coeffs[j][j][3]);

    for (int k = 0; k < 3; k++)
      tangential_coeffs[i][j][k] = mix_geom(tangential_coeffs[i][i][k], tangential_coeffs[j][j][k]);


    if (roll != ROLL_NONE){
      for (int k = 0; k < 3; k++)
        roll_coeffs[i][j][k] = mix_geom(roll_coeffs[i][i][k], roll_coeffs[j][j][k]);
    }

    if (twist != TWIST_NONE && twist != TWIST_MARSHALL){
      for (int k = 0; k < 3; k++)
        twist_coeffs[i][j][k] = mix_geom(twist_coeffs[i][i][k], twist_coeffs[j][j][k]);
    }
  }

  // It is possible that cut[i][j] at this point is still 0.0. This can happen when
  // there is a future fix_pour after the current run. A cut[i][j] = 0.0 creates
  // problems because neighbor.cpp uses min(cut[i][j]) to decide on the bin size
  // To avoid this issue, for cases involving  cut[i][j] = 0.0 (possible only
  // if there is no current information about radius/cutoff of type i and j).
  // we assign cutoff = max(cut[i][j]) for i,j such that cut[i][j] > 0.0.

  if (cutoff_global < 0){
    if (((maxrad_dynamic[i] > 0.0) && (maxrad_dynamic[j] > 0.0)) ||
        ((maxrad_dynamic[i] > 0.0) &&  (maxrad_frozen[j] > 0.0)) ||
        ((maxrad_frozen[i] > 0.0)  && (maxrad_dynamic[j] > 0.0))) { // radius info about both i and j exist
      cutoff = maxrad_dynamic[i]+maxrad_dynamic[j];
      cutoff = MAX(cutoff,maxrad_frozen[i]+maxrad_dynamic[j]);
      cutoff = MAX(cutoff,maxrad_dynamic[i]+maxrad_frozen[j]);
    }
    else { // radius info about either i or j does not exist (i.e. not present and not about to get poured; set to largest value to not interfere with neighbor list)
      double cutmax = 0.0;
      for (int k = 1; k <= atom->ntypes; k++) {
        cutmax = MAX(cutmax,2.0*maxrad_dynamic[k]);
        cutmax = MAX(cutmax,2.0*maxrad_frozen[k]);
      }
      cutoff = cutmax;
    }
  }
  else{
    cutoff = cutoff_global;
  }
  return cutoff;
}


/* ----------------------------------------------------------------------
	 proc 0 writes to restart file
	 ------------------------------------------------------------------------- */

void PairGranular::write_restart(FILE *fp)
{
  int i,j;
  fwrite(&normal,sizeof(int),1,fp);
  fwrite(&damping,sizeof(int),1,fp);
  fwrite(&tangential,sizeof(int),1,fp);
  fwrite(&roll,sizeof(int),1,fp);
  fwrite(&twist,sizeof(int),1,fp);
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&normal_coeffs[i][j],sizeof(double),4,fp);
        fwrite(&tangential_coeffs[i][j],sizeof(double),3,fp);
        fwrite(&roll_coeffs[i][j],sizeof(double),3,fp);
        fwrite(&twist_coeffs[i][j],sizeof(double),3,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
	 proc 0 reads from restart file, bcasts
	 ------------------------------------------------------------------------- */

void PairGranular::read_restart(FILE *fp)
{
  allocate();
  int i,j;
  int me = comm->me;
  if (me == 0){
    fread(&normal,sizeof(int),1,fp);
    fread(&damping,sizeof(int),1,fp);
    fread(&tangential,sizeof(int),1,fp);
    fread(&roll,sizeof(int),1,fp);
    fread(&twist,sizeof(int),1,fp);
  }
  MPI_Bcast(&normal,1,MPI_INT,0,world);
  MPI_Bcast(&damping,1,MPI_INT,0,world);
  MPI_Bcast(&tangential,1,MPI_INT,0,world);
  MPI_Bcast(&roll,1,MPI_INT,0,world);
  MPI_Bcast(&twist,1,MPI_INT,0,world);
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&normal_coeffs[i][j],sizeof(double),4,fp);
          fread(&tangential_coeffs[i][j],sizeof(double),3,fp);
          fread(&roll_coeffs[i][j],sizeof(double),3,fp);
          fread(&twist_coeffs[i][j],sizeof(double),3,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&normal_coeffs[i][j],4,MPI_DOUBLE,0,world);
        MPI_Bcast(&tangential_coeffs[i][j],3,MPI_DOUBLE,0,world);
        MPI_Bcast(&roll_coeffs[i][j],3,MPI_DOUBLE,0,world);
        MPI_Bcast(&twist_coeffs[i][j],3,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}


/* ---------------------------------------------------------------------- */

void PairGranular::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

double PairGranular::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce)
{
  double radi,radj,radsum;
  double r,rinv,rsqinv,delx,dely,delz, nx, ny, nz, Reff;
  double dR, dR2;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff,damp,ccel,tor1,tor2,tor3;
  double relrot1,relrot2,relrot3,vrl1,vrl2,vrl3,vrlmag,vrlmaginv;

  double knfac, damp_normal;
  double k_tangential, damp_tangential;
  double Fne, Ft, Fdamp, Fntot, Fcrit, Fscrit, Frcrit;
  double fs, fs1, fs2, fs3;

  //For JKR
  double R2, coh, F_pulloff, delta_pulloff, dist_pulloff, a, a2, E;
  double delta, t0, t1, t2, t3, t4, t5, t6;
  double sqrt1, sqrt2, sqrt3, sqrt4;


  //Rolling
  double k_roll, damp_roll;
  double roll1, roll2, roll3, torroll1, torroll2, torroll3;
  double rollmag, rolldotn, scalefac;
  double fr, fr1, fr2, fr3;

  //Twisting
  double k_twist, damp_twist, mu_twist;
  double signtwist, magtwist, magtortwist, Mtcrit;
  double tortwist1, tortwist2, tortwist3;

  double shrmag,rsht;
  int jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *history,*allhistory,**firsthistory;

  double *radius = atom->radius;
  radi = radius[i];
  radj = radius[j];
  radsum = radi + radj;
  Reff = radi*radj/(radi+radj);

  bool touchflag;
  if (normal == JKR){
    R2 = Reff*Reff;
    coh = normal_coeffs[itype][jtype][3];
    a = cbrt(9.0*M_PI*coh*R2/(4*E));
    delta_pulloff = a*a/Reff - 2*sqrt(M_PI*coh*a/E);
    dist_pulloff = radsum+delta_pulloff;
    touchflag = (rsq <= dist_pulloff*dist_pulloff);
  }
  else{
    touchflag = (rsq <= radsum*radsum);
  }

  if (touchflag){
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }

  double **x = atom->x;
  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];
  r = sqrt(rsq);
  rinv = 1.0/r;

  nx = delx*rinv;
  ny = dely*rinv;
  nz = delz*rinv;

  // relative translational velocity

  double **v = atom->v;
  vr1 = v[i][0] - v[j][0];
  vr2 = v[i][1] - v[j][1];
  vr3 = v[i][2] - v[j][2];

  // normal component

  vnnr = vr1*nx + vr2*ny + vr3*nz;
  vn1 = nx*vnnr;
  vn2 = ny*vnnr;
  vn3 = nz*vnnr;

  double *rmass = atom->rmass;
  int *mask = atom->mask;
  mi = rmass[i];
  mj = rmass[j];
  if (fix_rigid) {
    if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
    if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
  }

  meff = mi*mj / (mi+mj);
  if (mask[i] & freeze_group_bit) meff = mj;
  if (mask[j] & freeze_group_bit) meff = mi;

  delta = radsum - r;
  dR = delta*Reff;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  double **omega = atom->omega;
  wr1 = (radi*omega[i][0] + radj*omega[j][0]);
  wr2 = (radi*omega[i][1] + radj*omega[j][1]);
  wr3 = (radi*omega[i][2] + radj*omega[j][2]);

  // meff = effective mass of pair of particles
  // if I or J part of rigid body, use body mass
  // if I or J is frozen, meff is other particle

  int *type = atom->type;

  mi = rmass[i];
  mj = rmass[j];
  if (fix_rigid) {
    // NOTE: ensure mass_rigid is current for owned+ghost atoms?
    if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
    if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
  }

  meff = mi*mj / (mi+mj);
  if (mask[i] & freeze_group_bit) meff = mj;
  if (mask[j] & freeze_group_bit) meff = mi;

  delta = radsum - r;
  dR = delta*Reff;
  if (normal == JKR){
    dR2 = dR*dR;
    t0 = coh*coh*R2*R2*E;
    t1 = PI27SQ*t0;
    t2 = 8*dR*dR2*E*E*E;
    t3 = 4*dR2*E;
    sqrt1 = MAX(0, t0*(t1+2*t2)); //In case of sqrt(0) < 0 due to precision issues
    t4 = cbrt(t1+t2+THREEROOT3*M_PI*sqrt(sqrt1));
    t5 = t3/t4 + t4/E;
    sqrt2 = MAX(0, 2*dR + t5);
    t6 = sqrt(sqrt2);
    sqrt3 = MAX(0, 4*dR - t5 + SIXROOT6*coh*M_PI*R2/(E*t6));
    a = INVROOT6*(t6 + sqrt(sqrt3));
    a2 = a*a;
    knfac = FOURTHIRDS*E*a;
    Fne = knfac*a2/Reff - TWOPI*a2*sqrt(4*coh*E/(M_PI*a));
  }
  else{
    knfac = E;
    Fne = knfac*delta;
    if (normal != HOOKE)
      a = sqrt(dR);
      Fne *= a;
    if (normal == DMT)
      Fne -= 4*MY_PI*normal_coeffs[itype][jtype][3]*Reff;
  }

  //Consider restricting Hooke to only have 'velocity' as an option for damping?
  if (damping == VELOCITY){
    damp_normal = normal_coeffs[itype][jtype][1];
  }
  else if (damping == VISCOELASTIC){
    if (normal == HOOKE) a = sqrt(dR);
    damp_normal = normal_coeffs[itype][jtype][1]*a*meff;
  }
  else if (damping == TSUJI){
    damp_normal = normal_coeffs[itype][jtype][1]*sqrt(meff*knfac);
  }

  Fdamp = -damp_normal*vnnr;

  Fntot = Fne + Fdamp;

  jnum = list->numneigh[i];
  jlist = list->firstneigh[i];

  if (use_history){
    allhistory = fix_history->firstvalue[i];
    for (int jj = 0; jj < jnum; jj++) {
      neighprev++;
      if (neighprev >= jnum) neighprev = 0;
      if (jlist[neighprev] == j) break;
    }
    history = &allhistory[size_history*neighprev];
  }

  //****************************************
  //Tangential force, including history effects
  //****************************************

  // tangential component
  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity
  wr1 = (radi*omega[i][0] + radj*omega[j][0]);
  wr2 = (radi*omega[i][1] + radj*omega[j][1]);
  wr3 = (radi*omega[i][2] + radj*omega[j][2]);

  // relative tangential velocities
  vtr1 = vt1 - (nz*wr2-ny*wr3);
  vtr2 = vt2 - (nx*wr3-nz*wr1);
  vtr3 = vt3 - (ny*wr1-nx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  Fcrit = fabs(Fne);
  if (normal == JKR){
    F_pulloff = 3*M_PI*coh*Reff;
    Fcrit = fabs(Fne + 2*F_pulloff);
  }

  //------------------------------
  //Tangential forces
  //------------------------------
  k_tangential = tangential_coeffs[itype][jtype][0];
  if (normal != HOOKE){
    k_tangential *= a;
  }
  damp_tangential = tangential_coeffs[itype][jtype][1]*damp_normal;

  if (tangential_history){
    shrmag = sqrt(history[0]*history[0] + history[1]*history[1] +
        history[2]*history[2]);

    // tangential forces = history + tangential velocity damping
    fs1 = -k_tangential*history[0] - damp_tangential*vtr1;
    fs2 = -k_tangential*history[1] - damp_tangential*vtr2;
    fs3 = -k_tangential*history[2] - damp_tangential*vtr3;

    // rescale frictional displacements and forces if needed
    Fscrit = tangential_coeffs[itype][jtype][2] * Fcrit;
    fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
    if (fs > Fscrit) {
      if (shrmag != 0.0) {
        history[0] = -1.0/k_tangential*(Fscrit*fs1/fs + damp_tangential*vtr1);
        history[1] = -1.0/k_tangential*(Fscrit*fs2/fs + damp_tangential*vtr2);
        history[2] = -1.0/k_tangential*(Fscrit*fs3/fs + damp_tangential*vtr3);
        fs1 *= Fscrit/fs;
        fs2 *= Fscrit/fs;
        fs3 *= Fscrit/fs;
      } else fs1 = fs2 = fs3 = 0.0;
    }
  }
  else{ //Classic pair gran/hooke (no history)
    fs = meff*damp_tangential*vrel;
    if (vrel != 0.0) Ft = MIN(Fne,fs) / vrel;
    else Ft = 0.0;
    fs1 = -Ft*vtr1;
    fs2 = -Ft*vtr2;
    fs3 = -Ft*vtr3;
  }

  //****************************************
  // Rolling resistance
  //****************************************

  if (roll != ROLL_NONE){
    relrot1 = omega[i][0] - omega[j][0];
    relrot2 = omega[i][1] - omega[j][1];
    relrot3 = omega[i][2] - omega[j][2];

    // rolling velocity, see eq. 31 of Wang et al, Particuology v 23, p 49 (2015)
    // This is different from the Marshall papers, which use the Bagi/Kuhn formulation
    // for rolling velocity (see Wang et al for why the latter is wrong)
    vrl1 = Reff*(relrot2*nz - relrot3*ny); //- 0.5*((radj-radi)/radsum)*vtr1;
    vrl2 = Reff*(relrot3*nx - relrot1*nz); //- 0.5*((radj-radi)/radsum)*vtr2;
    vrl3 = Reff*(relrot1*ny - relrot2*nx); //- 0.5*((radj-radi)/radsum)*vtr3;
    vrlmag = sqrt(vrl1*vrl1+vrl2*vrl2+vrl3*vrl3);
    if (vrlmag != 0.0) vrlmaginv = 1.0/vrlmag;
    else vrlmaginv = 0.0;

    if (roll_history){
      int rhist0 = roll_history_index;
      int rhist1 = rhist0 + 1;
      int rhist2 = rhist1 + 1;

      // Rolling displacement
      rollmag = sqrt(history[rhist0]*history[rhist0] +
          history[rhist1]*history[rhist1] +
          history[rhist2]*history[rhist2]);

      rolldotn = history[rhist0]*nx + history[rhist1]*ny + history[rhist2]*nz;

      k_roll = roll_coeffs[itype][jtype][0];
      damp_roll = roll_coeffs[itype][jtype][1];
      fr1 = -k_roll*history[rhist0] - damp_roll*vrl1;
      fr2 = -k_roll*history[rhist1] - damp_roll*vrl2;
      fr3 = -k_roll*history[rhist2] - damp_roll*vrl3;

      // rescale frictional displacements and forces if needed
      Frcrit = roll_coeffs[itype][jtype][2] * Fcrit;

      fr = sqrt(fr1*fr1 + fr2*fr2 + fr3*fr3);
      if (fr > Frcrit) {
        if (rollmag != 0.0) {
          history[rhist0] = -1.0/k_roll*(Frcrit*fr1/fr + damp_roll*vrl1);
          history[rhist1] = -1.0/k_roll*(Frcrit*fr2/fr + damp_roll*vrl2);
          history[rhist2] = -1.0/k_roll*(Frcrit*fr3/fr + damp_roll*vrl3);
          fr1 *= Frcrit/fr;
          fr2 *= Frcrit/fr;
          fr3 *= Frcrit/fr;
        } else fr1 = fr2 = fr3 = 0.0;
      }
    }
    else{ //
      fr = meff*roll_coeffs[itype][jtype][1]*vrlmag;
      if (vrlmag != 0.0) fr = MIN(Fne, fr) / vrlmag;
      else fr = 0.0;
      fr1 = -fr*vrl1;
      fr2 = -fr*vrl2;
      fr3 = -fr*vrl3;
    }
  }

  //****************************************
  // Twisting torque, including history effects
  //****************************************
  if (twist != TWIST_NONE){
    magtwist = relrot1*nx + relrot2*ny + relrot3*nz; //Omega_T (eq 29 of Marshall)
    if (twist == TWIST_MARSHALL){
      k_twist = 0.5*k_tangential*a*a;; //eq 32
      damp_twist = 0.5*damp_tangential*a*a;
      mu_twist = TWOTHIRDS*a;
    }
    else{
      k_twist = twist_coeffs[itype][jtype][0];
      damp_twist = twist_coeffs[itype][jtype][1];
      mu_twist = twist_coeffs[itype][jtype][2];
    }
    if (twist_history){
      magtortwist = -k_twist*history[twist_history_index] - damp_twist*magtwist;//M_t torque (eq 30)
      signtwist = (magtwist > 0) - (magtwist < 0);
      Mtcrit = TWOTHIRDS*a*Fscrit;//critical torque (eq 44)
      if (fabs(magtortwist) > Mtcrit) {
        history[twist_history_index] = 1.0/k_twist*(Mtcrit*signtwist - damp_twist*magtwist);
        magtortwist = -Mtcrit * signtwist; //eq 34
      }
    }
    else{
      if (magtwist > 0) magtortwist = -damp_twist*magtwist;
      else magtortwist = 0;
    }
  }

  // set single_extra quantities

  svector[0] = fs1;
  svector[1] = fs2;
  svector[2] = fs3;
  svector[3] = fs;
  svector[4] = fr1;
  svector[5] = fr2;
  svector[6] = fr3;
  svector[7] = fr;
  svector[8] = magtortwist;
  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairGranular::pack_forward_comm(int n, int *list, double *buf,
    int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = mass_rigid[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairGranular::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    mass_rigid[i] = buf[m++];
}

/* ----------------------------------------------------------------------
	 memory usage of local atom-based arrays
	 ------------------------------------------------------------------------- */

double PairGranular::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
 mixing of Young's modulus (E)
------------------------------------------------------------------------- */

double PairGranular::mix_stiffnessE(double Eii, double Ejj, double Gii, double Gjj)
{
  double poisii = Eii/(2.0*Gii) - 1.0;
  double poisjj = Ejj/(2.0*Gjj) - 1.0;
  return 1/((1-poisii*poisjj)/Eii+(1-poisjj*poisjj)/Ejj);
}

/* ----------------------------------------------------------------------
	 mixing of shear modulus (G)
	 ------------------------------------------------------------------------- */

double PairGranular::mix_stiffnessG(double Eii, double Ejj, double Gii, double Gjj)
{
  double poisii = Eii/(2.0*Gii) - 1.0;
  double poisjj = Ejj/(2.0*Gjj) - 1.0;
  return 1/((2.0 -poisjj)/Gii+(2.0-poisjj)/Gjj);
}

/* ----------------------------------------------------------------------
	 mixing of everything else 
------------------------------------------------------------------------- */

double PairGranular::mix_geom(double valii, double valjj)
{
  return sqrt(valii*valjj);
}


/* ----------------------------------------------------------------------
     Compute pull-off distance (beyond contact) for a given radius and atom type
------------------------------------------------------------------------- */

double PairGranular::pulloff_distance(double radius, int itype)
{
  double E, coh, a, delta_pulloff;
  coh = normal_coeffs[itype][itype][3];
  E = mix_stiffnessE(normal_coeffs[itype][itype][0], normal_coeffs[itype][itype][0],
      normal_coeffs[itype][itype][2], normal_coeffs[itype][itype][2]);
  a = cbrt(9*M_PI*coh*radius*radius/(4*E));
  return a*a/radius - 2*sqrt(M_PI*coh*a/E);
}

