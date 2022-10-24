/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/*  ----------------------------------------------------------------------
   Contributing authors: Christopher Barrett (MSU) barrett@me.msstate.edu
                              Doyl Dickel (MSU) doyl@me.msstate.edu
    ----------------------------------------------------------------------*/
/*
“The research described and the resulting data presented herein, unless
otherwise noted, was funded under PE 0602784A, Project T53 "Military
Engineering Applied Research", Task 002 under Contract No. W56HZV-17-C-0095,
managed by the U.S. Army Combat Capabilities Development Command (CCDC) and
the Engineer Research and Development Center (ERDC).  The work described in
this document was conducted at CAVS, MSU.  Permission was granted by ERDC
to publish this information. Any opinions, findings and conclusions or
recommendations expressed in this material are those of the author(s) and
do not necessarily reflect the views of the United States Army.​”

DISTRIBUTION A. Approved for public release; distribution unlimited. OPSEC#4918
 */

#include "rann_fingerprint.h"
#include "pair_rann.h"

#include <cmath>

using namespace LAMMPS_NS::RANN;

Fingerprint::Fingerprint(PairRANN *_pair)
{
  spin = false;
  screen = false;
  empty = true;
  fullydefined = false;
  n_body_type = 0;
  style = "empty";
  pair = _pair;
}

// Smooth cutoff, goes from 1 to zero over the interval rc-dr to rc.
// Same as MEAM uses. Used by generateradialtable and generatexpcuttable.

double Fingerprint::cutofffunction(double r, double rc, double dr)
{
  double out;
  if (r < (rc - dr))
    out = 1;
  else if (r > rc)
    out = 0;
  else {
    out = 1 - (rc - r) / dr;
    out *= out;
    out *= out;
    out = 1 - out;
    out *= out;
  }
  return out;
}

void Fingerprint::generate_rinvssqrttable()
{
  int buf = 5;
  int m;
  double r1;
  double cutmax = pair->cutmax;
  int res = pair->res;
  rinvsqrttable = new double[res + buf];
  for (m = 0; m < (res + buf); m++) {
    r1 = cutmax * cutmax * (double) (m) / (double) (res);
    rinvsqrttable[m] = 1 / sqrt(r1);
  }
}
