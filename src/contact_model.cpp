// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
-------------------------------------------------------------------------

   This class contains a series of DEM contact models
   that can be defined and used to calculate forces
   and torques based on contact geometry
*/

#include <cmath>
#include "contact_model.h"

namespace LAMMPS_NS {
  namespace Contact_Model{

enum {HOOKE, HERTZ, HERTZ_MATERIAL, DMT, JKR};
enum {VELOCITY, MASS_VELOCITY, VISCOELASTIC, TSUJI};
enum {TANGENTIAL_NOHISTORY, TANGENTIAL_HISTORY,
      TANGENTIAL_MINDLIN, TANGENTIAL_MINDLIN_RESCALE,
      TANGENTIAL_MINDLIN_FORCE, TANGENTIAL_MINDLIN_RESCALE_FORCE};
enum {TWIST_NONE, TWIST_SDS, TWIST_MARSHALL};
enum {ROLL_NONE, ROLL_SDS};

#define PI27SQ 266.47931882941264802866    // 27*PI**2
#define THREEROOT3 5.19615242270663202362  // 3*sqrt(3)
#define SIXROOT6 14.69693845669906728801   // 6*sqrt(6)
#define INVROOT6 0.40824829046386307274    // 1/sqrt(6)
#define FOURTHIRDS (4.0/3.0)                 // 4/3
#define THREEQUARTERS 0.75                 // 3/4

#define EPSILON 1e-10

ContactModel::ContactModel()
{

}

/* ----------------------------------------------------------------------
   get volume-correct r basis in: basis*cbrt(vol) = q*r
------------------------------------------------------------------------- */

void ContactModel::()
{

}


}}
