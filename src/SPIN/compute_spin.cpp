// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "compute_spin.h"

#include "atom.h"
#include "error.h"
#include "fix_precession_spin.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "pair_hybrid.h"
#include "pair_spin.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeSpin::ComputeSpin(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), lockprecessionspin(nullptr), pair(nullptr), spin_pairs(nullptr)
{
  if ((narg != 3) && (narg != 4)) error->all(FLERR,"Illegal compute compute/spin command");

  vector_flag = 1;
  size_vector = 6;
  extvector = 0;

  // initialize the magnetic interaction flags

  pair_spin_flag = 0;
  long_spin_flag = 0;
  precession_spin_flag = 0;

  ComputeSpin::init();

  allocate();

}

/* ---------------------------------------------------------------------- */

ComputeSpin::~ComputeSpin()
{
  memory->destroy(vector);
  delete[] spin_pairs;
  delete[] lockprecessionspin;
}

/* ---------------------------------------------------------------------- */

void ComputeSpin::init()
{
  hbar = force->hplanck/MY_2PI;
  kb = force->boltz;
  npairs = npairspin = 0;
  precession_spin_flag = 0;

  // set ptrs on Pair/Spin styles

  // loop 1: obtain # of Pairs, and # of Pair/Spin styles

  PairHybrid *hybrid = dynamic_cast<PairHybrid *>(force->pair_match("^hybrid",0));
  if (force->pair_match("^spin",0,0)) {        // only one Pair/Spin style
    pair = force->pair_match("^spin",0,0);
    if (hybrid == nullptr) npairs = 1;
    else npairs = hybrid->nstyles;
    npairspin = 1;
  } else if (force->pair_match("^spin",0,1)) { // more than one Pair/Spin style
    pair = force->pair_match("^spin",0,1);
    if (hybrid == nullptr) npairs = 1;
    else npairs = hybrid->nstyles;
    for (int i = 0; i<npairs; i++) {
      if (force->pair_match("^spin",0,i)) {
        npairspin++;
      }
    }
  }

  // init length of vector of ptrs to Pair/Spin styles

  if (npairspin > 0) {
    spin_pairs = new PairSpin*[npairspin];
  }

  // loop 2: fill vector with ptrs to Pair/Spin styles

  int count = 0;
  if (npairspin == 1) {
    count = 1;
    spin_pairs[0] = dynamic_cast<PairSpin *>(force->pair_match("^spin",0,0));
  } else if (npairspin > 1) {
    for (int i = 0; i<npairs; i++) {
      if (force->pair_match("^spin",0,i)) {
        spin_pairs[count] = dynamic_cast<PairSpin *>(force->pair_match("^spin",0,i));
        count++;
      }
    }
  }

  if (count != npairspin)
    error->all(FLERR,"Incorrect number of spin pair styles");

  // set pair/spin and long/spin flags

  if (npairspin >= 1) pair_spin_flag = 1;

  for (int i = 0; i<npairs; i++) {
    if (force->pair_match("^spin/long",0,i)) {
      long_spin_flag = 1;
    }
  }

  // set ptrs for fix precession/spin styles

  auto precfixes = modify->get_fix_by_style("^precession/spin");
  nprecspin = precfixes.size();

  if (nprecspin > 0) {
    lockprecessionspin = new FixPrecessionSpin *[nprecspin];
    precession_spin_flag = 1;

    int i = 0;
    for (auto &ifix : precfixes)
      lockprecessionspin[i++] = dynamic_cast<FixPrecessionSpin *>(ifix);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpin::compute_vector()
{
  int i;
  int countsp, countsptot;
  double mag[4], magtot[4];
  double magenergy, magenergytot;
  double tempnum, tempnumtot;
  double tempdenom, tempdenomtot;
  double spintemperature;

  invoked_vector = update->ntimestep;

  countsp = countsptot = 0.0;
  mag[0] = mag[1] = mag[2] = mag[3] = 0.0;
  magtot[0] = magtot[1] = magtot[2] = magtot[3] = 0.0;
  magenergy = magenergytot = 0.0;
  tempnum = tempnumtot = 0.0;
  tempdenom = tempdenomtot = 0.0;
  spintemperature = 0.0;

  int *mask = atom->mask;
  double **sp = atom->sp;
  double **fm = atom->fm;
  double tx,ty,tz;

  int nlocal = atom->nlocal;

  // compute total magnetization and magnetic energy
  // compute spin temperature (Nurdin et al., Phys. Rev. E 61, 2000)

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (atom->sp_flag) {

        // compute first moment

        mag[0] += sp[i][0];
        mag[1] += sp[i][1];
        mag[2] += sp[i][2];

        // update magnetic precession energies

        if (precession_spin_flag) {
          for (int k = 0; k < nprecspin; k++) {
            magenergy += lockprecessionspin[k]->emag[i];
          }
        }

        // update magnetic pair interactions

        if (pair_spin_flag) {
          for (int k = 0; k < npairspin; k++) {
            magenergy += spin_pairs[k]->emag[i];
          }
        }

        tx = sp[i][1]*fm[i][2]-sp[i][2]*fm[i][1];
        ty = sp[i][2]*fm[i][0]-sp[i][0]*fm[i][2];
        tz = sp[i][0]*fm[i][1]-sp[i][1]*fm[i][0];
        tempnum += tx*tx+ty*ty+tz*tz;
        tempdenom += sp[i][0]*fm[i][0]+fm[i][1]*sp[i][1]+sp[i][2]*fm[i][2];
        countsp++;
      }
    }
    else error->all(FLERR,"Compute compute/spin requires atom/spin style");
  }

  MPI_Allreduce(mag,magtot,4,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&magenergy,&magenergytot,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&tempnum,&tempnumtot,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&tempdenom,&tempdenomtot,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&countsp,&countsptot,1,MPI_INT,MPI_SUM,world);

  // compute average magnetization

  double scale = 1.0/countsptot;
  magtot[0] *= scale;
  magtot[1] *= scale;
  magtot[2] *= scale;
  magtot[3] = sqrt((magtot[0]*magtot[0])+(magtot[1]*magtot[1])+(magtot[2]*magtot[2]));

  // compute spin temperature

  spintemperature = hbar*tempnumtot;
  spintemperature /= (2.0*kb*tempdenomtot);

  vector[0] = magtot[0];
  vector[1] = magtot[1];
  vector[2] = magtot[2];
  vector[3] = magtot[3];
  vector[4] = magenergytot;
  vector[5] = spintemperature;
}

/* ----------------------------------------------------------------------
   free and reallocate arrays
------------------------------------------------------------------------- */

void ComputeSpin::allocate()
{
  memory->create(vector,size_vector,"compute/spin:vector");
}

