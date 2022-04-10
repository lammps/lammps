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
------------------------------------------------------------------------- */

/** Fix Drude Transform ******************************************************/

#include "fix_drude_transform.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_drude.h"
#include "modify.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */
template <bool inverse>
FixDrudeTransform<inverse>::FixDrudeTransform(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), mcoeff(nullptr)
{
  if (narg != 3) error->all(FLERR,"Illegal fix drude/transform command");
  comm_forward = 9;
  fix_drude = nullptr;
}

/* ---------------------------------------------------------------------- */
template <bool inverse>
FixDrudeTransform<inverse>::~FixDrudeTransform()
{
  delete[] mcoeff;
}

/* ---------------------------------------------------------------------- */
template <bool inverse>
void FixDrudeTransform<inverse>::init()
{
  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(modify->fix[ifix]->style,"drude") == 0) break;
  if (ifix == modify->nfix) error->all(FLERR, "fix drude/transform requires fix drude");
  fix_drude = (FixDrude *) modify->fix[ifix];
}

/* ---------------------------------------------------------------------- */
template <bool inverse>
int FixDrudeTransform<inverse>::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */
template <bool inverse>
void FixDrudeTransform<inverse>::setup(int) {
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;
  int * type = atom->type;
  double * rmass = atom->rmass, * mass = atom->mass;
  tagint * drudeid = fix_drude->drudeid;
  int * drudetype = fix_drude->drudetype;

  if (!rmass) {
    if (!mcoeff) mcoeff = new double[ntypes+1];
    auto mcoeff_loc = new double[ntypes+1];
    for (int itype=0; itype<=ntypes; itype++) mcoeff_loc[itype] = 2.; // an impossible value: mcoeff is at most 1.
    for (int i=0; i<nlocal; i++) {
      if (drudetype[type[i]] == DRUDE_TYPE) {
        int j = atom->map(drudeid[i]);
        // i is drude, j is core
        if (mcoeff_loc[type[i]] < 1.5) { // already done
          if (mcoeff_loc[type[j]] > 1.5) { // not yet done ??
            error->all(FLERR,"There must be one Drude type per core type");}
          continue;
        }
        mcoeff_loc[type[i]] = mass[type[i]] / (mass[type[i]] + mass[type[j]]);
        mcoeff_loc[type[j]] = -mass[type[i]] / mass[type[j]];
      }
    }

    MPI_Allreduce(mcoeff_loc, mcoeff, ntypes+1, MPI_DOUBLE, MPI_MIN, world);
    // mcoeff is 2 for non polarizable
    // 0 < mcoeff < 1 for drude
    // mcoeff < 0 for core
    delete[] mcoeff_loc;
  }
}

/* ---------------------------------------------------------------------- */
namespace LAMMPS_NS { // required for specialization
template <>
void FixDrudeTransform<false>::initial_integrate(int) {
  comm->forward_comm(this);
  real_to_reduced();
  //comm->forward_comm(this); // Normally not needed
}

template <>
void FixDrudeTransform<false>::final_integrate() {
  comm->forward_comm(this);
  real_to_reduced();
  //comm->forward_comm(this); // Normally not needed
}

template <>
void FixDrudeTransform<true>::initial_integrate(int) {
  comm->forward_comm(this);
  reduced_to_real();
  //comm->forward_comm(this); // Normally not needed
}

template <>
void FixDrudeTransform<true>::final_integrate() {
  comm->forward_comm(this);
  reduced_to_real();
  //comm->forward_comm(this); // Normally not needed
}

} // end of namespace

/* ---------------------------------------------------------------------- */
template <bool inverse>
void FixDrudeTransform<inverse>::real_to_reduced()
{
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;
  int dim = domain->dimension;
  int * mask = atom->mask, * type = atom->type;
  double ** x = atom->x, ** v = atom->v, ** f = atom->f;
  double * rmass = atom->rmass, * mass = atom->mass;
  double mcore, mdrude, coeff;
  int icore, idrude;
  tagint * drudeid = fix_drude->drudeid;
  int * drudetype = fix_drude->drudetype;

  if (!rmass) { // TODO: maybe drudetype can be used instead?
    for (int itype=1; itype<=ntypes; itype++)
      if (mcoeff[itype] < 1.5) mass[itype] *= 1. - mcoeff[itype];
  }
  for (int i=0; i<nlocal; i++) {
    if (mask[i] & groupbit && drudetype[type[i]] != NOPOL_TYPE) {
      drudeid[i] = (tagint) domain->closest_image(i, atom->map(drudeid[i]));
    }
  }
  for (int i=0; i<nlocal; i++) {
    if (mask[i] & groupbit && drudetype[type[i]] != NOPOL_TYPE) {
      int j = (int) drudeid[i];
      if (drudetype[type[i]] == DRUDE_TYPE && j < nlocal) continue;

      if (drudetype[type[i]] == DRUDE_TYPE) {
        idrude = i;
        icore = j;
      } else {
        icore = i;
        idrude = j;
      }
      if (rmass) {
        mcore = rmass[icore];
        mdrude = rmass[idrude];
        rmass[icore] += mdrude;
        rmass[idrude] *= mcore / rmass[icore];
        coeff = mdrude / (mcore + mdrude);
      } else { // TODO check that all atoms of this types are in the group
        coeff = mcoeff[type[idrude]];
      }
      for (int k=0; k<dim; k++) {
        x[idrude][k] -= x[icore][k];
        x[icore][k] += coeff * x[idrude][k];
        v[idrude][k] -= v[icore][k];
        v[icore][k] += coeff * v[idrude][k];
        f[icore][k] += f[idrude][k];
        f[idrude][k] -= coeff * f[icore][k];
      }
    }
  }
  fix_drude->is_reduced = true;
}

/* ---------------------------------------------------------------------- */
template <bool inverse>
void FixDrudeTransform<inverse>::reduced_to_real()
{
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;
  int dim = domain->dimension;
  int * mask = atom->mask, * type = atom->type;
  double ** x = atom->x, ** v = atom->v, ** f = atom->f;
  double * rmass = atom->rmass, * mass = atom->mass;
  double mcore, mdrude, coeff;
  int icore, idrude;
  tagint * drudeid = fix_drude->drudeid;
  int * drudetype = fix_drude->drudetype;

  for (int i=0; i<nlocal; i++) {
    if (mask[i] & groupbit && drudetype[type[i]] != NOPOL_TYPE) {
      int j = (int) drudeid[i]; // local index of drude partner because drudeid is in reduced form
      if (drudetype[type[i]] == DRUDE_TYPE && j < nlocal) continue;

      if (drudetype[type[i]] == DRUDE_TYPE) {
        idrude = i;
        icore = j;
      } else {
        icore = i;
        idrude = j;
      }
      if (rmass) {
        double s = sqrt(1. - rmass[idrude]/rmass[icore]);
        rmass[idrude] = 0.5 * rmass[icore] * (1. - s);
        mdrude = rmass[idrude];
        rmass[icore] -= mdrude;
        mcore = rmass[icore];
        coeff = mdrude / (mcore + mdrude);
      } else {
        if (!mcoeff[type[icore]]) { // TODO: should it be > 1.5 ?
          double s = sqrt(1. - mass[type[idrude]] / mass[type[icore]]);
          mass[type[idrude]] = 0.5 * mass[type[icore]] * (1. - s);
          mdrude = mass[type[idrude]];
          mass[type[icore]] -= mdrude;
          mcore = mass[type[icore]];
          mcoeff[type[icore]] = mdrude / (mcore + mdrude);
        }
        coeff = mcoeff[type[idrude]];
      }
      for (int k=0; k<dim; k++) {
        x[icore][k] -= coeff * x[idrude][k];
        x[idrude][k] += x[icore][k];
        v[icore][k] -= coeff * v[idrude][k];
        v[idrude][k] += v[icore][k];
        f[idrude][k] += coeff * f[icore][k];
        f[icore][k] -= f[idrude][k];
      }
    }
  }
  for (int i=0; i<nlocal; i++) {
    if (mask[i] & groupbit && drudetype[type[i]] != NOPOL_TYPE) {
      drudeid[i] = atom->tag[(int) drudeid[i]];
    }
  }
  if (!rmass) {
    for (int itype=1; itype<=ntypes; itype++)
      if (mcoeff[itype] < 1.5) mass[itype] /= 1. - mcoeff[itype];
  }
  fix_drude->is_reduced = false;
}

/* ---------------------------------------------------------------------- */
template <bool inverse>
int FixDrudeTransform<inverse>::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  double ** x = atom->x, ** v = atom->v, ** f = atom->f;
  int * type = atom->type, * drudetype = fix_drude->drudetype;
  double dx,dy,dz;
  int dim = domain->dimension;
  int m = 0;
  for (int i=0; i<n; i++) {
    int j = list[i];
    if (pbc_flag == 0 ||
        (fix_drude->is_reduced && drudetype[type[j]] == DRUDE_TYPE)) {
        for (int k=0; k<dim; k++) buf[m++] = x[j][k];
    }
    else {
        if (domain->triclinic != 0) {
            dx = pbc[0]*domain->xprd + pbc[5]*domain->xy;
            dy = pbc[1]*domain->yprd;
            if (dim == 3) {
                dx += + pbc[4]*domain->xz;
                dy += pbc[3]*domain->yz;
                dz = pbc[2]*domain->zprd;
            }
        }
        else {
            dx = pbc[0]*domain->xprd;
            dy = pbc[1]*domain->yprd;
            if (dim == 3)
                dz = pbc[2]*domain->zprd;
        }
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        if (dim == 3)
            buf[m++] = x[j][2] + dz;
    }
    for (int k=0; k<dim; k++) buf[m++] = v[j][k];
    for (int k=0; k<dim; k++) buf[m++] = f[j][k];
  }
  return m;
}

/* ---------------------------------------------------------------------- */
template <bool inverse>
void FixDrudeTransform<inverse>::unpack_forward_comm(int n, int first, double *buf)
{
  double ** x = atom->x, ** v = atom->v, ** f = atom->f;
  int dim = domain->dimension;
  int m = 0;
  int last = first + n;
  for (int i=first; i<last; i++) {
    for (int k=0; k<dim; k++) x[i][k] = buf[m++];
    for (int k=0; k<dim; k++) v[i][k] = buf[m++];
    for (int k=0; k<dim; k++) f[i][k] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */
template class LAMMPS_NS::FixDrudeTransform<false>;
template class LAMMPS_NS::FixDrudeTransform<true>;

