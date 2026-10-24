/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meissner (TUHH)
------------------------------------------------------------------------- */

#include "electrode_vector_tip4p.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "neigh_list.h"
#include "pair.h"

#include <cmath>

using namespace LAMMPS_NS;

ElectrodeVectorTIP4P::ElectrodeVectorTIP4P(LAMMPS *lmp, int sensor_group, int source_group,
                                           double eta, bool invert_source) :
    ElectrodeVector(lmp, sensor_group, source_group, eta, invert_source),
    qdist(0.0), cut_coul(0.0), alpha(0.0), typeO(0), typeH(0)
{
}

/* ---------------------------------------------------------------------- */

void ElectrodeVectorTIP4P::setup(class Pair *fix_pair, class NeighList *fix_neighlist,
                                 bool timer_flag)
{
  ElectrodeVector::setup(fix_pair, fix_neighlist, timer_flag);

  int dim = 0;
  auto *p_qdist = (double *) fix_pair->extract("qdist", dim);
  dim = 0;
  auto *p_typeO = (int *) fix_pair->extract("typeO", dim);
  dim = 0;
  auto *p_typeH = (int *) fix_pair->extract("typeH", dim);
  dim = 0;
  auto *p_typeA = (int *) fix_pair->extract("typeA", dim);
  dim = 0;
  auto *p_typeB = (int *) fix_pair->extract("typeB", dim);
  dim = 0;
  auto *p_cut_coul = (double *) fix_pair->extract("cut_coul", dim);

  if (!p_qdist || !p_typeO || !p_typeH || !p_typeA || !p_typeB || !p_cut_coul)
    error->all(FLERR, "Pair style is incompatible with TIP4P ELECTRODE vector");

  qdist = *p_qdist;
  typeO = *p_typeO;
  typeH = *p_typeH;
  cut_coul = *p_cut_coul;

  if (force->angle == nullptr || force->bond == nullptr || force->angle->setflag == nullptr ||
      force->bond->setflag == nullptr)
    error->all(FLERR, "Bond and angle potentials must be defined for TIP4P");

  int typeA = *p_typeA;
  int typeB = *p_typeB;

  if (typeA < 1 || typeA > atom->nangletypes || force->angle->setflag[typeA] == 0)
    error->all(FLERR, "Bad TIP4P angle type for ELECTRODE");

  if (typeB < 1 || typeB > atom->nbondtypes || force->bond->setflag[typeB] == 0)
    error->all(FLERR, "Bad TIP4P bond type for ELECTRODE");

  double theta = force->angle->equilibrium_angle(typeA);
  double blen = force->bond->equilibrium_distance(typeB);
  alpha = qdist / (cos(0.5 * theta) * blen);
}

/* ---------------------------------------------------------------------- */

void ElectrodeVectorTIP4P::charge_position(int i, double *xsite)
{
  if (atom->type[i] != typeO) {
    xsite[0] = atom->x[i][0];
    xsite[1] = atom->x[i][1];
    xsite[2] = atom->x[i][2];
    return;
  }

  int iH1, iH2;
  find_M(i, iH1, iH2, xsite);
}

/* ---------------------------------------------------------------------- */

void ElectrodeVectorTIP4P::charge_force(int i, const double *fcharge)
{
  if (atom->type[i] != typeO) {
    ElectrodeVector::charge_force(i, fcharge);
    return;
  }

  int iH1, iH2;
  double xM[3];
  find_M(i, iH1, iH2, xM);

  atom->f[i][0] += (1.0 - alpha) * fcharge[0];
  atom->f[i][1] += (1.0 - alpha) * fcharge[1];
  atom->f[i][2] += (1.0 - alpha) * fcharge[2];

  atom->f[iH1][0] += 0.5 * alpha * fcharge[0];
  atom->f[iH1][1] += 0.5 * alpha * fcharge[1];
  atom->f[iH1][2] += 0.5 * alpha * fcharge[2];

  atom->f[iH2][0] += 0.5 * alpha * fcharge[0];
  atom->f[iH2][1] += 0.5 * alpha * fcharge[1];
  atom->f[iH2][2] += 0.5 * alpha * fcharge[2];
}

/* ---------------------------------------------------------------------- */

double ElectrodeVectorTIP4P::charge_force_alpha() const
{
  return alpha;
}

/* ---------------------------------------------------------------------- */
int ElectrodeVectorTIP4P::charge_force_virial(int i, const double *fcharge,
                                              double *v, int *vlist)
{
  if (atom->type[i] != typeO) {
    return ElectrodeVector::charge_force_virial(i, fcharge, v, vlist);
  }

  int iH1, iH2;
  double xM[3];
  find_M(i, iH1, iH2, xM);

  const double fO[3] = {
      (1.0 - alpha) * fcharge[0],
      (1.0 - alpha) * fcharge[1],
      (1.0 - alpha) * fcharge[2]};
  const double fH[3] = {
      0.5 * alpha * fcharge[0],
      0.5 * alpha * fcharge[1],
      0.5 * alpha * fcharge[2]};

  const double *xO = atom->x[i];
  const double *xH1 = atom->x[iH1];
  const double *xH2 = atom->x[iH2];

  vlist[0] = i;
  vlist[1] = iH1;
  vlist[2] = iH2;

  v[0] = xO[0] * fO[0] + xH1[0] * fH[0] + xH2[0] * fH[0];
  v[1] = xO[1] * fO[1] + xH1[1] * fH[1] + xH2[1] * fH[1];
  v[2] = xO[2] * fO[2] + xH1[2] * fH[2] + xH2[2] * fH[2];
  v[3] = xO[0] * fO[1] + xH1[0] * fH[1] + xH2[0] * fH[1];
  v[4] = xO[0] * fO[2] + xH1[0] * fH[2] + xH2[0] * fH[2];
  v[5] = xO[1] * fO[2] + xH1[1] * fH[2] + xH2[1] * fH[2];

  return 3;
}

/* ---------------------------------------------------------------------- */
double ElectrodeVectorTIP4P::pair_cutsq(int, int) const
{
  return cut_coul * cut_coul;
}

/* ---------------------------------------------------------------------- */

void ElectrodeVectorTIP4P::find_M(int i, int &iH1, int &iH2, double *xM)
{
  double **x = atom->x;

  iH1 = atom->map(atom->tag[i] + 1);
  iH2 = atom->map(atom->tag[i] + 2);

  if (iH1 == -1 || iH2 == -1)
    error->one(FLERR, "TIP4P hydrogen is missing");

  if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
    error->one(FLERR, "TIP4P hydrogen has incorrect atom type");

  iH1 = domain->closest_image(i, iH1);
  iH2 = domain->closest_image(i, iH2);

  double delx1 = x[iH1][0] - x[i][0];
  double dely1 = x[iH1][1] - x[i][1];
  double delz1 = x[iH1][2] - x[i][2];

  double delx2 = x[iH2][0] - x[i][0];
  double dely2 = x[iH2][1] - x[i][1];
  double delz2 = x[iH2][2] - x[i][2];

  xM[0] = x[i][0] + alpha * 0.5 * (delx1 + delx2);
  xM[1] = x[i][1] + alpha * 0.5 * (dely1 + dely2);
  xM[2] = x[i][2] + alpha * 0.5 * (delz1 + delz2);
}
