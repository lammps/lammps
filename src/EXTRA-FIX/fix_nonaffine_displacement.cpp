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

/* ----------------------------------------------------------------------
   Contributing authors: Joel Clemmer (SNL), Ishan Srivastava (LBNL)
------------------------------------------------------------------------- */

#include "fix_nonaffine_displacement.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_store_atom.h"
#include "force.h"
#include "group.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

enum { TYPE, RADIUS, CUSTOM };
enum { INTEGRATED, D2MIN };
enum { FIXED, OFFSET, UPDATE };

static const char cite_nonaffine_d2min[] =
  "@article{PhysRevE.57.7192,\n"
  " title = {Dynamics of viscoplastic deformation in amorphous solids},\n"
  " author = {Falk, M. L. and Langer, J. S.},\n"
  " journal = {Phys. Rev. E},\n"
  " volume = {57},\n"
  " issue = {6},\n"
  " pages = {7192--7205},\n"
  " numpages = {0},\n"
  " year = {1998},\n"
  " month = {Jun},\n"
  " publisher = {American Physical Society},\n"
  " doi = {10.1103/PhysRevE.57.7192},\n"
  "url = {https://link.aps.org/doi/10.1103/PhysRevE.57.7192}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixNonaffineDisplacement::FixNonaffineDisplacement(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), id_fix(nullptr), X(nullptr), Y(nullptr), F(nullptr), norm(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal fix nonaffine/displacement command");

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal nevery value {} in fix nonaffine/displacement", nevery);

  reference_timestep = update_timestep = offset_timestep = -1;
  int iarg = 4;
  if (strcmp(arg[iarg], "integrated") == 0) {
    nad_style = INTEGRATED;
    nevery = 1;
    iarg += 1;
  } else if (strcmp(arg[iarg], "d2min") == 0) {
    if (iarg + 1 > narg) error->all(FLERR,"Illegal fix nonaffine/displacement command");
    nad_style = D2MIN;
    if (strcmp(arg[iarg + 1], "type") == 0) {
      cut_style = TYPE;
    } else if (strcmp(arg[iarg + 1], "radius") == 0) {
      cut_style = RADIUS;
    } else if (strcmp(arg[iarg + 1], "custom") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Illegal fix nonaffine/displacement command");
      cut_style = CUSTOM;
      cutoff_custom = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      cutsq_custom = cutoff_custom * cutoff_custom;
      if (cutoff_custom <= 0)
        error->all(FLERR, "Illegal custom cutoff length {}", arg[iarg + 2]);
      iarg += 1;
    } else error->all(FLERR,"Illegal cutoff style {} in fix nonaffine/displacement", arg[iarg + 1]);
    iarg += 2;
  } else error->all(FLERR,"Illegal nonaffine displacement style {} in fix nonaffine/displacement", arg[iarg]);

  if (iarg + 2 > narg) error->all(FLERR,"Illegal fix nonaffine/displacement command");
  if (strcmp(arg[iarg], "fixed") == 0) {
    reference_style = FIXED;
    reference_timestep = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
    if (reference_timestep < 0)
      error->all(FLERR, "Illegal reference timestep {} in fix nonaffine/displacement", arg[iarg + 1]);
  } else if (strcmp(arg[iarg], "update") == 0) {
    reference_style = UPDATE;
    update_timestep = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
    if (update_timestep < 0)
      error->all(FLERR, "Illegal update timestep {} in fix nonaffine/displacement", arg[iarg + 1]);
  } else if (strcmp(arg[iarg], "offset") == 0) {
    reference_style = OFFSET;
    offset_timestep = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
    if ((offset_timestep <= 0) || (offset_timestep > nevery))
      error->all(FLERR, "Illegal offset timestep {} in fix nonaffine/displacement", arg[iarg + 1]);
  } else error->all(FLERR,"Illegal reference style {} in fix nonaffine/displacement", arg[iarg]);

  if (nad_style == D2MIN)
    if (cut_style == RADIUS && (!atom->radius_flag))
      error->all(FLERR, "Fix nonaffine/displacement radius style requires atom attribute radius");

  if (nad_style == INTEGRATED && reference_style == OFFSET)
    error->all(FLERR, "Fix nonaffine/displacement cannot use the integrated style with an offset reference state");

  peratom_flag = 1;
  peratom_freq = nevery;
  nmax = -1;
  reference_saved = 0;
  restart_global = 1;

  size_peratom_cols = 3;
  comm_reverse = 0;
  comm_forward = 0;
  if (nad_style == D2MIN) {
    comm_reverse = 18;
    comm_forward = 9;
  }

  if (nad_style == D2MIN && lmp->citeme) lmp->citeme->add(cite_nonaffine_d2min);
}

/* ---------------------------------------------------------------------- */

FixNonaffineDisplacement::~FixNonaffineDisplacement()
{
  if (id_fix && modify->nfix) modify->delete_fix(id_fix);
  delete[] id_fix;

  if (nad_style == D2MIN) {
    memory->destroy(X);
    memory->destroy(Y);
    memory->destroy(F);
    memory->destroy(norm);
    memory->destroy(array_atom);
  }
}

/* ---------------------------------------------------------------------- */

int FixNonaffineDisplacement::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::post_constructor()
{
  // Create persistent peratom storage for either an integrated velocity or reference position
  // Ghost atoms need reference coordinates for D2min
  std::string ghost_status = "0";
  if (nad_style == D2MIN) ghost_status = "1";

  id_fix = utils::strdup(id + std::string("_FIX_PA"));
  fix = dynamic_cast<FixStoreAtom *>(modify->add_fix(fmt::format("{} {} STORE/ATOM 3 0 {} 1", id_fix, group->names[igroup], ghost_status)));

  if (nad_style == INTEGRATED)
    array_atom = fix->astore;

  if (nad_style == D2MIN)
    grow_arrays(atom->nmax);

  for (int i = 0; i < atom->nlocal; i++)
    for (int j = 0; j < 3; j++) array_atom[i][j] = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::init()
{
  dtv = update->dt;

  if ((!reference_saved) && (reference_style == FIXED) && (update->ntimestep > reference_timestep))
    error->all(FLERR, "Initial timestep exceeds that of the reference state in fix nonaffine/displacement");

  if (nad_style == D2MIN) {
    if ((!force->pair) && (cut_style == TYPE))
    error->all(FLERR,"Fix nonaffine/displacement D2Min option requires a pair style be defined "
               "or cutoff specified");

    // need an occasional half neighbor list

    if (cut_style == RADIUS) {
      neighbor->add_request(this, NeighConst::REQ_SIZE | NeighConst::REQ_OCCASIONAL);
    } else {
      auto req = neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
      if (cut_style == CUSTOM) {
        double skin = neighbor->skin;
        mycutneigh = cutoff_custom + skin;

        double cutghost;            // as computed by Neighbor and Comm
        if (force->pair)
          cutghost = MAX(force->pair->cutforce + skin, comm->cutghostuser);
        else
          cutghost = comm->cutghostuser;

        if (mycutneigh > cutghost)
          error->all(FLERR,"Fix nonaffine/displacement D2Min option cutoff exceeds ghost atom range - use comm_modify cutoff command");

        req->set_cutoff(mycutneigh);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::setup(int /*vflag*/)
{
  post_force(0); // Save state if needed before starting the 1st timestep
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::post_force(int /*vflag*/)
{
  if (reference_saved && (!update->setupflag)) {
    if (nad_style == INTEGRATED) {
      integrate_velocity();
    } else {
      if ((update->ntimestep % nevery) == 0) calculate_D2Min();
    }
  }

  if (reference_style == FIXED)
    if (update->ntimestep == reference_timestep)
      save_reference_state();

  if (reference_style == UPDATE)
    if ((update->ntimestep % update_timestep) == 0)
      save_reference_state();

  if (reference_style == OFFSET)
    if (((update->ntimestep + offset_timestep) % nevery) == 0)
      save_reference_state();
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::write_restart(FILE *fp)
{
  if (comm->me == 0) {
    int size = sizeof(int);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(&reference_saved, sizeof(int), 1, fp);
  }
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::restart(char *buf)
{
  reference_saved = (int) ubuf(buf[0]).i;
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::integrate_velocity()
{
  dtv = update->dt;

  double **v = atom->v;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int m = 0; m < 3; m++) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        array_atom[i][m] += dtv * v[i][m];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::save_reference_state()
{
  double **x = atom->x;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (nad_style == D2MIN) {
    for (int m = 0; m < 3; m++) {
      for (int i = 0; i < nall; i++) {
        if (mask[i] & groupbit)  array_atom[i][m] = x[i][m];
      }
    }
  } else {
    for (int m = 0; m < 3; m++) {
      for (int i = 0; i < nall; i++) {
        if (mask[i] & groupbit)  array_atom[i][m] = 0.0;
      }
    }
  }

  if (nad_style == D2MIN) {
    xprd0 = domain->xprd;
    yprd0 = domain->yprd;
    zprd0 = domain->zprd;
    xprd0_half = domain->xprd_half;
    yprd0_half = domain->yprd_half;
    zprd0_half = domain->zprd_half;
    xy0 = domain->xy;
    xz0 = domain->xz;
    yz0 = domain->yz;
  }

  reference_saved = 1;
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::calculate_D2Min()
{
  // invoke half neighbor list (will copy or build if necessary)
  neighbor->build_one(list);

  if (atom->nmax > nmax)
    grow_arrays(atom->nmax);

  int i, j, k, l, ii, jj, inum, jnum, itype, jtype;
  double evol, j2, edev;
  double r[3], r0[3], rsq, radsum, temp[3];
  double X_tmp[3][3], Y_tmp[3][3], F_tmp[3][3], E[3][3];
  double Y_inv[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}}; // Zero for 2d since not all entries used
  int *ilist, *jlist, *numneigh, **firstneigh;

  double **x = atom->x;
  double **x0 = array_atom;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int dim = domain->dimension;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  Pair *pair = force->pair;
  double **cutsq;
  if (pair) cutsq = force->pair->cutsq;

  for (i = 0; i < nmax; i++) {
    for (k = 0; k < 3; k++) {
      for (l = 0; l < 3; l++) {
        X[i][k][l] = 0.0;
        Y[i][k][l] = 0.0;
      }
    }
    norm[i] = 0;
    array_atom[i][0] = 0;
  }

  // First loop through neighbors
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;

      jtype = type[j];
      r[0] = x[i][0] - x[j][0];
      r[1] = x[i][1] - x[j][1];
      r[2] = x[i][2] - x[j][2];
      rsq = lensq3(r);

      // Only include contributions from atoms that are CURRENTLY neighbors
      if (cut_style == TYPE) {
        if (rsq > cutsq[itype][jtype]) continue;
      } else if (cut_style == CUSTOM) {
        if (rsq > cutsq_custom) continue;
      } else {
        radsum = radius[i] + radius[j];
        if (rsq > (radsum * radsum)) continue;
      }

      r0[0] = x0[i][0] - x0[j][0];
      r0[1] = x0[i][1] - x0[j][1];
      r0[2] = x0[i][2] - x0[j][2];
      minimum_image0(r0);

      // Using notation from Falk & Langer 1998
      outer3(r, r0, X_tmp);
      outer3(r0, r0, Y_tmp);

      for (k = 0; k < 3; k++) {
        for (l = 0; l < 3; l++) {
          X[i][k][l] += X_tmp[k][l];
          Y[i][k][l] += Y_tmp[k][l];
        }
      }

      if (newton_pair || j < nlocal) {
        for (k = 0; k < 3; k++) {
          for (l = 0; l < 3; l++) {
            X[j][k][l] += X_tmp[k][l];
            Y[j][k][l] += Y_tmp[k][l];
          }
        }
      }
    }
  }

  comm_flag = 0;
  if (newton_pair) comm->reverse_comm(this, 18);

  // Calculate contributions to strain tensor
  double denom;
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        Y_tmp[j][k] = Y[i][j][k];
        X_tmp[j][k] = X[i][j][k];
      }
    }

    if (dim == 3) {
      invert3(Y_tmp, Y_inv);
    } else {
      denom = Y_tmp[0][0] * Y_tmp[1][1] - Y_tmp[0][1] * Y_tmp[1][0];
      if (denom != 0.0) denom = 1.0 / denom;
      Y_inv[0][0] = Y_tmp[1][1] * denom;
      Y_inv[0][1] = -Y_tmp[0][1] * denom;
      Y_inv[1][0] = -Y_tmp[1][0] * denom;
      Y_inv[1][1] = Y_tmp[0][0] * denom;
    }

    times3(X_tmp, Y_inv, F_tmp);

    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        F[i][j][k] = F_tmp[j][k];
      }
    }
  }

  comm->forward_comm(this);

  // Second loop through neighbors
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;

      jtype = type[j];
      r[0] = x[i][0] - x[j][0];
      r[1] = x[i][1] - x[j][1];
      r[2] = x[i][2] - x[j][2];
      rsq = lensq3(r);

      // Only include contributions from atoms that are CURRENTLY neighbors
      if (cut_style == TYPE) {
        if (rsq >= cutsq[itype][jtype]) continue;
      } else if (cut_style == CUSTOM) {
        if (rsq >= cutsq_custom) continue;
      } else {
        radsum = radius[i] + radius[j];
        if (rsq >= radsum * radsum) continue;
      }

      r0[0] = x0[i][0] - x0[j][0];
      r0[1] = x0[i][1] - x0[j][1];
      r0[2] = x0[i][2] - x0[j][2];
      minimum_image0(r0);

      // E * r0
      for (k = 0; k < 3; k++) {
        temp[k] = 0.0;
        for (l = 0; l < 3; l++)
          temp[k] += F[i][k][l] * r0[l];
      }

      sub3(r, temp, temp);
      array_atom[i][0] += lensq3(temp);
      norm[i] += 1;

      if (newton_pair || j < nlocal) {
        for (k = 0; k < 3; k++) {
          temp[k] = 0.0;
          for (l = 0; l < 3; l++)
            temp[k] += F[j][k][l] * r0[l];
        }

        sub3(r, temp, temp);
        array_atom[j][0] += lensq3(temp);
        norm[j] += 1;
      }
    }
  }

  comm_flag = 1;
  if (newton_pair) comm->reverse_comm(this, 2);

  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    if (norm[i] != 0)
      array_atom[i][0] /= norm[i];
    else
      array_atom[i][0] = 0.0;
    array_atom[i][0] = sqrt(array_atom[i][0]);

    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        F_tmp[j][k] = F[i][j][k];

    transpose_times3(F_tmp, F_tmp, E);
    for (j = 0; j < dim; j++) E[j][j] -= 1.0;

    evol = (E[0][0] + E[1][1] + E[2][2]) / dim;

    // Calculate deviatoric strain
    for (j = 0; j < dim; j++) E[j][j] -= evol;
    j2 = 0.0;
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        j2 += E[j][k] * E[j][k];

    edev = sqrt(0.5 * j2);

    array_atom[i][1] = evol;
    array_atom[i][2] = edev;
  }
}

/* ---------------------------------------------------------------------- */

int FixNonaffineDisplacement::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last, k, l;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (comm_flag == 0) {
      for (k = 0; k < 3; k++) {
        for (l = 0; l < 3; l++) {
          buf[m++] = X[i][k][l];
          buf[m++] = Y[i][k][l];
        }
      }
    } else {
      buf[m++] = array_atom[i][0];
      buf[m++] = ubuf(norm[i]).d;
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m, k, l;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (comm_flag == 0) {
      for (k = 0; k < 3; k++) {
        for (l = 0; l < 3; l++) {
          X[j][k][l] += buf[m++];
          Y[j][k][l] += buf[m++];
        }
      }
    } else {
      array_atom[j][0] += buf[m++];
      norm[j] += (int) ubuf(buf[m++]).i;
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixNonaffineDisplacement::pack_forward_comm(int n, int *list, double *buf,
                                          int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m, k, l;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 0; k < 3; k++) {
      for (l = 0; l < 3; l ++) {
        buf[m++] = F[j][k][l];
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last, k, l;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (k = 0; k < 3; k++) {
      for (l = 0; l < 3; l ++) {
        F[i][k][l] = buf[m++];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::minimum_image0(double *delta)
{
  if (domain->triclinic == 0) {
    if (domain->xperiodic) {
      while (fabs(delta[0]) > xprd0_half) {
        if (delta[0] < 0.0) delta[0] += xprd0;
        else delta[0] -= xprd0;
      }
    }
    if (domain->yperiodic) {
      while (fabs(delta[1]) > yprd0_half) {
        if (delta[1] < 0.0) delta[1] += yprd0;
        else delta[1] -= yprd0;
      }
    }
    if (domain->zperiodic) {
      while (fabs(delta[2]) > zprd0_half) {
        if (delta[2] < 0.0) delta[2] += zprd0;
        else delta[2] -= zprd0;
      }
    }

  } else {
    if (domain->zperiodic) {
      while (fabs(delta[2]) > zprd0_half) {
        if (delta[2] < 0.0) {
          delta[2] += zprd0;
          delta[1] += yz0;
          delta[0] += xz0;
        } else {
          delta[2] -= zprd0;
          delta[1] -= yz0;
          delta[0] -= xz0;
        }
      }
    }
    if (domain->yperiodic) {
      while (fabs(delta[1]) > yprd0_half) {
        if (delta[1] < 0.0) {
          delta[1] += yprd0;
          delta[0] += xy0;
        } else {
          delta[1] -= yprd0;
          delta[0] -= xy0;
        }
      }
    }
    if (domain->xperiodic) {
      while (fabs(delta[0]) > xprd0_half) {
        if (delta[0] < 0.0) delta[0] += xprd0;
        else delta[0] -= xprd0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNonaffineDisplacement::grow_arrays(int nmax_new)
{
  nmax = nmax_new;
  memory->destroy(X);
  memory->destroy(Y);
  memory->destroy(F);
  memory->destroy(norm);
  memory->create(X, nmax, 3, 3, "fix_nonaffine_displacement:X");
  memory->create(Y, nmax, 3, 3, "fix_nonaffine_displacement:Y");
  memory->create(F, nmax, 3, 3, "fix_nonaffine_displacement:F");
  memory->create(norm, nmax, "fix_nonaffine_displacement:norm");
}
