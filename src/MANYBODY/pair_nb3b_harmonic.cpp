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
   Contributing author: Todd R. Zeitler (SNL)
   (based on Stillinger-Weber pair style)
------------------------------------------------------------------------- */

#include "pair_nb3b_harmonic.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::MY_PI;

#define DELTA 4
#define SMALL 0.001

static const char *substyle[] = {"nb3n/harmonic", "nb3b/screened"};

/* ---------------------------------------------------------------------- */

PairNb3bHarmonic::PairNb3bHarmonic(LAMMPS *lmp) : Pair(lmp), params(nullptr)
{
  variant = HARMONIC;
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  unit_convert_flag = utils::get_supported_conversions(utils::ENERGY);
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairNb3bHarmonic::~PairNb3bHarmonic()
{
  memory->destroy(params);
  memory->destroy(elem3param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairNb3bHarmonic::compute(int eflag, int vflag)
{
  int i, j, k, ii, jj, kk, inum, jnum, jnumm1;
  int itype, jtype, ktype, ijparam, ikparam, ijkparam;
  double xtmp, ytmp, ztmp, evdwl;
  double rsq1, rsq2;
  double delr1[3], delr2[3], fj[3], fk[3];
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    jnumm1 = jnum - 1;

    for (jj = 0; jj < jnumm1; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2];
      if (rsq1 > params[ijparam].cutsq) continue;

      for (kk = jj + 1; kk < jnum; kk++) {
        k = jlist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
        if (rsq2 > params[ikparam].cutsq) continue;

        threebody(&params[ijparam], &params[ikparam], &params[ijkparam], rsq1, rsq2, delr1, delr2,
                  fj, fk, eflag, evdwl);

        f[i][0] -= fj[0] + fk[0];
        f[i][1] -= fj[1] + fk[1];
        f[i][2] -= fj[2] + fk[2];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        f[j][2] += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (evflag) ev_tally3(i, j, k, evdwl, 0.0, fj, fk, delr1, delr2);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairNb3bHarmonic::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  memory->create(cutsq, np1, np1, "pair:cutsq");

  map = new int[np1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairNb3bHarmonic::settings(int narg, char ** /*arg*/)
{
  if (narg != 0) error->all(FLERR, "Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairNb3bHarmonic::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  map_element2type(narg - 3, arg + 3);

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairNb3bHarmonic::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR, "Pair style {} requires atom IDs", substyle[variant]);
  if (force->newton_pair == 0)
    error->all(FLERR, "Pair style {} requires newton pair on", substyle[variant]);

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairNb3bHarmonic::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairNb3bHarmonic::read_file(char *file)
{
  memory->sfree(params);
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, substyle[variant], unit_convert_flag);
    char *line;

    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY, unit_convert);
    while ((line = reader.next_line(NPARAMS_PER_LINE + ((variant == SCREENED) ? 1 : 0)))) {
      try {
        ValueTokenizer values(line);

        std::string iname = values.next_string();
        std::string jname = values.next_string();
        std::string kname = values.next_string();

        // ielement,jelement,kelement = 1st args
        // if all 3 args are in element list, then parse this line
        // else skip to next entry in file
        int ielement, jelement, kelement;

        for (ielement = 0; ielement < nelements; ielement++)
          if (iname == elements[ielement]) break;
        if (ielement == nelements) continue;
        for (jelement = 0; jelement < nelements; jelement++)
          if (jname == elements[jelement]) break;
        if (jelement == nelements) continue;
        for (kelement = 0; kelement < nelements; kelement++)
          if (kname == elements[kelement]) break;
        if (kelement == nelements) continue;

        // load up parameter settings and error check their values

        if (nparams == maxparam) {
          maxparam += DELTA;
          params = (Param *) memory->srealloc(params, maxparam * sizeof(Param), "pair:params");

          // make certain all addional allocated storage is initialized
          // to avoid false positives when checking with valgrind

          memset(params + nparams, 0, DELTA * sizeof(Param));
        }

        params[nparams].ielement = ielement;
        params[nparams].jelement = jelement;
        params[nparams].kelement = kelement;
        params[nparams].k_theta = values.next_double();
        params[nparams].theta0 = values.next_double();
        params[nparams].invrho = 1.0;    // dummy value
        if (variant == SCREENED) {
          double rho = values.next_double();
          if (rho > 0.0)
            params[nparams].invrho = 1.0 / rho;
          else
            throw TokenizerException("Incorrect value for potential parameter", "rho");
        }
        params[nparams].cutoff = values.next_double();

        if (unit_convert) params[nparams].k_theta *= conversion_factor;
      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }

      if ((params[nparams].k_theta < 0.0) || (params[nparams].theta0 < 0.0) ||
          (params[nparams].cutoff < 0.0))
        error->one(FLERR, "Illegal {} parameter", substyle[variant]);

      nparams++;
    }
  }

  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxparam, 1, MPI_INT, 0, world);

  if (comm->me != 0) {
    params = (Param *) memory->srealloc(params, maxparam * sizeof(Param), "pair:params");
  }

  MPI_Bcast(params, maxparam * sizeof(Param), MPI_BYTE, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairNb3bHarmonic::setup_params()
{
  int i, j, k, m, n;
  double rtmp;

  // set elem3param for all triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem3param);
  memory->create(elem3param, nelements, nelements, nelements, "pair:elem3param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement && k == params[m].kelement) {
            if (n >= 0)
              error->all(FLERR, "Potential file has a duplicate entry for: {} {} {}", elements[i],
                         elements[j], elements[k]);
            n = m;
          }
        }
        if (n < 0)
          error->all(FLERR, "Potential file is missing an entry for: {} {} {}", elements[i],
                     elements[j], elements[k]);
        elem3param[i][j][k] = n;
      }

  // compute parameter values derived from inputs

  // set cutsq using shortcut to reduce neighbor list for accelerated
  // calculations. cut must remain unchanged as it is a potential parameter
  // (cut = a*sigma)

  for (m = 0; m < nparams; m++) {

    params[m].cut = params[m].cutoff;
    params[m].cutsq = params[m].cut * params[m].cut;

    params[m].theta0 = params[m].theta0 / 180.0 * MY_PI;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp = sqrt(params[m].cutsq);
    if (rtmp > cutmax) cutmax = rtmp;
  }
}

/* ---------------------------------------------------------------------- */

void PairNb3bHarmonic::threebody(Param * /*paramij*/, Param * /*paramik*/, Param *paramijk,
                                 double rsq1, double rsq2, double *delr1, double *delr2, double *fj,
                                 double *fk, int eflag, double &eng)
{
  double dtheta, tk;
  double r1, r2, c, s, a, a11, a12, a22;

  // angle (cos and sin)

  r1 = sqrt(rsq1);
  r2 = sqrt(rsq2);

  c = delr1[0] * delr2[0] + delr1[1] * delr2[1] + delr1[2] * delr2[2];
  c /= r1 * r2;

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  s = sqrt(1.0 - c * c);
  if (s < SMALL) s = SMALL;
  s = 1.0 / s;

  // force & energy

  dtheta = acos(c) - paramijk->theta0;
  tk = paramijk->k_theta * dtheta;

  if (eflag) eng = tk * dtheta;

  a = -2.0 * tk * s;
  a11 = a * c / rsq1;
  a12 = -a / (r1 * r2);
  a22 = a * c / rsq2;

  fj[0] = a11 * delr1[0] + a12 * delr2[0];
  fj[1] = a11 * delr1[1] + a12 * delr2[1];
  fj[2] = a11 * delr1[2] + a12 * delr2[2];
  fk[0] = a22 * delr2[0] + a12 * delr1[0];
  fk[1] = a22 * delr2[1] + a12 * delr1[1];
  fk[2] = a22 * delr2[2] + a12 * delr1[2];
}
