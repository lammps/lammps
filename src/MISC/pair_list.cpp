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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_list.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_special.h"
#include "memory.h"
#include "text_file_reader.h"

#include <cmath>
#include <cstring>
#include <exception>
#include <map>

using namespace LAMMPS_NS;
using MathSpecial::square;

enum { NONE = 0, HARM, MORSE, LJ126, QUARTIC };

// clang-format off
static std::map<std::string, int> stylename = {
    {"none",     NONE},
    {"harmonic", HARM},
    {"morse",    MORSE},
    {"lj126",    LJ126},
    {"quartic",  QUARTIC}
};
// clang-format on

// fast power function for integer exponent > 0
static double mypow(double x, int n)
{
  double yy;

  if (x == 0.0) return 0.0;

  for (yy = 1.0; n != 0; n >>= 1, x *= x)
    if (n & 1) yy *= x;

  return yy;
}

typedef struct {
  double x, y, z;
} dbl3_t;

/* ---------------------------------------------------------------------- */

PairList::PairList(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  respa_enable = 0;
  cut_global = 0.0;
  params = nullptr;
  check_flag = 1;
}

/* ---------------------------------------------------------------------- */

PairList::~PairList()
{
  memory->destroy(setflag);
  memory->destroy(cutsq);
  memory->destroy(params);
}

/* ----------------------------------------------------------------------
   in this pair style we don't use a neighbor list, but loop through
   a list of pairwise interactions, determines the corresponding local
   atom indices and compute those forces.
------------------------------------------------------------------------- */

void PairList::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  // get maximum allowed tag.

  bigint maxtag_one, maxtag;
  maxtag_one = maxtag = 0;
  const int nlocal = atom->nlocal;
  const tagint *_noalias const tag = atom->tag;
  for (int i = 0; i < nlocal; ++i) maxtag_one = MAX(maxtag_one, tag[i]);
  MPI_Allreduce(&maxtag_one, &maxtag, 1, MPI_LMP_TAGINT, MPI_MAX, world);

  const int newton_pair = force->newton_pair;
  const dbl3_t *_noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t *_noalias const f = (dbl3_t *) atom->f[0];    // NOLINT

  double fpair, epair;
  int i, j;

  int pc = 0;
  for (int n = 0; n < npairs; ++n) {
    const list_param &par = params[n];

    // can only use valid tags or else atom->map() below will segfault.
    if ((par.id1 < 1) || (par.id1 > maxtag)) {
      if (check_flag)
        error->all(FLERR, "Invalid pair list atom ID {}", par.id1);
      else
        continue;
    }
    if ((par.id2 < 1) || (par.id2 > maxtag)) {
      if (check_flag)
        error->all(FLERR, "Invalid pair list atom ID {}", par.id2);
      else
        continue;
    }

    i = atom->map(par.id1);
    j = atom->map(par.id2);

    // if one of the two atoms is missing on the node skip
    if ((i < 0) || (j < 0)) continue;

    // both atoms are ghosts -> skip
    if ((i >= nlocal) && (j >= nlocal)) continue;

    // with newton pair and one ghost we have to skip half the cases.
    // if id1 is a ghost, we skip if the sum of both ids is even.
    // if id2 is a ghost, we skip if the sum of both ids is odd.
    if (newton_pair) {
      if ((i >= nlocal) && ((par.id1 + par.id2) & 1) == 0) continue;
      if ((j >= nlocal) && ((par.id1 + par.id2) & 1) == 1) continue;
    }

    const double dx = x[i].x - x[j].x;
    const double dy = x[i].y - x[j].y;
    const double dz = x[i].z - x[j].z;
    const double rsq = dx * dx + dy * dy + dz * dz;

    fpair = epair = 0.0;
    if (check_flag) {
      if (newton_pair || i < nlocal) ++pc;
      if (newton_pair || j < nlocal) ++pc;
    }

    if (rsq < par.cutsq) {
      const double r2inv = 1.0 / rsq;

      if (par.style == HARM) {
        const double r = sqrt(rsq);
        const double dr = par.param.harm.r0 - r;
        fpair = 2.0 * par.param.harm.k * dr / r;

        if (eflag_either) epair = par.param.harm.k * dr * dr - par.offset;

      } else if (par.style == MORSE) {

        const double r = sqrt(rsq);
        const double dr = r - par.param.morse.r0;
        const double dexp = exp(-par.param.morse.alpha * dr);
        fpair = 2.0 * par.param.morse.d0 * par.param.morse.alpha * (dexp * dexp - dexp) / r;

        if (eflag_either)
          epair = par.param.morse.d0 * (dexp * dexp - 2.0 * dexp + 1.0) - par.offset;

      } else if (par.style == LJ126) {

        const double r6inv = r2inv * r2inv * r2inv;
        const double sig6 = mypow(par.param.lj126.sigma, 6);
        fpair = 24.0 * par.param.lj126.epsilon * r6inv * (2.0 * sig6 * sig6 * r6inv - sig6) * r2inv;

        if (eflag_either)
          epair = 4.0 * par.param.lj126.epsilon * r6inv * (sig6 * sig6 * r6inv - sig6) - par.offset;

      } else if (par.style == QUARTIC) {

        const double r = sqrt(rsq);
        double dr = r - par.param.quartic.r0;
        double ra = dr - par.param.quartic.b1;
        double rb = dr - par.param.quartic.b2;
        double r2 = dr * dr;
        fpair = -par.param.quartic.k / r * (r2 * (ra + rb) + 2.0 * dr * ra * rb);

        if (eflag_either) epair = par.param.quartic.k * r2 * ra * rb;
      }

      if (newton_pair || i < nlocal) {
        f[i].x += dx * fpair;
        f[i].y += dy * fpair;
        f[i].z += dz * fpair;
      }

      if (newton_pair || j < nlocal) {
        f[j].x -= dx * fpair;
        f[j].y -= dy * fpair;
        f[j].z -= dz * fpair;
      }

      if (evflag) ev_tally(i, j, nlocal, newton_pair, epair, 0.0, fpair, dx, dy, dz);
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();

  if (check_flag) {
    int tmp;
    MPI_Allreduce(&pc, &tmp, 1, MPI_INT, MPI_SUM, world);
    if (tmp != 2 * npairs)
      error->all(FLERR, "Not all pairs processed in pair_style list: {} vs {}", tmp, 2 * npairs);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairList::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");
}

/* ----------------------------------------------------------------------
   create one pair style for each arg in list
------------------------------------------------------------------------- */

void PairList::settings(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "pair_style list", error);

  cut_global = utils::numeric(FLERR, arg[1], false, lmp);
  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "nocheck") == 0) {
      check_flag = 0;
      ++iarg;
    } else if (strcmp(arg[2], "check") == 0) {
      check_flag = 1;
      ++iarg;
    } else
      error->all(FLERR, "Unknown pair style list keyword: {}", arg[iarg]);
  }

  std::vector<int> mystyles;
  std::vector<list_param> myparams;

  // read and parse potential file only on MPI rank 0.
  if (comm->me == 0) {
    int nharm, nmorse, nlj126, nquartic, nskipped;
    FILE *fp = utils::open_potential(arg[0], lmp, nullptr);
    if (!fp)
      error->one(FLERR, "Error opening pair list coeffs file {}: {}", arg[0], utils::getsyserror());
    TextFileReader reader(fp, "pair list coeffs");
    npairs = nharm = nmorse = nlj126 = nquartic = nskipped = 0;
    char *line;

    try {
      while ((line = reader.next_line())) {
        ValueTokenizer values(line);
        list_param oneparam;
        oneparam.offset = 0.0;
        oneparam.id1 = values.next_tagint();
        oneparam.id2 = values.next_tagint();
        oneparam.style = stylename[values.next_string()];
        ++npairs;

        switch (oneparam.style) {

          case HARM:
            oneparam.param.harm.k = values.next_double();
            oneparam.param.harm.r0 = values.next_double();
            ++nharm;
            break;

          case MORSE:
            oneparam.param.morse.d0 = values.next_double();
            oneparam.param.morse.alpha = values.next_double();
            oneparam.param.morse.r0 = values.next_double();
            ++nmorse;
            break;

          case LJ126:
            oneparam.param.lj126.epsilon = values.next_double();
            oneparam.param.lj126.sigma = values.next_double();
            ++nlj126;
            break;

          case QUARTIC:
            oneparam.param.quartic.k = values.next_double();
            oneparam.param.quartic.r0 = values.next_double();
            oneparam.param.quartic.b1 = values.next_double();
            oneparam.param.quartic.b2 = values.next_double();
            ++nquartic;
            break;

          case NONE:    // fallthrough
            error->warning(FLERR, "Skipping unrecognized pair list potential entry: {}",
                           utils::trim(line));
            ++nskipped;
            break;
        }
        if (values.has_next())
          oneparam.cutsq = square(values.next_double());
        else
          oneparam.cutsq = cut_global * cut_global;

        myparams.push_back(oneparam);
      }
    } catch (std::exception &e) {
      error->one(FLERR, "Error reading pair list coeffs file: {}\n{}", e.what(), line);
    }
    utils::logmesg(lmp,
                   "Read {} ({}/{}/{}/{}) interacting pair lines from {}. "
                   "{} skipped entries.\n",
                   npairs, nharm, nmorse, nlj126, nquartic, arg[0], nskipped);

    memory->create(params, npairs, "pair_list:params");
    memcpy(params, myparams.data(), npairs * sizeof(list_param));
    fclose(fp);
  }
  MPI_Bcast(&npairs, 1, MPI_INT, 0, world);
  if (comm->me != 0) memory->create(params, npairs, "pair_list:params");
  MPI_Bcast(params, npairs * sizeof(list_param), MPI_BYTE, 0, world);
}

/* ----------------------------------------------------------------------
   there are no coeffs to be set, but we need to update setflag and pretend there are
------------------------------------------------------------------------- */

void PairList::coeff(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "pair_coeff list", error);
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style: compute energy offset at cutoff
------------------------------------------------------------------------- */

void PairList::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style list requires atom IDs");

  if (atom->map_style == Atom::MAP_NONE) error->all(FLERR, "Pair style list requires an atom map");

  if (offset_flag) {
    for (int n = 0; n < npairs; ++n) {
      list_param &par = params[n];

      if (par.style == HARM) {
        const double dr = sqrt(par.cutsq) - par.param.harm.r0;
        par.offset = par.param.harm.k * dr * dr;

      } else if (par.style == MORSE) {
        const double dr = par.param.morse.r0 - sqrt(par.cutsq);
        const double dexp = exp(par.param.morse.alpha * dr);
        par.offset = par.param.morse.d0 * (dexp * dexp - 2.0 * dexp - 1.0);

      } else if (par.style == LJ126) {
        const double r6inv = par.cutsq * par.cutsq * par.cutsq;
        const double sig6 = mypow(par.param.lj126.sigma, 6);
        par.offset = 4.0 * par.param.lj126.epsilon * r6inv * (sig6 * sig6 * r6inv - sig6);

      } else if (par.style == QUARTIC) {
        // the offset is always 0 at rc
        par.offset = 0.0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   since we don't use atom types or neighbor lists, this is a NOP.
------------------------------------------------------------------------- */

double PairList::init_one(int, int)
{
  return cut_global;
}

/* ----------------------------------------------------------------------
   memory usage of each sub-style
------------------------------------------------------------------------- */

double PairList::memory_usage()
{
  double bytes = (double) npairs * sizeof(int);
  bytes += (double) npairs * sizeof(list_param);
  const int n = atom->ntypes + 1;
  bytes += (double) n * (n * sizeof(int) + sizeof(int *));
  bytes += (double) n * (n * sizeof(double) + sizeof(double *));
  return bytes;
}
