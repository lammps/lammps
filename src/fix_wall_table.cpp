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

#include "fix_wall_table.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "table_file_reader.h"
#include "tokenizer.h"

#include <cstring>

using namespace LAMMPS_NS;

enum { NONE, LINEAR, SPLINE };

static constexpr double BIGNUM = 1.0e300;

/* ---------------------------------------------------------------------- */

FixWallTable::FixWallTable(LAMMPS *lmp, int narg, char **arg) :
    FixWall(lmp, narg, arg), tables(nullptr)
{
  dynamic_group_allow = 1;
  tabstyle = NONE;
  if (strcmp(arg[3], "linear") == 0)
    tabstyle = LINEAR;
  else if (strcmp(arg[3], "spline") == 0)
    tabstyle = SPLINE;
  else
    error->all(FLERR, "Unknown table style {} in fix {}", arg[3], style);

  tablength = utils::inumeric(FLERR, arg[4], false, lmp);
  if (tablength < 2) error->all(FLERR, "Illegal number of fix {} table entries", style);
}

/* ---------------------------------------------------------------------- */

void FixWallTable::post_constructor()
{
  memory->sfree(tables);
  tables = (Table *) memory->smalloc(nwall * sizeof(Table), "wall:tables");

  for (int m = 0; m < nwall; ++m) {
    Table &tb = tables[m];
    null_table(tb);

    if (comm->me == 0) read_table(tb, fstr[m], kstr[m]);
    bcast_table(tb);

    // error check on table parameters
    if (tb.ninput <= 1) error->all(FLERR, "Invalid fix {} table length: {}", style, tb.ninput);
    tb.lo = tb.rfile[0];
    tb.hi = tb.rfile[tb.ninput - 1];
    if (tb.lo >= tb.hi) error->all(FLERR, "Fix {} table distance values do not increase", style);
    if (cutoff[m] > tb.hi)
      error->all(FLERR, "Fix {} wall cutoff {} is larger than table outer cutoff {}", style,
                 cutoff[m], tb.hi);

    // spline read-in data and compute r,e,f arrays within table
    spline_table(tb);
    compute_table(tb);
  }
}

/* ---------------------------------------------------------------------- */

FixWallTable::~FixWallTable()
{
  for (int m = 0; m < nwall; m++) free_table(tables[m]);
  memory->sfree(tables);
}

/* ----------------------------------------------------------------------
   compute the potential energy offset so it can be shifted to zero at the cutoff
------------------------------------------------------------------------- */

void FixWallTable::precompute(int m)
{
  double u, mdu;

  // always use linear interpolation to determine the offset, since we may be at the
  // outer table cutoff and spline interpolation would cause out-of-bounds array accesses
  int savedtabstyle = tabstyle;
  tabstyle = LINEAR;
  uf_lookup(m, cutoff[m], u, mdu);
  offset[m] = u;
  tabstyle = savedtabstyle;
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

void FixWallTable::wall_particle(int m, int which, double coord)
{
  double delta, fwall, vn, u, mdu;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = which / 2;
  int side = which % 2;
  if (side == 0) side = -1;

  int onflag = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (side < 0)
        delta = x[i][dim] - coord;
      else
        delta = coord - x[i][dim];
      if (delta >= cutoff[m]) continue;
      if (delta <= 0.0) {
        onflag = 1;
        continue;
      }

      uf_lookup(m, delta, u, mdu);
      fwall = side * mdu;
      f[i][dim] -= fwall;
      ewall[0] += u - offset[m];
      ewall[m + 1] += fwall;

      if (evflag) {
        if (side < 0)
          vn = -fwall * delta;
        else
          vn = fwall * delta;
        v_tally(dim, i, vn);
      }
    }

  if (onflag) error->one(FLERR, "Particle on or inside fix {} surface", style);
}

/* ---------------------------------------------------------------------- */

void FixWallTable::null_table(Table &tb)
{
  tb.rfile = tb.efile = tb.ffile = nullptr;
  tb.e2file = tb.f2file = nullptr;
  tb.r = tb.e = tb.de = nullptr;
  tb.f = tb.df = tb.e2 = tb.f2 = nullptr;
}

/* ---------------------------------------------------------------------- */

void FixWallTable::free_table(Table &tb)
{
  memory->destroy(tb.rfile);
  memory->destroy(tb.efile);
  memory->destroy(tb.ffile);
  memory->destroy(tb.e2file);
  memory->destroy(tb.f2file);

  memory->destroy(tb.r);
  memory->destroy(tb.e);
  memory->destroy(tb.de);
  memory->destroy(tb.f);
  memory->destroy(tb.df);
  memory->destroy(tb.e2);
  memory->destroy(tb.f2);
}

/* ---------------------------------------------------------------------- */

void FixWallTable::read_table(Table &tb, const std::string &tabfile, const std::string &keyword)
{
  TableFileReader reader(lmp, tabfile, "fix wall/table");

  char *line = reader.find_section_start(keyword);
  if (!line) error->one(FLERR, "Did not find keyword {} in table file", keyword);

  // read args on 2nd line of section
  // allocate table arrays for file values

  line = reader.next_line();
  param_extract(tb, line);
  memory->create(tb.rfile, tb.ninput, "wall:rfile");
  memory->create(tb.efile, tb.ninput, "wall:efile");
  memory->create(tb.ffile, tb.ninput, "wall:ffile");

  // read r,e,f table values from file

  reader.skip_line();
  for (int i = 0; i < tb.ninput; i++) {
    line = reader.next_line();
    if (!line)
      error->one(FLERR, "Data missing when parsing wall table '{}' line {} of {}.", keyword, i + 1,
                 tb.ninput);
    try {
      ValueTokenizer values(line);
      values.next_int();
      tb.rfile[i] = values.next_double();
      tb.efile[i] = values.next_double();
      tb.ffile[i] = values.next_double();
    } catch (TokenizerException &e) {
      error->one(FLERR, "Error parsing wall table '{}' line {} of {}. {}\nLine was: {}", keyword,
                 i + 1, tb.ninput, e.what(), line);
    }
  }

  // warn if force != dE/dr at any point that is not an inflection point
  // check via secant approximation to dE/dr
  // skip two end points since do not have surrounding secants
  // inflection point is where curvature changes sign

  double r, e, f, rprev, rnext, eprev, enext, fleft, fright;

  int ferror = 0;
  for (int i = 1; i < tb.ninput - 1; i++) {
    r = tb.rfile[i];
    rprev = tb.rfile[i - 1];
    rnext = tb.rfile[i + 1];
    e = tb.efile[i];
    eprev = tb.efile[i - 1];
    enext = tb.efile[i + 1];
    f = tb.ffile[i];
    fleft = -(e - eprev) / (r - rprev);
    fright = -(enext - e) / (rnext - r);
    if (f < fleft && f < fright) ferror++;
    if (f > fleft && f > fright) ferror++;
  }

  if (ferror)
    error->warning(FLERR,
                   "{} of {} force values in table are inconsistent with -dE/dr.\n"
                   "WARNING:  Should only be flagged at inflection points",
                   ferror, tb.ninput);
}

/* ----------------------------------------------------------------------
   build spline representation of e,f over entire range of read-in table
   this function sets these values in e2file,f2file
------------------------------------------------------------------------- */

void FixWallTable::spline_table(Table &tb)
{
  memory->create(tb.e2file, tb.ninput, "wall:e2file");
  memory->create(tb.f2file, tb.ninput, "wall:f2file");

  double ep0 = -tb.ffile[0];
  double epn = -tb.ffile[tb.ninput - 1];
  spline(tb.rfile, tb.efile, tb.ninput, ep0, epn, tb.e2file);

  if (tb.fpflag == 0) {
    tb.fplo = (tb.ffile[1] - tb.ffile[0]) / (tb.rfile[1] - tb.rfile[0]);
    tb.fphi = (tb.ffile[tb.ninput - 1] - tb.ffile[tb.ninput - 2]) /
        (tb.rfile[tb.ninput - 1] - tb.rfile[tb.ninput - 2]);
  }

  double fp0 = tb.fplo;
  double fpn = tb.fphi;
  spline(tb.rfile, tb.ffile, tb.ninput, fp0, fpn, tb.f2file);
}

/* ----------------------------------------------------------------------
   compute r,e,f vectors from splined values
------------------------------------------------------------------------- */

void FixWallTable::compute_table(Table &tb)
{
  // delta = table spacing for N-1 bins
  int tlm1 = tablength - 1;

  tb.delta = (tb.hi - tb.lo) / tlm1;
  tb.invdelta = 1.0 / tb.delta;
  tb.deltasq6 = tb.delta * tb.delta / 6.0;

  // N-1 evenly spaced bins in r from min to max
  // r,e,f = value at lower edge of bin
  // de,df values = delta values of e,f
  // r,e,f are N in length so de,df arrays can compute difference

  memory->create(tb.r, tablength, "wall:r");
  memory->create(tb.e, tablength, "wall:e");
  memory->create(tb.de, tablength, "wall:de");
  memory->create(tb.f, tablength, "wall:f");
  memory->create(tb.df, tablength, "wall:df");
  memory->create(tb.e2, tablength, "wall:e2");
  memory->create(tb.f2, tablength, "wall:f2");

  double a;
  for (int i = 0; i < tablength; i++) {
    a = tb.lo + i * tb.delta;
    tb.r[i] = a;
    tb.e[i] = splint(tb.rfile, tb.efile, tb.e2file, tb.ninput, a);
    tb.f[i] = splint(tb.rfile, tb.ffile, tb.f2file, tb.ninput, a);
  }

  for (int i = 0; i < tlm1; i++) {
    tb.de[i] = tb.e[i + 1] - tb.e[i];
    tb.df[i] = tb.f[i + 1] - tb.f[i];
  }
  // get final elements from linear extrapolation
  tb.de[tlm1] = 2.0 * tb.de[tlm1 - 1] - tb.de[tlm1 - 2];
  tb.df[tlm1] = 2.0 * tb.df[tlm1 - 1] - tb.df[tlm1 - 2];

  double ep0 = -tb.f[0];
  double epn = -tb.f[tlm1];
  spline(tb.r, tb.e, tablength, ep0, epn, tb.e2);
  spline(tb.r, tb.f, tablength, tb.fplo, tb.fphi, tb.f2);
}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value FP fplo fphi
   N is required, other params are optional
------------------------------------------------------------------------- */

void FixWallTable::param_extract(Table &tb, char *line)
{
  tb.ninput = 0;
  tb.fpflag = 0;

  try {
    ValueTokenizer values(line);

    while (values.has_next()) {
      std::string word = values.next_string();

      if (word == "N") {
        tb.ninput = values.next_int();
      } else if (word == "FP") {
        tb.fpflag = 1;
        tb.fplo = values.next_double();
        tb.fphi = values.next_double();
      } else {
        error->one(FLERR, "Invalid keyword {} in fix {} table parameters", word, style);
      }
    }
  } catch (TokenizerException &e) {
    error->one(FLERR, e.what());
  }

  if (tb.ninput == 0) error->one(FLERR, "Fix {} parameters did not set N", style);
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,rfile,efile,ffile,fpflag,fplo,fphi
------------------------------------------------------------------------- */

void FixWallTable::bcast_table(Table &tb)
{
  MPI_Bcast(&tb.ninput, 1, MPI_INT, 0, world);

  if (comm->me > 0) {
    memory->create(tb.rfile, tb.ninput, "wall:rfile");
    memory->create(tb.efile, tb.ninput, "wall:efile");
    memory->create(tb.ffile, tb.ninput, "wall:ffile");
  }

  MPI_Bcast(tb.rfile, tb.ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb.efile, tb.ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb.ffile, tb.ninput, MPI_DOUBLE, 0, world);

  MPI_Bcast(&tb.fpflag, 1, MPI_INT, 0, world);
  if (tb.fpflag) {
    MPI_Bcast(&tb.fplo, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&tb.fphi, 1, MPI_DOUBLE, 0, world);
  }
}

/* ----------------------------------------------------------------------
   spline and splint routines modified from Numerical Recipes
------------------------------------------------------------------------- */

void FixWallTable::spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
  int i, k;
  double p, qn, sig, un;
  auto u = new double[n];

  if (yp1 > BIGNUM)
    y2[0] = u[0] = 0.0;
  else {
    y2[0] = -0.5;
    u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }
  for (i = 1; i < n - 1; i++) {
    sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
  }
  if (ypn > BIGNUM)
    qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
  }
  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
  for (k = n - 2; k >= 0; k--) y2[k] = y2[k] * y2[k + 1] + u[k];

  delete[] u;
}

/* ---------------------------------------------------------------------- */

double FixWallTable::splint(double *xa, double *ya, double *y2a, int n, double x)
{
  int klo, khi, k;
  double h, b, a, y;

  klo = 0;
  khi = n - 1;
  while (khi - klo > 1) {
    k = (khi + klo) >> 1;
    if (xa[k] > x)
      khi = k;
    else
      klo = k;
  }
  h = xa[khi] - xa[klo];
  a = (xa[khi] - x) / h;
  b = (x - xa[klo]) / h;
  y = a * ya[klo] + b * ya[khi] +
      ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
  return y;
}

/* ----------------------------------------------------------------------
   calculate potential u and force f at distance x
   ensure x is between wall min/max, exit with error if not
------------------------------------------------------------------------- */

void FixWallTable::uf_lookup(int m, double x, double &u, double &f)
{
  double fraction, a, b;
  const Table &tb = tables[m];
  const int itable = static_cast<int>((x - tb.lo) * tb.invdelta);
  if (itable < 0)
    error->one(FLERR, "Particle / wall distance < table inner cutoff: {:.8}", x);
  else if (itable >= tablength)
    error->one(FLERR, "Particle / wall distance > table outer cutoff: {:.8}", x);

  if (tabstyle == LINEAR) {
    fraction = (x - tb.r[itable]) * tb.invdelta;
    u = tb.e[itable] + fraction * tb.de[itable];
    f = tb.f[itable] + fraction * tb.df[itable];
  } else if (tabstyle == SPLINE) {
    fraction = (x - tb.r[itable]) * tb.invdelta;

    b = (x - tb.r[itable]) * tb.invdelta;
    a = 1.0 - b;
    u = a * tb.e[itable] + b * tb.e[itable + 1] +
        ((a * a * a - a) * tb.e2[itable] + (b * b * b - b) * tb.e2[itable + 1]) * tb.deltasq6;
    f = a * tb.f[itable] + b * tb.f[itable + 1] +
        ((a * a * a - a) * tb.f2[itable] + (b * b * b - b) * tb.f2[itable + 1]) * tb.deltasq6;
  }
}
