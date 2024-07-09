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

#include "bond.h"

#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "suffix.h"
#include "update.h"

using namespace LAMMPS_NS;

enum { NONE, LINEAR, SPLINE };

// allocate space for static class instance variable and initialize it

int Bond::instance_total = 0;

/* -----------------------------------------------------------------------
   set bond contribution to Vdwl energy to 0.0
   a particular bond style can override this
------------------------------------------------------------------------- */

Bond::Bond(LAMMPS *_lmp) : Pointers(_lmp)
{
  instance_me = instance_total++;

  energy = 0.0;
  virial[0] = virial[1] = virial[2] = virial[3] = virial[4] = virial[5] = 0.0;
  writedata = 1;
  reinitflag = 1;

  comm_forward = comm_reverse = comm_reverse_off = 0;

  allocated = 0;
  suffix_flag = Suffix::NONE;
  born_matrix_enable = 0;
  partial_flag = 0;

  single_extra = 0;
  svector = nullptr;

  maxeatom = maxvatom = 0;
  eatom = nullptr;
  vatom = nullptr;
  setflag = nullptr;

  execution_space = Host;
  datamask_read = ALL_MASK;
  datamask_modify = ALL_MASK;
  copymode = kokkosable = 0;
}

/* ---------------------------------------------------------------------- */

Bond::~Bond()
{
  if (copymode) return;

  memory->destroy(eatom);
  memory->destroy(vatom);
}

/* ----------------------------------------------------------------------
   check if all coeffs are set
------------------------------------------------------------------------- */

void Bond::init()
{
  if (!allocated && atom->nbondtypes) error->all(FLERR, "Bond coeffs are not set");
  for (int i = 1; i <= atom->nbondtypes; i++)
    if (setflag[i] == 0) error->all(FLERR, "All bond coeffs are not set");
  init_style();
}

/* ----------------------------------------------------------------------
   check that there are no arguments
------------------------------------------------------------------------- */

void Bond::settings(int narg, char **args)
{
  if (narg > 0) error->all(FLERR, "Illegal bond_style {} argument: {}", force->bond_style, args[0]);
}

/* ----------------------------------------------------------------------
   setup for energy, virial computation
   see integrate::ev_set() for bitwise settings of eflag/vflag
   set the following flags, values are otherwise set to 0:
     evflag       != 0 if any bits of eflag or vflag are set
     eflag_global != 0 if ENERGY_GLOBAL bit of eflag set
     eflag_atom   != 0 if ENERGY_ATOM bit of eflag set
     eflag_either != 0 if eflag_global or eflag_atom is set
     vflag_global != 0 if VIRIAL_PAIR or VIRIAL_FDOTR bit of vflag set
     vflag_atom   != 0 if VIRIAL_ATOM or VIRIAL_CENTROID bit of vflag set
                       two-body and centroid stress are identical for bonds
     vflag_either != 0 if vflag_global or vflag_atom is set
------------------------------------------------------------------------- */

void Bond::ev_setup(int eflag, int vflag, int alloc)
{
  int i, n;

  evflag = 1;

  eflag_either = eflag;
  eflag_global = eflag & ENERGY_GLOBAL;
  eflag_atom = eflag & ENERGY_ATOM;

  vflag_either = vflag;
  vflag_global = vflag & (VIRIAL_PAIR | VIRIAL_FDOTR);
  vflag_atom = vflag & (VIRIAL_ATOM | VIRIAL_CENTROID);

  // reallocate per-atom arrays if necessary

  if (eflag_atom && atom->nmax > maxeatom) {
    maxeatom = atom->nmax;
    if (alloc) {
      memory->destroy(eatom);
      memory->create(eatom, comm->nthreads * maxeatom, "bond:eatom");
    }
  }
  if (vflag_atom && atom->nmax > maxvatom) {
    maxvatom = atom->nmax;
    if (alloc) {
      memory->destroy(vatom);
      memory->create(vatom, comm->nthreads * maxvatom, 6, "bond:vatom");
    }
  }

  // zero accumulators

  if (eflag_global) energy = 0.0;
  if (vflag_global)
    for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (eflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton_bond) n += atom->nghost;
    for (i = 0; i < n; i++) eatom[i] = 0.0;
  }
  if (vflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton_bond) n += atom->nghost;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
------------------------------------------------------------------------- */

void Bond::ev_tally(int i, int j, int nlocal, int newton_bond, double ebond, double fbond,
                    double delx, double dely, double delz)
{
  double ebondhalf, v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond)
        energy += ebond;
      else {
        ebondhalf = 0.5 * ebond;
        if (i < nlocal) energy += ebondhalf;
        if (j < nlocal) energy += ebondhalf;
      }
    }
    if (eflag_atom) {
      ebondhalf = 0.5 * ebond;
      if (newton_bond || i < nlocal) eatom[i] += ebondhalf;
      if (newton_bond || j < nlocal) eatom[j] += ebondhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx * delx * fbond;
    v[1] = dely * dely * fbond;
    v[2] = delz * delz * fbond;
    v[3] = delx * dely * fbond;
    v[4] = delx * delz * fbond;
    v[5] = dely * delz * fbond;

    if (vflag_global) {
      if (newton_bond) {
        virial[0] += v[0];
        virial[1] += v[1];
        virial[2] += v[2];
        virial[3] += v[3];
        virial[4] += v[4];
        virial[5] += v[5];
      } else {
        if (i < nlocal) {
          virial[0] += 0.5 * v[0];
          virial[1] += 0.5 * v[1];
          virial[2] += 0.5 * v[2];
          virial[3] += 0.5 * v[3];
          virial[4] += 0.5 * v[4];
          virial[5] += 0.5 * v[5];
        }
        if (j < nlocal) {
          virial[0] += 0.5 * v[0];
          virial[1] += 0.5 * v[1];
          virial[2] += 0.5 * v[2];
          virial[3] += 0.5 * v[3];
          virial[4] += 0.5 * v[4];
          virial[5] += 0.5 * v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        vatom[i][0] += 0.5 * v[0];
        vatom[i][1] += 0.5 * v[1];
        vatom[i][2] += 0.5 * v[2];
        vatom[i][3] += 0.5 * v[3];
        vatom[i][4] += 0.5 * v[4];
        vatom[i][5] += 0.5 * v[5];
      }
      if (newton_bond || j < nlocal) {
        vatom[j][0] += 0.5 * v[0];
        vatom[j][1] += 0.5 * v[1];
        vatom[j][2] += 0.5 * v[2];
        vatom[j][3] += 0.5 * v[3];
        vatom[j][4] += 0.5 * v[4];
        vatom[j][5] += 0.5 * v[5];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally energy and virial into global or per-atom accumulators
   for virial, have delx,dely,delz and fx,fy,fz
------------------------------------------------------------------------- */

void Bond::ev_tally_xyz(int i, int j, int nlocal, int newton_bond, double ebond, double fx,
                        double fy, double fz, double delx, double dely, double delz)
{
  double ebondhalf, v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) {
        energy += ebond;
      } else {
        ebondhalf = 0.5 * ebond;
        if (i < nlocal) energy += ebondhalf;
        if (j < nlocal) energy += ebondhalf;
      }
    }
    if (eflag_atom) {
      ebondhalf = 0.5 * ebond;
      if (newton_bond || i < nlocal) eatom[i] += ebondhalf;
      if (newton_bond || j < nlocal) eatom[j] += ebondhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx * fx;
    v[1] = dely * fy;
    v[2] = delz * fz;
    v[3] = delx * fy;
    v[4] = delx * fz;
    v[5] = dely * fz;

    if (vflag_global) {
      if (newton_bond) {
        virial[0] += v[0];
        virial[1] += v[1];
        virial[2] += v[2];
        virial[3] += v[3];
        virial[4] += v[4];
        virial[5] += v[5];
      } else {
        if (i < nlocal) {
          virial[0] += 0.5 * v[0];
          virial[1] += 0.5 * v[1];
          virial[2] += 0.5 * v[2];
          virial[3] += 0.5 * v[3];
          virial[4] += 0.5 * v[4];
          virial[5] += 0.5 * v[5];
        }
        if (j < nlocal) {
          virial[0] += 0.5 * v[0];
          virial[1] += 0.5 * v[1];
          virial[2] += 0.5 * v[2];
          virial[3] += 0.5 * v[3];
          virial[4] += 0.5 * v[4];
          virial[5] += 0.5 * v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        vatom[i][0] += 0.5 * v[0];
        vatom[i][1] += 0.5 * v[1];
        vatom[i][2] += 0.5 * v[2];
        vatom[i][3] += 0.5 * v[3];
        vatom[i][4] += 0.5 * v[4];
        vatom[i][5] += 0.5 * v[5];
      }
      if (newton_bond || j < nlocal) {
        vatom[j][0] += 0.5 * v[0];
        vatom[j][1] += 0.5 * v[1];
        vatom[j][2] += 0.5 * v[2];
        vatom[j][3] += 0.5 * v[3];
        vatom[j][4] += 0.5 * v[4];
        vatom[j][5] += 0.5 * v[5];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   write a table of bond potential energy/force vs distance to a file
------------------------------------------------------------------------- */

void Bond::write_file(int narg, char **arg)
{
  if (narg != 6 && narg != 8) error->all(FLERR, "Illegal bond_write command");

  // parse optional arguments

  int itype = 0;
  int jtype = 0;
  if (narg == 8) {
    itype = utils::inumeric(FLERR, arg[6], false, lmp);
    jtype = utils::inumeric(FLERR, arg[7], false, lmp);
    if ((itype < 1) || (itype > atom->ntypes))
      error->all(FLERR, "Invalid atom type {} in bond_write command", itype);
    if ((jtype < 1) || (jtype > atom->ntypes))
      error->all(FLERR, "Invalid atom type {} in bond_write command", jtype);
  }

  int btype = utils::inumeric(FLERR, arg[0], false, lmp);
  int n = utils::inumeric(FLERR, arg[1], false, lmp);
  double inner = utils::numeric(FLERR, arg[2], false, lmp);
  double outer = utils::numeric(FLERR, arg[3], false, lmp);
  if ((inner <= 0.0) || (inner >= outer))
    error->all(FLERR, "Invalid rlo={} / rhi={} values in bond_write command.", inner, outer);

  if (n < 2) error->all(FLERR, "Must have at least 2 table values");
  double r0 = equilibrium_distance(btype);

  // open file in append mode if exists
  // add line with DATE: and UNITS: tag when creating new file
  // print header in format used by bond_style table

  FILE *fp = nullptr;
  if (comm->me == 0) {
    std::string table_file = arg[4];

    // units sanity check:
    // - if this is the first time we write to this potential file,
    //   write out a line with "DATE:" and "UNITS:" tags
    // - if the file already exists, print a message about appending
    //   while printing the date and check that units are consistent.
    if (platform::file_is_readable(table_file)) {
      std::string units = utils::get_potential_units(table_file, "table");
      if (!units.empty() && (units != update->unit_style)) {
        error->one(FLERR, "Trying to append to a table file with UNITS: {} while units are {}",
                   units, update->unit_style);
      }
      std::string date = utils::get_potential_date(table_file, "table");
      utils::logmesg(lmp, "Appending to table file {} with DATE: {}\n", table_file, date);
      fp = fopen(table_file.c_str(), "a");
    } else {
      utils::logmesg(lmp, "Creating table file {} with DATE: {}\n", table_file,
                     utils::current_date());
      fp = fopen(table_file.c_str(), "w");
      if (fp)
        fmt::print(fp, "# DATE: {} UNITS: {} Created by bond_write\n", utils::current_date(),
                   update->unit_style);
    }
    if (fp == nullptr)
      error->one(FLERR, "Cannot open bond_write file {}: {}", arg[4], utils::getsyserror());
  }

  // initialize potentials before evaluating bond potential
  // ensures all bond coeffs are set and force constants
  // also initialize neighbor so that neighbor requests are processed
  // NOTE: might be safest to just do lmp->init()

  force->init();
  neighbor->init();

  if (comm->me == 0) {
    double r, e, f;

    // evaluate energy and force at each of N distances
    // note that Bond::single() takes r**2 and returns f/r.

    fprintf(fp, "# Bond potential %s for bond type %d: i,r,energy,force\n", force->bond_style,
            btype);
    fprintf(fp, "\n%s\nN %d EQ %.15g\n\n", arg[5], n, r0);

    const double dr = (outer - inner) / static_cast<double>(n - 1);
    for (int i = 0; i < n; i++) {
      r = inner + dr * static_cast<double>(i);
      e = single(btype, r * r, itype, jtype, f);
      fprintf(fp, "%8d %- 22.15g %- 22.15g %- 22.15g\n", i + 1, r, e, f * r);
    }
    fclose(fp);
  }
}

/* ---------------------------------------------------------------------- */

double Bond::memory_usage()
{
  double bytes = (double) comm->nthreads * maxeatom * sizeof(double);
  bytes += (double) comm->nthreads * maxvatom * 6 * sizeof(double);
  return bytes;
}

/* -----------------------------------------------------------------------
   reset all type-based bond params via init()
-------------------------------------------------------------------------- */

void Bond::reinit()
{
  if (!reinitflag) error->all(FLERR, "Fix adapt interface to this bond style not supported");

  init();
}
