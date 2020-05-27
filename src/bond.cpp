/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "bond.h"
#include <mpi.h>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "suffix.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NONE,LINEAR,SPLINE};

/* -----------------------------------------------------------------------
   set bond contribution to Vdwl energy to 0.0
   a particular bond style can override this
------------------------------------------------------------------------- */

Bond::Bond(LAMMPS *lmp) : Pointers(lmp)
{
  energy = 0.0;
  writedata = 1;

  allocated = 0;
  suffix_flag = Suffix::NONE;

  maxeatom = maxvatom = 0;
  eatom = NULL;
  vatom = NULL;
  setflag = NULL;

  execution_space = Host;
  datamask_read = ALL_MASK;
  datamask_modify = ALL_MASK;

  copymode = 0;
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
  if (!allocated && atom->nbondtypes)
    error->all(FLERR,"Bond coeffs are not set");
  for (int i = 1; i <= atom->nbondtypes; i++)
    if (setflag[i] == 0) error->all(FLERR,"All bond coeffs are not set");
  init_style();
}

/* ----------------------------------------------------------------------
   setup for energy, virial computation
   see integrate::ev_set() for values of eflag (0-3) and vflag (0-6)
------------------------------------------------------------------------- */

void Bond::ev_setup(int eflag, int vflag, int alloc)
{
  int i,n;

  evflag = 1;

  eflag_either = eflag;
  eflag_global = eflag % 2;
  eflag_atom = eflag / 2;

  vflag_either = vflag;
  vflag_global = vflag % 4;
  // per-atom virial and per-atom centroid virial are the same for bonds
  vflag_atom = vflag / 4;

  // reallocate per-atom arrays if necessary

  if (eflag_atom && atom->nmax > maxeatom) {
    maxeatom = atom->nmax;
    if (alloc) {
      memory->destroy(eatom);
      memory->create(eatom,comm->nthreads*maxeatom,"bond:eatom");
    }
  }
  if (vflag_atom && atom->nmax > maxvatom) {
    maxvatom = atom->nmax;
    if (alloc) {
      memory->destroy(vatom);
      memory->create(vatom,comm->nthreads*maxvatom,6,"bond:vatom");
    }
  }

  // zero accumulators

  if (eflag_global) energy = 0.0;
  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
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

void Bond::ev_tally(int i, int j, int nlocal, int newton_bond,
                    double ebond, double fbond,
                    double delx, double dely, double delz)
{
  double ebondhalf,v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) energy += ebond;
      else {
        ebondhalf = 0.5*ebond;
        if (i < nlocal) energy += ebondhalf;
        if (j < nlocal) energy += ebondhalf;
      }
    }
    if (eflag_atom) {
      ebondhalf = 0.5*ebond;
      if (newton_bond || i < nlocal) eatom[i] += ebondhalf;
      if (newton_bond || j < nlocal) eatom[j] += ebondhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx*delx*fbond;
    v[1] = dely*dely*fbond;
    v[2] = delz*delz*fbond;
    v[3] = delx*dely*fbond;
    v[4] = delx*delz*fbond;
    v[5] = dely*delz*fbond;

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
          virial[0] += 0.5*v[0];
          virial[1] += 0.5*v[1];
          virial[2] += 0.5*v[2];
          virial[3] += 0.5*v[3];
          virial[4] += 0.5*v[4];
          virial[5] += 0.5*v[5];
        }
        if (j < nlocal) {
          virial[0] += 0.5*v[0];
          virial[1] += 0.5*v[1];
          virial[2] += 0.5*v[2];
          virial[3] += 0.5*v[3];
          virial[4] += 0.5*v[4];
          virial[5] += 0.5*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        vatom[i][0] += 0.5*v[0];
        vatom[i][1] += 0.5*v[1];
        vatom[i][2] += 0.5*v[2];
        vatom[i][3] += 0.5*v[3];
        vatom[i][4] += 0.5*v[4];
        vatom[i][5] += 0.5*v[5];
      }
      if (newton_bond || j < nlocal) {
        vatom[j][0] += 0.5*v[0];
        vatom[j][1] += 0.5*v[1];
        vatom[j][2] += 0.5*v[2];
        vatom[j][3] += 0.5*v[3];
        vatom[j][4] += 0.5*v[4];
        vatom[j][5] += 0.5*v[5];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   write a table of bond potential energy/force vs distance to a file
------------------------------------------------------------------------- */

void Bond::write_file(int narg, char **arg)
{
  if (narg != 6 && narg !=8) error->all(FLERR,"Illegal bond_write command");

  // parse optional arguments

  int itype = 0;
  int jtype = 0;
  if (narg == 8) {
    itype = force->inumeric(FLERR,arg[6]);
    jtype = force->inumeric(FLERR,arg[7]);
    if (itype < 1 || itype > atom->ntypes || jtype < 1 || jtype > atom->ntypes)
    error->all(FLERR,"Invalid atom types in bond_write command");
  }

  int btype = force->inumeric(FLERR,arg[0]);
  int n = force->inumeric(FLERR,arg[1]);
  double inner = force->numeric(FLERR,arg[2]);
  double outer = force->numeric(FLERR,arg[3]);
  if (inner <= 0.0 || inner >= outer)
    error->all(FLERR,"Invalid rlo/rhi values in bond_write command");


  double r0 = equilibrium_distance(btype);

  // open file in append mode
  // print header in format used by bond_style table

  int me;
  MPI_Comm_rank(world,&me);
  FILE *fp;
  if (me == 0) {
    fp = fopen(arg[4],"a");
    if (fp == NULL) error->one(FLERR,"Cannot open bond_write file");
  }

  // initialize potentials before evaluating bond potential
  // insures all bond coeffs are set and force constants
  // also initialize neighbor so that neighbor requests are processed
  // NOTE: might be safest to just do lmp->init()

  force->init();
  neighbor->init();

  if (me == 0) {
    double r,e,f;

    // evaluate energy and force at each of N distances
    // note that Bond::single() takes r**2 and returns f/r.

    fprintf(fp,"# Bond potential %s for bond type %d: i,r,energy,force\n",
            force->bond_style,btype);
    fprintf(fp,"\n%s\nN %d EQ %.15g\n\n",arg[5],n,r0);

    const double dr = (outer-inner) / static_cast<double>(n-1);
    for (int i = 0; i < n; i++) {
      r = inner + dr * static_cast<double>(i);
      e = single(btype,r*r,itype,jtype,f);
      fprintf(fp,"%d %.15g %.15g %.15g\n",i+1,r,e,f*r);
    }
    fclose(fp);
  }
}

/* ---------------------------------------------------------------------- */

double Bond::memory_usage()
{
  double bytes = comm->nthreads*maxeatom * sizeof(double);
  bytes += comm->nthreads*maxvatom*6 * sizeof(double);
  return bytes;
}

/* -----------------------------------------------------------------------
   Reset all type-based bond params via init.
-------------------------------------------------------------------------- */
void Bond::reinit()
{
  if (!reinitflag)
    error->all(FLERR,"Fix adapt interface to this bond style not supported");

  init();
}
