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

#include "compute_snav_atom.h"

#include "sna.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

ComputeSNAVAtom::ComputeSNAVAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), cutsq(nullptr), list(nullptr), snav(nullptr),
  radelem(nullptr), wjelem(nullptr), sinnerelem(nullptr), dinnerelem(nullptr)
{

  // begin code common to all SNAP computes

  double rfac0, rmin0;
  int twojmax, switchflag, bzeroflag, bnormflag, wselfallflag;

  int ntypes = atom->ntypes;
  int nargmin = 6 + 2 * ntypes;

  if (narg < nargmin) error->all(FLERR, "Illegal compute {} command", style);

  // default values

  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  quadraticflag = 0;
  chemflag = 0;
  bnormflag = 0;
  wselfallflag = 0;
  switchinnerflag = 0;
  nelements = 1;

  // process required arguments

  memory->create(radelem, ntypes + 1, "sna/atom:radelem"); // offset by 1 to match up with types
  memory->create(wjelem, ntypes + 1, "sna/atom:wjelem");

  rcutfac = utils::numeric(FLERR, arg[3], false, lmp);
  rfac0 = utils::numeric(FLERR, arg[4], false, lmp);
  twojmax = utils::inumeric(FLERR, arg[5], false, lmp);

  for (int i = 0; i < ntypes; i++)
    radelem[i + 1] =
        utils::numeric(FLERR, arg[6 + i], false, lmp);
  for (int i = 0; i < ntypes; i++)
    wjelem[i + 1] =
        utils::numeric(FLERR, arg[6 + ntypes + i], false, lmp);

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq, ntypes + 1, ntypes + 1, "sna/atom:cutsq");
  for (int i = 1; i <= ntypes; i++) {
    cut = 2.0 * radelem[i] * rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[i][i] = cut * cut;
    for (int j = i + 1; j <= ntypes; j++) {
      cut = (radelem[i] + radelem[j]) * rcutfac;
      cutsq[i][j] = cutsq[j][i] = cut * cut;
    }
  }

  // set local input checks

  int sinnerflag = 0;
  int dinnerflag = 0;

  // process optional args

  int iarg = nargmin;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "rmin0") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      rmin0 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "switchflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      switchflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "bzeroflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      bzeroflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "quadraticflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      quadraticflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "chem") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      chemflag = 1;
      memory->create(map, ntypes + 1, "compute_sna_grid:map");
      nelements = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      for (int i = 0; i < ntypes; i++) {
        int jelem = utils::inumeric(FLERR, arg[iarg + 2 + i], false, lmp);
        if (jelem < 0 || jelem >= nelements) error->all(FLERR, "Illegal compute {} command", style);
        map[i + 1] = jelem;
      }
      iarg += 2 + ntypes;
    } else if (strcmp(arg[iarg], "bnormflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      bnormflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "wselfallflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      wselfallflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "switchinnerflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      switchinnerflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "sinner") == 0) {
      iarg++;
      if (iarg + ntypes > narg) error->all(FLERR, "Illegal compute {} command", style);
      memory->create(sinnerelem, ntypes + 1, "snap:sinnerelem");
      for (int i = 0; i < ntypes; i++)
        sinnerelem[i + 1] = utils::numeric(FLERR, arg[iarg + i], false, lmp);
      sinnerflag = 1;
      iarg += ntypes;
    } else if (strcmp(arg[iarg], "dinner") == 0) {
      iarg++;
      if (iarg + ntypes > narg) error->all(FLERR, "Illegal compute {} command", style);
      memory->create(dinnerelem, ntypes + 1, "snap:dinnerelem");
      for (int i = 0; i < ntypes; i++)
        dinnerelem[i + 1] = utils::numeric(FLERR, arg[iarg + i], false, lmp);
      dinnerflag = 1;
      iarg += ntypes;
    } else
      error->all(FLERR, "Illegal compute {} command", style);
  }

  if (switchinnerflag && !(sinnerflag && dinnerflag))
    error->all(
        FLERR,
        "Illegal compute {} command: switchinnerflag = 1, missing sinner/dinner keyword",
        style);

  if (!switchinnerflag && (sinnerflag || dinnerflag))
    error->all(
        FLERR,
        "Illegal compute {} command: switchinnerflag = 0, unexpected sinner/dinner keyword",
        style);

  snaptr = new SNA(lmp, rfac0, twojmax, rmin0, switchflag, bzeroflag, chemflag, bnormflag,
                   wselfallflag, nelements, switchinnerflag);

  ncoeff = snaptr->ncoeff;
  nvalues = ncoeff;
  if (quadraticflag) nvalues += (ncoeff * (ncoeff + 1)) / 2;

  // end code common to all SNAP computes

  size_peratom_cols = 6*nvalues*atom->ntypes;
  comm_reverse = size_peratom_cols;
  peratom_flag = 1;

  nmax = 0;
  snav = nullptr;

}

/* ---------------------------------------------------------------------- */

ComputeSNAVAtom::~ComputeSNAVAtom()
{
  memory->destroy(snav);
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);
  delete snaptr;

  if (chemflag) memory->destroy(map);

  if (switchinnerflag) {
    memory->destroy(sinnerelem);
    memory->destroy(dinnerelem);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSNAVAtom::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute snav/atom requires a pair style be defined");

  if (cutmax > force->pair->cutforce)
    error->all(FLERR,"Compute snav/atom cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if (modify->get_compute_by_style("snav/atom").size() > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute snav/atom");
  snaptr->init();
}

/* ---------------------------------------------------------------------- */

void ComputeSNAVAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSNAVAtom::compute_peratom()
{
  int ntotal = atom->nlocal + atom->nghost;

  invoked_peratom = update->ntimestep;

  // grow snav array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(snav);
    nmax = atom->nmax;
    memory->create(snav,nmax,size_peratom_cols,
                   "snav/atom:snav");
    array_atom = snav;
  }

  // clear local array

  for (int i = 0; i < ntotal; i++)
    for (int icoeff = 0; icoeff < size_peratom_cols; icoeff++) {
      snav[i][icoeff] = 0.0;
    }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int* const ilist = list->ilist;
  const int* const numneigh = list->numneigh;
  int** const firstneigh = list->firstneigh;
  int * const type = atom->type;
  // compute sna derivatives for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double** const x = atom->x;
  const int* const mask = atom->mask;

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {

      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const int itype = type[i];
      int ielem = 0;
      if (chemflag)
        ielem = map[itype];
      const double radi = radelem[itype];

      const int* const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      const int typeoffset = 6*nvalues*(atom->type[i]-1);

      // ensure rij, inside, and typej  are of size jnum

      snaptr->grow_rij(jnum);

      // rij[][3] = displacements between atom I and those neighbors
      // inside = indices of neighbors of I within cutoff
      // typej = types of neighbors of I within cutoff
      // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

      int ninside = 0;
      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        j &= NEIGHMASK;

        const double delx = x[j][0] - xtmp;
        const double dely = x[j][1] - ytmp;
        const double delz = x[j][2] - ztmp;
        const double rsq = delx*delx + dely*dely + delz*delz;
        int jtype = type[j];
        int jelem = 0;
        if (chemflag)
          jelem = map[jtype];
        if (rsq < cutsq[itype][jtype]&&rsq>1e-20) {
          snaptr->rij[ninside][0] = delx;
          snaptr->rij[ninside][1] = dely;
          snaptr->rij[ninside][2] = delz;
          snaptr->inside[ninside] = j;
          snaptr->wj[ninside] = wjelem[jtype];
          snaptr->rcutij[ninside] = (radi+radelem[jtype])*rcutfac;
          if (switchinnerflag) {
            snaptr->sinnerij[ninside] = 0.5*(sinnerelem[itype]+sinnerelem[jtype]);
            snaptr->dinnerij[ninside] = 0.5*(dinnerelem[itype]+dinnerelem[jtype]);
          }
          if (chemflag) snaptr->element[ninside] = jelem;
          ninside++;
        }
      }

      snaptr->compute_ui(ninside, ielem);
      snaptr->compute_zi();
      if (quadraticflag) {
        snaptr->compute_bi(ielem);
      }

      for (int jj = 0; jj < ninside; jj++) {
        const int j = snaptr->inside[jj];

        snaptr->compute_duidrj(jj);
        snaptr->compute_dbidrj();

        // Accumulate -dBi/dRi*Ri, -dBi/dRj*Rj

        double *snavi = snav[i]+typeoffset;
        double *snavj = snav[j]+typeoffset;

        for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
          snavi[icoeff]           += snaptr->dblist[icoeff][0]*xtmp;
          snavi[icoeff+nvalues]   += snaptr->dblist[icoeff][1]*ytmp;
          snavi[icoeff+2*nvalues] += snaptr->dblist[icoeff][2]*ztmp;
          snavi[icoeff+3*nvalues] += snaptr->dblist[icoeff][1]*ztmp;
          snavi[icoeff+4*nvalues] += snaptr->dblist[icoeff][0]*ztmp;
          snavi[icoeff+5*nvalues] += snaptr->dblist[icoeff][0]*ytmp;
          snavj[icoeff]           -= snaptr->dblist[icoeff][0]*x[j][0];
          snavj[icoeff+nvalues]   -= snaptr->dblist[icoeff][1]*x[j][1];
          snavj[icoeff+2*nvalues] -= snaptr->dblist[icoeff][2]*x[j][2];
          snavj[icoeff+3*nvalues] -= snaptr->dblist[icoeff][1]*x[j][2];
          snavj[icoeff+4*nvalues] -= snaptr->dblist[icoeff][0]*x[j][2];
          snavj[icoeff+5*nvalues] -= snaptr->dblist[icoeff][0]*x[j][1];
        }

        if (quadraticflag) {
          const int quadraticoffset = ncoeff;
          snavi += quadraticoffset;
          snavj += quadraticoffset;
          int ncount = 0;
          for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
            double bi = snaptr->blist[icoeff];
            double bix = snaptr->dblist[icoeff][0];
            double biy = snaptr->dblist[icoeff][1];
            double biz = snaptr->dblist[icoeff][2];

            // diagonal element of quadratic matrix

            double dbxtmp = bi*bix;
            double dbytmp = bi*biy;
            double dbztmp = bi*biz;
            snavi[ncount] +=           dbxtmp*xtmp;
            snavi[ncount+nvalues] +=   dbytmp*ytmp;
            snavi[ncount+2*nvalues] += dbztmp*ztmp;
            snavi[ncount+3*nvalues] += dbytmp*ztmp;
            snavi[ncount+4*nvalues] += dbxtmp*ztmp;
            snavi[ncount+5*nvalues] += dbxtmp*ytmp;
            snavj[ncount] -=            dbxtmp*x[j][0];
            snavj[ncount+nvalues] -=    dbytmp*x[j][1];
            snavj[ncount+2*nvalues] -=  dbztmp*x[j][2];
            snavj[ncount+3*nvalues] -=  dbytmp*x[j][2];
            snavj[ncount+4*nvalues] -=  dbxtmp*x[j][2];
            snavj[ncount+5*nvalues] -=  dbxtmp*x[j][1];
            ncount++;

            // upper-triangular elements of quadratic matrix

            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
              double dbxtmp = bi*snaptr->dblist[jcoeff][0]
                + bix*snaptr->blist[jcoeff];
              double dbytmp = bi*snaptr->dblist[jcoeff][1]
                + biy*snaptr->blist[jcoeff];
              double dbztmp = bi*snaptr->dblist[jcoeff][2]
                + biz*snaptr->blist[jcoeff];
              snavi[ncount] +=           dbxtmp*xtmp;
              snavi[ncount+nvalues] +=   dbytmp*ytmp;
              snavi[ncount+2*nvalues] += dbztmp*ztmp;
              snavi[ncount+3*nvalues] += dbytmp*ztmp;
              snavi[ncount+4*nvalues] += dbxtmp*ztmp;
              snavi[ncount+5*nvalues] += dbxtmp*ytmp;
              snavj[ncount] -=           dbxtmp*x[j][0];
              snavj[ncount+nvalues] -=   dbytmp*x[j][1];
              snavj[ncount+2*nvalues] -= dbztmp*x[j][2];
              snavj[ncount+3*nvalues] -= dbytmp*x[j][2];
              snavj[ncount+4*nvalues] -= dbxtmp*x[j][2];
              snavj[ncount+5*nvalues] -= dbxtmp*x[j][1];
              ncount++;
            }
          }
        }
      }
    }
  }

  // communicate snav contributions between neighbor procs

  comm->reverse_comm(this);

}

/* ---------------------------------------------------------------------- */

int ComputeSNAVAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last,icoeff;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    for (icoeff = 0; icoeff < size_peratom_cols; icoeff++)
      buf[m++] = snav[i][icoeff];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeSNAVAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m,icoeff;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (icoeff = 0; icoeff < size_peratom_cols; icoeff++)
      snav[j][icoeff] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeSNAVAtom::memory_usage()
{
  double bytes = (double)nmax*size_peratom_cols * sizeof(double); // snav
  bytes += snaptr->memory_usage();                        // SNA object

  return bytes;
}
