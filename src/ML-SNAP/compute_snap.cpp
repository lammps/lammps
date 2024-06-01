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

#include "compute_snap.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "sna.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

enum { SCALAR, VECTOR, ARRAY };

ComputeSnap::ComputeSnap(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), cutsq(nullptr), list(nullptr), snap(nullptr), snapall(nullptr),
    snap_peratom(nullptr), radelem(nullptr), wjelem(nullptr), map(nullptr), sinnerelem(nullptr),
    dinnerelem(nullptr), snaptr(nullptr), c_pe(nullptr), c_virial(nullptr)
{
  array_flag = 1;
  extarray = 0;

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
  bikflag = 0;
  dgradflag = 0;
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
    radelem[i + 1] = utils::numeric(FLERR, arg[6 + i], false, lmp);
  for (int i = 0; i < ntypes; i++)
    wjelem[i + 1] = utils::numeric(FLERR, arg[6 + ntypes + i], false, lmp);

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
    } else if (strcmp(arg[iarg], "bikflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute {} command", style);
      bikflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dgradflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Illegal compute snap command");
      dgradflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"switchinnerflag") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Illegal compute snap command");
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
    error->all(FLERR, "Illegal compute {} command: switchinnerflag = 1, "
               "missing sinner/dinner keyword", style);

  if (!switchinnerflag && (sinnerflag || dinnerflag))
    error->all(FLERR, "Illegal compute {} command: switchinnerflag = 0, "
               "unexpected sinner/dinner keyword", style);

  if (dgradflag && !bikflag)
    error->all(FLERR, "Illegal compute snap command: dgradflag=1 requires bikflag=1");

  if (dgradflag && quadraticflag)
    error->all(FLERR, "Illegal compute snap command: dgradflag=1 not implemented for quadratic SNAP");

  snaptr = new SNA(lmp, rfac0, twojmax, rmin0, switchflag, bzeroflag, chemflag, bnormflag,
                   wselfallflag, nelements, switchinnerflag);

  ncoeff = snaptr->ncoeff;
  nvalues = ncoeff;
  if (quadraticflag) nvalues += (ncoeff * (ncoeff + 1)) / 2;

  // end code common to all SNAP computes

  ndims_force = 3;
  ndims_virial = 6;
  yoffset = nvalues;
  zoffset = 2*nvalues;
  natoms = atom->natoms;
  bik_rows = 1;
  if (bikflag) bik_rows = natoms;
  dgrad_rows = ndims_force*natoms;
  size_array_rows = bik_rows+dgrad_rows + ndims_virial;
  if (dgradflag) {
    size_array_rows = bik_rows + 3*natoms*natoms + 1;
    size_array_cols = nvalues + 3;
    if (comm->me == 0)
      error->warning(FLERR, "dgradflag=1 creates a N^2 array, beware of large systems.");
  }
  else size_array_cols = nvalues*atom->ntypes + 1;
  lastcol = size_array_cols-1;

  ndims_peratom = ndims_force;
  size_peratom = ndims_peratom*nvalues*atom->ntypes;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeSnap::~ComputeSnap()
{
  memory->destroy(snap);
  memory->destroy(snapall);
  memory->destroy(snap_peratom);
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

void ComputeSnap::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute snap requires a pair style be defined");

  if (cutmax > force->pair->cutforce)
    error->all(FLERR,"Compute snap cutoff {} is longer than pairwise cutoff {}",
               cutmax, force->pair->cutforce);

  // need an occasional full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if (modify->get_compute_by_style("snap").size() > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute snap");
  snaptr->init();

  // allocate memory for global array

  memory->create(snap,size_array_rows,size_array_cols, "snap:snap");
  memory->create(snapall,size_array_rows,size_array_cols, "snap:snapall");
  array = snapall;

  // find compute for global reference potential energy

  c_pe = modify->get_compute_by_id("thermo_pe");
  if (!c_pe) error->all(FLERR,"compute thermo_pe does not exist.");

  // add compute for global reference virial tensor

  id_virial = id + std::string("_press");
  c_virial = modify->add_compute(id_virial + " all pressure NULL virial");
}


/* ---------------------------------------------------------------------- */

void ComputeSnap::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSnap::compute_array()
{
  int ntotal = atom->nlocal + atom->nghost;

  invoked_array = update->ntimestep;

  // grow snap_peratom array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(snap_peratom);
    nmax = atom->nmax;
    memory->create(snap_peratom,nmax,size_peratom, "snap:snap_peratom");
  }

  // clear global array

  for (int irow = 0; irow < size_array_rows; irow++)
    for (int icoeff = 0; icoeff < size_array_cols; icoeff++)
      snap[irow][icoeff] = 0.0;

  // clear local peratom array

  for (int i = 0; i < ntotal; i++)
    for (int icoeff = 0; icoeff < size_peratom; icoeff++)
      snap_peratom[i][icoeff] = 0.0;

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
    int irow = 0;
    if (bikflag) irow = atom->tag[ilist[ii] & NEIGHMASK]-1;
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
      const int typeoffset_local = ndims_peratom*nvalues*(itype-1);
      const int typeoffset_global = nvalues*(itype-1);

      if (dgradflag) {

        // dBi/dRi tags

        snap[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 0][0] = atom->tag[i]-1;
        snap[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 0][1] = atom->tag[i]-1;
        snap[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 0][2] = 0;
        snap[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 1][0] = atom->tag[i]-1;
        snap[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 1][1] = atom->tag[i]-1;
        snap[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 1][2] = 1;
        snap[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 2][0] = atom->tag[i]-1;
        snap[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 2][1] = atom->tag[i]-1;
        snap[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 2][2] = 2;

        // dBi/dRj tags

        for (int j=0; j<natoms; j++) {
          snap[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 0][0] = atom->tag[i]-1;
          snap[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 0][1] = j;
          snap[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 0][2] = 0;
          snap[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 1][0] = atom->tag[i]-1;
          snap[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 1][1] = j;
          snap[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 1][2] = 1;
          snap[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 2][0] = atom->tag[i]-1;
          snap[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 2][1] = j;
          snap[bik_rows + ((j)*3*natoms) + 3*(atom->tag[i]-1) + 2][2] = 2;
        }
      }

      // ensure rij, inside, and typej are of size jnum

      snaptr->grow_rij(jnum);

      // rij[][3] = displacements between atom I and those neighbors
      // inside = indices of neighbors of I within cutoff
      // typej = types of neighbors of I within cutoff
      // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

      // assign quantities in snaptr

      int ninside=0;
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

      // compute bispectrum for atom i

      snaptr->compute_ui(ninside, ielem);
      snaptr->compute_zi();
      snaptr->compute_bi(ielem);

      // loop over neighbors for descriptors derivatives

      for (int jj = 0; jj < ninside; jj++) {
        const int j = snaptr->inside[jj];

        snaptr->compute_duidrj(jj);
        snaptr->compute_dbidrj();

        // accumulate dBi/dRi, -dBi/dRj

        if (!dgradflag) {

          double *snadi = snap_peratom[i]+typeoffset_local;
          double *snadj = snap_peratom[j]+typeoffset_local;

          for (int icoeff = 0; icoeff < ncoeff; icoeff++) {

            snadi[icoeff] += snaptr->dblist[icoeff][0];
            snadi[icoeff+yoffset] += snaptr->dblist[icoeff][1];
            snadi[icoeff+zoffset] += snaptr->dblist[icoeff][2];

            snadj[icoeff] -= snaptr->dblist[icoeff][0];
            snadj[icoeff+yoffset] -= snaptr->dblist[icoeff][1];
            snadj[icoeff+zoffset] -= snaptr->dblist[icoeff][2];
          }

          if (quadraticflag) {
            const int quadraticoffset = ncoeff;
            snadi += quadraticoffset;
            snadj += quadraticoffset;
            int ncount = 0;
            for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
              double bi = snaptr->blist[icoeff];
              double bix = snaptr->dblist[icoeff][0];
              double biy = snaptr->dblist[icoeff][1];
              double biz = snaptr->dblist[icoeff][2];

              // diagonal elements of quadratic matrix

              double dbxtmp = bi*bix;
              double dbytmp = bi*biy;
              double dbztmp = bi*biz;

              snadi[ncount] +=         dbxtmp;
              snadi[ncount+yoffset] += dbytmp;
              snadi[ncount+zoffset] += dbztmp;
              snadj[ncount] -=         dbxtmp;
              snadj[ncount+yoffset] -= dbytmp;
              snadj[ncount+zoffset] -= dbztmp;

              ncount++;

              // upper-triangular elements of quadratic matrix

              for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                double dbxtmp = bi*snaptr->dblist[jcoeff][0]
                              + bix*snaptr->blist[jcoeff];
                double dbytmp = bi*snaptr->dblist[jcoeff][1]
                              + biy*snaptr->blist[jcoeff];
                double dbztmp = bi*snaptr->dblist[jcoeff][2]
                              + biz*snaptr->blist[jcoeff];

                snadi[ncount] +=         dbxtmp;
                snadi[ncount+yoffset] += dbytmp;
                snadi[ncount+zoffset] += dbztmp;
                snadj[ncount] -=         dbxtmp;
                snadj[ncount+yoffset] -= dbytmp;
                snadj[ncount+zoffset] -= dbztmp;

                ncount++;
              }
            }
          }
        } else {

            for (int icoeff = 0; icoeff < ncoeff; icoeff++) {

              // add to snap array for this proc

              // dBi/dRj

              snap[bik_rows + ((atom->tag[j]-1)*3*natoms) + 3*(atom->tag[i]-1) + 0][icoeff+3] -= snaptr->dblist[icoeff][0];
              snap[bik_rows + ((atom->tag[j]-1)*3*natoms) + 3*(atom->tag[i]-1) + 1][icoeff+3] -= snaptr->dblist[icoeff][1];
              snap[bik_rows + ((atom->tag[j]-1)*3*natoms) + 3*(atom->tag[i]-1) + 2][icoeff+3] -= snaptr->dblist[icoeff][2];

              // dBi/dRi

              snap[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 0][icoeff+3] += snaptr->dblist[icoeff][0];
              snap[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 1][icoeff+3] += snaptr->dblist[icoeff][1];
              snap[bik_rows + ((atom->tag[i]-1)*3*natoms) + 3*(atom->tag[i]-1) + 2][icoeff+3] += snaptr->dblist[icoeff][2];
            }
          }
      } // loop over jj inside

      // accumulate Bi

      if (!dgradflag) {

        // linear contributions

        int k = typeoffset_global;
        for (int icoeff = 0; icoeff < ncoeff; icoeff++)
          snap[irow][k++] += snaptr->blist[icoeff];

        // quadratic contributions

        if (quadraticflag) {
          for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
            double bveci = snaptr->blist[icoeff];
            snap[irow][k++] += 0.5*bveci*bveci;
            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
              double bvecj = snaptr->blist[jcoeff];
              snap[irow][k++] += bveci*bvecj;
            }
          }
        }

      } else {
          int k = 3;
          for (int icoeff = 0; icoeff < ncoeff; icoeff++)
          snap[irow][k++] += snaptr->blist[icoeff];
        }
    }
  } // for (int ii = 0; ii < inum; ii++) {

  // accumulate bispectrum force contributions to global array

  if (!dgradflag) {
    for (int itype = 0; itype < atom->ntypes; itype++) {
      const int typeoffset_local = ndims_peratom*nvalues*itype;
      const int typeoffset_global = nvalues*itype;
      for (int icoeff = 0; icoeff < nvalues; icoeff++) {
        for (int i = 0; i < ntotal; i++) {
          double *snadi = snap_peratom[i]+typeoffset_local;
          int iglobal = atom->tag[i];
          int irow = 3*(iglobal-1)+bik_rows;
          snap[irow++][icoeff+typeoffset_global] += snadi[icoeff];
          snap[irow++][icoeff+typeoffset_global] += snadi[icoeff+yoffset];
          snap[irow][icoeff+typeoffset_global] += snadi[icoeff+zoffset];
        }
      }
    }
  }

 // accumulate forces to global array

  if (!dgradflag) {
    for (int i = 0; i < atom->nlocal; i++) {
      int iglobal = atom->tag[i];
      int irow = 3*(iglobal-1)+bik_rows;
      snap[irow++][lastcol] = atom->f[i][0];
      snap[irow++][lastcol] = atom->f[i][1];
      snap[irow][lastcol] = atom->f[i][2];
    }
  } else {

    // for dgradflag=1, put forces at first 3 columns of bik rows

    for (int i=0; i<atom->nlocal; i++) {
      int iglobal = atom->tag[i];
      snap[iglobal-1][0+0] = atom->f[i][0];
      snap[iglobal-1][0+1] = atom->f[i][1];
      snap[iglobal-1][0+2] = atom->f[i][2];
    }
  }

  // accumulate bispectrum virial contributions to global array

  dbdotr_compute();

  // sum up over all processes

  MPI_Allreduce(&snap[0][0],&snapall[0][0],size_array_rows*size_array_cols,MPI_DOUBLE,MPI_SUM,world);

  // assign energy to last column

  if (!dgradflag) {
    for (int i = 0; i < bik_rows; i++) snapall[i][lastcol] = 0;
    int irow = 0;
    double reference_energy = c_pe->compute_scalar();
    snapall[irow][lastcol] = reference_energy;
  } else {

    // assign reference energy right after the dgrad rows, first column

    int irow = bik_rows + 3*natoms*natoms;
    double reference_energy = c_pe->compute_scalar();
    snapall[irow][0] = reference_energy;
  }

  // assign virial stress to last column
  // switch to Voigt notation

  if (!dgradflag) {
    c_virial->compute_vector();
    int irow = 3*natoms+bik_rows;
    snapall[irow++][lastcol] = c_virial->vector[0];
    snapall[irow++][lastcol] = c_virial->vector[1];
    snapall[irow++][lastcol] = c_virial->vector[2];
    snapall[irow++][lastcol] = c_virial->vector[5];
    snapall[irow++][lastcol] = c_virial->vector[4];
    snapall[irow][lastcol] = c_virial->vector[3];
  }

}

/* ----------------------------------------------------------------------
   compute global virial contributions via summing r_i.dB^j/dr_i over
   own & ghost atoms
------------------------------------------------------------------------- */

void ComputeSnap::dbdotr_compute()
{

  // no virial terms for dgrad yet

  if (dgradflag) return;

  double **x = atom->x;

  int irow0 = bik_rows+ndims_force*natoms;

  // sum over bispectrum contributions to forces
  // on all particles including ghosts

  int nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; i++)
    for (int itype = 0; itype < atom->ntypes; itype++) {
      const int typeoffset_local = ndims_peratom*nvalues*itype;
      const int typeoffset_global = nvalues*itype;
      double *snadi = snap_peratom[i]+typeoffset_local;
      for (int icoeff = 0; icoeff < nvalues; icoeff++) {
        double dbdx = snadi[icoeff];
        double dbdy = snadi[icoeff+yoffset];
        double dbdz = snadi[icoeff+zoffset];
        int irow = irow0;
        snap[irow++][icoeff+typeoffset_global] += dbdx*x[i][0];
        snap[irow++][icoeff+typeoffset_global] += dbdy*x[i][1];
        snap[irow++][icoeff+typeoffset_global] += dbdz*x[i][2];
        snap[irow++][icoeff+typeoffset_global] += dbdz*x[i][1];
        snap[irow++][icoeff+typeoffset_global] += dbdz*x[i][0];
        snap[irow][icoeff+typeoffset_global] += dbdy*x[i][0];
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeSnap::memory_usage()
{

  double bytes = (double)size_array_rows*size_array_cols *
    sizeof(double);                                     // snap
  bytes += (double)size_array_rows*size_array_cols *
    sizeof(double);                                     // snapall
  bytes += (double)nmax*size_peratom * sizeof(double);  // snap_peratom
  bytes += snaptr->memory_usage();                      // SNA object
  int n = atom->ntypes+1;
  bytes += (double)n*sizeof(int);        // map

  return bytes;
}
