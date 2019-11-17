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
/* IDEAS

-DONE: Need to define a local peratom array for snad and snad on local and ghost atoms
-DONE: Reverse communicate local peratom array
-DONE: Copy peratom array into output array
-DONE: size_array_cols = nperdim (ncoeff [+quadratic])
-DONE: size_array_rows = 1 + total number of atoms + 6
-DONE: size_peratom = (3+6)*nperdim*ntypes
INCOMPLETE: Mappy from local to global
INCOMPLETE: modify->find_compute() 
INCOMPLETE: eliminate local peratom array for viral, replace with fdotr 

 */
#include "compute_snap.h"
#include <cstring>
#include <cstdlib>
#include "sna.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{SCALAR,VECTOR,ARRAY};

ComputeSnap::ComputeSnap(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), cutsq(NULL), list(NULL), snap(NULL),
  radelem(NULL), wjelem(NULL), snap_peratom(NULL)
{

  array_flag = 1;
  extarray = 0;

  double rfac0, rmin0;
  int twojmax, switchflag, bzeroflag;
  radelem = NULL;
  wjelem = NULL;

  int ntypes = atom->ntypes;
  int nargmin = 6+2*ntypes;

  if (narg < nargmin) error->all(FLERR,"Illegal compute snap command");

  // default values

  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  quadraticflag = 0;

  // process required arguments

  memory->create(radelem,ntypes+1,"snap:radelem"); // offset by 1 to match up with types
  memory->create(wjelem,ntypes+1,"snap:wjelem");
  rcutfac = atof(arg[3]);
  rfac0 = atof(arg[4]);
  twojmax = atoi(arg[5]);
  for(int i = 0; i < ntypes; i++)
    radelem[i+1] = atof(arg[6+i]);
  for(int i = 0; i < ntypes; i++)
    wjelem[i+1] = atof(arg[6+ntypes+i]);

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq,ntypes+1,ntypes+1,"snap:cutsq");
  for(int i = 1; i <= ntypes; i++) {
    cut = 2.0*radelem[i]*rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[i][i] = cut*cut;
    for(int j = i+1; j <= ntypes; j++) {
      cut = (radelem[i]+radelem[j])*rcutfac;
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }
  }

  // process optional args

  int iarg = nargmin;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"rmin0") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      rmin0 = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"bzeroflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      bzeroflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"switchflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      switchflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"quadraticflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      quadraticflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal compute snap command");
  }

  snaptr = new SNA(lmp,rfac0,twojmax,
                   rmin0,switchflag,bzeroflag);

  ncoeff = snaptr->ncoeff;
  nperdim = ncoeff;
  if (quadraticflag) nperdim += (ncoeff*(ncoeff+1))/2;
  yoffset = nperdim;
  zoffset = 2*nperdim;
  virialoffset = 3*nperdim;
  natoms = atom->natoms;
  size_array_rows = 1+3*natoms+6;
  size_array_cols = nperdim*atom->ntypes+1; // extra col for reference potential
  ndims_peratom = 9;
  size_peratom = ndims_peratom*nperdim*atom->ntypes; // local atom force and virial data
  comm_reverse = size_peratom;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeSnap::~ComputeSnap()
{
  memory->destroy(snap);
  memory->destroy(snap_peratom);
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);
  delete snaptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSnap::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute snap requires a pair style be defined");

  if (cutmax > force->pair->cutforce)
    error->all(FLERR,"Compute snap cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"snap") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute snap");
  snaptr->init();

  // allocate memory for global array

  //  printf("allocate memory for global array rows = %d cols = %d\n",
  //         size_array_rows,size_array_cols);
  memory->create(snap,size_array_rows,size_array_cols,
                 "snap:snap");
  array = snap;

  // INCOMPLETE: modify->find_compute() 
  // was called 223960 times by snappy Ta example
  // that is over 600 times per config?
  // how is this possible???

  // find compute for reference energy

  char *id_pe = (char *) "thermo_pe";
  int ipe = modify->find_compute(id_pe);
  c_pe = modify->compute[ipe];

  // add compute for reference virial tensor

  char *id_virial = (char *) "snap_press";
  int ivirial = modify->find_compute(id_virial);
  if (ivirial == -1)
    error->all(FLERR,"compute snap requires that compute snap_press exists!");
  c_virial = modify->compute[ivirial];
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

  //  printf("Invoking compute snap on timestep %d\n",invoked_array);

  // grow snap_peratom array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(snap_peratom);
    nmax = atom->nmax;
    memory->create(snap_peratom,nmax,size_peratom,
                   "snap:snap_peratom");
  }

  // clear global array
  // only need to zero out first row

  for (int icoeff = 0; icoeff < size_array_cols; icoeff++)
    snap[0][icoeff] = 0.0;

  // clear peratom array

  for (int i = 0; i < ntotal; i++)
    for (int icoeff = 0; icoeff < size_peratom; icoeff++) {
      snap_peratom[i][icoeff] = 0.0;
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
      const double radi = radelem[itype];
      const int* const jlist = firstneigh[i];
      const int jnum = numneigh[i];
      const int typeoffset = ndims_peratom*nperdim*(atom->type[i]-1);

      // insure rij, inside, and typej  are of size jnum

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
        if (rsq < cutsq[itype][jtype]&&rsq>1e-20) {
          snaptr->rij[ninside][0] = delx;
          snaptr->rij[ninside][1] = dely;
          snaptr->rij[ninside][2] = delz;
          snaptr->inside[ninside] = j;
          snaptr->wj[ninside] = wjelem[jtype];
          snaptr->rcutij[ninside] = (radi+radelem[jtype])*rcutfac;
          ninside++;
        }
      }

      snaptr->compute_ui(ninside);
      snaptr->compute_zi();
      //      if (quadraticflag) {
      snaptr->compute_bi();
      //      }

      for (int jj = 0; jj < ninside; jj++) {
        const int j = snaptr->inside[jj];
        snaptr->compute_duidrj(snaptr->rij[jj],
                                    snaptr->wj[jj],
                                    snaptr->rcutij[jj],jj);
        snaptr->compute_dbidrj();

        // Accumulate -dBi/dRi, -dBi/dRj

        double *snadi = snap_peratom[i]+typeoffset;
        double *snadj = snap_peratom[j]+typeoffset;
        double *snavi = snadi+virialoffset;
        double *snavj = snadj+virialoffset;

        for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
          snadi[icoeff] += snaptr->dblist[icoeff][0];
          snadi[icoeff+yoffset] += snaptr->dblist[icoeff][1];
          snadi[icoeff+zoffset] += snaptr->dblist[icoeff][2];
          snadj[icoeff] -= snaptr->dblist[icoeff][0];
          snadj[icoeff+yoffset] -= snaptr->dblist[icoeff][1];
          snadj[icoeff+zoffset] -= snaptr->dblist[icoeff][2];

          snavi[icoeff]           += snaptr->dblist[icoeff][0]*xtmp;
          snavi[icoeff+nperdim]   += snaptr->dblist[icoeff][1]*ytmp;
          snavi[icoeff+2*nperdim] += snaptr->dblist[icoeff][2]*ztmp;
          snavi[icoeff+3*nperdim] += snaptr->dblist[icoeff][1]*ztmp;
          snavi[icoeff+4*nperdim] += snaptr->dblist[icoeff][0]*ztmp;
          snavi[icoeff+5*nperdim] += snaptr->dblist[icoeff][0]*ytmp;
          snavj[icoeff]           -= snaptr->dblist[icoeff][0]*x[j][0];
          snavj[icoeff+nperdim]   -= snaptr->dblist[icoeff][1]*x[j][1];
          snavj[icoeff+2*nperdim] -= snaptr->dblist[icoeff][2]*x[j][2];
          snavj[icoeff+3*nperdim] -= snaptr->dblist[icoeff][1]*x[j][2];
          snavj[icoeff+4*nperdim] -= snaptr->dblist[icoeff][0]*x[j][2];
          snavj[icoeff+5*nperdim] -= snaptr->dblist[icoeff][0]*x[j][1];
        }

        if (quadraticflag) {
          const int quadraticoffset = ncoeff;
          snadi += quadraticoffset;
          snadj += quadraticoffset;
          snavi += quadraticoffset;
          snavj += quadraticoffset;
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

            snavi[ncount] +=           dbxtmp*xtmp;
            snavi[ncount+nperdim] +=   dbytmp*ytmp;
            snavi[ncount+2*nperdim] += dbztmp*ztmp;
            snavi[ncount+3*nperdim] += dbytmp*ztmp;
            snavi[ncount+4*nperdim] += dbxtmp*ztmp;
            snavi[ncount+5*nperdim] += dbxtmp*ytmp;
            snavj[ncount] -=            dbxtmp*x[j][0];
            snavj[ncount+nperdim] -=    dbytmp*x[j][1];
            snavj[ncount+2*nperdim] -=  dbztmp*x[j][2];
            snavj[ncount+3*nperdim] -=  dbytmp*x[j][2];
            snavj[ncount+4*nperdim] -=  dbxtmp*x[j][2];
            snavj[ncount+5*nperdim] -=  dbxtmp*x[j][1];

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

              snavi[ncount] +=           dbxtmp*xtmp;
              snavi[ncount+nperdim] +=   dbytmp*ytmp;
              snavi[ncount+2*nperdim] += dbztmp*ztmp;
              snavi[ncount+3*nperdim] += dbytmp*ztmp;
              snavi[ncount+4*nperdim] += dbxtmp*ztmp;
              snavi[ncount+5*nperdim] += dbxtmp*ytmp;
              snavj[ncount] -=           dbxtmp*x[j][0];
              snavj[ncount+nperdim] -=   dbytmp*x[j][1];
              snavj[ncount+2*nperdim] -= dbztmp*x[j][2];
              snavj[ncount+3*nperdim] -= dbytmp*x[j][2];
              snavj[ncount+4*nperdim] -= dbxtmp*x[j][2];
              snavj[ncount+5*nperdim] -= dbxtmp*x[j][1];

              ncount++;
            }
          }
        }
      }

      // Accumulate Bi
      
      // linear contributions

      for (int icoeff = 0; icoeff < ncoeff; icoeff++)
        snap[0][icoeff] += snaptr->blist[icoeff];

      // quadratic contributions

      if (quadraticflag) {
        for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
          double bveci = snaptr->blist[icoeff];
          snap[0][icoeff] += 0.5*bveci*bveci;
          for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
            double bvecj = snaptr->blist[jcoeff];
            snap[0][icoeff] += bveci*bvecj;
          }
        }
      }
    }
  }

  // INCOMPLETE
  // can get rid of virial from snap_peratom by doing
  // equivalent of Pair::virial_fdotr_compute()
  // before reverse communicate of snap_peratom

  // communicate snap contributions between neighbor procs

  comm->reverse_comm_compute(this);

  // construct global array

  for (int itype = 0; itype < atom->ntypes; itype++) {
    const int typeoffset = 3*nperdim*itype;
    for (int icoeff = 0; icoeff < nperdim; icoeff++) {

      // assign force rows
      // INCOMPLETE ignore local-global mapping for now

      int irow = 1;
      for (int i = 0; i < atom->nlocal; i++) {
        double *snadi = snap_peratom[i]+typeoffset;
        snap[irow++][icoeff+typeoffset] = snadi[icoeff];
        snap[irow++][icoeff+typeoffset] = snadi[icoeff+yoffset];
        snap[irow++][icoeff+typeoffset] = snadi[icoeff+zoffset];

      }

      // assign virial row

      int irow0 = irow;
      snap[irow++][icoeff+typeoffset] = 0.0; 
      snap[irow++][icoeff+typeoffset] = 0.0; 
      snap[irow++][icoeff+typeoffset] = 0.0; 
      snap[irow++][icoeff+typeoffset] = 0.0; 
      snap[irow++][icoeff+typeoffset] = 0.0; 
      snap[irow++][icoeff+typeoffset] = 0.0;

      for (int i = 0; i < atom->nlocal; i++) {
        double *snavi = snap_peratom[i]+typeoffset+virialoffset;
        irow = irow0;
        snap[irow++][icoeff+typeoffset] += snavi[icoeff];
        snap[irow++][icoeff+typeoffset] += snavi[icoeff+1*nperdim];
        snap[irow++][icoeff+typeoffset] += snavi[icoeff+2*nperdim];
        snap[irow++][icoeff+typeoffset] += snavi[icoeff+3*nperdim];
        snap[irow++][icoeff+typeoffset] += snavi[icoeff+4*nperdim];
        snap[irow++][icoeff+typeoffset] += snavi[icoeff+5*nperdim];
      }

    }
  }

  // assign energy row

  int icol = size_array_cols-1;
  int irow = 0;
  double reference_energy = c_pe->compute_scalar();
  snap[irow++][icol] = reference_energy; 

  // assign force rows
  // INCOMPLETE ignore local-global mapping for now

  for (int i = 0; i < atom->nlocal; i++) {
    snap[irow++][icol] = atom->f[i][0];
    snap[irow++][icol] = atom->f[i][1];
    snap[irow++][icol] = atom->f[i][2];
  }

  // assign virial row
  // switch to Voigt notation

  c_virial->compute_vector();
  snap[irow++][icol] = c_virial->vector[0]; 
  snap[irow++][icol] = c_virial->vector[1]; 
  snap[irow++][icol] = c_virial->vector[2]; 
  snap[irow++][icol] = c_virial->vector[5]; 
  snap[irow++][icol] = c_virial->vector[4]; 
  snap[irow++][icol] = c_virial->vector[3]; 

}

/* ---------------------------------------------------------------------- */

int ComputeSnap::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last,icoeff;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    for (icoeff = 0; icoeff < size_peratom; icoeff++)
      buf[m++] = snap_peratom[i][icoeff];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeSnap::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m,icoeff;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (icoeff = 0; icoeff < size_peratom; icoeff++)
      snap_peratom[j][icoeff] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeSnap::memory_usage()
{

  double bytes = size_array_rows*size_array_cols * 
    sizeof(double);                                     // snap
  bytes += nmax*size_peratom * sizeof(double);          // snap_peratom
  bytes += snaptr->memory_usage();                      // SNA object

  return bytes;
}
