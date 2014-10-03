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

/* ----------------------------------------------------------------------
   Contributing authors: Paul Crozier (SNL), Jeff Greathouse (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "compute_rdf.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeRDF::ComputeRDF(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4 || (narg-4) % 2) error->all(FLERR,"Illegal compute rdf command");

  array_flag = 1;
  extarray = 0;

  nbin = force->inumeric(FLERR,arg[3]);
  if (nbin < 1) error->all(FLERR,"Illegal compute rdf command");
  if (narg == 4) npairs = 1;
  else npairs = (narg-4)/2;

  size_array_rows = nbin;
  size_array_cols = 1 + 2*npairs;

  int ntypes = atom->ntypes;
  memory->create(rdfpair,npairs,ntypes+1,ntypes+1,"rdf:rdfpair");
  memory->create(nrdfpair,ntypes+1,ntypes+1,"rdf:nrdfpair");
  ilo = new int[npairs];
  ihi = new int[npairs];
  jlo = new int[npairs];
  jhi = new int[npairs];

  if (narg == 4) {
    ilo[0] = 1; ihi[0] = ntypes;
    jlo[0] = 1; jhi[0] = ntypes;
    npairs = 1;

  } else {
    npairs = 0;
    int iarg = 4;
    while (iarg < narg) {
      force->bounds(arg[iarg],atom->ntypes,ilo[npairs],ihi[npairs]);
      force->bounds(arg[iarg+1],atom->ntypes,jlo[npairs],jhi[npairs]);
      if (ilo[npairs] > ihi[npairs] || jlo[npairs] > jhi[npairs])
        error->all(FLERR,"Illegal compute rdf command");
      npairs++;
      iarg += 2;
    }
  }

  int i,j;
  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      nrdfpair[i][j] = 0;

  for (int m = 0; m < npairs; m++)
    for (i = ilo[m]; i <= ihi[m]; i++)
      for (j = jlo[m]; j <= jhi[m]; j++)
        rdfpair[nrdfpair[i][j]++][i][j] = m;

  memory->create(hist,npairs,nbin,"rdf:hist");
  memory->create(histall,npairs,nbin,"rdf:histall");
  memory->create(array,nbin,1+2*npairs,"rdf:array");
  typecount = new int[ntypes+1];
  icount = new int[npairs];
  jcount = new int[npairs];
}

/* ---------------------------------------------------------------------- */

ComputeRDF::~ComputeRDF()
{
  memory->destroy(rdfpair);
  memory->destroy(nrdfpair);
  delete [] ilo;
  delete [] ihi;
  delete [] jlo;
  delete [] jhi;
  memory->destroy(hist);
  memory->destroy(histall);
  memory->destroy(array);
  delete [] typecount;
  delete [] icount;
  delete [] jcount;
}

/* ---------------------------------------------------------------------- */

void ComputeRDF::init()
{
  int i,m;

  if (force->pair) delr = force->pair->cutforce / nbin;
  else error->all(FLERR,"Compute rdf requires a pair style be defined");
  delrinv = 1.0/delr;

  // set 1st column of output array to bin coords

  for (int i = 0; i < nbin; i++)
    array[i][0] = (i+0.5) * delr;

  // count atoms of each type that are also in group

  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;

  for (i = 1; i <= ntypes; i++) typecount[i] = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) typecount[type[i]]++;

  // icount = # of I atoms participating in I,J pairs for each histogram
  // jcount = # of J atoms participating in I,J pairs for each histogram

  for (m = 0; m < npairs; m++) {
    icount[m] = 0;
    for (i = ilo[m]; i <= ihi[m]; i++) icount[m] += typecount[i];
    jcount[m] = 0;
    for (i = jlo[m]; i <= jhi[m]; i++) jcount[m] += typecount[i];
  }

  int *scratch = new int[npairs];
  MPI_Allreduce(icount,scratch,npairs,MPI_INT,MPI_SUM,world);
  for (i = 0; i < npairs; i++) icount[i] = scratch[i];
  MPI_Allreduce(jcount,scratch,npairs,MPI_INT,MPI_SUM,world);
  for (i = 0; i < npairs; i++) jcount[i] = scratch[i];
  delete [] scratch;

  // need an occasional half neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeRDF::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeRDF::compute_array()
{
  int i,j,m,ii,jj,inum,jnum,itype,jtype,ipair,jpair,ibin,ihisto;
  double xtmp,ytmp,ztmp,delx,dely,delz,r;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double factor_lj,factor_coul;

  invoked_array = update->ntimestep;

  // invoke half neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero the histogram counts

  for (i = 0; i < npairs; i++)
    for (j = 0; j < nbin; j++)
      hist[i][j] = 0;

  // tally the RDF
  // both atom i and j must be in fix group
  // itype,jtype must have been specified by user
  // consider I,J as one interaction even if neighbor pair is stored on 2 procs
  // tally I,J pair each time I is central atom, and each time J is central

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      // if both weighting factors are 0, skip this pair
      // could be 0 and still be in neigh list for long-range Coulombics
      // want consistency with non-charged pairs which wouldn't be in list

      if (factor_lj == 0.0 && factor_coul == 0.0) continue;

      if (!(mask[j] & groupbit)) continue;
      jtype = type[j];
      ipair = nrdfpair[itype][jtype];
      jpair = nrdfpair[jtype][itype];
      if (!ipair && !jpair) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      r = sqrt(delx*delx + dely*dely + delz*delz);
      ibin = static_cast<int> (r*delrinv);
      if (ibin >= nbin) continue;

      if (ipair)
        for (ihisto = 0; ihisto < ipair; ihisto++)
          hist[rdfpair[ihisto][itype][jtype]][ibin] += 1.0;
      if (newton_pair || j < nlocal) {
        if (jpair)
          for (ihisto = 0; ihisto < jpair; ihisto++)
            hist[rdfpair[ihisto][jtype][itype]][ibin] += 1.0;
      }
    }
  }

  // sum histograms across procs

  MPI_Allreduce(hist[0],histall[0],npairs*nbin,MPI_DOUBLE,MPI_SUM,world);

  // convert counts to g(r) and coord(r) and copy into output array
  // nideal = # of J atoms surrounding single I atom in a single bin
  //   assuming J atoms are at uniform density

  double constant,nideal,gr,ncoord,rlower,rupper;

  if (domain->dimension == 3) {
    constant = 4.0*MY_PI / (3.0*domain->xprd*domain->yprd*domain->zprd);

    for (m = 0; m < npairs; m++) {
      ncoord = 0.0;
      for (ibin = 0; ibin < nbin; ibin++) {
        rlower = ibin*delr;
        rupper = (ibin+1)*delr;
        nideal = constant *
          (rupper*rupper*rupper - rlower*rlower*rlower) * jcount[m];
        if (icount[m]*nideal != 0.0)
          gr = histall[m][ibin] / (icount[m]*nideal);
        else gr = 0.0;
        ncoord += gr*nideal;
        array[ibin][1+2*m] = gr;
        array[ibin][2+2*m] = ncoord;
      }
    }

  } else {
    constant = MY_PI / (domain->xprd*domain->yprd);

    for (m = 0; m < npairs; m++) {
      ncoord = 0.0;
      for (ibin = 0; ibin < nbin; ibin++) {
        rlower = ibin*delr;
        rupper = (ibin+1)*delr;
        nideal = constant * (rupper*rupper - rlower*rlower) * jcount[m];
        if (icount[m]*nideal != 0.0)
          gr = histall[m][ibin] / (icount[m]*nideal);
        else gr = 0.0;
        ncoord += gr*nideal;
        array[ibin][1+2*m] = gr;
        array[ibin][2+2*m] = ncoord;
      }
    }
  }
}
