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
   Contributing author: Xiaowang Zhou (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_eim.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairEIM::PairEIM(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  setfl = NULL;
  nmax = 0;
  rho = NULL;
  fp = NULL;

  nelements = 0;
  elements = NULL;

  negativity = NULL;
  q0 = NULL;
  cutforcesq = NULL;
  Fij = NULL;
  Gij = NULL;
  phiij = NULL;

  Fij_spline = NULL;
  Gij_spline = NULL;
  phiij_spline = NULL;

  // set comm size needed by this Pair

  comm_forward = 1;
  comm_reverse = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairEIM::~PairEIM()
{
  memory->destroy(rho);
  memory->destroy(fp);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
    memory->destroy(type2Fij);
    memory->destroy(type2Gij);
    memory->destroy(type2phiij);
  }

  for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;

  deallocate_setfl();

  delete [] negativity;
  delete [] q0;
  memory->destroy(cutforcesq);
  memory->destroy(Fij);
  memory->destroy(Gij);
  memory->destroy(phiij);

  memory->destroy(Fij_spline);
  memory->destroy(Gij_spline);
  memory->destroy(phiij_spline);
}

/* ---------------------------------------------------------------------- */

void PairEIM::compute(int eflag, int vflag)
{
  int i,j,ii,jj,m,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,phip,phi,coul,coulp,recip,psip;
  double *coeff;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  // grow energy array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(rho,nmax,"pair:rho");
    memory->create(fp,nmax,"pair:fp");
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density

  if (newton_pair) {
    m = nlocal + atom->nghost;
    for (i = 0; i < m; i++) {
      rho[i] = 0.0;
      fp[i] = 0.0;
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      rho[i] = 0.0;
      fp[i] = 0.0;
    }
  }

  // rho = density at each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq[itype][jtype]) {
        p = sqrt(rsq)*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
        coeff = Fij_spline[type2Fij[itype][jtype]][m];
        rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        if (newton_pair || j < nlocal) {
          coeff = Fij_spline[type2Fij[jtype][itype]][m];
          rho[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        }
      }
    }
  }

  // communicate and sum densities

  rhofp = 1;
  if (newton_pair) comm->reverse_comm_pair(this);
  comm->forward_comm_pair(this);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq[itype][jtype]) {
        p = sqrt(rsq)*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
        coeff = Gij_spline[type2Gij[itype][jtype]][m];
        fp[i] += rho[j]*(((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6]);
        if (newton_pair || j < nlocal) {
          fp[j] += rho[i]*(((coeff[3]*p + coeff[4])*p + coeff[5])*p +
                           coeff[6]);
        }
      }
    }
  }

  // communicate and sum modified densities

  rhofp = 2;
  if (newton_pair) comm->reverse_comm_pair(this);
  comm->forward_comm_pair(this);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    if (eflag) {
      phi = 0.5*rho[i]*fp[i];
      if (eflag_global) eng_vdwl += phi;
      if (eflag_atom) eatom[i] += phi;
    }
  }

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq[itype][jtype]) {
        r = sqrt(rsq);
        p = r*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'

        coeff = Fij_spline[type2Fij[jtype][itype]][m];
        rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = Fij_spline[type2Fij[itype][jtype]][m];
        rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = phiij_spline[type2phiij[itype][jtype]][m];
        phip = (coeff[0]*p + coeff[1])*p + coeff[2];
        phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = Gij_spline[type2Gij[itype][jtype]][m];
        coul = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coulp = (coeff[0]*p + coeff[1])*p + coeff[2];
        psip = phip + (rho[i]*rho[j]-q0[itype]*q0[jtype])*coulp +
               fp[i]*rhojp + fp[j]*rhoip;
        recip = 1.0/r;
        fpair = -psip*recip;
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) evdwl = phi-q0[itype]*q0[jtype]*coul;
        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairEIM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
  for (int i = 1; i <= n; i++) map[i] = -1;

  memory->create(type2Fij,n+1,n+1,"pair:type2Fij");
  memory->create(type2Gij,n+1,n+1,"pair:type2Gij");
  memory->create(type2phiij,n+1,n+1,"pair:type2phiij");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEIM::settings(int narg, char **arg)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs from set file
------------------------------------------------------------------------- */

void PairEIM::coeff(int narg, char **arg)
{
  int i,j,m,n;

  if (!allocated) allocate();

  if (narg < 5) error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read EIM element names before filename
  // nelements = # of EIM elements to read from file
  // elements = list of unique element names

  if (nelements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  nelements = narg - 3 - atom->ntypes;
  if (nelements < 1) error->all(FLERR,"Incorrect args for pair coefficients");
  elements = new char*[nelements];

  for (i = 0; i < nelements; i++) {
    n = strlen(arg[i+2]) + 1;
    elements[i] = new char[n];
    strcpy(elements[i],arg[i+2]);
  }

  // read EIM file

  deallocate_setfl();
  setfl = new Setfl();
  read_file(arg[2+nelements]);

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL

  for (i = 3 + nelements; i < narg; i++) {
    m = i - (3+nelements) + 1;
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    if (j < nelements) map[m] = j;
    else if (strcmp(arg[i],"NULL") == 0) map[m] = -1;
    else error->all(FLERR,"Incorrect args for pair coefficients");
  }

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass of atom type if i = j

  int count = 0;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(i,setfl->mass[map[i]]);
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEIM::init_style()
{
  // convert read-in file(s) to arrays and spline them

  file2array();
  array2spline();

  neighbor->request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEIM::init_one(int i, int j)
{
  cutmax = sqrt(cutforcesq[i][j]);
  return cutmax;
}

/* ----------------------------------------------------------------------
   read potential values from a set file
------------------------------------------------------------------------- */

void PairEIM::read_file(char *filename)
{
  // open potential file

  int me = comm->me;
  FILE *fptr;

  if (me == 0) {
    fptr = force->open_potential(filename);
    if (fptr == NULL) {
      char str[128];
      sprintf(str,"Cannot open EIM potential file %s",filename);
      error->one(FLERR,str);
    }
  }

  int npair = nelements*(nelements+1)/2;
  setfl->ielement = new int[nelements];
  setfl->mass = new double[nelements];
  setfl->negativity = new double[nelements];
  setfl->ra = new double[nelements];
  setfl->ri = new double[nelements];
  setfl->Ec = new double[nelements];
  setfl->q0 = new double[nelements];
  setfl->rcutphiA = new double[npair];
  setfl->rcutphiR = new double[npair];
  setfl->Eb = new double[npair];
  setfl->r0 = new double[npair];
  setfl->alpha = new double[npair];
  setfl->beta = new double[npair];
  setfl->rcutq = new double[npair];
  setfl->Asigma = new double[npair];
  setfl->rq = new double[npair];
  setfl->rcutsigma = new double[npair];
  setfl->Ac = new double[npair];
  setfl->zeta = new double[npair];
  setfl->rs = new double[npair];
  setfl->tp = new int[npair];

  if (me == 0)
    if (!grabglobal(fptr))
      error->one(FLERR,"Could not grab global entry from EIM potential file");
  MPI_Bcast(&setfl->division,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&setfl->rbig,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&setfl->rsmall,1,MPI_DOUBLE,0,world);

  for (int i = 0; i < nelements; i++) {
    if (me == 0)
      if (!grabsingle(fptr,i))
        error->one(FLERR,"Could not grab element entry from EIM potential file");
    MPI_Bcast(&setfl->ielement[i],1,MPI_INT,0,world);
    MPI_Bcast(&setfl->mass[i],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&setfl->negativity[i],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&setfl->ra[i],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&setfl->ri[i],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&setfl->Ec[i],1,MPI_DOUBLE,0,world);
    MPI_Bcast(&setfl->q0[i],1,MPI_DOUBLE,0,world);
  }

  for (int i = 0; i < nelements; i++) {
    for (int j = i; j < nelements; j++) {
      int ij;
      if (i == j) ij = i;
      else if (i < j) ij = nelements*(i+1) - (i+1)*(i+2)/2 + j;
      else ij = nelements*(j+1) - (j+1)*(j+2)/2 + i;
      if (me == 0)
        if (grabpair(fptr,i,j) == 0)
          error->one(FLERR,"Could not grab pair entry from EIM potential file");
      MPI_Bcast(&setfl->rcutphiA[ij],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&setfl->rcutphiR[ij],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&setfl->Eb[ij],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&setfl->r0[ij],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&setfl->alpha[ij],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&setfl->beta[ij],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&setfl->rcutq[ij],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&setfl->Asigma[ij],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&setfl->rq[ij],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&setfl->rcutsigma[ij],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&setfl->Ac[ij],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&setfl->zeta[ij],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&setfl->rs[ij],1,MPI_DOUBLE,0,world);
      MPI_Bcast(&setfl->tp[ij],1,MPI_INT,0,world);
    }
  }

  setfl->nr = 5000;
  setfl->cut = 0.0;
  for (int i = 0; i < npair; i++) {
    if (setfl->cut < setfl->rcutphiA[i]) setfl->cut = setfl->rcutphiA[i];
    if (setfl->cut < setfl->rcutphiR[i]) setfl->cut = setfl->rcutphiR[i];
    if (setfl->cut < setfl->rcutq[i]) setfl->cut = setfl->rcutq[i];
    if (setfl->cut < setfl->rcutsigma[i]) setfl->cut = setfl->rcutsigma[i];
  }
  setfl->dr = setfl->cut/(setfl->nr-1.0);

  memory->create(setfl->cuts,nelements,nelements,"pair:cuts");
  for (int i = 0; i < nelements; i++) {
    for (int j = 0; j < nelements; j++) {
      if (i > j) {
        setfl->cuts[i][j] = setfl->cuts[j][i];
      } else {
        int ij;
        if (i == j) {
          ij = i;
        } else {
          ij = nelements*(i+1) - (i+1)*(i+2)/2 + j;
        }
        setfl->cuts[i][j] = setfl->rcutphiA[ij];
        if (setfl->cuts[i][j] < setfl->rcutphiR[ij])
          setfl->cuts[i][j] = setfl->rcutphiR[ij];
        if (setfl->cuts[i][j] < setfl->rcutq[ij])
          setfl->cuts[i][j] = setfl->rcutq[ij];
        if (setfl->cuts[i][j] < setfl->rcutsigma[ij])
          setfl->cuts[i][j] = setfl->rcutsigma[ij];
      }
    }
  }

  memory->create(setfl->Fij,nelements,nelements,setfl->nr+1,"pair:Fij");
  memory->create(setfl->Gij,nelements,nelements,setfl->nr+1,"pair:Gij");
  memory->create(setfl->phiij,nelements,nelements,setfl->nr+1,"pair:phiij");

  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++) {
      for (int k = 0; k < setfl->nr; k++) {
        if (i > j) {
          setfl->phiij[i][j][k+1] = setfl->phiij[j][i][k+1];
        } else {
          double r = k*setfl->dr;
          setfl->phiij[i][j][k+1] = funcphi(i,j,r);
        }
      }
    }

  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++) {
      for (int k = 0; k < setfl->nr; k++) {
        double r = k*setfl->dr;
        setfl->Fij[i][j][k+1] = funcsigma(i,j,r);
      }
    }

  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++) {
      for (int k = 0; k < setfl->nr; k++) {
        if (i > j) {
          setfl->Gij[i][j][k+1] = setfl->Gij[j][i][k+1];
        } else {
          double r = k*setfl->dr;
          setfl->Gij[i][j][k+1] = funccoul(i,j,r);
        }
      }
    }

  // close the potential file

  if (me == 0) fclose(fptr);
}

/* ----------------------------------------------------------------------
   deallocate data associated with setfl file
------------------------------------------------------------------------- */

void PairEIM::deallocate_setfl()
{
  if (!setfl) return;
  delete [] setfl->ielement;
  delete [] setfl->mass;
  delete [] setfl->negativity;
  delete [] setfl->ra;
  delete [] setfl->ri;
  delete [] setfl->Ec;
  delete [] setfl->q0;
  delete [] setfl->rcutphiA;
  delete [] setfl->rcutphiR;
  delete [] setfl->Eb;
  delete [] setfl->r0;
  delete [] setfl->alpha;
  delete [] setfl->beta;
  delete [] setfl->rcutq;
  delete [] setfl->Asigma;
  delete [] setfl->rq;
  delete [] setfl->rcutsigma;
  delete [] setfl->Ac;
  delete [] setfl->zeta;
  delete [] setfl->rs;
  delete [] setfl->tp;
  memory->destroy(setfl->cuts);
  memory->destroy(setfl->Fij);
  memory->destroy(setfl->Gij);
  memory->destroy(setfl->phiij);
  delete setfl;
}

/* ----------------------------------------------------------------------
   convert read-in potentials to standard array format
   interpolate all file values to a single grid and cutoff
------------------------------------------------------------------------- */

void PairEIM::file2array()
{
  int i,j,m,n;
  int irow,icol;
  int ntypes = atom->ntypes;

  delete [] negativity;
  delete [] q0;
  delete [] cutforcesq;
  negativity = new double[ntypes+1];
  q0 = new double[ntypes+1];
  memory->create(cutforcesq,ntypes+1,ntypes+1,"pair:cutforcesq");
  for (i = 1; i <= ntypes; i++) {
    if (map[i] == -1) {
      negativity[i]=0.0;
      q0[i]=0.0;
    } else {
      negativity[i]=setfl->negativity[map[i]];
      q0[i]=setfl->q0[map[i]];
    }
  }

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++) {
      if (map[i] == -1 || map[j] == -1) {
        cutforcesq[i][j] = setfl->cut;
        cutforcesq[i][j] =  cutforcesq[i][j]*cutforcesq[i][j];
      } else {
        cutforcesq[i][j] = setfl->cuts[map[i]][map[j]];
        cutforcesq[i][j] =  cutforcesq[i][j]*cutforcesq[i][j];
      }
    }

  nr = setfl->nr;
  dr = setfl->dr;

  // ------------------------------------------------------------------
  // setup Fij arrays
  // ------------------------------------------------------------------

  nFij = nelements*nelements + 1;
  memory->destroy(Fij);
  memory->create(Fij,nFij,nr+1,"pair:Fij");

  // copy each element's Fij to global Fij

  n=0;
  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++) {
      for (m = 1; m <= nr; m++) Fij[n][m] = setfl->Fij[i][j][m];
      n++;
    }

  // add extra Fij of zeroes for non-EIM types to point to (pair hybrid)

  for (m = 1; m <= nr; m++) Fij[nFij-1][m] = 0.0;

  // type2Fij[i][j] = which Fij array (0 to nFij-1) each type pair maps to
  // setfl of Fij arrays
  // value = n = sum over rows of matrix until reach irow,icol
  // if atom type doesn't point to element (non-EIM atom in pair hybrid)
  // then map it to last Fij array of zeroes

  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2Fij[i][j] = nFij-1;
      } else {
        n = 0;
        for (m = 0; m < irow; m++) n += nelements;
        n += icol;
        type2Fij[i][j] = n;
      }
    }
  }

  // ------------------------------------------------------------------
  // setup Gij arrays
  // ------------------------------------------------------------------

  nGij = nelements * (nelements+1) / 2 + 1;
  memory->destroy(Gij);
  memory->create(Gij,nGij,nr+1,"pair:Gij");

  // copy each element's Gij to global Gij, only for I >= J

  n=0;
  for (i = 0; i < nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= nr; m++) Gij[n][m] = setfl->Gij[i][j][m];
      n++;
    }

  // add extra Gij of zeroes for non-EIM types to point to (pair hybrid)

  for (m = 1; m <= nr; m++) Gij[nGij-1][m] = 0.0;

  // type2Gij[i][j] = which Gij array (0 to nGij-1) each type pair maps to
  // setfl of Gij arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if atom type doesn't point to element (non-EIM atom in pair hybrid)
  // then map it to last Gij array of zeroes

  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2Gij[i][j] = nGij-1;
      } else {
        if (irow < icol) {
          irow = map[j];
          icol = map[i];
        }
        n = 0;
        for (m = 0; m < irow; m++) n += m + 1;
        n += icol;
        type2Gij[i][j] = n;
      }
    }
  }

  // ------------------------------------------------------------------
  // setup phiij arrays
  // ------------------------------------------------------------------

  nphiij = nelements * (nelements+1) / 2 + 1;
  memory->destroy(phiij);
  memory->create(phiij,nphiij,nr+1,"pair:phiij");

  // copy each element pair phiij to global phiij, only for I >= J

  n = 0;
  for (i = 0; i < nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= nr; m++) phiij[n][m] = setfl->phiij[i][j][m];
      n++;
    }

  // add extra phiij of zeroes for non-EIM types to point to (pair hybrid)

  for (m = 1; m <= nr; m++) phiij[nphiij-1][m] = 0.0;

  // type2phiij[i][j] = which phiij array (0 to nphiij-1)
  //                    each type pair maps to
  // setfl of phiij arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if atom type doesn't point to element (non-EIM atom in pair hybrid)
  // then map it to last phiij array of zeroes

  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2phiij[i][j] = nphiij-1;
      } else {
        if (irow < icol) {
          irow = map[j];
          icol = map[i];
        }
        n = 0;
        for (m = 0; m < irow; m++) n += m + 1;
        n += icol;
        type2phiij[i][j] = n;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairEIM::array2spline()
{
  rdr = 1.0/dr;

  memory->destroy(Fij_spline);
  memory->destroy(Gij_spline);
  memory->destroy(phiij_spline);

  memory->create(Fij_spline,nFij,nr+1,7,"pair:Fij");
  memory->create(Gij_spline,nGij,nr+1,7,"pair:Gij");
  memory->create(phiij_spline,nphiij,nr+1,7,"pair:phiij");

  for (int i = 0; i < nFij; i++)
    interpolate(nr,dr,Fij[i],Fij_spline[i],0.0);

  for (int i = 0; i < nGij; i++)
    interpolate(nr,dr,Gij[i],Gij_spline[i],0.0);

  for (int i = 0; i < nphiij; i++)
    interpolate(nr,dr,phiij[i],phiij_spline[i],0.0);
}

/* ---------------------------------------------------------------------- */

void PairEIM::interpolate(int n, double delta, double *f,
                          double **spline, double origin)
{
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
  spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
  spline[n][5] = 0.0;

  for (int m = 3; m <= n-2; m++)
    spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) +
                    8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;

  for (int m = 1; m <= n-1; m++) {
    spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) -
      2.0*spline[m][5] - spline[m+1][5];
    spline[m][3] = spline[m][5] + spline[m+1][5] -
      2.0*(spline[m+1][6]-spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5]/delta;
    spline[m][1] = 2.0*spline[m][4]/delta;
    spline[m][0] = 3.0*spline[m][3]/delta;
  }
}

/* ----------------------------------------------------------------------
   grab global line from file and store info in setfl
   return 0 if error
------------------------------------------------------------------------- */

int PairEIM::grabglobal(FILE *fptr)
{
  char line[MAXLINE];

  char *pch = NULL, *data = NULL;
  while (pch == NULL) {
    if (fgets(line,MAXLINE,fptr) == NULL) break;
    pch = strstr(line,"global");
    if (pch != NULL) {
      data = strtok (line," \t\n\r\f");
      data = strtok (NULL,"?");
      sscanf(data,"%lg %lg %lg",&setfl->division,&setfl->rbig,&setfl->rsmall);
    }
  }
  if (pch == NULL) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   grab elemental line from file and store info in setfl
   return 0 if error
------------------------------------------------------------------------- */

int PairEIM::grabsingle(FILE *fptr, int i)
{
  char line[MAXLINE];

  rewind(fptr);

  char *pch1 = NULL, *pch2 = NULL, *data = NULL;
  while (pch1 == NULL || pch2 == NULL) {
    if (fgets(line,MAXLINE,fptr) == NULL) break;
    pch1 = strtok (line," \t\n\r\f");
    pch1 = strstr(pch1,"element:");
    if (pch1 != NULL) {
      pch2 = strtok(NULL, " \t\n\r\f");
      if (pch2 != NULL) data = strtok (NULL, "?");
      if (strcmp(pch2,elements[i]) == 0) {
        sscanf(data,"%d %lg %lg %lg %lg %lg %lg",&setfl->ielement[i],
          &setfl->mass[i],&setfl->negativity[i],&setfl->ra[i],
          &setfl->ri[i],&setfl->Ec[i],&setfl->q0[i]);
      } else {
        pch2 = NULL;
      }
    }
  }
  if (pch1 == NULL || pch2 == NULL) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   grab pair line from file and store info in setfl
   return 0 if error
------------------------------------------------------------------------- */

int PairEIM::grabpair(FILE *fptr, int i, int j)
{
  char line[MAXLINE];

  rewind(fptr);

  int ij;
  if (i == j) ij = i;
  else if (i < j) ij = nelements*(i+1) - (i+1)*(i+2)/2 + j;
  else ij = nelements*(j+1) - (j+1)*(j+2)/2 + i;

  char *pch1 = NULL, *pch2 = NULL, *pch3 = NULL, *data = NULL;
  while (pch1 == NULL || pch2 == NULL || pch3 == NULL) {
    if (fgets(line,MAXLINE,fptr) == NULL) break;
    pch1 = strtok (line," \t\n\r\f");
    pch1 = strstr(pch1,"pair:");
    if (pch1 != NULL) {
      pch2 = strtok (NULL, " \t\n\r\f");
      if (pch2 != NULL) pch3 = strtok (NULL, " \t\n\r\f");
      if (pch3 != NULL) data = strtok (NULL, "?");
      if ((strcmp(pch2,elements[i]) == 0 &&
        strcmp(pch3,elements[j]) == 0) ||
        (strcmp(pch2,elements[j]) == 0 &&
        strcmp(pch3,elements[i]) == 0)) {
        sscanf(data,"%lg %lg %lg %lg %lg",
          &setfl->rcutphiA[ij],&setfl->rcutphiR[ij],
          &setfl->Eb[ij],&setfl->r0[ij],&setfl->alpha[ij]);
        fgets(line,MAXLINE,fptr);
        sscanf(line,"%lg %lg %lg %lg %lg",
          &setfl->beta[ij],&setfl->rcutq[ij],&setfl->Asigma[ij],
          &setfl->rq[ij],&setfl->rcutsigma[ij]);
        fgets(line,MAXLINE,fptr);
        sscanf(line,"%lg %lg %lg %d",
          &setfl->Ac[ij],&setfl->zeta[ij],&setfl->rs[ij],
          &setfl->tp[ij]);
      } else {
         pch1 = NULL;
         pch2 = NULL;
         pch3 = NULL;
      }
    }
  }
  if (pch1 == NULL || pch2 == NULL || pch3 == NULL) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   cutoff function
------------------------------------------------------------------------- */

double PairEIM::funccutoff(double rp, double rc, double r)
{
  double rbig = setfl->rbig;
  double rsmall = setfl->rsmall;

  double a = (rsmall-rbig)/(rc-rp)*(r-rp)+rbig;
  a = erfc(a);
  double b = erfc(rbig);
  double c = erfc(rsmall);
  return (a-c)/(b-c);
}

/* ----------------------------------------------------------------------
   pair interaction function phi
------------------------------------------------------------------------- */

double PairEIM::funcphi(int i, int j, double r)
{
  int ij;
  double value = 0.0;
  if (i == j) ij = i;
  else if (i < j) ij = nelements*(i+1) - (i+1)*(i+2)/2 + j;
  else ij = nelements*(j+1) - (j+1)*(j+2)/2 + i;
  if (r < 0.2) r = 0.2;
  if (setfl->tp[ij] == 1) {
    double a = setfl->Eb[ij]*setfl->alpha[ij] /
      (setfl->beta[ij]-setfl->alpha[ij]);
    double b = setfl->Eb[ij]*setfl->beta[ij] /
      (setfl->beta[ij]-setfl->alpha[ij]);
    if (r < setfl->rcutphiA[ij]) {
      value -= a*exp(-setfl->beta[ij]*(r/setfl->r0[ij]-1.0))*
        funccutoff(setfl->r0[ij],setfl->rcutphiA[ij],r);
    }
    if (r < setfl-> rcutphiR[ij]) {
      value += b*exp(-setfl->alpha[ij]*(r/setfl->r0[ij]-1.0))*
        funccutoff(setfl->r0[ij],setfl->rcutphiR[ij],r);
    }
  } else if (setfl->tp[ij] == 2) {
    double a=setfl->Eb[ij]*setfl->alpha[ij]*pow(setfl->r0[ij],setfl->beta[ij])/
      (setfl->beta[ij]-setfl->alpha[ij]);
    double b=a*setfl->beta[ij]/setfl->alpha[ij]*
      pow(setfl->r0[ij],setfl->alpha[ij]-setfl->beta[ij]);
    if (r < setfl->rcutphiA[ij]) {
      value -= a/pow(r,setfl->beta[ij])*
        funccutoff(setfl->r0[ij],setfl->rcutphiA[ij],r);
    }
    if (r < setfl-> rcutphiR[ij]) {
      value += b/pow(r,setfl->alpha[ij])*
        funccutoff(setfl->r0[ij],setfl->rcutphiR[ij],r);
    }
  }
  return value;
}

/* ----------------------------------------------------------------------
   ion propensity function sigma
------------------------------------------------------------------------- */

double PairEIM::funcsigma(int i, int j, double r)
{
  int ij;
  double value = 0.0;
  if (i == j) ij = i;
  else if (i < j) ij = nelements*(i+1) - (i+1)*(i+2)/2 + j;
  else ij = nelements*(j+1) - (j+1)*(j+2)/2 + i;
  if (r < 0.2) r = 0.2;
  if (r < setfl->rcutq[ij]) {
    value = setfl->Asigma[ij]*(setfl->negativity[j]-setfl->negativity[i]) *
      funccutoff(setfl->rq[ij],setfl->rcutq[ij],r);
  }
  return value;
}

/* ----------------------------------------------------------------------
   charge-charge interaction function sigma
------------------------------------------------------------------------- */

double PairEIM::funccoul(int i, int j, double r)
{
  int ij;
  double value = 0.0;
  if (i == j) ij = i;
  else if (i < j) ij = nelements*(i+1) - (i+1)*(i+2)/2 + j;
  else ij = nelements*(j+1) - (j+1)*(j+2)/2 + i;
  if (r < 0.2) r = 0.2;
  if (r < setfl->rcutsigma[ij]) {
    value = setfl->Ac[ij]*exp(-setfl->zeta[ij]*r)*
      funccutoff(setfl->rs[ij],setfl->rcutsigma[ij],r);
  }
  return value;
}

/* ---------------------------------------------------------------------- */

int PairEIM::pack_forward_comm(int n, int *list, double *buf, 
                               int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  if (rhofp == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = rho[j];
    }
  }
  if (rhofp == 2) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = fp[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairEIM::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (rhofp == 1) {
    for (i = first; i < last; i++) rho[i] = buf[m++];
  }
  if (rhofp == 2) {
    for (i = first; i < last; i++) fp[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int PairEIM::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (rhofp == 1) {
    for (i = first; i < last; i++) buf[m++] = rho[i];
  }
  if (rhofp == 2) {
    for (i = first; i < last; i++) buf[m++] = fp[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairEIM::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (rhofp == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      rho[j] += buf[m++];
    }
  }
  if (rhofp == 2) {
    for (i = 0; i < n; i++) {
      j = list[i];
      fp[j] += buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairEIM::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += 2 * nmax * sizeof(double);
  return bytes;
}
