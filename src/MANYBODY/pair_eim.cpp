// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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

#include "pair_eim.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairEIM::PairEIM(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  unit_convert_flag = utils::get_supported_conversions(utils::ENERGY);

  setfl = nullptr;
  nmax = 0;
  rho = nullptr;
  fp = nullptr;

  negativity = nullptr;
  q0 = nullptr;
  cutforcesq = nullptr;
  Fij = nullptr;
  Gij = nullptr;
  phiij = nullptr;

  Fij_spline = nullptr;
  Gij_spline = nullptr;
  phiij_spline = nullptr;

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
    memory->destroy(type2Fij);
    memory->destroy(type2Gij);
    memory->destroy(type2phiij);
  }

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
  ev_init(eflag,vflag);

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

void PairEIM::settings(int narg, char **/*arg*/)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs from set file
------------------------------------------------------------------------- */

void PairEIM::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg < 5) error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  const int ntypes = atom->ntypes;
  map_element2type(ntypes,arg+(narg-ntypes));

  // read EIM file

  deallocate_setfl();
  setfl = new Setfl();
  read_file(arg[2+nelements]);

  // set per-type atomic masses

  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++)
      if ((map[i] >= 0) && (map[j] >= 0))
        if (i == j) atom->set_mass(FLERR,i,setfl->mass[map[i]]);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEIM::init_style()
{
  // convert read-in file(s) to arrays and spline them

  file2array();
  array2spline();

  neighbor->request(this,instance_me);
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

  // read potential file
  if ( comm->me == 0) {
    EIMPotentialFileReader reader(lmp, filename, unit_convert_flag);

    reader.get_global(setfl);

    for (int i = 0; i < nelements; i++) {
      reader.get_element(setfl, i, elements[i]);
    }

    for (int i = 0; i < nelements; i++) {
      for (int j = i; j < nelements; j++) {
        int ij;
        if (i == j) ij = i;
        else if (i < j) ij = nelements*(i+1) - (i+1)*(i+2)/2 + j;
        else ij = nelements*(j+1) - (j+1)*(j+2)/2 + i;
        reader.get_pair(setfl, ij, elements[i], elements[j]);
      }
    }
  }

  // broadcast potential information to other procs
  MPI_Bcast(&setfl->division, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&setfl->rbig, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&setfl->rsmall, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(setfl->ielement, nelements, MPI_INT, 0, world);
  MPI_Bcast(setfl->mass, nelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->negativity, nelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->ra, nelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->ri, nelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->Ec, nelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->q0, nelements, MPI_DOUBLE, 0, world);

  MPI_Bcast(setfl->rcutphiA, npair, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->rcutphiR, npair, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->Eb, npair, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->r0, npair, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->alpha, npair, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->beta, npair, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->rcutq, npair, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->Asigma, npair, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->rq, npair, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->rcutsigma, npair, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->Ac, npair, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->zeta, npair, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->rs, npair, MPI_DOUBLE, 0, world);
  MPI_Bcast(setfl->tp, npair, MPI_INT, 0, world);

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
  memory->destroy(cutforcesq);
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
                          double **spline, double /*origin*/)
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
                               int /*pbc_flag*/, int * /*pbc*/)
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
  double bytes = (double)maxeatom * sizeof(double);
  bytes += (double)maxvatom*6 * sizeof(double);
  bytes += (double)2 * nmax * sizeof(double);
  return bytes;
}

EIMPotentialFileReader::EIMPotentialFileReader(LAMMPS *lmp,
                                               const std::string &filename,
                                               const int auto_convert) :
  Pointers(lmp), filename(filename)
{
  if (comm->me != 0) {
    error->one(FLERR, "EIMPotentialFileReader should only be called by proc 0!");
  }

  int unit_convert = auto_convert;
  FILE *fp = utils::open_potential(filename, lmp, &unit_convert);
  conversion_factor = utils::get_conversion_factor(utils::ENERGY,unit_convert);

  if (fp == nullptr) {
    error->one(FLERR,"cannot open eim potential file {}", filename);
  }

  parse(fp);

  fclose(fp);
}

std::pair<std::string, std::string> EIMPotentialFileReader::get_pair(const std::string &a, const std::string &b) {
  if (a < b) {
    return std::make_pair(a, b);
  }
  return std::make_pair(b, a);
}

char *EIMPotentialFileReader::next_line(FILE * fp) {
  // concatenate lines if they end with '&'
  // strip comments after '#'
  int n = 0;
  int nwords = 0;
  bool concat = false;

  char *ptr = fgets(line, MAXLINE, fp);

  if (ptr == nullptr) {
    // EOF
    return nullptr;
  }

  // strip comment
  if ((ptr = strchr(line, '#'))) *ptr = '\0';

  // strip ampersand
  if ((ptr = strrchr(line, '&'))) {
    concat = true;
    *ptr = '\0';
  }

  nwords = utils::count_words(line);

  if (nwords > 0) {
    n = strlen(line);
  }

  while (n == 0 || concat) {
    ptr = fgets(&line[n], MAXLINE - n, fp);

    if (ptr == nullptr) {
      // EOF
      return line;
    }

    // strip comment
    if ((ptr = strchr(line, '#'))) *ptr = '\0';

    // strip ampersand
    if ((ptr = strrchr(line, '&'))) {
      concat = true;
      *ptr = '\0';
    } else {
      concat = false;
    }

    nwords += utils::count_words(&line[n]);

    // skip line if blank
    if (nwords > 0) {
      n = strlen(line);
    }
  }

  return line;
}

void EIMPotentialFileReader::parse(FILE * fp)
{
  char *line = nullptr;
  bool found_global = false;

  while ((line = next_line(fp))) {
    ValueTokenizer values(line);
    std::string type = values.next_string();

    if (type == "global:") {
      if (values.count() != 4) {
        error->one(FLERR, "Invalid global line in EIM potential file");
      }

      division = values.next_double();
      rbig     = values.next_double();
      rsmall   = values.next_double();

      found_global = true;
    } else if (type == "element:") {
      if (values.count() != 9) {
        error->one(FLERR, "Invalid element line in EIM potential file");
      }

      std::string name = values.next_string();

      ElementData data;
      data.ielement   = values.next_int();
      data.mass       = values.next_double();
      data.negativity = values.next_double();
      data.ra         = values.next_double();
      data.ri         = values.next_double();
      data.Ec         = values.next_double();
      data.q0         = values.next_double();

      if (elements.find(name) == elements.end()) {
        elements[name] = data;
      } else {
        error->one(FLERR, "Duplicate pair line in EIM potential file");
      }

    } else if (type == "pair:") {
      if (values.count() != 17) {
        error->one(FLERR, "Invalid element line in EIM potential file");
      }

      std::string elementA = values.next_string();
      std::string elementB = values.next_string();

      PairData data;
      data.rcutphiA  = values.next_double();
      data.rcutphiR  = values.next_double();
      data.Eb        = values.next_double() * conversion_factor;
      data.r0        = values.next_double();
      data.alpha     = values.next_double();
      data.beta      = values.next_double();
      data.rcutq     = values.next_double();
      data.Asigma    = values.next_double();
      data.rq        = values.next_double();
      data.rcutsigma = values.next_double();
      data.Ac        = values.next_double() * conversion_factor;
      data.zeta      = values.next_double();
      data.rs        = values.next_double();

      // should be next_int, but since existing potential files have 1.0000e+00 format
      // we're doing this instead to keep compatibility
      data.tp       = (int)values.next_double();

      auto p = get_pair(elementA, elementB);

      if (pairs.find(p) == pairs.end()) {
        pairs[p] = data;
      } else {
        error->one(FLERR, "Duplicate pair line in EIM potential file");
      }
    }
  }

  if (!found_global) {
    error->one(FLERR, "Missing global line in EIM potential file");
  }
}

void EIMPotentialFileReader::get_global(PairEIM::Setfl *setfl) {
  setfl->division  = division;
  setfl->rbig      = rbig;
  setfl->rsmall    = rsmall;
}

void EIMPotentialFileReader::get_element(PairEIM::Setfl *setfl, int i,
                                         const std::string &name) {
  if (elements.find(name) == elements.end())
    error->one(FLERR,"Element " + name + " not defined in EIM potential file");

  ElementData &data = elements[name];
  setfl->ielement[i] = data.ielement;
  setfl->mass[i] = data.mass;
  setfl->negativity[i] = data.negativity;
  setfl->ra[i] = data.ra;
  setfl->ri[i] = data.ri;
  setfl->Ec[i] = data.Ec;
  setfl->q0[i] = data.q0;
}

void EIMPotentialFileReader::get_pair(PairEIM::Setfl *setfl, int ij,
                                      const std::string &elemA,
                                      const std::string &elemB) {
  auto p = get_pair(elemA, elemB);

  if (pairs.find(p) == pairs.end())
    error->one(FLERR,"Element pair (" + elemA + ", " + elemB
               + ") is not defined in EIM potential file");

  PairData &data = pairs[p];
  setfl->rcutphiA[ij] = data.rcutphiA;
  setfl->rcutphiR[ij] = data.rcutphiR;
  setfl->Eb[ij] = data.Eb;
  setfl->r0[ij] = data.r0;
  setfl->alpha[ij] = data.alpha;
  setfl->beta[ij] = data.beta;
  setfl->rcutq[ij] = data.rcutq;
  setfl->Asigma[ij] = data.Asigma;
  setfl->rq[ij] = data.rq;
  setfl->rcutsigma[ij] = data.rcutsigma;
  setfl->Ac[ij] = data.Ac;
  setfl->zeta[ij] = data.zeta;
  setfl->rs[ij] = data.rs;
  setfl->tp[ij] = data.tp;
}
