// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
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
   Contributing authors:
   Tanmoy Sanyal, M.Scott Shell, UC Santa Barbara
   David Rosenberger, TU Darmstadt
------------------------------------------------------------------------- */

#include "pair_local_density.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cstring>

using namespace LAMMPS_NS;

#define MAXLINE 1024

static const char cite_pair_local_density[] =
  "pair_style  local/density  command:\n\n"
  "@Article{Sanyal16,\n"
  " author =  {T.Sanyal and M.Scott Shell},\n"
  " title =   {Coarse-grained models using local-density potentials optimized with the relative entropy: Application to implicit solvation},\n"
  " journal = {J.~Chem.~Phys.},\n"
  " year =    2016,\n"
  " DOI = doi.org/10.1063/1.4958629"
  "}\n\n"
  "@Article{Sanyal18,\n"
  " author =  {T.Sanyal and M.Scott Shell},\n"
  " title =   {Transferable coarse-grained models of liquid-liquid equilibrium using local density potentials optimized with the relative entropy},\n"
  " journal = {J.~Phys.~Chem. B},\n"
  " year =    2018,\n"
  " DOI = doi.org/10.1021/acs.jpcb.7b12446"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairLocalDensity::PairLocalDensity(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  one_coeff = 1;
  single_enable = 1;

  // stuff read from tabulated file
  nLD = 0;
  nrho = 0;
  rho_min = nullptr;
  rho_max = nullptr;
  a = nullptr;
  b = nullptr;
  c0 = nullptr;
  c2 = nullptr;
  c4 = nullptr;
  c6 = nullptr;
  uppercut = nullptr;
  lowercut = nullptr;
  uppercutsq = nullptr;
  lowercutsq = nullptr;
  frho = nullptr;
  rho = nullptr;

  // splined arrays
  frho_spline = nullptr;

  // per-atom arrays
  nmax = 0;
  fp = nullptr;
  localrho = nullptr;

  // set comm size needed by this pair
  comm_forward = 1;
  comm_reverse = 1;

  // cite publication
  if (lmp->citeme) lmp->citeme->add(cite_pair_local_density);
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairLocalDensity::~PairLocalDensity()
{

  memory->destroy(localrho);
  memory->destroy(fp);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }

  memory->destroy(frho_spline);

  memory->destroy(rho_min);
  memory->destroy(rho_max);
  memory->destroy(delta_rho);
  memory->destroy(c0);
  memory->destroy(c2);
  memory->destroy(c4);
  memory->destroy(c6);
  memory->destroy(uppercut);
  memory->destroy(lowercut);
  memory->destroy(uppercutsq);
  memory->destroy(lowercutsq);
  memory->destroy(frho);
  memory->destroy(rho);

  memory->destroy(a);
  memory->destroy(b);
}

/* ---------------------------------------------------------------------- */

void PairLocalDensity::compute(int eflag, int vflag)
{

  int i,j,ii,jj,m,k,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double rsqinv, phi, uLD, dphi, evdwl,fpair;
  double p, *coeff;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);
  phi = uLD = evdwl = fpair = rsqinv = 0.0;

  /* localrho = LD at each atom
     fp = derivative of embedding energy at each atom for each LD potential
     uLD = embedding energy of each atom due to each LD potential*/

  // grow LD and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(localrho);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(localrho, nLD, nmax, "pairLD:localrho");
    memory->create(fp, nLD, nmax, "pairLD:fp");
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

  // zero out LD and fp

  if (newton_pair) {
    m = nlocal + atom->nghost;
    for (k = 0; k < nLD; k++) {
        for (i = 0; i < m; i++) {
            localrho[k][i] = 0.0;
            fp[k][i] = 0.0;
        }
    }
  }
  else {
    for (k = 0; k < nLD; k++) {
        for (i = 0; i < nlocal; i++) {
            localrho[k][i] = 0.0;
            fp[k][i] = 0.0;
        }
    }
   }

  // loop over neighs of central atoms and types of LDs

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

      // calculate distance-squared between i,j atom-types

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      // calculating LDs based on central and neigh filters

      for (k = 0; k < nLD; k++) {
        if (rsq < lowercutsq[k]) {
             phi = 1.0;
        }
          else if (rsq > uppercutsq[k]) {
             phi = 0.0;
        }
          else {
             phi = c0[k] + rsq * (c2[k] + rsq * (c4[k] + c6[k]*rsq));
        }
        localrho[k][i] += (phi * b[k][jtype]);

        /*checking for both i,j is necessary
        since a half neighbor list is processed.*/

        if (newton_pair || j<nlocal) {
            localrho[k][j] += (phi * b[k][itype]);
        }
      }
    }
  }

  // communicate and sum LDs over all procs
  if (newton_pair) comm->reverse_comm_pair(this);

  //

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    uLD = 0.0;

    for (k = 0; k < nLD; k++) {

        /*skip over this loop if the LD potential
          is not intendend for central atomtype <itype>*/
        if (!(a[k][itype])) continue;

        // linear extrapolation at rho_min and rho_max

        if (localrho[k][i] <= rho_min[k]) {
            coeff = frho_spline[k][0];
            fp[k][i] = coeff[2];
            uLD +=  a[k][itype] * ( coeff[6] + fp[k][i]*(localrho[k][i] - rho_min[k]) );
        }
        else if (localrho[k][i] >= rho_max[k]) {
            coeff = frho_spline[k][nrho-2];
            fp[k][i] = coeff[0] + coeff[1]  + coeff[2];
            uLD +=  a[k][itype] * ( (coeff[3] + coeff[4] + coeff[5] + coeff[6]) + fp[k][i]*(localrho[k][i] - rho_max[k]) );
        }
        else {
            p = (localrho[k][i] - rho_min[k]) / delta_rho[k];
            m = static_cast<int> (p);
            m = MAX(0, MIN(m, nrho-2));
            p -= m;
            p = MIN(p, 1.0);
            coeff = frho_spline[k][m];
            fp[k][i] = (coeff[0]*p + coeff[1])*p + coeff[2];
            uLD +=  a[k][itype] * (((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6]);
        }
    }

    if (eflag) {
        if (eflag_global) eng_vdwl += uLD;
        if (eflag_atom) eatom[i] += uLD;
    }
 }

  // communicate LD and fp to all procs

  comm->forward_comm_pair(this);

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

      // calculate square of distance between i,j atoms

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      // calculate force between two atoms
      fpair = 0.0;
      if (rsq < cutforcesq) {   // global cutoff check
        rsqinv = 1.0/rsq;
        for (k = 0; k < nLD; k++) {
            if (rsq >= lowercutsq[k] && rsq < uppercutsq[k]) {
               dphi = rsq * (2.0*c2[k] + rsq * (4.0*c4[k] + 6.0*c6[k]*rsq));
               fpair += -(a[k][itype]*b[k][jtype]*fp[k][i] + a[k][jtype]*b[k][itype]*fp[k][j]) * dphi;
            }
        }
        fpair *= rsqinv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
        }

      /*eng_vdwl has already been completely built,
        so no need to add anything here*/

      if (eflag) evdwl = 0.0;

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

void PairLocalDensity::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLocalDensity::settings(int narg, char ** /* arg */)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for all type pairs
   read tabulated LD input file
------------------------------------------------------------------------- */

void PairLocalDensity::coeff(int narg, char **arg)
{
  int i, j;
  if (!allocated) allocate();

  if (narg != 3) error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // parse LD file

  parse_file(arg[2]);


 // clear setflag since coeff() called once with I,J = * *

  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
      setflag[i][j] = 0;

  // set setflag for all i,j type pairs

  int count = 0;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
        setflag[i][j] = 1;
        count++;
      }
    }
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLocalDensity::init_style()
{
  // spline rho and frho arrays
  // request half neighbor list

  array2spline();

  // half neighbor request
  neighbor->request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLocalDensity::init_one(int /* i */, int /* j */)
{
  // single global cutoff = max of all uppercuts read in from LD file

  cutmax = 0.0;
  for (int k = 0; k < nLD; k++)
    cutmax = MAX(cutmax,uppercut[k]);

  cutforcesq = cutmax*cutmax;

  return cutmax;
}


/*--------------------------------------------------------------------------
  pair_write functionality for this pair style that gives just a snap-shot
  of the LD potential without doing an actual MD run
 ---------------------------------------------------------------------------*/

double PairLocalDensity::single(int /* i */, int /* j */, int itype, int jtype,
                                double rsq, double /* factor_coul */,
                                double /* factor_lj */, double &fforce)
{
    int m, k, index;
    double rsqinv, p, uLD;
    double *coeff, **LD;
    double dFdrho, phi, dphi;

    uLD = dFdrho = dphi = 0.0;

    memory->create(LD, nLD, 3, "pairLD:LD");
    for (k = 0; k < nLD; k++) {
        LD[k][1] = 0.0; // itype:- 1
        LD[k][2] = 0.0; // jtype:- 2
        }

    rsqinv = 1.0/rsq;
    for (k = 0; k < nLD; k++) {
        if (rsq < lowercutsq[k]) {
             phi = 1.0;
        }
        else if (rsq > uppercutsq[k]) {
             phi = 0.0;
        }
        else {
             phi = c0[k] + rsq * (c2[k] + rsq * (c4[k] + c6[k]*rsq));
        }
        LD[k][1] += (phi * b[k][jtype]);
        LD[k][2] += (phi * b[k][itype]);
    }

    for (k = 0; k < nLD; k++) {
        if (a[k][itype]) index = 1;
        if (a[k][jtype]) index = 2;

        if (LD[k][index] <= rho_min[k]) {
            coeff = frho_spline[k][0];
            dFdrho = coeff[2];
            uLD +=  a[k][itype] * ( coeff[6] + dFdrho*(LD[k][index] - rho_min[k]) );
        }
        else if (LD[k][index] >= rho_max[k]) {
            coeff = frho_spline[k][nrho-1];
            dFdrho = coeff[0] + coeff[1]  + coeff[2];
            uLD +=  a[k][itype] * ( (coeff[3] + coeff[4] + coeff[5] + coeff[6]) + dFdrho*(LD[k][index] - rho_max[k]) );
        }
        else {
            p = (LD[k][index] - rho_min[k]) / delta_rho[k];
            m = static_cast<int> (p);
            m = MAX(0, MIN(m, nrho-2));
            p -= m;
            p = MIN(p, 1.0);
            coeff = frho_spline[k][m];
            dFdrho = (coeff[0]*p + coeff[1])*p + coeff[2];
            uLD +=  a[k][itype] * (((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6]);
        }

        if (rsq < lowercutsq[k]) {
           dphi = 0.0;
        }
        else if (rsq > uppercutsq[k]) {
           dphi = 0.0;
        }
        else {
           dphi = rsq * (2.0*c2[k] + rsq * (4.0*c4[k] + 6.0*c6[k]*rsq));
        }
        fforce +=  -(a[k][itype]*b[k][jtype]*dFdrho + a[k][jtype]*b[k][itype]*dFdrho) * dphi *rsqinv;
    }
    memory->destroy(LD);

    return uLD;
}

/*--------------------------------------------------------------------
   spline the array frho read in from the file to create
   frho_spline
---------------------------------------------------------------------- */

void PairLocalDensity::array2spline() {
  memory->destroy(frho_spline);
  memory->create(frho_spline, nLD, nrho, 7, "pairLD:frho_spline");

  for (int k = 0; k < nLD; k++)
    interpolate_cbspl(nrho, delta_rho[k], frho[k], frho_spline[k]);

}

/* ----------------------------------------------------------------------
  (one-dimensional) cubic spline interpolation sub-routine,
  which determines the coeffs for a clamped cubic spline
  given tabulated data
 ------------------------------------------------------------------------*/

void PairLocalDensity::interpolate_cbspl(int n, double delta,
                                         double *f, double **spline)
{
/*   inputs:
          n         number of interpolating points

          f                 array containing function values to
                                be interpolated;  f[i] is the function
                                value corresponding to x[i]
                                    ('x' refers to the independent var)

              delta     difference in tabulated values of x

         outputs: (packaged as columns of the coeff matrix)
          coeff_b       coeffs of linear terms
              coeff_c   coeffs of quadratic terms
              coeff_d   coeffs of cubic terms
          spline    matrix that collects b,c,d


         other parameters:
          fpa           derivative of function at x=a
          fpb           derivative of function at x=b
*/

     double *dl, *dd, *du;
     double *coeff_b, *coeff_c, *coeff_d;
     double fpa, fpb;

     int i;

     coeff_b = new double [n];
     coeff_c = new double [n];
     coeff_d = new double [n];
     dl = new double [n];
     dd = new double [n];
     du = new double [n];

     // initialize values
     for ( i = 0; i<n; i++) {
         coeff_b[i] = coeff_c[i] = coeff_d[i] = 0.;
         dl[i] = dd[i] = du[i] = 0.;
     }

     // set slopes at beginning and end
     fpa = 0.;
     fpb = 0.;

     for (i = 0; i < n-1; i++) {
         dl[i] = du[i] = delta;
     }

     dd[0] = 2.0 * delta;
     dd[n-1] = 2.0 * delta;
     coeff_c[0] = ( 3.0 / delta ) * ( f[1] - f[0] ) - 3.0 * fpa;
     coeff_c[n-1] = 3.0 * fpb - ( 3.0 / delta ) * ( f[n-1] - f[n-2] );
     for (i = 0; i < n-2; i++) {
         dd[i+1] = 4.0 * delta;
         coeff_c[i+1] = ( 3.0 / delta ) * ( f[i+2] - f[i+1] ) -
                        ( 3.0 / delta ) * ( f[i+1] - f[i] );
     }

     // tridiagonal solver
     for (i = 0; i < n-1; i++) {
         du[i] /= dd[i];
         dd[i+1] -= dl[i]*du[i];
     }

     coeff_c[0] /= dd[0];
     for (i = 1; i < n; i++)
         coeff_c[i] = ( coeff_c[i] - dl[i-1] * coeff_c[i-1] ) / dd[i];

     for (i = n-2; i >= 0; i--)
         coeff_c[i] -= coeff_c[i+1] * du[i];

     for (i = 0; i < n-1; i++) {
         coeff_d[i] = ( coeff_c[i+1] - coeff_c[i] ) / ( 3.0 * delta );
         coeff_b[i] = ( f[i+1] - f[i] ) / delta - delta * ( coeff_c[i+1] + 2.0*coeff_c[i] ) / 3.0;
     }

     // normalize
     for (i = 0; i < n-1; i++) {
         coeff_b[i] = coeff_b[i] * delta ;
         coeff_c[i] = coeff_c[i] * delta*delta ;
         coeff_d[i] = coeff_d[i] * delta*delta*delta;
     }

     //copy to coefficient matrix
     for ( i = 0; i < n; i++) {
         spline[i][3] = coeff_d[i];
         spline[i][4] = coeff_c[i];
         spline[i][5] = coeff_b[i];
         spline[i][6] = f[i];
         spline[i][2] = spline[i][5]/delta;
         spline[i][1] = 2.0*spline[i][4]/delta;
         spline[i][0] = 3.0*spline[i][3]/delta;
     }

     delete [] coeff_b;
     delete [] coeff_c;
     delete [] coeff_d;
     delete [] du;
     delete [] dd;
     delete [] dl;
}

/* ----------------------------------------------------------------------
   read potential values from tabulated LD input file
------------------------------------------------------------------------- */

void PairLocalDensity::parse_file(char *filename) {

  int k, n;
  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];
  double ratio, lc2, uc2, denom;

  if (me == 0) {
    fptr = fopen(filename, "r");
    if (fptr == nullptr)
      error->one(FLERR,"Cannot open Local Density potential file {}: {}",filename,utils::getsyserror());
  }

  double *ftmp; // tmp var to extract the complete 2D frho array from file

  // broadcast number of LD potentials and number of (rho,frho) pairs
  if (me == 0) {

    // first 2 comment lines ignored
    utils::sfgets(FLERR,line,MAXLINE,fptr,filename,error);
    utils::sfgets(FLERR,line,MAXLINE,fptr,filename,error);

    // extract number of potentials and number of (frho, rho) points
    utils::sfgets(FLERR,line,MAXLINE,fptr,filename,error);
    sscanf(line, "%d %d", &nLD, &nrho);
    utils::sfgets(FLERR,line,MAXLINE,fptr,filename,error);
  }

  MPI_Bcast(&nLD,1,MPI_INT,0,world);
  MPI_Bcast(&nrho,1,MPI_INT,0,world);

  // setting up all arrays to be read from files and broadcasted
  memory->create(uppercut, nLD, "pairLD:uppercut");
  memory->create(lowercut, nLD, "pairLD:lowercut");
  memory->create(uppercutsq, nLD, "pairLD:uppercutsq");
  memory->create(lowercutsq, nLD, "pairLD:lowercutsq");
  memory->create(c0, nLD, "pairLD:c0");
  memory->create(c2, nLD, "pairLD:c2");
  memory->create(c4, nLD, "pairLD:c4");
  memory->create(c6, nLD, "pairLD:c6");
  memory->create(rho_min,  nLD, "pairLD:rho_min");
  memory->create(rho_max,  nLD, "pairLD:rho_max");
  memory->create(delta_rho, nLD,"pairLD:delta_rho");
  memory->create(ftmp, nrho*nLD, "pairLD:ftmp");

  // setting up central and neighbor atom filters
  memory->create(a, nLD, atom->ntypes+1 , "pairLD:a");
  memory->create(b, nLD, atom->ntypes+1, "pairLD:b");
  if (me == 0) {
    for (n = 1; n <= atom->ntypes; n++) {
        for (k = 0; k < nLD; k++) {
            a[k][n] = 0;
            b[k][n] = 0;
        }
    }
  }

 // read file block by block

  if (me == 0) {
    for (k = 0; k < nLD; k++) {

        // parse upper and lower cut values
        if (fgets(line,MAXLINE,fptr)==nullptr) break;
        sscanf(line, "%lf %lf", &lowercut[k], &uppercut[k]);

        // parse and broadcast central atom filter
        utils::sfgets(FLERR,line, MAXLINE, fptr,filename,error);
        char *tmp = strtok(line, " /t/n/r/f");
        while (tmp != nullptr) {
            a[k][atoi(tmp)] = 1;
            tmp = strtok(nullptr, " /t/n/r/f");
        }

        // parse neighbor atom filter
        utils::sfgets(FLERR,line, MAXLINE, fptr,filename,error);
        tmp = strtok(line, " /t/n/r/f");
        while (tmp != nullptr) {
            b[k][atoi(tmp)] = 1;
            tmp = strtok(nullptr, " /t/n/r/f");
        }

        // parse min, max and delta rho values
        utils::sfgets(FLERR,line, MAXLINE, fptr,filename,error);
        sscanf(line, "%lf %lf %lf", &rho_min[k], &rho_max[k], &delta_rho[k]);
        // recompute delta_rho from scratch for precision
        delta_rho[k] = (rho_max[k] - rho_min[k]) / (nrho - 1);

        // parse tabulated frho values from each line into temporary array
        for (n = 0; n < nrho; n++) {
          utils::sfgets(FLERR,line,MAXLINE,fptr,filename,error);
            sscanf(line, "%lf", &ftmp[k*nrho + n]);
        }

        // ignore blank line at the end of every block
        utils::sfgets(FLERR,line,MAXLINE,fptr,filename,error);

        // set coefficients for local density indicator function
        uc2 = uppercut[k] * uppercut[k];
        uppercutsq[k] = uc2;
        lc2 = lowercut[k] * lowercut[k];
        lowercutsq[k] = lc2;
        ratio = lc2/uc2;
        denom = 1.0 - ratio;
        denom = denom*denom*denom;
        c0[k] = (1 - 3.0 * ratio) / denom;
        c2[k] = (6.0 * ratio) / (uc2 * denom);
        c4[k] = -(3.0 + 3.0*ratio) / (uc2*uc2 * denom);
        c6[k] = 2.0 / (uc2*uc2*uc2 * denom);
      }
  }

  // Broadcast all parsed arrays
  MPI_Bcast(&lowercut[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&uppercut[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&lowercutsq[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&uppercutsq[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&c0[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&c2[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&c4[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&c6[0], nLD, MPI_DOUBLE, 0, world);
  for (k = 0; k < nLD; k++) {
      MPI_Bcast(&a[k][1], atom->ntypes, MPI_INT, 0, world);
      MPI_Bcast(&b[k][1], atom->ntypes, MPI_INT, 0, world);
  }
  MPI_Bcast(&rho_min[0],  nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&rho_max[0],  nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&delta_rho[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&ftmp[0], nLD*nrho, MPI_DOUBLE, 0, world);

  if (me == 0) fclose(fptr);

  // set up rho and frho arrays
  memory->create(rho, nLD, nrho, "pairLD:rho");
  memory->create(frho, nLD, nrho, "pairLD:frho");

  for (k = 0; k < nLD; k++) {
    for (n = 0; n < nrho; n++) {
        rho[k][n] = rho_min[k] + n*delta_rho[k];
        frho[k][n] = ftmp[k*nrho + n];
    }
 }

  // delete temporary array
  memory->destroy(ftmp);
}

/* ----------------------------------------------------------------------
   communication routines
------------------------------------------------------------------------- */

int PairLocalDensity::pack_comm(int n, int *list, double *buf,
                                int /* pbc_flag */, int * /* pbc */) {
  int i,j,k;
  int m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 0; k < nLD; k++) {
      buf[m++] = fp[k][j];
    }
  }

  return nLD;
}

/* ---------------------------------------------------------------------- */

void PairLocalDensity::unpack_comm(int n, int first, double *buf) {

  int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (k = 0; k < nLD; k++) {
      fp[k][i] = buf[m++];
    }
 }
}

/* ---------------------------------------------------------------------- */

int PairLocalDensity::pack_reverse_comm(int n, int first, double *buf) {

  int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (k = 0; k < nLD; k++) {
      buf[m++] = localrho[k][i];
    }
  }
  return nLD;
}

/* ---------------------------------------------------------------------- */

void PairLocalDensity::unpack_reverse_comm(int n, int *list, double *buf) {

  int i,j,k;
  int m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 0; k < nLD; k++) {
      localrho[k][j] += buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairLocalDensity::memory_usage()
{
  double bytes = (double)maxeatom * sizeof(double);
  bytes += (double)maxvatom*6 * sizeof(double);
  bytes += (double)2 * (nmax*nLD) * sizeof(double);
  return bytes;
}

