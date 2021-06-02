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

#include "pair_sph_lj.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHLJ::PairSPHLJ(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
}

/* ---------------------------------------------------------------------- */

PairSPHLJ::~PairSPHLJ() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(viscosity);
  }
}

/* ---------------------------------------------------------------------- */

void PairSPHLJ::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc, h, ih, ihsq, ihcub;
  double rsq, wfd, delVdotDelR, mu, deltaE, ci, cj, lrc;

  ev_init(eflag, vflag);

  double **v = atom->vest;
  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *desph = atom->desph;
  double *esph = atom->esph;
  double *cv = atom->cv;
  double *drho = atom->drho;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = mass[itype];

    // compute pressure of particle i with LJ EOS
    LJEOS2(rho[i], esph[i], cv[i], &fi, &ci);
    fi /= (rho[i] * rho[i]);
    //printf("fi = %f\n", fi);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jmass = mass[jtype];

      if (rsq < cutsq[itype][jtype]) {
        h = cut[itype][jtype];
        ih = 1.0 / h;
        ihsq = ih * ih;
        ihcub = ihsq * ih;

        wfd = h - sqrt(rsq);
        if (domain->dimension == 3) {
          // Lucy Kernel, 3d
          // Note that wfd, the derivative of the weight function with respect to r,
          // is lacking a factor of r.
          // The missing factor of r is recovered by
          // (1) using delV . delX instead of delV . (delX/r) and
          // (2) using f[i][0] += delx * fpair instead of f[i][0] += (delx/r) * fpair
          wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
        } else {
          // Lucy Kernel, 2d
          wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
        }

        // function call to LJ EOS
        LJEOS2(rho[j], esph[j], cv[j], &fj, &cj);
        fj /= (rho[j] * rho[j]);

        // apply long-range correction to model a LJ fluid with cutoff
        // this implies that the modelled LJ fluid has cutoff == SPH cutoff
        lrc = - 11.1701 * (ihcub * ihcub * ihcub - 1.5 * ihcub);
        fi += lrc;
        fj += lrc;

        // dot product of velocity delta and distance vector
        delVdotDelR = delx * (vxtmp - v[j][0]) + dely * (vytmp - v[j][1])
            + delz * (vztmp - v[j][2]);

        // artificial viscosity (Monaghan 1992)
        if (delVdotDelR < 0.) {
          mu = h * delVdotDelR / (rsq + 0.01 * h * h);
          fvisc = -viscosity[itype][jtype] * (ci + cj) * mu / (rho[i] + rho[j]);
        } else {
          fvisc = 0.;
        }

        // total pair force & thermal energy increment
        fpair = -imass * jmass * (fi + fj + fvisc) * wfd;
        deltaE = -0.5 * fpair * delVdotDelR;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        // and change in density
        drho[i] += jmass * delVdotDelR * wfd;

        // change in thermal energy
        desph[i] += deltaE;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
          desph[j] += deltaE;
          drho[j] += imass * delVdotDelR * wfd;
        }

        if (evflag)
          ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHLJ::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(viscosity, n + 1, n + 1, "pair:viscosity");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHLJ::settings(int narg, char **/*arg*/) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of arguments for pair_style sph/lj");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHLJ::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR,
        "Incorrect args for pair_style sph/lj coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR,arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR,arg[1], 1, atom->ntypes, jlo, jhi, error);

  double viscosity_one = utils::numeric(FLERR,arg[2],false,lmp);
  double cut_one = utils::numeric(FLERR,arg[3],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      viscosity[i][j] = viscosity_one;
      printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHLJ::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/lj coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  viscosity[j][i] = viscosity[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHLJ::single(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/,
    double /*rsq*/, double /*factor_coul*/, double /*factor_lj*/, double &fforce) {
  fforce = 0.0;

  return 0.0;
}


/*double PairSPHLJ::LJEOS2(double rho, double e, double cv) {


  double T = e / cv;
  if (T < 1.e-2) T = 1.e-2;
  //printf("%f %f\n", T, rho);
  double iT = 0.1e1 / T;
  //double itpow1_4 = exp(0.25 * log(iT)); //pow(iT, 0.1e1 / 0.4e1);
  double itpow1_4 = pow(iT, 0.1e1 / 0.4e1);
  double x = rho * itpow1_4;
  double xsq = x * x;
  double xpow3 = xsq * x;
  double xpow4 = xsq * xsq;
  double xpow9 = xpow3 * xpow3 * xpow3;


  return (0.1e1 + rho * (0.3629e1 + 0.7264e1 * x + 0.104925e2 * xsq + 0.11460e2
      * xpow3 + 0.21760e1 * xpow9 - itpow1_4 * itpow1_4 * (0.5369e1 + 0.13160e2
      * x + 0.18525e2 * xsq - 0.17076e2 * xpow3 + 0.9320e1 * xpow4) + iT
      * (-0.3492e1 + 0.18698e2 * x - 0.35505e2 * xsq + 0.31816e2 * xpow3
          - 0.11195e2 * xpow4)) * itpow1_4) * rho * T;
}*/


/* --------------------------------------------------------------------------------------------- */
/* Lennard-Jones EOS,
   Francis H. Ree
   "Analytic representation of thermodynamic data for the Lennard‐Jones fluid",
   Journal of Chemical Physics 73 pp. 5401-5403 (1980)
*/

void PairSPHLJ::LJEOS2(double rho, double e, double cv, double *p, double *c) {
  double T = e/cv;
  double beta = 1.0 / T;
  double beta_sqrt = sqrt(beta);
  double x = rho * sqrt(beta_sqrt);

  double xsq = x * x;
  double xpow3 = xsq * x;
  double xpow4 = xsq * xsq;

  /* differential of Helmholtz free energy w.r.t. x */
  double diff_A_NkT = 3.629 + 7.264*x - beta*(3.492 - 18.698*x + 35.505*xsq - 31.816*xpow3 + 11.195*xpow4)
                    - beta_sqrt*(5.369 + 13.16*x + 18.525*xsq - 17.076*xpow3 + 9.32*xpow4)
                    + 10.4925*xsq + 11.46*xpow3 + 2.176*xpow4*xpow4*x;

 /* differential of Helmholtz free energy w.r.t. x^2 */
  double d2A_dx2 = 7.264 + 20.985*x \
                 + beta*(18.698 - 71.01*x + 95.448*xsq - 44.78*xpow3)\
                 - beta_sqrt*(13.16 + 37.05*x - 51.228*xsq + 37.28*xpow3)\
                 + 34.38*xsq + 19.584*xpow4*xpow4;

  // p = rho k T * (1 + rho * d(A/(NkT))/drho)
  // dx/drho = rho/x
  *p = rho * T * (1.0 + diff_A_NkT * x); // pressure
  double csq = T * (1.0 + 2.0 * diff_A_NkT * x + d2A_dx2 * x * x); // soundspeed squared
  if (csq > 0.0) {
    *c = sqrt(csq); // soundspeed
  } else {
    *c = 0.0;
  }
}

/* ------------------------------------------------------------------------------ */

/* Jirí Kolafa, Ivo Nezbeda
 * "The Lennard-Jones fluid: an accurate analytic and theoretically-based equation of state",
 *  Fluid Phase Equilibria 100 pp. 1-34 (1994) */
/*double PairSPHLJ::LJEOS2(double rho, double e, double cv) {
 double T = e / cv;

 double sT = sqrt(T);
 double isT = 1.0 / sT;
 double dC = -0.063920968 * log(T) + 0.011117524 / T - 0.076383859 / sT
 + 1.080142248 + 0.000693129 * sT;
 double eta = 3.141592654 / 6. * rho * (dC * dC * dC);
 double zHS = (1 + eta * (1 + eta * (1 - eta / 1.5 * (1 + eta))))
 / ((1. - eta) * (1. - eta) * (1. - eta));
 double BC = (((((-0.58544978 * isT + 0.43102052) * isT + .87361369) * isT
 - 4.13749995) * isT + 2.90616279) * isT - 7.02181962) / T + 0.02459877;
 double gammaBH = 1.92907278;

 double sum = ((2.01546797 * 2 + rho * ((-28.17881636) * 3 + rho
 * (28.28313847 * 4 + rho * (-10.42402873) * 5))) + (-19.58371655 * 2
 + rho * (+75.62340289 * 3 + rho * ((-120.70586598) * 4 + rho
 * (+93.92740328 * 5 + rho * (-27.37737354) * 6)))) / sqrt(T)
 + ((29.34470520 * 2 + rho * ((-112.35356937) * 3 + rho * (+170.64908980
 * 4 + rho * ((-123.06669187) * 5 + rho * 34.42288969 * 6))))
 + ((-13.37031968) * 2 + rho * (65.38059570 * 3 + rho
 * ((-115.09233113) * 4 + rho * (88.91973082 * 5 + rho
 * (-25.62099890) * 6)))) / T) / T) * rho * rho;
 return ((zHS + BC / exp(gammaBH * rho * rho) * rho * (1 - 2 * gammaBH * rho
 * rho)) * T + sum) * rho;
 }
*/
