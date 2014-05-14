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
   Contributing author: Stan Moore (Sandia)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "pppm_stagger.h"
#include "atom.h"
#include "gridcomm.h"
#include "force.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

#include "math_const.h"
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define OFFSET 16384
#define EPS_HOC 1.0e-7

enum{REVERSE_RHO};
enum{FORWARD_IK,FORWARD_AD,FORWARD_IK_PERATOM,FORWARD_AD_PERATOM};

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ---------------------------------------------------------------------- */

PPPMStagger::PPPMStagger(LAMMPS *lmp, int narg, char **arg) : 
  PPPM(lmp, narg, arg)
{
  if (narg < 1) error->all(FLERR,"Illegal kspace_style pppm/stagger command");
  stagger_flag = 1;
  group_group_enable = 0;

  memory->create(gf_b2,8,7,"pppm_stagger:gf_b2");
  gf_b2[1][0] = 1.0;
  gf_b2[2][0] = 5.0 / 6.0;
  gf_b2[2][1] = 1.0 / 6.0;
  gf_b2[3][0] = 61.0 / 120.0;
  gf_b2[3][1] = 29.0 / 60.0;
  gf_b2[3][2] = 1.0 / 120.0;
  gf_b2[4][0] = 277.0 / 1008.0;
  gf_b2[4][1] = 1037.0 / 1680.0;
  gf_b2[4][2] = 181.0 / 1680.00;
  gf_b2[4][3] = 1.0 / 5040.0;
  gf_b2[5][0] = 50521.0 / 362880.0;
  gf_b2[5][1] = 7367.0 / 12960.0;
  gf_b2[5][2] = 16861.0 / 60480.0;
  gf_b2[5][3] = 1229.0 / 90720.0;
  gf_b2[5][4] = 1.0 / 362880.0;
  gf_b2[6][0] = 540553.0 / 7983360.0;
  gf_b2[6][1] = 17460701.0 / 39916800.0;
  gf_b2[6][2] = 8444893.0 / 19958400.0;
  gf_b2[6][3] = 1409633.0 / 19958400.0;
  gf_b2[6][4] = 44281.0 / 39916800.0;
  gf_b2[6][5] = 1.0 / 39916800.0;
  gf_b2[7][0] = 199360981.0 / 6227020800.0;
  gf_b2[7][1] = 103867703.0 / 345945600.0;
  gf_b2[7][2] = 66714163.0 / 138378240.0;
  gf_b2[7][3] = 54085121.0 / 311351040.0;
  gf_b2[7][4] = 1640063.0 / 138378240.0;
  gf_b2[7][5] = 671.0 / 10483200.0;
  gf_b2[7][6] = 1.0 / 6227020800.0;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

PPPMStagger::~PPPMStagger()
{
  memory->destroy(gf_b2);
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void PPPMStagger::init()
{
  // error check

  if (domain->triclinic)
   error->all(FLERR,"Cannot (yet) use kspace_style pppm/stagger "
              "with triclinic systems");

  PPPM::init();
}

/* ----------------------------------------------------------------------
   compute the PPPM long-range force, energy, virial
------------------------------------------------------------------------- */

void PPPMStagger::compute(int eflag, int vflag)
{
  int i,j;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  if (evflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
    cg_peratom->ghost_notify();
    cg_peratom->setup();
  }

  // convert atoms from box to lamda coords

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // extend size of per-atom arrays if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm:part2grid");
  }

  nstagger = 2;

  stagger = 0.0;
  for (int n=0; n<nstagger; n++) {

    // find grid points for all my particles
    // map my particle charge onto my local 3d density grid

    particle_map();
    make_rho();
    
    // all procs communicate density values from their ghost cells
    //   to fully sum contribution in their 3d bricks
    // remap from 3d decomposition to FFT decomposition

    cg->reverse_comm(this,REVERSE_RHO);
    brick2fft();

    // compute potential gradient on my FFT grid and
    //   portion of e_long on this proc's FFT grid
    // return gradients (electric fields) in 3d brick decomposition
    // also performs per-atom calculations via poisson_peratom()

    poisson();

    // all procs communicate E-field values
    // to fill ghost cells surrounding their 3d bricks
    
    if (differentiation_flag == 1) cg->forward_comm(this,FORWARD_AD);
    else cg->forward_comm(this,FORWARD_IK);

    // extra per-atom energy/virial communication

    if (evflag_atom) {
      if (differentiation_flag == 1 && vflag_atom) 
        cg_peratom->forward_comm(this,FORWARD_AD_PERATOM);
      else if (differentiation_flag == 0)
        cg_peratom->forward_comm(this,FORWARD_IK_PERATOM);
    }

    // calculate the force on my particles

    fieldforce();

    // extra per-atom energy/virial communication

    if (evflag_atom) fieldforce_peratom();

    stagger += 1.0/float(nstagger);
  }

  // sum global energy across procs and add in volume-dependent term
  // reset qsum and qsqsum if atom count has changed

  const double qscale = qqrd2e * scale;

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;

    if (atom->natoms != natoms_original) {
      qsum_qsq(0);
      natoms_original = atom->natoms;
    }

    energy *= 0.5*volume/float(nstagger);
    energy -= g_ewald*qsqsum/MY_PIS +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= qscale;
  }

  // sum global virial across procs

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) 
      virial[i] = 0.5*qscale*volume*virial_all[i]/float(nstagger);
  }

  // per-atom energy/virial
  // energy includes self-energy correction
  // notal accounts for TIP4P tallying eatom/vatom for ghost atoms

  if (evflag_atom) {
    double *q = atom->q;
    int nlocal = atom->nlocal;
    int ntotal = nlocal;
    if (tip4pflag) ntotal += atom->nghost;

    if (eflag_atom) {
      for (i = 0; i < nlocal; i++) {
        eatom[i] *= 0.5;
        eatom[i] -= g_ewald*q[i]*q[i]/MY_PIS + MY_PI2*q[i]*qsum /
          (g_ewald*g_ewald*volume);
        eatom[i] *= qscale;
      }
      for (i = nlocal; i < ntotal; i++) eatom[i] *= 0.5*qscale;
    }

    if (vflag_atom) {
      for (i = 0; i < ntotal; i++)
        for (j = 0; j < 6; j++) vatom[i][j] *= 0.5*qscale;
    }
  }

  // 2d slab correction

  if (slabflag == 1) slabcorr();

  // convert atoms back from lamda to box coords

  if (triclinic) domain->lamda2x(atom->nlocal);
}

/* ----------------------------------------------------------------------
   compute qopt
------------------------------------------------------------------------- */

double PPPMStagger::compute_qopt()
{
  if (differentiation_flag == 1)
    return compute_qopt_ad();

  double qopt = 0.0;
  const double * const prd = domain->prd;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double zprd_slab = zprd*slab_volfactor;
  const double unitkx = (MY_2PI/xprd);
  const double unitky = (MY_2PI/yprd);
  const double unitkz = (MY_2PI/zprd_slab);

  double snx,sny,snz;
  double cnx,cny,cnz;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,sum2,dot1,dot2;
  double numerator,denominator;
  double u1,u2,u3,sqk;

  int k,l,m,nx,ny,nz,kper,lper,mper;

  const int nbx = 2;
  const int nby = 2;
  const int nbz = 2;

  const int twoorder = 2*order;

  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);
    snz = square(sin(0.5*unitkz*mper*zprd_slab/nz_pppm));
    cnz = cos(0.5*unitkz*mper*zprd_slab/nz_pppm);

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);
      sny = square(sin(0.5*unitky*lper*yprd/ny_pppm));
      cny = cos(0.5*unitky*lper*yprd/ny_pppm);

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);
        snx = square(sin(0.5*unitkx*kper*xprd/nx_pppm));
        cnx = cos(0.5*unitkx*kper*xprd/nx_pppm);

        sqk = square(unitkx*kper) + square(unitky*lper) + square(unitkz*mper);

        if (sqk != 0.0) {
          numerator = MY_4PI/sqk;
          denominator = 0.5*(gf_denom(snx,sny,snz) + gf_denom2(cnx,cny,cnz));
          sum1 = 0.0;
          sum2 = 0.0;

          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm*nx);
            sx = exp(-0.25*square(qx/g_ewald));
            argx = 0.5*qx*xprd/nx_pppm;
            wx = powsinxx(argx,twoorder);

            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm*ny);
              sy = exp(-0.25*square(qy/g_ewald));
              argy = 0.5*qy*yprd/ny_pppm;
              wy = powsinxx(argy,twoorder);

              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm*nz);
                sz = exp(-0.25*square(qz/g_ewald));
                argz = 0.5*qz*zprd_slab/nz_pppm;
                wz = powsinxx(argz,twoorder);

                dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
                dot2 = qx*qx + qy*qy + qz*qz;
                u1   = sx*sy*sz;
                u2   = wx*wy*wz;
                u3   = numerator*u1*u2*dot1;
                sum1 += u1*u1*MY_4PI*MY_4PI/dot2;
                sum2 += u3*u3/dot2;
              }
            }
          }
          qopt += sum1 - sum2/denominator;
        }
      }
    }
  }
  double qopt_all;
  MPI_Allreduce(&qopt,&qopt_all,1,MPI_DOUBLE,MPI_SUM,world);
  return qopt_all;
}

/* ----------------------------------------------------------------------
   compute qopt_ad
------------------------------------------------------------------------- */

double PPPMStagger::compute_qopt_ad()
{
  double qopt = 0.0;
  const double * const prd = domain->prd;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double zprd_slab = zprd*slab_volfactor;
  const double unitkx = (MY_2PI/xprd);
  const double unitky = (MY_2PI/yprd);
  const double unitkz = (MY_2PI/zprd_slab);

  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,sum2,sum3,sum4,sum5,sum6,dot2;
  double u1,u2,sqk;

  int k,l,m,nx,ny,nz,kper,lper,mper;

  const int nbx = 2;
  const int nby = 2;
  const int nbz = 2;

  const int twoorder = 2*order;

  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);

        sqk = square(unitkx*kper) + square(unitky*lper) + square(unitkz*mper);

        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          sum4 = 0.0;
          sum5 = 0.0;
          sum6 = 0.0;

          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm*nx);
            sx = exp(-0.25*square(qx/g_ewald));
            argx = 0.5*qx*xprd/nx_pppm;
            wx = powsinxx(argx,twoorder);

            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm*ny);
              sy = exp(-0.25*square(qy/g_ewald));
              argy = 0.5*qy*yprd/ny_pppm;
              wy = powsinxx(argy,twoorder);

              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm*nz);
                sz = exp(-0.25*square(qz/g_ewald));
                argz = 0.5*qz*zprd_slab/nz_pppm;
                wz = powsinxx(argz,twoorder);

                dot2 = qx*qx + qy*qy + qz*qz;
                u1   = sx*sy*sz;
                u2   = wx*wy*wz;
                sum1 += u1*u1/dot2*MY_4PI*MY_4PI;
                sum2 += u1*u1*u2*u2*MY_4PI*MY_4PI;
                sum3 += u2;
                sum4 += dot2*u2;
                sum5 += u2*powint(-1.0,nx+ny+nz);
                sum6 += dot2*u2*powint(-1.0,nx+ny+nz);
              }
            }
          }
          qopt += sum1 - sum2/(0.5*(sum3*sum4 + sum5*sum6));
        }
      }
    }
  }
  double qopt_all;
  MPI_Allreduce(&qopt,&qopt_all,1,MPI_DOUBLE,MPI_SUM,world);
  return qopt_all;
}

/* ----------------------------------------------------------------------
   pre-compute Green's function denominator expansion coeffs, Gamma(2n)
------------------------------------------------------------------------- */

void PPPMStagger::compute_gf_denom()
{
  if (gf_b) memory->destroy(gf_b);
  memory->create(gf_b,order,"pppm:gf_b");

  int k,l,m;

  for (l = 1; l < order; l++) gf_b[l] = 0.0;
  gf_b[0] = 1.0;

  for (m = 1; m < order; m++) {
    for (l = m; l > 0; l--)
      gf_b[l] = 4.0 * (gf_b[l]*(l-m)*(l-m-0.5)-gf_b[l-1]*(l-m-1)*(l-m-1));
    gf_b[0] = 4.0 * (gf_b[0]*(l-m)*(l-m-0.5));
  }

  bigint ifact = 1;
  for (k = 1; k < 2*order; k++) ifact *= k;
  double gaminv = 1.0/ifact;
  for (l = 0; l < order; l++) gf_b[l] *= gaminv;
}

/* ----------------------------------------------------------------------
   pre-compute modified (Hockney-Eastwood) Coulomb Green's function
------------------------------------------------------------------------- */

void PPPMStagger::compute_gf_ik()
{
  const double * const prd = domain->prd;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double zprd_slab = zprd*slab_volfactor;
  const double unitkx = (MY_2PI/xprd);
  const double unitky = (MY_2PI/yprd);
  const double unitkz = (MY_2PI/zprd_slab);

  double snx,sny,snz;
  double cnx,cny,cnz;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,dot1,dot2;
  double numerator,denominator;
  double sqk;

  int k,l,m,n,nx,ny,nz,kper,lper,mper;

  const int nbx = static_cast<int> ((g_ewald*xprd/(MY_PI*nx_pppm)) *
                                    pow(-log(EPS_HOC),0.25));
  const int nby = static_cast<int> ((g_ewald*yprd/(MY_PI*ny_pppm)) *
                                    pow(-log(EPS_HOC),0.25));
  const int nbz = static_cast<int> ((g_ewald*zprd_slab/(MY_PI*nz_pppm)) *
                                    pow(-log(EPS_HOC),0.25));
  const int twoorder = 2*order;

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);
    snz = square(sin(0.5*unitkz*mper*zprd_slab/nz_pppm));
    cnz = cos(0.5*unitkz*mper*zprd_slab/nz_pppm);

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);
      sny = square(sin(0.5*unitky*lper*yprd/ny_pppm));
      cny = cos(0.5*unitky*lper*yprd/ny_pppm);

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);
        snx = square(sin(0.5*unitkx*kper*xprd/nx_pppm));
        cnx = cos(0.5*unitkx*kper*xprd/nx_pppm);

        sqk = square(unitkx*kper) + square(unitky*lper) + square(unitkz*mper);

        if (sqk != 0.0) {
          numerator = MY_4PI/sqk;
          denominator = 0.5*(gf_denom(snx,sny,snz) + gf_denom2(cnx,cny,cnz));
          sum1 = 0.0;

          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx*(kper+nx_pppm*nx);
            sx = exp(-0.25*square(qx/g_ewald));
            argx = 0.5*qx*xprd/nx_pppm;
            wx = powsinxx(argx,twoorder);

            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm*ny);
              sy = exp(-0.25*square(qy/g_ewald));
              argy = 0.5*qy*yprd/ny_pppm;
              wy = powsinxx(argy,twoorder);

              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm*nz);
                sz = exp(-0.25*square(qz/g_ewald));
                argz = 0.5*qz*zprd_slab/nz_pppm;
                wz = powsinxx(argz,twoorder);

                dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
                dot2 = qx*qx+qy*qy+qz*qz;
                sum1 += (dot1/dot2) * sx*sy*sz * wx*wy*wz;
              }
            }
          }
          greensfn[n++] = numerator*sum1/denominator;
        } else greensfn[n++] = 0.0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute optimized Green's function for energy calculation
------------------------------------------------------------------------- */

void PPPMStagger::compute_gf_ad()
{
  const double * const prd = domain->prd;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double zprd_slab = zprd*slab_volfactor;
  const double unitkx = (MY_2PI/xprd);
  const double unitky = (MY_2PI/yprd);
  const double unitkz = (MY_2PI/zprd_slab);

  double snx,sny,snz,sqk;
  double cnx,cny,cnz;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double numerator,denominator;
  int k,l,m,n,kper,lper,mper;

  const int twoorder = 2*order;

  for (int i = 0; i < 6; i++) sf_coeff[i] = 0.0;

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm*(2*m/nz_pppm);
    qz = unitkz*mper;
    snz = square(sin(0.5*qz*zprd_slab/nz_pppm));
    cnz = cos(0.5*qz*zprd_slab/nz_pppm);
    sz = exp(-0.25*square(qz/g_ewald));
    argz = 0.5*qz*zprd_slab/nz_pppm;
    wz = powsinxx(argz,twoorder);

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm*(2*l/ny_pppm);
      qy = unitky*lper;
      sny = square(sin(0.5*qy*yprd/ny_pppm));
      cny = cos(0.5*qy*yprd/ny_pppm);
      sy = exp(-0.25*square(qy/g_ewald));
      argy = 0.5*qy*yprd/ny_pppm;
      wy = powsinxx(argy,twoorder);

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm*(2*k/nx_pppm);
        qx = unitkx*kper;
        snx = square(sin(0.5*qx*xprd/nx_pppm));
        cnx = cos(0.5*qx*xprd/nx_pppm);
        sx = exp(-0.25*square(qx/g_ewald));
        argx = 0.5*qx*xprd/nx_pppm;
        wx = powsinxx(argx,twoorder);

        sqk = qx*qx + qy*qy + qz*qz;

        if (sqk != 0.0) {
          numerator = MY_4PI/sqk;
          denominator = 0.5*(gf_denom(snx,sny,snz) + gf_denom2(cnx,cny,cnz));
          greensfn[n] = numerator*sx*sy*sz*wx*wy*wz/denominator;
          sf_coeff[0] += sf_precoeff1[n]*greensfn[n];
          sf_coeff[1] += sf_precoeff2[n]*greensfn[n];
          sf_coeff[2] += sf_precoeff3[n]*greensfn[n];
          sf_coeff[3] += sf_precoeff4[n]*greensfn[n];
          sf_coeff[4] += sf_precoeff5[n]*greensfn[n];
          sf_coeff[5] += sf_precoeff6[n]*greensfn[n];
          n++;
        } else {
          greensfn[n] = 0.0;
          sf_coeff[0] += sf_precoeff1[n]*greensfn[n];
          sf_coeff[1] += sf_precoeff2[n]*greensfn[n];
          sf_coeff[2] += sf_precoeff3[n]*greensfn[n];
          sf_coeff[3] += sf_precoeff4[n]*greensfn[n];
          sf_coeff[4] += sf_precoeff5[n]*greensfn[n];
          sf_coeff[5] += sf_precoeff6[n]*greensfn[n];
          n++;
        }
      }
    }
  }

  // compute the coefficients for the self-force correction

  double prex, prey, prez;
  prex = prey = prez = MY_PI/volume;
  prex *= nx_pppm/xprd;
  prey *= ny_pppm/yprd;
  prez *= nz_pppm/zprd_slab;
  sf_coeff[0] *= prex;
  sf_coeff[1] *= prex*2;
  sf_coeff[2] *= prey;
  sf_coeff[3] *= prey*2;
  sf_coeff[4] *= prez;
  sf_coeff[5] *= prez*2;

  // communicate values with other procs

  double tmp[6];
  MPI_Allreduce(sf_coeff,tmp,6,MPI_DOUBLE,MPI_SUM,world);
  for (n = 0; n < 6; n++) sf_coeff[n] = tmp[n];
}

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array
------------------------------------------------------------------------- */

void PPPMStagger::particle_map()
{
  int nx,ny,nz;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {

    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    nx = static_cast<int> ((x[i][0]-boxlo[0])*delxinv+shift + stagger) - OFFSET;
    ny = static_cast<int> ((x[i][1]-boxlo[1])*delyinv+shift + stagger) - OFFSET;
    nz = static_cast<int> ((x[i][2]-boxlo[2])*delzinv+shift + stagger) - OFFSET;

    part2grid[i][0] = nx;
    part2grid[i][1] = ny;
    part2grid[i][2] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick

    if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
        ny+nlower < nylo_out || ny+nupper > nyhi_out ||
        nz+nlower < nzlo_out || nz+nupper > nzhi_out)
      flag = 1;
  }

  if (flag) error->one(FLERR,"Out of range atoms - cannot compute PPPM");
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void PPPMStagger::make_rho()
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;

  // clear 3d density array

  memset(&(density_brick[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv - stagger;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv - stagger;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv - stagger;

    compute_rho1d(dx,dy,dz);

    z0 = delvolinv * q[i];
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          density_brick[mz][my][mx] += x0*rho1d[0][l];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ik
------------------------------------------------------------------------- */

void PPPMStagger::fieldforce_ik()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ekx,eky,ekz;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv - stagger;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv - stagger;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv - stagger;

    compute_rho1d(dx,dy,dz);

    ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          ekx -= x0*vdx_brick[mz][my][mx];
          eky -= x0*vdy_brick[mz][my][mx];
          ekz -= x0*vdz_brick[mz][my][mx];
        }
      }
    }

    // convert E-field to force

    const double qfactor = qqrd2e * scale * q[i] / float(nstagger);
    f[i][0] += qfactor*ekx;
    f[i][1] += qfactor*eky;
    if (slabflag != 2) f[i][2] += qfactor*ekz;
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ad
------------------------------------------------------------------------- */

void PPPMStagger::fieldforce_ad()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz;
  FFT_SCALAR ekx,eky,ekz;
  double s1,s2,s3;
  double sf = 0.0;
  double *prd;

  prd = domain->prd;
  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];

  double hx_inv = nx_pppm/xprd;
  double hy_inv = ny_pppm/yprd;
  double hz_inv = nz_pppm/zprd;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv - stagger;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv - stagger;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv - stagger;

    compute_rho1d(dx,dy,dz);
    compute_drho1d(dx,dy,dz);

    ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          ekx += drho1d[0][l]*rho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          eky += rho1d[0][l]*drho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          ekz += rho1d[0][l]*rho1d[1][m]*drho1d[2][n]*u_brick[mz][my][mx];
        }
      }
    }
    ekx *= hx_inv;
    eky *= hy_inv;
    ekz *= hz_inv;

    // convert E-field to force and substract self forces

    const double qfactor = qqrd2e * scale / float(nstagger);

    s1 = x[i][0]*hx_inv + stagger;
    s2 = x[i][1]*hy_inv + stagger;
    s3 = x[i][2]*hz_inv + stagger;
    sf = sf_coeff[0]*sin(2*MY_PI*s1);
    sf += sf_coeff[1]*sin(4*MY_PI*s1);
    sf *= 2*q[i]*q[i];
    f[i][0] += qfactor*(ekx*q[i] - sf);

    sf = sf_coeff[2]*sin(2*MY_PI*s2);
    sf += sf_coeff[3]*sin(4*MY_PI*s2);
    sf *= 2*q[i]*q[i];
    f[i][1] += qfactor*(eky*q[i] - sf);

    sf = sf_coeff[4]*sin(2*MY_PI*s3);
    sf += sf_coeff[5]*sin(4*MY_PI*s3);
    sf *= 2*q[i]*q[i];
    if (slabflag != 2) f[i][2] += qfactor*(ekz*q[i] - sf);
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get per-atom energy/virial
------------------------------------------------------------------------- */

void PPPMStagger::fieldforce_peratom()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR u,v0,v1,v2,v3,v4,v5;

  // loop over my charges, interpolate from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv - stagger;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv - stagger;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv - stagger;

    compute_rho1d(dx,dy,dz);

    u = v0 = v1 = v2 = v3 = v4 = v5 = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          if (eflag_atom) u += x0*u_brick[mz][my][mx];
          if (vflag_atom) {
            v0 += x0*v0_brick[mz][my][mx];
            v1 += x0*v1_brick[mz][my][mx];
            v2 += x0*v2_brick[mz][my][mx];
            v3 += x0*v3_brick[mz][my][mx];
            v4 += x0*v4_brick[mz][my][mx];
            v5 += x0*v5_brick[mz][my][mx];
          }
        }
      }
    }

    if (eflag_atom) eatom[i] += q[i]*u/float(nstagger);
    if (vflag_atom) {
      vatom[i][0] += q[i]*v0/float(nstagger);
      vatom[i][1] += q[i]*v1/float(nstagger);
      vatom[i][2] += q[i]*v2/float(nstagger);
      vatom[i][3] += q[i]*v3/float(nstagger);
      vatom[i][4] += q[i]*v4/float(nstagger);
      vatom[i][5] += q[i]*v5/float(nstagger);
    }
  }
}

/* ----------------------------------------------------------------------
   perform and time the 1d FFTs required for N timesteps
------------------------------------------------------------------------- */

int PPPMStagger::timing_1d(int n, double &time1d)
{
  PPPM::timing_1d(n,time1d);
  time1d *= 2.0;

  if (differentiation_flag) return 2;
  return 4;
}

/* ----------------------------------------------------------------------
   perform and time the 3d FFTs required for N timesteps
------------------------------------------------------------------------- */

int PPPMStagger::timing_3d(int n, double &time3d)
{
  PPPM::timing_3d(n,time3d);
  time3d *= 2.0;

  if (differentiation_flag) return 2;
  return 4;
}
