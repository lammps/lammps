// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#include "pppm_dielectric.h"

#include "atom.h"
#include "atom_vec_dielectric.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fft3d_wrap.h"
#include "force.h"
#include "grid3d.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define SMALL 0.00001

enum {REVERSE_RHO};
enum {FORWARD_IK,FORWARD_AD,FORWARD_IK_PERATOM,FORWARD_AD_PERATOM};

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ---------------------------------------------------------------------- */

PPPMDielectric::PPPMDielectric(LAMMPS *_lmp) : PPPM(_lmp)
{
  group_group_enable = 0;

  efield = nullptr;
  phi = nullptr;
  potflag = 0;
  use_qscaled = true;

  // no warnings about non-neutral systems from qsum_qsq()
  warn_nonneutral = 2;

  avec = dynamic_cast<AtomVecDielectric *>(atom->style_match("dielectric"));
  if (!avec) error->all(FLERR,"pppm/dielectric requires atom style dielectric");
}

/* ---------------------------------------------------------------------- */

PPPMDielectric::~PPPMDielectric()
{
  memory->destroy(efield);
  memory->destroy(phi);
}

/* ----------------------------------------------------------------------
   compute the PPPMDielectric long-range force, energy, virial
------------------------------------------------------------------------- */

void PPPMDielectric::compute(int eflag, int vflag)
{
  int i,j;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  ev_init(eflag,vflag);

  if (evflag_atom && !peratom_allocate_flag) allocate_peratom();

  // if atom count has changed, update qsum and qsqsum

  if (atom->natoms != natoms_original) {
    qsum_qsq();
    natoms_original = atom->natoms;
  }

  // recompute the average epsilon of all the atoms

  compute_ave_epsilon();

  // return if there are no charges or dipoles

  if (qsqsum == 0.0) return;

  // convert atoms from box to lamda coords

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(part2grid);
    memory->destroy(efield);
    memory->destroy(phi);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm/dielectric:part2grid");
    memory->create(efield,nmax,3,"pppm/dielectric:efield");
    memory->create(phi,nmax,"pppm/dielectric:phi");
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid

  particle_map();
  make_rho();

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  gc->reverse_comm(Grid3d::KSPACE,this,REVERSE_RHO,1,sizeof(FFT_SCALAR),
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  brick2fft();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  // also performs per-atom calculations via poisson_peratom()

  double energy_before_poisson = energy;
  poisson();

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  if (differentiation_flag == 1)
    gc->forward_comm(Grid3d::KSPACE,this,FORWARD_AD,1,sizeof(FFT_SCALAR),
                     gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  else
    gc->forward_comm(Grid3d::KSPACE,this,FORWARD_IK,3,sizeof(FFT_SCALAR),
                     gc_buf1,gc_buf2,MPI_FFT_SCALAR);

  // extra per-atom energy/virial communication

  if (evflag_atom) {
    if (differentiation_flag == 1 && vflag_atom)
      gc->forward_comm(Grid3d::KSPACE,this,FORWARD_AD_PERATOM,6,sizeof(FFT_SCALAR),
                       gc_buf1,gc_buf2,MPI_FFT_SCALAR);
    else if (differentiation_flag == 0)
      gc->forward_comm(Grid3d::KSPACE,this,FORWARD_IK_PERATOM,7,sizeof(FFT_SCALAR),
                       gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  }

  // calculate the force on my particles

  fieldforce();

  // extra per-atom energy/virial communication

  if (evflag_atom) fieldforce_peratom();

  // sum global energy across procs and add in volume-dependent term
  // NOTE: electrostatic energy is not linearly dependent on charge density (unlike forces)
  //       recall that we are using atom->q_scaled for make_rho()
  //       need to switch to atom->q to compute energy (elong)
  //       also need to use average epsilon (assuming that global dielectric = 1, not set in the input script)

  const double qscale = qqrd2e * scale;

  if (eflag_global) {

    energy = energy_before_poisson;

    // switch to unscaled charges to find charge density

    use_qscaled = false;

    // redo the charge density

    make_rho();

    // communicate for charge density

    gc->reverse_comm(Grid3d::KSPACE,this,REVERSE_RHO,1,sizeof(FFT_SCALAR),
                     gc_buf1,gc_buf2,MPI_FFT_SCALAR);

    brick2fft();

    // compute electrostatic energy with the unscaled charges and average epsilon

    poisson();

    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;
    energy *= 0.5*volume;
      energy -= g_ewald*qsqsum/MY_PIS +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= (qscale/epsilon_ave);

    // revert to qscaled charges (for force in the next time step)

    use_qscaled = true;
  }

  // sum global virial across procs

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*volume*virial_all[i];
  }

  // per-atom energy/virial
  // energy includes self-energy correction
  // ntotal accounts for TIP4P tallying eatom/vatom for ghost atoms

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
   compute the average dielectric constant of all the atoms
   NOTE: for dielectric use cases
------------------------------------------------------------------------- */

void PPPMDielectric::compute_ave_epsilon()
{
  const double * const epsilon = atom->epsilon;
  const int nlocal = atom->nlocal;
  double epsilon_local(0.0);

#if defined(_OPENMP)
#pragma omp parallel for default(shared) reduction(+:epsilon_local)
#endif
  for (int i = 0; i < nlocal; i++) {
    epsilon_local += epsilon[i];
  }

  MPI_Allreduce(&epsilon_local,&epsilon_ave,1,MPI_DOUBLE,MPI_SUM,world);
  epsilon_ave /= (double)atom->natoms;
}

/* ----------------------------------------------------------------------
   compute qsum,qsqsum,q2 and give error/warning if not charge neutral
   called initially, when particle count changes, when charges are changed
------------------------------------------------------------------------- */

void PPPMDielectric::qsum_qsq(int warning_flag)
{
  const double * const q = atom->q;
  const double * const epsilon = atom->epsilon;
  const int nlocal = atom->nlocal;
  double qsum_local(0.0), qsqsum_local(0.0), qsqsume_local(0.0);
  double qsqsume;

#if defined(_OPENMP)
#pragma omp parallel for default(shared) reduction(+:qsum_local,qsqsum_local)
#endif
  for (int i = 0; i < nlocal; i++) {
    qsum_local += q[i];
    qsqsum_local += q[i]*q[i];
    qsqsume_local += q[i]*q[i]/epsilon[i];
  }

  MPI_Allreduce(&qsum_local,&qsum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&qsqsum_local,&qsqsum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&qsqsume_local,&qsqsume,1,MPI_DOUBLE,MPI_SUM,world);

  if ((qsqsum == 0.0) && (comm->me == 0) && warn_nocharge && warning_flag) {
    error->warning(FLERR,"Using kspace solver on system with no charge");
    warn_nocharge = 0;
  }

  // q2 is used to compute the mesh spacing, here using qsqsume to match with regular pppm
  q2 = qsqsume * force->qqrd2e; //q2 = qsqsum * force->qqrd2e;

  // not yet sure of the correction needed for non-neutral systems
  // so issue warning or error

  if (fabs(qsum) > SMALL) {
    std::string message = fmt::format("System is not charge neutral, net "
                                      "charge = {:.8}",qsum);
    if (!warn_nonneutral) error->all(FLERR,message);
    if (warn_nonneutral == 1 && comm->me == 0) error->warning(FLERR,message);
    warn_nonneutral = 2;
  }
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
   NOTE: compute charge density using q_scaled if use_qscaled==true
                        else using unscaled charge values
------------------------------------------------------------------------- */

void PPPMDielectric::make_rho()
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

  double *q = atom->q_scaled;
  if (!use_qscaled) q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

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

void PPPMDielectric::fieldforce_ik()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ekx,eky,ekz,u;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q_scaled;
  double **x = atom->x;
  double **f = atom->f;
  double *eps = atom->epsilon;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    u = ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          if (potflag) u += x0*u_brick[mz][my][mx];
          ekx -= x0*vdx_brick[mz][my][mx];
          eky -= x0*vdy_brick[mz][my][mx];
          ekz -= x0*vdz_brick[mz][my][mx];
        }
      }
    }

    // electrostatic potential

    if (potflag) phi[i] = u;

    // convert E-field to force
    const double efactor = scale * eps[i];
    efield[i][0] = efactor*ekx;
    efield[i][1] = efactor*eky;
    efield[i][2] = efactor*ekz;

    const double qfactor = qqrd2e * efactor * q[i];
    f[i][0] += qfactor*ekx;
    f[i][1] += qfactor*eky;
    if (slabflag != 2) f[i][2] += qfactor*ekz;
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ad
------------------------------------------------------------------------- */

void PPPMDielectric::fieldforce_ad()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz;
  FFT_SCALAR ekx,eky,ekz,u;

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

  double *q = atom->q_scaled;
  double **x = atom->x;
  double **f = atom->f;
  double *eps = atom->epsilon;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);
    compute_drho1d(dx,dy,dz);

    u = ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          u += rho1d[0][l]*rho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          ekx += drho1d[0][l]*rho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          eky += rho1d[0][l]*drho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          ekz += rho1d[0][l]*rho1d[1][m]*drho1d[2][n]*u_brick[mz][my][mx];
        }
      }
    }
    ekx *= hx_inv;
    eky *= hy_inv;
    ekz *= hz_inv;

    // electrical potential

    if (potflag) phi[i] = u;

    // convert E-field to force and substract self forces

    const double qfactor = qqrd2e * scale * eps[i];
    double qtmp = q[i];

    s1 = x[i][0]*hx_inv;
    s2 = x[i][1]*hy_inv;
    s3 = x[i][2]*hz_inv;
    sf = sf_coeff[0]*sin(2*MY_PI*s1);
    sf += sf_coeff[1]*sin(4*MY_PI*s1);
    sf *= 2*qtmp*qtmp;
    f[i][0] += qfactor*(ekx*qtmp - sf);
    if (qtmp != 0) efield[i][0] = qfactor*(ekx - sf/qtmp);
    else efield[i][0] = qfactor*ekx;

    sf = sf_coeff[2]*sin(2*MY_PI*s2);
    sf += sf_coeff[3]*sin(4*MY_PI*s2);
    sf *= 2*qtmp*qtmp;
    f[i][1] += qfactor*(eky*qtmp - sf);
    if (qtmp != 0) efield[i][1] = qfactor*(eky - sf/qtmp);
    else efield[i][1] = qfactor*eky;

    sf = sf_coeff[4]*sin(2*MY_PI*s3);
    sf += sf_coeff[5]*sin(4*MY_PI*s3);
    sf *= 2*qtmp*qtmp;
    if (slabflag != 2) {
      f[i][2] += qfactor*(ekz*qtmp - sf);
      if (qtmp != 0) efield[i][2] = qfactor*(ekz - sf/qtmp);
      else efield[i][2] = qfactor*ekz;
    }

  }
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane. Also
   extended to non-neutral systems (J. Chem. Phys. 131, 094107).
------------------------------------------------------------------------- */

void PPPMDielectric::slabcorr()
{
  // compute local contribution to global dipole moment

  double *q = atom->q_scaled;
  double **x = atom->x;
  double *eps = atom->epsilon;
  double zprd_slab = domain->zprd*slab_volfactor;
  int nlocal = atom->nlocal;

  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) dipole += q[i]*x[i][2];

  // sum local contributions to get global dipole moment

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // need to make non-neutral systems and/or
  //  per-atom energy translationally invariant

  double dipole_r2 = 0.0;
  if (eflag_atom || fabs(qsum) > SMALL) {
    for (int i = 0; i < nlocal; i++)
      dipole_r2 += q[i]*x[i][2]*x[i][2];

    // sum local contributions

    double tmp;
    MPI_Allreduce(&dipole_r2,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    dipole_r2 = tmp;
  }

  // compute corrections

  const double e_slabcorr = MY_2PI*(dipole_all*dipole_all -
    qsum*dipole_r2 - qsum*qsum*zprd_slab*zprd_slab/12.0)/volume;
  const double qscale = qqrd2e * scale;

  if (eflag_global) energy += qscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    double efact = qscale * MY_2PI/volume;
    for (int i = 0; i < nlocal; i++)
      eatom[i] += efact * eps[i]*q[i]*(x[i][2]*dipole_all - 0.5*(dipole_r2 +
        qsum*x[i][2]*x[i][2]) - qsum*zprd_slab*zprd_slab/12.0);
  }

  // add on force corrections

  double ffact = qscale * (-4.0*MY_PI/volume);
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) {
    f[i][2] += ffact * eps[i]*q[i]*(dipole_all - qsum*x[i][2]);
    efield[i][2] += ffact * eps[i]*(dipole_all - qsum*x[i][2]);
  }
}
