/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS

KSpaceStyle(pppm/old,PPPMOld)

#else

#ifndef LMP_PPPM_OLD_H
#define LMP_PPPM_OLD_H

#include "lmptype.h"
#include "mpi.h"

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_CFLOAT
#else
typedef double FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

#include "kspace.h"

namespace LAMMPS_NS {

class PPPMOld : public KSpace {
 public:
  PPPMOld(class LAMMPS *, int, char **);
  virtual ~PPPMOld();
  virtual void init();
  virtual void setup();
  virtual void compute(int, int);
  virtual int timing_1d(int, double &);
  virtual int timing_3d(int, double &);
  virtual double memory_usage();

  virtual void compute_group_group(int, int, int);

 protected:
  int me,nprocs;
  int nfactors;
  int *factors;
  double cutoff;
  double volume;
  double delxinv,delyinv,delzinv,delvolinv;
  double shift,shiftone;
  int peratom_allocate_flag;

  int nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in;
  int nxlo_out,nylo_out,nzlo_out,nxhi_out,nyhi_out,nzhi_out;
  int nxlo_ghost,nxhi_ghost,nylo_ghost,nyhi_ghost,nzlo_ghost,nzhi_ghost;
  int nxlo_fft,nylo_fft,nzlo_fft,nxhi_fft,nyhi_fft,nzhi_fft;
  int nlower,nupper;
  int ngrid,nfft,nfft_both;
  int nbuf,nbuf_peratom;

  FFT_SCALAR ***density_brick;
  FFT_SCALAR ***vdx_brick,***vdy_brick,***vdz_brick;
  FFT_SCALAR ***u_brick;
  FFT_SCALAR ***v0_brick,***v1_brick,***v2_brick;
  FFT_SCALAR ***v3_brick,***v4_brick,***v5_brick;
  double *greensfn;
  double **vg;
  double *fkx,*fky,*fkz;
  FFT_SCALAR *density_fft;
  FFT_SCALAR *work1,*work2;
  FFT_SCALAR *buf1,*buf2,*buf3,*buf4;

  double *gf_b;
  FFT_SCALAR **rho1d,**rho_coeff;

  // group-group interactions

  int group_allocate_flag;
  FFT_SCALAR ***density_A_brick,***density_B_brick;
  FFT_SCALAR *density_A_fft,*density_B_fft;


  class FFT3d *fft1,*fft2;
  class Remap *remap;

  int **part2grid;             // storage for particle -> grid mapping
  int nmax;

  int triclinic;               // domain settings, orthog or triclinic
  double *boxlo;
                               // TIP4P settings
  int typeH,typeO;             // atom types of TIP4P water H and O atoms
  double qdist;                // distance from O site to negative charge
  double alpha;                // geometric factor

  void set_grid();
  virtual void allocate();
  virtual void allocate_peratom();
  virtual void deallocate();
  virtual void deallocate_peratom();
  int factorable(int);
  double rms(double, double, bigint, double, double **);
  double diffpr(double, double, double, double, double **);
  void compute_gf_denom();

  virtual void particle_map();
  virtual void make_rho();
  virtual void brick2fft();
  virtual void fillbrick();
  virtual void fillbrick_peratom();
  virtual void poisson(int,int);
  virtual void poisson_peratom();
  virtual void fieldforce();
  virtual void fieldforce_peratom();
  void procs2grid2d(int,int,int,int *, int*);
  void compute_rho1d(const FFT_SCALAR &, const FFT_SCALAR &,
                     const FFT_SCALAR &);
  void compute_rho_coeff();
  void slabcorr();

  // group-group interactions

  virtual void allocate_groups();
  virtual void deallocate_groups();
  virtual void make_rho_groups(int, int, int);
  virtual void poisson_groups(int);

/* ----------------------------------------------------------------------
   denominator for Hockney-Eastwood Green's function
     of x,y,z = sin(kx*deltax/2), etc

            inf                 n-1
   S(n,k) = Sum  W(k+pi*j)**2 = Sum b(l)*(z*z)**l
           j=-inf               l=0

          = -(z*z)**n /(2n-1)! * (d/dx)**(2n-1) cot(x)  at z = sin(x)
   gf_b = denominator expansion coeffs
------------------------------------------------------------------------- */

  inline double gf_denom(const double &x, const double &y,
                         const double &z) const {
    double sx,sy,sz;
    sz = sy = sx = 0.0;
    for (int l = order-1; l >= 0; l--) {
      sx = gf_b[l] + sx*x;
      sy = gf_b[l] + sy*y;
      sz = gf_b[l] + sz*z;
    }
    double s = sx*sy*sz;
    return s*s;
  };
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use PPPM with 2d simulation

The kspace style pppm cannot be used in 2d simulations.  You can use
2d PPPM in a 3d simulation; see the kspace_modify command.

E: Kspace style requires atom attribute q

The atom style defined does not have these attributes.

E: Cannot use nonperiodic boundaries with PPPM

For kspace style pppm, all 3 dimensions must have periodic boundaries
unless you use the kspace_modify command to define a 2d slab with a
non-periodic z dimension.

E: Incorrect boundaries with slab PPPM

Must have periodic x,y dimensions and non-periodic z dimension to use
2d slab option with PPPM.

E: PPPM order cannot be < 2 or > than %d

This is a limitation of the PPPM implementation in LAMMPS.

E: KSpace style is incompatible with Pair style

Setting a kspace style requires that a pair style with a long-range
Coulombic or dispersion component be used.

E: Bond and angle potentials must be defined for TIP4P

Cannot use TIP4P pair potential unless bond and angle potentials
are defined.

E: Bad TIP4P angle type for PPPM/TIP4P

Specified angle type is not valid.

E: Bad TIP4P bond type for PPPM/TIP4P

Specified bond type is not valid.

E: Cannot use kspace solver on system with no charge

No atoms in system have a non-zero charge.

W: System is not charge neutral, net charge = %g

The total charge on all atoms on the system is not 0.0, which
is not valid for the long-range Coulombic solvers.

W: Reducing PPPM order b/c stencil extends beyond neighbor processor

This may lead to a larger grid than desired.  See the kspace_modify overlap
command to prevent changing of the PPPM order.

E: PPPM grid is too large

The global PPPM grid is larger than OFFSET in one or more dimensions.
OFFSET is currently set to 4096.  You likely need to decrease the
requested accuracy.

E: PPPM order has been reduced to 0

The auto-adjust of the order failed.  You will need to
set the grid size and order directly via kspace_modify.

E: KSpace accuracy must be > 0

The kspace accuracy designated in the input must be greater than zero.

E: Cannot compute PPPM G

The Ewald factor could not be computed for the current choice of
grid size, cutoff, accuracy.

E: Out of range atoms - cannot compute PPPM

One or more atoms are attempting to map their charge to a PPPM grid
point that is not owned by a processor.  This is likely for one of two
reasons, both of them bad.  First, it may mean that an atom near the
boundary of a processor's sub-domain has moved more than 1/2 the
"neighbor skin distance"_neighbor.html without neighbor lists being
rebuilt and atoms being migrated to new processors.  This also means
you may be missing pairwise interactions that need to be computed.
The solution is to change the re-neighboring criteria via the
"neigh_modify"_neigh_modify command.  The safest settings are "delay 0
every 1 check yes".  Second, it may mean that an atom has moved far
outside a processor's sub-domain or even the entire simulation box.
This indicates bad physics, e.g. due to highly overlapping atoms, too
large a timestep, etc.

E: Cannot (yet) use K-space slab correction with compute group/group

This option is not yet supported.

*/
