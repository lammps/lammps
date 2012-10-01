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

KSpaceStyle(msm,MSM)

#else

#ifndef LMP_MSM_H
#define LMP_MSM_H

#include "lmptype.h"
#include "mpi.h"

#include "kspace.h"

namespace LAMMPS_NS {

class MSM : public KSpace {
 public:
  MSM(class LAMMPS *, int, char **);
  ~MSM();
  void init();
  void setup();
  void compute(int, int);

 protected:
  int me,nprocs;
  double precision;
  int nfactors;
  int *factors;
  double qsum,qsqsum,q2;
  double qqrd2e;
  double cutoff;
  double volume;
  double *delxinv,*delyinv,*delzinv,*delvolinv;

  int *nx_msm,*ny_msm,*nz_msm;
  int *nxlo_in,*nylo_in,*nzlo_in;
  int *nxhi_in,*nyhi_in,*nzhi_in;
  int *nxlo_in_d,*nylo_in_d,*nzlo_in_d;
  int *nxhi_in_d,*nyhi_in_d,*nzhi_in_d;
  int *nxlo_out,*nylo_out,*nzlo_out;
  int *nxhi_out,*nyhi_out,*nzhi_out;
  int *nxlo_ghost,*nxhi_ghost,*nylo_ghost;
  int *nyhi_ghost,*nzlo_ghost,*nzhi_ghost;
  int *ngrid;
  int nxlo_direct,nxhi_direct,nylo_direct;
  int nyhi_direct,nzlo_direct,nzhi_direct;
  int nmax_direct;
  int nlower,nupper;
  int nbuf,nbuf_peratom;
  int peratom_allocate_flag;
  int levels;


  double ****qgrid;
  double ****egrid;
  double ****fxgrid,****fygrid,****fzgrid;
  double ****v0grid,****v1grid,****v2grid;
  double ****v3grid,****v4grid,****v5grid;
  double **g_direct;
  double **dgx_direct,**dgy_direct,**dgz_direct;
  double **v0_direct,**v1_direct,**v2_direct;
  double **v3_direct,**v4_direct,**v5_direct;

  double *buf1,*buf2,*buf3,*buf4;

  double **phi1d,**dphi1d;

  int **part2grid;             // storage for particle -> grid mapping
  int nmax;

  int triclinic;               // domain settings, orthog or triclinic
  double *boxlo;

  void set_grid();
  double estimate_1d_error(double,double);
  double estimate_3d_error();
  double estimate_total_error();
  void allocate();
  void allocate_peratom();
  void deallocate();
  void deallocate_peratom();
  void allocate_levels();
  void deallocate_levels();
  int factorable(int,int&,int&);
  void compute_gf_denom();
  double gf_denom(double, double, double);
  void particle_map();
  void make_rho();
  void ghost_swap(int);
  void grid_swap(int,double*** &);
  void direct_ad(int);
  void direct(int);
  void direct_peratom(int);
  void restriction(int);
  void prolongation(int);
  void fillbrick_ad_peratom(int);
  void fillbrick(int);
  void fieldforce_ad();
  void fieldforce();
  void fieldforce_peratom();
  void procs2grid2d(int,int,int,int*,int*);
  void compute_phis_and_dphis(const double &, const double &, const double &);
  double compute_phi(const double &);
  double compute_dphi(const double &);
  void get_g_direct();
  void get_dg_direct();
  void get_virial_direct();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot (yet) use MSM with triclinic box

This feature is not yet supported.

E: Cannot (yet) use MSM with 2d simulation

This feature is not yet supported.

E: Kspace style requires atom attribute q

The atom style defined does not have these attributes.

E: Cannot (yet) use nonperiodic boundaries with MSM

This feature is not yet supported.

E: Cannot use slab correction with MSM

Slab correction can only be used with Ewald and PPPM, not MSM

E: MSM order cannot be < 4 or > than 10

This is a limitation of the MSM implementation in LAMMPS.

Currently the order may only range from 4 to 10

E: MSM order must be even

Currently the MSM order must be an even number

E: KSpace style is incompatible with Pair style

Setting a kspace style requires that a pair style with a long-range
Coulombic component be selected that is compatible with MSM.  Note
that TIP4P is not (yet) supported by MSM

E: Cannot use kspace solver on system with no charge

No atoms in system have a non-zero charge.

E: System is not charge neutral, net charge = %g

The total charge on all atoms on the system is not 0.0, which
is not valid for MSM.

W: MSM parallel communication error, try reducing accuracy or number of procs

Currently only nearest neighbor communication between processors is implemented in MSM.
If charge from an atom spans more than one processor domain this error will result.

E: MSM grid is too large

The global MSM grid is larger than OFFSET in one or more dimensions.
OFFSET is currently set to 16384.  You likely need to decrease the
requested accuracy.

E: KSpace accuracy must be > 0

The kspace accuracy designated in the input must be greater than zero.

E: Out of range atoms - cannot compute MSM

One or more atoms are attempting to map their charge to a MSM grid
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

*/
