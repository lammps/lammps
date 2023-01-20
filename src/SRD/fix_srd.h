/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(srd,FixSRD);
// clang-format on
#else

#ifndef LMP_FIX_SRD_H
#define LMP_FIX_SRD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSRD : public Fix {
 public:
  FixSRD(class LAMMPS *, int, char **);
  ~FixSRD() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void pre_neighbor() override;
  void post_force(int) override;
  double compute_vector(int) override;

  double memory_usage() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 private:
  int me, nprocs;
  int bigexist, biggroup, biggroupbit;
  int collidestyle, lamdaflag, overlap, insideflag, exactflag, maxbounceallow;
  int cubicflag, shiftuser, shiftseed, shiftflag, tstat;
  int rescale_rotate, rescale_collide;
  double gridsrd, gridsearch, lamda, radfactor, cubictol;
  int triclinic, change_size, change_shape, deformflag;

  double dt_big, dt_srd;
  double mass_big, mass_srd;
  double temperature_srd;
  double sigma;
  double srd_per_cell;
  double dmax, vmax, vmaxsq;
  double maxbigdiam, minbigdiam;
  double dist_ghost, dist_srd, dist_srd_reneigh;    // explained in code

  int wallexist, nwall, wallvarflag;
  class FixWallSRD *wallfix;
  int *wallwhich;
  double *xwall, *xwallhold, *vwall;
  double **fwall;
  double walltrigger;

  class AtomVecEllipsoid *avec_ellipsoid;
  class AtomVecLine *avec_line;
  class AtomVecTri *avec_tri;

  // for orthogonal box, these are in box units
  // for triclinic box, these are in lamda units

  double srdlo[3], srdhi[3];                    // SRDs must stay inside
  double srdlo_reneigh[3], srdhi_reneigh[3];    // SRDs trigger a reneigh

  int dimension;
  int initflag, setupflag, reneighflag;
  class RanMars *random;
  class RanPark *randomshift;

  // stats

  int ncheck, ncollide, ninside, nrescale, reneighcount;
  int nbounce, bouncemaxnum, bouncemax;
  int stats_flag;
  int srd_bin_count;
  double srd_bin_temp;
  double stats[12], stats_all[12];

  double **flocal;    // local ptrs to atom force and torque
  double **tlocal;

  // info to store for each owned and ghost big particle and wall

  struct Big {
    int index;               // local index of particle/wall
    int type;                // SPHERE or ELLIPSOID or LINE or TRI or WALL
    double radius, radsq;    // radius of sphere
    double aradsqinv;        // 3 ellipsoid radii
    double bradsqinv;
    double cradsqinv;
    double length;                 // length of line segment
    double normbody[3];            // normal of tri in body-frame
    double cutbinsq;               // add big to bin if within this distance
    double omega[3];               // current omega for sphere/ellipsoid/tri/line
    double ex[3], ey[3], ez[3];    // current orientation vecs for ellipsoid/tri
    double norm[3];                // current unit normal of tri in space-frame
    double theta;                  // current orientation of line
  };

  Big *biglist;      // list of info for each owned & ghost big and wall
  int torqueflag;    // 1 if any big particle is torqued

  // current size of particle-based arrays

  int nbig;      // # of owned/ghost big particles and walls
  int maxbig;    // max number of owned/ghost big particles and walls
  int nmax;      // max number of SRD particles

  // bins for SRD velocity remap, shifting and communication
  // binsize and inv are in lamda units for triclinic

  int nbins1, nbin1x, nbin1y, nbin1z;
  double binsize1x, binsize1y, binsize1z;
  double bininv1x, bininv1y, bininv1z;

  struct BinAve {
    int owner;           // 1 if I am owner of this bin, 0 if not
    int n;               // # of SRD particles in bin
    double xctr[3];      // center point of bin, only used for triclinic
    double vsum[3];      // sum of v components for SRD particles in bin
    double random;       // random value if I am owner
    double value[12];    // extra per-bin values
  };

  struct BinComm {
    int nsend, nrecv;            // # of bins to send/recv
    int sendproc, recvproc;      // who to send/recv to/from
    int *sendlist, *recvlist;    // list of bins to send/recv
  };

  struct BinShift {
    int commflag;                      // 1 if this shift requires any comm
    int nbins, nbinx, nbiny, nbinz;    // extent of my bins
    int maxbinsq, maxvbin;
    int binlo[3], binhi[3];    // extent of my bins in global array
    double corner[3];          // lower,left corner to offset from
                               // corner is in lamda units for triclinic
    BinAve *vbin;              // my bins
    BinComm bcomm[6];          // bin communication pattern for overlaps
  };
  BinShift shifts[2];    // 0 = no shift, 1 = shift

  int maxbin1;
  int *binhead;    // 1st SRD particle in each bin
  int *binnext;    // next SRD particle in same bin
  int maxbuf;
  double *sbuf1, *sbuf2;    // buffers for send/recv of velocity bin data
  double *rbuf1, *rbuf2;

  // bins and stencil for collision searching for SRDs & big particles

  int nbins2, nbin2x, nbin2y, nbin2z;
  int maxbin2;
  double binsize2x, binsize2y, binsize2z;
  double bininv2x, bininv2y, bininv2z;
  double xblo2, yblo2, zblo2;

  int *nbinbig;      // # of big particles overlapping each bin
  int **binbig;      // indices of big particles overlapping each bin
  int *binsrd;       // which bin each SRD particle is in
  int nstencil;      // # of bins in stencil
  int maxstencil;    // max # of bins stencil array can hold
  int **stencil;     // list of 3d bin offsets a big particle can overlap

  // persistent data for line/tri collision calculations

  double tfraction, theta0, theta1;
  double xs0[3], xs1[3], xsc[3];
  double xb0[3], xb1[3], xbc[3];
  double nbc[3];

  // shared data for triangle collision calculations

  // private functions

  void reset_velocities();
  void vbin_comm(int);
  void vbin_pack(BinAve *, int, int *, double *);
  void vbin_unpack(double *, BinAve *, int, int *);

  void xbin_comm(int, int);
  void xbin_pack(BinAve *, int, int *, double *, int);
  void xbin_unpack(double *, BinAve *, int, int *, int);

  void collisions_single();
  void collisions_multi();

  int inside_sphere(double *, double *, Big *);
  int inside_ellipsoid(double *, double *, Big *);
  int inside_line(double *, double *, double *, double *, Big *, double);
  int inside_tri(double *, double *, double *, double *, Big *, double);
  int inside_wall(double *, int);

  double collision_sphere_exact(double *, double *, double *, double *, Big *, double *, double *,
                                double *);
  void collision_sphere_inexact(double *, double *, Big *, double *, double *, double *);
  double collision_ellipsoid_exact(double *, double *, double *, double *, Big *, double *,
                                   double *, double *);
  void collision_ellipsoid_inexact(double *, double *, Big *, double *, double *, double *);
  double collision_line_exact(double *, double *, double *, double *, Big *, double, double *,
                              double *, double *);
  double collision_tri_exact(double *, double *, double *, double *, Big *, double, double *,
                             double *, double *);
  double collision_wall_exact(double *, int, double *, double *, double *, double *);
  void collision_wall_inexact(double *, int, double *, double *, double *);

  void slip(double *, double *, double *, Big *, double *, double *, double *);
  void slip_wall(double *, int, double *, double *);
  void noslip(double *, double *, double *, Big *, int, double *, double *, double *);

  void force_torque(double *, double *, double *, double *, double *, double *);
  void force_wall(double *, double *, int);

  int update_srd(int, double, double *, double *, double *, double *);

  void parameterize();
  void setup_bounds();
  void setup_velocity_bins();
  void setup_velocity_shift(int, int);
  void setup_search_bins();
  void setup_search_stencil();
  void big_static();
  void big_dynamic();

  double point_bin_distance(double *, int, int, int);
  double bin_bin_distance(int, int, int);
  void velocity_stats(int);

  double newton_raphson(double, double);
  void lineside(double, double &, double &);
  void triside(double, double &, double &);

  double distance(int, int);
  void print_collision(int, int, int, double, double, double *, double *, double *, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
