/* -*- c++ -*- -----------------------------------------------------------
   LAMMPS 2003 (July 31) - Molecular Dynamics Simulator
   Sandia National Laboratories, www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   For more info, see the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------ */

#ifdef FIX_CLASS
// clang-format off
FixStyle(lb/fluid,FixLbFluid)
// clang-format on
#else

#ifndef LMP_FIX_LB_FLUID_H
#define LMP_FIX_LB_FLUID_H

#include "fix.h"

#if defined(MPI_STUBS)
#error "The LATBOLTZ package cannot be compiled in serial with MPI STUBS"
#endif

namespace LAMMPS_NS {

static const double kappa_lb = 0.0;

// 15-velocity lattice propogation vectors
static const int e15[15][3] = {{0, 0, 0},  {1, 0, 0},  {0, 1, 0},   {-1, 0, 0},   {0, -1, 0},
                               {0, 0, 1},  {0, 0, -1}, {1, 1, 1},   {-1, 1, 1},   {-1, -1, 1},
                               {1, -1, 1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, -1}, {1, -1, -1}};

// 15-velocity weights
static const double w_lb15[15] = {2. / 9.,  1. / 9.,  1. / 9.,  1. / 9.,  1. / 9.,
                                  1. / 9.,  1. / 9.,  1. / 72., 1. / 72., 1. / 72.,
                                  1. / 72., 1. / 72., 1. / 72., 1. / 72., 1. / 72.};

// 15-velocity normalizations
static const double Ng_lb15[15] = {1., 3., 3.,       3.,       9. / 2.,  9. / 2., 9. / 2., 9.,
                                   9., 9., 27. / 2., 27. / 2., 27. / 2., 9.,      1.};

// 15-velcity transformation matrix for f_i to moments
// clang-format off
static const double mg_lb15[15][15] = {
  {     1.,        1.,        1.,        1.,        1.,        1.,        1.,     1.,     1.,     1.,     1.,     1.,     1.,     1.,     1.},
  {     0.,        1.,        0.,       -1.,        0.,        0.,        0.,     1.,    -1.,    -1.,     1.,     1.,    -1.,    -1.,     1.},
  {     0.,        0.,        1.,        0.,       -1.,        0.,        0.,     1.,     1.,    -1.,    -1.,     1.,     1.,    -1.,    -1.},
  {     0.,        0.,        0.,        0.,        0.,        1.,       -1.,     1.,     1.,     1.,     1.,    -1.,    -1.,    -1.,    -1.},
  { -1./3.,     2./3.,    -1./3.,     2./3.,    -1./3.,    -1./3.,    -1./3.,  2./3.,  2./3.,  2./3.,  2./3.,  2./3.,  2./3.,  2./3.,  2./3.},
  { -1./3.,    -1./3.,     2./3.,    -1./3.,     2./3.,    -1./3.,    -1./3.,  2./3.,  2./3.,  2./3.,  2./3.,  2./3.,  2./3.,  2./3.,  2./3.},
  { -1./3.,    -1./3.,    -1./3.,    -1./3.,    -1./3.,     2./3.,     2./3.,  2./3.,  2./3.,  2./3.,  2./3.,  2./3.,  2./3.,  2./3.,  2./3.},
  {     0.,        0.,        0.,        0.,        0.,        0.,        0.,     1.,    -1.,     1.,    -1.,     1.,    -1.,     1.,    -1.},
  {     0.,        0.,        0.,        0.,        0.,        0.,        0.,     1.,     1.,    -1.,    -1.,    -1.,    -1.,     1.,     1.},
  {     0.,        0.,        0.,        0.,        0.,        0.,        0.,     1.,    -1.,    -1.,     1.,    -1.,     1.,     1.,    -1.},
  {     0.,        0.,    -1./3.,        0.,     1./3.,        0.,        0.,  2./3.,  2./3., -2./3., -2./3.,  2./3.,  2./3., -2./3., -2./3.},
  {     0.,        0.,        0.,        0.,        0.,    -1./3.,     1./3.,  2./3.,  2./3.,  2./3.,  2./3., -2./3., -2./3., -2./3., -2./3.},
  {     0.,    -1./3.,        0.,     1./3.,        0.,        0.,        0.,  2./3., -2./3., -2./3.,  2./3.,  2./3., -2./3., -2./3.,  2./3.},
  {     0.,        0.,        0.,        0.,        0.,        0.,        0.,     1.,    -1.,     1.,    -1.,    -1.,     1.,    -1.,     1.},
  {M_SQRT2,-M_SQRT1_2,-M_SQRT1_2,-M_SQRT1_2,-M_SQRT1_2,-M_SQRT1_2,-M_SQRT1_2,M_SQRT2,M_SQRT2,M_SQRT2,M_SQRT2,M_SQRT2,M_SQRT2,M_SQRT2,M_SQRT2}
};
// clang-format on

// 15-velocity opposite lattice directions for bounce-back, i.e. od[i] = j such that e15[j]=-e15[i]
static const int od[15] = {0, 3, 4, 1, 2, 6, 5, 13, 14, 11, 12, 9, 10, 7, 8};

// 15-velocity bounce-back list
// bbl[i][0] = number of bounce-back directions for orientation i
// bbl[i][j]...bbl[i][bbl[i][0]] directions that would be coming from inside the wall so need to come from bounce-back
// bbl[i][[bbl[i][0]+1]...bbl[i][16] directions where standard propogation can proceed (pointing into or along wall)
// inside edge has 1/4 inside domain, 3/4 outside domain
// outside edge has 3/4 outside domain, 1/4 inside domain
// Note: 1. This list is not exhaustive (eg. there should be 12 inside and 12 outside edges possible, it just covers cases
//          accessible in the pit routines.  Could be generalized to include other geometries
//       2. Need better labelling for corners (particularly in-out) that distinguishes different cases (e.g. 10 and 29 are NOT same, also 11,31)
// ori   wall normals (point into domain)
//  0    not relevent, ori==0 only for lattice type 0 (standard bulk fluid) and 2 (outside domain)
//  1    wall +x
//  2    wall +y
//  3    wall +z
//  4    outside edge +x,+z
//  5    inside edge  +x,+z
//  6    inside edge +y,+z
//  7    outside edge  +y,+z
//  8    inside edge -x,-y
//  9    inside edge -x,+y
// 10    in-out corner +x,+y,+z
// 11    in-out corner +x,-y,+z
// 12    inside corner -x,+y,+z
// 13    inside corner -x,-y,+z
// 14    wall -x
// 15    wall -y
// 16    wall -z
// 17    outside edge -x,+z
// 18    inside edge  -x,+z
// 19    inside edge -y,+z
// 20    outside edge  -y,+z
// 21    inside edge +x,-y
// 22    inside edge +x,+y
// 23    in-out corner -x,+y,+z
// 24    in-out corner -x,-y,+z
// 25    inside corner +x,+y,+z
// 26    inside corner +x,-y,+z
// 27    inside edge +y,-z
// 28    inside edge -y,-z
// 29    in-out corner +x,+y,+z
// 30    in-out corner +x,-y,+z
// 31    in-out corner -x,+y,+z
// 32    in-out corner -x,-y,+z
// clang-format off
static const int bbl[33][16] = {
  { 0,      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  { 5,  1,  7, 10, 11, 14,      0,  2,  3,  4,  5,  6,  8,  9, 12, 13},
  { 5,  2,  7,  8, 11, 12,      0,  1,  3,  4,  5,  6,  9, 10, 13, 14},
  { 5,  5,  7,  8,  9, 10,      0,  1,  2,  3,  4,  6, 11, 12, 13, 14},
  { 4,  7, 10,  1,  5,      0,  2,  3,  4,  6,  8,  9, 11, 12, 13, 14},
  { 8,  1,  5,  7, 10,  8,  9, 11, 14,      0,  2,  3,  4,  6, 12, 13},
  { 8,  2,  5,  7,  8,  9, 10, 11, 12,      0,  1,  3,  4,  6, 13, 14},
  { 4,  2,  5,  7,  8,      0,  1,  3,  4,  6,  9, 10, 11, 12, 13, 14},
  { 8,  3,  4,  9, 13,  8, 10, 12, 14,      0,  1,  2,  5,  6,  7, 11},
  { 8,  2,  3,  8, 12,  7,  9, 11, 13,      0,  1,  4,  5,  6, 10, 14},
  { 3,  7,  8, 10,      0,  1,  2,  3,  4,  5,  6,  9, 11, 12, 13, 14},
  { 3,  7,  9, 10,      0,  1,  2,  3,  4,  5,  6,  8, 11, 12, 13, 14},
  {10,  2,  3,  5,  8,  7,  9, 10, 11, 12, 13,      0,  1,  4,  6, 14},
  {10,  3,  4,  5,  9,  7,  8, 10, 12, 13, 14,      0,  1,  2,  6, 11},
  { 5,  3,  8,  9, 12, 13,      0,  1,  2,  4,  5,  6,  7, 10, 11, 14},
  { 5,  4,  9, 10, 13, 14,      0,  1,  2,  3,  5,  6,  7,  8, 11, 12},
  { 5,  6, 11, 12, 13, 14,      0,  1,  2,  3,  4,  5,  7,  8,  9, 10},
  { 4,  8,  9,  3,  5,      0,  1,  2,  4,  6,  7, 10, 11, 12, 13, 14},
  { 8,  3,  5,  8,  9,  7, 10, 12, 13,      0,  1,  2,  4,  6, 11, 14},
  { 8,  4,  5,  9, 10,  7,  8, 13, 14,      0,  1,  2,  3,  6, 11, 12},
  { 4,  4,  5,  9, 10,      0,  1,  2,  3,  6,  7,  8, 11, 12, 13, 14},
  { 8,  1,  4, 10, 14,  7,  9, 11, 13,      0,  2,  3,  5,  6,  8, 12},
  { 8,  1,  2,  7, 11,  8, 10, 12, 14,      0,  3,  4,  5,  6,  9, 13},
  { 3,  7,  8,  9,      0,  1,  2,  3,  4,  5,  6, 10, 11, 12, 13, 14},
  { 3,  8,  9, 10,      0,  1,  2,  3,  4,  5,  6,  7, 11, 12, 13, 14},
  {10,  1,  2,  5,  7,  8,  9, 10, 11, 12, 14,      0,  3,  4,  6, 13},
  {10,  1,  4,  5, 10,  7,  8,  9, 11, 13, 14,      0,  2,  3,  6, 12},
  { 8,  2,  6, 11, 12,  7,  8, 13, 14,      0,  1,  3,  4,  5,  9, 10},
  { 8,  4,  6, 13, 14,  9, 10, 11, 12,      0,  1,  2,  3,  5,  7,  8},
  { 6,  2,  7,  8, 11, 10, 12,      0,  1,  3,  4,  5,  6,  9, 13, 14},
  { 6,  4,  9, 10, 14,  7, 13,      0,  1,  2,  3,  5,  6,  8, 11, 12},
  { 6,  2,  7,  8, 12,  9, 11,      0,  1,  3,  4,  5,  6, 10, 13, 14},
  { 6,  4,  9, 10, 13,  8, 14,      0,  1,  2,  3,  5,  6,  7, 11, 12}
};
//clang-format on

// 19-velocity lattice propogation vectors
static const int e19[19][3] = {{0, 0, 0},   {1, 0, 0},  {0, 1, 0},  {-1, 0, 0}, {0, -1, 0},
                               {0, 0, 1},   {0, 0, -1}, {1, 1, 0},  {1, -1, 0}, {-1, 1, 0},
                               {-1, -1, 0}, {1, 0, 1},  {1, 0, -1}, {-1, 0, 1}, {-1, 0, -1},
                               {0, 1, 1},   {0, 1, -1}, {0, -1, 1}, {0, -1, -1}};

static const double w_lb19[19] = {1. / 3.,  1. / 18., 1. / 18., 1. / 18., 1. / 18.,
                                  1. / 18., 1. / 18., 1. / 36., 1. / 36., 1. / 36.,
                                  1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36.,
                                  1. / 36., 1. / 36., 1. / 36., 1. / 36.};

static const double Ng_lb19[19] = {1.,  3.,  3.,        3.,        9. / 2.,  9. / 2.,  9. / 2.,
                                   9.,  9.,  9.,        27. / 2.,  27. / 2., 27. / 2., 18.,
                                   18., 18., 162. / 7., 126. / 5., 30.};

// clang-format off
static const double mg_lb19[19][19] = {
  {    1.,     1.,     1.,     1.,     1.,     1.,     1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.},
  {    0.,     1.,     0.,    -1.,     0.,     0.,     0.,    1.,    1.,   -1.,   -1.,    1.,    1.,   -1.,   -1.,    0.,    0.,    0.,    0.},
  {    0.,     0.,     1.,     0.,    -1.,     0.,     0.,    1.,   -1.,    1.,   -1.,    0.,    0.,    0.,    0.,    1.,    1.,   -1.,   -1.},
  {    0.,     0.,     0.,     0.,     0.,     1.,    -1.,    0.,    0.,    0.,    0.,    1.,   -1.,    1.,   -1.,    1.,   -1.,    1.,   -1.},
  {-1./3.,  2./3., -1./3.,  2./3., -1./3., -1./3., -1./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3.,-1./3.,-1./3.,-1./3.,-1./3.},
  {-1./3., -1./3.,  2./3., -1./3.,  2./3., -1./3., -1./3., 2./3., 2./3., 2./3., 2./3.,-1./3.,-1./3.,-1./3.,-1./3., 2./3., 2./3., 2./3., 2./3.},
  {-1./3., -1./3., -1./3., -1./3., -1./3.,  2./3.,  2./3.,-1./3.,-1./3.,-1./3.,-1./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3.},
  {    0.,     0.,     0.,     0.,     0.,     0.,     0.,    1.,   -1.,   -1.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.},
  {    0.,     0.,     0.,     0.,     0.,     0.,     0.,    0.,    0.,    0.,    0.,    1.,   -1.,   -1.,    1.,    0.,    0.,    0.,    0.},
  {    0.,     0.,     0.,     0.,     0.,     0.,     0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,   -1.,   -1.,    1.},
  {    0., -1./3.,     0.,  1./3.,     0.,     0.,     0., 2./3., 2./3.,-2./3.,-2./3.,-1./3.,-1./3., 1./3., 1./3.,    0.,    0.,    0.,    0.},
  {    0.,     0., -1./3.,     0.,  1./3.,     0.,     0., 2./3.,-2./3., 2./3.,-2./3.,    0.,    0.,    0.,    0.,-1./3.,-1./3., 1./3., 1./3.},
  {    0.,     0.,     0.,     0.,     0., -1./3.,  1./3.,    0.,    0.,    0.,    0., 2./3.,-2./3., 2./3.,-2./3.,-1./3., 1./3.,-1./3., 1./3.},
  {    0.,   -0.5,     0.,    0.5,     0.,     0.,     0.,    0.,    0.,    0.,    0.,   0.5,   0.5,  -0.5,  -0.5,    0.,    0.,    0.,    0.},
  {    0.,     0.,     0.,     0.,     0.,   -0.5,    0.5,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,   0.5,  -0.5,   0.5,  -0.5},
  {    0.,     0.,   -0.5,     0.,    0.5,     0.,     0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,   0.5,   0.5,  -0.5,  -0.5},
  {1./18.,-5./18.,-5./18.,-5./18.,-5./18.,  2./9.,  2./9.,7./18.,7./18.,7./18.,7./18.,-1./9.,-1./9.,-1./9.,-1./9.,-1./9.,-1./9.,-1./9.,-1./9.},
  {1./14.,-5./14.,  1./7.,-5./14.,  1./7.,-3./14.,-3./14.,    0.,    0.,    0.,    0.,5./14.,5./14.,5./14.,5./14.,-1./7.,-1./7.,-1./7.,-1./7.},
  {1./10.,     0.,-3./10.,     0.,-3./10.,-3./10.,-3./10.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,3./10.,3./10.,3./10.,3./10.}
};
// clang-format on

class Site {
 public:
  int type;
  int orientation;
};

class FixLbFluid : public Fix {
  friend class FixLbMomentum;
  friend class FixLbViscous;

 public:
  FixLbFluid(class LAMMPS *, int, char **);
  ~FixLbFluid() override;
  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void setup(int) override;
  void pre_force(int) override;
  void post_force(int) override;
  void final_integrate() override;
  void end_of_step() override;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  double compute_scalar() override;
  double compute_vector(int) override;

  void dump(int);

 private:
  double viscosity, densityinit_real, a_0_real, T;
  int setdx, seta0, setdof;
  int numvel;

  double dm_lb, dx_lb, dt_lb;    // Lattice units for mass, distance, time.

  int Nbx, Nby, Nbz;             // Total # of x,y,z grid points.
  int subNbx, subNby, subNbz;    // # of x,y,z, grid points (including buffer)
                                 //   on local processor.
  int fluid_local_n0[3];         // Size of local including both lower and upper end points
  int fluid_global_o0[3];        // Offset of local in global from lower end point
  int fluid_global_n0[3];        // Size of global including both lower and upper end points

  int me, nprocs;           // MPI variables: processor ID, # of processors
  MPI_Datatype oneslice;    // MPI datatypes to pass arrays.
  MPI_Datatype passxu, passyu, passzu;
  MPI_Datatype passxf, passyf, passzf;
  MPI_Datatype passxrho, passyrho, passzrho;
  MPI_Datatype passxtemp, passytemp, passztemp;
  MPI_Datatype passxWtemp, passyWtemp, passzWtemp;
  MPI_Datatype fluid_density_2_mpitype;
  MPI_Datatype fluid_velocity_2_mpitype;

  MPI_Datatype realType3_mpitype;    // MPI type for a 3-vector

  double kB, densityinit, a_0;    // Boltzmann constant, initial density,
                                  //   and a_0 all in lattice units.
  double *Gamma;
  double **hydroF;
  double *massp;

  int dump_interval, dump_time_index;
  std::string dump_file_name_xdmf;
  std::string dump_file_name_raw;
  FILE *dump_file_handle_xdmf;
  MPI_File dump_file_handle_raw;
  MPI_Datatype dump_file_mpitype;

  int groupbit_viscouslb;

  double ***density_lb;    // fluid density
  double ****u_lb;         // fluid velocity
  double ****f_lb;         // distributions
  double ****fnew;         // used in the calculation of the new
                           //   distributions.
  double ****feq;          // equilibrium distributions

  double ****Ff, ***Wf;    // Force, total weight, from the MD particles on the fluid.
  double ****Fftempx, ***Wftempx;
  double ****Fftempy, ***Wftempy;
  double ****Fftempz, ***Wftempz;

  double tau;    // Lattice Boltzmann variables.
  double K_0;

  int step;

  int n_stencil = 2;    // Number of points for spread/interpolate stencil

  double bodyforcex, bodyforcey, bodyforcez;    // Body Forces acting on the fluid (default=0)
  double vwtp, vwbt;                            // Velocities of the z walls in the y
                                                //   direction. (must have fixed boundary
                                                //   conditions in z)

  int pressure;                  // pressure boundary conditions on/off
  double rhofactor;              // factor for density/pressure jump at boundary
  double rhoH, rhoL;             // density target for high/low side of jump
  double meanrho1, meanrhoLx;    // current densities at boundary on high/low side of jump

  int noisestress;    // 1 to include noise in the system,
                      //   0 otherwise.
  int lin_init;       // 1 to initialize with linear interpolation
                      //   between boundaries.
                      // 0 initialize to uniform density, 0.0 velocities.

  double namp, noisefactor;
  int seed;
  class RanMars *random;

  int readrestart;    // 1 to read in data from a restart file.
  MPI_File pFileRead;

  int printrestart;    // 1 to write data to a restart file.
  MPI_File pFileWrite;

  double timeEqb, timeUpdate, timePCalc, timefluidForce, timeCorrectU;

  double dof_lb = 0;

  int fixviscouslb;

  void rescale(void);
  void SetupBuffers(void);
  void InitializeFirstRun(void);
  void initializeLB(void);
  void initialize_feq(void);
  void (FixLbFluid::*equilibriumdist)(int, int, int, int, int, int);
  void equilibriumdist15(int, int, int, int, int, int);
  void equilibriumdist19(int, int, int, int, int, int);
  void parametercalc_part(int, int, int, int, int, int);
  void correctu_part(int, int, int, int, int, int);
  void parametercalc_full(void);
  void update_periodic(int, int, int, int, int, int);
  void correctu_full(void);

  void calc_mass_momentum(double &totalmass, double totalmomentum[3]);
  void calc_MPT(double &totalmass, double totalmomentum[3], double &Tave);

  void (FixLbFluid::*update_full)(void);
  void update_full15(void);
  void update_full19(void);

  void read_restartfile(void);
  void write_restartfile(void);

  void (FixLbFluid::*interpolate)(int, int);
  void (FixLbFluid::*interpolationweight)(int);
  void keys_interpolation(int, int);
  void keys_interpolationweight(int);
  void trilinear_interpolation(int, int);
  void trilinear_interpolationweight(int);
  void IBM3_interpolation(int, int);
  void IBM3_interpolationweight(int);

  void calc_fluidforceI(void);
  void calc_fluidforceII(void);
  void calc_fluidforceweight(void);

  int adjust_dof_fix();
  double dof_compute();

  /* nanopit parameters */
  int npits;           // number of nanopits
  int h_s;             // slit height
  int h_p;             // pit height
  int w_p;             // pit width
  int l_p;             // pit length
  int l_pp;            // distance between consecutive pits
  int l_e;             // length of end segments
  int sw;              // side walls on/off
  int openingsites;    // Number of active fluid sites at x=0

  Site ***sublattice, ***wholelattice;    // lattice geometry

  /* nanopit routines */
  void addslit(int &, int, int, int, int);
  void addpit(int &, int, int, int, int, int);
  void initializeGeometry(void);
  void initializeGlobalGeometry(void);
};
}    // namespace LAMMPS_NS
#endif
#endif
