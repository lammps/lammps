/* -*- c++ -*- -----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------ */

#ifdef FIX_CLASS
// clang-format off
FixStyle(lb/fluid,FixLbFluid);
// clang-format on
#else

#ifndef LMP_FIX_LB_FLUID_H
#define LMP_FIX_LB_FLUID_H

#include "fix.h"

#if defined(MPI_STUBS)
#error "The LATBOLTZ package cannot be compiled in serial with MPI STUBS"
#endif

namespace LAMMPS_NS {

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
  void post_force(int) override;
  void final_integrate() override;
  void end_of_step() override;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  double compute_scalar() override;
  double compute_vector(int) override;

  void dump(const bigint);

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

  int n_stencil;    // Number of points for spread/interpolate stencil

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

  double dof_lb;

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

  class Site ***sublattice, ***wholelattice;    // lattice geometry

  /* nanopit routines */
  void addslit(int &, int, int, int, int);
  void addpit(int &, int, int, int, int, int);
  void initializeGeometry(void);
  void initializeGlobalGeometry(void);
};
}    // namespace LAMMPS_NS
#endif
#endif
