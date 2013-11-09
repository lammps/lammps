/* -----------------------------------------------------------------------
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

FixStyle(lb/fluid,FixLbFluid)

#else 

#ifndef LMP_FIX_LB_FLUID_H
#define LMP_FIX_LB_FLUID_H

#include "fix.h"

#if defined(MPI_STUBS)
#error "The USER-LB package cannot be compiled in serial with MPI STUBS"
#endif

namespace LAMMPS_NS {

  class FixLbFluid : public Fix {
    friend class FixLbMomentum;
    friend class FixLbRigidPCSphere;
    friend class FixLbPC;
    friend class FixLbViscous;

public:
    FixLbFluid(class LAMMPS *, int, char **);
    ~FixLbFluid();
    int setmask();
    void init();
    void initial_integrate(int);
    void setup(int);
    void post_force(int);
    void end_of_step();

    void grow_arrays(int);
    void copy_arrays(int, int, int);
    int pack_exchange(int, double *);
    int unpack_exchange(int, double *);

  private:
    double viscosity,densityinit_real,a_0_real,T;
    int setdx,seta0;
    int numvel;

    double dm_lb,dx_lb,dt_lb;                        // Lattice units for mass, distance, time.

    int Nbx,Nby,Nbz;                                 // Total # of x,y,z grid points.
    int subNbx,subNby,subNbz;                        // # of x,y,z, grid points (including buffer) 
                                                     //   on local processor.
    int me, nprocs;                                  // MPI variables: processor ID, # of processors
    MPI_Datatype oneslice;                           // MPI datatypes to pass arrays.
    MPI_Datatype passxu,passyu,passzu;
    MPI_Datatype passxf,passyf,passzf;
    MPI_Datatype passxrho,passyrho,passzrho;
    MPI_Datatype passxtemp,passytemp,passztemp;

    double kB,densityinit,a_0;                       // Boltzmann constant, initial density, 
                                                     //   and a_0 all in lattice units.
    double *Gamma;      
    double *NodeArea;
    int setGamma,setArea;
    double **hydroF;

    int groupbit_viscouslb, groupbit_pc, groupbit_rigid_pc_sphere;

    double ***density_lb;                            // fluid density
    double ****u_lb;                                 // fluid velocity
    double ****f_lb;                                 // distributions
    double ****fnew;                                 // used in the calculation of the new 
                                                     //   distributions.
    double ****feq;                                  // equilibrium distributions
    double ****feqold;                               // equilibrium distributions from previous 
                                                     //   timestep

    double ****feqn;                                 // equilibrium distributions without noise.
    double ****feqoldn;                              // equilibrium distributions from previous 
                                                     //   timestep without noise.
    double ****Ff;                                   // Force from the MD particles on the fluid.
    double ****Fftempx;
    double ****Fftempy;
    double ****Fftempz;

    double *Ng_lb;                                   // Lattice Boltzmann variables.  
    double *w_lb;
    double **mg_lb;
    int **e;
    double tau;
    double expminusdtovertau;
    double Dcoeff;
    double K_0;
    double dtoverdtcollision;

    int step;

    double ****buf;                                  // arrays used to output data.
    double ****buf2;
    double ****altogether;
    double ****altogether2;

    double bodyforcex,bodyforcey,bodyforcez;         // Body Forces acting on the fluid (default=0)
    double vwtp,vwbt;                                // Velocities of the z walls in the y 
                                                     //   direction. (must have fixed boundary 
                                                     //   conditions in z)

    int noisestress;                                 // 1 to include noise in the system, 
                                                     //   0 otherwise.
    double namp,noisefactor;
    int seed;
    class RanMars *random;

    int force_diagnostic;                            // 1 to print out the force action on a group 
                                                     //   of particles, 0 otherwise.
    int igroupforce;                                 // the group for which the force is to be 
                                                     //   printed.

    int typeLB;                                      
 
    int trilinear_stencil;                           // 1 to use the trilinear stencil, 0 to use the 
                                                     //   peskin stencil.

    int readrestart;                                 // 1 to read in data from a restart file.
    MPI_File pFileRead;                                

    int printrestart;                                // 1 to write data to a restart file.
    MPI_File pFileWrite;

    int printfluid;
    int fixviscouslb;

    void rescale(void);
    void (FixLbFluid::*initializeLB)(void);
    void initializeLB15(void);
    void initializeLB19(void);
    void initialize_feq(void);
    void (FixLbFluid::*equilibriumdist)(int,int,int,int,int,int);
    void equilibriumdist15(int,int,int,int,int,int);
    void equilibriumdist19(int,int,int,int,int,int);
    void parametercalc_part(int,int,int,int,int,int);
    void parametercalc_full(void);
    void update_periodic(int,int,int,int,int,int);
    void (FixLbFluid::*update_full)(void);
    void update_full15(void);
    void update_full19(void);
    void streamout(void);
    void read_restartfile(void);
    void write_restartfile(void);
    void peskin_interpolation(int);
    void trilinear_interpolation(int);
    void calc_fluidforce(void);

  };
}
#endif
#endif
    
