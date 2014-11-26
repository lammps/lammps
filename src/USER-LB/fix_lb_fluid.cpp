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
   Contributing authors: Frances Mackay, Santtu Ollila, Colin Denniston (UWO)
------------------------------------------------------------------------- */

#include "fix_lb_fluid.h"
#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "atom.h"
#include <iostream>
#include <iomanip>
#include "group.h"
#include "random_mars.h"
#include "update.h"
#include "force.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static const double kappa_lb=0.0;

FixLbFluid::FixLbFluid(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  //=====================================================================================================
  //  Sample inputfile call:
  // fix # group lb/fluid nevery typeLB viscosity densityinit_real
  //  
  //  where: nevery:            call this fix every nevery timesteps. 
  //		                 (nevery generally set to 1).
  //         typeLB:            there are two different integrators 
  //                             in the code labelled "1" and "2".
  //         viscosity:         the viscosity of the fluid. 
  //         densityinit_real:  the density of the fluid.
  //
  // optional arguments:
  //  "setArea" type node_area:                       set the surface area per node associated with a
  //                                                   given atom type.  By default the surface area 
  //                                                   is set at 1.0*dx_lb^2.
  //  "setGamma" gamma:                               specify a user-defined value for the force
  //                                                   coupling constant, instead of using the default
  //                                                   value.
  //  "scaleGamma" type scale_factor:                 scale the user provided force coupling constant
  //                                                   by the factor, scale_factor, for the given atom
  //                                                   type.
  //  "dx" dx_lb:                                     the lattice-Boltzmann grid spacing.
  //  "dm" dm_lb:                                     the lattice-Boltzmann mass unit.
  //  "a0" a_0_real:                                  the square of the sound speed in the fluid.
  //  "noise" Temperature seed:                       include noise in the system.  
  //                                                   Temperature is the temperature for the fluid.
  //                                                   seed is the seed for the random number generator.  
  //  "calcforce" N group:                            print the force acting on a given group every 
  //                                                   N timesteps.
  //  "trilinear":                                    use the trilinear interpolation stencil.
  //  "read_restart" restart_file:                    restart a fluid run from restart_file.
  //  "write_restart" N:                              write a fluid restart file every N timesteps.
  //  "zwall_velocity" velocity_bottom velocity_top:  assign velocities to the z-walls
  //                                                   in the system.
  //  "bodyforce" bodyforcex bodyforcey bodyforcez:   add a constant body force to the
  //                                                   fluid.
  //  "printfluid" N:                                 print the fluid density and velocity at each
  //                                                   grid point every N timesteps.
  //  "D3Q19":                                        use the 19 velocity D3Q19 model.  By default,
  //                                                   the 15 velocity D3Q15 model is used.
  //=====================================================================================================

  if(narg <7) error->all(FLERR,"Illegal fix lb/fluid command");

  if (comm->style != 0) 
    error->universe_all(FLERR,"Fix lb/fluid can only currently be used with "
                        "comm_style brick");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  nevery = atoi(arg[3]);
  typeLB = atoi(arg[4]);
  viscosity = atof(arg[5]);
  densityinit_real = atof(arg[6]);
  
  // Default values for optional arguments:
  force_diagnostic=0;
  noisestress = 0;
  trilinear_stencil = 0;
  readrestart = 0;
  printrestart = 0;
  bodyforcex = bodyforcey = bodyforcez = 0.0;
  vwtp = vwbt = 0.0;
  printfluid = 0;
  T = 300.0;
  dm_lb = 1.0;
  fixviscouslb = 0;
  setdx = 1;
  seta0 = 1;
  setGamma = 0;
  setArea = 0;
  numvel = 15;

  Gamma = NULL;
  NodeArea = NULL;

  int iarg = 7;
  while (iarg < narg){
    if(strcmp(arg[iarg],"setArea")==0){
      if(setGamma == 1)
	error->all(FLERR,"Illegal fix lb/fluid command: cannot use a combination of default and user-specified gamma values");
      setArea = 1;
      int itype = atoi(arg[iarg+1]);
      double areafactor = atof(arg[iarg+2]);
      if(itype <= 0 || itype > atom->ntypes || areafactor < 0.0)
	error->all(FLERR,"Illegal fix lb/fluid command: setArea");
      if(NodeArea == NULL){
	NodeArea = new double[atom->ntypes+1];
	for(int i=0; i<=atom->ntypes; i++) NodeArea[i] = -1.0;
      }
      NodeArea[itype] = areafactor;
      iarg += 3;
    }
    else if(strcmp(arg[iarg],"setGamma")==0){
      if(setArea == 1)
	error->all(FLERR,"Illegal fix lb/fluid command: cannot use a combination of default and user-specified gamma values");
      setGamma = 1;
      double Gammaone;
      Gammaone = atof(arg[iarg+1]);
      if(Gamma == NULL)
	Gamma = new double[atom->ntypes+1];
      for(int i=0; i<=atom->ntypes; i++) Gamma[i] = Gammaone;
      iarg += 2;
    }
    else if(strcmp(arg[iarg],"scaleGamma")==0){
      if(setGamma == 0)
	error->all(FLERR,"Illegal fix lb/fluid command: must set a value for Gamma before scaling it");
      int itype = atoi(arg[iarg+1]);
      double scalefactor = atof(arg[iarg+2]);
      if(itype <= 0 || itype > atom->ntypes || scalefactor < 0.0)
	error->all(FLERR,"Illegal fix lb/fluid command: scaleGamma");
      Gamma[itype] *= scalefactor;
      iarg += 3;
    }     
    else if(strcmp(arg[iarg],"dx")==0){
      dx_lb = atof(arg[iarg+1]);
      iarg += 2;
      setdx = 0;
    }
    else if(strcmp(arg[iarg],"dm")==0){
      dm_lb = atof(arg[iarg+1]);
      iarg += 2;
    }
    else if(strcmp(arg[iarg],"a0")==0){
      a_0_real = atof(arg[iarg+1]);
      iarg += 2;
      seta0 = 0;
    }
    else if(strcmp(arg[iarg],"noise")== 0){
      noisestress = 1;
      T = atof(arg[iarg+1]);
      seed = atoi(arg[iarg+2]);
      iarg += 3;
    }
    else if(strcmp(arg[iarg],"calcforce")==0){
      force_diagnostic = atoi(arg[iarg+1]);
      if(force_diagnostic % nevery != 0){
	char str[200];
	sprintf(str,"Requesting calcforce output every %i timesteps. Will only print output for those timesteps that are a multiple of nevery.",force_diagnostic);
	error->warning(FLERR,str);	
      }
      igroupforce=group->find(arg[iarg+2]);
      iarg += 3;
    }
    else if(strcmp(arg[iarg],"trilinear")==0){
      trilinear_stencil = 1;
      iarg += 1;
    }
    else if(strcmp(arg[iarg],"read_restart")==0){
      readrestart = 1;  
      int nlength = strlen(arg[iarg+1]) + 16;
      char *filename = new char[nlength];
      strcpy(filename,arg[iarg+1]); 
      MPI_File_open(world,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,&pFileRead);
      delete [] filename;
      iarg += 2;
    }
    else if(strcmp(arg[iarg],"write_restart")==0){
      printrestart = atoi(arg[iarg+1]);
      if(printrestart % nevery != 0){
	char str[200];
	sprintf(str,"Requesting restart files every %i timesteps. Will only print restart files for those timesteps that are a multiple of nevery.",printrestart);
	error->warning(FLERR,str);	
      }
      iarg += 2;
    }
    else if(strcmp(arg[iarg],"zwall_velocity")==0){
      if(domain->periodicity[2]!=0) error->all(FLERR,"fix lb/fluid error: setting \
a z wall velocity without implementing fixed BCs in z");
      vwbt = atof(arg[iarg+1]);
      vwtp = atof(arg[iarg+2]);
      iarg += 3;
    }
    else if(strcmp(arg[iarg],"bodyforce")==0){
      bodyforcex = atof(arg[iarg+1]);
      bodyforcey = atof(arg[iarg+2]);
      bodyforcez = atof(arg[iarg+3]);
      iarg += 4;
    }
    else if(strcmp(arg[iarg],"printfluid")==0){
      printfluid = atoi(arg[iarg+1]);
      iarg += 2;
    }
    else if(strcmp(arg[iarg],"D3Q19")==0){
      numvel = 19;
      iarg += 1;
    }
    else error->all(FLERR,"Illegal fix lb/fluid command");
  }

  //--------------------------------------------------------------------------
  //Choose between D3Q15 and D3Q19 functions:
  //--------------------------------------------------------------------------
  if(numvel == 15){
    initializeLB = &FixLbFluid::initializeLB15;
    equilibriumdist = &FixLbFluid::equilibriumdist15;
    update_full = &FixLbFluid::update_full15;
  }else{
    initializeLB = &FixLbFluid::initializeLB19;
    equilibriumdist = &FixLbFluid::equilibriumdist19;
    update_full = &FixLbFluid::update_full19; 
  }  
  
  //--------------------------------------------------------------------------
  // perform initial allocation of atom-based array register
  // with Atom class
  //--------------------------------------------------------------------------
  hydroF = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  for(int i=0; i<atom->nmax; i++)
    for(int j=0; j<3; j++)
    hydroF[i][j] = 0.0;
 
  Ng_lb = NULL;
  w_lb = NULL;
  mg_lb = NULL;
  e = NULL;
  feq = NULL;
  feqold = NULL;
  feqn = NULL;
  feqoldn = NULL;
  f_lb = NULL;
  fnew = NULL;
  density_lb = NULL;
  u_lb = NULL;
  altogether = NULL;
  buf = NULL;
  Ff = NULL;
  Fftempx = NULL;
  Fftempy = NULL;
  Fftempz = NULL;

  //--------------------------------------------------------------------------
  // Set the lattice Boltzmann dt. 
  //--------------------------------------------------------------------------
  dt_lb=nevery*(update->dt);

  //--------------------------------------------------------------------------
  // Set the lattice Boltzmann dx if it wasn't specified in the 
  // input.
  //--------------------------------------------------------------------------
  if(setdx == 1){
    double dx_lb1 = sqrt(3.0*viscosity*dt_lb/densityinit_real);
    double mindomain = std::min(std::min(domain->xprd/comm->procgrid[0],domain->yprd/comm->procgrid[1]),domain->zprd/comm->procgrid[2]);
    dx_lb = mindomain/floor(mindomain/dx_lb1);

    if(comm->me==0){
      char str[128];
      sprintf(str,"Setting the lattice-Boltzmann dx to %10.6f",dx_lb);
      error->message(FLERR,str);
    }   
  }
  //--------------------------------------------------------------------------
  // If the area per node has not been set by the user, set to the 
  // default value of dx_lb*dx_lb.
  //--------------------------------------------------------------------------
  if(setGamma == 0){
    if(setArea == 0){ 
      if(comm->me==0){
	error->message(FLERR,"Assuming an area per node of dx*dx for all of the MD particles.  This should only be used if these all correspond to point particles; otherwise, change using the setArea keyword");
      }
      NodeArea = new double[atom->ntypes+1];
      for(int i=0; i<=atom->ntypes; i++) NodeArea[i] = -1.0;
    }
    for(int i=0; i<=atom->ntypes; i++)
      if(NodeArea[i] < 0.0) NodeArea[i] = dx_lb*dx_lb;
  }
  //--------------------------------------------------------------------------
  // Set a0 if it wasn't specified in the input
  //--------------------------------------------------------------------------
  if(seta0 == 1)
    a_0_real = 0.33333333*dx_lb*dx_lb/dt_lb/dt_lb;

  //--------------------------------------------------------------------------
  // Check to make sure that the total number of grid points in each direction
  // divides evenly among the processors in that direction.
  // Shrink-wrapped boundary conditions (which are not permitted by this fix)
  // might cause a problem, so check for this.  A full check of the boundary 
  // conditions is performed in the init routine, rather than here, as it is
  // possible to change the BCs between runs.
  //--------------------------------------------------------------------------
  double aa;
  double eps=1.0e-8;
  aa = (domain->xprd/comm->procgrid[0])/dx_lb;
  if(fabs(aa - floor(aa+0.5)) > eps){
    if(domain->boundary[0][0] != 0){
      error->all(FLERR,"the x-direction must be periodic");
    }
    char errormessage[200];
    sprintf(errormessage,"With dx= %f, and the simulation domain divided by %i processors in the x direction, the simulation domain in the x direction must be a multiple of %f",dx_lb,comm->procgrid[0],comm->procgrid[0]*dx_lb);
    error->all(FLERR,errormessage);
  }
  aa = (domain->yprd/comm->procgrid[1])/dx_lb;
  if(fabs(aa - floor(aa+0.5)) > eps){
    if(domain->boundary[1][0] != 0){
      error->all(FLERR,"the y-direction must be periodic");
    }
    char errormessage[200];
    sprintf(errormessage,"With dx= %f, and the simulation domain divided by %i processors in the y direction, the simulation domain in the y direction must be a multiple of %f",dx_lb,comm->procgrid[1],comm->procgrid[1]*dx_lb);
    error->all(FLERR,errormessage);
  }
  aa = (domain->zprd/comm->procgrid[2])/dx_lb;
  if(fabs(aa - floor(aa+0.5)) > eps){
    if(domain->boundary[2][0] == 2 || domain->boundary[2][0] == 3){
      error->all(FLERR,"the z-direction can not have shrink-wrap boundary conditions");
    }
    char errormessage[200];
    sprintf(errormessage,"With dx= %f, and the simulation domain divided by %i processors in the z direction, the simulation domain in the z direction must be a multiple of %f",dx_lb,comm->procgrid[2],comm->procgrid[2]*dx_lb);
    error->all(FLERR,errormessage);
  }
  
  //--------------------------------------------------------------------------
  // Set the total number of grid points in each direction.
  //--------------------------------------------------------------------------
  Nbx = (int)(domain->xprd/dx_lb + 0.5);
  Nby = (int)(domain->yprd/dx_lb + 0.5);
  Nbz = (int)(domain->zprd/dx_lb + 0.5);

  //--------------------------------------------------------------------------
  // Set the number of grid points in each dimension for the local subgrids.
  //--------------------------------------------------------------------------
  subNbx= Nbx/comm->procgrid[0] + 2;
  subNby= Nby/comm->procgrid[1] + 2;
  subNbz= Nbz/comm->procgrid[2] + 2;
 
  //--------------------------------------------------------------------------
  // In order to calculate the fluid forces correctly, need to have atleast
  // 5 grid points in each direction per processor.
  //--------------------------------------------------------------------------
  if(subNbx<7 || subNby < 7 || subNbz<7)
    error->all(FLERR,"Need at least 5 grid points in each direction per processor");

  // If there are walls in the z-direction add an extra grid point.
  if(domain->periodicity[2]==0){
   Nbz += 1;
   if(comm->myloc[2]==comm->procgrid[2]-1)
     subNbz += 1;
  }

  if(comm->me==0){
    char str[128];
    if(setdx == 1){
      sprintf(str,"Using a lattice-Boltzmann grid of %i by %i by %i total grid points.  To change, use the dx keyword",Nbx,Nby,Nbz);
    }else{
      sprintf(str,"Using a lattice-Boltzmann grid of %i by %i by %i total grid points.",Nbx,Nby,Nbz);
    }
    error->message(FLERR,str);   
  }

  //--------------------------------------------------------------------------
  // Store the largest value of subNbz, which is needed for allocating the
  // buf array (since a processor with comm->myloc[2] == comm->procgrid[2]-1
  // may have an additional subNbz point as compared with the rest).
  //--------------------------------------------------------------------------
  int subNbzmax;
  MPI_Allreduce(&subNbz,&subNbzmax,1,MPI_INT,MPI_MAX,world);

  //--------------------------------------------------------------------------
  // Create the MPI datatypes used to pass portions of arrays:
  // datatypes to pass the f and feq arrays.
  //--------------------------------------------------------------------------
  MPI_Aint lb, sizeofdouble;
  MPI_Type_get_extent(MPI_DOUBLE,&lb,&sizeofdouble);
  
  MPI_Type_vector(subNbz-2,numvel,numvel,MPI_DOUBLE,&oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNby-2,1,numvel*subNbz*sizeofdouble,oneslice,&passxf);
  MPI_Type_commit(&passxf);
 
  MPI_Type_create_hvector(subNbx,1,numvel*subNbz*subNby*sizeofdouble,oneslice,&passyf);
  MPI_Type_commit(&passyf);
  
  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNby,numvel,numvel*subNbz,MPI_DOUBLE,&oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNbx,1,numvel*subNbz*subNby*sizeofdouble,oneslice,&passzf);
  MPI_Type_commit(&passzf);

  // datatypes to pass the u array, and the Ff array.
  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNbz+3,3,3,MPI_DOUBLE,&oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNby+3,1,3*(subNbz+3)*sizeofdouble,oneslice,&passxu);
  MPI_Type_commit(&passxu);
  
  MPI_Type_create_hvector(subNbx+3,1,3*(subNbz+3)*(subNby+3)*sizeofdouble,oneslice,&passyu);
  MPI_Type_commit(&passyu);
  
  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNby+3,3,3*(subNbz+3),MPI_DOUBLE,&oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNbx+3,1,3*(subNbz+3)*(subNby+3)*sizeofdouble,oneslice,&passzu);
  MPI_Type_commit(&passzu);

  // datatypes to pass the density array.
  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNbz+3,1,1,MPI_DOUBLE,&oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNby+3,1,1*(subNbz+3)*sizeofdouble,oneslice,&passxrho);
  MPI_Type_commit(&passxrho);
  
  MPI_Type_create_hvector(subNbx+3,1,1*(subNbz+3)*(subNby+3)*sizeofdouble,oneslice,&passyrho);
  MPI_Type_commit(&passyrho);
  
  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNby+3,1,1*(subNbz+3),MPI_DOUBLE,&oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNbx+3,1,1*(subNbz+3)*(subNby+3)*sizeofdouble,oneslice,&passzrho);
  MPI_Type_commit(&passzrho);

  // datatypes to receive a portion of the Ff array.
  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNbz+3,3,3,MPI_DOUBLE,&oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNby+3,1,3*(subNbz+3)*sizeofdouble,oneslice,&passxtemp);
  MPI_Type_commit(&passxtemp);
  
  MPI_Type_create_hvector(subNbx+3,1,3*(subNbz+3)*5*sizeofdouble,oneslice,&passytemp);
  MPI_Type_commit(&passytemp);
  
  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNby+3,3,3*5,MPI_DOUBLE,&oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNbx+3,1,3*5*(subNby+3)*sizeofdouble,oneslice,&passztemp);
  MPI_Type_commit(&passztemp);

  MPI_Type_free(&oneslice);

  //--------------------------------------------------------------------------
  // Allocate the necessary arrays.
  //--------------------------------------------------------------------------
  memory->create(Ng_lb,numvel,"FixLbFluid:Ng_lb");
  memory->create(w_lb,numvel,"FixLbFluid:w_lb");
  memory->create(mg_lb,numvel,numvel,"FixLbFluid:mg_lb");
  memory->create(e,numvel,3,"FixLbFluid:e");
  memory->create(feq,subNbx,subNby,subNbz,numvel,"FixLbFluid:feq");
  if(typeLB == 2){
    memory->create(feqold,subNbx,subNby,subNbz,numvel,"FixLbFluid:feqold");
    memory->create(feqn,subNbx,subNby,subNbz,numvel,"FixLbFluid:feqn");
    memory->create(feqoldn,subNbx,subNby,subNbz,numvel,"FixLbFluid:feqoldn");
  }
  memory->create(f_lb,subNbx,subNby,subNbz,numvel,"FixLbFluid:f_lb");
  memory->create(fnew,subNbx,subNby,subNbz,numvel,"FixLbFluid:fnew");
  memory->create(density_lb,subNbx+3,subNby+3,subNbz+3,"FixLbFluid:density_lb");
  memory->create(u_lb,subNbx+3,subNby+3,subNbz+3,3,"FixLbFluid:u_lb");
  if(printfluid > 0){
    memory->create(buf,subNbx,subNby,subNbzmax,4,"FixLbFluid:buf");
    if(me==0)
      memory->create(altogether,Nbx,Nby,Nbz,4,"FixLbFluid:altogether");
  }
  memory->create(Ff,subNbx+3,subNby+3,subNbz+3,3,"FixLbFluid:Ff");
  memory->create(Fftempx,5,subNby+3,subNbz+3,3,"FixLbFluid:Fftempx");
  memory->create(Fftempy,subNbx+3,5,subNbz+3,3,"FixLbFluid:Fftempy");
  memory->create(Fftempz,subNbx+3,subNby+3,5,3,"FixLbFluid:Fftempz");

  if(noisestress==1){
    random = new RanMars(lmp,seed + comm->me);
  }

  //--------------------------------------------------------------------------
  // Rescale the variables to Lattice Boltzmann dimensionless units.
  //--------------------------------------------------------------------------
  rescale();
  
  //--------------------------------------------------------------------------
  // Initialize the arrays.
  //--------------------------------------------------------------------------
  (*this.*initializeLB)();
  initialize_feq();

}

FixLbFluid::~FixLbFluid()
{

  atom->delete_callback(id,0);
  memory->destroy(hydroF);
 
  memory->destroy(Ng_lb);
  memory->destroy(w_lb);
  memory->destroy(mg_lb);
  memory->destroy(e);
  memory->destroy(feq);
  if(typeLB == 2){
    memory->destroy(feqold);
    memory->destroy(feqn);
    memory->destroy(feqoldn);
  }
  memory->destroy(f_lb);
  memory->destroy(fnew);
  memory->destroy(density_lb);
  memory->destroy(u_lb);
  if(printfluid>0){
    if(me==0)
      memory->destroy(altogether);
    memory->destroy(buf);
  }
  memory->destroy(Ff);
  memory->destroy(Fftempx);
  memory->destroy(Fftempy);
  memory->destroy(Fftempz);
  
  if(noisestress==1){
    delete random;
  }

  if(setGamma == 1){
    delete [] Gamma;
  }else{
    delete [] NodeArea;
  }
}

int FixLbFluid::setmask()
{
  int mask =0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

void FixLbFluid::init(void)
{
  int i,j;

  if (comm->style != 0) 
    error->universe_all(FLERR,"Fix lb/fluid can only currently be used with "
                        "comm_style brick");

  //--------------------------------------------------------------------------
  // Check to see if the MD timestep has changed between runs.
  //--------------------------------------------------------------------------
  double dt_lb_now;
  dt_lb_now=nevery*(update->dt);

  if(fabs(dt_lb_now - dt_lb) > 1.0e-12){
    error->warning(FLERR,"Timestep has changed between runs with the same lb/fluid.  Unphysical results may occur");
  }   
  
  //--------------------------------------------------------------------------
  // Make sure the size of the simulation domain has not changed
  // between runs.
  //--------------------------------------------------------------------------
  int Nbx_now,Nby_now,Nbz_now;
  Nbx_now = (int)(domain->xprd/dx_lb + 0.5);
  Nby_now = (int)(domain->yprd/dx_lb + 0.5);
  Nbz_now = (int)(domain->zprd/dx_lb + 0.5);
  // If there are walls in the z-direction add an extra grid point.
  if(domain->periodicity[2]==0){
   Nbz_now += 1;
  }
  
  if(Nbx_now != Nbx || Nby_now != Nby || Nbz_now != Nbz){
    error->all(FLERR,"the simulation domain can not change shape between runs with the same lb/fluid");
  }
  
  //--------------------------------------------------------------------------
  // Check to make sure that the chosen LAMMPS boundary types are compatible
  // with this fix.
  //    shrink-wrap is not compatible in any dimension.
  //    fixed only works in the z-direction.
  //--------------------------------------------------------------------------
  if(domain->boundary[0][0] != 0){
    error->all(FLERR,"the x-direction must be periodic");
  }
  if(domain->boundary[1][0] != 0){
    error->all(FLERR,"the y-direction must be periodic");
  }
  if(domain->boundary[2][0] == 2 || domain->boundary[2][0] == 3){
    error->all(FLERR,"the z-direction can not have shrink-wrap boundary conditions");
  }
  
  //--------------------------------------------------------------------------
  // Check if the lb/viscous fix is also called:
  //--------------------------------------------------------------------------
  groupbit_viscouslb = groupbit_pc = groupbit_rigid_pc_sphere = 0;
  for (i = 0; i < modify->nfix; i++){
    if (strcmp(modify->fix[i]->style,"lb/viscous") == 0){
      fixviscouslb = 1;
      groupbit_viscouslb = group->bitmask[modify->fix[i]->igroup];
    }
    if(strcmp(modify->fix[i]->style,"lb/pc")==0){
      groupbit_pc = group->bitmask[modify->fix[i]->igroup];
    }
    if(strcmp(modify->fix[i]->style,"lb/rigid/pc/sphere")==0){
      groupbit_rigid_pc_sphere = group->bitmask[modify->fix[i]->igroup];
    }
  }

  // Warn if the fluid force is not applied to any of the particles.
  if(!(groupbit_viscouslb || groupbit_pc || groupbit_rigid_pc_sphere) && comm->me==0){
    error->message(FLERR,"Not adding the fluid force to any of the MD particles.  To add this force use one of the lb/viscous, lb/pc, or lb/rigid/pc/sphere fixes");
  }
  
  // If fix lb/viscous is called for a particular atom, make sure 
  // lb/pc or lb/rigid/pc/sphere are not:
  if(fixviscouslb == 1){
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
      for(j=0; j<nlocal; j++){
	if((mask[j] & groupbit) && (mask[j] & groupbit_viscouslb) && (mask[j] & groupbit_pc))
	  error->one(FLERR,"should not use the lb/viscous command when integrating with the lb/pc fix");
	if((mask[j] & groupbit) && (mask[j] & groupbit_viscouslb) && (mask[j] & groupbit_rigid_pc_sphere))
	  error->one(FLERR,"should not use the lb/viscous command when integrating with the lb/rigid/pc/sphere fix");
      }
   }
 
}

void FixLbFluid::setup(int vflag)
{
  //--------------------------------------------------------------------------
  // Need to calculate the force on the fluid for a restart run.
  //--------------------------------------------------------------------------
  if(step > 0)
    calc_fluidforce();
}  

void FixLbFluid::initial_integrate(int vflag)
{
  // only call every nevery timesteps (by default nevery only affects how
  // often end_of_step is called.
  if(update->ntimestep % nevery == 0){
    //--------------------------------------------------------------------------
    // Print a header labelling any output printed to the screen.
    //--------------------------------------------------------------------------
    static int printheader = 1;
    
    if(printheader == 1){
      if(force_diagnostic > 0 && me == 0){
	printf("-------------------------------------------------------------------------------\n");
	printf("     F_x          F_y          F_z          T_x          T_y          T_z\n");
	printf("-------------------------------------------------------------------------------\n");
      }
      
      if(printfluid > 0 && me == 0){
	printf("---------------------------------------------------------------------\n");
	printf("     density            u_x              u_y              u_z \n");
	printf("---------------------------------------------------------------------\n");
      }
      printheader = 0;
    }
    
    //--------------------------------------------------------------------------
    // Determine the equilibrium distribution on the local subgrid.
    //--------------------------------------------------------------------------
    (*this.*equilibriumdist)(1,subNbx-1,1,subNby-1,1,subNbz-1);
    
    //--------------------------------------------------------------------------
    // Using the equilibrium distribution, calculate the new
    // distribution function.
    //--------------------------------------------------------------------------
    (*this.*update_full)();
    
    std::swap(f_lb,fnew); 
    
    //--------------------------------------------------------------------------
    // Calculate moments of the distribution function.
    //--------------------------------------------------------------------------
    parametercalc_full();
    
    //--------------------------------------------------------------------------
    // Store the equilibrium distribution function, it is needed in
    // the next time step by the update routine.
    //--------------------------------------------------------------------------
    if(typeLB == 2){
      std::swap(feqold,feq);
      std::swap(feqoldn,feqn);
    }
  
  }  
 
  //--------------------------------------------------------------------------
  // Perform diagnostics, and print output for the graphics program
  //--------------------------------------------------------------------------
  if(printfluid > 0 && update->ntimestep > 0 && (update->ntimestep % printfluid == 0))
    streamout();
   
}
void FixLbFluid::post_force(int vflag)
{
  // only call every nevery timesteps (by default nevery only affects how
  // often end_of_step is called.
  if(update->ntimestep % nevery == 0){
    if(fixviscouslb==1)
      calc_fluidforce();
  }
}

void FixLbFluid::end_of_step()
{  
  // end_of_step is only called every nevery timesteps
  if(fixviscouslb==0)
    calc_fluidforce();
 
  if(printrestart>0){
    if((update->ntimestep)%printrestart == 0){
      write_restartfile();
    }
  }

}

//==========================================================================
//   allocate atom-based array
//==========================================================================
void FixLbFluid::grow_arrays(int nmax)
{
  memory->grow(hydroF,nmax,3,"FixLbFluid:hydroF");
}

//==========================================================================
//   copy values within local atom-based array
//==========================================================================
void FixLbFluid::copy_arrays(int i, int j, int delflag)
{
  hydroF[j][0] = hydroF[i][0];
  hydroF[j][1] = hydroF[i][1];
  hydroF[j][2] = hydroF[i][2];
}

//==========================================================================
//   pack values in local atom-based array for exchange with another proc
//==========================================================================
int FixLbFluid::pack_exchange(int i, double *buf)
{
  buf[0] = hydroF[i][0];
  buf[1] = hydroF[i][1];
  buf[2] = hydroF[i][2];

  return 3;
}

//==========================================================================
//   unpack values in local atom-based array from exchange with another proc
//==========================================================================
int FixLbFluid::unpack_exchange(int nlocal, double *buf)
{
  hydroF[nlocal][0] = buf[0];
  hydroF[nlocal][1] = buf[1];
  hydroF[nlocal][2] = buf[2];

  return 3;
}

//==========================================================================
//   calculate the force from the local atoms acting on the fluid.
//==========================================================================
void FixLbFluid::calc_fluidforce(void)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  int i,j,k,m;
  MPI_Request requests[20];
  double forceloc[3],force[3];
  double torqueloc[3],torque[3];
  
  //--------------------------------------------------------------------------
  // Zero out arrays
  //--------------------------------------------------------------------------
  std::fill(&Ff[0][0][0][0],&Ff[0][0][0][0] + (subNbx+3)*(subNby+3)*(subNbz+3)*3,0.0);
  std::fill(&Fftempx[0][0][0][0],&Fftempx[0][0][0][0] + 5*(subNby+3)*(subNbz+3)*3,0.0);
  std::fill(&Fftempy[0][0][0][0],&Fftempy[0][0][0][0] + (subNbx+3)*5*(subNbz+3)*3,0.0);
  std::fill(&Fftempz[0][0][0][0],&Fftempz[0][0][0][0] + (subNbx+3)*(subNby+3)*5*3,0.0);

  forceloc[0] = forceloc[1] = forceloc[2] = 0.0;
  torqueloc[0] = torqueloc[1] = torqueloc[2] = 0.0;

  for(i=0; i<atom->nmax; i++)
    for(j=0; j<3; j++)
      hydroF[i][j] = 0.0;
  
  
  double unwrap[3];
  double dx,dy,dz;
  double massone;
  imageint *image = atom->image;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  double sum[4],xcm[4];

  if(force_diagnostic > 0 && update->ntimestep > 0 && (update->ntimestep % force_diagnostic == 0)){
    //Calculate the center of mass of the particle group
    //(needed to calculate the torque).
    sum[0] = sum[1] = sum[2] = sum[3] = 0.0;
    for(i=0; i<nlocal; i++){
      if(mask[i] & group->bitmask[igroupforce]){

	domain->unmap(x[i],image[i],unwrap);

	if(rmass) massone = rmass[i];
	else massone = mass[type[i]];
	
	sum[0] += unwrap[0]*massone;
	sum[1] += unwrap[1]*massone;
	sum[2] += unwrap[2]*massone;
	sum[3] += massone;
      }
    }
    MPI_Allreduce(&sum[0],&xcm[0],4,MPI_DOUBLE,MPI_SUM,world);
    xcm[0] = xcm[0]/xcm[3];
    xcm[1] = xcm[1]/xcm[3];
    xcm[2] = xcm[2]/xcm[3];
  }

  //--------------------------------------------------------------------------
  //Calculate the contribution to the force on the fluid.
  //--------------------------------------------------------------------------
  for(i=0; i<nlocal; i++){
    if(mask[i] & groupbit){
      if(trilinear_stencil==1) {
	trilinear_interpolation(i);
      }else{
	peskin_interpolation(i);
      }
      
      if(force_diagnostic > 0 && update->ntimestep > 0 && (update->ntimestep % force_diagnostic == 0)){
	if(mask[i] & group->bitmask[igroupforce]){
	  
	  domain->unmap(x[i],image[i],unwrap);
	  dx = unwrap[0] - xcm[0];
	  dy = unwrap[1] - xcm[1];
	  dz = unwrap[2] - xcm[2];
	  
	  forceloc[0] += hydroF[i][0];
	  forceloc[1] += hydroF[i][1];
	  forceloc[2] += hydroF[i][2];
	  torqueloc[0] += dy*hydroF[i][2] - dz*hydroF[i][1];
	  torqueloc[1] += dz*hydroF[i][0] - dx*hydroF[i][2];
	  torqueloc[2] += dx*hydroF[i][1] - dy*hydroF[i][0];
	}	    
      } 
    }
  }

  //--------------------------------------------------------------------------
  //Communicate the force contributions which lie outside the local processor
  //sub domain.
  //--------------------------------------------------------------------------
  for(i=0; i<10; i++)
    requests[i]=MPI_REQUEST_NULL;
  MPI_Isend(&Ff[0][0][0][0],1,passxu,comm->procneigh[0][0],10,world,&requests[0]);
  MPI_Isend(&Ff[subNbx+2][0][0][0],1,passxu,comm->procneigh[0][0],20,world,&requests[1]);
  MPI_Isend(&Ff[subNbx-1][0][0][0],1,passxu,comm->procneigh[0][1],30,world,&requests[2]);
  MPI_Isend(&Ff[subNbx][0][0][0],1,passxu,comm->procneigh[0][1],40,world,&requests[3]);
  MPI_Isend(&Ff[subNbx+1][0][0][0],1,passxu,comm->procneigh[0][1],50,world,&requests[4]);
  MPI_Irecv(&Fftempx[0][0][0][0],1,passxtemp,comm->procneigh[0][1],10,world,&requests[5]);
  MPI_Irecv(&Fftempx[1][0][0][0],1,passxtemp,comm->procneigh[0][1],20,world,&requests[6]);
  MPI_Irecv(&Fftempx[2][0][0][0],1,passxtemp,comm->procneigh[0][0],30,world,&requests[7]);
  MPI_Irecv(&Fftempx[3][0][0][0],1,passxtemp,comm->procneigh[0][0],40,world,&requests[8]);
  MPI_Irecv(&Fftempx[4][0][0][0],1,passxtemp,comm->procneigh[0][0],50,world,&requests[9]);
  MPI_Waitall(10,requests,MPI_STATUS_IGNORE);
  
  for(j=0; j<subNby+3; j++){
    for(k=0; k<subNbz+3; k++){
      for(m=0; m<3; m++){
	Ff[subNbx-2][j][k][m] += Fftempx[0][j][k][m];
	Ff[subNbx-3][j][k][m] += Fftempx[1][j][k][m];
	Ff[1][j][k][m] += Fftempx[2][j][k][m];
	Ff[2][j][k][m] += Fftempx[3][j][k][m];
	Ff[3][j][k][m] += Fftempx[4][j][k][m];
      }
    }
  }

  for(i=0; i<10; i++)
    requests[i]=MPI_REQUEST_NULL;
  MPI_Isend(&Ff[0][0][0][0],1,passyu,comm->procneigh[1][0],10,world,&requests[0]);
  MPI_Isend(&Ff[0][subNby+2][0][0],1,passyu,comm->procneigh[1][0],20,world,&requests[1]);
  MPI_Isend(&Ff[0][subNby-1][0][0],1,passyu,comm->procneigh[1][1],30,world,&requests[2]);
  MPI_Isend(&Ff[0][subNby][0][0],1,passyu,comm->procneigh[1][1],40,world,&requests[3]);
  MPI_Isend(&Ff[0][subNby+1][0][0],1,passyu,comm->procneigh[1][1],50,world,&requests[4]);
  MPI_Irecv(&Fftempy[0][0][0][0],1,passytemp,comm->procneigh[1][1],10,world,&requests[5]);
  MPI_Irecv(&Fftempy[0][1][0][0],1,passytemp,comm->procneigh[1][1],20,world,&requests[6]);
  MPI_Irecv(&Fftempy[0][2][0][0],1,passytemp,comm->procneigh[1][0],30,world,&requests[7]);
  MPI_Irecv(&Fftempy[0][3][0][0],1,passytemp,comm->procneigh[1][0],40,world,&requests[8]);
  MPI_Irecv(&Fftempy[0][4][0][0],1,passytemp,comm->procneigh[1][0],50,world,&requests[9]);
  MPI_Waitall(10,requests,MPI_STATUS_IGNORE);

  for(i=0; i<subNbx+3; i++){
    for(k=0; k<subNbz+3; k++){
      for(m=0; m<3; m++){
	Ff[i][subNby-2][k][m] += Fftempy[i][0][k][m];
	Ff[i][subNby-3][k][m] += Fftempy[i][1][k][m];
	Ff[i][1][k][m] += Fftempy[i][2][k][m];
	Ff[i][2][k][m] += Fftempy[i][3][k][m];
	Ff[i][3][k][m] += Fftempy[i][4][k][m];
      }
    }
  }

  for(i=0; i<10; i++)
    requests[i]=MPI_REQUEST_NULL;
  MPI_Isend(&Ff[0][0][0][0],1,passzu,comm->procneigh[2][0],10,world,&requests[0]);
  MPI_Isend(&Ff[0][0][subNbz+2][0],1,passzu,comm->procneigh[2][0],20,world,&requests[1]);
  MPI_Isend(&Ff[0][0][subNbz-1][0],1,passzu,comm->procneigh[2][1],30,world,&requests[2]);
  MPI_Isend(&Ff[0][0][subNbz][0],1,passzu,comm->procneigh[2][1],40,world,&requests[3]);
  MPI_Isend(&Ff[0][0][subNbz+1][0],1,passzu,comm->procneigh[2][1],50,world,&requests[4]);
  MPI_Irecv(&Fftempz[0][0][0][0],1,passztemp,comm->procneigh[2][1],10,world,&requests[5]);
  MPI_Irecv(&Fftempz[0][0][1][0],1,passztemp,comm->procneigh[2][1],20,world,&requests[6]);
  MPI_Irecv(&Fftempz[0][0][2][0],1,passztemp,comm->procneigh[2][0],30,world,&requests[7]);
  MPI_Irecv(&Fftempz[0][0][3][0],1,passztemp,comm->procneigh[2][0],40,world,&requests[8]);
  MPI_Irecv(&Fftempz[0][0][4][0],1,passztemp,comm->procneigh[2][0],50,world,&requests[9]);
  MPI_Waitall(10,requests,MPI_STATUS_IGNORE);  

  for(i=0; i<subNbx+3; i++){
    for(j=0; j<subNby+3; j++){
      for(m=0; m<3; m++){
	Ff[i][j][subNbz-2][m] += Fftempz[i][j][0][m];
	Ff[i][j][subNbz-3][m] += Fftempz[i][j][1][m];
	Ff[i][j][1][m] += Fftempz[i][j][2][m];
	Ff[i][j][2][m] += Fftempz[i][j][3][m];
	Ff[i][j][3][m] += Fftempz[i][j][4][m];
      }
    }
  }

  if(force_diagnostic > 0 && update->ntimestep > 0 && (update->ntimestep % force_diagnostic == 0)){
    force[0] = force[1] = force[2] = 0.0;
    torque[0] = torque[1] = torque[2] =0.0;
    
    MPI_Allreduce(&forceloc[0],&force[0],3,MPI_DOUBLE,MPI_SUM,world); 
    MPI_Allreduce(&torqueloc[0],&torque[0],3,MPI_DOUBLE,MPI_SUM,world);
    
    if(me==0){
      printf("%E %E %E %E %E %E\n",force[0],force[1],force[2],
 	     torque[0],torque[1],torque[2]);

    }
  }
  
}
//==========================================================================
// uses the Peskin stencil to perform the velocity, density and
// force interpolations.
//==========================================================================
void FixLbFluid::peskin_interpolation(int i)
{
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double massone;
  int ix,iy,iz;
  int ixp,iyp,izp;
  double dx1,dy1,dz1;
  int isten,ii,jj,kk;
  double r,rsq,weightx,weighty,weightz;
  double FfP[64];
  int k;
  double unode[3];
  double mnode;
  double gammavalue;

  //--------------------------------------------------------------------------
  //Calculate nearest leftmost grid point.
  //Since array indices from 1 to subNb-2 correspond to the
  // local subprocessor domain (not indices from 0), use the 
  // ceiling value.
  //--------------------------------------------------------------------------
  ix = (int)ceil((x[i][0]-domain->sublo[0])/dx_lb);
  iy = (int)ceil((x[i][1]-domain->sublo[1])/dx_lb);
  iz = (int)ceil((x[i][2]-domain->sublo[2])/dx_lb);
	
  //--------------------------------------------------------------------------
  //Calculate distances to the nearest points.
  //--------------------------------------------------------------------------
  dx1 = x[i][0] - (domain->sublo[0] + (ix-1)*dx_lb);
  dy1 = x[i][1] - (domain->sublo[1] + (iy-1)*dx_lb);
  dz1 = x[i][2] - (domain->sublo[2] + (iz-1)*dx_lb);
  
  // Need to convert these to lattice units:
  dx1 = dx1/dx_lb;
  dy1 = dy1/dx_lb;
  dz1 = dz1/dx_lb;
  
  unode[0]=0.0; unode[1]=0.0; unode[2]=0.0;
  mnode = 0.0;
  isten=0;
    
  //--------------------------------------------------------------------------
  // Calculate the interpolation weights, and interpolated values of
  // the fluid velocity, and density.
  //--------------------------------------------------------------------------
  for(ii=-1; ii<3; ii++){
    rsq=(-dx1+ii)*(-dx1+ii);
    
    if(rsq>=4)
      weightx=0.0;
    else{
      r=sqrt(rsq);
      if(rsq>1){
	weightx=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
      } else{
	weightx=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
      }
    }
    for(jj=-1; jj<3; jj++){
      rsq=(-dy1+jj)*(-dy1+jj);
      if(rsq>=4)
	weighty=0.0;
      else{
	r=sqrt(rsq);
	if(rsq>1){
	  weighty=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
	} else{
	  weighty=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
	}
      }
      for(kk=-1; kk<3; kk++){
	rsq=(-dz1+kk)*(-dz1+kk);
	if(rsq>=4)
	  weightz=0.0;
	else{
	  r=sqrt(rsq);
	  if(rsq>1){
	    weightz=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
	  } else{
	    weightz=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
	  }
	}
	ixp = ix+ii;
	iyp = iy+jj;
	izp = iz+kk;

	//The atom is allowed to be within one lattice grid point outside the
	//local processor sub-domain.  
	if(ixp < -1 || ixp > (subNbx+1) || iyp < -1 || iyp > (subNby+1) || izp < -1 || izp > (subNbz+1))
	  error->one(FLERR,"Atom outside local processor simulation domain.  Either unstable fluid pararmeters, or \
require more frequent neighborlist rebuilds");

	if(domain->periodicity[2] == 0 && comm->myloc[2] == 0 && izp < 1)
	  error->warning(FLERR,"Atom too close to lower z wall.  Unphysical results may occur");
	if(domain->periodicity[2] == 0 && comm->myloc[2] == (comm->procgrid[2]-1) && (izp > (subNbz-2) ))
	  error->warning(FLERR,"Atom too close to upper z wall.  Unphysical results may occur");
	
	if(ixp==-1) ixp=subNbx+2;
	if(iyp==-1) iyp=subNby+2;
	if(izp==-1) izp=subNbz+2;
	
	FfP[isten] = weightx*weighty*weightz;
	// interpolated velocity based on delta function.
	for(k=0; k<3; k++){
	  unode[k] += u_lb[ixp][iyp][izp][k]*FfP[isten];
	}
	if(setGamma==0)
	  mnode += density_lb[ixp][iyp][izp]*FfP[isten];
	
	isten++;
      }
    }
  }
  if(setGamma==0){
    mnode *= NodeArea[type[i]];

    if(rmass) massone = rmass[i];
    else massone = mass[type[i]];
    massone = massone/dm_lb; 

    gammavalue = 2.0*(mnode*massone)*dtoverdtcollision/(mnode+massone);
  }    
  else{
    gammavalue = Gamma[type[i]];
  }
  
  isten=0;
  for(ii=-1; ii<3; ii++)
    for(jj=-1; jj<3; jj++)
      for(kk=-1; kk<3; kk++){
	ixp = ix+ii;
	iyp = iy+jj;
	izp = iz+kk;
	
	if(ixp==-1) ixp=subNbx+2;
	if(iyp==-1) iyp=subNby+2;
	if(izp==-1) izp=subNbz+2;
	// Compute the force on the fluid.  Need to convert the velocity from
	// LAMMPS units to LB units.
	for(k=0; k<3; k++){
	  Ff[ixp][iyp][izp][k] += gammavalue*((v[i][k]*dt_lb/dx_lb)-unode[k])*FfP[isten];
	}
	
	isten++;
      }
  for(k=0; k<3; k++)
    hydroF[i][k] = -1.0*gammavalue*((v[i][k]*dt_lb/dx_lb)-unode[k])*dm_lb*dx_lb/dt_lb/dt_lb;
}

//==========================================================================
// uses the trilinear stencil to perform the velocity, density and
// force interpolations.
//==========================================================================
void FixLbFluid::trilinear_interpolation(int i)
{
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double massone;
  int ix,iy,iz;
  int ixp,iyp,izp;
  double dx1,dy1,dz1;
  double FfP[8];
  int k;
  double unode[3];
  double mnode;
  double gammavalue;

  //--------------------------------------------------------------------------
  // Calculate nearest leftmost grid point.
  // Since array indices from 1 to subNb-2 correspond to the
  // local subprocessor domain (not indices from 0), use the 
  // ceiling value.
  //--------------------------------------------------------------------------
  ix = (int)ceil((x[i][0]-domain->sublo[0])/dx_lb);
  iy = (int)ceil((x[i][1]-domain->sublo[1])/dx_lb);
  iz = (int)ceil((x[i][2]-domain->sublo[2])/dx_lb);
	
  //--------------------------------------------------------------------------
  //Calculate distances to the nearest points.
  //--------------------------------------------------------------------------
  dx1 = x[i][0] - (domain->sublo[0] + (ix-1)*dx_lb);
  dy1 = x[i][1] - (domain->sublo[1] + (iy-1)*dx_lb);
  dz1 = x[i][2] - (domain->sublo[2] + (iz-1)*dx_lb);
 
  //--------------------------------------------------------------------------
  // Need to convert these to lattice units:
  //--------------------------------------------------------------------------
  dx1 = dx1/dx_lb;
  dy1 = dy1/dx_lb;
  dz1 = dz1/dx_lb;

  //--------------------------------------------------------------------------
  // Calculate the interpolation weights
  //--------------------------------------------------------------------------
  FfP[0] = (1.-dx1)*(1.-dy1)*(1.-dz1);
  FfP[1] = (1.-dx1)*(1.-dy1)*dz1;
  FfP[2] = (1.-dx1)*dy1*(1.-dz1);
  FfP[3] = (1.-dx1)*dy1*dz1;
  FfP[4] = dx1*(1.-dy1)*(1.-dz1);
  FfP[5] = dx1*(1.-dy1)*dz1;
  FfP[6] = dx1*dy1*(1.-dz1);
  FfP[7] = dx1*dy1*dz1;
  
  ixp = (ix+1);
  iyp = (iy+1);
  izp = (iz+1);

  //The atom is allowed to be within one lattice grid point outside the
  //local processor sub-domain.  
  if(ix < 0 || ixp > (subNbx+1) || iy < 0 || iyp > (subNby+1) || iz < 0 || izp > (subNbz+1))
    error->one(FLERR,"Atom outside local processor simulation domain.  Either unstable fluid pararmeters, or \
require more frequent neighborlist rebuilds");

  if(domain->periodicity[2] == 0 && comm->myloc[2] == 0 && (iz < 1 || izp < 1))
    error->warning(FLERR,"Atom too close to lower z wall.  Unphysical results may occur");
  if(domain->periodicity[2] == 0 && comm->myloc[2] == (comm->procgrid[2]-1) && (izp > (subNbz-2) || iz > (subNbz-2)))
    error->warning(FLERR,"Atom too close to upper z wall.  Unphysical results may occur");  
  
   
  for (k=0; k<3; k++) { 	// tri-linearly interpolated velocity at node
    unode[k] = u_lb[ix][iy][iz][k]*FfP[0]
      + u_lb[ix][iy][izp][k]*FfP[1]
      + u_lb[ix][iyp][iz][k]*FfP[2]
      + u_lb[ix][iyp][izp][k]*FfP[3]
      + u_lb[ixp][iy][iz][k]*FfP[4]
      + u_lb[ixp][iy][izp][k]*FfP[5]
      + u_lb[ixp][iyp][iz][k]*FfP[6]
      + u_lb[ixp][iyp][izp][k]*FfP[7];
  }

  if(setGamma==0){
    mnode = density_lb[ix][iy][iz]*FfP[0]
      + density_lb[ix][iy][izp]*FfP[1]
      + density_lb[ix][iyp][iz]*FfP[2]
      + density_lb[ix][iyp][izp]*FfP[3]
      + density_lb[ixp][iy][iz]*FfP[4]
      + density_lb[ixp][iy][izp]*FfP[5]
      + density_lb[ixp][iyp][iz]*FfP[6]
      + density_lb[ixp][iyp][izp]*FfP[7];

    mnode *= NodeArea[type[i]];

    if(rmass) massone = rmass[i];
    else massone = mass[type[i]];
    massone = massone/dm_lb; 

    gammavalue = 2.0*(mnode*massone)*dtoverdtcollision/(mnode+massone);
  }else{
    gammavalue = Gamma[type[i]];
  }
  

  for(k=0; k<3; k++){
    Ff[ix][iy][iz][k]    += gammavalue*((v[i][k]*dt_lb/dx_lb)-unode[k])*FfP[0];
    Ff[ix][iy][izp][k]   += gammavalue*((v[i][k]*dt_lb/dx_lb)-unode[k])*FfP[1];
    Ff[ix][iyp][iz][k]   += gammavalue*((v[i][k]*dt_lb/dx_lb)-unode[k])*FfP[2];
    Ff[ix][iyp][izp][k]  += gammavalue*((v[i][k]*dt_lb/dx_lb)-unode[k])*FfP[3];
    Ff[ixp][iy][iz][k]   += gammavalue*((v[i][k]*dt_lb/dx_lb)-unode[k])*FfP[4];
    Ff[ixp][iy][izp][k]  += gammavalue*((v[i][k]*dt_lb/dx_lb)-unode[k])*FfP[5];
    Ff[ixp][iyp][iz][k]  += gammavalue*((v[i][k]*dt_lb/dx_lb)-unode[k])*FfP[6];
    Ff[ixp][iyp][izp][k] += gammavalue*((v[i][k]*dt_lb/dx_lb)-unode[k])*FfP[7];
  }

  for(k=0; k<3; k++)
    hydroF[i][k] = -1.0*gammavalue*((v[i][k]*dt_lb/dx_lb)-unode[k])*dm_lb*dx_lb/dt_lb/dt_lb;

}

//==========================================================================
// read in a fluid restart file.  This is only used to restart the
// fluid portion of a LAMMPS simulation.  
//==========================================================================
void FixLbFluid::read_restartfile(void)
{
  MPI_Status status;
  MPI_Datatype realtype;
  MPI_Datatype filetype;


  int realsizes[4] = {subNbx,subNby,subNbz,numvel};
  int realstarts[4] = {1,1,1,0};
  int gsizes[4] = {Nbx,Nby,Nbz,numvel};
  int lsizes[4] = {subNbx-2,subNby-2,subNbz-2,numvel};
  int starts[4] = {comm->myloc[0]*(subNbx-2),comm->myloc[1]*(subNby-2),comm->myloc[2]*(subNbz-2),0};
  if(domain->periodicity[2]==0 && comm->myloc[2]==comm->procgrid[2]-1){
    starts[2] = comm->myloc[2]*(subNbz-3);
  }

  MPI_Type_create_subarray(4,realsizes,lsizes,realstarts,MPI_ORDER_C,MPI_DOUBLE,&realtype);
  MPI_Type_commit(&realtype);

  MPI_Type_create_subarray(4,gsizes,lsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&filetype);
  MPI_Type_commit(&filetype);

  MPI_File_set_view(pFileRead,0,MPI_DOUBLE,filetype,(char *) "native",
                    MPI_INFO_NULL);
  MPI_File_seek(pFileRead,0,MPI_SEEK_SET);
  MPI_File_read_all(pFileRead,&f_lb[0][0][0][0],1,realtype,&status);
  if(typeLB == 2){
    MPI_File_read_all(pFileRead,&feqold[0][0][0][0],1,realtype,&status);
    MPI_File_read_all(pFileRead,&feqoldn[0][0][0][0],1,realtype,&status);  
  }

  MPI_Type_free(&realtype);
  MPI_Type_free(&filetype);
  MPI_File_close(&pFileRead);

}

//==========================================================================
// write a fluid restart file.   
//==========================================================================
void FixLbFluid::write_restartfile(void)
{

  MPI_File fh;
  MPI_Status status;
  MPI_Datatype realtype;
  MPI_Datatype filetype;

  char *hfile;
  hfile = new char[32];
  sprintf(hfile,"FluidRestart_" BIGINT_FORMAT ".dat",update->ntimestep);
  
  MPI_File_open(world,hfile,MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,&fh);

  int realsizes[4] = {subNbx,subNby,subNbz,numvel};
  int realstarts[4] = {1,1,1,0};
  int gsizes[4] = {Nbx,Nby,Nbz,numvel};
  int lsizes[4] = {subNbx-2,subNby-2,subNbz-2,numvel};
  int starts[4] = {comm->myloc[0]*(subNbx-2),comm->myloc[1]*(subNby-2),comm->myloc[2]*(subNbz-2),0};
  if(domain->periodicity[2]==0 && comm->myloc[2]==comm->procgrid[2]-1){
    starts[2] = comm->myloc[2]*(subNbz-3);
  }

  MPI_Type_create_subarray(4,realsizes,lsizes,realstarts,MPI_ORDER_C,MPI_DOUBLE,&realtype);
  MPI_Type_commit(&realtype);

  MPI_Type_create_subarray(4,gsizes,lsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&filetype);
  MPI_Type_commit(&filetype);

  MPI_File_set_view(fh,0,MPI_DOUBLE,filetype,(char *) "native",MPI_INFO_NULL);
  MPI_File_write_all(fh,&f_lb[0][0][0][0],1,realtype,&status);
  if(typeLB == 2){
    MPI_File_write_all(fh,&feqold[0][0][0][0],1,realtype,&status);
    MPI_File_write_all(fh,&feqoldn[0][0][0][0],1,realtype,&status);  
  }

  MPI_Type_free(&realtype);
  MPI_Type_free(&filetype);
  MPI_File_close(&fh);
  delete [] hfile;

}

//==========================================================================
// rescale the simulation parameters so that dx_lb=dt_lb=dm_lb=1.
// This assumes that all the simulation parameters have been given in
// terms of distance, time and mass units. 
//==========================================================================
void FixLbFluid::rescale(void)
{
  vwtp = vwtp*dt_lb/dx_lb;
  vwbt = vwbt*dt_lb/dx_lb;
  
  bodyforcex = bodyforcex*dt_lb*dt_lb/dx_lb;
  bodyforcey = bodyforcey*dt_lb*dt_lb/dx_lb;
  bodyforcez = bodyforcez*dt_lb*dt_lb/dx_lb;
  
  tau=(3.0*viscosity/densityinit_real)*dt_lb*dt_lb/dx_lb/dx_lb;
  tau /= dt_lb;
  if(typeLB==1)
    tau = tau + 0.5;
   
  if(setGamma == 0){
    for(int i=0; i<= atom->ntypes; i++){
      NodeArea[i] = NodeArea[i]/dx_lb/dx_lb;
    }
  }else{
    for(int i=0; i<= atom->ntypes; i++){
      Gamma[i] = Gamma[i]*dt_lb/dm_lb;
    }
  }
   
  densityinit = densityinit_real*dx_lb*dx_lb*dx_lb/dm_lb;

  a_0 = a_0_real*dt_lb*dt_lb/(dx_lb*dx_lb);

  // Warn if using the D3Q19 model with noise, and a0 is too small.
  if(numvel==19 && noisestress==1 && a_0 < 0.2){
    error->warning(FLERR,"Fix lb/fluid WARNING: Chosen value for a0 may be too small. \
Check temperature reproduction.\n");
  }

  if(noisestress==1){
    if(a_0>0.5555555){
      error->all(FLERR,"Fix lb/fluid ERROR: the Lattice Boltzmann dx and dt need \
to be chosen such that the scaled a_0 < 5/9\n");
    }
  }

  // Courant Condition:
  if(a_0 >= 1.0){
    error->all(FLERR,"Fix lb/fluid ERROR: the lattice Boltzmann dx and dt do not \
satisfy the Courant condition.\n");
  }

  kB = (force->boltz/force->mvv2e)*dt_lb*dt_lb/dx_lb/dx_lb/dm_lb;

  if(typeLB==1){
    expminusdtovertau = 0.0;
    Dcoeff = 0.0;
    namp = 2.0*kB*T*(tau-0.5)/3.0;
    noisefactor = 1.0;
    if(a_0 <= 0.333333333333333){
      K_0 = 5.17*(0.333333333333333 - a_0);
    }else{
      K_0 = 2.57*(a_0 - 0.333333333333333);
    }
     dtoverdtcollision = dt_lb*6.0*viscosity/densityinit_real/dx_lb/dx_lb;
  }else if(typeLB==2){
    expminusdtovertau=exp(-1.0/tau);
    Dcoeff=(1.0-(1.0-expminusdtovertau)*tau);
    namp = 2.0*kB*T/3.;
    noisefactor=sqrt((1.0-expminusdtovertau*expminusdtovertau)/
		     (2.0))/(1.0-expminusdtovertau);
    K_0 = 4.5*(1.0/3.0-a_0);
    dtoverdtcollision = dt_lb*3.0*viscosity/densityinit_real/dx_lb/dx_lb;
  }

}

//==========================================================================
// Set the lattice-Boltzmann velocity vectors and weights for the D3Q15
// model.  Initialize the fluid velocity and density.
//==========================================================================
void FixLbFluid::initializeLB15(void)
{
  int i,j,k,m;

  //velocity vectors.
  e[0][0]= 0;
  e[0][1]= 0;
  e[0][2]= 0;

  e[1][0]= 1;
  e[1][1]= 0;
  e[1][2]= 0;

  e[2][0]= 0;
  e[2][1]= 1;
  e[2][2]= 0;

  e[3][0]= -1;
  e[3][1]= 0;
  e[3][2]= 0;

  e[4][0]= 0;
  e[4][1]= -1;
  e[4][2]= 0;

  e[5][0]= 0;
  e[5][1]= 0;
  e[5][2]= 1;

  e[6][0]= 0;
  e[6][1]= 0;
  e[6][2]= -1;

  e[7][0]= 1;
  e[7][1]= 1;
  e[7][2]= 1;

  e[8][0]= -1;
  e[8][1]= 1;
  e[8][2]= 1;

  e[9][0]= -1;
  e[9][1]= -1;
  e[9][2]= 1;

  e[10][0]= 1;
  e[10][1]= -1;
  e[10][2]= 1;

  e[11][0]= 1;
  e[11][1]= 1;
  e[11][2]= -1;

  e[12][0]= -1;
  e[12][1]= 1;
  e[12][2]= -1;

  e[13][0]= -1;
  e[13][1]= -1;
  e[13][2]= -1;

  e[14][0]= 1;
  e[14][1]= -1;
  e[14][2]= -1;

  //weights.
  w_lb[0]=2./9.;
  w_lb[1]=1./9.;
  w_lb[2]=1./9.;
  w_lb[3]=1./9.;
  w_lb[4]=1./9.;
  w_lb[5]=1./9.;
  w_lb[6]=1./9.;
  w_lb[7]=1./72.;
  w_lb[8]=1./72.;
  w_lb[9]=1./72.;
  w_lb[10]=1./72.;
  w_lb[11]=1./72.;
  w_lb[12]=1./72.;
  w_lb[13]=1./72.;
  w_lb[14]=1./72.;

  Ng_lb[0]=1.;
  Ng_lb[1]=3.;
  Ng_lb[2]=3.;
  Ng_lb[3]=3.;
  Ng_lb[4]=9./2.;
  Ng_lb[5]=9./2.;
  Ng_lb[6]=9./2.;
  Ng_lb[7]=9.;
  Ng_lb[8]=9.;
  Ng_lb[9]=9.;
  Ng_lb[10]=27./2.;
  Ng_lb[11]=27./2.;
  Ng_lb[12]=27./2.;
  Ng_lb[13]=9.;
  Ng_lb[14]=1.;

  mg_lb[0][0]=1.;mg_lb[0][1]=1.;mg_lb[0][2]=1.;mg_lb[0][3]=1.;mg_lb[0][4]=1.;
  mg_lb[0][5]=1.;mg_lb[0][6]=1.;mg_lb[0][7]=1.;mg_lb[0][8]=1.;mg_lb[0][9]=1.;
  mg_lb[0][10]=1.;mg_lb[0][11]=1.;mg_lb[0][12]=1.;mg_lb[0][13]=1.;mg_lb[0][14]=1.;
  mg_lb[1][0]=0;mg_lb[1][1]=1.;mg_lb[1][2]=0;mg_lb[1][3]=-1.;mg_lb[1][4]=0;
  mg_lb[1][5]=0;mg_lb[1][6]=0;mg_lb[1][7]=1.;mg_lb[1][8]=-1.;mg_lb[1][9]=-1.;
  mg_lb[1][10]=1.;mg_lb[1][11]=1.;mg_lb[1][12]=-1.;mg_lb[1][13]=-1.;mg_lb[1][14]=1.;
  mg_lb[2][0]=0;mg_lb[2][1]=0;mg_lb[2][2]=1.;mg_lb[2][3]=0;mg_lb[2][4]=-1.;
  mg_lb[2][5]=0;mg_lb[2][6]=0;mg_lb[2][7]=1.;mg_lb[2][8]=1.;mg_lb[2][9]=-1.;
  mg_lb[2][10]=-1.;mg_lb[2][11]=1.;mg_lb[2][12]=1.;mg_lb[2][13]=-1.;mg_lb[2][14]=-1.;
  mg_lb[3][0]=0;mg_lb[3][1]=0;mg_lb[3][2]=0;mg_lb[3][3]=0;mg_lb[3][4]=0;
  mg_lb[3][5]=1.;mg_lb[3][6]=-1.;mg_lb[3][7]=1.;mg_lb[3][8]=1.;mg_lb[3][9]=1.;
  mg_lb[3][10]=1.;mg_lb[3][11]=-1.;mg_lb[3][12]=-1.;mg_lb[3][13]=-1.;mg_lb[3][14]=-1.;
  mg_lb[4][0]=-1./3.;mg_lb[4][1]=2./3.;mg_lb[4][2]=-1./3.;mg_lb[4][3]=2./3.;mg_lb[4][4]=-1./3.;
  mg_lb[4][5]=-1./3.;mg_lb[4][6]=-1./3.;mg_lb[4][7]=2./3.;mg_lb[4][8]=2./3.;mg_lb[4][9]=2./3.;
  mg_lb[4][10]=2./3.;mg_lb[4][11]=2./3.;mg_lb[4][12]=2./3.;mg_lb[4][13]=2./3.;mg_lb[4][14]=2./3.;
  mg_lb[5][0]=-1./3.;mg_lb[5][1]=-1./3.;mg_lb[5][2]=2./3.;mg_lb[5][3]=-1./3.;mg_lb[5][4]=2./3.;
  mg_lb[5][5]=-1./3.;mg_lb[5][6]=-1./3.;mg_lb[5][7]=2./3.;mg_lb[5][8]=2./3.;mg_lb[5][9]=2./3.;
  mg_lb[5][10]=2./3.;mg_lb[5][11]=2./3.;mg_lb[5][12]=2./3.;mg_lb[5][13]=2./3.;mg_lb[5][14]=2./3.;
  mg_lb[6][0]=-1./3.;mg_lb[6][1]=-1./3.;mg_lb[6][2]=-1./3.;mg_lb[6][3]=-1./3.;mg_lb[6][4]=-1./3.;
  mg_lb[6][5]=2./3.;mg_lb[6][6]=2./3.;mg_lb[6][7]=2./3.;mg_lb[6][8]=2./3.;mg_lb[6][9]=2./3.;
  mg_lb[6][10]=2./3.;mg_lb[6][11]=2./3.;mg_lb[6][12]=2./3.;mg_lb[6][13]=2./3.;mg_lb[6][14]=2./3.;
  mg_lb[7][0]=0;mg_lb[7][1]=0;mg_lb[7][2]=0;mg_lb[7][3]=0;mg_lb[7][4]=0;
  mg_lb[7][5]=0;mg_lb[7][6]=0;mg_lb[7][7]=1;mg_lb[7][8]=-1;mg_lb[7][9]=1;
  mg_lb[7][10]=-1;mg_lb[7][11]=1;mg_lb[7][12]=-1;mg_lb[7][13]=1;mg_lb[7][14]=-1; 
  mg_lb[8][0]=0;mg_lb[8][1]=0;mg_lb[8][2]=0;mg_lb[8][3]=0;mg_lb[8][4]=0;
  mg_lb[8][5]=0;mg_lb[8][6]=0;mg_lb[8][7]=1;mg_lb[8][8]=1;mg_lb[8][9]=-1;
  mg_lb[8][10]=-1;mg_lb[8][11]=-1;mg_lb[8][12]=-1;mg_lb[8][13]=1;mg_lb[8][14]=1; 
  mg_lb[9][0]=0;mg_lb[9][1]=0;mg_lb[9][2]=0;mg_lb[9][3]=0;mg_lb[9][4]=0;
  mg_lb[9][5]=0;mg_lb[9][6]=0;mg_lb[9][7]=1;mg_lb[9][8]=-1;mg_lb[9][9]=-1;
  mg_lb[9][10]=1;mg_lb[9][11]=-1;mg_lb[9][12]=1;mg_lb[9][13]=1;mg_lb[9][14]=-1; 
  mg_lb[10][0]=0;mg_lb[10][1]=0;mg_lb[10][2]=-1./3.;mg_lb[10][3]=0;mg_lb[10][4]=1./3.;
  mg_lb[10][5]=0;mg_lb[10][6]=0;mg_lb[10][7]=2./3.;mg_lb[10][8]=2./3.;mg_lb[10][9]=-2./3.;
  mg_lb[10][10]=-2./3.;mg_lb[10][11]=2./3.;mg_lb[10][12]=2./3.;mg_lb[10][13]=-2./3.;mg_lb[10][14]=-2./3.;
  mg_lb[11][0]=0;mg_lb[11][1]=0;mg_lb[11][2]=0;mg_lb[11][3]=0;mg_lb[11][4]=0;
  mg_lb[11][5]=-1./3.;mg_lb[11][6]=1./3.;mg_lb[11][7]=2./3.;mg_lb[11][8]=2./3.;mg_lb[11][9]=2./3.;
  mg_lb[11][10]=2./3.;mg_lb[11][11]=-2./3.;mg_lb[11][12]=-2./3.;mg_lb[11][13]=-2./3.;mg_lb[11][14]=-2./3.;
  mg_lb[12][0]=0;mg_lb[12][1]=-1./3.;mg_lb[12][2]=0;mg_lb[12][3]=1./3.;mg_lb[12][4]=0;
  mg_lb[12][5]=0;mg_lb[12][6]=0;mg_lb[12][7]=2./3.;mg_lb[12][8]=-2./3.;mg_lb[12][9]=-2./3.;
  mg_lb[12][10]=2./3.;mg_lb[12][11]=2./3.;mg_lb[12][12]=-2./3.;mg_lb[12][13]=-2./3.;mg_lb[12][14]=2./3.;
  mg_lb[13][0]=0;mg_lb[13][1]=0;mg_lb[13][2]=0;mg_lb[13][3]=0;mg_lb[13][4]=0;
  mg_lb[13][5]=0;mg_lb[13][6]=0;mg_lb[13][7]=1;mg_lb[13][8]=-1;mg_lb[13][9]=1;
  mg_lb[13][10]=-1;mg_lb[13][11]=-1;mg_lb[13][12]=1;mg_lb[13][13]=-1;mg_lb[13][14]=1;
  mg_lb[14][0]=sqrt(2.);mg_lb[14][1]=-1./sqrt(2.);mg_lb[14][2]=-1./sqrt(2.);
  mg_lb[14][3]=-1./sqrt(2.);mg_lb[14][4]=-1./sqrt(2.);
  mg_lb[14][5]=-1./sqrt(2.);mg_lb[14][6]=-1./sqrt(2.);mg_lb[14][7]=sqrt(2.);
  mg_lb[14][8]=sqrt(2.);mg_lb[14][9]=sqrt(2.);
  mg_lb[14][10]=sqrt(2.);mg_lb[14][11]=sqrt(2.);mg_lb[14][12]=sqrt(2.);
  mg_lb[14][13]=sqrt(2.);mg_lb[14][14]=sqrt(2.);

  for(i=0; i<subNbx+3; i++)
    for(j=0; j<subNby+3; j++)
      for(k=0; k<subNbz+3; k++){
	u_lb[i][j][k][0]=0.0;
	u_lb[i][j][k][1]=0.0;
	u_lb[i][j][k][2]=0.0;
	density_lb[i][j][k] = densityinit;
  }
  for(i=0; i<subNbx; i++)
    for(j=0; j<subNby; j++)
      for(k=0; k<subNbz; k++)
	for(m=0; m<15; m++)
	  f_lb[i][j][k][m] = density_lb[i][j][k]/15.0;

}

//==========================================================================
// Set the lattice-Boltzmann velocity vectors and weights for the D3Q19
// model.  Initialize the fluid velocity and density.
//==========================================================================
void FixLbFluid::initializeLB19(void)
{
  int i,j,k,m;

  //velocity vectors.
  e[0][0]= 0;
  e[0][1]= 0;
  e[0][2]= 0;

  e[1][0]= 1;
  e[1][1]= 0;
  e[1][2]= 0;

  e[2][0]= 0;
  e[2][1]= 1;
  e[2][2]= 0;

  e[3][0]= -1;
  e[3][1]= 0;
  e[3][2]= 0;

  e[4][0]= 0;
  e[4][1]= -1;
  e[4][2]= 0;

  e[5][0]= 0;
  e[5][1]= 0;
  e[5][2]= 1;

  e[6][0]= 0;
  e[6][1]= 0;
  e[6][2]= -1;

  e[7][0] = 1;
  e[7][1] = 1;
  e[7][2] = 0;

  e[8][0] = 1;
  e[8][1] = -1;
  e[8][2] = 0;

  e[9][0] = -1;
  e[9][1] = 1;
  e[9][2] = 0;  

  e[10][0] = -1;
  e[10][1] = -1;
  e[10][2] = 0;

  e[11][0] = 1;
  e[11][1] = 0;
  e[11][2] = 1;

  e[12][0] = 1;
  e[12][1] = 0;
  e[12][2] = -1;

  e[13][0] = -1;
  e[13][1] = 0;
  e[13][2] = 1;

  e[14][0] = -1;
  e[14][1] = 0;
  e[14][2] = -1;

  e[15][0] = 0;
  e[15][1] = 1;
  e[15][2] = 1;

  e[16][0] = 0;
  e[16][1] = 1;
  e[16][2] = -1;

  e[17][0] = 0;
  e[17][1] = -1;
  e[17][2] = 1;

  e[18][0] = 0;
  e[18][1] = -1;
  e[18][2] = -1;

  //weights.
  w_lb[0]=1./3.;
  w_lb[1]=1./18.;
  w_lb[2]=1./18.;
  w_lb[3]=1./18.;
  w_lb[4]=1./18.;
  w_lb[5]=1./18.;
  w_lb[6]=1./18.;
  w_lb[7]=1./36.;
  w_lb[8]=1./36.;
  w_lb[9]=1./36.;
  w_lb[10]=1./36.;
  w_lb[11]=1./36.;
  w_lb[12]=1./36.;
  w_lb[13]=1./36.;
  w_lb[14]=1./36.;
  w_lb[15]=1./36.;
  w_lb[16]=1./36.;
  w_lb[17]=1./36.;
  w_lb[18]=1./36.;

  Ng_lb[0]=1.;
  Ng_lb[1]=3.;
  Ng_lb[2]=3.;
  Ng_lb[3]=3.;
  Ng_lb[4]=9./2.;
  Ng_lb[5]=9./2.;
  Ng_lb[6]=9./2.;
  Ng_lb[7]=9.;
  Ng_lb[8]=9.;
  Ng_lb[9]=9.;
  Ng_lb[10]=27./2.;
  Ng_lb[11]=27./2.;
  Ng_lb[12]=27./2.;
  Ng_lb[13]=18.;
  Ng_lb[14]=18.;
  Ng_lb[15]=18.;
  Ng_lb[16]=162./7.;
  Ng_lb[17]=126./5.;
  Ng_lb[18]=30.;

  mg_lb[0][0] = 1.; mg_lb[0][1] = 1.; mg_lb[0][2] = 1.; mg_lb[0][3] = 1.; mg_lb[0][4] = 1.;
  mg_lb[0][5] = 1.; mg_lb[0][6] = 1.; mg_lb[0][7] = 1.; mg_lb[0][8] = 1.; mg_lb[0][9] = 1.;
  mg_lb[0][10]= 1.; mg_lb[0][11]= 1.; mg_lb[0][12]= 1.; mg_lb[0][13]= 1.; mg_lb[0][14]= 1.;
  mg_lb[0][15]= 1.; mg_lb[0][16]= 1.; mg_lb[0][17]= 1.; mg_lb[0][18]= 1.;

  mg_lb[1][0] = 0.; mg_lb[1][1] = 1.; mg_lb[1][2] = 0.; mg_lb[1][3] =-1.; mg_lb[1][4] = 0.;
  mg_lb[1][5] = 0.; mg_lb[1][6] = 0.; mg_lb[1][7] = 1.; mg_lb[1][8] = 1.; mg_lb[1][9] =-1.;
  mg_lb[1][10]=-1.; mg_lb[1][11]= 1.; mg_lb[1][12]= 1.; mg_lb[1][13]=-1.; mg_lb[1][14]=-1.;
  mg_lb[1][15]= 0.; mg_lb[1][16]= 0.; mg_lb[1][17]= 0.; mg_lb[1][18]= 0.;

  mg_lb[2][0] = 0.; mg_lb[2][1] = 0.; mg_lb[2][2] = 1.; mg_lb[2][3] = 0.; mg_lb[2][4] =-1.;
  mg_lb[2][5] = 0.; mg_lb[2][6] = 0.; mg_lb[2][7] = 1.; mg_lb[2][8] =-1.; mg_lb[2][9] = 1.;
  mg_lb[2][10]=-1.; mg_lb[2][11]= 0.; mg_lb[2][12]= 0.; mg_lb[2][13]= 0.; mg_lb[2][14]= 0.;
  mg_lb[2][15]= 1.; mg_lb[2][16]= 1.; mg_lb[2][17]=-1.; mg_lb[2][18]=-1.;  

  mg_lb[3][0] = 0.; mg_lb[3][1] = 0.; mg_lb[3][2] = 0.; mg_lb[3][3] = 0.; mg_lb[3][4] = 0.;
  mg_lb[3][5] = 1.; mg_lb[3][6] =-1.; mg_lb[3][7] = 0.; mg_lb[3][8] = 0.; mg_lb[3][9] = 0.;
  mg_lb[3][10]= 0.; mg_lb[3][11]= 1.; mg_lb[3][12]=-1.; mg_lb[3][13]= 1.; mg_lb[3][14]=-1.;
  mg_lb[3][15]= 1.; mg_lb[3][16]=-1.; mg_lb[3][17]= 1.; mg_lb[3][18]=-1.;

  mg_lb[4][0] =-1./3.; mg_lb[4][1] = 2./3.; mg_lb[4][2] =-1./3.; mg_lb[4][3] = 2./3.; mg_lb[4][4] =-1./3.;
  mg_lb[4][5] =-1./3.; mg_lb[4][6] =-1./3.; mg_lb[4][7] = 2./3.; mg_lb[4][8] = 2./3.; mg_lb[4][9] = 2./3.;
  mg_lb[4][10]= 2./3.; mg_lb[4][11]= 2./3.; mg_lb[4][12]= 2./3.; mg_lb[4][13]= 2./3.; mg_lb[4][14]= 2./3.;
  mg_lb[4][15]=-1./3.; mg_lb[4][16]=-1./3.; mg_lb[4][17]=-1./3.; mg_lb[4][18]=-1./3.;

  mg_lb[5][0] =-1./3.; mg_lb[5][1] =-1./3.; mg_lb[5][2] = 2./3.; mg_lb[5][3] =-1./3.; mg_lb[5][4] = 2./3.;
  mg_lb[5][5] =-1./3.; mg_lb[5][6] =-1./3.; mg_lb[5][7] = 2./3.; mg_lb[5][8] = 2./3.; mg_lb[5][9] = 2./3.;
  mg_lb[5][10]= 2./3.; mg_lb[5][11]=-1./3.; mg_lb[5][12]=-1./3.; mg_lb[5][13]=-1./3.; mg_lb[5][14]=-1./3.;
  mg_lb[5][15]= 2./3.; mg_lb[5][16]= 2./3.; mg_lb[5][17]= 2./3.; mg_lb[5][18]= 2./3.;

  mg_lb[6][0] =-1./3.; mg_lb[6][1] =-1./3.; mg_lb[6][2] =-1./3.; mg_lb[6][3] =-1./3.; mg_lb[6][4] =-1./3.;
  mg_lb[6][5] = 2./3.; mg_lb[6][6] = 2./3.; mg_lb[6][7] =-1./3.; mg_lb[6][8] =-1./3.; mg_lb[6][9] =-1./3.;
  mg_lb[6][10]=-1./3.; mg_lb[6][11]= 2./3.; mg_lb[6][12]= 2./3.; mg_lb[6][13]= 2./3.; mg_lb[6][14]= 2./3.;
  mg_lb[6][15]= 2./3.; mg_lb[6][16]= 2./3.; mg_lb[6][17]= 2./3.; mg_lb[6][18]= 2./3.;

  mg_lb[7][0] = 0.; mg_lb[7][1] = 0.; mg_lb[7][2] = 0.; mg_lb[7][3] = 0.; mg_lb[7][4] = 0.;
  mg_lb[7][5] = 0.; mg_lb[7][6] = 0.; mg_lb[7][7] = 1.; mg_lb[7][8] =-1.; mg_lb[7][9] =-1.;
  mg_lb[7][10]= 1.; mg_lb[7][11]= 0.; mg_lb[7][12]= 0.; mg_lb[7][13]= 0.; mg_lb[7][14]= 0.;
  mg_lb[7][15]= 0.; mg_lb[7][16]= 0.; mg_lb[7][17]= 0.; mg_lb[7][18]= 0.;

  mg_lb[8][0] = 0.; mg_lb[8][1] = 0.; mg_lb[8][2] = 0.; mg_lb[8][3] = 0.; mg_lb[8][4] = 0.;
  mg_lb[8][5] = 0.; mg_lb[8][6] = 0.; mg_lb[8][7] = 0.; mg_lb[8][8] = 0.; mg_lb[8][9] = 0.;
  mg_lb[8][10]= 0.; mg_lb[8][11]= 1.; mg_lb[8][12]=-1.; mg_lb[8][13]=-1.; mg_lb[8][14]= 1.;
  mg_lb[8][15]= 0.; mg_lb[8][16]= 0.; mg_lb[8][17]= 0.; mg_lb[8][18]= 0.;

  mg_lb[9][0] = 0.; mg_lb[9][1] = 0.; mg_lb[9][2] = 0.; mg_lb[9][3] = 0.; mg_lb[9][4] = 0.;
  mg_lb[9][5] = 0.; mg_lb[9][6] = 0.; mg_lb[9][7] = 0.; mg_lb[9][8] = 0.; mg_lb[9][9] = 0.;
  mg_lb[9][10]= 0.; mg_lb[9][11]= 0.; mg_lb[9][12]= 0.; mg_lb[9][13]= 0.; mg_lb[9][14]= 0.;
  mg_lb[9][15]= 1.; mg_lb[9][16]=-1.; mg_lb[9][17]=-1.; mg_lb[9][18]= 1.;

  mg_lb[10][0] = 0.;    mg_lb[10][1] =-1./3.; mg_lb[10][2] = 0.;    mg_lb[10][3] = 1./3.; mg_lb[10][4] = 0.;
  mg_lb[10][5] = 0.;    mg_lb[10][6] = 0.;    mg_lb[10][7] = 2./3.; mg_lb[10][8] = 2./3.; mg_lb[10][9] =-2./3.;
  mg_lb[10][10]=-2./3.; mg_lb[10][11]=-1./3.; mg_lb[10][12]=-1./3.; mg_lb[10][13]= 1./3.; mg_lb[10][14]= 1./3.;
  mg_lb[10][15]= 0.;    mg_lb[10][16]= 0.;    mg_lb[10][17]= 0.;    mg_lb[10][18]= 0.;

  mg_lb[11][0] = 0.;    mg_lb[11][1] = 0.;    mg_lb[11][2] =-1./3.; mg_lb[11][3] = 0.;    mg_lb[11][4] = 1./3.;
  mg_lb[11][5] = 0.;    mg_lb[11][6] = 0.;    mg_lb[11][7] = 2./3.; mg_lb[11][8] =-2./3.; mg_lb[11][9] = 2./3.;
  mg_lb[11][10]=-2./3.; mg_lb[11][11]= 0.;    mg_lb[11][12]= 0.;    mg_lb[11][13]= 0.;    mg_lb[11][14]= 0.;
  mg_lb[11][15]=-1./3.; mg_lb[11][16]=-1./3.; mg_lb[11][17]= 1./3.; mg_lb[11][18]= 1./3.;

  mg_lb[12][0] = 0.;    mg_lb[12][1] = 0.;    mg_lb[12][2] = 0.;    mg_lb[12][3] = 0.;    mg_lb[12][4] = 0.;
  mg_lb[12][5] =-1./3.; mg_lb[12][6] = 1./3.; mg_lb[12][7] = 0.;    mg_lb[12][8] = 0.;    mg_lb[12][9] = 0.;
  mg_lb[12][10]= 0.;    mg_lb[12][11]= 2./3.; mg_lb[12][12]=-2./3.; mg_lb[12][13]= 2./3.; mg_lb[12][14]=-2./3.;
  mg_lb[12][15]=-1./3.; mg_lb[12][16]= 1./3.; mg_lb[12][17]=-1./3.; mg_lb[12][18]= 1./3.;

  mg_lb[13][0] = 0.; mg_lb[13][1] =-0.5; mg_lb[13][2] = 0.;  mg_lb[13][3] = 0.5; mg_lb[13][4] = 0.;
  mg_lb[13][5] = 0.; mg_lb[13][6] = 0.;  mg_lb[13][7] = 0.;  mg_lb[13][8] = 0.;  mg_lb[13][9] = 0.;
  mg_lb[13][10]= 0.; mg_lb[13][11]= 0.5; mg_lb[13][12]= 0.5; mg_lb[13][13]=-0.5; mg_lb[13][14]=-0.5;
  mg_lb[13][15]= 0.; mg_lb[13][16]= 0.;  mg_lb[13][17]= 0.;  mg_lb[13][18]= 0.;

  mg_lb[14][0] = 0.;  mg_lb[14][1] = 0.;  mg_lb[14][2] = 0.;  mg_lb[14][3] = 0.;  mg_lb[14][4] = 0.;
  mg_lb[14][5] =-0.5; mg_lb[14][6] = 0.5; mg_lb[14][7] = 0.;  mg_lb[14][8] = 0.;  mg_lb[14][9] = 0.;
  mg_lb[14][10]= 0.;  mg_lb[14][11]= 0.;  mg_lb[14][12]= 0.;  mg_lb[14][13]= 0.;  mg_lb[14][14]= 0.;
  mg_lb[14][15]= 0.5; mg_lb[14][16]=-0.5; mg_lb[14][17]= 0.5; mg_lb[14][18]=-0.5;

  mg_lb[15][0] = 0.;  mg_lb[15][1] = 0.;  mg_lb[15][2] =-0.5; mg_lb[15][3] = 0.;  mg_lb[15][4] = 0.5;
  mg_lb[15][5] = 0.;  mg_lb[15][6] = 0.;  mg_lb[15][7] = 0.;  mg_lb[15][8] = 0.;  mg_lb[15][9] = 0.;
  mg_lb[15][10]= 0.;  mg_lb[15][11]= 0.;  mg_lb[15][12]= 0.;  mg_lb[15][13]= 0.;  mg_lb[15][14]= 0.;
  mg_lb[15][15]= 0.5; mg_lb[15][16]= 0.5; mg_lb[15][17]=-0.5; mg_lb[15][18]=-0.5;

  mg_lb[16][0] = 1./18.; mg_lb[16][1] =-5./18.; mg_lb[16][2] =-5./18.; mg_lb[16][3] =-5./18.; mg_lb[16][4] =-5./18.;
  mg_lb[16][5] = 2./9.;  mg_lb[16][6] = 2./9.;  mg_lb[16][7] = 7./18.; mg_lb[16][8] = 7./18.; mg_lb[16][9] = 7./18.;
  mg_lb[16][10]= 7./18.; mg_lb[16][11]=-1./9.;  mg_lb[16][12]=-1./9.;  mg_lb[16][13]=-1./9.;  mg_lb[16][14]=-1./9.;
  mg_lb[16][15]=-1./9.;  mg_lb[16][16]=-1./9.;  mg_lb[16][17]=-1./9.;  mg_lb[16][18]=-1./9.;

  mg_lb[17][0] = 1./14.; mg_lb[17][1] =-5./14.; mg_lb[17][2] = 1./7.;  mg_lb[17][3] =-5./14.; mg_lb[17][4] = 1./7.;
  mg_lb[17][5] =-3./14.; mg_lb[17][6] =-3./14.; mg_lb[17][7] = 0.;     mg_lb[17][8] = 0.;     mg_lb[17][9] = 0.;
  mg_lb[17][10]= 0.;     mg_lb[17][11]= 5./14.; mg_lb[17][12]= 5./14.; mg_lb[17][13]= 5./14.; mg_lb[17][14]= 5./14.;
  mg_lb[17][15]=-1./7.;  mg_lb[17][16]=-1./7.;  mg_lb[17][17]=-1./7.;  mg_lb[17][18]=-1./7.;

  mg_lb[18][0] = 1./10.; mg_lb[18][1] = 0.;     mg_lb[18][2] =-3./10.; mg_lb[18][3] = 0.;    mg_lb[18][4] =-3./10.;
  mg_lb[18][5] =-3./10.; mg_lb[18][6] =-3./10.; mg_lb[18][7] = 0.;     mg_lb[18][8] = 0.;    mg_lb[18][9] = 0.;
  mg_lb[18][10]= 0.;     mg_lb[18][11]= 0.;     mg_lb[18][12]= 0.;     mg_lb[18][13]= 0.;    mg_lb[18][14]= 0.;
  mg_lb[18][15]= 3./10.; mg_lb[18][16]= 3./10.; mg_lb[18][17]= 3./10.; mg_lb[18][18]= 3./10.;

  for(i=0; i<subNbx+3; i++)
    for(j=0; j<subNby+3; j++)
      for(k=0; k<subNbz+3; k++){
	u_lb[i][j][k][0]=0.0;
	u_lb[i][j][k][1]=0.0;
	u_lb[i][j][k][2]=0.0;
	density_lb[i][j][k] = densityinit;
  }
  for(i=0; i<subNbx; i++)
    for(j=0; j<subNby; j++)
      for(k=0; k<subNbz; k++)
	for(m=0; m<19; m++)
	  f_lb[i][j][k][m] = density_lb[i][j][k]/19.0;

}

//==========================================================================
// Initialize the equilibrium distribution functions 
// (this just uses the initial fluid parameters, and assumes no forces).
//==========================================================================
void FixLbFluid::initialize_feq(void)
{
  int i,j,k,p;
  MPI_Request requests[8];
  int numrequests;

  // If using the standary LB integrator, do not need to send feqn.
  if(typeLB == 1){
    numrequests = 4;
  }else{
    numrequests = 8;
  }

  std::fill(&Ff[0][0][0][0],&Ff[0][0][0][0] + (subNbx+3)*(subNby+3)*(subNbz+3)*3,0.0);
  std::fill(&Fftempx[0][0][0][0],&Fftempx[0][0][0][0] + 5*(subNby+3)*(subNbz+3)*3,0.0);
  std::fill(&Fftempy[0][0][0][0],&Fftempy[0][0][0][0] + (subNbx+3)*5*(subNbz+3)*3,0.0);
  std::fill(&Fftempz[0][0][0][0],&Fftempz[0][0][0][0] + (subNbx+3)*(subNby+3)*5*3,0.0);  

  if(readrestart == 0){
    step=0;

    parametercalc_full();
    (*this.*equilibriumdist)(1,subNbx-1,1,subNby-1,1,subNbz-1);  

    for(i=0; i<numrequests; i++)
      requests[i]=MPI_REQUEST_NULL;
    MPI_Isend(&feq[1][1][1][0],1,passxf,comm->procneigh[0][0],15,world,&requests[0]);
    MPI_Irecv(&feq[0][1][1][0],1,passxf,comm->procneigh[0][0],25,world,&requests[1]);
    MPI_Isend(&feq[subNbx-2][1][1][0],1,passxf,comm->procneigh[0][1],25,world,&requests[2]);
    MPI_Irecv(&feq[subNbx-1][1][1][0],1,passxf,comm->procneigh[0][1],15,world,&requests[3]);
    if(typeLB == 2){
      MPI_Isend(&feqn[1][1][1][0],1,passxf,comm->procneigh[0][0],10,world,&requests[4]);
      MPI_Irecv(&feqn[0][1][1][0],1,passxf,comm->procneigh[0][0],20,world,&requests[5]);
      MPI_Isend(&feqn[subNbx-2][1][1][0],1,passxf,comm->procneigh[0][1],20,world,&requests[6]);
      MPI_Irecv(&feqn[subNbx-1][1][1][0],1,passxf,comm->procneigh[0][1],10,world,&requests[7]);
    }  
    MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);
    
    for(i=0; i<numrequests; i++)
      requests[i]=MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][1][1][0],1,passyf,comm->procneigh[1][0],15,world,&requests[0]);
    MPI_Irecv(&feq[0][0][1][0],1,passyf,comm->procneigh[1][0],25,world,&requests[1]);   
    MPI_Isend(&feq[0][subNby-2][1][0],1,passyf,comm->procneigh[1][1],25,world,&requests[2]);
    MPI_Irecv(&feq[0][subNby-1][1][0],1,passyf,comm->procneigh[1][1],15,world,&requests[3]);
    if(typeLB == 2){
      MPI_Isend(&feqn[0][1][1][0],1,passyf,comm->procneigh[1][0],10,world,&requests[4]);
      MPI_Irecv(&feqn[0][0][1][0],1,passyf,comm->procneigh[1][0],20,world,&requests[5]);   
      MPI_Isend(&feqn[0][subNby-2][1][0],1,passyf,comm->procneigh[1][1],20,world,&requests[6]);
      MPI_Irecv(&feqn[0][subNby-1][1][0],1,passyf,comm->procneigh[1][1],10,world,&requests[7]);
    }
    MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);
    
    for(i=0; i<numrequests; i++)
      requests[i]=MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][0][1][0],1,passzf,comm->procneigh[2][0],15,world,&requests[0]);
    MPI_Irecv(&feq[0][0][0][0],1,passzf,comm->procneigh[2][0],25,world,&requests[1]);
    MPI_Isend(&feq[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],25,world,&requests[2]);
    MPI_Irecv(&feq[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],15,world,&requests[3]);
    if(typeLB == 2){
      MPI_Isend(&feqn[0][0][1][0],1,passzf,comm->procneigh[2][0],10,world,&requests[4]);
      MPI_Irecv(&feqn[0][0][0][0],1,passzf,comm->procneigh[2][0],20,world,&requests[5]);
      MPI_Isend(&feqn[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],20,world,&requests[6]);
      MPI_Irecv(&feqn[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],10,world,&requests[7]);
    } 
    MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);
    
    //Save feqold.
    if(typeLB == 2){
      for(i=0; i<subNbx; i++)
	for(j=0; j<subNby; j++)
	  for(k=0; k<subNbz; k++)
	    for(p=0; p<numvel; p++){
	      feqold[i][j][k][p] = feq[i][j][k][p];
	      feqoldn[i][j][k][p] = feqn[i][j][k][p];
	    }
    }
  }else{
    step = 1;
    
    read_restartfile();
    
    if(typeLB == 2){
      for(i=0; i<8; i++)
	requests[i]=MPI_REQUEST_NULL;
      MPI_Isend(&feqold[1][1][1][0],1,passxf,comm->procneigh[0][0],15,world,&requests[0]);
      MPI_Irecv(&feqold[0][1][1][0],1,passxf,comm->procneigh[0][0],25,world,&requests[1]);
      MPI_Isend(&feqold[subNbx-2][1][1][0],1,passxf,comm->procneigh[0][1],25,world,&requests[2]);
      MPI_Irecv(&feqold[subNbx-1][1][1][0],1,passxf,comm->procneigh[0][1],15,world,&requests[3]);
      MPI_Isend(&feqoldn[1][1][1][0],1,passxf,comm->procneigh[0][0],10,world,&requests[4]);
      MPI_Irecv(&feqoldn[0][1][1][0],1,passxf,comm->procneigh[0][0],20,world,&requests[5]);
      MPI_Isend(&feqoldn[subNbx-2][1][1][0],1,passxf,comm->procneigh[0][1],20,world,&requests[6]);
      MPI_Irecv(&feqoldn[subNbx-1][1][1][0],1,passxf,comm->procneigh[0][1],10,world,&requests[7]);  
      MPI_Waitall(8,requests,MPI_STATUS_IGNORE);
      
      for(i=0; i<8; i++)
	requests[i]=MPI_REQUEST_NULL;
      MPI_Isend(&feqold[0][1][1][0],1,passyf,comm->procneigh[1][0],15,world,&requests[0]);
      MPI_Irecv(&feqold[0][0][1][0],1,passyf,comm->procneigh[1][0],25,world,&requests[1]);   
      MPI_Isend(&feqold[0][subNby-2][1][0],1,passyf,comm->procneigh[1][1],25,world,&requests[2]);
      MPI_Irecv(&feqold[0][subNby-1][1][0],1,passyf,comm->procneigh[1][1],15,world,&requests[3]);
      MPI_Isend(&feqoldn[0][1][1][0],1,passyf,comm->procneigh[1][0],10,world,&requests[4]);
      MPI_Irecv(&feqoldn[0][0][1][0],1,passyf,comm->procneigh[1][0],20,world,&requests[5]);   
      MPI_Isend(&feqoldn[0][subNby-2][1][0],1,passyf,comm->procneigh[1][1],20,world,&requests[6]);
      MPI_Irecv(&feqoldn[0][subNby-1][1][0],1,passyf,comm->procneigh[1][1],10,world,&requests[7]);
      MPI_Waitall(8,requests,MPI_STATUS_IGNORE);
      
      for(i=0; i<8; i++)
	requests[i]=MPI_REQUEST_NULL;
      MPI_Isend(&feqold[0][0][1][0],1,passzf,comm->procneigh[2][0],15,world,&requests[0]);
      MPI_Irecv(&feqold[0][0][0][0],1,passzf,comm->procneigh[2][0],25,world,&requests[1]);
      MPI_Isend(&feqold[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],25,world,&requests[2]);
      MPI_Irecv(&feqold[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],15,world,&requests[3]);   
      MPI_Isend(&feqoldn[0][0][1][0],1,passzf,comm->procneigh[2][0],10,world,&requests[4]);
      MPI_Irecv(&feqoldn[0][0][0][0],1,passzf,comm->procneigh[2][0],20,world,&requests[5]);
      MPI_Isend(&feqoldn[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],20,world,&requests[6]);
      MPI_Irecv(&feqoldn[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],10,world,&requests[7]); 
      MPI_Waitall(8,requests,MPI_STATUS_IGNORE);
    }
    parametercalc_full();
  }
}

//==========================================================================
// Compute the lattice Boltzmann equilibrium distribution functions for
// the D3Q15 model.
//==========================================================================
void FixLbFluid::equilibriumdist15(int xstart, int xend, int ystart, int yend, int zstart, int zend) {

  double rho;
  int i, j, k, l, iup, idwn, jup, jdwn, kup, kdwn;
  double Fx_w, Fy_w, Fz_w;

  double total_density(0.0);
  double drhox, drhoy, drhoz, drhoxx, drhoyy, drhozz;
  double Pxx, Pyy, Pzz, Pxy, Pxz, Pyz;
  double grs, p0;
  double dPdrho;

  double S[2][3],std;
  int jj;
 
  double etacov[15],ghostnoise;


  for (i=xstart; i<xend; i++) {
    iup=i+1;
    idwn=i-1;
    for (j=ystart; j<yend; j++) {
      jup=j+1;
      jdwn=j-1;
      for (k=zstart; k<zend; k++) {
	kup=k+1;
	kdwn=k-1;

	rho=density_lb[i][j][k];
	total_density += rho;

	// Derivatives.
	drhox = (density_lb[iup][j][k] - density_lb[idwn][j][k])/2.0;
	drhoxx = (density_lb[iup][j][k] - 2.0*density_lb[i][j][k] + 
		  density_lb[idwn][j][k]);

	drhoy = (density_lb[i][jup][k] - density_lb[i][jdwn][k])/2.0;
	drhoyy = (density_lb[i][jup][k] - 2.0*density_lb[i][j][k] + 
		  density_lb[i][jdwn][k]);

	drhoz = (density_lb[i][j][kup] - density_lb[i][j][kdwn])/2.0;
	drhozz = (density_lb[i][j][kup] - 2.0*density_lb[i][j][k] + 
		  density_lb[i][j][kdwn]);

	// Need one-sided derivatives for the boundary of the domain, if fixed boundary
	// conditions are used.
	if(domain->periodicity[2]==0){
	  if(comm->myloc[2]==0 && k==1){
	    drhoz = (-3.0*density_lb[i][j][k] + 4.0*density_lb[i][j][k+1] - 
		     density_lb[i][j][k+2])/2.0;
	    drhozz = (-density_lb[i][j][k+3] + 4.0*density_lb[i][j][k+2] - 
		      5.0*density_lb[i][j][k+1] + 2.0*rho);
	  }
	  if(comm->myloc[2]==comm->procgrid[2]-1 && k==subNbz-2){
	    drhoz = -(-3.0*density_lb[i][j][k] + 4.0*density_lb[i][j][k-1] - 
		      density_lb[i][j][k-2])/2.0;
	    drhozz = (-density_lb[i][j][k-3] + 4.0*density_lb[i][j][k-2] - 
		      5.0*density_lb[i][j][k-1] + 2.0*rho);
	  }
	}

	grs = drhox*drhox + drhoy*drhoy + drhoz*drhoz;	

	p0 = rho*a_0-kappa_lb*rho*(drhoxx + drhoyy + drhozz);
//                   kappa_lb is the square gradient coeff in the pressure tensor

	dPdrho = a_0; //assuming here that kappa_lb = 0.


	if(typeLB==1){
	  Pxx = p0 + kappa_lb*(drhox*drhox - 0.5*grs)+(tau-0.5)*(1.0/3.0-dPdrho)*
	    (3.0*u_lb[i][j][k][0]*drhox+u_lb[i][j][k][1]*drhoy+u_lb[i][j][k][2]*drhoz);
	  Pyy = p0 + kappa_lb*(drhoy*drhoy - 0.5*grs)+(tau-0.5)*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhox+3.0*u_lb[i][j][k][1]*drhoy+u_lb[i][j][k][2]*drhoz);
	  Pzz = p0 + kappa_lb*(drhoz*drhoz - 0.5*grs)+(tau-0.5)*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhox+u_lb[i][j][k][1]*drhoy+3.0*u_lb[i][j][k][2]*drhoz);
	  Pxy = kappa_lb*drhox*drhoy+(tau-0.5)*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhoy+u_lb[i][j][k][1]*drhox);
	  Pxz = kappa_lb*drhox*drhoz+(tau-0.5)*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhoz+u_lb[i][j][k][2]*drhox);
	  Pyz = kappa_lb*drhoy*drhoz+(tau-0.5)*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][1]*drhoz+u_lb[i][j][k][2]*drhoy);
	}else if(typeLB==2){
	  Pxx = p0 + kappa_lb*(drhox*drhox - 0.5*grs)+tau*(1.0/3.0-dPdrho)*
	    (3.0*u_lb[i][j][k][0]*drhox+u_lb[i][j][k][1]*drhoy+u_lb[i][j][k][2]*drhoz);
	  Pyy = p0 + kappa_lb*(drhoy*drhoy - 0.5*grs)+tau*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhox+3.0*u_lb[i][j][k][1]*drhoy+u_lb[i][j][k][2]*drhoz);
	  Pzz = p0 + kappa_lb*(drhoz*drhoz - 0.5*grs)+tau*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhox+u_lb[i][j][k][1]*drhoy+3.0*u_lb[i][j][k][2]*drhoz);
	  Pxy = kappa_lb*drhox*drhoy+tau*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhoy+u_lb[i][j][k][1]*drhox);
	  Pxz = kappa_lb*drhox*drhoz+tau*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhoz+u_lb[i][j][k][2]*drhox);
	  Pyz = kappa_lb*drhoy*drhoz+tau*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][1]*drhoz+u_lb[i][j][k][2]*drhoy);
	}	  

 	Fx_w = Ff[i][j][k][0];
 	Fy_w = Ff[i][j][k][1];
 	Fz_w = Ff[i][j][k][2];

	etacov[0] = rho;
	etacov[1] = rho*u_lb[i][j][k][0] + Fx_w*tau + rho*bodyforcex*tau;
	etacov[2] = rho*u_lb[i][j][k][1] + Fy_w*tau + rho*bodyforcey*tau;
	etacov[3] = rho*u_lb[i][j][k][2] + Fz_w*tau + rho*bodyforcez*tau;

	etacov[4] = Pxx + rho*u_lb[i][j][k][0]*u_lb[i][j][k][0] -rho/3. + 
	  tau*(2.0*u_lb[i][j][k][0]*(Fx_w+rho*bodyforcex));
	etacov[5] = Pyy + rho*u_lb[i][j][k][1]*u_lb[i][j][k][1] -rho/3. + 
	  tau*(2.0*u_lb[i][j][k][1]*(Fy_w+rho*bodyforcey));
	etacov[6] = Pzz + rho*u_lb[i][j][k][2]*u_lb[i][j][k][2] -rho/3. + 
	  tau*(2.0*u_lb[i][j][k][2]*(Fz_w+rho*bodyforcez));
	etacov[7] = Pxy + rho*u_lb[i][j][k][0]*u_lb[i][j][k][1] + 
	  tau*(u_lb[i][j][k][0]*(Fy_w+rho*bodyforcey) + (Fx_w+rho*bodyforcex)*u_lb[i][j][k][1]);
	etacov[8] = Pyz + rho*u_lb[i][j][k][1]*u_lb[i][j][k][2] + 
	  tau*(u_lb[i][j][k][1]*(Fz_w+rho*bodyforcez) + (Fy_w+rho*bodyforcey)*u_lb[i][j][k][2]);
	etacov[9] = Pxz + rho*u_lb[i][j][k][0]*u_lb[i][j][k][2] + 
	  tau*(u_lb[i][j][k][0]*(Fz_w+rho*bodyforcez) + (Fx_w+rho*bodyforcex)*u_lb[i][j][k][2]);
	etacov[10] = 0.0; 
	etacov[11] = 0.0; 
	etacov[12] = 0.0;
	etacov[13] = rho*u_lb[i][j][k][0]*u_lb[i][j][k][1]*u_lb[i][j][k][2];
	const double TrP = Pxx+Pyy+Pzz;
	etacov[14] = K_0*(rho-TrP);
       
	for (l=0; l<15; l++) {

	  feq[i][j][k][l] = 0.0;
 	  for (int ii=0; ii<15; ii++) 
 	    feq[i][j][k][l] += w_lb[l]*mg_lb[ii][l]*etacov[ii]*Ng_lb[ii];

	  if(typeLB == 2){
	    feqn[i][j][k][l] = feq[i][j][k][l];
	  }
	}

	if(noisestress==1){
	  std = sqrt(namp*rho);

	  for(jj=0; jj<3; jj++)
	    S[0][jj] = std*random->gaussian();
	  for(jj=0; jj<3; jj++)
	    S[1][jj] = std*random->gaussian(); 

	  etacov[4] = (S[0][0]*sqrt(3.0-3.0*a_0));
	  etacov[5] = ((1.0-3.0*a_0)*S[0][0]/sqrt(3.0-3.0*a_0)+
		       sqrt((8.0-12.0*a_0)/(3.0-3.0*a_0))*S[0][1]);
	  etacov[6] = ((1.0-3.0*a_0)*S[0][0]/sqrt(3.0-3.0*a_0)+
		       (2.0-6.0*a_0)*S[0][1]/sqrt((8.0-12.0*a_0)*(3.0-3.0*a_0))+
		       sqrt((5.0-9.0*a_0)/(2.0-3.0*a_0))*S[0][2]);
	  etacov[7] = S[1][0];
	  etacov[8] = S[1][1];
	  etacov[9] = S[1][2];

	  for (l=10; l<15; l++) {
	    etacov[l] = sqrt(9.0*namp*rho/Ng_lb[l])*random->gaussian();
	  }
	  etacov[14] += -K_0*(etacov[4]+etacov[5]+etacov[6]);  //correction from noise to TrP

	  for (l=0; l<15; l++) {
	    ghostnoise = w_lb[l]*
	      (mg_lb[4][l]*etacov[4]*Ng_lb[4] + mg_lb[5][l]*etacov[5]*Ng_lb[5] + 
	       mg_lb[6][l]*etacov[6]*Ng_lb[6] + mg_lb[7][l]*etacov[7]*Ng_lb[7] + 
	       mg_lb[8][l]*etacov[8]*Ng_lb[8] + mg_lb[9][l]*etacov[9]*Ng_lb[9] + 
	       mg_lb[10][l]*etacov[10]*Ng_lb[10] + mg_lb[11][l]*etacov[11]*Ng_lb[11]
	       + mg_lb[12][l]*etacov[12]*Ng_lb[12] + mg_lb[13][l]*etacov[13]*Ng_lb[13]
	       + mg_lb[14][l]*etacov[14]*Ng_lb[14]);
	    feq[i][j][k][l] += ghostnoise*noisefactor;
	  }
	}	
      }
    }
  }
}

//==========================================================================
// Compute the lattice Boltzmann equilibrium distribution functions for
// the D3Q19 model.
//==========================================================================
void FixLbFluid::equilibriumdist19(int xstart, int xend, int ystart, int yend, int zstart, int zend) {

  double rho;
  int i, j, k, l, iup, idwn, jup, jdwn, kup, kdwn;
  double Fx_w, Fy_w, Fz_w;

  double total_density(0.0);
  double drhox, drhoy, drhoz, drhoxx, drhoyy, drhozz;
  double Pxx, Pyy, Pzz, Pxy, Pxz, Pyz;
  double grs, p0;
  double dPdrho;

  double S[2][3],std;
  int jj;
 
  double etacov[19],ghostnoise;

  for (i=xstart; i<xend; i++) {
    iup=i+1;
    idwn=i-1;
    for (j=ystart; j<yend; j++) {
      jup=j+1;
      jdwn=j-1;
      for (k=zstart; k<zend; k++) {
	kup=k+1;
	kdwn=k-1;

	rho=density_lb[i][j][k];
	total_density += rho;

	// Derivatives.
	drhox = (density_lb[iup][j][k] - density_lb[idwn][j][k])/2.0;
	drhoxx = (density_lb[iup][j][k] - 2.0*density_lb[i][j][k] + 
		  density_lb[idwn][j][k]);

	drhoy = (density_lb[i][jup][k] - density_lb[i][jdwn][k])/2.0;
	drhoyy = (density_lb[i][jup][k] - 2.0*density_lb[i][j][k] + 
		  density_lb[i][jdwn][k]);

	drhoz = (density_lb[i][j][kup] - density_lb[i][j][kdwn])/2.0;
	drhozz = (density_lb[i][j][kup] - 2.0*density_lb[i][j][k] + 
		  density_lb[i][j][kdwn]);

	// Need one-sided derivatives for the boundary of the domain, if fixed boundary
	// conditions are used.
	if(domain->periodicity[2]==0){
	  if(comm->myloc[2]==0 && k==1){
	    drhoz = (-3.0*density_lb[i][j][k] + 4.0*density_lb[i][j][k+1] - 
		     density_lb[i][j][k+2])/2.0;
	    drhozz = (-density_lb[i][j][k+3] + 4.0*density_lb[i][j][k+2] - 
		      5.0*density_lb[i][j][k+1] + 2.0*rho);
	  }
	  if(comm->myloc[2]==comm->procgrid[2]-1 && k==subNbz-2){
	    drhoz = -(-3.0*density_lb[i][j][k] + 4.0*density_lb[i][j][k-1] - 
		      density_lb[i][j][k-2])/2.0;
	    drhozz = (-density_lb[i][j][k-3] + 4.0*density_lb[i][j][k-2] - 
		      5.0*density_lb[i][j][k-1] + 2.0*rho);
	  }
	}

	grs = drhox*drhox + drhoy*drhoy + drhoz*drhoz;	

	p0 = rho*a_0-kappa_lb*rho*(drhoxx + drhoyy + drhozz);
//                   kappa_lb is the square gradient coeff in the pressure tensor

	dPdrho = a_0; //assuming here that kappa_lb = 0.


	if(typeLB==1){
	  Pxx = p0 + kappa_lb*(drhox*drhox - 0.5*grs)+(tau-0.5)*(1.0/3.0-dPdrho)*
	    (3.0*u_lb[i][j][k][0]*drhox+u_lb[i][j][k][1]*drhoy+u_lb[i][j][k][2]*drhoz);
	  Pyy = p0 + kappa_lb*(drhoy*drhoy - 0.5*grs)+(tau-0.5)*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhox+3.0*u_lb[i][j][k][1]*drhoy+u_lb[i][j][k][2]*drhoz);
	  Pzz = p0 + kappa_lb*(drhoz*drhoz - 0.5*grs)+(tau-0.5)*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhox+u_lb[i][j][k][1]*drhoy+3.0*u_lb[i][j][k][2]*drhoz);
	  Pxy = kappa_lb*drhox*drhoy+(tau-0.5)*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhoy+u_lb[i][j][k][1]*drhox);
	  Pxz = kappa_lb*drhox*drhoz+(tau-0.5)*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhoz+u_lb[i][j][k][2]*drhox);
	  Pyz = kappa_lb*drhoy*drhoz+(tau-0.5)*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][1]*drhoz+u_lb[i][j][k][2]*drhoy);
	}else if(typeLB==2){
	  Pxx = p0 + kappa_lb*(drhox*drhox - 0.5*grs)+tau*(1.0/3.0-dPdrho)*
	    (3.0*u_lb[i][j][k][0]*drhox+u_lb[i][j][k][1]*drhoy+u_lb[i][j][k][2]*drhoz);
	  Pyy = p0 + kappa_lb*(drhoy*drhoy - 0.5*grs)+tau*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhox+3.0*u_lb[i][j][k][1]*drhoy+u_lb[i][j][k][2]*drhoz);
	  Pzz = p0 + kappa_lb*(drhoz*drhoz - 0.5*grs)+tau*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhox+u_lb[i][j][k][1]*drhoy+3.0*u_lb[i][j][k][2]*drhoz);
	  Pxy = kappa_lb*drhox*drhoy+tau*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhoy+u_lb[i][j][k][1]*drhox);
	  Pxz = kappa_lb*drhox*drhoz+tau*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][0]*drhoz+u_lb[i][j][k][2]*drhox);
	  Pyz = kappa_lb*drhoy*drhoz+tau*(1.0/3.0-dPdrho)*
	    (u_lb[i][j][k][1]*drhoz+u_lb[i][j][k][2]*drhoy);
	}	  

 	Fx_w = Ff[i][j][k][0];
 	Fy_w = Ff[i][j][k][1];
 	Fz_w = Ff[i][j][k][2];

	etacov[0] = rho;
	etacov[1] = rho*u_lb[i][j][k][0] + Fx_w*tau + rho*bodyforcex*tau;
	etacov[2] = rho*u_lb[i][j][k][1] + Fy_w*tau + rho*bodyforcey*tau;
	etacov[3] = rho*u_lb[i][j][k][2] + Fz_w*tau + rho*bodyforcez*tau;

	etacov[4] = Pxx + rho*u_lb[i][j][k][0]*u_lb[i][j][k][0] -rho/3. + 
	  tau*(2.0*u_lb[i][j][k][0]*(Fx_w+rho*bodyforcex));
	etacov[5] = Pyy + rho*u_lb[i][j][k][1]*u_lb[i][j][k][1] -rho/3. + 
	  tau*(2.0*u_lb[i][j][k][1]*(Fy_w+rho*bodyforcey));
	etacov[6] = Pzz + rho*u_lb[i][j][k][2]*u_lb[i][j][k][2] -rho/3. + 
	  tau*(2.0*u_lb[i][j][k][2]*(Fz_w+rho*bodyforcez));
	etacov[7] = Pxy + rho*u_lb[i][j][k][0]*u_lb[i][j][k][1] + 
	  tau*(u_lb[i][j][k][0]*(Fy_w+rho*bodyforcey) + (Fx_w+rho*bodyforcex)*u_lb[i][j][k][1]);
	etacov[8] = Pxz + rho*u_lb[i][j][k][0]*u_lb[i][j][k][2] + 
	  tau*(u_lb[i][j][k][0]*(Fz_w+rho*bodyforcez) + (Fx_w+rho*bodyforcex)*u_lb[i][j][k][2]);
	etacov[9] = Pyz + rho*u_lb[i][j][k][1]*u_lb[i][j][k][2] + 
	  tau*(u_lb[i][j][k][1]*(Fz_w+rho*bodyforcez) + (Fy_w+rho*bodyforcey)*u_lb[i][j][k][2]);
	etacov[10] = 0.0; 
	etacov[11] = 0.0; 
	etacov[12] = 0.0;
	etacov[13] = 0.0;
	etacov[14] = 0.0;
	etacov[15] = 0.0;
	etacov[16] = 0.0;
	etacov[17] = 0.0;
	etacov[18] = 0.0;
	
	for (l=0; l<19; l++) {

	  feq[i][j][k][l] = 0.0;
 	  for (int ii=0; ii<19; ii++) 
 	    feq[i][j][k][l] += w_lb[l]*mg_lb[ii][l]*etacov[ii]*Ng_lb[ii];

	  if(typeLB == 2){
	    feqn[i][j][k][l] = feq[i][j][k][l];
	  }
	}

	if(noisestress==1){
	  std = sqrt(namp*rho);

	  for(jj=0; jj<3; jj++)
	    S[0][jj] = std*random->gaussian();
	  for(jj=0; jj<3; jj++)
	    S[1][jj] = std*random->gaussian(); 

	  etacov[4] = (S[0][0]*sqrt(3.0-3.0*a_0));
	  etacov[5] = ((1.0-3.0*a_0)*S[0][0]/sqrt(3.0-3.0*a_0)+
		       sqrt((8.0-12.0*a_0)/(3.0-3.0*a_0))*S[0][1]);
	  etacov[6] = ((1.0-3.0*a_0)*S[0][0]/sqrt(3.0-3.0*a_0)+
		       (2.0-6.0*a_0)*S[0][1]/sqrt((8.0-12.0*a_0)*(3.0-3.0*a_0))+
		       sqrt((5.0-9.0*a_0)/(2.0-3.0*a_0))*S[0][2]);
	  etacov[7] = S[1][0];
	  etacov[8] = S[1][1];
	  etacov[9] = S[1][2];

	  for (l=10; l<19; l++) {
	    etacov[l] = sqrt(9.0*namp*rho/Ng_lb[l])*random->gaussian();
	  }
	  
	  for (l=0; l<19; l++) {
	    ghostnoise = w_lb[l]*
	      (mg_lb[4][l]*etacov[4]*Ng_lb[4] + mg_lb[5][l]*etacov[5]*Ng_lb[5] + 
	       mg_lb[6][l]*etacov[6]*Ng_lb[6] + mg_lb[7][l]*etacov[7]*Ng_lb[7] + 
	       mg_lb[8][l]*etacov[8]*Ng_lb[8] + mg_lb[9][l]*etacov[9]*Ng_lb[9] + 
	       mg_lb[10][l]*etacov[10]*Ng_lb[10] + mg_lb[11][l]*etacov[11]*Ng_lb[11]
	       + mg_lb[12][l]*etacov[12]*Ng_lb[12] + mg_lb[13][l]*etacov[13]*Ng_lb[13]
	       + mg_lb[14][l]*etacov[14]*Ng_lb[14] + mg_lb[15][l]*etacov[15]*Ng_lb[15]
	       + mg_lb[16][l]*etacov[16]*Ng_lb[16] + mg_lb[17][l]*etacov[17]*Ng_lb[17]
	       + mg_lb[18][l]*etacov[18]*Ng_lb[18]);
	    feq[i][j][k][l] += ghostnoise*noisefactor;
	  }
	}
	
      }
    }
  }
}

//==========================================================================
// Calculate the fluid density and velocity over the entire simulation
// domain.
//==========================================================================
void FixLbFluid::parametercalc_full(void)
{
  MPI_Request requests[4];
  MPI_Request requests2[12];
  int numrequests;
  int i;

  //--------------------------------------------------------------------------
  // send the boundaries of f_lb, as they will be needed later by the update
  // routine, and use these to calculate the density and velocity on the 
  // boundary.
  //--------------------------------------------------------------------------
  for(i=0; i<4; i++)
    requests[i]=MPI_REQUEST_NULL;
  MPI_Isend(&f_lb[1][1][1][0],1,passxf,comm->procneigh[0][0],10,world,&requests[0]);
  MPI_Irecv(&f_lb[0][1][1][0],1,passxf,comm->procneigh[0][0],20,world,&requests[1]);
  MPI_Isend(&f_lb[subNbx-2][1][1][0],1,passxf,comm->procneigh[0][1],20,world,&requests[2]);
  MPI_Irecv(&f_lb[subNbx-1][1][1][0],1,passxf,comm->procneigh[0][1],10,world,&requests[3]);
  parametercalc_part(1,subNbx-1,1,subNby-1,1,subNbz-1);
  MPI_Waitall(4,requests,MPI_STATUS_IGNORE);

  for(i=0; i<4; i++)
    requests[i]=MPI_REQUEST_NULL;
  MPI_Isend(&f_lb[0][1][1][0],1,passyf,comm->procneigh[1][0],10,world,&requests[0]);
  MPI_Irecv(&f_lb[0][0][1][0],1,passyf,comm->procneigh[1][0],20,world,&requests[1]);   
  MPI_Isend(&f_lb[0][subNby-2][1][0],1,passyf,comm->procneigh[1][1],20,world,&requests[2]);
  MPI_Irecv(&f_lb[0][subNby-1][1][0],1,passyf,comm->procneigh[1][1],10,world,&requests[3]);
  parametercalc_part(0,1,1,subNby-1,1,subNbz-1);
  parametercalc_part(subNbx-1,subNbx,1,subNby-1,1,subNbz-1);
  MPI_Waitall(4,requests,MPI_STATUS_IGNORE);

  for(i=0; i<4; i++)
    requests[i]=MPI_REQUEST_NULL;
  MPI_Isend(&f_lb[0][0][1][0],1,passzf,comm->procneigh[2][0],10,world,&requests[0]);
  MPI_Irecv(&f_lb[0][0][0][0],1,passzf,comm->procneigh[2][0],20,world,&requests[1]);
  MPI_Isend(&f_lb[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],20,world,&requests[2]);
  MPI_Irecv(&f_lb[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],10,world,&requests[3]);
  parametercalc_part(0,subNbx,0,1,1,subNbz-1);
  parametercalc_part(0,subNbx,subNby-1,subNby,1,subNbz-1);
  MPI_Waitall(4,requests,MPI_STATUS_IGNORE);
  
  parametercalc_part(0,subNbx,0,subNby,0,1);
  parametercalc_part(0,subNbx,0,subNby,subNbz-1,subNbz);

  //--------------------------------------------------------------------------
  // Send the remaining portions of the u array (and density array if Gamma 
  // is set the default way).
  //--------------------------------------------------------------------------
  if(setGamma == 0) numrequests = 12;
  else numrequests = 6;

  for(i=0; i<numrequests; i++)
    requests2[i]=MPI_REQUEST_NULL;
  MPI_Isend(&u_lb[2][0][0][0],1,passxu,comm->procneigh[0][0],10,world,&requests2[0]);
  MPI_Isend(&u_lb[3][0][0][0],1,passxu,comm->procneigh[0][0],20,world,&requests2[1]);
  MPI_Isend(&u_lb[subNbx-3][0][0][0],1,passxu,comm->procneigh[0][1],30,world,&requests2[2]);
  MPI_Irecv(&u_lb[subNbx][0][0][0],1,passxu,comm->procneigh[0][1],10,world,&requests2[3]);
  MPI_Irecv(&u_lb[subNbx+1][0][0][0],1,passxu,comm->procneigh[0][1],20,world,&requests2[4]);
  MPI_Irecv(&u_lb[subNbx+2][0][0][0],1,passxu,comm->procneigh[0][0],30,world,&requests2[5]);
  if(setGamma==0){
    MPI_Isend(&density_lb[2][0][0],1,passxrho,comm->procneigh[0][0],40,world,&requests2[6]);
    MPI_Isend(&density_lb[3][0][0],1,passxrho,comm->procneigh[0][0],50,world,&requests2[7]);
    MPI_Isend(&density_lb[subNbx-3][0][0],1,passxrho,comm->procneigh[0][1],60,world,&requests2[8]);
    MPI_Irecv(&density_lb[subNbx][0][0],1,passxrho,comm->procneigh[0][1],40,world,&requests2[9]);
    MPI_Irecv(&density_lb[subNbx+1][0][0],1,passxrho,comm->procneigh[0][1],50,world,&requests2[10]);
    MPI_Irecv(&density_lb[subNbx+2][0][0],1,passxrho,comm->procneigh[0][0],60,world,&requests2[11]);
  }
  MPI_Waitall(numrequests,requests2,MPI_STATUS_IGNORE);

  for(i=0; i<numrequests; i++)
    requests2[i]=MPI_REQUEST_NULL;
  MPI_Isend(&u_lb[0][2][0][0],1,passyu,comm->procneigh[1][0],10,world,&requests2[0]);
  MPI_Isend(&u_lb[0][3][0][0],1,passyu,comm->procneigh[1][0],20,world,&requests2[1]);
  MPI_Isend(&u_lb[0][subNby-3][0][0],1,passyu,comm->procneigh[1][1],30,world,&requests2[2]);
  MPI_Irecv(&u_lb[0][subNby][0][0],1,passyu,comm->procneigh[1][1],10,world,&requests2[3]);
  MPI_Irecv(&u_lb[0][subNby+1][0][0],1,passyu,comm->procneigh[1][1],20,world,&requests2[4]);
  MPI_Irecv(&u_lb[0][subNby+2][0][0],1,passyu,comm->procneigh[1][0],30,world,&requests2[5]);
  if(setGamma==0){
    MPI_Isend(&density_lb[0][2][0],1,passyrho,comm->procneigh[1][0],40,world,&requests2[6]);
    MPI_Isend(&density_lb[0][3][0],1,passyrho,comm->procneigh[1][0],50,world,&requests2[7]);
    MPI_Isend(&density_lb[0][subNby-3][0],1,passyrho,comm->procneigh[1][1],60,world,&requests2[8]);
    MPI_Irecv(&density_lb[0][subNby][0],1,passyrho,comm->procneigh[1][1],40,world,&requests2[9]);
    MPI_Irecv(&density_lb[0][subNby+1][0],1,passyrho,comm->procneigh[1][1],50,world,&requests2[10]);
    MPI_Irecv(&density_lb[0][subNby+2][0],1,passyrho,comm->procneigh[1][0],60,world,&requests2[11]);
  }
  MPI_Waitall(numrequests,requests2,MPI_STATUS_IGNORE);

  for(i=0; i<12; i++)
    requests2[i]=MPI_REQUEST_NULL;
  int requestcount=0;
  if(domain->periodicity[2]!=0 || comm->myloc[2] != 0){
    MPI_Isend(&u_lb[0][0][2][0],1,passzu,comm->procneigh[2][0],10,world,&requests2[requestcount]);
    MPI_Isend(&u_lb[0][0][3][0],1,passzu,comm->procneigh[2][0],20,world,&requests2[requestcount+1]);
    MPI_Irecv(&u_lb[0][0][subNbz+2][0],1,passzu,comm->procneigh[2][0],30,world,&requests2[requestcount+2]);
    requestcount=requestcount+3;
    if(setGamma==0){
      MPI_Isend(&density_lb[0][0][2],1,passzrho,comm->procneigh[2][0],40,world,&requests2[requestcount]);
      MPI_Isend(&density_lb[0][0][3],1,passzrho,comm->procneigh[2][0],50,world,&requests2[requestcount+1]);
      MPI_Irecv(&density_lb[0][0][subNbz+2],1,passzrho,comm->procneigh[2][0],60,world,&requests2[requestcount+2]);
      requestcount=requestcount+3;
    }
  }
  if(domain->periodicity[2]!=0 || comm->myloc[2] != (comm->procgrid[2]-1)){
    MPI_Isend(&u_lb[0][0][subNbz-3][0],1,passzu,comm->procneigh[2][1],30,world,&requests2[requestcount]);
    MPI_Irecv(&u_lb[0][0][subNbz][0],1,passzu,comm->procneigh[2][1],10,world,&requests2[requestcount+1]);
    MPI_Irecv(&u_lb[0][0][subNbz+1][0],1,passzu,comm->procneigh[2][1],20,world,&requests2[requestcount+2]);
    requestcount=requestcount+3;
    if(setGamma==0){
      MPI_Isend(&density_lb[0][0][subNbz-3],1,passzrho,comm->procneigh[2][1],60,world,&requests2[requestcount]);
      MPI_Irecv(&density_lb[0][0][subNbz],1,passzrho,comm->procneigh[2][1],40,world,&requests2[requestcount+1]);
      MPI_Irecv(&density_lb[0][0][subNbz+1],1,passzrho,comm->procneigh[2][1],50,world,&requests2[requestcount+2]);
      requestcount=requestcount+3;
    }
  }    
  MPI_Waitall(requestcount,requests2,MPI_STATUS_IGNORE); 

}

//==========================================================================
// Calculate the fluid density and velocity over a simulation volume
// specified by xstart,xend; ystart,yend; zstart,zend.
//==========================================================================
void FixLbFluid::parametercalc_part(int xstart, int xend, int ystart, int yend, int zstart, int zend)
{
  int i,j,k,m;

  for(i=xstart; i<xend; i++){
    for(j=ystart; j<yend; j++){
      for(k=zstart; k<zend; k++){

	density_lb[i][j][k]=0.0;
	u_lb[i][j][k][0]=0.0;
	u_lb[i][j][k][1]=0.0;
	u_lb[i][j][k][2]=0.0;
	for (m=0; m<numvel; m++) {
	 
	  density_lb[i][j][k] += f_lb[i][j][k][m];
	  
	  u_lb[i][j][k][0] += f_lb[i][j][k][m]*e[m][0];
	  u_lb[i][j][k][1] += f_lb[i][j][k][m]*e[m][1];
	  u_lb[i][j][k][2] += f_lb[i][j][k][m]*e[m][2];
	  
	}
	
	//For the on-lattice wall scheme, need to set this velocity to zero.
	if(domain->periodicity[2]==0){
	  if(comm->myloc[2]==0){
	    if(k==1){
	      u_lb[i][j][k][2]=0.0;
	    }
	  }
	  if(comm->myloc[2]==comm->procgrid[2]-1){
	    if(k==subNbz-2){
	      u_lb[i][j][k][2]=0.0;
	    }
	  }
	  
	}
      
	u_lb[i][j][k][0]=u_lb[i][j][k][0]/density_lb[i][j][k];
	u_lb[i][j][k][1]=u_lb[i][j][k][1]/density_lb[i][j][k];
	u_lb[i][j][k][2]=u_lb[i][j][k][2]/density_lb[i][j][k];
      }
    }
  }
  
}

//==========================================================================
// Update the distribution function over a simulation volume specified 
// by xstart,xend; ystart,yend; zstart,zend.
//==========================================================================
void FixLbFluid::update_periodic(int xstart, int xend, int ystart, int yend, int zstart, int zend)
{
  int i,j,k,m;
  int imod,jmod,kmod,imodm,jmodm,kmodm;

  for(i=xstart; i<xend; i++)
    for(j=ystart; j<yend; j++)
      for(k=zstart; k<zend; k++){
	
	if(typeLB==1){
	  for(m=0; m<numvel; m++){
	    imod = i-e[m][0];
	    jmod = j-e[m][1];
	    kmod = k-e[m][2];

	    fnew[i][j][k][m] = f_lb[imod][jmod][kmod][m] + (feq[imod][jmod][kmod][m]-f_lb[imod][jmod][kmod][m])/tau;
	  }	    
	}else if(typeLB==2){
	  for(m=0; m<numvel; m++){
	    imod = i-e[m][0];
	    jmod = j-e[m][1];
	    kmod = k-e[m][2];
	    
	    fnew[i][j][k][m] = feq[imod][jmod][kmod][m] + (f_lb[imod][jmod][kmod][m] - feq[imod][jmod][kmod][m])*expminusdtovertau;
	  }
	  
	  fnew[i][j][k][0]+=Dcoeff*(feq[i][j][k][0]-feqold[i][j][k][0]);
	  for(m=1; m<numvel; m++){
	    imod = i-e[m][0];
	    jmod = j-e[m][1];
	    kmod = k-e[m][2];
	    imodm = i+e[m][0];
	    jmodm = j+e[m][1];
	    kmodm = k+e[m][2];
	    
	     fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) + (0.5-Dcoeff*(tau+0.5))*
	       (feqn[imodm][jmodm][kmodm][m] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);

	  }
	}		   
      }  
}

//==========================================================================
//   Print the fluid properties to the screen.
//==========================================================================
void FixLbFluid::streamout(void)
{
  int i,j,k;
  int istart,jstart,kstart;
  int iend,jend,kend;
  int w,iproc;
  int size,sizeloc;
  MPI_Request request_send,request_recv;
  MPI_Status status;

  //--------------------------------------------------------------------------
  // **Uncomment in order to test conservation of mass and momentum.
  //--------------------------------------------------------------------------
  // massloc=0.0;
  // momentumloc[0]=momentumloc[1]=momentumloc[2]=0.0;
  // for(i=1; i<subNbx-1; i++){
  //   for(j=1; j<subNby-1; j++){
  //     for(k=1; k<subNbz-1; k++){
  // 	massloc += density_lb[i][j][k];
  // 	momentumloc[0] += density_lb[i][j][k]*u_lb[i][j][k][0];
  // 	momentumloc[1] += density_lb[i][j][k]*u_lb[i][j][k][1];
  // 	momentumloc[2] += density_lb[i][j][k]*u_lb[i][j][k][2];
  //     }
  //   }
  // }

  // MPI_Allreduce(&massloc,&mass,1,MPI_DOUBLE,MPI_SUM,world);
  // MPI_Allreduce(&momentumloc[0],&momentum[0],3,MPI_DOUBLE,MPI_SUM,world);

  // if(comm->me==0){
  //   printf("%16.12f %16.12f %16.12f %16.12f\n",mass*dm_lb,momentum[0]*dm_lb*dx_lb/dt_lb,momentum[1]*dm_lb*dx_lb/dt_lb,momentum[2]*dm_lb*dx_lb/dt_lb);
  //  }

  sizeloc=(subNbx*subNby*subNbz*4);
  MPI_Allreduce(&sizeloc,&size,1,MPI_INT,MPI_MAX,world);

  if(me==0){
    for(iproc=0; iproc < comm->nprocs; iproc++){
      if(iproc){
	MPI_Irecv(&buf[0][0][0][0],size,MPI_DOUBLE,iproc,0,world,&request_recv);
	MPI_Wait(&request_recv,&status);

	istart=static_cast<int> (buf[0][0][0][0]);
	jstart=static_cast<int> (buf[0][0][0][1]);
	kstart=static_cast<int> (buf[0][0][0][2]);
	iend=static_cast<int> (buf[0][0][1][0]);
	jend=static_cast<int> (buf[0][0][1][1]);
	kend=static_cast<int> (buf[0][0][1][2]);

	for(i=istart; i<iend; i++){
	  for(j=jstart; j<jend; j++){
	    for(k=kstart; k<kend; k++){
	      for(w=0; w<4; w++){
		altogether[i][j][k][w]=buf[i-istart+1][j-jstart+1][k-kstart+1][w];
	      }
	    }
	  }
	}
      }else{
	for(i=1; i<subNbx-1; i++){
	  for(j=1; j<subNby-1; j++){
	    for(k=1; k<subNbz-1; k++){
	      altogether[i-1][j-1][k-1][0]=density_lb[i][j][k];
	      altogether[i-1][j-1][k-1][1]=u_lb[i][j][k][0];
	      altogether[i-1][j-1][k-1][2]=u_lb[i][j][k][1];
	      altogether[i-1][j-1][k-1][3]=u_lb[i][j][k][2];
	    }
	  }
	}
      }
    }
    //i = Nbx/2;
    //j = Nby/2;
    for(i=0; i<Nbx; i++)
      for(j=0; j<Nby; j++)
	for(k=0; k<Nbz; k++){
	  printf("%16.12f %16.12f %16.12f %16.12f\n",altogether[i][j][k][0]*dm_lb/dx_lb/dx_lb/dx_lb,altogether[i][j][k][1]*dx_lb/dt_lb,altogether[i][j][k][2]*dx_lb/dt_lb,altogether[i][j][k][3]*dx_lb/dt_lb);
	}
    
    
  } else {
    istart=comm->myloc[0]*(subNbx-2);
    jstart=comm->myloc[1]*(subNby-2);
    if(domain->periodicity[2]==0){
      if(comm->myloc[2]==comm->procgrid[2]-1){
	kstart=comm->myloc[2]*(subNbz-3);
      }else{
	kstart=comm->myloc[2]*(subNbz-2);
      }
    }else{
      kstart=comm->myloc[2]*(subNbz-2);
    }
    iend=istart+subNbx-2;
    jend=jstart+subNby-2;
    kend=kstart+subNbz-2;
    for(i=0; i<subNbx; i++){
      for(j=0; j<subNby; j++){
	for(k=0; k<subNbz; k++){
	  buf[i][j][k][0]=density_lb[i][j][k];
	  buf[i][j][k][1]=u_lb[i][j][k][0];
	  buf[i][j][k][2]=u_lb[i][j][k][1];
	  buf[i][j][k][3]=u_lb[i][j][k][2];
	}
      }
    }
    buf[0][0][0][0]=istart;
    buf[0][0][0][1]=jstart;
    buf[0][0][0][2]=kstart;
    buf[0][0][1][0]=iend;
    buf[0][0][1][1]=jend;
    buf[0][0][1][2]=kend;

    MPI_Isend(&buf[0][0][0][0],size,MPI_DOUBLE,0,0,world,&request_send);
    MPI_Wait(&request_send,&status); 
  }

}

//==========================================================================
// Update the distribution functions over the entire simulation domain for
// the D3Q15 model.
//==========================================================================
void FixLbFluid::update_full15(void)
{
  
   MPI_Request req_send15,req_recv15;
   MPI_Request req_send25,req_recv25;
   MPI_Request requests[8];
   int numrequests;
   double tmp1;
   MPI_Status status;
   double rb;
   int i,j,k,m;
   int imod,jmod,kmod;
   int imodm,jmodm,kmodm;

   //--------------------------------------------------------------------------
   // If using the standard LB integrator, do not need to send info about feqn.
   //--------------------------------------------------------------------------
   if(typeLB == 1){
     numrequests = 4;
   }else{
     numrequests = 8;
   }
 
   //--------------------------------------------------------------------------
   // Fixed z boundary conditions.
   //--------------------------------------------------------------------------
   if(domain->periodicity[2]==0){  
 
     for(i=0; i<numrequests; i++)
       requests[i]=MPI_REQUEST_NULL;
     MPI_Isend(&feq[1][1][1][0],1,passxf,comm->procneigh[0][0],15,world,&requests[0]);
     MPI_Irecv(&feq[0][1][1][0],1,passxf,comm->procneigh[0][0],25,world,&requests[1]);
     MPI_Isend(&feq[subNbx-2][1][1][0],1,passxf,comm->procneigh[0][1],25,world,&requests[2]);
     MPI_Irecv(&feq[subNbx-1][1][1][0],1,passxf,comm->procneigh[0][1],15,world,&requests[3]);
     if(typeLB == 2){
       MPI_Isend(&feqn[1][1][1][0],1,passxf,comm->procneigh[0][0],10,world,&requests[4]);
       MPI_Irecv(&feqn[0][1][1][0],1,passxf,comm->procneigh[0][0],20,world,&requests[5]);
       MPI_Isend(&feqn[subNbx-2][1][1][0],1,passxf,comm->procneigh[0][1],20,world,&requests[6]);
       MPI_Irecv(&feqn[subNbx-1][1][1][0],1,passxf,comm->procneigh[0][1],10,world,&requests[7]);
     }
     update_periodic(2,subNbx-2,2,subNby-2,2,subNbz-2);
     MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);


     for(i=0; i<numrequests; i++)
       requests[i]=MPI_REQUEST_NULL;
     MPI_Isend(&feq[0][1][1][0],1,passyf,comm->procneigh[1][0],15,world,&requests[0]);
     MPI_Irecv(&feq[0][0][1][0],1,passyf,comm->procneigh[1][0],25,world,&requests[1]);   
     MPI_Isend(&feq[0][subNby-2][1][0],1,passyf,comm->procneigh[1][1],25,world,&requests[2]);
     MPI_Irecv(&feq[0][subNby-1][1][0],1,passyf,comm->procneigh[1][1],15,world,&requests[3]);
     if(typeLB == 2){
       MPI_Isend(&feqn[0][1][1][0],1,passyf,comm->procneigh[1][0],10,world,&requests[4]);
       MPI_Irecv(&feqn[0][0][1][0],1,passyf,comm->procneigh[1][0],20,world,&requests[5]);   
       MPI_Isend(&feqn[0][subNby-2][1][0],1,passyf,comm->procneigh[1][1],20,world,&requests[6]);
       MPI_Irecv(&feqn[0][subNby-1][1][0],1,passyf,comm->procneigh[1][1],10,world,&requests[7]);
     }
     update_periodic(1,2,2,subNby-2,2,subNbz-2);
     update_periodic(subNbx-2,subNbx-1,2,subNby-2,2,subNbz-2);
     MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);
     
     for(i=0; i<numrequests; i++)
       requests[i]=MPI_REQUEST_NULL;
     MPI_Isend(&feq[0][0][1][0],1,passzf,comm->procneigh[2][0],15,world,&requests[0]);
     MPI_Irecv(&feq[0][0][0][0],1,passzf,comm->procneigh[2][0],25,world,&requests[1]);
     MPI_Isend(&feq[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],25,world,&requests[2]);
     MPI_Irecv(&feq[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],15,world,&requests[3]);
     if(typeLB == 2){
       MPI_Isend(&feqn[0][0][1][0],1,passzf,comm->procneigh[2][0],10,world,&requests[4]);
       MPI_Irecv(&feqn[0][0][0][0],1,passzf,comm->procneigh[2][0],20,world,&requests[5]);
       MPI_Isend(&feqn[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],20,world,&requests[6]);
       MPI_Irecv(&feqn[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],10,world,&requests[7]);
     } 
     update_periodic(1,subNbx-1,1,2,2,subNbz-2);
     update_periodic(1,subNbx-1,subNby-2,subNby-1,2,subNbz-2);
     MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);

     if(typeLB==1){
       update_periodic(1,subNbx-1,1,subNby-1,1,2);
       update_periodic(1,subNbx-1,1,subNby-1,subNbz-2,subNbz-1);
     }else if(typeLB==2){
       if(comm->myloc[2]==0){
     	 for(i=1; i<subNbx-1; i++){
     	   for(j=1;j<subNby-1;j++){
     	     k=1;
     	     for(m=0; m<15; m++){
     	       imod = i-e[m][0];
     	       jmod = j-e[m][1];
     	       kmod = k-e[m][2];
	     
     	       fnew[i][j][k][m] = feq[imod][jmod][kmod][m] + (f_lb[imod][jmod][kmod][m]-feq[imod][jmod][kmod][m])*expminusdtovertau;
     	     }

     	     for(m=0; m<15; m++){
     	       imod = i-e[m][0];
     	       jmod = j-e[m][1];
     	       kmod = k-e[m][2];
	       imodm = i+e[m][0];
	       jmodm = j+e[m][1];
	       kmodm = k+e[m][2];
	       
     	       if(m==5)
     	       	 fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][6] - feqold[imod][jmod][kmod][m]) + 
     	       	   (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][6] - feqn[imod][jmod][kmod][6]);
     	       else if(m==7)
     	       	 fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][11] - feqold[imod][jmod][kmod][m]) +
     	       	   (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][11] - feqn[imod][jmod][kmod][11]);
     	       else if(m==8)
     	       	 fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][12] - feqold[imod][jmod][kmod][m]) +
     	       	   (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][12] - feqn[imod][jmod][kmod][12]);
     	       else if(m==9)
     	       	 fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][13] - feqold[imod][jmod][kmod][m]) +
     	       	   (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][13] - feqn[imod][jmod][kmod][13]);
     	       else if(m==10)
     	       	 fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][14] - feqold[imod][jmod][kmod][m]) +
     	       	   (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][14] - feqn[imod][jmod][kmod][14]);
	       else if(m==6)
		 fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m]-feqold[imod][jmod][kmod][m]) +
		   (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][5] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);
	       else if(m==11)
		 fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m]-feqold[imod][jmod][kmod][m]) +
		   (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][7] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);	
	       else if(m==12)
		 fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m]-feqold[imod][jmod][kmod][m]) +
		   (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][8] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);	 
	       else if(m==13)
		 fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m]-feqold[imod][jmod][kmod][m]) +
		   (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][9] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);	 
	       else if(m==14)
		 fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m]-feqold[imod][jmod][kmod][m]) +
		   (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][10] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);	
     	       else
     	       	 fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) + 
     	       	   (0.5-Dcoeff*(tau+0.5))*(feqn[imodm][jmodm][kmodm][m] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);

     	     }
     	   }
     	 }      
       }else{
     	 update_periodic(1,subNbx-1,1,subNby-1,1,2);
       }
       if(comm->myloc[2]==comm->procgrid[2]-1){
     	 for(i=1;i<subNbx-1;i++){
     	   for(j=1;j<subNby-1;j++){
     	     k=subNbz-2;
     	     for(m=0; m<15; m++){
     	       imod = i-e[m][0];
     	       jmod = j-e[m][1];
     	       kmod = k-e[m][2];

     	       fnew[i][j][k][m] = feq[imod][jmod][kmod][m] + (f_lb[imod][jmod][kmod][m]-feq[imod][jmod][kmod][m])*expminusdtovertau;
     	     }	    
     	     for(m=0; m<15; m++){
     	       imod = i-e[m][0];
     	       jmod = j-e[m][1];
     	       kmod = k-e[m][2];
	       imodm = i+e[m][0];
	       jmodm = j+e[m][1];
	       kmodm = k+e[m][2];
	       
     	       if(m==6)
     	       	 fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][5] - feqold[imod][jmod][kmod][m]) + 
     	       	   (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][5] - feqn[imod][jmod][kmod][5]);
     	       else if(m==11)
     	       	 fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][7] - feqold[imod][jmod][kmod][m]) +
     	       	   (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][7] - feqn[imod][jmod][kmod][7]);
     	       else if(m==12)
     	       	 fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][8] - feqold[imod][jmod][kmod][m]) +
     	       	   (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][8] - feqn[imod][jmod][kmod][8]);
     	       else if(m==13)
     	       	 fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][9] - feqold[imod][jmod][kmod][m]) +
     	       	   (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][9] - feqn[imod][jmod][kmod][9]);
     	       else if(m==14)
     	       	 fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][10] - feqold[imod][jmod][kmod][m]) +
     	       	   (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][10] - feqn[imod][jmod][kmod][10]);
	       else if(m==5)
		 fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		   (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][6] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);
	       else if(m==7)
		 fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		   (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][11] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);	   
	       else if(m==8)
		 fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		   (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][12] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);	
 	       else if(m==9)
		 fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		   (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][13] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);	   
	       else if(m==10)
		 fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		   (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][14] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);	
     	       else
     	       	 fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) + 
     	       	   (0.5-Dcoeff*(tau+0.5))*(feqn[imodm][jmodm][kmodm][m] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);
	       
     	     }
     	   }
     	 }
       }
       else{
     	 update_periodic(1,subNbx-1,1,subNby-1,subNbz-2,subNbz-1);
       }
     }
     
     req_send15=MPI_REQUEST_NULL;
     req_recv25=MPI_REQUEST_NULL;
     req_send25=MPI_REQUEST_NULL;
     req_recv15=MPI_REQUEST_NULL;
     
     if(comm->myloc[2]==0){
       MPI_Isend(&fnew[0][0][1][0],1,passzf,comm->procneigh[2][0],15,world,&req_send15);
       MPI_Irecv(&fnew[0][0][0][0],1,passzf,comm->procneigh[2][0],25,world,&req_recv25);
     }
     
     if(comm->myloc[2]==comm->procgrid[2]-1){
       MPI_Isend(&fnew[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],25,world,&req_send25);
       MPI_Irecv(&fnew[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],15,world,&req_recv15);
     }
     if(comm->myloc[2]==0){
       MPI_Wait(&req_send15,&status);
       MPI_Wait(&req_recv25,&status);
       
       for(i=1;i<subNbx-1;i++){
	 for(j=1;j<subNby-1;j++){
	   k=1;
	   if(typeLB == 1){
	     fnew[i][j][k][5]=fnew[i][j][k-1][6];
	     tmp1=fnew[i][j][k-1][11]+fnew[i][j][k-1][12]+fnew[i][j][k-1][13]+fnew[i][j][k-1][14]; 
	   }
	   else{
	     fnew[i][j][k][5]=fnew[i][j][k-1][6] + (0.5-Dcoeff*(tau+0.5))*feqn[i][j][k+1][5];
	     tmp1=fnew[i][j][k-1][11]+fnew[i][j][k-1][12]+fnew[i][j][k-1][13]+fnew[i][j][k-1][14] + 
	       (0.5-Dcoeff*(tau+0.5))*(feqn[i-1][j-1][k+1][7] + feqn[i+1][j-1][k+1][8] + 
				       feqn[i+1][j+1][k+1][9] + feqn[i-1][j+1][k+1][10]);
	   }

	   fnew[i][j][k][7]=-0.25*(fnew[i][j][k][1]+fnew[i][j][k][2]-fnew[i][j][k][3]-
	   			   fnew[i][j][k][4]+2.0*fnew[i][j][k][11]-2.0*fnew[i][j][k][13]-tmp1);
	   fnew[i][j][k][8]=0.25*(fnew[i][j][k][1]-fnew[i][j][k][2]-fnew[i][j][k][3]+
	   			  fnew[i][j][k][4]+2.0*fnew[i][j][k][14]-2.0*fnew[i][j][k][12]+tmp1);
	   fnew[i][j][k][9]=0.25*(fnew[i][j][k][1]+fnew[i][j][k][2]-fnew[i][j][k][3]-
	   			  fnew[i][j][k][4]+2.0*fnew[i][j][k][11]-2.0*fnew[i][j][k][13]+tmp1);
	   fnew[i][j][k][10]=-0.25*(fnew[i][j][k][1]-fnew[i][j][k][2]-fnew[i][j][k][3]+
	   			    fnew[i][j][k][4]+2.0*fnew[i][j][k][14]-2.0*fnew[i][j][k][12]-tmp1);



	   rb=fnew[i][j][k][0]+fnew[i][j][k][1]+fnew[i][j][k][2]+fnew[i][j][k][3]+fnew[i][j][k][4]+
	     fnew[i][j][k][5]+fnew[i][j][k][6]+tmp1+fnew[i][j][k][11]+fnew[i][j][k][12]+
	     fnew[i][j][k][13]+fnew[i][j][k][14];
	   
	   fnew[i][j][k][7] += 0.25*rb*vwbt;
	   fnew[i][j][k][8] += 0.25*rb*vwbt;
	   fnew[i][j][k][9] += -0.25*rb*vwbt;
	   fnew[i][j][k][10] += -0.25*rb*vwbt;
	 }
       }

     }
     if(comm->myloc[2]==comm->procgrid[2]-1){
       MPI_Wait(&req_send25,&status);
       MPI_Wait(&req_recv15,&status);
       
       for(i=1;i<subNbx-1;i++){
	 for(j=1;j<subNby-1;j++){
	   k=subNbz-2;
	   
	   if(typeLB == 1){
	     fnew[i][j][k][6]=fnew[i][j][k+1][5];
	     tmp1=fnew[i][j][k+1][7]+fnew[i][j][k+1][8]+fnew[i][j][k+1][9]+fnew[i][j][k+1][10];
	   }
	   else{
	     fnew[i][j][k][6]=fnew[i][j][k+1][5] + (0.5-Dcoeff*(tau+0.5))*feqn[i][j][k-1][6];
	     tmp1=fnew[i][j][k+1][7]+fnew[i][j][k+1][8]+fnew[i][j][k+1][9]+fnew[i][j][k+1][10] + 
	       (0.5-Dcoeff*(tau+0.5))*(feqn[i-1][j-1][k-1][11] + feqn[i+1][j-1][k-1][12] + 
				       feqn[i+1][j+1][k-1][13] + feqn[i-1][j+1][k-1][14]);
	   }
  
	   fnew[i][j][k][11]=-0.25*(fnew[i][j][k][1]+fnew[i][j][k][2]-fnew[i][j][k][3]-
	   			    fnew[i][j][k][4]+2.0*fnew[i][j][k][7]-2.0*fnew[i][j][k][9]-tmp1);
	   fnew[i][j][k][12]=0.25*(fnew[i][j][k][1]-fnew[i][j][k][2]-fnew[i][j][k][3]+
	   			   fnew[i][j][k][4]-2.0*fnew[i][j][k][8]+2.0*fnew[i][j][k][10]+tmp1);
	   fnew[i][j][k][13]=0.25*(fnew[i][j][k][1]+fnew[i][j][k][2]-fnew[i][j][k][3]-
	   			   fnew[i][j][k][4]+2.0*fnew[i][j][k][7]-2.0*fnew[i][j][k][9]+tmp1);
	   fnew[i][j][k][14]=-0.25*(fnew[i][j][k][1]-fnew[i][j][k][2]-fnew[i][j][k][3]+
	   			    fnew[i][j][k][4]-2.0*fnew[i][j][k][8]+2.0*fnew[i][j][k][10]-tmp1);

	   
	   rb=fnew[i][j][k][0]+fnew[i][j][k][1]+fnew[i][j][k][2]+fnew[i][j][k][3]+fnew[i][j][k][4]+
	     fnew[i][j][k][5]+fnew[i][j][k][6]+fnew[i][j][k][7]+fnew[i][j][k][8]+fnew[i][j][k][9]+
	     fnew[i][j][k][10]+tmp1;
	   
	   fnew[i][j][k][11] += 0.25*rb*vwtp;
	   fnew[i][j][k][12] += 0.25*rb*vwtp;
	   fnew[i][j][k][13] += -0.25*rb*vwtp;
	   fnew[i][j][k][14] += -0.25*rb*vwtp;
	 }
       }	
     }
     
     //--------------------------------------------------------------------------
     // Periodic z boundary conditions.
     //--------------------------------------------------------------------------
   }else {

     for(i=0; i<numrequests; i++)
       requests[i]=MPI_REQUEST_NULL;
     MPI_Isend(&feq[1][1][1][0],1,passxf,comm->procneigh[0][0],15,world,&requests[0]);
     MPI_Irecv(&feq[0][1][1][0],1,passxf,comm->procneigh[0][0],25,world,&requests[1]);
     MPI_Isend(&feq[subNbx-2][1][1][0],1,passxf,comm->procneigh[0][1],25,world,&requests[2]);
     MPI_Irecv(&feq[subNbx-1][1][1][0],1,passxf,comm->procneigh[0][1],15,world,&requests[3]);
     if(typeLB == 2){
       MPI_Isend(&feqn[1][1][1][0],1,passxf,comm->procneigh[0][0],10,world,&requests[4]);
       MPI_Irecv(&feqn[0][1][1][0],1,passxf,comm->procneigh[0][0],20,world,&requests[5]);
       MPI_Isend(&feqn[subNbx-2][1][1][0],1,passxf,comm->procneigh[0][1],20,world,&requests[6]);
       MPI_Irecv(&feqn[subNbx-1][1][1][0],1,passxf,comm->procneigh[0][1],10,world,&requests[7]);
     }
     update_periodic(2,subNbx-2,2,subNby-2,2,subNbz-2);
     MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);

     for(i=0; i<numrequests; i++)
       requests[i]=MPI_REQUEST_NULL;
     MPI_Isend(&feq[0][1][1][0],1,passyf,comm->procneigh[1][0],15,world,&requests[0]);
     MPI_Irecv(&feq[0][0][1][0],1,passyf,comm->procneigh[1][0],25,world,&requests[1]);   
     MPI_Isend(&feq[0][subNby-2][1][0],1,passyf,comm->procneigh[1][1],25,world,&requests[2]);
     MPI_Irecv(&feq[0][subNby-1][1][0],1,passyf,comm->procneigh[1][1],15,world,&requests[3]);
     if(typeLB == 2){
       MPI_Isend(&feqn[0][1][1][0],1,passyf,comm->procneigh[1][0],10,world,&requests[4]);
       MPI_Irecv(&feqn[0][0][1][0],1,passyf,comm->procneigh[1][0],20,world,&requests[5]);   
       MPI_Isend(&feqn[0][subNby-2][1][0],1,passyf,comm->procneigh[1][1],20,world,&requests[6]);
       MPI_Irecv(&feqn[0][subNby-1][1][0],1,passyf,comm->procneigh[1][1],10,world,&requests[7]);
     }
     update_periodic(1,2,2,subNby-2,2,subNbz-2);
     update_periodic(subNbx-2,subNbx-1,2,subNby-2,2,subNbz-2);
     MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);

     for(i=0; i<numrequests; i++)
       requests[i]=MPI_REQUEST_NULL;
     MPI_Isend(&feq[0][0][1][0],1,passzf,comm->procneigh[2][0],15,world,&requests[0]);
     MPI_Irecv(&feq[0][0][0][0],1,passzf,comm->procneigh[2][0],25,world,&requests[1]);
     MPI_Isend(&feq[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],25,world,&requests[2]);
     MPI_Irecv(&feq[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],15,world,&requests[3]);
     if(typeLB == 2){
       MPI_Isend(&feqn[0][0][1][0],1,passzf,comm->procneigh[2][0],10,world,&requests[4]);
       MPI_Irecv(&feqn[0][0][0][0],1,passzf,comm->procneigh[2][0],20,world,&requests[5]);
       MPI_Isend(&feqn[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],20,world,&requests[6]);
       MPI_Irecv(&feqn[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],10,world,&requests[7]);
     }  
     update_periodic(1,subNbx-1,1,2,2,subNbz-2);
     update_periodic(1,subNbx-1,subNby-2,subNby-1,2,subNbz-2);
     MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);
     
     update_periodic(1,subNbx-1,1,subNby-1,1,2);
     update_periodic(1,subNbx-1,1,subNby-1,subNbz-2,subNbz-1);
   }
 
}

//==========================================================================
// Update the distribution functions over the entire simulation domain for
// the D3Q19 model.
//==========================================================================
void FixLbFluid::update_full19(void)
{
  
  MPI_Request req_send15,req_recv15;
  MPI_Request req_send25,req_recv25;
  MPI_Request requests[8];
  int numrequests;
  double tmp1,tmp2,tmp3;
  MPI_Status status;
  double rb;
  int i,j,k,m;
  int imod,jmod,kmod;
  int imodm,jmodm,kmodm;
  
  //--------------------------------------------------------------------------
  // If using the standard LB integrator, do not need to send info about feqn.
  //--------------------------------------------------------------------------
  if(typeLB == 1){
    numrequests = 4;
  }else{
    numrequests = 8;
  }
  
  //--------------------------------------------------------------------------
  // Fixed z boundary conditions.
  //--------------------------------------------------------------------------
  if(domain->periodicity[2]==0){  
    
    for(i=0; i<numrequests; i++)
      requests[i]=MPI_REQUEST_NULL;
    MPI_Isend(&feq[1][1][1][0],1,passxf,comm->procneigh[0][0],15,world,&requests[0]);
    MPI_Irecv(&feq[0][1][1][0],1,passxf,comm->procneigh[0][0],25,world,&requests[1]);
    MPI_Isend(&feq[subNbx-2][1][1][0],1,passxf,comm->procneigh[0][1],25,world,&requests[2]);
    MPI_Irecv(&feq[subNbx-1][1][1][0],1,passxf,comm->procneigh[0][1],15,world,&requests[3]);
    if(typeLB == 2){
      MPI_Isend(&feqn[1][1][1][0],1,passxf,comm->procneigh[0][0],10,world,&requests[4]);
      MPI_Irecv(&feqn[0][1][1][0],1,passxf,comm->procneigh[0][0],20,world,&requests[5]);
      MPI_Isend(&feqn[subNbx-2][1][1][0],1,passxf,comm->procneigh[0][1],20,world,&requests[6]);
      MPI_Irecv(&feqn[subNbx-1][1][1][0],1,passxf,comm->procneigh[0][1],10,world,&requests[7]);
    }
    update_periodic(2,subNbx-2,2,subNby-2,2,subNbz-2);
    MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);
    
    for(i=0; i<numrequests; i++)
      requests[i]=MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][1][1][0],1,passyf,comm->procneigh[1][0],15,world,&requests[0]);
    MPI_Irecv(&feq[0][0][1][0],1,passyf,comm->procneigh[1][0],25,world,&requests[1]);   
    MPI_Isend(&feq[0][subNby-2][1][0],1,passyf,comm->procneigh[1][1],25,world,&requests[2]);
    MPI_Irecv(&feq[0][subNby-1][1][0],1,passyf,comm->procneigh[1][1],15,world,&requests[3]);
    if(typeLB == 2){
      MPI_Isend(&feqn[0][1][1][0],1,passyf,comm->procneigh[1][0],10,world,&requests[4]);
      MPI_Irecv(&feqn[0][0][1][0],1,passyf,comm->procneigh[1][0],20,world,&requests[5]);   
      MPI_Isend(&feqn[0][subNby-2][1][0],1,passyf,comm->procneigh[1][1],20,world,&requests[6]);
      MPI_Irecv(&feqn[0][subNby-1][1][0],1,passyf,comm->procneigh[1][1],10,world,&requests[7]);
    }
    update_periodic(1,2,2,subNby-2,2,subNbz-2);
    update_periodic(subNbx-2,subNbx-1,2,subNby-2,2,subNbz-2);
    MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);
    
    for(i=0; i<numrequests; i++)
      requests[i]=MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][0][1][0],1,passzf,comm->procneigh[2][0],15,world,&requests[0]);
    MPI_Irecv(&feq[0][0][0][0],1,passzf,comm->procneigh[2][0],25,world,&requests[1]);
    MPI_Isend(&feq[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],25,world,&requests[2]);
    MPI_Irecv(&feq[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],15,world,&requests[3]);
    if(typeLB == 2){
      MPI_Isend(&feqn[0][0][1][0],1,passzf,comm->procneigh[2][0],10,world,&requests[4]);
      MPI_Irecv(&feqn[0][0][0][0],1,passzf,comm->procneigh[2][0],20,world,&requests[5]);
      MPI_Isend(&feqn[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],20,world,&requests[6]);
      MPI_Irecv(&feqn[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],10,world,&requests[7]);
    } 
    update_periodic(1,subNbx-1,1,2,2,subNbz-2);
    update_periodic(1,subNbx-1,subNby-2,subNby-1,2,subNbz-2);
    MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);
    
    if(typeLB==1){
      update_periodic(1,subNbx-1,1,subNby-1,1,2);
      update_periodic(1,subNbx-1,1,subNby-1,subNbz-2,subNbz-1);
    }else if(typeLB==2){
      if(comm->myloc[2]==0){
	for(i=1; i<subNbx-1; i++){
	  for(j=1;j<subNby-1;j++){
	    k=1;
	    for(m=0; m<19; m++){
	      imod = i-e[m][0];
	      jmod = j-e[m][1];
	      kmod = k-e[m][2];
	      
	      fnew[i][j][k][m] = feq[imod][jmod][kmod][m] + (f_lb[imod][jmod][kmod][m]-feq[imod][jmod][kmod][m])*expminusdtovertau;
	    }
	    
	    for(m=0; m<19; m++){
	      imod = i-e[m][0];
	      jmod = j-e[m][1];
	      kmod = k-e[m][2];
	      imodm = i+e[m][0];
	      jmodm = j+e[m][1];
	      kmodm = k+e[m][2];
	      
	      if(m==5)
		fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][6] - feqold[imod][jmod][kmod][m]) + 
		  (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][6] - feqn[imod][jmod][kmod][6]);
	      else if(m==11)
		fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][12] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][12] - feqn[imod][jmod][kmod][12]);
	      else if(m==13)
		fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][14] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][14] - feqn[imod][jmod][kmod][14]);
	      else if(m==15)
		fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][16] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][16] - feqn[imod][jmod][kmod][16]);
	      else if(m==17)
		fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][18] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][18] - feqn[imod][jmod][kmod][18]);
	      else if(m==6)
		fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][5] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);
	      else if(m==12)
		fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][11] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);
	      else if(m==14)
		fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][13] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);
	      else if(m==16)
		fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][15] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);
	      else if(m==18)
		fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][17] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);
	      else
		fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) + 
		  (0.5-Dcoeff*(tau+0.5))*(feqn[imodm][jmodm][kmodm][m] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);
	    }
	  }
	}      
      }else{
	update_periodic(1,subNbx-1,1,subNby-1,1,2);
      }
      if(comm->myloc[2]==comm->procgrid[2]-1){
	for(i=1;i<subNbx-1;i++){
	  for(j=1;j<subNby-1;j++){
	    k=subNbz-2;
	    for(m=0; m<19; m++){
	      imod = i-e[m][0];
	      jmod = j-e[m][1];
	      kmod = k-e[m][2];
	      
	      fnew[i][j][k][m] = feq[imod][jmod][kmod][m] + (f_lb[imod][jmod][kmod][m]-feq[imod][jmod][kmod][m])*expminusdtovertau;
	    }	    
	    for(m=0; m<19; m++){
	      imod = i-e[m][0];
	      jmod = j-e[m][1];
	      kmod = k-e[m][2];
	      imodm = i+e[m][0];
	      jmodm = j+e[m][1];
	      kmodm = k+e[m][2];
	      
	      if(m==6)
		fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][5] - feqold[imod][jmod][kmod][m]) + 
		  (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][5] - feqn[imod][jmod][kmod][5]);
	      else if(m==12)
		fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][11] - feqold[imod][jmod][kmod][m]) + 
		  (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][11] - feqn[imod][jmod][kmod][11]);
	      else if(m==14)
		fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][13] - feqold[imod][jmod][kmod][m]) + 
		  (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][13] - feqn[imod][jmod][kmod][13]);
	      else if(m==16)
		fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][15] - feqold[imod][jmod][kmod][m]) + 
		  (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][15] - feqn[imod][jmod][kmod][15]);
	      else if(m==18)
		fnew[i][j][k][m] += Dcoeff*(feq[imod][jmod][kmod][17] - feqold[imod][jmod][kmod][m]) + 
		  (0.5-Dcoeff*(tau+0.5))*(feqoldn[imod][jmod][kmod][m] - feqoldn[imod][jmod][kmod][17] - feqn[imod][jmod][kmod][17]);
	      else if(m==5)
		fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][6] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);	  
	      else if(m==11)
		fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][12] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);	  
	      else if(m==13)
		fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][14] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);	
	      else if(m==15)
		fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][16] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);
	      else if(m==17)
		fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) +
		  (0.5-Dcoeff*(tau+0.5))*(feqn[i][j][k][18] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);	  	  
	      else
		fnew[i][j][k][m] += Dcoeff*(feq[i][j][k][m] - feqold[imod][jmod][kmod][m]) + 
		  (0.5-Dcoeff*(tau+0.5))*(feqn[imodm][jmodm][kmodm][m] - feqoldn[i][j][k][m] - feqn[i][j][k][m] + feqoldn[imod][jmod][kmod][m]);	       
	    }
	  }
	}
      }
      else{
	update_periodic(1,subNbx-1,1,subNby-1,subNbz-2,subNbz-1);
      }
    }
    
    req_send15=MPI_REQUEST_NULL;
    req_recv25=MPI_REQUEST_NULL;
    req_send25=MPI_REQUEST_NULL;
    req_recv15=MPI_REQUEST_NULL;
    
    if(comm->myloc[2]==0){
      MPI_Isend(&fnew[0][0][1][0],1,passzf,comm->procneigh[2][0],15,world,&req_send15);
      MPI_Irecv(&fnew[0][0][0][0],1,passzf,comm->procneigh[2][0],25,world,&req_recv25);
    }
    
    if(comm->myloc[2]==comm->procgrid[2]-1){
      MPI_Isend(&fnew[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],25,world,&req_send25);
      MPI_Irecv(&fnew[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],15,world,&req_recv15);
    }
    if(comm->myloc[2]==0){
      MPI_Wait(&req_send15,&status);
      MPI_Wait(&req_recv25,&status);
      
      for(i=1;i<subNbx-1;i++){
	for(j=1;j<subNby-1;j++){
	  k=1;
	  
	  if(typeLB == 1){
	    fnew[i][j][k][5]=fnew[i][j][k-1][6];
	    tmp1=fnew[i][j][k-1][12]+fnew[i][j][k-1][14]+fnew[i][j][k-1][16]+fnew[i][j][k-1][18];
	  }
	  else{
	    fnew[i][j][k][5]=fnew[i][j][k-1][6] + (0.5-Dcoeff*(tau+0.5))*feqn[i][j][k+1][5];
	    tmp1=fnew[i][j][k-1][12]+fnew[i][j][k-1][14]+fnew[i][j][k-1][16]+fnew[i][j][k-1][18] +
	      (0.5-Dcoeff*(tau+0.5))*(feqn[i-1][j][k+1][11] + feqn[i+1][j][k+1][13] + 
				      feqn[i][j-1][k+1][15] + feqn[i][j+1][k+1][17]);	     
	  }
	  
	  tmp2=fnew[i][j][k][3]+fnew[i][j][k][9]+fnew[i][j][k][10]+fnew[i][j][k][14]-
	    fnew[i][j][k][1]-fnew[i][j][k][7]-fnew[i][j][k][8]-fnew[i][j][k][12];
	  
	  rb=fnew[i][j][k][0]+fnew[i][j][k][1]+fnew[i][j][k][2]+fnew[i][j][k][3]+fnew[i][j][k][4]+
	    fnew[i][j][k][5]+fnew[i][j][k][6]+fnew[i][j][k][7]+fnew[i][j][k][8]+fnew[i][j][k][9]+
	    fnew[i][j][k][10]+fnew[i][j][k][12]+fnew[i][j][k][14]+fnew[i][j][k][16]+fnew[i][j][k][18]+tmp1;
	  
	  tmp3=rb*vwbt-fnew[i][j][k][2]+fnew[i][j][k][4]-fnew[i][j][k][7]+fnew[i][j][k][8]-fnew[i][j][k][9]+
	    fnew[i][j][k][10]-fnew[i][j][k][16]+fnew[i][j][k][18];
	  
	  fnew[i][j][k][11] = 0.25*(tmp1+2.0*tmp2);
	  fnew[i][j][k][13] = 0.25*(tmp1-2.0*tmp2);
	  fnew[i][j][k][15] = 0.25*(tmp1+2.0*tmp3);
	  fnew[i][j][k][17] = 0.25*(tmp1-2.0*tmp3);
	}
      }
      
    }
    if(comm->myloc[2]==comm->procgrid[2]-1){
      MPI_Wait(&req_send25,&status);
      MPI_Wait(&req_recv15,&status);
      
      for(i=1;i<subNbx-1;i++){
	for(j=1;j<subNby-1;j++){
	  k=subNbz-2;
	  
	  if(typeLB == 1){
	    fnew[i][j][k][6]=fnew[i][j][k+1][5];
	    tmp1=fnew[i][j][k+1][11]+fnew[i][j][k+1][13]+fnew[i][j][k+1][15]+fnew[i][j][k+1][17];
	  }
	  else{
	    fnew[i][j][k][6]=fnew[i][j][k+1][5] + (0.5-Dcoeff*(tau+0.5))*feqn[i][j][k-1][6];
	    tmp1=fnew[i][j][k+1][11]+fnew[i][j][k+1][13]+fnew[i][j][k+1][15]+fnew[i][j][k+1][17] +
	      (0.5-Dcoeff*(tau+0.5))*(feqn[i-1][j][k-1][12] + feqn[i+1][j][k-1][14] + 
				      feqn[i][j-1][k-1][16] + feqn[i][j+1][k-1][18]);
	  }
	  
	  tmp2=fnew[i][j][k][3]+fnew[i][j][k][9]+fnew[i][j][k][10]+fnew[i][j][k][13]-fnew[i][j][k][1]-
	    fnew[i][j][k][7]-fnew[i][j][k][8]-fnew[i][j][k][11];
	  
	  rb=fnew[i][j][k][0]+fnew[i][j][k][1]+fnew[i][j][k][2]+fnew[i][j][k][3]+fnew[i][j][k][4]+
	    fnew[i][j][k][5]+fnew[i][j][k][6]+fnew[i][j][k][7]+fnew[i][j][k][8]+fnew[i][j][k][9]+
	    fnew[i][j][k][10]+fnew[i][j][k][11]+fnew[i][j][k][13]+fnew[i][j][k][15]+fnew[i][j][k][17]+tmp1;
	  
	  tmp3=rb*vwtp-fnew[i][j][k][2]+fnew[i][j][k][4]-fnew[i][j][k][7]+fnew[i][j][k][8]-fnew[i][j][k][9]+
	    fnew[i][j][k][10]-fnew[i][j][k][15]+fnew[i][j][k][17];
	  
	  fnew[i][j][k][12] = 0.25*(tmp1+2.0*tmp2);
	  fnew[i][j][k][14] = 0.25*(tmp1-2.0*tmp2);
	  fnew[i][j][k][16] = 0.25*(tmp1+2.0*tmp3);
	  fnew[i][j][k][18] = 0.25*(tmp1-2.0*tmp3);
	}
      }	
    }
    
    //--------------------------------------------------------------------------
    // Periodic z boundary conditions.
    //--------------------------------------------------------------------------
  }else {
    
    for(i=0; i<numrequests; i++)
      requests[i]=MPI_REQUEST_NULL;
    MPI_Isend(&feq[1][1][1][0],1,passxf,comm->procneigh[0][0],15,world,&requests[0]);
    MPI_Irecv(&feq[0][1][1][0],1,passxf,comm->procneigh[0][0],25,world,&requests[1]);
    MPI_Isend(&feq[subNbx-2][1][1][0],1,passxf,comm->procneigh[0][1],25,world,&requests[2]);
    MPI_Irecv(&feq[subNbx-1][1][1][0],1,passxf,comm->procneigh[0][1],15,world,&requests[3]);
    if(typeLB == 2){
      MPI_Isend(&feqn[1][1][1][0],1,passxf,comm->procneigh[0][0],10,world,&requests[4]);
      MPI_Irecv(&feqn[0][1][1][0],1,passxf,comm->procneigh[0][0],20,world,&requests[5]);
      MPI_Isend(&feqn[subNbx-2][1][1][0],1,passxf,comm->procneigh[0][1],20,world,&requests[6]);
      MPI_Irecv(&feqn[subNbx-1][1][1][0],1,passxf,comm->procneigh[0][1],10,world,&requests[7]);
    }
    update_periodic(2,subNbx-2,2,subNby-2,2,subNbz-2);
    MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);
    
    for(i=0; i<numrequests; i++)
      requests[i]=MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][1][1][0],1,passyf,comm->procneigh[1][0],15,world,&requests[0]);
    MPI_Irecv(&feq[0][0][1][0],1,passyf,comm->procneigh[1][0],25,world,&requests[1]);   
    MPI_Isend(&feq[0][subNby-2][1][0],1,passyf,comm->procneigh[1][1],25,world,&requests[2]);
    MPI_Irecv(&feq[0][subNby-1][1][0],1,passyf,comm->procneigh[1][1],15,world,&requests[3]);
    if(typeLB == 2){
      MPI_Isend(&feqn[0][1][1][0],1,passyf,comm->procneigh[1][0],10,world,&requests[4]);
      MPI_Irecv(&feqn[0][0][1][0],1,passyf,comm->procneigh[1][0],20,world,&requests[5]);   
      MPI_Isend(&feqn[0][subNby-2][1][0],1,passyf,comm->procneigh[1][1],20,world,&requests[6]);
      MPI_Irecv(&feqn[0][subNby-1][1][0],1,passyf,comm->procneigh[1][1],10,world,&requests[7]);
    }
    update_periodic(1,2,2,subNby-2,2,subNbz-2);
    update_periodic(subNbx-2,subNbx-1,2,subNby-2,2,subNbz-2);
    MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);
    
    for(i=0; i<numrequests; i++)
      requests[i]=MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][0][1][0],1,passzf,comm->procneigh[2][0],15,world,&requests[0]);
    MPI_Irecv(&feq[0][0][0][0],1,passzf,comm->procneigh[2][0],25,world,&requests[1]);
    MPI_Isend(&feq[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],25,world,&requests[2]);
    MPI_Irecv(&feq[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],15,world,&requests[3]);
    if(typeLB == 2){
      MPI_Isend(&feqn[0][0][1][0],1,passzf,comm->procneigh[2][0],10,world,&requests[4]);
      MPI_Irecv(&feqn[0][0][0][0],1,passzf,comm->procneigh[2][0],20,world,&requests[5]);
      MPI_Isend(&feqn[0][0][subNbz-2][0],1,passzf,comm->procneigh[2][1],20,world,&requests[6]);
      MPI_Irecv(&feqn[0][0][subNbz-1][0],1,passzf,comm->procneigh[2][1],10,world,&requests[7]);
    }  
    update_periodic(1,subNbx-1,1,2,2,subNbz-2);
    update_periodic(1,subNbx-1,subNby-2,subNby-1,2,subNbz-2);
    MPI_Waitall(numrequests,requests,MPI_STATUS_IGNORE);
    
    update_periodic(1,subNbx-1,1,subNby-1,1,2);
    update_periodic(1,subNbx-1,1,subNby-1,subNbz-2,subNbz-1);
  }
  
}


