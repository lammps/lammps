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

/* ----------------------------------------------------------------------
   Contributing authors: Frances Mackay,
                         Santtu Ollila,
                         Tyson Whitehead,
                         Colin Denniston (UWO)
------------------------------------------------------------------------ */

#include "fix_lb_fluid.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "update.h"

#include <cfloat>
#include <cmath>
#include <cstring>
#include <vector>

#include "latboltz_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_fix_lbfluid[] =
    "fix lb/fluid command: doi:10.1016/j.cpc.2022.108318\n\n"
    "@Article{Denniston et al.,\n"
    "  author = {C. Denniston and N. Afrasiabian and M. G. Cole-Andre,"
    "    F. E. Mackay and S. T. T. Ollila and T. Whitehead},\n"
    " title =   {{LAMMPS} lb/fluid fix version 2: Improved Hydrodynamic "
    "    Forces Implemented into {LAMMPS} Through a Lattice-{B}oltzmann Fluid},"
    " journal = {Comput.\\ Phys.\\ Commun.},\n"
    " year =    2022,\n"
    " volume =  275,\n"
    " pages =   {108318}\n"
    "}\n\n";

/* ------------------------------------------------------------------------ */

namespace LAMMPS_NS {
class Site {
 public:
  int type;
  int orientation;
};
}    // namespace LAMMPS_NS

/* ------------------------------------------------------------------------ */

FixLbFluid::FixLbFluid(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), Gamma(nullptr), hydroF(nullptr), massp(nullptr), density_lb(nullptr),
    u_lb(nullptr), f_lb(nullptr), fnew(nullptr), feq(nullptr), Ff(nullptr), Wf(nullptr),
    Fftempx(nullptr), Wftempx(nullptr), Fftempy(nullptr), Wftempy(nullptr), Fftempz(nullptr),
    Wftempz(nullptr), n_stencil(2), random(nullptr), dof_lb(0), sublattice(nullptr),
    wholelattice(nullptr)
{
  //=====================================================================================================
  //  Sample inputfile call:
  // fix # group lb/fluid nevery viscosity densityinit_real
  //
  //  where: nevery:            call this fix every nevery timesteps.
  //                             (keep this set to 1 for now).
  //         viscosity:         the viscosity of the fluid.
  //         densityinit_real:  the density of the fluid.
  //
  // optional arguments:
  //  "dx" dx_lb:                                    the lattice-Boltzmann grid spacing.
  //  "dm" dm_lb:                                    the lattice-Boltzmann mass unit.
  //  "noise" Temperature seed:                      include noise in the system.
  //                                                  Temperature is the temperature for the fluid.
  //                                                  seed is the seed for the random number generator.
  //  "stencil" n_stencil                            stencil for spread/interpolation of particles
  //                                                  n_stencil = 2 is trilinear stencil
  //                                                  n_stencil = 3 is 3-point standard IBM stencil
  //                                                  n_stencil = 4 is 4-point Keys' interpolation stencil
  //  "read_restart" restart_file:                   restart a fluid run from restart_file.
  //  "write_restart" N:                             write a fluid restart file every N timesteps.
  //  "zwall_velocity" velocity_bottom velocity_top: assign velocities to the z-walls
  //                                                  in the system.
  //  "pressurebcx" pgradav:                         pressure boundary conditions in x-direction. The
  //                                                  pressure jump at the x-boundary is pgradav*Lx*1000
  //  "bodyforce" bodyforcex bodyforcey bodyforcez:  add a constant body force acceleration to the fluid.
  //  "D3Q19":                                       use the 19 velocity D3Q19 model.  By default,
  //                                                  the 15 velocity D3Q15 model is used.
  //  "dumpxdmf" N file timeI:                       dump the fluid density and velocity at each grid
  //                                                  point every N timesteps to xdmf file.  Time will be
  //                                                  indexed by the simulation time if timeI=1 and by the
  //                                                  output frame index (0,1,2,...) if timeI=0.
  //  "linearInit"                                   Initialize density and velocity using linear interpolation
  //                                                  between boundaries. (default is uniform density and
  //                                                  0 velocities.
  //  "dof" dof:                                     specify the number of degrees of freedom to use when
  //                                                  computing temperature (useful when using rigid fix).
  //
  // advanced arguments: (don't use these unless you are an expert)
  //  "scaleGamma" type scale_factor:                scale the default force coupling constant by the
  //                                                  factor, scale_factor, for the given atom type.
  //                                                  Will accept '*' wildtype for type.
  //                                                  Negative value results in IBM (prescribed motion
  //                                                  or stationary node)
  //  "a0" a_0_real:                                 the square of the sound speed in the fluid.
  //
  // Pit Geometry arguments:                         Use these arguments ONLY if you actually want more
  // ----------------------                          complex than rectangular/cubic geometry.
  //                                                  Pit geometry must fill system in x-direction
  //                                                  but can be longer and then be truncated (which
  //                                                  enables asymmetric entrance/exit end sections)
  //  "npits" npits h_p l_p l_pp l_e:                setup to use pit geometry.  Size arguments are measured
  //                                                 in integer multiples of dx_lb
  //                                                  npits = number of pit regions
  //                                                  h_p = z-height of pit regions (floor to bottom of slit)
  //                                                  l_p = x-length of pit regions
  //                                                  l_pp = x-length of slit regions between consecutive pits
  //                                                  l_e = x-length of slit regions at ends
  //  "wp" w_p:                                      y-width of slit regions (defaults to full width)
  //  "sw"                                           y-sidewalls (in xz plane) on/off
  //  Sideview (in xz plane) of pit geometry:
  //    ______________________________________________________________________
  //      slit                          slit                          slit     ^
  //                                                                           |
  //    <---le---><---------lp-------><---lpp---><-------lp--------><---le---> hs = (Nbz-1) - hp
  //                                                                           |
  //    __________                    __________                    __________ v
  //              |                  |          |                  |           ^       z
  //              |                  |          |                  |           |       |
  //              |       pit        |          |       pit        |           hp      +-x
  //              |                  |          |                  |           |
  //              |__________________|          |__________________|           v
  //
  // Endview (in yz plane) of pit geometry (no sw so wp is active):
  //   _____________________
  //                          ^
  //                          |
  //                          hs
  //                          |
  //   _____________________  v
  //       |          |       ^
  //       |          |       |          z
  //       |<---wp--->|       hp         |
  //       |          |       |          +-y
  //       |__________|       v
  //
  // FIX REFERENCE VARIABLES
  //
  // Scalar variable accessible via f_ID in the command script:
  //      Particle temperature (for particles in lb/fluid fix, specify dof using the option above if using
  //                            the rigid fix)
  //
  // Array variables accessible via f_ID[I] in the command script:
  //  I   Variable
  //
  //  1   Temperature of fluid (does NOT account for particles/forces causing regions to move coherently)
  //  2   Total mass of fluid + particles
  //  3   Total x momentem of fluid + particles
  //  4   Total y momentem of fluid + particles
  //  5   Total z momentem of fluid + particles
  //=====================================================================================================

  if (lmp->citeme) lmp->citeme->add(cite_fix_lbfluid);

  // we require continuous time stepping
  time_depend = 1;

  if (narg < 6) error->all(FLERR, "Illegal fix lb/fluid command");

  if (comm->style != 0)
    error->universe_all(FLERR, "Fix lb/fluid can only currently be used with comm_style brick");

  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nprocs);

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  viscosity = utils::numeric(FLERR, arg[4], false, lmp);
  densityinit_real = utils::numeric(FLERR, arg[5], false, lmp);

  // Default values for optional arguments:
  noisestress = 0;
  readrestart = 0;
  printrestart = 0;
  bodyforcex = bodyforcey = bodyforcez = 0.0;
  vwtp = vwbt = 0.0;
  dump_interval = 0;
  T = 300.0;
  dm_lb = -1.0;
  fixviscouslb = 0;
  setdx = 1;
  seta0 = 1;
  setdof = 0;
  numvel = 15;
  pressure = 0;
  rhofactor = 0.0;
  npits = -1;
  sw = 0;
  h_s = h_p = w_p = l_pp = l_e = l_p = 0;
  lin_init = 0;

  const int ntypes = atom->ntypes;
  Gamma = new double[ntypes + 1];
  for (int i = 0; i <= ntypes; i++) Gamma[i] = 1.0;

  // Flags for fix references (i.e. quantities accessible via f_ID[n]
  vector_flag = 1;
  size_vector = 5;

  scalar_flag = 1;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "scaleGamma") == 0) {
      int itype;
      if (strcmp(arg[iarg + 1], "*") == 0)
        itype = 0;
      else
        itype = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      double scalefactor = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if ((itype < 0) || (itype > ntypes))
        error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      if (itype)
        Gamma[itype] = scalefactor;
      else
        for (int it = 1; it <= atom->ntypes; it++) Gamma[it] = scalefactor;
      iarg += 3;
    } else if (strcmp(arg[iarg], "dx") == 0) {
      if ((iarg + 2) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      dx_lb = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
      setdx = 0;
    } else if (strcmp(arg[iarg], "dm") == 0) {
      if ((iarg + 2) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      dm_lb = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "a0") == 0) {
      if ((iarg + 2) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      a_0_real = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
      seta0 = 0;
    } else if (strcmp(arg[iarg], "noise") == 0) {
      if ((iarg + 3) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      noisestress = 1;
      T = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      seed = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      iarg += 3;
    } else if (strcmp(arg[iarg], "stencil") == 0) {
      if ((iarg + 2) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      n_stencil = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "read_restart") == 0) {
      if ((iarg + 2) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      readrestart = 1;
      char *filename = utils::strdup(arg[iarg + 1]);
      MPI_File_open(world, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &pFileRead);
      delete[] filename;
      iarg += 2;
    } else if (strcmp(arg[iarg], "write_restart") == 0) {
      if ((iarg + 2) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      printrestart = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "zwall_velocity") == 0) {
      if ((iarg + 3) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      if (domain->periodicity[2] != 0)
        error->all(FLERR, "setting a z wall velocity without implementing fixed BCs in z");
      vwbt = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      vwtp = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      iarg += 3;
    } else if (strcmp(arg[iarg], "pressurebcx") == 0) {
      if ((iarg + 2) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      pressure = 1;
      rhofactor = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "bodyforce") == 0) {
      if ((iarg + 4) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      bodyforcex = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      bodyforcey = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      bodyforcez = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg], "D3Q19") == 0) {
      numvel = 19;
      iarg += 1;
    } else if (strcmp(arg[iarg], "dumpxdmf") == 0) {
      if ((iarg + 4) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      dump_interval = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      dump_file_name_xdmf = std::string(arg[iarg + 2]) + std::string(".xdmf");
      dump_file_name_raw = std::string(arg[iarg + 2]) + std::string(".raw");
      dump_time_index = utils::inumeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg], "linearInit") == 0) {
      lin_init = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg], "dof") == 0) {
      if ((iarg + 2) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      setdof = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "npits") == 0) {
      if ((iarg + 6) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      npits = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      h_p = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      l_p = utils::inumeric(FLERR, arg[iarg + 3], false, lmp);
      l_pp = utils::inumeric(FLERR, arg[iarg + 4], false, lmp);
      l_e = utils::inumeric(FLERR, arg[iarg + 5], false, lmp);
      iarg += 6;
    } else if (strcmp(arg[iarg], "sw") == 0) {
      sw = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg], "wp") == 0) {
      if ((iarg + 2) > narg) error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
      w_p = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix lb/fluid command: {}", arg[iarg]);
  }

  //-------------------------------------------------------------------------
  //Choose stencil
  //-------------------------------------------------------------------------
  if (n_stencil == 4) {
    interpolate = &FixLbFluid::keys_interpolation;
    interpolationweight = &FixLbFluid::keys_interpolationweight;
  } else if (n_stencil == 3) {
    interpolate = &FixLbFluid::IBM3_interpolation;
    interpolationweight = &FixLbFluid::IBM3_interpolationweight;
  } else {
    interpolate = &FixLbFluid::trilinear_interpolation;
    interpolationweight = &FixLbFluid::trilinear_interpolationweight;
  }

  //--------------------------------------------------------------------------
  //Choose between D3Q15 and D3Q19 functions:
  //--------------------------------------------------------------------------
  if (numvel == 19) {
    equilibriumdist = &FixLbFluid::equilibriumdist19;
    update_full = &FixLbFluid::update_full19;
  } else {
    equilibriumdist = &FixLbFluid::equilibriumdist15;
    update_full = &FixLbFluid::update_full15;
  }

  //--------------------------------------------------------------------------
  // perform initial allocation of atom-based array register
  // with Atom class
  //--------------------------------------------------------------------------
  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  for (int i = 0; i < atom->nmax; i++) {
    massp[i] = 0.0;
    for (int j = 0; j < 3; j++) { hydroF[i][j] = 0.0; }
  }

  //--------------------------------------------------------------------------
  // Set the lattice Boltzmann dt.
  //--------------------------------------------------------------------------
  dt_lb = nevery * (update->dt);

  //--------------------------------------------------------------------------
  // Set the lattice Boltzmann dx if it wasn't specified in the
  // input.
  //--------------------------------------------------------------------------
  if (setdx == 1) {
    double dx_lb1 = sqrt(3.0 * viscosity * dt_lb / densityinit_real);
    double mindomain =
        std::min(std::min(domain->xprd / comm->procgrid[0], domain->yprd / comm->procgrid[1]),
                 domain->zprd / comm->procgrid[2]);
    dx_lb = mindomain / floor(mindomain / dx_lb1);

    if (comm->me == 0) utils::logmesg(lmp, "Setting lattice-Boltzmann dx to {:10.6f}", dx_lb);
  }
  //--------------------------------------------------------------------------
  // Set a0 if it wasn't specified in the input
  //--------------------------------------------------------------------------
  if (seta0 == 1) a_0_real = 0.333333333333333 * dx_lb * dx_lb / dt_lb / dt_lb;

  // convert average pressure gradient into fractional change in density
  rhofactor = rhofactor / 1000 * domain->xprd / densityinit_real / a_0_real;
  if (fabs(rhofactor) > 2) error->all(FLERR, "Illegal pressure jump");
  if (comm->me == 0 && fabs(rhofactor) > 0.15)
    error->warning(FLERR, "Huge pressure jump requested.");

  //--------------------------------------------------------------------------
  // Check to make sure that the total number of grid points in each direction
  // divides evenly among the processors in that direction.
  // Shrink-wrapped boundary conditions (which are not permitted by this fix)
  // might cause a problem, so check for this.  A full check of the boundary
  // conditions is performed in the init routine, rather than here, as it is
  // possible to change the BCs between runs.
  //--------------------------------------------------------------------------
  double aa;
  double eps = 1.0e-8;
  aa = (domain->xprd / comm->procgrid[0]) / dx_lb;
  if (fabs(aa - floor(aa + 0.5)) > eps) {
    if (domain->boundary[0][0] != 0) error->all(FLERR, "the x-direction must be periodic");
    error->all(FLERR,
               "With dx= {}, and the simulation domain divided by {} processors in x "
               "direction, the simulation domain in x direction must be a multiple of {}",
               dx_lb, comm->procgrid[0], comm->procgrid[0] * dx_lb);
  }
  aa = (domain->yprd / comm->procgrid[1]) / dx_lb;
  if (fabs(aa - floor(aa + 0.5)) > eps) {
    if (domain->boundary[1][0] == 2 || domain->boundary[1][0] == 3)
      error->all(FLERR, "the y-direction can not have shrink-wrap boundary conditions");
    error->all(FLERR,
               "With dx= {}, and the simulation domain divided by {} processors in y "
               "direction, the simulation domain in y direction must be a multiple of {}",
               dx_lb, comm->procgrid[1], comm->procgrid[1] * dx_lb);
  }
  aa = (domain->zprd / comm->procgrid[2]) / dx_lb;
  if (fabs(aa - floor(aa + 0.5)) > eps) {
    if (domain->boundary[2][0] == 2 || domain->boundary[2][0] == 3) {
      error->all(FLERR, "the z-direction can not have shrink-wrap boundary conditions");
    }
    error->all(FLERR,
               "With dx= {}, and the simulation domain divided by {} processors in z "
               "direction, the simulation domain in z direction must be a multiple of {}",
               dx_lb, comm->procgrid[2], comm->procgrid[2] * dx_lb);
  }

  //--------------------------------------------------------------------------
  // Set the total number of grid points in each direction.
  //--------------------------------------------------------------------------
  Nbx = (int) (domain->xprd / dx_lb + 0.5);
  Nby = (int) (domain->yprd / dx_lb + 0.5);
  Nbz = (int) (domain->zprd / dx_lb + 0.5);

  //--------------------------------------------------------------------------
  // Set the number of grid points in each dimension for the local subgrids.
  //--------------------------------------------------------------------------
  subNbx = Nbx / comm->procgrid[0] + 2;
  subNby = Nby / comm->procgrid[1] + 2;
  subNbz = Nbz / comm->procgrid[2] + 2;

  //--------------------------------------------------------------------------
  // In order to calculate the fluid forces correctly, need to have atleast
  // 5 grid points in each direction per processor.
  //--------------------------------------------------------------------------
  if (subNbx < 7 || subNby < 7 || subNbz < 7)
    error->all(FLERR, "Need at least 5 grid points in each direction per processor");

  // If there are walls in the z-direction add an extra grid point.
  if (domain->periodicity[2] == 0) {
    Nbz += 1;
    if (comm->myloc[2] == comm->procgrid[2] - 1) subNbz += 1;
  }

  if (comm->me == 0) {
    utils::logmesg(lmp, "Using a lattice-Boltzmann grid of {} by {} by {} points. ", Nbx, Nby, Nbz);
    if (setdx == 1)
      utils::logmesg(lmp, "To change, use the dx keyword\n");
    else
      utils::logmesg(lmp, "\n");
  }

  //--------------------------------------------------------------------------
  // Grid checks for pits
  //--------------------------------------------------------------------------
  if (npits != -1) {
    h_s = (Nbz - 1) - h_p;
    if (w_p == 0) w_p = Nby;
    if (comm->me == 0) {
      if (h_s <= 1) error->all(FLERR, "hp too big for system in z-direction");
      if (w_p > Nby) error->all(FLERR, "wp bigger than system in y-direction");
      if (w_p && sw)
        utils::logmesg(lmp, "wp ignored with side walls, pits are full width in y-direction\n");
      if (2 * l_e + npits * l_p + (npits - 1) * l_pp < Nbx)
        error->all(FLERR,
                   "length of pits and end segments too small to fill system in x-direction");
      else if (2 * l_e + npits * l_p + (npits - 1) * l_pp > Nbx)
        utils::logmesg(lmp,
                       "length of pits and end segments larger than system "
                       "in x-direction: truncation will occur\n");
      if (numvel == 19)
        error->all(FLERR, "Pit geometry options not available for D3Q19, use D3Q15 instead");
      if (vwtp != 0.0 || vwbt != 0.0)
        error->all(FLERR, "Moving walls not compatible with pit geometry options");
    }
  } else {
    if (sw && comm->me == 0)
      error->all(FLERR,
                 "Sidewalls require pit geometry (you can enable by setting npits to 0 or higher)");
  }

  //--------------------------------------------------------------------------
  // Store the largest value of subNbz, which is needed for allocating the
  // buf array (since a processor with comm->myloc[2] == comm->procgrid[2]-1
  // may have an additional subNbz point as compared with the rest).
  //--------------------------------------------------------------------------
  int subNbzmax;
  MPI_Allreduce(&subNbz, &subNbzmax, 1, MPI_INT, MPI_MAX, world);

  //--------------------------------------------------------------------------
  // Create the MPI datatypes used to pass portions of arrays:
  //--------------------------------------------------------------------------
  // MPI 3-vector
  MPI_Type_contiguous(3, MPI_DOUBLE, &realType3_mpitype);
  MPI_Type_commit(&realType3_mpitype);

  // datatypes to pass the f and feq arrays.
  MPI_Aint lb, sizeofdouble;
  MPI_Type_get_extent(MPI_DOUBLE, &lb, &sizeofdouble);

  MPI_Type_vector(subNbz - 2, numvel, numvel, MPI_DOUBLE, &oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNby - 2, 1, numvel * subNbz * sizeofdouble, oneslice, &passxf);
  MPI_Type_commit(&passxf);

  MPI_Type_create_hvector(subNbx, 1, numvel * subNbz * subNby * sizeofdouble, oneslice, &passyf);
  MPI_Type_commit(&passyf);

  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNby, numvel, numvel * subNbz, MPI_DOUBLE, &oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNbx, 1, numvel * subNbz * subNby * sizeofdouble, oneslice, &passzf);
  MPI_Type_commit(&passzf);

  // datatypes to pass the u array, and the Ff array.
  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNbz + 3, 3, 3, MPI_DOUBLE, &oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNby + 3, 1, 3 * (subNbz + 3) * sizeofdouble, oneslice, &passxu);
  MPI_Type_commit(&passxu);

  MPI_Type_create_hvector(subNbx + 3, 1, 3 * (subNbz + 3) * (subNby + 3) * sizeofdouble, oneslice,
                          &passyu);
  MPI_Type_commit(&passyu);

  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNby + 3, 3, 3 * (subNbz + 3), MPI_DOUBLE, &oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNbx + 3, 1, 3 * (subNbz + 3) * (subNby + 3) * sizeofdouble, oneslice,
                          &passzu);
  MPI_Type_commit(&passzu);

  // datatypes to pass the density array, the Wf array.
  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNbz + 3, 1, 1, MPI_DOUBLE, &oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNby + 3, 1, 1 * (subNbz + 3) * sizeofdouble, oneslice, &passxrho);
  MPI_Type_commit(&passxrho);

  MPI_Type_create_hvector(subNbx + 3, 1, 1 * (subNbz + 3) * (subNby + 3) * sizeofdouble, oneslice,
                          &passyrho);
  MPI_Type_commit(&passyrho);

  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNby + 3, 1, 1 * (subNbz + 3), MPI_DOUBLE, &oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNbx + 3, 1, 1 * (subNbz + 3) * (subNby + 3) * sizeofdouble, oneslice,
                          &passzrho);
  MPI_Type_commit(&passzrho);

  // datatypes to receive a portion of the Ff array.
  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNbz + 3, 3, 3, MPI_DOUBLE, &oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNby + 3, 1, 3 * (subNbz + 3) * sizeofdouble, oneslice, &passxtemp);
  MPI_Type_commit(&passxtemp);

  MPI_Type_create_hvector(subNbx + 3, 1, 3 * (subNbz + 3) * 5 * sizeofdouble, oneslice, &passytemp);
  MPI_Type_commit(&passytemp);

  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNby + 3, 3, 3 * 5, MPI_DOUBLE, &oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNbx + 3, 1, 3 * 5 * (subNby + 3) * sizeofdouble, oneslice, &passztemp);
  MPI_Type_commit(&passztemp);

  // datatypes to receive a portion of the Wf and Mf array.
  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNbz + 3, 1, 1, MPI_DOUBLE, &oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNby + 3, 1, 1 * (subNbz + 3) * sizeofdouble, oneslice, &passxWtemp);
  MPI_Type_commit(&passxWtemp);

  MPI_Type_create_hvector(subNbx + 3, 1, 1 * (subNbz + 3) * 5 * sizeofdouble, oneslice,
                          &passyWtemp);
  MPI_Type_commit(&passyWtemp);

  MPI_Type_free(&oneslice);
  MPI_Type_vector(subNby + 3, 1, 1 * 5, MPI_DOUBLE, &oneslice);
  MPI_Type_commit(&oneslice);
  MPI_Type_create_hvector(subNbx + 3, 1, 1 * 5 * (subNby + 3) * sizeofdouble, oneslice,
                          &passzWtemp);
  MPI_Type_commit(&passzWtemp);

  MPI_Type_free(&oneslice);

  //--------------------------------------------------------------------------
  // Allocate the necessary arrays.
  //--------------------------------------------------------------------------
  memory->create(feq, subNbx, subNby, subNbz, numvel, "FixLbFluid:feq");
  memory->create(f_lb, subNbx, subNby, subNbz, numvel, "FixLbFluid:f_lb");
  memory->create(fnew, subNbx, subNby, subNbz, numvel, "FixLbFluid:fnew");
  memory->create(density_lb, subNbx + 3, subNby + 3, subNbz + 3, "FixLbFluid:density_lb");
  memory->create(u_lb, subNbx + 3, subNby + 3, subNbz + 3, 3, "FixLbFluid:u_lb");
  memory->create(Ff, subNbx + 3, subNby + 3, subNbz + 3, 3, "FixLbFluid:Ff");
  std::fill(&Ff[0][0][0][0], &Ff[0][0][0][0] + (subNbx + 3) * (subNby + 3) * (subNbz + 3) * 3, 0.0);
  memory->create(Fftempx, 5, subNby + 3, subNbz + 3, 3, "FixLbFluid:Fftempx");
  memory->create(Fftempy, subNbx + 3, 5, subNbz + 3, 3, "FixLbFluid:Fftempy");
  memory->create(Fftempz, subNbx + 3, subNby + 3, 5, 3, "FixLbFluid:Fftempz");
  memory->create(Wf, subNbx + 3, subNby + 3, subNbz + 3, "FixLbFluid:Wf");
  std::fill(&Wf[0][0][0], &Wf[0][0][0] + (subNbx + 3) * (subNby + 3) * (subNbz + 3), 0.0);
  memory->create(Wftempx, 5, subNby + 3, subNbz + 3, "FixLbFluid:Wftempx");
  memory->create(Wftempy, subNbx + 3, 5, subNbz + 3, "FixLbFluid:Wftempy");
  memory->create(Wftempz, subNbx + 3, subNby + 3, 5, "FixLbFluid:Wftempz");

  if (noisestress == 1) { random = new RanMars(lmp, seed + comm->me); }

  //--------------------------------------------------------------------------
  // Rescale the variables to Lattice Boltzmann dimensionless units.
  //--------------------------------------------------------------------------
  rescale();

  //Intialize Timers
  timeEqb = timeUpdate = timePCalc = timefluidForce = timeCorrectU = 0.0;

  // Setup buffers for output.
  SetupBuffers();
  InitializeFirstRun();
}

FixLbFluid::~FixLbFluid()
{
  // Timer output
  double timeAll;
  MPI_Reduce(&timeEqb, &timeAll, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  if (me == 0) printf("\nLB equilibriumDist time: %f\n", timeAll / nprocs);
  MPI_Reduce(&timeUpdate, &timeAll, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  if (me == 0) printf("LB update time: %f\n", timeAll / nprocs);
  MPI_Reduce(&timePCalc, &timeAll, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  if (me == 0) printf("LB PCalc time: %f\n", timeAll / nprocs);
  MPI_Reduce(&timefluidForce, &timeAll, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  if (me == 0) printf("LB fluidForce time: %f\n", timeAll / nprocs);
  MPI_Reduce(&timeCorrectU, &timeAll, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  if (me == 0) printf("LB CorrectU time: %f\n", timeAll / nprocs);

  // Close off output
  if (dump_interval) {
    if (me == 0) {
      fprintf(dump_file_handle_xdmf, "    </Grid>\n  </Domain>\n</Xdmf>\n");
      if (fclose(dump_file_handle_xdmf))
        error->one(FLERR, "Unable to close \"{}\": {}", dump_file_name_xdmf, utils::getsyserror());
    }
    MPI_File_close(&dump_file_handle_raw);
  }
  MPI_Type_free(&realType3_mpitype);
  MPI_Type_free(&passxf);
  MPI_Type_free(&passyf);
  MPI_Type_free(&passzf);
  MPI_Type_free(&passxu);
  MPI_Type_free(&passyu);
  MPI_Type_free(&passzu);
  MPI_Type_free(&passxrho);
  MPI_Type_free(&passyrho);
  MPI_Type_free(&passzrho);
  MPI_Type_free(&passxtemp);
  MPI_Type_free(&passytemp);
  MPI_Type_free(&passztemp);
  MPI_Type_free(&passxWtemp);
  MPI_Type_free(&passyWtemp);
  MPI_Type_free(&passzWtemp);
  MPI_Type_free(&fluid_density_2_mpitype);
  MPI_Type_free(&fluid_velocity_2_mpitype);
  MPI_Type_free(&dump_file_mpitype);

  // Unregister fix callback
  atom->delete_callback(id, Atom::GROW);

  // Clean up memory
  memory->destroy(hydroF);
  memory->destroy(massp);

  memory->destroy(feq);
  memory->destroy(f_lb);
  memory->destroy(fnew);
  memory->destroy(density_lb);
  memory->destroy(u_lb);

  memory->destroy(Ff);
  memory->destroy(Fftempx);
  memory->destroy(Fftempy);
  memory->destroy(Fftempz);
  memory->destroy(Wf);
  memory->destroy(Wftempx);
  memory->destroy(Wftempy);
  memory->destroy(Wftempz);

  if (noisestress == 1) { delete random; }

  delete[] Gamma;

  memory->destroy(sublattice); /* nanopit arrays */
}

int FixLbFluid::setmask()
{
  return FixConst::INITIAL_INTEGRATE | FixConst::POST_FORCE | FixConst::FINAL_INTEGRATE |
      FixConst::END_OF_STEP;
}

void FixLbFluid::init()
{
  if (modify->get_fix_by_style("lb/fluid").size() > 1)
    error->all(FLERR, "Only one fix lb/fluid at a time is supported");

  if (modify->get_fix_by_style("dt/reset").size() > 1)
    error->all(FLERR, "Fix lb/fluid is not compatible with fix dt/reset");

  //--------------------------------------------------------------------------
  // Check to see if the MD timestep has changed between runs.
  //--------------------------------------------------------------------------
  double dt_lb_now;
  dt_lb_now = nevery * (update->dt);

  if (fabs(dt_lb_now - dt_lb) > 1.0e-12)
    error->warning(FLERR, "Timestep changed between runs.  Must re-issue fix lb/fluid command.");

  //--------------------------------------------------------------------------
  // Make sure the size of the simulation domain has not changed
  // between runs.
  //--------------------------------------------------------------------------
  int Nbx_now, Nby_now, Nbz_now;
  Nbx_now = (int) (domain->xprd / dx_lb + 0.5);
  Nby_now = (int) (domain->yprd / dx_lb + 0.5);
  Nbz_now = (int) (domain->zprd / dx_lb + 0.5);
  // If there are walls in the z-direction add an extra grid point.
  if (domain->periodicity[2] == 0) { Nbz_now += 1; }

  if (Nbx_now != Nbx || Nby_now != Nby || Nbz_now != Nbz) {
    error->all(FLERR, "Simulation domain can not change shape between runs with the same lb/fluid");
  }

  //--------------------------------------------------------------------------
  // Check to make sure that the chosen LAMMPS boundary types are compatible
  // with this fix.
  //    shrink-wrap is not compatible in any dimension.
  //    fixed only works in the y (pits sw) or z-direction (pits or zwall_vel.)
  //--------------------------------------------------------------------------
  if (domain->boundary[0][0] != 0) { error->all(FLERR, "the x-direction must be periodic"); }
  if (domain->boundary[1][0] == 2 || domain->boundary[1][0] == 3) {
    error->all(FLERR, "the y-direction cannot have shrink-wrap boundary conditions");
  }
  if (comm->me == 0 && domain->boundary[1][0] == 0 && sw) {
    error->warning(FLERR,
                   "Particle-Particle interactions are periodic in y-direction but there are fixed "
                   "boundaries in y-direction for lb/fluid.");
  }
  if (comm->me == 0 && domain->boundary[1][0] == 1 && sw == 0) {
    error->warning(FLERR,
                   "lb/fluid boundaries are periodic in y-direction but there are fixed boundies "
                   "in y-direction for particle-particle interactions.");
  }
  if (domain->boundary[2][0] == 2 || domain->boundary[2][0] == 3) {
    error->all(FLERR, "the z-direction cannot have shrink-wrap boundary conditions");
  }
  if (domain->boundary[2][0] == 0 && npits != -1) {
    error->all(FLERR,
               "the z-direction cannot have periodic boundary conditions when there are pits.");
  }

  //--------------------------------------------------------------------------
  // Check if the lb/viscous fix is also called:
  //--------------------------------------------------------------------------
  groupbit_viscouslb = 0;
  auto fixlbv = modify->get_fix_by_style("lb/viscous");
  if (fixlbv.size() > 0) {
    if (fixlbv.size() > 1)
      error->all(FLERR, "More than one fix lb/viscous at a time is not supported");
    fixviscouslb = 1;
    groupbit_viscouslb = group->bitmask[fixlbv[0]->igroup];
  }

  // Warn if the fluid force is not applied to any of the particles.
  if (!(groupbit_viscouslb) && comm->me == 0)
    utils::logmesg(lmp,
                   "Not adding the fluid force to any of the MD particles.  "
                   "To add this force use lb/viscous fix\n");
}

void FixLbFluid::setup(int /* vflag */)
{
  //--------------------------------------------------------------------------
  // We could calculate the force on the fluid for a restart run here but that
  // would not match the original continuation.  In fact, the recalcuated
  // forces would actually be zero if the algorithm were exact whereas the
  // original forces would be nonzero.  As such, seems better to just leave
  // the hydro forces as zero for the first 1/2 step on restart.
  //--------------------------------------------------------------------------
}

void FixLbFluid::initial_integrate(int /* vflag */)
{
  //--------------------------------------------------------------------------
  // Determine the equilibrium distribution on the local subgrid.
  //--------------------------------------------------------------------------
  double st = MPI_Wtime();
  (*this.*equilibriumdist)(1, subNbx - 1, 1, subNby - 1, 1, subNbz - 1);
  timeEqb += MPI_Wtime() - st;

  //--------------------------------------------------------------------------
  // Using the equilibrium distribution, calculate the new
  // distribution function.
  //--------------------------------------------------------------------------
  st = MPI_Wtime();
  (*this.*update_full)();
  timeUpdate += MPI_Wtime() - st;

  std::swap(f_lb, fnew);

  //--------------------------------------------------------------------------
  // Calculate moments of the distribution function.
  //--------------------------------------------------------------------------
  st = MPI_Wtime();
  parametercalc_full();
  timePCalc += MPI_Wtime() - st;
}

void FixLbFluid::post_force(int /*vflag*/)
{
  if (fixviscouslb == 1) {
    double st = MPI_Wtime();
    calc_fluidforceweight();
    calc_fluidforceI();
    timefluidForce += MPI_Wtime() - st;
  }
}

void FixLbFluid::final_integrate()
{
  // Correct u by adding in the force contribution, but only on local subgrid
  double st = MPI_Wtime();
  correctu_full();
  timeCorrectU += MPI_Wtime() - st;
}

void FixLbFluid::end_of_step()
{
  if (printrestart > 0) {
    if ((update->ntimestep) % printrestart == 0) { write_restartfile(); }
  }

  // Output fluid to dumpfile
  dump(update->ntimestep);
}

//==========================================================================
//   allocate atom-based array
//==========================================================================
void FixLbFluid::grow_arrays(int nmax)
{
  memory->grow(hydroF, nmax, 3, "FixLbFluid:hydroF");
  memory->grow(massp, nmax, "FixLbFluid:massp");
}

//==========================================================================
//   copy values within local atom-based array
//==========================================================================
void FixLbFluid::copy_arrays(int i, int j, int /* delflag */)
{
  hydroF[j][0] = hydroF[i][0];
  hydroF[j][1] = hydroF[i][1];
  hydroF[j][2] = hydroF[i][2];
  massp[j] = massp[i];
}

//==========================================================================
//   pack values in local atom-based array for exchange with another proc
//==========================================================================
int FixLbFluid::pack_exchange(int i, double *buf)
{
  buf[0] = hydroF[i][0];
  buf[1] = hydroF[i][1];
  buf[2] = hydroF[i][2];
  buf[3] = massp[i];
  return 4;
}

//==========================================================================
//   unpack values in local atom-based array from exchange with another proc
//==========================================================================
int FixLbFluid::unpack_exchange(int nlocal, double *buf)
{
  hydroF[nlocal][0] = buf[0];
  hydroF[nlocal][1] = buf[1];
  hydroF[nlocal][2] = buf[2];
  massp[nlocal] = buf[3];
  return 4;
}

//==========================================================================
// rescale the simulation parameters so that dx_lb=dt_lb=dm_lb=1.
// This assumes that all the simulation parameters have been given in
// terms of distance, time and mass units.
//==========================================================================
void FixLbFluid::rescale()
{

  densityinit = densityinit_real * dx_lb * dx_lb * dx_lb;
  if (dm_lb < 0) {    // indicating it was not set by user
    dm_lb = densityinit;
    densityinit = 1.0;
  } else
    densityinit /= dm_lb;

  vwtp = vwtp * dt_lb / dx_lb;
  vwbt = vwbt * dt_lb / dx_lb;

  bodyforcex = bodyforcex * dt_lb * dt_lb / dx_lb;
  bodyforcey = bodyforcey * dt_lb * dt_lb / dx_lb;
  bodyforcez = bodyforcez * dt_lb * dt_lb / dx_lb;

  if (pressure) {
    rhoH = (1.0 + rhofactor * 0.5) * densityinit_real;
    rhoL = (1.0 - rhofactor * 0.5) * densityinit_real;
    rhoH = rhoH * dx_lb * dx_lb * dx_lb / dm_lb;
    rhoL = rhoL * dx_lb * dx_lb * dx_lb / dm_lb;
  } else
    rhoH = rhoL = densityinit;

  tau = (3.0 * viscosity / densityinit_real) * dt_lb * dt_lb / dx_lb / dx_lb;
  tau /= dt_lb;
  tau = tau + 0.5;

  a_0 = a_0_real * dt_lb * dt_lb / (dx_lb * dx_lb);

  // Warn if using the D3Q19 model with noise, and a0 is too small.
  if (numvel == 19 && noisestress == 1 && a_0 < 0.2) {
    error->warning(FLERR, "Fix lb/fluid WARNING: Chosen value for a0 may be too small. \
Check temperature reproduction.\n");
  }

  if (noisestress == 1) {
    if (a_0 > 0.5555555) {
      error->all(FLERR, "Fix lb/fluid ERROR: the Lattice Boltzmann dx and dt need \
to be chosen such that the scaled a_0 < 5/9\n");
    }
  }

  // Courant Condition:
  if (a_0 >= 1.0) {
    error->all(FLERR, "Fix lb/fluid ERROR: the lattice Boltzmann dx and dt do not \
satisfy the Courant condition.\n");
  }

  kB = (force->boltz / force->mvv2e) * dt_lb * dt_lb / dx_lb / dx_lb / dm_lb;

  namp = 2.0 * kB * T * (tau - 0.5) / 3.0;
  noisefactor = 1.0;
  if (a_0 <= 0.333333333333333) {
    K_0 = 5.17 * (0.333333333333333 - a_0);
  } else {
    K_0 = 2.57 * (a_0 - 0.333333333333333);
  }
}

void FixLbFluid::InitializeFirstRun()
{

  // Define geometry for pits

  // Allocate the necessary arrays
  memory->create(wholelattice, Nbx, Nby, Nbz, "FixLBFluid:lattice");
  memory->create(sublattice, subNbx, subNby, subNbz, "FixLBFluid:sublattice");

  // Initialize global lattice geometry.
  initializeGlobalGeometry();

  // Initialize local lattice geometry based on global geometry.
  initializeGeometry();
  if (comm->me == 0) utils::logmesg(lmp, "Local Grid Geometry created.\n");

  // Destroy redundant global lattice.
  memory->destroy(wholelattice);
  MPI_Barrier(world);

  //--------------------------------------------------------------------------
  // Initialize the arrays.
  //--------------------------------------------------------------------------

  if (readrestart == 0) {
    step = 0;

    initializeLB();    // initialize density, fluid velocity, f_lb

  } else {
    step = 1;

    read_restartfile();    // get f_lb from restart file
  }
  parametercalc_full();    // necessary after restart, consistently fill arrays for regular init

  // When we calculate feq in initial_integrate we will need forces, as hydro forces are
  // velocity dependent and use the 1/2 step velocities which are not saved, we leave forces
  // at their initilized zero for first half-step

  // Output for t=0
  dump(update->ntimestep);

  if (me == 0) utils::logmesg(lmp, "First Run initialized\n");
}

//==========================================================================
// Initialize the fluid velocity and density.
//==========================================================================

void FixLbFluid::initializeLB()
{
  for (int i = 0; i < subNbx + 3; i++) {

    int ix;
    if (i == subNbx + 2)
      ix = round(((domain->sublo[0] - domain->boxlo[0]) / dx_lb)) - 2;
    else
      ix = round(((domain->sublo[0] - domain->boxlo[0]) / dx_lb)) + (i - 1);

    for (int j = 0; j < subNby + 3; j++)
      for (int k = 0; k < subNbz + 3; k++) {

        if (lin_init)
          density_lb[i][j][k] = rhoH +
              (rhoL - rhoH) * ix / (Nbx - 1);    //linear interpolation between boundary values
        else
          density_lb[i][j][k] = densityinit;

        int iz;
        if (k == subNbz + 2)
          iz = round(((domain->sublo[2] - domain->boxlo[2]) / dx_lb)) - 2;
        else
          iz = round(((domain->sublo[2] - domain->boxlo[2]) / dx_lb)) + (k - 1);

        u_lb[i][j][k][0] = 0.0;
        if (lin_init)
          u_lb[i][j][k][1] = vwbt +
              (vwtp - vwbt) * iz / (Nbz - 1);    //linear interpolation between boundary velocities
        else
          u_lb[i][j][k][1] = 0.0;
        u_lb[i][j][k][2] = 0.0;
      }
  }

  for (int i = 0; i < subNbx; i++)
    for (int j = 0; j < subNby; j++)
      for (int k = 0; k < subNbz; k++)
        if (numvel == 15)
          for (int m = 0; m < 15; m++)
            f_lb[i][j][k][m] = density_lb[i][j][k] * w_lb15[m] *
                (1 +
                 (u_lb[i][j][k][0] * e15[m][0] + u_lb[i][j][k][1] * e15[m][1] +
                  u_lb[i][j][k][2] * e15[m][2]) /
                     a_0);
        else
          for (int m = 0; m < 19; m++)
            f_lb[i][j][k][m] = density_lb[i][j][k] * w_lb19[m] *
                (1 +
                 (u_lb[i][j][k][0] * e19[m][0] + u_lb[i][j][k][1] * e19[m][1] +
                  u_lb[i][j][k][2] * e19[m][2]) /
                     a_0);
}

//==========================================================================
//   calculate the force from the local atoms acting on the fluid.
//==========================================================================
void FixLbFluid::calc_fluidforceI()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i, j, k, m;

  //--------------------------------------------------------------------------
  // Zero out arrays
  //--------------------------------------------------------------------------
  std::fill(&Ff[0][0][0][0], &Ff[0][0][0][0] + (subNbx + 3) * (subNby + 3) * (subNbz + 3) * 3, 0.0);
  std::fill(&Fftempx[0][0][0][0], &Fftempx[0][0][0][0] + 5 * (subNby + 3) * (subNbz + 3) * 3, 0.0);
  std::fill(&Fftempy[0][0][0][0], &Fftempy[0][0][0][0] + (subNbx + 3) * 5 * (subNbz + 3) * 3, 0.0);
  std::fill(&Fftempz[0][0][0][0], &Fftempz[0][0][0][0] + (subNbx + 3) * (subNby + 3) * 5 * 3, 0.0);

  std::fill(&Wftempx[0][0][0], &Wftempx[0][0][0] + 5 * (subNby + 3) * (subNbz + 3), 0.0);
  std::fill(&Wftempy[0][0][0], &Wftempy[0][0][0] + (subNbx + 3) * 5 * (subNbz + 3), 0.0);
  std::fill(&Wftempz[0][0][0], &Wftempz[0][0][0] + (subNbx + 3) * (subNby + 3) * 5, 0.0);

  //--------------------------------------------------------------------------
  //Calculate the contribution to the force on the fluid.
  //--------------------------------------------------------------------------
  dof_lb = 0.0;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) { (*this.*interpolate)(i, 0); }
  }

  //--------------------------------------------------------------------------
  //Communicate the force contributions which lie outside the local processor
  //sub domain.
  //--------------------------------------------------------------------------
  MPI_Request requests[10];
  for (i = 0; i < 10; i++) requests[i] = MPI_REQUEST_NULL;
  MPI_Isend(&Ff[0][0][0][0], 1, passxu, comm->procneigh[0][0], 10, world, &requests[0]);
  MPI_Isend(&Ff[subNbx + 2][0][0][0], 1, passxu, comm->procneigh[0][0], 20, world, &requests[1]);
  MPI_Isend(&Ff[subNbx - 1][0][0][0], 1, passxu, comm->procneigh[0][1], 30, world, &requests[2]);
  MPI_Isend(&Ff[subNbx][0][0][0], 1, passxu, comm->procneigh[0][1], 40, world, &requests[3]);
  MPI_Isend(&Ff[subNbx + 1][0][0][0], 1, passxu, comm->procneigh[0][1], 50, world, &requests[4]);
  MPI_Irecv(&Fftempx[0][0][0][0], 1, passxtemp, comm->procneigh[0][1], 10, world, &requests[5]);
  MPI_Irecv(&Fftempx[1][0][0][0], 1, passxtemp, comm->procneigh[0][1], 20, world, &requests[6]);
  MPI_Irecv(&Fftempx[2][0][0][0], 1, passxtemp, comm->procneigh[0][0], 30, world, &requests[7]);
  MPI_Irecv(&Fftempx[3][0][0][0], 1, passxtemp, comm->procneigh[0][0], 40, world, &requests[8]);
  MPI_Irecv(&Fftempx[4][0][0][0], 1, passxtemp, comm->procneigh[0][0], 50, world, &requests[9]);

  MPI_Waitall(10, requests, MPI_STATUS_IGNORE);
  for (j = 0; j < subNby + 3; j++) {
    for (k = 0; k < subNbz + 3; k++) {
      for (m = 0; m < 3; m++) {
        Ff[subNbx - 2][j][k][m] += Fftempx[0][j][k][m];
        Ff[subNbx - 3][j][k][m] += Fftempx[1][j][k][m];
        Ff[1][j][k][m] += Fftempx[2][j][k][m];
        Ff[2][j][k][m] += Fftempx[3][j][k][m];
        Ff[3][j][k][m] += Fftempx[4][j][k][m];
      }
    }
  }

  for (i = 0; i < 10; i++) requests[i] = MPI_REQUEST_NULL;
  MPI_Isend(&Ff[0][0][0][0], 1, passyu, comm->procneigh[1][0], 10, world, &requests[0]);
  MPI_Isend(&Ff[0][subNby + 2][0][0], 1, passyu, comm->procneigh[1][0], 20, world, &requests[1]);
  MPI_Isend(&Ff[0][subNby - 1][0][0], 1, passyu, comm->procneigh[1][1], 30, world, &requests[2]);
  MPI_Isend(&Ff[0][subNby][0][0], 1, passyu, comm->procneigh[1][1], 40, world, &requests[3]);
  MPI_Isend(&Ff[0][subNby + 1][0][0], 1, passyu, comm->procneigh[1][1], 50, world, &requests[4]);
  MPI_Irecv(&Fftempy[0][0][0][0], 1, passytemp, comm->procneigh[1][1], 10, world, &requests[5]);
  MPI_Irecv(&Fftempy[0][1][0][0], 1, passytemp, comm->procneigh[1][1], 20, world, &requests[6]);
  MPI_Irecv(&Fftempy[0][2][0][0], 1, passytemp, comm->procneigh[1][0], 30, world, &requests[7]);
  MPI_Irecv(&Fftempy[0][3][0][0], 1, passytemp, comm->procneigh[1][0], 40, world, &requests[8]);
  MPI_Irecv(&Fftempy[0][4][0][0], 1, passytemp, comm->procneigh[1][0], 50, world, &requests[9]);

  MPI_Waitall(10, requests, MPI_STATUS_IGNORE);
  for (i = 0; i < subNbx + 3; i++) {
    for (k = 0; k < subNbz + 3; k++) {
      for (m = 0; m < 3; m++) {
        Ff[i][subNby - 2][k][m] += Fftempy[i][0][k][m];
        Ff[i][subNby - 3][k][m] += Fftempy[i][1][k][m];
        Ff[i][1][k][m] += Fftempy[i][2][k][m];
        Ff[i][2][k][m] += Fftempy[i][3][k][m];
        Ff[i][3][k][m] += Fftempy[i][4][k][m];
      }
    }
  }

  for (i = 0; i < 10; i++) requests[i] = MPI_REQUEST_NULL;
  MPI_Isend(&Ff[0][0][0][0], 1, passzu, comm->procneigh[2][0], 10, world, &requests[0]);
  MPI_Isend(&Ff[0][0][subNbz + 2][0], 1, passzu, comm->procneigh[2][0], 20, world, &requests[1]);
  MPI_Isend(&Ff[0][0][subNbz - 1][0], 1, passzu, comm->procneigh[2][1], 30, world, &requests[2]);
  MPI_Isend(&Ff[0][0][subNbz][0], 1, passzu, comm->procneigh[2][1], 40, world, &requests[3]);
  MPI_Isend(&Ff[0][0][subNbz + 1][0], 1, passzu, comm->procneigh[2][1], 50, world, &requests[4]);
  MPI_Irecv(&Fftempz[0][0][0][0], 1, passztemp, comm->procneigh[2][1], 10, world, &requests[5]);
  MPI_Irecv(&Fftempz[0][0][1][0], 1, passztemp, comm->procneigh[2][1], 20, world, &requests[6]);
  MPI_Irecv(&Fftempz[0][0][2][0], 1, passztemp, comm->procneigh[2][0], 30, world, &requests[7]);
  MPI_Irecv(&Fftempz[0][0][3][0], 1, passztemp, comm->procneigh[2][0], 40, world, &requests[8]);
  MPI_Irecv(&Fftempz[0][0][4][0], 1, passztemp, comm->procneigh[2][0], 50, world, &requests[9]);
  MPI_Waitall(10, requests, MPI_STATUS_IGNORE);

  for (i = 0; i < subNbx + 3; i++) {
    for (j = 0; j < subNby + 3; j++) {
      for (m = 0; m < 3; m++) {
        Ff[i][j][subNbz - 2][m] += Fftempz[i][j][0][m];
        Ff[i][j][subNbz - 3][m] += Fftempz[i][j][1][m];
        Ff[i][j][1][m] += Fftempz[i][j][2][m];
        Ff[i][j][2][m] += Fftempz[i][j][3][m];
        Ff[i][j][3][m] += Fftempz[i][j][4][m];
      }
    }
  }
}

void FixLbFluid::calc_fluidforceII()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i;

  //--------------------------------------------------------------------------
  //Calculate the contribution to the force on the fluid.
  //--------------------------------------------------------------------------
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) { (*this.*interpolate)(i, 1); }
  }
}

//==========================================================================
//   calculate the force weight from the local atoms acting on the fluid.
//==========================================================================
void FixLbFluid::calc_fluidforceweight()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  MPI_Request requests[20];

  //--------------------------------------------------------------------------
  // Zero out arrays
  //--------------------------------------------------------------------------
  std::fill(&Wf[0][0][0], &Wf[0][0][0] + (subNbx + 3) * (subNby + 3) * (subNbz + 3), 0.0);
  std::fill(&Wftempx[0][0][0], &Wftempx[0][0][0] + 5 * (subNby + 3) * (subNbz + 3), 0.0);
  std::fill(&Wftempy[0][0][0], &Wftempy[0][0][0] + (subNbx + 3) * 5 * (subNbz + 3), 0.0);
  std::fill(&Wftempz[0][0][0], &Wftempz[0][0][0] + (subNbx + 3) * (subNby + 3) * 5, 0.0);

  //--------------------------------------------------------------------------
  //Calculate the contribution to the interpolation weight on the fluid.
  //--------------------------------------------------------------------------
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) { (*this.*interpolationweight)(i); }
  }

  //--------------------------------------------------------------------------
  //Communicate the force weight contributions which lie outside the local
  //processor sub domain.
  //--------------------------------------------------------------------------
  for (int i = 0; i < 10; i++) requests[i] = MPI_REQUEST_NULL;
  MPI_Isend(&Wf[0][0][0], 1, passxrho, comm->procneigh[0][0], 10, world, &requests[0]);
  MPI_Isend(&Wf[subNbx + 2][0][0], 1, passxrho, comm->procneigh[0][0], 20, world, &requests[1]);
  MPI_Isend(&Wf[subNbx - 1][0][0], 1, passxrho, comm->procneigh[0][1], 30, world, &requests[2]);
  MPI_Isend(&Wf[subNbx][0][0], 1, passxrho, comm->procneigh[0][1], 40, world, &requests[3]);
  MPI_Isend(&Wf[subNbx + 1][0][0], 1, passxrho, comm->procneigh[0][1], 50, world, &requests[4]);
  MPI_Irecv(&Wftempx[0][0][0], 1, passxWtemp, comm->procneigh[0][1], 10, world, &requests[5]);
  MPI_Irecv(&Wftempx[1][0][0], 1, passxWtemp, comm->procneigh[0][1], 20, world, &requests[6]);
  MPI_Irecv(&Wftempx[2][0][0], 1, passxWtemp, comm->procneigh[0][0], 30, world, &requests[7]);
  MPI_Irecv(&Wftempx[3][0][0], 1, passxWtemp, comm->procneigh[0][0], 40, world, &requests[8]);
  MPI_Irecv(&Wftempx[4][0][0], 1, passxWtemp, comm->procneigh[0][0], 50, world, &requests[9]);
  MPI_Waitall(10, requests, MPI_STATUS_IGNORE);

  for (int j = 0; j < subNby + 3; j++) {
    for (int k = 0; k < subNbz + 3; k++) {
      Wf[subNbx - 2][j][k] += Wftempx[0][j][k];
      Wf[subNbx - 3][j][k] += Wftempx[1][j][k];
      Wf[1][j][k] += Wftempx[2][j][k];
      Wf[2][j][k] += Wftempx[3][j][k];
      Wf[3][j][k] += Wftempx[4][j][k];
    }
  }

  for (int i = 0; i < 10; i++) requests[i] = MPI_REQUEST_NULL;
  MPI_Isend(&Wf[0][0][0], 1, passyrho, comm->procneigh[1][0], 10, world, &requests[0]);
  MPI_Isend(&Wf[0][subNby + 2][0], 1, passyrho, comm->procneigh[1][0], 20, world, &requests[1]);
  MPI_Isend(&Wf[0][subNby - 1][0], 1, passyrho, comm->procneigh[1][1], 30, world, &requests[2]);
  MPI_Isend(&Wf[0][subNby][0], 1, passyrho, comm->procneigh[1][1], 40, world, &requests[3]);
  MPI_Isend(&Wf[0][subNby + 1][0], 1, passyrho, comm->procneigh[1][1], 50, world, &requests[4]);
  MPI_Irecv(&Wftempy[0][0][0], 1, passyWtemp, comm->procneigh[1][1], 10, world, &requests[5]);
  MPI_Irecv(&Wftempy[0][1][0], 1, passyWtemp, comm->procneigh[1][1], 20, world, &requests[6]);
  MPI_Irecv(&Wftempy[0][2][0], 1, passyWtemp, comm->procneigh[1][0], 30, world, &requests[7]);
  MPI_Irecv(&Wftempy[0][3][0], 1, passyWtemp, comm->procneigh[1][0], 40, world, &requests[8]);
  MPI_Irecv(&Wftempy[0][4][0], 1, passyWtemp, comm->procneigh[1][0], 50, world, &requests[9]);
  MPI_Waitall(10, requests, MPI_STATUS_IGNORE);

  for (int i = 0; i < subNbx + 3; i++) {
    for (int k = 0; k < subNbz + 3; k++) {
      Wf[i][subNby - 2][k] += Wftempy[i][0][k];
      Wf[i][subNby - 3][k] += Wftempy[i][1][k];
      Wf[i][1][k] += Wftempy[i][2][k];
      Wf[i][2][k] += Wftempy[i][3][k];
      Wf[i][3][k] += Wftempy[i][4][k];
    }
  }

  for (int i = 0; i < 10; i++) requests[i] = MPI_REQUEST_NULL;
  MPI_Isend(&Wf[0][0][0], 1, passzrho, comm->procneigh[2][0], 10, world, &requests[0]);
  MPI_Isend(&Wf[0][0][subNbz + 2], 1, passzrho, comm->procneigh[2][0], 20, world, &requests[1]);
  MPI_Isend(&Wf[0][0][subNbz - 1], 1, passzrho, comm->procneigh[2][1], 30, world, &requests[2]);
  MPI_Isend(&Wf[0][0][subNbz], 1, passzrho, comm->procneigh[2][1], 40, world, &requests[3]);
  MPI_Isend(&Wf[0][0][subNbz + 1], 1, passzrho, comm->procneigh[2][1], 50, world, &requests[4]);
  MPI_Irecv(&Wftempz[0][0][0], 1, passzWtemp, comm->procneigh[2][1], 10, world, &requests[5]);
  MPI_Irecv(&Wftempz[0][0][1], 1, passzWtemp, comm->procneigh[2][1], 20, world, &requests[6]);
  MPI_Irecv(&Wftempz[0][0][2], 1, passzWtemp, comm->procneigh[2][0], 30, world, &requests[7]);
  MPI_Irecv(&Wftempz[0][0][3], 1, passzWtemp, comm->procneigh[2][0], 40, world, &requests[8]);
  MPI_Irecv(&Wftempz[0][0][4], 1, passzWtemp, comm->procneigh[2][0], 50, world, &requests[9]);
  MPI_Waitall(10, requests, MPI_STATUS_IGNORE);

  for (int i = 0; i < subNbx + 3; i++) {
    for (int j = 0; j < subNby + 3; j++) {
      Wf[i][j][subNbz - 2] += Wftempz[i][j][0];
      Wf[i][j][subNbz - 3] += Wftempz[i][j][1];
      Wf[i][j][1] += Wftempz[i][j][2];
      Wf[i][j][2] += Wftempz[i][j][3];
      Wf[i][j][3] += Wftempz[i][j][4];
    }
  }

  // Wf will potentially be needed on all ghost sites so pass back completed Wf
  int numrequests = 10;
  for (int i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
  MPI_Isend(&Wf[1][0][0], 1, passxrho, comm->procneigh[0][0], 30, world, &requests[0]);
  MPI_Isend(&Wf[2][0][0], 1, passxrho, comm->procneigh[0][0], 40, world, &requests[1]);
  MPI_Isend(&Wf[3][0][0], 1, passxrho, comm->procneigh[0][0], 50, world, &requests[2]);
  MPI_Isend(&Wf[subNbx - 3][0][0], 1, passxrho, comm->procneigh[0][1], 60, world, &requests[3]);
  MPI_Isend(&Wf[subNbx - 2][0][0], 1, passxrho, comm->procneigh[0][1], 70, world, &requests[4]);

  MPI_Irecv(&Wf[subNbx - 1][0][0], 1, passxrho, comm->procneigh[0][1], 30, world, &requests[5]);
  MPI_Irecv(&Wf[subNbx][0][0], 1, passxrho, comm->procneigh[0][1], 40, world, &requests[6]);
  MPI_Irecv(&Wf[subNbx + 1][0][0], 1, passxrho, comm->procneigh[0][1], 50, world, &requests[7]);
  MPI_Irecv(&Wf[subNbx + 2][0][0], 1, passxrho, comm->procneigh[0][0], 60, world, &requests[8]);
  MPI_Irecv(&Wf[0][0][0], 1, passxrho, comm->procneigh[0][0], 70, world, &requests[9]);
  MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

  for (int i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
  MPI_Isend(&Wf[0][1][0], 1, passyrho, comm->procneigh[1][0], 30, world, &requests[0]);
  MPI_Isend(&Wf[0][2][0], 1, passyrho, comm->procneigh[1][0], 40, world, &requests[1]);
  MPI_Isend(&Wf[0][3][0], 1, passyrho, comm->procneigh[1][0], 50, world, &requests[2]);
  MPI_Isend(&Wf[0][subNby - 3][0], 1, passyrho, comm->procneigh[1][1], 60, world, &requests[3]);
  MPI_Isend(&Wf[0][subNby - 2][0], 1, passyrho, comm->procneigh[1][1], 70, world, &requests[4]);

  MPI_Irecv(&Wf[0][subNby - 1][0], 1, passyrho, comm->procneigh[1][1], 30, world, &requests[5]);
  MPI_Irecv(&Wf[0][subNby][0], 1, passyrho, comm->procneigh[1][1], 40, world, &requests[6]);
  MPI_Irecv(&Wf[0][subNby + 1][0], 1, passyrho, comm->procneigh[1][1], 50, world, &requests[7]);
  MPI_Irecv(&Wf[0][subNby + 2][0], 1, passyrho, comm->procneigh[1][0], 60, world, &requests[8]);
  MPI_Irecv(&Wf[0][0][0], 1, passyrho, comm->procneigh[1][0], 70, world, &requests[9]);
  MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

  for (int i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
  int requestcount = 0;
  if (domain->periodicity[2] != 0 || comm->myloc[2] != 0) {
    MPI_Isend(&Wf[0][0][1], 1, passzrho, comm->procneigh[2][0], 30, world, &requests[requestcount]);
    MPI_Isend(&Wf[0][0][2], 1, passzrho, comm->procneigh[2][0], 40, world,
              &requests[requestcount + 1]);
    MPI_Isend(&Wf[0][0][3], 1, passzrho, comm->procneigh[2][0], 50, world,
              &requests[requestcount + 2]);

    MPI_Irecv(&Wf[0][0][subNbz + 2], 1, passzrho, comm->procneigh[2][0], 60, world,
              &requests[requestcount + 3]);
    MPI_Irecv(&Wf[0][0][0], 1, passzrho, comm->procneigh[2][0], 70, world,
              &requests[requestcount + 4]);
    requestcount = requestcount + 5;
  }
  if (domain->periodicity[2] != 0 || comm->myloc[2] != (comm->procgrid[2] - 1)) {
    MPI_Isend(&Wf[0][0][subNbz - 3], 1, passzrho, comm->procneigh[2][1], 60, world,
              &requests[requestcount]);
    MPI_Isend(&Wf[0][0][subNbz - 2], 1, passzrho, comm->procneigh[2][1], 70, world,
              &requests[requestcount + 1]);

    MPI_Irecv(&Wf[0][0][subNbz - 1], 1, passzrho, comm->procneigh[2][1], 30, world,
              &requests[requestcount + 2]);
    MPI_Irecv(&Wf[0][0][subNbz], 1, passzrho, comm->procneigh[2][1], 40, world,
              &requests[requestcount + 3]);
    MPI_Irecv(&Wf[0][0][subNbz + 1], 1, passzrho, comm->procneigh[2][1], 50, world,
              &requests[requestcount + 4]);
    requestcount = requestcount + 5;
  }
  MPI_Waitall(requestcount, requests, MPI_STATUS_IGNORE);
}

//==========================================================================
// Adds Keys cubic Hermite spline stencil weight for particle i to Wf array
//==========================================================================
void FixLbFluid::keys_interpolationweight(int i)
{
  double **x = atom->x;
  int ix, iy, iz;
  int ixp, iyp, izp;
  double dx1, dy1, dz1;
  double r, weightx, weighty, weightz;

  //--------------------------------------------------------------------------
  //Calculate nearest leftmost grid point.
  //Since array indices from 1 to subNb-2 correspond to the
  // local subprocessor domain (not indices from 0), use the
  // ceiling value.
  //--------------------------------------------------------------------------
  ix = (int) ceil((x[i][0] - domain->sublo[0]) / dx_lb);
  iy = (int) ceil((x[i][1] - domain->sublo[1]) / dx_lb);
  iz = (int) ceil((x[i][2] - domain->sublo[2]) / dx_lb);

  //The atom is allowed to be within one lattice grid point outside the
  //local processor sub-domain.
  if (ix < 0 || ix > (subNbx - 1) || iy < 0 || iy > (subNby - 1) || iz < 0 || iz > (subNbz - 1)) {
    error->one(FLERR,
               "Atom outside local processor simulation domain!"
               " Either unstable fluid pararmeters or require more frequent neighborlist rebuilds");
  }
  if (domain->periodicity[2] == 0 && comm->myloc[2] == 0 && iz < 2)
    error->warning(FLERR, "Atom too close to lower z wall.  Unphysical results may occur");
  if (domain->periodicity[2] == 0 && comm->myloc[2] == (comm->procgrid[2] - 1) &&
      (iz > (subNbz - 4)))
    error->warning(FLERR, "Atom too close to upper z wall.  Unphysical results may occur");

  //--------------------------------------------------------------------------
  //Calculate distances to the nearest points.
  //--------------------------------------------------------------------------
  dx1 = x[i][0] - (domain->sublo[0] + (ix - 1) * dx_lb);
  dy1 = x[i][1] - (domain->sublo[1] + (iy - 1) * dx_lb);
  dz1 = x[i][2] - (domain->sublo[2] + (iz - 1) * dx_lb);

  // Need to convert these to lattice units:
  dx1 = dx1 / dx_lb;
  dy1 = dy1 / dx_lb;
  dz1 = dz1 / dx_lb;

  //--------------------------------------------------------------------------
  // Calculate the interpolation weights
  //--------------------------------------------------------------------------
  for (int ii = -1; ii < 3; ii++) {
    r = fabs(-dx1 + ii);
    if (r >= 2)
      weightx = 0.0;
    else {
      if (r > 1) {
        weightx = 2.0 + (-4.0 + (2.5 - 0.5 * r) * r) * r;
      } else {
        weightx = 1.0 + r * r * (-2.5 + 1.5 * r);
      }
    }
    for (int jj = -1; jj < 3; jj++) {
      r = fabs(-dy1 + jj);
      if (r >= 2)
        weighty = 0.0;
      else {
        if (r > 1) {
          weighty = 2.0 + (-4.0 + (2.5 - 0.5 * r) * r) * r;
        } else {
          weighty = 1.0 + r * r * (-2.5 + 1.5 * r);
        }
      }
      for (int kk = -1; kk < 3; kk++) {
        r = fabs(-dz1 + kk);
        if (r >= 2)
          weightz = 0.0;
        else {
          if (r > 1) {
            weightz = 2.0 + (-4.0 + (2.5 - 0.5 * r) * r) * r;
          } else {
            weightz = 1.0 + r * r * (-2.5 + 1.5 * r);
          }
        }
        ixp = ix + ii;
        iyp = iy + jj;
        izp = iz + kk;

        if (ixp == -1) ixp = subNbx + 2;
        if (iyp == -1) iyp = subNby + 2;
        if (izp == -1) izp = subNbz + 2;

        Wf[ixp][iyp][izp] += fabs(weightx * weighty * weightz);
      }
    }
  }
}

//==========================================================================
// uses Keys cubic Hermite spline stencil to perform the velocity, density and
// force interpolations.
//==========================================================================
void FixLbFluid::keys_interpolation(int i, int uonly)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double massone;
  int ix, iy, iz;
  int ixp, iyp, izp;
  double dx1, dy1, dz1;
  int isten, ii, jj, kk;
  double r, weightx, weighty, weightz;
  double FfP[64], TfP[64];
  int k;
  double unew[3];
  double mnode;
  double gammavalue;

  //--------------------------------------------------------------------------
  //Calculate nearest leftmost grid point.
  //Since array indices from 1 to subNb-2 correspond to the
  // local subprocessor domain (not indices from 0), use the
  // ceiling value.
  //--------------------------------------------------------------------------
  ix = (int) ceil((x[i][0] - domain->sublo[0]) / dx_lb);
  iy = (int) ceil((x[i][1] - domain->sublo[1]) / dx_lb);
  iz = (int) ceil((x[i][2] - domain->sublo[2]) / dx_lb);

  //--------------------------------------------------------------------------
  //Calculate distances to the nearest points.
  //--------------------------------------------------------------------------
  dx1 = x[i][0] - (domain->sublo[0] + (ix - 1) * dx_lb);
  dy1 = x[i][1] - (domain->sublo[1] + (iy - 1) * dx_lb);
  dz1 = x[i][2] - (domain->sublo[2] + (iz - 1) * dx_lb);

  // Need to convert these to lattice units:
  dx1 = dx1 / dx_lb;
  dy1 = dy1 / dx_lb;
  dz1 = dz1 / dx_lb;

  unew[0] = 0.0;
  unew[1] = 0.0;
  unew[2] = 0.0;
  mnode = 0.0;
  double myweight = 0.0;
  isten = 0;
  double myweightsq = 0.0;

  //--------------------------------------------------------------------------
  // Calculate the interpolation weights, and interpolated values of
  // the fluid velocity, and density.
  //--------------------------------------------------------------------------
  for (int ii = -1; ii < 3; ii++) {
    r = fabs(-dx1 + ii);
    if (r >= 2)
      weightx = 0.0;
    else {
      if (r > 1) {
        weightx = 2.0 + (-4.0 + (2.5 - 0.5 * r) * r) * r;
      } else {
        weightx = 1.0 + r * r * (-2.5 + 1.5 * r);
      }
    }
    for (int jj = -1; jj < 3; jj++) {
      r = fabs(-dy1 + jj);
      if (r >= 2)
        weighty = 0.0;
      else {
        if (r > 1) {
          weighty = 2.0 + (-4.0 + (2.5 - 0.5 * r) * r) * r;
        } else {
          weighty = 1.0 + r * r * (-2.5 + 1.5 * r);
        }
      }
      for (int kk = -1; kk < 3; kk++) {
        r = fabs(-dz1 + kk);
        if (r >= 2)
          weightz = 0.0;
        else {
          if (r > 1) {
            weightz = 2.0 + (-4.0 + (2.5 - 0.5 * r) * r) * r;
          } else {
            weightz = 1.0 + r * r * (-2.5 + 1.5 * r);
          }
        }

        ixp = ix + ii;
        iyp = iy + jj;
        izp = iz + kk;

        //The atom is assumed to be within one lattice grid point outside the
        //local processor sub-domain.  This was checked when weights were constructed.

        if (ixp == -1) ixp = subNbx + 2;
        if (iyp == -1) iyp = subNby + 2;
        if (izp == -1) izp = subNbz + 2;

        FfP[isten] = weightx * weighty * weightz;
        if (fabs(FfP[isten]) > 1e-14)
          TfP[isten] = FfP[isten] * fabs(FfP[isten]) / Wf[ixp][iyp][izp];
        else
          TfP[isten] = 0.0;
        // interpolated mass and momentum.
        for (k = 0; k < 3; k++) {
          unew[k] += density_lb[ixp][iyp][izp] * u_lb[ixp][iyp][izp][k] * FfP[isten];
        }
        mnode += density_lb[ixp][iyp][izp] * FfP[isten];
        myweight += TfP[isten];
        myweightsq += TfP[isten] * FfP[isten];

        isten++;
      }
    }
  }

  //myweight=1.0;
  for (k = 0; k < 3; k++)    //mass weighted velocity
    unew[k] /= mnode;

  if (uonly) {
    // for(k=0; k<3; k++)
    //   unode[i][k] = unew[k]*dx_lb/dt_lb;
    // note this in LB units
    return;
  }

  dof_lb += myweight;
  mnode = mnode * myweight * myweight / myweightsq;

  massp[i] = mnode * dm_lb;    // convert to LAMMPS units for external use

  if (rmass)
    massone = rmass[i];
  else
    massone = mass[type[i]];
  massone = massone / dm_lb;

  // Compute the force on the fluid/particle.  Need to convert the velocity and f between
  // LAMMPS units and LB units.
  double hF[3], Ffull[3];    //, Fp[3];
  if (Gamma[type[i]] <
      0) {    // stationary or prescribed motion for particles, massone->infinity limit
    gammavalue = fabs(Gamma[type[i]]) * 2.0 * mnode;
    for (k = 0; k < 3; k++) {
      hF[k] = gammavalue * ((v[i][k] * dt_lb / dx_lb) - unew[k]);

      hydroF[i][k] = 0.0;    //force on particle in LAMMPS units

      //Full force on fluid:
      Ffull[k] = hF[k] / myweight;
    }
  } else {    // Normal case where particles respond to fluid forces
    gammavalue = Gamma[type[i]] * 2.0 * (mnode * massone) / (mnode + massone);
    for (k = 0; k < 3; k++) {
      hF[k] = gammavalue * ((v[i][k] * dt_lb / dx_lb) - unew[k]);

      hydroF[i][k] = -hF[k] * dm_lb * dx_lb / dt_lb / dt_lb;    //force on particle in LAMMPS units

      //Full force on fluid:
      Ffull[k] =
          (hF[k] + f[i][k] * dt_lb * dt_lb / dm_lb / dx_lb * mnode / (mnode + massone)) / myweight;
    }
  }

  // Distribute force on fluid to fluid mesh:
  isten = 0;
  for (ii = -1; ii < 3; ii++)
    for (jj = -1; jj < 3; jj++)
      for (kk = -1; kk < 3; kk++) {
        ixp = ix + ii;
        iyp = iy + jj;
        izp = iz + kk;

        if (ixp == -1) ixp = subNbx + 2;
        if (iyp == -1) iyp = subNby + 2;
        if (izp == -1) izp = subNbz + 2;

        for (k = 0; k < 3; k++) { Ff[ixp][iyp][izp][k] += Ffull[k] * TfP[isten]; }

        isten++;
      }
}

//==========================================================================
// Adds the Immersed Boundary 3-pt stencil weight for particle i to Wf array
//==========================================================================
void FixLbFluid::IBM3_interpolationweight(int i)
{
  double **x = atom->x;
  int ix, iy, iz;
  int ixp, iyp, izp;
  double dx1, dy1, dz1;
  double r, weightx, weighty, weightz;

  //--------------------------------------------------------------------------
  //Calculate nearest leftmost grid point.
  //Since array indices from 1 to subNb-2 correspond to the
  // local subprocessor domain (not indices from 0), use the
  // ceiling value.
  //--------------------------------------------------------------------------
  ix = (int) ceil((x[i][0] - domain->sublo[0]) / dx_lb);
  iy = (int) ceil((x[i][1] - domain->sublo[1]) / dx_lb);
  iz = (int) ceil((x[i][2] - domain->sublo[2]) / dx_lb);

  //--------------------------------------------------------------------------
  //Calculate distances to the nearest points.
  //--------------------------------------------------------------------------
  dx1 = x[i][0] - (domain->sublo[0] + (ix - 1) * dx_lb);
  dy1 = x[i][1] - (domain->sublo[1] + (iy - 1) * dx_lb);
  dz1 = x[i][2] - (domain->sublo[2] + (iz - 1) * dx_lb);

  // Need to convert these to lattice units:
  dx1 = dx1 / dx_lb;
  dy1 = dy1 / dx_lb;
  dz1 = dz1 / dx_lb;

  //The atom is allowed to be within one lattice grid point outside the
  //local processor sub-domain.  ??? Not sure these all make sense ???
  if (ix < 0 || ix > (subNbx - 1) || iy < 0 || iy > (subNby - 1) || iz < 0 || iz > (subNbz - 1)) {
    error->one(FLERR,
               "Atom outside local processor simulation domain!"
               " Either unstable fluid pararmeters or require more frequent neighborlist rebuilds");
  }
  if (domain->periodicity[2] == 0 && comm->myloc[2] == 0 && iz < 1)
    error->warning(FLERR, "Atom too close to lower z wall.  Unphysical results may occur");
  if (domain->periodicity[2] == 0 && comm->myloc[2] == (comm->procgrid[2] - 1) &&
      (iz > (subNbz - 2)))
    error->warning(FLERR, "Atom too close to upper z wall.  Unphysical results may occur");

  int imin = 0, jmin = 0, kmin = 0;
  if (dx1 < 0.5) imin = -1;
  if (dy1 < 0.5) jmin = -1;
  if (dz1 < 0.5) kmin = -1;

  //--------------------------------------------------------------------------
  // Calculate the interpolation weights
  //--------------------------------------------------------------------------
  for (int ii = imin; ii < 3 + imin; ii++) {
    r = fabs(-dx1 + ii);
    if (r >= 1.5)
      weightx = 0.0;
    else if (r > 0.5)
      weightx = (5.0 - 3.0 * r - sqrt(1.0 - 3.0 * (1.0 - r) * (1.0 - r))) / 6.0;
    else
      weightx = (1.0 + sqrt(1.0 - 3.0 * r * r)) / 3;

    for (int jj = jmin; jj < 3 + jmin; jj++) {
      r = fabs(-dy1 + jj);
      if (r >= 1.5)
        weighty = 0.0;
      else if (r > 0.5)
        weighty = (5.0 - 3.0 * r - sqrt(1.0 - 3.0 * (1.0 - r) * (1.0 - r))) / 6.0;
      else
        weighty = (1.0 + sqrt(1.0 - 3.0 * r * r)) / 3;

      for (int kk = kmin; kk < 3 + kmin; kk++) {
        r = fabs(-dz1 + kk);
        if (r >= 1.5)
          weightz = 0.0;
        else if (r > 0.5)
          weightz = (5.0 - 3.0 * r - sqrt(1.0 - 3.0 * (1.0 - r) * (1.0 - r))) / 6.0;
        else
          weightz = (1.0 + sqrt(1.0 - 3.0 * r * r)) / 3;

        ixp = ix + ii;
        iyp = iy + jj;
        izp = iz + kk;

        if (ixp == -1) ixp = subNbx + 2;
        if (iyp == -1) iyp = subNby + 2;
        if (izp == -1) izp = subNbz + 2;

        Wf[ixp][iyp][izp] += weightx * weighty * weightz;
      }
    }
  }
}

//==========================================================================
// uses the Immersed Boundary 3-point stencil to perform the velocity, density and
// force interpolations.
//==========================================================================
void FixLbFluid::IBM3_interpolation(int i, int uonly)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double massone;
  int ix, iy, iz;
  int ixp, iyp, izp;
  double dx1, dy1, dz1;
  int isten, ii, jj, kk;
  double r, weightx, weighty, weightz;
  double FfP[27], TfP[27];
  int k;
  double unew[3];
  double mnode;
  double gammavalue;

  //--------------------------------------------------------------------------
  //Calculate nearest leftmost grid point.
  //Since array indices from 1 to subNb-2 correspond to the
  // local subprocessor domain (not indices from 0), use the
  // ceiling value.
  //--------------------------------------------------------------------------
  ix = (int) ceil((x[i][0] - domain->sublo[0]) / dx_lb);
  iy = (int) ceil((x[i][1] - domain->sublo[1]) / dx_lb);
  iz = (int) ceil((x[i][2] - domain->sublo[2]) / dx_lb);

  //--------------------------------------------------------------------------
  //Calculate distances to the nearest points.
  //--------------------------------------------------------------------------
  dx1 = x[i][0] - (domain->sublo[0] + (ix - 1) * dx_lb);
  dy1 = x[i][1] - (domain->sublo[1] + (iy - 1) * dx_lb);
  dz1 = x[i][2] - (domain->sublo[2] + (iz - 1) * dx_lb);

  // Need to convert these to lattice units:
  dx1 = dx1 / dx_lb;
  dy1 = dy1 / dx_lb;
  dz1 = dz1 / dx_lb;

  unew[0] = 0.0;
  unew[1] = 0.0;
  unew[2] = 0.0;
  mnode = 0.0;
  double myweight = 0.0;
  isten = 0;
  double myweightsq = 0.0;

  //--------------------------------------------------------------------------
  // Calculate the interpolation weights, and interpolated values of
  // the fluid velocity, and density.
  //--------------------------------------------------------------------------
  int imin = 0, jmin = 0, kmin = 0;
  if (dx1 < 0.5) imin = -1;
  if (dy1 < 0.5) jmin = -1;
  if (dz1 < 0.5) kmin = -1;

  for (int ii = imin; ii < 3 + imin; ii++) {
    r = fabs(-dx1 + ii);
    if (r >= 1.5)
      weightx = 0.0;
    else if (r > 0.5)
      weightx = (5.0 - 3.0 * r - sqrt(1.0 - 3.0 * (1.0 - r) * (1.0 - r))) / 6.0;
    else
      weightx = (1.0 + sqrt(1.0 - 3.0 * r * r)) / 3;

    for (int jj = jmin; jj < 3 + jmin; jj++) {
      r = fabs(-dy1 + jj);
      if (r >= 1.5)
        weighty = 0.0;
      else if (r > 0.5)
        weighty = (5.0 - 3.0 * r - sqrt(1.0 - 3.0 * (1.0 - r) * (1.0 - r))) / 6.0;
      else
        weighty = (1.0 + sqrt(1.0 - 3.0 * r * r)) / 3;

      for (int kk = kmin; kk < 3 + kmin; kk++) {
        r = fabs(-dz1 + kk);
        if (r >= 1.5)
          weightz = 0.0;
        else if (r > 0.5)
          weightz = (5.0 - 3.0 * r - sqrt(1.0 - 3.0 * (1.0 - r) * (1.0 - r))) / 6.0;
        else
          weightz = (1.0 + sqrt(1.0 - 3.0 * r * r)) / 3;

        ixp = ix + ii;
        iyp = iy + jj;
        izp = iz + kk;

        //The atom is assumed to be within one lattice grid point outside the
        //local processor sub-domain.  This was checked when weights were constructed.

        if (ixp == -1) ixp = subNbx + 2;
        if (iyp == -1) iyp = subNby + 2;
        if (izp == -1) izp = subNbz + 2;

        FfP[isten] = weightx * weighty * weightz;

        if (fabs(FfP[isten]) > 1e-14)
          TfP[isten] = FfP[isten] * FfP[isten] / Wf[ixp][iyp][izp];
        else
          TfP[isten] = 0.0;
        // interpolated mass and momentum.
        for (k = 0; k < 3; k++) {
          unew[k] += density_lb[ixp][iyp][izp] * u_lb[ixp][iyp][izp][k] * FfP[isten];
        }
        mnode += density_lb[ixp][iyp][izp] * FfP[isten];
        myweight += TfP[isten];
        myweightsq += TfP[isten] * FfP[isten];

        isten++;
      }
    }
  }

  //myweight=1.0;
  for (k = 0; k < 3; k++)    //mass weighted velocity
    unew[k] /= mnode;

  if (uonly) {
    // for(k=0; k<3; k++)
    //   unode[i][k] = unew[k]*dx_lb/dt_lb;
    return;
  }
  dof_lb += myweight;
  mnode = mnode * myweight * myweight / myweightsq;

  massp[i] = mnode * dm_lb;    // convert to LAMMPS units for external use

  if (rmass)
    massone = rmass[i];
  else
    massone = mass[type[i]];
  massone = massone / dm_lb;

  // Compute the force on the fluid/particle.  Need to convert the velocity and f between
  // LAMMPS units and LB units.
  double hF[3], Ffull[3];    //, Fp[3];
  if (Gamma[type[i]] <
      0) {    // stationary or prescribed motion for particles, massone->infinity limit
    gammavalue = fabs(Gamma[type[i]]) * 2.0 * mnode;
    for (k = 0; k < 3; k++) {
      hF[k] = gammavalue * ((v[i][k] * dt_lb / dx_lb) - unew[k]);

      hydroF[i][k] = 0.0;    //force on particle in LAMMPS units

      //Full force on fluid:
      Ffull[k] = hF[k] / myweight;
    }
  } else {    // Normal case where particles respond to fluid forces
    gammavalue = Gamma[type[i]] * 2.0 * (mnode * massone) / (mnode + massone);
    for (k = 0; k < 3; k++) {
      hF[k] = gammavalue * ((v[i][k] * dt_lb / dx_lb) - unew[k]);

      hydroF[i][k] = -hF[k] * dm_lb * dx_lb / dt_lb / dt_lb;    //force on particle in LAMMPS units

      //Full force on fluid:
      Ffull[k] =
          (hF[k] + f[i][k] * dt_lb * dt_lb / dm_lb / dx_lb * mnode / (mnode + massone)) / myweight;
    }
  }

  // Distribute force on fluid to fluid mesh:
  isten = 0;
  for (ii = imin; ii < 3 + imin; ii++)
    for (jj = jmin; jj < 3 + jmin; jj++)
      for (kk = kmin; kk < 3 + kmin; kk++) {
        ixp = ix + ii;
        iyp = iy + jj;
        izp = iz + kk;

        if (ixp == -1) ixp = subNbx + 2;
        if (iyp == -1) iyp = subNby + 2;
        if (izp == -1) izp = subNbz + 2;

        for (k = 0; k < 3; k++) { Ff[ixp][iyp][izp][k] += Ffull[k] * TfP[isten]; }

        isten++;
      }
}

//==========================================================================
// uses the trilinear stencil to perform the velocity, density and
// force interpolations.
//==========================================================================
void FixLbFluid::trilinear_interpolationweight(int i)
{
  double **x = atom->x;
  int ix, iy, iz;
  int ixp, iyp, izp;
  double dx1, dy1, dz1;
  double FfP[8];

  //--------------------------------------------------------------------------
  // Calculate nearest leftmost grid point.
  // Since array indices from 1 to subNb-2 correspond to the
  // local subprocessor domain (not indices from 0), use the
  // ceiling value.
  //--------------------------------------------------------------------------
  ix = (int) ceil((x[i][0] - domain->sublo[0]) / dx_lb);
  iy = (int) ceil((x[i][1] - domain->sublo[1]) / dx_lb);
  iz = (int) ceil((x[i][2] - domain->sublo[2]) / dx_lb);

  //--------------------------------------------------------------------------
  //Calculate distances to the nearest points.
  //--------------------------------------------------------------------------
  dx1 = x[i][0] - (domain->sublo[0] + (ix - 1) * dx_lb);
  dy1 = x[i][1] - (domain->sublo[1] + (iy - 1) * dx_lb);
  dz1 = x[i][2] - (domain->sublo[2] + (iz - 1) * dx_lb);

  //--------------------------------------------------------------------------
  // Need to convert these to lattice units:
  //--------------------------------------------------------------------------
  dx1 = dx1 / dx_lb;
  dy1 = dy1 / dx_lb;
  dz1 = dz1 / dx_lb;

  //--------------------------------------------------------------------------
  // Calculate the interpolation weights
  //--------------------------------------------------------------------------
  FfP[0] = (1.0 - dx1) * (1.0 - dy1) * (1.0 - dz1);
  FfP[1] = (1.0 - dx1) * (1.0 - dy1) * dz1;
  FfP[2] = (1.0 - dx1) * dy1 * (1.0 - dz1);
  FfP[3] = (1.0 - dx1) * dy1 * dz1;
  FfP[4] = dx1 * (1.0 - dy1) * (1.0 - dz1);
  FfP[5] = dx1 * (1.0 - dy1) * dz1;
  FfP[6] = dx1 * dy1 * (1.0 - dz1);
  FfP[7] = dx1 * dy1 * dz1;

  ixp = (ix + 1);
  iyp = (iy + 1);
  izp = (iz + 1);

  //The atom is allowed to be within one lattice grid point outside the
  //local processor sub-domain.
  if (ix < 0 || ixp > (subNbx + 1) || iy < 0 || iyp > (subNby + 1) || iz < 0 || izp > (subNbz + 1))
    error->one(
        FLERR,
        "Atom outside local processor simulation domain.  Either unstable fluid pararmeters, or \
require more frequent neighborlist rebuilds");

  if (domain->periodicity[2] == 0 && comm->myloc[2] == 0 && (iz < 1 || izp < 1))
    error->warning(FLERR, "Atom too close to lower z wall.  Unphysical results may occur");
  if (domain->periodicity[2] == 0 && comm->myloc[2] == (comm->procgrid[2] - 1) &&
      (izp > (subNbz - 2) || iz > (subNbz - 2)))
    error->warning(FLERR, "Atom too close to upper z wall.  Unphysical results may occur");

  Wf[ix][iy][iz] += FfP[0];
  Wf[ix][iy][izp] += FfP[1];
  Wf[ix][iyp][iz] += FfP[2];
  Wf[ix][iyp][izp] += FfP[3];
  Wf[ixp][iy][iz] += FfP[4];
  Wf[ixp][iy][izp] += FfP[5];
  Wf[ixp][iyp][iz] += FfP[6];
  Wf[ixp][iyp][izp] += FfP[7];
}

void FixLbFluid::trilinear_interpolation(int i, int uonly)
{    // usual call had uonly=0
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double massone;
  int ix, iy, iz;
  int ixp, iyp, izp;
  double dx1, dy1, dz1;
  double FfP[8];
  int k;
  double unew[3];
  double mnode;
  double gammavalue;

  //--------------------------------------------------------------------------
  // Calculate nearest leftmost grid point.
  // Since array indices from 1 to subNb-2 correspond to the
  // local subprocessor domain (not indices from 0), use the
  // ceiling value.
  //--------------------------------------------------------------------------
  ix = (int) ceil((x[i][0] - domain->sublo[0]) / dx_lb);
  iy = (int) ceil((x[i][1] - domain->sublo[1]) / dx_lb);
  iz = (int) ceil((x[i][2] - domain->sublo[2]) / dx_lb);

  //--------------------------------------------------------------------------
  //Calculate distances to the nearest points.
  //--------------------------------------------------------------------------
  dx1 = x[i][0] - (domain->sublo[0] + (ix - 1) * dx_lb);
  dy1 = x[i][1] - (domain->sublo[1] + (iy - 1) * dx_lb);
  dz1 = x[i][2] - (domain->sublo[2] + (iz - 1) * dx_lb);

  //--------------------------------------------------------------------------
  // Need to convert these to lattice units:
  //--------------------------------------------------------------------------
  dx1 = dx1 / dx_lb;
  dy1 = dy1 / dx_lb;
  dz1 = dz1 / dx_lb;

  //--------------------------------------------------------------------------
  // Calculate the interpolation weights
  //--------------------------------------------------------------------------
  FfP[0] = (1.0 - dx1) * (1.0 - dy1) * (1.0 - dz1);
  FfP[1] = (1.0 - dx1) * (1.0 - dy1) * dz1;
  FfP[2] = (1.0 - dx1) * dy1 * (1.0 - dz1);
  FfP[3] = (1.0 - dx1) * dy1 * dz1;
  FfP[4] = dx1 * (1.0 - dy1) * (1.0 - dz1);
  FfP[5] = dx1 * (1.0 - dy1) * dz1;
  FfP[6] = dx1 * dy1 * (1.0 - dz1);
  FfP[7] = dx1 * dy1 * dz1;

  ixp = (ix + 1);
  iyp = (iy + 1);
  izp = (iz + 1);

  mnode = density_lb[ix][iy][iz] * FfP[0] + density_lb[ix][iy][izp] * FfP[1] +
      density_lb[ix][iyp][iz] * FfP[2] + density_lb[ix][iyp][izp] * FfP[3] +
      density_lb[ixp][iy][iz] * FfP[4] + density_lb[ixp][iy][izp] * FfP[5] +
      density_lb[ixp][iyp][iz] * FfP[6] + density_lb[ixp][iyp][izp] * FfP[7];

  for (k = 0; k < 3; k++) {    // tri-linearly interpolated mass-weighted velocity at node
    unew[k] = (density_lb[ix][iy][iz] * u_lb[ix][iy][iz][k] * FfP[0] +
               density_lb[ix][iy][izp] * u_lb[ix][iy][izp][k] * FfP[1] +
               density_lb[ix][iyp][iz] * u_lb[ix][iyp][iz][k] * FfP[2] +
               density_lb[ix][iyp][izp] * u_lb[ix][iyp][izp][k] * FfP[3] +
               density_lb[ixp][iy][iz] * u_lb[ixp][iy][iz][k] * FfP[4] +
               density_lb[ixp][iy][izp] * u_lb[ixp][iy][izp][k] * FfP[5] +
               density_lb[ixp][iyp][iz] * u_lb[ixp][iyp][iz][k] * FfP[6] +
               density_lb[ixp][iyp][izp] * u_lb[ixp][iyp][izp][k] * FfP[7]) /
        mnode;
  }

  if (uonly) {
    // for(k=0; k<3; k++)
    //   unode[i][k] = unew[k]*dx_lb/dt_lb;
    return;
  }

  double TfP[8] = {0};
  if (FfP[0] > 0) TfP[0] = FfP[0] * FfP[0] / Wf[ix][iy][iz];
  if (FfP[1] > 0) TfP[1] = FfP[1] * FfP[1] / Wf[ix][iy][izp];
  if (FfP[2] > 0) TfP[2] = FfP[2] * FfP[2] / Wf[ix][iyp][iz];
  if (FfP[3] > 0) TfP[3] = FfP[3] * FfP[3] / Wf[ix][iyp][izp];
  if (FfP[4] > 0) TfP[4] = FfP[4] * FfP[4] / Wf[ixp][iy][iz];
  if (FfP[5] > 0) TfP[5] = FfP[5] * FfP[5] / Wf[ixp][iy][izp];
  if (FfP[6] > 0) TfP[6] = FfP[6] * FfP[6] / Wf[ixp][iyp][iz];
  if (FfP[7] > 0) TfP[7] = FfP[7] * FfP[7] / Wf[ixp][iyp][izp];

  double myweight = TfP[0] + TfP[1] + TfP[2] + TfP[3] + TfP[4] + TfP[5] + TfP[6] + TfP[7];
  double myweightsq = TfP[0] * FfP[0] + TfP[1] * FfP[1] + TfP[2] * FfP[2] + TfP[3] * FfP[3] +
      TfP[4] * FfP[4] + TfP[5] * FfP[5] + TfP[6] * FfP[6] + TfP[7] * FfP[7];

  dof_lb += myweight;
  mnode = mnode * myweight * myweight / myweightsq;

  massp[i] = mnode * dm_lb;    // convert to LAMMPS units for external use

  if (rmass)
    massone = rmass[i];
  else
    massone = mass[type[i]];
  massone = massone / dm_lb;    //convert to lattice units

  if (Gamma[type[i]] < 0)    // node has prescribed motion, equivalent to massone->infinity
    gammavalue = fabs(Gamma[type[i]]) * 2.0 * mnode;
  else    // normal case
    gammavalue = Gamma[type[i]] * 2.0 * (mnode * massone) / (mnode + massone);

  // Compute The force on the fluid/particle.  Need to convert the velocity and f from
  // LAMMPS units to LB units.
  for (k = 0; k < 3; k++) {
    double hF = gammavalue * ((v[i][k] * dt_lb / dx_lb) - unew[k]);

    double Fp;
    if (Gamma[type[i]] < 0) {    // massone->infinity means no effect on particles, i.e. F->0
      hydroF[i][k] = 0.0;
      Fp = 0;
    } else {    // normal case
      hydroF[i][k] = -hF * dm_lb * dx_lb / dt_lb / dt_lb;
      Fp = f[i][k] * dt_lb * dt_lb / dm_lb / dx_lb * mnode / (mnode + massone);
    }
    double Ffull = (hF + Fp) / myweight;

    Ff[ix][iy][iz][k] += Ffull * TfP[0];
    Ff[ix][iy][izp][k] += Ffull * TfP[1];
    Ff[ix][iyp][iz][k] += Ffull * TfP[2];
    Ff[ix][iyp][izp][k] += Ffull * TfP[3];
    Ff[ixp][iy][iz][k] += Ffull * TfP[4];
    Ff[ixp][iy][izp][k] += Ffull * TfP[5];
    Ff[ixp][iyp][iz][k] += Ffull * TfP[6];
    Ff[ixp][iyp][izp][k] += Ffull * TfP[7];
  }
}

/// Create MPI_Datatype for writing to the local buffer portion of a global buffer with specified number of local
/// and global ghost points (this is the local buffer embedded in the global buffer).
///
/// This routine creates a MPI_Datatype for a global buffer (e.g., a file) from which the local buffer plus the
/// global ghost points (which may be less than the actual number of local ghost points) are to be written.  The
/// MPI_Datatype is only suitable for writing as the the local buffers do not overlap in the global buffer.  Use
/// mpiTypeGlobalRead() to create a MPI_Datatype suitable for reading.
///
/// \param[in] local_ghost   Ghost points per side in local buffer (>= global_ghost)
/// \param[in] local_size    Size of local buffer
/// \param[in] global_offset  Offset of local buffer in global buffer excluding ghost points
/// \param[in] global_ghost  Ghost points per side in local buffer (<= local_ghost)
/// \param[in] global_size   Size of global buffer excluding ghost points
/// \param[in] mpitype       MPI_Datatype of an element of the buffers
///
/// The following is a 2D slice taken from the centre of a buffer showing a local buffer, denoted by o and O,
/// emedded in a global buffer, denoted by X.  Areas outside the lines are ghost points.  The capital vs small
/// letters denote regions incuded and excluded, respectively, in the created type.  Note that lower case x is not
/// part of the global buffer while lower-case o is part of the local buffer.
///
///     oooooooxx                                xxxxxxxxx
///     oOOOoooXx  local_ghost  = 2              xXXXXXXXx  local_ghost  = 2
///     oO+-+-+Xx  global_ghost = 1              xX+---+oo  global_ghost = 1
///     oO|O|o|Xx                                xX|ooo|oo
///     oo+-+o|Xx  global_offset = { 2, 2, 2 }   xX|o+-+Oo  global_offset = { 3, 3, 3 }
///     oo|ooo|Xx                                xX|o|O|Oo
///     oo+---+Xx  local_size  = { 3, 3, 3 }     xX+-+-+Oo  local_size  = { 3, 3, 3 }
///     xXXXXXXXx  global_size = { 5, 5, 5 }     xXooOOOOo  global_size = { 5, 5, 5 }
///     xxxxxxxxx                                xxooooooo
///
static MPI_Datatype mpiTypeGlobalWrite(const int /*local_ghost*/, const int *local_size,
                                       const int *global_offset, const int global_ghost,
                                       const int *global_size, const MPI_Datatype mpitype)
{
  MPI_Datatype global_mpitype;

  {
    bool endpoint_lower[] = {global_offset[0] == 0, global_offset[1] == 0, global_offset[2] == 0};
    bool endpoint_upper[] = {global_offset[0] + local_size[0] == global_size[0],
                             global_offset[1] + local_size[1] == global_size[1],
                             global_offset[2] + local_size[2] == global_size[2]};

    int sizes[] = {global_ghost + global_size[0] + global_ghost,
                   global_ghost + global_size[1] + global_ghost,
                   global_ghost + global_size[2] + global_ghost};
    int subsizes[] = {(global_ghost * endpoint_lower[0] + local_size[0] - 1 * !endpoint_upper[0] +
                       global_ghost * endpoint_upper[0]),
                      (global_ghost * endpoint_lower[1] + local_size[1] - 1 * !endpoint_upper[1] +
                       global_ghost * endpoint_upper[1]),
                      (global_ghost * endpoint_lower[2] + local_size[2] - 1 * !endpoint_upper[2] +
                       global_ghost * endpoint_upper[2])};
    int starts[] = {global_ghost * !endpoint_lower[0] + global_offset[0],
                    global_ghost * !endpoint_lower[1] + global_offset[1],
                    global_ghost * !endpoint_lower[2] + global_offset[2]};
    // Note Fortran ordering as we switch order for paraview output
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, mpitype,
                             &global_mpitype);
  }

  return global_mpitype;
}

/// Create MPI_Datatype for writing to the local buffer portion of a global buffer with specified number of local
/// and global ghost points (this is the local buffer appart from the global buffer).
///
/// This routine creates a MPI_Datatype for a local buffer to which a porition of a global buffer plus global ghost
/// points (which may be less than the actual number of local ghost points) are to be written.  The MPI_Datatype is
/// only suitable for writing as the the local buffers do not overlap in the global buffer.  Use mpiTypeLocalRead()
/// to create a MPI_Datatype suitable for writing.
///
/// \param[in] local_ghost   Ghost points per side in local buffer (>= global_ghost)
/// \param[in] local_size    Size of local buffer
/// \param[in] global_offset Offset of local buffer in global buffer excluding ghost points
/// \param[in] global_ghost  Ghost points per side in local buffer (<= local_ghost)
/// \param[in] global_size   Size of global buffer excluding ghost points
/// \param[in] mpitype       MPI_Datatype of an element of the buffer
///
/// Reference the diagram under mpiTypeGlobalRead().
///
static MPI_Datatype mpiTypeLocalWrite(const int local_ghost, const int *local_size,
                                      const int *global_offset, const int global_ghost,
                                      const int *global_size, const MPI_Datatype mpitype)
{
  MPI_Datatype local_mpitype;

  {
    bool endpoint_lower[] = {global_offset[0] == 0, global_offset[1] == 0, global_offset[2] == 0};
    bool endpoint_upper[] = {global_offset[0] + local_size[0] == global_size[0],
                             global_offset[1] + local_size[1] == global_size[1],
                             global_offset[2] + local_size[2] == global_size[2]};

    int sizes[] = {local_ghost + local_size[0] + local_ghost,
                   local_ghost + local_size[1] + local_ghost,
                   local_ghost + local_size[2] + local_ghost};
    int subsizes[] = {(global_ghost * endpoint_lower[0] + local_size[0] - 1 * !endpoint_upper[0] +
                       global_ghost * endpoint_upper[0]),
                      (global_ghost * endpoint_lower[1] + local_size[1] - 1 * !endpoint_upper[1] +
                       global_ghost * endpoint_upper[1]),
                      (global_ghost * endpoint_lower[2] + local_size[2] - 1 * !endpoint_upper[2] +
                       global_ghost * endpoint_upper[2])};
    // Note: odd shift by -1 as bottom is actually wrapped to top in MPI version of package
    int starts[] = {local_ghost - 1 - global_ghost * endpoint_lower[0],
                    local_ghost - 1 - global_ghost * endpoint_lower[1],
                    local_ghost - 1 - global_ghost * endpoint_lower[2]};
    // Note Fortran ordering as we switch order for paraview output
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, mpitype,
                             &local_mpitype);
  }

  return local_mpitype;
}

/// Create MPI_Datatype for writing to the local portion of the global dump buffer.
///
/// \param[in] local_size    Size of local buffer
/// \param[in] global_offset Offset of local buffer in global buffer excluding ghost points
/// \param[in] global_size   Size of global buffer excluding ghost points
///
static MPI_Datatype mpiTypeDumpGlobal(const int *local_size, const int *global_offset,
                                      const int *global_size)
{
  MPI_Datatype dump;

  // Global MPI types for our porition of the global dump file
  {
    MPI_Datatype realType3_mpitype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &realType3_mpitype);

    // Density and velocity types for our chunk of the global file
    MPI_Datatype density_mpitype =
        mpiTypeGlobalWrite(2, local_size, global_offset, 0, global_size, MPI_DOUBLE);
    MPI_Datatype velocity_mpitype =
        mpiTypeGlobalWrite(2, local_size, global_offset, 0, global_size, realType3_mpitype);

    // Density followed by velocity type for our chunk of the global dump file
    {
      MPI_Aint density_lb, density_extent;
      MPI_Type_get_extent(density_mpitype, &density_lb, &density_extent);

      int blocklengths[] = {1, 1};
      MPI_Aint displacements[] = {0, density_lb + density_extent};
      MPI_Datatype datatypes[] = {density_mpitype, velocity_mpitype};

      MPI_Type_create_struct(2, blocklengths, displacements, datatypes, &dump);
    }

    // Release local use types
    MPI_Type_free(&velocity_mpitype);
    MPI_Type_free(&density_mpitype);
    MPI_Type_free(&realType3_mpitype);
  }

  return dump;
}

// ------------------------------------------------------------------------------------------------------------- //
//==========================================================================
// Creates various buffers for output.
//==========================================================================
void FixLbFluid::SetupBuffers()
{
  fluid_global_n0[0] = Nbx + 1;    // Both domain boundaries -- this unifies the seperate periodic
  fluid_global_n0[1] = Nby + 1;    // and non-periodic sizes of the original code
  fluid_global_n0[2] = Nbz + 1 - (domain->periodicity[2] == 0);

  fluid_local_n0[0] = subNbx - 1;    // These include both domain boundaries as with fluid_global_n0
  fluid_local_n0[1] = subNby - 1;    // but no additional points beyond these
  fluid_local_n0[2] =
      subNbz - 1 - (domain->periodicity[2] == 0 && (comm->myloc[2] == comm->procgrid[2] - 1));

  fluid_global_o0[0] = (fluid_local_n0[0] - 1) * comm->myloc[0];
  fluid_global_o0[1] = (fluid_local_n0[1] - 1) * comm->myloc[1];
  fluid_global_o0[2] = (fluid_local_n0[2] - 1) * comm->myloc[2];

  // Local write MPI types for our portion of the global dump file
  fluid_density_2_mpitype =
      mpiTypeLocalWrite(2, fluid_local_n0, fluid_global_o0, 0, fluid_global_n0, MPI_DOUBLE);
  fluid_velocity_2_mpitype =
      mpiTypeLocalWrite(2, fluid_local_n0, fluid_global_o0, 0, fluid_global_n0, realType3_mpitype);

  MPI_Type_commit(&fluid_density_2_mpitype);
  MPI_Type_commit(&fluid_velocity_2_mpitype);

  // Global write MPI types for our porition of the global dump file
  dump_file_mpitype = mpiTypeDumpGlobal(fluid_local_n0, fluid_global_o0, fluid_global_n0);
  MPI_Type_commit(&dump_file_mpitype);

  // Output
  if (dump_interval) {
    if (me == 0) {
      dump_file_handle_xdmf = fopen(dump_file_name_xdmf.c_str(), "w");
      if (!dump_file_handle_xdmf)
        error->one(FLERR, "Unable to truncate/create \"{}\": {}", dump_file_name_xdmf,
                   utils::getsyserror());
      fprintf(dump_file_handle_xdmf,
              "<?xml version=\"1.0\" ?>\n"
              "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
              "<Xdmf Version=\"2.0\">\n"
              "  <Domain>\n"
              "    <Grid Name=\"fluid\" GridType=\"Collection\" CollectionType=\"Temporal\">\n\n");
    }

    MPI_File_open(world, const_cast<char *>(dump_file_name_raw.c_str()),
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &dump_file_handle_raw);

    MPI_File_set_size(dump_file_handle_raw, 0);
    MPI_File_set_view(dump_file_handle_raw, 0, MPI_DOUBLE, dump_file_mpitype, "native",
                      MPI_INFO_NULL);
  }
}

//==========================================================================
// Output fluid density and velocity to file in XDMF format
//==========================================================================
void FixLbFluid::dump(const bigint step)
{
  static bigint frameindex = 0;

  if (dump_interval && step % dump_interval == 0) {
    // Write XDMF grid entry for time step
    if (me == 0) {
      bigint block =
          (bigint) fluid_global_n0[0] * fluid_global_n0[1] * fluid_global_n0[2] * sizeof(double);
      bigint offset = frameindex * block * (1 + 3);
      double time = dump_time_index ? update->ntimestep * dt_lb : frameindex;

      fmt::print(dump_file_handle_xdmf,
                 "      <Grid Name=\"{}\">\n"
                 "        <Time Value=\"{:f}\"/>\n\n"
                 "        <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"{} {} {}\"/>\n"
                 "        <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"
                 "          <DataItem Dimensions=\"3\">\n"
                 "            {:f} {:f} {:f}\n"
                 "          </DataItem>\n"
                 "          <DataItem Dimensions=\"3\">\n"
                 "            {:f} {:f} {:f}\n"
                 "          </DataItem>\n"
                 "        </Geometry>\n\n",
                 step, time, fluid_global_n0[2], fluid_global_n0[1], fluid_global_n0[0],
                 domain->boxlo[2], domain->boxlo[1], domain->boxlo[0], dx_lb, dx_lb, dx_lb);
      fmt::print(dump_file_handle_xdmf,
                 "        <Attribute Name=\"density\">\n"
                 "          <DataItem ItemType=\"Function\" Function=\"$0 * {:f}\" "
                 "Dimensions=\"{} {} {}\">\n"
                 "            <DataItem Precision=\"{}\" Format=\"Binary\" Seek=\"{}\" "
                 "Dimensions=\"{} {} {}\">\n"
                 "              {}\n"
                 "            </DataItem>\n"
                 "          </DataItem>\n"
                 "        </Attribute>\n\n",
                 dm_lb / (dx_lb * dx_lb * dx_lb), fluid_global_n0[2], fluid_global_n0[1],
                 fluid_global_n0[0], sizeof(double), offset, fluid_global_n0[2], fluid_global_n0[1],
                 fluid_global_n0[0], dump_file_name_raw.c_str());
      fmt::print(dump_file_handle_xdmf,
                 "        <Attribute Name=\"velocity\" AttributeType=\"Vector\">\n"
                 "          <DataItem ItemType=\"Function\" Function=\"$0 * {:f}\" "
                 "Dimensions=\"{} {} {} 3\">\n"
                 "            <DataItem Precision=\"{}\" Format=\"Binary\" Seek=\"{}\" "
                 "Dimensions=\"{} {} {} 3\">\n"
                 "              {}\n"
                 "            </DataItem>\n"
                 "          </DataItem>\n"
                 "        </Attribute>\n\n",
                 dx_lb / dt_lb, fluid_global_n0[2], fluid_global_n0[1], fluid_global_n0[0],
                 sizeof(double), offset + block * 1, fluid_global_n0[2], fluid_global_n0[1],
                 fluid_global_n0[0], dump_file_name_raw.c_str());
      fmt::print(dump_file_handle_xdmf, "      </Grid>\n\n");

      frameindex++;
    }

    // Write raw data
    {
      const size_t size2 = (subNbx + 3) * (subNby + 3) * (subNbz + 3);

      // Transpose local arrays to fortran-order for paraview output
      std::vector<double> density_2_fort(size2);
      std::vector<double> velocity_2_fort(size2 * 3);
      int indexc = 0;
      for (int i = 0; i < subNbx + 3; i++)
        for (int j = 0; j < subNby + 3; j++)
          for (int k = 0; k < subNbz + 3; k++) {
            density_2_fort[i + (subNbx + 3) * (j + (subNby + 3) * k)] = density_lb[i][j][k];
            velocity_2_fort[0 + 3 * (i + (subNbx + 3) * (j + (subNby + 3) * k))] = u_lb[i][j][k][0];
            velocity_2_fort[1 + 3 * (i + (subNbx + 3) * (j + (subNby + 3) * k))] = u_lb[i][j][k][1];
            velocity_2_fort[2 + 3 * (i + (subNbx + 3) * (j + (subNby + 3) * k))] = u_lb[i][j][k][2];
            indexc++;
          }

      MPI_File_write_all(dump_file_handle_raw, &density_2_fort[0], 1, fluid_density_2_mpitype,
                         MPI_STATUS_IGNORE);
      MPI_File_write_all(dump_file_handle_raw, &velocity_2_fort[0], 1, fluid_velocity_2_mpitype,
                         MPI_STATUS_IGNORE);

      // For C output use the following but switch to MPI_ORDER_C in mpiTypeXXXWrite
      //MPI_File_write_all(dump_file_handle_raw, &density_lb[0][0][0], 1, fluid_density_2_mpitype, MPI_STATUS_IGNORE);
      //MPI_File_write_all(dump_file_handle_raw, &u_lb[0][0][0][0], 1, fluid_velocity_2_mpitype, MPI_STATUS_IGNORE);
    }
  }
}

//==========================================================================
// read in a fluid restart file.  This is only used to restart the
// fluid portion of a LAMMPS simulation.
//==========================================================================
void FixLbFluid::read_restartfile()
{
  MPI_Status status;
  MPI_Datatype realtype;
  MPI_Datatype filetype;

  int realsizes[4] = {subNbx, subNby, subNbz, numvel};
  int realstarts[4] = {1, 1, 1, 0};
  int gsizes[4] = {Nbx, Nby, Nbz, numvel};
  int lsizes[4] = {subNbx - 2, subNby - 2, subNbz - 2, numvel};
  int starts[4] = {comm->myloc[0] * (subNbx - 2), comm->myloc[1] * (subNby - 2),
                   comm->myloc[2] * (subNbz - 2), 0};
  if (domain->periodicity[2] == 0 && comm->myloc[2] == comm->procgrid[2] - 1) {
    starts[2] = comm->myloc[2] * (subNbz - 3);
  }

  MPI_Type_create_subarray(4, realsizes, lsizes, realstarts, MPI_ORDER_C, MPI_DOUBLE, &realtype);
  MPI_Type_commit(&realtype);

  MPI_Type_create_subarray(4, gsizes, lsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &filetype);
  MPI_Type_commit(&filetype);

  MPI_File_set_view(pFileRead, 0, MPI_DOUBLE, filetype, (char *) "native", MPI_INFO_NULL);
  MPI_File_seek(pFileRead, 0, MPI_SEEK_SET);
  MPI_File_read_all(pFileRead, &f_lb[0][0][0][0], 1, realtype, &status);

  MPI_Type_free(&realtype);
  MPI_Type_free(&filetype);
  MPI_File_close(&pFileRead);
}

//==========================================================================
// write a fluid restart file.
//==========================================================================
void FixLbFluid::write_restartfile()
{
  // We first create the distribution with the correct momentum on the full step instead
  // of the 1/2 step that is used in the algorithm.  The main difference is the velocity
  // shift from a 1/2 collision.  As the restart will have zero forces for first 1/2 step
  // we only take 1/2 force here.  We can use fnew as it will be overwritten in initial_integrate.
  // This ensures total momentum is conserved after a restart.

  double etacov[numvel];
  for (int i = 0; i < subNbx; i++)
    for (int j = 0; j < subNby; j++)
      for (int k = 0; k < subNbz; k++) {
        if (numvel == 15) {
          etacov[0] = 0.0;
          etacov[1] = 0.5 * (Ff[i][j][k][0] + density_lb[i][j][k] * bodyforcex);
          etacov[2] = 0.5 * (Ff[i][j][k][1] + density_lb[i][j][k] * bodyforcey);
          etacov[3] = 0.5 * (Ff[i][j][k][2] + density_lb[i][j][k] * bodyforcez);
          etacov[4] = 0.0;    // elements 4-9 also have a correction that we are ignoring
          etacov[5] = 0.0;
          etacov[6] = 0.0;
          etacov[7] = 0.0;
          etacov[8] = 0.0;
          etacov[9] = 0.0;
          etacov[10] = 0.0;
          etacov[11] = 0.0;
          etacov[12] = 0.0;
          etacov[13] = 0.0;
          etacov[14] = 0.0;
        } else {
          etacov[0] = 0.0;
          etacov[1] = 0.5 * (Ff[i][j][k][0] + density_lb[i][j][k] * bodyforcex);
          etacov[2] = 0.5 * (Ff[i][j][k][1] + density_lb[i][j][k] * bodyforcey);
          etacov[3] = 0.5 * (Ff[i][j][k][2] + density_lb[i][j][k] * bodyforcez);
          etacov[4] = 0.0;    // elements 4-9 also have a correction that we are ignoring
          etacov[5] = 0.0;
          etacov[6] = 0.0;
          etacov[7] = 0.0;
          etacov[8] = 0.0;
          etacov[9] = 0.0;
          etacov[10] = 0.0;
          etacov[11] = 0.0;
          etacov[12] = 0.0;
          etacov[13] = 0.0;
          etacov[14] = 0.0;
          etacov[15] = 0.0;
          etacov[16] = 0.0;
          etacov[17] = 0.0;
          etacov[18] = 0.0;
        }

        if (numvel == 15)
          for (int l = 0; l < 15; l++) {
            fnew[i][j][k][l] = f_lb[i][j][k][l];
            for (int ii = 0; ii < 15; ii++) {
              fnew[i][j][k][l] += w_lb15[l] * mg_lb15[ii][l] * etacov[ii] * Ng_lb15[ii];
            }
          }
        else
          for (int l = 0; l < 19; l++) {
            fnew[i][j][k][l] = f_lb[i][j][k][l];
            for (int ii = 0; ii < 19; ii++) {
              fnew[i][j][k][l] += w_lb19[l] * mg_lb19[ii][l] * etacov[ii] * Ng_lb19[ii];
            }
          }
      }

  MPI_File fh;
  MPI_Status status;
  MPI_Datatype realtype;
  MPI_Datatype filetype;

  char *hfile = utils::strdup(fmt::format("FluidRestart_{}.dat", update->ntimestep));
  MPI_File_open(world, hfile, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

  int realsizes[4] = {subNbx, subNby, subNbz, numvel};
  int realstarts[4] = {1, 1, 1, 0};
  int gsizes[4] = {Nbx, Nby, Nbz, numvel};
  int lsizes[4] = {subNbx - 2, subNby - 2, subNbz - 2, numvel};
  int starts[4] = {comm->myloc[0] * (subNbx - 2), comm->myloc[1] * (subNby - 2),
                   comm->myloc[2] * (subNbz - 2), 0};
  if (domain->periodicity[2] == 0 && comm->myloc[2] == comm->procgrid[2] - 1)
    starts[2] = comm->myloc[2] * (subNbz - 3);

  MPI_Type_create_subarray(4, realsizes, lsizes, realstarts, MPI_ORDER_C, MPI_DOUBLE, &realtype);
  MPI_Type_commit(&realtype);

  MPI_Type_create_subarray(4, gsizes, lsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &filetype);
  MPI_Type_commit(&filetype);

  MPI_File_set_view(fh, 0, MPI_DOUBLE, filetype, (char *) "native", MPI_INFO_NULL);
  MPI_File_write_all(fh, &fnew[0][0][0][0], 1, realtype, &status);

  MPI_Type_free(&realtype);
  MPI_Type_free(&filetype);
  MPI_File_close(&fh);
  delete[] hfile;
}

//==========================================================================
// Compute the lattice Boltzmann equilibrium distribution functions for
// the D3Q15 model.
//==========================================================================
void FixLbFluid::equilibriumdist15(int xstart, int xend, int ystart, int yend, int zstart, int zend)
{
  for (int i = xstart; i < xend; ++i) {
    for (int j = ystart; j < yend; ++j) {
      for (int k = zstart; k < zend; ++k) {

        if (sublattice[i][j][k].type != 2) {
          const double rho = density_lb[i][j][k];
          const double p0 = rho * a_0;
          const double tauR = tau - 0.5;

          const double Fx_w = Ff[i][j][k][0];
          const double Fy_w = Ff[i][j][k][1];
          const double Fz_w = Ff[i][j][k][2];

          double etacov[15];
          etacov[0] = rho;
          etacov[1] = rho * u_lb[i][j][k][0] + (Fx_w + rho * bodyforcex) * tauR;
          etacov[2] = rho * u_lb[i][j][k][1] + (Fy_w + rho * bodyforcey) * tauR;
          etacov[3] = rho * u_lb[i][j][k][2] + (Fz_w + rho * bodyforcez) * tauR;

          etacov[4] = p0 + rho * u_lb[i][j][k][0] * u_lb[i][j][k][0] - rho / 3.0 +
              tauR * (2.0 * u_lb[i][j][k][0] * (Fx_w + rho * bodyforcex));
          etacov[5] = p0 + rho * u_lb[i][j][k][1] * u_lb[i][j][k][1] - rho / 3.0 +
              tauR * (2.0 * u_lb[i][j][k][1] * (Fy_w + rho * bodyforcey));
          etacov[6] = p0 + rho * u_lb[i][j][k][2] * u_lb[i][j][k][2] - rho / 3.0 +
              tauR * (2.0 * u_lb[i][j][k][2] * (Fz_w + rho * bodyforcez));
          etacov[7] = rho * u_lb[i][j][k][0] * u_lb[i][j][k][1] +
              tauR *
                  (u_lb[i][j][k][0] * (Fy_w + rho * bodyforcey) +
                   (Fx_w + rho * bodyforcex) * u_lb[i][j][k][1]);
          etacov[8] = rho * u_lb[i][j][k][1] * u_lb[i][j][k][2] +
              tauR *
                  (u_lb[i][j][k][1] * (Fz_w + rho * bodyforcez) +
                   (Fy_w + rho * bodyforcey) * u_lb[i][j][k][2]);
          etacov[9] = rho * u_lb[i][j][k][0] * u_lb[i][j][k][2] +
              tauR *
                  (u_lb[i][j][k][0] * (Fz_w + rho * bodyforcez) +
                   (Fx_w + rho * bodyforcex) * u_lb[i][j][k][2]);
          etacov[10] = 0.0;
          etacov[11] = 0.0;
          etacov[12] = 0.0;
          etacov[13] = rho * u_lb[i][j][k][0] * u_lb[i][j][k][1] * u_lb[i][j][k][2];
          etacov[14] = K_0 * rho * (1.0 - 3.0 * a_0);    // should be looked at if kappa != 0;

          for (int l = 0; l < 15; l++) {
            feq[i][j][k][l] = 0.0;
            for (int ii = 0; ii < 15; ii++)
              feq[i][j][k][l] += w_lb15[l] * mg_lb15[ii][l] * etacov[ii] * Ng_lb15[ii];
          }

          if (noisestress == 1) {
            const double stdv = sqrt(namp * rho);
            double S[2][3];

            for (int jj = 0; jj < 3; jj++) S[0][jj] = stdv * random->gaussian();
            for (int jj = 0; jj < 3; jj++) S[1][jj] = stdv * random->gaussian();

            etacov[4] = (S[0][0] * sqrt(3.0 - 3.0 * a_0));
            etacov[5] = ((1.0 - 3.0 * a_0) * S[0][0] / sqrt(3.0 - 3.0 * a_0) +
                         sqrt((8.0 - 12.0 * a_0) / (3.0 - 3.0 * a_0)) * S[0][1]);
            etacov[6] =
                ((1.0 - 3.0 * a_0) * S[0][0] / sqrt(3.0 - 3.0 * a_0) +
                 (2.0 - 6.0 * a_0) * S[0][1] / sqrt((8.0 - 12.0 * a_0) * (3.0 - 3.0 * a_0)) +
                 sqrt((5.0 - 9.0 * a_0) / (2.0 - 3.0 * a_0)) * S[0][2]);
            etacov[7] = S[1][0];
            etacov[8] = S[1][1];
            etacov[9] = S[1][2];

            for (int l = 10; l < 15; l++) {
              etacov[l] = sqrt(9.0 * namp * rho / Ng_lb15[l]) * random->gaussian();
            }
            etacov[14] +=
                -K_0 * (etacov[4] + etacov[5] + etacov[6]);    //correction from noise to TrP

            for (int l = 0; l < 15; l++) {
              const double ghostnoise = w_lb15[l] *
                  (mg_lb15[4][l] * etacov[4] * Ng_lb15[4] + mg_lb15[5][l] * etacov[5] * Ng_lb15[5] +
                   mg_lb15[6][l] * etacov[6] * Ng_lb15[6] + mg_lb15[7][l] * etacov[7] * Ng_lb15[7] +
                   mg_lb15[8][l] * etacov[8] * Ng_lb15[8] + mg_lb15[9][l] * etacov[9] * Ng_lb15[9] +
                   mg_lb15[10][l] * etacov[10] * Ng_lb15[10] +
                   mg_lb15[11][l] * etacov[11] * Ng_lb15[11] +
                   mg_lb15[12][l] * etacov[12] * Ng_lb15[12] +
                   mg_lb15[13][l] * etacov[13] * Ng_lb15[13] +
                   mg_lb15[14][l] * etacov[14] * Ng_lb15[14]);
              feq[i][j][k][l] += ghostnoise * noisefactor;
            }
          }
        } else {    // non-active site, this should be redundant
          memset(feq[i][j][k], 0, 15 * sizeof(double));
        }
      }
    }
  }
}

//==========================================================================
// Compute the lattice Boltzmann equilibrium distribution functions for
// the D3Q19 model.
//==========================================================================
void FixLbFluid::equilibriumdist19(int xstart, int xend, int ystart, int yend, int zstart, int zend)
{
  for (int i = xstart; i < xend; ++i) {
    for (int j = ystart; j < yend; ++j) {
      for (int k = zstart; k < zend; ++k) {
        const int iup = i + 1;
        const int idwn = i - 1;
        const int jup = j + 1;
        const int jdwn = j - 1;
        const int kup = k + 1;
        const int kdwn = k - 1;

        const double rho = density_lb[i][j][k];
        double drhox, drhoy, drhoz, drhoxx, drhoyy, drhozz;
        if (kappa_lb > 0) {    // kappa_lb is the square gradient coeff in the pressure tensor
          // Derivatives.
          drhox = (density_lb[iup][j][k] - density_lb[idwn][j][k]) / 2.0;
          drhoxx = (density_lb[iup][j][k] - 2.0 * density_lb[i][j][k] + density_lb[idwn][j][k]);

          drhoy = (density_lb[i][jup][k] - density_lb[i][jdwn][k]) / 2.0;
          drhoyy = (density_lb[i][jup][k] - 2.0 * density_lb[i][j][k] + density_lb[i][jdwn][k]);

          drhoz = (density_lb[i][j][kup] - density_lb[i][j][kdwn]) / 2.0;
          drhozz = (density_lb[i][j][kup] - 2.0 * density_lb[i][j][k] + density_lb[i][j][kdwn]);

          // Need one-sided derivatives for the boundary of the domain, if fixed boundary
          // conditions are used.

          if (domain->periodicity[2] == 0) {
            if (comm->myloc[2] == 0 && k == 1) {
              drhoz = (-3.0 * density_lb[i][j][k] + 4.0 * density_lb[i][j][k + 1] -
                       density_lb[i][j][k + 2]) /
                  2.0;
              drhozz = (-density_lb[i][j][k + 3] + 4.0 * density_lb[i][j][k + 2] -
                        5.0 * density_lb[i][j][k + 1] + 2.0 * rho);
            }
            if (comm->myloc[2] == comm->procgrid[2] - 1 && k == subNbz - 2) {
              drhoz = -(-3.0 * density_lb[i][j][k] + 4.0 * density_lb[i][j][k - 1] -
                        density_lb[i][j][k - 2]) /
                  2.0;
              drhozz = (-density_lb[i][j][k - 3] + 4.0 * density_lb[i][j][k - 2] -
                        5.0 * density_lb[i][j][k - 1] + 2.0 * rho);
            }
          }
          // clang-format on
        } else {
          drhox = drhoy = drhoz = 0.0;
          drhoxx = drhoyy = drhozz = 0.0;
        }
        const double grs = drhox * drhox + drhoy * drhoy + drhoz * drhoz;
        const double p0 = rho * a_0 - kappa_lb * rho * (drhoxx + drhoyy + drhozz);
        const double dPdrho = a_0;    //assuming here that kappa_lb = 0.
        const double tauR = tau - 0.5;

        const double Pxx = p0 + kappa_lb * (drhox * drhox - 0.5 * grs) +
            tauR * (1.0 / 3.0 - dPdrho) *
                (3.0 * u_lb[i][j][k][0] * drhox + u_lb[i][j][k][1] * drhoy +
                 u_lb[i][j][k][2] * drhoz);
        const double Pyy = p0 + kappa_lb * (drhoy * drhoy - 0.5 * grs) +
            tauR * (1.0 / 3.0 - dPdrho) *
                (u_lb[i][j][k][0] * drhox + 3.0 * u_lb[i][j][k][1] * drhoy +
                 u_lb[i][j][k][2] * drhoz);
        const double Pzz = p0 + kappa_lb * (drhoz * drhoz - 0.5 * grs) +
            tauR * (1.0 / 3.0 - dPdrho) *
                (u_lb[i][j][k][0] * drhox + u_lb[i][j][k][1] * drhoy +
                 3.0 * u_lb[i][j][k][2] * drhoz);
        const double Pxy = kappa_lb * drhox * drhoy +
            tauR * (1.0 / 3.0 - dPdrho) * (u_lb[i][j][k][0] * drhoy + u_lb[i][j][k][1] * drhox);
        const double Pxz = kappa_lb * drhox * drhoz +
            tauR * (1.0 / 3.0 - dPdrho) * (u_lb[i][j][k][0] * drhoz + u_lb[i][j][k][2] * drhox);
        const double Pyz = kappa_lb * drhoy * drhoz +
            tauR * (1.0 / 3.0 - dPdrho) * (u_lb[i][j][k][1] * drhoz + u_lb[i][j][k][2] * drhoy);

        const double Fx_w = Ff[i][j][k][0];
        const double Fy_w = Ff[i][j][k][1];
        const double Fz_w = Ff[i][j][k][2];

        double etacov[19];
        etacov[0] = rho;
        etacov[1] = rho * u_lb[i][j][k][0] + (Fx_w + rho * bodyforcex) * tauR;
        etacov[2] = rho * u_lb[i][j][k][1] + (Fy_w + rho * bodyforcey) * tauR;
        etacov[3] = rho * u_lb[i][j][k][2] + (Fz_w + rho * bodyforcez) * tauR;

        etacov[4] = Pxx + rho * u_lb[i][j][k][0] * u_lb[i][j][k][0] - rho / 3.0 +
            tauR * (2.0 * u_lb[i][j][k][0] * (Fx_w + rho * bodyforcex));
        etacov[5] = Pyy + rho * u_lb[i][j][k][1] * u_lb[i][j][k][1] - rho / 3.0 +
            tauR * (2.0 * u_lb[i][j][k][1] * (Fy_w + rho * bodyforcey));
        etacov[6] = Pzz + rho * u_lb[i][j][k][2] * u_lb[i][j][k][2] - rho / 3.0 +
            tauR * (2.0 * u_lb[i][j][k][2] * (Fz_w + rho * bodyforcez));
        etacov[7] = Pxy + rho * u_lb[i][j][k][0] * u_lb[i][j][k][1] +
            tauR *
                (u_lb[i][j][k][0] * (Fy_w + rho * bodyforcey) +
                 (Fx_w + rho * bodyforcex) * u_lb[i][j][k][1]);
        etacov[8] = Pxz + rho * u_lb[i][j][k][0] * u_lb[i][j][k][2] +
            tauR *
                (u_lb[i][j][k][0] * (Fz_w + rho * bodyforcez) +
                 (Fx_w + rho * bodyforcex) * u_lb[i][j][k][2]);
        etacov[9] = Pyz + rho * u_lb[i][j][k][1] * u_lb[i][j][k][2] +
            tauR *
                (u_lb[i][j][k][1] * (Fz_w + rho * bodyforcez) +
                 (Fy_w + rho * bodyforcey) * u_lb[i][j][k][2]);
        etacov[10] = 0.0;
        etacov[11] = 0.0;
        etacov[12] = 0.0;
        etacov[13] = 0.0;
        etacov[14] = 0.0;
        etacov[15] = 0.0;
        etacov[16] = 0.0;
        etacov[17] = 0.0;
        etacov[18] = 0.0;

        for (int l = 0; l < 19; l++) {
          feq[i][j][k][l] = 0.0;
          for (int ii = 0; ii < 19; ii++)
            feq[i][j][k][l] += w_lb19[l] * mg_lb19[ii][l] * etacov[ii] * Ng_lb19[ii];
        }

        if (noisestress == 1) {
          const double std = sqrt(namp * rho);
          double S[2][3];

          for (int jj = 0; jj < 3; jj++) S[0][jj] = std * random->gaussian();
          for (int jj = 0; jj < 3; jj++) S[1][jj] = std * random->gaussian();

          etacov[4] = (S[0][0] * sqrt(3.0 - 3.0 * a_0));
          etacov[5] = ((1.0 - 3.0 * a_0) * S[0][0] / sqrt(3.0 - 3.0 * a_0) +
                       sqrt((8.0 - 12.0 * a_0) / (3.0 - 3.0 * a_0)) * S[0][1]);
          etacov[6] = ((1.0 - 3.0 * a_0) * S[0][0] / sqrt(3.0 - 3.0 * a_0) +
                       (2.0 - 6.0 * a_0) * S[0][1] / sqrt((8.0 - 12.0 * a_0) * (3.0 - 3.0 * a_0)) +
                       sqrt((5.0 - 9.0 * a_0) / (2.0 - 3.0 * a_0)) * S[0][2]);
          etacov[7] = S[1][0];
          etacov[8] = S[1][1];
          etacov[9] = S[1][2];

          for (int l = 10; l < 19; l++) {
            etacov[l] = sqrt(9.0 * namp * rho / Ng_lb19[l]) * random->gaussian();
          }

          for (int l = 0; l < 19; l++) {
            const double ghostnoise = w_lb19[l] *
                (mg_lb19[4][l] * etacov[4] * Ng_lb19[4] + mg_lb19[5][l] * etacov[5] * Ng_lb19[5] +
                 mg_lb19[6][l] * etacov[6] * Ng_lb19[6] + mg_lb19[7][l] * etacov[7] * Ng_lb19[7] +
                 mg_lb19[8][l] * etacov[8] * Ng_lb19[8] + mg_lb19[9][l] * etacov[9] * Ng_lb19[9] +
                 mg_lb19[10][l] * etacov[10] * Ng_lb19[10] +
                 mg_lb19[11][l] * etacov[11] * Ng_lb19[11] +
                 mg_lb19[12][l] * etacov[12] * Ng_lb19[12] +
                 mg_lb19[13][l] * etacov[13] * Ng_lb19[13] +
                 mg_lb19[14][l] * etacov[14] * Ng_lb19[14] +
                 mg_lb19[15][l] * etacov[15] * Ng_lb19[15] +
                 mg_lb19[16][l] * etacov[16] * Ng_lb19[16] +
                 mg_lb19[17][l] * etacov[17] * Ng_lb19[17] +
                 mg_lb19[18][l] * etacov[18] * Ng_lb19[18]);
            feq[i][j][k][l] += ghostnoise * noisefactor;
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
void FixLbFluid::parametercalc_full()
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
  for (i = 0; i < 4; i++) requests[i] = MPI_REQUEST_NULL;
  MPI_Isend(&f_lb[1][1][1][0], 1, passxf, comm->procneigh[0][0], 10, world, &requests[0]);
  MPI_Irecv(&f_lb[0][1][1][0], 1, passxf, comm->procneigh[0][0], 20, world, &requests[1]);
  MPI_Isend(&f_lb[subNbx - 2][1][1][0], 1, passxf, comm->procneigh[0][1], 20, world, &requests[2]);
  MPI_Irecv(&f_lb[subNbx - 1][1][1][0], 1, passxf, comm->procneigh[0][1], 10, world, &requests[3]);
  parametercalc_part(1, subNbx - 1, 1, subNby - 1, 1, subNbz - 1);
  MPI_Waitall(4, requests, MPI_STATUS_IGNORE);

  for (i = 0; i < 4; i++) requests[i] = MPI_REQUEST_NULL;
  MPI_Isend(&f_lb[0][1][1][0], 1, passyf, comm->procneigh[1][0], 10, world, &requests[0]);
  MPI_Irecv(&f_lb[0][0][1][0], 1, passyf, comm->procneigh[1][0], 20, world, &requests[1]);
  MPI_Isend(&f_lb[0][subNby - 2][1][0], 1, passyf, comm->procneigh[1][1], 20, world, &requests[2]);
  MPI_Irecv(&f_lb[0][subNby - 1][1][0], 1, passyf, comm->procneigh[1][1], 10, world, &requests[3]);
  parametercalc_part(0, 1, 1, subNby - 1, 1, subNbz - 1);
  parametercalc_part(subNbx - 1, subNbx, 1, subNby - 1, 1, subNbz - 1);
  MPI_Waitall(4, requests, MPI_STATUS_IGNORE);

  for (i = 0; i < 4; i++) requests[i] = MPI_REQUEST_NULL;
  MPI_Isend(&f_lb[0][0][1][0], 1, passzf, comm->procneigh[2][0], 10, world, &requests[0]);
  MPI_Irecv(&f_lb[0][0][0][0], 1, passzf, comm->procneigh[2][0], 20, world, &requests[1]);
  MPI_Isend(&f_lb[0][0][subNbz - 2][0], 1, passzf, comm->procneigh[2][1], 20, world, &requests[2]);
  MPI_Irecv(&f_lb[0][0][subNbz - 1][0], 1, passzf, comm->procneigh[2][1], 10, world, &requests[3]);
  parametercalc_part(0, subNbx, 0, 1, 1, subNbz - 1);
  parametercalc_part(0, subNbx, subNby - 1, subNby, 1, subNbz - 1);
  MPI_Waitall(4, requests, MPI_STATUS_IGNORE);

  parametercalc_part(0, subNbx, 0, subNby, 0, 1);
  parametercalc_part(0, subNbx, 0, subNby, subNbz - 1, subNbz);

  //--------------------------------------------------------------------------
  // Send the remaining portions of the u array and density array
  //--------------------------------------------------------------------------
  numrequests = 12;

  for (i = 0; i < numrequests; i++) requests2[i] = MPI_REQUEST_NULL;
  MPI_Isend(&u_lb[2][0][0][0], 1, passxu, comm->procneigh[0][0], 10, world, &requests2[0]);
  MPI_Isend(&u_lb[3][0][0][0], 1, passxu, comm->procneigh[0][0], 20, world, &requests2[1]);
  MPI_Isend(&u_lb[subNbx - 3][0][0][0], 1, passxu, comm->procneigh[0][1], 30, world, &requests2[2]);
  MPI_Irecv(&u_lb[subNbx][0][0][0], 1, passxu, comm->procneigh[0][1], 10, world, &requests2[3]);
  MPI_Irecv(&u_lb[subNbx + 1][0][0][0], 1, passxu, comm->procneigh[0][1], 20, world, &requests2[4]);
  MPI_Irecv(&u_lb[subNbx + 2][0][0][0], 1, passxu, comm->procneigh[0][0], 30, world, &requests2[5]);

  MPI_Isend(&density_lb[2][0][0], 1, passxrho, comm->procneigh[0][0], 40, world, &requests2[6]);
  MPI_Isend(&density_lb[3][0][0], 1, passxrho, comm->procneigh[0][0], 50, world, &requests2[7]);
  MPI_Isend(&density_lb[subNbx - 3][0][0], 1, passxrho, comm->procneigh[0][1], 60, world,
            &requests2[8]);
  MPI_Irecv(&density_lb[subNbx][0][0], 1, passxrho, comm->procneigh[0][1], 40, world,
            &requests2[9]);
  MPI_Irecv(&density_lb[subNbx + 1][0][0], 1, passxrho, comm->procneigh[0][1], 50, world,
            &requests2[10]);
  MPI_Irecv(&density_lb[subNbx + 2][0][0], 1, passxrho, comm->procneigh[0][0], 60, world,
            &requests2[11]);

  MPI_Waitall(numrequests, requests2, MPI_STATUS_IGNORE);

  for (i = 0; i < numrequests; i++) requests2[i] = MPI_REQUEST_NULL;
  MPI_Isend(&u_lb[0][2][0][0], 1, passyu, comm->procneigh[1][0], 10, world, &requests2[0]);
  MPI_Isend(&u_lb[0][3][0][0], 1, passyu, comm->procneigh[1][0], 20, world, &requests2[1]);
  MPI_Isend(&u_lb[0][subNby - 3][0][0], 1, passyu, comm->procneigh[1][1], 30, world, &requests2[2]);
  MPI_Irecv(&u_lb[0][subNby][0][0], 1, passyu, comm->procneigh[1][1], 10, world, &requests2[3]);
  MPI_Irecv(&u_lb[0][subNby + 1][0][0], 1, passyu, comm->procneigh[1][1], 20, world, &requests2[4]);
  MPI_Irecv(&u_lb[0][subNby + 2][0][0], 1, passyu, comm->procneigh[1][0], 30, world, &requests2[5]);

  MPI_Isend(&density_lb[0][2][0], 1, passyrho, comm->procneigh[1][0], 40, world, &requests2[6]);
  MPI_Isend(&density_lb[0][3][0], 1, passyrho, comm->procneigh[1][0], 50, world, &requests2[7]);
  MPI_Isend(&density_lb[0][subNby - 3][0], 1, passyrho, comm->procneigh[1][1], 60, world,
            &requests2[8]);
  MPI_Irecv(&density_lb[0][subNby][0], 1, passyrho, comm->procneigh[1][1], 40, world,
            &requests2[9]);
  MPI_Irecv(&density_lb[0][subNby + 1][0], 1, passyrho, comm->procneigh[1][1], 50, world,
            &requests2[10]);
  MPI_Irecv(&density_lb[0][subNby + 2][0], 1, passyrho, comm->procneigh[1][0], 60, world,
            &requests2[11]);

  MPI_Waitall(numrequests, requests2, MPI_STATUS_IGNORE);

  for (i = 0; i < 12; i++) requests2[i] = MPI_REQUEST_NULL;
  int requestcount = 0;
  if (domain->periodicity[2] != 0 || comm->myloc[2] != 0) {
    MPI_Isend(&u_lb[0][0][2][0], 1, passzu, comm->procneigh[2][0], 10, world,
              &requests2[requestcount]);
    MPI_Isend(&u_lb[0][0][3][0], 1, passzu, comm->procneigh[2][0], 20, world,
              &requests2[requestcount + 1]);
    MPI_Irecv(&u_lb[0][0][subNbz + 2][0], 1, passzu, comm->procneigh[2][0], 30, world,
              &requests2[requestcount + 2]);
    requestcount = requestcount + 3;

    MPI_Isend(&density_lb[0][0][2], 1, passzrho, comm->procneigh[2][0], 40, world,
              &requests2[requestcount]);
    MPI_Isend(&density_lb[0][0][3], 1, passzrho, comm->procneigh[2][0], 50, world,
              &requests2[requestcount + 1]);
    MPI_Irecv(&density_lb[0][0][subNbz + 2], 1, passzrho, comm->procneigh[2][0], 60, world,
              &requests2[requestcount + 2]);
    requestcount = requestcount + 3;
  }
  if (domain->periodicity[2] != 0 || comm->myloc[2] != (comm->procgrid[2] - 1)) {
    MPI_Isend(&u_lb[0][0][subNbz - 3][0], 1, passzu, comm->procneigh[2][1], 30, world,
              &requests2[requestcount]);
    MPI_Irecv(&u_lb[0][0][subNbz][0], 1, passzu, comm->procneigh[2][1], 10, world,
              &requests2[requestcount + 1]);
    MPI_Irecv(&u_lb[0][0][subNbz + 1][0], 1, passzu, comm->procneigh[2][1], 20, world,
              &requests2[requestcount + 2]);
    requestcount = requestcount + 3;

    MPI_Isend(&density_lb[0][0][subNbz - 3], 1, passzrho, comm->procneigh[2][1], 60, world,
              &requests2[requestcount]);
    MPI_Irecv(&density_lb[0][0][subNbz], 1, passzrho, comm->procneigh[2][1], 40, world,
              &requests2[requestcount + 1]);
    MPI_Irecv(&density_lb[0][0][subNbz + 1], 1, passzrho, comm->procneigh[2][1], 50, world,
              &requests2[requestcount + 2]);
    requestcount = requestcount + 3;
  }
  MPI_Waitall(requestcount, requests2, MPI_STATUS_IGNORE);
}

//==========================================================================
// Calculate the fluid density and velocity over a simulation volume
// specified by xstart,xend; ystart,yend; zstart,zend.
//==========================================================================
void FixLbFluid::parametercalc_part(int xstart, int xend, int ystart, int yend, int zstart,
                                    int zend)
{
  for (int i = xstart; i < xend; i++) {
    for (int j = ystart; j < yend; j++) {
      for (int k = zstart; k < zend; k++) {

        density_lb[i][j][k] = 0.0;
        u_lb[i][j][k][0] = 0.0;
        u_lb[i][j][k][1] = 0.0;
        u_lb[i][j][k][2] = 0.0;

        if (sublattice[i][j][k].type != 2) {    // type 2 is outside domain
          if (numvel == 15) {
            for (int m = 0; m < 15; m++) {
              density_lb[i][j][k] += f_lb[i][j][k][m];

              u_lb[i][j][k][0] += f_lb[i][j][k][m] * e15[m][0];
              u_lb[i][j][k][1] += f_lb[i][j][k][m] * e15[m][1];
              u_lb[i][j][k][2] += f_lb[i][j][k][m] * e15[m][2];
            }
          } else {
            for (int m = 0; m < 19; m++) {
              density_lb[i][j][k] += f_lb[i][j][k][m];

              u_lb[i][j][k][0] += f_lb[i][j][k][m] * e19[m][0];
              u_lb[i][j][k][1] += f_lb[i][j][k][m] * e19[m][1];
              u_lb[i][j][k][2] += f_lb[i][j][k][m] * e19[m][2];
            }
          }
          u_lb[i][j][k][0] = u_lb[i][j][k][0] / density_lb[i][j][k];
          u_lb[i][j][k][1] = u_lb[i][j][k][1] / density_lb[i][j][k];
          u_lb[i][j][k][2] = u_lb[i][j][k][2] / density_lb[i][j][k];
        }
      }
    }
  }
}

//==========================================================================
// Add force contribution to fluid velocity over the entire simulation
// domain.
//==========================================================================
void FixLbFluid::correctu_full()
{
  MPI_Request requests2[12];
  int numrequests;
  int i;

  correctu_part(1, subNbx - 1, 1, subNby - 1, 1, subNbz - 1);

  //--------------------------------------------------------------------------
  // Send the remaining portions of the u array
  //--------------------------------------------------------------------------
  numrequests = 10;

  for (i = 0; i < numrequests; i++) requests2[i] = MPI_REQUEST_NULL;
  MPI_Isend(&u_lb[1][0][0][0], 1, passxu, comm->procneigh[0][0], 5, world, &requests2[0]);
  MPI_Isend(&u_lb[2][0][0][0], 1, passxu, comm->procneigh[0][0], 10, world, &requests2[1]);
  MPI_Isend(&u_lb[3][0][0][0], 1, passxu, comm->procneigh[0][0], 20, world, &requests2[2]);
  MPI_Isend(&u_lb[subNbx - 3][0][0][0], 1, passxu, comm->procneigh[0][1], 30, world, &requests2[3]);
  MPI_Isend(&u_lb[subNbx - 2][0][0][0], 1, passxu, comm->procneigh[0][1], 35, world, &requests2[4]);
  MPI_Irecv(&u_lb[subNbx - 1][0][0][0], 1, passxu, comm->procneigh[0][1], 5, world, &requests2[5]);
  MPI_Irecv(&u_lb[subNbx][0][0][0], 1, passxu, comm->procneigh[0][1], 10, world, &requests2[6]);
  MPI_Irecv(&u_lb[subNbx + 1][0][0][0], 1, passxu, comm->procneigh[0][1], 20, world, &requests2[7]);
  MPI_Irecv(&u_lb[subNbx + 2][0][0][0], 1, passxu, comm->procneigh[0][0], 30, world, &requests2[8]);
  MPI_Irecv(&u_lb[0][0][0][0], 1, passxu, comm->procneigh[0][0], 35, world, &requests2[9]);

  MPI_Waitall(numrequests, requests2, MPI_STATUS_IGNORE);

  for (i = 0; i < numrequests; i++) requests2[i] = MPI_REQUEST_NULL;
  MPI_Isend(&u_lb[0][1][0][0], 1, passyu, comm->procneigh[1][0], 5, world, &requests2[0]);
  MPI_Isend(&u_lb[0][2][0][0], 1, passyu, comm->procneigh[1][0], 10, world, &requests2[1]);
  MPI_Isend(&u_lb[0][3][0][0], 1, passyu, comm->procneigh[1][0], 20, world, &requests2[2]);
  MPI_Isend(&u_lb[0][subNby - 3][0][0], 1, passyu, comm->procneigh[1][1], 30, world, &requests2[3]);
  MPI_Isend(&u_lb[0][subNby - 2][0][0], 1, passyu, comm->procneigh[1][1], 35, world, &requests2[4]);
  MPI_Irecv(&u_lb[0][subNby - 1][0][0], 1, passyu, comm->procneigh[1][1], 5, world, &requests2[5]);
  MPI_Irecv(&u_lb[0][subNby][0][0], 1, passyu, comm->procneigh[1][1], 10, world, &requests2[6]);
  MPI_Irecv(&u_lb[0][subNby + 1][0][0], 1, passyu, comm->procneigh[1][1], 20, world, &requests2[7]);
  MPI_Irecv(&u_lb[0][subNby + 2][0][0], 1, passyu, comm->procneigh[1][0], 30, world, &requests2[8]);
  MPI_Irecv(&u_lb[0][0][0][0], 1, passyu, comm->procneigh[1][0], 35, world, &requests2[9]);

  MPI_Waitall(numrequests, requests2, MPI_STATUS_IGNORE);

  for (i = 0; i < 12; i++) requests2[i] = MPI_REQUEST_NULL;
  int requestcount = 0;
  if (domain->periodicity[2] != 0 || comm->myloc[2] != 0) {
    MPI_Isend(&u_lb[0][0][1][0], 1, passzu, comm->procneigh[2][0], 5, world,
              &requests2[requestcount]);
    MPI_Isend(&u_lb[0][0][2][0], 1, passzu, comm->procneigh[2][0], 10, world,
              &requests2[requestcount + 1]);
    MPI_Isend(&u_lb[0][0][3][0], 1, passzu, comm->procneigh[2][0], 20, world,
              &requests2[requestcount + 2]);
    MPI_Irecv(&u_lb[0][0][subNbz + 2][0], 1, passzu, comm->procneigh[2][0], 30, world,
              &requests2[requestcount + 3]);
    MPI_Irecv(&u_lb[0][0][0][0], 1, passzu, comm->procneigh[2][0], 35, world,
              &requests2[requestcount + 4]);
    requestcount = requestcount + 5;
  }
  if (domain->periodicity[2] != 0 || comm->myloc[2] != (comm->procgrid[2] - 1)) {
    MPI_Isend(&u_lb[0][0][subNbz - 3][0], 1, passzu, comm->procneigh[2][1], 30, world,
              &requests2[requestcount]);
    MPI_Isend(&u_lb[0][0][subNbz - 2][0], 1, passzu, comm->procneigh[2][1], 35, world,
              &requests2[requestcount + 1]);
    MPI_Irecv(&u_lb[0][0][subNbz - 1][0], 1, passzu, comm->procneigh[2][1], 5, world,
              &requests2[requestcount + 2]);
    MPI_Irecv(&u_lb[0][0][subNbz][0], 1, passzu, comm->procneigh[2][1], 10, world,
              &requests2[requestcount + 3]);
    MPI_Irecv(&u_lb[0][0][subNbz + 1][0], 1, passzu, comm->procneigh[2][1], 20, world,
              &requests2[requestcount + 4]);
    requestcount = requestcount + 5;
  }
  MPI_Waitall(requestcount, requests2, MPI_STATUS_IGNORE);
}

//==========================================================================
// Add force contribution to fluid velocity
//==========================================================================
void FixLbFluid::correctu_part(int xstart, int xend, int ystart, int yend, int zstart, int zend)
{
  for (int i = xstart; i < xend; i++) {
    for (int j = ystart; j < yend; j++) {
      for (int k = zstart; k < zend; k++) {

        if (sublattice[i][j][k].type ==
            0) {    // update bulk sites, not boundaries (with set u), or exterior
          u_lb[i][j][k][0] += 0.5 * (Ff[i][j][k][0] / density_lb[i][j][k] + bodyforcex);
          u_lb[i][j][k][1] += 0.5 * (Ff[i][j][k][1] / density_lb[i][j][k] + bodyforcey);
          u_lb[i][j][k][2] += 0.5 * (Ff[i][j][k][2] / density_lb[i][j][k] + bodyforcez);
        }
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
  for (int i = xstart; i < xend; i++)
    for (int j = ystart; j < yend; j++)
      for (int k = zstart; k < zend; k++) {

        int type = sublattice[i][j][k].type;
        int ori = sublattice[i][j][k].orientation;

        if (type == 0) {    // bulk fluid nodes
          if (numvel == 15) {
            for (int m = 0; m < 15; m++) {
              int imod = i - e15[m][0];
              int jmod = j - e15[m][1];
              int kmod = k - e15[m][2];

              fnew[i][j][k][m] = f_lb[imod][jmod][kmod][m] +
                  (feq[imod][jmod][kmod][m] - f_lb[imod][jmod][kmod][m]) / tau;
            }
          } else {
            for (int m = 0; m < 19; m++) {
              int imod = i - e19[m][0];
              int jmod = j - e19[m][1];
              int kmod = k - e19[m][2];

              fnew[i][j][k][m] = f_lb[imod][jmod][kmod][m] +
                  (feq[imod][jmod][kmod][m] - f_lb[imod][jmod][kmod][m]) / tau;
            }
          }
        } else if (type == 1) {    // pit geometry boundary fluid nodes
          // bounce back
          for (int nbb = 1; nbb <= bbl[ori][0]; nbb++) {
            int ll = bbl[ori][nbb];

            fnew[i][j][k][ll] =
                f_lb[i][j][k][od[ll]] - (f_lb[i][j][k][od[ll]] - feq[i][j][k][od[ll]]) / tau;
          }
          // normal propagation
          fnew[i][j][k][0] = f_lb[i][j][k][0] - (f_lb[i][j][k][0] - feq[i][j][k][0]) / tau;
          for (int nbb = bbl[ori][0] + 1; nbb <= 15; nbb++) {
            int ll = bbl[ori][nbb];
            int imod = (i - e15[ll][0]);
            int jmod = (j - e15[ll][1]);
            int kmod = (k - e15[ll][2]);

            fnew[i][j][k][ll] = f_lb[imod][jmod][kmod][ll] -
                (f_lb[imod][jmod][kmod][ll] - feq[imod][jmod][kmod][ll]) / tau;
          }
        } else if (type == 3) {    // on-wall z boundaries
          if (numvel == 15) {
            // normal propagation directions (includes on-wall which will be modified)
            fnew[i][j][k][0] = f_lb[i][j][k][0] - (f_lb[i][j][k][0] - feq[i][j][k][0]) / tau;
            for (int nbb = bbl[ori][0] + 1; nbb <= 15; nbb++) {
              int ll = bbl[ori][nbb];
              int imod = (i - e15[ll][0]);
              int jmod = (j - e15[ll][1]);
              int kmod = (k - e15[ll][2]);

              fnew[i][j][k][ll] = f_lb[imod][jmod][kmod][ll] -
                  (f_lb[imod][jmod][kmod][ll] - feq[imod][jmod][kmod][ll]) / tau;
            }
            if (ori == 3) {    // bottom wall, modified return
              double rb, rbvw, rvw2, seqyy, r;
              rb = fnew[i][j][k][0] + fnew[i][j][k][1] + fnew[i][j][k][2] + fnew[i][j][k][3] +
                  fnew[i][j][k][4] + fnew[i][j][k][6] + fnew[i][j][k][11] + fnew[i][j][k][12] +
                  fnew[i][j][k][13] + fnew[i][j][k][14] + f_lb[i][j][k][6] +
                  (feq[i][j][k][6] - f_lb[i][j][k][6]) / tau + f_lb[i][j][k][11] +
                  (feq[i][j][k][11] - f_lb[i][j][k][11]) / tau + f_lb[i][j][k][12] +
                  (feq[i][j][k][12] - f_lb[i][j][k][12]) / tau + f_lb[i][j][k][13] +
                  (feq[i][j][k][13] - f_lb[i][j][k][13]) / tau + f_lb[i][j][k][14] +
                  (feq[i][j][k][14] - f_lb[i][j][k][14]) / tau;
              rbvw = rb * vwbt;
              rvw2 = rbvw * vwbt;
              seqyy = rb * a_0 + rvw2;
              r = rb / density_lb[i][j][k];

              fnew[i][j][k][2] = r * feq[i][j][k][2];
              fnew[i][j][k][4] = r * feq[i][j][k][4];
              fnew[i][j][k][0] = rb -
                  2.0 *
                      (fnew[i][j][k][6] + fnew[i][j][k][11] + fnew[i][j][k][12] +
                       fnew[i][j][k][13] + fnew[i][j][k][14] + fnew[i][j][k][2] +
                       fnew[i][j][k][4]) +
                  rvw2;
              fnew[i][j][k][1] = 0.5 * (fnew[i][j][k][2] + fnew[i][j][k][4] - rvw2);
              fnew[i][j][k][3] = fnew[i][j][k][1];
              fnew[i][j][k][5] = fnew[i][j][k][6] +
                  2.0 *
                      (fnew[i][j][k][11] + fnew[i][j][k][12] + fnew[i][j][k][13] +
                       fnew[i][j][k][14]) +
                  fnew[i][j][k][2] + fnew[i][j][k][4] - seqyy;
              fnew[i][j][k][7] = 0.25 * (rbvw + seqyy - 2.0 * fnew[i][j][k][2]) - fnew[i][j][k][11];
              fnew[i][j][k][8] = 0.25 * (rbvw + seqyy - 2.0 * fnew[i][j][k][2]) - fnew[i][j][k][12];
              fnew[i][j][k][9] =
                  0.25 * (-rbvw + seqyy - 2.0 * fnew[i][j][k][4]) - fnew[i][j][k][13];
              fnew[i][j][k][10] =
                  0.25 * (-rbvw + seqyy - 2.0 * fnew[i][j][k][4]) - fnew[i][j][k][14];

            } else {    // top wall modified return (ori should now be 16 or something is messed up)
              double rb, rbvw, rvw2, seqyy, r;
              rb = fnew[i][j][k][0] + fnew[i][j][k][1] + fnew[i][j][k][2] + fnew[i][j][k][3] +
                  fnew[i][j][k][4] + fnew[i][j][k][5] + fnew[i][j][k][7] + fnew[i][j][k][8] +
                  fnew[i][j][k][9] + fnew[i][j][k][10] + f_lb[i][j][k][5] +
                  (feq[i][j][k][5] - f_lb[i][j][k][5]) / tau + f_lb[i][j][k][7] +
                  (feq[i][j][k][7] - f_lb[i][j][k][7]) / tau + f_lb[i][j][k][8] +
                  (feq[i][j][k][8] - f_lb[i][j][k][8]) / tau + f_lb[i][j][k][9] +
                  (feq[i][j][k][9] - f_lb[i][j][k][9]) / tau + f_lb[i][j][k][10] +
                  (feq[i][j][k][10] - f_lb[i][j][k][10]) / tau;
              rbvw = rb * vwtp;
              rvw2 = rbvw * vwtp;
              seqyy = rb * a_0 + rvw2;
              r = rb / density_lb[i][j][k];

              fnew[i][j][k][2] = r * feq[i][j][k][2];
              fnew[i][j][k][4] = r * feq[i][j][k][4];
              fnew[i][j][k][0] = rb -
                  2.0 *
                      (fnew[i][j][k][5] + fnew[i][j][k][7] + fnew[i][j][k][8] + fnew[i][j][k][9] +
                       fnew[i][j][k][10] + fnew[i][j][k][2] + fnew[i][j][k][4]) +
                  rvw2;
              fnew[i][j][k][1] = 0.5 * (fnew[i][j][k][2] + fnew[i][j][k][4] - rvw2);
              fnew[i][j][k][3] = fnew[i][j][k][1];
              fnew[i][j][k][6] = fnew[i][j][k][5] +
                  2.0 *
                      (fnew[i][j][k][7] + fnew[i][j][k][8] + fnew[i][j][k][9] + fnew[i][j][k][10]) +
                  fnew[i][j][k][2] + fnew[i][j][k][4] - seqyy;
              fnew[i][j][k][11] = 0.25 * (rbvw + seqyy - 2.0 * fnew[i][j][k][2]) - fnew[i][j][k][7];
              fnew[i][j][k][12] = 0.25 * (rbvw + seqyy - 2.0 * fnew[i][j][k][2]) - fnew[i][j][k][8];
              fnew[i][j][k][13] =
                  0.25 * (-rbvw + seqyy - 2.0 * fnew[i][j][k][4]) - fnew[i][j][k][9];
              fnew[i][j][k][14] =
                  0.25 * (-rbvw + seqyy - 2.0 * fnew[i][j][k][4]) - fnew[i][j][k][10];
            }
          } else {
            for (int m = 0; m < 19; m++) {
              int imod = i - e19[m][0];
              int jmod = j - e19[m][1];
              int kmod = k - e19[m][2];

              fnew[i][j][k][m] = f_lb[imod][jmod][kmod][m] +
                  (feq[imod][jmod][kmod][m] - f_lb[imod][jmod][kmod][m]) / tau;
            }
          }
        }
      }
}

//==========================================================================
// Update the distribution functions over the entire simulation domain for
// the D3Q15 model.
//==========================================================================
void FixLbFluid::update_full15()
{
  MPI_Request requests[8];
  int numrequests = 4;

  //--------------------------------------------------------------------------
  // Fixed z boundary conditions.
  //--------------------------------------------------------------------------
  if (domain->periodicity[2] == 0) {

    for (int i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
    MPI_Isend(&feq[1][1][1][0], 1, passxf, comm->procneigh[0][0], 15, world, &requests[0]);
    MPI_Irecv(&feq[0][1][1][0], 1, passxf, comm->procneigh[0][0], 25, world, &requests[1]);
    MPI_Isend(&feq[subNbx - 2][1][1][0], 1, passxf, comm->procneigh[0][1], 25, world, &requests[2]);
    MPI_Irecv(&feq[subNbx - 1][1][1][0], 1, passxf, comm->procneigh[0][1], 15, world, &requests[3]);

    MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

    for (int i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][1][1][0], 1, passyf, comm->procneigh[1][0], 15, world, &requests[0]);
    MPI_Irecv(&feq[0][0][1][0], 1, passyf, comm->procneigh[1][0], 25, world, &requests[1]);
    MPI_Isend(&feq[0][subNby - 2][1][0], 1, passyf, comm->procneigh[1][1], 25, world, &requests[2]);
    MPI_Irecv(&feq[0][subNby - 1][1][0], 1, passyf, comm->procneigh[1][1], 15, world, &requests[3]);

    MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

    for (int i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][0][1][0], 1, passzf, comm->procneigh[2][0], 15, world, &requests[0]);
    MPI_Irecv(&feq[0][0][0][0], 1, passzf, comm->procneigh[2][0], 25, world, &requests[1]);
    MPI_Isend(&feq[0][0][subNbz - 2][0], 1, passzf, comm->procneigh[2][1], 25, world, &requests[2]);
    MPI_Irecv(&feq[0][0][subNbz - 1][0], 1, passzf, comm->procneigh[2][1], 15, world, &requests[3]);

    MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

    update_periodic(2, subNbx - 2, 2, subNby - 2, 2, subNbz - 2);
    update_periodic(1, 2, 2, subNby - 2, 2, subNbz - 2);
    update_periodic(subNbx - 2, subNbx - 1, 2, subNby - 2, 2, subNbz - 2);
    update_periodic(1, subNbx - 1, 1, 2, 2, subNbz - 2);
    update_periodic(1, subNbx - 1, subNby - 2, subNby - 1, 2, subNbz - 2);
    update_periodic(1, subNbx - 1, 1, subNby - 1, 1, 2);
    update_periodic(1, subNbx - 1, 1, subNby - 1, subNbz - 2, subNbz - 1);

    //--------------------------------------------------------------------------
    // Periodic z boundary conditions.
    //--------------------------------------------------------------------------
  } else {

    for (int i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
    MPI_Isend(&feq[1][1][1][0], 1, passxf, comm->procneigh[0][0], 15, world, &requests[0]);
    MPI_Irecv(&feq[0][1][1][0], 1, passxf, comm->procneigh[0][0], 25, world, &requests[1]);
    MPI_Isend(&feq[subNbx - 2][1][1][0], 1, passxf, comm->procneigh[0][1], 25, world, &requests[2]);
    MPI_Irecv(&feq[subNbx - 1][1][1][0], 1, passxf, comm->procneigh[0][1], 15, world, &requests[3]);

    update_periodic(2, subNbx - 2, 2, subNby - 2, 2, subNbz - 2);
    MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

    for (int i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][1][1][0], 1, passyf, comm->procneigh[1][0], 15, world, &requests[0]);
    MPI_Irecv(&feq[0][0][1][0], 1, passyf, comm->procneigh[1][0], 25, world, &requests[1]);
    MPI_Isend(&feq[0][subNby - 2][1][0], 1, passyf, comm->procneigh[1][1], 25, world, &requests[2]);
    MPI_Irecv(&feq[0][subNby - 1][1][0], 1, passyf, comm->procneigh[1][1], 15, world, &requests[3]);

    update_periodic(1, 2, 2, subNby - 2, 2, subNbz - 2);
    update_periodic(subNbx - 2, subNbx - 1, 2, subNby - 2, 2, subNbz - 2);
    MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

    for (int i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][0][1][0], 1, passzf, comm->procneigh[2][0], 15, world, &requests[0]);
    MPI_Irecv(&feq[0][0][0][0], 1, passzf, comm->procneigh[2][0], 25, world, &requests[1]);
    MPI_Isend(&feq[0][0][subNbz - 2][0], 1, passzf, comm->procneigh[2][1], 25, world, &requests[2]);
    MPI_Irecv(&feq[0][0][subNbz - 1][0], 1, passzf, comm->procneigh[2][1], 15, world, &requests[3]);

    update_periodic(1, subNbx - 1, 1, 2, 2, subNbz - 2);
    update_periodic(1, subNbx - 1, subNby - 2, subNby - 1, 2, subNbz - 2);
    MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

    update_periodic(1, subNbx - 1, 1, subNby - 1, 1, 2);
    update_periodic(1, subNbx - 1, 1, subNby - 1, subNbz - 2, subNbz - 1);
  }

  // Rescale for pressure boundary conditions in x-direction
  // modify populations leaving after collision
  if (pressure == 1) {
    if (comm->myloc[0] == 0) {
      for (int j = 1; j < subNby - 1; j++) {
        int jup = j + 1;
        int jdwn = j - 1;
        for (int k = 1; k < subNbz - 1; k++) {
          int kup = k + 1;
          int kdwn = k - 1;
          if (sublattice[1][j][k].type ==
              0) {    // this works ok, but should have some case for walls and fnew from walls.

            fnew[1][j][k][1] = fnew[1][j][k][1] * rhoH / density_lb[0][j][k];
            fnew[1][j][k][7] = fnew[1][j][k][7] * rhoH / density_lb[0][jdwn][kdwn];
            fnew[1][j][k][10] = fnew[1][j][k][10] * rhoH / density_lb[0][jup][kdwn];
            fnew[1][j][k][11] = fnew[1][j][k][11] * rhoH / density_lb[0][jdwn][kup];
            fnew[1][j][k][14] = fnew[1][j][k][14] * rhoH / density_lb[0][jup][kup];
          }
        }
      }
    }
    if (comm->myloc[0] == comm->procgrid[0] - 1) {
      for (int j = 1; j < subNby - 1; j++) {
        int jup = j + 1;
        int jdwn = j - 1;
        for (int k = 1; k < subNbz - 1; k++) {
          int kup = k + 1;
          int kdwn = k - 1;
          if (sublattice[subNbx - 2][j][k].type == 0) {

            fnew[subNbx - 2][j][k][3] =
                fnew[subNbx - 2][j][k][3] * rhoL / density_lb[subNbx - 1][j][k];
            fnew[subNbx - 2][j][k][8] =
                fnew[subNbx - 2][j][k][8] * rhoL / density_lb[subNbx - 1][jdwn][kdwn];
            fnew[subNbx - 2][j][k][9] =
                fnew[subNbx - 2][j][k][9] * rhoL / density_lb[subNbx - 1][jup][kdwn];
            fnew[subNbx - 2][j][k][12] =
                fnew[subNbx - 2][j][k][12] * rhoL / density_lb[subNbx - 1][jdwn][kup];
            fnew[subNbx - 2][j][k][13] =
                fnew[subNbx - 2][j][k][13] * rhoL / density_lb[subNbx - 1][jup][kup];
          }
        }
      }
    }
  }
}

//==========================================================================
// Update the distribution functions over the entire simulation domain for
// the D3Q19 model.
//==========================================================================
void FixLbFluid::update_full19()
{

  MPI_Request req_send15, req_recv15;
  MPI_Request req_send25, req_recv25;
  MPI_Request requests[8];
  int numrequests;
  double tmp1, tmp2, tmp3;
  MPI_Status status;
  double rb;
  int i, j, k;

  numrequests = 4;

  //--------------------------------------------------------------------------
  // Fixed z boundary conditions.
  //--------------------------------------------------------------------------
  if (domain->periodicity[2] == 0) {

    for (i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
    MPI_Isend(&feq[1][1][1][0], 1, passxf, comm->procneigh[0][0], 15, world, &requests[0]);
    MPI_Irecv(&feq[0][1][1][0], 1, passxf, comm->procneigh[0][0], 25, world, &requests[1]);
    MPI_Isend(&feq[subNbx - 2][1][1][0], 1, passxf, comm->procneigh[0][1], 25, world, &requests[2]);
    MPI_Irecv(&feq[subNbx - 1][1][1][0], 1, passxf, comm->procneigh[0][1], 15, world, &requests[3]);

    update_periodic(2, subNbx - 2, 2, subNby - 2, 2, subNbz - 2);
    MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

    for (i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][1][1][0], 1, passyf, comm->procneigh[1][0], 15, world, &requests[0]);
    MPI_Irecv(&feq[0][0][1][0], 1, passyf, comm->procneigh[1][0], 25, world, &requests[1]);
    MPI_Isend(&feq[0][subNby - 2][1][0], 1, passyf, comm->procneigh[1][1], 25, world, &requests[2]);
    MPI_Irecv(&feq[0][subNby - 1][1][0], 1, passyf, comm->procneigh[1][1], 15, world, &requests[3]);

    update_periodic(1, 2, 2, subNby - 2, 2, subNbz - 2);
    update_periodic(subNbx - 2, subNbx - 1, 2, subNby - 2, 2, subNbz - 2);
    MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

    for (i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][0][1][0], 1, passzf, comm->procneigh[2][0], 15, world, &requests[0]);
    MPI_Irecv(&feq[0][0][0][0], 1, passzf, comm->procneigh[2][0], 25, world, &requests[1]);
    MPI_Isend(&feq[0][0][subNbz - 2][0], 1, passzf, comm->procneigh[2][1], 25, world, &requests[2]);
    MPI_Irecv(&feq[0][0][subNbz - 1][0], 1, passzf, comm->procneigh[2][1], 15, world, &requests[3]);

    update_periodic(1, subNbx - 1, 1, 2, 2, subNbz - 2);
    update_periodic(1, subNbx - 1, subNby - 2, subNby - 1, 2, subNbz - 2);
    MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

    update_periodic(1, subNbx - 1, 1, subNby - 1, 1, 2);
    update_periodic(1, subNbx - 1, 1, subNby - 1, subNbz - 2, subNbz - 1);

    req_send15 = MPI_REQUEST_NULL;
    req_recv25 = MPI_REQUEST_NULL;
    req_send25 = MPI_REQUEST_NULL;
    req_recv15 = MPI_REQUEST_NULL;

    if (comm->myloc[2] == 0) {
      MPI_Isend(&fnew[0][0][1][0], 1, passzf, comm->procneigh[2][0], 15, world, &req_send15);
      MPI_Irecv(&fnew[0][0][0][0], 1, passzf, comm->procneigh[2][0], 25, world, &req_recv25);
    }

    if (comm->myloc[2] == comm->procgrid[2] - 1) {
      MPI_Isend(&fnew[0][0][subNbz - 2][0], 1, passzf, comm->procneigh[2][1], 25, world,
                &req_send25);
      MPI_Irecv(&fnew[0][0][subNbz - 1][0], 1, passzf, comm->procneigh[2][1], 15, world,
                &req_recv15);
    }
    if (comm->myloc[2] == 0) {
      MPI_Wait(&req_send15, &status);
      MPI_Wait(&req_recv25, &status);

      for (i = 1; i < subNbx - 1; i++) {
        for (j = 1; j < subNby - 1; j++) {
          k = 1;

          fnew[i][j][k][5] = fnew[i][j][k - 1][6];
          tmp1 = fnew[i][j][k - 1][12] + fnew[i][j][k - 1][14] + fnew[i][j][k - 1][16] +
              fnew[i][j][k - 1][18];

          tmp2 = fnew[i][j][k][3] + fnew[i][j][k][9] + fnew[i][j][k][10] + fnew[i][j][k][14] -
              fnew[i][j][k][1] - fnew[i][j][k][7] - fnew[i][j][k][8] - fnew[i][j][k][12];

          rb = fnew[i][j][k][0] + fnew[i][j][k][1] + fnew[i][j][k][2] + fnew[i][j][k][3] +
              fnew[i][j][k][4] + fnew[i][j][k][5] + fnew[i][j][k][6] + fnew[i][j][k][7] +
              fnew[i][j][k][8] + fnew[i][j][k][9] + fnew[i][j][k][10] + fnew[i][j][k][12] +
              fnew[i][j][k][14] + fnew[i][j][k][16] + fnew[i][j][k][18] + tmp1;

          tmp3 = rb * vwbt - fnew[i][j][k][2] + fnew[i][j][k][4] - fnew[i][j][k][7] +
              fnew[i][j][k][8] - fnew[i][j][k][9] + fnew[i][j][k][10] - fnew[i][j][k][16] +
              fnew[i][j][k][18];

          fnew[i][j][k][11] = 0.25 * (tmp1 + 2.0 * tmp2);
          fnew[i][j][k][13] = 0.25 * (tmp1 - 2.0 * tmp2);
          fnew[i][j][k][15] = 0.25 * (tmp1 + 2.0 * tmp3);
          fnew[i][j][k][17] = 0.25 * (tmp1 - 2.0 * tmp3);
        }
      }
    }
    if (comm->myloc[2] == comm->procgrid[2] - 1) {
      MPI_Wait(&req_send25, &status);
      MPI_Wait(&req_recv15, &status);

      for (i = 1; i < subNbx - 1; i++) {
        for (j = 1; j < subNby - 1; j++) {
          k = subNbz - 2;

          fnew[i][j][k][6] = fnew[i][j][k + 1][5];
          tmp1 = fnew[i][j][k + 1][11] + fnew[i][j][k + 1][13] + fnew[i][j][k + 1][15] +
              fnew[i][j][k + 1][17];

          tmp2 = fnew[i][j][k][3] + fnew[i][j][k][9] + fnew[i][j][k][10] + fnew[i][j][k][13] -
              fnew[i][j][k][1] - fnew[i][j][k][7] - fnew[i][j][k][8] - fnew[i][j][k][11];

          rb = fnew[i][j][k][0] + fnew[i][j][k][1] + fnew[i][j][k][2] + fnew[i][j][k][3] +
              fnew[i][j][k][4] + fnew[i][j][k][5] + fnew[i][j][k][6] + fnew[i][j][k][7] +
              fnew[i][j][k][8] + fnew[i][j][k][9] + fnew[i][j][k][10] + fnew[i][j][k][11] +
              fnew[i][j][k][13] + fnew[i][j][k][15] + fnew[i][j][k][17] + tmp1;

          tmp3 = rb * vwtp - fnew[i][j][k][2] + fnew[i][j][k][4] - fnew[i][j][k][7] +
              fnew[i][j][k][8] - fnew[i][j][k][9] + fnew[i][j][k][10] - fnew[i][j][k][15] +
              fnew[i][j][k][17];

          fnew[i][j][k][12] = 0.25 * (tmp1 + 2.0 * tmp2);
          fnew[i][j][k][14] = 0.25 * (tmp1 - 2.0 * tmp2);
          fnew[i][j][k][16] = 0.25 * (tmp1 + 2.0 * tmp3);
          fnew[i][j][k][18] = 0.25 * (tmp1 - 2.0 * tmp3);
        }
      }
    }

    //--------------------------------------------------------------------------
    // Periodic z boundary conditions.
    //--------------------------------------------------------------------------
  } else {

    for (i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
    MPI_Isend(&feq[1][1][1][0], 1, passxf, comm->procneigh[0][0], 15, world, &requests[0]);
    MPI_Irecv(&feq[0][1][1][0], 1, passxf, comm->procneigh[0][0], 25, world, &requests[1]);
    MPI_Isend(&feq[subNbx - 2][1][1][0], 1, passxf, comm->procneigh[0][1], 25, world, &requests[2]);
    MPI_Irecv(&feq[subNbx - 1][1][1][0], 1, passxf, comm->procneigh[0][1], 15, world, &requests[3]);

    update_periodic(2, subNbx - 2, 2, subNby - 2, 2, subNbz - 2);
    MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

    for (i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][1][1][0], 1, passyf, comm->procneigh[1][0], 15, world, &requests[0]);
    MPI_Irecv(&feq[0][0][1][0], 1, passyf, comm->procneigh[1][0], 25, world, &requests[1]);
    MPI_Isend(&feq[0][subNby - 2][1][0], 1, passyf, comm->procneigh[1][1], 25, world, &requests[2]);
    MPI_Irecv(&feq[0][subNby - 1][1][0], 1, passyf, comm->procneigh[1][1], 15, world, &requests[3]);

    update_periodic(1, 2, 2, subNby - 2, 2, subNbz - 2);
    update_periodic(subNbx - 2, subNbx - 1, 2, subNby - 2, 2, subNbz - 2);
    MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

    for (i = 0; i < numrequests; i++) requests[i] = MPI_REQUEST_NULL;
    MPI_Isend(&feq[0][0][1][0], 1, passzf, comm->procneigh[2][0], 15, world, &requests[0]);
    MPI_Irecv(&feq[0][0][0][0], 1, passzf, comm->procneigh[2][0], 25, world, &requests[1]);
    MPI_Isend(&feq[0][0][subNbz - 2][0], 1, passzf, comm->procneigh[2][1], 25, world, &requests[2]);
    MPI_Irecv(&feq[0][0][subNbz - 1][0], 1, passzf, comm->procneigh[2][1], 15, world, &requests[3]);

    update_periodic(1, subNbx - 1, 1, 2, 2, subNbz - 2);
    update_periodic(1, subNbx - 1, subNby - 2, subNby - 1, 2, subNbz - 2);
    MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);

    update_periodic(1, subNbx - 1, 1, subNby - 1, 1, 2);
    update_periodic(1, subNbx - 1, 1, subNby - 1, subNbz - 2, subNbz - 1);
  }

  // Rescale for pressure boundary conditions in x-direction
  // modify populations leaving after collision
  if (pressure == 1) {
    if (comm->myloc[0] == 0) {
      for (j = 1; j < subNby - 1; j++) {
        int jup = j + 1;
        int jdwn = j - 1;
        for (k = 1; k < subNbz - 1; k++) {
          int kup = k + 1;
          int kdwn = k - 1;

          fnew[1][j][k][1] = fnew[1][j][k][1] * rhoH / density_lb[0][j][k];
          fnew[1][j][k][7] = fnew[1][j][k][7] * rhoH / density_lb[0][jdwn][k];
          fnew[1][j][k][8] = fnew[1][j][k][8] * rhoH / density_lb[0][jup][k];
          fnew[1][j][k][11] = fnew[1][j][k][11] * rhoH / density_lb[0][j][kdwn];
          fnew[1][j][k][12] = fnew[1][j][k][12] * rhoH / density_lb[0][j][kup];
        }
      }
    }

    if (comm->myloc[0] == comm->procgrid[0] - 1) {
      for (j = 1; j < subNby - 1; j++) {
        int jup = j + 1;
        int jdwn = j - 1;
        for (k = 1; k < subNbz - 1; k++) {
          int kup = k + 1;
          int kdwn = k - 1;

          fnew[subNbx - 2][j][k][3] =
              fnew[subNbx - 2][j][k][3] * rhoL / density_lb[subNbx - 1][j][k];
          fnew[subNbx - 2][j][k][9] =
              fnew[subNbx - 2][j][k][9] * rhoL / density_lb[subNbx - 1][jdwn][k];
          fnew[subNbx - 2][j][k][10] =
              fnew[subNbx - 2][j][k][10] * rhoL / density_lb[subNbx - 1][jup][k];
          fnew[subNbx - 2][j][k][13] =
              fnew[subNbx - 2][j][k][13] * rhoL / density_lb[subNbx - 1][j][kdwn];
          fnew[subNbx - 2][j][k][14] =
              fnew[subNbx - 2][j][k][14] * rhoL / density_lb[subNbx - 1][j][kup];
        }
      }
    }
  }
}

/* nanopit routines */

void FixLbFluid::initializeGlobalGeometry()
{
  // --------------------------------------------
  // initializeGlobalGeometry defines the global
  // lattice structure through the lattice array
  // --------------------------------------------

  if (npits > -1) {

    // Add topography block by block along x-direction.
    int ix = 0;
    if (npits == 0) {
      addslit(ix, h_s, h_p, Nbx - 1, sw);    // two parallel walls
    } else
      addslit(ix, h_s, h_p, l_e, sw);    // if l_e = 0, does not add anything

    for (int nn = 0; nn < npits; nn++) {
      addpit(ix, h_s, h_p, w_p, l_p, sw);    // if l_p=Nbx, no corners/edges along z axis
      if (nn < npits - 1) { addslit(ix, h_s, h_p, l_pp - 1, sw); }
    }

    if (npits > 0 && ix < Nbx - 1) {
      addslit(ix, h_s, h_p, l_e - 1, sw);    //***should probably make this one go to Nbx
    }

  } else {    // no pit geometry so all sites active if periodic or zwalls.
    for (int i = 0; i < Nbx; i++)
      for (int j = 0; j < Nby; j++) {
        if (domain->periodicity[2] == 0) {    // If there are walls in the z-direction
          wholelattice[i][j][0].type = wholelattice[i][j][Nbz - 1].type = 3;
          wholelattice[i][j][0].orientation = 3;
          wholelattice[i][j][Nbz - 1].orientation = 16;
        } else {
          wholelattice[i][j][0].type = wholelattice[i][j][Nbz - 1].type = 0;
          wholelattice[i][j][0].orientation = wholelattice[i][j][Nbz - 1].orientation = 0;
        }
        for (int k = 1; k < Nbz - 1; k++) {
          wholelattice[i][j][k].type = 0;
          wholelattice[i][j][k].orientation = 0;
        }
      }
  }

  int inletsites = 0, outletsites = 0;
  for (int j = 0; j < Nby; j++)
    for (int k = 0; k < Nbz; k++) {
      if (wholelattice[0][j][k].type != 2) inletsites++;
      if (wholelattice[Nbx - 1][j][k].type != 2) outletsites++;
    }
  if (comm->me == 0)
    if (inletsites != outletsites)
      error->all(FLERR, "inlet and outlet geometry do not match in x-direction");
  openingsites = inletsites;
}

void FixLbFluid::initializeGeometry()
{
  // --------------------------------------------
  //  initializeGeometry defines the local lattice
  //  structure through the sublattice array
  // --------------------------------------------
#define GEODUMP 0

  int i, j, k;
  int boxlims[3][2];    // absolute subvolume coordinates

  // determine absolute subvolume coordinates
  boxlims[0][0] = comm->myloc[0] * (subNbx - 2);                       // inclusive limits!
  boxlims[0][1] = comm->myloc[0] * (subNbx - 2) + (subNbx - 2 - 1);    // inclusive limits!
  boxlims[1][0] = comm->myloc[1] * (subNby - 2);                       // inclusive limits!
  boxlims[1][1] = comm->myloc[1] * (subNby - 2) + (subNby - 2 - 1);    // inclusive limits!
  if (domain->periodicity[2] == 0 && comm->myloc[2] == comm->procgrid[2] - 1) {
    boxlims[2][0] = comm->myloc[2] * (subNbz - 3);                       // inclusive limits!
    boxlims[2][1] = comm->myloc[2] * (subNbz - 3) + (subNbz - 2 - 1);    // inclusive limits!
  } else {
    boxlims[2][0] = comm->myloc[2] * (subNbz - 2);                       // inclusive limits!
    boxlims[2][1] = comm->myloc[2] * (subNbz - 2) + (subNbz - 2 - 1);    // inclusive limits!
  }

  for (i = 0; i < subNbx; i++)    // Sweep over the sublattice including the buffer zone.
    for (j = 0; j < subNby; j++)
      for (k = 0; k < subNbz; k++) {    // Buffers set assuming we have periodic boundary conditions
        sublattice[i][j][k] =
            wholelattice[(boxlims[0][0] + i - 1 + Nbx) % Nbx][(boxlims[1][0] + j - 1 + Nby) % Nby]
                        [(boxlims[2][0] + k - 1 + Nbz) % Nbz];
        if (k == 0 &&
            comm->myloc[2] == 0) {    // true boundary ghost rather than neighboring processor ghost
          Site just_inside =
              wholelattice[(boxlims[0][0] + i - 1 + Nbx) % Nbx][(boxlims[1][0] + j - 1 + Nby) % Nby]
                          [(boxlims[2][0] + k + 1 - 1 + Nbz) % Nbz];
          if (sublattice[i][j][k].type == 1 && just_inside.type == 2) {
            //we have incorrectly wrapped around a fixed boundary that doesn't bound any actual fluid on our side
            sublattice[i][j][k].type = 2;
            sublattice[i][j][k].orientation = 0;
            ;
          }
        } else if (k == subNbz - 1 && comm->myloc[2] == comm->procgrid[2] - 1) {
          Site just_inside =
              wholelattice[(boxlims[0][0] + i - 1 + Nbx) % Nbx][(boxlims[1][0] + j - 1 + Nby) % Nby]
                          [(boxlims[2][0] + k - 1 - 1 + Nbz) % Nbz];
          if (sublattice[i][j][k].type == 1 && just_inside.type == 2) {
            //we have incorrectly wrapped around a fixed boundary that doesn't bound any actual fluid on our side
            sublattice[i][j][k].type = 2;
            sublattice[i][j][k].orientation = 0;
            ;
          }
        }
      }

      // Output local geometry to data files labeled with processor number
      // Type dump
#if GEODUMP
  auto datfile = fmt::format("subgeom_{}_end_type.dmp", me);
  FILE *outfile = fopen(datfile.c_str(), "w");
  if (!outfile)
    error->one(FLERR, " file {} could not be opened: {}", datfile, utils::getsyserror());

  fmt::print(outfile, "\n me: {} px: {} py: {} pz: {}\n", me, comm->myloc[0], comm->myloc[1],
             comm->myloc[2]);

  for (i = 0; i < subNbx; i++) {
    fmt::print(outfile, "i={}\n", i);
    for (k = subNbz - 1; k > -1; k--) {
      if (k == subNbz - 2 || k == 0) {
        for (j = 0; j < subNby + 2; j++) fputs("---", outfile);
        fputs("\n", outfile);
      }
      for (j = 0; j < subNby; j++) {
        fmt::print(outfile, " {} ", sublattice[i][j][k].type);
        if (j == 0 || j == subNby - 2) fputs(" | ", outfile);
        if (j == subNby - 1) fputs("\n", outfile);
      }
    }
    fputs(" \n \n", outfile);
  }
  fputs("\n", outfile);
  fclose(outfile);

  // Orientation dump
  datfile = fmt::format("subgeom_{}_end_ori.dmp", me);
  outfile = fopen(datfile.c_str(), "w");

  if (!outfile)
    error->one(FLERR, " file {} could not be opened: {}", datfile, utils::getsyserror());

  fmt::print("\nme: {}\n", me);
  for (i = 0; i < subNbx; i++) {
    fmt::print("i={}\n", i);
    for (k = subNbz - 1; k > -1; k--) {
      if (k == subNbz - 2 || k == 0) {
        for (j = 0; j < subNby + 2; j++) fputs("---", outfile);
        fputs("\bn", outfile);
      }
      for (j = 0; j < subNby; j++) {
        fmt::print(outfile, " {} ", sublattice[i][j][k].orientation);
        if (j == 0 || j == subNby - 2) fputs(" | ", outfile);
        if (j == subNby - 1) fputs("\n", outfile);
      }
    }
    fputs(" \n \n", outfile);
  }
  fputs("\n", outfile);
  fclose(outfile);
#endif
}

void FixLbFluid::addslit(int &x0, const int HS, const int HP, const int LE, const int SW)
{
  // preconditions: if LE = 0, don't add any slit region
  // x0 : leftmost coordinate
  // HS : height of slit
  // HP : height of pit
  // LE : length of end segment
  // SW : side walls on/off

  Site temp;
  int xend;

  if (x0 + LE < Nbx)
    xend = LE;
  else
    xend = Nbx - x0;

  for (int ii = 0; ii < xend; ii++) {
    for (int jj = 0; jj < Nby; jj++) {
      for (int kk = 0; kk < (HS + HP + 1); kk++) {
        temp.type = 2;
        temp.orientation = 0;
        wholelattice[x0 + ii][jj][kk] = temp;
      }
    }
  }

  for (int ii = 0; ii < xend; ii++) {
    for (int jj = 0; jj < Nby; jj++) {
      for (int kk = HP; kk < (HS + HP + 1); kk++) {
        if (kk == HP) {    // bottom of slit & bottom edges
          temp.type = 1;
          temp.orientation = 3;
          if (SW == 1 && jj == 0) {
            temp.orientation = 6;
          } else if (SW == 1 && jj == Nby - 1) {
            temp.orientation = 19;
          }
        } else if (kk == (HP + HS)) {    // slit ceiling & top edges
          temp.type = 1;
          temp.orientation = 16;
          if (SW == 1 && jj == 0) {
            temp.orientation = 27;
          } else if (SW == 1 && jj == Nby - 1) {
            temp.orientation = 28;
          }
        } else {    // bulk nodes & side walls
          temp.type = 0;
          temp.orientation = 0;
          if (SW == 1 && jj == 0) {
            temp.orientation = 2;
            temp.type = 1;
          } else if (SW == 1 && jj == Nby - 1) {
            temp.orientation = 15;
            temp.type = 1;
          }
        }
        wholelattice[x0 + ii][jj][kk] = temp;
      }
      temp.type = 2;
      temp.orientation = 0;
      wholelattice[x0 + ii][jj][HP + HS + 1] = temp;
    }
  }
  x0 += xend;
  //x0 += LE;
}

void FixLbFluid::addpit(int &x0, const int HS, const int HP, const int WP, const int LP,
                        const int SW)
{
  // This routine works but is in bad need of a complete revamp to clean it up.
  // preconditions: if x0 == 0 and x0+l_p >= Nbx, then do not add walls at x0 and x0+LP
  // x0 : leftmost coordinate
  // HS : height of slit
  // HP : height of pit
  // WP : width of pit
  // LP : length of pit

  Site temp;
  int ii, jj, kk;
  int xend;

  int p_edge1y = (Nby - WP) / 2;
  int p_edge2y = p_edge1y + WP - 1;

  if (SW == 1) {
    p_edge1y = 0;
    p_edge2y = Nby - 1;
  }

  if (SW == 0) {
    if (p_edge1y <= 0) p_edge1y = -1;
    if (p_edge2y >= Nby - 1) p_edge2y = Nby;
  }

  if (x0 + LP < Nbx)
    xend = x0 + LP;
  else
    xend = Nbx - 1;

  for (int ii = x0; ii <= xend; ii++) {    //Initialize whole block to "not in fluid"
    for (int jj = 0; jj < Nby; jj++) {
      for (int kk = 0; kk < (HS + HP + 1); kk++) {
        temp.type = 2;
        temp.orientation = 0;
        wholelattice[ii][jj][kk] = temp;
      }
    }
  }
  // There is some redundancy in the checks, should be cleaned up at some point, but get "vertical chute" working first
  if (SW == 0) {
    for (int ii = x0; ii <= xend; ii++) {
      for (int jj = 0; jj < Nby; jj++) {
        for (int kk = 0; kk < (HS + HP + 1); kk++) {

          if ((jj < p_edge1y || jj > p_edge2y) && kk < HP) {
            temp.orientation = 0;
            temp.type = 2;
          }

          // top plate
          else if (kk == HP + HS && HP < (Nbz - 1)) {
            temp.orientation = 16;
            temp.type = 1;
          }

          // slit regions: bottom plate
          else if ((jj < p_edge1y || jj > p_edge2y) && kk == HP) {
            temp.orientation = 3;
            temp.type = 1;
          }

          // bottom of pit, ok for troughs, too
          else if ((ii > x0 && ii < LP + x0) && (jj > p_edge1y && jj < p_edge2y) && kk == 0 &&
                   HP < (Nbz - 1)) {
            temp.orientation = 3;
            temp.type = 1;
          }

          // top corners at the level of slit floor WHEN THERE ARE NO SIDE WALLS
          else if (ii == x0 && jj == p_edge1y && p_edge1y > 0 && kk == HP) {
            if (x0 > 0 && LP < Nbx) {
              temp.orientation = 10;
              temp.type = 1;
            } else {
              temp.orientation = 7;
              temp.type = 1;
            }
          } else if (ii == (x0 + LP) && jj == p_edge1y && p_edge1y > 0 && kk == HP) {
            if ((x0 + LP) < (Nbx - 1) && LP < Nbx) {
              temp.orientation = 23;
              temp.type = 1;
            } else {
              temp.orientation = 7;
              temp.type = 1;
            }
          } else if (ii == x0 && jj == p_edge2y && (p_edge2y < Nby) && kk == HP) {
            if (x0 > 0 && LP < Nbx) {
              temp.orientation = 11;
              temp.type = 1;
            } else {
              temp.orientation = 20;
              temp.type = 1;
            }
          } else if (ii == (x0 + LP) && jj == p_edge2y && (p_edge2y < Nby) && kk == HP) {
            if ((x0 + LP) < (Nbx - 1) && LP < Nbx) {
              temp.orientation = 24;
              temp.type = 1;
            } else {
              temp.orientation = 20;
              temp.type = 1;
            }
          }

          // UPSTREAM corners at the level of pit floor
          else if (ii == x0 && jj == p_edge1y && (p_edge1y > 0) && kk == 0) {
            if (x0 > 0 && LP < Nbx && HP < (Nbz - 1))
              temp.orientation = 25;    // corner whose "normal" pointing to fluid is  (+1,+1,+1)
            else if (x0 > 0 && LP < Nbx && HP >= (Nbz - 1))
              temp.orientation =
                  22;    // VERTICAL SHUTE inside edge along z axis "normal" pointing to fluid is  (+1,+1,0)
            else
              temp.orientation = 6;
            temp.type = 1;
          }
          // DOWNSTREAM corners at the level of pit floor
          else if (ii == (x0 + LP) && jj == p_edge1y && (p_edge1y > 0) && kk == 0) {
            if ((x0 + LP) < (Nbx - 1) && LP < Nbx && HP < (Nbz - 1))
              temp.orientation = 12;    // ***redundant check
            else if ((x0 + LP) < (Nbx - 1) && LP < Nbx && HP >= (Nbz - 1))
              temp.orientation = 9;    // VERTICAL SHUTE: HP >= Nbz-1 i.e. height of pit
            else
              temp.orientation = 6;
            temp.type = 1;
          }
          // UPSTREAM corners at the level of pit floor
          else if (ii == x0 && jj == p_edge2y && (p_edge2y < Nby) && kk == 0) {
            if (x0 > 0 && LP < Nbx && HP < (Nbz - 1))
              temp.orientation = 26;
            else if (x0 > 0 && LP < Nbx && HP >= (Nbz - 1))
              temp.orientation = 21;
            else
              temp.orientation = 19;
            temp.type = 1;
          }
          // DOWNSTREAM corners at the level of pit floor
          else if (ii == (x0 + LP) && jj == p_edge2y && (p_edge2y < Nby) && kk == 0) {
            if ((x0 + LP) < (Nbx - 1) && LP < Nbx && HP < (Nbz - 1))
              temp.orientation = 13;
            else if ((x0 + LP) < (Nbx - 1) && LP < Nbx && HP >= (Nbz - 1))
              temp.orientation = 8;
            else
              temp.orientation = 19;
            temp.type = 1;
          }

          // side walls + edges
          else if (ii == x0 && (jj > p_edge1y) && (jj < p_edge2y)) {
            if (x0 > 0 && (HP > (Nbz - 1)) && LP < Nbx && kk < HP) {
              temp.orientation = 1;
              temp.type = 1;
            }    // wall of VERTICAL SHUTE with normal (+1,0,0)
            else if (x0 > 0 && (kk < HP) && (LP < Nbx) && kk > 0) {
              temp.orientation = 1;
              temp.type = 1;
            }                      // UPSTREAM WALL with normal (+1,0,0)
            else if (kk == 0) {    // bottom edge along y
              if (x0 > 0 && LP < Nbx) {
                temp.orientation = 5;
                temp.type = 1;
              } else {
                temp.orientation = 3;
                temp.type = 1;
              }
            } else if (kk == HP) {    // bottom edge along y
              if (x0 > 0 && LP < Nbx) {
                temp.orientation = 4;
                temp.type = 1;
              } else {
                temp.orientation = 0;
                temp.type = 0;
              }    // ***redundant
            } else {
              temp.orientation = 0;
              temp.type = 0;
            }
          } else if (ii == (x0 + LP) && (jj > p_edge1y) && (jj < p_edge2y)) {
            if ((x0 + LP) < (Nbx - 1) && LP < Nbx && kk < HP && HP > (Nbz - 1)) {
              temp.orientation = 14;
              temp.type = 1;
            } else if ((x0 + LP) < (Nbx - 1) && LP < Nbx && kk < HP && kk > 0) {
              temp.orientation = 14;
              temp.type = 1;
            } else if (kk == 0) {    // bottom edge along y
              if ((x0 + LP) < (Nbx - 1) && LP < Nbx) {
                temp.orientation = 18;
                temp.type = 1;
              } else {
                temp.orientation = 3;
                temp.type = 1;
              }
            } else if (kk == HP) {
              if (LP < Nbx) {
                temp.orientation = 17;
                temp.type = 1;
              } else {
                temp.orientation = 0;
                temp.type = 0;
              }
            } else {
              temp.orientation = 0;
              temp.type = 0;
            }
          } else if (ii > x0 && (ii < (x0 + LP)) && jj == p_edge1y && (p_edge1y > 0) && kk < HP) {
            if (HP >= (Nbz - 1)) {
              temp.orientation = 2;
              temp.type = 1;
            } else if (kk > 0 && HP < (Nbz - 1)) {
              temp.orientation = 2;
              temp.type = 1;
            } else if (kk == 0) {
              temp.orientation = 6;
              temp.type = 1;
            }    // edge along x axis
            else {
              temp.orientation = 0;
              temp.type = 0;
            }
          } else if (ii > x0 && (ii < (x0 + LP)) && jj == p_edge2y && (p_edge2y < Nby) && kk < HP) {
            if (HP >= (Nbz - 1)) {
              temp.orientation = 15;
              temp.type = 1;
            } else if (kk > 0 && HP < (Nbz - 1)) {
              temp.orientation = 15;
              temp.type = 1;
            } else if (kk == 0) {
              temp.orientation = 19;
              temp.type = 1;
            }    // edge along x axis
            else {
              temp.orientation = 0;
              temp.type = 0;
            }
          }

          // top edges along x axis. If sidewalls => they become part of flat wall with (0,+/-1,0) surface normals
          else if (ii > x0 && ii < (x0 + LP) && jj == p_edge1y && kk == HP) {
            temp.orientation = 7;
            temp.type = 1;
          } else if (ii > x0 && ii < (x0 + LP) && jj == p_edge2y && kk == HP) {
            temp.orientation = 20;
            temp.type = 1;
          }

          // edges along z axis
          else if (ii == x0 && jj == p_edge1y && kk < HP) {
            if (x0 > 0 && LP < Nbx && HP >= (Nbz - 1)) {
              temp.orientation = 22;
              temp.type = 1;
            } else if (x0 > 0 && LP < Nbx && kk > 0 && HP < (Nbz - 1)) {
              temp.orientation = 22;
              temp.type = 1;
            } else {
              temp.orientation = 2;
              temp.type = 1;
            }
          } else if (ii == x0 && jj == p_edge2y && kk < HP) {
            if (x0 > 0 && LP < Nbx && HP >= (Nbz - 1)) {
              temp.orientation = 21;
              temp.type = 1;
            } else if (x0 > 0 && LP < Nbx && kk > 0 && HP < (Nbz - 1)) {
              temp.orientation = 21;
              temp.type = 1;
            } else {
              temp.orientation = 15;
              temp.type = 1;
            }
          } else if (ii == (x0 + LP) && jj == p_edge1y && kk < HP) {
            if ((x0 + LP) < (Nbx - 1) && LP < Nbx && HP >= (Nbz - 1)) {
              temp.orientation = 9;
              temp.type = 1;
            } else if ((x0 + LP) < (Nbx - 1) && LP < Nbx && kk > 0 && HP < (Nbz - 1)) {
              temp.orientation = 9;
              temp.type = 1;
            } else {
              temp.orientation = 2;
              temp.type = 1;
            }
          } else if (ii == (x0 + LP) && jj == p_edge2y && kk < HP) {
            if ((x0 + LP) < (Nbx - 1) && LP < Nbx && HP >= (Nbz - 1)) {
              temp.orientation = 8;
              temp.type = 1;
            } else if ((x0 + LP) < (Nbx - 1) && LP < Nbx && kk > 0 && HP < (Nbz - 1)) {
              temp.orientation = 8;
              temp.type = 1;
            } else {
              temp.orientation = 15;
              temp.type = 1;
            }
          } else {
            temp.orientation = 0;
            temp.type = 0;
          }
          wholelattice[ii][jj][kk] = temp;
        }
      }
    }
  }

  // PIT REGION IN PRESENCE OF SIDE WALLS IN SYSTEM
  else if (SW == 1) {
    for (ii = x0; ii <= xend; ii++) {
      for (jj = 0; jj < Nby; jj++) {
        for (kk = 0; kk < (HS + HP + 1); kk++) {

          // ceiling & edges
          if (kk == HP + HS && HP < (Nbz - 1)) {
            if (jj == p_edge1y) {
              temp.orientation = 27;
              temp.type = 1;
            } else if (jj == p_edge2y) {
              temp.orientation = 28;
              temp.type = 1;
            } else {
              temp.orientation = 16;
              temp.type = 1;
            }
          } else if (kk > HP && kk < HP + HS) {
            if (jj == p_edge1y) {
              temp.orientation = 2;
              temp.type = 1;
            } else if (jj == p_edge2y) {
              temp.orientation = 15;
              temp.type = 1;
            } else {
              temp.orientation = 0;
              temp.type = 0;
            }
          }

          // bottom of pit, ok for troughs, too
          else if ((ii > x0 && ii < LP + x0) && (jj > p_edge1y && jj < p_edge2y) && kk == 0 &&
                   HP < (Nbz - 1)) {
            temp.orientation = 3;
            temp.type = 1;
          }

          // top corners at the level of slit floor WHEN THERE ARE SIDE WALLS
          else if (ii == x0 && jj == p_edge1y && kk == HP) {
            if (x0 > 0 && LP < Nbx) {
              temp.orientation = 29;
              temp.type = 1;
            } else {
              temp.orientation = 2;
              temp.type = 1;
            }
          } else if (ii == (x0 + LP) && jj == p_edge1y && kk == HP) {
            if ((x0 + LP - 1) < (Nbx - 1) && LP < Nbx) {
              temp.orientation = 31;
              temp.type = 1;
            }
            //	    else { temp.orientation = 6; temp.type = 1; }
          } else if (ii == x0 && jj == p_edge2y && kk == HP) {
            if (x0 > 0 && LP < Nbx) {
              temp.orientation = 30;
              temp.type = 1;
            } else {
              temp.orientation = 15;
              temp.type = 1;
            }
          } else if (ii == (x0 + LP) && jj == p_edge2y && kk == HP) {
            if ((x0 + LP - 1) < (Nbx - 1) && LP < Nbx) {
              temp.orientation = 32;
              temp.type = 1;
            }
            //	    else { temp.orientation = 19; temp.type = 1; }
          }

          // UPSTREAM corners at the bottom of the pit
          else if (ii == x0 && jj == p_edge1y && kk == 0) {
            if (x0 > 0 && LP < Nbx)
              temp.orientation = 25;    // corner whose "normal" pointing to fluid is  (+1,+1,+1)
            else
              temp.orientation = 6;
            temp.type = 1;
          }
          // DOWNSTREAM corners at the bottom of the pit
          else if (ii == (x0 + LP) && jj == p_edge1y && kk == 0) {
            if ((x0 + LP - 1) < (Nbx - 1) && LP < Nbx) temp.orientation = 12;
            temp.type = 1;
          }
          // UPSTREAM corners at the bottom of the pit
          else if (ii == x0 && jj == p_edge2y && kk == 0) {
            if (x0 > 0 && LP < Nbx)
              temp.orientation = 26;
            else
              temp.orientation = 19;
            temp.type = 1;
          }
          // DOWNSTREAM corners at the bottom of the pit
          else if (ii == (x0 + LP) && jj == p_edge2y && kk == 0) {
            if ((x0 + LP - 1) < (Nbx - 1) && LP < Nbx) temp.orientation = 13;
            temp.type = 1;
          }

          // UPSTREAM and DOWNSTREAM pit walls + edges
          else if (ii == x0 && (jj > p_edge1y) && (jj < p_edge2y)) {
            if (x0 > 0 && (kk < HP) && (LP < Nbx) && kk > 0) {
              temp.orientation = 1;
              temp.type = 1;
            }                      // UPSTREAM WALL with normal (+1,0,0)
            else if (kk == 0) {    // bottom edge along y
              if (x0 > 0 && LP < Nbx) {
                temp.orientation = 5;
                temp.type = 1;
              } else {
                temp.orientation = 3;
                temp.type = 1;
              }
            } else if (kk == HP) {    // bottom edge along y
              if (x0 > 0 && LP < Nbx) {
                temp.orientation = 4;
                temp.type = 1;
              } else {
                temp.orientation = 0;
                temp.type = 0;
              }
            } else {
              temp.orientation = 0;
              temp.type = 0;
            }
          } else if (ii == (x0 + LP) && (jj > p_edge1y) && (jj < p_edge2y)) {
            if ((x0 + LP) < (Nbx - 1) && LP < Nbx && kk < HP && kk > 0) {
              temp.orientation = 14;
              temp.type = 1;
            } else if (kk == 0) {    // bottom edge along y
              if ((x0 + LP) < (Nbx - 1) && LP < Nbx) {
                temp.orientation = 18;
                temp.type = 1;
              } else {
                temp.orientation = 3;
                temp.type = 1;
              }
            } else if (kk == HP) {
              if (LP < Nbx) {
                temp.orientation = 17;
                temp.type = 1;
              } else {
                temp.orientation = 0;
                temp.type = 0;
              }
            } else {
              temp.orientation = 0;
              temp.type = 0;
            }
          }

          // SIDE WALLS INSIDE PIT
          else if (ii > x0 && (ii < (x0 + LP)) && jj == p_edge1y) {
            if (kk > 0 && kk < HP + HS) {
              temp.orientation = 2;
              temp.type = 1;
            }    // wall, surface normal (0,+1,0)
            else if (kk == 0) {
              temp.orientation = 6;
              temp.type = 1;
            }    // edge along x axis
            else {
              temp.orientation = 0;
              temp.type = 0;
            }
          } else if (ii > x0 && (ii < (x0 + LP)) && jj == p_edge2y) {
            if (kk > 0 && kk < HP + HS) {
              temp.orientation = 15;
              temp.type = 1;
            }    // wall, surface normal (0,-1,0)
            else if (kk == 0) {
              temp.orientation = 19;
              temp.type = 1;
            }    // edge along x axis
            else {
              temp.orientation = 0;
              temp.type = 0;
            }
          }

          // edges along z axis
          else if (ii == x0 && jj == p_edge1y && kk < HP) {
            if (x0 > 0 && LP < Nbx && kk > 0) {
              temp.orientation = 22;
              temp.type = 1;
            } else {
              temp.orientation = 2;
              temp.type = 1;
            }
          } else if (ii == x0 && jj == p_edge2y && kk < HP) {
            if (x0 > 0 && LP < Nbx && kk > 0) {
              temp.orientation = 21;
              temp.type = 1;
            } else {
              temp.orientation = 15;
              temp.type = 1;
            }
          } else if (ii == (x0 + LP) && jj == p_edge1y && kk < HP) {
            if ((x0 + LP) < (Nbx - 1) && LP < Nbx && kk > 0) {
              temp.orientation = 9;
              temp.type = 1;
            }
            //	    else { temp.orientation = 2; temp.type = 1;}
          } else if (ii == (x0 + LP) && jj == p_edge2y && kk < HP) {
            if ((x0 + LP) < (Nbx - 1) && LP < Nbx && kk > 0) {
              temp.orientation = 8;
              temp.type = 1;
            }
            //	    else { temp.orientation = 15; temp.type = 1; }
          } else {
            temp.orientation = 0;
            temp.type = 0;
          }
          wholelattice[ii][jj][kk] = temp;
        }
      }
    }
  }
  x0 += (xend - x0) + 1;
}

/* nanopit routines end */

//==========================================================================
//   Compute total mass and momentum of the fluid
//==========================================================================
void FixLbFluid::calc_mass_momentum(double &totalmass, double totalmomentum[3])
{
  double localmass = 0.0;
  double localmomentum[3] = {0., 0., 0.};

  for (int i = 1; i < subNbx - 1; i++)
    for (int j = 1; j < subNby - 1; j++)
      for (int k = 1; k < subNbz - 1; k++) {
        localmass += density_lb[i][j][k];
        localmomentum[0] += density_lb[i][j][k] * u_lb[i][j][k][0];
        localmomentum[1] += density_lb[i][j][k] * u_lb[i][j][k][1];
        localmomentum[2] += density_lb[i][j][k] * u_lb[i][j][k][2];
      }
  //Need to do an allreduce if there is more than one processor.
  MPI_Allreduce(&localmass, &totalmass, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&localmomentum[0], &totalmomentum[0], 3, MPI_DOUBLE, MPI_SUM, world);

  totalmomentum[0] *= dm_lb * dx_lb / dt_lb;
  totalmomentum[1] *= dm_lb * dx_lb / dt_lb;
  totalmomentum[2] *= dm_lb * dx_lb / dt_lb;
  totalmass *= dm_lb;
}

void FixLbFluid::calc_MPT(double &totalmass, double totalmomentum[3], double &Tave)
{
  double localmass = 0.0;
  double localmomentum[3] = {0., 0., 0.};
  double localTave = 0.0;

  for (int i = 1; i < subNbx - 1; i++)
    for (int j = 1; j < subNby - 1; j++)
      for (int k = 1; k < subNbz - 1; k++) {
        localmass += density_lb[i][j][k];
        localmomentum[0] += density_lb[i][j][k] * u_lb[i][j][k][0];
        localmomentum[1] += density_lb[i][j][k] * u_lb[i][j][k][1];
        localmomentum[2] += density_lb[i][j][k] * u_lb[i][j][k][2];

        localTave += (u_lb[i][j][k][0] * u_lb[i][j][k][0] + u_lb[i][j][k][1] * u_lb[i][j][k][1] +
                      u_lb[i][j][k][2] * u_lb[i][j][k][2]) *
            density_lb[i][j][k];
      }

  //Need to do an allreduce if there is more than one processor.
  MPI_Allreduce(&localmass, &totalmass, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&localmomentum[0], &totalmomentum[0], 3, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&localTave, &Tave, 1, MPI_DOUBLE, MPI_SUM, world);

  totalmomentum[0] *= dm_lb * dx_lb / dt_lb;
  totalmomentum[1] *= dm_lb * dx_lb / dt_lb;
  totalmomentum[2] *= dm_lb * dx_lb / dt_lb;
  totalmass *= dm_lb;

  double kB = (force->boltz / force->mvv2e) * dt_lb * dt_lb / dx_lb / dx_lb / dm_lb;
  Tave = Tave / (3.0 * Nbx * Nby * Nbz * kB);
}

/* ----------------------------------------------------------------------
   lb/fluid variables that can be accessed from outside:

------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

int FixLbFluid::adjust_dof_fix() /* Based on same private method in compute class */
{                                /* altered to return fix_dof */
  int fix_dof = 0;
  for (auto &ifix : modify->get_fix_list())
    if (ifix->dof_flag) fix_dof += ifix->dof(igroup);
  return fix_dof;
}

double FixLbFluid::dof_compute() /* Based on same protected member of compute_temp class */
{ /* with extra_dof variable replaced by its default domain->dimension */

  double dof;

  if (setdof)
    dof = setdof;
  else {
    MPI_Allreduce(&dof_lb, &dof, 1, MPI_DOUBLE, MPI_SUM, world);
    dof = 3.0 * dof;
  }

  double tfactor;
  if (dof > FLT_EPSILON)
    tfactor = force->mvv2e / (dof * force->boltz);
  else
    tfactor = 0.0;
  return tfactor;
}

/* ---------------------------------------------------------------------- */

double FixLbFluid::compute_scalar() /* Based on same member of compute_temp class */
{
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        t += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * (rmass[i] + massp[i]);
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        t += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) *
            (mass[type[i]] + massp[i]);
      }
  }

  double scalar;    // extra declaration needed as not in compute
  MPI_Allreduce(&t, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

  scalar *= dof_compute();
  return scalar;
}

double FixLbFluid::compute_vector(int n)
{
  static int laststep = -1;
  static double totalmass, totalmomentum[3], Tav;

  if (laststep != update->ntimestep) {
    laststep = update->ntimestep;

    double fluidmass, fluidmomentum[3];
    calc_MPT(fluidmass, fluidmomentum, Tav);

    double particlemass = 0, particlevcm[3] = {0, 0, 0};
    particlemass = group->mass(0);    // use igroup=0 here and for vcm to get group 'all' particles,
                                      // not just ones in this fix
    group->vcm(0, particlemass, particlevcm);

    totalmass = fluidmass + particlemass;

    totalmomentum[0] = particlemass * particlevcm[0] + fluidmomentum[0];
    totalmomentum[1] = particlemass * particlevcm[1] + fluidmomentum[1];
    totalmomentum[2] = particlemass * particlevcm[2] + fluidmomentum[2];
  }

  switch (n) {
    case 0:
      return Tav;
      break;
    case 1:
      return totalmass;
      break;
    case 2:
      return totalmomentum[0];
      break;
    case 3:
      return totalmomentum[1];
      break;
    case 4:
      return totalmomentum[2];
      break;
    default:
      error->all(FLERR, "No fix lb/fluid variable available for that index.");
      return -1;
  }
}
