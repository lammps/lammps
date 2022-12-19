// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Paul Crozier (SNL)
                         Carolyn Phillips (University of Michigan)
------------------------------------------------------------------------- */

#include "fix_ttm.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "random_mars.h"
#include "respa.h"
#include "potential_file_reader.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using namespace FixConst;

// OFFSET avoids outside-of-box atoms being rounded to grid pts incorrectly
// SHIFT = 0.0 assigns atoms to lower-left grid pt
// SHIFT = 0.5 assigns atoms to nearest grid pt
// use SHIFT = 0.0 for now since it allows fix ave/chunk
//   to spatially average consistent with the TTM grid

static constexpr int OFFSET = 16384;
static constexpr double SHIFT = 0.0;

/* ---------------------------------------------------------------------- */

FixTTM::FixTTM(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  random(nullptr),
  gfactor1(nullptr), gfactor2(nullptr), ratio(nullptr), flangevin(nullptr),
  T_electron(nullptr), T_electron_old(nullptr),
  net_energy_transfer(nullptr), net_energy_transfer_all(nullptr)
{
  if (narg < 13) error->all(FLERR,"Illegal fix ttm command");

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 1;
  nevery = 1;
  restart_peratom = 1;
  restart_global = 1;

  seed = utils::inumeric(FLERR,arg[3],false,lmp);
  electronic_specific_heat = utils::numeric(FLERR,arg[4],false,lmp);
  electronic_density = utils::numeric(FLERR,arg[5],false,lmp);
  electronic_thermal_conductivity = utils::numeric(FLERR,arg[6],false,lmp);
  gamma_p = utils::numeric(FLERR,arg[7],false,lmp);
  gamma_s = utils::numeric(FLERR,arg[8],false,lmp);
  v_0 = utils::numeric(FLERR,arg[9],false,lmp);
  nxgrid = utils::inumeric(FLERR,arg[10],false,lmp);
  nygrid = utils::inumeric(FLERR,arg[11],false,lmp);
  nzgrid = utils::inumeric(FLERR,arg[12],false,lmp);

  tinit = 0.0;
  infile = outfile = nullptr;

  int iarg = 13;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"set") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ttm command");
      tinit = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (tinit <= 0.0)
        error->all(FLERR,"Fix ttm initial temperature must be > 0.0");
      iarg += 2;
    } else if (strcmp(arg[iarg],"infile") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ttm command");
      infile = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"outfile") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix ttm command");
      outevery = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      outfile = utils::strdup(arg[iarg+2]);
      iarg += 3;
    } else error->all(FLERR,"Illegal fix ttm command");
  }

  // error check

  if (seed <= 0)
    error->all(FLERR,"Invalid random number seed in fix ttm command");
  if (electronic_specific_heat <= 0.0)
    error->all(FLERR,"Fix ttm electronic_specific_heat must be > 0.0");
  if (electronic_density <= 0.0)
    error->all(FLERR,"Fix ttm electronic_density must be > 0.0");
  if (electronic_thermal_conductivity < 0.0)
    error->all(FLERR,"Fix ttm electronic_thermal_conductivity must be >= 0.0");
  if (gamma_p <= 0.0) error->all(FLERR,"Fix ttm gamma_p must be > 0.0");
  if (gamma_s < 0.0) error->all(FLERR,"Fix ttm gamma_s must be >= 0.0");
  if (v_0 < 0.0) error->all(FLERR,"Fix ttm v_0 must be >= 0.0");
  if (nxgrid <= 0 || nygrid <= 0 || nzgrid <= 0)
    error->all(FLERR,"Fix ttm grid sizes must be > 0");

  v_0_sq = v_0*v_0;

  // grid OFFSET to perform
  // SHIFT to map atom to nearest or lower-left grid point

  shift = OFFSET + SHIFT;

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  // allocate per-type arrays for force prefactors

  gfactor1 = new double[atom->ntypes+1];
  gfactor2 = new double[atom->ntypes+1];

  // check for allowed maximum number of total grid points

  bigint totalgrid = (bigint) nxgrid * nygrid * nzgrid;
  if (totalgrid > MAXSMALLINT)
    error->all(FLERR,"Too many grid points in fix ttm");
  ngridtotal = totalgrid;

  // allocate per-atom flangevin and zero it

  flangevin = nullptr;
  FixTTM::grow_arrays(atom->nmax);

  for (int i = 0; i < atom->nmax; i++) {
    flangevin[i][0] = 0.0;
    flangevin[i][1] = 0.0;
    flangevin[i][2] = 0.0;
  }

  // set 2 callbacks

  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  // determines which class deallocate_grid() is called from

  deallocate_flag = 0;
}

/* ---------------------------------------------------------------------- */

FixTTM::~FixTTM()
{
  delete [] infile;
  delete [] outfile;

  delete random;

  delete [] gfactor1;
  delete [] gfactor2;

  memory->destroy(flangevin);

  if (!deallocate_flag) FixTTM::deallocate_grid();
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_constructor()
{
  // allocate global grid on each proc
  // needs to be done in post_contructor() beccause is virtual method

  allocate_grid();

  // initialize electron temperatures on grid

  int ix,iy,iz;
  for (iz = 0; iz < nzgrid; iz++)
    for (iy = 0; iy < nygrid; iy++)
      for (ix = 0; ix < nxgrid; ix++)
        T_electron[iz][iy][ix] = tinit;

  // zero net_energy_transfer_all
  // in case compute_vector accesses it on timestep 0

  outflag = 0;
  memset(&net_energy_transfer_all[0][0][0],0,ngridtotal*sizeof(double));

  // set initial electron temperatures from user input file

  if (infile) read_electron_temperatures(infile);
}

/* ---------------------------------------------------------------------- */

int FixTTM::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTTM::init()
{
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix ttm with 2d simulation");
  if (domain->nonperiodic != 0)
    error->all(FLERR,"Cannot use non-periodic boundares with fix ttm");
  if (domain->triclinic)
    error->all(FLERR,"Cannot use fix ttm with triclinic box");

  // set force prefactors

  for (int i = 1; i <= atom->ntypes; i++) {
    gfactor1[i] = - gamma_p / force->ftm2v;
    gfactor2[i] =
      sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
  }

  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixTTM::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet")) {
    post_force_setup(vflag);
  } else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(nlevels_respa-1);
    post_force_respa_setup(vflag,nlevels_respa-1,0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_setup(int /*vflag*/)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // apply langevin forces that have been stored from previous run

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force(int /*vflag*/)
{
  int ix,iy,iz;
  double gamma1,gamma2;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double dxinv = nxgrid/domain->xprd;
  double dyinv = nygrid/domain->yprd;
  double dzinv = nzgrid/domain->zprd;

  // apply damping and thermostat to all atoms in fix group

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + shift) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + shift) - OFFSET;
      if (ix < 0) ix += nxgrid;
      if (iy < 0) iy += nygrid;
      if (iz < 0) iz += nzgrid;
      if (ix >= nxgrid) ix -= nxgrid;
      if (iy >= nygrid) iy -= nygrid;
      if (iz >= nzgrid) iz -= nzgrid;

      if (T_electron[iz][iy][ix] < 0)
        error->one(FLERR,"Electronic temperature dropped below zero");

      double tsqrt = sqrt(T_electron[iz][iy][ix]);

      gamma1 = gfactor1[type[i]];
      double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      if (vsq > v_0_sq) gamma1 *= (gamma_p + gamma_s)/gamma_p;
      gamma2 = gfactor2[type[i]] * tsqrt;

      flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
      flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
      flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);

      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_respa_setup(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force_setup(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTM::end_of_step()
{
  int ix,iy,iz;

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double dxinv = nxgrid/domain->xprd;
  double dyinv = nygrid/domain->yprd;
  double dzinv = nzgrid/domain->zprd;

  for (iz = 0; iz < nzgrid; iz++)
    for (iy = 0; iy < nygrid; iy++)
      for (ix = 0; ix < nxgrid; ix++)
        net_energy_transfer[iz][iy][ix] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + shift) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + shift) - OFFSET;
      if (ix < 0) ix += nxgrid;
      if (iy < 0) iy += nygrid;
      if (iz < 0) iz += nzgrid;
      if (ix >= nxgrid) ix -= nxgrid;
      if (iy >= nygrid) iy -= nygrid;
      if (iz >= nzgrid) iz -= nzgrid;

      net_energy_transfer[iz][iy][ix] +=
        (flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] +
         flangevin[i][2]*v[i][2]);
    }

  outflag = 0;
  MPI_Allreduce(&net_energy_transfer[0][0][0],&net_energy_transfer_all[0][0][0],
                ngridtotal,MPI_DOUBLE,MPI_SUM,world);

  double dx = domain->xprd/nxgrid;
  double dy = domain->yprd/nygrid;
  double dz = domain->zprd/nzgrid;
  double del_vol = dx*dy*dz;

  // num_inner_timesteps = # of inner steps (thermal solves)
  // required this MD step to maintain a stable explicit solve

  int num_inner_timesteps = 1;
  double inner_dt = update->dt;

  double stability_criterion = 1.0 -
    2.0*inner_dt/(electronic_specific_heat*electronic_density) *
    (electronic_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));

  if (stability_criterion < 0.0) {
    inner_dt = 0.5*(electronic_specific_heat*electronic_density) /
      (electronic_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
    num_inner_timesteps = static_cast<int>(update->dt/inner_dt) + 1;
    inner_dt = update->dt/double(num_inner_timesteps);
    if (num_inner_timesteps > 1000000)
      error->warning(FLERR,"Too many inner timesteps in fix ttm");
  }

  // finite difference iterations to update T_electron

  for (int istep = 0; istep < num_inner_timesteps; istep++) {

    for (iz = 0; iz < nzgrid; iz++)
      for (iy = 0; iy < nygrid; iy++)
        for (ix = 0; ix < nxgrid; ix++)
          T_electron_old[iz][iy][ix] = T_electron[iz][iy][ix];

    // compute new electron T profile

    for (iz = 0; iz < nzgrid; iz++)
      for (iy = 0; iy < nygrid; iy++)
        for (ix = 0; ix < nxgrid; ix++) {
          int xright = ix + 1;
          int yright = iy + 1;
          int zright = iz + 1;
          if (xright == nxgrid) xright = 0;
          if (yright == nygrid) yright = 0;
          if (zright == nzgrid) zright = 0;
          int xleft = ix - 1;
          int yleft = iy - 1;
          int zleft = iz - 1;
          if (xleft == -1) xleft = nxgrid - 1;
          if (yleft == -1) yleft = nygrid - 1;
          if (zleft == -1) zleft = nzgrid - 1;

          T_electron[iz][iy][ix] =
            T_electron_old[iz][iy][ix] +
            inner_dt/(electronic_specific_heat*electronic_density) *
            (electronic_thermal_conductivity *

             ((T_electron_old[iz][iy][xright] + T_electron_old[iz][iy][xleft] -
               2.0*T_electron_old[iz][iy][ix])/dx/dx +
              (T_electron_old[iz][yright][ix] + T_electron_old[iz][yleft][ix] -
               2.0*T_electron_old[iz][iy][ix])/dy/dy +
              (T_electron_old[zright][iy][ix] + T_electron_old[zleft][iy][ix] -
               2.0*T_electron_old[iz][iy][ix])/dz/dz) -

             (net_energy_transfer_all[iz][iy][ix])/del_vol);
        }
  }

  // output of grid electron temperatures to file

  if (outfile && (update->ntimestep % outevery == 0))
    write_electron_temperatures(fmt::format("{}.{}",outfile,update->ntimestep));
}

/* ----------------------------------------------------------------------
   read in initial electron temperatures from a user-specified file
   only read by proc 0, grid values are Bcast to other procs
------------------------------------------------------------------------- */

void FixTTM::read_electron_temperatures(const std::string &filename)
{
  if (comm->me == 0) {

    int ***T_initial_set;
    memory->create(T_initial_set,nxgrid,nygrid,nzgrid,"ttm:T_initial_set");
    memset(&T_initial_set[0][0][0],0,ngridtotal*sizeof(int));

    // read initial electron temperature values from file
    bigint nread = 0;

    try {
      PotentialFileReader reader(lmp, filename, "electron temperature grid");

      while (nread < ngridtotal) {
        // reader will skip over comment-only lines
        auto values = reader.next_values(4);
        ++nread;

        int ix = values.next_int() - 1;
        int iy = values.next_int() - 1;
        int iz = values.next_int() - 1;
        double T_tmp  = values.next_double();

        // check correctness of input data

        if ((ix < 0) || (ix >= nxgrid) || (iy < 0) || (iy >= nygrid) || (iz < 0) || (iz >= nzgrid))
          throw TokenizerException("Fix ttm invalid grid index in fix ttm grid file","");

        if (T_tmp < 0.0)
          throw TokenizerException("Fix ttm electron temperatures must be > 0.0","");

        T_electron[iz][iy][ix] = T_tmp;
        T_initial_set[iz][iy][ix] = 1;
      }
    } catch (std::exception &e) {
      error->one(FLERR, e.what());
    }

    // check completeness of input data

    for (int iz = 0; iz < nzgrid; iz++)
      for (int iy = 0; iy < nygrid; iy++)
        for (int ix = 0; ix < nxgrid; ix++)
          if (T_initial_set[iz][iy][ix] == 0)
            error->all(FLERR,"Fix ttm infile did not set all temperatures");

    memory->destroy(T_initial_set);
  }

  MPI_Bcast(&T_electron[0][0][0],ngridtotal,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   write out current electron temperatures to user-specified file
   only written by proc 0
------------------------------------------------------------------------- */

void FixTTM::write_electron_temperatures(const std::string &filename)
{
  if (comm->me) return;

  FILE *fp = fopen(filename.c_str(),"w");
  if (!fp) error->one(FLERR,"Fix ttm could not open output file {}: {}",
                      filename,utils::getsyserror());
  fmt::print(fp,"# DATE: {} UNITS: {} COMMENT: Electron temperature on "
             "{}x{}x{} grid at step {} - created by fix {}\n", utils::current_date(),
             update->unit_style, nxgrid, nygrid, nzgrid, update->ntimestep, style);

  int ix,iy,iz;

  for (iz = 0; iz < nzgrid; iz++)
    for (iy = 0; iy < nygrid; iy++)
      for (ix = 0; ix < nxgrid; ix++)
        fprintf(fp,"%d %d %d %20.16g\n",ix+1,iy+1,iz+1,T_electron[iz][iy][ix]);

  fclose(fp);
}

/* ---------------------------------------------------------------------- */

void FixTTM::reset_dt()
{
  for (int i = 1; i <= atom->ntypes; i++)
    gfactor2[i] =
      sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
}

/* ---------------------------------------------------------------------- */

void FixTTM::grow_arrays(int ngrow)
{
  memory->grow(flangevin,ngrow,3,"ttm:flangevin");
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTTM::write_restart(FILE *fp)
{
  double *rlist;
  memory->create(rlist,nxgrid*nygrid*nzgrid+4,"ttm:rlist");

  int n = 0;
  rlist[n++] = nxgrid;
  rlist[n++] = nygrid;
  rlist[n++] = nzgrid;
  rlist[n++] = seed;

  // store global grid values

  for (int iz = 0; iz < nzgrid; iz++)
    for (int iy = 0; iy < nygrid; iy++)
      for (int ix = 0; ix < nxgrid; ix++)
        rlist[n++] =  T_electron[iz][iy][ix];

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(rlist,sizeof(double),n,fp);
  }

  memory->destroy(rlist);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTTM::restart(char *buf)
{
  int n = 0;
  auto rlist = (double *) buf;

  // check that restart grid size is same as current grid size

  int nxgrid_old = static_cast<int> (rlist[n++]);
  int nygrid_old = static_cast<int> (rlist[n++]);
  int nzgrid_old = static_cast<int> (rlist[n++]);

  if (nxgrid_old != nxgrid || nygrid_old != nygrid || nzgrid_old != nzgrid)
    error->all(FLERR,"Must restart fix ttm with same grid size");

  // change RN seed from initial seed, to avoid same Langevin factors
  // just increment by 1, since for RanMars that is a new RN stream

  seed = static_cast<int> (rlist[n++]) + 1;
  delete random;
  random = new RanMars(lmp,seed+comm->me);

  // restore global grid values

  for (int iz = 0; iz < nzgrid; iz++)
    for (int iy = 0; iy < nygrid; iy++)
      for (int ix = 0; ix < nxgrid; ix++)
        T_electron[iz][iy][ix] = rlist[n++];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixTTM::pack_restart(int i, double *buf)
{
  // pack buf[0] this way because other fixes unpack it

  buf[0] = 4;
  buf[1] = flangevin[i][0];
  buf[2] = flangevin[i][1];
  buf[3] = flangevin[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixTTM::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  flangevin[nlocal][0] = extra[nlocal][m++];
  flangevin[nlocal][1] = extra[nlocal][m++];
  flangevin[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixTTM::size_restart(int /*nlocal*/)
{
  return 4;
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixTTM::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   return the energy of the electronic subsystem or the net_energy transfer
   between the subsystems
------------------------------------------------------------------------- */

double FixTTM::compute_vector(int n)
{
  if (outflag == 0) {
    e_energy = 0.0;
    transfer_energy = 0.0;

    int ix,iy,iz;

    double dx = domain->xprd/nxgrid;
    double dy = domain->yprd/nygrid;
    double dz = domain->zprd/nzgrid;
    double del_vol = dx*dy*dz;

    for (iz = 0; iz < nzgrid; iz++)
      for (iy = 0; iy < nygrid; iy++)
        for (ix = 0; ix < nxgrid; ix++) {
          e_energy +=
            T_electron[iz][iy][ix]*electronic_specific_heat*
            electronic_density*del_vol;
          transfer_energy +=
            net_energy_transfer_all[iz][iy][ix]*update->dt;
          //printf("TRANSFER %d %d %d %g\n",ix,iy,iz,transfer_energy);
        }

    //printf("TRANSFER %g\n",transfer_energy);

    outflag = 1;
  }

  if (n == 0) return e_energy;
  if (n == 1) return transfer_energy;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage for flangevin and 3d grids
------------------------------------------------------------------------- */

double FixTTM::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) atom->nmax * 3 * sizeof(double);
  bytes += (double) 4*ngridtotal * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTM::allocate_grid()
{
  memory->create(T_electron_old,nzgrid,nygrid,nxgrid,"ttm:T_electron_old");
  memory->create(T_electron,nzgrid,nygrid,nxgrid,"ttm:T_electron");
  memory->create(net_energy_transfer,nzgrid,nygrid,nxgrid,
                 "ttm:net_energy_transfer");
  memory->create(net_energy_transfer_all,nzgrid,nygrid,nxgrid,
                 "ttm:net_energy_transfer_all");
}

/* ----------------------------------------------------------------------
   deallocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTM::deallocate_grid()
{
  memory->destroy(T_electron_old);
  memory->destroy(T_electron);
  memory->destroy(net_energy_transfer);
  memory->destroy(net_energy_transfer_all);
}
