// clang-format off
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
   Contributing authors: (in addition to authors of original fix ttm)
   Sergey Starikov (Joint Institute for High Temperatures of RAS)
   Vasily Pisarev (Joint Institute for High Temperatures of RAS)
------------------------------------------------------------------------- */

#include "fix_ttm_mod.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
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
using namespace MathConst;

// OFFSET avoids outside-of-box atoms being rounded to grid pts incorrectly
// SHIFT = 0.0 assigns atoms to lower-left grid pt
// SHIFT = 0.5 assigns atoms to nearest grid pt
// use SHIFT = 0.0 for now since it allows fix ave/chunk
//   to spatially average consistent with the TTM grid

static const char cite_fix_ttm_mod[] =
  "fix ttm/mod command: doi:10.1088/0953-8984/26/47/475401, doi:10.1002/ctpp.201310025\n\n"
  "@article{Pisarev2014,\n"
  "author = {Pisarev, V. V. and Starikov, S. V.},\n"
  "title = {Atomistic Simulation of Ion Track Formation in {UO$_2$}.},\n"
  "journal = {J.~Phys.\\ Condens.\\ Matter},\n"
  "volume = {26},\n"
  "number = {47},\n"
  "pages = {475401},\n"
  "year = {2014}\n"
  "}\n\n"
  "@article{Norman2013,\n"
  "author = {Norman, G. E. and Starikov, S. V. and Stegailov, V. V. and Saitov, I. M. and Zhilyaev, P. A.},\n"
  "title = {Atomistic Modeling of Warm Dense Matter in the Two-Temperature State},\n"
  "journal = {Contrib.\\ Plasma Phys.},\n"
  "number = {2},\n"
  "volume = {53},\n"
  "pages = {129--139},\n"
  "year = {2013}\n"
  "}\n\n";

static constexpr int OFFSET = 16384;
static constexpr double SHIFT = 0.0;

/* ---------------------------------------------------------------------- */

FixTTMMod::FixTTMMod(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  random(nullptr), nsum(nullptr), nsum_all(nullptr),
  gfactor1(nullptr), gfactor2(nullptr), ratio(nullptr), flangevin(nullptr),
  T_electron(nullptr), T_electron_old(nullptr), sum_vsq(nullptr), sum_mass_vsq(nullptr),
  sum_vsq_all(nullptr), sum_mass_vsq_all(nullptr), net_energy_transfer(nullptr),
  net_energy_transfer_all(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_ttm_mod);

  if (narg < 8) error->all(FLERR,"Illegal fix ttm/mod command");

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 1;
  nevery = 1;
  restart_peratom = 1;
  restart_global = 1;

  seed = utils::inumeric(FLERR,arg[3],false,lmp);

  nxgrid = utils::inumeric(FLERR,arg[5],false,lmp);
  nygrid = utils::inumeric(FLERR,arg[6],false,lmp);
  nzgrid = utils::inumeric(FLERR,arg[7],false,lmp);

  double tinit = 0.0;
  infile = outfile = nullptr;

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"set") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ttm/mod command");
      tinit = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (tinit <= 0.0)
        error->all(FLERR,"Fix ttm/mod initial temperature must be > 0.0");
      iarg += 2;
    } else if (strcmp(arg[iarg],"infile") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ttm/mod command");
      infile = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"outfile") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix ttm/mod command");
      outevery = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      outfile = utils::strdup(arg[iarg+2]);
      iarg += 3;
    } else error->all(FLERR,"Illegal fix ttm/mod command");
  }

  // error check

  if (seed <= 0)
    error->all(FLERR,"Invalid random number seed in fix ttm/mod command");
  if (nxgrid <= 0 || nygrid <= 0 || nzgrid <= 0)
    error->all(FLERR,"Fix ttm/mod grid sizes must be > 0");

  // check for allowed maximum number of total grid points

  bigint total_ngrid = (bigint) nxgrid * nygrid * nzgrid;
  if (total_ngrid > MAXSMALLINT)
    error->all(FLERR,"Too many grid points in fix ttm/mod");
  ngridtotal = total_ngrid;

  // t_surface is determined by electronic temperature (not constant)

  read_parameters(arg[4]);

  t_surface_l = surface_l;
  mult_factor = intensity;
  duration = 0.0;
  v_0_sq = v_0*v_0;
  surface_double = double(t_surface_l)*(domain->xprd/nxgrid);
  if ((C_limit+esheat_0) < 0.0)
    error->all(FLERR,"Fix ttm/mod electronic_specific_heat must be >= 0.0");
  if (electronic_density <= 0.0)
    error->all(FLERR,"Fix ttm/mod electronic_density must be > 0.0");
  if (gamma_p < 0.0) error->all(FLERR,"Fix ttm/mod gamma_p must be >= 0.0");
  if (gamma_s < 0.0) error->all(FLERR,"Fix ttm/mod gamma_s must be >= 0.0");
  if (v_0 < 0.0) error->all(FLERR,"Fix ttm/mod v_0 must be >= 0.0");
  if (ionic_density <= 0.0) error->all(FLERR,"Fix ttm/mod ionic_density must be > 0.0");
  if (surface_l < 0) error->all(FLERR,"Surface coordinates must be >= 0");
  if (surface_l >= surface_r) error->all(FLERR, "Left surface coordinate must be less than right surface coordinate");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  // allocate per-type arrays for force prefactors

  gfactor1 = new double[atom->ntypes+1];
  gfactor2 = new double[atom->ntypes+1];

  // allocate 3d grid variables

  memory->create(nsum,nxgrid,nygrid,nzgrid,"ttm/mod:nsum");
  memory->create(nsum_all,nxgrid,nygrid,nzgrid,"ttm/mod:nsum_all");
  memory->create(sum_vsq,nxgrid,nygrid,nzgrid,"ttm/mod:sum_vsq");
  memory->create(sum_mass_vsq,nxgrid,nygrid,nzgrid,"ttm/mod:sum_mass_vsq");
  memory->create(sum_vsq_all,nxgrid,nygrid,nzgrid,"ttm/mod:sum_vsq_all");
  memory->create(sum_mass_vsq_all,nxgrid,nygrid,nzgrid,
                 "ttm/mod:sum_mass_vsq_all");
  memory->create(T_electron_old,nxgrid,nygrid,nzgrid,"ttm/mod:T_electron_old");
  memory->create(T_electron_first,nxgrid,nygrid,nzgrid,"ttm/mod:T_electron_first");
  memory->create(T_electron,nxgrid,nygrid,nzgrid,"ttm/mod:T_electron");
  memory->create(net_energy_transfer,nxgrid,nygrid,nzgrid,
                 "ttm/mod:net_energy_transfer");
  memory->create(net_energy_transfer_all,nxgrid,nygrid,nzgrid,
                 "ttm/mod:net_energy_transfer_all");
  flangevin = nullptr;
  grow_arrays(atom->nmax);

  // grid OFFSET to perform
  // SHIFT to map atom to nearest or lower-left grid point

  shift = OFFSET + SHIFT;

  // zero out the flangevin array

  for (int i = 0; i < atom->nmax; i++) {
    flangevin[i][0] = 0.0;
    flangevin[i][1] = 0.0;
    flangevin[i][2] = 0.0;
  }

  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  // initialize electron temperatures on grid

  int ix,iy,iz;
  for (ix = 0; ix < nxgrid; ix++)
    for (iy = 0; iy < nygrid; iy++)
      for (iz = 0; iz < nzgrid; iz++)
        T_electron[ix][iy][iz] = tinit;

  // if specified, read initial electron temperatures from file

  if (infile) read_electron_temperatures(infile);
}

/* ---------------------------------------------------------------------- */

FixTTMMod::~FixTTMMod()
{
  delete random;
  delete[] gfactor1;
  delete[] gfactor2;

  memory->destroy(nsum);
  memory->destroy(nsum_all);
  memory->destroy(sum_vsq);
  memory->destroy(sum_mass_vsq);
  memory->destroy(sum_vsq_all);
  memory->destroy(sum_mass_vsq_all);
  memory->destroy(T_electron_first);
  memory->destroy(T_electron_old);
  memory->destroy(T_electron);
  memory->destroy(flangevin);
  memory->destroy(net_energy_transfer);
  memory->destroy(net_energy_transfer_all);
}

/* ---------------------------------------------------------------------- */

int FixTTMMod::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::init()
{
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix ttm/mod with 2d simulation");
  if (domain->nonperiodic != 0)
    error->all(FLERR,"Cannot use non-periodic boundares with fix ttm/mod");
  if (domain->triclinic)
    error->all(FLERR,"Cannot use fix ttm/mod with triclinic box");

  // set force prefactors

  for (int i = 1; i <= atom->ntypes; i++) {
    gfactor1[i] = - gamma_p / force->ftm2v;
    gfactor2[i] =
      sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
  }

  for (int ix = 0; ix < nxgrid; ix++)
    for (int iy = 0; iy < nygrid; iy++)
      for (int iz = 0; iz < nzgrid; iz++)
        net_energy_transfer_all[ix][iy][iz] = 0;

  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::setup(int vflag)
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

void FixTTMMod::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double dx = domain->xprd/nxgrid;
  double dy = domain->yprd/nygrid;
  double dz = domain->zprd/nzgrid;
  double gamma1,gamma2;

  // apply damping and thermostat to all atoms in fix group

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double xscale = (x[i][0] - domain->boxlo[0])/domain->xprd;
      double yscale = (x[i][1] - domain->boxlo[1])/domain->yprd;
      double zscale = (x[i][2] - domain->boxlo[2])/domain->zprd;
      int ix = static_cast<int>(xscale*nxgrid + shift) - OFFSET;
      int iy = static_cast<int>(yscale*nygrid + shift) - OFFSET;
      int iz = static_cast<int>(zscale*nzgrid + shift) - OFFSET;
      while (ix > nxgrid-1) ix -= nxgrid;
      while (iy > nygrid-1) iy -= nygrid;
      while (iz > nzgrid-1) iz -= nzgrid;
      while (ix < 0) ix += nxgrid;
      while (iy < 0) iy += nygrid;
      while (iz < 0) iz += nzgrid;

      if (T_electron[ix][iy][iz] < 0)
        error->all(FLERR,"Electronic temperature dropped below zero");

      double tsqrt = sqrt(T_electron[ix][iy][iz]);

      gamma1 = gfactor1[type[i]];
      double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      if (vsq > v_0_sq) gamma1 *= (gamma_p + gamma_s)/gamma_p;
      gamma2 = gfactor2[type[i]] * tsqrt;
      if (ix >= surface_l) {
        if (ix < surface_r) {
          flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
          flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
          flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
          double x_surf = dx*double(surface_l)+dx;
          double x_at = x[i][0] - domain->boxlo[0];
          int right_x = ix + 1;
          int right_y = iy + 1;
          int right_z = iz + 1;
          if (right_x == nxgrid) right_x = 0;
          if (right_y == nygrid) right_y = 0;
          if (right_z == nzgrid) right_z = 0;
          int left_x = ix - 1;
          int left_y = iy - 1;
          int left_z = iz - 1;
          if (left_x == -1) left_x = nxgrid - 1;
          if (left_y == -1) left_y = nygrid - 1;
          if (left_z == -1) left_z = nzgrid - 1;
          double T_i = T_electron[ix][iy][iz];
          double T_ir = T_electron[right_x][iy][iz];
          double T_iu = T_electron[ix][right_y][iz];
          double T_if = T_electron[ix][iy][right_z];
          double C_i = el_properties(T_electron[ix][iy][iz]).el_heat_capacity;
          double C_ir = el_properties(T_electron[right_x][iy][iz]).el_heat_capacity;
          double C_iu = el_properties(T_electron[ix][right_y][iz]).el_heat_capacity;
          double C_if = el_properties(T_electron[ix][iy][right_z]).el_heat_capacity;
          double diff_x = (x_at - x_surf)*(x_at - x_surf);
          diff_x = pow(diff_x,0.5);
          double len_factor = diff_x/(diff_x+free_path);
          if (movsur == 1) {
            if (x_at >= x_surf) {
              flangevin[i][0] -= pres_factor/ionic_density*((C_ir*T_ir*free_path/(diff_x+free_path)/(diff_x+free_path)) +
                                                            (len_factor/dx)*(C_ir*T_ir-C_i*T_i));
              flangevin[i][1] -= pres_factor/ionic_density/dy*(C_iu*T_iu-C_i*T_i);
              flangevin[i][2] -= pres_factor/ionic_density/dz*(C_if*T_if-C_i*T_i);
            }
          } else {
            flangevin[i][0] -= pres_factor/ionic_density/dx*(C_ir*T_ir-C_i*T_i);
            flangevin[i][1] -= pres_factor/ionic_density/dy*(C_iu*T_iu-C_i*T_i);
            flangevin[i][2] -= pres_factor/ionic_density/dz*(C_if*T_if-C_i*T_i);
          }
          f[i][0] += flangevin[i][0];
          f[i][1] += flangevin[i][1];
          f[i][2] += flangevin[i][2];
        }
      }
      if (movsur == 1) {
        if (ix < surface_l) {
          t_surface_l = ix;
        }
      }
    }
  }
  MPI_Allreduce(&t_surface_l,&surface_l,1,MPI_INT,MPI_MIN,world);
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::post_force_setup(int /*vflag*/)
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

void FixTTMMod::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::post_force_respa_setup(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force_setup(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::reset_dt()
{
  for (int i = 1; i <= atom->ntypes; i++)
    gfactor2[i] =
      sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
}

/* ----------------------------------------------------------------------
   read in ttm/mod parameters from a user-specified file
   only called by proc 0
------------------------------------------------------------------------- */

void FixTTMMod::read_parameters(const std::string &filename)
{
  if (comm->me == 0) {

    try {
      PotentialFileReader reader(lmp, filename, "ttm/mod parameter");

      // C0 (metal)

      reader.next_line();
      esheat_0 = reader.next_values(1).next_double();

      // C1 (metal*10^3)

      reader.next_line();
      esheat_1 = reader.next_values(1).next_double();

      // C2 (metal*10^6)

      reader.next_line();
      esheat_2 = reader.next_values(1).next_double();

      // C3 (metal*10^9)

      reader.next_line();
      esheat_3 = reader.next_values(1).next_double();

      // C4 (metal*10^12)

      reader.next_line();
      esheat_4 = reader.next_values(1).next_double();

      // C_limit

      reader.next_line();
      C_limit = reader.next_values(1).next_double();

      // Temperature damping factor

      reader.next_line();
      T_damp = reader.next_values(1).next_double();

      // rho_e

      reader.next_line();
      electronic_density = reader.next_values(1).next_double();

      // thermal_diffusion

      reader.next_line();
      el_th_diff = reader.next_values(1).next_double();

      // gamma_p

      reader.next_line();
      gamma_p = reader.next_values(1).next_double();

      // gamma_s

      reader.next_line();
      gamma_s = reader.next_values(1).next_double();

      // v0

      reader.next_line();
      v_0 = reader.next_values(1).next_double();

      // average intensity of pulse (source of energy) (metal units)

      reader.next_line();
      intensity = reader.next_values(1).next_double();

      // coordinate of 1st surface in x-direction (in box units) - constant

      reader.next_line();
      surface_l = reader.next_values(1).next_int();

      // coordinate of 2nd surface in x-direction (in box units) - constant

      reader.next_line();
      surface_r = reader.next_values(1).next_int();

      // skin_layer = intensity is reduced (I=I0*exp[-x/skin_layer])

      reader.next_line();
      skin_layer =  reader.next_values(1).next_int();

      // width of pulse (picoseconds)

      reader.next_line();
      width = reader.next_values(1).next_double();

      // factor of electronic pressure (PF) Pe = PF*Ce*Te

      reader.next_line();
      pres_factor = reader.next_values(1).next_double();

      // effective free path of electrons (angstrom)

      reader.next_line();
      free_path = reader.next_values(1).next_double();

      // ionic density (ions*angstrom^{-3})

      reader.next_line();
      ionic_density = reader.next_values(1).next_double();

      // if movsur = 0: surface is frozen

      reader.next_line();
      movsur = reader.next_values(1).next_int();

      // electron_temperature_min

      reader.next_line();
      electron_temperature_min = reader.next_values(1).next_double();
    } catch (std::exception &e) {
      error->one(FLERR,e.what());
    }
  }
  MPI_Bcast(&esheat_0, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&esheat_1, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&esheat_2, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&esheat_3, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&esheat_4, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&C_limit, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&T_damp, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&electronic_density, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&el_th_diff, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gamma_p, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gamma_s, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&v_0, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&intensity, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&surface_l, 1, MPI_INT, 0, world);
  MPI_Bcast(&surface_r, 1, MPI_INT, 0, world);
  MPI_Bcast(&skin_layer, 1, MPI_INT, 0, world);
  MPI_Bcast(&width, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&pres_factor, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&free_path, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&ionic_density, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&movsur, 1, MPI_INT, 0, world);
  MPI_Bcast(&electron_temperature_min, 1, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   read in initial electron temperatures from a user-specified file
   only read by proc 0, grid values are Bcast to other procs
------------------------------------------------------------------------- */

void FixTTMMod::read_electron_temperatures(const std::string &filename)
{
  if (comm->me == 0) {

    int ***T_initial_set;
    memory->create(T_initial_set,nxgrid,nygrid,nzgrid,"ttm/mod:T_initial_set");
    memset(&T_initial_set[0][0][0],0,ngridtotal*sizeof(int));

    // read initial electron temperature values from file
    bigint nread = 0;

    try {
      PotentialFileReader reader(lmp, filename, "electron temperature grid");

      while (nread < ngridtotal) {
        // reader will skip over comment-only lines
        auto values = reader.next_values(4);
        ++nread;

        int ix = values.next_int();
        int iy = values.next_int();
        int iz = values.next_int();
        double T_tmp  = values.next_double();

        // check correctness of input data

        if ((ix < 0) || (ix >= nxgrid) || (iy < 0) || (iy >= nygrid) || (iz < 0) || (iz >= nzgrid))
          throw TokenizerException("Fix ttm invalid grid index in fix ttm/mod grid file","");

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
            error->all(FLERR,"Fix ttm/mod infile did not set all temperatures");

    memory->destroy(T_initial_set);
  }

  MPI_Bcast(&T_electron[0][0][0],ngridtotal,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   write out current electron temperatures to user-specified file
   only written by proc 0
------------------------------------------------------------------------- */

void FixTTMMod::write_electron_temperatures(const std::string &filename)
{
  if (comm->me) return;

  FILE *fp = fopen(filename.c_str(),"w");
  if (!fp) error->one(FLERR,"Fix ttm/mod could not open output file {}: {}",
                      filename, utils::getsyserror());
  fmt::print(fp,"# DATE: {} UNITS: {} COMMENT: Electron temperature "
             "{}x{}x{} grid at step {}. Created by fix {}\n", utils::current_date(),
             update->unit_style, nxgrid, nygrid, nzgrid, update->ntimestep, style);

  int ix,iy,iz;

  for (ix = 0; ix < nxgrid; ix++)
    for (iy = 0; iy < nygrid; iy++)
      for (iz = 0; iz < nzgrid; iz++) {
        if (movsur == 1 && T_electron[ix][iy][iz] == 0.0)
          T_electron[ix][iy][iz] = electron_temperature_min;
        fprintf(fp,"%d %d %d %20.16g\n",ix,iy,iz,T_electron[ix][iy][iz]);
      }

  fclose(fp);
}

/* ---------------------------------------------------------------------- */

el_heat_capacity_thermal_conductivity FixTTMMod::el_properties(double T_e)
{
  el_heat_capacity_thermal_conductivity properties;
  double T_temp = T_e/1000.0, T_reduced = T_damp*T_temp;
  double T2 = T_temp*T_temp;
  double T3 = T2*T_temp;
  double T4 = T3*T_temp;
  double poly = esheat_0 + esheat_1*T_temp + esheat_2*T2 + esheat_3*T3 + esheat_4*T4;
  properties.el_heat_capacity = electronic_density*(poly*exp(-T_reduced*T_reduced) + C_limit); // heat capacity
  properties.el_thermal_conductivity = el_th_diff*properties.el_heat_capacity; // thermal conductivity
  return properties;
}
double FixTTMMod::el_sp_heat_integral(double T_e)
{
  double T_temp = T_e/1000.0, T_reduced = T_damp*T_temp;
  if (T_damp != 0)
    return electronic_density*(MY_PIS*(3*esheat_4/pow(T_damp,5)+2*esheat_2/pow(T_damp,3)+4*esheat_0/T_damp)*erf(T_reduced)+
                               4*esheat_3/pow(T_damp,4)+4*esheat_1/T_damp/T_damp-
                               ((6*esheat_4*T_temp+4*esheat_3)/pow(T_damp,4)+
                                (4*esheat_1+4*esheat_4*pow(T_temp,3)+4*esheat_3*T_temp*T_temp+4*esheat_2*T_temp)/T_damp/T_damp)*exp(-T_reduced*T_reduced))*125.0+electronic_density*C_limit*T_e;
  else
    return electronic_density*((esheat_0 + C_limit)*T_e + esheat_1*T_temp*T_e/2.0 + esheat_2*T_temp*T_temp*T_e/3.0 + esheat_3*pow(T_temp,3)*T_e/4.0 + esheat_4*pow(T_temp,4)*T_e/5.0);
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::end_of_step()
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (movsur == 1) {
    for (int ix = 0; ix < nxgrid; ix++)
      for (int iy = 0; iy < nygrid; iy++)
        for (int iz = 0; iz < nzgrid; iz++) {
          double TTT = T_electron[ix][iy][iz];
          if (TTT > 0) {
            if (ix < t_surface_l)
              t_surface_l = ix;
          }
        }
  }
  for (int ix = 0; ix < nxgrid; ix++)
    for (int iy = 0; iy < nygrid; iy++)
      for (int iz = 0; iz < nzgrid; iz++)
        net_energy_transfer[ix][iy][iz] = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      double xscale = (x[i][0] - domain->boxlo[0])/domain->xprd;
      double yscale = (x[i][1] - domain->boxlo[1])/domain->yprd;
      double zscale = (x[i][2] - domain->boxlo[2])/domain->zprd;
      int ix = static_cast<int>(xscale*nxgrid + shift) - OFFSET;
      int iy = static_cast<int>(yscale*nygrid + shift) - OFFSET;
      int iz = static_cast<int>(zscale*nzgrid + shift) - OFFSET;
      while (ix > nxgrid-1) ix -= nxgrid;
      while (iy > nygrid-1) iy -= nygrid;
      while (iz > nzgrid-1) iz -= nzgrid;
      while (ix < 0) ix += nxgrid;
      while (iy < 0) iy += nygrid;
      while (iz < 0) iz += nzgrid;
      if (ix >= t_surface_l) {
        if (ix < surface_r)
          net_energy_transfer[ix][iy][iz] +=
            (flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] +
             flangevin[i][2]*v[i][2]);
      }
    }

  MPI_Allreduce(&net_energy_transfer[0][0][0],
                &net_energy_transfer_all[0][0][0],
                ngridtotal,MPI_DOUBLE,MPI_SUM,world);

  double dx = domain->xprd/nxgrid;
  double dy = domain->yprd/nygrid;
  double dz = domain->zprd/nzgrid;
  double del_vol = dx*dy*dz;
  double el_specific_heat = 0.0;
  double el_thermal_conductivity = el_properties(electron_temperature_min).el_thermal_conductivity;
  for (int ix = 0; ix < nxgrid; ix++)
    for (int iy = 0; iy < nygrid; iy++)
      for (int iz = 0; iz < nzgrid; iz++)
        {
          if (el_properties(T_electron[ix][iy][iz]).el_thermal_conductivity > el_thermal_conductivity)
            el_thermal_conductivity = el_properties(T_electron[ix][iy][iz]).el_thermal_conductivity;
          if (el_specific_heat > 0.0)
            {
              if ((T_electron[ix][iy][iz] > 0.0) && (el_properties(T_electron[ix][iy][iz]).el_heat_capacity < el_specific_heat))
                el_specific_heat = el_properties(T_electron[ix][iy][iz]).el_heat_capacity;
            }
          else if (T_electron[ix][iy][iz] > 0.0) el_specific_heat = el_properties(T_electron[ix][iy][iz]).el_heat_capacity;
        }
  // num_inner_timesteps = # of inner steps (thermal solves)
  // required this MD step to maintain a stable explicit solve

  int num_inner_timesteps = 1;
  double inner_dt = update->dt;
  double stability_criterion = 0.0;

  for (int ix = 0; ix < nxgrid; ix++)
    for (int iy = 0; iy < nygrid; iy++)
      for (int iz = 0; iz < nzgrid; iz++)
        T_electron_first[ix][iy][iz] =
          T_electron[ix][iy][iz];
  do {
    for (int ix = 0; ix < nxgrid; ix++)
      for (int iy = 0; iy < nygrid; iy++)
        for (int iz = 0; iz < nzgrid; iz++)
          T_electron[ix][iy][iz] =
            T_electron_first[ix][iy][iz];

    stability_criterion = 1.0 -
      2.0*inner_dt/el_specific_heat *
      (el_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
    if (stability_criterion < 0.0) {
      inner_dt = 0.25*el_specific_heat /
        (el_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
    }
    num_inner_timesteps = static_cast<unsigned int>(update->dt/inner_dt) + 1;
    inner_dt = update->dt/double(num_inner_timesteps);
    if (num_inner_timesteps > 1000000)
      error->warning(FLERR,"Too many inner timesteps in fix ttm/mod");
    for (int ith_inner_timestep = 0; ith_inner_timestep < num_inner_timesteps;
         ith_inner_timestep++) {
      for (int ix = 0; ix < nxgrid; ix++)
        for (int iy = 0; iy < nygrid; iy++)
          for (int iz = 0; iz < nzgrid; iz++)
            T_electron_old[ix][iy][iz] =
              T_electron[ix][iy][iz];
      // compute new electron T profile
      duration = duration + inner_dt;
      for (int ix = 0; ix < nxgrid; ix++)
        for (int iy = 0; iy < nygrid; iy++)
          for (int iz = 0; iz < nzgrid; iz++) {
            int right_x = ix + 1;
            int right_y = iy + 1;
            int right_z = iz + 1;
            if (right_x == nxgrid) right_x = 0;
            if (right_y == nygrid) right_y = 0;
            if (right_z == nzgrid) right_z = 0;
            int left_x = ix - 1;
            int left_y = iy - 1;
            int left_z = iz - 1;
            if (left_x == -1) left_x = nxgrid - 1;
            if (left_y == -1) left_y = nygrid - 1;
            if (left_z == -1) left_z = nzgrid - 1;
            auto  skin_layer_d = double(skin_layer);
            auto  ix_d = double(ix);
            auto  surface_d = double(t_surface_l);
            mult_factor = 0.0;
            if (duration < width) {
              if (ix >= t_surface_l) mult_factor = (intensity/(dx*skin_layer_d))*exp((-1.0)*(ix_d - surface_d)/skin_layer_d);
            }
            if (ix < t_surface_l) net_energy_transfer_all[ix][iy][iz] = 0.0;
            double cr_vac = 1;
            if (T_electron_old[ix][iy][iz] == 0) cr_vac = 0;
            double cr_v_l_x = 1;
            if (T_electron_old[left_x][iy][iz] == 0) cr_v_l_x = 0;
            double cr_v_r_x = 1;
            if (T_electron_old[right_x][iy][iz] == 0) cr_v_r_x = 0;
            double cr_v_l_y = 1;
            if (T_electron_old[ix][left_y][iz] == 0) cr_v_l_y = 0;
            double cr_v_r_y = 1;
            if (T_electron_old[ix][right_y][iz] == 0) cr_v_r_y = 0;
            double cr_v_l_z = 1;
            if (T_electron_old[ix][iy][left_z] == 0) cr_v_l_z = 0;
            double cr_v_r_z = 1;
            if (T_electron_old[ix][iy][right_z] == 0) cr_v_r_z = 0;
            if (cr_vac != 0) {
              T_electron[ix][iy][iz] =
                T_electron_old[ix][iy][iz] +
                inner_dt/el_properties(T_electron_old[ix][iy][iz]).el_heat_capacity *
                ((cr_v_r_x*el_properties(T_electron_old[ix][iy][iz]/2.0+T_electron_old[right_x][iy][iz]/2.0).el_thermal_conductivity*
                  (T_electron_old[right_x][iy][iz]-T_electron_old[ix][iy][iz])/dx -
                  cr_v_l_x*el_properties(T_electron_old[ix][iy][iz]/2.0+T_electron_old[left_x][iy][iz]/2.0).el_thermal_conductivity*
                  (T_electron_old[ix][iy][iz]-T_electron_old[left_x][iy][iz])/dx)/dx +
                 (cr_v_r_y*el_properties(T_electron_old[ix][iy][iz]/2.0+T_electron_old[ix][right_y][iz]/2.0).el_thermal_conductivity*
                  (T_electron_old[ix][right_y][iz]-T_electron_old[ix][iy][iz])/dy -
                  cr_v_l_y*el_properties(T_electron_old[ix][iy][iz]/2.0+T_electron_old[ix][left_y][iz]/2.0).el_thermal_conductivity*
                  (T_electron_old[ix][iy][iz]-T_electron_old[ix][left_y][iz])/dy)/dy +
                 (cr_v_r_z*el_properties(T_electron_old[ix][iy][iz]/2.0+T_electron_old[ix][iy][right_z]/2.0).el_thermal_conductivity*
                  (T_electron_old[ix][iy][right_z]-T_electron_old[ix][iy][iz])/dz -
                  cr_v_l_z*el_properties(T_electron_old[ix][iy][iz]/2.0+T_electron_old[ix][iy][left_z]/2.0).el_thermal_conductivity*
                  (T_electron_old[ix][iy][iz]-T_electron_old[ix][iy][left_z])/dz)/dz);
              T_electron[ix][iy][iz]+=inner_dt/el_properties(T_electron[ix][iy][iz]).el_heat_capacity*
                (mult_factor -
                 net_energy_transfer_all[ix][iy][iz]/del_vol);
            }
            else T_electron[ix][iy][iz] =
                   T_electron_old[ix][iy][iz];
            if ((T_electron[ix][iy][iz] > 0.0) && (T_electron[ix][iy][iz] < electron_temperature_min))
              T_electron[ix][iy][iz] = T_electron[ix][iy][iz] + 0.5*(electron_temperature_min - T_electron[ix][iy][iz]);

            if (el_properties(T_electron[ix][iy][iz]).el_thermal_conductivity > el_thermal_conductivity)
              el_thermal_conductivity = el_properties(T_electron[ix][iy][iz]).el_thermal_conductivity;
            if ((T_electron[ix][iy][iz] > 0.0) && (el_properties(T_electron[ix][iy][iz]).el_heat_capacity < el_specific_heat))
              el_specific_heat = el_properties(T_electron[ix][iy][iz]).el_heat_capacity;
          }
    }
    stability_criterion = 1.0 -
      2.0*inner_dt/el_specific_heat *
      (el_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));

  } while (stability_criterion < 0.0);

  // output of grid electron temperatures to file

  if (outfile && (update->ntimestep % outevery == 0))
    write_electron_temperatures(fmt::format("{}.{}", outfile, update->ntimestep));
}

/* ----------------------------------------------------------------------
   memory usage of 3d grid
------------------------------------------------------------------------- */

double FixTTMMod::memory_usage()
{
  double bytes = 0.0;
  bytes += (double)5*ngridtotal * sizeof(int);
  bytes += (double)14*ngridtotal * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::grow_arrays(int ngrow)
{
  memory->grow(flangevin,ngrow,3,"ttm/mod:flangevin");
}

/* ----------------------------------------------------------------------
   return the energy of the electronic subsystem or the net_energy transfer
   between the subsystems
------------------------------------------------------------------------- */

double FixTTMMod::compute_vector(int n)
{
  double e_energy = 0.0;
  double transfer_energy = 0.0;

  double dx = domain->xprd/nxgrid;
  double dy = domain->yprd/nygrid;
  double dz = domain->zprd/nzgrid;
  double del_vol = dx*dy*dz;

  for (int ix = 0; ix < nxgrid; ix++)
    for (int iy = 0; iy < nygrid; iy++)
      for (int iz = 0; iz < nzgrid; iz++) {
        e_energy += el_sp_heat_integral(T_electron[ix][iy][iz])*del_vol;
        transfer_energy +=
          net_energy_transfer_all[ix][iy][iz]*update->dt;
      }

  if (n == 0) return e_energy;
  if (n == 1) return transfer_energy;
  return 0.0;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTTMMod::write_restart(FILE *fp)
{
  double *rlist;
  memory->create(rlist,nxgrid*nygrid*nzgrid+4,"ttm/mod:rlist");

  int n = 0;
  rlist[n++] = nxgrid;
  rlist[n++] = nygrid;
  rlist[n++] = nzgrid;
  rlist[n++] = seed;

  for (int ix = 0; ix < nxgrid; ix++)
    for (int iy = 0; iy < nygrid; iy++)
      for (int iz = 0; iz < nzgrid; iz++)
        rlist[n++] =  T_electron[ix][iy][iz];

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

void FixTTMMod::restart(char *buf)
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

  // restore global frid values

  for (int ix = 0; ix < nxgrid; ix++)
    for (int iy = 0; iy < nygrid; iy++)
      for (int iz = 0; iz < nzgrid; iz++)
        T_electron[ix][iy][iz] = rlist[n++];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixTTMMod::pack_restart(int i, double *buf)
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

void FixTTMMod::unpack_restart(int nlocal, int nth)
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
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixTTMMod::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixTTMMod::size_restart(int /*nlocal*/)
{
  return 4;
}
