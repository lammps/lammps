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
   Contributing authors: Paul Crozier (SNL)
                         Carolyn Phillips (University of Michigan)
------------------------------------------------------------------------- */

#include "fix_ttm_grid.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "gridcomm.h"
#include "neighbor.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

#include "tokenizer.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAXLINE 1024
#define OFFSET 16384

enum{T_ELECTRON,NET_ENERGY_TRANSFER,DSUM,SUM_VSQ,SUM_MASS_VSQ};

/* ---------------------------------------------------------------------- */

FixTTMGrid::FixTTMGrid(LAMMPS *lmp, int narg, char **arg) : 
  FixTTM(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixTTMGrid::post_force(int /*vflag*/)
{
  int ix,iy,iz;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double gamma1,gamma2;

  double *boxlo = domain->boxlo;

  // apply damping and thermostat to all atoms in fix group

  int flag = 0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      ix = static_cast<int> ((x[i][0]-boxlo[0])*delxinv + shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*delyinv + shift) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*delzinv + shift) - OFFSET;

      if (T_electron[ix][iy][iz] < 0)
        error->all(FLERR,"Electronic temperature dropped below zero");

      // check that ix,iy,iz is within my ghost cell range

      if (ix < nxlo_out || ix > nxhi_out || iy < nylo_out || iy > nyhi_out ||
          iz < nzlo_out || iz > nzhi_out) flag = 1;

      double tsqrt = sqrt(T_electron[ix][iy][iz]);

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

  if (flag) error->one(FLERR,"Out of range fix ttm/grid atoms");
}

/* ---------------------------------------------------------------------- */

void FixTTMGrid::end_of_step()
{
  int ix,iy,iz;

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;

  memset(&net_energy_transfer[nzlo_out][nylo_out][nxlo_out],0,
         ngridout*sizeof(double));

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ix = static_cast<int> ((x[i][0]-boxlo[0])*delxinv+shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*delyinv+shift) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*delzinv+shift) - OFFSET;
      net_energy_transfer[ix][iy][iz] +=
        (flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] +
         flangevin[i][2]*v[i][2]);
    }

  gc->reverse_comm(1,this,1,sizeof(double),NET_ENERGY_TRANSFER,
                   gc_buf1,gc_buf2,MPI_DOUBLE);

  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nznodes;
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
      error->warning(FLERR,"Too many inner timesteps in fix ttm/grid");
  }

  for (int ith_inner_timestep = 0; ith_inner_timestep < num_inner_timesteps;
       ith_inner_timestep++) {

    memcpy(&T_electron_old[nzlo_out][nylo_out][nxlo_out],
           &T_electron[nzlo_out][nylo_out][nxlo_out],ngridout*sizeof(double));

    // compute new electron T profile
 
    for (iz = nzlo_in; iz <= nzhi_in; iz++)
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++) {
          T_electron[ix][iy][iz] =
            T_electron_old[ix][iy][iz] +
            inner_dt/(electronic_specific_heat*electronic_density) *
            (electronic_thermal_conductivity *

             ((T_electron_old[ix-1][iy][iz] + T_electron_old[ix+1][iy][iz] -
               2*T_electron_old[ix][iy][iz])/dx/dx +
              (T_electron_old[ix][iy-1][iz] + T_electron_old[ix][iy+1][iz] -
               2*T_electron_old[ix][iy][iz])/dy/dy +
              (T_electron_old[ix][iy][iz-1] + T_electron_old[ix][iy][iz+1] -
               2*T_electron_old[ix][iy][iz])/dz/dz) -
             
             net_energy_transfer_all[ix][iy][iz]/del_vol);
        }
  }

  // comm new T_electron to ghost grid points

  gc->forward_comm(1,this,1,sizeof(double),T_ELECTRON,
                   gc_buf1,gc_buf2,MPI_DOUBLE);

  // output nodal temperatures for current timestep

  if ((nfileevery) && !(update->ntimestep % nfileevery)) {

    // compute atomic Ta for each grid point

    memset(&dsum[nzlo_out][nylo_out][nxlo_out],0,
           ngridout*sizeof(double));
    memset(&sum_vsq[nzlo_out][nylo_out][nxlo_out],0,
           ngridout*sizeof(double));
    memset(&sum_mass_vsq[nzlo_out][nylo_out][nxlo_out],0,
           ngridout*sizeof(double));

    double massone;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];

        ix = static_cast<int> ((x[i][0]-boxlo[0])*delxinv+shift) - OFFSET;
        iy = static_cast<int> ((x[i][1]-boxlo[1])*delyinv+shift) - OFFSET;
        iz = static_cast<int> ((x[i][2]-boxlo[2])*delzinv+shift) - OFFSET;

        dsum[ix][iy][iz] += 1.0;
        sum_vsq[ix][iy][iz] += vsq;
        sum_mass_vsq[ix][iy][iz] += massone*vsq;
      }

    // can GridComm work on a larger struct of several values?

    gc->reverse_comm(1,this,1,sizeof(double),DSUM,
                     gc_buf1,gc_buf2,MPI_DOUBLE);
    gc->reverse_comm(1,this,1,sizeof(double),SUM_VSQ,
                     gc_buf1,gc_buf2,MPI_DOUBLE);
    gc->reverse_comm(1,this,1,sizeof(double),SUM_MASS_VSQ,
                     gc_buf1,gc_buf2,MPI_DOUBLE);
    
    // NOTE: should this be a write function, work with parallel pings

    if (comm->me == 0) {
      fmt::print(fp,"{}",update->ntimestep);

      double T_a;
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++) {
            T_a = 0;
            if (dsum[ix][iy][iz] > 0)
              T_a = sum_mass_vsq[ix][iy][iz] /
                (3.0*force->boltz*dsum[ix][iy][iz]/force->mvv2e);
            fmt::print(fp," {}",T_a);
          }

      fputs("\t",fp);
      for (iz = nzlo_in; iz <= nzhi_in; iz++)
        for (iy = nylo_in; iy <= nyhi_in; iy++)
          for (ix = nxlo_in; ix <= nxhi_in; ix++)
            fmt::print(fp," {}",T_electron[ix][iy][iz]);
      fputs("\n",fp);
    }
  }
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMGrid::pack_forward_grid(int flag, void *vbuf, int nlist, int *list)
{
  double *buf = (double *) vbuf;

  double *src;
  if (flag == T_ELECTRON)
    src = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++)
    buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_forward_grid(int flag, void *vbuf, int nlist, int *list)
{
  double *buf = (double *) vbuf;

  double *dest;
  if (flag == T_ELECTRON)
    dest = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++)
    dest[list[i]] = buf[i];
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMGrid::pack_reverse_grid(int flag, void *vbuf, int nlist, int *list)
{
  double *buf = (double *) vbuf;

  double *src;
  if (flag == NET_ENERGY_TRANSFER)
    src = &net_energy_transfer[nzlo_out][nylo_out][nxlo_out];
  else if (flag == DSUM)
    src = &dsum[nzlo_out][nylo_out][nxlo_out];
  else if (flag == SUM_VSQ)
    src = &sum_vsq[nzlo_out][nylo_out][nxlo_out];
  else if (flag == SUM_MASS_VSQ)
    src = &sum_mass_vsq[nzlo_out][nylo_out][nxlo_out];
  
  for (int i = 0; i < nlist; i++)
    buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_reverse_grid(int flag, void *vbuf, int nlist, int *list)
{
  double *buf = (double *) vbuf;

  double *dest;
  if (flag == NET_ENERGY_TRANSFER)
    dest = &net_energy_transfer[nzlo_out][nylo_out][nxlo_out];
  else if (flag == DSUM)
    dest = &dsum[nzlo_out][nylo_out][nxlo_out];
  else if (flag == SUM_VSQ)
    dest = &sum_vsq[nzlo_out][nylo_out][nxlo_out];
  else if (flag == SUM_MASS_VSQ)
    dest = &sum_mass_vsq[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++)
    dest[list[i]] += buf[i];
}

/* ----------------------------------------------------------------------
   return the energy of the electronic subsystem or the net_energy transfer
   between the subsystems
------------------------------------------------------------------------- */

double FixTTMGrid::compute_vector(int n)
{
  int ix,iy,iz;

  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nznodes;
  double del_vol = dx*dy*dz;

  double e_energy_me = 0.0;
  double transfer_energy_me = 0.0;

  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        e_energy_me +=
          T_electron[ix][iy][iz]*electronic_specific_heat*
          electronic_density*del_vol;
        transfer_energy_me +=
          net_energy_transfer_all[ix][iy][iz]*update->dt;
      }

  // NOTE: only do allreduce once ?
  double e_energy,transfer_energy;
  MPI_Allreduce(&e_energy_me,&e_energy,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&transfer_energy_me,&transfer_energy,1,MPI_DOUBLE,MPI_SUM,world);

  if (n == 0) return e_energy;
  if (n == 1) return transfer_energy;
  return 0.0;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTTMGrid::write_restart(FILE *fp)
{
  int ix,iy,iz;

  double *rlist;
  memory->create(rlist,nxnodes*nynodes*nznodes+1,"TTM:rlist");

  int n = 0;
  rlist[n++] = seed;

  // NOTE: this is a parallel grid now, not a serial grid

  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
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

void FixTTMGrid::restart(char *buf)
{
  int ix,iy,iz;

  int n = 0;
  double *rlist = (double *) buf;

  // the seed must be changed from the initial seed
  // NOTE: 0.5 is whacky, could go to zero
  // NOTE: maybe GridComm should provide a method to pack a grid with bounds

  seed = static_cast<int> (0.5*rlist[n++]);

  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        T_electron[ix][iy][iz] = rlist[n++];

  delete random;
  random = new RanMars(lmp,seed+comm->me);
}

/* ----------------------------------------------------------------------
   memory usage for flangevin and 3d grid
------------------------------------------------------------------------- */

double FixTTMGrid::memory_usage()
{
  double bytes = 0.0;
  bytes += (double)3*atom->nmax * sizeof(double);
  bytes += (double)6*ngridout * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTMGrid::allocate_grid()
{
  // global indices of grid range from 0 to N-1
  // nlo_in,nhi_in = lower/upper limits of the 3d sub-brick of
  //   global grid that I own without ghost cells
  // both non-tiled and tiled proc layouts use 0-1 fractional subdomain info

  if (comm->layout != Comm::LAYOUT_TILED) {
    nxlo_in = static_cast<int> (comm->xsplit[comm->myloc[0]] * nxnodes);
    nxhi_in = static_cast<int> (comm->xsplit[comm->myloc[0]+1] * nxnodes) - 1;

    nylo_in = static_cast<int> (comm->ysplit[comm->myloc[1]] * nynodes);
    nyhi_in = static_cast<int> (comm->ysplit[comm->myloc[1]+1] * nynodes) - 1;

    nzlo_in = static_cast<int> (comm->zsplit[comm->myloc[2]] * nznodes);
    nzhi_in = static_cast<int> (comm->zsplit[comm->myloc[2]+1] * nznodes) - 1;

  } else {
    nxlo_in = static_cast<int> (comm->mysplit[0][0] * nxnodes);
    nxhi_in = static_cast<int> (comm->mysplit[0][1] * nxnodes) - 1;

    nylo_in = static_cast<int> (comm->mysplit[1][0] * nynodes);
    nyhi_in = static_cast<int> (comm->mysplit[1][1] * nynodes) - 1;

    nzlo_in = static_cast<int> (comm->mysplit[1][0] * nznodes);
    nzhi_in = static_cast<int> (comm->mysplit[1][1] * nznodes) - 1;
  }

  // nlo,nhi = min/max index of global grid pt my owned atoms can be mapped to
  // finite difference stencil requires extra grid pt around my owned grid pts
  // max of these 2 quantities is the ghost cells needed in each di
  // nlo_out,nhi_out = nlo_in,nhi_in + ghost cells

  // NOTE: need to wait until init() when neighskin is defined?

  double *boxlo = domain->boxlo;
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  shift = OFFSET + 0.0;    // change this to 0.5 for nearest grid pt
  delxinv = nxnodes/domain->xprd;
  delyinv = nxnodes/domain->yprd;
  delzinv = nxnodes/domain->zprd;

  int nlo,nhi;
  double cuthalf = 0.5*neighbor->skin;

  nlo = static_cast<int> ((sublo[0]-cuthalf-boxlo[0])*delxinv + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[0]+cuthalf-boxlo[0])*delxinv + shift) - OFFSET;
  nxlo_out = MIN(nlo,nxlo_in-1);
  nxhi_out = MAX(nhi,nxhi_in+1);

  nlo = static_cast<int> ((sublo[1]-cuthalf-boxlo[1])*delyinv + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[1]+cuthalf-boxlo[1])*delyinv + shift) - OFFSET;
  nylo_out = MIN(nlo,nylo_in-1);
  nyhi_out = MAX(nhi,nyhi_in+1);

  nlo = static_cast<int> ((sublo[2]-cuthalf-boxlo[2])*delzinv + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[2]+cuthalf-boxlo[2])*delzinv + shift) - OFFSET;
  nzlo_out = MIN(nlo,nzlo_in-1);
  nzhi_out = MAX(nhi,nzhi_in+1);

  ngridout = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) * 
    (nzhi_out-nzlo_out+1);

  gc = new GridComm(lmp,world,nxnodes,nynodes,nznodes,
                    nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                    nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out);

  gc->setup(ngc_buf1,ngc_buf2);

  memory->create(gc_buf1,ngc_buf1,"ttm/grid:gc_buf1");
  memory->create(gc_buf2,ngc_buf2,"ttm/grid:gc_buf2");

  memory->create3d_offset(nsum,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"ttm:nsum");
  memory->create3d_offset(sum_vsq,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"ttm:sum_vsq");
  memory->create3d_offset(sum_mass_vsq,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"ttm:sum_mass_vsq");
  memory->create3d_offset(T_electron_old,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"ttm:T_electron_old");
  memory->create3d_offset(T_electron,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"ttm:T_electron");
  memory->create3d_offset(net_energy_transfer,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"ttm:net_energy_transfer");
}

/* ----------------------------------------------------------------------
   initialize 3d grid quantities
------------------------------------------------------------------------- */

void FixTTMGrid::init_grid()
{
  memset(&net_energy_transfer[nzlo_out][nylo_out][nxlo_out],0,
         ngridout*sizeof(double));
}

/* ----------------------------------------------------------------------
   deallocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTMGrid::deallocate_grid()
{
  memory->destroy3d_offset(nsum,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(sum_vsq,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(sum_mass_vsq,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(T_electron_old,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(T_electron,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(net_energy_transfer,nzlo_out,nylo_out,nxlo_out);
}

/* ----------------------------------------------------------------------
   read in initial electron temperatures from a user-specified file
   only called by proc 0
------------------------------------------------------------------------- */

void FixTTMGrid::read_initial_electron_temperatures(const char *filename)
{
  int ***T_initial_set;
  memory->create(T_initial_set,nxnodes,nynodes,nznodes,"ttm:T_initial_set");
  memset(&T_initial_set[0][0][0],0,total_nnodes*sizeof(int));

  std::string name = utils::get_potential_file_path(filename);
  if (name.empty())
    error->one(FLERR,"Cannot open input file: {}",filename);
  FILE *fpr = fopen(name.c_str(),"r");

  // read initial electron temperature values from file

  char line[MAXLINE];
  int ixnode,iynode,iznode;
  double T_tmp;
  while (1) {
    if (fgets(line,MAXLINE,fpr) == nullptr) break;
    ValueTokenizer values(line);
    if (values.has_next()) ixnode = values.next_int();
    if (values.has_next()) iynode = values.next_int();
    if (values.has_next()) iznode = values.next_int();
    if (values.has_next()) T_tmp  = values.next_double();
    else error->one(FLERR,"Incorrect format in fix ttm/grid input file");

    // check correctness of input data

    if ((ixnode < 0) || (ixnode >= nxnodes)
        || (iynode < 0) || (iynode >= nynodes)
        || (iznode < 0) || (iznode >= nznodes))
      error->one(FLERR,"Fix ttm/grid invalid node index in input");

    if (T_tmp < 0.0)
      error->one(FLERR,"Fix ttm/grid electron temperatures must be > 0.0");

    T_electron[ixnode][iynode][iznode] = T_tmp;
    T_initial_set[ixnode][iynode][iznode] = 1;
  }
  fclose(fpr);

  // check completeness of input data

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        if (T_initial_set[ixnode][iynode][iznode] == 0)
          error->one(FLERR,"Initial temperatures not all set in fix ttm/grid");

  memory->destroy(T_initial_set);
}
