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

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr int MAXLINE = 256;
static constexpr int CHUNK = 1024;

#define OFFSET 16384

/* ---------------------------------------------------------------------- */

FixTTMGrid::FixTTMGrid(LAMMPS *lmp, int narg, char **arg) : 
  FixTTM(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

FixTTMGrid::~FixTTMGrid()
{
  deallocate_grid();
  deallocate_flag = 1;
}

/* ---------------------------------------------------------------------- */

void FixTTMGrid::post_constructor()
{
  // allocate 3d grid variables

  allocate_grid();

  // zero net_energy_transfer
  // in case compute_vector accesses it on timestep 0

  memset(&net_energy_transfer[nzlo_out][nylo_out][nxlo_out],0,
         ngridout*sizeof(double));

  // set initial electron temperatures from user input file

  read_electron_temperatures(infile);
}

/* ---------------------------------------------------------------------- */

void FixTTMGrid::post_force(int /*vflag*/)
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
  double dxinv = nxnodes/domain->xprd;
  double dyinv = nxnodes/domain->yprd;
  double dzinv = nxnodes/domain->zprd;

  // apply damping and thermostat to all atoms in fix group

  int flag = 0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv + shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv + shift) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv + shift) - OFFSET;

      // flag if ix,iy,iz is not within my ghost cell range

      if (ix < nxlo_out || ix > nxhi_out || 
          iy < nylo_out || iy > nyhi_out ||
          iz < nzlo_out || iz > nzhi_out) {
        flag = 1;
        continue;
      }

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

  if (flag) error->one(FLERR,"Out of range fix ttm/grid atoms");
}

/* ---------------------------------------------------------------------- */

void FixTTMGrid::end_of_step()
{
  int ix,iy,iz;

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *boxlo = domain->boxlo;
  double dxinv = nxnodes/domain->xprd;
  double dyinv = nxnodes/domain->yprd;
  double dzinv = nxnodes/domain->zprd;
  double volgrid = 1.0 / (dxinv*dyinv*dzinv);

  memset(&net_energy_transfer[nzlo_out][nylo_out][nxlo_out],0,
         ngridout*sizeof(double));

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ix = static_cast<int> ((x[i][0]-boxlo[0])*dxinv+shift) - OFFSET;
      iy = static_cast<int> ((x[i][1]-boxlo[1])*dyinv+shift) - OFFSET;
      iz = static_cast<int> ((x[i][2]-boxlo[2])*dzinv+shift) - OFFSET;
      net_energy_transfer[iz][iy][ix] +=
        (flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] +
         flangevin[i][2]*v[i][2]);
    }

  gc->reverse_comm(1,this,1,sizeof(double),0,gc_buf1,gc_buf2,MPI_DOUBLE);

  // num_inner_timesteps = # of inner steps (thermal solves)
  // required this MD step to maintain a stable explicit solve

  int num_inner_timesteps = 1;
  double inner_dt = update->dt;

  double stability_criterion = 1.0 -
    2.0*inner_dt/(electronic_specific_heat*electronic_density) *
    (electronic_thermal_conductivity*(dxinv*dxinv + dyinv*dyinv + dzinv*dzinv));

  if (stability_criterion < 0.0) {
    inner_dt = 0.5*(electronic_specific_heat*electronic_density) /
      (electronic_thermal_conductivity*(dxinv*dxinv + dyinv*dyinv + dzinv*dzinv));
    num_inner_timesteps = static_cast<int>(update->dt/inner_dt) + 1;
    inner_dt = update->dt/double(num_inner_timesteps);
    if (num_inner_timesteps > 1000000)
      error->warning(FLERR,"Too many inner timesteps in fix ttm/grid");
  }

  // finite difference iterations to update T_electron

  for (int ith_inner_timestep = 0; ith_inner_timestep < num_inner_timesteps;
       ith_inner_timestep++) {

    memcpy(&T_electron_old[nzlo_out][nylo_out][nxlo_out],
           &T_electron[nzlo_out][nylo_out][nxlo_out],ngridout*sizeof(double));

    // compute new electron T profile
 
    for (iz = nzlo_in; iz <= nzhi_in; iz++)
      for (iy = nylo_in; iy <= nyhi_in; iy++)
        for (ix = nxlo_in; ix <= nxhi_in; ix++) {
          T_electron[iz][iy][ix] =
            T_electron_old[iz][iy][ix] +
            inner_dt/(electronic_specific_heat*electronic_density) *
            (electronic_thermal_conductivity *

             ((T_electron_old[iz][iy][ix-1] + T_electron_old[iz][iy][ix+1] -
               2*T_electron_old[iz][iy][ix])*dxinv*dxinv +
              (T_electron_old[iz][iy-1][ix] + T_electron_old[iz][iy+1][ix] -
               2*T_electron_old[iz][iy][ix])*dyinv*dyinv +
              (T_electron_old[iz-1][iy][ix] + T_electron_old[iz+1][iy][ix] -
               2*T_electron_old[iz][iy][ix])*dzinv*dzinv) -
             
             net_energy_transfer[iz][iy][ix]/volgrid);
        }
  }

  // check that all temperatures are >= 0.0

  int flag = 0;

  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        if (T_electron[iz][iy][ix] < 0.0) flag = 1;

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) 
    error->all(FLERR,"Fix ttm electron temperature became negative");

  // communicate new T_electron values to ghost grid points

  gc->forward_comm(1,this,1,sizeof(double),0,gc_buf1,gc_buf2,MPI_DOUBLE);

  // assign electron temperature to each atom for fix output


}

/* ----------------------------------------------------------------------
   read electron temperatures on grid from a user-specified file
   proc 0 reads one chunk at a time, broadcasts to other procs
   each proc stores values for grid points it owns
------------------------------------------------------------------------- */

void FixTTMGrid::read_electron_temperatures(const char *filename)
{
  int i,j,ix,iy,iz,nchunk,eof;

  int me = comm->me;

  // initialize my own+ghost grid values to zero

  memset(&T_electron[nzlo_out][nylo_out][nxlo_out],0,ngridout*sizeof(double));

  // proc 0 opens file

  FILE *fp = nullptr;

  if (me == 0) {
    std::string name = utils::get_potential_file_path(filename);
    if (name.empty()) error->one(FLERR,"Cannot open input file: {}",filename);
    fp = fopen(name.c_str(),"r");
  }

  // read electron temperature values from file, one chunk at a time
 
  char **values = new char*[4];
  char *buffer = new char[CHUNK*MAXLINE];
  bigint ntotal = (bigint) nxnodes * nynodes * nznodes;
  bigint nread = 0;
  char *buf,*next;

  while (nread < ntotal) {
    nchunk = MIN(ntotal-nread,CHUNK);
    eof = utils::read_lines_from_file(fp,nchunk,MAXLINE,buffer,me,world);
    if (eof) error->all(FLERR,"Unexpected end of data file");

    buf = buffer;
    next = strchr(buf,'\n');
    *next = '\0';
    int nwords = utils::trim_and_count_words(buf);
    *next = '\n';

    if (nwords != 4) error->all(FLERR,"Incorrect format in fix ttm data file");
    
    // loop over lines of grid point values
    // tokenize the line into ix,iy,iz grid index plus temperature value
    // if I own grid point, store the value
    
    for (i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');
      
      for (j = 0; j < nwords; j++) {
        buf += strspn(buf," \t\n\r\f");
        buf[strcspn(buf," \t\n\r\f")] = '\0';
        values[j] = buf;
        buf += strlen(buf)+1;
      }
    
      ix = utils::inumeric(FLERR,values[0],false,lmp);
      iy = utils::inumeric(FLERR,values[1],false,lmp);
      iz = utils::inumeric(FLERR,values[2],false,lmp);

      if (ix < 0 || ix >= nxnodes || iy < 0 || iy >= nynodes ||
          iz < 0 || iz >= nznodes)
        error->all(FLERR,"Fix ttm/grid invalid node index in input");
      
      if (ix >= nxlo_in && ix <= nxhi_in &&
          iy >= nylo_in && iy <= nyhi_in &&
          iz >= nzlo_in && iz <= nzhi_in)
        T_electron[iz][iy][ix] = utils::numeric(FLERR,values[3],true,lmp);
      
      buf = next + 1;
    }

    nread += nchunk;
  }

  // close file

  if (me == 0) fclose(fp);

  // clean up

  delete [] values;
  delete [] buffer;

  // check that all owned temperature values are > 0.0

  int flag = 0;

  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++)
        if (T_electron[iz][iy][ix] <= 0.0) flag = 1;

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) 
    error->all(FLERR,
               "Fix ttm infile did not set all temperatures or some are <= 0.0");

  // communicate new T_electron values to ghost grid points

  gc->forward_comm(1,this,1,sizeof(double),0,gc_buf1,gc_buf2,MPI_DOUBLE);
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMGrid::pack_forward_grid(int flag, void *vbuf, int nlist, int *list)
{
  double *buf = (double *) vbuf;
  double *src = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++)
    buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_forward_grid(int flag, void *vbuf, int nlist, int *list)
{
  double *buf = (double *) vbuf;
  double *dest = &T_electron[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++)
    dest[list[i]] = buf[i];
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void FixTTMGrid::pack_reverse_grid(int flag, void *vbuf, int nlist, int *list)
{
  double *buf = (double *) vbuf;
  double *src = &net_energy_transfer[nzlo_out][nylo_out][nxlo_out];
  
  for (int i = 0; i < nlist; i++)
    buf[i] = src[list[i]];
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void FixTTMGrid::unpack_reverse_grid(int flag, void *vbuf, int nlist, int *list)
{
  double *buf = (double *) vbuf;
  double *dest = &net_energy_transfer[nzlo_out][nylo_out][nxlo_out];

  for (int i = 0; i < nlist; i++)
    dest[list[i]] += buf[i];
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

  double *boxlo = domain->boxlo;
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;
  double dxinv = nxnodes/domain->xprd;
  double dyinv = nxnodes/domain->yprd;
  double dzinv = nxnodes/domain->zprd;

  shift = OFFSET + 0.0;    // change this to 0.5 for nearest grid pt

  int nlo,nhi;
  double cuthalf = 0.5*neighbor->skin;

  nlo = static_cast<int> ((sublo[0]-cuthalf-boxlo[0])*dxinv + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[0]+cuthalf-boxlo[0])*dxinv + shift) - OFFSET;
  nxlo_out = MIN(nlo,nxlo_in-1);
  nxhi_out = MAX(nhi,nxhi_in+1);

  nlo = static_cast<int> ((sublo[1]-cuthalf-boxlo[1])*dyinv + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[1]+cuthalf-boxlo[1])*dyinv + shift) - OFFSET;
  nylo_out = MIN(nlo,nylo_in-1);
  nyhi_out = MAX(nhi,nyhi_in+1);

  nlo = static_cast<int> ((sublo[2]-cuthalf-boxlo[2])*dzinv + shift) - OFFSET;
  nhi = static_cast<int> ((subhi[2]+cuthalf-boxlo[2])*dzinv + shift) - OFFSET;
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

  memory->create3d_offset(T_electron_old,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"ttm:T_electron_old");
  memory->create3d_offset(T_electron,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"ttm:T_electron");
  memory->create3d_offset(net_energy_transfer,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"ttm:net_energy_transfer");
}

/* ----------------------------------------------------------------------
   deallocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTMGrid::deallocate_grid()
{
  delete gc;
  memory->destroy(gc_buf1);
  memory->destroy(gc_buf2);

  memory->destroy3d_offset(T_electron_old,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(T_electron,nzlo_out,nylo_out,nxlo_out);
  memory->destroy3d_offset(net_energy_transfer,nzlo_out,nylo_out,nxlo_out);
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
        T_electron[iz][iy][ix] = rlist[n++];

  delete random;
  random = new RanMars(lmp,seed+comm->me);
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
        rlist[n++] =  T_electron[iz][iy][ix];

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(rlist,sizeof(double),n,fp);
  }

  memory->destroy(rlist);
}

/* ----------------------------------------------------------------------
   return the energy of the electronic subsystem 
   or the net_energy transfer between the subsystems
------------------------------------------------------------------------- */

double FixTTMGrid::compute_vector(int n)
{
  int ix,iy,iz;

  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nznodes;
  double volgrid = dx*dy*dz;

  double e_energy_me = 0.0;
  double transfer_energy_me = 0.0;

  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        e_energy_me +=
          T_electron[iz][iy][ix]*electronic_specific_heat*
          electronic_density*volgrid;
        transfer_energy_me +=
          net_energy_transfer[iz][iy][ix]*update->dt;
      }

  double e_energy,transfer_energy;
  MPI_Allreduce(&e_energy_me,&e_energy,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&transfer_energy_me,&transfer_energy,1,MPI_DOUBLE,MPI_SUM,world);

  if (n == 0) return e_energy;
  if (n == 1) return transfer_energy;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage for flangevin and 3d grid
------------------------------------------------------------------------- */

double FixTTMGrid::memory_usage()
{
  double bytes = 0.0;
  bytes += (double)3*atom->nmax * sizeof(double);
  bytes += (double)3*ngridout * sizeof(double);
  return bytes;
}

