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
   Contributing authors: Trung Dac Nguyen, W. Michael Brown (ORNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_eam_gpu.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "neigh_request.h"
#include "gpu_extra.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MAXLINE 1024

// External functions from cuda library for atom decomposition

int eam_gpu_init(const int ntypes, double host_cutforcesq,
                 int **host_type2rhor, int **host_type2z2r,
                 int *host_type2frho, 
                 double ***host_rhor_spline, double ***host_z2r_spline,
                 double ***host_frho_spline,
                 double rdr, double rdrho, int nrhor, int nrho, 
                 int nz2r, int nfrho, int nr,
                 const int nlocal, const int nall, const int max_nbors, 
                 const int maxspecial, const double cell_size, 
                 int &gpu_mode, FILE *screen);
void eam_gpu_clear();
int** eam_gpu_compute_energy_n(const int ago, const int inum_full,
                         const int nall, double **host_x, int *host_type,
                         double *sublo, double *subhi, int *tag, int **nspecial, 
                         int **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum,  const double cpu_time,
                         bool &success, double *host_fp, double *boxlo,
                         double *prd, int &inum);
void eam_gpu_compute_energy(const int ago, const int inum_full, const int nall,
                      double **host_x, int *host_type, int *ilist, int *numj,
                      int **firstneigh, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom, int &host_start,
                      const double cpu_time, bool &success, double *host_fp,
                      const int nlocal, double *boxlo, double *prd);
void eam_gpu_compute_n(const int ago, const int inum_full,
                         const int nall, double **host_x, int *host_type,
                         double *sublo, double *subhi, int *tag, int **nspecial, 
                         int **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum,  const double cpu_time,
                         bool &success, double *host_fp, double *boxlo,
                         double *prd, int inum);
void eam_gpu_compute(const int ago, const int inum_full, const int nall,
                      double **host_x, int *host_type, int *ilist, int *numj,
                      int **firstneigh, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom, int &host_start,
                      const double cpu_time, bool &success, double *host_fp,
                      const int nlocal, double *boxlo, double *prd);
double eam_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairEAMGPU::PairEAMGPU(LAMMPS *lmp) : PairEAM(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  cpu_time = 0.0;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairEAMGPU::~PairEAMGPU()
{
  eam_gpu_clear();
}

/* ---------------------------------------------------------------------- */

double PairEAMGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + eam_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairEAMGPU::compute(int eflag, int vflag)
{
  int i,j,ii,jj,m,jnum,itype,jtype;
  double evdwl,*coeff;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;
 
  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length
    
  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(rho,nmax,"pair:rho");
    memory->create(fp,nmax,"pair:fp");
  }

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // zero out density

  if (newton_pair) {
    m = nlocal + atom->nghost;
    for (i = 0; i < m; i++) rho[i] = 0.0; 
  } else for (i = 0; i < nlocal; i++) rho[i] = 0.0; 

  
  // compute density on each atom on GPU

  int nall = atom->nlocal + atom->nghost;  
  int inum, host_start, inum_dev;
  
  bool success = true;
  int *ilist, *numneigh, **firstneigh; 
  if (gpu_mode != GPU_FORCE) { 
    inum = atom->nlocal;
    
    firstneigh = eam_gpu_compute_energy_n(neighbor->ago, inum, nall, atom->x,
             atom->type, domain->sublo, domain->subhi,
             atom->tag, atom->nspecial, atom->special,
             eflag, vflag, eflag_atom, vflag_atom,
             host_start, &ilist, &numneigh, cpu_time,
             success, fp, domain->boxlo, 
             domain->prd, inum_dev);
  } else { // gpu_mode == GPU_FORCE
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    eam_gpu_compute_energy(neighbor->ago, inum, nall, atom->x, atom->type,
		    ilist, numneigh, firstneigh, eflag, vflag, eflag_atom,
		    vflag_atom, host_start, cpu_time, success, fp,
		    atom->nlocal, domain->boxlo, domain->prd);
  }
    
  if (!success)
    error->one(FLERR,"Out of memory on GPGPU");

  if (host_start<inum) {
    cpu_time = MPI_Wtime();
    cpu_compute_energy(host_start, inum, eflag, vflag, ilist, numneigh, firstneigh);
    cpu_time = MPI_Wtime() - cpu_time;
  }
  
  // communicate derivative of embedding function

  comm->forward_comm_pair(this);
    
  // compute forces on each atom on GPU

  
  if (gpu_mode != GPU_FORCE) {
    eam_gpu_compute_n(neighbor->ago, inum, nall, atom->x,
				   atom->type, domain->sublo, domain->subhi,
				   atom->tag, atom->nspecial, atom->special,
				   eflag, vflag, eflag_atom, vflag_atom,
				   host_start, &ilist, &numneigh, cpu_time,
				   success, fp, domain->boxlo, 
				   domain->prd, inum_dev);
  } else { // gpu_mode == GPU_FORCE
    eam_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
		    ilist, numneigh, firstneigh, eflag, vflag, eflag_atom,
		    vflag_atom, host_start, cpu_time, success, fp,
		    atom->nlocal, domain->boxlo, domain->prd);
  }
  
  if (host_start<inum) {
    double cpu_time2 = MPI_Wtime();
    cpu_compute(host_start, inum, eflag, vflag, ilist, numneigh, firstneigh);
    cpu_time += MPI_Wtime() - cpu_time2;
  }
  
}

void PairEAMGPU::cpu_compute_energy(int start, int inum, int eflag, int vflag,
				      int *ilist, int *numneigh,
				      int **firstneigh) 
{
  int i,j,ii,jj,m,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,phi;
  double *coeff;
  int *jlist;
  
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  
  // rho = density at each atom
  // loop over neighbors of my atoms

  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
	jtype = type[j];
	p = sqrt(rsq)*rdr + 1.0;
	m = static_cast<int> (p);
	m = MIN(m,nr-1);
	p -= m;
	p = MIN(p,1.0);
	coeff = rhor_spline[type2rhor[jtype][itype]][m];
	rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
	    }
    }
  }

  // communicate and sum densities

  if (newton_pair) comm->reverse_comm_pair(this);
  
  
  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  
  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    p = rho[i]*rdrho + 1.0;    
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    if (eflag) {
      phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      if (eflag_global) eng_vdwl += phi;
      if (eflag_atom) eatom[i] += phi;
    
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairEAMGPU::cpu_compute(int start, int inum, int eflag, int vflag,
				      int *ilist, int *numneigh,
				      int **firstneigh)
{
  int i,j,ii,jj,m,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
  double *coeff;
  int *jlist;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  
  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
	jtype = type[j];
	r = sqrt(rsq);
	p = r*rdr + 1.0;
	m = static_cast<int> (p);
	m = MIN(m,nr-1);
	p -= m;
	p = MIN(p,1.0);

	// rhoip = derivative of (density at atom j due to atom i)
	// rhojp = derivative of (density at atom i due to atom j)
	// phi = pair potential energy
	// phip = phi'
	// z2 = phi * r
	// z2p = (phi * r)' = (phi' r) + phi
	// psip needs both fp[i] and fp[j] terms since r_ij appears in two
	//   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
	//   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

	coeff = rhor_spline[type2rhor[itype][jtype]][m];
	rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
	coeff = rhor_spline[type2rhor[jtype][itype]][m];
	rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
	coeff = z2r_spline[type2z2r[itype][jtype]][m];
	z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
	z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

	recip = 1.0/r;
	phi = z2*recip;
	phip = z2p*recip - phi*recip;
	psip = fp[i]*rhojp + fp[j]*rhoip + phip;
	fpair = -psip*recip;

	f[i][0] += delx*fpair;
	f[i][1] += dely*fpair;
	f[i][2] += delz*fpair;

	if (eflag) evdwl = phi;
	if (evflag) ev_tally_full(i,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEAMGPU::init_style()
{
  if (force->newton_pair) 
    error->all(FLERR,"Cannot use newton pair with eam/gpu pair style");
  
  if (!allocated) error->all(FLERR,"Not allocate memory eam/gpu pair style");
  
  // convert read-in file(s) to arrays and spline them

  file2array();
  array2spline();
  
  // Repeat cutsq calculation because done after call to init_style
  double maxcut = -1.0;
  double cut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i,j);
        cut *= cut;
        if (cut > maxcut)
          maxcut = cut;
        cutsq[i][j] = cutsq[j][i] = cut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }
  double cell_size = sqrt(maxcut) + neighbor->skin;
  
  int maxspecial=0;
  if (atom->molecular)
    maxspecial=atom->maxspecial;
  int success = eam_gpu_init(atom->ntypes+1, cutforcesq,
          type2rhor, type2z2r, type2frho,
          rhor_spline, z2r_spline, frho_spline,
          rdr, rdrho, nrhor, nrho, nz2r, nfrho, nr, atom->nlocal, 
          atom->nlocal+atom->nghost, 300, maxspecial,
          cell_size, gpu_mode, screen);
  GPU_EXTRA::check_flag(success,error,world);
  
  if (gpu_mode == GPU_FORCE) {
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  }
}




