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
#include "pair_eam_lj_gpu.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "math_const.h"
#include "error.h"
#include "neigh_request.h"
#include "gpu_extra.h"

using namespace LAMMPS_NS;
using namespace MathConst;

// External functions from cuda library for atom decomposition

int eam_lj_gpu_init(const int ntypes, double host_cutforcesq, 
                 int **host_type2rhor, int **host_type2z2r, int *host_type2frho,
                 double ***host_rhor_spline, double ***host_z2r_spline,
                 double ***host_frho_spline, double **host_cutljsq, 
                 double **host_lj1, double **host_lj2, double **host_lj3, 
                 double **host_lj4, double **offset, double *special_lj,
                 double rdr, double rdrho, int nrhor, 
                 int nrho, int nz2r, int nfrho, int nr, 
                 const int nlocal, const int nall, const int max_nbors, 
                 const int maxspecial, const double cell_size, 
                 int &gpu_mode, FILE *screen, int &fp_size);
void eam_lj_gpu_clear();
int** eam_lj_gpu_compute_n(const int ago, const int inum_full,
                 const int nall, double **host_x, int *host_type,
                 double *sublo, double *subhi, int *tag, int **nspecial, 
                 int **special, const bool eflag, const bool vflag,
                 const bool eatom, const bool vatom, int &host_start,
                 int **ilist, int **jnum,  const double cpu_time,
                 bool &success, int &inum, void **fp_ptr);
void eam_lj_gpu_compute(const int ago, const int inum_full, const int nlocal, 
                 const int nall,double **host_x, int *host_type, 
                 int *ilist, int *numj, int **firstneigh, 
                 const bool eflag, const bool vflag,
                 const bool eatom, const bool vatom, int &host_start,
                 const double cpu_time, bool &success, void **fp_ptr);
void eam_lj_gpu_compute_force(int *ilist, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom);

double eam_lj_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairEAMLJGPU::PairEAMLJGPU(LAMMPS *lmp) : PairEAM(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  cpu_time = 0.0;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairEAMLJGPU::~PairEAMLJGPU()
{
  eam_lj_gpu_clear();

  memory->destroy(cut);
  memory->destroy(epsilon);
  memory->destroy(sigma);
  memory->destroy(lj1);
  memory->destroy(lj2);
  memory->destroy(lj3);
  memory->destroy(lj4);
  memory->destroy(offset);
}

/* ---------------------------------------------------------------------- */

double PairEAMLJGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + eam_lj_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairEAMLJGPU::compute(int eflag, int vflag)
{
  int i,j,ii,jj,m,jnum,itype,jtype;
  double evdwl,*coeff;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;
 
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // compute density on each atom on GPU

  int nall = atom->nlocal + atom->nghost;  
  int inum, host_start, inum_dev;
  
  bool success = true;
  int *ilist, *numneigh, **firstneigh; 
  if (gpu_mode != GPU_FORCE) { 
    inum = atom->nlocal;
    firstneigh = eam_lj_gpu_compute_n(neighbor->ago, inum, nall, 
             atom->x,atom->type, domain->sublo, domain->subhi,
             atom->tag, atom->nspecial, atom->special,
             eflag, vflag, eflag_atom, vflag_atom,
             host_start, &ilist, &numneigh, cpu_time,
             success, inum_dev, &fp_pinned);
  } else { // gpu_mode == GPU_FORCE
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    eam_lj_gpu_compute(neighbor->ago, inum, nlocal, nall, atom->x, atom->type,
		    ilist, numneigh, firstneigh, eflag, vflag, eflag_atom,
		    vflag_atom, host_start, cpu_time, success, &fp_pinned);
  }
    
  if (!success)
    error->one(FLERR,"Out of memory on GPGPU");

  if (host_start<inum) {
    cpu_time = MPI_Wtime();
    cpu_compute_energy(host_start, inum, eflag, vflag, 
            ilist, numneigh, firstneigh);
    cpu_time = MPI_Wtime() - cpu_time;
  }
  
  // communicate derivative of embedding function

  comm->forward_comm_pair(this);
    
  // compute forces on each atom on GPU
  if (gpu_mode != GPU_FORCE) 
    eam_lj_gpu_compute_force(NULL, eflag, vflag, eflag_atom, vflag_atom);
  else
    eam_lj_gpu_compute_force(ilist, eflag, vflag, eflag_atom, vflag_atom);

  if (host_start<inum) {
    double cpu_time2 = MPI_Wtime();
    cpu_compute(host_start, inum, eflag, vflag, ilist, numneigh, firstneigh);
    cpu_time += MPI_Wtime() - cpu_time2;
  }
  
}

void PairEAMLJGPU::cpu_compute_energy(int start, int inum, int eflag, int vflag,
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

void PairEAMLJGPU::cpu_compute(int start, int inum, int eflag, int vflag,
				      int *ilist, int *numneigh,
				      int **firstneigh)
{
  int i,j,ii,jj,m,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
  double r2inv,r6inv,factor_lj,force_eam,force_lj;
  double *coeff;
  int *jlist;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  
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
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      
      if (rsq < cutsq[itype][jtype]) {
        jtype = type[j];
        if (rsq < cutforcesq  && (itype ==2 && jtype ==2)) {
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
          force_eam = -psip*recip;
        } else {
          force_eam=0.0;
          phi = 0.0;
        }
        
      if (rsq < cutsq[itype][jtype] && (itype !=2 || jtype !=2)) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        force_lj = factor_lj*r2inv*r6inv * (lj1[itype][jtype]*r6inv 
          - lj2[itype][jtype]);
      } else force_lj=0.0;
  
  fpair = force_eam + force_lj;
	f[i][0] += delx*fpair;
	f[i][1] += dely*fpair;
	f[i][2] += delz*fpair;

	if (eflag) {
    evdwl = phi;
    if (rsq < cutsq[itype][jtype] && (itype !=2 || jtype !=2)) {
      double e = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
        offset[itype][jtype];
      evdwl += factor_lj*e;  
    }
    }
	if (evflag) ev_tally_full(i,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
    
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairEAMLJGPU::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
  for (int i = 1; i <= n; i++) map[i] = -1;

  type2frho = new int[n+1];
  memory->create(type2rhor,n+1,n+1,"pair:type2rhor");
  memory->create(type2z2r,n+1,n+1,"pair:type2z2r");
  
  // LJ
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEAMLJGPU::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style eam/lj/gpu command");

  cut_global = force->numeric(arg[0]);

  // reset cutoffs that have been explicitly set
  
  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read DYNAMO funcfl file
------------------------------------------------------------------------- */

void PairEAMLJGPU::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  
  if (narg > 5) 
    error->all(FLERR,"Incorrect args for pair coefficients");

  // parse pair of atom types

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);
  
  int count = 0;
  if (narg == 3) { // eam
 
    // read funcfl file if hasn't already been read
    // store filename in Funcfl data struct

    int ifuncfl;
    for (ifuncfl = 0; ifuncfl < nfuncfl; ifuncfl++)
      if (strcmp(arg[2],funcfl[ifuncfl].file) == 0) break;

    if (ifuncfl == nfuncfl) {
      nfuncfl++;
      funcfl = (Funcfl *) 
        memory->srealloc(funcfl,nfuncfl*sizeof(Funcfl),"pair:funcfl");
      read_file(arg[2]);
      int n = strlen(arg[2]) + 1;
      funcfl[ifuncfl].file = new char[n];
      strcpy(funcfl[ifuncfl].file,arg[2]);
    }

    // set setflag and map only for i,i type pairs
    // set mass of atom type if i = j

    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        if (i == j) {
	  setflag[i][i] = 1;
	  map[i] = ifuncfl;
	  atom->set_mass(i,funcfl[ifuncfl].mass);
	  count++;
         }
      }
    }
  } else if (narg >= 4 || narg <= 5) { // LJ
    double epsilon_one = force->numeric(arg[2]);
    double sigma_one = force->numeric(arg[3]);

    double cut_one = cut_global;
    if (narg == 5) cut_one = force->numeric(arg[4]);
  
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        epsilon[i][j] = epsilon_one;
        sigma[i][j] = sigma_one;
        cut[i][j] = cut_one;
        setflag[i][j] = 1;
        count++;
      }
    }
  }
  
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEAMLJGPU::init_style()
{
  if (force->newton_pair) 
    error->all(FLERR,"Cannot use newton pair with eam/lj/gpu pair style");
  
  if (!allocated) error->all(FLERR,"Not allocate memory eam/lj/gpu pair style");
  
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
  int fp_size;
  int success = eam_lj_gpu_init(atom->ntypes+1, cutforcesq,
          type2rhor, type2z2r, type2frho,
          rhor_spline, z2r_spline, frho_spline, cutsq,
          lj1, lj2, lj3, lj4, offset, force->special_lj,
          rdr, rdrho, nrhor, nrho, nz2r, nfrho, nr, atom->nlocal, 
          atom->nlocal+atom->nghost, 300, maxspecial,
			     cell_size, gpu_mode, screen, fp_size);
  GPU_EXTRA::check_flag(success,error,world);
  
  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }

  if (fp_size == sizeof(double))
    fp_single = false;
  else
    fp_single = true;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEAMLJGPU::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
			       sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  if (i==2 && j==2) {
    // single global cutoff = max of cut from all files read in
    // for funcfl could be multiple files
    // for setfl or fs, just one file

    if (funcfl) {
      cutmax = 0.0;
      for (int m = 0; m < nfuncfl; m++)
        cutmax = MAX(cutmax,funcfl[m].cut);
    } else if (setfl) cutmax = setfl->cut;
    else if (fs) cutmax = fs->cut;

    cutforcesq = cutmax*cutmax;
    cut[i][j] = cutforcesq;
    return cutmax;
  }
  
  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  if (offset_flag) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);
        
    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] * 
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9); 
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] * 
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9); 
  } 
  
  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

int PairEAMLJGPU::pack_comm(int n, int *list, double *buf, int pbc_flag, 
			  int *pbc)
{
  int i,j,m;

  m = 0;

  if (fp_single) {
    float *fp_ptr = (float *)fp_pinned;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = static_cast<double>(fp_ptr[j]);
    }
  } else {
    double *fp_ptr = (double *)fp_pinned;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = fp_ptr[j];
    }
  }

  return 1;
}

/* ---------------------------------------------------------------------- */

void PairEAMLJGPU::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (fp_single) {
    float *fp_ptr = (float *)fp_pinned;
    for (i = first; i < last; i++) fp_ptr[i] = buf[m++];
  } else {
    double *fp_ptr = (double *)fp_pinned;
    for (i = first; i < last; i++) fp_ptr[i] = buf[m++];
  }
}
