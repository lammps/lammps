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
   Contributing authors: Trung Dac Nguyen (ORNL), W. Michael Brown (ORNL)
------------------------------------------------------------------------- */

#include "pair_eam_fs_gpu.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "neigh_request.h"
#include "gpu_extra.h"
#include "domain.h"
#include "utils.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

// External functions from cuda library for atom decomposition

int eam_fs_gpu_init(const int ntypes, double host_cutforcesq,
                 int **host_type2rhor, int **host_type2z2r,
                 int *host_type2frho, double ***host_rhor_spline,
                 double ***host_z2r_spline, double ***host_frho_spline,
                 double rdr, double rdrho, double rhomax,
                 int nrhor, int nrho, int nz2r, int nfrho, int nr,
                 const int nlocal, const int nall, const int max_nbors,
                 const int maxspecial, const double cell_size, int &gpu_mode,
                 FILE *screen, int &fp_size);
void eam_fs_gpu_clear();
int** eam_fs_gpu_compute_n(const int ago, const int inum_full, const int nall,
                        double **host_x, int *host_type, double *sublo,
                        double *subhi, tagint *tag, int **nspecial, tagint **special,
                        const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int &host_start, int **ilist,
                        int **jnum,  const double cpu_time, bool &success,
                        int &inum, void **fp_ptr);
void eam_fs_gpu_compute(const int ago, const int inum_full, const int nlocal,
                     const int nall,double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh,
                     const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success, void **fp_ptr);
void eam_fs_gpu_compute_force(int *ilist, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom);
double eam_fs_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairEAMFSGPU::PairEAMFSGPU(LAMMPS *lmp) : PairEAM(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ---------------------------------------------------------------------- */

PairEAMFSGPU::~PairEAMFSGPU()
{
  eam_fs_gpu_clear();
}

/* ---------------------------------------------------------------------- */

double PairEAMFSGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + eam_fs_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairEAMFSGPU::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  // compute density on each atom on GPU

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int inum, host_start, inum_dev;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = eam_fs_gpu_compute_n(neighbor->ago, inum, nall, atom->x,
                                   atom->type, domain->sublo, domain->subhi,
                                   atom->tag, atom->nspecial, atom->special,
                                   eflag, vflag, eflag_atom, vflag_atom,
                                   host_start, &ilist, &numneigh, cpu_time,
                                   success, inum_dev, &fp_pinned);
  } else { // gpu_mode == GPU_FORCE
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    eam_fs_gpu_compute(neighbor->ago, inum, nlocal, nall, atom->x, atom->type,
                    ilist, numneigh, firstneigh, eflag, vflag, eflag_atom,
                    vflag_atom, host_start, cpu_time, success, &fp_pinned);
  }

  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");

  // communicate derivative of embedding function

  comm->forward_comm_pair(this);

  // compute forces on each atom on GPU
  if (gpu_mode != GPU_FORCE)
    eam_fs_gpu_compute_force(NULL, eflag, vflag, eflag_atom, vflag_atom);
  else
    eam_fs_gpu_compute_force(ilist, eflag, vflag, eflag_atom, vflag_atom);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEAMFSGPU::init_style()
{
  if (force->newton_pair)
    error->all(FLERR,"Cannot use newton pair with eam/fs/gpu pair style");

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
  int success = eam_fs_gpu_init(atom->ntypes+1, cutforcesq, type2rhor, type2z2r,
                             type2frho, rhor_spline, z2r_spline, frho_spline,
                             rdr, rdrho, rhomax, nrhor, nrho, nz2r, nfrho, nr,
                             atom->nlocal, atom->nlocal+atom->nghost, 300,
                             maxspecial, cell_size, gpu_mode, screen, fp_size);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }

  if (fp_size == sizeof(double))
    fp_single = false;
  else
    fp_single = true;
}

/* ---------------------------------------------------------------------- */

double PairEAMFSGPU::single(int i, int j, int itype, int jtype,
                            double rsq, double /* factor_coul */,
                            double /* factor_lj */, double &fforce)
{
  int m;
  double r,p,rhoip,rhojp,z2,z2p,recip,phi,phip,psip;
  double *coeff;

  r = sqrt(rsq);
  p = r*rdr + 1.0;
  m = static_cast<int> (p);
  m = MIN(m,nr-1);
  p -= m;
  p = MIN(p,1.0);

  coeff = rhor_spline[type2rhor[itype][jtype]][m];
  rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = rhor_spline[type2rhor[jtype][itype]][m];
  rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = z2r_spline[type2z2r[itype][jtype]][m];
  z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
  z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

  double fp_i,fp_j;
  if (fp_single == false) {
    fp_i = ((double*)fp_pinned)[i];
    fp_j = ((double*)fp_pinned)[j];
  } else {
    fp_i = ((float*)fp_pinned)[i];
    fp_j = ((float*)fp_pinned)[j];
  }

  recip = 1.0/r;
  phi = z2*recip;
  phip = z2p*recip - phi*recip;
  psip = fp_i*rhojp + fp_j*rhoip + phip;
  fforce = -psip*recip;

  return phi;
}

/* ---------------------------------------------------------------------- */

int PairEAMFSGPU::pack_forward_comm(int n, int *list, double *buf,
                                    int /* pbc_flag */, int * /* pbc */)
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

  return m;
}

/* ---------------------------------------------------------------------- */

void PairEAMFSGPU::unpack_forward_comm(int n, int first, double *buf)
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

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read EAM Finnis-Sinclair file
------------------------------------------------------------------------- */

void PairEAMFSGPU::coeff(int narg, char **arg)
{
  int i,j;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read EAM Finnis-Sinclair file

  if (fs) {
    for (i = 0; i < fs->nelements; i++) delete [] fs->elements[i];
    delete [] fs->elements;
    delete [] fs->mass;
    memory->destroy(fs->frho);
    memory->destroy(fs->rhor);
    memory->destroy(fs->z2r);
    delete fs;
  }
  fs = new Fs();
  read_file(arg[2]);

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL

  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < fs->nelements; j++)
      if (strcmp(arg[i],fs->elements[j]) == 0) break;
    if (j < fs->nelements) map[i-2] = j;
    else error->all(FLERR,"No matching element in EAM potential file");
  }

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass of atom type if i = j

  int count = 0;
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(FLERR,i,fs->mass[map[i]]);
        count++;
      }
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   read a multi-element DYNAMO setfl file
------------------------------------------------------------------------- */

void PairEAMFSGPU::read_file(char *filename)
{
  Fs *file = fs;

  // open potential file

  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];
  char *r_token;

  if (me == 0) {
    fptr = force->open_potential(filename);
    if (fptr == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open EAM potential file %s",filename);
      error->one(FLERR,str);
    }
  }

  // read and broadcast header
  // extract element names from nelements line

  int n;
  if (me == 0) {
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);

  sscanf(line,"%d",&file->nelements);
  int nwords = atom->count_words(line);
  if (nwords != file->nelements + 1)
    error->all(FLERR,"Incorrect element names in EAM potential file");

  char **words = new char*[file->nelements+1];
  r_token = line;
  nwords = 0;
  utils::strtok_r(r_token," \t\n\r\f",&r_token);
  while ((words[nwords++] = utils::strtok_r(NULL," \t\n\r\f",&r_token))) continue;

  file->elements = new char*[file->nelements];
  for (int i = 0; i < file->nelements; i++) {
    n = strlen(words[i]) + 1;
    file->elements[i] = new char[n];
    strcpy(file->elements[i],words[i]);
  }
  delete [] words;

  if (me == 0) {
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%d %lg %d %lg %lg",
           &file->nrho,&file->drho,&file->nr,&file->dr,&file->cut);
  }

  MPI_Bcast(&file->nrho,1,MPI_INT,0,world);
  MPI_Bcast(&file->drho,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->nr,1,MPI_INT,0,world);
  MPI_Bcast(&file->dr,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->cut,1,MPI_DOUBLE,0,world);

  file->mass = new double[file->nelements];
  memory->create(file->frho,file->nelements,file->nrho+1,
                                              "pair:frho");
  memory->create(file->rhor,file->nelements,file->nelements,
                 file->nr+1,"pair:rhor");
  memory->create(file->z2r,file->nelements,file->nelements,
                 file->nr+1,"pair:z2r");

  int i,j,tmp;
  for (i = 0; i < file->nelements; i++) {
    if (me == 0) {
      fgets(line,MAXLINE,fptr);
      sscanf(line,"%d %lg",&tmp,&file->mass[i]);
    }
    MPI_Bcast(&file->mass[i],1,MPI_DOUBLE,0,world);

    if (me == 0) grab(fptr,file->nrho,&file->frho[i][1]);
    MPI_Bcast(&file->frho[i][1],file->nrho,MPI_DOUBLE,0,world);

    for (j = 0; j < file->nelements; j++) {
      if (me == 0) grab(fptr,file->nr,&file->rhor[i][j][1]);
      MPI_Bcast(&file->rhor[i][j][1],file->nr,MPI_DOUBLE,0,world);
    }
  }

  for (i = 0; i < file->nelements; i++)
    for (j = 0; j <= i; j++) {
      if (me == 0) grab(fptr,file->nr,&file->z2r[i][j][1]);
      MPI_Bcast(&file->z2r[i][j][1],file->nr,MPI_DOUBLE,0,world);
    }

  // close the potential file

  if (me == 0) fclose(fptr);
}

/* ----------------------------------------------------------------------
   copy read-in setfl potential to standard array format
------------------------------------------------------------------------- */

void PairEAMFSGPU::file2array()
{
  int i,j,m,n;
  int ntypes = atom->ntypes;

  // set function params directly from fs file

  nrho = fs->nrho;
  nr = fs->nr;
  drho = fs->drho;
  dr = fs->dr;
  rhomax = (nrho-1) * drho;

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of fs elements + 1 for zero array

  nfrho = fs->nelements + 1;
  memory->destroy(frho);
  memory->create(frho,nfrho,nrho+1,"pair:frho");

  // copy each element's frho to global frho

  for (i = 0; i < fs->nelements; i++)
    for (m = 1; m <= nrho; m++) frho[i][m] = fs->frho[i][m];

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho-1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to element (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0) type2frho[i] = map[i];
    else type2frho[i] = nfrho-1;

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = square of # of fs elements

  nrhor = fs->nelements * fs->nelements;
  memory->destroy(rhor);
  memory->create(rhor,nrhor,nr+1,"pair:rhor");

  // copy each element pair rhor to global rhor

  n = 0;
  for (i = 0; i < fs->nelements; i++)
    for (j = 0; j < fs->nelements; j++) {
      for (m = 1; m <= nr; m++) rhor[n][m] = fs->rhor[i][j][m];
      n++;
    }

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for fs files, there is a full NxN set of rhor arrays
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      type2rhor[i][j] = map[i] * fs->nelements + map[j];

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of fs elements

  nz2r = fs->nelements * (fs->nelements+1) / 2;
  memory->destroy(z2r);
  memory->create(z2r,nz2r,nr+1,"pair:z2r");

  // copy each element pair z2r to global z2r, only for I >= J

  n = 0;
  for (i = 0; i < fs->nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= nr; m++) z2r[n][m] = fs->z2r[i][j][m];
      n++;
    }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow,icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}
