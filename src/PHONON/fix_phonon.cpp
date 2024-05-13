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
   Contributing author: Ling-Ti Kong

   Contact:
     School of Materials Science and Engineering,
     Shanghai Jiao Tong University,
     800 Dongchuan Road, Minhang,
     Shanghai 200240, CHINA

     konglt@sjtu.edu.cn; konglt@gmail.com
------------------------------------------------------------------------- */

#include "fix_phonon.h"

#include "atom.h"
#include "citeme.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fft3d_wrap.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "tokenizer.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr int MAXLINE = 512;

enum{ FORWARD=-1, BACKWARD=1 };

static const char cite_fix_phonon[] =
  "fix phonon command: doi:10.1016/j.cpc.2011.04.019\n\n"
  "@Article{Kong11,\n"
  " author = {L. T. Kong},\n"
  " title = {Phonon Dispersion Measured Directly from Molecular Dynamics Simulations},\n"
  " journal = {Comput.\\ Phys.\\ Commun.},\n"
  " year =    2011,\n"
  " volume =  182,\n"
  " pages =   {2201--2207}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixPhonon::FixPhonon(LAMMPS *lmp,  int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_phonon);

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  if (narg < 8) error->all(FLERR,"Illegal fix phonon command: number of arguments < 8");

  nevery = utils::inumeric(FLERR, arg[3],false,lmp);   // Calculate this fix every n steps!
  if (nevery < 1) error->all(FLERR,"Illegal fix phonon command");

  nfreq  = utils::inumeric(FLERR, arg[4],false,lmp);   // frequency to output result
  if (nfreq < 1) error->all(FLERR,"Illegal fix phonon command");

  waitsteps = utils::bnumeric(FLERR,arg[5],false,lmp); // Wait this many timesteps before actually measuring
  if (waitsteps < 0) error->all(FLERR,"Illegal fix phonon command: waitsteps < 0 !");

  mapfile = utils::strdup(arg[6]);
  prefix = utils::strdup(arg[7]);
  logfile = utils::strdup(std::string(prefix)+".log");

  int sdim = sysdim = domain->dimension;
  int iarg = 8;
  nasr = 20;

  // other command line options
  while (iarg < narg) {
    if (strcmp(arg[iarg],"sysdim") == 0) {
      if (++iarg >= narg) error->all(FLERR,"Illegal fix phonon command: incomplete command line options.");
      sdim = utils::inumeric(FLERR, arg[iarg],false,lmp);
      if (sdim < 1) error->all(FLERR,"Illegal fix phonon command: sysdim should not be less than 1.");

    } else if (strcmp(arg[iarg],"nasr") == 0) {
      if (++iarg >= narg) error->all(FLERR,"Illegal fix phonon command: incomplete command line options.");
      nasr = utils::inumeric(FLERR, arg[iarg],false,lmp);

    } else {
      error->all(FLERR,"Illegal fix phonon command: unknown option read!");
    }

    ++iarg;
  }

  // get the dimension of the simulation; 1D is possible by specifying the option of "sysdim 1"
  if (sdim < sysdim) sysdim = sdim;
  nasr = MAX(0, nasr);

  // get the total number of atoms in group and run min/max checks
  bigint ng = group->count(igroup);
  if (ng > MAXSMALLINT) error->all(FLERR,"Too many atoms for fix phonon");
  if (ng < 1) error->all(FLERR,"No atom found for fix phonon!");
  ngroup = static_cast<int>(ng);


  // MPI gatherv related variables
  recvcnts = new int[nprocs];
  displs   = new int[nprocs];

  // mapping index
  tag2surf.clear(); // clear map info
  surf2tag.clear();

  // get the mapping between lattice indices and atom IDs

  atom->map_init();
  readmap();
  delete[] mapfile;
  if (nucell == 1) nasr = MIN(1,nasr);

  // get the mass matrix for dynamic matrix
  getmass();

  // create FFT and allocate memory for FFT
  // here the parallization is done on the x direction only
  nxlo = 0;
  int *nx_loc = new int [nprocs];
  for (int i = 0; i < nprocs; ++i) {
    nx_loc[i] = nx / nprocs;
    if (i < nx%nprocs) ++nx_loc[i];
  }
  for (int i = 0; i < me; ++i) nxlo += nx_loc[i];
  nxhi  = nxlo + nx_loc[me] - 1;
  mynpt = nx_loc[me] * ny * nz;
  mynq  = mynpt;

  fft_dim   = nucell  * sysdim;
  fft_dim2  = fft_dim * fft_dim;
  fft_nsend = mynpt   * fft_dim;

  fft_cnts  = new int[nprocs];
  fft_disp  = new int[nprocs];
  fft_disp[0] = 0;
  for (int i = 0; i < nprocs; ++i) fft_cnts[i] = nx_loc[i] * ny * nz * fft_dim;
  for (int i = 1; i < nprocs; ++i) fft_disp[i] = fft_disp[i-1] + fft_cnts[i-1];
  delete []nx_loc;

  fft = new FFT3d(lmp,world,nz,ny,nx,0,nz-1,0,ny-1,nxlo,nxhi,0,nz-1,0,ny-1,nxlo,nxhi,0,0,&mysize,0);
  memory->create(fft_data, MAX(1,mynq)*2, "fix_phonon:fft_data");

  // allocate variables; MAX(1,... is used because a null buffer will result in error for MPI
  memory->create(RIloc,ngroup,(sysdim+1),"fix_phonon:RIloc");
  memory->create(RIall,ngroup,(sysdim+1),"fix_phonon:RIall");
  memory->create(Rsort,ngroup, sysdim, "fix_phonon:Rsort");

  memory->create(Rnow, MAX(1,mynpt),fft_dim,"fix_phonon:Rnow");
  memory->create(Rsum, MAX(1,mynpt),fft_dim,"fix_phonon:Rsum");

  memory->create(basis,nucell, sysdim, "fix_phonon:basis");

  // because of hermit, only nearly half of q points are stored
  memory->create(Rqnow,MAX(1,mynq),fft_dim, "fix_phonon:Rqnow");
  memory->create(Rqsum,MAX(1,mynq),fft_dim2,"fix_phonon:Rqsum");
  memory->create(Phi_q,MAX(1,mynq),fft_dim2,"fix_phonon:Phi_q");

  // variable to collect all local Phi to root
  if (me == 0) memory->create(Phi_all,ntotal,fft_dim2,"fix_phonon:Phi_all");
  else memory->create(Phi_all,1,1,"fix_phonon:Phi_all");

  // output some information on the system to log file
  if (me == 0) {
    flog = fopen(logfile, "w");
    if (flog == nullptr)
      error->one(FLERR,"Can not open output file {}: {}", logfile,utils::getsyserror());
    fmt::print(flog,"############################################################\n");
    fmt::print(flog,"# group name of the atoms under study      : {}\n", group->names[igroup]);
    fmt::print(flog,"# total number of atoms in the group       : {}\n", ngroup);
    fmt::print(flog,"# dimension of the system                  : {} D\n", sysdim);
    fmt::print(flog,"# number of atoms per unit cell            : {}\n", nucell);
    fmt::print(flog,"# dimension of the FFT mesh                : {} x {} x {}\n", nx, ny, nz);
    fmt::print(flog,"# number of wait steps before measurement  : {}\n", waitsteps);
    fmt::print(flog,"# frequency of the measurement             : {}\n", nevery);
    fmt::print(flog,"# output result after this many measurement: {}\n", nfreq);
    fmt::print(flog,"# number of processors used by this run    : {}\n", nprocs);
    fmt::print(flog,"############################################################\n");
    fmt::print(flog,"# mapping information between lattice indices and atom id\n");
    fmt::print(flog,"# nx ny nz nucell\n");
    fmt::print(flog,"{} {} {} {}\n", nx, ny, nz, nucell);
    fmt::print(flog,"# l1 l2 l3 k atom_id\n");
    int ix, iy, iz, iu;
    for (idx = 0; idx < ngroup; ++idx) {
      itag = surf2tag[idx];
      iu   = idx%nucell;
      iz   = (idx/nucell)%nz;
      iy   = (idx/(nucell*nz))%ny;
      ix   = (idx/(nucell*nz*ny))%nx;
      fmt::print(flog,"{} {} {} {} {}\n", ix, iy, iz, iu, itag);
    }
    fmt::print(flog,"############################################################\n");
    fflush(flog);
  }
  surf2tag.clear();

  // default temperature is from thermo
  TempSum = new double[sysdim];
  id_temp = utils::strdup("thermo_temp");
  int icompute = modify->find_compute(id_temp);
  temperature = modify->compute[icompute];
  inv_nTemp = 1.0/group->count(temperature->igroup);

} // end of constructor

/* ---------------------------------------------------------------------- */

void FixPhonon::post_run()
{
  // compute and output final results
  if (ifreq > 0 && ifreq != nfreq) postprocess();
  if (me == 0) fclose(flog);
}

/* ---------------------------------------------------------------------- */

FixPhonon::~FixPhonon()
{
  // delete locally stored array
  memory->destroy(RIloc);
  memory->destroy(RIall);
  memory->destroy(Rsort);
  memory->destroy(Rnow);
  memory->destroy(Rsum);

  memory->destroy(basis);

  memory->destroy(Rqnow);
  memory->destroy(Rqsum);
  memory->destroy(Phi_q);
  memory->destroy(Phi_all);

  delete []recvcnts;
  delete []displs;
  delete []prefix;
  delete []logfile;
  delete []fft_cnts;
  delete []fft_disp;
  delete []id_temp;
  delete []TempSum;
  delete []M_inv_sqrt;
  delete []basetype;

  // destroy FFT
  delete fft;
  memory->sfree(fft_data);

  // clear map info
  tag2surf.clear();
  surf2tag.clear();

}

/* ---------------------------------------------------------------------- */

int FixPhonon::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;

return mask;
}

/* ---------------------------------------------------------------------- */

void FixPhonon::init()
{
  // warn if more than one fix-phonon
  int count = 0;
  for (int i = 0; i < modify->nfix; ++i) if (strcmp(modify->fix[i]->style,"phonon") == 0) ++count;
  if (count > 1 && me == 0) error->warning(FLERR,"More than one fix phonon defined"); // just warn, but allowed.
}

/* ---------------------------------------------------------------------- */

void FixPhonon::setup(int /*flag*/)
{
  // initialize accumulating variables
  for (int i = 0; i < sysdim; ++i) TempSum[i] = 0.;

  for (int i = 0; i < mynpt; ++i)
  for (int j = 0; j < fft_dim;  ++j) Rsum[i][j] = 0.;

  for (int i =0; i < mynq; ++i)
  for (int j =0; j < fft_dim2; ++j) Rqsum[i][j] = std::complex<double> (0.,0.);

  for (int i = 0; i < 6; ++i) hsum[i] = 0.;

  for (int i = 0; i < nucell; ++i)
  for (int j = 0; j < sysdim; ++j) basis[i][j] = 0.;

  neval = ifreq = 0;
  prev_nstep = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixPhonon::end_of_step()
{
  if ( (update->ntimestep-prev_nstep) <= waitsteps) return;

  double **x = atom->x;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double *h = domain->h;

  int i,idim,jdim,ndim;
  double xcur[3];

  // to get the current temperature
  if (!(temperature->invoked_flag & Compute::INVOKED_VECTOR)) temperature->compute_vector();
  for (idim = 0; idim < sysdim; ++idim) TempSum[idim] += temperature->vector[idim];

  // evaluate R(r) on local proc
  nfind = 0;
  for (i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      itag = tag[i];
      idx  = tag2surf[itag];

      domain->unmap(x[i], image[i], xcur);

      for (idim = 0; idim < sysdim; ++idim) RIloc[nfind][idim] = xcur[idim];
      RIloc[nfind++][sysdim] = static_cast<double>(idx);
    }
  }

  // gather R(r) on local proc, then sort and redistribute to all procs for FFT
  nfind *= (sysdim+1);
  displs[0] = 0;
  for (i = 0; i < nprocs; ++i) recvcnts[i] = 0;
  MPI_Gather(&nfind,1,MPI_INT,recvcnts,1,MPI_INT,0,world);
  for (i = 1; i < nprocs; ++i) displs[i] = displs[i-1] + recvcnts[i-1];

  MPI_Gatherv(RIloc[0],nfind,MPI_DOUBLE,RIall[0],recvcnts,displs,MPI_DOUBLE,0,world);
  if (me == 0) {
    for (i = 0; i < ngroup; ++i) {
      idx = static_cast<int>(RIall[i][sysdim]);
      for (idim = 0; idim < sysdim; ++idim) Rsort[idx][idim] = RIall[i][idim];
    }
  }
  MPI_Scatterv(Rsort[0],fft_cnts,fft_disp, MPI_DOUBLE, Rnow[0], fft_nsend, MPI_DOUBLE,0,world);

  // get Rsum
  for (idx = 0; idx < mynpt; ++idx)
  for (idim = 0; idim < fft_dim; ++idim) Rsum[idx][idim] += Rnow[idx][idim];

  // FFT R(r) to get R(q)
  for (idim = 0; idim < fft_dim; ++idim) {
    int m = 0;
    for (idx = 0; idx < mynpt; ++idx) {
      fft_data[m++] = static_cast<FFT_SCALAR>(Rnow[idx][idim]);
      fft_data[m++] = static_cast<FFT_SCALAR>(0.);
    }

    fft->compute(fft_data,fft_data,FORWARD);

    m = 0;
    for (idq = 0; idq < mynq; ++idq) {
      Rqnow[idq][idim] = std::complex<double>(static_cast<double>(fft_data[m]), static_cast<double>(fft_data[m+1]));
      m += 2;
    }
  }

  // to get sum(R(q).R(q)*)
  for (idq = 0; idq < mynq; ++idq) {
    ndim = 0;
    for (idim = 0; idim < fft_dim; ++idim)
    for (jdim = 0; jdim < fft_dim; ++jdim) Rqsum[idq][ndim++] += Rqnow[idq][idim] * std::conj(Rqnow[idq][jdim]);
  }

  // get basis info
  if (fft_dim > sysdim) {
    double dist2orig[3];
    for (idx = 0; idx < mynpt; ++idx) {
      ndim = sysdim;
      for (i = 1; i < nucell; ++i) {
        for (idim = 0; idim < sysdim; ++idim) dist2orig[idim] = Rnow[idx][ndim++] - Rnow[idx][idim];
        domain->minimum_image(dist2orig);
        for (idim = 0; idim < sysdim; ++idim) basis[i][idim] += dist2orig[idim];
      }
    }
  }
  // get lattice vector info
  for (int i = 0; i < 6; ++i) hsum[i] += h[i];

  // increment counter
  ++neval;

  // compute and output Phi_q after every nfreq evaluations
  if (++ifreq == nfreq) postprocess();

}   // end of end_of_step()

/* ---------------------------------------------------------------------- */

double FixPhonon::memory_usage()
{
  double bytes = (double)sizeof(double)*2*mynq
               + sizeof(std::map<int,int>)*2*ngroup
               + sizeof(double)*(ngroup*(3*sysdim+2)+mynpt*fft_dim*2)
               + sizeof(std::complex<double>)*MAX(1,mynq)*fft_dim *(1+2*fft_dim)
               + sizeof(std::complex<double>)*ntotal*fft_dim2
               + sizeof(int) * nprocs * 4;
  return bytes;
}

/* ---------------------------------------------------------------------- */

int FixPhonon::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    delete [] id_temp;
    id_temp = utils::strdup(arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify temp ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,"Fix_modify temp ID does not compute temperature");
    inv_nTemp = 1.0/group->count(temperature->igroup);

    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
 * private method, to get the mass matrix for dynamic matrix
 * --------------------------------------------------------------------*/
void FixPhonon::getmass()
{
  int nlocal = atom->nlocal;
  int *mask  = atom->mask;
  tagint *tag   = atom->tag;
  int *type  = atom->type;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double *mass_one, *mass_all;
  double *type_one, *type_all;

  mass_one = new double[nucell];
  mass_all = new double[nucell];
  type_one = new double[nucell];
  type_all = new double[nucell];
  for (int i = 0; i < nucell; ++i)  mass_one[i] = type_one[i] = 0.;

  if (rmass) {
    for (int i = 0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        itag = tag[i];
        idx  = tag2surf[itag];
        int iu = idx%nucell;
        mass_one[iu] += rmass[i];
        type_one[iu] += double(type[i]);
      }
    }
  } else {
    for (int i = 0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        itag = tag[i];
        idx  = tag2surf[itag];
        int iu = idx%nucell;
        mass_one[iu] += mass[type[i]];
        type_one[iu] += double(type[i]);
      }
    }
  }

  MPI_Allreduce(mass_one,mass_all,nucell,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(type_one,type_all,nucell,MPI_DOUBLE,MPI_SUM,world);

  M_inv_sqrt = new double[nucell];
  basetype   = new int[nucell];

  double inv_total = 1./double(ntotal);
  for (int i = 0; i < nucell; ++i) {
    mass_all[i] *= inv_total;
    M_inv_sqrt[i] = sqrt(1./mass_all[i]);

    basetype[i] = int(type_all[i]*inv_total);
  }
  delete []mass_one;
  delete []mass_all;
  delete []type_one;
  delete []type_all;
}


/* ----------------------------------------------------------------------
 * private method, to read the mapping info from file
 * --------------------------------------------------------------------*/

void FixPhonon::readmap()
{
  int info = 0;

  // auto-generate mapfile for "cluster" (gamma only system)
  if (strcmp(mapfile, "GAMMA") == 0) {
    nx = ny = nz = ntotal = 1;
    nucell = ngroup;

    tagint *tag_loc, *tag_all;
    memory->create(tag_loc,ngroup,"fix_phonon:tag_loc");
    memory->create(tag_all,ngroup,"fix_phonon:tag_all");

    // get atom IDs on local proc
    int nfind = 0;
    for (int i = 0; i < atom->nlocal; ++i) {
      if (atom->mask[i] & groupbit) tag_loc[nfind++] = atom->tag[i];
    }

    // gather IDs on local proc
    displs[0] = 0;
    for (int i = 0; i < nprocs; ++i) recvcnts[i] = 0;
    MPI_Allgather(&nfind,1,MPI_INT,recvcnts,1,MPI_INT,world);
    for (int i = 1; i < nprocs; ++i) displs[i] = displs[i-1] + recvcnts[i-1];

    MPI_Allgatherv(tag_loc,nfind,MPI_LMP_TAGINT,tag_all,recvcnts,displs,MPI_LMP_TAGINT,world);
    for (int i = 0; i < ngroup; ++i) {
      itag = tag_all[i];
      tag2surf[itag] = i;
      surf2tag[i] = itag;
    }

    memory->destroy(tag_loc);
    memory->destroy(tag_all);
    return;
  }

  // read from map file for others
  char line[MAXLINE] = {'\0'};
  FILE *fp = fopen(mapfile, "r");
  if (fp == nullptr)
    error->all(FLERR,"Cannot open input map file {}: {}", mapfile, utils::getsyserror());

  if (fgets(line,MAXLINE,fp) == nullptr)
    error->all(FLERR,"Error while reading header of mapping file!");
  try {
    ValueTokenizer values(line);

    nx = values.next_int();
    ny = values.next_int();
    nz = values.next_int();
    nucell = values.next_int();
  } catch (TokenizerException &e) {
    error->all(FLERR, "Incorrect header format: {}", e.what());
  }

  ntotal = nx*ny*nz;
  if (ntotal*nucell != ngroup)
    error->all(FLERR,"FFT mesh and number of atoms in group mismatch!");

  // second line of mapfile is comment
  if (fgets(line,MAXLINE,fp) == nullptr)
    error->all(FLERR,"Error while reading comment of mapping file!");

  try {
    int ix, iy, iz, iu;
    // the remaining lines carry the mapping info
    for (int i = 0; i < ngroup; ++i) {
      if (fgets(line,MAXLINE,fp) == nullptr) {info = 1; break;}
      ValueTokenizer values(line);
      ix   = values.next_int();
      iy   = values.next_int();
      iz   = values.next_int();
      iu   = values.next_int();
      itag   = values.next_tagint();

      // check if index is in correct range
      if (ix < 0 || ix >= nx || iy < 0 || iy >= ny ||
          iz < 0 || iz >= nz || iu < 0 || iu >= nucell) {info = 2; break;}
      // 1 <= itag <= natoms
      if (itag < 1 || itag > atom->map_tag_max) {info = 3; break;}
      idx = ((ix*ny+iy)*nz+iz)*nucell + iu;
      tag2surf[itag] = idx;
      surf2tag[idx]  = itag;
    }
  } catch (TokenizerException &e) {
    error->all(FLERR, "Incorrect map file format: {}", e.what());
  }
  fclose(fp);

  if (tag2surf.size() != surf2tag.size() ||
      tag2surf.size() != static_cast<std::size_t>(ngroup) )
    error->all(FLERR,"The mapping is incomplete!");
  if (info) error->all(FLERR,"Error while reading mapping file!");

  // check the correctness of mapping
  int *mask  = atom->mask;
  tagint *tag   = atom->tag;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      itag = tag[i];
      idx  = tag2surf[itag];
      if (itag != surf2tag[idx])
        error->one(FLERR,"The mapping info read is incorrect!");
    }
  }
}

/* ----------------------------------------------------------------------
 * private method, to output the force constant matrix
 * --------------------------------------------------------------------*/
void FixPhonon::postprocess( )
{
  if (neval < 1) return;

  ifreq = 0;
  int idim, jdim, ndim;
  double inv_neval = 1. /double(neval);

  // to get <Rq.Rq*>
  for (idq = 0; idq < mynq; ++idq)
  for (idim = 0; idim < fft_dim2; ++idim) Phi_q[idq][idim] = Rqsum[idq][idim] * inv_neval;

  // to get <R>
  for (idx = 0; idx < mynpt; ++idx)
  for (idim = 0; idim < fft_dim; ++idim) Rnow[idx][idim] = Rsum[idx][idim] * inv_neval;

  // to get <R>q
  for (idim = 0; idim < fft_dim; ++idim) {
    int m = 0;
    for (idx = 0; idx < mynpt; ++idx) {
      fft_data[m++] = static_cast<FFT_SCALAR>(Rnow[idx][idim]);
      fft_data[m++] = static_cast<FFT_SCALAR>(0.);
    }

    fft->compute(fft_data,fft_data,FORWARD);

    m = 0;
    for (idq = 0; idq < mynq; ++idq) {
      Rqnow[idq][idim]  = std::complex<double>(static_cast<double>(fft_data[m]), static_cast<double>(fft_data[m+1]));
      m += 2;
    }
  }

  // to get G(q) = <Rq.Rq*> - <R>q.<R*>q
  for (idq = 0; idq < mynq; ++idq) {
    ndim = 0;
    for (idim = 0; idim < fft_dim; ++idim)
    for (jdim = 0; jdim < fft_dim; ++jdim) Phi_q[idq][ndim++] -= Rqnow[idq][idim] * std::conj(Rqnow[idq][jdim]);
  }

  // to get Phi = KT.G^-1; normalization of FFTW data is done here
  double boltz = force->boltz, TempAve = 0.;
  auto kbtsqrt = new double[sysdim];
  double TempFac = inv_neval * inv_nTemp;
  double NormFac = TempFac * double(ntotal);

  for (idim = 0; idim < sysdim; ++idim) {
    kbtsqrt[idim] = sqrt(TempSum[idim] * NormFac);
    TempAve += TempSum[idim] * TempFac;
  }
  TempAve /= sysdim*boltz;

  for (idq = 0; idq < mynq; ++idq) {
    GaussJordan(fft_dim, Phi_q[idq]);
    ndim =0;
    for (idim = 0; idim < fft_dim; ++idim)
    for (jdim = 0; jdim < fft_dim; ++jdim) Phi_q[idq][ndim++] *= kbtsqrt[idim%sysdim]*kbtsqrt[jdim%sysdim];
  }

  // to collect all local Phi_q to root
  displs[0]=0;
  for (int i = 0; i < nprocs; ++i) recvcnts[i] = fft_cnts[i]*fft_dim*2;
  for (int i = 1; i < nprocs; ++i) displs[i] = displs[i-1] + recvcnts[i-1];
  MPI_Gatherv(Phi_q[0],mynq*fft_dim2*2,MPI_DOUBLE,Phi_all[0],recvcnts,displs,MPI_DOUBLE,0,world);

  // to collect all basis info and averaged it on root
  auto basis_root = new double[fft_dim];
  if (fft_dim > sysdim) MPI_Reduce(&basis[1][0], &basis_root[sysdim], fft_dim-sysdim, MPI_DOUBLE, MPI_SUM, 0, world);

  if (me == 0) { // output dynamic matrix by root

    // get basis info
    for (idim = 0;      idim < sysdim;  ++idim) basis_root[idim]  = 0.;
    for (idim = sysdim; idim < fft_dim; ++idim) basis_root[idim] /= double(ntotal)*double(neval);
    // get unit cell base vector info; might be incorrect if MD pbc and FixPhonon pbc mismatch.
    double basevec[9];
    basevec[1] = basevec[2] = basevec[5] = 0.;
    basevec[0] = hsum[0] * inv_neval / double(nx);
    basevec[4] = hsum[1] * inv_neval / double(ny);
    basevec[8] = hsum[2] * inv_neval / double(nz);
    basevec[7] = hsum[3] * inv_neval / double(nz);
    basevec[6] = hsum[4] * inv_neval / double(nz);
    basevec[3] = hsum[5] * inv_neval / double(ny);

    // write binary file, in fact, it is the force constants matrix that is written
    // Enforcement of ASR and the conversion of dynamical matrix is done in the postprocessing code
    auto fname = fmt::format("{}.bin.{}",prefix,update->ntimestep);
    FILE *fp_bin = fopen(fname.c_str(),"wb");

    fwrite(&sysdim, sizeof(int),    1, fp_bin);
    fwrite(&nx,     sizeof(int),    1, fp_bin);
    fwrite(&ny,     sizeof(int),    1, fp_bin);
    fwrite(&nz,     sizeof(int),    1, fp_bin);
    fwrite(&nucell, sizeof(int),    1, fp_bin);
    fwrite(&boltz,  sizeof(double), 1, fp_bin);

    fwrite(Phi_all[0],sizeof(double),(bigint)ntotal*fft_dim2*2,fp_bin);

    fwrite(&TempAve,      sizeof(double),1,      fp_bin);
    fwrite(&basevec[0],   sizeof(double),9,      fp_bin);
    fwrite(&basis_root[0],sizeof(double),fft_dim,fp_bin);
    fwrite(basetype,      sizeof(int),   nucell, fp_bin);
    fwrite(M_inv_sqrt,    sizeof(double),nucell, fp_bin);

    fclose(fp_bin);

    // write log file, here however, it is the dynamical matrix that is written
    fmt::print(flog,"############################################################\n");
    fmt::print(flog,"# Current time step                      : {}\n", update->ntimestep);
    fmt::print(flog,"# Total number of measurements           : {}\n", neval);
    fmt::print(flog,"# Average temperature of the measurement : {}\n", TempAve);
    fmt::print(flog,"# Boltzmann constant under current units : {}\n", boltz);
    fmt::print(flog,"# basis vector A1 = [{} {} {}]\n", basevec[0], basevec[1], basevec[2]);
    fmt::print(flog,"# basis vector A2 = [{} {} {}]\n", basevec[3], basevec[4], basevec[5]);
    fmt::print(flog,"# basis vector A3 = [{} {} {}]\n", basevec[6], basevec[7], basevec[8]);
    fmt::print(flog,"############################################################\n");
    fmt::print(flog,"# qx\t qy \t qz \t\t Phi(q)\n");

    EnforceASR();

    // to get D = 1/M x Phi
    for (idq = 0; idq < ntotal; ++idq) {
      ndim =0;
      for (idim = 0; idim < fft_dim; ++idim)
        for (jdim = 0; jdim < fft_dim; ++jdim)
          Phi_all[idq][ndim++] *= M_inv_sqrt[idim/sysdim]*M_inv_sqrt[jdim/sysdim];
    }

    idq =0;
    for (int ix = 0; ix < nx; ++ix) {
      double qx = double(ix)/double(nx);
      for (int iy = 0; iy < ny; ++iy) {
        double qy = double(iy)/double(ny);
        for (int iz = 0; iz < nz; ++iz) {
          double qz = double(iz)/double(nz);
          fmt::print(flog,"{} {} {}", qx, qy, qz);
          for (idim = 0; idim < fft_dim2; ++idim)
            fmt::print(flog, " {} {}", std::real(Phi_all[idq][idim]), std::imag(Phi_all[idq][idim]));
          fmt::print(flog, "\n");
          ++idq;
        }
      }
    }
    fflush(flog);
  }
  delete[] kbtsqrt;
  delete[] basis_root;
}   // end of postprocess

/* ----------------------------------------------------------------------
 * private method, to get the inverse of a complex matrix by means of
 * Gaussian-Jordan Elimination with full pivoting; square matrix required.
 *
 * Adapted from the Numerical Recipes in Fortran.
 * --------------------------------------------------------------------*/
void FixPhonon::GaussJordan(int n, std::complex<double> *Mat)
{
  int i,icol,irow,j,k,l,ll,idr,idc;
  int *indxc,*indxr,*ipiv;
  double big, nmjk;
  std::complex<double> dum, pivinv;

  indxc = new int[n];
  indxr = new int[n];
  ipiv  = new int[n];

  for (i = 0; i < n; ++i) ipiv[i] = 0;
  for (i = 0; i < n; ++i) {
    big = 0.;
    for (j = 0; j < n; ++j) {
      if (ipiv[j] != 1) {
        for (k = 0; k < n; ++k) {
          if (ipiv[k] == 0) {
            idr = j*n+k;
            nmjk = norm(Mat[idr]);
            if (nmjk >= big) {
              big  = nmjk;
              irow = j;
              icol = k;
            }
          } else if (ipiv[k] > 1) error->one(FLERR,"Singular matrix in complex GaussJordan!");
        }
      }
    }
    ipiv[icol] += 1;
    if (irow != icol) {
      for (l = 0; l < n; ++l) {
        idr  = irow*n+l;
        idc  = icol*n+l;
        dum  = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    idr = icol*n+icol;
    if (Mat[idr] == std::complex<double>(0.,0.)) error->one(FLERR,"Singular matrix in complex GaussJordan!");

    pivinv = 1./ Mat[idr];
    Mat[idr] = std::complex<double>(1.,0.);
    idr = icol*n;
    for (l = 0; l < n; ++l) Mat[idr+l] *= pivinv;
    for (ll = 0; ll < n; ++ll) {
      if (ll != icol) {
        idc = ll*n + icol;
        dum = Mat[idc];
        Mat[idc] = 0.;
        idc -= icol;
        for (l = 0; l < n; ++l) Mat[idc+l] -= Mat[idr+l]*dum;
      }
    }
  }

  for (l = n-1; l >= 0; --l) {
    int rl = indxr[l];
    int cl = indxc[l];
    if (rl != cl) {
      for (k = 0; k < n; ++k) {
        idr = k*n + rl;
        idc = k*n + cl;
        dum = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
  }
  delete []indxr;
  delete []indxc;
  delete []ipiv;
}

/* ----------------------------------------------------------------------
 * private method, to apply the acoustic sum rule on force constant matrix
 * at gamma point. Should be executed on root only.
 * --------------------------------------------------------------------*/
void FixPhonon::EnforceASR()
{
  if (nasr < 1) return;

  for (int iit = 0; iit < nasr; ++iit) {
    // simple ASR; the resultant matrix might not be symmetric
    for (int a = 0; a < sysdim; ++a)
    for (int b = 0; b < sysdim; ++b) {
      for (int k = 0; k < nucell; ++k) {
        double sum = 0.;
        for (int kp = 0; kp < nucell; ++kp) {
          int idx = (k*sysdim+a)*fft_dim + kp*sysdim + b;
          sum += std::real(Phi_all[0][idx]);
        }
        sum /= double(nucell);
        for (int kp = 0; kp < nucell; ++kp) {
          int idx = (k*sysdim+a)*fft_dim + kp*sysdim + b;
          Phi_all[0][idx] -= sum;
        }
      }
    }

    // symmetrize
    for (int k = 0; k < nucell; ++k)
    for (int kp = k; kp < nucell; ++kp) {
      double csum = 0.;
      for (int a = 0; a < sysdim; ++a)
      for (int b = 0; b < sysdim; ++b) {
        int idx = (k*sysdim+a)*fft_dim + kp*sysdim + b;
        int jdx = (kp*sysdim+b)*fft_dim + k*sysdim + a;
        csum = (std::real(Phi_all[0][idx])+std::real(Phi_all[0][jdx]))*0.5;
        Phi_all[0][idx] = std::complex<double>(csum, std::imag(Phi_all[0][idx]));
        Phi_all[0][jdx] = std::complex<double>(csum, std::imag(Phi_all[0][jdx]));
      }
    }
  }

  // symmetric ASR
  for (int a = 0; a < sysdim; ++a)
  for (int b = 0; b < sysdim; ++b) {
    for (int k = 0; k < nucell; ++k) {
      double sum = 0.;
      for (int kp = 0; kp < nucell; ++kp) {
        int idx = (k*sysdim+a)*fft_dim + kp*sysdim + b;
        sum += std::real(Phi_all[0][idx]);
      }
      sum /= double(nucell-k);
      for (int kp = k; kp < nucell; ++kp) {
        int idx = (k*sysdim+a)*fft_dim + kp*sysdim + b;
        int jdx = (kp*sysdim+b)*fft_dim + k*sysdim + a;
        Phi_all[0][idx] -= sum;
        Phi_all[0][jdx] = std::complex<double>(std::real(Phi_all[0][idx]),
                                               std::imag(Phi_all[0][jdx]));
      }
    }
  }
}
/* --------------------------------------------------------------------*/
