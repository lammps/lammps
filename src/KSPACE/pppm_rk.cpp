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
 * For the original PPPM class
   Contributing authors: Roy Pollock (LLNL), Paul Crozier (SNL)
     per-atom energy/virial & group/group energy/force added by Stan Moore (BYU)
     analytic diff (2 FFT) option added by Rolf Isele-Holder (Aachen University)
     triclinic added by Stan Moore (SNL)
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
 * For this new derivative class PPPM_RK of PPPM:
   Author: Brian Dandurand (Queen's University Belfast)
   This class PPPM_RK is needed to implement the enhanced baseline of:
      Brian Dandurand, Hans Vandierendonck, and Bronis de Supinski.
      "Improving Parallel Scalability for Molecular Dynamics Simulations in the Exascale Era".
      in Proceedings of the IPDPS Conference. 2025.
   The enhanced baseline in turn was inspired by the earlier contribution of
      D. F. Richards, J. N. Glosli, B. Chan, M. R. Dorr, E. W. Draeger, J.-
      L. Fattebert, W. D. Krauss, T. Spelce, F. H. Streitz, M. P. Surh, and
      J. A. Gunnels,
      "Beyond homogeneous decomposition: scaling long-range forces
      on massively parallel systems,"
      in Proceedings of the Conference on High Performance Computing Networking,
      Storage and Analysis, ser. SC '09. New York, NY, USA:
      Association for Computing Machinery, 2009.
------------------------------------------------------------------------- */

#include "pppm_rk.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "grid3d.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "universe.h"

using namespace LAMMPS_NS;
using namespace MathConst;

static const char cite_pppm_rk[] =
    "kspace_style pppm/rk command: "
    "https://doi.org/10.1109/IPDPS64566.2025.00077\n\n"
    "@inproceedings{11078563,\n"
    "author={Dandurand, Brian and Vandierendonck, Hans and De Supinski, Bronis R.},\n"
    "booktitle={2025 IEEE International Parallel and Distributed Processing Symposium (IPDPS)},\n"
    "title={Improving Parallel Scalability for Molecular Dynamics Simulations in the Exascale Era},\n"
    "year={2025},\n"
    "volume={},\n"
    "number={},\n"
    "pages={813-823},\n"
    "keywords={Distributed processing;Accuracy;Scalability;Computational modeling;Parallel processing;"
    "Trajectory;Timing;Maintenance;History;Optimization;molecular dynamics simulation;parallelism;"
    "asynchrony;scalability;LAMMPS;MPI},\n"
    "doi={10.1109/IPDPS64566.2025.00077}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

PPPM_RK::PPPM_RK(LAMMPS *lmp) :
    PPPM(lmp), density_brick_buf(nullptr), vdx_brick_buf(nullptr), vdy_brick_buf(nullptr),
    vdz_brick_buf(nullptr), u_brick_buf(nullptr), density_sizes(nullptr), density_disps(nullptr),
    partitionInfoK(nullptr), block_size(0)
{
  rk_flag = 1;
  group_group_enable = 0;
  if (lmp->citeme) {
      lmp->citeme->add(cite_pppm_rk);
  }
}

/* ---------------------------------------------------------------------- */

void PPPM_RK::init()
{
  PPPM::init();
  if (me == 0) utils::logmesg(lmp,"PPPM_RK initialization ...\n");
  memory->create3d_offset(density_brick_buf,
                          nzlo_in,nzhi_in,nylo_in,nyhi_in,
                          nxlo_in,nxhi_in,"pppm_split:density_brick_buf");

  if (differentiation_flag == 1) {
    memory->create3d_offset(u_brick_buf,
                            nzlo_in,nzhi_in,nylo_in,nyhi_in,
                            nxlo_in,nxhi_in,"pppm_split:u_brick_buf");
  } else {
    memory->create3d_offset(vdx_brick_buf,
                            nzlo_in,nzhi_in,nylo_in,nyhi_in,
                            nxlo_in,nxhi_in,"pppm_split:vdx_brick_buf");
    memory->create3d_offset(vdy_brick_buf,
                            nzlo_in,nzhi_in,nylo_in,nyhi_in,
                            nxlo_in,nxhi_in,"pppm_split:vdy_brick_buf");
    memory->create3d_offset(vdz_brick_buf,
                            nzlo_in,nzhi_in,nylo_in,nyhi_in,
                            nxlo_in,nxhi_in,"pppm_split:vdz_brick_buf");
  }
  setupRKBlock();
}

/* ---------------------------------------------------------------------- */

//Adapted from the constructor of VerletSplit
/*Main purpose is to instantiate the inter-RK block communicator
  and to assign values to ratio, rproc, me_block, and block_size*/
void PPPM_RK::setupRKBlock()
{
  // error checks on partitions

  if (universe->nworlds != 2)
    error->universe_all(FLERR,"pppm/rk requires 2 partitions");
  if (universe->procs_per_world[0] % universe->procs_per_world[1])
    error->universe_all(FLERR,"pppm/rk requires Rspace partition "
                        "size be multiple of Kspace partition size");
  if (comm->style != Comm::BRICK)
    error->universe_all(FLERR,"pppm/rk can only currently be used with comm_style brick");

  // rproc = 1 for Rspace procs, 0 for Kspace procs

  if (universe->iworld == 0) rproc = 1;
  else rproc = 0;

  ratio = universe->procs_per_world[0] / universe->procs_per_world[1];

  // Kspace root proc broadcasts info about Kspace proc layout to Rspace procs

  int kspace_procgrid[3];

  if (universe->me == universe->root_proc[1]) {
    kspace_procgrid[0] = comm->procgrid[0];
    kspace_procgrid[1] = comm->procgrid[1];
    kspace_procgrid[2] = comm->procgrid[2];
  }
  MPI_Bcast(kspace_procgrid,3,MPI_INT,universe->root_proc[1],universe->uworld);

  int ***kspace_grid2proc;
  memory->create(kspace_grid2proc,kspace_procgrid[0],
                 kspace_procgrid[1],kspace_procgrid[2],
                 "pppm/rk:kspace_grid2proc");

  if (universe->me == universe->root_proc[1]) {
    for (int i = 0; i < comm->procgrid[0]; i++)
      for (int j = 0; j < comm->procgrid[1]; j++)
        for (int k = 0; k < comm->procgrid[2]; k++)
          kspace_grid2proc[i][j][k] = comm->grid2proc[i][j][k];
  }
  MPI_Bcast(&kspace_grid2proc[0][0][0],
            kspace_procgrid[0]*kspace_procgrid[1]*kspace_procgrid[2],MPI_INT,
            universe->root_proc[1],universe->uworld);

  // Rspace partition must be multiple of Kspace partition in each dim
  // so atoms of one Kspace proc coincide with atoms of several Rspace procs

  if (rproc) {
    int flag = 0;
    if (comm->procgrid[0] % kspace_procgrid[0]) flag = 1;
    if (comm->procgrid[1] % kspace_procgrid[1]) flag = 1;
    if (comm->procgrid[2] % kspace_procgrid[2]) flag = 1;
    if (flag)
      error->one(FLERR,
                 "pppm/rk requires Rspace partition layout be "
                 "multiple of Kspace partition layout in each dim");
  }

  // block = 1 Kspace proc with set of Rspace procs it overlays
  // me_block = 0 for Kspace proc
  // me_block = 1 to ratio for Rspace procs
  // block = MPI communicator for that set of procs

  int iblock,key;

  if (!rproc) {
    iblock = comm->me;
    key = 0;
  } else {
    int kpx = comm->myloc[0] / (comm->procgrid[0]/kspace_procgrid[0]);
    int kpy = comm->myloc[1] / (comm->procgrid[1]/kspace_procgrid[1]);
    int kpz = comm->myloc[2] / (comm->procgrid[2]/kspace_procgrid[2]);
    iblock = kspace_grid2proc[kpx][kpy][kpz];
    key = 1;
  }

  MPI_Comm_split(universe->uworld,iblock,key,&block);
  MPI_Comm_rank(block, &me_block);
  MPI_Comm_size(block, &block_size);

  // output block groupings to universe screen/logfile
  // bmap is ordered by block and then by proc within block

  int *bmap = new int[universe->nprocs];
  for (int i = 0; i < universe->nprocs; i++) bmap[i] = -1;
  bmap[iblock*(ratio+1)+me_block] = universe->me;

  int *bmapall = new int[universe->nprocs];
  MPI_Allreduce(bmap,bmapall,universe->nprocs,MPI_INT,MPI_MAX,universe->uworld);

  if (universe->me == 0) {
    if (universe->uscreen) {
      fprintf(universe->uscreen,
              "Per-block Rspace/Kspace proc IDs (original proc IDs):\n");
      int m = 0;
      for (int i = 0; i < universe->nprocs/(ratio+1); i++) {
        fprintf(universe->uscreen,"  block %d:",i);
        int kspace_proc = bmapall[m];
        for (int j = 1; j <= ratio; j++)
          fprintf(universe->uscreen," %d",bmapall[m+j]);
        fprintf(universe->uscreen," %d",kspace_proc);
        kspace_proc = bmapall[m];
        for (int j = 1; j <= ratio; j++) {
          if (j == 1) fprintf(universe->uscreen," (");
          else fprintf(universe->uscreen," ");
          fprintf(universe->uscreen,"%d",
                  universe->uni2orig[bmapall[m+j]]);
        }
        fprintf(universe->uscreen," %d)\n",universe->uni2orig[kspace_proc]);
        m += ratio + 1;
      }
    }
    if (universe->ulogfile) {
      fprintf(universe->ulogfile,
              "Per-block Rspace/Kspace proc IDs (original proc IDs):\n");
      int m = 0;
      for (int i = 0; i < universe->nprocs/(ratio+1); i++) {
        fprintf(universe->ulogfile,"  block %d:",i);
        int kspace_proc = bmapall[m];
        for (int j = 1; j <= ratio; j++)
          fprintf(universe->ulogfile," %d",bmapall[m+j]);

        fprintf(universe->ulogfile," %d",kspace_proc);
        kspace_proc = bmapall[m];
        for (int j = 1; j <= ratio; j++) {
          if (j == 1) fprintf(universe->ulogfile," (");
          else fprintf(universe->ulogfile," ");
          fprintf(universe->ulogfile,"%d",
                  universe->uni2orig[bmapall[m+j]]);
        }
        fprintf(universe->ulogfile," %d)\n",universe->uni2orig[kspace_proc]);
        m += ratio + 1;
      }
    }
  }

  memory->destroy(kspace_grid2proc);
  delete[] bmap;
  delete[] bmapall;

  /*Copy communicators, assign structures helpful for inter-RK communication*/
  MPI_Comm_dup(block, &block_density_brick);
  MPI_Comm_dup(block, &block_kforce_x);
  MPI_Comm_dup(block, &block_kforce_y);
  MPI_Comm_dup(block, &block_kforce_z);
  MPI_Comm_dup(block, &block_kforce_u);
  MPI_Comm_dup(block, &block_energy);
  MPI_Comm_dup(block, &block_virial);

  density_sizes = new int[block_size];
  density_disps = new int[block_size];

  memory->create(partitionInfoK,block_size,7,"pppm:partitionInfoK");

  int partitionInfo[7];
  partitionInfo[0] = (nzhi_in - nzlo_in + 1)*(nyhi_in - nylo_in + 1)*(nxhi_in - nxlo_in + 1);
  partitionInfo[1] = nzlo_in;
  partitionInfo[2] = nzhi_in;
  partitionInfo[3] = nylo_in;
  partitionInfo[4] = nyhi_in;
  partitionInfo[5] = nxlo_in;
  partitionInfo[6] = nxhi_in;
  MPI_Gather(partitionInfo,7,MPI_INT,&partitionInfoK[0][0],7,MPI_INT,0,block);
  density_sizes[0] = 0;
  density_disps[0] = 0;
  for (int i = 1; i < block_size; i++) {
    density_sizes[i] = partitionInfoK[i][0];
    density_disps[i] = density_disps[i-1]+density_sizes[i-1];
  }
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

PPPM_RK::~PPPM_RK()
{
  if (copymode) return;

  memory->destroy3d_offset(density_brick_buf,nzlo_in,nylo_in,nxlo_in);
  if (differentiation_flag == 1) {
    memory->destroy3d_offset(u_brick_buf,nzlo_in,nylo_in,nxlo_in);
  } else {
    memory->destroy3d_offset(vdx_brick_buf,nzlo_in,nylo_in,nxlo_in);
    memory->destroy3d_offset(vdy_brick_buf,nzlo_in,nylo_in,nxlo_in);
    memory->destroy3d_offset(vdz_brick_buf,nzlo_in,nylo_in,nxlo_in);
  }
  MPI_Comm_free(&block);

  MPI_Comm_free(&block_density_brick);
  MPI_Comm_free(&block_kforce_x);
  MPI_Comm_free(&block_kforce_y);
  MPI_Comm_free(&block_kforce_z);
  MPI_Comm_free(&block_kforce_u);
  MPI_Comm_free(&block_energy);
  MPI_Comm_free(&block_virial);

  delete[] density_sizes;
  delete[] density_disps;

  memory->destroy(partitionInfoK);
}


/* ----------------------------------------------------------------------
   communicate charge densities to Kspace processes
   also eflag,vflag and box bounds if needed
------------------------------------------------------------------------- */
void PPPM_RK::r2k_comm(int &eflag, int &vflag)
{
  if (rproc) {
    compute_charge_densities(eflag,vflag);
    //Post sending of density brick gathering
    for (int k = nzlo_in; k <= nzhi_in; k++)
      for (int j = nylo_in; j <= nyhi_in; j++)
        memcpy(&density_brick_buf[k][j][nxlo_in],&density_brick[k][j][nxlo_in],(nxhi_in-nxlo_in+1)*sizeof(FFT_SCALAR));

    MPI_Igatherv(&density_brick_buf[nzlo_in][nylo_in][nxlo_in],nfft_brick,MPI_FFT_SCALAR,
                 nullptr,nullptr,nullptr,MPI_FFT_SCALAR,0,block_density_brick,&mpi_requests_density);

    // Rspace sends eflag,vflag (and possibly domain box_change) to Kspace
    if (me_block == 1) {
      flags_buf[0] = eflag; flags_buf[1] = vflag;
      MPI_Isend(flags_buf,2,MPI_INT,0,TAG_FLAGS,block,&mpi_requests_flags);
      // send box bounds from Rspace to Kspace if simulation box is dynamic
      if (domain->box_change) {
        memcpy(domain_boxbuf,domain->boxlo,3*sizeof(double));
        memcpy(domain_boxbuf+3,domain->boxhi,3*sizeof(double));
        MPI_Isend(domain_boxbuf,6,MPI_DOUBLE,0,TAG_BOX,block,&mpi_requests_domain_box);
      }
    }

    //Posting R-process recv communications
    if (eflag) MPI_Ibcast(&energy_buf,1,MPI_DOUBLE,0,block_energy,&mpi_requests_energy);
    if (vflag) MPI_Ibcast(virial_mpi_buf,6,MPI_DOUBLE,0,block_virial,&mpi_requests_virial);

    int buf_len = (nzhi_in - nzlo_in + 1)*(nyhi_in - nylo_in + 1)*(nxhi_in - nxlo_in + 1);
    if (differentiation_flag == 1) {
      FFT_SCALAR *usrc = &u_brick_buf[nzlo_in][nylo_in][nxlo_in];
      MPI_Iscatterv(nullptr,nullptr,nullptr,MPI_FFT_SCALAR,
                    usrc,buf_len,MPI_FFT_SCALAR,0,block_kforce_u,&mpi_requests_grid_u);
    } else {
      FFT_SCALAR *xsrc = &vdx_brick_buf[nzlo_in][nylo_in][nxlo_in];
      FFT_SCALAR *ysrc = &vdy_brick_buf[nzlo_in][nylo_in][nxlo_in];
      FFT_SCALAR *zsrc = &vdz_brick_buf[nzlo_in][nylo_in][nxlo_in];
      MPI_Iscatterv(nullptr,nullptr,nullptr,MPI_FFT_SCALAR,
                    xsrc,buf_len,MPI_FFT_SCALAR,0,block_kforce_x,&mpi_requests_grid_x);
      MPI_Iscatterv(nullptr,nullptr,nullptr,MPI_FFT_SCALAR,
                    ysrc,buf_len,MPI_FFT_SCALAR,0,block_kforce_y,&mpi_requests_grid_y);
      MPI_Iscatterv(nullptr,nullptr,nullptr,MPI_FFT_SCALAR,
                    zsrc,buf_len,MPI_FFT_SCALAR,0,block_kforce_z,&mpi_requests_grid_z);
    }
  } else { //kproc
    // Kspace receives eflag,vflag from Rspace
    MPI_Recv(flags_buf,2,MPI_INT,1,TAG_FLAGS,block,MPI_STATUS_IGNORE);
    eflag = flags_buf[0]; vflag = flags_buf[1];
    // send box bounds from Rspace to Kspace if simulation box is dynamic
    if (domain->box_change) {
      MPI_Recv(domain_boxbuf,6,MPI_DOUBLE,1,TAG_BOX,block,MPI_STATUS_IGNORE);
      memcpy(domain->boxlo,domain_boxbuf  ,3*sizeof(double));
      memcpy(domain->boxhi,domain_boxbuf+3,3*sizeof(double));
      domain->set_global_box();
      domain->set_local_box();
      force->kspace->setup();
    }
    //K-processes to receive brick charge densities
    //this rank is the root of the gather and contributes no data itself:
    //use MPI_IN_PLACE instead of a null send buffer, which is not standard
    //conforming and rejected by some MPI implementations.  the data type
    //must be MPI_FFT_SCALAR to match the senders and single precision FFTs
    MPI_Igatherv(MPI_IN_PLACE,0,MPI_FFT_SCALAR,&density_brick_buf[nzlo_in][nylo_in][nxlo_in],
                 density_sizes,density_disps,MPI_FFT_SCALAR,0,block_density_brick,&mpi_requests_density);
    MPI_Wait(&mpi_requests_density,MPI_STATUS_IGNORE);

    FFT_SCALAR *buf_loc = &density_brick_buf[nzlo_in][nylo_in][nxlo_in];
    int idx = 0;
    int zlo, zhi, ylo, yhi, xlo, xhi;
    for (int k = 1; k<block_size; k++) {
      zlo = partitionInfoK[k][1];
      zhi = partitionInfoK[k][2];
      ylo = partitionInfoK[k][3];
      yhi = partitionInfoK[k][4];
      xlo = partitionInfoK[k][5];
      xhi = partitionInfoK[k][6];
      for(int zz=zlo; zz<=zhi; zz++)
        for(int yy=ylo; yy<=yhi; yy++)
          for(int xx=xlo; xx<=xhi; xx++)
            density_brick[zz][yy][xx] = buf_loc[idx++];
    }
  }
}

/* ----------------------------------------------------------------------
   communicate grid potentials to Rspace processes
------------------------------------------------------------------------- */
void PPPM_RK::k2r_comm(int eflag, int vflag)
{
  if (rproc) {
    //R-side waits on its send of flags, domain_box, and charge densities
    if (me_block==1) {
      MPI_Wait(&mpi_requests_flags,MPI_STATUS_IGNORE);
      if (domain->box_change)
        MPI_Wait(&mpi_requests_domain_box,MPI_STATUS_IGNORE);
    }
    MPI_Wait(&mpi_requests_density,MPI_STATUS_IGNORE);
    /***Done waiting on the R-process sends***/

    /***R-process wait on the receives***/
    //R-processes wait to receive the grid potentials
    waitReceiptGridPotentialsEV(eflag, vflag);
    //R-processes interpolate long-range forces from grid potentials
    compute_interpolate_forces(eflag, vflag);
  }
  if (!rproc) {
    //K-side communicates grid values to R-side
    post_sending_scatter_grid_potentials_ev(eflag, vflag);
    wait_sending_scatter_grid_potentials_ev(eflag, vflag);
  }
}

void PPPM_RK::compute_charge_densities(int eflag, int vflag)
{
  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  ev_init(eflag,vflag);

  if (evflag_atom && !peratom_allocate_flag) allocate_peratom();

  // if atom count has changed, update qsum and qsqsum

  if (atom->natoms != natoms_original) {
    qsum_qsq();
    natoms_original = atom->natoms;
  }

  // return if there are no charges

  if (qsqsum == 0.0) return;

  // convert atoms from box to lamda coords

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm:part2grid");
  }

  // find grid points for all my particles
  particle_map(); //Needed for poisson() and after
  // map my particle charge onto my local 3d density grid
  make_rho();  //Needed for poisson()

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  gc->reverse_comm(Grid3d::KSPACE,this,REVERSE_RHO,1,sizeof(FFT_SCALAR),
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);

}

// Assumes partition of processes as enforced in init through setupRKBlock
void PPPM_RK::compute(int eflag, int vflag)
{
  r2k_comm(eflag,vflag);
  if(!rproc) compute_grid_potentials(eflag, vflag);
  k2r_comm(eflag,vflag);
}

//Computes only the grid potentials from K-space processes
void PPPM_RK::compute_grid_potentials(int eflag, int vflag)
{
  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  ev_init(eflag,vflag);

  if (evflag_atom && !peratom_allocate_flag) allocate_peratom();

  // remap from 3d decomposition to FFT decomposition
  brick2fft();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  // also performs per-atom calculations via poisson_peratom()

  poisson();

  if (eflag_global) MPI_Allreduce(&energy,&energy_buf,1,MPI_DOUBLE,MPI_SUM,world);
  if (vflag_global) MPI_Allreduce(virial,virial_mpi_buf,6,MPI_DOUBLE,MPI_SUM,world);
}

/*Given updated grid potentials, finish computing the long-range forces
  from the R-space processes*/
void PPPM_RK::compute_interpolate_forces(int /*eflag*/, int /*vflag*/)
{
  int i,j;
  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  if (differentiation_flag == 1)
    gc->forward_comm(Grid3d::KSPACE,this,FORWARD_AD,1,sizeof(FFT_SCALAR),
                     gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  else
    gc->forward_comm(Grid3d::KSPACE,this,FORWARD_IK,3,sizeof(FFT_SCALAR),
                     gc_buf1,gc_buf2,MPI_FFT_SCALAR);

  // extra per-atom energy/virial communication

  if (evflag_atom) {
    if (differentiation_flag == 1 && vflag_atom)
      gc->forward_comm(Grid3d::KSPACE,this,FORWARD_AD_PERATOM,6,sizeof(FFT_SCALAR),
                       gc_buf1,gc_buf2,MPI_FFT_SCALAR);
    else if (differentiation_flag == 0)
      gc->forward_comm(Grid3d::KSPACE,this,FORWARD_IK_PERATOM,7,sizeof(FFT_SCALAR),
                       gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  }

  // calculate the force on my particles

  fieldforce();

  // extra per-atom energy/virial communication

  if (evflag_atom) fieldforce_peratom();

  // sum global energy across procs and add in volume-dependent term

  const double qscale = qqrd2e * scale;

  if (eflag_global) {
    energy = energy_buf;

    energy *= 0.5*volume;
    energy -= g_ewald*qsqsum/MY_PIS +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= qscale;
  }

  // sum global virial across procs

  if (vflag_global) {
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*volume*virial_mpi_buf[i];
  }

  // per-atom energy/virial
  // energy includes self-energy correction
  // ntotal accounts for TIP4P tallying eatom/vatom for ghost atoms

  if (evflag_atom) {
    double *q = atom->q;
    int nlocal = atom->nlocal;
    int ntotal = nlocal;
    if (tip4pflag) ntotal += atom->nghost;

    if (eflag_atom) {
      for (i = 0; i < nlocal; i++) {
        eatom[i] *= 0.5;
        eatom[i] -= g_ewald*q[i]*q[i]/MY_PIS + MY_PI2*q[i]*qsum /
          (g_ewald*g_ewald*volume);
        eatom[i] *= qscale;
      }
      for (i = nlocal; i < ntotal; i++) eatom[i] *= 0.5*qscale;
    }

    if (vflag_atom) {
      for (i = 0; i < ntotal; i++)
        for (j = 0; j < 6; j++) vatom[i][j] *= 0.5*qscale;
    }
  }

  // 2d slab correction

  if (slabflag == 1) slabcorr();

  // convert atoms back from lamda to box coords

  if (triclinic) domain->lamda2x(atom->nlocal);
}

