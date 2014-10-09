/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#include <stdio.h>

#define _type2frho MY_AP(coeff1)
#define _type2rhor MY_AP(coeff2)
#define _type2z2r MY_AP(coeff3)
#define _rdr MY_AP(rdr)
#define _rdrho MY_AP(rdrho)
#define _nr MY_AP(nr)
#define _nrho MY_AP(nrho)
#define _nfrho MY_AP(nfrho)
#define _nrhor MY_AP(nrhor)
#define _nz2r MY_AP(nz2r)
#define _frho_spline MY_AP(frho_spline)
#define _rhor_spline MY_AP(rhor_spline)
#define _z2r_spline MY_AP(z2r_spline)
#define _rho MY_AP(rho)
#define _fp MY_AP(fp)

__device__ __constant__ F_CFLOAT MY_AP(rdr);
__device__ __constant__ F_CFLOAT MY_AP(rdrho);
__device__ __constant__ int MY_AP(nr);
__device__ __constant__ int MY_AP(nrho);
__device__ __constant__ int MY_AP(nfrho);
__device__ __constant__ int MY_AP(nrhor);
__device__ __constant__ int MY_AP(nz2r);
__device__ __constant__ F_CFLOAT* MY_AP(frho_spline);
__device__ __constant__ F_CFLOAT* MY_AP(rhor_spline);
__device__ __constant__ F_CFLOAT* MY_AP(z2r_spline);
__device__ __constant__ F_CFLOAT* MY_AP(rho);
__device__ __constant__ F_CFLOAT* MY_AP(fp);

#define _rhor_spline_tex         MY_AP(rhor_spline_tex)
#if F_PRECISION == 1
texture<float4, 1> _rhor_spline_tex;
#else
texture<int4, 1> _rhor_spline_tex;
#endif


#define _z2r_spline_tex         MY_AP(z2r_spline_tex)
#if F_PRECISION == 1
texture<float4, 1> _z2r_spline_tex;
#else
texture<int4, 1> _z2r_spline_tex;
#endif



#include "pair_eam_cuda_cu.h"
#include "pair_eam_cuda_kernel_nc.cu"
#include <time.h>

int eam_buff_offset;
int rhor_spline_size;
void* rhor_spline_pointer;
int z2r_spline_size;
void* z2r_spline_pointer;


inline void BindEAMTextures(cuda_shared_data* sdata)
{
  _rhor_spline_tex.normalized = false;                      // access with normalized texture coordinates
  _rhor_spline_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
  _rhor_spline_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates

  const textureReference* rhor_spline_texture_ptr = &MY_AP(rhor_spline_tex);

#if F_PRECISION == 1
  cudaChannelFormatDesc channelDescRhor = cudaCreateChannelDesc<float4>();
  cudaBindTexture(0, rhor_spline_texture_ptr, rhor_spline_pointer, &channelDescRhor, rhor_spline_size);
#else
  cudaChannelFormatDesc channelDescRhor = cudaCreateChannelDesc<int4>();
  cudaBindTexture(0, rhor_spline_texture_ptr, rhor_spline_pointer, &channelDescRhor, rhor_spline_size);
#endif

  _z2r_spline_tex.normalized = false;                      // access with normalized texture coordinates
  _z2r_spline_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
  _z2r_spline_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates

  const textureReference* z2r_spline_texture_ptr = &MY_AP(z2r_spline_tex);

#if F_PRECISION == 1
  cudaChannelFormatDesc channelDescZ2r = cudaCreateChannelDesc<float4>();
  cudaBindTexture(0, z2r_spline_texture_ptr, z2r_spline_pointer, &channelDescZ2r, z2r_spline_size);
#else
  cudaChannelFormatDesc channelDescZ2r = cudaCreateChannelDesc<int4>();
  cudaBindTexture(0, z2r_spline_texture_ptr, z2r_spline_pointer, &channelDescZ2r, z2r_spline_size);
#endif

}

void Cuda_PairEAMCuda_UpdateBuffer(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist)
{
  CUT_CHECK_ERROR("Cuda_PairEAMCuda: before updateBuffer failed");
  int3 layout = getgrid(sneighlist->inum, 7 * sizeof(F_CFLOAT));
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);
  int size = (unsigned)(layout.y * layout.x) * 7 * sizeof(F_CFLOAT);

  if(sdata->buffersize < size) {
    MYDBG(printf("Cuda_PairEAMCuda Resizing Buffer at %p with %i kB to\n", sdata->buffer, sdata->buffersize);)

    if(sdata->buffer != NULL) cudaFree(sdata->buffer);

    cudaMalloc((void**)&sdata->buffer, size);
    sdata->buffersize = size;
    sdata->buffer_new++;
    MYDBG(printf("New buffer at %p with %i kB\n", sdata->buffer, sdata->buffersize);)
  }

  cudaMemcpyToSymbol(MY_AP(buffer), & sdata->buffer, sizeof(int*));
  CUT_CHECK_ERROR("Cuda_PairEAMCuda: updateBuffer failed");
}

void Cuda_PairEAMCuda_UpdateNeighbor(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist)
{
  cudaMemcpyToSymbol(MY_AP(neighbor_maxlocal) , & sneighlist->firstneigh.dim[0]  , sizeof(unsigned));
  cudaMemcpyToSymbol(MY_AP(firstneigh), & sneighlist->firstneigh.dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(ilist)     , & sneighlist->ilist     .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(inum)      , & sneighlist->inum               , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal               , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nmax)      , & sdata->atom.nmax               , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(numneigh)  , & sneighlist->numneigh  .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(neighbors)      , & sneighlist->neighbors  .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(maxneighbors)       , & sneighlist->maxneighbors     , sizeof(int));
}

void Cuda_PairEAMCuda_UpdateNmax(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist)
{
  CUT_CHECK_ERROR("Cuda_PairEAMCuda: before updateNmax failed");
  cudaMemcpyToSymbol(MY_AP(x)         , & sdata->atom.x         .dev_data, sizeof(X_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(x_type)         	, & sdata->atom.x_type    .dev_data, sizeof(X_CFLOAT4*));
  cudaMemcpyToSymbol(MY_AP(f)         			, & sdata->atom.f         .dev_data, sizeof(F_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(type)      			, & sdata->atom.type      .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(tag)      			, & sdata->atom.tag       .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(eatom)     			, & sdata->atom.eatom     .dev_data, sizeof(ENERGY_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(vatom)     			, & sdata->atom.vatom     .dev_data, sizeof(ENERGY_CFLOAT*));
  CUT_CHECK_ERROR("Cuda_PairEAMCuda: updateNmax failed");
}


void Cuda_PairEAMCuda_Init(cuda_shared_data* sdata, double rdr, double rdrho, int nfrho, int nrhor, int nr, int nrho, int nz2r,
                           void* frho_spline, void* rhor_spline, void* z2r_spline, void* rho, void* fp,
                           int* type2frho, int** type2z2r, int** type2rhor)
{
  // !! LAMMPS indexes atom types starting with 1 !!

  unsigned cuda_ntypes = sdata->atom.ntypes + 1;

  if(cuda_ntypes * cuda_ntypes > CUDA_MAX_TYPES2)
    printf("# CUDA: Cuda_PairEAMCuda_Init: you need %u types. this is more than %u "
           "(assumed at compile time). re-compile with -DCUDA_MAX_TYPES_PLUS_ONE=99 "
           "or ajust this in cuda_common.h\n", cuda_ntypes, CUDA_MAX_TYPES2);

  unsigned nI = sizeof(F_CFLOAT) * cuda_ntypes * cuda_ntypes;

  X_CFLOAT cutsq_global;
  cutsq_global = (X_CFLOAT)(sdata->pair.cut_global);
  cudaMemcpyToSymbol(MY_AP(cutsq_global)	, &cutsq_global  				, sizeof(X_CFLOAT));


  F_CFLOAT* coeff_buf = new F_CFLOAT[cuda_ntypes * cuda_ntypes];

  for(int i = 0; i < cuda_ntypes; i++) coeff_buf[i] = type2frho[i];

  cudaMemcpyToSymbol(MY_AP(coeff1)        , coeff_buf             , cuda_ntypes * sizeof(F_CFLOAT));

  for(int i = 0; i < cuda_ntypes * cuda_ntypes; i++) coeff_buf[i] = (&type2rhor[0][0])[i];

  cudaMemcpyToSymbol(MY_AP(coeff2)        , coeff_buf             , nI);

  for(int i = 0; i < cuda_ntypes * cuda_ntypes; i++) coeff_buf[i] = (&type2z2r[0][0])[i];

  cudaMemcpyToSymbol(MY_AP(coeff3)        , coeff_buf             , nI);

  delete [] coeff_buf;
  X_CFLOAT box_size[3] = {
    sdata->domain.subhi[0] - sdata->domain.sublo[0],
    sdata->domain.subhi[1] - sdata->domain.sublo[1],
    sdata->domain.subhi[2] - sdata->domain.sublo[2]
  };
  F_CFLOAT rdr_F = rdr;
  F_CFLOAT rdrho_F = rdrho;
  cudaMemcpyToSymbol(MY_AP(box_size)   , box_size                 , sizeof(X_CFLOAT) * 3);
  cudaMemcpyToSymbol(MY_AP(cuda_ntypes), & cuda_ntypes            , sizeof(unsigned));
  cudaMemcpyToSymbol(MY_AP(virial)     , &sdata->pair.virial.dev_data   , sizeof(ENERGY_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(eng_vdwl)     , &sdata->pair.eng_vdwl.dev_data   , sizeof(ENERGY_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(periodicity), sdata->domain.periodicity, sizeof(int) * 3);
  cudaMemcpyToSymbol(MY_AP(collect_forces_later), &sdata->pair.collect_forces_later  , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(rdr), &rdr_F, sizeof(F_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(rdrho), &rdrho_F, sizeof(F_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(nr), &nr, sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nrho), &nrho, sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nfrho), &nfrho, sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nrhor), &nrhor, sizeof(int));
  cudaMemcpyToSymbol(MY_AP(rho), &rho, sizeof(F_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(fp), &fp, sizeof(F_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(frho_spline), &frho_spline, sizeof(F_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(rhor_spline), &rhor_spline, sizeof(F_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(z2r_spline), &z2r_spline, sizeof(F_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(nrhor), &nrhor, sizeof(int));

  rhor_spline_size = nrhor * (nr + 1) * EAM_COEFF_LENGTH * sizeof(F_CFLOAT);
  z2r_spline_size = nz2r * (nr + 1) * EAM_COEFF_LENGTH * sizeof(F_CFLOAT);
  rhor_spline_pointer = rhor_spline;
  z2r_spline_pointer = z2r_spline;

  CUT_CHECK_ERROR("Cuda_PairEAMCuda: init failed");

}



void Cuda_PairEAM1Cuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag, int eflag_atom, int vflag_atom)
{

  if(sdata->atom.update_nmax)
    Cuda_PairEAMCuda_UpdateNmax(sdata, sneighlist);

  if(sdata->atom.update_neigh)
    Cuda_PairEAMCuda_UpdateNeighbor(sdata, sneighlist);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  if(sdata->buffer_new)
    Cuda_PairEAMCuda_UpdateBuffer(sdata, sneighlist);

  cudaMemcpyToSymbol(MY_AP(eatom)     			, & sdata->atom.eatom     .dev_data, sizeof(ENERGY_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(vatom)     			, & sdata->atom.vatom     .dev_data, sizeof(ENERGY_CFLOAT*));

  int sharedperproc = 0;

  if(eflag || eflag_atom) sharedperproc = 1;

  if(vflag || vflag_atom) sharedperproc = 7;

  int3 layout = getgrid(sneighlist->inum, sharedperproc * sizeof(ENERGY_CFLOAT));
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  eam_buff_offset = grid.x * grid.y;

  BindXTypeTexture(sdata);
  BindEAMTextures(sdata); // initialize only on first call


  MYDBG(printf("# CUDA: Cuda_PairEAMCuda: kernel start eflag: %i vflag: %i\n", eflag, vflag);)
  CUT_CHECK_ERROR("Cuda_PairEAMCuda: pre pair Kernel 1 problems before kernel invocation");
  PairEAMCuda_Kernel1 <<< grid, threads, sharedperproc* sizeof(ENERGY_CFLOAT)*threads.x>>> (eflag, vflag, eflag_atom, vflag_atom);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_PairEAMCuda: pair Kernel 1 execution failed");



  MYDBG(printf("# CUDA: Cuda_PairEAMCoulLongCuda: kernel done\n");)

}

void Cuda_PairEAM2Cuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag, int eflag_atom, int vflag_atom)
{
  int sharedperproc = 0;

  if(eflag || eflag_atom) sharedperproc = 1;

  if(vflag || vflag_atom) sharedperproc = 7;

  int3 layout = getgrid(sneighlist->inum, sharedperproc * sizeof(ENERGY_CFLOAT));
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  BindXTypeTexture(sdata);
  BindEAMTextures(sdata); // initialize only on first call
  // initialize only on first call
  sdata->pair.lastgridsize = grid.x * grid.y;
  sdata->pair.n_energy_virial = sharedperproc;

  MYDBG(printf("# CUDA: Cuda_PairEAMCuda: kernel start eflag: %i vflag: %i\n", eflag, vflag);)
  CUT_CHECK_ERROR("Cuda_PairEAMCuda: pre pair Kernel 2 problems before kernel invocation");
  PairEAMCuda_Kernel2 <<< grid, threads, sharedperproc* sizeof(ENERGY_CFLOAT)*threads.x>>> (eflag, vflag, eflag_atom, vflag_atom);
  CUT_CHECK_ERROR("Cuda_PairEAMCuda: pair Kernel 2 start failed");
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_PairEAMCuda: pair Kernel 2 execution failed");

  if(eflag || vflag) {
    int n = grid.x * grid.y;
    grid.x = sharedperproc;
    grid.y = 1;
    threads.x = 256;
    MY_AP(PairVirialCompute_reduce) <<< grid, threads, threads.x* sizeof(ENERGY_CFLOAT)*sharedperproc>>>(n);
    cudaThreadSynchronize();
    CUT_CHECK_ERROR("Cuda_PairEAMCuda: virial compute Kernel execution failed");
  }

  MYDBG(printf("# CUDA: Cuda_PairEAMCoulLongCuda: kernel done\n");)

}

void Cuda_PairEAMCuda_PackComm(cuda_shared_data* sdata, int n, int iswap, void* buf_send)
{
  int3 layout = getgrid(n, 0);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);
  F_CFLOAT* buf = (F_CFLOAT*)(& ((double*)sdata->buffer)[eam_buff_offset]);

  PairEAMCuda_PackComm_Kernel <<< grid, threads, 0>>> ((int*) sdata->comm.sendlist.dev_data, n
      , sdata->comm.maxlistlength, iswap, buf);
  cudaThreadSynchronize();
  cudaMemcpy(buf_send, buf, n* sizeof(F_CFLOAT), cudaMemcpyDeviceToHost);
  cudaThreadSynchronize();
}

void Cuda_PairEAMCuda_UnpackComm(cuda_shared_data* sdata, int n, int first, void* buf_recv, void* fp)
{
  F_CFLOAT* fp_first = &(((F_CFLOAT*) fp)[first]);
  cudaMemcpy(fp_first, buf_recv, n * sizeof(F_CFLOAT), cudaMemcpyHostToDevice);
}

#undef _type2frho
#undef _type2rhor
#undef _type2z2r


/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   need i < nlocal test since called by bond_quartic and dihedral_charmm
------------------------------------------------------------------------- */

