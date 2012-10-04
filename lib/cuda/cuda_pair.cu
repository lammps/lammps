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

enum PAIR_FORCES {PAIR_NONE, PAIR_BORN, PAIR_BUCK, PAIR_CG_CMM, PAIR_LJ_CHARMM, PAIR_LJ_CLASS2, PAIR_LJ_CUT, PAIR_LJ_EXPAND, PAIR_LJ_GROMACS, PAIR_LJ_SMOOTH, PAIR_LJ96_CUT, PAIR_MORSE, PAIR_MORSE_R6};
enum COUL_FORCES {COUL_NONE, COUL_CHARMM, COUL_CHARMM_IMPLICIT, COUL_CUT, COUL_LONG, COUL_DEBYE, COUL_GROMACS, COUL_SPECIAL};
#define DATA_NONE 0
#define DATA_V 1
#define DATA_TAG 2
#define DATA_RMASS 4
#define DATA_MASS 8
#define DATA_TORQUE 16
#define DATA_OMEGA 32
#define DATA_RADIUS 64
#define DATA_DENSITY 128
#define DATA_MASK 256
#define DATA_V_RADIUS 512
#define DATA_OMEGA_RMASS 1024

#define NEIGHMASK 0x3FFFFFFF

#define MY_PREFIX cuda_pair
#define IncludeCommonNeigh
#include "cuda_shared.h"
#include "cuda_common.h"
#include "cuda_wrapper_cu.h"
#include "crm_cuda_utils.cu"

//constants used by multiple forces

//general
#define _cutsq MY_AP(cutsq)
#define _offset MY_AP(offset)
#define _special_lj MY_AP(special_lj)
#define _special_coul MY_AP(special_coul)
#define _cutsq_global MY_AP(cutsq_global)
#define _collect_forces_later MY_AP(collect_forces_later)

__device__ __constant__ X_FLOAT _cutsq[CUDA_MAX_TYPES2];
__device__ __constant__ ENERGY_FLOAT _offset[CUDA_MAX_TYPES2];
__device__ __constant__ F_FLOAT _special_lj[4];
__device__ __constant__ F_FLOAT _special_coul[4];
__device__ __constant__ X_FLOAT _cutsq_global;
__device__ __constant__ int _collect_forces_later;

__device__ __constant__ F_FLOAT MY_AP(coeff1)[CUDA_MAX_TYPES2]; //pair force coefficients in case ntypes < CUDA_MAX_TYPES (coeffs fit into constant space)
__device__ __constant__ F_FLOAT MY_AP(coeff2)[CUDA_MAX_TYPES2];
__device__ __constant__ F_FLOAT MY_AP(coeff3)[CUDA_MAX_TYPES2];
__device__ __constant__ F_FLOAT MY_AP(coeff4)[CUDA_MAX_TYPES2];
__device__ __constant__ F_FLOAT MY_AP(coeff5)[CUDA_MAX_TYPES2];


__device__ __constant__ F_FLOAT* MY_AP(coeff1_gm); //pair force coefficients in case ntypes > CUDA_MAX_TYPES (coeffs do not fit into constant space)
__device__ __constant__ F_FLOAT* MY_AP(coeff2_gm);
__device__ __constant__ F_FLOAT* MY_AP(coeff3_gm);
__device__ __constant__ F_FLOAT* MY_AP(coeff4_gm);
__device__ __constant__ F_FLOAT* MY_AP(coeff5_gm);
__device__ __constant__ F_FLOAT* MY_AP(coeff6_gm);
__device__ __constant__ F_FLOAT* MY_AP(coeff7_gm);
__device__ __constant__ F_FLOAT* MY_AP(coeff8_gm);
__device__ __constant__ F_FLOAT* MY_AP(coeff9_gm);
__device__ __constant__ F_FLOAT* MY_AP(coeff10_gm);

#define _coeff1_gm_tex         MY_AP(coeff1_gm_tex)
#if F_PRECISION == 1
texture<float> _coeff1_gm_tex;
#else
texture<int2, 1> _coeff1_gm_tex;
#endif

#define _coeff2_gm_tex         MY_AP(coeff2_gm_tex)
#if F_PRECISION == 1
texture<float> _coeff2_gm_tex;
#else
texture<int2, 1> _coeff2_gm_tex;
#endif

#define _coeff3_gm_tex         MY_AP(coeff3_gm_tex)
#if F_PRECISION == 1
texture<float> _coeff3_gm_tex;
#else
texture<int2, 1> _coeff3_gm_tex;
#endif

#define _coeff4_gm_tex         MY_AP(coeff4_gm_tex)
#if F_PRECISION == 1
texture<float> _coeff4_gm_tex;
#else
texture<int2, 1> _coeff4_gm_tex;
#endif

#define _coeff5_gm_tex         MY_AP(coeff5_gm_tex)
#if F_PRECISION == 1
texture<float> _coeff5_gm_tex;
#else
texture<int2, 1> _coeff5_gm_tex;
#endif

#define _coeff6_gm_tex         MY_AP(coeff6_gm_tex)
#if F_PRECISION == 1
texture<float> _coeff6_gm_tex;
#else
texture<int2, 1> _coeff6_gm_tex;
#endif

#define _coeff7_gm_tex         MY_AP(coeff7_gm_tex)
#if F_PRECISION == 1
texture<float> _coeff7_gm_tex;
#else
texture<int2, 1> _coeff7_gm_tex;
#endif

#define _coeff8_gm_tex         MY_AP(coeff8_gm_tex)
#if F_PRECISION == 1
texture<float> _coeff8_gm_tex;
#else
texture<int2, 1> _coeff8_gm_tex;
#endif

#define _coeff9_gm_tex         MY_AP(coeff9_gm_tex)
#if F_PRECISION == 1
texture<float> _coeff9_gm_tex;
#else
texture<int2, 1> _coeff9_gm_tex;
#endif

#define _coeff10_gm_tex         MY_AP(coeff10_gm_tex)
#if F_PRECISION == 1
texture<float> _coeff10_gm_tex;
#else
texture<int2, 1> _coeff10_gm_tex;
#endif

//if more than 5 coefficients are needed for a pair potential add them here


//coulomb
#define _cut_coulsq MY_AP(cut_coulsq)
#define _cut_coulsq_global MY_AP(cut_coulsq_global)
#define _g_ewald MY_AP(g_ewald)
#define _qqrd2e MY_AP(qqrd2e)
#define _kappa MY_AP(kappa)
__device__ __constant__ X_FLOAT _cut_coulsq[CUDA_MAX_TYPES2];
__device__ __constant__ X_FLOAT _cut_coulsq_global;
__device__ __constant__ F_FLOAT _g_ewald;
__device__ __constant__ F_FLOAT _qqrd2e;
__device__ __constant__ F_FLOAT _kappa;

//inner cutoff
#define _cut_innersq MY_AP(cut_innersq)
#define _cut_innersq_global MY_AP(cut_innersq_global)
__device__ __constant__ X_FLOAT _cut_innersq[CUDA_MAX_TYPES2];
__device__ __constant__ X_FLOAT _cut_innersq_global;


template <const PAIR_FORCES pair_type, const COUL_FORCES coul_type, const unsigned int extended_data>
__global__ void Pair_Kernel_TpA(int eflag, int vflag, int eflag_atom, int vflag_atom);

template <const PAIR_FORCES pair_type, const COUL_FORCES coul_type, const unsigned int extended_data>
__global__ void Pair_Kernel_BpA(int eflag, int vflag, int eflag_atom, int vflag_atom);

template <const PAIR_FORCES pair_type, const COUL_FORCES coul_type, const unsigned int extended_data>
__global__ void Pair_Kernel_TpA_opt(int eflag, int vflag, int eflag_atom, int vflag_atom, int comm_phase);

template <const PAIR_FORCES pair_type, const COUL_FORCES coul_type, const unsigned int extended_data>
__global__ void Pair_Kernel_BpA_opt(int eflag, int vflag, int eflag_atom, int vflag_atom, int comm_phase);

#include <stdio.h>
#include "cuda_pair_cu.h"
#include "cuda_pair_virial_kernel_nc.cu"

//Functions which are shared by pair styles

//Update Buffersize
void Cuda_UpdateBuffer(cuda_shared_data* sdata, int size)
{
  CUT_CHECK_ERROR("Cuda_Pair_UpdateBuffer_AllStyles: before updateBuffer failed");

  if(sdata->buffersize < size) {
    MYDBG(printf("Resizing Buffer at %p with %i kB to\n", sdata->buffer, sdata->buffersize);)
    CudaWrapper_FreeCudaData(sdata->buffer, sdata->buffersize);
    sdata->buffer = CudaWrapper_AllocCudaData(size);
    sdata->buffersize = size;
    sdata->buffer_new++;
    MYDBG(printf("New buffer at %p with %i kB\n", sdata->buffer, sdata->buffersize);)
  }

  cudaMemcpyToSymbol(MY_AP(buffer), & sdata->buffer, sizeof(int*));
  CUT_CHECK_ERROR("Cuda_Pair_UpdateBuffer_AllStyles failed");
}

void Cuda_Pair_UpdateNeighbor_AllStyles(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist)
{
  //Neighbor
  cudaMemcpyToSymbol(MY_AP(neighbor_maxlocal)  , & sneighlist->firstneigh.dim[0]  , sizeof(unsigned));
  cudaMemcpyToSymbol(MY_AP(firstneigh)     , & sneighlist->firstneigh.dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(ilist)          , & sneighlist->ilist     .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(inum)           , & sneighlist->inum               , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(numneigh)       , & sneighlist->numneigh  .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(neighbors)      , & sneighlist->neighbors  .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(maxneighbors)       , & sneighlist->maxneighbors     , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(overlap_comm)       , & sdata->overlap_comm, sizeof(int));

  if(sdata->overlap_comm) {
    cudaMemcpyToSymbol(MY_AP(numneigh_border)  , & sneighlist->numneigh_border .dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(numneigh_inner)   , & sneighlist->numneigh_inner  .dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(neighbors_border) , & sneighlist->neighbors_border.dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(neighbors_inner)  , & sneighlist->neighbors_inner .dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(ilist_border)     , & sneighlist->ilist_border    .dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(inum_border)      , & sneighlist->inum_border     .dev_data, sizeof(int*));
  }

}
//Update constants after nmax change which are generally needed by all pair styles
void Cuda_Pair_UpdateNmax_AllStyles(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist)
{
  CUT_CHECK_ERROR("Cuda_Pair_UpdateNmax_AllStyles: Begin");

  //System
  cudaMemcpyToSymbol(MY_AP(nlocal)    			, & sdata->atom.nlocal             , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nall)      			, & sdata->atom.nall               , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nmax)      			, & sdata->atom.nmax               , sizeof(int));

  //Atom
  cudaMemcpyToSymbol(MY_AP(x)         			, & sdata->atom.x         .dev_data, sizeof(X_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(x_type)         	, & sdata->atom.x_type    .dev_data, sizeof(X_FLOAT4*));
  cudaMemcpyToSymbol(MY_AP(f)         			, & sdata->atom.f         .dev_data, sizeof(F_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(type)      			, & sdata->atom.type      .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(q)         			, & sdata->atom.q         .dev_data, sizeof(F_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(tag)      			, & sdata->atom.tag       .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(eatom)     			, & sdata->atom.eatom     .dev_data, sizeof(ENERGY_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(vatom)     			, & sdata->atom.vatom     .dev_data, sizeof(ENERGY_FLOAT*));


  //Other
  cudaMemcpyToSymbol(MY_AP(debugdata)      , & sdata->debugdata      , sizeof(int*));
  CUT_CHECK_ERROR("Cuda_Pair_UpdateNmax_AllStyles: End");
}

//Initialisation of GPU Constants which rarely change
void Cuda_Pair_Init_AllStyles(cuda_shared_data* sdata, int ncoeff, bool need_q = false, bool use_global_params = false, bool need_innercut = false, bool need_cut = true)
{
  unsigned cuda_ntypes = sdata->atom.ntypes + 1;
  unsigned cuda_ntypes2 = cuda_ntypes * cuda_ntypes;
  unsigned n = sizeof(F_FLOAT) * cuda_ntypes2;
  unsigned nx = sizeof(X_FLOAT) * cuda_ntypes2;

  //check if enough constant memory is available
  if((cuda_ntypes2 > CUDA_MAX_TYPES2) && !use_global_params)
    printf("# CUDA: Cuda_Pair_Init: you need %u types. this is more than %u "
           "(assumed at compile time). re-compile with -DCUDA_MAX_TYPES_PLUS_ONE=32 "
           "or ajust this in cuda_common.h\n", cuda_ntypes, CUDA_MAX_TYPES_PLUS_ONE - 1);

  if((cuda_ntypes2 > CUDA_MAX_TYPES2) && !use_global_params)
    exit(0);

  //type conversion of cutoffs and parameters
  if(need_cut) {
    X_FLOAT cutsq[cuda_ntypes2];

    for(int i = 1; i <= sdata->atom.ntypes; ++i) {
      for(int j = 1; j <= sdata->atom.ntypes; ++j) {
        cutsq[i * cuda_ntypes + j] = (X_FLOAT)(sdata->pair.cut_global * sdata->pair.cut_global);
      }
    }

    int cutsqdiffer = 0;
    X_FLOAT cutsq_global;
    cutsq_global = (X_FLOAT)(sdata->pair.cut_global * sdata->pair.cut_global);

    if(sdata->pair.cut) {
      for(int i = 1; i <= sdata->atom.ntypes; ++i) {
        for(int j = i; j <= sdata->atom.ntypes; ++j) {
          if(sdata->pair.cut[i][j] > 1e-6) {
            cutsq[i * cuda_ntypes + j] = (X_FLOAT)(sdata->pair.cut[i][j] * sdata->pair.cut[i][j]);
            cutsq[j * cuda_ntypes + i] = (X_FLOAT)(sdata->pair.cut[i][j] * sdata->pair.cut[i][j]);
          }

          if(i == 1 && j == 1) cutsq_global = cutsq[i * cuda_ntypes + j];

          if((cutsq_global - cutsq[i * cuda_ntypes + j]) * (cutsq_global - cutsq[i * cuda_ntypes + j]) > 1e-6)
            cutsqdiffer++;
        }
      }
    }

    if(sdata->pair.cutsq) {
      for(int i = 1; i <= sdata->atom.ntypes; ++i) {
        for(int j = i; j <= sdata->atom.ntypes; ++j) {
          if(sdata->pair.cut[i][j] > 1e-6) {
            cutsq[i * cuda_ntypes + j] = (X_FLOAT)(sdata->pair.cutsq[i][j]);
            cutsq[j * cuda_ntypes + i] = (X_FLOAT)(sdata->pair.cutsq[i][j]);
          }

          if(i == 1 && j == 1) cutsq_global = cutsq[i * cuda_ntypes + j];

          if((cutsq_global - cutsq[i * cuda_ntypes + j]) * (cutsq_global - cutsq[i * cuda_ntypes + j]) > 1e-6)
            cutsqdiffer++;
        }
      }
    }

    //printf("CUTSQGLOB: %i %e\n",cutsqdiffer,cutsq_global);
    if(cutsqdiffer) {

      cutsq_global = -1.0;
      cudaMemcpyToSymbol(MY_AP(cutsq)      	, cutsq                    		, nx);
    }

    cudaMemcpyToSymbol(MY_AP(cutsq_global)	, &cutsq_global  				, sizeof(X_FLOAT));
  }

  if(need_innercut) {
    X_FLOAT cut_innersq[cuda_ntypes2];

    for(int i = 1; i <= sdata->atom.ntypes; ++i) {
      for(int j = 1; j <= sdata->atom.ntypes; ++j) {
        cut_innersq[i * cuda_ntypes + j] = (X_FLOAT)(sdata->pair.cut_inner_global * sdata->pair.cut_inner_global);
      }
    }

    int cutsqdiffer = 0;
    X_FLOAT cut_innersq_global;
    cut_innersq_global = (X_FLOAT)(sdata->pair.cut_inner_global * sdata->pair.cut_inner_global);

    if(sdata->pair.cut_inner) {
      for(int i = 1; i <= sdata->atom.ntypes; ++i) {
        for(int j = i; j <= sdata->atom.ntypes; ++j) {
          if(sdata->pair.cut_inner[i][j] > 1e-6) {
            cut_innersq[i * cuda_ntypes + j] = (X_FLOAT)(sdata->pair.cut_inner[i][j] * sdata->pair.cut_inner[i][j]);
            cut_innersq[j * cuda_ntypes + i] = (X_FLOAT)(sdata->pair.cut_inner[i][j] * sdata->pair.cut_inner[i][j]);
          }

          if(i == 1 && j == 1) cut_innersq_global = cut_innersq[i * cuda_ntypes + j];

          if((cut_innersq_global - cut_innersq[i * cuda_ntypes + j]) * (cut_innersq_global - cut_innersq[i * cuda_ntypes + j]) > 1e-6)
            cutsqdiffer++;
        }
      }
    }

    if(cutsqdiffer) {
      cut_innersq_global = -1.0;
      cudaMemcpyToSymbol(MY_AP(cut_innersq)      	, cut_innersq                    		, nx);
    }

    cudaMemcpyToSymbol(MY_AP(cut_innersq_global)	, &cut_innersq_global  				, sizeof(X_FLOAT));
  }

  if(need_q) {
    X_FLOAT cut_coulsq[cuda_ntypes2];

    for(int i = 1; i <= sdata->atom.ntypes; ++i) {
      for(int j = 1; j <= sdata->atom.ntypes; ++j) {
        cut_coulsq[i * cuda_ntypes + j] = (X_FLOAT)(sdata->pair.cut_coul_global * sdata->pair.cut_coul_global);
      }
    }

    int cutsqdiffer = 0;
    X_FLOAT cut_coulsq_global;
    cut_coulsq_global = (X_FLOAT)(sdata->pair.cut_coul_global * sdata->pair.cut_coul_global);

    if(sdata->pair.cut_coulsq_global > cut_coulsq_global)  cut_coulsq_global = (X_FLOAT) sdata->pair.cut_coulsq_global;

    if(sdata->pair.cut_coul) {
      for(int i = 1; i <= sdata->atom.ntypes; ++i) {
        for(int j = i; j <= sdata->atom.ntypes; ++j) {
          if(sdata->pair.cut_coul[i][j] > 1e-6) {
            cut_coulsq[i * cuda_ntypes + j] = (X_FLOAT)(sdata->pair.cut_coul[i][j] * sdata->pair.cut_coul[i][j]);
            cut_coulsq[j * cuda_ntypes + i] = (X_FLOAT)(sdata->pair.cut_coul[i][j] * sdata->pair.cut_coul[i][j]);
          }

          if(i == 1 && j == 1) cut_coulsq_global = cut_coulsq[i * cuda_ntypes + j];

          if((cut_coulsq_global - cut_coulsq[i * cuda_ntypes + j]) * (cut_coulsq_global - cut_coulsq[i * cuda_ntypes + j]) > 1e-6)
            cutsqdiffer++;
        }
      }
    }

    if(cutsqdiffer) {
      cut_coulsq_global = -1.0;
      cudaMemcpyToSymbol(MY_AP(cut_coulsq)      	, cut_coulsq                    		, nx);
    }

    cudaMemcpyToSymbol(MY_AP(cut_coulsq_global), &cut_coulsq_global  					, sizeof(X_FLOAT));
  }

  CUT_CHECK_ERROR("Cuda_Pair: init pre Coeff failed");

  if(ncoeff > 0) {
    F_FLOAT coeff1[cuda_ntypes2];

    for(int i = 1; i <= sdata->atom.ntypes; ++i) {
      for(int j = 1; j <= sdata->atom.ntypes; ++j) {
        coeff1[i * cuda_ntypes + j] = (F_FLOAT) sdata->pair.coeff1[i][j];
      }
    }

    if(use_global_params) {
      cudaMemcpyToSymbol(MY_AP(coeff1_gm)  , &sdata->pair.coeff1_gm.dev_data   , sizeof(F_FLOAT*));
      cudaMemcpy((sdata->pair.coeff1_gm.dev_data), coeff1, n, cudaMemcpyHostToDevice);

      _coeff1_gm_tex.normalized = false;                      // access with normalized texture coordinates
      _coeff1_gm_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
      _coeff1_gm_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
      const textureReference* coeff1_gm_texture_ptr = &MY_AP(coeff1_gm_tex);
      CUT_CHECK_ERROR("Cuda_Pair: init Coeff0 a failed");

#if F_PRECISION == 1
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<float>();
      CUT_CHECK_ERROR("Cuda_Pair: init Coeff0 b failed");
      cudaBindTexture(0, coeff1_gm_texture_ptr, sdata->pair.coeff1_gm.dev_data, &channelDescXType, sdata->atom.nmax * sizeof(F_FLOAT));
      CUT_CHECK_ERROR("Cuda_Pair: init Coeff0 c failed");
#else
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<int2>();
      CUT_CHECK_ERROR("Cuda_Pair: init Coeff0 b-d failed");
      cudaBindTexture(0, coeff1_gm_texture_ptr, sdata->pair.coeff1_gm.dev_data, &channelDescXType, sdata->atom.nmax * 2 * sizeof(int2));
      CUT_CHECK_ERROR("Cuda_Pair: init Coeff0 c-d failed");
#endif

    } else
      cudaMemcpyToSymbol(MY_AP(coeff1), coeff1 , n);
  }

  CUT_CHECK_ERROR("Cuda_Pair: init Coeff0 failed");

  if(ncoeff > 1) {
    F_FLOAT coeff2[cuda_ntypes2];

    for(int i = 1; i <= sdata->atom.ntypes; ++i) {
      for(int j = 1; j <= sdata->atom.ntypes; ++j) {
        coeff2[i * cuda_ntypes + j] = (F_FLOAT) sdata->pair.coeff2[i][j];
      }
    }

    if(use_global_params) {
      cudaMemcpyToSymbol(MY_AP(coeff2_gm)  , &sdata->pair.coeff2_gm.dev_data   , sizeof(F_FLOAT*));
      cudaMemcpy(sdata->pair.coeff2_gm.dev_data, coeff2, n, cudaMemcpyHostToDevice);

      _coeff2_gm_tex.normalized = false;                      // access with normalized texture coordinates
      _coeff2_gm_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
      _coeff2_gm_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
      const textureReference* coeff2_gm_texture_ptr = &MY_AP(coeff2_gm_tex);

#if F_PRECISION == 1
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<float>();
      cudaBindTexture(0, coeff2_gm_texture_ptr, sdata->pair.coeff2_gm.dev_data, &channelDescXType, sdata->atom.nmax * sizeof(F_FLOAT));
#else
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<int2>();
      cudaBindTexture(0, coeff2_gm_texture_ptr, sdata->pair.coeff2_gm.dev_data, &channelDescXType, sdata->atom.nmax * 2 * sizeof(int2));
#endif

    } else
      cudaMemcpyToSymbol(MY_AP(coeff2), coeff2 , n);
  }

  CUT_CHECK_ERROR("Cuda_Pair: init Coeff1 failed");

  if(ncoeff > 2) {
    F_FLOAT coeff3[cuda_ntypes2];

    for(int i = 1; i <= sdata->atom.ntypes; ++i) {
      for(int j = 1; j <= sdata->atom.ntypes; ++j) {
        coeff3[i * cuda_ntypes + j] = (F_FLOAT) sdata->pair.coeff3[i][j];
      }
    }

    if(use_global_params) {
      cudaMemcpyToSymbol(MY_AP(coeff3_gm)  , &sdata->pair.coeff3_gm.dev_data   , sizeof(F_FLOAT*));
      cudaMemcpy(sdata->pair.coeff3_gm.dev_data, coeff3, n, cudaMemcpyHostToDevice);
      _coeff3_gm_tex.normalized = false;                      // access with normalized texture coordinates
      _coeff3_gm_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
      _coeff3_gm_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
      const textureReference* coeff3_gm_texture_ptr = &MY_AP(coeff3_gm_tex);

#if F_PRECISION == 1
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<float>();
      cudaBindTexture(0, coeff3_gm_texture_ptr, sdata->pair.coeff3_gm.dev_data, &channelDescXType, sdata->atom.nmax * sizeof(F_FLOAT));
#else
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<int2>();
      cudaBindTexture(0, coeff3_gm_texture_ptr, sdata->pair.coeff3_gm.dev_data, &channelDescXType, sdata->atom.nmax * 2 * sizeof(int2));
#endif
    } else
      cudaMemcpyToSymbol(MY_AP(coeff3), coeff3 , n);
  }

  CUT_CHECK_ERROR("Cuda_Pair: init Coeff3 failed");

  if(ncoeff > 3) {
    F_FLOAT coeff4[cuda_ntypes2];

    for(int i = 1; i <= sdata->atom.ntypes; ++i) {
      for(int j = 1; j <= sdata->atom.ntypes; ++j) {
        coeff4[i * cuda_ntypes + j] = (F_FLOAT) sdata->pair.coeff4[i][j];
      }
    }

    if(use_global_params) {
      cudaMemcpyToSymbol(MY_AP(coeff4_gm)  , &sdata->pair.coeff4_gm.dev_data   , sizeof(F_FLOAT*));
      cudaMemcpy(sdata->pair.coeff4_gm.dev_data, coeff4, n, cudaMemcpyHostToDevice);
      _coeff4_gm_tex.normalized = false;                      // access with normalized texture coordinates
      _coeff4_gm_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
      _coeff4_gm_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
      const textureReference* coeff4_gm_texture_ptr = &MY_AP(coeff4_gm_tex);

#if F_PRECISION == 1
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<float>();
      cudaBindTexture(0, coeff4_gm_texture_ptr, sdata->pair.coeff4_gm.dev_data, &channelDescXType, sdata->atom.nmax * sizeof(F_FLOAT));
#else
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<int2>();
      cudaBindTexture(0, coeff4_gm_texture_ptr, sdata->pair.coeff4_gm.dev_data, &channelDescXType, sdata->atom.nmax * 2 * sizeof(int2));
#endif
    } else
      cudaMemcpyToSymbol(MY_AP(coeff4), coeff4 , n);
  }

  CUT_CHECK_ERROR("Cuda_Pair: init Coeff4 failed");

  if(ncoeff > 4) {
    F_FLOAT coeff5[cuda_ntypes2];

    for(int i = 1; i <= sdata->atom.ntypes; ++i) {
      for(int j = 1; j <= sdata->atom.ntypes; ++j) {
        coeff5[i * cuda_ntypes + j] = (F_FLOAT) sdata->pair.coeff5[i][j];
      }
    }

    if(use_global_params) {
      cudaMemcpyToSymbol(MY_AP(coeff5_gm)  , &sdata->pair.coeff5_gm.dev_data   , sizeof(F_FLOAT*));
      cudaMemcpy(sdata->pair.coeff5_gm.dev_data, coeff5, n, cudaMemcpyHostToDevice);
      _coeff5_gm_tex.normalized = false;                      // access with normalized texture coordinates
      _coeff5_gm_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
      _coeff5_gm_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
      const textureReference* coeff5_gm_texture_ptr = &MY_AP(coeff5_gm_tex);

#if F_PRECISION == 1
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<float>();
      cudaBindTexture(0, coeff5_gm_texture_ptr, sdata->pair.coeff5_gm.dev_data, &channelDescXType, sdata->atom.nmax * sizeof(F_FLOAT));
#else
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<int2>();
      cudaBindTexture(0, coeff5_gm_texture_ptr, sdata->pair.coeff5_gm.dev_data, &channelDescXType, sdata->atom.nmax * 2 * sizeof(int2));
#endif
    } else
      cudaMemcpyToSymbol(MY_AP(coeff5), coeff5 , n);
  }

  CUT_CHECK_ERROR("Cuda_Pair: init Coeff5 failed");

  if(ncoeff > 5) {
    F_FLOAT coeff6[cuda_ntypes2];

    for(int i = 1; i <= sdata->atom.ntypes; ++i) {
      for(int j = 1; j <= sdata->atom.ntypes; ++j) {
        coeff6[i * cuda_ntypes + j] = (F_FLOAT) sdata->pair.coeff6[i][j];
      }
    }

    if(use_global_params) {
      cudaMemcpyToSymbol(MY_AP(coeff6_gm)  , &sdata->pair.coeff6_gm.dev_data   , sizeof(F_FLOAT*));
      cudaMemcpy(sdata->pair.coeff6_gm.dev_data, coeff6, n, cudaMemcpyHostToDevice);
      _coeff6_gm_tex.normalized = false;                      // access with normalized texture coordinates
      _coeff6_gm_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
      _coeff6_gm_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
      const textureReference* coeff6_gm_texture_ptr = &MY_AP(coeff6_gm_tex);

#if F_PRECISION == 1
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<float>();
      cudaBindTexture(0, coeff6_gm_texture_ptr, sdata->pair.coeff6_gm.dev_data, &channelDescXType, sdata->atom.nmax * sizeof(F_FLOAT));
#else
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<int2>();
      cudaBindTexture(0, coeff6_gm_texture_ptr, sdata->pair.coeff6_gm.dev_data, &channelDescXType, sdata->atom.nmax * 2 * sizeof(int2));
#endif
    }
  }

  CUT_CHECK_ERROR("Cuda_Pair: init Coeff6 failed");

  if(ncoeff > 6) {
    F_FLOAT coeff7[cuda_ntypes2];

    for(int i = 1; i <= sdata->atom.ntypes; ++i) {
      for(int j = 1; j <= sdata->atom.ntypes; ++j) {
        coeff7[i * cuda_ntypes + j] = (F_FLOAT) sdata->pair.coeff7[i][j];
      }
    }

    if(use_global_params) {
      cudaMemcpyToSymbol(MY_AP(coeff7_gm)  , &sdata->pair.coeff7_gm.dev_data   , sizeof(F_FLOAT*));
      cudaMemcpy(sdata->pair.coeff7_gm.dev_data, coeff7, n, cudaMemcpyHostToDevice);
      _coeff7_gm_tex.normalized = false;                      // access with normalized texture coordinates
      _coeff7_gm_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
      _coeff7_gm_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
      const textureReference* coeff7_gm_texture_ptr = &MY_AP(coeff7_gm_tex);

#if F_PRECISION == 1
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<float>();
      cudaBindTexture(0, coeff7_gm_texture_ptr, sdata->pair.coeff7_gm.dev_data, &channelDescXType, sdata->atom.nmax * sizeof(F_FLOAT));
#else
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<int2>();
      cudaBindTexture(0, coeff7_gm_texture_ptr, sdata->pair.coeff7_gm.dev_data, &channelDescXType, sdata->atom.nmax * 2 * sizeof(int2));
#endif
    }
  }

  CUT_CHECK_ERROR("Cuda_Pair: init Coeff7 failed");

  if(ncoeff > 7) {
    F_FLOAT coeff8[cuda_ntypes2];

    for(int i = 1; i <= sdata->atom.ntypes; ++i) {
      for(int j = 1; j <= sdata->atom.ntypes; ++j) {
        coeff8[i * cuda_ntypes + j] = (F_FLOAT) sdata->pair.coeff8[i][j];
      }
    }

    if(use_global_params) {
      cudaMemcpyToSymbol(MY_AP(coeff8_gm)  , &sdata->pair.coeff8_gm.dev_data   , sizeof(F_FLOAT*));
      cudaMemcpy(sdata->pair.coeff8_gm.dev_data, coeff8, n, cudaMemcpyHostToDevice);
      _coeff8_gm_tex.normalized = false;                      // access with normalized texture coordinates
      _coeff8_gm_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
      _coeff8_gm_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
      const textureReference* coeff8_gm_texture_ptr = &MY_AP(coeff8_gm_tex);

#if F_PRECISION == 1
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<float>();
      cudaBindTexture(0, coeff8_gm_texture_ptr, sdata->pair.coeff8_gm.dev_data, &channelDescXType, sdata->atom.nmax * sizeof(F_FLOAT));
#else
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<int2>();
      cudaBindTexture(0, coeff8_gm_texture_ptr, sdata->pair.coeff8_gm.dev_data, &channelDescXType, sdata->atom.nmax * 2 * sizeof(int2));
#endif
    }
  }

  CUT_CHECK_ERROR("Cuda_Pair: init Coeff8 failed");

  if(ncoeff > 8) {
    F_FLOAT coeff9[cuda_ntypes2];

    for(int i = 1; i <= sdata->atom.ntypes; ++i) {
      for(int j = 1; j <= sdata->atom.ntypes; ++j) {
        coeff9[i * cuda_ntypes + j] = (F_FLOAT) sdata->pair.coeff9[i][j];
      }
    }

    if(use_global_params) {
      cudaMemcpyToSymbol(MY_AP(coeff9_gm)  , &sdata->pair.coeff9_gm.dev_data   , sizeof(F_FLOAT*));
      cudaMemcpy(sdata->pair.coeff9_gm.dev_data, coeff9, n, cudaMemcpyHostToDevice);
      _coeff9_gm_tex.normalized = false;                      // access with normalized texture coordinates
      _coeff9_gm_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
      _coeff9_gm_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
      const textureReference* coeff9_gm_texture_ptr = &MY_AP(coeff9_gm_tex);

#if F_PRECISION == 1
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<float>();
      cudaBindTexture(0, coeff9_gm_texture_ptr, sdata->pair.coeff9_gm.dev_data, &channelDescXType, sdata->atom.nmax * sizeof(F_FLOAT));
#else
      cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<int2>();
      cudaBindTexture(0, coeff9_gm_texture_ptr, sdata->pair.coeff9_gm.dev_data, &channelDescXType, sdata->atom.nmax * 2 * sizeof(int2));
#endif
    }
  }

  CUT_CHECK_ERROR("Cuda_Pair: init Coeff9 failed");

  F_FLOAT special_lj[4];
  special_lj[0] = sdata->pair.special_lj[0];
  special_lj[1] = sdata->pair.special_lj[1];
  special_lj[2] = sdata->pair.special_lj[2];
  special_lj[3] = sdata->pair.special_lj[3];


  X_FLOAT box_size[3] = {
    sdata->domain.subhi[0] - sdata->domain.sublo[0],
    sdata->domain.subhi[1] - sdata->domain.sublo[1],
    sdata->domain.subhi[2] - sdata->domain.sublo[2]
  };

  cudaMemcpyToSymbol(MY_AP(box_size)   	, box_size                 		, sizeof(X_FLOAT) * 3);
  cudaMemcpyToSymbol(MY_AP(cuda_ntypes)	, &cuda_ntypes            		, sizeof(unsigned));
  cudaMemcpyToSymbol(MY_AP(special_lj) 	, special_lj               		, sizeof(F_FLOAT) * 4);
  cudaMemcpyToSymbol(MY_AP(virial)     	, &sdata->pair.virial.dev_data   , sizeof(ENERGY_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(eng_vdwl)     	, &sdata->pair.eng_vdwl.dev_data , sizeof(ENERGY_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(periodicity)	, sdata->domain.periodicity		, sizeof(int) * 3);
  cudaMemcpyToSymbol(MY_AP(collect_forces_later), &sdata->pair.collect_forces_later  , sizeof(int));

  if(need_q) {
    F_FLOAT qqrd2e_tmp = sdata->pppm.qqrd2e;
    F_FLOAT special_coul[4];
    special_coul[0] = sdata->pair.special_coul[0];
    special_coul[1] = sdata->pair.special_coul[1];
    special_coul[2] = sdata->pair.special_coul[2];
    special_coul[3] = sdata->pair.special_coul[3];

    cudaMemcpyToSymbol(MY_AP(special_coul)	, special_coul             		, sizeof(F_FLOAT) * 4);
    cudaMemcpyToSymbol(MY_AP(g_ewald)    	, &sdata->pair.g_ewald	   		, sizeof(F_FLOAT));
    cudaMemcpyToSymbol(MY_AP(qqrd2e)     	, &qqrd2e_tmp	   				, sizeof(F_FLOAT));
    cudaMemcpyToSymbol(MY_AP(kappa)     	, &sdata->pair.kappa				, sizeof(F_FLOAT));
    cudaMemcpyToSymbol(MY_AP(eng_coul)     , &sdata->pair.eng_coul.dev_data , sizeof(ENERGY_FLOAT*));
  }

  CUT_CHECK_ERROR("Cuda_Pair: init failed");
}
timespec startpairtime, endpairtime;
//Function which is called prior to kernel invocation, determins grid, Binds Textures, updates constant memory if necessary
void Cuda_Pair_PreKernel_AllStyles(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag, dim3 &grid, dim3 &threads, int &sharedperproc, bool need_q = false, int maxthreads = 256)
{
  if(sdata->atom.nlocal == 0) return;

  if(sdata->atom.update_neigh)
    Cuda_Pair_UpdateNeighbor_AllStyles(sdata, sneighlist);

  if(sdata->atom.update_nmax)
    Cuda_Pair_UpdateNmax_AllStyles(sdata, sneighlist);

  if(sdata->atom.update_nlocal) {
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
    cudaMemcpyToSymbol(MY_AP(nall)    , & sdata->atom.nall          , sizeof(int));
  }



  BindXTypeTexture(sdata);

  if(need_q) BindQTexture(sdata);


  sharedperproc = 0;

  if(sdata->pair.use_block_per_atom) sharedperproc += 3;

  if(eflag) sharedperproc += 1;

  if(need_q && eflag) sharedperproc += 1;

  if(vflag) sharedperproc += 6;

  int threadnum = sneighlist->inum;

  if(sdata->comm.comm_phase == 2)threadnum = sneighlist->inum_border2;

  if(sdata->pair.use_block_per_atom) {
    threadnum *= 64;
    maxthreads = 64;
  }

  int3 layout = getgrid(threadnum, sharedperproc * sizeof(ENERGY_FLOAT), maxthreads, true); //need to limit to 192 threads due to register limit
  threads.x = layout.z;
  threads.y = 1;
  threads.z = 1;
  grid.x = layout.x;
  grid.y = layout.y;
  grid.z = 1;

  int size = (unsigned)(layout.y * layout.x) * sharedperproc * sizeof(ENERGY_FLOAT);

  if(sdata->pair.collect_forces_later) size += (unsigned)(sdata->atom.nmax * 3 * sizeof(F_FLOAT));

  Cuda_UpdateBuffer(sdata, size);

  if(sdata->pair.use_block_per_atom)
    cudaMemset(sdata->buffer, 0, size);

  sdata->pair.lastgridsize = grid.x * grid.y;
  sdata->pair.n_energy_virial = sharedperproc;

  if(sdata->pair.use_block_per_atom) sdata->pair.n_energy_virial -= 3;

  clock_gettime(CLOCK_REALTIME, &startpairtime);

  MYDBG(printf("# CUDA: Cuda_Pair: kernel start eflag: %i vflag: %i config: %i %i %i %i\n", eflag, vflag, grid.x, grid.y, threads.x, sharedperproc * sizeof(ENERGY_FLOAT)*threads.x);)
}

//Function which is called after the kernel invocation, collects energy and virial
void Cuda_Pair_PostKernel_AllStyles(cuda_shared_data* sdata, dim3 &grid, int &sharedperproc, int eflag, int vflag)
{
  if((not sdata->pair.collect_forces_later) && (eflag || vflag)) { //not sdata->comm.comm_phase==2))
    cudaThreadSynchronize();
    clock_gettime(CLOCK_REALTIME, &endpairtime);
    sdata->cuda_timings.pair_kernel +=
      endpairtime.tv_sec - startpairtime.tv_sec + 1.0 * (endpairtime.tv_nsec - startpairtime.tv_nsec) / 1000000000;
    CUT_CHECK_ERROR("Cuda_Pair: Kernel execution failed");

    if(eflag || vflag) {
      int n = grid.x * grid.y;

      if(sdata->pair.use_block_per_atom)
        grid.x = sharedperproc - 3;
      else
        grid.x = sharedperproc;

      grid.y = 1;
      dim3 threads(128, 1, 1);
      MYDBG(printf("# CUDA: Cuda_Pair: virial compute kernel start eflag: %i vflag: %i config: %i %i %i %i\n", eflag, vflag, grid.x, grid.y, threads.x, sharedperproc * sizeof(ENERGY_FLOAT)*threads.x);)
      MY_AP(PairVirialCompute_reduce) <<< grid, threads, threads.x* sizeof(ENERGY_FLOAT)>>>(n);
      cudaThreadSynchronize();
      CUT_CHECK_ERROR("Cuda_Pair: virial compute Kernel execution failed");
    }

    MYDBG(printf("# CUDA: Cuda_Pair: kernel done\n");)
  }
}


#include "pair_born_coul_long_cuda.cu"
#include "pair_buck_coul_cut_cuda.cu"
#include "pair_buck_coul_long_cuda.cu"
#include "pair_buck_cuda.cu"
#include "pair_lj_sdk_cuda.cu"
#include "pair_lj_sdk_coul_cut_cuda.cu"
#include "pair_lj_sdk_coul_debye_cuda.cu"
#include "pair_lj_sdk_coul_long_cuda.cu"
#include "pair_gran_hooke_cuda.cu"
#include "pair_lj_charmm_coul_charmm_implicit_cuda.cu"
#include "pair_lj_charmm_coul_charmm_cuda.cu"
#include "pair_lj_charmm_coul_long_cuda.cu"
#include "pair_lj_class2_coul_cut_cuda.cu"
#include "pair_lj_class2_coul_long_cuda.cu"
#include "pair_lj_class2_cuda.cu"
#include "pair_lj_cut_coul_cut_cuda.cu"
#include "pair_lj_cut_coul_debye_cuda.cu"
#include "pair_lj_cut_coul_long_cuda.cu"
#include "pair_lj_cut_cuda.cu"
#include "pair_lj_cut_experimental_cuda.cu"
#include "pair_lj_expand_cuda.cu"
#include "pair_lj_gromacs_cuda.cu"
#include "pair_lj_gromacs_coul_gromacs_cuda.cu"
#include "pair_lj_smooth_cuda.cu"
#include "pair_lj96_cut_cuda.cu"
#include "pair_morse_coul_long_cuda.cu"
#include "pair_morse_cuda.cu"
#include "pair_eam_cuda.cu"

#include "cuda_pair_kernel.cu"

#include "pair_manybody_const.h"
#include "pair_tersoff_cuda.cu"
#include "pair_sw_cuda.cu"

void Cuda_Pair_UpdateNmax(cuda_shared_data* sdata)
{
  CUT_CHECK_ERROR("Cuda_Pair: before updateNmax failed");
  cudaMemcpyToSymbol(MY_AP(nlocal)    , & sdata->atom.nlocal             , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nall)      , & sdata->atom.nall               , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nmax)      , & sdata->atom.nmax               , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(type)      , & sdata->atom.type       .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(x)         , & sdata->atom.x          .dev_data, sizeof(X_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(x_type)    , & sdata->atom.x_type     .dev_data, sizeof(X_FLOAT4*));
  cudaMemcpyToSymbol(MY_AP(xhold)     , & sdata->atom.xhold      .dev_data, sizeof(X_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(v)         , & sdata->atom.v          .dev_data, sizeof(V_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(radius)    , & sdata->atom.radius     .dev_data, sizeof(X_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(v_radius)  , & sdata->atom.v_radius   .dev_data, sizeof(V_FLOAT4*));
  cudaMemcpyToSymbol(MY_AP(omega)     , & sdata->atom.omega      .dev_data, sizeof(V_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(rmass)     , & sdata->atom.rmass      .dev_data, sizeof(V_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(omega_rmass), & sdata->atom.omega_rmass.dev_data, sizeof(V_FLOAT4*));
  cudaMemcpyToSymbol(MY_AP(map_array), & sdata->atom.map_array .dev_data, sizeof(int*));
  CUT_CHECK_ERROR("Cuda_Pair: updateNmax failed");
}


void Cuda_Pair_GenerateXType(cuda_shared_data* sdata)
{
  MYDBG(printf(" # CUDA: GenerateXType ... start %i %i %i %p %p %p %p\n", sdata->atom.nlocal, sdata->atom.nall, sdata->atom.nmax, sdata->atom.x.dev_data, sdata->atom.x_type.dev_data, sdata->atom.xhold.dev_data, sdata->atom.type.dev_data);)

  if(sdata->atom.update_nmax)
    Cuda_Pair_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal) {
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
    cudaMemcpyToSymbol(MY_AP(nall)    , & sdata->atom.nall          , sizeof(int));
  }

  MYDBG(printf(" # CUDA: GenerateXType ... getgrid\n"); fflush(stdout);)

  int3 layout = getgrid(sdata->atom.nall);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  MYDBG(printf(" # CUDA: GenerateXType ... kernel start test\n");  fflush(stdout);)
  Pair_GenerateXType_Kernel <<< grid, threads, 0>>>();
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_Pair GenerateXType: Kernel failed");
  MYDBG(printf(" # CUDA: GenerateXType ... end\n");  fflush(stdout);)
}

void Cuda_Pair_RevertXType(cuda_shared_data* sdata)
{
  MYDBG(printf(" # CUDA: RevertXType ... start\n");)

  if(sdata->atom.update_nmax)
    Cuda_Pair_UpdateNmax(sdata);

  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nall)    , & sdata->atom.nall          , sizeof(int));

  int3 layout = getgrid(sdata->atom.nall);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  Pair_RevertXType_Kernel <<< grid, threads, 0>>>();
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_Pair GenerateXType: Kernel failed");
  MYDBG(printf(" # CUDA: RevertXType ... end\n");)
}

void Cuda_Pair_GenerateVRadius(cuda_shared_data* sdata)
{
  MYDBG(printf(" # CUDA: GenerateVRadius ... start %i %i %i %p %p %p %p\n", sdata->atom.nlocal, sdata->atom.nall, sdata->atom.nmax, sdata->atom.x.dev_data, sdata->atom.x_type.dev_data, sdata->atom.xhold.dev_data, sdata->atom.type.dev_data);)

  if(sdata->atom.update_nmax)
    Cuda_Pair_UpdateNmax(sdata);

  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nall)    , & sdata->atom.nall          , sizeof(int));
  MYDBG(printf(" # CUDA: GenerateVRadius ... getgrid\n"); fflush(stdout);)

  int3 layout = getgrid(sdata->atom.nall);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  MYDBG(printf(" # CUDA: GenerateVRadius ... kernel start test\n");  fflush(stdout);)
  Pair_GenerateVRadius_Kernel <<< grid, threads, 0>>>();
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_Pair GenerateVRadius: Kernel failed");
  MYDBG(printf(" # CUDA: GenerateVRadius ... end\n");  fflush(stdout);)
}

void Cuda_Pair_GenerateOmegaRmass(cuda_shared_data* sdata)
{
  MYDBG(printf(" # CUDA: GenerateOmegaRmass ... start %i %i %i %p %p %p %p\n", sdata->atom.nlocal, sdata->atom.nall, sdata->atom.nmax, sdata->atom.x.dev_data, sdata->atom.x_type.dev_data, sdata->atom.xhold.dev_data, sdata->atom.type.dev_data);)

  if(sdata->atom.update_nmax)
    Cuda_Pair_UpdateNmax(sdata);

  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nall)    , & sdata->atom.nall          , sizeof(int));
  MYDBG(printf(" # CUDA: GenerateOmegaRmass ... getgrid\n"); fflush(stdout);)

  int3 layout = getgrid(sdata->atom.nall);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  MYDBG(printf(" # CUDA: GenerateOmegaRmass ... kernel start test\n");  fflush(stdout);)
  Pair_GenerateOmegaRmass_Kernel <<< grid, threads, 0>>>();
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_Pair GenerateOmegaRmass: Kernel failed");
  MYDBG(printf(" # CUDA: GenerateOmegaRmass ... end\n");  fflush(stdout);)
}

void Cuda_Pair_BuildXHold(cuda_shared_data* sdata)
{
  if(sdata->atom.update_nmax)
    Cuda_Pair_UpdateNmax(sdata);

  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nall)    , & sdata->atom.nall          , sizeof(int));

  int3 layout = getgrid(sdata->atom.nall);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  Pair_BuildXHold_Kernel <<< grid, threads, 0>>>();
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_Pair GenerateXType: Kernel failed");
}

void Cuda_Pair_CollectForces(cuda_shared_data* sdata, int eflag, int vflag)
{
  cudaThreadSynchronize();
  clock_gettime(CLOCK_REALTIME, &endpairtime);
  sdata->cuda_timings.pair_kernel +=
    endpairtime.tv_sec - startpairtime.tv_sec + 1.0 * (endpairtime.tv_nsec - startpairtime.tv_nsec) / 1000000000;
  CUT_CHECK_ERROR("Cuda_Pair: Kernel execution failed");
  dim3 threads;
  dim3 grid;

  if(eflag || vflag) {
    int n = sdata->pair.lastgridsize;
    grid.x = sdata->pair.n_energy_virial;
    grid.y = 1;
    threads.x = 128;
    //printf("A grid.x: %i\n",grid.x);
    MY_AP(PairVirialCompute_reduce) <<< grid, threads, threads.x* sizeof(ENERGY_FLOAT)>>>(n);
    cudaThreadSynchronize();
    CUT_CHECK_ERROR("Cuda_Pair_CollectForces: virial compute Kernel execution failed");
  }

  int3 layout = getgrid(sdata->atom.nlocal);
  threads.x = layout.z;
  grid.x = layout.x;
  grid.y = layout.y;
  Pair_CollectForces_Kernel <<< grid, threads, 0>>>(sdata->pair.n_energy_virial, sdata->pair.lastgridsize);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_Pair_CollectForces: Force Summation Kernel execution failed");

}
