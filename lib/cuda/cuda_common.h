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

#ifndef _CUDA_COMMON_H_
#define _CUDA_COMMON_H_

//#include "cutil.h"
#include "cuda_precision.h"
#include "cuda_wrapper_cu.h"

#define CUDA_MAX_TYPES_PLUS_ONE 12 //for pair styles which use constant space for parameters, this needs to be one larger than the number of atom types
//this can not be arbitrarly large, since constant space is limited.
//in principle one could alter potentials to use global memory for parameters, some du that already since the first examples I encountered had a high number (20+) of atom types
//Christian
#define CUDA_MAX_TYPES2 (CUDA_MAX_TYPES_PLUS_ONE * CUDA_MAX_TYPES_PLUS_ONE)
#define CUDA_MAX_NSPECIAL 25

// define some easy-to-use debug and emulation macros
#ifdef _DEBUG
#define MYDBG(a) a
#else
#define MYDBG(a)
#endif

#if __DEVICE_EMULATION__
#define MYEMU(a) a
#else
#define MYEMU(a)
#endif

#define MYEMUDBG(a) MYEMU(MYDBG(a))

// Add Prefix (needed as workaround, same constant's names in different files causes conflict)
#define MY_ADD_PREFIX(prefix, var) prefix##_##var
#define MY_ADD_PREFIX2(prefix, var) MY_ADD_PREFIX(prefix, var)
#define MY_AP(var) MY_ADD_PREFIX2(MY_PREFIX, var)

#define MY_VAR_TO_STR(var) #var
#define MY_VAR_TO_STR2(var) MY_VAR_TO_STR(var)
//#define &MY_AP(var) (MY_VAR_TO_STR2(MY_PREFIX) "_" MY_VAR_TO_STR2(var))
//#define &MY_AP(var) &(MY_AP(var))
#define CUDA_USE_TEXTURE
#define CUDA_USE_CFLOAT4

//constants used by many classes

//domain
#define _boxhi       MY_AP(boxhi)
#define _boxlo       MY_AP(boxlo)
#define _subhi       MY_AP(subhi)
#define _sublo       MY_AP(sublo)
#define _box_size    MY_AP(box_size)
#define _prd         MY_AP(prd)
#define _periodicity MY_AP(periodicity)
#define _triclinic	 MY_AP(triclinic)
#define _boxhi_lamda MY_AP(boxhi_lamda)
#define _boxlo_lamda MY_AP(boxlo_lamda)
#define _prd_lamda   MY_AP(prd_lamda)
#define _h		 	 MY_AP(h)
#define _h_inv	 	 MY_AP(h_inv)
#define _h_rate		 MY_AP(h_rate)
__device__ __constant__ X_CFLOAT _boxhi[3];
__device__ __constant__ X_CFLOAT _boxlo[3];
__device__ __constant__ X_CFLOAT _subhi[3];
__device__ __constant__ X_CFLOAT _sublo[3];
__device__ __constant__ X_CFLOAT _box_size[3];
__device__ __constant__ X_CFLOAT _prd[3];
__device__ __constant__ int _periodicity[3];
__device__ __constant__ int _triclinic;
__device__ __constant__ X_CFLOAT _boxhi_lamda[3];
__device__ __constant__ X_CFLOAT _boxlo_lamda[3];
__device__ __constant__ X_CFLOAT _prd_lamda[3];
__device__ __constant__ X_CFLOAT _h[6];
__device__ __constant__ X_CFLOAT _h_inv[6];
__device__ __constant__ V_CFLOAT _h_rate[6];


//atom properties
#define _x           MY_AP(x)
#define _v           MY_AP(v)
#define _f           MY_AP(f)
#define _tag         MY_AP(tag)
#define _type        MY_AP(type)
#define _mask        MY_AP(mask)
#define _image       MY_AP(image)
#define _q           MY_AP(q)
#define _mass        MY_AP(mass)
#define _rmass       MY_AP(rmass)
#define _rmass_flag  MY_AP(rmass_flag)
#define _eatom       MY_AP(eatom)
#define _vatom       MY_AP(vatom)
#define _x_type      MY_AP(x_type)
#define _radius      MY_AP(radius)
#define _density     MY_AP(density)
#define _omega       MY_AP(omega)
#define _torque      MY_AP(torque)
#define _special     MY_AP(special)
#define _maxspecial  MY_AP(maxspecial)
#define _nspecial    MY_AP(nspecial)
#define _special_flag  MY_AP(special_flag)
#define _molecule    MY_AP(molecule)
#define _v_radius    MY_AP(v_radius)
#define _omega_rmass MY_AP(omega_rmass)
#define _freeze_group_bit MY_AP(freeze_group_bit)
#define _map_array   MY_AP(map_array)
__device__ __constant__ X_CFLOAT* _x;  //holds pointer to positions
__device__ __constant__ V_CFLOAT* _v;
__device__ __constant__ F_CFLOAT* _f;
__device__ __constant__ int* _tag;
__device__ __constant__ int* _type;
__device__ __constant__ int* _mask;
__device__ __constant__ int* _image;
__device__ __constant__ V_CFLOAT* _mass;
__device__ __constant__ F_CFLOAT* _q;
__device__ __constant__ V_CFLOAT* _rmass;
__device__ __constant__ int _rmass_flag;
__device__ __constant__ ENERGY_CFLOAT* _eatom;
__device__ __constant__ ENERGY_CFLOAT* _vatom;
__device__ __constant__ X_CFLOAT4* _x_type;  //holds pointer to positions
__device__ __constant__ X_CFLOAT* _radius;
__device__ __constant__ F_CFLOAT* _density;
__device__ __constant__ V_CFLOAT* _omega;
__device__ __constant__ F_CFLOAT* _torque;
__device__ __constant__ int* _special;
__device__ __constant__ int _maxspecial;
__device__ __constant__ int* _nspecial;
__device__ __constant__ int _special_flag[4];
__device__ __constant__ int* _molecule;
__device__ __constant__ V_CFLOAT4* _v_radius;  //holds pointer to positions
__device__ __constant__ V_CFLOAT4* _omega_rmass;  //holds pointer to positions
__device__ __constant__ int _freeze_group_bit;
__device__ __constant__ int* _map_array;

#ifdef CUDA_USE_TEXTURE

#define _x_tex         MY_AP(x_tex)
#if X_PRECISION == 1
texture<float> _x_tex;
#else
texture<int2, 1> _x_tex;
#endif

#define _type_tex         MY_AP(type_tex)
texture<int> _type_tex;

#define _x_type_tex         MY_AP(x_type_tex)
#if X_PRECISION == 1
texture<float4, 1> _x_type_tex;
#else
texture<int4, 1> _x_type_tex;
#endif

#define _v_radius_tex         MY_AP(v_radius_tex)
#if V_PRECISION == 1
texture<float4, 1> _v_radius_tex;
#else
texture<int4, 1> _v_radius_tex;
#endif

#define _omega_rmass_tex         MY_AP(omega_rmass_tex)
#if V_PRECISION == 1
texture<float4, 1> _omega_rmass_tex;
#else
texture<int4, 1> _omega_rmass_tex;
#endif

#define _q_tex         MY_AP(q_tex)
#if F_PRECISION == 1
texture<float> _q_tex;
#else
texture<int2, 1> _q_tex;
#endif

#endif

//neighbor
#ifdef IncludeCommonNeigh
#define _inum        	MY_AP(inum)
#define _inum_border    MY_AP(inum_border)
#define _ilist       	MY_AP(ilist)
#define _ilist_border 	MY_AP(ilist_border)
#define _numneigh    	MY_AP(numneigh)
#define _numneigh_border 	MY_AP(numneigh_border)
#define _numneigh_inner		MY_AP(numneigh_inner)
#define _firstneigh  	MY_AP(firstneigh)
#define _neighbors 	MY_AP(neighbors)
#define _neighbors_border 	MY_AP(neighbors_border)
#define _neighbors_inner  	MY_AP(neighbors_inner)
#define _reneigh_flag 	MY_AP(reneigh_flag)
#define _triggerneighsq MY_AP(triggerneighsq)
#define _xhold       	MY_AP(xhold)
#define _maxhold     	MY_AP(maxhold)
#define _dist_check     MY_AP(dist_check)
#define _neighbor_maxlocal MY_AP(neighbor_maxlocal)
#define _maxneighbors   MY_AP(maxneighbors)
#define _overlap_comm   MY_AP(overlap_comm)
__device__ __constant__ int _inum;
__device__ __constant__ int* _inum_border;
__device__ __constant__ int* _ilist;
__device__ __constant__ int* _ilist_border;
__device__ __constant__ int* _numneigh;
__device__ __constant__ int* _numneigh_border;
__device__ __constant__ int* _numneigh_inner;
__device__ __constant__ int** _firstneigh;
__device__ __constant__ int* _neighbors;
__device__ __constant__ int* _neighbors_border;
__device__ __constant__ int* _neighbors_inner;
__device__ __constant__ int* _reneigh_flag;
__device__ __constant__ X_CFLOAT _triggerneighsq;
__device__ __constant__ X_CFLOAT* _xhold;  //holds pointer to positions
__device__ __constant__ int _maxhold;
__device__ __constant__ int _dist_check;
__device__ __constant__ int _neighbor_maxlocal;
__device__ __constant__ int _maxneighbors;
__device__ __constant__ int _overlap_comm;
#endif

//system properties
#define _nall        MY_AP(nall)
#define _nghost      MY_AP(nghost)
#define _nlocal      MY_AP(nlocal)
#define _nmax        MY_AP(nmax)
#define _cuda_ntypes MY_AP(cuda_ntypes)
#define _dtf         MY_AP(dtf)
#define _dtv         MY_AP(dtv)
#define _factor      MY_AP(factor)
#define _virial      MY_AP(virial)
#define _eng_vdwl    MY_AP(eng_vdwl)
#define _eng_coul    MY_AP(eng_coul)
#define _molecular   MY_AP(molecular)
__device__ __constant__ unsigned _nall;
__device__ __constant__ unsigned _nghost;
__device__ __constant__ unsigned _nlocal;
__device__ __constant__ unsigned _nmax;
__device__ __constant__ unsigned _cuda_ntypes;
__device__ __constant__ V_CFLOAT _dtf;
__device__ __constant__ X_CFLOAT _dtv;
__device__ __constant__ V_CFLOAT _factor;
__device__ __constant__ ENERGY_CFLOAT* _virial;
__device__ __constant__ ENERGY_CFLOAT* _eng_vdwl;
__device__ __constant__ ENERGY_CFLOAT* _eng_coul;
__device__ __constant__ int _molecular;

//other general constants
#define _buffer      MY_AP(buffer)
#define _flag		 MY_AP(flag)
#define _debugdata   MY_AP(debugdata)
__device__ __constant__ void* _buffer;
__device__ __constant__ int* _flag;
__device__ __constant__ int* _debugdata;

// pointers to data fields on GPU are hold in constant space
// -> reduces register usage and number of parameters for kernelcalls
// will be variables of file scope in cuda files




// maybe used to output cudaError_t
#define MY_OUTPUT_RESULT(result) \
  switch(result) \
  { \
  case cudaSuccess: printf(" => cudaSuccess\n"); break; \
  case cudaErrorInvalidValue: printf(" => cudaErrorInvalidValue\n"); break; \
  case cudaErrorInvalidSymbol: printf(" => cudaErrorInvalidSymbol\n"); break; \
  case cudaErrorInvalidDevicePointer: printf(" => cudaErrorInvalidDevicePointer\n"); break; \
  case cudaErrorInvalidMemcpyDirection: printf(" => cudaErrorInvalidMemcpyDirection\n"); break; \
  default: printf(" => unknown\n"); break; \
  }

#ifdef _DEBUG
#  define CUT_CHECK_ERROR(errorMessage) {                                    \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
      fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
              errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
      exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    err = cudaThreadSynchronize();                                           \
    if( cudaSuccess != err) {                                                \
      fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
              errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
      exit(EXIT_FAILURE);                                                  \
    }                                                                        \
  }
#else
#  define CUT_CHECK_ERROR(errorMessage) {                                    \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
      fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
              errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
      exit(EXIT_FAILURE);                                                  \
    }                                                                        \
  }
#endif

#  define CUDA_SAFE_CALL_NO_SYNC( call) {                                    \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
              __FILE__, __LINE__, cudaGetErrorString( err) );              \
      exit(EXIT_FAILURE);                                                  \
    } }

#  define CUDA_SAFE_CALL( call)     CUDA_SAFE_CALL_NO_SYNC(call);

#define X_MASK 1
#define V_MASK 2
#define F_MASK 4
#define TAG_MASK 8
#define TYPE_MASK 16
#define MASK_MASK 32
#define IMAGE_MASK 64
#define Q_MASK 128
#define MOLECULE_MASK 256
#define RMASS_MASK 512
#define RADIUS_MASK 1024
#define DENSITY_MASK 2048
#define OMEGA_MASK 4096
#define TORQUE_MASK 8192



#endif // #ifdef _CUDA_COMMON_H_
