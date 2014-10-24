// **************************************************************************
//                                   dpd.cu
//                             -------------------
//                           Trung Dac Nguyen (ORNL)
//
//  Device code for acceleration of the dpd pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : Jan 15, 2014
//    email                : nguyentd@ornl.gov
// ***************************************************************************/

#ifdef NV_KERNEL
#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
texture<float4> pos_tex;
texture<float4> vel_tex;
#else
texture<int4,1> pos_tex;
texture<int4,1> vel_tex;
#endif
#else
#define pos_tex x_
#define vel_tex v_
#endif

#define EPSILON (numtyp)1.0e-10

//#define _USE_UNIFORM_SARU_LCG
//#define _USE_UNIFORM_SARU_TEA8
//#define _USE_GAUSSIAN_SARU_LCG

#if !defined(_USE_UNIFORM_SARU_LCG) && !defined(_USE_UNIFORM_SARU_TEA8) && !defined(_USE_GAUSSIAN_SARU_LCG)
#define _USE_UNIFORM_SARU_LCG
#endif

// References: 
// 1. Y. Afshar, F. Schmid, A. Pishevar, S. Worley, Comput. Phys. Comm. 184 (2013), 1119–1128.
// 2. C. L. Phillips, J. A. Anderson, S. C. Glotzer, Comput. Phys. Comm. 230 (2011), 7191-7201.
// PRNG period = 3666320093*2^32 ~ 2^64 ~ 10^19

#define LCGA 0x4beb5d59 // Full period 32 bit LCG
#define LCGC 0x2600e1f7
#define oWeylPeriod 0xda879add // Prime period 3666320093
#define oWeylOffset 0x8009d14b
#define TWO_N32 0.232830643653869628906250e-9f /* 2^-32 */

// specifically implemented for steps = 1; high = 1.0; low = -1.0
// returns uniformly distributed random numbers u in [-1.0;1.0] 
// using the inherent LCG, then multiply u with sqrt(3) to "match" 
// with a normal random distribution. 
// Afshar et al. mutlplies u in [-0.5;0.5] with sqrt(12)
// Curly brackets to make variables local to the scope.
#ifdef _USE_UNIFORM_SARU_LCG
#define SQRT3 (numtyp)1.7320508075688772935274463
#define saru(seed1, seed2, seed, timestep, randnum) {                         \
  unsigned int seed3 = seed + timestep;                                       \
  seed3^=(seed1<<7)^(seed2>>6);                                               \
  seed2+=(seed1>>4)^(seed3>>15);                                              \
  seed1^=(seed2<<9)+(seed3<<8);                                               \
  seed3^=0xA5366B4D*((seed2>>11) ^ (seed1<<1));                               \
  seed2+=0x72BE1579*((seed1<<4)  ^ (seed3>>16));                              \
  seed1^=0x3F38A6ED*((seed3>>5)  ^ (((signed int)seed2)>>22));                \
  seed2+=seed1*seed3;                                                         \
  seed1+=seed3 ^ (seed2>>2);                                                  \
  seed2^=((signed int)seed2)>>17;                                             \
  unsigned int state  = 0x79dedea3*(seed1^(((signed int)seed1)>>14));         \
  unsigned int wstate = (state + seed2) ^ (((signed int)state)>>8);           \
  state  = state + (wstate*(wstate^0xdddf97f5));                              \
  wstate = 0xABCB96F7 + (wstate>>1);                                          \
  state = LCGA*state + LCGC;                                                  \
  wstate = wstate + oWeylOffset+((((signed int)wstate)>>31) & oWeylPeriod);   \
  unsigned int v = (state ^ (state>>26)) + wstate;                            \
  unsigned int s = (signed int)((v^(v>>20))*0x6957f5a7);                      \
  randnum = SQRT3*(s*TWO_N32*(numtyp)2.0-(numtyp)1.0);                        \
}
#endif

// specifically implemented for steps = 1; high = 1.0; low = -1.0
// returns uniformly distributed random numbers u in [-1.0;1.0] using TEA8 
// then multiply u with sqrt(3) to "match" with a normal random distribution 
// Afshar et al. mutlplies u in [-0.5;0.5] with sqrt(12)
#ifdef _USE_UNIFORM_SARU_TEA8
#define SQRT3 (numtyp)1.7320508075688772935274463
#define k0 0xA341316C
#define k1 0xC8013EA4
#define k2 0xAD90777D
#define k3 0x7E95761E
#define delta 0x9e3779b9
#define rounds 8
#define saru(seed1, seed2, seed, timestep, randnum) {                         \
  unsigned int seed3 = seed + timestep;                                       \
  seed3^=(seed1<<7)^(seed2>>6);                                               \
  seed2+=(seed1>>4)^(seed3>>15);                                              \
  seed1^=(seed2<<9)+(seed3<<8);                                               \
  seed3^=0xA5366B4D*((seed2>>11) ^ (seed1<<1));                               \
  seed2+=0x72BE1579*((seed1<<4)  ^ (seed3>>16));                              \
  seed1^=0x3F38A6ED*((seed3>>5)  ^ (((signed int)seed2)>>22));                \
  seed2+=seed1*seed3;                                                         \
  seed1+=seed3 ^ (seed2>>2);                                                  \
  seed2^=((signed int)seed2)>>17;                                             \
  unsigned int state  = 0x79dedea3*(seed1^(((signed int)seed1)>>14));         \
  unsigned int wstate = (state + seed2) ^ (((signed int)state)>>8);           \
  state  = state + (wstate*(wstate^0xdddf97f5));                              \
  wstate = 0xABCB96F7 + (wstate>>1);                                          \
  unsigned int sum = 0;                                                       \
  for (int i=0; i < rounds; i++) {                                            \
    sum += delta;                                                             \
    state += ((wstate<<4) + k0)^(wstate + sum)^((wstate>>5) + k1);            \
    wstate += ((state<<4) + k2)^(state + sum)^((state>>5) + k3);              \
  }                                                                           \
  unsigned int v = (state ^ (state>>26)) + wstate;                            \
  unsigned int s = (signed int)((v^(v>>20))*0x6957f5a7);                      \
  randnum = SQRT3*(s*TWO_N32*(numtyp)2.0-(numtyp)1.0);                        \
}
#endif

// specifically implemented for steps = 1; high = 1.0; low = -1.0
// returns two uniformly distributed random numbers r1 and r2 in [-1.0;1.0], 
// and uses the polar method (Marsaglia's) to transform to a normal random value
// This is used to compared with CPU DPD using RandMars::gaussian()
#ifdef _USE_GAUSSIAN_SARU_LCG
#define saru(seed1, seed2, seed, timestep, randnum) {                         \
  unsigned int seed3 = seed + timestep;                                       \
  seed3^=(seed1<<7)^(seed2>>6);                                               \
  seed2+=(seed1>>4)^(seed3>>15);                                              \
  seed1^=(seed2<<9)+(seed3<<8);                                               \
  seed3^=0xA5366B4D*((seed2>>11) ^ (seed1<<1));                               \
  seed2+=0x72BE1579*((seed1<<4)  ^ (seed3>>16));                              \
  seed1^=0x3F38A6ED*((seed3>>5)  ^ (((signed int)seed2)>>22));                \
  seed2+=seed1*seed3;                                                         \
  seed1+=seed3 ^ (seed2>>2);                                                  \
  seed2^=((signed int)seed2)>>17;                                             \
  unsigned int state=0x12345678;                                              \
  unsigned int wstate=12345678;                                               \
  state  = 0x79dedea3*(seed1^(((signed int)seed1)>>14));                      \
  wstate = (state + seed2) ^ (((signed int)state)>>8);                        \
  state  = state + (wstate*(wstate^0xdddf97f5));                              \
  wstate = 0xABCB96F7 + (wstate>>1);                                          \
  unsigned int v, s;                                                          \
  numtyp r1, r2, rsq;                                                         \
  while (1) {                                                                 \
    state = LCGA*state + LCGC;                                                \
    wstate = wstate + oWeylOffset+((((signed int)wstate)>>31) & oWeylPeriod); \
    v = (state ^ (state>>26)) + wstate;                                       \
    s = (signed int)((v^(v>>20))*0x6957f5a7);                                 \
    r1 = s*TWO_N32*(numtyp)2.0-(numtyp)1.0;                                   \
    state = LCGA*state + LCGC;                                                \
    wstate = wstate + oWeylOffset+((((signed int)wstate)>>31) & oWeylPeriod); \
    v = (state ^ (state>>26)) + wstate;                                       \
    s = (signed int)((v^(v>>20))*0x6957f5a7);                                 \
    r2 = s*TWO_N32*(numtyp)2.0-(numtyp)1.0;                                   \
    rsq = r1 * r1 + r2 * r2;                                                  \
    if (rsq < (numtyp)1.0) break;                                             \
  }                                                                           \
  numtyp fac = ucl_sqrt((numtyp)-2.0*log(rsq)/rsq);                           \
  randnum = r2*fac;                                                           \
}
#endif
                                                                             
__kernel void k_dpd(const __global numtyp4 *restrict x_, 
                    const __global numtyp4 *restrict coeff,
                    const int lj_types, 
                    const __global numtyp *restrict sp_lj, 
                    const __global int * dev_nbor, 
                    const __global int * dev_packed, 
                    __global acctyp4 *restrict ans, 
                    __global acctyp *restrict engv, 
                    const int eflag, const int vflag, const int inum,
                    const int nbor_pitch, 
                    const __global numtyp4 *restrict v_,
                    const __global numtyp *restrict cutsq,
                    const numtyp dtinvsqrt, const int seed, 
                    const int timestep, const int tstat_only,
                    const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  if (ii<inum) {
    int i, numj, nbor, nbor_end;
    __local int n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);
  
    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    numtyp4 iv; fetch4(iv,i,vel_tex); //v_[i];
    int itag=iv.w;

    numtyp factor_dpd;
    for ( ; nbor<nbor_end; nbor+=n_stride) {
  
      int j=dev_packed[nbor];
      factor_dpd = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      numtyp4 jv; fetch4(jv,j,vel_tex); //v_[j];
      int jtag=jv.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
        
      int mtype=itype*lj_types+jtype;
      if (rsq<cutsq[mtype]) {
        numtyp r=ucl_sqrt(rsq);
        if (r < EPSILON) continue;

        numtyp rinv=ucl_recip(r);
        numtyp delvx = iv.x - jv.x;
        numtyp delvy = iv.y - jv.y;
        numtyp delvz = iv.z - jv.z;
        numtyp dot = delx*delvx + dely*delvy + delz*delvz;
        numtyp wd = (numtyp)1.0 - r/coeff[mtype].w;

        unsigned int tag1=itag, tag2=jtag;
        if (tag1 > tag2) {
          tag1 = jtag; tag2 = itag;
        }
        
        numtyp randnum = (numtyp)0.0;  
        saru(tag1, tag2, seed, timestep, randnum);

        // conservative force = a0 * wd, or 0 if tstat only
        // drag force = -gamma * wd^2 * (delx dot delv) / r
        // random force = sigma * wd * rnd * dtinvsqrt;

        numtyp force = (numtyp)0.0;
        if (!tstat_only) force = coeff[mtype].x*wd;
        force -= coeff[mtype].y*wd*wd*dot*rinv;
        force += coeff[mtype].z*wd*randnum*dtinvsqrt;
        force*=factor_dpd*rinv;
      
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          // unshifted eng of conservative term:
          // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
          // eng shifted to 0.0 at cutoff
          numtyp e = (numtyp)0.5*coeff[mtype].x*coeff[mtype].w * wd*wd;
          energy+=factor_dpd*e; 
        }
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

__kernel void k_dpd_fast(const __global numtyp4 *restrict x_, 
                         const __global numtyp4 *restrict coeff_in,
                         const __global numtyp *restrict sp_lj_in, 
                         const __global int * dev_nbor, 
                         const __global int * dev_packed, 
                         __global acctyp4 *restrict ans, 
                         __global acctyp *restrict engv, 
                         const int eflag, const int vflag, const int inum, 
                         const int nbor_pitch,
                         const __global numtyp4 *restrict v_,
                         const __global numtyp *restrict cutsq,
                         const numtyp dtinvsqrt, const int seed, 
                         const int timestep, const int tstat_only,
                         const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  __local numtyp4 coeff[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  if (tid<4)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    coeff[tid]=coeff_in[tid];
  }
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __syncthreads();
  
  if (ii<inum) {
    int i, numj, nbor, nbor_end;
    __local int n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);
    numtyp4 iv; fetch4(iv,i,vel_tex); //v_[i];
    int itag=iv.w;

    numtyp factor_dpd;
    for ( ; nbor<nbor_end; nbor+=n_stride) {
  
      int j=dev_packed[nbor];
      factor_dpd = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype+jx.w;
      numtyp4 jv; fetch4(jv,j,vel_tex); //v_[j];
      int jtag=jv.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
        
      if (rsq<cutsq[mtype]) {
        numtyp r=ucl_sqrt(rsq);
        if (r < EPSILON) continue;

        numtyp rinv=ucl_recip(r);
        numtyp delvx = iv.x - jv.x;
        numtyp delvy = iv.y - jv.y;
        numtyp delvz = iv.z - jv.z;
        numtyp dot = delx*delvx + dely*delvy + delz*delvz;
        numtyp wd = (numtyp)1.0 - r/coeff[mtype].w;

        unsigned int tag1=itag, tag2=jtag;
        if (tag1 > tag2) {
          tag1 = jtag; tag2 = itag;
        }
        
        numtyp randnum = (numtyp)0.0;  
        saru(tag1, tag2, seed, timestep, randnum);

        // conservative force = a0 * wd, or 0 if tstat only
        // drag force = -gamma * wd^2 * (delx dot delv) / r
        // random force = sigma * wd * rnd * dtinvsqrt;

        numtyp force = (numtyp)0.0;
        if (!tstat_only) force = coeff[mtype].x*wd;
        force -= coeff[mtype].y*wd*wd*dot*rinv;
        force += coeff[mtype].z*wd*randnum*dtinvsqrt;
        force*=factor_dpd*rinv;
      
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          // unshifted eng of conservative term:
          // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
          // eng shifted to 0.0 at cutoff
          numtyp e = (numtyp)0.5*coeff[mtype].x*coeff[mtype].w * wd*wd;
          energy+=factor_dpd*e; 
        }
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

