// **************************************************************************
//                                   mdpd.cu
//                             -------------------
//                           Trung Dac Nguyen (ORNL)
//
//  Device code for acceleration of the mdpd pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : December 2023
//    email                : ndactrung@gmail.com
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
_texture( vel_tex,float4);
#else
_texture_2d( pos_tex,int4);
_texture_2d( vel_tex,int4);
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

#define LCGA 0x4beb5d59 /* Full period 32 bit LCG */
#define LCGC 0x2600e1f7
#define oWeylPeriod 0xda879add /* Prime period 3666320093 */
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

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) < (B) ? (B) : (A))

// coeff.x = A_att, coeff.y = B_rep, coeff.z = gamma, coeff.w = sigma
// coeff2.x = cut, coeff2.y = cut_r, coeff2.z = cutsq

__kernel void k_mdpd(const __global numtyp4 *restrict x_,
                    const __global numtyp4 *restrict extra,
                    const __global numtyp4 *restrict coeff,
                    const __global numtyp4 *restrict coeff2,
                    const int lj_types,
                    const __global numtyp *restrict sp_lj,
                    const __global numtyp *restrict sp_sqrt,
                    const __global int * dev_nbor,
                    const __global int * dev_packed,
                    __global acctyp3 *restrict ans,
                    __global acctyp *restrict engv,
                    const int eflag, const int vflag, const int inum,
                    const int nbor_pitch,
                    const __global numtyp4 *restrict v_,
                    const __global numtyp *restrict cutsq,
                    const numtyp dtinvsqrt, const int seed,
                    const int timestep, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  int n_stride;
  local_allocate_store_pair();

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int i, numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    numtyp4 iv; fetch4(iv,i,vel_tex); //v_[i];
    int itag=iv.w;

    const numtyp rhoi = extra[i].x;

    numtyp factor_dpd, factor_sqrt;
    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
      factor_dpd = sp_lj[sbmask(j)];
      factor_sqrt = sp_sqrt[sbmask(j)];
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
      if (rsq<coeff2[mtype].z) {  // cutsq[itype][jtype]
        numtyp r=ucl_sqrt(rsq);
        if (r < EPSILON) continue;

        numtyp rinv=ucl_recip(r);
        numtyp delvx = iv.x - jv.x;
        numtyp delvy = iv.y - jv.y;
        numtyp delvz = iv.z - jv.z;
        numtyp dot = delx*delvx + dely*delvy + delz*delvz;
        
        numtyp A_attij = coeff[mtype].x;
        numtyp B_repij = coeff[mtype].y;
        numtyp gammaij = coeff[mtype].z;
        numtyp sigmaij = coeff[mtype].w;
        numtyp cutij =   coeff2[mtype].x;
        numtyp cut_rij = coeff2[mtype].y;
        numtyp cutsqij = coeff2[mtype].z;
        
        numtyp wc = (numtyp)1.0 - r/cutij;
        numtyp wc_r = (numtyp)1.0 - r/cut_rij;
        wc_r = MAX(wc_r,(numtyp)0.0);
        numtyp wr = wc;

        const numtyp rhoj = extra[j].x;
        
        unsigned int tag1=itag, tag2=jtag;
        if (tag1 > tag2) {
          tag1 = jtag; tag2 = itag;
        }

        numtyp randnum = (numtyp)0.0;
        saru(tag1, tag2, seed, timestep, randnum);

        // conservative force = A_att * wc + B_rep*(rhoi+rhoj)*wc_r
        // drag force = -gamma * wr^2 * (delx dot delv) / r
        // random force = sigma * wr * rnd * dtinvsqrt;

        numtyp force = A_attij*wc + B_repij*(rhoi+rhoj)*wc_r;
        force -= gammaij*wr*wr*dot*rinv;
        force += sigmaij*wr*randnum*dtinvsqrt;
        force *= factor_dpd*rinv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          // unshifted eng of conservative term:
          // eng shifted to 0.0 at cutoff
          numtyp e = (numtyp)0.5*A_attij*cutij * wr*wr + (numtyp)0.5*B_repij*cut_rij*(rhoi+rhoj)*wc_r*wc_r;
          energy+=factor_dpd*e;
        }
        if (EVFLAG && vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
  } // if ii
  store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                ans,engv);
}

__kernel void k_mdpd_fast(const __global numtyp4 *restrict x_,
                          const __global numtyp4 *restrict extra,
                          const __global numtyp4 *restrict coeff_in,
                          const __global numtyp4 *restrict coeff2_in,
                          const __global numtyp *restrict sp_lj_in,
                          const __global numtyp *restrict sp_sqrt_in,
                          const __global int * dev_nbor,
                          const __global int * dev_packed,
                          __global acctyp3 *restrict ans,
                          __global acctyp *restrict engv,
                          const int eflag, const int vflag, const int inum,
                          const int nbor_pitch,
                          const __global numtyp4 *restrict v_,
                          const __global numtyp *restrict cutsq,
                          const numtyp dtinvsqrt, const int seed,
                          const int timestep, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  #ifndef ONETYPE
  __local numtyp4 coeff[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 coeff2[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  __local numtyp sp_sqrt[4];
  if (tid<4) {
    sp_lj[tid]=sp_lj_in[tid];
    sp_sqrt[tid]=sp_sqrt_in[tid];
  }
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    coeff[tid]=coeff_in[tid];
    coeff2[tid]=coeff2_in[tid];
  }
  __syncthreads();
  #else
  const numtyp A_attij=coeff_in[ONETYPE].x;
  const numtyp B_repij=coeff_in[ONETYPE].y;
  const numtyp gammaij=coeff_in[ONETYPE].z;
  const numtyp sigmaij=coeff_in[ONETYPE].w;
  const numtyp cutij=coeff2_in[ONETYPE].x;
  const numtyp cut_rij=coeff2_in[ONETYPE].y;
  const numtyp cutsqij=coeff2_in[ONETYPE].z;
  const numtyp cutsq_p=cutsq[ONETYPE];
  #endif

  int n_stride;
  local_allocate_store_pair();

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int i, numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    #ifndef ONETYPE
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);
    #endif
    numtyp4 iv; fetch4(iv,i,vel_tex); //v_[i];
    int itag=iv.w;

    const numtyp rhoi = extra[i].x;

    #ifndef ONETYPE
    numtyp factor_dpd;
    #endif
    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
      #ifndef ONETYPE
      factor_dpd = sp_lj[sbmask(j)];
      j &= NEIGHMASK;
      #endif

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      #ifndef ONETYPE
      int mtype=itype+jx.w;
      const numtyp cutsq_p=cutsq[mtype];
      #endif
      numtyp4 jv; fetch4(jv,j,vel_tex); //v_[j];
      int jtag=jv.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<cutsq_p) {
        numtyp r=ucl_sqrt(rsq);
        if (r < EPSILON) continue;

        numtyp rinv=ucl_recip(r);
        numtyp delvx = iv.x - jv.x;
        numtyp delvy = iv.y - jv.y;
        numtyp delvz = iv.z - jv.z;
        numtyp dot = delx*delvx + dely*delvy + delz*delvz;

        #ifndef ONETYPE
        numtyp A_attij = coeff[mtype].x;
        numtyp B_repij = coeff[mtype].y;
        numtyp gammaij = coeff[mtype].z;
        numtyp sigmaij = coeff[mtype].w;
        numtyp cutij =   coeff2[mtype].x;
        numtyp cut_rij = coeff2[mtype].y;
        numtyp cutsqij = coeff2[mtype].z;
        #endif

        numtyp wc = (numtyp)1.0 - r/cutij;
        numtyp wc_r = (numtyp)1.0 - r/cut_rij;
        wc_r = MAX(wc_r,(numtyp)0.0);
        numtyp wr = wc;

        const numtyp rhoj = extra[j].x;
        
        unsigned int tag1=itag, tag2=jtag;
        if (tag1 > tag2) {
          tag1 = jtag; tag2 = itag;
        }

        numtyp randnum = (numtyp)0.0;
        saru(tag1, tag2, seed, timestep, randnum);

        // conservative force = A_att * wc + B_rep*(rhoi+rhoj)*wc_r
        // drag force = -gamma * wr^2 * (delx dot delv) / r
        // random force = sigma * wr * rnd * dtinvsqrt;

        numtyp force = A_attij*wc + B_repij*(rhoi+rhoj)*wc_r;
        force -= gammaij*wr*wr*dot*rinv;
        force += sigmaij*wr*randnum*dtinvsqrt;
        #ifndef ONETYPE
        force *= factor_dpd*rinv;
        #else
        force*=rinv;
        #endif

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          // unshifted eng of conservative term:
          // eng shifted to 0.0 at cutoff
          numtyp e = (numtyp)0.5*A_attij*cutij * wr*wr + (numtyp)0.5*B_repij*cut_rij*(rhoi+rhoj)*wc_r*wc_r;
          #ifndef ONETYPE
          energy+=factor_dpd*e;
          #else
          energy+=e;
          #endif
        }
        if (EVFLAG && vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
  } // if ii
  store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                ans,engv);
}

