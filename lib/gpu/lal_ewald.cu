// **************************************************************************
//                                 ewald.cu
//                             -------------------
//                            W. Michael Brown (ORNL)
//                            Axel Kohlmeyer (Temple)
//
//  Device code for Ewald (reciprocal-space) acceleration
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : developers@lammps.org
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
_texture( q_tex,float);
#else
_texture_2d( pos_tex,int4);
_texture( q_tex,int2);
#endif
#else
#define pos_tex x_
#define q_tex q_
#endif

// Index into the cs/sn arrays.  The arrays store the cosine and sine of
// m*unitk[ic]*r for m in [-kmax,kmax] (full range, negatives included), all
// three Cartesian directions ic, and all local atoms i.  The atom index is
// the fastest-varying dimension so reads in the reduction kernel are coalesced.
#define CS_INDEX(m,ic,i,nlocal,kmax) ((((m)+(kmax))*3+(ic))*(nlocal)+(i))

// ---------------------------------------------------------------------------
// Compute cos(m*unitk[ic]*r_i) and sin(...) for every atom, direction, and
// |m| up to kmax using the standard cosine/sine recurrence.  One thread per
// atom.  Negative m are stored using cos(-m)=cos(m), sin(-m)=-sin(m).
// ---------------------------------------------------------------------------
__kernel void k_ewald_cssn(const __global numtyp4 *restrict x_,
                           __global numtyp *restrict cs,
                           __global numtyp *restrict sn,
                           const numtyp unitkx, const numtyp unitky,
                           const numtyp unitkz, const int kmax,
                           const int nlocal) {
  int i=GLOBAL_ID_X;
  if (i<nlocal) {
    numtyp4 p;
    fetch4(p,i,pos_tex);
    numtyp xx[3], unitk[3];
    xx[0]=p.x; xx[1]=p.y; xx[2]=p.z;
    unitk[0]=unitkx; unitk[1]=unitky; unitk[2]=unitkz;

    for (int ic=0; ic<3; ic++) {
      const numtyp arg=unitk[ic]*xx[ic];
      const numtyp c1=cos(arg);
      const numtyp s1=sin(arg);

      // m = 0
      cs[CS_INDEX(0,ic,i,nlocal,kmax)]=(numtyp)1.0;
      sn[CS_INDEX(0,ic,i,nlocal,kmax)]=(numtyp)0.0;

      // m = 1 and m = -1
      cs[CS_INDEX(1,ic,i,nlocal,kmax)]=c1;
      sn[CS_INDEX(1,ic,i,nlocal,kmax)]=s1;
      cs[CS_INDEX(-1,ic,i,nlocal,kmax)]=c1;
      sn[CS_INDEX(-1,ic,i,nlocal,kmax)]=-s1;

      numtyp csm1=c1, snm1=s1;
      for (int m=2; m<=kmax; m++) {
        const numtyp csm=csm1*c1-snm1*s1;
        const numtyp snm=snm1*c1+csm1*s1;
        cs[CS_INDEX(m,ic,i,nlocal,kmax)]=csm;
        sn[CS_INDEX(m,ic,i,nlocal,kmax)]=snm;
        cs[CS_INDEX(-m,ic,i,nlocal,kmax)]=csm;
        sn[CS_INDEX(-m,ic,i,nlocal,kmax)]=-snm;
        csm1=csm;
        snm1=snm;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Compute the local (per-MPI-rank) structure factors
//   sfacrl[k] = sum_i q_i cos(k.r_i),  sfacim[k] = sum_i q_i sin(k.r_i)
// One thread-block per k-vector reduces over all local atoms, so the result
// is deterministic (fixed-order tree reduction, no floating-point atomics).
// ---------------------------------------------------------------------------
__kernel void k_ewald_structure(const __global numtyp *restrict q_,
                                const __global numtyp *restrict cs,
                                const __global numtyp *restrict sn,
                                const __global int *restrict kxvecs,
                                const __global int *restrict kyvecs,
                                const __global int *restrict kzvecs,
                                __global acctyp *restrict sfacrl,
                                __global acctyp *restrict sfacim,
                                const int kmax, const int nlocal,
                                const int kcount) {
  __local acctyp red0[BLOCK_PAIR];
  __local acctyp red1[BLOCK_PAIR];

  const int tid=THREAD_ID_X;
  const int k=BLOCK_ID_X;

  acctyp sr=(acctyp)0.0;
  acctyp si=(acctyp)0.0;

  if (k<kcount) {
    const int kx=kxvecs[k];
    const int ky=kyvecs[k];
    const int kz=kzvecs[k];
    const int basex=CS_INDEX(kx,0,0,nlocal,kmax);
    const int basey=CS_INDEX(ky,1,0,nlocal,kmax);
    const int basez=CS_INDEX(kz,2,0,nlocal,kmax);

    for (int i=tid; i<nlocal; i+=BLOCK_SIZE_X) {
      const numtyp csx=cs[basex+i];
      const numtyp snx=sn[basex+i];
      const numtyp csy=cs[basey+i];
      const numtyp sny=sn[basey+i];
      const numtyp csz=cs[basez+i];
      const numtyp snz=sn[basez+i];
      const numtyp cypz=csy*csz-sny*snz;
      const numtyp sypz=sny*csz+csy*snz;
      const numtyp exprl=csx*cypz-snx*sypz;
      const numtyp expim=snx*cypz+csx*sypz;
      numtyp qi;
      fetch(qi,i,q_tex);
      sr+=(acctyp)(qi*exprl);
      si+=(acctyp)(qi*expim);
    }
  }

  // deterministic shared-memory tree reduction across the block.  This is
  // written out explicitly (rather than using block_reduce_add2 from
  // lal_aux_fun1.h, which only exists for SHUFFLE_AVAIL==0) so it compiles on
  // every backend regardless of whether warp-shuffle reductions are available.
  red0[tid]=sr;
  red1[tid]=si;
  __syncthreads();
  for (int s=BLOCK_SIZE_X/2; s>0; s>>=1) {
    if (tid<s) {
      red0[tid]+=red0[tid+s];
      red1[tid]+=red1[tid+s];
    }
    __syncthreads();
  }

  if (tid==0 && k<kcount) {
    sfacrl[k]=red0[0];
    sfacim[k]=red1[0];
  }
}

// ---------------------------------------------------------------------------
// K-space contribution to the per-atom electric field and force.
// One thread per atom sums over all k-vectors, reusing the resident cs/sn and
// the global structure factors sfacrl_all/sfacim_all.  Writes the k-space
// force qscale*q_i*ek_i into the answer array (merged with the pair force by
// fix gpu).
// ---------------------------------------------------------------------------
__kernel void k_ewald_field(const __global numtyp *restrict q_,
                            const __global numtyp *restrict cs,
                            const __global numtyp *restrict sn,
                            const __global int *restrict kxvecs,
                            const __global int *restrict kyvecs,
                            const __global int *restrict kzvecs,
                            const __global numtyp4 *restrict eg,
                            const __global numtyp *restrict ug,
                            const __global numtyp *restrict vg,
                            const __global acctyp *restrict sfacrl_all,
                            const __global acctyp *restrict sfacim_all,
                            __global acctyp3 *restrict ans,
                            __global acctyp *restrict eatom,
                            __global acctyp *restrict vatom,
                            const numtyp qscale, const int slabflag,
                            const int eflag_atom, const int vflag_atom,
                            const int kmax, const int nlocal, const int kcount) {
  int i=GLOBAL_ID_X;
  if (i<nlocal) {
    acctyp ekx=(acctyp)0.0;
    acctyp eky=(acctyp)0.0;
    acctyp ekz=(acctyp)0.0;
    acctyp ea=(acctyp)0.0;
    acctyp va0=(acctyp)0.0, va1=(acctyp)0.0, va2=(acctyp)0.0;
    acctyp va3=(acctyp)0.0, va4=(acctyp)0.0, va5=(acctyp)0.0;

    for (int k=0; k<kcount; k++) {
      const int kx=kxvecs[k];
      const int ky=kyvecs[k];
      const int kz=kzvecs[k];
      const numtyp csx=cs[CS_INDEX(kx,0,i,nlocal,kmax)];
      const numtyp snx=sn[CS_INDEX(kx,0,i,nlocal,kmax)];
      const numtyp csy=cs[CS_INDEX(ky,1,i,nlocal,kmax)];
      const numtyp sny=sn[CS_INDEX(ky,1,i,nlocal,kmax)];
      const numtyp csz=cs[CS_INDEX(kz,2,i,nlocal,kmax)];
      const numtyp snz=sn[CS_INDEX(kz,2,i,nlocal,kmax)];
      const numtyp cypz=csy*csz-sny*snz;
      const numtyp sypz=sny*csz+csy*snz;
      const numtyp exprl=csx*cypz-snx*sypz;
      const numtyp expim=snx*cypz+csx*sypz;
      const acctyp partial=expim*sfacrl_all[k]-exprl*sfacim_all[k];
      const numtyp4 egk=eg[k];
      ekx+=partial*egk.x;
      eky+=partial*egk.y;
      ekz+=partial*egk.z;

      if (eflag_atom || vflag_atom) {
        const acctyp pp=exprl*sfacrl_all[k]+expim*sfacim_all[k];
        const acctyp ugpp=ug[k]*pp;
        if (eflag_atom) ea+=ugpp;
        if (vflag_atom) {
          va0+=vg[k*6  ]*ugpp;
          va1+=vg[k*6+1]*ugpp;
          va2+=vg[k*6+2]*ugpp;
          va3+=vg[k*6+3]*ugpp;
          va4+=vg[k*6+4]*ugpp;
          va5+=vg[k*6+5]*ugpp;
        }
      }
    }

    numtyp qi;
    fetch(qi,i,q_tex);
    const numtyp qfac=qscale*qi;
    acctyp3 f;
    f.x=qfac*ekx;
    f.y=qfac*eky;
    f.z=(slabflag!=2) ? qfac*ekz : (acctyp)0.0;
    ans[i]=f;

    // raw per-atom contributions; the host applies the q_i factor, the
    // self-energy correction, and qscale to match the CPU exactly
    if (eflag_atom) eatom[i]=ea;
    if (vflag_atom) {
      vatom[i*6  ]=va0;
      vatom[i*6+1]=va1;
      vatom[i*6+2]=va2;
      vatom[i*6+3]=va3;
      vatom[i*6+4]=va4;
      vatom[i*6+5]=va5;
    }
  }
}
