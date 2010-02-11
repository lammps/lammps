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
   Contributing authors: Mike Brown (SNL), wmbrown@sandia.gov
                         Peng Wang (Nvidia), penwang@nvidia.com
                         Paul Crozier (SNL), pscrozi@sandia.gov
------------------------------------------------------------------------- */

#ifndef LJ_GPU_KERNEL
#define LJ_GPU_KERNEL

/* Cell list version of LJ kernel */
template<bool eflag, bool vflag, int blockSize>
__global__ void kernel_lj_cell(float3 *force3,
			       float *energy, float3 *virial, 
			       float3 *cell_list, unsigned int *cell_idx, 
			       int *cell_type, int *cell_atom,
			       const int inum, const int nall, const int ncell, 
			       const int ncellx, const int ncelly, const int ncellz)
{
	
  
	
  // calculate 3D block idx from 2d block
  int bx = blockIdx.x;
  int by = blockIdx.y % ncelly;
  int bz = blockIdx.y / ncelly;

  int tid = threadIdx.x;
  
  // compute cell idx from 3D block idx
  int cid = bx + INT_MUL(by, ncellx) + INT_MUL(bz, INT_MUL(ncellx,ncelly));
  
  __shared__ int typeSh[blockSize];
  __shared__ float posSh[blockSize*3];
  __shared__ float cutsqSh[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ float lj1Sh[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ float lj2Sh[MAX_SHARED_TYPES*MAX_SHARED_TYPES];

  extern __shared__ float smem[];

  __shared__ float *lj3Sh;
  __shared__ float *lj4Sh;
  __shared__ float *offsetSh;

  // load force parameters into shared memory
  for (int i = tid; i < MAX_SHARED_TYPES*MAX_SHARED_TYPES; i += blockSize) {
    int itype = i/MAX_SHARED_TYPES;
    int jtype = i%MAX_SHARED_TYPES;
    cutsqSh[i] = _cutsq_<float>(itype,jtype);
    lj1Sh[i]   = _lj1_<float>(itype,jtype).x;
    lj2Sh[i]   = _lj1_<float>(itype,jtype).y;
  }

  // Only allocate shared memory when needed, 
  // this reduces shared memory limitation on occupancy
  if (eflag || vflag) {
    lj3Sh = smem;
    lj4Sh = lj3Sh + MAX_SHARED_TYPES*MAX_SHARED_TYPES;
    offsetSh = lj4Sh + MAX_SHARED_TYPES*MAX_SHARED_TYPES;
    for (int i = tid; i < MAX_SHARED_TYPES*MAX_SHARED_TYPES; i += blockSize) {
      int itype = i/MAX_SHARED_TYPES;
      int jtype = i%MAX_SHARED_TYPES;
      lj3Sh[i]   = _lj3_<float>(itype,jtype).x+0.01;
      lj4Sh[i]   = _lj3_<float>(itype,jtype).y;
      offsetSh[i]= _offset_<float>(itype,jtype);
    }
  }

  __syncthreads();

  int nborz0 = max(bz-1,0), nborz1 = min(bz+1, ncellz-1),
      nbory0 = max(by-1,0), nbory1 = min(by+1, ncelly-1),
      nborx0 = max(bx-1,0), nborx1 = min(bx+1, ncellx-1);

  for (int ii = 0; ii < ceil((float)(cell_atom[cid])/blockSize); ii++) {
    float3 f = {0.0f, 0.0f, 0.0f};
    float ener = 0.0f;
    float3 v0 = {0.0f, 0.0f, 0.0f}, v1 = {0.0f, 0.0f, 0.0f};
    int itype;
    float ix, iy, iz;
    int i = tid + ii*blockSize;
    unsigned int answer_pos = cell_idx[cid*blockSize+i];

    // load current cell atom position and type into sMem
    for (int j = tid; j < cell_atom[cid]; j += blockSize) {
      int pid = cid*blockSize + j;
      float3 pos = cell_list[pid];
      posSh[j            ] = pos.x;
      posSh[j+  blockSize] = pos.y;
      posSh[j+2*blockSize] = pos.z;
      typeSh[j]            = cell_type[pid];
    }
    __syncthreads();
    if (answer_pos < inum) {
      itype = typeSh[i];
      ix = posSh[i            ];
      iy = posSh[i+  blockSize];
      iz = posSh[i+2*blockSize];

      // compute force from current cell
      for (int j = 0; j < cell_atom[cid]; j++) {
	if (j == i) continue;
	float delx = ix - posSh[j            ];
	float dely = iy - posSh[j+  blockSize];
	float delz = iz - posSh[j+2*blockSize];
	int jtype = typeSh[j];
	int mtype = itype + jtype*MAX_SHARED_TYPES;
	float r2inv = delx*delx + dely*dely + delz*delz;
	
	if (r2inv < cutsqSh[mtype]) {
	  r2inv = 1.0f/r2inv;
	  float r6inv = r2inv * r2inv * r2inv;
	  float force = r2inv*r6inv*(lj1Sh[mtype]*r6inv - lj2Sh[mtype]);
	  f.x += delx * force;
	  f.y += dely * force;
	  f.z += delz * force;

	  if (eflag) {
	    float e = r6inv*(lj3Sh[mtype]*r6inv - lj4Sh[mtype]);
	    ener += (e - offsetSh[mtype]); 
	  }
	  
	  if (vflag) {
	    v0.x += delx*delx*force;
	    v0.y += dely*dely*force;
	    v0.z += delz*delz*force;
	    v1.x += delx*dely*force;
	    v1.y += delx*delz*force;
	    v1.z += dely*delz*force;
	  }

	} 
      }
    }
    __syncthreads();

    // compute force from neigboring cells
    for (int nborz = nborz0; nborz <= nborz1; nborz++) {
      for (int nbory = nbory0; nbory <= nbory1; nbory++) {
	for (int nborx = nborx0; nborx <= nborx1; nborx++) {
	  if (nborz == bz && nbory == by && nborx == bx) continue;
	  
	  // compute cell id
	  int cid_nbor = nborx + INT_MUL(nbory,ncellx) + 
	    INT_MUL(nborz,INT_MUL(ncellx,ncelly));
	
	  // load neighbor cell position and type into smem
	  for (int j = tid; j < cell_atom[cid_nbor]; j += blockSize) {
	    int pid = INT_MUL(cid_nbor,blockSize) + j;
	    float3 pos = cell_list[pid];
	    posSh[j            ] = pos.x;
	    posSh[j+  blockSize] = pos.y;
	    posSh[j+2*blockSize] = pos.z;
	    typeSh[j]           = cell_type[pid];
	  }
	  __syncthreads();
	  // compute force
	  if (answer_pos < inum) {
	    for (int j = 0; j < cell_atom[cid_nbor]; j++) {
	      float delx = ix - posSh[j           ];
	      float dely = iy - posSh[j+  blockSize];
	      float delz = iz - posSh[j+2*blockSize];
	      int jtype = typeSh[j];
	      int mtype = itype + jtype*MAX_SHARED_TYPES;
	      float r2inv = delx*delx + dely*dely + delz*delz;
	      
	      if (r2inv < cutsqSh[mtype]) {
		r2inv = 1.0f/r2inv;
		float r6inv = r2inv * r2inv * r2inv;
		float force = r2inv*r6inv*(lj1Sh[mtype]*r6inv - lj2Sh[mtype]);
		f.x += delx * force;
		f.y += dely * force;
		f.z += delz * force;

		if (eflag) {
		  float e=r6inv*(lj3Sh[mtype]*r6inv - lj4Sh[mtype]);				
		  ener += (e-offsetSh[mtype]); 
		}
		if (vflag) {
		  v0.x += delx*delx*force;
		  v0.y += dely*dely*force;
		  v0.z += delz*delz*force;
		  v1.x += delx*dely*force;
		  v1.y += delx*delz*force;
		  v1.z += dely*delz*force;
		}
	      }
	    }
	  }
	  __syncthreads();
	}
      }
    }

    if (answer_pos < inum) {
      force3[answer_pos] = f;
      if (eflag)
	energy[answer_pos] = ener;
      if (vflag) {
	virial[2*answer_pos] = v0;
	virial[2*answer_pos+1] = v1;
      }
    }
  }

}


/* Neigbhor list version of LJ kernel */
template<class numtyp, class acctyp>
__global__ void kernel_lj(const numtyp *special_lj, const int *dev_nbor, 
                          const int *dev_ij, const int nbor_pitch, acctyp *ans, 
                          size_t ans_pitch, const bool eflag, 
                          const bool vflag, const int inum, const int nall) {
  __shared__ numtyp sp_lj[4];
                            
  // ii indexes the two interacting particles in gi
  int ii=threadIdx.x;
  if (ii<4)
    sp_lj[ii]=special_lj[ii];    
  ii+=INT_MUL(blockIdx.x,blockDim.x);

  if (ii<inum) {
  
    acctyp energy=(numtyp)0;
    acctyp fx=(numtyp)0;
    acctyp fy=(numtyp)0;
    acctyp fz=(numtyp)0;
    acctyp virial[6];
    for (int i=0; i<6; i++)
      virial[i]=(numtyp)0;
  
    const int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    const int *list=dev_ij+*nbor;
    const int *list_end=list+numj;
  
    numtyp ix=_x_<numtyp>(i,0);
    numtyp iy=_x_<numtyp>(i,1);
    numtyp iz=_x_<numtyp>(i,2);
    int itype=_x_<numtyp>(i,3);

    numtyp factor_lj;
    for ( ; list<list_end; list++) {
  
      int j=*list;
      if (j < nall) 
        factor_lj = 1.0;
      else {
        factor_lj = sp_lj[j/nall];
        j %= nall;
      }
      int jtype=_x_<numtyp>(j,3);

      // Compute r12
      numtyp delx = ix-_x_<numtyp>(j,0);
      numtyp dely = iy-_x_<numtyp>(j,1);
      numtyp delz = iz-_x_<numtyp>(j,2);
      numtyp r2inv = delx*delx+dely*dely+delz*delz;
        
      if (r2inv<_cutsq_<numtyp>(itype,jtype)) {
        r2inv=(numtyp)1.0/r2inv;
        numtyp r6inv =r2inv*r2inv*r2inv;
        numtyp force =factor_lj*r2inv*r6inv*(_lj1_<numtyp>(itype,jtype).x*r6inv-
                                             _lj1_<numtyp>(itype,jtype).y);
      
        fx+=delx*force;
        fy+=dely*force;
        fz+=delz*force;

        if (eflag) {
          numtyp e=r6inv*(_lj3_<numtyp>(itype,jtype).x*r6inv-
                          _lj3_<numtyp>(itype,jtype).y);
          energy+=factor_lj*(e-_offset_<numtyp>(1,1)); 
        }
        if (vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor

    // Store answers
    acctyp *ap1=ans+ii;
    if (eflag) {
      *ap1=energy;
      ap1+=ans_pitch;
    }
    if (vflag) {
      for (int i=0; i<6; i++) {
        *ap1=virial[i];
        ap1+=ans_pitch;
      }
    }
    *ap1=fx;
    ap1+=ans_pitch;
    *ap1=fy;
    ap1+=ans_pitch;
    *ap1=fz;

  } // if ii
}

template<class numtyp, class acctyp>
__global__ void kernel_lj_fast(const numtyp *special_lj, const int *dev_nbor, 
                               const int *dev_ij, const int nbor_pitch, 
                               acctyp *ans, size_t ans_pitch,const bool eflag, 
                               const bool vflag, const int inum, 
                               const int nall) {
                                
  // ii indexes the two interacting particles in gi
  int ii=threadIdx.x;
  __shared__ numtyp sp_lj[4];
  __shared__ numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp lj2[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp lj4[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp offset[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  if (ii<4)
    sp_lj[ii]=special_lj[ii];    
  if (ii<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    int itype=ii/MAX_SHARED_TYPES;
    int jtype=ii%MAX_SHARED_TYPES;
    cutsq[ii]=_cutsq_<numtyp>(itype,jtype);
    lj1[ii]=_lj1_<numtyp>(itype,jtype).x;
    lj2[ii]=_lj1_<numtyp>(itype,jtype).y;
    if (eflag) {
      lj3[ii]=_lj3_<numtyp>(itype,jtype).x;
      lj4[ii]=_lj3_<numtyp>(itype,jtype).y;
      offset[ii]=_offset_<numtyp>(itype,jtype);
    }
  }
  ii+=INT_MUL(blockIdx.x,blockDim.x);
  
  if (ii<inum) {
  
    acctyp energy=(numtyp)0;
    acctyp fx=(numtyp)0;
    acctyp fy=(numtyp)0;
    acctyp fz=(numtyp)0;
    acctyp virial[6];
    for (int i=0; i<6; i++)
      virial[i]=(numtyp)0;
  
    const int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    const int *list=dev_ij+*nbor;
    const int *list_end=list+numj;
  
    numtyp ix=_x_<numtyp>(i,0);
    numtyp iy=_x_<numtyp>(i,1);
    numtyp iz=_x_<numtyp>(i,2);
    int itype=INT_MUL(MAX_SHARED_TYPES,_x_<numtyp>(i,3));

    numtyp factor_lj;

    for ( ; list<list_end; list++) {

      int j= *list;
      
      if (j < nall) 
        factor_lj = 1.0;
      else {
        factor_lj = sp_lj[j/nall];
        j %= nall;
      }
      int mtype=itype+_x_<numtyp>(j,3);

      // Compute r12
      numtyp delx = ix-_x_<numtyp>(j,0);
      numtyp dely = iy-_x_<numtyp>(j,1);
      numtyp delz = iz-_x_<numtyp>(j,2);
      numtyp r2inv = delx*delx+dely*dely+delz*delz;

      if (r2inv<cutsq[mtype]) {
        r2inv=(numtyp)1.0/r2inv;
        numtyp r6inv = r2inv*r2inv*r2inv;
        numtyp force = factor_lj*r2inv*r6inv*(lj1[mtype]*r6inv-lj2[mtype]);
      
        fx+=delx*force;
        fy+=dely*force;
        fz+=delz*force;

        if (eflag) {
          numtyp e=r6inv*(lj3[mtype]*r6inv-lj4[mtype]);
          energy+=factor_lj*(e-offset[mtype]); 
        }
        if (vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor

    // Store answers
    acctyp *ap1=ans+ii;
    if (eflag) {
      *ap1=energy;
      ap1+=ans_pitch;
    }
    if (vflag) {
      for (int i=0; i<6; i++) {
        *ap1=virial[i];
        ap1+=ans_pitch;
      }
    }
    *ap1=fx;
    ap1+=ans_pitch;
    *ap1=fy;
    ap1+=ans_pitch;
    *ap1=fz;

  } // if ii
}

#endif
