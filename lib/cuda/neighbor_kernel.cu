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

#define SBBITS 30

__global__ void Binning_Kernel(int* binned_id, int bin_nmax, int bin_dim_x, int bin_dim_y, int bin_dim_z,
                               CUDA_CFLOAT rez_bin_size_x, CUDA_CFLOAT rez_bin_size_y, CUDA_CFLOAT rez_bin_size_z)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  /*int* bin_count=(int*) _buffer;
  bin_count=bin_count+20;
  CUDA_CFLOAT* binned_x=(CUDA_CFLOAT*)(bin_count+bin_dim_x*bin_dim_y*bin_dim_z);*/
  CUDA_CFLOAT* binned_x = (CUDA_CFLOAT*) _buffer;
  binned_x = &binned_x[2];
  int* bin_count = (int*) &binned_x[3 * bin_dim_x * bin_dim_y * bin_dim_z * bin_nmax];

  if(i < _nall) {
    // copy atom position from global device memory to local register
    // in this 3 steps to get as much coalesced access as possible
    X_CFLOAT* my_x = _x + i;
    CUDA_CFLOAT x_i = *my_x;
    my_x += _nmax;
    CUDA_CFLOAT y_i = *my_x;
    my_x += _nmax;
    CUDA_CFLOAT z_i = *my_x;


    // calculate flat bin index
    int bx = __float2int_rd(rez_bin_size_x * (x_i - _sublo[0])) + 2;
    int by = __float2int_rd(rez_bin_size_y * (y_i - _sublo[1])) + 2;
    int bz = __float2int_rd(rez_bin_size_z * (z_i - _sublo[2])) + 2;

    bx -= bx * negativCUDA(1.0f * bx);
    bx -= (bx - bin_dim_x + 1) * negativCUDA(1.0f * bin_dim_x - 1.0f - 1.0f * bx);
    by -= by * negativCUDA(1.0f * by);
    by -= (by - bin_dim_y + 1) * negativCUDA(1.0f * bin_dim_y - 1.0f - 1.0f * by);
    bz -= bz * negativCUDA(1.0f * bz);
    bz -= (bz - bin_dim_z + 1) * negativCUDA(1.0f * bin_dim_z - 1.0f - 1.0f * bz);


    const unsigned j = bin_dim_z * (bin_dim_y * bx + by) + bz;

    // add new atom to bin, get bin-array position
    const unsigned k = atomicAdd(& bin_count[j], 1);

    if(k < bin_nmax) {
      binned_id [bin_nmax * j + k] = i;
      binned_x [3 * bin_nmax * j + k] = x_i;
      binned_x [3 * bin_nmax * j + k + bin_nmax] = y_i;
      binned_x [3 * bin_nmax * j + k + 2 * bin_nmax] = z_i;
    } else {
      // normally, this should not happen:
      int errorn = atomicAdd((int*) _buffer, 1);
      MYEMUDBG(printf("# CUDA: Binning_Kernel: WARNING: atom %i ignored, no place left in bin %u\n", i, j);)
    }
  }
}


__device__ inline int exclusion(int &i, int &j, int &itype, int &jtype)
{
  int m;

  if(_nex_type)
    if(_ex_type[itype * _cuda_ntypes + jtype]) return 1;

  if(_nex_group) {
    for(m = 0; m < _nex_group; m++) {
      if(_mask[i] & _ex1_bit[m] && _mask[j] & _ex2_bit[m]) return 1;

      if(_mask[i] & _ex2_bit[m] && _mask[j] & _ex1_bit[m]) return 1;
    }
  }

  if(_nex_mol) {
    if(_molecule[i] == _molecule[j])
      for(m = 0; m < _nex_mol; m++)
        if(_mask[i] & _ex_mol_bit[m] && _mask[j] & _ex_mol_bit[m]) return 1;
  }

  return 0;
}

extern __shared__ CUDA_CFLOAT shared[];

__device__ inline int find_special(int3 &n, int* list, int &tag, int3 flag)
{
  int k = n.z;

  for(int l = 0; l < n.z; l++) k = ((list[l] == tag) ? l : k);

  return k < n.x ? flag.x : (k < n.y ? flag.y : (k < n.z ? flag.z : 0));
}

template <const unsigned int exclude>
__global__ void NeighborBuildFullBin_Kernel(int* binned_id, int bin_nmax, int bin_dim_x, int bin_dim_y, CUDA_CFLOAT globcutoff, int block_style, bool neighall)
{
  int natoms = neighall ? _nall : _nlocal;
  //const bool domol=false;
  int bin_dim_z = gridDim.y;
  CUDA_CFLOAT* binned_x = (CUDA_CFLOAT*) _buffer;
  binned_x = &binned_x[2];
  int* bin_count = (int*) &binned_x[3 * bin_dim_x * bin_dim_y * bin_dim_z * bin_nmax];
  int bin = __mul24(gridDim.y, blockIdx.x) + blockIdx.y;
  int bin_x = blockIdx.x / bin_dim_y;
  int bin_y = blockIdx.x - bin_x * bin_dim_y;
  int bin_z = blockIdx.y;
  int bin_c = bin_count[bin];


  CUDA_CFLOAT cut;

  if(globcutoff > 0)
    cut = globcutoff;

  int i = _nall;
  CUDA_CFLOAT* my_x;
  CUDA_CFLOAT x_i, y_i, z_i;

  for(int actOffset = 0; actOffset < bin_c; actOffset += blockDim.x) {

    int actIdx = threadIdx.x + actOffset;
    CUDA_CFLOAT* other_x = shared;
    int* other_id = (int*) &other_x[3 * blockDim.x];

    if(actIdx < bin_c) {
      i = binned_id[__mul24(bin, bin_nmax) + actIdx];
      my_x = binned_x + __mul24(__mul24(bin, 3), bin_nmax) + actIdx;
      x_i = *my_x;
      my_x += bin_nmax;
      y_i = *my_x;
      my_x += bin_nmax;
      z_i = *my_x;
    } else
      i = 2 * _nall;

    __syncthreads();

    int jnum = 0;
    int itype;

    if(i < natoms) {
      jnum = 0;
      _ilist[i] = i;
      itype = _type[i];
    }

    //__syncthreads();


    for(int otherActOffset = 0; otherActOffset < bin_c; otherActOffset += blockDim.x) {
      int otherActIdx = threadIdx.x + otherActOffset;

      if(otherActIdx < bin_c) {
        if(otherActOffset == actOffset) {
          other_id[threadIdx.x] = i;
          other_x[threadIdx.x] = x_i;
          other_x[threadIdx.x + blockDim.x] = y_i;
          other_x[threadIdx.x + 2 * blockDim.x] = z_i;
        } else {
          other_id[threadIdx.x] = binned_id[__mul24(bin, bin_nmax) + otherActIdx];
          my_x = binned_x + __mul24(__mul24(bin, 3), bin_nmax) + otherActIdx;
          other_x[threadIdx.x] = *my_x;
          my_x += bin_nmax;
          other_x[threadIdx.x + blockDim.x] = *my_x;
          my_x += bin_nmax;
          other_x[threadIdx.x + __mul24(2, blockDim.x)] = *my_x;

        }
      }

      __syncthreads();
      int kk = threadIdx.x;

      for(int k = 0; k < MIN(bin_c - otherActOffset, blockDim.x); ++k) {
        if(i < natoms) {
          kk++;
          kk = kk < MIN(bin_c - otherActOffset, blockDim.x) ? kk : 0;
          int j = other_id[kk];

          if(exclude && exclusion(i, j, itype, _type[j])) continue;

          if(globcutoff < 0) {
            int jtype = _type[j];
            cut = _cutneighsq[itype * _cuda_ntypes + jtype];
          }

          CUDA_CFLOAT delx = x_i - other_x[kk];
          CUDA_CFLOAT dely = y_i - other_x[kk + blockDim.x];
          CUDA_CFLOAT delz = z_i - other_x[kk + 2 * blockDim.x];
          CUDA_CFLOAT rsq = delx * delx + dely * dely + delz * delz;


          if(rsq <= cut && i != j) {
            if(jnum < _maxneighbors) {
              if(block_style)
                _neighbors[i * _maxneighbors + jnum] = j;
              else
                _neighbors[i + jnum * natoms] = j;
            }

            ++jnum;
          }
        }
      }

      __syncthreads();

    }

    for(int obin_x = bin_x - 1; obin_x < bin_x + 2; obin_x++)
      for(int obin_y = bin_y - 1; obin_y < bin_y + 2; obin_y++)
        for(int obin_z = bin_z - 1; obin_z < bin_z + 2; obin_z++) {
          if(obin_x < 0 || obin_y < 0 || obin_z < 0) continue;

          if(obin_x >= bin_dim_x || obin_y >= bin_dim_y || obin_z >= bin_dim_z) continue;

          int other_bin = bin_dim_z * (bin_dim_y * obin_x + obin_y) + obin_z;

          if(other_bin == bin) continue;

          int obin_c = bin_count[other_bin];

          for(int otherActOffset = 0; otherActOffset < obin_c; otherActOffset += blockDim.x) {
            int otherActIdx = otherActOffset + threadIdx.x;

            if(threadIdx.x < MIN(blockDim.x, obin_c - otherActOffset)) {
              other_id[threadIdx.x] = binned_id[__mul24(other_bin, bin_nmax) + otherActIdx];
              my_x = binned_x + __mul24(__mul24(other_bin, 3), bin_nmax) + otherActIdx;
              other_x[threadIdx.x] = *my_x;
              my_x += bin_nmax;
              other_x[threadIdx.x + blockDim.x] = *my_x;
              my_x += bin_nmax;
              other_x[threadIdx.x + 2 * blockDim.x] = *my_x;
            }

            __syncthreads();

            for(int k = 0; k < MIN(blockDim.x, obin_c - otherActOffset); ++k) {
              if(i < natoms) {
                int j = other_id[k];

                if(exclude && exclusion(i, j, itype, _type[j])) continue;

                if(globcutoff < 0) {
                  int jtype = _type[j];
                  cut = _cutneighsq[itype * _cuda_ntypes + jtype];
                }

                CUDA_CFLOAT delx = x_i - other_x[k];
                CUDA_CFLOAT dely = y_i - other_x[k + blockDim.x];
                CUDA_CFLOAT delz = z_i - other_x[k + 2 * blockDim.x];
                CUDA_CFLOAT rsq = delx * delx + dely * dely + delz * delz;

                if(rsq <= cut && i != j) {
                  if(jnum < _maxneighbors) {
                    if(block_style)
                      _neighbors[i * _maxneighbors + jnum] = j;
                    else
                      _neighbors[i + jnum * natoms] = j;
                  }

                  ++jnum;
                }
              }
            }

            __syncthreads();

          }
        }

    if(jnum > _maxneighbors)((int*)_buffer)[0] = -jnum;

    if(i < natoms)
      _numneigh[i] = jnum;
  }
}


__global__ void FindSpecial(int block_style)
{
  int ii = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  int which;
  int tag_mask = 0;
  int3 spec_flag;

  int3 mynspecial = {0, 0, 1};

  if(ii >= _nlocal) return;

  int special_id[CUDA_MAX_NSPECIAL];

  int i = _ilist[ii];

  if(i >= _nlocal) return;

  int jnum = _numneigh[i];

  if(_special_flag[1] == 0) spec_flag.x = -1;
  else if(_special_flag[1] == 1) spec_flag.x = 0;
  else spec_flag.x = 1;

  if(_special_flag[2] == 0) spec_flag.y = -1;
  else if(_special_flag[2] == 1) spec_flag.y = 0;
  else spec_flag.y = 2;

  if(_special_flag[3] == 0) spec_flag.z = -1;
  else if(_special_flag[3] == 1) spec_flag.z = 0;
  else spec_flag.z = 3;

  mynspecial.x = _nspecial[i];
  mynspecial.y = _nspecial[i + _nmax];
  mynspecial.z = _nspecial[i + 2 * _nmax];

  if(i < _nlocal) {
    int* list = &_special[i];

    for(int k = 0; k < mynspecial.z; k++) {
      special_id[k] = list[k * _nmax];
      tag_mask = tag_mask | special_id[k];
    }
  }


  for(int k = 0; k < MIN(jnum, _maxneighbors); k++) {
    int j;

    if(block_style)
      j = _neighbors[i * _maxneighbors + k];
    else
      j = _neighbors[i + k * _nlocal];

    int tag_j = _tag[j];
    which = 0;

    if((tag_mask & tag_j) == tag_j) {
      which = find_special(mynspecial, special_id, tag_j, spec_flag);

      if(which > 0) {
        if(block_style)
          _neighbors[i * _maxneighbors + k] = j ^ (which << SBBITS);
        else
          _neighbors[i + k * _nlocal] = j ^ (which << SBBITS);
      } else if(which < 0) {
        if(block_style)
          _neighbors[i * _maxneighbors + k] = _neighbors[i * _maxneighbors + jnum - 1];
        else
          _neighbors[i + k * _nlocal] = _neighbors[i + (jnum - 1) * _nlocal];

        jnum--;
        k--;
      }
    }
  }

  _numneigh[i] = jnum;
}

__global__ void NeighborBuildFullBin_OverlapComm_Kernel(int* binned_id, int bin_nmax, int bin_dim_x, int bin_dim_y, CUDA_CFLOAT globcutoff, int block_style)
{
  int bin_dim_z = gridDim.y;
  CUDA_CFLOAT* binned_x = (CUDA_CFLOAT*) _buffer;
  binned_x = &binned_x[2];
  int* bin_count = (int*) &binned_x[3 * bin_dim_x * bin_dim_y * bin_dim_z * bin_nmax];
  int bin = __mul24(gridDim.y, blockIdx.x) + blockIdx.y;
  int bin_x = blockIdx.x / bin_dim_y;
  int bin_y = blockIdx.x - bin_x * bin_dim_y;
  int bin_z = blockIdx.y;
  int bin_c = bin_count[bin];


  CUDA_CFLOAT cut;

  if(globcutoff > 0)
    cut = globcutoff;

  int i = _nall;
  CUDA_CFLOAT* my_x;
  CUDA_CFLOAT x_i, y_i, z_i;

  for(int actOffset = 0; actOffset < bin_c; actOffset += blockDim.x) {

    int actIdx = threadIdx.x + actOffset;
    CUDA_CFLOAT* other_x = shared;
    int* other_id = (int*) &other_x[3 * blockDim.x];

    if(actIdx < bin_c) {
      i = binned_id[__mul24(bin, bin_nmax) + actIdx];
      my_x = binned_x + __mul24(__mul24(bin, 3), bin_nmax) + actIdx;
      x_i = *my_x;
      my_x += bin_nmax;
      y_i = *my_x;
      my_x += bin_nmax;
      z_i = *my_x;
    } else
      i = 2 * _nall;

    __syncthreads();

    int jnum = 0;
    int jnum_border = 0;
    int jnum_inner = 0;
    int i_border = -1;
    int itype;

    if(i < _nlocal) {
      jnum = 0;
      _ilist[i] = i;
      itype = _type[i];
    }

    __syncthreads();


    for(int otherActOffset = 0; otherActOffset < bin_c; otherActOffset += blockDim.x) {
      int otherActIdx = threadIdx.x + otherActOffset;

      if(otherActIdx < bin_c) {
        if(otherActOffset == actOffset) {
          other_id[threadIdx.x] = i;
          other_x[threadIdx.x] = x_i;
          other_x[threadIdx.x + blockDim.x] = y_i;
          other_x[threadIdx.x + 2 * blockDim.x] = z_i;
        } else {
          other_id[threadIdx.x] = binned_id[__mul24(bin, bin_nmax) + otherActIdx];
          my_x = binned_x + __mul24(__mul24(bin, 3), bin_nmax) + otherActIdx;
          other_x[threadIdx.x] = *my_x;
          my_x += bin_nmax;
          other_x[threadIdx.x + blockDim.x] = *my_x;
          my_x += bin_nmax;
          other_x[threadIdx.x + __mul24(2, blockDim.x)] = *my_x;

        }
      }

      __syncthreads();
      int kk = threadIdx.x;

      for(int k = 0; k < MIN(bin_c - otherActOffset, blockDim.x); ++k) {
        if(i < _nlocal) {
          kk++;
          kk = kk < MIN(bin_c - otherActOffset, blockDim.x) ? kk : 0;
          int j = other_id[kk];

          if(globcutoff < 0) {
            int jtype = _type[j];
            cut = _cutneighsq[itype * _cuda_ntypes + jtype];
          }

          CUDA_CFLOAT delx = x_i - other_x[kk];
          CUDA_CFLOAT dely = y_i - other_x[kk + blockDim.x];
          CUDA_CFLOAT delz = z_i - other_x[kk + 2 * blockDim.x];
          CUDA_CFLOAT rsq = delx * delx + dely * dely + delz * delz;


          if(rsq <= cut && i != j) {
            if((j >= _nlocal) && (i_border < 0))
              i_border = atomicAdd(_inum_border, 1);

            if(jnum < _maxneighbors) {
              if(block_style) {
                _neighbors[i * _maxneighbors + jnum] = j;

                if(j >= _nlocal) {
                  _neighbors_border[i_border * _maxneighbors + jnum_border] = j;
                } else {
                  _neighbors_inner[i * _maxneighbors + jnum_inner] = j;
                }
              } else {
                _neighbors[i + jnum * _nlocal] = j;

                if(j >= _nlocal) {
                  _neighbors_border[i_border + jnum_border * _nlocal] = j;
                } else {
                  _neighbors_inner[i + jnum_inner * _nlocal] = j;
                }
              }
            }

            ++jnum;

            if(j >= _nlocal)
              jnum_border++;
            else
              jnum_inner++;
          }
        }
      }

      __syncthreads();
    }

    for(int obin_x = bin_x - 1; obin_x < bin_x + 2; obin_x++)
      for(int obin_y = bin_y - 1; obin_y < bin_y + 2; obin_y++)
        for(int obin_z = bin_z - 1; obin_z < bin_z + 2; obin_z++) {
          if(obin_x < 0 || obin_y < 0 || obin_z < 0) continue;

          if(obin_x >= bin_dim_x || obin_y >= bin_dim_y || obin_z >= bin_dim_z) continue;

          int other_bin = bin_dim_z * (bin_dim_y * obin_x + obin_y) + obin_z;

          if(other_bin == bin) continue;

          int obin_c = bin_count[other_bin];

          for(int otherActOffset = 0; otherActOffset < obin_c; otherActOffset += blockDim.x) {
            int otherActIdx = otherActOffset + threadIdx.x;

            if(threadIdx.x < MIN(blockDim.x, obin_c - otherActOffset)) {
              other_id[threadIdx.x] = binned_id[__mul24(other_bin, bin_nmax) + otherActIdx];
              my_x = binned_x + __mul24(__mul24(other_bin, 3), bin_nmax) + otherActIdx;
              other_x[threadIdx.x] = *my_x;
              my_x += bin_nmax;
              other_x[threadIdx.x + blockDim.x] = *my_x;
              my_x += bin_nmax;
              other_x[threadIdx.x + 2 * blockDim.x] = *my_x;
            }

            __syncthreads();

            for(int k = 0; k < MIN(blockDim.x, obin_c - otherActOffset); ++k) {
              if(i < _nlocal) {
                int j = other_id[k];

                if(globcutoff < 0) {
                  int jtype = _type[j];
                  cut = _cutneighsq[itype * _cuda_ntypes + jtype];
                }

                CUDA_CFLOAT delx = x_i - other_x[k];
                CUDA_CFLOAT dely = y_i - other_x[k + blockDim.x];
                CUDA_CFLOAT delz = z_i - other_x[k + 2 * blockDim.x];
                CUDA_CFLOAT rsq = delx * delx + dely * dely + delz * delz;

                if(rsq <= cut && i != j) {
                  if((j >= _nlocal) && (i_border < 0))
                    i_border = atomicAdd(_inum_border, 1);

                  if(jnum < _maxneighbors) {
                    if(block_style) {
                      _neighbors[i * _maxneighbors + jnum] = j;

                      if(j >= _nlocal) {
                        _neighbors_border[i_border * _maxneighbors + jnum_border] = j;
                      } else {
                        _neighbors_inner[i * _maxneighbors + jnum_inner] = j;
                      }
                    } else {
                      _neighbors[i + jnum * _nlocal] = j;

                      if(j >= _nlocal) {
                        _neighbors_border[i_border + jnum_border * _nlocal] = j;
                      } else {
                        _neighbors_inner[i + jnum_inner * _nlocal] = j;
                      }
                    }
                  }

                  ++jnum;

                  if(j >= _nlocal)
                    jnum_border++;
                  else
                    jnum_inner++;
                }
              }
            }

            __syncthreads();
          }
        }

    if(jnum > _maxneighbors)((int*)_buffer)[0] = -jnum;

    if(i < _nlocal) {
      _numneigh[i] = jnum;
      _numneigh_inner[i] = jnum_inner;

      if(i_border >= 0) _numneigh_border[i_border] = jnum_border;

      if(i_border >= 0) _ilist_border[i_border] = i;

    }
  }
}

__global__ void NeighborBuildFullNsq_Kernel()
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  int* buffer = (int*) _buffer;

  if(i < _nlocal) {
    X_CFLOAT* my_x = _x + i;
    CUDA_CFLOAT x_i = *my_x;
    my_x += _nmax;
    CUDA_CFLOAT y_i = *my_x;
    my_x += _nmax;
    CUDA_CFLOAT z_i = *my_x;
    int jnum = 0;
    int* jlist = _firstneigh[i];
    _ilist[i] = i;

    int itype = _type[i];
    __syncthreads();

    for(int j = 0; j < _nall; ++j) {
      my_x = _x + j;
      CUDA_CFLOAT x_j = *my_x;
      my_x += _nmax;
      CUDA_CFLOAT y_j = *my_x;
      my_x += _nmax;
      CUDA_CFLOAT z_j = *my_x;
      CUDA_CFLOAT delx = x_i - x_j;
      CUDA_CFLOAT dely = y_i - y_j;
      CUDA_CFLOAT delz = z_i - z_j;
      CUDA_CFLOAT rsq = delx * delx + dely * dely + delz * delz;
      int jtype = _type[j];

      if(rsq <= _cutneighsq[itype * _cuda_ntypes + jtype] && i != j) {
        if(jnum < _maxneighbors)
          jlist[jnum] = j;

        if(i == 151)((int*)_buffer)[jnum + 2] = j;

        ++jnum;
      }

      __syncthreads();
    }

    if(jnum > _maxneighbors) buffer[0] = 0;

    _numneigh[i] = jnum;

    if(i == 151)((int*)_buffer)[1] = jnum;
  }
}

