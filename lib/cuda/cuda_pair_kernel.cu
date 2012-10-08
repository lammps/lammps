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
#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429


template <const PAIR_FORCES pair_type, const COUL_FORCES coul_type, const unsigned int extended_data>
__global__ void Pair_Kernel_TpA(int eflag, int vflag, int eflag_atom, int vflag_atom)
{
  ENERGY_FLOAT evdwl = ENERGY_F(0.0);
  ENERGY_FLOAT ecoul = ENERGY_F(0.0);

  ENERGY_FLOAT* sharedE;
  ENERGY_FLOAT* sharedECoul;
  ENERGY_FLOAT* sharedV = &sharedmem[threadIdx.x];

  if(eflag || eflag_atom) {
    sharedE = &sharedmem[threadIdx.x];
    sharedE[0] = ENERGY_F(0.0);
    sharedV += blockDim.x;

    if(coul_type != COUL_NONE) {
      sharedECoul = sharedE + blockDim.x;
      sharedECoul[0] = ENERGY_F(0.0);
      sharedV += blockDim.x;
    }
  }

  if(vflag || vflag_atom) {
    sharedV[0 * blockDim.x] = ENERGY_F(0.0);
    sharedV[1 * blockDim.x] = ENERGY_F(0.0);
    sharedV[2 * blockDim.x] = ENERGY_F(0.0);
    sharedV[3 * blockDim.x] = ENERGY_F(0.0);
    sharedV[4 * blockDim.x] = ENERGY_F(0.0);
    sharedV[5 * blockDim.x] = ENERGY_F(0.0);
  }

  int ii = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  X_FLOAT xtmp, ytmp, ztmp;
  X_FLOAT4 myxtype;
  F_FLOAT fxtmp, fytmp, fztmp, fpair;
  F_FLOAT delx, dely, delz;
  F_FLOAT factor_lj, factor_coul;
  F_FLOAT qtmp;
  int itype, i, j;
  int jnum = 0;
  int* jlist;

  if(ii < _inum) {
    i = _ilist[ii];

    myxtype = fetchXType(i);
    xtmp = myxtype.x;
    ytmp = myxtype.y;
    ztmp = myxtype.z;
    itype = static_cast <int>(myxtype.w);


    fxtmp = F_F(0.0);
    fytmp = F_F(0.0);
    fztmp = F_F(0.0);

    if(coul_type != COUL_NONE)
      qtmp = fetchQ(i);

    jnum = _numneigh[i];
    jlist = &_neighbors[i];
  }

  __syncthreads();

  for(int jj = 0; jj < jnum; jj++) {
    if(ii < _inum)
      if(jj < jnum) {
        fpair = F_F(0.0);
        j = jlist[jj * _nlocal];
        factor_lj =  _special_lj[sbmask(j)];

        if(coul_type != COUL_NONE)
          factor_coul = _special_coul[sbmask(j)];

        j &= NEIGHMASK;

        myxtype = fetchXType(j);
        delx = xtmp - myxtype.x;
        dely = ytmp - myxtype.y;
        delz = ztmp - myxtype.z;
        int jtype = static_cast <int>(myxtype.w);


        const F_FLOAT rsq = delx * delx + dely * dely + delz * delz;

        bool in_cutoff = rsq < (_cutsq_global > X_F(0.0) ? _cutsq_global : _cutsq[itype * _cuda_ntypes + jtype]);

        if(in_cutoff) {
          switch(pair_type) {
            case PAIR_BORN:
              fpair += PairBornCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_BUCK:
              fpair += PairBuckCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_CG_CMM:
              fpair += PairLJSDKCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ_CHARMM:
              fpair += PairLJCharmmCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ_CLASS2:
              fpair += PairLJClass2Cuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ_CUT:
              fpair += PairLJCutCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ_EXPAND:
              fpair += PairLJExpandCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ_GROMACS:
              fpair += PairLJGromacsCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ_SMOOTH:
              fpair += PairLJSmoothCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ96_CUT:
              fpair += PairLJ96CutCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_MORSE_R6:
              fpair += PairMorseR6Cuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_MORSE:
              fpair += PairMorseCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;
          }
        }

        if(coul_type != COUL_NONE) {
          const F_FLOAT qiqj = qtmp * fetchQ(j);

          if(qiqj * qiqj > 1e-8) {
            const bool in_coul_cutoff =
              rsq < (_cut_coulsq_global > X_F(0.0) ? _cut_coulsq_global : _cut_coulsq[itype * _cuda_ntypes + jtype]);

            if(in_coul_cutoff) {
              switch(coul_type) {
                case COUL_CHARMM:
                  fpair += CoulCharmmCuda_Eval(rsq, factor_coul, eflag, ecoul, qiqj);
                  break;

                case COUL_CHARMM_IMPLICIT:
                  fpair += CoulCharmmImplicitCuda_Eval(rsq, factor_coul, eflag, ecoul, qiqj);
                  break;

                case COUL_CUT: {
                  const F_FLOAT forcecoul = factor_coul * _qqrd2e * qiqj * _RSQRT_(rsq);

                  if(eflag) {
                    ecoul += forcecoul;
                  }

                  fpair += forcecoul * (F_F(1.0) / rsq);
                }
                break;

                case COUL_DEBYE: {
                  const F_FLOAT r2inv = F_F(1.0) / rsq;
                  const X_FLOAT r = _RSQRT_(r2inv);
                  const X_FLOAT rinv = F_F(1.0) / r;
                  const F_FLOAT screening = _EXP_(-_kappa * r);
                  F_FLOAT forcecoul = factor_coul * _qqrd2e * qiqj * screening ;

                  if(eflag) {
                    ecoul += forcecoul * rinv;
                  }

                  forcecoul *= (_kappa + rinv);
                  fpair += forcecoul * r2inv;
                }
                break;

                case COUL_GROMACS:
                  fpair += CoulGromacsCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_coul, eflag, ecoul, qiqj);
                  break;

                case COUL_LONG: {
                  const F_FLOAT r2inv = F_F(1.0) / rsq;
                  const F_FLOAT r = _RSQRT_(r2inv);
                  const F_FLOAT grij = _g_ewald * r;
                  const F_FLOAT expm2 = _EXP_(-grij * grij);
                  const F_FLOAT t = F_F(1.0) / (F_F(1.0) + EWALD_P * grij);
                  const F_FLOAT erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
                  const F_FLOAT prefactor = _qqrd2e * qiqj * (F_F(1.0) / r);
                  F_FLOAT forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);

                  if(factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;

                  if(eflag) {
                    ecoul += prefactor * erfc;

                    if(factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
                  }

                  fpair += forcecoul * r2inv;
                }
                break;
              }
            }

            in_cutoff = in_cutoff || in_coul_cutoff;
          }
        }


        if(in_cutoff) {
          F_FLOAT dxfp, dyfp, dzfp;
          fxtmp += dxfp = delx * fpair;
          fytmp += dyfp = dely * fpair;
          fztmp += dzfp = delz * fpair;

          if(vflag) {
            sharedV[0 * blockDim.x] += delx * dxfp;
            sharedV[1 * blockDim.x] += dely * dyfp;
            sharedV[2 * blockDim.x] += delz * dzfp;
            sharedV[3 * blockDim.x] += delx * dyfp;
            sharedV[4 * blockDim.x] += delx * dzfp;
            sharedV[5 * blockDim.x] += dely * dzfp;
          }
        }
      }
  }

  __syncthreads();

  if(ii < _inum) {
    F_FLOAT* my_f;

    if(_collect_forces_later) {
      ENERGY_FLOAT* buffer = (ENERGY_FLOAT*) _buffer;

      if(eflag) {
        buffer = &buffer[1 * gridDim.x * gridDim.y];

        if(coul_type != COUL_NONE)
          buffer = &buffer[1 * gridDim.x * gridDim.y];
      }

      if(vflag) {
        buffer = &buffer[6 * gridDim.x * gridDim.y];
      }

      my_f = (F_FLOAT*) buffer;
      my_f += i;
      *my_f = fxtmp;
      my_f += _nmax;
      *my_f = fytmp;
      my_f += _nmax;
      *my_f = fztmp;
    } else {
      my_f = _f + i;
      *my_f += fxtmp;
      my_f += _nmax;
      *my_f += fytmp;
      my_f += _nmax;
      *my_f += fztmp;
    }
  }

  __syncthreads();

  if(eflag) {
    sharedE[0] = evdwl;

    if(coul_type != COUL_NONE)
      sharedECoul[0] = ecoul;
  }

  if(eflag_atom && i < _nlocal) {
    if(coul_type != COUL_NONE)
      _eatom[i] += evdwl + ecoul;
    else
      _eatom[i] += evdwl;
  }

  if(vflag_atom && i < _nlocal) {
    _vatom[i]         += ENERGY_F(0.5) * sharedV[0 * blockDim.x];
    _vatom[i + _nmax]   += ENERGY_F(0.5) * sharedV[1 * blockDim.x];
    _vatom[i + 2 * _nmax] += ENERGY_F(0.5) * sharedV[2 * blockDim.x];
    _vatom[i + 3 * _nmax] += ENERGY_F(0.5) * sharedV[3 * blockDim.x];
    _vatom[i + 4 * _nmax] += ENERGY_F(0.5) * sharedV[4 * blockDim.x];
    _vatom[i + 5 * _nmax] += ENERGY_F(0.5) * sharedV[5 * blockDim.x];
  }

  if(vflag || eflag) PairVirialCompute_A_Kernel(eflag, vflag, coul_type != COUL_NONE ? 1 : 0);
}

template <const PAIR_FORCES pair_type, const COUL_FORCES coul_type, const unsigned int extended_data>
__global__ void Pair_Kernel_BpA(int eflag, int vflag, int eflag_atom, int vflag_atom)
{
  int ii = (blockIdx.x * gridDim.y + blockIdx.y);

  if(ii >= _inum)
    return;

  ENERGY_FLOAT evdwl = ENERGY_F(0.0);
  ENERGY_FLOAT ecoul = ENERGY_F(0.0);
  F_FLOAT3* sharedVirial1;
  F_FLOAT3* sharedVirial2;
  F_FLOAT* sharedEnergy;
  F_FLOAT* sharedEnergyCoul;

  F_FLOAT3* sharedForce = (F_FLOAT3*) &sharedmem[0];

  if(vflag) {
    sharedVirial1 = &sharedForce[64];
    sharedVirial2 = &sharedVirial1[64];
  } else {
    sharedVirial1 = &sharedForce[0];
    sharedVirial2 = &sharedVirial1[0];
  }

  if(eflag) {
    if(vflag || vflag_atom)
      sharedEnergy = (F_FLOAT*) &sharedVirial2[64];
    else
      sharedEnergy = (F_FLOAT*) &sharedForce[64];

    if(coul_type != COUL_NONE)
      sharedEnergyCoul = (F_FLOAT*) &sharedEnergy[64];

  }

  F_FLOAT3 partialForce = { F_F(0.0),  F_F(0.0),  F_F(0.0) };
  F_FLOAT3 partialVirial1 = {  F_F(0.0),  F_F(0.0),  F_F(0.0) };
  F_FLOAT3 partialVirial2 = {  F_F(0.0),  F_F(0.0),  F_F(0.0) };

  X_FLOAT xtmp, ytmp, ztmp;
  X_FLOAT4 myxtype;
  F_FLOAT delx, dely, delz;
  F_FLOAT factor_lj, factor_coul;
  F_FLOAT fpair;
  F_FLOAT qtmp;
  int itype, jnum, i, j;
  int* jlist;

  i = _ilist[ii];

  myxtype = fetchXType(i);

  xtmp = myxtype.x;
  ytmp = myxtype.y;
  ztmp = myxtype.z;
  itype = static_cast <int>(myxtype.w);

  if(coul_type != COUL_NONE)
    qtmp = fetchQ(i);

  jnum = _numneigh[i];

  jlist = &_neighbors[i * _maxneighbors];
  __syncthreads();

  for(int jj = threadIdx.x; jj < jnum + blockDim.x; jj += blockDim.x) {
    if(jj < jnum) {
      fpair = F_F(0.0);
      j = jlist[jj];
      factor_lj =  _special_lj[sbmask(j)];

      if(coul_type != COUL_NONE)
        factor_coul = _special_coul[sbmask(j)];

      j &= NEIGHMASK;

      myxtype = fetchXType(j);

      delx = xtmp - myxtype.x;
      dely = ytmp - myxtype.y;
      delz = ztmp - myxtype.z;
      int jtype = static_cast <int>(myxtype.w);

      const F_FLOAT rsq = delx * delx + dely * dely + delz * delz;

      bool in_cutoff = rsq < (_cutsq_global > X_F(0.0) ? _cutsq_global : _cutsq[itype * _cuda_ntypes + jtype]);
      bool in_coul_cutoff;

      if(in_cutoff) {
        switch(pair_type) {
          case PAIR_BORN:
            fpair += PairBornCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_BUCK:
            fpair += PairBuckCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_CG_CMM:
            fpair += PairLJSDKCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ_CHARMM:
            fpair += PairLJCharmmCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ_CLASS2:
            fpair += PairLJClass2Cuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ_CUT:
            fpair += PairLJCutCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ_EXPAND:
            fpair += PairLJExpandCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ_GROMACS:
            fpair += PairLJGromacsCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ_SMOOTH:
            fpair += PairLJSmoothCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ96_CUT:
            fpair += PairLJ96CutCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_MORSE_R6:
            fpair += PairMorseR6Cuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_MORSE:
            fpair += PairMorseCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;
        }
      }

      if(coul_type != COUL_NONE) {
        const F_FLOAT qiqj = qtmp * fetchQ(j);

        if(qiqj * qiqj > (1e-8f)) {
          in_coul_cutoff =
            rsq < (_cut_coulsq_global > X_F(0.0) ? _cut_coulsq_global : _cut_coulsq[itype * _cuda_ntypes + jtype]);

          if(in_coul_cutoff) {
            switch(coul_type) {
              case COUL_CHARMM:
                fpair += CoulCharmmCuda_Eval(rsq, factor_coul, eflag, ecoul, qiqj);
                break;

              case COUL_CHARMM_IMPLICIT:
                fpair += CoulCharmmImplicitCuda_Eval(rsq, factor_coul, eflag, ecoul, qiqj);
                break;

              case COUL_GROMACS:
                fpair += CoulGromacsCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_coul, eflag, ecoul, qiqj);
                break;

              case COUL_LONG: {
                const F_FLOAT r2inv = F_F(1.0) / rsq;
                const F_FLOAT r = _RSQRT_(r2inv);
                const F_FLOAT grij = _g_ewald * r;
                const F_FLOAT expm2 = _EXP_(-grij * grij);
                const F_FLOAT t = F_F(1.0) / (F_F(1.0) + EWALD_P * grij);
                const F_FLOAT erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
                const F_FLOAT prefactor = _qqrd2e * qiqj * (F_F(1.0) / r);
                F_FLOAT forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);

                if(factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;

                if(eflag) {
                  ecoul += prefactor * erfc;

                  if(factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
                }

                fpair += forcecoul * r2inv;
              }
              break;

              case COUL_DEBYE: {
                const F_FLOAT r2inv = F_F(1.0) / rsq;
                const X_FLOAT r = _RSQRT_(r2inv);
                const X_FLOAT rinv = F_F(1.0) / r;
                const F_FLOAT screening = _EXP_(-_kappa * r);
                F_FLOAT forcecoul = factor_coul * _qqrd2e * qiqj * screening ;

                if(eflag) {
                  ecoul += forcecoul * rinv;
                }

                forcecoul *= (_kappa + rinv);
                fpair += forcecoul * r2inv;
              }
              break;

              case COUL_CUT: {
                const F_FLOAT forcecoul = factor_coul * _qqrd2e * qiqj * _RSQRT_(rsq);

                if(eflag) {
                  ecoul += forcecoul;
                }

                fpair += forcecoul * (F_F(1.0) / rsq);
              }
              break;


            }
          }
        }
      }



      if(in_cutoff || in_coul_cutoff) {
        F_FLOAT dxfp, dyfp, dzfp;
        partialForce.x += dxfp = delx * fpair;
        partialForce.y += dyfp = dely * fpair;
        partialForce.z += dzfp = delz * fpair;

        if(vflag) {
          partialVirial1.x += delx * dxfp;
          partialVirial1.y += dely * dyfp;
          partialVirial1.z += delz * dzfp;
          partialVirial2.x += delx * dyfp;
          partialVirial2.y += delx * dzfp;
          partialVirial2.z += dely * dzfp;
        }
      }
    }
  }

  if(eflag) {
    sharedEnergy[threadIdx.x] = evdwl;

    if(coul_type != COUL_NONE)
      sharedEnergyCoul[threadIdx.x] = ecoul;
  }

  sharedForce[threadIdx.x] = partialForce;

  if(vflag) {
    sharedVirial1[threadIdx.x] = partialVirial1;
    sharedVirial2[threadIdx.x] = partialVirial2;
  }

  __syncthreads();


  for(unsigned int s = blockDim.x >> 1; s > 0; s >>= 1) {

    if(threadIdx.x < s) {
      sharedForce[ threadIdx.x ].x += sharedForce[ threadIdx.x + s ].x;
      sharedForce[ threadIdx.x ].y += sharedForce[ threadIdx.x + s ].y;
      sharedForce[ threadIdx.x ].z += sharedForce[ threadIdx.x + s ].z;

      if(vflag) {
        sharedVirial1[ threadIdx.x ].x += sharedVirial1[ threadIdx.x + s ].x;
        sharedVirial1[ threadIdx.x ].y += sharedVirial1[ threadIdx.x + s ].y;
        sharedVirial1[ threadIdx.x ].z += sharedVirial1[ threadIdx.x + s ].z;

        sharedVirial2[ threadIdx.x ].x += sharedVirial2[ threadIdx.x + s ].x;
        sharedVirial2[ threadIdx.x ].y += sharedVirial2[ threadIdx.x + s ].y;
        sharedVirial2[ threadIdx.x ].z += sharedVirial2[ threadIdx.x + s ].z;
      }

      if(eflag) {
        sharedEnergy[ threadIdx.x ] += sharedEnergy[ threadIdx.x + s ];

        if(coul_type != COUL_NONE)
          sharedEnergyCoul[ threadIdx.x ] += sharedEnergyCoul[ threadIdx.x + s ];
      }
    }

    __syncthreads();
  }

  if(threadIdx.x == 0) {

    ENERGY_FLOAT* buffer = (ENERGY_FLOAT*) _buffer;

    if(eflag) {
      ENERGY_FLOAT tmp_evdwl;
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 0 * gridDim.x * gridDim.y] = tmp_evdwl = ENERGY_F(0.5) * sharedEnergy[0];

      if(eflag_atom)
        _eatom[i] = tmp_evdwl;

      buffer = &buffer[gridDim.x * gridDim.y];

      if(coul_type != COUL_NONE) {
        buffer[blockIdx.x * gridDim.y + blockIdx.y + 0 * gridDim.x * gridDim.y] = tmp_evdwl = ENERGY_F(0.5) * sharedEnergyCoul[0];

        if(eflag_atom)
          _eatom[i] += tmp_evdwl;

        buffer = &buffer[gridDim.x * gridDim.y];
      }
    }

    if(vflag) {
      ENERGY_FLOAT tmp;
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 0 * gridDim.x * gridDim.y] = tmp = ENERGY_F(0.5) * sharedVirial1[0].x;

      if(vflag_atom) _vatom[i + 0 * _nmax] = tmp;

      buffer[blockIdx.x * gridDim.y + blockIdx.y + 1 * gridDim.x * gridDim.y] = tmp = ENERGY_F(0.5) * sharedVirial1[0].y;

      if(vflag_atom) _vatom[i + 1 * _nmax] = tmp;

      buffer[blockIdx.x * gridDim.y + blockIdx.y + 2 * gridDim.x * gridDim.y] = tmp = ENERGY_F(0.5) * sharedVirial1[0].z;

      if(vflag_atom) _vatom[i + 2 * _nmax] = tmp;

      buffer[blockIdx.x * gridDim.y + blockIdx.y + 3 * gridDim.x * gridDim.y] = tmp = ENERGY_F(0.5) * sharedVirial2[0].x;

      if(vflag_atom) _vatom[i + 3 * _nmax] = tmp;

      buffer[blockIdx.x * gridDim.y + blockIdx.y + 4 * gridDim.x * gridDim.y] = tmp = ENERGY_F(0.5) * sharedVirial2[0].y;

      if(vflag_atom) _vatom[i + 4 * _nmax] = tmp;

      buffer[blockIdx.x * gridDim.y + blockIdx.y + 5 * gridDim.x * gridDim.y] = tmp = ENERGY_F(0.5) * sharedVirial2[0].z;

      if(vflag_atom) _vatom[i + 5 * _nmax] = tmp;

      buffer = &buffer[6 * gridDim.x * gridDim.y];
    }

    F_FLOAT* my_f;

    if(_collect_forces_later) {
      my_f = (F_FLOAT*) buffer;
      my_f += i;
      *my_f = sharedForce[0].x;
      my_f += _nmax;
      *my_f = sharedForce[0].y;
      my_f += _nmax;
      *my_f = sharedForce[0].z;
    } else {
      my_f = _f + i;
      *my_f += sharedForce[0].x;
      my_f += _nmax;
      *my_f += sharedForce[0].y;
      my_f += _nmax;
      *my_f += sharedForce[0].z;
    }
  }
}


template <const PAIR_FORCES pair_type, const COUL_FORCES coul_type, const unsigned int extended_data>
__global__ void Pair_Kernel_TpA_opt(int eflag, int vflag, int eflag_atom, int vflag_atom, int comm_phase)
{
  ENERGY_FLOAT evdwl = ENERGY_F(0.0);
  ENERGY_FLOAT ecoul = ENERGY_F(0.0);

  ENERGY_FLOAT* sharedE;
  ENERGY_FLOAT* sharedECoul;
  ENERGY_FLOAT* sharedV = &sharedmem[threadIdx.x];

  if(eflag || eflag_atom) {
    sharedE = &sharedmem[threadIdx.x];
    sharedE[0] = ENERGY_F(0.0);
    sharedV += blockDim.x;

    if(coul_type != COUL_NONE) {
      sharedECoul = sharedE + blockDim.x;
      sharedECoul[0] = ENERGY_F(0.0);
      sharedV += blockDim.x;
    }
  }

  if(vflag || vflag_atom) {
    sharedV[0 * blockDim.x] = ENERGY_F(0.0);
    sharedV[1 * blockDim.x] = ENERGY_F(0.0);
    sharedV[2 * blockDim.x] = ENERGY_F(0.0);
    sharedV[3 * blockDim.x] = ENERGY_F(0.0);
    sharedV[4 * blockDim.x] = ENERGY_F(0.0);
    sharedV[5 * blockDim.x] = ENERGY_F(0.0);
  }

  int ii = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  X_FLOAT xtmp, ytmp, ztmp;
  X_FLOAT4 myxtype;
  F_FLOAT fxtmp, fytmp, fztmp, fpair;
  F_FLOAT delx, dely, delz;
  F_FLOAT factor_lj, factor_coul;
  F_FLOAT qtmp;
  int itype, i, j;
  int jnum = 0;
  int* jlist;

  if(ii < (comm_phase < 2 ? _inum : _inum_border[0])) {
    i = comm_phase < 2 ? _ilist[ii] : _ilist_border[ii] ;

    myxtype = fetchXType(i);
    myxtype = _x_type[i];
    xtmp = myxtype.x;
    ytmp = myxtype.y;
    ztmp = myxtype.z;
    itype = static_cast <int>(myxtype.w);


    fxtmp = F_F(0.0);
    fytmp = F_F(0.0);
    fztmp = F_F(0.0);

    if(coul_type != COUL_NONE)
      qtmp = fetchQ(i);

    jnum = comm_phase == 0 ? _numneigh[i] : (comm_phase == 1 ? _numneigh_inner[i] : _numneigh_border[ii]);


    jlist = comm_phase == 0 ? &_neighbors[i] : (comm_phase == 1 ? &_neighbors_inner[i] : &_neighbors_border[ii]);
  }

  __syncthreads();

  for(int jj = 0; jj < jnum; jj++) {
    if(ii < (comm_phase < 2 ? _inum : _inum_border[0]))
      if(jj < jnum) {
        fpair = F_F(0.0);
        j = jlist[jj * _nlocal];

        factor_lj = j < _nall ? F_F(1.0) : _special_lj[j / _nall];

        if(coul_type != COUL_NONE)
          factor_coul = j < _nall ? F_F(1.0) : _special_coul[j / _nall];

        j = j < _nall ? j : j % _nall;

        myxtype = fetchXType(j);
        delx = xtmp - myxtype.x;
        dely = ytmp - myxtype.y;
        delz = ztmp - myxtype.z;
        int jtype = static_cast <int>(myxtype.w);


        const F_FLOAT rsq = delx * delx + dely * dely + delz * delz;

        bool in_cutoff = rsq < (_cutsq_global > X_F(0.0) ? _cutsq_global : _cutsq[itype * _cuda_ntypes + jtype]);

        if(in_cutoff) {
          switch(pair_type) {
            case PAIR_BORN:
              fpair += PairBornCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_BUCK:
              fpair += PairBuckCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_CG_CMM:
              fpair += PairLJSDKCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ_CHARMM:
              fpair += PairLJCharmmCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ_CLASS2:
              fpair += PairLJClass2Cuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ_CUT:
              fpair += PairLJCutCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ_EXPAND:
              fpair += PairLJExpandCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ_GROMACS:
              fpair += PairLJGromacsCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ_SMOOTH:
              fpair += PairLJSmoothCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_LJ96_CUT:
              fpair += PairLJ96CutCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_MORSE_R6:
              fpair += PairMorseR6Cuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;

            case PAIR_MORSE:
              fpair += PairMorseCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
              break;
          }
        }

        if(coul_type != COUL_NONE) {
          const F_FLOAT qiqj = qtmp * fetchQ(j);

          if(qiqj * qiqj > 1e-8) {
            const bool in_coul_cutoff =
              rsq < (_cut_coulsq_global > X_F(0.0) ? _cut_coulsq_global : _cut_coulsq[itype * _cuda_ntypes + jtype]);

            if(in_coul_cutoff) {
              switch(coul_type) {
                case COUL_CHARMM:
                  fpair += CoulCharmmCuda_Eval(rsq, factor_coul, eflag, ecoul, qiqj);
                  break;

                case COUL_CHARMM_IMPLICIT:
                  fpair += CoulCharmmImplicitCuda_Eval(rsq, factor_coul, eflag, ecoul, qiqj);
                  break;

                case COUL_CUT: {
                  const F_FLOAT forcecoul = factor_coul * _qqrd2e * qiqj * _RSQRT_(rsq);

                  if(eflag) {
                    ecoul += forcecoul;
                  }

                  fpair += forcecoul * (F_F(1.0) / rsq);
                }
                break;

                case COUL_DEBYE: {
                  const F_FLOAT r2inv = F_F(1.0) / rsq;
                  const X_FLOAT r = _RSQRT_(r2inv);
                  const X_FLOAT rinv = F_F(1.0) / r;
                  const F_FLOAT screening = _EXP_(-_kappa * r);
                  F_FLOAT forcecoul = factor_coul * _qqrd2e * qiqj * screening ;

                  if(eflag) {
                    ecoul += forcecoul * rinv;
                  }

                  forcecoul *= (_kappa + rinv);
                  fpair += forcecoul * r2inv;
                }
                break;

                case COUL_GROMACS:
                  fpair += CoulGromacsCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_coul, eflag, ecoul, qiqj);
                  break;

                case COUL_LONG: {
                  const F_FLOAT r2inv = F_F(1.0) / rsq;
                  const F_FLOAT r = _RSQRT_(r2inv);
                  const F_FLOAT grij = _g_ewald * r;
                  const F_FLOAT expm2 = _EXP_(-grij * grij);
                  const F_FLOAT t = F_F(1.0) / (F_F(1.0) + EWALD_P * grij);
                  const F_FLOAT erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
                  const F_FLOAT prefactor = _qqrd2e * qiqj * (F_F(1.0) / r);
                  F_FLOAT forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);

                  if(factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;

                  if(eflag) {
                    ecoul += prefactor * erfc;

                    if(factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
                  }

                  fpair += forcecoul * r2inv;
                }
                break;

              }
            }

            in_cutoff = in_cutoff || in_coul_cutoff;
          }
        }


        if(in_cutoff) {
          F_FLOAT dxfp, dyfp, dzfp;
          fxtmp += dxfp = delx * fpair;
          fytmp += dyfp = dely * fpair;
          fztmp += dzfp = delz * fpair;

          if(vflag) {
            sharedV[0 * blockDim.x] += delx * dxfp;
            sharedV[1 * blockDim.x] += dely * dyfp;
            sharedV[2 * blockDim.x] += delz * dzfp;
            sharedV[3 * blockDim.x] += delx * dyfp;
            sharedV[4 * blockDim.x] += delx * dzfp;
            sharedV[5 * blockDim.x] += dely * dzfp;
          }
        }
      }
  }

  __syncthreads();

  if(ii < (comm_phase < 2 ? _inum : _inum_border[0])) {
    F_FLOAT* my_f;

    if(_collect_forces_later) {
      ENERGY_FLOAT* buffer = (ENERGY_FLOAT*) _buffer;

      if(eflag) {
        buffer = &buffer[1 * gridDim.x * gridDim.y];

        if(coul_type != COUL_NONE)
          buffer = &buffer[1 * gridDim.x * gridDim.y];
      }

      if(vflag) {
        buffer = &buffer[6 * gridDim.x * gridDim.y];
      }

      my_f = (F_FLOAT*) buffer;
      my_f += i;
      *my_f = fxtmp;
      my_f += _nmax;
      *my_f = fytmp;
      my_f += _nmax;
      *my_f = fztmp;
    } else {
      my_f = _f + i;
      *my_f += fxtmp;
      my_f += _nmax;
      *my_f += fytmp;
      my_f += _nmax;
      *my_f += fztmp;
    }
  }

  __syncthreads();

  if(eflag) {
    sharedE[0] = evdwl;

    if(coul_type != COUL_NONE)
      sharedECoul[0] = ecoul;
  }

  if(eflag_atom && i < _nlocal) {
    if(coul_type != COUL_NONE)
      _eatom[i] += evdwl + ecoul;
    else
      _eatom[i] += evdwl;
  }

  if(vflag_atom && i < _nlocal) {
    _vatom[i]         += ENERGY_F(0.5) * sharedV[0 * blockDim.x];
    _vatom[i + _nmax]   += ENERGY_F(0.5) * sharedV[1 * blockDim.x];
    _vatom[i + 2 * _nmax] += ENERGY_F(0.5) * sharedV[2 * blockDim.x];
    _vatom[i + 3 * _nmax] += ENERGY_F(0.5) * sharedV[3 * blockDim.x];
    _vatom[i + 4 * _nmax] += ENERGY_F(0.5) * sharedV[4 * blockDim.x];
    _vatom[i + 5 * _nmax] += ENERGY_F(0.5) * sharedV[5 * blockDim.x];
  }

  if(vflag || eflag) PairVirialCompute_A_Kernel(eflag, vflag, coul_type != COUL_NONE ? 1 : 0);
}

template <const PAIR_FORCES pair_type, const COUL_FORCES coul_type, const unsigned int extended_data>
__global__ void Pair_Kernel_BpA_opt(int eflag, int vflag, int eflag_atom, int vflag_atom, int comm_phase)
{
  int ii = (blockIdx.x * gridDim.y + blockIdx.y);

  if(ii >= (comm_phase < 2 ? _inum : _inum_border[0]))
    return;

  ENERGY_FLOAT evdwl = ENERGY_F(0.0);
  ENERGY_FLOAT ecoul = ENERGY_F(0.0);
  F_FLOAT3* sharedVirial1;
  F_FLOAT3* sharedVirial2;
  F_FLOAT* sharedEnergy;
  F_FLOAT* sharedEnergyCoul;

  F_FLOAT3* sharedForce = (F_FLOAT3*) &sharedmem[0];

  if(vflag) {
    sharedVirial1 = &sharedForce[64];
    sharedVirial2 = &sharedVirial1[64];
  } else {
    sharedVirial1 = &sharedForce[0];
    sharedVirial2 = &sharedVirial1[0];
  }

  if(eflag) {
    if(vflag || vflag_atom)
      sharedEnergy = (F_FLOAT*) &sharedVirial2[64];
    else
      sharedEnergy = (F_FLOAT*) &sharedForce[64];

    if(coul_type != COUL_NONE)
      sharedEnergyCoul = (F_FLOAT*) &sharedEnergy[64];

  }

  F_FLOAT3 partialForce = { F_F(0.0),  F_F(0.0),  F_F(0.0) };
  F_FLOAT3 partialVirial1 = {  F_F(0.0),  F_F(0.0),  F_F(0.0) };
  F_FLOAT3 partialVirial2 = {  F_F(0.0),  F_F(0.0),  F_F(0.0) };

  X_FLOAT xtmp, ytmp, ztmp;
  X_FLOAT4 myxtype;
  F_FLOAT delx, dely, delz;
  F_FLOAT factor_lj, factor_coul;
  F_FLOAT fpair;
  F_FLOAT qtmp;
  int itype, jnum, i, j;
  int* jlist;

  i = comm_phase < 2 ? _ilist[ii] : _ilist_border[ii];

  myxtype = fetchXType(i);

  xtmp = myxtype.x;
  ytmp = myxtype.y;
  ztmp = myxtype.z;
  itype = static_cast <int>(myxtype.w);

  if(coul_type != COUL_NONE)
    qtmp = fetchQ(i);

  jnum = comm_phase == 0 ? _numneigh[i] : (comm_phase == 1 ? _numneigh_inner[i] : _numneigh_border[ii]);

  jlist = comm_phase == 0 ? &_neighbors[i * _maxneighbors] : (comm_phase == 1 ? &_neighbors_inner[i * _maxneighbors] : &_neighbors_border[ii * _maxneighbors]);
  __syncthreads();

  for(int jj = threadIdx.x; jj < jnum + blockDim.x; jj += blockDim.x) {
    if(jj < jnum) {
      fpair = F_F(0.0);
      j = jlist[jj];
      factor_lj   = j < _nall ? F_F(1.0) : _special_lj[j / _nall];

      if(coul_type != COUL_NONE)
        factor_coul = j < _nall ? F_F(1.0) : _special_coul[j / _nall];

      j 			= j < _nall ? j : j % _nall;

      myxtype = fetchXType(j);

      delx = xtmp - myxtype.x;
      dely = ytmp - myxtype.y;
      delz = ztmp - myxtype.z;
      int jtype = static_cast <int>(myxtype.w);

      const F_FLOAT rsq = delx * delx + dely * dely + delz * delz;

      bool in_cutoff = rsq < (_cutsq_global > X_F(0.0) ? _cutsq_global : _cutsq[itype * _cuda_ntypes + jtype]);
      bool in_coul_cutoff;

      if(in_cutoff) {
        switch(pair_type) {
          case PAIR_BORN:
            fpair += PairBornCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_BUCK:
            fpair += PairBuckCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_CG_CMM:
            fpair += PairLJSDKCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ_CHARMM:
            fpair += PairLJCharmmCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ_CLASS2:
            fpair += PairLJClass2Cuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ_CUT:
            fpair += PairLJCutCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ_EXPAND:
            fpair += PairLJExpandCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ_GROMACS:
            fpair += PairLJGromacsCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ_SMOOTH:
            fpair += PairLJSmoothCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_LJ96_CUT:
            fpair += PairLJ96CutCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_MORSE_R6:
            fpair += PairMorseR6Cuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;

          case PAIR_MORSE:
            fpair += PairMorseCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_lj, eflag, evdwl);
            break;
        }
      }

      if(coul_type != COUL_NONE) {
        const F_FLOAT qiqj = qtmp * fetchQ(j);

        if(qiqj * qiqj > (1e-8f)) {
          in_coul_cutoff =
            rsq < (_cut_coulsq_global > X_F(0.0) ? _cut_coulsq_global : _cut_coulsq[itype * _cuda_ntypes + jtype]);

          if(in_coul_cutoff) {
            switch(coul_type) {
              case COUL_CHARMM:
                fpair += CoulCharmmCuda_Eval(rsq, factor_coul, eflag, ecoul, qiqj);
                break;

              case COUL_CHARMM_IMPLICIT:
                fpair += CoulCharmmImplicitCuda_Eval(rsq, factor_coul, eflag, ecoul, qiqj);
                break;

              case COUL_GROMACS:
                fpair += CoulGromacsCuda_Eval(rsq, itype * _cuda_ntypes + jtype, factor_coul, eflag, ecoul, qiqj);
                break;

              case COUL_LONG: {
                const F_FLOAT r2inv = F_F(1.0) / rsq;
                const F_FLOAT r = _RSQRT_(r2inv);
                const F_FLOAT grij = _g_ewald * r;
                const F_FLOAT expm2 = _EXP_(-grij * grij);
                const F_FLOAT t = F_F(1.0) / (F_F(1.0) + EWALD_P * grij);
                const F_FLOAT erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
                const F_FLOAT prefactor = _qqrd2e * qiqj * (F_F(1.0) / r);
                F_FLOAT forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);

                if(factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;

                if(eflag) {
                  ecoul += prefactor * erfc;

                  if(factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
                }

                fpair += forcecoul * r2inv;
              }
              break;

              case COUL_DEBYE: {
                const F_FLOAT r2inv = F_F(1.0) / rsq;
                const X_FLOAT r = _RSQRT_(r2inv);
                const X_FLOAT rinv = F_F(1.0) / r;
                const F_FLOAT screening = _EXP_(-_kappa * r);
                F_FLOAT forcecoul = factor_coul * _qqrd2e * qiqj * screening ;

                if(eflag) {
                  ecoul += forcecoul * rinv;
                }

                forcecoul *= (_kappa + rinv);
                fpair += forcecoul * r2inv;
              }
              break;

              case COUL_CUT: {
                const F_FLOAT forcecoul = factor_coul * _qqrd2e * qiqj * _RSQRT_(rsq);

                if(eflag) {
                  ecoul += forcecoul;
                }

                fpair += forcecoul * (F_F(1.0) / rsq);
              }
              break;


            }
          }
        }
      }



      if(in_cutoff || in_coul_cutoff) {
        F_FLOAT dxfp, dyfp, dzfp;
        partialForce.x += dxfp = delx * fpair;
        partialForce.y += dyfp = dely * fpair;
        partialForce.z += dzfp = delz * fpair;

        if(vflag) {
          partialVirial1.x += delx * dxfp;
          partialVirial1.y += dely * dyfp;
          partialVirial1.z += delz * dzfp;
          partialVirial2.x += delx * dyfp;
          partialVirial2.y += delx * dzfp;
          partialVirial2.z += dely * dzfp;
        }
      }
    }
  }

  if(eflag) {
    sharedEnergy[threadIdx.x] = evdwl;

    if(coul_type != COUL_NONE)
      sharedEnergyCoul[threadIdx.x] = ecoul;
  }

  sharedForce[threadIdx.x] = partialForce;

  if(vflag) {
    sharedVirial1[threadIdx.x] = partialVirial1;
    sharedVirial2[threadIdx.x] = partialVirial2;
  }

  __syncthreads();


  for(unsigned int s = blockDim.x >> 1; s > 0; s >>= 1) {

    if(threadIdx.x < s) {
      sharedForce[ threadIdx.x ].x += sharedForce[ threadIdx.x + s ].x;
      sharedForce[ threadIdx.x ].y += sharedForce[ threadIdx.x + s ].y;
      sharedForce[ threadIdx.x ].z += sharedForce[ threadIdx.x + s ].z;

      if(vflag) {
        sharedVirial1[ threadIdx.x ].x += sharedVirial1[ threadIdx.x + s ].x;
        sharedVirial1[ threadIdx.x ].y += sharedVirial1[ threadIdx.x + s ].y;
        sharedVirial1[ threadIdx.x ].z += sharedVirial1[ threadIdx.x + s ].z;

        sharedVirial2[ threadIdx.x ].x += sharedVirial2[ threadIdx.x + s ].x;
        sharedVirial2[ threadIdx.x ].y += sharedVirial2[ threadIdx.x + s ].y;
        sharedVirial2[ threadIdx.x ].z += sharedVirial2[ threadIdx.x + s ].z;
      }

      if(eflag) {
        sharedEnergy[ threadIdx.x ] += sharedEnergy[ threadIdx.x + s ];

        if(coul_type != COUL_NONE)
          sharedEnergyCoul[ threadIdx.x ] += sharedEnergyCoul[ threadIdx.x + s ];
      }
    }

    __syncthreads();
  }

  if(threadIdx.x == 0) {

    ENERGY_FLOAT* buffer = (ENERGY_FLOAT*) _buffer;

    if(eflag) {
      ENERGY_FLOAT tmp_evdwl;
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 0 * gridDim.x * gridDim.y] = tmp_evdwl = ENERGY_F(0.5) * sharedEnergy[0];

      if(eflag_atom)
        _eatom[i] = tmp_evdwl;

      buffer = &buffer[gridDim.x * gridDim.y];

      if(coul_type != COUL_NONE) {
        buffer[blockIdx.x * gridDim.y + blockIdx.y + 0 * gridDim.x * gridDim.y] = tmp_evdwl = ENERGY_F(0.5) * sharedEnergyCoul[0];

        if(eflag_atom)
          _eatom[i] += tmp_evdwl;

        buffer = &buffer[gridDim.x * gridDim.y];
      }
    }

    if(vflag) {
      ENERGY_FLOAT tmp;
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 0 * gridDim.x * gridDim.y] = tmp = ENERGY_F(0.5) * sharedVirial1[0].x;

      if(vflag_atom) _vatom[i + 0 * _nmax] = tmp;

      buffer[blockIdx.x * gridDim.y + blockIdx.y + 1 * gridDim.x * gridDim.y] = tmp = ENERGY_F(0.5) * sharedVirial1[0].y;

      if(vflag_atom) _vatom[i + 1 * _nmax] = tmp;

      buffer[blockIdx.x * gridDim.y + blockIdx.y + 2 * gridDim.x * gridDim.y] = tmp = ENERGY_F(0.5) * sharedVirial1[0].z;

      if(vflag_atom) _vatom[i + 2 * _nmax] = tmp;

      buffer[blockIdx.x * gridDim.y + blockIdx.y + 3 * gridDim.x * gridDim.y] = tmp = ENERGY_F(0.5) * sharedVirial2[0].x;

      if(vflag_atom) _vatom[i + 3 * _nmax] = tmp;

      buffer[blockIdx.x * gridDim.y + blockIdx.y + 4 * gridDim.x * gridDim.y] = tmp = ENERGY_F(0.5) * sharedVirial2[0].y;

      if(vflag_atom) _vatom[i + 4 * _nmax] = tmp;

      buffer[blockIdx.x * gridDim.y + blockIdx.y + 5 * gridDim.x * gridDim.y] = tmp = ENERGY_F(0.5) * sharedVirial2[0].z;

      if(vflag_atom) _vatom[i + 5 * _nmax] = tmp;

      buffer = &buffer[6 * gridDim.x * gridDim.y];
    }

    F_FLOAT* my_f;

    if(_collect_forces_later) {
      my_f = (F_FLOAT*) buffer;
      my_f += i;
      *my_f = sharedForce[0].x;
      my_f += _nmax;
      *my_f = sharedForce[0].y;
      my_f += _nmax;
      *my_f = sharedForce[0].z;
    } else {
      my_f = _f + i;
      *my_f += sharedForce[0].x;
      my_f += _nmax;
      *my_f += sharedForce[0].y;
      my_f += _nmax;
      *my_f += sharedForce[0].z;
    }
  }
}

__global__ void Pair_GenerateXType_Kernel()
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nall) {
    X_FLOAT4 xtype;
    xtype.x = _x[i];
    xtype.y = _x[i + _nmax];
    xtype.z = _x[i + 2 * _nmax];
    xtype.w = _type[i];
    _x_type[i] = xtype;
  }

}

__global__ void Pair_GenerateVRadius_Kernel()
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nall) {
    V_FLOAT4 vradius;
    vradius.x = _v[i];
    vradius.y = _v[i + _nmax];
    vradius.z = _v[i + 2 * _nmax];
    vradius.w = _radius[i];
    _v_radius[i] = vradius;
  }
}

__global__ void Pair_GenerateOmegaRmass_Kernel()
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nall) {
    V_FLOAT4 omegarmass;
    omegarmass.x = _omega[i];
    omegarmass.y = _omega[i + _nmax];
    omegarmass.z = _omega[i + 2 * _nmax];
    omegarmass.w = _rmass[i];
    _omega_rmass[i] = omegarmass;
  }
}

__global__ void Pair_RevertXType_Kernel()
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nall) {
    X_FLOAT4 xtype = _x_type[i];
    _x[i] = xtype.x;
    _x[i + _nmax] = xtype.y;
    _x[i + 2 * _nmax] = xtype.z;
    _type[i] = static_cast <int>(xtype.w);
  }

}

__global__ void Pair_BuildXHold_Kernel()
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < _nall) {
    X_FLOAT4 xtype = _x_type[i];
    _xhold[i] = xtype.x;
    _xhold[i + _nmax] = xtype.y;
    _xhold[i + 2 * _nmax] = xtype.z;
  }

}

__global__ void Pair_CollectForces_Kernel(int nperblock, int n)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i >= _nlocal) return;

  ENERGY_FLOAT* buf = (ENERGY_FLOAT*) _buffer;

  F_FLOAT* buf_f = (F_FLOAT*) &buf[nperblock * n];
  F_FLOAT* my_f = _f + i;
  buf_f += i;
  *my_f += * buf_f;
  my_f += _nmax;
  buf_f += _nmax;
  *my_f += * buf_f;
  my_f += _nmax;
  buf_f += _nmax;
  *my_f += * buf_f;
  my_f += _nmax;
}
