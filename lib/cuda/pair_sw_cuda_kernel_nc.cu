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
#define Pi F_F(3.1415926535897932384626433832795)
#define PI Pi
#define PI2 F_F(0.5)*Pi
#define PI4 F_F(0.25)*Pi



__device__ void twobody(int iparam, F_CFLOAT rsq, F_CFLOAT &fforce,
                        int eflag, ENERGY_CFLOAT &eng)
{
  F_CFLOAT r, rp, rq, rainv, expsrainv;

  r = sqrt(rsq);
  rp = pow(r, -params_sw[iparam].powerp);
  rq = pow(r, -params_sw[iparam].powerq);
  rainv = 1.0 / (r - params_sw[iparam].cut);
  expsrainv = exp(params_sw[iparam].sigma * rainv);
  fforce = (params_sw[iparam].c1 * rp - params_sw[iparam].c2 * rq +
            (params_sw[iparam].c3 * rp - params_sw[iparam].c4 * rq) * rainv * rainv * r) * expsrainv / rsq;

  if(eflag) eng += (params_sw[iparam].c5 * rp - params_sw[iparam].c6 * rq) * expsrainv;
}

__device__ void threebody(int paramij, int paramik, int paramijk,
                          F_CFLOAT4 &delr1,
                          F_CFLOAT4 &delr2,
                          F_CFLOAT3 &fj, F_CFLOAT3 &fk, int eflag, ENERGY_CFLOAT &eng)
{
  F_CFLOAT r1, rinvsq1, rainv1, gsrainv1, gsrainvsq1, expgsrainv1;
  F_CFLOAT r2, rinvsq2, rainv2, gsrainv2, gsrainvsq2, expgsrainv2;
  F_CFLOAT rinv12, cs, delcs, delcssq, facexp, facrad, frad1, frad2;
  F_CFLOAT facang, facang12, csfacang, csfac1, csfac2;

  r1 = sqrt(delr1.w);
  rinvsq1 = F_F(1.0) / delr1.w;
  rainv1 = F_F(1.0) / (r1 - params_sw[paramij].cut);
  gsrainv1 = params_sw[paramij].sigma_gamma * rainv1;
  gsrainvsq1 = gsrainv1 * rainv1 / r1;
  expgsrainv1 = exp(gsrainv1);

  r2 = sqrt(delr2.w);
  rinvsq2 = F_F(1.0) / delr2.w;
  rainv2 = F_F(1.0) / (r2 - params_sw[paramik].cut);
  gsrainv2 = params_sw[paramik].sigma_gamma * rainv2;
  gsrainvsq2 = gsrainv2 * rainv2 / r2;
  expgsrainv2 = exp(gsrainv2);

  rinv12 = F_F(1.0) / (r1 * r2);
  cs = (delr1.x * delr2.x + delr1.y * delr2.y + delr1.z * delr2.z) * rinv12;
  delcs = cs - params_sw[paramijk].costheta;
  delcssq = delcs * delcs;

  facexp = expgsrainv1 * expgsrainv2;

  // facrad = sqrt(paramij->lambda_epsilon*paramik->lambda_epsilon) *
  //          facexp*delcssq;

  facrad = params_sw[paramijk].lambda_epsilon * facexp * delcssq;
  frad1 = facrad * gsrainvsq1;
  frad2 = facrad * gsrainvsq2;
  facang = params_sw[paramijk].lambda_epsilon2 * facexp * delcs;
  facang12 = rinv12 * facang;
  csfacang = cs * facang;
  csfac1 = rinvsq1 * csfacang;

  fj.x = delr1.x * (frad1 + csfac1) - delr2.x * facang12;
  fj.y = delr1.y * (frad1 + csfac1) - delr2.y * facang12;
  fj.z = delr1.z * (frad1 + csfac1) - delr2.z * facang12;

  csfac2 = rinvsq2 * csfacang;

  fk.x = delr2.x * (frad2 + csfac2) - delr1.x * facang12;
  fk.y = delr2.y * (frad2 + csfac2) - delr1.y * facang12;
  fk.z = delr2.z * (frad2 + csfac2) - delr1.z * facang12;

  if(eflag) eng += F_F(2.0) * facrad;
}

__device__ void threebody_fj(int paramij, int paramik, int paramijk,
                             F_CFLOAT4 &delr1,
                             F_CFLOAT4 &delr2,
                             F_CFLOAT3 &fj)
{
  F_CFLOAT r1, rinvsq1, rainv1, gsrainv1, gsrainvsq1, expgsrainv1;
  F_CFLOAT r2, rainv2, gsrainv2, expgsrainv2;
  F_CFLOAT rinv12, cs, delcs, delcssq, facexp, facrad, frad1;
  F_CFLOAT facang, facang12, csfacang, csfac1;

  r1 = sqrt(delr1.w);
  rinvsq1 = F_F(1.0) / delr1.w;
  rainv1 = F_F(1.0) / (r1 - params_sw[paramij].cut);
  gsrainv1 = params_sw[paramij].sigma_gamma * rainv1;
  gsrainvsq1 = gsrainv1 * rainv1 / r1;
  expgsrainv1 = exp(gsrainv1);

  r2 = sqrt(delr2.w);
  rainv2 = F_F(1.0) / (r2 - params_sw[paramik].cut);
  gsrainv2 = params_sw[paramik].sigma_gamma * rainv2;
  expgsrainv2 = exp(gsrainv2);

  rinv12 = F_F(1.0) / (r1 * r2);
  cs = (delr1.x * delr2.x + delr1.y * delr2.y + delr1.z * delr2.z) * rinv12;
  delcs = cs - params_sw[paramijk].costheta;
  delcssq = delcs * delcs;

  facexp = expgsrainv1 * expgsrainv2;

  // facrad = sqrt(paramij->lambda_epsilon*paramik->lambda_epsilon) *
  //          facexp*delcssq;

  facrad = params_sw[paramijk].lambda_epsilon * facexp * delcssq;
  frad1 = facrad * gsrainvsq1;
  facang = params_sw[paramijk].lambda_epsilon2 * facexp * delcs;
  facang12 = rinv12 * facang;
  csfacang = cs * facang;
  csfac1 = rinvsq1 * csfacang;

  fj.x = delr1.x * (frad1 + csfac1) - delr2.x * facang12;
  fj.y = delr1.y * (frad1 + csfac1) - delr2.y * facang12;
  fj.z = delr1.z * (frad1 + csfac1) - delr2.z * facang12;
}


__global__ void Pair_SW_Kernel_TpA_RIJ()//F_CFLOAT4* _glob_r_ij,int* _glob_numneigh_red,int* _glob_neighbors_red,int* _glob_neightype_red)
{
  int ii = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(ii >= _nall) return;

  X_CFLOAT4 myxtype;
  F_CFLOAT4 delij;
  F_CFLOAT xtmp, ytmp, ztmp;
  int itype, jnum, i, j;
  int* jlist;
  int neigh_red = 0;
  i = ii;//_ilist[ii];
  myxtype = fetchXType(i);

  xtmp = myxtype.x;
  ytmp = myxtype.y;
  ztmp = myxtype.z;
  itype = map[(static_cast <int>(myxtype.w))];

  jnum = _numneigh[i];
  jlist = &_neighbors[i];

  __syncthreads();

  for(int jj = 0; jj < jnum; jj++) {
    if(jj < jnum) {

      j = jlist[jj * _nall];
      j &= NEIGHMASK;
      myxtype = fetchXType(j);
      delij.x = xtmp - myxtype.x;
      delij.y = ytmp - myxtype.y;
      delij.z = ztmp - myxtype.z;
      int jtype = map[(static_cast <int>(myxtype.w))];
      int iparam_ij = elem2param[(itype * nelements + jtype) * nelements + jtype];
      delij.w = vec3_dot(delij, delij);

      if(delij.w < params_sw[iparam_ij].cutsq) {
        _glob_neighbors_red[i + neigh_red * _nall] = j;
        _glob_neightype_red[i + neigh_red * _nall] = jtype;
        _glob_r_ij[i + neigh_red * _nall] = delij;
        neigh_red++;
      }
    }
  }

  _glob_numneigh_red[i] = neigh_red;
}


template <int eflag, int vflagm>
__global__ void Pair_SW_Kernel_TpA(int eflag_atom, int vflag_atom) //,F_CFLOAT* _glob_zeta_ij,F_CFLOAT4* _glob_r_ij,int* _glob_numneigh_red,int* _glob_neighbors_red,int* _glob_neightype_red)
{
  ENERGY_CFLOAT evdwl = ENERGY_F(0.0);

  ENERGY_CFLOAT* sharedE = &sharedmem[threadIdx.x];
  ENERGY_CFLOAT* sharedV = &sharedmem[threadIdx.x];

  F_CFLOAT* shared_F_F = (F_CFLOAT*) sharedmem;

  if((eflag || eflag_atom) && (vflagm || vflag_atom)) shared_F_F = (F_CFLOAT*) &sharedmem[7 * blockDim.x];
  else if(eflag) shared_F_F = (F_CFLOAT*) &sharedmem[blockDim.x];
  else if(vflagm) shared_F_F = (F_CFLOAT*) &sharedmem[6 * blockDim.x];

  shared_F_F += threadIdx.x;

  if(eflag_atom || eflag) {
    sharedE[0] = ENERGY_F(0.0);
    sharedV += blockDim.x;
  }

  if(vflagm || vflag_atom) {
    sharedV[0 * blockDim.x] = ENERGY_F(0.0);
    sharedV[1 * blockDim.x] = ENERGY_F(0.0);
    sharedV[2 * blockDim.x] = ENERGY_F(0.0);
    sharedV[3 * blockDim.x] = ENERGY_F(0.0);
    sharedV[4 * blockDim.x] = ENERGY_F(0.0);
    sharedV[5 * blockDim.x] = ENERGY_F(0.0);
  }

  int jnum_red = 0;
#define fxtmp shared_F_F[0]
#define fytmp shared_F_F[blockDim.x]
#define fztmp shared_F_F[2*blockDim.x]
  //#define jnum_red (static_cast <int> (shared_F_F[3*blockDim.x]))

  int ii = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  X_CFLOAT4 myxtype_i, myxtype_j, myxtype_k;
  F_CFLOAT4 delij, delik, deljk;
  F_CFLOAT fpair;

  int itype, i, j;
  int* jlist_red;

  if(ii < _inum) {
    i = _ilist[ii];

    if(vflagm)
      myxtype_i = fetchXType(i);

    //itype=map[(static_cast <int> (myxtype_i.w))];
    itype = map[_type[i]];


    fxtmp = F_F(0.0);
    fytmp = F_F(0.0);
    fztmp = F_F(0.0);


    //shared_F_F[3*blockDim.x] = _glob_numneigh_red[i];
    jnum_red = _glob_numneigh_red[i];
    jlist_red = &_glob_neighbors_red[i];
  }

  __syncthreads();
#pragma unroll 1

  for(int jj = 0; jj < jnum_red; jj++) {
    if(i < _nlocal) {
      fpair = F_F(0.0);
      j = jlist_red[jj * _nall];
      j &= NEIGHMASK;

      if(vflagm)
        myxtype_j = fetchXType(j);

      int jtype = _glob_neightype_red[i + jj * _nall];
      delij = _glob_r_ij[i + jj * _nall];

      volatile int iparam_ij = elem2param[(itype * nelements + jtype) * nelements + jtype];
      volatile int iparam_ji = elem2param[(jtype * nelements + itype) * nelements + itype];

      if(delij.w < params_sw[iparam_ij].cutsq) {
        F_CFLOAT dxfp, dyfp, dzfp;
        twobody(iparam_ij, delij.w, fpair, eflag, evdwl);
        fxtmp += dxfp = delij.x * fpair;
        fytmp += dyfp = delij.y * fpair;
        fztmp += dzfp = delij.z * fpair;

        if(vflagm) {
          sharedV[0 * blockDim.x] += delij.x * dxfp;
          sharedV[1 * blockDim.x] += delij.y * dyfp;
          sharedV[2 * blockDim.x] += delij.z * dzfp;
          sharedV[3 * blockDim.x] += delij.x * dyfp;
          sharedV[4 * blockDim.x] += delij.x * dzfp;
          sharedV[5 * blockDim.x] += delij.y * dzfp;
        }






        vec3_scale(F_F(-1.0), delij, delij);

#pragma unroll 1

        for(int kk = jj + 1; kk < jnum_red; kk++) {
          int k = jlist_red[kk * _nall];
          k &= NEIGHMASK;

          if(vflagm)
            myxtype_k = fetchXType(k);

          delik = _glob_r_ij[i + kk * _nall];

          int ktype = _glob_neightype_red[i + kk * _nall];
          int iparam_ik = elem2param[(itype * nelements + ktype) * nelements + ktype];
          int iparam_ijk = elem2param[(itype * nelements + jtype) * nelements + ktype];
          vec3_scale(F_F(-1.0), delik, delik);

          if(delik.w <= params_sw[iparam_ijk].cutsq) {
            F_CFLOAT3 fj, fk;
            threebody(iparam_ij, iparam_ik, iparam_ijk,
                      delij, delik, fj, fk, eflag, evdwl);
            fxtmp -= fj.x + fk.x;
            fytmp -= fj.y + fk.y;
            fztmp -= fj.z + fk.z;

            if(vflagm) {
              sharedV[0 * blockDim.x] -= ENERGY_F(2.0) * myxtype_i.x * (fj.x + fk.x);
              sharedV[1 * blockDim.x] -= ENERGY_F(2.0) * myxtype_i.y * (fj.y + fk.y);
              sharedV[2 * blockDim.x] -= ENERGY_F(2.0) * myxtype_i.z * (fj.z + fk.z);
              sharedV[3 * blockDim.x] -= ENERGY_F(2.0) * myxtype_i.x * (fj.y + fk.y);
              sharedV[4 * blockDim.x] -= ENERGY_F(2.0) * myxtype_i.x * (fj.z + fk.z);
              sharedV[5 * blockDim.x] -= ENERGY_F(2.0) * myxtype_i.y * (fj.z + fk.z);

              sharedV[0 * blockDim.x] += ENERGY_F(2.0) * myxtype_j.x * fj.x;
              sharedV[1 * blockDim.x] += ENERGY_F(2.0) * myxtype_j.y * fj.y;
              sharedV[2 * blockDim.x] += ENERGY_F(2.0) * myxtype_j.z * fj.z;
              sharedV[3 * blockDim.x] += ENERGY_F(2.0) * myxtype_j.x * fj.y;
              sharedV[4 * blockDim.x] += ENERGY_F(2.0) * myxtype_j.x * fj.z;
              sharedV[5 * blockDim.x] += ENERGY_F(2.0) * myxtype_j.y * fj.z;

              sharedV[0 * blockDim.x] += ENERGY_F(2.0) * myxtype_k.x * fk.x;
              sharedV[1 * blockDim.x] += ENERGY_F(2.0) * myxtype_k.y * fk.y;
              sharedV[2 * blockDim.x] += ENERGY_F(2.0) * myxtype_k.z * fk.z;
              sharedV[3 * blockDim.x] += ENERGY_F(2.0) * myxtype_k.x * fk.y;
              sharedV[4 * blockDim.x] += ENERGY_F(2.0) * myxtype_k.x * fk.z;
              sharedV[5 * blockDim.x] += ENERGY_F(2.0) * myxtype_k.y * fk.z;
            }
          }
        }

        int j_jnum_red = _glob_numneigh_red[j];
        int* j_jlist_red = &_glob_neighbors_red[j];

        int j_ii = 0;

        //#pragma unroll 1
        for(int j_kk = 0; j_kk < j_jnum_red; j_kk++) {
          if(j_jlist_red[j_kk * _nall] == i) j_ii = j_kk;
        }

#pragma unroll 1

        for(int kk = 0; kk < j_jnum_red; kk++) {
          if(j_ii == kk) continue;

          int k = j_jlist_red[kk * _nall];
          k &= NEIGHMASK;
          deljk = _glob_r_ij[j + kk * _nall];
          vec3_scale(F_F(-1.0), deljk, deljk);
          int ktype = _glob_neightype_red[j + kk * _nall];

          int iparam_ji = elem2param[(jtype * nelements + itype) * nelements + itype];
          int iparam_jk = elem2param[(jtype * nelements + ktype) * nelements + ktype];
          int iparam_jik = elem2param[(jtype * nelements + itype) * nelements + ktype];


          vec3_scale(F_F(-1.0), delij, delij);

          if(deljk.w <= params_sw[iparam_jik].cutsq) {
            F_CFLOAT3 fj;

            threebody_fj(iparam_ji, iparam_jk, iparam_jik,
                         delij, deljk, fj);
            fxtmp += fj.x;
            fytmp += fj.y;
            fztmp += fj.z;

          }

          vec3_scale(F_F(-1.0), delij, delij);
        }
      }
    }

  }

  __syncthreads();

  if(ii < _inum) {
    F_CFLOAT* my_f;

    if(_collect_forces_later) {
      ENERGY_CFLOAT* buffer = (ENERGY_CFLOAT*) _buffer;

      if(eflag) {
        buffer = &buffer[1 * gridDim.x * gridDim.y];
      }

      if(vflagm) {
        buffer = &buffer[6 * gridDim.x * gridDim.y];
      }

      my_f = (F_CFLOAT*) buffer;
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
  }

  if(eflag_atom && i < _nlocal) {
    _eatom[i] = ENERGY_F(0.5) * evdwl;
  }

  if(vflag_atom && i < _nlocal) {
    _vatom[i]         = ENERGY_F(0.5) * sharedV[0 * blockDim.x];
    _vatom[i + _nmax]   = ENERGY_F(0.5) * sharedV[1 * blockDim.x];
    _vatom[i + 2 * _nmax] = ENERGY_F(0.5) * sharedV[2 * blockDim.x];
    _vatom[i + 3 * _nmax] = ENERGY_F(0.5) * sharedV[3 * blockDim.x];
    _vatom[i + 4 * _nmax] = ENERGY_F(0.5) * sharedV[4 * blockDim.x];
    _vatom[i + 5 * _nmax] = ENERGY_F(0.5) * sharedV[5 * blockDim.x];
  }

  if(vflagm && eflag) PairVirialCompute_A_Kernel_Template<1, 1>();
  else if(eflag) PairVirialCompute_A_Kernel_Template<1, 0>();
  else if(vflagm) PairVirialCompute_A_Kernel_Template<0, 1>();

#undef fxtmp
#undef fytmp
#undef fztmp
  //#undef jnum_red
}
