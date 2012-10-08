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


__global__ void PairGranHookeCuda_Kernel(int eflag, int vflag, int eflag_atom, int vflag_atom, int** firstneight, int* binned_id
    , F_FLOAT kn, F_FLOAT gamman, F_FLOAT gammat, F_FLOAT xmu)
{
  ENERGY_FLOAT evdwl = ENERGY_F(0.0);

  ENERGY_FLOAT* sharedE;
  ENERGY_FLOAT* sharedV;

  if(eflag || eflag_atom) {
    sharedE = &sharedmem[threadIdx.x];
    sharedV = &sharedmem[0];
    sharedE[0] = ENERGY_F(0.0);
    sharedV += blockDim.x;
  }

  if(vflag || vflag_atom) {
    sharedV += threadIdx.x;
    sharedV[0 * blockDim.x] = ENERGY_F(0.0);
    sharedV[1 * blockDim.x] = ENERGY_F(0.0);
    sharedV[2 * blockDim.x] = ENERGY_F(0.0);
    sharedV[3 * blockDim.x] = ENERGY_F(0.0);
    sharedV[4 * blockDim.x] = ENERGY_F(0.0);
    sharedV[5 * blockDim.x] = ENERGY_F(0.0);
  }

  int ii = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  MYEMUDBG(if(ii == 0) printf("# CUDA: PairGranHookeCuda_Kernel: -- no binning --\n");)

    X_FLOAT xtmp, ytmp, ztmp;

  X_FLOAT4 myxtype;
  V_FLOAT4 myvradius, ovradius;
  F_FLOAT fxtmp, fytmp, fztmp, torquextmp, torqueytmp, torqueztmp;
  F_FLOAT delx, dely, delz;
  F_FLOAT radi, radj, radsum, r, rsqinv;
  F_FLOAT vr1, vr2, vr3, vnnr, vn1, vn2, vn3, vt1, vt2, vt3;
  F_FLOAT wr1, wr2, wr3;
  F_FLOAT vtr1, vtr2, vtr3, vrel;
  F_FLOAT meff, damp, ccel, tor1, tor2, tor3;
  F_FLOAT fn, fs, ft, fs1, fs2, fs3;

  int jnum = 0;
  int i, j;
  int* jlist;

  if(ii < _inum) {
    i = _ilist[ii];

    myxtype = fetchXType(i);
    myvradius = fetchVRadius(i);

    xtmp = myxtype.x;
    ytmp = myxtype.y;
    ztmp = myxtype.z;
    radi = myvradius.w;

    fxtmp = F_F(0.0);
    fytmp = F_F(0.0);
    fztmp = F_F(0.0);
    torquextmp = F_F(0.0);
    torqueytmp = F_F(0.0);
    torqueztmp = F_F(0.0);

    jnum = _numneigh[i];

    jlist = &_neighbors[i];
  }

  __syncthreads();

  for(int jj = 0; jj < jnum; jj++) {
    if(ii < _inum)
      if(jj < jnum) {
        j = jlist[jj * _nlocal];

        myxtype = fetchXType(j);
        ovradius = fetchVRadius(j);

        delx = xtmp - myxtype.x;
        dely = ytmp - myxtype.y;
        delz = ztmp - myxtype.z;

        radj = ovradius.w;
        radsum = radi + radj;

        const F_FLOAT rsq = delx * delx + dely * dely + delz * delz;

        if(rsq < radsum * radsum) {
          const F_FLOAT rinv = _RSQRT_(rsq);
          r = F_F(1.0) / rinv;
          rsqinv = F_F(1.0) / rsq;

          // relative translational velocity

          vr1 = myvradius.x - ovradius.x;
          vr2 = myvradius.y - ovradius.y;
          vr3 = myvradius.z - ovradius.z;

          // normal component

          vnnr = vr1 * delx + vr2 * dely + vr3 * delz;
          vn1 = delx * vnnr * rsqinv;
          vn2 = dely * vnnr * rsqinv;
          vn3 = delz * vnnr * rsqinv;

          // tangential component

          vt1 = vr1 - vn1;
          vt2 = vr2 - vn2;
          vt3 = vr3 - vn3;

          // relative rotational velocity
          V_FLOAT4 omegarmass_i = fetchOmegaRmass(i);
          V_FLOAT4 omegarmass_j = fetchOmegaRmass(j);

          wr1 = (radi * omegarmass_i.x + radj * omegarmass_j.x) * rinv;
          wr2 = (radi * omegarmass_i.y + radj * omegarmass_j.y) * rinv;
          wr3 = (radi * omegarmass_i.z + radj * omegarmass_j.z) * rinv;

          meff = omegarmass_i.w * omegarmass_j.w / (omegarmass_i.w + omegarmass_j.w);

          if(_mask[i] & _freeze_group_bit) meff = omegarmass_j.w;

          if(_mask[j] & _freeze_group_bit) meff = omegarmass_i.w;

          damp = meff * gamman * vnnr * rsqinv;
          ccel = kn * (radsum - r) * rinv - damp;

          vtr1 = vt1 - (delz * wr2 - dely * wr3);
          vtr2 = vt2 - (delx * wr3 - delz * wr1);
          vtr3 = vt3 - (dely * wr1 - delx * wr2);
          vrel = vtr1 * vtr1 + vtr2 * vtr2 + vtr3 * vtr3;
          vrel = _SQRT_(vrel);

          fn = xmu * fabs(ccel * r);
          fs = meff * gammat * vrel;
          ft = (vrel != F_F(0.0)) ? MIN(fn, fs) / vrel : F_F(0.0);

          fs1 = -ft * vtr1;
          fs2 = -ft * vtr2;
          fs3 = -ft * vtr3;

          F_FLOAT dxfp, dyfp, dzfp;
          fxtmp += dxfp = delx * ccel + fs1;
          fytmp += dyfp = dely * ccel + fs2;
          fztmp += dzfp = delz * ccel + fs3;

          tor1 = rinv * (dely * fs3 - delz * fs2);
          tor2 = rinv * (delz * fs1 - delx * fs3);
          tor3 = rinv * (delx * fs2 - dely * fs1);

          torquextmp -= radi * tor1;
          torqueytmp -= radi * tor2;
          torqueztmp -= radi * tor3;

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
    F_FLOAT* my_f = _f + i;
    *my_f += fxtmp;
    my_f += _nmax;
    *my_f += fytmp;
    my_f += _nmax;
    *my_f += fztmp;
    F_FLOAT* my_torque = _torque + i;
    *my_torque += torquextmp;
    my_torque += _nmax;
    *my_torque += torqueytmp;
    my_torque += _nmax;
    *my_torque += torqueztmp;
  }

  __syncthreads();

  if(eflag) sharedE[0] = evdwl;

  if(eflag_atom && i < _nlocal) _eatom[i] += evdwl;

  if(vflag_atom && i < _nlocal) {
    _vatom[i]         += ENERGY_F(0.5) * sharedV[0 * blockDim.x];
    _vatom[i + _nmax]   += ENERGY_F(0.5) * sharedV[1 * blockDim.x];
    _vatom[i + 2 * _nmax] += ENERGY_F(0.5) * sharedV[2 * blockDim.x];
    _vatom[i + 3 * _nmax] += ENERGY_F(0.5) * sharedV[3 * blockDim.x];
    _vatom[i + 4 * _nmax] += ENERGY_F(0.5) * sharedV[4 * blockDim.x];
    _vatom[i + 5 * _nmax] += ENERGY_F(0.5) * sharedV[5 * blockDim.x];
  }

  if(vflag || eflag) PairVirialCompute_A_Kernel(eflag, vflag, 0);
}
