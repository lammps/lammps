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
template <const int eflag, const int vflag>
static inline __device__ void PairVirialCompute_A_Kernel_Template()
{
  __syncthreads();
  ENERGY_CFLOAT* shared = sharedmem;

  if(eflag) {
    reduceBlock(shared);
    shared += blockDim.x;
  }

  if(vflag) {
    reduceBlock(shared + 0 * blockDim.x);
    reduceBlock(shared + 1 * blockDim.x);
    reduceBlock(shared + 2 * blockDim.x);
    reduceBlock(shared + 3 * blockDim.x);
    reduceBlock(shared + 4 * blockDim.x);
    reduceBlock(shared + 5 * blockDim.x);
  }

  if(threadIdx.x == 0) {
    shared = sharedmem;
    ENERGY_CFLOAT* buffer = (ENERGY_CFLOAT*) _buffer;

    if(eflag) {
      buffer[blockIdx.x * gridDim.y + blockIdx.y] = ENERGY_F(0.5) * shared[0];
      shared += blockDim.x;
      buffer += gridDim.x * gridDim.y;
    }

    if(vflag) {
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 0 * gridDim.x * gridDim.y] = ENERGY_F(0.5) * shared[0 * blockDim.x];
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 1 * gridDim.x * gridDim.y] = ENERGY_F(0.5) * shared[1 * blockDim.x];
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 2 * gridDim.x * gridDim.y] = ENERGY_F(0.5) * shared[2 * blockDim.x];
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 3 * gridDim.x * gridDim.y] = ENERGY_F(0.5) * shared[3 * blockDim.x];
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 4 * gridDim.x * gridDim.y] = ENERGY_F(0.5) * shared[4 * blockDim.x];
      buffer[blockIdx.x * gridDim.y + blockIdx.y + 5 * gridDim.x * gridDim.y] = ENERGY_F(0.5) * shared[5 * blockDim.x];
    }
  }

  __syncthreads();
}

__global__ void virial_fdotr_compute_kernel(int eflag)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  ENERGY_CFLOAT* sharedE = (ENERGY_CFLOAT*) &sharedmem[0];
  ENERGY_CFLOAT* sharedVirial = (ENERGY_CFLOAT*) &sharedE[blockDim.x];
  sharedE += threadIdx.x;
  sharedVirial += threadIdx.x;

  if(i < _nlocal) {

    F_CFLOAT x = _x[i];
    F_CFLOAT y = _x[i + _nmax];
    F_CFLOAT z = _x[i + 2 * _nmax];
    F_CFLOAT fx = _f[i];
    F_CFLOAT fy = _f[i + _nmax];
    F_CFLOAT fz = _f[i + 2 * _nmax];
    //if(fz*z*fz*z>1e-5) printf("V %i %i %e %e %e %e %e %e\n",i,_tag[i],x,y,z,fx,fy,fz);
    sharedVirial[0] = fx * x;
    sharedVirial[1 * blockDim.x] = fy * y;
    sharedVirial[2 * blockDim.x] = fz * z;
    sharedVirial[3 * blockDim.x] = fy * x;
    sharedVirial[4 * blockDim.x] = fz * x;
    sharedVirial[5 * blockDim.x] = fz * y;
  } else {
    sharedVirial[0] = 0;
    sharedVirial[1 * blockDim.x] = 0;
    sharedVirial[2 * blockDim.x] = 0;
    sharedVirial[3 * blockDim.x] = 0;
    sharedVirial[4 * blockDim.x] = 0;
    sharedVirial[5 * blockDim.x] = 0;
  }

  sharedVirial = (ENERGY_CFLOAT*) &sharedmem[0];
  sharedVirial += blockDim.x;
  reduceBlockP2(sharedVirial);
  reduceBlockP2(&sharedVirial[1 * blockDim.x]);
  reduceBlockP2(&sharedVirial[2 * blockDim.x]);
  reduceBlockP2(&sharedVirial[3 * blockDim.x]);
  reduceBlockP2(&sharedVirial[4 * blockDim.x]);
  reduceBlockP2(&sharedVirial[5 * blockDim.x]);

  if(threadIdx.x < 6) {
    ENERGY_CFLOAT* buffer = (ENERGY_CFLOAT*) _buffer;

    if(eflag) buffer = &buffer[gridDim.x * gridDim.y];

    buffer[blockIdx.x * gridDim.y + blockIdx.y + threadIdx.x * gridDim.x * gridDim.y] = sharedVirial[threadIdx.x * blockDim.x];
  }
}

/*#define vec3_scale(K,X,Y) Y.x = K*X.x;  Y.y = K*X.y;  Y.z = K*X.z;
#define vec3_scaleadd(K,X,Y,Z) Z.x = K*X.x+Y.x;  Z.y = K*X.y+Y.y;  Z.z = K*X.z+Y.z;
#define vec3_add(X,Y,Z) Z.x = X.x+Y.x;  Z.y = X.y+Y.y;  Z.z = X.z+Y.z;
#define vec3_dot(X,Y) (X.x*Y.x + X.y*Y.y + X.z*Y.z)*/

__device__ inline void vec3_scale(F_CFLOAT k, F_CFLOAT3 &x, F_CFLOAT3 &y)
{
  y.x = k * x.x;
  y.y = k * x.y;
  y.z = k * x.z;
}

__device__ inline void vec3_scale(F_CFLOAT k, F_CFLOAT4 &x, F_CFLOAT3 &y)
{
  y.x = k * x.x;
  y.y = k * x.y;
  y.z = k * x.z;
}

__device__ inline void vec3_scale(F_CFLOAT k, F_CFLOAT4 &x, F_CFLOAT4 &y)
{
  y.x = k * x.x;
  y.y = k * x.y;
  y.z = k * x.z;
}

__device__ inline void vec3_scaleadd(F_CFLOAT k, F_CFLOAT3 &x, F_CFLOAT3 &y, F_CFLOAT3 &z)
{
  z.x = k * x.x + y.x;
  z.y = k * x.y + y.y;
  z.z = k * x.z + y.z;
}

__device__ inline void vec3_add(F_CFLOAT3 &x, F_CFLOAT3 &y, F_CFLOAT3 &z)
{
  z.x = x.x + y.x;
  z.y = x.y + y.y;
  z.z = x.z + y.z;
}

__device__ inline F_CFLOAT vec3_dot(F_CFLOAT3 x, F_CFLOAT3 y)
{
  return x.x * y.x + x.y * y.y + x.z * y.z;
}

__device__ inline F_CFLOAT vec3_dot(F_CFLOAT4 x, F_CFLOAT4 y)
{
  return x.x * y.x + x.y * y.y + x.z * y.z;
}

/* ----------------------------------------------------------------------
   Fermi-like smoothing function
------------------------------------------------------------------------- */

__device__ inline F_CFLOAT F_fermi(F_CFLOAT &r, int &iparam)
{
  return F_F(1.0) / (F_F(1.0) + exp(-params[iparam].ZBLexpscale * (r - params[iparam].ZBLcut)));
}

/* ----------------------------------------------------------------------
   Fermi-like smoothing function derivative with respect to r
------------------------------------------------------------------------- */

__device__ inline F_CFLOAT F_fermi_d(F_CFLOAT &r, int &iparam)
{
  volatile const F_CFLOAT tmp =  exp(-params[iparam].ZBLexpscale * (r - params[iparam].ZBLcut));
  return params[iparam].ZBLexpscale * tmp /
         ((F_F(1.0) + tmp) * (F_F(1.0) + tmp));
}

__device__ inline F_CFLOAT ters_fc(F_CFLOAT r, F_CFLOAT ters_R, F_CFLOAT ters_D)
{
  return (r < ters_R - ters_D) ? F_F(1.0) : ((r > ters_R + ters_D) ?
         F_F(0.0) : F_F(0.5) * (F_F(1.0) - sin(PI2 * (r - ters_R) / ters_D)));
}

__device__ inline F_CFLOAT ters_fc_d(F_CFLOAT r, F_CFLOAT ters_R, F_CFLOAT ters_D)
{
  return ((r < ters_R - ters_D) || (r > ters_R + ters_D)) ?
         F_F(0.0) : -(PI4 / ters_D) * cos(PI2 * (r - ters_R) / ters_D);
}


__device__ inline F_CFLOAT ters_gijk(F_CFLOAT &cos_theta, int iparam)
{
  F_CFLOAT ters_c = params[iparam].c;
  F_CFLOAT ters_d = params[iparam].d;

  return params[iparam].gamma * (F_F(1.0) + pow(params[iparam].c / params[iparam].d, F_F(2.0)) -
                                 pow(ters_c, F_F(2.0)) / (pow(ters_d, F_F(2.0)) + pow(params[iparam].h - cos_theta, F_F(2.0))));
}

__device__ F_CFLOAT ters_gijk2(F_CFLOAT &cos_theta, int iparam)
{
  F_CFLOAT ters_c = params[iparam].c;
  F_CFLOAT ters_d = params[iparam].d;

  return params[iparam].gamma * (F_F(1.0) + pow(ters_c / ters_d, F_F(2.0)) -
                                 pow(ters_c, F_F(2.0)) / (pow(ters_d, F_F(2.0)) + pow(params[iparam].h - cos_theta, F_F(2.0))));
}

__device__ inline F_CFLOAT ters_gijk_d(F_CFLOAT costheta, int iparam)
{
  F_CFLOAT numerator = -F_F(2.0) * pow(params[iparam].c, F_F(2.0)) * (params[iparam].h - costheta);
  F_CFLOAT denominator = pow(pow(params[iparam].d, F_F(2.0)) +
                            pow(params[iparam].h - costheta, F_F(2.0)), F_F(2.0));
  return params[iparam].gamma * numerator / denominator;
}

__device__ inline F_CFLOAT zeta(int iparam, const F_CFLOAT rsqij, const F_CFLOAT rsqik,
                               F_CFLOAT3 &delij, F_CFLOAT3 &delik)
{
  F_CFLOAT rij, rik, costheta, arg, ex_delr;

  rij = sqrt(rsqij);
  rik = sqrt(rsqik);
  costheta = vec3_dot(delij, delik) / (rij * rik);

  arg = (params[iparam].powermint == 3) ? (params[iparam].lam3 * (rij - rik) * params[iparam].lam3 * (rij - rik) * params[iparam].lam3 * (rij - rik)) : params[iparam].lam3 * (rij - rik);

  if(arg > F_F(69.0776)) ex_delr = F_F(1.e30);
  else if(arg < -F_F(69.0776)) ex_delr = F_F(0.0);
  else ex_delr = exp(arg);

  return ters_fc(rik, params[iparam].bigr, params[iparam].bigd) * ex_delr * params[iparam].gamma * (F_F(1.0) + (params[iparam].c * params[iparam].c / (params[iparam].d * params[iparam].d)) -
         (params[iparam].c * params[iparam].c) / ((params[iparam].d * params[iparam].d) + (params[iparam].h - costheta) * (params[iparam].h - costheta)));
}

__device__ void repulsive(int iparam, F_CFLOAT rsq, F_CFLOAT &fforce,
                          int eflag, ENERGY_CFLOAT &eng)
{
  F_CFLOAT r, tmp_fc, tmp_fc_d, tmp_exp;

  F_CFLOAT ters_R = params[iparam].bigr;
  F_CFLOAT ters_D = params[iparam].bigd;
  r = sqrt(rsq);
  tmp_fc = ters_fc(r, ters_R, ters_D);
  tmp_fc_d = ters_fc_d(r, ters_R, ters_D);
  tmp_exp = exp(-params[iparam].lam1 * r);

  if(!_zbl) {
    fforce = -params[iparam].biga * tmp_exp * (tmp_fc_d - tmp_fc * params[iparam].lam1) / r;

    if(eflag) eng += tmp_fc * params[iparam].biga * tmp_exp;
  } else {
    F_CFLOAT const fforce_ters = params[iparam].biga * tmp_exp * (tmp_fc_d - tmp_fc * params[iparam].lam1);
    ENERGY_CFLOAT eng_ters = tmp_fc * params[iparam].biga * tmp_exp;

    F_CFLOAT r_ov_a = r / params[iparam].a_ij;
    F_CFLOAT phi = F_F(0.1818) * exp(-F_F(3.2) * r_ov_a) + F_F(0.5099) * exp(-F_F(0.9423) * r_ov_a) +
                  F_F(0.2802) * exp(-F_F(0.4029) * r_ov_a) + F_F(0.02817) * exp(-F_F(0.2016) * r_ov_a);
    F_CFLOAT dphi = (F_F(1.0) / params[iparam].a_ij) * (-F_F(3.2) * F_F(0.1818) * exp(-F_F(3.2) * r_ov_a) -
                   F_F(0.9423) * F_F(0.5099) * exp(-F_F(0.9423) * r_ov_a) -
                   F_F(0.4029) * F_F(0.2802) * exp(-F_F(0.4029) * r_ov_a) -
                   F_F(0.2016) * F_F(0.02817) * exp(-F_F(0.2016) * r_ov_a));
    F_CFLOAT fforce_ZBL = params[iparam].premult / (-r * r) * phi + params[iparam].premult / r * dphi;
    ENERGY_CFLOAT eng_ZBL = params[iparam].premult * (F_F(1.0) / r) * phi;

    fforce = -(-F_fermi_d(r, iparam) * (eng_ZBL - eng_ters) + fforce_ZBL + F_fermi(r, iparam) * (fforce_ters - fforce_ZBL)) / r;

    if(eflag)
      eng += eng_ZBL + F_fermi(r, iparam) * (eng_ters - eng_ZBL);
  }


}

/* ---------------------------------------------------------------------- */

__device__ inline F_CFLOAT ters_fa(F_CFLOAT r, int iparam, F_CFLOAT ters_R, F_CFLOAT ters_D)
{
  if(r > ters_R + ters_D) return F_F(0.0);

  if(_zbl)
    return -params[iparam].bigb * exp(-params[iparam].lam2 * r) * ters_fc(r, ters_R, ters_D) * F_fermi(r, iparam);
  else
    return -params[iparam].bigb * exp(-params[iparam].lam2 * r) * ters_fc(r, ters_R, ters_D);
}

/* ---------------------------------------------------------------------- */

__device__ inline F_CFLOAT ters_fa_d(F_CFLOAT r, int iparam, F_CFLOAT ters_R, F_CFLOAT ters_D)
{
  if(r > ters_R + ters_D) return F_F(0.0);

  if(_zbl)
    return params[iparam].bigb * exp(-params[iparam].lam2 * r) *
           ((params[iparam].lam2 * ters_fc(r, ters_R, ters_D) - ters_fc_d(r, ters_R, ters_D)) * F_fermi(r, iparam)
            - ters_fc(r, ters_R, ters_D) * F_fermi_d(r, iparam));
  else
    return params[iparam].bigb * exp(-params[iparam].lam2 * r) *
           (params[iparam].lam2 * ters_fc(r, ters_R, ters_D) - ters_fc_d(r, ters_R, ters_D));
}

/* ---------------------------------------------------------------------- */

__device__ inline F_CFLOAT ters_bij(F_CFLOAT zeta, int iparam)
{
  F_CFLOAT tmp = params[iparam].beta * zeta;

  if(tmp > params[iparam].c1) return F_F(1.0) / sqrt(tmp);

  if(tmp > params[iparam].c2)
    return (F_F(1.0) - pow(tmp, -params[iparam].powern) / (F_F(2.0) * params[iparam].powern)) / sqrt(tmp);

  if(tmp < params[iparam].c4) return F_F(1.0);

  if(tmp < params[iparam].c3)
    return F_F(1.0) - pow(tmp, params[iparam].powern) / (F_F(2.0) * params[iparam].powern);

  return pow(F_F(1.0) + pow(tmp, params[iparam].powern), -F_F(1.0) / (F_F(2.0) * params[iparam].powern));
}

/* ---------------------------------------------------------------------- */

__device__ inline F_CFLOAT ters_bij_d(F_CFLOAT zeta, int iparam)
{
  F_CFLOAT tmp = params[iparam].beta * zeta;

  if(tmp > params[iparam].c1) return params[iparam].beta * -F_F(0.5) * pow(tmp, -F_F(1.5));

  if(tmp > params[iparam].c2)
    return params[iparam].beta * (-F_F(0.5) * pow(tmp, -F_F(1.5)) *
                                  (F_F(1.0) - F_F(0.5) * (F_F(1.0) +  F_F(1.0) / (F_F(2.0) * params[iparam].powern)) *
                                   pow(tmp, -params[iparam].powern)));

  if(tmp < params[iparam].c4) return F_F(0.0);

  if(tmp < params[iparam].c3)
    return -F_F(0.5) * params[iparam].beta * pow(tmp, params[iparam].powern - F_F(1.0));

  F_CFLOAT tmp_n = pow(tmp, params[iparam].powern);
  return -F_F(0.5) * pow(F_F(1.0) + tmp_n, -F_F(1.0) - (F_F(1.0) / (F_F(2.0) * params[iparam].powern))) * tmp_n / zeta;
}

__device__ void force_zeta(int iparam, F_CFLOAT rsq, F_CFLOAT zeta_ij,
                           F_CFLOAT &fforce, F_CFLOAT &prefactor,
                           int eflag, F_CFLOAT &eng)
{
  F_CFLOAT r, fa, fa_d, bij;
  F_CFLOAT ters_R = params[iparam].bigr;
  F_CFLOAT ters_D = params[iparam].bigd;
  r = sqrt(rsq);
  fa = ters_fa(r, iparam, ters_R, ters_D);
  fa_d = ters_fa_d(r, iparam, ters_R, ters_D);
  bij = ters_bij(zeta_ij, iparam);
  fforce = F_F(0.5) * bij * fa_d / r;
  prefactor = -F_F(0.5) * fa * ters_bij_d(zeta_ij, iparam);

  if(eflag) eng += bij * fa;
}

__device__ void force_zeta_prefactor_force(int iparam, F_CFLOAT rsq, F_CFLOAT zeta_ij,
    F_CFLOAT &fforce, F_CFLOAT &prefactor)
{
  F_CFLOAT r, fa, fa_d, bij;
  F_CFLOAT ters_R = params[iparam].bigr;
  F_CFLOAT ters_D = params[iparam].bigd;
  r = sqrt(rsq);
  fa = ters_fa(r, iparam, ters_R, ters_D);
  fa_d = ters_fa_d(r, iparam, ters_R, ters_D);
  bij = ters_bij(zeta_ij, iparam);
  fforce = F_F(0.5) * bij * fa_d / r;
  prefactor = -F_F(0.5) * fa * ters_bij_d(zeta_ij, iparam);
}

__device__ void force_zeta_prefactor(int iparam, F_CFLOAT rsq, F_CFLOAT zeta_ij,
                                     F_CFLOAT &prefactor)
{
  F_CFLOAT r, fa;
  r = sqrt(rsq);
  fa = ters_fa(r, iparam, params[iparam].bigr, params[iparam].bigd);
  prefactor = -F_F(0.5) * fa * ters_bij_d(zeta_ij, iparam);
}


__device__ void costheta_d(F_CFLOAT3 &rij_hat, F_CFLOAT &rij,
                           F_CFLOAT3 &rik_hat, F_CFLOAT &rik,
                           F_CFLOAT3 &dri, F_CFLOAT3 &drj, F_CFLOAT3 &drk)
{
  // first element is derivative wrt Ri, second wrt Rj, third wrt Rk

  F_CFLOAT cos_theta = vec3_dot(rij_hat, rik_hat);

  vec3_scaleadd(-cos_theta, rij_hat, rik_hat, drj);
  vec3_scale(F_F(1.0) / rij, drj, drj);
  vec3_scaleadd(-cos_theta, rik_hat, rij_hat, drk);
  vec3_scale(F_F(1.0) / rik, drk, drk);
  vec3_add(drj, drk, dri);
  vec3_scale(-F_F(1.0), dri, dri);
}

__device__ void ters_zetaterm_d(F_CFLOAT prefactor,
                                F_CFLOAT3 &rij_hat, F_CFLOAT rij,
                                F_CFLOAT3 &rik_hat, F_CFLOAT rik,
                                F_CFLOAT3 &dri, F_CFLOAT3 &drj, F_CFLOAT3 &drk,
                                int iparam)
{
  F_CFLOAT ex_delr, ex_delr_d, tmp;
  F_CFLOAT3 dcosdri, dcosdrj, dcosdrk;

  if(params[iparam].powermint == 3) tmp = (params[iparam].lam3 * (rij - rik) * params[iparam].lam3 * (rij - rik) * params[iparam].lam3 * (rij - rik));
  else tmp = params[iparam].lam3 * (rij - rik);

  if(tmp > F_F(69.0776)) ex_delr = F_F(1.e30);
  else if(tmp < -F_F(69.0776)) ex_delr = F_F(0.0);
  else ex_delr = exp(tmp);

  if(params[iparam].powermint == 3)
    ex_delr_d = F_F(3.0) * (params[iparam].lam3 * params[iparam].lam3 * params[iparam].lam3) * (rij - rik) * (rij - rik) * ex_delr;
  else ex_delr_d = params[iparam].lam3 * ex_delr;


  const F_CFLOAT cos_theta = vec3_dot(rij_hat, rik_hat);
  costheta_d(rij_hat, rij, rik_hat, rik, dcosdri, dcosdrj, dcosdrk);

  const F_CFLOAT gijk = params[iparam].gamma * (F_F(1.0) + (params[iparam].c * params[iparam].c) / (params[iparam].d * params[iparam].d) -
                       (params[iparam].c * params[iparam].c) / (params[iparam].d * params[iparam].d + (params[iparam].h - cos_theta) * (params[iparam].h - cos_theta)));
  const F_CFLOAT numerator = -F_F(2.0) * params[iparam].c * params[iparam].c * (params[iparam].h - cos_theta);
  const F_CFLOAT denominator = (params[iparam].d * params[iparam].d) +
                              (params[iparam].h - cos_theta) * (params[iparam].h - cos_theta);
  const F_CFLOAT gijk_d = params[iparam].gamma * numerator / (denominator * denominator); // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);
  const F_CFLOAT fc = ters_fc(rik, params[iparam].bigr, params[iparam].bigd);
  const F_CFLOAT dfc = ters_fc_d(rik, params[iparam].bigr, params[iparam].bigd);


  vec3_scale(-dfc * gijk * ex_delr, rik_hat, dri);
  vec3_scaleadd(fc * gijk_d * ex_delr, dcosdri, dri, dri);
  vec3_scaleadd(fc * gijk * ex_delr_d, rik_hat, dri, dri);
  vec3_scaleadd(-fc * gijk * ex_delr_d, rij_hat, dri, dri);
  vec3_scale(prefactor, dri, dri);
  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc * gijk_d * ex_delr, dcosdrj, drj);
  vec3_scaleadd(fc * gijk * ex_delr_d, rij_hat, drj, drj);
  vec3_scale(prefactor, drj, drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc * gijk * ex_delr, rik_hat, drk);
  vec3_scaleadd(fc * gijk_d * ex_delr, dcosdrk, drk, drk);
  vec3_scaleadd(-fc * gijk * ex_delr_d, rik_hat, drk, drk);
  vec3_scale(prefactor, drk, drk);
}

__device__ void ters_zetaterm_d_fi(F_CFLOAT &prefactor,
                                   F_CFLOAT3 &rij_hat, F_CFLOAT &rij,
                                   F_CFLOAT3 &rik_hat, F_CFLOAT &rik,
                                   F_CFLOAT3 &dri,  int &iparam)
{
  F_CFLOAT ex_delr, ex_delr_d, tmp;

  if(params[iparam].powermint == 3) tmp = (params[iparam].lam3 * (rij - rik) * params[iparam].lam3 * (rij - rik) * params[iparam].lam3 * (rij - rik));
  else tmp = params[iparam].lam3 * (rij - rik);

  if(tmp > F_F(69.0776)) ex_delr = F_F(1.e30);
  else if(tmp < -F_F(69.0776)) ex_delr = F_F(0.0);
  else ex_delr = exp(tmp);

  if(params[iparam].powermint == 3)
    ex_delr_d = F_F(3.0) * (params[iparam].lam3 * params[iparam].lam3 * params[iparam].lam3) * (rij - rik) * (rij - rik) * ex_delr;
  else ex_delr_d = params[iparam].lam3 * ex_delr;

  const F_CFLOAT cos_theta = vec3_dot(rij_hat, rik_hat);
  //costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);


  F_CFLOAT3 dcosdri;
  vec3_scaleadd(-cos_theta, rij_hat, rik_hat, dri);
  vec3_scale(F_F(1.0) / rij, dri, dri);
  vec3_scaleadd(-cos_theta, rik_hat, rij_hat, dcosdri);
  vec3_scale(F_F(1.0) / rik, dcosdri, dcosdri);
  vec3_add(dri, dcosdri, dcosdri);
  vec3_scale(-F_F(1.0), dcosdri, dcosdri);

  const F_CFLOAT gijk = params[iparam].gamma * (F_F(1.0) + (params[iparam].c * params[iparam].c) / (params[iparam].d * params[iparam].d) -
                       (params[iparam].c * params[iparam].c) / (params[iparam].d * params[iparam].d + (params[iparam].h - cos_theta) * (params[iparam].h - cos_theta)));
  const F_CFLOAT numerator = -F_F(2.0) * params[iparam].c * params[iparam].c * (params[iparam].h - cos_theta);
  const F_CFLOAT denominator = (params[iparam].d * params[iparam].d) +
                              (params[iparam].h - cos_theta) * (params[iparam].h - cos_theta);
  const F_CFLOAT gijk_d = params[iparam].gamma * numerator / (denominator * denominator); // compute the derivative wrt Ri
  //
  const F_CFLOAT fc = ters_fc(rik, params[iparam].bigr, params[iparam].bigd);
  const F_CFLOAT dfc = ters_fc_d(rik, params[iparam].bigr, params[iparam].bigd);

  vec3_scale(-dfc * gijk * ex_delr, rik_hat, dri);
  vec3_scaleadd(fc * gijk_d * ex_delr, dcosdri, dri, dri);
  vec3_scaleadd(fc * gijk * ex_delr_d, rik_hat, dri, dri);
  vec3_scaleadd(-fc * gijk * ex_delr_d, rij_hat, dri, dri);
  vec3_scale(prefactor, dri, dri);

}

__device__ void ters_zetaterm_d_fj(F_CFLOAT &prefactor,
                                   F_CFLOAT3 &rij_hat, F_CFLOAT &rij,
                                   F_CFLOAT3 &rik_hat, F_CFLOAT &rik,
                                   F_CFLOAT3 &drj, int &iparam)
{
  F_CFLOAT ex_delr, ex_delr_d, tmp;

  if(params[iparam].powermint == 3) tmp = (params[iparam].lam3 * (rij - rik) * params[iparam].lam3 * (rij - rik) * params[iparam].lam3 * (rij - rik));
  else tmp = params[iparam].lam3 * (rij - rik);

  if(tmp > F_F(69.0776)) ex_delr = F_F(1.e30);
  else if(tmp < -F_F(69.0776)) ex_delr = F_F(0.0);
  else ex_delr = exp(tmp);

  if(params[iparam].powermint == 3)
    ex_delr_d = F_F(3.0) * (params[iparam].lam3 * params[iparam].lam3 * params[iparam].lam3) * (rij - rik) * (rij - rik) * ex_delr;
  else ex_delr_d = params[iparam].lam3 * ex_delr;

  const F_CFLOAT cos_theta = vec3_dot(rij_hat, rik_hat);
  vec3_scaleadd(-cos_theta, rij_hat, rik_hat, drj);
  vec3_scale(F_F(1.0) / rij, drj, drj);

  const F_CFLOAT gijk = params[iparam].gamma * (F_F(1.0) + (params[iparam].c * params[iparam].c) / (params[iparam].d * params[iparam].d) -
                       (params[iparam].c * params[iparam].c) / (params[iparam].d * params[iparam].d + (params[iparam].h - cos_theta) * (params[iparam].h - cos_theta)));
  const F_CFLOAT numerator = -F_F(2.0) * params[iparam].c * params[iparam].c * (params[iparam].h - cos_theta);
  const F_CFLOAT denominator = (params[iparam].d * params[iparam].d) +
                              (params[iparam].h - cos_theta) * (params[iparam].h - cos_theta);
  const F_CFLOAT gijk_d = params[iparam].gamma * numerator / (denominator * denominator); // compute the derivative wrt Ri

  const F_CFLOAT fc = ters_fc(rik, params[iparam].bigr, params[iparam].bigd);

  vec3_scale(fc * gijk_d * ex_delr, drj, drj);
  vec3_scaleadd(fc * gijk * ex_delr_d, rij_hat, drj, drj);
  vec3_scale(prefactor, drj, drj);
}

__device__ void ters_zetaterm_d_fk(F_CFLOAT &prefactor,
                                   F_CFLOAT3 &rij_hat, F_CFLOAT &rij,
                                   F_CFLOAT3 &rik_hat, F_CFLOAT &rik,
                                   F_CFLOAT3 &drk, int &iparam)
{
  F_CFLOAT ex_delr, ex_delr_d, tmp;

  if(params[iparam].powermint == 3) tmp = (params[iparam].lam3 * (rij - rik) * params[iparam].lam3 * (rij - rik) * params[iparam].lam3 * (rij - rik));
  else tmp = params[iparam].lam3 * (rij - rik);

  if(tmp > F_F(69.0776)) ex_delr = F_F(1.e30);
  else if(tmp < -F_F(69.0776)) ex_delr = F_F(0.0);
  else ex_delr = exp(tmp);

  if(params[iparam].powermint == 3)
    ex_delr_d = F_F(3.0) * (params[iparam].lam3 * params[iparam].lam3 * params[iparam].lam3) * (rij - rik) * (rij - rik) * ex_delr;
  else ex_delr_d = params[iparam].lam3 * ex_delr;

  const F_CFLOAT cos_theta = vec3_dot(rij_hat, rik_hat);
  vec3_scaleadd(-cos_theta, rik_hat, rij_hat, drk);
  vec3_scale(F_F(1.0) / rik, drk, drk);

  const F_CFLOAT gijk = params[iparam].gamma * (F_F(1.0) + (params[iparam].c * params[iparam].c) / (params[iparam].d * params[iparam].d) -
                       (params[iparam].c * params[iparam].c) / (params[iparam].d * params[iparam].d + (params[iparam].h - cos_theta) * (params[iparam].h - cos_theta)));
  const F_CFLOAT numerator = -F_F(2.0) * params[iparam].c * params[iparam].c * (params[iparam].h - cos_theta);
  const F_CFLOAT denominator = (params[iparam].d * params[iparam].d) +
                              (params[iparam].h - cos_theta) * (params[iparam].h - cos_theta);
  const F_CFLOAT gijk_d = params[iparam].gamma * numerator / (denominator * denominator); // compute the derivative wrt Ri

  const F_CFLOAT fc = ters_fc(rik, params[iparam].bigr, params[iparam].bigd);
  const F_CFLOAT dfc = ters_fc_d(rik, params[iparam].bigr, params[iparam].bigd);

  vec3_scale(fc * gijk_d * ex_delr, drk, drk);
  vec3_scaleadd(dfc * gijk * ex_delr, rik_hat, drk, drk);
  vec3_scaleadd(-fc * gijk * ex_delr_d, rik_hat, drk, drk);
  vec3_scale(prefactor, drk, drk);
}

__device__ void attractive(int iparam, F_CFLOAT prefactor,
                           F_CFLOAT4 &delij,
                           F_CFLOAT4 &delik,
                           F_CFLOAT3 &fi, F_CFLOAT3 &fj, F_CFLOAT3 &fk)
{
  F_CFLOAT3 rij_hat, rik_hat;
  F_CFLOAT rij, rijinv, rik, rikinv;

  rij = sqrt(delij.w);
  rijinv = F_F(1.0) / rij;
  vec3_scale(rijinv, delij, rij_hat);

  rik = sqrt(delik.w);
  rikinv = F_F(1.0) / rik;
  vec3_scale(rikinv, delik, rik_hat);

  ters_zetaterm_d(prefactor, rij_hat, rij, rik_hat, rik, fi, fj, fk, iparam);
}

__device__ void attractive_fi(int &iparam, F_CFLOAT &prefactor,
                              F_CFLOAT4 &delij,
                              F_CFLOAT4 &delik,
                              F_CFLOAT3 &f)
{
  F_CFLOAT3 rij_hat, rik_hat;
  F_CFLOAT rij, rijinv, rik, rikinv;

  rij = sqrt(delij.w);
  rijinv = F_F(1.0) / rij;
  vec3_scale(rijinv, delij, rij_hat);

  rik = sqrt(delik.w);
  rikinv = F_F(1.0) / rik;
  vec3_scale(rikinv, delik, rik_hat);

  ters_zetaterm_d_fi(prefactor, rij_hat, rij, rik_hat, rik, f, iparam);
}

__device__ void attractive_fj(int iparam, F_CFLOAT prefactor,
                              F_CFLOAT4 &delij,
                              F_CFLOAT4 &delik,
                              F_CFLOAT3 &f)
{
  F_CFLOAT3 rij_hat, rik_hat;
  F_CFLOAT rij, rijinv, rik, rikinv;

  rij = sqrt(delij.w);
  rijinv = F_F(1.0) / rij;
  vec3_scale(rijinv, delij, rij_hat);

  rik = sqrt(delik.w);
  rikinv = F_F(1.0) / rik;
  vec3_scale(rikinv, delik, rik_hat);

  ters_zetaterm_d_fj(prefactor, rij_hat, rij, rik_hat, rik, f, iparam);
}

__device__ void attractive_fk(int iparam, F_CFLOAT prefactor,
                              F_CFLOAT4 &delij,
                              F_CFLOAT4 &delik,
                              F_CFLOAT3 &f)
{
  F_CFLOAT3 rij_hat, rik_hat;
  F_CFLOAT rij, rijinv, rik, rikinv;

  rij = sqrt(delij.w);
  rijinv = F_F(1.0) / rij;
  vec3_scale(rijinv, delij, rij_hat);

  rik = sqrt(delik.w);
  rikinv = F_F(1.0) / rik;
  vec3_scale(rikinv, delik, rik_hat);

  ters_zetaterm_d_fk(prefactor, rij_hat, rij, rik_hat, rik, f, iparam);
}

__global__ void Pair_Tersoff_Kernel_TpA_RIJ()//F_CFLOAT4* _glob_r_ij,int* _glob_numneigh_red,int* _glob_neighbors_red,int* _glob_neightype_red)
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

      if(delij.w < params[iparam_ij].cutsq) {
        _glob_neighbors_red[i + neigh_red * _nall] = j;
        _glob_neightype_red[i + neigh_red * _nall] = jtype;
        _glob_r_ij[i + neigh_red * _nall] = delij;
        neigh_red++;
      }
    }
  }

  _glob_numneigh_red[i] = neigh_red;
}


__global__ void Pair_Tersoff_Kernel_TpA_ZetaIJ()//F_CFLOAT* _glob_zeta_ij,F_CFLOAT4* _glob_r_ij,int* _glob_numneigh_red,int* _glob_neighbors_red,int* _glob_neightype_red)
{

  int ii = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(ii >= _nall) return;


  F_CFLOAT4 delij;
  F_CFLOAT4 delik;

  int itype, jnum, i, j;
  int* jlist;
  i = ii;
  itype = map[(static_cast <int>(_type[i]))];

  jnum = _glob_numneigh_red[i];
  jlist = &_glob_neighbors_red[i];

  __syncthreads();

  for(int jj = 0; jj < jnum; jj++) {
    if(jj < jnum) {

      j = jlist[jj * _nall];
      j &= NEIGHMASK;
      int jtype = _glob_neightype_red[i + jj * _nall];
      delij = _glob_r_ij[i + jj * _nall];

      int iparam_ij = elem2param[(itype * nelements + jtype) * nelements + jtype];

      if(delij.w < params[iparam_ij].cutsq) {
        F_CFLOAT zeta_ij = 0.0;
        F_CFLOAT3 delij3 = {delij.x, delij.y, delij.z};

        for(int kk = 0; kk < jnum; kk++) {
          if(jj == kk) continue;

          int k = jlist[kk * _nall];
          k &= NEIGHMASK;

          int ktype = _glob_neightype_red[i + kk * _nall];
          delik = _glob_r_ij[i + kk * _nall];
          F_CFLOAT3 delik3 = {delik.x, delik.y, delik.z};
          int iparam_ijk = elem2param[(itype * nelements + jtype) * nelements + ktype];
          const F_CFLOAT rsqki = delik.w;

          if(rsqki <= params[iparam_ijk].cutsq)
            zeta_ij += zeta(iparam_ijk, delij.w, rsqki, delij3, delik3);
        }

        _glob_zeta_ij[i + jj * _nall] = zeta_ij;
      }
    }
  }
}

//back3: num 12 steps 10: ZetaIJ/TPA 0.255/0.106
//back5: num 12 steps 10: ZetaIJ/TPA 0.257/0.098
//back6: num 12 steps 10: ZetaIJ/TPA 0.027/0.097 /rij berechnung extra
//back12: num 12 steps 10: ZetaIJ/TPA 0.026/0.070
//back15: num 12 steps 10: ZetaIJ/TPA 0.0137/0.0287 //pow beseitigt
//        num 12 steps 10: ZetaIJ/TPA 0.0137/0.027
template <int eflag, int vflagm>
__global__ void Pair_Tersoff_Kernel_TpA(int eflag_atom, int vflag_atom) //,F_CFLOAT* _glob_zeta_ij,F_CFLOAT4* _glob_r_ij,int* _glob_numneigh_red,int* _glob_neighbors_red,int* _glob_neightype_red)
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
  F_CFLOAT prefactor_ij, prefactor_ji;

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

      if(delij.w < params[iparam_ij].cutsq) {
        F_CFLOAT dxfp, dyfp, dzfp;
        repulsive(iparam_ij, delij.w, fpair, eflag, evdwl);
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



        force_zeta(iparam_ij, delij.w, _glob_zeta_ij[i + jj * _nall], fpair, prefactor_ij, eflag, evdwl);
        fxtmp -=
          dxfp = delij.x * fpair;
        fytmp -=
          dyfp = delij.y * fpair;
        fztmp -=
          dzfp = delij.z * fpair;

        if(vflagm) {
          sharedV[0 * blockDim.x] -= ENERGY_F(2.0) * delij.x * dxfp;
          sharedV[1 * blockDim.x] -= ENERGY_F(2.0) * delij.y * dyfp;
          sharedV[2 * blockDim.x] -= ENERGY_F(2.0) * delij.z * dzfp;
          sharedV[3 * blockDim.x] -= ENERGY_F(2.0) * delij.x * dyfp;
          sharedV[4 * blockDim.x] -= ENERGY_F(2.0) * delij.x * dzfp;
          sharedV[5 * blockDim.x] -= ENERGY_F(2.0) * delij.y * dzfp;
        }

        int j_jj = 0;

        //#pragma unroll 1
        for(int kk = 0; kk < _glob_numneigh_red[j]; kk++) {
          if(_glob_neighbors_red[j + kk * _nall] == i) j_jj = kk;
        }

        force_zeta_prefactor_force(iparam_ji, delij.w, _glob_zeta_ij[j + j_jj * _nall], fpair, prefactor_ji);

        fxtmp -=
          dxfp = delij.x * fpair;
        fytmp -=
          dyfp = delij.y * fpair;
        fztmp -=
          dzfp = delij.z * fpair;



        vec3_scale(F_F(-1.0), delij, delij);

#pragma unroll 1

        for(int kk = 0; kk < jnum_red; kk++) {
          if(jj == kk) continue;

          int k = jlist_red[kk * _nall];
          k &= NEIGHMASK;

          if(vflagm)
            myxtype_k = fetchXType(k);

          delik = _glob_r_ij[i + kk * _nall];

          int ktype = _glob_neightype_red[i + kk * _nall];
          int iparam_ijk = elem2param[(itype * nelements + jtype) * nelements + ktype];
          vec3_scale(F_F(-1.0), delik, delik);

          if(delik.w <= params[iparam_ijk].cutsq) {
            if(vflagm) {
              F_CFLOAT3 fi, fj, fk;
              attractive(iparam_ijk, prefactor_ij,
                         delij, delik, fi, fj, fk);
              fxtmp += fi.x;
              fytmp += fi.y;
              fztmp += fi.z;

              sharedV[0 * blockDim.x] += ENERGY_F(2.0) * myxtype_i.x * fi.x;
              sharedV[1 * blockDim.x] += ENERGY_F(2.0) * myxtype_i.y * fi.y;
              sharedV[2 * blockDim.x] += ENERGY_F(2.0) * myxtype_i.z * fi.z;
              sharedV[3 * blockDim.x] += ENERGY_F(2.0) * myxtype_i.x * fi.y;
              sharedV[4 * blockDim.x] += ENERGY_F(2.0) * myxtype_i.x * fi.z;
              sharedV[5 * blockDim.x] += ENERGY_F(2.0) * myxtype_i.y * fi.z;

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
            } else {
              F_CFLOAT3 fi; //local variable
              attractive_fi(iparam_ijk, prefactor_ij,
                            delij, delik, fi);
              fxtmp += fi.x;
              fytmp += fi.y;
              fztmp += fi.z;

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

          int iparam_jik = elem2param[(jtype * nelements + itype) * nelements + ktype];
          int iparam_jki = elem2param[(jtype * nelements + ktype) * nelements + itype];


          vec3_scale(F_F(-1.0), delij, delij);

          if(deljk.w <= params[iparam_jik].cutsq) {
            F_CFLOAT3 ftmp; //local variable

            attractive_fj(iparam_jik, prefactor_ji,
                          delij, deljk, ftmp);
            fxtmp += ftmp.x;
            fytmp += ftmp.y;
            fztmp += ftmp.z;
            int iparam_jk = elem2param[(jtype * nelements + ktype) * nelements + ktype];
            F_CFLOAT prefactor_jk;
            force_zeta_prefactor(iparam_jk, deljk.w, _glob_zeta_ij[j + kk * _nall], prefactor_jk);

            attractive_fk(iparam_jki, prefactor_jk,
                          deljk, delij, ftmp);
            fxtmp += ftmp.x;
            fytmp += ftmp.y;
            fztmp += ftmp.z;

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
