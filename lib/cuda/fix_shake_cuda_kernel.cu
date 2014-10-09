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

__device__ void v_tally(int &vflag_global, int &vflag_atom, int &n, int* list, ENERGY_CFLOAT total, ENERGY_CFLOAT* v)
{
  /*if(vflag_global)
  {
    ENERGY_CFLOAT fraction = n/total;
  ENERGY_CFLOAT* shared = &sharedmem[threadIdx.x];
    *shared += fraction*v[0]; shared+=blockDim.x;
    *shared += fraction*v[1]; shared+=blockDim.x;
    *shared += fraction*v[2]; shared+=blockDim.x;
    *shared += fraction*v[3]; shared+=blockDim.x;
    *shared += fraction*v[4]; shared+=blockDim.x;
    *shared += fraction*v[5];
  }*/
  if(vflag_atom) {
    ENERGY_CFLOAT fraction = ENERGY_F(1.0) / total;

    for(int i = 0; i < n; i++) {
      int m = list[i];
      ENERGY_CFLOAT* myvatom = &_vatom[m];

      *myvatom += fraction * v[0];
      myvatom += _nmax;
      *myvatom += fraction * v[1];
      myvatom += _nmax;
      *myvatom += fraction * v[2];
      myvatom += _nmax;
      *myvatom += fraction * v[3];
      myvatom += _nmax;
      *myvatom += fraction * v[4];
      myvatom += _nmax;
      *myvatom += fraction * v[5];
    }
  }
}

inline __device__ void minimum_image(X_CFLOAT3 &delta)
{
  if(_triclinic == 0) {
    if(_periodicity[0]) {
      delta.x += delta.x < -X_F(0.5) * _prd[0] ? _prd[0] :
                 (delta.x >  X_F(0.5) * _prd[0] ? -_prd[0] : X_F(0.0));
    }

    if(_periodicity[1]) {
      delta.y += delta.y < -X_F(0.5) * _prd[1] ? _prd[1] :
                 (delta.y >  X_F(0.5) * _prd[1] ? -_prd[1] : X_F(0.0));
    }

    if(_periodicity[2]) {
      delta.z += delta.z < -X_F(0.5) * _prd[2] ? _prd[2] :
                 (delta.z >  X_F(0.5) * _prd[2] ? -_prd[2] : X_F(0.0));
    }

  } else {
    if(_periodicity[1]) {
      delta.z += delta.z < -X_F(0.5) * _prd[2] ? _prd[2] :
                 (delta.z >  X_F(0.5) * _prd[2] ? -_prd[2] : X_F(0.0));
      delta.y += delta.z < -X_F(0.5) * _prd[2] ? _h[3] :
                 (delta.z >  X_F(0.5) * _prd[2] ? -_h[3] : X_F(0.0));
      delta.x += delta.z < -X_F(0.5) * _prd[2] ? _h[4] :
                 (delta.z >  X_F(0.5) * _prd[2] ? -_h[4] : X_F(0.0));

    }

    if(_periodicity[1]) {
      delta.y += delta.y < -X_F(0.5) * _prd[1] ? _prd[1] :
                 (delta.y >  X_F(0.5) * _prd[1] ? -_prd[1] : X_F(0.0));
      delta.x += delta.y < -X_F(0.5) * _prd[1] ? _h[5] :
                 (delta.y >  X_F(0.5) * _prd[1] ? -_h[5] : X_F(0.0));

    }

    if(_periodicity[0]) {
      delta.x += delta.x < -X_F(0.5) * _prd[0] ? _prd[0] :
                 (delta.x >  X_F(0.5) * _prd[0] ? -_prd[0] : X_F(0.0));
    }
  }
}

__global__ void FixShakeCuda_UnconstrainedUpdate_Kernel()
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i >= _nlocal) return;

  X_CFLOAT3 my_xshake = {X_F(0.0), X_F(0.0), X_F(0.0)};

  if(_shake_flag[i]) {
    F_CFLOAT* my_f = _f + i;
    V_CFLOAT* my_v = _v + i;
    X_CFLOAT* my_x = _x + i;

    V_CFLOAT 		dtfmsq = _dtfsq;

    if(_rmass_flag) dtfmsq *= V_F(1.0) / _rmass[i];
    else 			dtfmsq *= V_F(1.0) / _mass[_type[i]];

    my_xshake.x =  *my_x + _dtv* *my_v + dtfmsq* *my_f;
    my_f += _nmax;
    my_v += _nmax;
    my_x += _nmax;
    my_xshake.y =  *my_x + _dtv* *my_v + dtfmsq* *my_f;
    my_f += _nmax;
    my_v += _nmax;
    my_x += _nmax;
    my_xshake.z =  *my_x + _dtv* *my_v + dtfmsq* *my_f;
  }

  _xshake[i] = my_xshake;
}




__device__ void FixShakeCuda_Shake2(int &vflag, int &vflag_atom, int &m)
{
  int nlist, list[2];
  ENERGY_CFLOAT v[6];
  X_CFLOAT invmass0, invmass1;

  // local atom IDs and constraint distances

  int i0 = _map_array[_shake_atom[m]];
  int i1 = _map_array[_shake_atom[m + _nmax]];
  X_CFLOAT bond1 = _bond_distance[_shake_type[m]];

  // r01 = distance vec between atoms, with PBC

  X_CFLOAT3 r01;

  X_CFLOAT4 x_i0, x_i1;
  x_i0 = fetchXType(i0);
  x_i1 = fetchXType(i1);

  r01.x = x_i0.x - x_i1.x;
  r01.y = x_i0.y - x_i1.y;
  r01.z = x_i0.z - x_i1.z;
  minimum_image(r01);

  // s01 = distance vec after unconstrained update, with PBC

  X_CFLOAT3 s01;
  X_CFLOAT3 xs_i0 = _xshake[i0];
  X_CFLOAT3 xs_i1 = _xshake[i1];

  s01.x = xs_i0.x - xs_i1.x;
  s01.y = xs_i0.y - xs_i1.y;
  s01.z = xs_i0.z - xs_i1.z;
  minimum_image(s01);

  // scalar distances between atoms

  X_CFLOAT r01sq = r01.x * r01.x + r01.y * r01.y + r01.z * r01.z;
  X_CFLOAT s01sq = s01.x * s01.x + s01.y * s01.y + s01.z * s01.z;

  // a,b,c = coeffs in quadratic equation for lamda

  if(_rmass_flag) {
    invmass0 = X_F(1.0) / _rmass[i0];
    invmass1 = X_F(1.0) / _rmass[i1];
  } else {
    invmass0 = X_F(1.0) / _mass[static_cast <int>(x_i0.w)];
    invmass1 = X_F(1.0) / _mass[static_cast <int>(x_i1.w)];
  }

  X_CFLOAT a = (invmass0 + invmass1) * (invmass0 + invmass1) * r01sq;
  X_CFLOAT b = X_F(2.0) * (invmass0 + invmass1) *
              (s01.x * r01.x + s01.y * r01.y + s01.z * r01.z);
  X_CFLOAT c = s01sq - bond1 * bond1;

  // error check

  X_CFLOAT determ = b * b - X_F(4.0) * a * c;

  if(determ < X_F(0.0)) {
    _flag[0]++;
    determ = X_F(0.0);
  }

  // exact quadratic solution for lamda

  X_CFLOAT lamda, lamda1, lamda2;
  lamda1 = -b + _SQRT_(determ);
  lamda2 = -lamda1 - X_F(2.0) * b;
  lamda1 *= X_F(1.0) / (X_F(2.0) * a);
  lamda2 *= X_F(1.0) / (X_F(2.0) * a);

  lamda = (fabs(lamda1) <= fabs(lamda2)) ? lamda1 : lamda2;

  // update forces if atom is owned by this processor

  lamda *= X_F(1.0) / _dtfsq;


  //attenion: are shake clusters <-> atom unique?
  nlist = 0;

  if(i0 < _nlocal) {
    _f[i0]         += lamda * r01.x;
    _f[i0 + _nmax]   += lamda * r01.y;
    _f[i0 + 2 * _nmax] += lamda * r01.z;
    list[nlist++] = i0;
  }

  if(i1 < _nlocal) {
    _f[i1]         -= lamda * r01.x;
    _f[i1 + _nmax]   -= lamda * r01.y;
    _f[i1 + 2 * _nmax] -= lamda * r01.z;
    list[nlist++] = i1;
  }

  if(vflag || vflag_atom) {
    ENERGY_CFLOAT* shared = &sharedmem[threadIdx.x];
    X_CFLOAT factor = nlist;
    v[0] = lamda * r01.x * r01.x;
    *shared = factor * v[0];
    shared += blockDim.x; //times 2.0 since the reducing function is the same as in force calculations, which adds a factor 0.5
    v[1] = lamda * r01.y * r01.y;
    *shared = factor * v[1];
    shared += blockDim.x;
    v[2] = lamda * r01.z * r01.z;
    *shared = factor * v[2];
    shared += blockDim.x;
    v[3] = lamda * r01.x * r01.y;
    *shared = factor * v[3];
    shared += blockDim.x;
    v[4] = lamda * r01.x * r01.z;
    *shared = factor * v[4];
    shared += blockDim.x;
    v[5] = lamda * r01.y * r01.z;
    *shared = factor * v[5];
    shared += blockDim.x;

    v_tally(vflag, vflag_atom, nlist, list, 2.0, v);
  }
}


__device__ void FixShakeCuda_Shake3(int &vflag, int &vflag_atom, int &m)
{
  int nlist, list[3];
  ENERGY_CFLOAT v[6];
  X_CFLOAT invmass0, invmass1, invmass2;

  // local atom IDs and constraint distances

  int i0 = _map_array[_shake_atom[m]];
  int i1 = _map_array[_shake_atom[m + _nmax]];
  int i2 = _map_array[_shake_atom[m + 2 * _nmax]];
  X_CFLOAT bond1 = _bond_distance[_shake_type[m]];
  X_CFLOAT bond2 = _bond_distance[_shake_type[m + _nmax]];

  // r01 = distance vec between atoms, with PBC

  X_CFLOAT3 r01, r02;

  X_CFLOAT4 x_i0, x_i1, x_i2;
  x_i0 = fetchXType(i0);
  x_i1 = fetchXType(i1);
  x_i2 = fetchXType(i2);

  r01.x = x_i0.x - x_i1.x;
  r01.y = x_i0.y - x_i1.y;
  r01.z = x_i0.z - x_i1.z;
  minimum_image(r01);

  r02.x = x_i0.x - x_i2.x;
  r02.y = x_i0.y - x_i2.y;
  r02.z = x_i0.z - x_i2.z;
  minimum_image(r02);

  // s01 = distance vec after unconstrained update, with PBC

  X_CFLOAT3 s01, s02;
  X_CFLOAT3 xs_i0 = _xshake[i0];
  X_CFLOAT3 xs_i1 = _xshake[i1];
  X_CFLOAT3 xs_i2 = _xshake[i2];

  s01.x = xs_i0.x - xs_i1.x;
  s01.y = xs_i0.y - xs_i1.y;
  s01.z = xs_i0.z - xs_i1.z;
  minimum_image(s01);

  s02.x = xs_i0.x - xs_i2.x;
  s02.y = xs_i0.y - xs_i2.y;
  s02.z = xs_i0.z - xs_i2.z;
  minimum_image(s02);

  // scalar distances between atoms

  X_CFLOAT r01sq = r01.x * r01.x + r01.y * r01.y + r01.z * r01.z;
  X_CFLOAT r02sq = r02.x * r02.x + r02.y * r02.y + r02.z * r02.z;
  X_CFLOAT s01sq = s01.x * s01.x + s01.y * s01.y + s01.z * s01.z;
  X_CFLOAT s02sq = s02.x * s02.x + s02.y * s02.y + s02.z * s02.z;

  // a,b,c = coeffs in quadratic equation for lamda

  if(_rmass_flag) {
    invmass0 = X_F(1.0) / _rmass[i0];
    invmass1 = X_F(1.0) / _rmass[i1];
    invmass2 = X_F(1.0) / _rmass[i2];
  } else {
    invmass0 = X_F(1.0) / _mass[static_cast <int>(x_i0.w)];
    invmass1 = X_F(1.0) / _mass[static_cast <int>(x_i1.w)];
    invmass2 = X_F(1.0) / _mass[static_cast <int>(x_i2.w)];
  }

  X_CFLOAT a11 = X_F(2.0) * (invmass0 + invmass1) *
                (s01.x * r01.x + s01.y * r01.y + s01.z * r01.z);
  X_CFLOAT a12 = X_F(2.0) * invmass0 *
                (s01.x * r02.x + s01.y * r02.y + s01.z * r02.z);
  X_CFLOAT a21 = X_F(2.0) * invmass0 *
                (s02.x * r01.x + s02.y * r01.y + s02.z * r01.z);
  X_CFLOAT a22 = X_F(2.0) * (invmass0 + invmass2) *
                (s02.x * r02.x + s02.y * r02.y + s02.z * r02.z);

  // error check

  X_CFLOAT determ = a11 * a22 - a12 * a21;

  if(determ == X_F(0.0)) _flag[0]++;

  X_CFLOAT determinv = X_F(1.0) / determ;

  X_CFLOAT a11inv = a22 * determinv;
  X_CFLOAT a12inv = -a12 * determinv;
  X_CFLOAT a21inv = -a21 * determinv;
  X_CFLOAT a22inv = a11 * determinv;

  // quadratic correction coeffs

  X_CFLOAT r0102 = (r01.x * r02.x + r01.y * r02.y + r01.z * r02.z);

  X_CFLOAT quad1_0101 = (invmass0 + invmass1) * (invmass0 + invmass1) * r01sq;
  X_CFLOAT quad1_0202 = invmass0 * invmass0 * r02sq;
  X_CFLOAT quad1_0102 = X_F(2.0) * (invmass0 + invmass1) * invmass0 * r0102;

  X_CFLOAT quad2_0202 = (invmass0 + invmass2) * (invmass0 + invmass2) * r02sq;
  X_CFLOAT quad2_0101 = invmass0 * invmass0 * r01sq;
  X_CFLOAT quad2_0102 = X_F(2.0) * (invmass0 + invmass2) * invmass0 * r0102;

  // iterate until converged

  X_CFLOAT lamda01 = X_F(0.0);
  X_CFLOAT lamda02 = X_F(0.0);
  int niter = 0;
  int done = 0;

  X_CFLOAT quad1, quad2, b1, b2, lamda01_new, lamda02_new;

  //maybe all running full loop?
  while(__any(!done) && niter < _max_iter) {
    quad1 = quad1_0101 * lamda01 * lamda01 + quad1_0202 * lamda02 * lamda02 +
            quad1_0102 * lamda01 * lamda02;
    quad2 = quad2_0101 * lamda01 * lamda01 + quad2_0202 * lamda02 * lamda02 +
            quad2_0102 * lamda01 * lamda02;

    b1 = bond1 * bond1 - s01sq - quad1;
    b2 = bond2 * bond2 - s02sq - quad2;

    lamda01_new = a11inv * b1 + a12inv * b2;
    lamda02_new = a21inv * b1 + a22inv * b2;

    done++;
    done = (fabs(lamda01_new - lamda01) > _tolerance) ? 0 : done;
    done = (fabs(lamda02_new - lamda02) > _tolerance) ? 0 : done;


    lamda01 = done < 2 ? lamda01_new : lamda01;
    lamda02 = done < 2 ? lamda02_new : lamda02;
    niter++;
  }

  // update forces if atom is owned by this processor

  lamda01 *= X_F(1.0) / _dtfsq;
  lamda02 *= X_F(1.0) / _dtfsq;


  //attenion: are shake clusters <-> atom unique?
  nlist = 0;

  if(i0 < _nlocal) {
    _f[i0] += lamda01 * r01.x + lamda02 * r02.x;
    _f[i0 + _nmax] += lamda01 * r01.y + lamda02 * r02.y;
    _f[i0 + 2 * _nmax] += lamda01 * r01.z + lamda02 * r02.z;
    list[nlist++] = i0;
  }

  if(i1 < _nlocal) {
    _f[i1] -= lamda01 * r01.x;
    _f[i1 + _nmax] -= lamda01 * r01.y;
    _f[i1 + 2 * _nmax] -= lamda01 * r01.z;
    list[nlist++] = i1;
  }

  if(i2 < _nlocal) {
    _f[i2] -= lamda02 * r02.x;
    _f[i2 + _nmax] -= lamda02 * r02.y;
    _f[i2 + 2 * _nmax] -= lamda02 * r02.z;
    list[nlist++] = i2;
  }

  if(vflag || vflag_atom) {
    ENERGY_CFLOAT* shared = &sharedmem[threadIdx.x];
    X_CFLOAT factor = X_F(2.0) / X_F(3.0) * nlist;
    v[0] = lamda01 * r01.x * r01.x + lamda02 * r02.x * r02.x;
    *shared = factor * v[0];
    shared += blockDim.x; //times 2.0 since the reducing function is the same as in force calculations, which adds a factor 0.5
    v[1] = lamda01 * r01.y * r01.y + lamda02 * r02.y * r02.y;
    *shared = factor * v[1];
    shared += blockDim.x;
    v[2] = lamda01 * r01.z * r01.z + lamda02 * r02.z * r02.z;
    *shared = factor * v[2];
    shared += blockDim.x;
    v[3] = lamda01 * r01.x * r01.y + lamda02 * r02.x * r02.y;
    *shared = factor * v[3];
    shared += blockDim.x;
    v[4] = lamda01 * r01.x * r01.z + lamda02 * r02.x * r02.z;
    *shared = factor * v[4];
    shared += blockDim.x;
    v[5] = lamda01 * r01.y * r01.z + lamda02 * r02.y * r02.z;
    *shared = factor * v[5];
    shared += blockDim.x;

    v_tally(vflag, vflag_atom, nlist, list, 3.0, v);
  }
}

__device__ void FixShakeCuda_Shake4(int &vflag, int &vflag_atom, int &m)
{
  int nlist, list[4];
  ENERGY_CFLOAT v[6];
  X_CFLOAT invmass0, invmass1, invmass2, invmass3;

  // local atom IDs and constraint distances

  int i0 = _map_array[_shake_atom[m]];
  int i1 = _map_array[_shake_atom[m + _nmax]];
  int i2 = _map_array[_shake_atom[m + 2 * _nmax]];
  int i3 = _map_array[_shake_atom[m + 3 * _nmax]];
  X_CFLOAT bond1 = _bond_distance[_shake_type[m]];
  X_CFLOAT bond2 = _bond_distance[_shake_type[m + _nmax]];
  X_CFLOAT bond3 = _bond_distance[_shake_type[m + 2 * _nmax]];

  // r01 = distance vec between atoms, with PBC

  X_CFLOAT3 r01, r02, r03;

  X_CFLOAT4 x_i0, x_i1, x_i2, x_i3;
  x_i0 = fetchXType(i0);
  x_i1 = fetchXType(i1);
  x_i2 = fetchXType(i2);
  x_i3 = fetchXType(i3);

  r01.x = x_i0.x - x_i1.x;
  r01.y = x_i0.y - x_i1.y;
  r01.z = x_i0.z - x_i1.z;
  minimum_image(r01);

  r02.x = x_i0.x - x_i2.x;
  r02.y = x_i0.y - x_i2.y;
  r02.z = x_i0.z - x_i2.z;
  minimum_image(r02);

  r03.x = x_i0.x - x_i3.x;
  r03.y = x_i0.y - x_i3.y;
  r03.z = x_i0.z - x_i3.z;
  minimum_image(r03);

  // s01 = distance vec after unconstrained update, with PBC

  X_CFLOAT3 s01, s02, s03;
  X_CFLOAT3 xs_i0 = _xshake[i0];
  X_CFLOAT3 xs_i1 = _xshake[i1];
  X_CFLOAT3 xs_i2 = _xshake[i2];
  X_CFLOAT3 xs_i3 = _xshake[i3];

  s01.x = xs_i0.x - xs_i1.x;
  s01.y = xs_i0.y - xs_i1.y;
  s01.z = xs_i0.z - xs_i1.z;
  minimum_image(s01);

  s02.x = xs_i0.x - xs_i2.x;
  s02.y = xs_i0.y - xs_i2.y;
  s02.z = xs_i0.z - xs_i2.z;
  minimum_image(s02);

  s03.x = xs_i0.x - xs_i3.x;
  s03.y = xs_i0.y - xs_i3.y;
  s03.z = xs_i0.z - xs_i3.z;
  minimum_image(s03);

  // scalar distances between atoms

  X_CFLOAT r01sq = r01.x * r01.x + r01.y * r01.y + r01.z * r01.z;
  X_CFLOAT r02sq = r02.x * r02.x + r02.y * r02.y + r02.z * r02.z;
  X_CFLOAT r03sq = r03.x * r03.x + r03.y * r03.y + r03.z * r03.z;
  X_CFLOAT s01sq = s01.x * s01.x + s01.y * s01.y + s01.z * s01.z;
  X_CFLOAT s02sq = s02.x * s02.x + s02.y * s02.y + s02.z * s02.z;
  X_CFLOAT s03sq = s03.x * s03.x + s03.y * s03.y + s03.z * s03.z;

  // a,b,c = coeffs in quadratic equation for lamda

  if(_rmass_flag) {
    invmass0 = X_F(1.0) / _rmass[i0];
    invmass1 = X_F(1.0) / _rmass[i1];
    invmass2 = X_F(1.0) / _rmass[i2];
    invmass3 = X_F(1.0) / _rmass[i3];
  } else {
    invmass0 = X_F(1.0) / _mass[static_cast <int>(x_i0.w)];
    invmass1 = X_F(1.0) / _mass[static_cast <int>(x_i1.w)];
    invmass2 = X_F(1.0) / _mass[static_cast <int>(x_i2.w)];
    invmass3 = X_F(1.0) / _mass[static_cast <int>(x_i3.w)];
  }

  X_CFLOAT a11 = X_F(2.0) * (invmass0 + invmass1) *
                (s01.x * r01.x + s01.y * r01.y + s01.z * r01.z);
  X_CFLOAT a12 = X_F(2.0) * invmass0 *
                (s01.x * r02.x + s01.y * r02.y + s01.z * r02.z);
  X_CFLOAT a13 = X_F(2.0) * invmass0 *
                (s01.x * r03.x + s01.y * r03.y + s01.z * r03.z);
  X_CFLOAT a21 = X_F(2.0) * invmass0 *
                (s02.x * r01.x + s02.y * r01.y + s02.z * r01.z);
  X_CFLOAT a22 = X_F(2.0) * (invmass0 + invmass2) *
                (s02.x * r02.x + s02.y * r02.y + s02.z * r02.z);
  X_CFLOAT a23 = X_F(2.0) * (invmass0) *
                (s02.x * r03.x + s02.y * r03.y + s02.z * r03.z);
  X_CFLOAT a31 = X_F(2.0) * (invmass0) *
                (s03.x * r01.x + s03.y * r01.y + s03.z * r01.z);
  X_CFLOAT a32 = X_F(2.0) * (invmass0) *
                (s03.x * r02.x + s03.y * r02.y + s03.z * r02.z);
  X_CFLOAT a33 = X_F(2.0) * (invmass0 + invmass3) *
                (s03.x * r03.x + s03.y * r03.y + s03.z * r03.z);

  // error check

  X_CFLOAT determ = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 -
                   a11 * a23 * a32 - a12 * a21 * a33 - a13 * a22 * a31;

  if(determ == X_F(0.0)) _flag[0]++;

  X_CFLOAT determinv = X_F(1.0) / determ;

  X_CFLOAT a11inv = determinv * (a22 * a33 - a23 * a32);
  X_CFLOAT a12inv = -determinv * (a12 * a33 - a13 * a32);
  X_CFLOAT a13inv = determinv * (a12 * a23 - a13 * a22);
  X_CFLOAT a21inv = -determinv * (a21 * a33 - a23 * a31);
  X_CFLOAT a22inv = determinv * (a11 * a33 - a13 * a31);
  X_CFLOAT a23inv = -determinv * (a11 * a23 - a13 * a21);
  X_CFLOAT a31inv = determinv * (a21 * a32 - a22 * a31);
  X_CFLOAT a32inv = -determinv * (a11 * a32 - a12 * a31);
  X_CFLOAT a33inv = determinv * (a11 * a22 - a12 * a21);

  // quadratic correction coeffs

  X_CFLOAT r0102 = (r01.x * r02.x + r01.y * r02.y + r01.z * r02.z);
  X_CFLOAT r0103 = (r01.x * r03.x + r01.y * r03.y + r01.z * r03.z);
  X_CFLOAT r0203 = (r02.x * r03.x + r02.y * r03.y + r02.z * r03.z);

  X_CFLOAT quad1_0101 = (invmass0 + invmass1) * (invmass0 + invmass1) * r01sq;
  X_CFLOAT quad1_0202 = invmass0 * invmass0 * r02sq;
  X_CFLOAT quad1_0303 = invmass0 * invmass0 * r03sq;
  X_CFLOAT quad1_0102 = X_F(2.0) * (invmass0 + invmass1) * invmass0 * r0102;
  X_CFLOAT quad1_0103 = X_F(2.0) * (invmass0 + invmass1) * invmass0 * r0103;
  X_CFLOAT quad1_0203 = X_F(2.0) * invmass0 * invmass0 * r0203;

  X_CFLOAT quad2_0101 = invmass0 * invmass0 * r01sq;
  X_CFLOAT quad2_0202 = (invmass0 + invmass2) * (invmass0 + invmass2) * r02sq;
  X_CFLOAT quad2_0303 = invmass0 * invmass0 * r03sq;
  X_CFLOAT quad2_0102 = X_F(2.0) * (invmass0 + invmass2) * invmass0 * r0102;
  X_CFLOAT quad2_0103 = X_F(2.0) * invmass0 * invmass0 * r0103;
  X_CFLOAT quad2_0203 = X_F(2.0) * (invmass0 + invmass2) * invmass0 * r0203;

  X_CFLOAT quad3_0101 = invmass0 * invmass0 * r01sq;
  X_CFLOAT quad3_0202 = invmass0 * invmass0 * r02sq;
  X_CFLOAT quad3_0303 = (invmass0 + invmass3) * (invmass0 + invmass3) * r03sq;
  X_CFLOAT quad3_0102 = X_F(2.0) * invmass0 * invmass0 * r0102;
  X_CFLOAT quad3_0103 = X_F(2.0) * (invmass0 + invmass3) * invmass0 * r0103;
  X_CFLOAT quad3_0203 = X_F(2.0) * (invmass0 + invmass3) * invmass0 * r0203;
  // iterate until converged

  X_CFLOAT lamda01 = X_F(0.0);
  X_CFLOAT lamda02 = X_F(0.0);
  X_CFLOAT lamda03 = X_F(0.0);
  int niter = 0;
  int done = 0;

  X_CFLOAT quad1, quad2, quad3, b1, b2, b3, lamda01_new, lamda02_new, lamda03_new;

  //maybe all running full loop?
  while(__any(!done) && niter < _max_iter) {
    quad1 = quad1_0101 * lamda01 * lamda01 +
            quad1_0202 * lamda02 * lamda02 +
            quad1_0303 * lamda03 * lamda03 +
            quad1_0102 * lamda01 * lamda02 +
            quad1_0103 * lamda01 * lamda03 +
            quad1_0203 * lamda02 * lamda03;

    quad2 = quad2_0101 * lamda01 * lamda01 +
            quad2_0202 * lamda02 * lamda02 +
            quad2_0303 * lamda03 * lamda03 +
            quad2_0102 * lamda01 * lamda02 +
            quad2_0103 * lamda01 * lamda03 +
            quad2_0203 * lamda02 * lamda03;

    quad3 = quad3_0101 * lamda01 * lamda01 +
            quad3_0202 * lamda02 * lamda02 +
            quad3_0303 * lamda03 * lamda03 +
            quad3_0102 * lamda01 * lamda02 +
            quad3_0103 * lamda01 * lamda03 +
            quad3_0203 * lamda02 * lamda03;

    b1 = bond1 * bond1 - s01sq - quad1;
    b2 = bond2 * bond2 - s02sq - quad2;
    b3 = bond3 * bond3 - s03sq - quad3;

    lamda01_new = a11inv * b1 + a12inv * b2 + a13inv * b3;
    lamda02_new = a21inv * b1 + a22inv * b2 + a23inv * b3;
    lamda03_new = a31inv * b1 + a32inv * b2 + a33inv * b3;

    done++;
    done = (fabs(lamda01_new - lamda01) > _tolerance) ? 0 : done;
    done = (fabs(lamda02_new - lamda02) > _tolerance) ? 0 : done;
    done = (fabs(lamda03_new - lamda03) > _tolerance) ? 0 : done;

    lamda01 = done < 2 ? lamda01_new : lamda01;
    lamda02 = done < 2 ? lamda02_new : lamda02;
    lamda03 = done < 2 ? lamda03_new : lamda03;
    niter++;
  }

  // update forces if atom is owned by this processor

  lamda01 *= X_F(1.0) / _dtfsq;
  lamda02 *= X_F(1.0) / _dtfsq;
  lamda03 *= X_F(1.0) / _dtfsq;


  //attenion: are shake clusters <-> atom unique?
  nlist = 0;

  if(i0 < _nlocal) {
    _f[i0] 			+= lamda01 * r01.x + lamda02 * r02.x + lamda03 * r03.x;
    _f[i0 + _nmax] 	+= lamda01 * r01.y + lamda02 * r02.y + lamda03 * r03.y;
    _f[i0 + 2 * _nmax] 	+= lamda01 * r01.z + lamda02 * r02.z + lamda03 * r03.z;
    list[nlist++] = i0;
  }

  if(i1 < _nlocal) {
    _f[i1] -= lamda01 * r01.x;
    _f[i1 + _nmax] -= lamda01 * r01.y;
    _f[i1 + 2 * _nmax] -= lamda01 * r01.z;
    list[nlist++] = i1;
  }

  if(i2 < _nlocal) {
    _f[i2] -= lamda02 * r02.x;
    _f[i2 + _nmax] -= lamda02 * r02.y;
    _f[i2 + 2 * _nmax] -= lamda02 * r02.z;
    list[nlist++] = i2;
  }

  if(i3 < _nlocal) {
    _f[i3] -= lamda03 * r03.x;
    _f[i3 + _nmax] -= lamda03 * r03.y;
    _f[i3 + 2 * _nmax] -= lamda03 * r03.z;
    list[nlist++] = i3;
  }

  if(vflag || vflag_atom) {
    ENERGY_CFLOAT* shared = &sharedmem[threadIdx.x];
    X_CFLOAT factor = X_F(2.0) / X_F(4.0) * nlist;
    v[0] = lamda01 * r01.x * r01.x + lamda02 * r02.x * r02.x + lamda03 * r03.x * r03.x;
    *shared = factor * v[0];
    shared += blockDim.x; //times 2.0 since the reducing function is the same as in force calculations, which adds a factor 0.5
    v[1] = lamda01 * r01.y * r01.y + lamda02 * r02.y * r02.y + lamda03 * r03.y * r03.y;
    *shared = factor * v[1];
    shared += blockDim.x;
    v[2] = lamda01 * r01.z * r01.z + lamda02 * r02.z * r02.z + lamda03 * r03.z * r03.z;
    *shared = factor * v[2];
    shared += blockDim.x;
    v[3] = lamda01 * r01.x * r01.y + lamda02 * r02.x * r02.y + lamda03 * r03.x * r03.y;
    *shared = factor * v[3];
    shared += blockDim.x;
    v[4] = lamda01 * r01.x * r01.z + lamda02 * r02.x * r02.z + lamda03 * r03.x * r03.z;
    *shared = factor * v[4];
    shared += blockDim.x;
    v[5] = lamda01 * r01.y * r01.z + lamda02 * r02.y * r02.z + lamda03 * r03.y * r03.z;
    *shared = factor * v[5];
    shared += blockDim.x;

    v_tally(vflag, vflag_atom, nlist, list, 4.0, v);
  }
}

__device__ void FixShakeCuda_Shake3Angle(int &vflag, int &vflag_atom, int &m)
{
  int nlist, list[3];
  ENERGY_CFLOAT v[6];
  X_CFLOAT invmass0, invmass1, invmass2;

  // local atom IDs and constraint distances

  int i0 = _map_array[_shake_atom[m]];
  int i1 = _map_array[_shake_atom[m + _nmax]];
  int i2 = _map_array[_shake_atom[m + 2 * _nmax]];
  X_CFLOAT bond1 = _bond_distance[_shake_type[m]];
  X_CFLOAT bond2 = _bond_distance[_shake_type[m + _nmax]];
  X_CFLOAT bond12 = _angle_distance[_shake_type[m + 2 * _nmax]];

  // r01 = distance vec between atoms, with PBC

  X_CFLOAT3 r01, r02, r12;

  X_CFLOAT4 x_i0, x_i1, x_i2;
  x_i0 = fetchXType(i0);
  x_i1 = fetchXType(i1);
  x_i2 = fetchXType(i2);

  r01.x = x_i0.x - x_i1.x;
  r01.y = x_i0.y - x_i1.y;
  r01.z = x_i0.z - x_i1.z;
  minimum_image(r01);

  r02.x = x_i0.x - x_i2.x;
  r02.y = x_i0.y - x_i2.y;
  r02.z = x_i0.z - x_i2.z;
  minimum_image(r02);

  r12.x = x_i1.x - x_i2.x;
  r12.y = x_i1.y - x_i2.y;
  r12.z = x_i1.z - x_i2.z;
  minimum_image(r12);

  // s01 = distance vec after unconstrained update, with PBC

  X_CFLOAT3 s01, s02, s12;
  X_CFLOAT3 xs_i0 = _xshake[i0];
  X_CFLOAT3 xs_i1 = _xshake[i1];
  X_CFLOAT3 xs_i2 = _xshake[i2];

  s01.x = xs_i0.x - xs_i1.x;
  s01.y = xs_i0.y - xs_i1.y;
  s01.z = xs_i0.z - xs_i1.z;
  minimum_image(s01);

  s02.x = xs_i0.x - xs_i2.x;
  s02.y = xs_i0.y - xs_i2.y;
  s02.z = xs_i0.z - xs_i2.z;
  minimum_image(s02);

  s12.x = xs_i1.x - xs_i2.x;
  s12.y = xs_i1.y - xs_i2.y;
  s12.z = xs_i1.z - xs_i2.z;
  minimum_image(s12);

  // scalar distances between atoms

  X_CFLOAT r01sq = r01.x * r01.x + r01.y * r01.y + r01.z * r01.z;
  X_CFLOAT r02sq = r02.x * r02.x + r02.y * r02.y + r02.z * r02.z;
  X_CFLOAT r12sq = r12.x * r12.x + r12.y * r12.y + r12.z * r12.z;
  X_CFLOAT s01sq = s01.x * s01.x + s01.y * s01.y + s01.z * s01.z;
  X_CFLOAT s02sq = s02.x * s02.x + s02.y * s02.y + s02.z * s02.z;
  X_CFLOAT s12sq = s12.x * s12.x + s12.y * s12.y + s12.z * s12.z;

  // a,b,c = coeffs in quadratic equation for lamda

  if(_rmass_flag) {
    invmass0 = X_F(1.0) / _rmass[i0];
    invmass1 = X_F(1.0) / _rmass[i1];
    invmass2 = X_F(1.0) / _rmass[i2];
  } else {
    invmass0 = X_F(1.0) / _mass[static_cast <int>(x_i0.w)];
    invmass1 = X_F(1.0) / _mass[static_cast <int>(x_i1.w)];
    invmass2 = X_F(1.0) / _mass[static_cast <int>(x_i2.w)];
  }

  X_CFLOAT a11 = X_F(2.0) * (invmass0 + invmass1) *
                (s01.x * r01.x + s01.y * r01.y + s01.z * r01.z);
  X_CFLOAT a12 = X_F(2.0) * invmass0 *
                (s01.x * r02.x + s01.y * r02.y + s01.z * r02.z);
  X_CFLOAT a13 = - X_F(2.0) * invmass1 *
                (s01.x * r12.x + s01.y * r12.y + s01.z * r12.z);
  X_CFLOAT a21 = X_F(2.0) * invmass0 *
                (s02.x * r01.x + s02.y * r01.y + s02.z * r01.z);
  X_CFLOAT a22 = X_F(2.0) * (invmass0 + invmass2) *
                (s02.x * r02.x + s02.y * r02.y + s02.z * r02.z);
  X_CFLOAT a23 = X_F(2.0) * invmass2 *
                (s02.x * r12.x + s02.y * r12.y + s02.z * r12.z);
  X_CFLOAT a31 = - X_F(2.0) * invmass1 *
                (s12.x * r01.x + s12.y * r01.y + s12.z * r01.z);
  X_CFLOAT a32 = X_F(2.0) * invmass2 *
                (s12.x * r02.x + s12.y * r02.y + s12.z * r02.z);
  X_CFLOAT a33 = X_F(2.0) * (invmass1 + invmass2) *
                (s12.x * r12.x + s12.y * r12.y + s12.z * r12.z);

  // inverse of matrix

  X_CFLOAT determ = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 -
                   a11 * a23 * a32 - a12 * a21 * a33 - a13 * a22 * a31;

  if(determ == X_F(0.0)) _flag[0]++;

  X_CFLOAT determinv = X_F(1.0) / determ;

  X_CFLOAT a11inv = determinv * (a22 * a33 - a23 * a32);
  X_CFLOAT a12inv = -determinv * (a12 * a33 - a13 * a32);
  X_CFLOAT a13inv = determinv * (a12 * a23 - a13 * a22);
  X_CFLOAT a21inv = -determinv * (a21 * a33 - a23 * a31);
  X_CFLOAT a22inv = determinv * (a11 * a33 - a13 * a31);
  X_CFLOAT a23inv = -determinv * (a11 * a23 - a13 * a21);
  X_CFLOAT a31inv = determinv * (a21 * a32 - a22 * a31);
  X_CFLOAT a32inv = -determinv * (a11 * a32 - a12 * a31);
  X_CFLOAT a33inv = determinv * (a11 * a22 - a12 * a21);

  // quadratic correction coeffs

  X_CFLOAT r0102 = (r01.x * r02.x + r01.y * r02.y + r01.z * r02.z);
  X_CFLOAT r0112 = (r01.x * r12.x + r01.y * r12.y + r01.z * r12.z);
  X_CFLOAT r0212 = (r02.x * r12.x + r02.y * r12.y + r02.z * r12.z);

  X_CFLOAT quad1_0101 = (invmass0 + invmass1) * (invmass0 + invmass1) * r01sq;
  X_CFLOAT quad1_0202 = invmass0 * invmass0 * r02sq;
  X_CFLOAT quad1_1212 = invmass1 * invmass1 * r12sq;
  X_CFLOAT quad1_0102 = X_F(2.0) * (invmass0 + invmass1) * invmass0 * r0102;
  X_CFLOAT quad1_0112 = - X_F(2.0) * (invmass0 + invmass1) * invmass1 * r0112;
  X_CFLOAT quad1_0212 = - X_F(2.0) * invmass0 * invmass1 * r0212;

  X_CFLOAT quad2_0101 = invmass0 * invmass0 * r01sq;
  X_CFLOAT quad2_0202 = (invmass0 + invmass2) * (invmass0 + invmass2) * r02sq;
  X_CFLOAT quad2_1212 = invmass2 * invmass2 * r12sq;
  X_CFLOAT quad2_0102 = X_F(2.0) * (invmass0 + invmass2) * invmass0 * r0102;
  X_CFLOAT quad2_0112 = X_F(2.0) * invmass0 * invmass2 * r0112;
  X_CFLOAT quad2_0212 = X_F(2.0) * (invmass0 + invmass2) * invmass2 * r0212;

  X_CFLOAT quad3_0101 = invmass1 * invmass1 * r01sq;
  X_CFLOAT quad3_0202 = invmass2 * invmass2 * r02sq;
  X_CFLOAT quad3_1212 = (invmass1 + invmass2) * (invmass1 + invmass2) * r12sq;
  X_CFLOAT quad3_0102 = - X_F(2.0) * invmass1 * invmass2 * r0102;
  X_CFLOAT quad3_0112 = - X_F(2.0) * (invmass1 + invmass2) * invmass1 * r0112;
  X_CFLOAT quad3_0212 = X_F(2.0) * (invmass1 + invmass2) * invmass2 * r0212;
  // iterate until converged

  X_CFLOAT lamda01 = X_F(0.0);
  X_CFLOAT lamda02 = X_F(0.0);
  X_CFLOAT lamda12 = X_F(0.0);
  int niter = 0;
  int done = 0;

  X_CFLOAT quad1, quad2, quad3, b1, b2, b3, lamda01_new, lamda02_new, lamda12_new;

  //maybe all running full loop?
  while(__any(!done) && niter < _max_iter) {
    quad1 = quad1_0101 * lamda01 * lamda01 +
            quad1_0202 * lamda02 * lamda02 +
            quad1_1212 * lamda12 * lamda12 +
            quad1_0102 * lamda01 * lamda02 +
            quad1_0112 * lamda01 * lamda12 +
            quad1_0212 * lamda02 * lamda12;

    quad2 = quad2_0101 * lamda01 * lamda01 +
            quad2_0202 * lamda02 * lamda02 +
            quad2_1212 * lamda12 * lamda12 +
            quad2_0102 * lamda01 * lamda02 +
            quad2_0112 * lamda01 * lamda12 +
            quad2_0212 * lamda02 * lamda12;

    quad3 = quad3_0101 * lamda01 * lamda01 +
            quad3_0202 * lamda02 * lamda02 +
            quad3_1212 * lamda12 * lamda12 +
            quad3_0102 * lamda01 * lamda02 +
            quad3_0112 * lamda01 * lamda12 +
            quad3_0212 * lamda02 * lamda12;

    b1 = bond1 * bond1 - s01sq - quad1;
    b2 = bond2 * bond2 - s02sq - quad2;
    b3 = bond12 * bond12 - s12sq - quad3;

    lamda01_new = a11inv * b1 + a12inv * b2 + a13inv * b3;
    lamda02_new = a21inv * b1 + a22inv * b2 + a23inv * b3;
    lamda12_new = a31inv * b1 + a32inv * b2 + a33inv * b3;

    done++;
    done = (fabs(lamda01_new - lamda01) > _tolerance) ? 0 : done;
    done = (fabs(lamda02_new - lamda02) > _tolerance) ? 0 : done;
    done = (fabs(lamda12_new - lamda12) > _tolerance) ? 0 : done;

    lamda01 = done < 2 ? lamda01_new : lamda01;
    lamda02 = done < 2 ? lamda02_new : lamda02;
    lamda12 = done < 2 ? lamda12_new : lamda12;
    niter++;
  }

  // update forces if atom is owned by this processor

  lamda01 *= X_F(1.0) / _dtfsq;
  lamda02 *= X_F(1.0) / _dtfsq;
  lamda12 *= X_F(1.0) / _dtfsq;


  //attenion: are shake clusters <-> atom unique?
  nlist = 0;

  if(i0 < _nlocal) {
    _f[i0] 			+= lamda01 * r01.x + lamda02 * r02.x;
    _f[i0 + _nmax] 	+= lamda01 * r01.y + lamda02 * r02.y;
    _f[i0 + 2 * _nmax] 	+= lamda01 * r01.z + lamda02 * r02.z;
    list[nlist++] = i0;
  }

  if(i1 < _nlocal) {
    _f[i1] 			-= lamda01 * r01.x - lamda12 * r12.x;
    _f[i1 + _nmax] 	-= lamda01 * r01.y - lamda12 * r12.y;
    _f[i1 + 2 * _nmax] 	-= lamda01 * r01.z - lamda12 * r12.z;
    list[nlist++] = i1;
  }

  if(i2 < _nlocal) {
    _f[i2] 			-= lamda02 * r02.x + lamda12 * r12.x;
    _f[i2 + _nmax] 	-= lamda02 * r02.y + lamda12 * r12.y;
    _f[i2 + 2 * _nmax] 	-= lamda02 * r02.z + lamda12 * r12.z;
    list[nlist++] = i2;
  }

  if(vflag || vflag_atom) {
    ENERGY_CFLOAT* shared = &sharedmem[threadIdx.x];
    X_CFLOAT factor = X_F(2.0) / X_F(3.0) * nlist;
    v[0] = lamda01 * r01.x * r01.x + lamda02 * r02.x * r02.x + lamda12 * r12.x * r12.x;
    *shared = factor * v[0];
    shared += blockDim.x; //times 2.0 since the reducing function is the same as in force calculations, which adds a factor 0.5
    v[1] = lamda01 * r01.y * r01.y + lamda02 * r02.y * r02.y + lamda12 * r12.y * r12.y;
    *shared = factor * v[1];
    shared += blockDim.x;
    v[2] = lamda01 * r01.z * r01.z + lamda02 * r02.z * r02.z + lamda12 * r12.z * r12.z;
    *shared = factor * v[2];
    shared += blockDim.x;
    v[3] = lamda01 * r01.x * r01.y + lamda02 * r02.x * r02.y + lamda12 * r12.x * r12.y;
    *shared = factor * v[3];
    shared += blockDim.x;
    v[4] = lamda01 * r01.x * r01.z + lamda02 * r02.x * r02.z + lamda12 * r12.x * r12.z;
    *shared = factor * v[4];
    shared += blockDim.x;
    v[5] = lamda01 * r01.y * r01.z + lamda02 * r02.y * r02.z + lamda12 * r12.y * r12.z;
    *shared = factor * v[5];
    shared += blockDim.x;

    v_tally(vflag, vflag_atom, nlist, list, 3.0, v);
  }
}

__global__ void FixShakeCuda_Shake_Kernel(int vflag, int vflag_atom, int* list, int nlist)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < nlist) {

    int m = list[i];
    int sflag = _shake_flag[m];

    if(sflag == 2) FixShakeCuda_Shake2(vflag, vflag_atom, m);
    else if(sflag == 3) FixShakeCuda_Shake3(vflag, vflag_atom, m);
    else if(sflag == 4) FixShakeCuda_Shake4(vflag, vflag_atom, m);
    else FixShakeCuda_Shake3Angle(vflag, vflag_atom, m);
  } else {
    ENERGY_CFLOAT* shared = &sharedmem[threadIdx.x];
    *shared = ENERGY_F(0.0);
    shared += blockDim.x;
    *shared = ENERGY_F(0.0);
    shared += blockDim.x;
    *shared = ENERGY_F(0.0);
    shared += blockDim.x;
    *shared = ENERGY_F(0.0);
    shared += blockDim.x;
    *shared = ENERGY_F(0.0);
    shared += blockDim.x;
    *shared = ENERGY_F(0.0);
  }

  if(vflag) {
    __syncthreads();
    int eflag = 0;
    PairVirialCompute_A_Kernel(eflag, vflag);
  }

}

__global__ void FixShakeCuda_PackComm_Kernel(int* sendlist, int n, int maxlistlength, int iswap, X_CFLOAT dx, X_CFLOAT dy, X_CFLOAT dz)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  int* list = sendlist + iswap * maxlistlength;

  if(i < n) {
    int j = list[i];

    if(j > _nmax) _flag[0] = 1;

    X_CFLOAT3 xs = _xshake[j];
    ((X_CFLOAT*) _buffer)[i] = xs.x + dx;
    ((X_CFLOAT*) _buffer)[i + 1 * n] = xs.y + dy;
    ((X_CFLOAT*) _buffer)[i + 2 * n] = xs.z + dz;
  }

}

__global__ void FixShakeCuda_PackComm_Self_Kernel(int* sendlist, int n, int maxlistlength, int iswap, X_CFLOAT dx, X_CFLOAT dy, X_CFLOAT dz, int first)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  int* list = sendlist + iswap * maxlistlength;

  if(i < n) {
    int j = list[i];

    if(j > _nmax) _flag[0] = 1;

    X_CFLOAT3 xs = _xshake[j];
    xs.x += dx;
    xs.y += dy;
    xs.z += dz;
    _xshake[i + first] = xs;
  }

}

__global__ void FixShakeCuda_UnpackComm_Kernel(int n, int first)
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if(i < n) {
    X_CFLOAT3 xs;
    xs.x = ((X_CFLOAT*) _buffer)[i];
    xs.y = ((X_CFLOAT*) _buffer)[i + 1 * n];
    xs.z = ((X_CFLOAT*) _buffer)[i + 2 * n];
    _xshake[i + first] = xs;
  }
}

