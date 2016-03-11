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

// Mersenne Twister (MT19937) pseudo random number generator:
// M. Matsumoto & T. Nishimura,
// ACM Transactions on Modeling and Computer Simulation,
// vol. 8, no. 1, 1998, pp. 3-30.
//
// Uses the Marsaglia RNG in RanMars to generate the initial seeds

#include <math.h>
#include "random_mt.h"
#include "random_mars.h"
#include "math_inline.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MT_A 0x9908B0DF
#define MT_B 0x9D2C5680
#define MT_C 0xEFC60000

static const double xn = 3.6541528853610088;
static const double A  = 0.00492867323399;

/* ---------------------------------------------------------------------- */

RanMT::RanMT(LAMMPS *lmp, int seed) : Pointers(lmp)
{
  int i;
  const uint32_t f = 1812433253UL;
  _save = _second = 0;

  if (seed <= 0 || seed > 900000000)
    error->one(FLERR,"Invalid seed for Mersenne Twister random # generator");

  // start minimal initialization
  _m[0] = seed;
  _idx = MT_N-1;
  for (i=1; i < MT_N; ++i)
    _m[i] = (f * (_m[i-1] ^ (_m[i-1] >> 30)) + i);

  // to seed the RNG some more using a second RNG
  RanMars rng(lmp,seed);

  for (i=0; i < MT_N-1; ++i)
    _m[i+1] = (_m[i+1] ^ ((_m[i] ^ (_m[i] >> 30)) * 1664525UL))
      + (uint32_t) (rng.uniform()* (1U<<31)) + i;

  _m[0] = _m[MT_N-1];

  for (i=0; i < MT_N-1; ++i)
    _m[i+1] = (_m[i+1] ^ ((_m[i] ^ (_m[i] >> 30)) * 1566083941UL))-i-1;

  _m[0] = 0x80000000UL;

  // randomize one more turn
  _idx = 0;
  for (i=0; i < MT_N-1; ++i) _randomize();

  // build ziggurath tables
  _nlayers = 256;
  _xlayers = new double[_nlayers];
  _ylayers = new double[_nlayers];
  _xlayers[0] = xn;
  _ylayers[0] = exp(-0.5*xn*xn);
  for (int i=1; i < _nlayers; ++i) {
    _ylayers[i] = A/_xlayers[i-1] + _ylayers[i-1];
    _xlayers[i] = sqrt(-2.0*log(_ylayers[i]));
  }
  _xlayers[_nlayers-1] = 0.0;
  _ylayers[_nlayers-1] = 1.0;
}

/* ---------------------------------------------------------------------- */

RanMT::~RanMT()
{
  delete[] _xlayers;
  delete[] _ylayers;
}


/* ----------------------------------------------------------------------
   grab 32bits of randomness
------------------------------------------------------------------------- */

uint32_t RanMT::_randomize() {
  uint32_t r;

  if (_idx >= MT_N) {
    // fill the entire status array with new data in one sweep
    const uint32_t LMASK = (1LU << MT_R) - 1;  // Lower MT_R bits
    const uint32_t UMASK = 0xFFFFFFFF << MT_R; // Upper (32 - MT_R) bits
    static const uint32_t magic[2] = {0, MT_A};
    const int diff = MT_N-MT_M;
    int i;

    for (i=0; i < diff; ++i) {
      r = (_m[i] & UMASK) | (_m[i+1] & LMASK);
      _m[i] = _m[i+MT_M] ^ (r >> 1) ^ magic[r & 1];}

    for (i=diff; i < MT_N-1; ++i) {
      r = (_m[i] & UMASK) | (_m[i+1] & LMASK);
      _m[i] = _m[i-diff] ^ (r >> 1) ^ magic[r & 1];}

    r = (_m[MT_N-1] & UMASK) | (_m[0] & LMASK);
    _m[MT_N-1] = _m[MT_M-1] ^ (r >> 1) ^ magic[r & 1];
    _idx = 0;
  }
  r = _m[_idx++];

  r ^=  r >> MT_U;
  r ^= (r << MT_S) & MT_B;
  r ^= (r << MT_T) & MT_C;
  r ^=  r >> MT_L;

  return r;
}

/* ----------------------------------------------------------------------
   uniform distributed RN. just grab a 32bit integer and convert to double
------------------------------------------------------------------------- */
static const double conv_u32int = 1.0 / (256.0*256.0*256.0*256.0);
double RanMT::uniform()
{
  double uni = (double) _randomize();
  return uni*conv_u32int;
}

/* ----------------------------------------------------------------------
   gaussian distributed RNG
------------------------------------------------------------------------- */

double RanMT::gaussian()
{
  double first,v1,v2,rsq,fac;

  if (!_save) {
    int again = 1;
    while (again) {
      v1 = 2.0*uniform()-1.0;
      v2 = 2.0*uniform()-1.0;
      rsq = v1*v1 + v2*v2;
      if (rsq < 1.0 && rsq != 0.0) again = 0;
    }
    // fac = sqrt(-2.0*log(rsq)/rsq);
    fac = MathInline::sqrtlgx2byx(rsq);
    _second = v1*fac;
    first = v2*fac;
    _save = 1;
  } else {
    first = _second;
    _save = 0;
  }
  return first;
}

/* ----------------------------------------------------------------------
   gaussian distributed RNG via ziggurat method
------------------------------------------------------------------------- */

double RanMT::gaussian_z()
{
  double x,y;

  int n = static_cast<int>(static_cast<double>(_nlayers)*uniform());
  if (n == 0) {                 // special case: tail
    const double xbound = A / _ylayers[0];
    x = 2.0*xbound * uniform() - xbound;
    if (fabs(x) < _xlayers[0]) {
      return x;
    } else {
      double rv;
      do {
        const double x0 = _xlayers[0];
        x = -log(uniform()) / x0;
        y = -log(uniform());
        rv = (uniform() < 0.5) ? -x0-x : x0+x;
      } while (2*y < x*x);
      return rv;
    }
  } else if (n == _nlayers-1) {  // special case: top
    const double xbound = _xlayers[n-1];
    x = 2.0*xbound * uniform() - xbound;
    const double delta = _ylayers[n]-_ylayers[n-1];
    y = delta*uniform() + _ylayers[n-1];
    
    if (y < exp(-0.5*x*x))
      return x;
    else
      return gaussian_z();
  } else {                      // normal case
    const double xbound = _xlayers[n];
    x = 2.0*xbound * uniform() - xbound;
    if (fabs(x) < _xlayers[n+1])
      return x;

    const double delta = _ylayers[n+1]-_ylayers[n];
    y = delta*uniform() + _ylayers[n];
    if (y < exp(-0.5*x*x))
      return x;
    else
      return gaussian_z();    
  }
}
