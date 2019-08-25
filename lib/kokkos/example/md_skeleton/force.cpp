/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <system.h>
#include <cstdio>


/* Simple Lennard Jones Force Kernel using neighborlists
 * Calculates for every pair of atoms (i,j) with distance smaller r_cut
 * f_ij = 4*epsilon * ( (sigma/r_ij)^12 - (sigma/r_ij)^6 )
 * where r_ij is the distance of atoms (i,j).
 * The force on atom i is the sum over f_ij:
 * f_i = sum_j (f_ij)
 * Neighborlists are used in order to pre calculate which atoms j are
 * close enough to i to be able to contribute. By choosing a larger neighbor
 * cutoff then the force cutoff, the neighbor list can be reused several times
 * (typically 10 - 100).
 */

struct ForceFunctor {

  typedef t_x_array::execution_space execution_space; //Device Type for running the kernel
  typedef double2 value_type; // When energy calculation is requested return energy, and virial

  t_x_array_randomread x;       //atom positions
  t_f_array f;                  //atom forces
  t_int_1d_const numneigh;      //number of neighbors per atom
  t_neighbors_const neighbors;  //neighborlist
  double cutforcesq;            //force cutoff
  double epsilon;               //Potential parameter
  double sigma6;                //Potential parameter


  ForceFunctor(System s) {
    x = s.d_x;
    f = s.f;
    numneigh = s.numneigh;
    neighbors = s.neighbors;
    cutforcesq = s.force_cutsq;
    epsilon = 1.0;
    sigma6 = 1.0;
  }

  /* Operator for not calculating energy and virial */

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &i) const {
    force<0>(i);
  }

  /* Operator for calculating energy and virial */

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &i, double2 &energy_virial) const {
    double2 ev = force<1>(i);
    energy_virial.x += ev.x;
    energy_virial.y += ev.y;
  }

  template<int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  double2 force(const int &i) const
  {
    const int numneighs = numneigh[i];
    const double xtmp = x(i, 0);
    const double ytmp = x(i, 1);
    const double ztmp = x(i, 2);
    double fix = 0;
    double fiy = 0;
    double fiz = 0;
    double energy = 0;
    double virial = 0;

    //pragma simd forces vectorization (ignoring the performance objections of the compiler)
    //give hint to compiler that fix, fiy and fiz are used for reduction only

  #ifdef USE_SIMD
    #pragma simd reduction (+: fix,fiy,fiz,energy,virial)
  #endif
    for(int k = 0; k < numneighs; k++) {
      const int j = neighbors(i, k);
      const double delx = xtmp - x(j, 0);
      const double dely = ytmp - x(j, 1);
      const double delz = ztmp - x(j, 2);
      const double rsq = delx * delx + dely * dely + delz * delz;

      //if(i==0) printf("%i %i %lf %lf\n",i,j,rsq,cutforcesq);
      if(rsq < cutforcesq) {
        const double sr2 = 1.0 / rsq;
        const double sr6 = sr2 * sr2 * sr2  * sigma6;
        const double force = 48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon;
        fix += delx * force;
        fiy += dely * force;
        fiz += delz * force;

        if(EVFLAG) {
          energy += sr6 * (sr6 - 1.0) * epsilon;
          virial += delx * delx * force + dely * dely * force + delz * delz * force;
        }
      }
    }

    f(i, 0) += fix;
    f(i, 1) += fiy;
    f(i, 2) += fiz;

    double2 energy_virial ;
    energy_virial.x = 4.0 * energy ;
    energy_virial.y = 0.5 * virial ;
    return energy_virial;
  }

  /* init and join functions when doing the reduction to obtain energy and virial */

  KOKKOS_FUNCTION
  static void init(volatile value_type &update) {
    update.x = update.y = 0;
  }
  KOKKOS_FUNCTION
  static void join(volatile value_type &update ,
                   const volatile value_type &source) {
    update.x += source.x ;
    update.y += source.y ;
  }

};


/* Calling function */

double2 force(System &s,int evflag) {

  ForceFunctor f(s);

  double2 ev ; ev.x = 0 ; ev.y = 0 ;
  if(!evflag)
    Kokkos::parallel_for(s.nlocal,f);
  else
    Kokkos::parallel_reduce(s.nlocal,f,ev);

  execution_space().fence();
  return ev;
}

