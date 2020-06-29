/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#include <Kokkos_UnorderedMap.hpp>

namespace Kokkos {
namespace Impl {

uint32_t find_hash_size(uint32_t size) {
  if (size == 0u) return 0u;

  // these primes try to preserve randomness of hash
  static const uint32_t primes[] = {
      3,         7,         13,        23,        53,        97,
      193,       389,       769,       1543,      2237,      2423,
      2617,      2797,      2999,      3167,      3359,      3539,
      3727,      3911,      4441,      4787,      5119,      5471,
      5801,      6143,      6521,      6827,      7177,      7517,
      7853,      8887,      9587,      10243,     10937,     11617,
      12289,     12967,     13649,     14341,     15013,     15727,
      17749,     19121,     20479,     21859,     23209,     24593,
      25939,     27329,     28669,     30047,     31469,     35507,
      38231,     40961,     43711,     46439,     49157,     51893,
      54617,     57347,     60077,     62801,     70583,     75619,
      80669,     85703,     90749,     95783,     100823,    105871,
      110909,    115963,    120997,    126031,    141157,    151237,
      161323,    171401,    181499,    191579,    201653,    211741,
      221813,    231893,    241979,    252079,    282311,    302483,
      322649,    342803,    362969,    383143,    403301,    423457,
      443629,    463787,    483953,    504121,    564617,    604949,
      645313,    685609,    725939,    766273,    806609,    846931,
      887261,    927587,    967919,    1008239,   1123477,   1198397,
      1273289,   1348177,   1423067,   1497983,   1572869,   1647761,
      1722667,   1797581,   1872461,   1947359,   2022253,   2246953,
      2396759,   2546543,   2696363,   2846161,   2995973,   3145739,
      3295541,   3445357,   3595117,   3744941,   3894707,   4044503,
      4493921,   4793501,   5093089,   5392679,   5692279,   5991883,
      6291469,   6591059,   6890641,   7190243,   7489829,   7789447,
      8089033,   8987807,   9586981,   10186177,  10785371,  11384539,
      11983729,  12582917,  13182109,  13781291,  14380469,  14979667,
      15578861,  16178053,  17895707,  19014187,  20132683,  21251141,
      22369661,  23488103,  24606583,  25725083,  26843549,  27962027,
      29080529,  30198989,  31317469,  32435981,  35791397,  38028379,
      40265327,  42502283,  44739259,  46976221,  49213237,  51450131,
      53687099,  55924061,  58161041,  60397993,  62634959,  64871921,
      71582857,  76056727,  80530643,  85004567,  89478503,  93952427,
      98426347,  102900263, 107374217, 111848111, 116322053, 120795971,
      125269877, 129743807, 143165587, 152113427, 161061283, 170009141,
      178956983, 187904819, 196852693, 205800547, 214748383, 223696237,
      232644089, 241591943, 250539763, 259487603, 268435399};

  const uint32_t num_primes = sizeof(primes) / sizeof(uint32_t);

  uint32_t hsize = primes[num_primes - 1];
  for (uint32_t i = 0; i < num_primes; ++i) {
    if (size <= primes[i]) {
      hsize = primes[i];
      break;
    }
  }
  return hsize;
}

}  // namespace Impl
}  // namespace Kokkos
