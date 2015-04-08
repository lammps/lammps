/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/



#include <Kokkos_Core.hpp>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#ifndef KOKKOS_RANDOM_HPP
#define KOKKOS_RANDOM_HPP

// These generators are based on Vigna, Sebastiano (2014). "An experimental exploration of Marsaglia's xorshift generators, scrambled"
// See: http://arxiv.org/abs/1402.6246

namespace Kokkos {

  /*Template functions to get equidistributed random numbers from a generator for a specific Scalar type

       template<class Generator,Scalar>
       struct rand{

         //Max value returned by draw(Generator& gen)
         KOKKOS_INLINE_FUNCTION
         static Scalar max();

         //Returns a value between zero and max()
         KOKKOS_INLINE_FUNCTION
         static Scalar draw(Generator& gen);

         //Returns a value between zero and range()
         //Note: for floating point values range can be larger than max()
         KOKKOS_INLINE_FUNCTION
         static Scalar draw(Generator& gen, const Scalar& range){}

         //Return value between start and end
         KOKKOS_INLINE_FUNCTION
         static Scalar draw(Generator& gen, const Scalar& start, const Scalar& end);
      };

    The Random number generators themselves have two components a state-pool and the actual generator
    A state-pool manages a number of generators, so that each active thread is able to grep its own.
    This allows the generation of random numbers which are independent between threads. Note that
    in contrast to CuRand none of the functions of the pool (or the generator) are collectives,
    i.e. all functions can be called inside conditionals.

    template<class Device>
    class Pool {
     public:
      //The Kokkos device type
      typedef Device device_type;
      //The actual generator type
      typedef Generator<Device> generator_type;

      //Default constructor: does not initialize a pool
      Pool();

      //Initializing constructor: calls init(seed,Device_Specific_Number);
      Pool(unsigned int seed);

      //Intialize Pool with seed as a starting seed with a pool_size of num_states
      //The Random_XorShift64 generator is used in serial to initialize all states,
      //thus the intialization process is platform independent and deterministic.
      void init(unsigned int seed, int num_states);

      //Get a generator. This will lock one of the states, guaranteeing that each thread
      //will have its private generator. Note: on Cuda getting a state involves atomics,
      //and is thus not deterministic!
      generator_type get_state();

      //Give a state back to the pool. This unlocks the state, and writes the modified
      //state of the generator back to the pool.
      void free_state(generator_type gen);

    }

    template<class Device>
    class Generator {
     public:
     //The Kokkos device type
    typedef DeviceType device_type;

    //Max return values of respective [X]rand[S]() functions
    enum {MAX_URAND = 0xffffffffU};
    enum {MAX_URAND64 = 0xffffffffffffffffULL-1};
    enum {MAX_RAND = static_cast<int>(0xffffffffU/2)};
    enum {MAX_RAND64 = static_cast<int64_t>(0xffffffffffffffffULL/2-1)};


    //Init with a state and the idx with respect to pool. Note: in serial the
    //Generator can be used by just giving it the necessary state arguments
    KOKKOS_INLINE_FUNCTION
    Generator (STATE_ARGUMENTS, int state_idx = 0);

    //Draw a equidistributed uint32_t in the range (0,MAX_URAND]
    KOKKOS_INLINE_FUNCTION
    uint32_t urand();

    //Draw a equidistributed uint64_t in the range (0,MAX_URAND64]
    KOKKOS_INLINE_FUNCTION
    uint64_t urand64();

    //Draw a equidistributed uint32_t in the range (0,range]
    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& range);

    //Draw a equidistributed uint32_t in the range (start,end]
    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& start, const uint32_t& end );

    //Draw a equidistributed uint64_t in the range (0,range]
    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& range);

    //Draw a equidistributed uint64_t in the range (start,end]
    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& start, const uint64_t& end );

    //Draw a equidistributed int in the range (0,MAX_RAND]
    KOKKOS_INLINE_FUNCTION
    int rand();

    //Draw a equidistributed int in the range (0,range]
    KOKKOS_INLINE_FUNCTION
    int rand(const int& range);

    //Draw a equidistributed int in the range (start,end]
    KOKKOS_INLINE_FUNCTION
    int rand(const int& start, const int& end );

    //Draw a equidistributed int64_t in the range (0,MAX_RAND64]
    KOKKOS_INLINE_FUNCTION
    int64_t rand64();

    //Draw a equidistributed int64_t in the range (0,range]
    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& range);

    //Draw a equidistributed int64_t in the range (start,end]
    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& start, const int64_t& end );

    //Draw a equidistributed float in the range (0,1.0]
    KOKKOS_INLINE_FUNCTION
    float frand();

    //Draw a equidistributed float in the range (0,range]
    KOKKOS_INLINE_FUNCTION
    float frand(const float& range);

    //Draw a equidistributed float in the range (start,end]
    KOKKOS_INLINE_FUNCTION
    float frand(const float& start, const float& end );

    //Draw a equidistributed double in the range (0,1.0]
    KOKKOS_INLINE_FUNCTION
    double drand();

    //Draw a equidistributed double in the range (0,range]
    KOKKOS_INLINE_FUNCTION
    double drand(const double& range);

    //Draw a equidistributed double in the range (start,end]
    KOKKOS_INLINE_FUNCTION
    double drand(const double& start, const double& end );

    //Draw a standard normal distributed double
    KOKKOS_INLINE_FUNCTION
    double normal() ;

    //Draw a normal distributed double with given mean and standard deviation
    KOKKOS_INLINE_FUNCTION
    double normal(const double& mean, const double& std_dev=1.0);
    }

    //Additional Functions:

    //Fills view with random numbers in the range (0,range]
    template<class ViewType, class PoolType>
    void fill_random(ViewType view, PoolType pool, ViewType::value_type range);

    //Fills view with random numbers in the range (start,end]
    template<class ViewType, class PoolType>
    void fill_random(ViewType view, PoolType pool,
                     ViewType::value_type start, ViewType::value_type end);

*/

  template<class Generator, class Scalar>
  struct rand;


  template<class Generator>
  struct rand<Generator,char> {

    KOKKOS_INLINE_FUNCTION
    static short max(){return 127;}
    KOKKOS_INLINE_FUNCTION
    static short draw(Generator& gen)
                          {return short((gen.rand()&0xff+256)%256);}
    KOKKOS_INLINE_FUNCTION
    static short draw(Generator& gen, const char& range)
                          {return char(gen.rand(range));}
    KOKKOS_INLINE_FUNCTION
    static short draw(Generator& gen, const char& start, const char& end)
                          {return char(gen.rand(start,end));}

  };

  template<class Generator>
  struct rand<Generator,short> {
    KOKKOS_INLINE_FUNCTION
    static short max(){return 32767;}
    KOKKOS_INLINE_FUNCTION
    static short draw(Generator& gen)
                          {return short((gen.rand()&0xffff+65536)%32768);}
    KOKKOS_INLINE_FUNCTION
    static short draw(Generator& gen, const short& range)
                          {return short(gen.rand(range));}
    KOKKOS_INLINE_FUNCTION
    static short draw(Generator& gen, const short& start, const short& end)
                          {return short(gen.rand(start,end));}

  };

  template<class Generator>
  struct rand<Generator,int> {
    KOKKOS_INLINE_FUNCTION
    static int max(){return Generator::MAX_RAND;}
    KOKKOS_INLINE_FUNCTION
    static int draw(Generator& gen)
                          {return gen.rand();}
    KOKKOS_INLINE_FUNCTION
    static int draw(Generator& gen, const int& range)
                          {return gen.rand(range);}
    KOKKOS_INLINE_FUNCTION
    static int draw(Generator& gen, const int& start, const int& end)
                          {return gen.rand(start,end);}

  };

  template<class Generator>
  struct rand<Generator,unsigned int> {
    KOKKOS_INLINE_FUNCTION
    static unsigned int max () {
      return Generator::MAX_URAND;
    }
    KOKKOS_INLINE_FUNCTION
    static unsigned int draw (Generator& gen) {
      return gen.urand ();
    }
    KOKKOS_INLINE_FUNCTION
    static unsigned int draw(Generator& gen, const unsigned int& range) {
      return gen.urand (range);
    }
    KOKKOS_INLINE_FUNCTION
    static unsigned int
    draw (Generator& gen, const unsigned int& start, const unsigned int& end) {
      return gen.urand (start, end);
    }
  };

  template<class Generator>
  struct rand<Generator,long> {
    KOKKOS_INLINE_FUNCTION
    static long max () {
      // FIXME (mfh 26 Oct 2014) It would be better to select the
      // return value at compile time, using something like enable_if.
      return sizeof (long) == 4 ?
        static_cast<long> (Generator::MAX_RAND) :
        static_cast<long> (Generator::MAX_RAND64);
    }
    KOKKOS_INLINE_FUNCTION
    static long draw (Generator& gen) {
      // FIXME (mfh 26 Oct 2014) It would be better to select the
      // return value at compile time, using something like enable_if.
      return sizeof (long) == 4 ?
        static_cast<long> (gen.rand ()) :
        static_cast<long> (gen.rand64 ());
    }
    KOKKOS_INLINE_FUNCTION
    static long draw (Generator& gen, const long& range) {
      // FIXME (mfh 26 Oct 2014) It would be better to select the
      // return value at compile time, using something like enable_if.
      return sizeof (long) == 4 ?
        static_cast<long> (gen.rand (static_cast<int> (range))) :
        static_cast<long> (gen.rand64 (range));
    }
    KOKKOS_INLINE_FUNCTION
    static long draw (Generator& gen, const long& start, const long& end) {
      // FIXME (mfh 26 Oct 2014) It would be better to select the
      // return value at compile time, using something like enable_if.
      return sizeof (long) == 4 ?
        static_cast<long> (gen.rand (static_cast<int> (start),
                                     static_cast<int> (end))) :
        static_cast<long> (gen.rand64 (start, end));
    }
  };

  template<class Generator>
  struct rand<Generator,unsigned long> {
    KOKKOS_INLINE_FUNCTION
    static unsigned long max () {
      // FIXME (mfh 26 Oct 2014) It would be better to select the
      // return value at compile time, using something like enable_if.
      return sizeof (unsigned long) == 4 ?
        static_cast<unsigned long> (Generator::MAX_URAND) :
        static_cast<unsigned long> (Generator::MAX_URAND64);
    }
    KOKKOS_INLINE_FUNCTION
    static unsigned long draw (Generator& gen) {
      // FIXME (mfh 26 Oct 2014) It would be better to select the
      // return value at compile time, using something like enable_if.
      return sizeof (unsigned long) == 4 ?
        static_cast<unsigned long> (gen.urand ()) :
        static_cast<unsigned long> (gen.urand64 ());
    }
    KOKKOS_INLINE_FUNCTION
    static unsigned long draw(Generator& gen, const unsigned long& range) {
      // FIXME (mfh 26 Oct 2014) It would be better to select the
      // return value at compile time, using something like enable_if.
      return sizeof (unsigned long) == 4 ?
        static_cast<unsigned long> (gen.urand (static_cast<unsigned int> (range))) :
        static_cast<unsigned long> (gen.urand64 (range));
    }
    KOKKOS_INLINE_FUNCTION
    static unsigned long
    draw (Generator& gen, const unsigned long& start, const unsigned long& end) {
      // FIXME (mfh 26 Oct 2014) It would be better to select the
      // return value at compile time, using something like enable_if.
      return sizeof (unsigned long) == 4 ?
        static_cast<unsigned long> (gen.urand (static_cast<unsigned int> (start),
                                               static_cast<unsigned int> (end))) :
        static_cast<unsigned long> (gen.urand64 (start, end));
    }
  };

  // NOTE (mfh 26 oct 2014) This is a partial specialization for long
  // long, a C99 / C++11 signed type which is guaranteed to be at
  // least 64 bits.  Do NOT write a partial specialization for
  // int64_t!!!  This is just a typedef!  It could be either long or
  // long long.  We don't know which a priori, and I've seen both.
  // The types long and long long are guaranteed to differ, so it's
  // always safe to specialize for both.
  template<class Generator>
  struct rand<Generator, long long> {
    KOKKOS_INLINE_FUNCTION
    static long long max () {
      // FIXME (mfh 26 Oct 2014) It's legal for long long to be > 64 bits.
      return Generator::MAX_RAND64;
    }
    KOKKOS_INLINE_FUNCTION
    static long long draw (Generator& gen) {
      // FIXME (mfh 26 Oct 2014) It's legal for long long to be > 64 bits.
      return gen.rand64 ();
    }
    KOKKOS_INLINE_FUNCTION
    static long long draw (Generator& gen, const long long& range) {
      // FIXME (mfh 26 Oct 2014) It's legal for long long to be > 64 bits.
      return gen.rand64 (range);
    }
    KOKKOS_INLINE_FUNCTION
    static long long draw (Generator& gen, const long long& start, const long long& end) {
      // FIXME (mfh 26 Oct 2014) It's legal for long long to be > 64 bits.
      return gen.rand64 (start, end);
    }
  };

  // NOTE (mfh 26 oct 2014) This is a partial specialization for
  // unsigned long long, a C99 / C++11 unsigned type which is
  // guaranteed to be at least 64 bits.  Do NOT write a partial
  // specialization for uint64_t!!!  This is just a typedef!  It could
  // be either unsigned long or unsigned long long.  We don't know
  // which a priori, and I've seen both.  The types unsigned long and
  // unsigned long long are guaranteed to differ, so it's always safe
  // to specialize for both.
  template<class Generator>
  struct rand<Generator,unsigned long long> {
    KOKKOS_INLINE_FUNCTION
    static unsigned long long max () {
      // FIXME (mfh 26 Oct 2014) It's legal for unsigned long long to be > 64 bits.
      return Generator::MAX_URAND64;
    }
    KOKKOS_INLINE_FUNCTION
    static unsigned long long draw (Generator& gen) {
      // FIXME (mfh 26 Oct 2014) It's legal for unsigned long long to be > 64 bits.
      return gen.urand64 ();
    }
    KOKKOS_INLINE_FUNCTION
    static unsigned long long draw (Generator& gen, const unsigned long long& range) {
      // FIXME (mfh 26 Oct 2014) It's legal for long long to be > 64 bits.
      return gen.urand64 (range);
    }
    KOKKOS_INLINE_FUNCTION
    static unsigned long long
    draw (Generator& gen, const unsigned long long& start, const unsigned long long& end) {
      // FIXME (mfh 26 Oct 2014) It's legal for long long to be > 64 bits.
      return gen.urand64 (start, end);
    }
  };

  template<class Generator>
  struct rand<Generator,float> {
    KOKKOS_INLINE_FUNCTION
    static float max(){return 1.0f;}
    KOKKOS_INLINE_FUNCTION
    static float draw(Generator& gen)
                          {return gen.frand();}
    KOKKOS_INLINE_FUNCTION
    static float draw(Generator& gen, const float& range)
                          {return gen.frand(range);}
    KOKKOS_INLINE_FUNCTION
    static float draw(Generator& gen, const float& start, const float& end)
                          {return gen.frand(start,end);}

  };

  template<class Generator>
  struct rand<Generator,double> {
    KOKKOS_INLINE_FUNCTION
    static double max(){return 1.0;}
    KOKKOS_INLINE_FUNCTION
    static double draw(Generator& gen)
                          {return gen.drand();}
    KOKKOS_INLINE_FUNCTION
    static double draw(Generator& gen, const double& range)
                          {return gen.drand(range);}
    KOKKOS_INLINE_FUNCTION
    static double draw(Generator& gen, const double& start, const double& end)
                          {return gen.drand(start,end);}

  };

  template<class DeviceType>
  class Random_XorShift64_Pool;

  template<class DeviceType>
  class Random_XorShift64 {
  private:
    uint64_t state_;
    const int state_idx_;
    friend class Random_XorShift64_Pool<DeviceType>;
  public:

    typedef DeviceType device_type;

    enum {MAX_URAND = 0xffffffffU};
    enum {MAX_URAND64 = 0xffffffffffffffffULL-1};
    enum {MAX_RAND = static_cast<int>(0xffffffff/2)};
    enum {MAX_RAND64 = static_cast<int64_t>(0xffffffffffffffffLL/2-1)};

    KOKKOS_INLINE_FUNCTION
    Random_XorShift64 (uint64_t state, int state_idx = 0)
     : state_(state),state_idx_(state_idx){}

    KOKKOS_INLINE_FUNCTION
    uint32_t urand() {
      state_ ^= state_ >> 12;
      state_ ^= state_ << 25;
      state_ ^= state_ >> 27;

      uint64_t tmp = state_ * 2685821657736338717ULL;
      tmp = tmp>>16;
      return static_cast<uint32_t>(tmp&MAX_URAND);
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64() {
      state_ ^= state_ >> 12;
      state_ ^= state_ << 25;
      state_ ^= state_ >> 27;
      return (state_ * 2685821657736338717ULL) - 1;
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& range) {
      const uint32_t max_val = (MAX_URAND/range)*range;
      uint32_t tmp = urand();
      while(tmp>=max_val)
        tmp = urand();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& start, const uint32_t& end ) {
      return urand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& range) {
      const uint64_t max_val = (MAX_URAND64/range)*range;
      uint64_t tmp = urand64();
      while(tmp>=max_val)
        tmp = urand64();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& start, const uint64_t& end ) {
      return urand64(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    int rand() {
      return static_cast<int>(urand()/2);
    }

    KOKKOS_INLINE_FUNCTION
    int rand(const int& range) {
      const int max_val = (MAX_RAND/range)*range;
      int tmp = rand();
      while(tmp>=max_val)
        tmp = rand();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    int rand(const int& start, const int& end ) {
      return rand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64() {
      return static_cast<int64_t>(urand64()/2);
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& range) {
      const int64_t max_val = (MAX_RAND64/range)*range;
      int64_t tmp = rand64();
      while(tmp>=max_val)
        tmp = rand64();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& start, const int64_t& end ) {
      return rand64(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    float frand() {
      return 1.0f * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    float frand(const float& range) {
      return range * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    float frand(const float& start, const float& end ) {
      return frand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    double drand() {
      return 1.0 * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    double drand(const double& range) {
      return range * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    double drand(const double& start, const double& end ) {
      return drand(end-start)+start;
    }

    //Marsaglia polar method for drawing a standard normal distributed random number
    KOKKOS_INLINE_FUNCTION
    double normal() {
      double S = 2.0;
      double U;
      while(S>=1.0) {
        U = drand();
        const double V = drand();
        S = U*U+V*V;
      }
      return U*sqrt(-2.0*log(S)/S);
    }

    KOKKOS_INLINE_FUNCTION
    double normal(const double& mean, const double& std_dev=1.0) {
      return mean + normal()*std_dev;
    }

  };

  template<class DeviceType = Kokkos::DefaultExecutionSpace>
  class Random_XorShift64_Pool {
  private:
    typedef View<int*,DeviceType> lock_type;
    typedef View<uint64_t*,DeviceType> state_data_type;
    lock_type locks_;
    state_data_type state_;
    int num_states_;

  public:
    typedef Random_XorShift64<DeviceType> generator_type;
    typedef DeviceType device_type;

    Random_XorShift64_Pool() {
      num_states_ = 0;
    }
    Random_XorShift64_Pool(unsigned int seed) {
      num_states_ = 0;
      init(seed,DeviceType::max_hardware_threads());
    }

    Random_XorShift64_Pool(const Random_XorShift64_Pool& src):
      locks_(src.locks_),
      state_(src.state_),
      num_states_(src.num_states_)
    {}

    Random_XorShift64_Pool operator = (const Random_XorShift64_Pool& src) {
      locks_ = src.locks_;
      state_ = src.state_;
      num_states_ = src.num_states_;
      return *this;
    }

    void init(unsigned int seed, int num_states) {
      num_states_ = num_states;

      locks_ = lock_type("Kokkos::Random_XorShift64::locks",num_states_);
      state_ = state_data_type("Kokkos::Random_XorShift64::state",num_states_);

      typename state_data_type::HostMirror h_state = create_mirror_view(state_);
      typename lock_type::HostMirror h_lock = create_mirror_view(locks_);

      // Execute on the HostMirror's default execution space.
      Random_XorShift64<typename state_data_type::HostMirror::execution_space> gen(seed,0);
      for(int i = 0; i < 17; i++)
        gen.rand();
      for(int i = 0; i < num_states_; i++) {
        int n1 = gen.rand();
        int n2 = gen.rand();
        int n3 = gen.rand();
        int n4 = gen.rand();
        h_state(i) = (((static_cast<uint64_t>(n1)) & 0xffff)<<00) |
                     (((static_cast<uint64_t>(n2)) & 0xffff)<<16) |
                     (((static_cast<uint64_t>(n3)) & 0xffff)<<32) |
                     (((static_cast<uint64_t>(n4)) & 0xffff)<<48);
        h_lock(i) = 0;
      }
      deep_copy(state_,h_state);
      deep_copy(locks_,h_lock);
    }

    KOKKOS_INLINE_FUNCTION
    Random_XorShift64<DeviceType> get_state() const {
      const int i = DeviceType::hardware_thread_id();;
      return Random_XorShift64<DeviceType>(state_(i),i);
    }

    KOKKOS_INLINE_FUNCTION
    void free_state(const Random_XorShift64<DeviceType>& state) const {
      state_(state.state_idx_) = state.state_;
    }
  };


  template<class DeviceType>
  class Random_XorShift1024_Pool;

  template<class DeviceType>
  class Random_XorShift1024 {
  private:
    int p_;
    const int state_idx_;
    uint64_t state_[16];
    friend class Random_XorShift1024_Pool<DeviceType>;
  public:

    typedef DeviceType device_type;

    enum {MAX_URAND = 0xffffffffU};
    enum {MAX_URAND64 = 0xffffffffffffffffULL-1};
    enum {MAX_RAND = static_cast<int>(0xffffffffU/2)};
    enum {MAX_RAND64 = static_cast<int64_t>(0xffffffffffffffffULL/2-1)};

    KOKKOS_INLINE_FUNCTION
    Random_XorShift1024 (uint64_t* state, int p, int state_idx = 0):
      p_(p),state_idx_(state_idx){
      for(int i=0 ; i<16; i++)
        state_[i] = state[i];
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand() {
      uint64_t state_0 = state_[ p_ ];
      uint64_t state_1 = state_[ p_ = ( p_ + 1 ) & 15 ];
      state_1 ^= state_1 << 31;
      state_1 ^= state_1 >> 11;
      state_0 ^= state_0 >> 30;
      uint64_t tmp = ( state_[ p_ ] = state_0 ^ state_1 ) * 1181783497276652981ULL;
      tmp = tmp>>16;
      return static_cast<uint32_t>(tmp&MAX_URAND);
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64() {
      uint64_t state_0 = state_[ p_ ];
      uint64_t state_1 = state_[ p_ = ( p_ + 1 ) & 15 ];
      state_1 ^= state_1 << 31;
      state_1 ^= state_1 >> 11;
      state_0 ^= state_0 >> 30;
      return (( state_[ p_ ] = state_0 ^ state_1 ) * 1181783497276652981LL) - 1;
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& range) {
      const uint32_t max_val = (MAX_URAND/range)*range;
      uint32_t tmp = urand();
      while(tmp>=max_val)
        tmp = urand();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& start, const uint32_t& end ) {
      return urand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& range) {
      const uint64_t max_val = (MAX_URAND64/range)*range;
      uint64_t tmp = urand64();
      while(tmp>=max_val)
        tmp = urand64();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& start, const uint64_t& end ) {
      return urand64(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    int rand() {
      return static_cast<int>(urand()/2);
    }

    KOKKOS_INLINE_FUNCTION
    int rand(const int& range) {
      const int max_val = (MAX_RAND/range)*range;
      int tmp = rand();
      while(tmp>=max_val)
        tmp = rand();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    int rand(const int& start, const int& end ) {
      return rand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64() {
      return static_cast<int64_t>(urand64()/2);
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& range) {
      const int64_t max_val = (MAX_RAND64/range)*range;
      int64_t tmp = rand64();
      while(tmp>=max_val)
        tmp = rand64();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& start, const int64_t& end ) {
      return rand64(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    float frand() {
      return 1.0f * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    float frand(const float& range) {
      return range * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    float frand(const float& start, const float& end ) {
      return frand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    double drand() {
      return 1.0 * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    double drand(const double& range) {
      return range * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    double drand(const double& start, const double& end ) {
      return frand(end-start)+start;
    }

    //Marsaglia polar method for drawing a standard normal distributed random number
    KOKKOS_INLINE_FUNCTION
    double normal() {
      double S = 2.0;
      double U;
      while(S>=1.0) {
        U = drand();
        const double V = drand();
        S = U*U+V*V;
      }
      return U*sqrt(-2.0*log(S)/S);
    }

    KOKKOS_INLINE_FUNCTION
    double normal(const double& mean, const double& std_dev=1.0) {
      return mean + normal()*std_dev;
    }
  };


  template<class DeviceType = Kokkos::DefaultExecutionSpace>
  class Random_XorShift1024_Pool {
  private:
    typedef View<int*,DeviceType> int_view_type;
    typedef View<uint64_t*[16],DeviceType> state_data_type;

    int_view_type locks_;
    state_data_type state_;
    int_view_type p_;
    int num_states_;

  public:
    typedef Random_XorShift1024<DeviceType> generator_type;

    typedef DeviceType device_type;

    Random_XorShift1024_Pool() {
      num_states_ = 0;
    }

    inline
    Random_XorShift1024_Pool(unsigned int seed){
      num_states_ = 0;
      init(seed,DeviceType::max_hardware_threads());
    }

    Random_XorShift1024_Pool(const Random_XorShift1024_Pool& src):
      locks_(src.locks_),
      state_(src.state_),
      p_(src.p_),
      num_states_(src.num_states_)
    {}

    Random_XorShift1024_Pool operator = (const Random_XorShift1024_Pool& src) {
      locks_ = src.locks_;
      state_ = src.state_;
      p_ = src.p_;
      num_states_ = src.num_states_;
      return *this;
    }

    inline
    void init(unsigned int seed, int num_states) {
      num_states_ = num_states;

      locks_ = int_view_type("Kokkos::Random_XorShift1024::locks",num_states_);
      state_ = state_data_type("Kokkos::Random_XorShift1024::state",num_states_);
      p_ = int_view_type("Kokkos::Random_XorShift1024::p",num_states_);

      typename state_data_type::HostMirror h_state = create_mirror_view(state_);
      typename int_view_type::HostMirror h_lock = create_mirror_view(locks_);
      typename int_view_type::HostMirror h_p = create_mirror_view(p_);

      // Execute on the HostMirror's default execution space.
      Random_XorShift64<typename state_data_type::HostMirror::execution_space> gen(seed,0);
      for(int i = 0; i < 17; i++)
        gen.rand();
      for(int i = 0; i < num_states_; i++) {
        for(int j = 0; j < 16 ; j++) {
          int n1 = gen.rand();
          int n2 = gen.rand();
          int n3 = gen.rand();
          int n4 = gen.rand();
          h_state(i,j) = (((static_cast<uint64_t>(n1)) & 0xffff)<<00) |
                         (((static_cast<uint64_t>(n2)) & 0xffff)<<16) |
                         (((static_cast<uint64_t>(n3)) & 0xffff)<<32) |
                         (((static_cast<uint64_t>(n4)) & 0xffff)<<48);
        }
        h_p(i) = 0;
        h_lock(i) = 0;
      }
      deep_copy(state_,h_state);
      deep_copy(locks_,h_lock);
    }

    KOKKOS_INLINE_FUNCTION
    Random_XorShift1024<DeviceType> get_state() const {
      const int i = DeviceType::hardware_thread_id();
      return Random_XorShift1024<DeviceType>(&state_(i,0),p_(i),i);
    };

    KOKKOS_INLINE_FUNCTION
    void free_state(const Random_XorShift1024<DeviceType>& state) const {
      for(int i = 0; i<16; i++)
        state_(state.state_idx_,i) = state.state_[i];
      p_(state.state_idx_) = state.p_;
    }
  };

#if defined(KOKKOS_HAVE_CUDA) && defined(__CUDACC__)

  template<>
  class Random_XorShift1024<Kokkos::Cuda> {
  private:
    int p_;
    const int state_idx_;
    uint64_t* state_;
    friend class Random_XorShift1024_Pool<Kokkos::Cuda>;
  public:

    typedef Kokkos::Cuda device_type;

    enum {MAX_URAND = 0xffffffffU};
    enum {MAX_URAND64 = 0xffffffffffffffffULL-1};
    enum {MAX_RAND = static_cast<int>(0xffffffffU/2)};
    enum {MAX_RAND64 = static_cast<int64_t>(0xffffffffffffffffULL/2-1)};

    KOKKOS_INLINE_FUNCTION
    Random_XorShift1024 (uint64_t* state, int p, int state_idx = 0):
      p_(p),state_idx_(state_idx),state_(state){
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand() {
      uint64_t state_0 = state_[ p_ ];
      uint64_t state_1 = state_[ p_ = ( p_ + 1 ) & 15 ];
      state_1 ^= state_1 << 31;
      state_1 ^= state_1 >> 11;
      state_0 ^= state_0 >> 30;
      uint64_t tmp = ( state_[ p_ ] = state_0 ^ state_1 ) * 1181783497276652981ULL;
      tmp = tmp>>16;
      return static_cast<uint32_t>(tmp&MAX_URAND);
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64() {
      uint64_t state_0 = state_[ p_ ];
      uint64_t state_1 = state_[ p_ = ( p_ + 1 ) & 15 ];
      state_1 ^= state_1 << 31;
      state_1 ^= state_1 >> 11;
      state_0 ^= state_0 >> 30;
      return (( state_[ p_ ] = state_0 ^ state_1 ) * 1181783497276652981LL) - 1;
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& range) {
      const uint32_t max_val = (MAX_URAND/range)*range;
      uint32_t tmp = urand();
      while(tmp>=max_val)
        urand();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& start, const uint32_t& end ) {
      return urand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& range) {
      const uint64_t max_val = (MAX_URAND64/range)*range;
      uint64_t tmp = urand64();
      while(tmp>=max_val)
        urand64();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& start, const uint64_t& end ) {
      return urand64(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    int rand() {
      return static_cast<int>(urand()/2);
    }

    KOKKOS_INLINE_FUNCTION
    int rand(const int& range) {
      const int max_val = (MAX_RAND/range)*range;
      int tmp = rand();
      while(tmp>=max_val)
        rand();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    int rand(const int& start, const int& end ) {
      return rand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64() {
      return static_cast<int64_t>(urand64()/2);
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& range) {
      const int64_t max_val = (MAX_RAND64/range)*range;
      int64_t tmp = rand64();
      while(tmp>=max_val)
        rand64();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& start, const int64_t& end ) {
      return rand64(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    float frand() {
      return 1.0f * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    float frand(const float& range) {
      return range * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    float frand(const float& start, const float& end ) {
      return frand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    double drand() {
      return 1.0 * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    double drand(const double& range) {
      return range * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    double drand(const double& start, const double& end ) {
      return frand(end-start)+start;
    }

    //Marsaglia polar method for drawing a standard normal distributed random number
    KOKKOS_INLINE_FUNCTION
    double normal() {
      double S = 2.0;
      double U;
      while(S>=1.0) {
        U = drand();
        const double V = drand();
        S = U*U+V*V;
      }
      return U*sqrt(-2.0*log(S)/S);
    }

    KOKKOS_INLINE_FUNCTION
    double normal(const double& mean, const double& std_dev=1.0) {
      return mean + normal()*std_dev;
    }
  };

template<>
inline
Random_XorShift64_Pool<Kokkos::Cuda>::Random_XorShift64_Pool(unsigned int seed) {
  num_states_ = 0;
  init(seed,4*32768);
}

template<>
KOKKOS_INLINE_FUNCTION
Random_XorShift64<Kokkos::Cuda> Random_XorShift64_Pool<Kokkos::Cuda>::get_state() const {
#ifdef __CUDA_ARCH__
  const int i_offset = (threadIdx.x*blockDim.y + threadIdx.y)*blockDim.z+threadIdx.z;
  int i = ((blockIdx.x*gridDim.y+blockIdx.y)*gridDim.z + blockIdx.z) *
           blockDim.x*blockDim.y*blockDim.z + i_offset;
  while(Kokkos::atomic_compare_exchange(&locks_(i),0,1)) {
      i+=blockDim.x*blockDim.y*blockDim.z;
      if(i>=num_states_) {i = i_offset;}
  }

  return Random_XorShift64<Kokkos::Cuda>(state_(i),i);
#else
  return Random_XorShift64<Kokkos::Cuda>(state_(0),0);
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
void Random_XorShift64_Pool<Kokkos::Cuda>::free_state(const Random_XorShift64<Kokkos::Cuda> &state) const {
#ifdef __CUDA_ARCH__
  state_(state.state_idx_) = state.state_;
  locks_(state.state_idx_) = 0;
  return;
#endif
}


template<>
inline
Random_XorShift1024_Pool<Kokkos::Cuda>::Random_XorShift1024_Pool(unsigned int seed) {
  num_states_ = 0;
  init(seed,4*32768);
}

template<>
KOKKOS_INLINE_FUNCTION
Random_XorShift1024<Kokkos::Cuda> Random_XorShift1024_Pool<Kokkos::Cuda>::get_state() const {
#ifdef __CUDA_ARCH__
  const int i_offset = (threadIdx.x*blockDim.y + threadIdx.y)*blockDim.z+threadIdx.z;
  int i = ((blockIdx.x*gridDim.y+blockIdx.y)*gridDim.z + blockIdx.z) *
           blockDim.x*blockDim.y*blockDim.z + i_offset;
  while(Kokkos::atomic_compare_exchange(&locks_(i),0,1)) {
      i+=blockDim.x*blockDim.y*blockDim.z;
      if(i>=num_states_) {i = i_offset;}
  }

  return Random_XorShift1024<Kokkos::Cuda>(&state_(i,0), p_(i), i);
#else
  return Random_XorShift1024<Kokkos::Cuda>(&state_(0,0), p_(0), 0);
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
void Random_XorShift1024_Pool<Kokkos::Cuda>::free_state(const Random_XorShift1024<Kokkos::Cuda> &state) const {
#ifdef __CUDA_ARCH__
  for(int i=0; i<16; i++)
    state_(state.state_idx_,i) = state.state_[i];
  locks_(state.state_idx_) = 0;
  return;
#endif
}


#endif



template<class ViewType, class RandomPool, int loops, int rank>
struct fill_random_functor_range;
template<class ViewType, class RandomPool, int loops, int rank>
struct fill_random_functor_begin_end;

template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_range<ViewType,RandomPool,loops,1>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type range;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_range(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type range_):
    a(a_),rand_pool(rand_pool_),range(range_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0())
        a(idx) = Rand::draw(gen,range);
    }
    rand_pool.free_state(gen);
  }
};

template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_range<ViewType,RandomPool,loops,2>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type range;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_range(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type range_):
    a(a_),rand_pool(rand_pool_),range(range_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()) {
        for(unsigned int k=0;k<a.dimension_1();k++)
          a(idx,k) = Rand::draw(gen,range);
      }
    }
    rand_pool.free_state(gen);
  }
};


template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_range<ViewType,RandomPool,loops,3>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type range;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_range(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type range_):
    a(a_),rand_pool(rand_pool_),range(range_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()) {
        for(unsigned int k=0;k<a.dimension_1();k++)
          for(unsigned int l=0;l<a.dimension_2();l++)
            a(idx,k,l) = Rand::draw(gen,range);
      }
    }
    rand_pool.free_state(gen);
  }
};

template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_range<ViewType,RandomPool,loops,4>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type range;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_range(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type range_):
    a(a_),rand_pool(rand_pool_),range(range_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()) {
        for(unsigned int k=0;k<a.dimension_1();k++)
          for(unsigned int l=0;l<a.dimension_2();l++)
            for(unsigned int m=0;m<a.dimension_3();m++)
              a(idx,k,l,m) = Rand::draw(gen,range);
      }
    }
    rand_pool.free_state(gen);
  }
};

template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_range<ViewType,RandomPool,loops,5>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type range;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_range(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type range_):
    a(a_),rand_pool(rand_pool_),range(range_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()) {
        for(unsigned int k=0;k<a.dimension_1();k++)
          for(unsigned int l=0;l<a.dimension_2();l++)
            for(unsigned int m=0;m<a.dimension_3();m++)
              for(unsigned int n=0;n<a.dimension_4();n++)
              a(idx,k,l,m,n) = Rand::draw(gen,range);
      }
    }
    rand_pool.free_state(gen);
  }
};

template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_range<ViewType,RandomPool,loops,6>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type range;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_range(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type range_):
    a(a_),rand_pool(rand_pool_),range(range_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()) {
        for(unsigned int k=0;k<a.dimension_1();k++)
          for(unsigned int l=0;l<a.dimension_2();l++)
            for(unsigned int m=0;m<a.dimension_3();m++)
              for(unsigned int n=0;n<a.dimension_4();n++)
                for(unsigned int o=0;o<a.dimension_5();o++)
              a(idx,k,l,m,n,o) = Rand::draw(gen,range);
      }
    }
    rand_pool.free_state(gen);
  }
};

template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_range<ViewType,RandomPool,loops,7>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type range;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_range(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type range_):
    a(a_),rand_pool(rand_pool_),range(range_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()) {
        for(unsigned int k=0;k<a.dimension_1();k++)
          for(unsigned int l=0;l<a.dimension_2();l++)
            for(unsigned int m=0;m<a.dimension_3();m++)
              for(unsigned int n=0;n<a.dimension_4();n++)
                for(unsigned int o=0;o<a.dimension_5();o++)
                  for(unsigned int p=0;p<a.dimension_6();p++)
              a(idx,k,l,m,n,o,p) = Rand::draw(gen,range);
      }
    }
    rand_pool.free_state(gen);
  }
};

template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_range<ViewType,RandomPool,loops,8>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type range;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_range(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type range_):
    a(a_),rand_pool(rand_pool_),range(range_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()) {
        for(unsigned int k=0;k<a.dimension_1();k++)
          for(unsigned int l=0;l<a.dimension_2();l++)
            for(unsigned int m=0;m<a.dimension_3();m++)
              for(unsigned int n=0;n<a.dimension_4();n++)
                for(unsigned int o=0;o<a.dimension_5();o++)
                  for(unsigned int p=0;p<a.dimension_6();p++)
                    for(unsigned int q=0;q<a.dimension_7();q++)
              a(idx,k,l,m,n,o,p,q) = Rand::draw(gen,range);
      }
    }
    rand_pool.free_state(gen);
  }
};
template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_begin_end<ViewType,RandomPool,loops,1>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type begin,end;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_begin_end(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type begin_, typename ViewType::const_value_type end_):
    a(a_),rand_pool(rand_pool_),begin(begin_),end(end_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0())
        a(idx) = Rand::draw(gen,begin,end);
    }
    rand_pool.free_state(gen);
  }
};

template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_begin_end<ViewType,RandomPool,loops,2>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type begin,end;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_begin_end(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type begin_, typename ViewType::const_value_type end_):
    a(a_),rand_pool(rand_pool_),begin(begin_),end(end_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()) {
        for(unsigned int k=0;k<a.dimension_1();k++)
          a(idx,k) = Rand::draw(gen,begin,end);
      }
    }
    rand_pool.free_state(gen);
  }
};


template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_begin_end<ViewType,RandomPool,loops,3>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type begin,end;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_begin_end(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type begin_, typename ViewType::const_value_type end_):
    a(a_),rand_pool(rand_pool_),begin(begin_),end(end_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()) {
        for(unsigned int k=0;k<a.dimension_1();k++)
          for(unsigned int l=0;l<a.dimension_2();l++)
            a(idx,k,l) = Rand::draw(gen,begin,end);
      }
    }
    rand_pool.free_state(gen);
  }
};

template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_begin_end<ViewType,RandomPool,loops,4>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type begin,end;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_begin_end(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type begin_, typename ViewType::const_value_type end_):
    a(a_),rand_pool(rand_pool_),begin(begin_),end(end_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()) {
        for(unsigned int k=0;k<a.dimension_1();k++)
          for(unsigned int l=0;l<a.dimension_2();l++)
            for(unsigned int m=0;m<a.dimension_3();m++)
              a(idx,k,l,m) = Rand::draw(gen,begin,end);
      }
    }
    rand_pool.free_state(gen);
  }
};

template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_begin_end<ViewType,RandomPool,loops,5>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type begin,end;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_begin_end(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type begin_, typename ViewType::const_value_type end_):
    a(a_),rand_pool(rand_pool_),begin(begin_),end(end_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()){
        for(unsigned int l=0;l<a.dimension_1();l++)
          for(unsigned int m=0;m<a.dimension_2();m++)
            for(unsigned int n=0;n<a.dimension_3();n++)
              for(unsigned int o=0;o<a.dimension_4();o++)
          a(idx,l,m,n,o) = Rand::draw(gen,begin,end);
      }
    }
    rand_pool.free_state(gen);
  }
};

template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_begin_end<ViewType,RandomPool,loops,6>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type begin,end;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_begin_end(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type begin_, typename ViewType::const_value_type end_):
    a(a_),rand_pool(rand_pool_),begin(begin_),end(end_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()) {
        for(unsigned int k=0;k<a.dimension_1();k++)
          for(unsigned int l=0;l<a.dimension_2();l++)
            for(unsigned int m=0;m<a.dimension_3();m++)
              for(unsigned int n=0;n<a.dimension_4();n++)
                for(unsigned int o=0;o<a.dimension_5();o++)
          a(idx,k,l,m,n,o) = Rand::draw(gen,begin,end);
      }
    }
    rand_pool.free_state(gen);
  }
};


template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_begin_end<ViewType,RandomPool,loops,7>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type begin,end;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_begin_end(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type begin_, typename ViewType::const_value_type end_):
    a(a_),rand_pool(rand_pool_),begin(begin_),end(end_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()) {
        for(unsigned int k=0;k<a.dimension_1();k++)
          for(unsigned int l=0;l<a.dimension_2();l++)
            for(unsigned int m=0;m<a.dimension_3();m++)
              for(unsigned int n=0;n<a.dimension_4();n++)
                for(unsigned int o=0;o<a.dimension_5();o++)
                  for(unsigned int p=0;p<a.dimension_6();p++)
            a(idx,k,l,m,n,o,p) = Rand::draw(gen,begin,end);
      }
    }
    rand_pool.free_state(gen);
  }
};

template<class ViewType, class RandomPool, int loops>
struct fill_random_functor_begin_end<ViewType,RandomPool,loops,8>{
  typedef typename ViewType::device_type device_type;
  ViewType a;
  RandomPool rand_pool;
  typename ViewType::const_value_type begin,end;

  typedef rand<typename RandomPool::generator_type, typename ViewType::non_const_value_type> Rand;

  fill_random_functor_begin_end(ViewType a_, RandomPool rand_pool_,
      typename ViewType::const_value_type begin_, typename ViewType::const_value_type end_):
    a(a_),rand_pool(rand_pool_),begin(begin_),end(end_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned int i) const {
    typename RandomPool::generator_type gen = rand_pool.get_state();
    for(unsigned int j=0;j<loops;j++) {
      const uint64_t idx = i*loops+j;
      if(idx<a.dimension_0()) {
        for(unsigned int k=0;k<a.dimension_1();k++)
          for(unsigned int l=0;l<a.dimension_2();l++)
            for(unsigned int m=0;m<a.dimension_3();m++)
              for(unsigned int n=0;n<a.dimension_4();n++)
                for(unsigned int o=0;o<a.dimension_5();o++)
                  for(unsigned int p=0;p<a.dimension_6();p++)
                    for(unsigned int q=0;q<a.dimension_7();q++)
              a(idx,k,l,m,n,o,p,q) = Rand::draw(gen,begin,end);
      }
    }
    rand_pool.free_state(gen);
  }
};

template<class ViewType, class RandomPool>
void fill_random(ViewType a, RandomPool g, typename ViewType::const_value_type range) {
  int64_t LDA = a.dimension_0();
  if(LDA>0)
    parallel_for((LDA+127)/128,fill_random_functor_range<ViewType,RandomPool,128,ViewType::Rank>(a,g,range));
}

template<class ViewType, class RandomPool>
void fill_random(ViewType a, RandomPool g, typename ViewType::const_value_type begin,typename ViewType::const_value_type end ) {
  int64_t LDA = a.dimension_0();
  if(LDA>0)
    parallel_for((LDA+127)/128,fill_random_functor_begin_end<ViewType,RandomPool,128,ViewType::Rank>(a,g,begin,end));
}
}

#endif
