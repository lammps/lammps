/*

Copyright (c) 2014, NVIDIA Corporation
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef KOKKOS_SYNCHRONIC_N3998_HPP
#define KOKKOS_SYNCHRONIC_N3998_HPP

#include <impl/Kokkos_Synchronic.hpp>
#include <functional>

/*
In the section below, a synchronization point represents a point at which a
thread may block until a given synchronization condition has been reached or
at which it may notify other threads that a synchronization condition has
been achieved.
*/
namespace Kokkos { namespace Impl {

    /*
    A latch maintains an internal counter that is initialized when the latch
    is created. The synchronization condition is reached when the counter is
    decremented to 0. Threads may block at a synchronization point waiting
    for the condition to be reached. When the condition is reached, any such
    blocked threads will be released.
    */
    struct latch {
        latch(int val) : count(val), released(false) { }
        latch(const latch&) = delete;
        latch& operator=(const latch&) = delete;
        ~latch( ) { }
        void arrive( ) {
            __arrive( );
        }
        void arrive_and_wait( ) {
            if(!__arrive( ))
                wait( );
        }
        void wait( ) {
            while(!released.load_when_not_equal(false,std::memory_order_acquire))
                ;
        }
        bool try_wait( ) {
            return released.load(std::memory_order_acquire);
        }
    private:
        bool __arrive( ) {
            if(count.fetch_add(-1,std::memory_order_release)!=1)
                return false;
            released.store(true,std::memory_order_release);
            return true;
        }
        std::atomic<int> count;
        synchronic<bool> released;
    };

    /*
    A barrier is created with an initial value representing the number of threads
    that can arrive at the synchronization point. When that many threads have
    arrived, the  synchronization condition is reached and the threads are
    released. The barrier will then reset, and may be reused for a new cycle, in
    which the same set of threads may arrive again at the synchronization point.
    The same set of threads shall arrive at the barrier in each cycle, otherwise
    the behaviour is undefined.
    */
    struct barrier {
        barrier(int val) : expected(val), arrived(0), nexpected(val), epoch(0) { }
        barrier(const barrier&) = delete;
        barrier& operator=(const barrier&) = delete;
        ~barrier() { }
        void arrive_and_wait() {
            int const myepoch = epoch.load(std::memory_order_relaxed);
            if(!__arrive(myepoch))
                while(epoch.load_when_not_equal(myepoch,std::memory_order_acquire) == myepoch)
                    ;
        }
        void arrive_and_drop() {
            nexpected.fetch_add(-1,std::memory_order_relaxed);
            __arrive(epoch.load(std::memory_order_relaxed));
        }
    private:
        bool __arrive(int const myepoch) {
            int const myresult = arrived.fetch_add(1,std::memory_order_acq_rel) + 1;
            if(__builtin_expect(myresult == expected,0)) {
                expected = nexpected.load(std::memory_order_relaxed);
                arrived.store(0,std::memory_order_relaxed);
                epoch.store(myepoch+1,std::memory_order_release);
                return true;
            }
            return false;
        }
        int expected;
        std::atomic<int> arrived, nexpected;
        synchronic<int> epoch;
    };

    /*
    A notifying barrier behaves as a barrier, but is constructed with a callable
    completion function that is invoked after all threads have arrived at the
    synchronization point, and before the synchronization condition is reached.
    The completion may modify the set of threads that arrives at the barrier in
    each cycle.
    */
    struct notifying_barrier {
        template <typename T>
        notifying_barrier(int val, T && f) : expected(val), arrived(0), nexpected(val), epoch(0), completion(std::forward<T>(f)) { }
        notifying_barrier(const notifying_barrier&) = delete;
        notifying_barrier& operator=(const notifying_barrier&) = delete;
        ~notifying_barrier( ) { }
        void arrive_and_wait() {
            int const myepoch = epoch.load(std::memory_order_relaxed);
            if(!__arrive(myepoch))
                while(epoch.load_when_not_equal(myepoch,std::memory_order_acquire) == myepoch)
                    ;
        }
        void arrive_and_drop() {
            nexpected.fetch_add(-1,std::memory_order_relaxed);
            __arrive(epoch.load(std::memory_order_relaxed));
        }
    private:
        bool __arrive(int const myepoch) {
            int const myresult = arrived.fetch_add(1,std::memory_order_acq_rel) + 1;
            if(__builtin_expect(myresult == expected,0)) {
                int const newexpected = completion();
                expected = newexpected ? newexpected : nexpected.load(std::memory_order_relaxed);
                arrived.store(0,std::memory_order_relaxed);
                epoch.store(myepoch+1,std::memory_order_release);
                return true;
            }
            return false;
        }
        int expected;
        std::atomic<int> arrived, nexpected;
        synchronic<int> epoch;
        std::function<int()> completion;
    };
}}

#endif //__N3998_H
