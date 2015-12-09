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

#ifndef TEST_SYNCHRONIC_HPP
#define TEST_SYNCHRONIC_HPP

#include <impl/Kokkos_Synchronic.hpp>
#include <mutex>

namespace Test {

template <bool truly>
struct dumb_mutex {

    dumb_mutex () : locked(0) {
    }

    void lock() {
        while(1) {
            bool state = false;
            if (locked.compare_exchange_weak(state,true,std::memory_order_acquire)) {
                break;
            }
            while (locked.load(std::memory_order_relaxed)) {
              if (!truly) {
                Kokkos::Impl::portable_yield();
              }
            }
        }
    }

    void unlock() {
        locked.store(false,std::memory_order_release);
    }

private :
    std::atomic<bool> locked;
};

#ifdef WIN32
#include <winsock2.h>
#include <windows.h>
#include <synchapi.h>
struct srw_mutex {

    srw_mutex () {
        InitializeSRWLock(&_lock);
    }

    void lock() {
        AcquireSRWLockExclusive(&_lock);
    }
    void unlock() {
        ReleaseSRWLockExclusive(&_lock);
    }

private :
    SRWLOCK _lock;
};
#endif

struct ttas_mutex {

    ttas_mutex() : locked(false) {
    }

	ttas_mutex(const ttas_mutex&) = delete;
	ttas_mutex& operator=(const ttas_mutex&) = delete;

    void lock() {
        for(int i = 0;; ++i) {
            bool state = false;
            if(locked.compare_exchange_weak(state,true,std::memory_order_relaxed,Kokkos::Impl::notify_none))
                break;
            locked.expect_update(true);
        }
        std::atomic_thread_fence(std::memory_order_acquire);
    }
    void unlock() {
        locked.store(false,std::memory_order_release);
    }

private :
    Kokkos::Impl::synchronic<bool> locked;
};

struct ticket_mutex {

    ticket_mutex() : active(0), queue(0) {
    }

	ticket_mutex(const ticket_mutex&) = delete;
	ticket_mutex& operator=(const ticket_mutex&) = delete;

    void lock() {
        int const me = queue.fetch_add(1, std::memory_order_relaxed);
        while(me != active.load_when_equal(me, std::memory_order_acquire))
            ;
    }

    void unlock() {
        active.fetch_add(1,std::memory_order_release);
    }
private :
    Kokkos::Impl::synchronic<int> active;
    std::atomic<int> queue;
};

struct mcs_mutex {

    mcs_mutex() : head(nullptr) {
    }

	mcs_mutex(const mcs_mutex&) = delete;
	mcs_mutex& operator=(const mcs_mutex&) = delete;

    struct unique_lock {

        unique_lock(mcs_mutex & arg_m) : m(arg_m), next(nullptr), ready(false) {

            unique_lock * const h = m.head.exchange(this,std::memory_order_acquire);
            if(__builtin_expect(h != nullptr,0)) {
                h->next.store(this,std::memory_order_seq_cst,Kokkos::Impl::notify_one);
                while(!ready.load_when_not_equal(false,std::memory_order_acquire))
                    ;
            }
        }

	    unique_lock(const unique_lock&) = delete;
	    unique_lock& operator=(const unique_lock&) = delete;

        ~unique_lock() {
            unique_lock * h = this;
            if(__builtin_expect(!m.head.compare_exchange_strong(h,nullptr,std::memory_order_release, std::memory_order_relaxed),0)) {
                unique_lock * n = next.load(std::memory_order_relaxed);
                while(!n)
                    n = next.load_when_not_equal(n,std::memory_order_relaxed);
                n->ready.store(true,std::memory_order_release,Kokkos::Impl::notify_one);
            }
        }

    private:
        mcs_mutex & m;
        Kokkos::Impl::synchronic<unique_lock*> next;
        Kokkos::Impl::synchronic<bool> ready;
    };

private :
    std::atomic<unique_lock*> head;
};

}

namespace std {
template<>
struct unique_lock<Test::mcs_mutex> : Test::mcs_mutex::unique_lock {
  unique_lock(Test::mcs_mutex & arg_m) : Test::mcs_mutex::unique_lock(arg_m) {
  }
  unique_lock(const unique_lock&) = delete;
  unique_lock& operator=(const unique_lock&) = delete;
};

}

/* #include <cmath> */
#include <stdlib.h>

namespace Test {

//-------------------------------------
//  MersenneTwister
//-------------------------------------
#define MT_IA  397
#define MT_LEN 624

class MersenneTwister
{
    volatile unsigned long m_buffer[MT_LEN][64/sizeof(unsigned long)];
    volatile int m_index;

public:
    MersenneTwister() {
        for (int i = 0; i < MT_LEN; i++)
            m_buffer[i][0] = rand();
        m_index = 0;
        for (int i = 0; i < MT_LEN * 100; i++)
            integer();
    }
    unsigned long integer() {
        // Indices
        int i = m_index;
        int i2 = m_index + 1; if (i2 >= MT_LEN) i2 = 0; // wrap-around
        int j = m_index + MT_IA; if (j >= MT_LEN) j -= MT_LEN; // wrap-around

        // Twist
        unsigned long s = (m_buffer[i][0] & 0x80000000) | (m_buffer[i2][0] & 0x7fffffff);
        unsigned long r = m_buffer[j][0] ^ (s >> 1) ^ ((s & 1) * 0x9908B0DF);
        m_buffer[m_index][0] = r;
        m_index = i2;

        // Swizzle
        r ^= (r >> 11);
        r ^= (r << 7) & 0x9d2c5680UL;
        r ^= (r << 15) & 0xefc60000UL;
        r ^= (r >> 18);
        return r;
    }
    float poissonInterval(float ooLambda) {
        return -logf(1.0f - integer() * 2.3283e-10f) * ooLambda;
    }
};

} // namespace Test

#endif //TEST_HPP
