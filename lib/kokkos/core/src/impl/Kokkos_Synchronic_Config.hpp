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

#ifndef KOKKOS_SYNCHRONIC_CONFIG_H
#define KOKKOS_SYNCHRONIC_CONFIG_H

#include <thread>
#include <chrono>

namespace Kokkos {
namespace Impl {

//the default yield function used inside the implementation is the Standard one
#define __synchronic_yield std::this_thread::yield
#define __synchronic_relax __synchronic_yield

#if defined(_MSC_VER)
    //this is a handy GCC optimization that I use inside the implementation
    #define __builtin_expect(condition,common) condition
    #if _MSC_VER <= 1800
        //using certain keywords that VC++ temporarily doesn't support
        #define _ALLOW_KEYWORD_MACROS
        #define noexcept
        #define constexpr
    #endif
    //yes, I define multiple assignment operators
    #pragma warning(disable:4522)
    //I don't understand how Windows is so bad at timing functions, but is OK
    //with straight-up yield loops
    #define __do_backoff(b) __synchronic_yield()
#else
#define __do_backoff(b) b.sleep_for_step()
#endif

//certain platforms have efficient support for spin-waiting built into the operating system
#if defined(__linux__) || (defined(_WIN32_WINNT) && _WIN32_WINNT >= 0x0602)
#if defined(_WIN32_WINNT)
#include <winsock2.h>
#include <Windows.h>
    //the combination of WaitOnAddress and WakeByAddressAll is supported on Windows 8.1+
    #define __synchronic_wait(x,v) WaitOnAddress((PVOID)x,(PVOID)&v,sizeof(v),-1)
    #define __synchronic_wait_timed(x,v,t) WaitOnAddress((PVOID)x,(PVOID)&v,sizeof(v),std::chrono::duration_cast<std::chrono::milliseconds>(t).count())
    #define __synchronic_wake_one(x) WakeByAddressSingle((PVOID)x)
    #define __synchronic_wake_all(x) WakeByAddressAll((PVOID)x)
    #define __synchronic_wait_volatile(x,v) WaitOnAddress((PVOID)x,(PVOID)&v,sizeof(v),-1)
    #define __synchronic_wait_timed_volatile(x,v,t) WaitOnAddress((PVOID)x,(PVOID)&v,sizeof(v),std::chrono::duration_cast<std::chrono::milliseconds>(t).count())
    #define __synchronic_wake_one_volatile(x) WakeByAddressSingle((PVOID)x)
    #define __synchronic_wake_all_volatile(x) WakeByAddressAll((PVOID)x)
    #define __SYNCHRONIC_COMPATIBLE(x) (std::is_pod<x>::value && (sizeof(x) <= 8))

    inline void native_sleep(unsigned long microseconds)
    {
      // What to do if microseconds is < 1000?
      Sleep(microseconds / 1000);
    }

    inline void native_yield()
    {
      SwitchToThread();
    }
#elif defined(__linux__)
    #include <chrono>
    #include <time.h>
    #include <unistd.h>
    #include <pthread.h>
    #include <linux/futex.h>
    #include <sys/syscall.h>
    #include <climits>
    #include <cassert>
    template < class Rep, class Period>
    inline timespec to_timespec(std::chrono::duration<Rep,Period> const& delta) {
      struct timespec ts;
      ts.tv_sec = static_cast<long>(std::chrono::duration_cast<std::chrono::seconds>(delta).count());
      assert(!ts.tv_sec);
      ts.tv_nsec = static_cast<long>(std::chrono::duration_cast<std::chrono::nanoseconds>(delta).count());
      return ts;
    }
    inline long futex(void const* addr1, int op, int val1) {
        return syscall(SYS_futex, addr1, op, val1, 0, 0, 0);
    }
    inline long futex(void const* addr1, int op, int val1, struct timespec timeout) {
        return syscall(SYS_futex, addr1, op, val1, &timeout, 0, 0);
    }
    inline void native_sleep(unsigned long microseconds)
    {
      usleep(microseconds);
    }
    inline void native_yield()
    {
      pthread_yield();
    }

    //the combination of SYS_futex(WAIT) and SYS_futex(WAKE) is supported on all recent Linux distributions
    #define __synchronic_wait(x,v) futex(x, FUTEX_WAIT_PRIVATE, v)
    #define __synchronic_wait_timed(x,v,t) futex(x, FUTEX_WAIT_PRIVATE, v, to_timespec(t))
    #define __synchronic_wake_one(x) futex(x, FUTEX_WAKE_PRIVATE, 1)
    #define __synchronic_wake_all(x) futex(x, FUTEX_WAKE_PRIVATE, INT_MAX)
    #define __synchronic_wait_volatile(x,v) futex(x, FUTEX_WAIT, v)
    #define __synchronic_wait_volatile_timed(x,v,t) futex(x, FUTEX_WAIT, v, to_timespec(t))
    #define __synchronic_wake_one_volatile(x) futex(x, FUTEX_WAKE, 1)
    #define __synchronic_wake_all_volatile(x) futex(x, FUTEX_WAKE, INT_MAX)
    #define __SYNCHRONIC_COMPATIBLE(x) (std::is_integral<x>::value && (sizeof(x) <= 4))

    //the yield function on Linux is better replaced by sched_yield, which is tuned for spin-waiting
    #undef __synchronic_yield
    #define __synchronic_yield sched_yield

    //for extremely short wait times, just let another hyper-thread run
    #undef __synchronic_relax
    #define __synchronic_relax() asm volatile("rep; nop" ::: "memory")

#endif
#endif

#ifdef _GLIBCXX_USE_NANOSLEEP
inline void portable_sleep(std::chrono::microseconds const& time)
{ std::this_thread::sleep_for(time); }
#else
inline void portable_sleep(std::chrono::microseconds const& time)
{ native_sleep(time.count()); }
#endif

#ifdef _GLIBCXX_USE_SCHED_YIELD
inline void portable_yield()
{ std::this_thread::yield(); }
#else
inline void portable_yield()
{ native_yield(); }
#endif

//this is the number of times we initially spin, on the first wait attempt
#define __SYNCHRONIC_SPIN_COUNT_A 16

//this is how decide to yield instead of just spinning, 'c' is the current trip count
//#define __SYNCHRONIC_SPIN_YIELD(c) true
#define __SYNCHRONIC_SPIN_RELAX(c) (c>>3)

//this is the number of times we normally spin, on every subsequent wait attempt
#define __SYNCHRONIC_SPIN_COUNT_B 8

}
}

#endif //__SYNCHRONIC_CONFIG_H
