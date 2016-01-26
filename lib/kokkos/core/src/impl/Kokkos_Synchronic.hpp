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

#ifndef KOKKOS_SYNCHRONIC_HPP
#define KOKKOS_SYNCHRONIC_HPP

#include <impl/Kokkos_Synchronic_Config.hpp>

#include <atomic>
#include <chrono>
#include <thread>
#include <functional>
#include <algorithm>

namespace Kokkos {
namespace Impl {

enum notify_hint {
  notify_all,
  notify_one,
  notify_none
};
enum expect_hint {
  expect_urgent,
  expect_delay
};

namespace Details {

template <class S, class T>
bool __synchronic_spin_wait_for_update(S const& arg, T const& nval, int attempts) noexcept {
  int i = 0;
  for(;i < __SYNCHRONIC_SPIN_RELAX(attempts); ++i)
    if(__builtin_expect(arg.load(std::memory_order_relaxed) != nval,1))
      return true;
    else
      __synchronic_relax();
  for(;i < attempts; ++i)
    if(__builtin_expect(arg.load(std::memory_order_relaxed) != nval,1))
      return true;
    else
      __synchronic_yield();
  return false;
}

struct __exponential_backoff {
  __exponential_backoff(int arg_maximum=512) : maximum(arg_maximum), microseconds(8), x(123456789), y(362436069), z(521288629) {
  }
  static inline void sleep_for(std::chrono::microseconds const& time) {
    auto t = time.count();
    if(__builtin_expect(t > 75,0)) {
      portable_sleep(time);
    }
    else if(__builtin_expect(t > 25,0))
      __synchronic_yield();
    else
      __synchronic_relax();
  }
  void sleep_for_step() {
    sleep_for(step());
  }
  std::chrono::microseconds step() {
    float const f = ranfu();
    int const t = int(microseconds * f);
    if(__builtin_expect(f >= 0.95f,0))
      microseconds = 8;
    else
      microseconds = (std::min)(microseconds>>1,maximum);
    return std::chrono::microseconds(t);
  }
private :
  int maximum, microseconds, x, y, z;
  int xorshf96() {
    int t;
    x ^= x << 16; x ^= x >> 5; x ^= x << 1;
    t = x; x = y; y = z; z = t ^ x ^ y;
    return z;
  }
  float ranfu() {
    return (float)(xorshf96()&(~0UL>>1)) / (float)(~0UL>>1);
  }
};

template <class T, class Enable = void>
struct __synchronic_base {

protected:
  std::atomic<T> atom;

  void notify(notify_hint = notify_all) noexcept {
  }
  void notify(notify_hint = notify_all) volatile noexcept {
  }

public :
  __synchronic_base() noexcept = default;
  constexpr __synchronic_base(T v) noexcept : atom(v) { }
  __synchronic_base(const __synchronic_base&) = delete;
  ~__synchronic_base() { }
  __synchronic_base& operator=(const __synchronic_base&) = delete;
  __synchronic_base& operator=(const __synchronic_base&) volatile = delete;

  void expect_update(T val, expect_hint = expect_urgent) const noexcept {
    if(__synchronic_spin_wait_for_update(atom, val, __SYNCHRONIC_SPIN_COUNT_A))
      return;
    __exponential_backoff b;
    while(atom.load(std::memory_order_relaxed) == val) {
      __do_backoff(b);
      if(__synchronic_spin_wait_for_update(atom, val, __SYNCHRONIC_SPIN_COUNT_B))
        return;
    }
  }
  void expect_update(T val, expect_hint = expect_urgent) const volatile noexcept {
    if(__synchronic_spin_wait_for_update(atom, val, __SYNCHRONIC_SPIN_COUNT_A))
      return;
    __exponential_backoff b;
    while(atom.load(std::memory_order_relaxed) == val) {
      __do_backoff(b);
      if(__synchronic_spin_wait_for_update(atom, val, __SYNCHRONIC_SPIN_COUNT_B))
        return;
    }
  }

  template <class Clock, class Duration>
  void expect_update_until(T val, std::chrono::time_point<Clock,Duration> const& then, expect_hint = expect_urgent) const {
    if(__synchronic_spin_wait_for_update(atom, val, __SYNCHRONIC_SPIN_COUNT_A))
      return;
    __exponential_backoff b;
    std::chrono::milliseconds remains = then - std::chrono::high_resolution_clock::now();
    while(remains > std::chrono::milliseconds::zero() && atom.load(std::memory_order_relaxed) == val) {
      __do_backoff(b);
      if(__synchronic_spin_wait_for_update(atom, val, __SYNCHRONIC_SPIN_COUNT_B))
        return;
      remains = then - std::chrono::high_resolution_clock::now();
    }
  }
  template <class Clock, class Duration>
  void expect_update_until(T val, std::chrono::time_point<Clock,Duration> const& then, expect_hint = expect_urgent) const volatile {
    if(__synchronic_spin_wait_for_update(atom, val, __SYNCHRONIC_SPIN_COUNT_A))
      return;
    __exponential_backoff b;
    std::chrono::milliseconds remains = then - std::chrono::high_resolution_clock::now();
    while(remains > std::chrono::milliseconds::zero() && atom.load(std::memory_order_relaxed) == val) {
      __do_backoff(b);
      if(__synchronic_spin_wait_for_update(atom, val, __SYNCHRONIC_SPIN_COUNT_B))
        return;
      remains = then - std::chrono::high_resolution_clock::now();
    }
  }
};

#ifdef __SYNCHRONIC_COMPATIBLE
template <class T>
struct __synchronic_base<T, typename std::enable_if<__SYNCHRONIC_COMPATIBLE(T)>::type> {

public:
  std::atomic<T> atom;

  void notify(notify_hint hint = notify_all) noexcept {
    if(__builtin_expect(hint == notify_none,1))
      return;
    auto const x = count.fetch_add(0,std::memory_order_acq_rel);
    if(__builtin_expect(x,0)) {
      if(__builtin_expect(hint == notify_all,1))
        __synchronic_wake_all(&atom);
      else
        __synchronic_wake_one(&atom);
    }
  }
  void notify(notify_hint hint = notify_all) volatile noexcept {
    if(__builtin_expect(hint == notify_none,1))
      return;
    auto const x = count.fetch_add(0,std::memory_order_acq_rel);
    if(__builtin_expect(x,0)) {
      if(__builtin_expect(hint == notify_all,1))
        __synchronic_wake_all_volatile(&atom);
      else
        __synchronic_wake_one_volatile(&atom);
    }
  }

public :
  __synchronic_base() noexcept : count(0) { }
  constexpr __synchronic_base(T v) noexcept : atom(v), count(0) { }
  __synchronic_base(const __synchronic_base&) = delete;
  ~__synchronic_base() { }
  __synchronic_base& operator=(const __synchronic_base&) = delete;
  __synchronic_base& operator=(const __synchronic_base&) volatile = delete;

  void expect_update(T val, expect_hint = expect_urgent) const noexcept {
    if(__builtin_expect(__synchronic_spin_wait_for_update(atom, val,__SYNCHRONIC_SPIN_COUNT_A),1))
      return;
    while(__builtin_expect(atom.load(std::memory_order_relaxed) == val,1)) {
      count.fetch_add(1,std::memory_order_release);
      __synchronic_wait(&atom,val);
      count.fetch_add(-1,std::memory_order_acquire);
    }
  }
  void expect_update(T val, expect_hint = expect_urgent) const volatile noexcept {
    if(__builtin_expect(__synchronic_spin_wait_for_update(atom, val,__SYNCHRONIC_SPIN_COUNT_A),1))
      return;
    while(__builtin_expect(atom.load(std::memory_order_relaxed) == val,1)) {
      count.fetch_add(1,std::memory_order_release);
      __synchronic_wait_volatile(&atom,val);
      count.fetch_add(-1,std::memory_order_acquire);
    }
  }

  template <class Clock, class Duration>
  void expect_update_until(T val, std::chrono::time_point<Clock,Duration> const& then, expect_hint = expect_urgent) const {
    if(__builtin_expect(__synchronic_spin_wait_for_update(atom, val,__SYNCHRONIC_SPIN_COUNT_A),1))
      return;
    std::chrono::milliseconds remains = then - std::chrono::high_resolution_clock::now();
    while(__builtin_expect(remains > std::chrono::milliseconds::zero() && atom.load(std::memory_order_relaxed) == val,1)) {
      count.fetch_add(1,std::memory_order_release);
      __synchronic_wait_timed(&atom,val,remains);
      count.fetch_add(-1,std::memory_order_acquire);
      remains = then - std::chrono::high_resolution_clock::now();
    }
  }
  template <class Clock, class Duration>
  void expect_update_until(T val, std::chrono::time_point<Clock,Duration> const& then, expect_hint = expect_urgent) const volatile {
    if(__builtin_expect(__synchronic_spin_wait_for_update(atom, val,__SYNCHRONIC_SPIN_COUNT_A),1))
      return;
    std::chrono::milliseconds remains = then - std::chrono::high_resolution_clock::now();
    while(__builtin_expect(remains > std::chrono::milliseconds::zero() && atom.load(std::memory_order_relaxed) == val,1)) {
      count.fetch_add(1,std::memory_order_release);
      __synchronic_wait_timed_volatile(&atom,val,remains);
      count.fetch_add(-1,std::memory_order_acquire);
      remains = then - std::chrono::high_resolution_clock::now();
    }
  }
private:
  mutable std::atomic<int> count;
};
#endif

template <class T, class Enable = void>
struct __synchronic : public __synchronic_base<T> {

  __synchronic() noexcept = default;
  constexpr __synchronic(T v) noexcept : __synchronic_base<T>(v) { }
  __synchronic(const __synchronic&) = delete;
  __synchronic& operator=(const __synchronic&) = delete;
  __synchronic& operator=(const __synchronic&) volatile = delete;
};

template <class T>
struct __synchronic<T,typename std::enable_if<std::is_integral<T>::value>::type> : public __synchronic_base<T> {

  T fetch_add(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) volatile noexcept {
    auto const t = this->atom.fetch_add(v,m);
    this->notify(n);
    return t;
  }
  T fetch_add(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) noexcept {
    auto const t = this->atom.fetch_add(v,m);
    this->notify(n);
    return t;
  }
  T fetch_sub(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) volatile noexcept {
    auto const t = this->atom.fetch_sub(v,m);
    this->notify(n);
    return t;
  }
  T fetch_sub(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) noexcept {
    auto const t = this->atom.fetch_sub(v,m);
    this->notify(n);
    return t;
  }
  T fetch_and(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) volatile noexcept {
    auto const t = this->atom.fetch_and(v,m);
    this->notify(n);
    return t;
  }
  T fetch_and(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) noexcept {
    auto const t = this->atom.fetch_and(v,m);
    this->notify(n);
    return t;
  }
  T fetch_or(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) volatile noexcept {
    auto const t = this->atom.fetch_or(v,m);
    this->notify(n);
    return t;
  }
  T fetch_or(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) noexcept {
    auto const t = this->atom.fetch_or(v,m);
    this->notify(n);
    return t;
  }
  T fetch_xor(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) volatile noexcept {
    auto const t = this->atom.fetch_xor(v,m);
    this->notify(n);
    return t;
  }
  T fetch_xor(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) noexcept {
    auto const t = this->atom.fetch_xor(v,m);
    this->notify(n);
    return t;
  }

  __synchronic() noexcept = default;
  constexpr __synchronic(T v) noexcept : __synchronic_base<T>(v) { }
  __synchronic(const __synchronic&) = delete;
  __synchronic& operator=(const __synchronic&) = delete;
  __synchronic& operator=(const __synchronic&) volatile = delete;

  T operator=(T v) volatile noexcept {
    auto const t = this->atom = v;
    this->notify();
    return t;
  }
  T operator=(T v) noexcept {
    auto const t = this->atom = v;
    this->notify();
    return t;
  }
  T operator++(int) volatile noexcept {
    auto const t = ++this->atom;
    this->notify();
    return t;
  }
  T operator++(int) noexcept {
    auto const t = ++this->atom;
    this->notify();
    return t;
  }
  T operator--(int) volatile noexcept {
    auto const t = --this->atom;
    this->notify();
    return t;
  }
  T operator--(int) noexcept {
    auto const t = --this->atom;
    this->notify();
    return t;
  }
  T operator++() volatile noexcept {
    auto const t = this->atom++;
    this->notify();
    return t;
  }
  T operator++() noexcept {
    auto const t = this->atom++;
    this->notify();
    return t;
  }
  T operator--() volatile noexcept {
    auto const t = this->atom--;
    this->notify();
    return t;
  }
  T operator--() noexcept {
    auto const t = this->atom--;
    this->notify();
    return t;
  }
  T operator+=(T v) volatile noexcept {
    auto const t = this->atom += v;
    this->notify();
    return t;
  }
  T operator+=(T v) noexcept {
    auto const t = this->atom += v;
    this->notify();
    return t;
  }
  T operator-=(T v) volatile noexcept {
    auto const t = this->atom -= v;
    this->notify();
    return t;
  }
  T operator-=(T v) noexcept {
    auto const t = this->atom -= v;
    this->notify();
    return t;
  }
  T operator&=(T v) volatile noexcept {
    auto const t = this->atom &= v;
    this->notify();
    return t;
  }
  T operator&=(T v) noexcept {
    auto const t = this->atom &= v;
    this->notify();
    return t;
  }
  T operator|=(T v) volatile noexcept {
    auto const t = this->atom |= v;
    this->notify();
    return t;
  }
  T operator|=(T v) noexcept {
    auto const t = this->atom |= v;
    this->notify();
    return t;
  }
  T operator^=(T v) volatile noexcept {
    auto const t = this->atom ^= v;
    this->notify();
    return t;
  }
  T operator^=(T v) noexcept {
    auto const t = this->atom ^= v;
    this->notify();
    return t;
  }
};

template <class T>
struct __synchronic<T*> : public __synchronic_base<T*> {

  T* fetch_add(ptrdiff_t v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) volatile noexcept {
    auto const t = this->atom.fetch_add(v,m);
    this->notify(n);
    return t;
  }
  T* fetch_add(ptrdiff_t v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) noexcept {
    auto const t = this->atom.fetch_add(v,m);
    this->notify(n);
    return t;
  }
  T* fetch_sub(ptrdiff_t v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) volatile noexcept {
    auto const t = this->atom.fetch_sub(v,m);
    this->notify(n);
    return t;
  }
  T* fetch_sub(ptrdiff_t v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) noexcept {
    auto const t = this->atom.fetch_sub(v,m);
    this->notify(n);
    return t;
  }

  __synchronic() noexcept = default;
  constexpr __synchronic(T* v) noexcept : __synchronic_base<T*>(v) { }
  __synchronic(const __synchronic&) = delete;
  __synchronic& operator=(const __synchronic&) = delete;
  __synchronic& operator=(const __synchronic&) volatile = delete;

  T* operator=(T* v) volatile noexcept {
    auto const t = this->atom = v;
    this->notify();
    return t;
  }
  T* operator=(T* v) noexcept {
    auto const t = this->atom = v;
    this->notify();
    return t;
  }
  T* operator++(int) volatile noexcept {
    auto const t = ++this->atom;
    this->notify();
    return t;
  }
  T* operator++(int) noexcept {
    auto const t = ++this->atom;
    this->notify();
    return t;
  }
  T* operator--(int) volatile noexcept {
    auto const t = --this->atom;
    this->notify();
    return t;
  }
  T* operator--(int) noexcept {
    auto const t = --this->atom;
    this->notify();
    return t;
  }
  T* operator++() volatile noexcept {
    auto const t = this->atom++;
    this->notify();
    return t;
  }
  T* operator++() noexcept {
    auto const t = this->atom++;
    this->notify();
    return t;
  }
  T* operator--() volatile noexcept {
    auto const t = this->atom--;
    this->notify();
    return t;
  }
  T* operator--() noexcept {
    auto const t = this->atom--;
    this->notify();
    return t;
  }
  T* operator+=(ptrdiff_t v) volatile noexcept {
    auto const t = this->atom += v;
    this->notify();
    return t;
  }
  T* operator+=(ptrdiff_t v) noexcept {
    auto const t = this->atom += v;
    this->notify();
    return t;
  }
  T* operator-=(ptrdiff_t v) volatile noexcept {
    auto const t = this->atom -= v;
    this->notify();
    return t;
  }
  T* operator-=(ptrdiff_t v) noexcept {
    auto const t = this->atom -= v;
    this->notify();
    return t;
  }
};

} //namespace Details

template <class T>
struct synchronic : public Details::__synchronic<T> {

  bool is_lock_free() const volatile noexcept { return this->atom.is_lock_free(); }
  bool is_lock_free() const noexcept { return this->atom.is_lock_free(); }
  void store(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) volatile noexcept {
    this->atom.store(v,m);
    this->notify(n);
  }
  void store(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) noexcept {
    this->atom.store(v,m);
    this->notify(n);
  }
  T load(std::memory_order m = std::memory_order_seq_cst) const volatile noexcept { return this->atom.load(m); }
  T load(std::memory_order m = std::memory_order_seq_cst) const noexcept { return this->atom.load(m); }

  operator T() const volatile noexcept { return (T)this->atom; }
  operator T() const noexcept { return (T)this->atom; }

  T exchange(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) volatile noexcept {
    auto const t = this->atom.exchange(v,m);
    this->notify(n);
    return t;
  }
  T exchange(T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) noexcept {
    auto const t = this->atom.exchange(v,m);
    this->notify(n);
    return t;
  }
  bool compare_exchange_weak(T& r, T v, std::memory_order m1, std::memory_order m2, notify_hint n = notify_all) volatile noexcept {
    auto const t = this->atom.compare_exchange_weak(r,v,m1,m2);
    this->notify(n);
    return t;
  }
  bool compare_exchange_weak(T& r, T v, std::memory_order m1, std::memory_order m2, notify_hint n = notify_all) noexcept {
    auto const t = this->atom.compare_exchange_weak(r,v,m1, m2);
    this->notify(n);
    return t;
  }
  bool compare_exchange_strong(T& r, T v, std::memory_order m1, std::memory_order m2, notify_hint n = notify_all) volatile noexcept {
    auto const t = this->atom.compare_exchange_strong(r,v,m1,m2);
    this->notify(n);
    return t;
  }
  bool compare_exchange_strong(T& r, T v, std::memory_order m1, std::memory_order m2, notify_hint n = notify_all) noexcept {
    auto const t = this->atom.compare_exchange_strong(r,v,m1,m2);
    this->notify(n);
    return t;
  }
  bool compare_exchange_weak(T& r, T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) volatile noexcept {
    auto const t = this->atom.compare_exchange_weak(r,v,m);
    this->notify(n);
    return t;
  }
  bool compare_exchange_weak(T& r, T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) noexcept {
    auto const t = this->atom.compare_exchange_weak(r,v,m);
    this->notify(n);
    return t;
  }
  bool compare_exchange_strong(T& r, T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) volatile noexcept {
    auto const t = this->atom.compare_exchange_strong(r,v,m);
    this->notify(n);
    return t;
  }
  bool compare_exchange_strong(T& r, T v, std::memory_order m = std::memory_order_seq_cst, notify_hint n = notify_all) noexcept {
    auto const t = this->atom.compare_exchange_strong(r,v,m);
    this->notify(n);
    return t;
  }

  synchronic() noexcept = default;
  constexpr synchronic(T val) noexcept : Details::__synchronic<T>(val) { }
  synchronic(const synchronic&) = delete;
  ~synchronic() { }
  synchronic& operator=(const synchronic&) = delete;
  synchronic& operator=(const synchronic&) volatile = delete;
  T operator=(T val) noexcept {
    return Details::__synchronic<T>::operator=(val);
  }
  T operator=(T val) volatile noexcept {
    return Details::__synchronic<T>::operator=(val);
  }

  T load_when_not_equal(T val, std::memory_order order = std::memory_order_seq_cst, expect_hint h = expect_urgent) const noexcept {
    Details::__synchronic<T>::expect_update(val,h);
    return load(order);
  }
  T load_when_not_equal(T val, std::memory_order order = std::memory_order_seq_cst, expect_hint h = expect_urgent) const volatile noexcept {
    Details::__synchronic<T>::expect_update(val,h);
    return load(order);
  }
  T load_when_equal(T val, std::memory_order order = std::memory_order_seq_cst, expect_hint h = expect_urgent) const noexcept {
    for(T nval = load(std::memory_order_relaxed); nval != val; nval = load(std::memory_order_relaxed))
      Details::__synchronic<T>::expect_update(nval,h);
    return load(order);
  }
  T load_when_equal(T val, std::memory_order order = std::memory_order_seq_cst, expect_hint h = expect_urgent) const volatile noexcept {
    for(T nval = load(std::memory_order_relaxed); nval != val; nval = load(std::memory_order_relaxed))
      expect_update(nval,h);
    return load(order);
  }
  template <class Rep, class Period>
  void expect_update_for(T val, std::chrono::duration<Rep,Period> const& delta, expect_hint h = expect_urgent) const {
    Details::__synchronic<T>::expect_update_until(val, std::chrono::high_resolution_clock::now() + delta,h);
  }
  template < class Rep, class Period>
  void expect_update_for(T val, std::chrono::duration<Rep,Period> const& delta, expect_hint h = expect_urgent) const volatile {
    Details::__synchronic<T>::expect_update_until(val, std::chrono::high_resolution_clock::now() + delta,h);
  }
};

#include <inttypes.h>

typedef synchronic<char> synchronic_char;
typedef synchronic<char> synchronic_schar;
typedef synchronic<unsigned char> synchronic_uchar;
typedef synchronic<short> synchronic_short;
typedef synchronic<unsigned short> synchronic_ushort;
typedef synchronic<int> synchronic_int;
typedef synchronic<unsigned int> synchronic_uint;
typedef synchronic<long> synchronic_long;
typedef synchronic<unsigned long> synchronic_ulong;
typedef synchronic<long long> synchronic_llong;
typedef synchronic<unsigned long long> synchronic_ullong;
//typedef synchronic<char16_t> synchronic_char16_t;
//typedef synchronic<char32_t> synchronic_char32_t;
typedef synchronic<wchar_t> synchronic_wchar_t;

typedef synchronic<int_least8_t> synchronic_int_least8_t;
typedef synchronic<uint_least8_t> synchronic_uint_least8_t;
typedef synchronic<int_least16_t> synchronic_int_least16_t;
typedef synchronic<uint_least16_t> synchronic_uint_least16_t;
typedef synchronic<int_least32_t> synchronic_int_least32_t;
typedef synchronic<uint_least32_t> synchronic_uint_least32_t;
//typedef synchronic<int_least_64_t> synchronic_int_least_64_t;
typedef synchronic<uint_least64_t> synchronic_uint_least64_t;
typedef synchronic<int_fast8_t> synchronic_int_fast8_t;
typedef synchronic<uint_fast8_t> synchronic_uint_fast8_t;
typedef synchronic<int_fast16_t> synchronic_int_fast16_t;
typedef synchronic<uint_fast16_t> synchronic_uint_fast16_t;
typedef synchronic<int_fast32_t> synchronic_int_fast32_t;
typedef synchronic<uint_fast32_t> synchronic_uint_fast32_t;
typedef synchronic<int_fast64_t> synchronic_int_fast64_t;
typedef synchronic<uint_fast64_t> synchronic_uint_fast64_t;
typedef synchronic<intptr_t> synchronic_intptr_t;
typedef synchronic<uintptr_t> synchronic_uintptr_t;
typedef synchronic<size_t> synchronic_size_t;
typedef synchronic<ptrdiff_t> synchronic_ptrdiff_t;
typedef synchronic<intmax_t> synchronic_intmax_t;
typedef synchronic<uintmax_t> synchronic_uintmax_t;

}
}

#endif //__SYNCHRONIC_H
