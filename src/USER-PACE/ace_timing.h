//
// Created by  Yury Lysogorskiy  on 19.02.20.
//

#ifndef ACE_TIMING_H
#define ACE_TIMING_H
//keywords for macros
#define __start_moment __start_moment
#define __duration __duration

//////////////////////////////////////////
#ifdef FINE_TIMING

#include <chrono>

using namespace std::chrono;
using Clock = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::time_point<Clock>;
using Duration = Clock::duration;


//timer definition macros
#define ACE_DEFINE_TIMER(timer_name) TimePoint timer_name ## __start_moment;\
Duration timer_name ## __duration;

//timer init macros
#define ACE_TIMER_INIT(timer_name) timer_name ## __duration  = std::chrono::nanoseconds(0);

//timer start macros
#define ACE_TIMER_START(timer_name) timer_name ## __start_moment = Clock::now();

//timer stop macros
#define ACE_TIMER_STOP(timer_name) timer_name##__duration += Clock::now() - timer_name##__start_moment;

#define ACE_TIMER_MICROSECONDS(timer_name)  (std::chrono::duration_cast<std::chrono::microseconds>(timer_name##__duration).count())
#define ACE_TIMER_NANOSECONDS(timer_name)  (std::chrono::duration_cast<std::chrono::nanoseconds>(timer_name##__duration).count())

#define TIMER_MICROSECONDS_FROM(timer_name, obj)  (std::chrono::duration_cast<std::chrono::microseconds>(obj.timer_name##__duration).count())
#define TIMER_NANOSECONDS_FROM(timer_name, obj)  (std::chrono::duration_cast<std::chrono::nanoseconds>(obj.timer_name##__duration).count())

#else  /// EMPTY Definitions
//timer definition macros
#define ACE_DEFINE_TIMER(timer_name)
//timer init macros
#define ACE_TIMER_INIT(timer_name)

//timer start macros
#define ACE_TIMER_START(timer_name)

//timer stop macros
#define ACE_TIMER_STOP(timer_name)

#define ACE_TIMER_MICROSECONDS(timer_name)
#define ACE_TIMER_NANOSECONDS(timer_name)

#endif
//////////////////////////////////////////


#endif //ACE_TIMING_H
