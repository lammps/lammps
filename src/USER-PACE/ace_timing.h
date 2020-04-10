//
// Created by  Yury Lysogorskiy  on 19.02.20.
//

#ifndef ACE_TIMING_H
#define ACE_TIMING_H
//keywords for macros
#define __start_moment __start_moment
#define __duration __duration


#include <chrono>

using namespace std::chrono;
using Clock = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::time_point<Clock>;
using Duration = Clock::duration;

//////////////////////////////////////////
#ifdef FINE_TIMING

struct ACETimer {
    Duration duration;
    TimePoint start_moment;

    ACETimer() { init(); };

    void init() { duration = std::chrono::nanoseconds(0); }

    void start() { start_moment = Clock::now(); }

    void stop() { duration += Clock::now() - start_moment; }

    long as_microseconds() { return std::chrono::duration_cast<std::chrono::microseconds>(duration).count(); }

    long as_nanoseconds() { return std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count(); }

};

#else  /// EMPTY Definitions

struct ACETimer {
    Duration duration;
    TimePoint start_moment;

    ACETimer() {};

    void init() {}

    void start() {}

    void stop() {}

    long as_microseconds() { return 0; }

    long as_nanoseconds() { return 0; }

};

#endif
//////////////////////////////////////////


#endif //ACE_TIMING_H
