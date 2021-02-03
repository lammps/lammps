/*
 * Performant implementation of atomic cluster expansion and interface to LAMMPS
 *
 * Copyright 2021  (c) Yury Lysogorskiy^1, Cas van der Oord^2, Anton Bochkarev^1,
 * Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
 * Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1
 *
 * ^1: Ruhr-University Bochum, Bochum, Germany
 * ^2: University of Cambridge, Cambridge, United Kingdom
 * ^3: Sandia National Laboratories, Albuquerque, New Mexico, USA
 * ^4: University of British Columbia, Vancouver, BC, Canada
 *
 *
 * See the LICENSE file.
 */


// Created by  Yury Lysogorskiy  on 19.02.20.

#ifndef ACE_TIMING_H
#define ACE_TIMING_H

#include <chrono>

using namespace std::chrono;
using Clock = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::time_point<Clock>;
using Duration = Clock::duration;

//////////////////////////////////////////
#ifdef FINE_TIMING
/**
 * Helper class for timing the code.
 * The timer should be initialized to reset measured time and
 * then call "start" and "stop" before and after measured code.
 * The measured time is stored in "duration" variable
 */
struct ACETimer {
    Duration duration; ///< measured duration
    TimePoint start_moment; ///< start moment of current measurement

    ACETimer() { init(); };

    /**
     * Reset timer
     */
    void init() { duration = std::chrono::nanoseconds(0); }

    /**
     * Start timer
     */
    void start() { start_moment = Clock::now(); }

    /**
     * Stop timer, update measured "duration"
     */
    void stop() { duration += Clock::now() - start_moment; }

    /**
     * Get duration in microseconds
     */
    long as_microseconds() { return std::chrono::duration_cast<std::chrono::microseconds>(duration).count(); }

    /**
     * Get duration in nanoseconds
     */
    long as_nanoseconds() { return std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count(); }

};

#else  // EMPTY Definitions
/**
 * Helper class for timing the code.
 * The timer should be initialized to reset measured time and
 * then call "start" and "stop" before and after measured code.
 * The measured time is stored in "duration" variable
 */
struct ACETimer {
    Duration duration; ///< measured duration
    TimePoint start_moment; ///< start moment of current measurement

    ACETimer() {};

    /**
    * Reset timer
    */
    void init() {}

    /**
    * Start timer
    */
    void start() {}

    /**
     * Stop timer, update measured "duration"
     */
    void stop() {}

    /**
    * Get duration in microseconds
    */
    long as_microseconds() {return 0; }

    /**
    * Get duration in nanoseconds
    */
    long as_nanoseconds() {return 0; }

};

#endif
//////////////////////////////////////////


#endif //ACE_TIMING_H
