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
 * This FILENAME is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
