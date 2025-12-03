//
// Created by Yury Lysogorskiy on 06.12.23.
//

#ifndef LAMMPS_PACE_UTILS_H
#define LAMMPS_PACE_UTILS_H

#include <cstring>
#include <chrono>

#include "ace-evaluator/ace_arraynd.h"

namespace PACE {
    static char const *const elements_pace[] = {
            "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
            "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
            "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
            "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr",
            "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
            "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
            "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};
    static constexpr int elements_num_pace = sizeof(elements_pace) / sizeof(const char *);

    static int AtomicNumberByName(char *elname) {
        for (int i = 1; i < elements_num_pace; i++)
            if (strcmp(elname, elements_pace[i]) == 0) return i;
        return -1;
    }

    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    using Duration = Clock::duration;

//////////////////////////////////////////
/**
 * Helper class for timing the code.
 * The timer should be initialized to reset measured time and
 * then call "start" and "stop" before and after measured code.
 * The measured time is stored in "duration" variable
 */
    struct ACETimer {
        Duration duration{}; ///< measured duration
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
        long as_microseconds() const { return std::chrono::duration_cast<std::chrono::microseconds>(duration).count(); }

        /**
         * Get duration in nanoseconds
         */
        long as_nanoseconds() const { return std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count(); }

    };


    inline Array2D<double> inv3x3(Array2D<double>& a) {
        DOUBLE_TYPE det = a(0,0) * (a(2,2) * a(1,1) - a(2,1) * a(1,2))
                          - a(1,0) * (a(2,2) * a(0,1) - a(2,1) * a(0,2))
                          + a(2,0) * (a(1,2) * a(0,1) - a(1,1) * a(0,2));
        if (fabs(det) < 1e-9)
            throw std::overflow_error("Couldn't invert matrix - determinant is almost zero");

        Array2D<double> ainv(3,3);
        ainv.fill(0);

        ainv(0,0) = (a(2,2) * a(1,1) - a(2,1) * a(1,2)) / det;
        ainv(0,1) = -(a(2,2) * a(0,1) - a(2,1) * a(0,2)) / det;
        ainv(0,2) = (a(1,2) * a(0,1) - a(1,1) * a(0,2)) / det;

        ainv(1,0) = -(a(2,2) * a(1,0) - a(2,0) * a(1,2)) / det;
        ainv(1,1) = (a(2,2) * a(0,0) - a(2,0) * a(0,2)) / det;
        ainv(1,2) = -(a(1,2) * a(0,0) - a(1,0) * a(0,2)) / det;

        ainv(2,0) = (a(2,1) * a(1,0) - a(2,0) * a(1,1)) / det;
        ainv(2,1) = -(a(2,1) * a(0,0) - a(2,0) * a(0,1)) / det;
        ainv(2,2) = (a(1,1) * a(0,0) - a(1,0) * a(0,1)) / det;

        return ainv;
    };

}

#endif //LAMMPS_PACE_UTILS_H
