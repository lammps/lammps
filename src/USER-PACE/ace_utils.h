//
// Created by lysogy36 on 27.02.20.
//

#ifndef ACE_ACE_UTILS_H
#define ACE_ACE_UTILS_H

const DOUBLE_TYPE pi = 4.0 * atan(1.0);


inline int sign(DOUBLE_TYPE x) {
    if (x < 0) return -1;
    else if (x > 0) return +1;
    else return 0;
}

inline double absolute_relative_error(double x, double y, double zero_threshold = 5e-6) {
    if (x == 0 && y == 0) return 0;
    else if (x == 0 || y == 0) return (abs(x + y) < zero_threshold ? 0 : 2);
    else return 2 * abs(x - y) / (abs(x) + abs(y));
}

#endif //ACE_ACE_UTILS_H
