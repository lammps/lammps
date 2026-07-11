
#include <cctype>
#include <limits>

extern "C" {

#include "lmp_f2c.h"

// undefine conflicting f2c macros
#undef min
#undef max

doublereal dlamch_(const char *cmach)
{
    if (!cmach) return 0.0;
    char select = toupper(*cmach);

    // BLAS assumes rounding not truncation => epsilon is half
    const double eps = 0.5 * std::numeric_limits<double>::epsilon();
    if (select == 'E') return eps;

    double min = std::numeric_limits<double>::min();
    const double max = std::numeric_limits<double>::max();
    double small = 1.0 / max;
    if (small >= min) min = small * (1.0 + eps);
    if (select == 'S') return min;

    const double radix = std::numeric_limits<double>::radix;
    if (select == 'B') return radix;

    if (select == 'P') return eps * radix;

    if (select == 'N') return std::numeric_limits<double>::digits;

    if (select == 'M') return std::numeric_limits<double>::min_exponent;

    if (select == 'U') return min;

    if (select == 'L') return std::numeric_limits<double>::max_exponent;

    if (select == 'O') return max;

    return 0.0;
}
}
