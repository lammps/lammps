
#include <cmath>

extern "C" {

#include "lmp_f2c.h"

logical disnan_(const doublereal *din)
{
    if (!din) return TRUE_;

    return std::isnan(*din) ? TRUE_ : FALSE_;
}
}
