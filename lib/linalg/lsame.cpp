
#include <cctype>

extern "C" {

#include "lmp_f2c.h"

logical lsame_(const char *a, const char *b)
{
    char ua, ub;
    if (!a || !b) return FALSE_;

    ua = toupper(*a);
    ub = toupper(*b);
    return (ua == ub) ? TRUE_ : FALSE_;
}
}
