
#include "error_stats.h"
#include <iostream>
#include <string>
#include <cmath>


void ErrorStats::reset() {
    num = 0;
    maxidx = -1;
    sum = sumsq = maxerr =0.0;
}

void ErrorStats::add(const double &val) {
    ++num;
    if (val > maxerr) {
        maxidx = num;
        maxerr = val;
    }
    sum += val;
    sumsq += val*val;
}

double ErrorStats::avg() const {
    return (num > 0) ? sum/num : 0.0;
}

double ErrorStats::dev() const {
    return (num > 0) ? sqrt(sumsq/num - sum/num*sum/num) : 0.0;
}

std::ostream &operator<<(std::ostream &out, const ErrorStats &stats)
{
    const std::ios_base::fmtflags flags = out.flags();
    const std::streamsize width = out.width(10);
    const std::streamsize prec = out.precision(3);

    out << std::scientific
        << "Average: " << stats.avg()
        << " StdDev: " << stats.dev()
        << " MaxErr: " << stats.max();

    out.precision(prec);
    out.width(width);
    out.flags(flags);

    return out << " @ item: " << stats.idx();
}

