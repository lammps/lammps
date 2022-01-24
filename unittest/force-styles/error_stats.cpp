/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "error_stats.h"
#include "fmt/format.h"
#include <cmath>
#include <iostream>
#include <string>

void ErrorStats::reset()
{
    num    = 0;
    maxidx = -1;
    sum = sumsq = maxerr = 0.0;
}

void ErrorStats::add(const double &val)
{
    ++num;
    if (val > maxerr) {
        maxidx = num;
        maxerr = val;
    }
    sum += val;
    sumsq += val * val;
}

double ErrorStats::avg() const
{
    return (num > 0) ? sum / num : 0.0;
}

double ErrorStats::dev() const
{
    return (num > 0) ? sqrt(sumsq / num - sum / num * sum / num) : 0.0;
}

std::ostream &operator<<(std::ostream &out, const ErrorStats &stats)
{
    out << fmt::format("Average: {:10.3e} StdDev: {:10.3e} MaxErr: {:10.3e} @ item: {}",
                       stats.avg(), stats.dev(), stats.max(), stats.idx());
    return out;
}
