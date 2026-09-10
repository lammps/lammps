/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef ERROR_STATS_H
#define ERROR_STATS_H

#include <iostream>

class ErrorStats {
public:
    friend std::ostream &operator<<(std::ostream &out, const ErrorStats &stats);

    ErrorStats() { reset(); }
    virtual ~ErrorStats() {}

    void reset();
    void add(const double &val);
    [[nodiscard]] double avg() const;
    [[nodiscard]] double dev() const;
    [[nodiscard]] double max() const { return maxerr; }
    [[nodiscard]] double idx() const { return maxidx; }
    [[nodiscard]] bool has_data() const { return num > 0; }

private:
    double sum, sumsq, maxerr;
    int num, maxidx;
};

extern std::ostream &operator<<(std::ostream &out, const ErrorStats &stats);

#endif
