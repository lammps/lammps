/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#ifndef LMP_TOKENIZER_H
#define LMP_TOKENIZER_H

#include <string>
#include <vector>

namespace LAMMPS_NS {

class Tokenizer {
    std::vector<std::string> tokens;
public:
    typedef std::vector<std::string>::iterator iterator;

    Tokenizer(const std::string & str, const std::string & seperators = " \t\r\n\f");

    iterator begin();
    iterator end();

    const std::string & operator[](size_t index);
    const size_t count() const;
};

}

#endif
