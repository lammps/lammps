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

#include "tokenizer.h"

using namespace LAMMPS_NS;

Tokenizer::Tokenizer(const std::string & str, const std::string & seperators) {
    size_t end = -1;

    do {
        size_t start = str.find_first_not_of(seperators, end + 1);
        if(start == std::string::npos) break;

        end = str.find_first_of(seperators, start);

        if(end == std::string::npos) {
            tokens.push_back(str.substr(start));
        } else {
            tokens.push_back(str.substr(start, end-start));
        }
    } while(end != std::string::npos);
}

Tokenizer::iterator Tokenizer::begin() {
    return tokens.begin();
}

Tokenizer::iterator Tokenizer::end() {
    return tokens.end();
}

const std::string & Tokenizer::operator[](size_t index) {
    return tokens[index];
}

const size_t Tokenizer::count() const {
    return tokens.size();
}