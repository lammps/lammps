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
#include "utils.h"

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

Tokenizer::const_iterator Tokenizer::cbegin() const {
    return tokens.cbegin();
}

Tokenizer::const_iterator Tokenizer::cend() const {
    return tokens.cend();
}

std::string & Tokenizer::operator[](size_t index) {
    return tokens[index];
}

size_t Tokenizer::count() const {
    return tokens.size();
}


ValueTokenizer::ValueTokenizer(const std::string & str, const std::string & seperators) : tokens(str, seperators) {
    current  = tokens.begin();
}

bool ValueTokenizer::has_next() const {
    return current != tokens.cend();
}

std::string ValueTokenizer::next_string() {
    if (has_next()) {
        std::string value = *current;
        ++current;
        return value;
    }
    return "";
}

int ValueTokenizer::next_int() {
    if (has_next()) {
        if(!utils::is_integer(*current)) {
            throw InvalidIntegerException(*current);
        }
        int value = atoi(current->c_str());
        ++current;
        return value;
    }
    return 0;
}

bigint ValueTokenizer::next_bigint() {
    if (has_next()) {
        if(!utils::is_integer(*current)) {
            throw InvalidIntegerException(*current);
        }
        bigint value = ATOBIGINT(current->c_str());
        ++current;
        return value;
    }
    return 0;
}

tagint ValueTokenizer::next_tagint() {
    if (current != tokens.end()) {
        if(!utils::is_integer(*current)) {
            throw InvalidIntegerException(*current);
        }
        tagint value = ATOTAGINT(current->c_str());
        ++current;
        return value;
    }
    return 0;
}

double ValueTokenizer::next_double() {
    if (current != tokens.end()) {
        if(!utils::is_double(*current)) {
            throw InvalidFloatException(*current);
        }

        double value = atof(current->c_str());
        ++current;
        return value;
    }
    return 0.0;
}

void ValueTokenizer::skip(int ntokens) {
    current = std::next(current, ntokens);
}

size_t ValueTokenizer::count() const {
    return tokens.count();
}
