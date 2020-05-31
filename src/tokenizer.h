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
#include "lmptype.h"
#include <exception>

namespace LAMMPS_NS {

class Tokenizer {
    std::vector<std::string> tokens;
public:
    typedef std::vector<std::string>::iterator iterator;
    typedef std::vector<std::string>::const_iterator const_iterator;

    Tokenizer(const std::string & str, const std::string & seperators = " \t\r\n\f");

    iterator begin();
    iterator end();
    const_iterator cbegin() const;
    const_iterator cend() const;

    std::string & operator[](size_t index);
    size_t count() const;
};

class TokenizerException : public std::exception {
  std::string message;
public:
  TokenizerException(const std::string & msg, const std::string & token) : message(msg + ": '" + token + "'") {
  }

  ~TokenizerException() throw() {
  }

  virtual const char * what() const throw() {
    return message.c_str();
  }
};

class InvalidIntegerException : public TokenizerException {
public:
    InvalidIntegerException(const std::string & token) : TokenizerException("Not a valid integer number", token) {
    }
};

class InvalidFloatException : public TokenizerException {
public:
    InvalidFloatException(const std::string & token) : TokenizerException("Not a valid floating-point number", token) {
    }
};

class ValueTokenizer {
    Tokenizer tokens;
    Tokenizer::const_iterator current;
public:
    ValueTokenizer(const std::string & str, const std::string & seperators = " \t\r\n\f");

    std::string next_string();
    tagint next_tagint();
    bigint next_bigint();
    int    next_int();
    double next_double();

    bool has_next() const;
    void skip(int ntokens);

    size_t count() const;
};


}

#endif
