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
#include "fmt/format.h"

#include <exception>
#include <utility>

using namespace LAMMPS_NS;

TokenizerException::TokenizerException(const std::string &msg, const std::string &token){
    if(token.empty()) {
        message = msg;
    } else {
        message = fmt::format("{}: '{}'", msg, token);
    }
}

/** Class for splitting text into words
 *
 * This tokenizer will break down a string into sub-strings (i.e words)
 * separated by the given separator characters.
 *
\verbatim embed:rst

*See also*
   :cpp:class:`ValueTokenizer`, :cpp:func:`utils::split_words`

\endverbatim
 *
 * \param str         string to be processed
 * \param separators  string with separator characters (default: " \t\r\n\f") */

Tokenizer::Tokenizer(const std::string &str, const std::string &separators) :
    text(str), separators(separators), start(0), ntokens(std::string::npos)
{
    reset();
}

Tokenizer::Tokenizer(const Tokenizer &rhs) :
    text(rhs.text), separators(rhs.separators), ntokens(rhs.ntokens)
{
    reset();
}

Tokenizer::Tokenizer(Tokenizer && rhs) :
    text(std::move(rhs.text)), separators(std::move(rhs.separators)), ntokens(rhs.ntokens)
{
    reset();
}

/*! Re-position the tokenizer state to the first word,
 * i.e. the first non-separator character */
void Tokenizer::reset() {
    start = text.find_first_not_of(separators);
}

/*! Search the text to be processed for a sub-string.
 *
 * \param  str  string to be searched for
 * \return      true if string was found, false if not */
bool Tokenizer::contains(const std::string &str) const {
    return text.find(str) != std::string::npos;
}

/*! Skip over a given number of tokens
 *
 * \param  n  number of tokens to skip over */
void Tokenizer::skip(int n) {
    for(int i = 0; i < n; ++i) {
        if(!has_next()) throw TokenizerException("No more tokens", "");

        size_t end = text.find_first_of(separators, start);

        if(end == std::string::npos) {
            start = end;
        } else {
            start = text.find_first_not_of(separators, end+1);
        }
    }
}

/*! Indicate whether more tokens are available
 *
 * \return   true if there are more tokens, false if not */
bool Tokenizer::has_next() const {
    return start != std::string::npos;
}

/*! Retrieve next token.
 *
 * \return   string with the next token */
std::string Tokenizer::next() {
    if(!has_next()) throw TokenizerException("No more tokens", "");

    size_t end = text.find_first_of(separators, start);

    if(end == std::string::npos) {
        std::string token = text.substr(start);
        start = end;
        return token;
    }

    std::string token = text.substr(start, end-start);
    start = text.find_first_not_of(separators, end+1);
    return token;
}

/*! Count number of tokens in text.
 *
 * \return   number of counted tokens */
size_t Tokenizer::count() {
    // lazy evaluation
    if (ntokens == std::string::npos) {
      ntokens = utils::count_words(text, separators);
    }
    return ntokens;
}

/*! Retrieve the entire text converted to an STL vector of tokens.
 *
 * \return   The STL vector */
std::vector<std::string> Tokenizer::as_vector() {
  // store current state
  size_t current = start;

  reset();

  // generate vector
  std::vector<std::string> tokens;

  while(has_next()) {
    tokens.emplace_back(next());
  }

  // restore state
  start = current;

  return tokens;
}

/*! Class for reading text with numbers
 *
\verbatim embed:rst

*See also*
   :cpp:class:`Tokenizer`

\endverbatim
 *
 * \param str         String to be processed
 * \param separators  String with separator characters (default: " \t\r\n\f")
 *
 * \see Tokenizer InvalidIntegerException InvalidFloatException */

ValueTokenizer::ValueTokenizer(const std::string &str, const std::string &separators) : tokens(str, separators) {
}

ValueTokenizer::ValueTokenizer(const ValueTokenizer &rhs) : tokens(rhs.tokens) {
}

ValueTokenizer::ValueTokenizer(ValueTokenizer &&rhs) : tokens(std::move(rhs.tokens)) {
}

/*! Indicate whether more tokens are available
 *
 * \return   true if there are more tokens, false if not */
bool ValueTokenizer::has_next() const {
    return tokens.has_next();
}

/*! Search the text to be processed for a sub-string.
 *
 * \param  value  string with value to be searched for
 * \return        true if string was found, false if not */
bool ValueTokenizer::contains(const std::string &value) const {
    return tokens.contains(value);
}

/*! Retrieve next token
 *
 * \return   string with next token */
std::string ValueTokenizer::next_string() {
    if (has_next()) {
        std::string value = tokens.next();
        return value;
    }
    return "";
}

/*! Retrieve next token and convert to int
 *
 * \return   value of next token */
int ValueTokenizer::next_int() {
    if (has_next()) {
        std::string current = tokens.next();
        if(!utils::is_integer(current)) {
            throw InvalidIntegerException(current);
        }
        int value = atoi(current.c_str());
        return value;
    }
    return 0;
}

/*! Retrieve next token and convert to bigint
 *
 * \return   value of next token */
bigint ValueTokenizer::next_bigint() {
    if (has_next()) {
        std::string current = tokens.next();
        if(!utils::is_integer(current)) {
            throw InvalidIntegerException(current);
        }
        bigint value = ATOBIGINT(current.c_str());
        return value;
    }
    return 0;
}

/*! Retrieve next token and convert to tagint
 *
 * \return   value of next token */
tagint ValueTokenizer::next_tagint() {
    if (has_next()) {
        std::string current = tokens.next();
        if(!utils::is_integer(current)) {
            throw InvalidIntegerException(current);
        }
        tagint value = ATOTAGINT(current.c_str());
        return value;
    }
    return 0;
}

/*! Retrieve next token and convert to double
 *
 * \return   value of next token */
double ValueTokenizer::next_double() {
    if (has_next()) {
        std::string current = tokens.next();
        if(!utils::is_double(current)) {
            throw InvalidFloatException(current);
        }
        double value = atof(current.c_str());
        return value;
    }
    return 0.0;
}

/*! Skip over a given number of tokens
 *
 * \param  n  number of tokens to skip over */
void ValueTokenizer::skip(int n) {
    tokens.skip(n);
}

/*! Count number of tokens in text.
 *
 * \return   number of counted tokens */
size_t ValueTokenizer::count() {
    return tokens.count();
}
