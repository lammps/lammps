/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef YAML_WRITER_H
#define YAML_WRITER_H

#include "yaml.h"
#include <cstdio>
#include <string>

class YamlWriter {
public:
    YamlWriter(const char *outfile);
    virtual ~YamlWriter();

    // emitters
    void emit(const std::string &key, const double value);
    void emit(const std::string &key, const long value);
    void emit(const std::string &key, const int value);
    void emit(const std::string &key, const std::string &value);
    void emit_block(const std::string &key, const std::string &value);

private:
    FILE *fp;
    yaml_emitter_t emitter;
    yaml_event_t event;

private:
    YamlWriter(){};
    YamlWriter(const YamlWriter &){};
};

#endif
