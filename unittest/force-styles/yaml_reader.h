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

#ifndef YAML_READER_H
#define YAML_READER_H

#include "yaml.h"

#include <cstdio>
#include <iostream>
#include <map>
#include <string>

template <typename ConsumerClass> class YamlReader {
private:
    enum StateValue {
        START,
        ACCEPT_KEY,
        ACCEPT_VALUE,
        STOP,
        ERROR,
    };

    StateValue state;
    bool accepted;
    std::string key;
    std::string basename;

protected:
    typedef void (ConsumerClass::*EventConsumer)(const yaml_event_t &event);
    std::map<std::string, EventConsumer> consumers;

public:
    YamlReader() {}
    virtual ~YamlReader() {}

    std::string get_basename() const { return basename; }

    int parse_file(const std::string &infile)
    {
        basename          = infile;
        std::size_t found = basename.rfind(".yaml");
        if (found > 0) basename = basename.substr(0, found);
        found = basename.find_last_of("/\\");
        if (found != std::string::npos) basename = basename.substr(found + 1);

        FILE *fp = fopen(infile.c_str(), "r");
        yaml_parser_t parser;
        yaml_event_t event;

        if (!fp) {
            std::cerr << "Cannot open yaml file '" << infile << "': " << strerror(errno)
                      << std::endl;
            return 1;
        }
        yaml_parser_initialize(&parser);
        yaml_parser_set_input_file(&parser, fp);
        state = START;
        do {
            if (!yaml_parser_parse(&parser, &event)) {
                state = STOP;
            }

            if (!consume_event(event)) {
                state = STOP;
            }

            if (accepted) {
                if (!consume_key_value(key, event)) {
                    std::cerr << "Ignoring unknown key/value pair: " << key << " = "
                              << event.data.scalar.value << std::endl;
                }
            }
            yaml_event_delete(&event);
        } while (state != STOP);

        yaml_parser_delete(&parser);
        fclose(fp);
        return 0;
    }

protected:
    bool consume_key_value(const std::string &key, const yaml_event_t &event)
    {
        auto it                 = consumers.find(key);
        ConsumerClass *consumer = dynamic_cast<ConsumerClass *>(this);

        if (consumer) {
            if (it != consumers.end()) {
                // std::cerr << "Loading: " << key << std::endl;
                (consumer->*(it->second))(event);
                return true;
            }
            std::cerr << "UNKNOWN" << std::endl;
        } else {
            std::cerr << "ConsumerClass is not valid" << std::endl;
        }
        return false;
    }

    bool consume_event(yaml_event_t &event)
    {
        accepted = false;
        switch (state) {
            case START:
                switch (event.type) {
                    case YAML_MAPPING_START_EVENT:
                        state = ACCEPT_KEY;
                        break;
                    case YAML_SCALAR_EVENT:
                    case YAML_SEQUENCE_START_EVENT:
                        state = ERROR;
                        break;
                    case YAML_STREAM_END_EVENT:
                        state = STOP;
                        break;
                    case YAML_STREAM_START_EVENT:
                    case YAML_DOCUMENT_START_EVENT:
                    case YAML_DOCUMENT_END_EVENT:
                        // ignore
                        break;
                    default:
                        std::cerr << "UNHANDLED YAML EVENT:  " << event.type << std::endl;
                        state = ERROR;
                        break;
                }
                break;

            case ACCEPT_KEY:
                switch (event.type) {
                    case YAML_SCALAR_EVENT:
                        key   = (char *)event.data.scalar.value;
                        state = ACCEPT_VALUE;
                        break;
                    case YAML_MAPPING_END_EVENT:
                        state = STOP;
                        break;
                    default:
                        std::cerr << "UNHANDLED YAML EVENT (key): " << event.type
                                  << "\nVALUE: " << event.data.scalar.value << std::endl;
                        state = ERROR;
                        break;
                }
                break;

            case ACCEPT_VALUE:
                switch (event.type) {
                    case YAML_SCALAR_EVENT:
                        accepted = true;
                        state    = ACCEPT_KEY;
                        break;
                    default:
                        std::cerr << "UNHANDLED YAML EVENT (value): " << event.type
                                  << "\nVALUE: " << event.data.scalar.value << std::endl;
                        state = ERROR;
                        break;
                }
                break;

            case ERROR:
            case STOP:
                break;
        }
        return (state != ERROR);
    }
};

#endif
