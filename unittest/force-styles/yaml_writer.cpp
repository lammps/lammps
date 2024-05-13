/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "yaml_writer.h"
#include "fmt/format.h"
#include "yaml.h"

#include <cstdio>
#include <string>

YamlWriter::YamlWriter(const char *outfile)
{
    yaml_emitter_initialize(&emitter);
    fp = fopen(outfile, "w");
    if (!fp) {
        perror(__FILE__);
        return;
    }

    yaml_emitter_set_output_file(&emitter, fp);

    yaml_stream_start_event_initialize(&event, YAML_UTF8_ENCODING);
    yaml_emitter_emit(&emitter, &event);
    yaml_document_start_event_initialize(&event, nullptr, nullptr, nullptr, 0);
    yaml_emitter_emit(&emitter, &event);
    yaml_mapping_start_event_initialize(&event, nullptr, (yaml_char_t *)YAML_MAP_TAG, 1,
                                        YAML_ANY_MAPPING_STYLE);
    yaml_emitter_emit(&emitter, &event);
}

YamlWriter::~YamlWriter()
{
    yaml_mapping_end_event_initialize(&event);
    yaml_emitter_emit(&emitter, &event);
    yaml_document_end_event_initialize(&event, 0);
    yaml_emitter_emit(&emitter, &event);
    yaml_stream_end_event_initialize(&event);
    yaml_emitter_emit(&emitter, &event);
    yaml_emitter_delete(&emitter);
    fclose(fp);
}

void YamlWriter::emit(const std::string &key, const double value)
{
    emit(key, fmt::format("{}", value));
}

void YamlWriter::emit(const std::string &key, const long value)
{
    emit(key, fmt::format("{}", value));
}

void YamlWriter::emit(const std::string &key, const int value)
{
    emit(key, fmt::format("{}", value));
}

void YamlWriter::emit(const std::string &key, const std::string &value)
{
    yaml_scalar_event_initialize(&event, nullptr, (yaml_char_t *)YAML_STR_TAG,
                                 (yaml_char_t *)key.c_str(), key.size(), 1, 0,
                                 YAML_PLAIN_SCALAR_STYLE);
    yaml_emitter_emit(&emitter, &event);
    yaml_scalar_event_initialize(&event, nullptr, (yaml_char_t *)YAML_STR_TAG,
                                 (yaml_char_t *)value.c_str(), value.size(), 1, 0,
                                 YAML_PLAIN_SCALAR_STYLE);
    yaml_emitter_emit(&emitter, &event);
}

void YamlWriter::emit_block(const std::string &key, const std::string &value)
{
    yaml_scalar_event_initialize(&event, nullptr, (yaml_char_t *)YAML_STR_TAG,
                                 (yaml_char_t *)key.c_str(), key.size(), 1, 0,
                                 YAML_PLAIN_SCALAR_STYLE);
    yaml_emitter_emit(&emitter, &event);
    yaml_scalar_event_initialize(&event, nullptr, (yaml_char_t *)YAML_STR_TAG,
                                 (yaml_char_t *)value.c_str(), value.size(), 1, 0,
                                 YAML_LITERAL_SCALAR_STYLE);
    yaml_emitter_emit(&emitter, &event);
}
