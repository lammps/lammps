/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#ifndef TESTCASE_COMPRESSED_DUMP__H
#define TESTCASE_COMPRESSED_DUMP__H

#include "../testing/core.h"
#include "../testing/systems/melt.h"
#include <string>

extern const char *COMPRESS_SUFFIX;
extern const char *COMPRESS_EXTENSION;
extern char *COMPRESS_EXECUTABLE;

class CompressedDumpTest : public MeltTest {
protected:
    std::string dump_style;
    std::string compression_style;

public:
    CompressedDumpTest(const std::string &dump_style) : MeltTest(), dump_style(dump_style)
    {
        compression_style = fmt::format("{}/{}", dump_style, COMPRESS_SUFFIX);
    }

    std::string text_dump_filename(std::string ident)
    {
        return fmt::format("dump_{}_text_{}", COMPRESS_SUFFIX, ident);
    }

    std::string compressed_dump_filename(std::string ident)
    {
        return fmt::format("dump_{}_compressed_{}.{}", COMPRESS_SUFFIX, ident, COMPRESS_EXTENSION);
    }

    std::string converted_dump_filename(std::string ident)
    {
        return fmt::format("dump_{}_compressed_{}", COMPRESS_SUFFIX, ident);
    }

    void enable_triclinic()
    {
        BEGIN_HIDE_OUTPUT();
        command("change_box all triclinic");
        END_HIDE_OUTPUT();
    }

    void generate_dump(std::string dump_file, std::string dump_modify_options, int ntimesteps)
    {
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("dump id all {} 1 {}", dump_style, dump_file));

        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id {}", dump_modify_options));
        }

        command(fmt::format("run {}", ntimesteps));
        END_HIDE_OUTPUT();
    }

    void generate_text_and_compressed_dump(std::string text_file, std::string compressed_file,
                                           std::string dump_options,
                                           std::string dump_modify_options, int ntimesteps)
    {
        generate_text_and_compressed_dump(text_file, compressed_file, dump_options, dump_options,
                                          dump_modify_options, dump_modify_options, ntimesteps);
    }

    void generate_text_and_compressed_dump(std::string text_file, std::string compressed_file,
                                           std::string text_options, std::string compressed_options,
                                           std::string text_modify_options,
                                           std::string compressed_modify_options, int ntimesteps)
    {
        BEGIN_HIDE_OUTPUT();
        command(fmt::format("dump id0 all {} 1 {} {}", dump_style, text_file, text_options));
        command(fmt::format("dump id1 all {} 1 {} {}", compression_style, compressed_file,
                            compressed_options));

        if (!text_modify_options.empty()) {
            command(fmt::format("dump_modify id0 {}", text_modify_options));
        }

        if (!compressed_modify_options.empty()) {
            command(fmt::format("dump_modify id1 {}", compressed_modify_options));
        }

        command(fmt::format("run {}", ntimesteps));
        END_HIDE_OUTPUT();
    }

    std::string convert_compressed_to_text(std::string compressed_file)
    {
        BEGIN_HIDE_OUTPUT();
        std::string converted_file = compressed_file.substr(0, compressed_file.find_last_of('.'));
        std::string cmdline =
            fmt::format("\"{}\" -d -c {} > {}", COMPRESS_EXECUTABLE, compressed_file, converted_file);
        system(cmdline.c_str());
        END_HIDE_OUTPUT();
        return converted_file;
    }
};

#endif
