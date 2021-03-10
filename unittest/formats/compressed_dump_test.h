/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

const char * COMPRESS_SUFFIX = nullptr;
const char * COMPRESS_EXTENSION = nullptr;
char * COMPRESS_BINARY = nullptr;

class CompressedDumpTest : public MeltTest {
protected:
    std::string dump_style;
    std::string compression_style;

public:
    CompressedDumpTest(const std::string & dump_style) : MeltTest(), dump_style(dump_style) {
        compression_style = fmt::format("{}/{}", dump_style, COMPRESS_SUFFIX);
    }

    std::string text_dump_filename(std::string ident) {
        return fmt::format("dump_{}_text_{}", COMPRESS_SUFFIX, ident);
    }

    std::string compressed_dump_filename(std::string ident) {
        return fmt::format("dump_{}_compressed_{}.{}", COMPRESS_SUFFIX, ident, COMPRESS_EXTENSION);
    }

    std::string converted_dump_filename(std::string ident) {
        return fmt::format("dump_{}_compressed_{}", COMPRESS_SUFFIX, ident);
    }

    void enable_triclinic()
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        command("change_box all triclinic");
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void generate_dump(std::string dump_file, std::string dump_modify_options, int ntimesteps)
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        command(fmt::format("dump id all {} 1 {}", dump_style, dump_file));

        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id {}", dump_modify_options));
        }

        command(fmt::format("run {}", ntimesteps));
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void generate_text_and_compressed_dump(std::string text_file, std::string compressed_file,
                                           std::string dump_modify_options, int ntimesteps)
    {
        generate_text_and_compressed_dump(text_file, compressed_file,
                                          dump_modify_options, dump_modify_options, ntimesteps);
    }

    void generate_text_and_compressed_dump(std::string text_file, std::string compressed_file,
                                           std::string text_options, std::string compressed_options, int ntimesteps)
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        command(fmt::format("dump id0 all {} 1 {}", dump_style, text_file));
        command(fmt::format("dump id1 all {} 1 {}", compression_style, compressed_file));

        if (!text_options.empty()) {
            command(fmt::format("dump_modify id0 {}", text_options));
        }

        if (!compressed_options.empty()) {
            command(fmt::format("dump_modify id1 {}", compressed_options));
        }

        command(fmt::format("run {}", ntimesteps));
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    std::string convert_compressed_to_text(std::string compressed_file)
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        std::string converted_file = compressed_file.substr(0, compressed_file.find_last_of('.'));
        std::string cmdline =
            fmt::format("{} -d -c {} > {}", COMPRESS_BINARY, compressed_file, converted_file);
        system(cmdline.c_str());
        if (!verbose) ::testing::internal::GetCapturedStdout();
        return converted_file;
    }
};

#endif
