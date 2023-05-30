/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(info,Info);
// clang-format on
#else

#ifndef LMP_INFO_H
#define LMP_INFO_H

#include "command.h"

#include <vector>

namespace LAMMPS_NS {

class Info : public Command {
 public:
  Info(class LAMMPS *lmp) : Command(lmp){};
  void command(int, char **) override;

  bool is_active(const char *, const char *);
  bool is_defined(const char *, const char *);
  bool is_available(const char *, const char *);

  bool has_style(const std::string &category, const std::string &name);
  std::vector<std::string> get_available_styles(const std::string &category);

  static bool has_gzip_support();
  static bool has_png_support();
  static bool has_jpeg_support();
  static bool has_ffmpeg_support();
  static bool has_fft_single_support();
  static bool has_exceptions();
  static bool has_package(const std::string &);
  static bool has_accelerator_feature(const std::string &, const std::string &,
                                      const std::string &);
  static bool has_gpu_device();
  static std::string get_gpu_device_info();
  static std::string get_accelerator_info(const std::string &pkg = "");

  void get_memory_info(double *);
  char **get_variable_names(int &num);

 private:
  void available_styles(FILE *out, int flags);

  void atom_styles(FILE *out);
  void integrate_styles(FILE *out);
  void minimize_styles(FILE *out);
  void pair_styles(FILE *out);
  void bond_styles(FILE *out);
  void angle_styles(FILE *out);
  void dihedral_styles(FILE *out);
  void improper_styles(FILE *out);
  void kspace_styles(FILE *out);
  void fix_styles(FILE *out);
  void compute_styles(FILE *out);
  void region_styles(FILE *out);
  void dump_styles(FILE *out);
  void command_styles(FILE *out);
};

}    // namespace LAMMPS_NS

#endif
#endif
