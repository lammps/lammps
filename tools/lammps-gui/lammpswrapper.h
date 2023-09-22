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

#ifndef LAMMPSWRAPPER_H
#define LAMMPSWRAPPER_H

class LammpsWrapper {
public:
    LammpsWrapper();

public:
    void open(int nargs, char **args);
    void close();
    void finalize();

    void file(const char *);
    void command(const char *);
    void commands_string(const char *);

    void force_timeout();

    int extract_setting(const char *keyword);
    void *extract_global(const char *keyword);
    void *extract_atom(const char *keyword);

    int id_count(const char *idtype);
    int id_name(const char *idtype, int idx, char *buf, int buflen);
    int style_count(const char *keyword);
    int style_name(const char *keyword, int idx, char *buf, int buflen);
    int variable_info(int idx, char *buf, int buflen);

    double get_thermo(const char *keyword);
    void *last_thermo(const char *keyword, int idx);

    bool is_open() const { return lammps_handle != nullptr; }
    bool is_running();

    bool has_error() const;
    int get_last_error_message(char *errorbuf, int buflen);

    bool config_accelerator(const char *package, const char *category, const char *setting) const;
    bool config_has_package(const char *pkg) const;
    bool has_gpu_device() const;

    bool load_lib(const char *lammpslib);
    bool has_plugin() const;

private:
    void *lammps_handle;
    void *plugin_handle;
};
#endif

// Local Variables:
// c-basic-offset: 4
// End:
