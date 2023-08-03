/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "lammpswrapper.h"

#if defined(LAMMPS_GUI_USE_PLUGIN)
#include "liblammpsplugin.h"
#else
#include "library.h"
#endif

LammpsWrapper::LammpsWrapper() : lammps_handle(nullptr), plugin_handle(nullptr) {}

void LammpsWrapper::open(int narg, char **args)
{
    if (!lammps_handle) {
#if defined(LAMMPS_GUI_USE_PLUGIN)
        lammps_handle = ((liblammpsplugin_t *)plugin_handle)->open_no_mpi(narg, args, nullptr);
#else
        lammps_handle = lammps_open_no_mpi(narg, args, nullptr);
#endif
    }
}

int LammpsWrapper::extract_setting(const char *keyword)
{
    int val = 0;
    if (lammps_handle) {
#if defined(LAMMPS_GUI_USE_PLUGIN)
        val = ((liblammpsplugin_t *)plugin_handle)->extract_setting(lammps_handle, keyword);
#else
        val           = lammps_extract_setting(lammps_handle, keyword);
#endif
    }
    return val;
}

double LammpsWrapper::get_thermo(const char *keyword)
{
    double val = 0.0;
    if (lammps_handle) {
#if defined(LAMMPS_GUI_USE_PLUGIN)
        val = ((liblammpsplugin_t *)plugin_handle)->get_thermo(lammps_handle, keyword);
#else
        val           = lammps_get_thermo(lammps_handle, keyword);
#endif
    }
    return val;
}

bool LammpsWrapper::is_running()
{
    int val = 0;
    if (lammps_handle) {
#if defined(LAMMPS_GUI_USE_PLUGIN)
        val = ((liblammpsplugin_t *)plugin_handle)->is_running(lammps_handle);
#else
        val           = lammps_is_running(lammps_handle);
#endif
    }
    return val != 0;
}

void LammpsWrapper::command(const char *input)
{
    if (lammps_handle) {
#if defined(LAMMPS_GUI_USE_PLUGIN)
        ((liblammpsplugin_t *)plugin_handle)->command(lammps_handle, input);
#else
        lammps_command(lammps_handle, input);
#endif
    }
}

void LammpsWrapper::commands_string(const char *input)
{
    if (lammps_handle) {
#if defined(LAMMPS_GUI_USE_PLUGIN)
        ((liblammpsplugin_t *)plugin_handle)->commands_string(lammps_handle, input);
#else
        lammps_commands_string(lammps_handle, input);
#endif
    }
}

// may be called with null handle. returns global error then.
bool LammpsWrapper::has_error() const
{
#if defined(LAMMPS_GUI_USE_PLUGIN)
    return ((liblammpsplugin_t *)plugin_handle)->has_error(lammps_handle) != 0;
#else
    return lammps_has_error(lammps_handle) != 0;
#endif
}

int LammpsWrapper::get_last_error_message(char *buf, int buflen)
{
#if defined(LAMMPS_GUI_USE_PLUGIN)
    return ((liblammpsplugin_t *)plugin_handle)->get_last_error_message(lammps_handle, buf, buflen);
#else
    return lammps_get_last_error_message(lammps_handle, buf, buflen);
#endif
}

void LammpsWrapper::force_timeout()
{
#if defined(LAMMPS_GUI_USE_PLUGIN)
    ((liblammpsplugin_t *)plugin_handle)->force_timeout(lammps_handle);
#else
    lammps_force_timeout(lammps_handle);
#endif
}

void LammpsWrapper::close()
{
#if defined(LAMMPS_GUI_USE_PLUGIN)
    if (lammps_handle) ((liblammpsplugin_t *)plugin_handle)->close(lammps_handle);
#else
    if (lammps_handle) lammps_close(lammps_handle);
#endif
    lammps_handle = nullptr;
}

void LammpsWrapper::finalize()
{
#if defined(LAMMPS_GUI_USE_PLUGIN)
    if (lammps_handle) {
        liblammpsplugin_t *lammps = (liblammpsplugin_t *)plugin_handle;
        lammps->close(lammps_handle);
        lammps->mpi_finalize();
        lammps->kokkos_finalize();
        lammps->python_finalize();
    }
#else
    if (lammps_handle) {
        lammps_close(lammps_handle);
        lammps_mpi_finalize();
        lammps_kokkos_finalize();
        lammps_python_finalize();
    }
#endif
}

bool LammpsWrapper::config_has_package(const char *package) const
{
#if defined(LAMMPS_GUI_USE_PLUGIN)
    return ((liblammpsplugin_t *)plugin_handle)->config_has_package(package) != 0;
#else
    return lammps_config_has_package(package) != 0;
#endif
}

bool LammpsWrapper::config_accelerator(const char *package, const char *category,
                                       const char *setting) const
{
#if defined(LAMMPS_GUI_USE_PLUGIN)
    return ((liblammpsplugin_t *)plugin_handle)->config_accelerator(package, category, setting) !=
           0;
#else
    return lammps_config_accelerator(package, category, setting) != 0;
#endif
}

bool LammpsWrapper::has_gpu_device() const
{
#if defined(LAMMPS_GUI_USE_PLUGIN)
    return ((liblammpsplugin_t *)plugin_handle)->has_gpu_device() != 0;
#else
    return lammps_has_gpu_device() != 0;
#endif
}

#if defined(LAMMPS_GUI_USE_PLUGIN)
bool LammpsWrapper::has_plugin() const
{
    return true;
}

bool LammpsWrapper::load_lib(const char *libfile)
{
    if (plugin_handle) liblammpsplugin_release((liblammpsplugin_t *)plugin_handle);
    plugin_handle = liblammpsplugin_load(libfile);
    if (!plugin_handle) return false;
    if (((liblammpsplugin_t *)plugin_handle)->abiversion != LAMMPSPLUGIN_ABI_VERSION) {
        liblammpsplugin_release((liblammpsplugin_t *)plugin_handle);
        plugin_handle = nullptr;
        return false;
    }
    return true;
}
#else
bool LammpsWrapper::has_plugin() const
{
    return false;
}

bool LammpsWrapper::load_lib(const char *)
{
    return true;
}
#endif

// Local Variables:
// c-basic-offset: 4
// End:
