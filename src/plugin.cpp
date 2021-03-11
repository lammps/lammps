/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "plugin.h"

#include "error.h"
#include "force.h"
#include "lammps.h"
#include "modify.h"

#include <vector>

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace LAMMPS_NS
{
  static std::vector<lammpsplugin_t> pluginlist;

  void plugin_load(const char *file, void *ptr)
  {
    LAMMPS *lmp = (LAMMPS *)ptr;
#if defined(WIN32)
    utils::logmesg(lmp,"Loading of plugins not supported on Windows yet\n");
#else
    void *dso = dlopen(file,RTLD_NOW|RTLD_GLOBAL);
    if (dso == nullptr) {
      utils::logmesg(lmp,fmt::format("Loading of plugin from file {} failed: "
                                     "{}", file, utils::getsyserror()));
      return;
    }

    void *initfunc = dlsym(dso,"lammpsplugin_init");
    if (initfunc == nullptr) {
      dlclose(dso);
      utils::logmesg(lmp,fmt::format("Symbol lookup failure in file {}\n",file));
      return;
    }
    ((lammpsplugin_initfunc)(initfunc))(ptr,dso,&plugin_register);
//    dlclose(dso);
#endif
  }

  void plugin_register(lammpsplugin_t *plugin, void *ptr)
  {
    LAMMPS *lmp = (LAMMPS *)ptr;

    if (plugin == nullptr) return;
    if (plugin->version)
      utils::logmesg(lmp,fmt::format("Loading: {}\n  compiled for LAMMPS "
                                     "version {} loaded into LAMMPS "
                                     "version {}\n", plugin->info,
                                     plugin->version, lmp->version));

    pluginlist.push_back(*plugin);

    std::string pstyle = plugin->style;
    if (pstyle == "pair") {
      (*(lmp->force->pair_map))[plugin->name] =
        (LAMMPS_NS::Force::PairCreator)plugin->creator;
    } else {
      utils::logmesg(lmp,fmt::format("Loading plugin for {} styles not "
                                     "yet implemented\n", pstyle));
      pluginlist.pop_back();
    }
  }

  int plugin_get_num_plugins() 
  {
    return pluginlist.size();
  }

  const lammpsplugin_t *plugin_info(int idx)
  {
    if ((idx < 0) || idx >= pluginlist.size()) return nullptr;

    return &pluginlist[idx];
  }

  int plugin_find(const char *style, const char *name)
  {
    for (int i=0; i < pluginlist.size(); ++i) {
      if ((strcmp(style,pluginlist[i].style) == 0)
          && (strcmp(name,pluginlist[i].name) == 0))
        return i;
    }
    return -1;
  }
  
  void plugin_unload(const char *style, const char *name, void *ptr)
  {
    LAMMPS *lmp = (LAMMPS *)ptr;
    std::string pstyle = style;
    if (pstyle == "pair") {
      auto found = lmp->force->pair_map->find(name);
      if (found != lmp->force->pair_map->end()) {
        lmp->force->pair_map->erase(found);
      }
    }
  }
}
