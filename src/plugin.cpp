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
#include "pair.h"

#include <map>
#include <list>

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace LAMMPS_NS
{
  // store list of plugin information data for loaded styles
  static std::list<lammpsplugin_t> pluginlist;
  // map of dso handles 
  static std::map<void *, int> dso_refcounter;

  // load DSO and call included registration function
  void plugin_load(const char *file, LAMMPS *lmp)
  {
#if defined(WIN32)
    utils::logmesg(lmp,"Loading of plugins not supported on Windows yet\n");
#else

    // open DSO from given path
    
    void *dso = dlopen(file,RTLD_NOW|RTLD_GLOBAL);
    if (dso == nullptr) {
      utils::logmesg(lmp,fmt::format("Loading of plugin from file {} failed: "
                                     "{}", file, utils::getsyserror()));
      return;
    }

    // look up lammpsplugin_init() function in DSO. must have C bindings.

    void *initfunc = dlsym(dso,"lammpsplugin_init");
    if (initfunc == nullptr) {
      dlclose(dso);
      utils::logmesg(lmp,fmt::format("Symbol lookup failure in file {}\n",file));
      return;
    }

    // call initializer function loaded from DSO and pass pointer to LAMMPS instance,
    // the DSO handle (for reference counting) and plugin registration function pointer

    ((lammpsplugin_initfunc)(initfunc))((void *)lmp, dso, (void *)&plugin_register);
#endif
  }

  // register new style from plugin with LAMMPS
  void plugin_register(lammpsplugin_t *plugin, void *ptr)
  {
    LAMMPS *lmp = (LAMMPS *)ptr;

    if (plugin == nullptr) return;
    utils::logmesg(lmp,fmt::format("Loading plugin: {} by {}\n",
                                   plugin->info, plugin->author));
    if ((plugin->version) && (strcmp(plugin->version,lmp->version) != 0))
      utils::logmesg(lmp,fmt::format("  compiled for LAMMPS version {} "
                                     "loaded into LAMMPS version {}\n",
                                     plugin->version, lmp->version));
    pluginlist.push_back(*plugin);
    if (dso_refcounter.find(plugin->handle) != dso_refcounter.end()) {
      ++ dso_refcounter[plugin->handle];
    } else {
      dso_refcounter[plugin->handle] = 1;
    }

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

  const lammpsplugin_t *plugin_get_info(int idx)
  {
    int i=0;
    for (auto p=pluginlist.begin(); p != pluginlist.end(); ++p) {
      if (i == idx) return &(*p);
      ++i;
    }
    return nullptr;
  }

  int plugin_find(const char *style, const char *name)
  {
    int i=0;
    for (auto entry : pluginlist) {
      if ((strcmp(style,entry.style) == 0)
          && (strcmp(name,entry.name) == 0))
        return i;
      ++i;
    }
    return -1;
  }

  void plugin_erase(const char *style, const char *name)
  {
    fmt::print("Erasing {} style {} from list with {} plugins\n",
               style, name, pluginlist.size());
    for (auto p=pluginlist.begin(); p != pluginlist.end(); ++p) {
      if ((strcmp(style,p->style) == 0)
          && (strcmp(name,p->name) == 0)) {
        pluginlist.erase(p);
        return;
      }
    }
  }

  // remove plugin from given style table and plugin list.
  // close dso handle if this is the last plugin from that dso.
  void plugin_unload(const char *style, const char *name, LAMMPS *lmp)
  {
    int idx = plugin_find(style,name);
    if (idx < 0)
      lmp->error->all(FLERR,fmt::format("{} style {} is not loaded from "
                                        "a plugin", style, name));
    auto plugin = plugin_get_info(idx);
    void *handle = plugin->handle;
    utils::logmesg(lmp,fmt::format("Unloading {} style {}\n",style,name));
    plugin_erase(style,name);

    std::string pstyle = style;
    if (pstyle == "pair") {
      auto found = lmp->force->pair_map->find(name);
      if (found != lmp->force->pair_map->end()) {
        lmp->force->pair_map->erase(found);
      }

      // must delete pair style instance if in use
      if (lmp->force->pair_style) {
        if (utils::strmatch(lmp->force->pair_style,"^hybrid")) {
          if (lmp->force->pair_match(name,1,1) != nullptr)
            lmp->force->create_pair("none",0);          
        } else {
          if (strcmp(lmp->force->pair_style,name) == 0)
            lmp->force->create_pair("none",0);
        }
      }
    }

    -- dso_refcounter[handle];
    if (dso_refcounter[handle] == 0) {
      dlclose(handle);
    }
  }
}
