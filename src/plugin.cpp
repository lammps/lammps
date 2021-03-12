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

#include "comm.h"
#include "error.h"
#include "input.h"
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
  // list of plugin information data for loaded styles
  static std::list<lammpsplugin_t> pluginlist;

  // map for counting references to dso handles
  static std::map<void *, int> dso_refcounter;

  // load DSO and call included registration function
  void plugin_load(const char *file, LAMMPS *lmp)
  {
    int me = lmp->comm->me;
#if defined(WIN32)
    lmp->error->all(FLERR,"Loading of plugins on Windows not yet supported\n");
#else

    // open DSO file from given path; load symbols globally

    void *dso = dlopen(file,RTLD_NOW|RTLD_GLOBAL);
    if (dso == nullptr) {
      if (me == 0)
        utils::logmesg(lmp,fmt::format("Open of file {} failed\n",file));
      return;
    }

    // look up lammpsplugin_init() function in DSO
    // function must have C bindings so there is no name mangling

    void *initfunc = dlsym(dso,"lammpsplugin_init");
    if (initfunc == nullptr) {
      dlclose(dso);

      if (me == 0)
        utils::logmesg(lmp,fmt::format("Plugin symbol lookup failure in "
                                       "file {}\n",file));
      return;
    }

    // call initializer function loaded from DSO and pass a pointer
    // to the LAMMPS instance, the DSO handle (for reference counting)
    // and plugin registration function pointer

    (*(lammpsplugin_initfunc)(initfunc))((void *)lmp, dso,
                                         (void *)&plugin_register);
#endif
  }

  /* --------------------------------------------------------------------
     register a new style from a plugin with LAMMPS
     this is the callback function that is called from within
     the plugin initializer function. all plugin information
     is taken from the lammpsplugin_t struct.
     -------------------------------------------------------------------- */

  void plugin_register(lammpsplugin_t *plugin, void *ptr)
  {
    LAMMPS *lmp = (LAMMPS *)ptr;
    int me = lmp->comm->me;

    if (plugin == nullptr) return;

    // ignore load request if same plugin already loaded
    int idx = plugin_find(plugin->style,plugin->name);
    if (idx >= 0) {
      if (me == 0)
        utils::logmesg(lmp,fmt::format("Ignoring load of {} style {}: must "
                                       "unload existing {} plugin first\n",
                                       plugin->style,plugin->name,plugin->name));
      return;
    }

    if (me == 0) {
      utils::logmesg(lmp,fmt::format("Loading plugin: {} by {}\n",
                                     plugin->info,plugin->author));
      // print version info only if the versions of host and plugin don't match
      if ((plugin->version) && (strcmp(plugin->version,lmp->version) != 0))
        utils::logmesg(lmp,fmt::format("  compiled for LAMMPS version {} "
                                       "loaded into LAMMPS version {}\n",
                                       plugin->version,lmp->version));
    }

    pluginlist.push_back(*plugin);

    if (dso_refcounter.find(plugin->handle) != dso_refcounter.end()) {
      ++ dso_refcounter[plugin->handle];
    } else {
      dso_refcounter[plugin->handle] = 1;
    }

    std::string pstyle = plugin->style;
    if (pstyle == "pair") {
      auto pair_map = lmp->force->pair_map;
      if (pair_map->find(plugin->name) != pair_map->end()) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,fmt::format("Overriding built-in pair "
                                                "style {} from plugin",
                                                plugin->name));
      }
      (*pair_map)[plugin->name] = (Force::PairCreator)plugin->creator.v1;

    } else if (pstyle == "fix") {
      auto fix_map = lmp->modify->fix_map;
      if (fix_map->find(plugin->name) != fix_map->end()) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,fmt::format("Overriding built-in fix "
                                                "style {} from plugin",
                                                plugin->name));
      }
      (*fix_map)[plugin->name] = (Modify::FixCreator)plugin->creator.v2;

    } else if (pstyle == "command") {
      auto command_map = lmp->input->command_map;
      if (command_map->find(plugin->name) != command_map->end()) {
        if (lmp->comm->me == 0)
          lmp->error->warning(FLERR,fmt::format("Overriding built-in command "
                                                "style {} from plugin",
                                                plugin->name));
      }
      (*command_map)[plugin->name] = (Input::CommandCreator)plugin->creator.v3;

    } else {
      utils::logmesg(lmp,fmt::format("Loading plugin for {} styles not "
                                     "yet implemented\n",pstyle));
      pluginlist.pop_back();
    }
  }

  /* --------------------------------------------------------------------
     remove plugin from given style table and plugin list
     optionally close the DSO handle if it is the last plugin from that DSO
     must also delete style instances if style is currently active
     -------------------------------------------------------------------- */

  void plugin_unload(const char *style, const char *name, LAMMPS *lmp)
  {
    int me = lmp->comm->me;

    // ignore unload request if not loaded from a plugin
    int idx = plugin_find(style,name);
    if (idx < 0) {
      if (me == 0)
        utils::logmesg(lmp,fmt::format("Ignoring unload of {} style {}: not "
                                       "loaded from a plugin\n",style,name));
      return;
    }

    // copy of DSO handle for later
    void *handle = plugin_get_info(idx)->handle;

    // remove selected plugin from list of plugins

    if (me == 0)
      utils::logmesg(lmp,fmt::format("Unloading {} style {}\n",style,name));
    plugin_erase(style,name);

    // remove style of given name from corresponding map
    // must delete style instance if currently active so
    // we can close the DSO handle if the last reference is gone.

    std::string pstyle = style;
    if (pstyle == "pair") {
      auto found = lmp->force->pair_map->find(name);
      if (found != lmp->force->pair_map->end())
        lmp->force->pair_map->erase(found);

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

    } else if (pstyle == "fix") {
      for (int ifix = lmp->modify->find_fix_by_style(name);
           ifix >= 0; ifix = lmp->modify->find_fix_by_style(name))
        lmp->modify->delete_fix(ifix);

    } else if (pstyle == "command") {
      auto command_map = lmp->input->command_map;
      auto cmd = command_map->find(name);
      if (cmd != command_map->end()) command_map->erase(name);
    }

    // if reference count is down to zero, close DSO handle.

    -- dso_refcounter[handle];
    if (dso_refcounter[handle] == 0) {
#ifndef WIN32
      dlclose(handle);
#endif
    }
  }

  /* --------------------------------------------------------------------
     remove plugin of given name and style from internal lists
     -------------------------------------------------------------------- */

  void plugin_erase(const char *style, const char *name)
  {
    for (auto p=pluginlist.begin(); p != pluginlist.end(); ++p) {
      if ((strcmp(style,p->style) == 0)
          && (strcmp(name,p->name) == 0)) {
        pluginlist.erase(p);
        return;
      }
    }
  }

  /* --------------------------------------------------------------------
     number of styles loaded from plugin files
     -------------------------------------------------------------------- */

  int plugin_get_num_plugins()
  {
    return pluginlist.size();
  }

  /* --------------------------------------------------------------------
     return position index in list of given plugin of given style
     -------------------------------------------------------------------- */

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

  /* --------------------------------------------------------------------
     get pointer to plugin initializer struct at position idx
     -------------------------------------------------------------------- */

  const lammpsplugin_t *plugin_get_info(int idx)
  {
    int i=0;
    for (auto p=pluginlist.begin(); p != pluginlist.end(); ++p) {
      if (i == idx) return &(*p);
      ++i;
    }
    return nullptr;
  }
}
