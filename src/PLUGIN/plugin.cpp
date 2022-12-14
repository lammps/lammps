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

#include "plugin.h"

#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "input.h"
#include "modify.h"

#include <cstring>
#include <list>
#include <map>

namespace LAMMPS_NS {
// list of plugin information data for loaded styles
static std::list<lammpsplugin_t> pluginlist;

// map for counting references to dso handles
static std::map<void *, int> dso_refcounter;

static bool verbose = true;

/* ---------------------------------------------------------------------- */

Plugin::Plugin(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void Plugin::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal plugin command");

  std::string cmd = arg[0];
  if (cmd == "load") {
    if (narg < 2) error->all(FLERR, "Illegal plugin load command");
    for (int i = 1; i < narg; ++i) plugin_load(arg[i], lmp);

  } else if (cmd == "unload") {
    if (narg != 3) error->all(FLERR, "Illegal plugin unload command");
    plugin_unload(arg[1], arg[2], lmp);

  } else if (cmd == "clear") {
    plugin_clear(lmp);

  } else if (cmd == "list") {
    if (comm->me == 0) {
      int num = plugin_get_num_plugins();
      utils::logmesg(lmp, "Currently loaded plugins\n");
      for (int i = 0; i < num; ++i) {
        auto entry = plugin_get_info(i);
        utils::logmesg(lmp, "{:4}: {} style plugin {}\n", i + 1, entry->style, entry->name);
      }
    }
  } else
    error->all(FLERR, "Illegal plugin command");
}

// auto-load DSOs from designated folder(s)
void plugin_auto_load(LAMMPS *lmp)
{
#if defined(LMP_PLUGIN)
  for (const auto &plugin_dir : platform::list_pathenv("LAMMPS_PLUGIN_PATH")) {
    verbose = false;
    int count = 0;
    for (const auto &file : platform::list_directory(plugin_dir)) {
      if (utils::strmatch(file, "\\plugin.so$"))
        count += plugin_load(platform::path_join(plugin_dir, file).c_str(), lmp);
    }
    if (lmp->comm->me == 0) utils::logmesg(lmp, "Loaded {} plugins from {}\n", count, plugin_dir);
  }
#endif
}

// load DSO and call included registration function
int plugin_load(const char *file, LAMMPS *lmp)
{
#if defined(LMP_PLUGIN)
  int me = lmp->comm->me;

  // open DSO file from given path; load symbols globally

  platform::dlerror();
  void *dso = platform::dlopen(file);
  if (dso == nullptr) {
    if (me == 0) utils::logmesg(lmp, "Open of file {} failed: {}\n", file, platform::dlerror());
    return 0;
  }

  // look up lammpsplugin_init() function in DSO
  // function must have C bindings so there is no name mangling

  platform::dlerror();
  void *initfunc = platform::dlsym(dso, "lammpsplugin_init");
  if (initfunc == nullptr) {
    platform::dlclose(dso);

    if (me == 0)
      utils::logmesg(lmp, "Plugin symbol lookup failure in file {}: {}\n", file,
                     platform::dlerror());
    return 0;
  }

  // call initializer function loaded from DSO and pass a pointer
  // to the LAMMPS instance, the DSO handle (for reference counting)
  // and plugin registration function pointer

  (*(lammpsplugin_initfunc) (initfunc))((void *) lmp, dso, (void *) &plugin_register);
  return 1;
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
#if defined(LMP_PLUGIN)
  auto lmp = (LAMMPS *) ptr;
  int me = lmp->comm->me;

  if (plugin == nullptr) return;

  // ignore load request if same plugin already loaded
  int idx = plugin_find(plugin->style, plugin->name);
  if (idx >= 0) {
    if (verbose && (me == 0))
      utils::logmesg(lmp, "Ignoring load of {} style {}: must unload existing {} plugin first\n",
                     plugin->style, plugin->name, plugin->name);
    return;
  }

  if (verbose && (me == 0)) {
    utils::logmesg(lmp, "Loading plugin: {} by {}\n", plugin->info, plugin->author);
    // print version info only if the versions of host and plugin don't match
    if ((plugin->version) && (strcmp(plugin->version, lmp->version) != 0))
      utils::logmesg(lmp, "  compiled for LAMMPS version {}, loaded into LAMMPS version {}\n",
                     plugin->version, lmp->version);
  }

  pluginlist.push_back(*plugin);

  if (dso_refcounter.find(plugin->handle) != dso_refcounter.end()) {
    ++dso_refcounter[plugin->handle];
  } else {
    dso_refcounter[plugin->handle] = 1;
  }

  std::string pstyle = plugin->style;
  if (pstyle == "pair") {
    auto pair_map = lmp->force->pair_map;
    if (pair_map->find(plugin->name) != pair_map->end()) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in pair style {} from plugin", plugin->name);
    }
    (*pair_map)[plugin->name] = (Force::PairCreator) plugin->creator.v1;

  } else if (pstyle == "bond") {
    auto bond_map = lmp->force->bond_map;
    if (bond_map->find(plugin->name) != bond_map->end()) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in bond style {} from plugin", plugin->name);
    }
    (*bond_map)[plugin->name] = (Force::BondCreator) plugin->creator.v1;

  } else if (pstyle == "angle") {
    auto angle_map = lmp->force->angle_map;
    if (angle_map->find(plugin->name) != angle_map->end()) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in angle style {} from plugin", plugin->name);
    }
    (*angle_map)[plugin->name] = (Force::AngleCreator) plugin->creator.v1;

  } else if (pstyle == "dihedral") {
    auto dihedral_map = lmp->force->dihedral_map;
    if (dihedral_map->find(plugin->name) != dihedral_map->end()) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in dihedral style {} from plugin",
                            plugin->name);
    }
    (*dihedral_map)[plugin->name] = (Force::DihedralCreator) plugin->creator.v1;

  } else if (pstyle == "improper") {
    auto improper_map = lmp->force->improper_map;
    if (improper_map->find(plugin->name) != improper_map->end()) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in improper style {} from plugin",
                            plugin->name);
    }
    (*improper_map)[plugin->name] = (Force::ImproperCreator) plugin->creator.v1;

  } else if (pstyle == "kspace") {
    auto kspace_map = lmp->force->kspace_map;
    if (kspace_map->find(plugin->name) != kspace_map->end()) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in kspace style {} from plugin", plugin->name);
    }
    (*kspace_map)[plugin->name] = (Force::KSpaceCreator) plugin->creator.v1;

  } else if (pstyle == "compute") {
    auto compute_map = lmp->modify->compute_map;
    if (compute_map->find(plugin->name) != compute_map->end()) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in compute style {} from plugin",
                            plugin->name);
    }
    (*compute_map)[plugin->name] = (Modify::ComputeCreator) plugin->creator.v2;

  } else if (pstyle == "fix") {
    auto fix_map = lmp->modify->fix_map;
    if (fix_map->find(plugin->name) != fix_map->end()) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in fix style {} from plugin", plugin->name);
    }
    (*fix_map)[plugin->name] = (Modify::FixCreator) plugin->creator.v2;

  } else if (pstyle == "region") {
    auto region_map = lmp->domain->region_map;
    if (region_map->find(plugin->name) != region_map->end()) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in region style {} from plugin", plugin->name);
    }
    (*region_map)[plugin->name] = (Domain::RegionCreator) plugin->creator.v2;

  } else if (pstyle == "command") {
    auto command_map = lmp->input->command_map;
    if (command_map->find(plugin->name) != command_map->end()) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in command style {} from plugin",
                            plugin->name);
    }
    (*command_map)[plugin->name] = (Input::CommandCreator) plugin->creator.v1;

  } else {
    utils::logmesg(lmp, "Loading plugins for {} styles not yet implemented\n", pstyle);
    pluginlist.pop_back();
  }
#endif
}

/* --------------------------------------------------------------------
     remove plugin from given style table and plugin list
     optionally close the DSO handle if it is the last plugin from that DSO
     must also delete style instances if style is currently active
     -------------------------------------------------------------------- */

void plugin_unload(const char *style, const char *name, LAMMPS *lmp)
{
#if defined(LMP_PLUGIN)
  int me = lmp->comm->me;

  // ignore unload request from unsupported style categories
  if ((strcmp(style, "pair") != 0) && (strcmp(style, "bond") != 0) &&
      (strcmp(style, "angle") != 0) && (strcmp(style, "dihedral") != 0) &&
      (strcmp(style, "improper") != 0) && (strcmp(style, "kspace") != 0) &&
      (strcmp(style, "compute") != 0) && (strcmp(style, "fix") != 0) &&
      (strcmp(style, "region") != 0) && (strcmp(style, "command") != 0)) {
    if (me == 0)
      utils::logmesg(lmp, "Ignoring unload: {} is not a supported plugin style\n", style);
    return;
  }

  // ignore unload request if not loaded from a plugin
  int idx = plugin_find(style, name);
  if (idx < 0) {
    if (me == 0)
      utils::logmesg(lmp, "Ignoring unload of {} style {}: not from a plugin\n", style, name);
    return;
  }

  // make copy of DSO handle for later use
  void *handle = plugin_get_info(idx)->handle;

  // remove selected plugin from list of plugins

  if (verbose && (me == 0)) utils::logmesg(lmp, "Unloading {} style {}\n", style, name);
  plugin_erase(style, name);

  // remove style of given name from corresponding map
  // must delete style instance if currently active so
  // we can close the DSO handle if the last reference is gone.

  std::string pstyle = style;
  if (pstyle == "pair") {

    auto found = lmp->force->pair_map->find(name);
    if (found != lmp->force->pair_map->end()) lmp->force->pair_map->erase(found);

    // must delete pair style instance if in use

    if (lmp->force->pair_style) {
      if (utils::strmatch(lmp->force->pair_style, "^hybrid")) {
        if (lmp->force->pair_match(name, 1, 1) != nullptr) lmp->force->create_pair("none", 0);
      } else {
        if (strcmp(lmp->force->pair_style, name) == 0) lmp->force->create_pair("none", 0);
      }
    }

  } else if (pstyle == "bond") {

    auto found = lmp->force->bond_map->find(name);
    if (found != lmp->force->bond_map->end()) lmp->force->bond_map->erase(found);

    // must delete bond style instance if in use

    if ((lmp->force->bond_style != nullptr) && (lmp->force->bond_match(name) != nullptr))
      lmp->force->create_bond("none", 0);

  } else if (pstyle == "angle") {

    auto found = lmp->force->angle_map->find(name);
    if (found != lmp->force->angle_map->end()) lmp->force->angle_map->erase(found);

    // must delete angle style instance if in use

    if ((lmp->force->angle_style != nullptr) && (lmp->force->angle_match(name) != nullptr))
      lmp->force->create_angle("none", 0);

  } else if (pstyle == "dihedral") {

    auto found = lmp->force->dihedral_map->find(name);
    if (found != lmp->force->dihedral_map->end()) lmp->force->dihedral_map->erase(found);

    // must delete dihedral style instance if in use

    if ((lmp->force->dihedral_style) && (lmp->force->dihedral_match(name) != nullptr))
      lmp->force->create_dihedral("none", 0);

  } else if (pstyle == "improper") {

    auto found = lmp->force->improper_map->find(name);
    if (found != lmp->force->improper_map->end()) lmp->force->improper_map->erase(found);

    // must delete improper style instance if in use

    if ((lmp->force->improper_style != nullptr) && (lmp->force->improper_match(name) != nullptr))
      lmp->force->create_improper("none", 0);

  } else if (pstyle == "kspace") {

    auto kspace_map = lmp->force->kspace_map;
    auto found = kspace_map->find(name);
    if (found != kspace_map->end()) kspace_map->erase(name);

  } else if (pstyle == "compute") {

    auto compute_map = lmp->modify->compute_map;
    auto found = compute_map->find(name);
    if (found != compute_map->end()) compute_map->erase(name);

    // must delete all compute instances using this compute style

    for (auto &icompute : lmp->modify->get_compute_by_style(name))
      lmp->modify->delete_compute(icompute->id);

  } else if (pstyle == "fix") {

    auto fix_map = lmp->modify->fix_map;
    auto found = fix_map->find(name);
    if (found != fix_map->end()) fix_map->erase(name);

    // must delete all fix instances using this fix style

    for (auto &ifix : lmp->modify->get_fix_by_style(name)) lmp->modify->delete_fix(ifix->id);

  } else if (pstyle == "region") {

    auto region_map = lmp->domain->region_map;
    auto found = region_map->find(name);
    if (found != region_map->end()) region_map->erase(name);

    for (auto &iregion : lmp->domain->get_region_by_style(name))
      lmp->domain->delete_region(iregion);

  } else if (pstyle == "command") {

    auto command_map = lmp->input->command_map;
    auto found = command_map->find(name);
    if (found != command_map->end()) command_map->erase(name);
  }

  // if reference count is down to zero, close DSO handle.

  --dso_refcounter[handle];
  if (dso_refcounter[handle] == 0) { platform::dlclose(handle); }
#endif
}

/* --------------------------------------------------------------------
     unload all loaded plugins
     -------------------------------------------------------------------- */

void plugin_clear(LAMMPS *lmp)
{
  verbose = false;
  while (pluginlist.size() > 0) {
    auto p = pluginlist.begin();
    plugin_unload(p->style, p->name, lmp);
  }
  verbose = true;
}

/* --------------------------------------------------------------------
     remove plugin of given name and style from internal lists
     -------------------------------------------------------------------- */

void plugin_erase(const char *style, const char *name)
{
  for (auto p = pluginlist.begin(); p != pluginlist.end(); ++p) {
    if ((strcmp(style, p->style) == 0) && (strcmp(name, p->name) == 0)) {
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
  int i = 0;
  for (const auto &entry : pluginlist) {
    if ((strcmp(style, entry.style) == 0) && (strcmp(name, entry.name) == 0)) return i;
    ++i;
  }
  return -1;
}

/* --------------------------------------------------------------------
     get pointer to plugin initializer struct at position idx
     -------------------------------------------------------------------- */

const lammpsplugin_t *plugin_get_info(int idx)
{
  int i = 0;
  for (const auto &p : pluginlist) {
    if (i == idx) return &p;
    ++i;
  }
  return nullptr;
}
}    // namespace LAMMPS_NS
