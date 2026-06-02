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
#include "update.h"

#include <cstring>
#include <list>
#include <map>
#include <mutex>

namespace LAMMPS_NS {
namespace {
  // list of plugin information data for loaded styles
  std::list<lammpsplugin_t> pluginlist;

  // map for counting references to dso handles
  std::map<void *, int> dso_refcounter;

  // serializes all access to pluginlist / dso_refcounter so that concurrent
  // LAMMPS instances can load, unload, or list plugins at the same time.  Every
  // entry point that touches these containers (plugin_register, plugin_unload,
  // plugin_clear, plugin_finalize, and the "plugin list" command) takes this
  // lock; the helpers they call while holding it (plugin_find, plugin_erase,
  // plugin_get_info, plugin_get_num_plugins, plugin_unload_locked) assume it is
  // already held and must not re-acquire it.  RAII (std::lock_guard) guarantees
  // the lock is released even if an error handler throws out of the critical
  // section.
  std::mutex plugin_mutex;

  bool verbose = true;
}    // namespace

/* ---------------------------------------------------------------------- */

Plugin::Plugin(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void Plugin::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal plugin command");

  std::string cmd = arg[0];
  if (cmd == "load") {
    if (narg < 2) utils::missing_cmd_args(FLERR, "plugin load", error);
    for (int i = 1; i < narg; ++i) plugin_load(arg[i], lmp);

  } else if (cmd == "unload") {
    if (narg != 3) error->all(FLERR, Error::ARGZERO, "Plugin unload command requires exactly 2 arguments");
    plugin_unload(arg[1], arg[2], lmp);

  } else if (cmd == "clear") {
    plugin_clear(lmp);

  } else if (cmd == "restore") {
    // loaded plugins now live in the process-global style registry and persist
    // across the "clear" command, so there is nothing left to restore.
    if (comm->me == 0)
      utils::logmesg(lmp, "Loaded plugins persist across 'clear'; "
                          "'plugin restore' is no longer needed\n");

  } else if (cmd == "list") {
    if (comm->me == 0) {
      std::lock_guard<std::mutex> guard(plugin_mutex);
      int num = plugin_get_num_plugins();
      utils::logmesg(lmp, "Currently loaded plugins: {}\n", num);
      for (int i = 0; i < num; ++i) {
        const auto *entry = plugin_get_info(i);
        utils::logmesg(lmp, "{:4}: {} style plugin {}\n", i + 1, entry->style, entry->name);
      }
    }
  } else
    error->all(FLERR, Error::ARGZERO, "Unknown plugin command {}", cmd);
}

// auto-load DSOs from designated folder(s)
void plugin_auto_load(LAMMPS *lmp)
{
#if defined(LMP_PLUGIN)
  bool oldverbose = verbose;
  verbose = false;
  for (const auto &plugin_dir : platform::list_pathenv("LAMMPS_PLUGIN_PATH")) {
    int count = 0;
    for (const auto &file : platform::list_directory(plugin_dir)) {
      if (utils::strmatch(file, "\\plugin.so$"))
        count += plugin_load(platform::path_join(plugin_dir, file).c_str(), lmp);
    }
    if (lmp->comm->me == 0) utils::logmesg(lmp, "Loaded {} plugins from {}\n", count, plugin_dir);
  }
  verbose = oldverbose;
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

  (*(lammpsplugin_initfunc) initfunc)((void *) lmp, dso, (void *) &plugin_register);
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
  std::lock_guard<std::mutex> guard(plugin_mutex);
  auto *lmp = (LAMMPS *) ptr;
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
    if (Force::pair_styles().has_builtin(plugin->name)) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in pair style {} from plugin", plugin->name);
    }
    Force::pair_styles().set_plugin(plugin->name, (Force::PairCreator) plugin->creator.v1);

  } else if (pstyle == "bond") {
    if (Force::bond_styles().has_builtin(plugin->name)) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in bond style {} from plugin", plugin->name);
    }
    Force::bond_styles().set_plugin(plugin->name, (Force::BondCreator) plugin->creator.v1);

  } else if (pstyle == "angle") {
    if (Force::angle_styles().has_builtin(plugin->name)) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in angle style {} from plugin", plugin->name);
    }
    Force::angle_styles().set_plugin(plugin->name, (Force::AngleCreator) plugin->creator.v1);

  } else if (pstyle == "dihedral") {
    if (Force::dihedral_styles().has_builtin(plugin->name)) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in dihedral style {} from plugin",
                            plugin->name);
    }
    Force::dihedral_styles().set_plugin(plugin->name, (Force::DihedralCreator) plugin->creator.v1);

  } else if (pstyle == "improper") {
    if (Force::improper_styles().has_builtin(plugin->name)) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in improper style {} from plugin",
                            plugin->name);
    }
    Force::improper_styles().set_plugin(plugin->name, (Force::ImproperCreator) plugin->creator.v1);

  } else if (pstyle == "kspace") {
    if (Force::kspace_styles().has_builtin(plugin->name)) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in kspace style {} from plugin", plugin->name);
    }
    Force::kspace_styles().set_plugin(plugin->name, (Force::KSpaceCreator) plugin->creator.v1);

  } else if (pstyle == "compute") {
    if (Modify::compute_styles().has_builtin(plugin->name)) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in compute style {} from plugin",
                            plugin->name);
    }
    Modify::compute_styles().set_plugin(plugin->name, (Modify::ComputeCreator) plugin->creator.v2);

  } else if (pstyle == "fix") {
    if (Modify::fix_styles().has_builtin(plugin->name)) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in fix style {} from plugin", plugin->name);
    }
    Modify::fix_styles().set_plugin(plugin->name, (Modify::FixCreator) plugin->creator.v2);

  } else if (pstyle == "region") {
    if (Domain::region_styles().has_builtin(plugin->name)) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in region style {} from plugin", plugin->name);
    }
    Domain::region_styles().set_plugin(plugin->name, (Domain::RegionCreator) plugin->creator.v2);

  } else if (pstyle == "command") {
    if (Input::command_styles().has_builtin(plugin->name)) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in command style {} from plugin",
                            plugin->name);
    }
    Input::command_styles().set_plugin(plugin->name, (Input::CommandCreator) plugin->creator.v1);

  } else if (pstyle == "run") {
    if (Update::integrate_styles().has_builtin(plugin->name)) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in run style {} from plugin", plugin->name);
    }
    Update::integrate_styles().set_plugin(plugin->name, (Update::IntegrateCreator) plugin->creator.v2);

  } else if (pstyle == "min") {
    if (Update::minimize_styles().has_builtin(plugin->name)) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR, "Overriding built-in minimize style {} from plugin",
                            plugin->name);
    }
    Update::minimize_styles().set_plugin(plugin->name, (Update::MinimizeCreator) plugin->creator.v1);

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

#if defined(LMP_PLUGIN)
// unload a single plugin.  the caller must already hold plugin_mutex.

static void plugin_unload_locked(const char *style, const char *name, LAMMPS *lmp)
{
  int me = lmp->comm->me;

  // ignore unload request from unsupported style categories
  if ((strcmp(style, "pair") != 0) && (strcmp(style, "bond") != 0) &&
      (strcmp(style, "angle") != 0) && (strcmp(style, "dihedral") != 0) &&
      (strcmp(style, "improper") != 0) && (strcmp(style, "kspace") != 0) &&
      (strcmp(style, "compute") != 0) && (strcmp(style, "fix") != 0) &&
      (strcmp(style, "region") != 0) && (strcmp(style, "command") != 0) &&
      (strcmp(style, "run") != 0) && (strcmp(style, "min") != 0)) {
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

    // must delete pair style instance if in use

    if (lmp->force->pair_style) {
      if (utils::strmatch(lmp->force->pair_style, "^hybrid")) {
        if (lmp->force->pair_match(name, 1, 1) != nullptr) lmp->force->create_pair("none", 0);
      } else {
        if (strcmp(lmp->force->pair_style, name) == 0) lmp->force->create_pair("none", 0);
      }
    }

    Force::pair_styles().clear_plugin(name);

  } else if (pstyle == "bond") {

    // must delete bond style instance if in use

    if ((lmp->force->bond_style != nullptr) && (lmp->force->bond_match(name) != nullptr))
      lmp->force->create_bond("none", 0);

    Force::bond_styles().clear_plugin(name);

  } else if (pstyle == "angle") {

    // must delete angle style instance if in use

    if ((lmp->force->angle_style != nullptr) && (lmp->force->angle_match(name) != nullptr))
      lmp->force->create_angle("none", 0);

    Force::angle_styles().clear_plugin(name);

  } else if (pstyle == "dihedral") {

    // must delete dihedral style instance if in use

    if ((lmp->force->dihedral_style) && (lmp->force->dihedral_match(name) != nullptr))
      lmp->force->create_dihedral("none", 0);

    Force::dihedral_styles().clear_plugin(name);

  } else if (pstyle == "improper") {

    // must delete improper style instance if in use

    if ((lmp->force->improper_style != nullptr) && (lmp->force->improper_match(name) != nullptr))
      lmp->force->create_improper("none", 0);

    Force::improper_styles().clear_plugin(name);

  } else if (pstyle == "kspace") {

    // must delete kspace style instance if in use

    if ((lmp->force->kspace_style != nullptr) && (lmp->force->kspace_match(name, 1) != nullptr))
      lmp->force->create_kspace("none", 0);

    Force::kspace_styles().clear_plugin(name);

  } else if (pstyle == "compute") {

    // must delete all compute instances using this compute style

    for (const auto &icompute : lmp->modify->get_compute_by_style(name))
      lmp->modify->delete_compute(icompute->id);

    Modify::compute_styles().clear_plugin(name);

  } else if (pstyle == "fix") {

    // must delete all fix instances using this fix style

    for (const auto &ifix : lmp->modify->get_fix_by_style(name)) lmp->modify->delete_fix(ifix->id);

    Modify::fix_styles().clear_plugin(name);

  } else if (pstyle == "region") {

    // must delete all region instances using this region style

    for (const auto &iregion : lmp->domain->get_region_by_style(name))
      lmp->domain->delete_region(iregion);

    Domain::region_styles().clear_plugin(name);

  } else if (pstyle == "command") {

    Input::command_styles().clear_plugin(name);

  } else if (pstyle == "run") {

    // must restore default run style if plugin style is in use

    if (strcmp(name, lmp->update->integrate_style) == 0) {
      char *str = (char *) "verlet";
      lmp->update->create_integrate(1, &str, 1);
    }
    Update::integrate_styles().clear_plugin(name);

  } else if (pstyle == "min") {

    // must restore default minimize style if plugin style is in use

    if (strcmp(name, lmp->update->minimize_style) == 0) {
      char *str = (char *) "cg";
      lmp->update->create_minimize(1, &str, 1);
    }
    Update::minimize_styles().clear_plugin(name);
  }

  // if reference count is down to zero, close DSO handle.

  --dso_refcounter[handle];
  if (dso_refcounter[handle] == 0) { platform::dlclose(handle); }
}
#endif

void plugin_unload(const char *style, const char *name, LAMMPS *lmp)
{
#if defined(LMP_PLUGIN)
  std::lock_guard<std::mutex> guard(plugin_mutex);
  plugin_unload_locked(style, name, lmp);
#endif
}

/* --------------------------------------------------------------------
     unload all loaded plugins
     -------------------------------------------------------------------- */

void plugin_clear(LAMMPS *lmp)
{
#if defined(LMP_PLUGIN)
  std::lock_guard<std::mutex> guard(plugin_mutex);
  bool oldverbose = verbose;
  verbose = true;
  while (pluginlist.size() > 0) {
    auto p = pluginlist.begin();
    plugin_unload_locked(p->style, p->name, lmp);
  }
  verbose = oldverbose;
#endif
}

/* --------------------------------------------------------------------
     unload all shared objects
     -------------------------------------------------------------------- */

void plugin_finalize()
{
#if defined(LMP_PLUGIN)
  std::lock_guard<std::mutex> guard(plugin_mutex);
  while (pluginlist.size() > 0) {
    auto p = pluginlist.begin();

    void *handle = p->handle;
    plugin_erase(p->style, p->name);

    // if reference count is down to zero, close DSO handle.

    --dso_refcounter[handle];
    if (dso_refcounter[handle] == 0) { platform::dlclose(handle); }
  }
#endif
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
