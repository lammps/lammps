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

#ifndef LMP_CREATOR_REGISTRY_H
#define LMP_CREATOR_REGISTRY_H

#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace LAMMPS_NS {

// Process-global registry mapping a style keyword to a factory function.
//
// One CreatorRegistry instance exists per style category (pair, fix, ...) as a
// function-local static (Meyers singleton), so it is shared by all LAMMPS
// instances in a process and persists across the "clear" command.  Each entry
// has two slots: a built-in creator (registered once from the generated
// style_*.cpp) and an optional plugin creator that overrides it.  Removing a
// plugin override transparently restores the built-in.
//
// All operations take a per-registry mutex via std::lock_guard so that
// concurrent LAMMPS instances loading plugins cannot race.  RAII guarantees the
// lock is released even if an exception (e.g. from Error::one) unwinds through
// a caller; the mutex is never held across a factory invocation.

template <typename Creator> class CreatorRegistry {
 public:
  // register a built-in style (called from generated style_*.cpp)
  void add_builtin(const std::string &name, Creator fn)
  {
    std::lock_guard<std::mutex> guard(mutex);
    entries[name].builtin = fn;
  }

  // register or override a style from a plugin
  void set_plugin(const std::string &name, Creator fn)
  {
    std::lock_guard<std::mutex> guard(mutex);
    entries[name].plugin = fn;
  }

  // drop a plugin override; the built-in (if any) becomes active again
  void clear_plugin(const std::string &name)
  {
    std::lock_guard<std::mutex> guard(mutex);
    auto it = entries.find(name);
    if (it == entries.end()) return;
    it->second.plugin = nullptr;
    if (it->second.builtin == nullptr) entries.erase(it);
  }

  // look up a creator; a plugin override takes precedence over the built-in.
  // returns nullptr if no style of that name is registered.
  Creator find(const std::string &name) const
  {
    std::lock_guard<std::mutex> guard(mutex);
    auto it = entries.find(name);
    if (it == entries.end()) return nullptr;
    return it->second.plugin ? it->second.plugin : it->second.builtin;
  }

  bool contains(const std::string &name) const
  {
    std::lock_guard<std::mutex> guard(mutex);
    return entries.find(name) != entries.end();
  }

  bool has_builtin(const std::string &name) const
  {
    std::lock_guard<std::mutex> guard(mutex);
    auto it = entries.find(name);
    return (it != entries.end()) && (it->second.builtin != nullptr);
  }

  bool has_plugin(const std::string &name) const
  {
    std::lock_guard<std::mutex> guard(mutex);
    auto it = entries.find(name);
    return (it != entries.end()) && (it->second.plugin != nullptr);
  }

  // list of all registered style keywords (for "info styles", help, ...)
  std::vector<std::string> keys() const
  {
    std::lock_guard<std::mutex> guard(mutex);
    std::vector<std::string> result;
    result.reserve(entries.size());
    for (const auto &item : entries) result.push_back(item.first);
    return result;
  }

 private:
  struct Entry {
    Creator builtin = nullptr;
    Creator plugin = nullptr;
  };
  std::unordered_map<std::string, Entry> entries;
  mutable std::mutex mutex;
};

// register all built-in styles into their global registries.  Idempotent and
// thread-safe (guarded internally by std::call_once); call early from the
// LAMMPS constructor so registration happens exactly once per process and the
// generated registration translation units are pulled into a static library.

void register_builtin_styles();

}    // namespace LAMMPS_NS

#endif
