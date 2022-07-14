
#include "lammpsplugin.h"
#include "version.h"

#include "pair_pace.h"

using namespace LAMMPS_NS;

static Pair *pair_pace_creator(LAMMPS *lmp)
{
  return new PairPACE(lmp);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  // register pace pair style
  plugin.version = LAMMPS_VERSION;
  plugin.style = "pair";
  plugin.name = "pace";
  plugin.info = "PACE plugin pair style v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &pair_pace_creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
