
#include "lammpsplugin.h"

#include "version.h"

#include <cstring>

#include "fix_nve2.h"

using namespace LAMMPS_NS;

static Fix *nve2creator(LAMMPS *lmp, int argc, char **argv)
{
  return new FixNVE2(lmp, argc, argv);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  plugin.version = LAMMPS_VERSION;
  plugin.style = "fix";
  plugin.name = "nve2";
  plugin.info = "NVE2 variant fix style v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v2 = (lammpsplugin_factory2 *) &nve2creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
