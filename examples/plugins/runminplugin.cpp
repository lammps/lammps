
#include "lammpsplugin.h"

#include "version.h"

#include <cstring>

#include "min_cg2.h"
#include "verlet2.h"

using namespace LAMMPS_NS;

static Min *min_cg2creator(LAMMPS *lmp)
{
  return new MinCG2(lmp);
}

static Integrate *verlet2creator(LAMMPS *lmp, int argc, char **argv)
{
  return new Verlet2(lmp, argc, argv);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  plugin.version = LAMMPS_VERSION;

  plugin.style = "min";
  plugin.name = "cg2";
  plugin.info = "CG2 minimize style v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &min_cg2creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);

  plugin.style = "run";
  plugin.name = "verlet2";
  plugin.info = "Verlet2 run style v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v2 = (lammpsplugin_factory2 *) &verlet2creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
