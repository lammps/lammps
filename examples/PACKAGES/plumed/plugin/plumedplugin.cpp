
#include "lammpsplugin.h"
#include "version.h"

#include "fix_plumed.h"

using namespace LAMMPS_NS;

static Fix *fix_plumed_creator(LAMMPS *lmp, int argc, char **argv)
{
  return new FixPlumed(lmp, argc, argv);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  // register plumed fix style
  plugin.version = LAMMPS_VERSION;
  plugin.style = "fix";
  plugin.name = "plumed";
  plugin.info = "Plumed2 plugin fix style v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &fix_plumed_creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
