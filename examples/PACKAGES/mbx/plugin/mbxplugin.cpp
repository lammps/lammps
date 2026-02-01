
#include "lammpsplugin.h"
#include "version.h"

#include "pair_mbx.h"
#include "fix_mbx.h"

using namespace LAMMPS_NS;

static Fix *fix_mbx_creator(LAMMPS *lmp, int argc, char **argv)
{
  return new FixMBX(lmp, argc, argv);
}

static Pair *pair_mbx_creator(LAMMPS *lmp)
{
  return new PairMBX(lmp);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  // register MBX pair style
  plugin.version = LAMMPS_VERSION;
  plugin.style = "pair";
  plugin.name = "mbx";
  plugin.info = "MBX plugin pair style v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &pair_mbx_creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);

  // register MBX fix style
  plugin.version = LAMMPS_VERSION;
  plugin.style = "fix";
  plugin.name = "MBX";
  plugin.info = "MBX plugin fix style v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &fix_mbx_creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
