
#include "lammpsplugin.h"
#include "version.h"

#include "pair_kim.h"
#include "fix_store_kim.h"
#include "kim_command.h"

using namespace LAMMPS_NS;

static Pair *pair_kim_creator(LAMMPS *lmp)
{
  return new PairKIM(lmp);
}

static Fix *fix_store_kim_creator(LAMMPS *lmp, int argc, char **argv)
{
  return new FixStoreKIM(lmp, argc, argv);
}

static Command *kim_command_creator(LAMMPS *lmp)
{
  return new KimCommand(lmp);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  // register kim pair style
  plugin.version = LAMMPS_VERSION;
  plugin.style = "pair";
  plugin.name = "kim";
  plugin.info = "KIM plugin pair style v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &pair_kim_creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);

  // register fix STORE/KIM only need to update changed fields
  plugin.style = "fix";
  plugin.name = "STORE/KIM";
  plugin.info = "Internal settings storage for KIM fix style v1.0";
  plugin.creator.v2 = (lammpsplugin_factory2 *) &fix_store_kim_creator;
  (*register_plugin)(&plugin, lmp);

  // register KIM command
  plugin.style = "command";
  plugin.name = "kim";
  plugin.info = "kim command style v1.0";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &kim_command_creator;
  (*register_plugin)(&plugin, lmp);
}
