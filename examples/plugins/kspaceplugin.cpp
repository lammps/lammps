
#include "lammpsplugin.h"

#include "comm.h"
#include "command.h"
#include "error.h"
#include "version.h"

#include <cstring>

#include "kspace_zero2.h"

using namespace LAMMPS_NS;

static KSpace *zero2creator(LAMMPS *lmp)
{
  KSpace *ptr = (KSpace *) new KSpaceZero2(lmp);
  return ptr;
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  plugin.version = LAMMPS_VERSION;
  plugin.style = "kspace";
  plugin.name = "zero2";
  plugin.info = "zero2 KSpace style v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &zero2creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
