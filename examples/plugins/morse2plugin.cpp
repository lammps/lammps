
#include "lammpsplugin.h"

#include "version.h"

#include <cstring>

#include "pair_morse2.h"
#include "pair_morse2_omp.h"

using namespace LAMMPS_NS;

static Pair *morse2creator(LAMMPS *lmp)
{
  return new PairMorse2(lmp);
}

static Pair *morse2ompcreator(LAMMPS *lmp)
{
  return new PairMorse2OMP(lmp);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  // register plain morse2 pair style
  plugin.version = LAMMPS_VERSION;
  plugin.style = "pair";
  plugin.name = "morse2";
  plugin.info = "Morse2 variant pair style v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &morse2creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);

  // also register morse2/omp pair style. only need to update changed fields
  plugin.name = "morse2/omp";
  plugin.info = "Morse2 variant pair style for OpenMP v1.0";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &morse2ompcreator;
  (*register_plugin)(&plugin, lmp);
}
