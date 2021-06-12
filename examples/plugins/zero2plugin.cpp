
#include "lammpsplugin.h"
#include "version.h"

#include <cstring>

#include "angle_zero2.h"
#include "bond_zero2.h"
#include "dihedral_zero2.h"
#include "improper_zero2.h"
#include "pair_zero2.h"

using namespace LAMMPS_NS;

static Pair *pairzerocreator(LAMMPS *lmp)
{
  return new PairZero2(lmp);
}

static Bond *bondzerocreator(LAMMPS *lmp)
{
  return new BondZero2(lmp);
}

static Angle *anglezerocreator(LAMMPS *lmp)
{
  return new AngleZero2(lmp);
}

static Dihedral *dihedralzerocreator(LAMMPS *lmp)
{
  return new DihedralZero2(lmp);
}

static Improper *improperzerocreator(LAMMPS *lmp)
{
  return new ImproperZero2(lmp);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  // register zero2 pair style
  plugin.version = LAMMPS_VERSION;
  plugin.style = "pair";
  plugin.name = "zero2";
  plugin.info = "Zero2 variant pair style v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &pairzerocreator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);

  // register zero2 bond style
  plugin.style = "bond";
  plugin.info = "Zero2 variant bond style v1.0";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &bondzerocreator;
  (*register_plugin)(&plugin, lmp);

  // register zero2 angle style
  plugin.style = "angle";
  plugin.info = "Zero2 variant angle style v1.0";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &anglezerocreator;
  (*register_plugin)(&plugin, lmp);

  // register zero2 dihedral style
  plugin.style = "dihedral";
  plugin.info = "Zero2 variant dihedral style v1.0";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &dihedralzerocreator;
  (*register_plugin)(&plugin, lmp);

  // register zero2 improper style
  plugin.style = "improper";
  plugin.info = "Zero2 variant improper style v1.0";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &improperzerocreator;
  (*register_plugin)(&plugin, lmp);
}
