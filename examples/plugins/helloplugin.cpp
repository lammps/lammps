
#include "lammpsplugin.h"

#include "comm.h"
#include "command.h"
#include "error.h"
#include "version.h"

#include <cstring>

namespace LAMMPS_NS {
class Hello : public Command {
 public:
  Hello(class LAMMPS *lmp) : Command(lmp){};
  void command(int, char **) override;
};
}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;

void Hello::command(int argc, char **argv)
{
  if (argc != 1) error->all(FLERR, "Illegal hello command");
  if (comm->me == 0) utils::logmesg(lmp, fmt::format("Hello, {}!\n", argv[0]));
}

static Command *hellocreator(LAMMPS *lmp)
{
  return new Hello(lmp);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  plugin.version = LAMMPS_VERSION;
  plugin.style = "command";
  plugin.name = "hello";
  plugin.info = "Hello world command v1.1";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &hellocreator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
