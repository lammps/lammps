
#include "lammpsplugin.h"

#include "lammps.h"
#include "version.h"

#include <cstring>

#include "pair_morse2.h"

using namespace LAMMPS_NS;

static Pair *morse2creator(LAMMPS *lmp)
{
    return new PairMorse2(lmp);
}

static lammpsplugin_t plugin;
void lammpsplugin_init(void *lmp)
{
    memset(&plugin,0,sizeof(lammpsplugin_t));
    plugin.version = LAMMPS_VERSION;
    plugin.style   = "pair";
    plugin.name    = "morse2";
    plugin.info    = "Morse2 variant pair style v1.0 by Axel Kohlmeyer";
    plugin.creator = (lammpsplugin_factory *) &morse2creator;
    
    lammpsplugin_register(&plugin, (LAMMPS *)lmp);
}
