/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

#include "lammpsrunner.h"

#if defined(LAMMPS_GUI_USE_PLUGIN)
#include "liblammpsplugin.h"
#else
#include "library.h"
#endif

#include <cstdio>

LammpsRunner::LammpsRunner(QObject *parent) :
    QThread(parent), handle(nullptr), plugin(nullptr), input(nullptr)
{
}

void LammpsRunner::run()
{
    if (handle) {
#if defined(LAMMPS_GUI_USE_PLUGIN)
        liblammpsplugin_t *lammps = (liblammpsplugin_t *)plugin;
        lammps->commands_string(handle, input);
#else
        lammps_commands_string(handle, input);
#endif
    }
    emit resultReady();
}

// Local Variables:
// c-basic-offset: 4
// End:
