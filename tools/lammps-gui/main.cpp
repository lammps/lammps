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

#include "lammpsgui.h"

#include <QApplication>
#include <QFileInfo>

#include <cstdio>
#include <cstring>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    const char *infile = nullptr;
    if (argc > 1) {
        infile = argv[1];
        if ((strcmp(infile, "-help") == 0) || (strcmp(infile, "-h") == 0)) {
            printf("This is LAMMPS-GUI version " LAMMPS_GUI_VERSION
                   " using Qt version " QT_VERSION_STR "\n");
            printf("Usage: %s [-h|-help|<inputfile>]\n", argv[0]);
            return 1;
        }
    }

    LammpsGui w(nullptr, infile);
    w.show();
    return a.exec();
}

// Local Variables:
// c-basic-offset: 4
// End:
