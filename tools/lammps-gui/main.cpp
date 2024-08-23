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
#include <QCommandLineOption>
#include <QCommandLineParser>
#include <QLocale>

#include <cstdio>
#include <cstring>

#define stringify(x) myxstr(x)
#define myxstr(x) #x

int main(int argc, char *argv[])
{
    Q_INIT_RESOURCE(lammpsgui);
    QApplication app(argc, argv);
    // enforce using the plain ASCII C locale within the GUI.
    QLocale::setDefault(QLocale::c());
    QCoreApplication::setOrganizationName("The LAMMPS Developers");
    QCoreApplication::setOrganizationDomain("lammps.org");
    QCoreApplication::setApplicationName("LAMMPS-GUI - QT" stringify(QT_VERSION_MAJOR));
    QCoreApplication::setApplicationVersion(LAMMPS_GUI_VERSION);
    QCommandLineParser parser;
    parser.setApplicationDescription(
        "\nThis is LAMMPS-GUI v" LAMMPS_GUI_VERSION "\n"
        "\nA graphical editor for LAMMPS input files with syntax highlighting and\n"
        "auto-completion that can run LAMMPS directly. It has built-in capabilities\n"
        "for monitoring, visualization, plotting, and capturing console output.");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument("file", "The LAMMPS input file to open (optional).");
    parser.process(app); // this removes known arguments

    const char *infile = nullptr;
    if (argc > 1) infile = argv[1];
    LammpsGui w(nullptr, infile);
    w.show();
    return app.exec();
}

// Local Variables:
// c-basic-offset: 4
// End:
