/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef PREFERENCES_H
#define PREFERENCES_H

#include <QDialog>

class QDialogButtonBox;
class QSettings;
class QTabWidget;
class LammpsWrapper;

class Preferences : public QDialog {
    Q_OBJECT

public:
    explicit Preferences(LammpsWrapper *lammps, QWidget *parent = nullptr);
    ~Preferences() override;

private:
    QTabWidget *tabWidget;
    QDialogButtonBox *buttonBox;
    QSettings *settings;
    LammpsWrapper *lammps;
};

// individual tabs

class AcceleratorTab : public QWidget {
    Q_OBJECT

public:
    explicit AcceleratorTab(QSettings *settings, LammpsWrapper *lammps, QWidget *parent = nullptr);
    enum { None, Intel, Kokkos, OpenMP, Opt };

private:
    QSettings *settings;
    LammpsWrapper *lammps;
};

class SnapshotTab : public QWidget {
    Q_OBJECT

public:
    explicit SnapshotTab(QSettings *settings, QWidget *parent = nullptr);

private:
    QSettings *settings;
};

#endif

// Local Variables:
// c-basic-offset: 4
// End:
