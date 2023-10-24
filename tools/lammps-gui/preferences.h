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
class QFont;
class QSettings;
class QTabWidget;
class LammpsWrapper;

class Preferences : public QDialog {
    Q_OBJECT

public:
    explicit Preferences(LammpsWrapper *lammps, QWidget *parent = nullptr);
    ~Preferences() override;

private slots:
    void accept() override;

public:
    bool need_relaunch;

private:
    QTabWidget *tabWidget;
    QDialogButtonBox *buttonBox;
    QSettings *settings;
    LammpsWrapper *lammps;
};

// individual tabs

class GeneralTab : public QWidget {
    Q_OBJECT

public:
    explicit GeneralTab(QSettings *settings, LammpsWrapper *lammps, QWidget *parent = nullptr);

private slots:
    void pluginpath();
    void newallfont();
    void newtextfont();

private:
    void updatefonts(const QFont &all, const QFont &text);
    QSettings *settings;
    LammpsWrapper *lammps;
};

class AcceleratorTab : public QWidget {
    Q_OBJECT

public:
    explicit AcceleratorTab(QSettings *settings, LammpsWrapper *lammps, QWidget *parent = nullptr);
    enum { None, Opt, OpenMP, Intel, Kokkos, Gpu };

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

class EditorTab : public QWidget {
    Q_OBJECT

public:
    explicit EditorTab(QSettings *settings, QWidget *parent = nullptr);

private:
    QSettings *settings;
};

#endif

// Local Variables:
// c-basic-offset: 4
// End:
