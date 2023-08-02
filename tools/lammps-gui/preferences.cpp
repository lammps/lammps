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

#include "preferences.h"

#include <QDialogButtonBox>
#include <QDoubleValidator>
#include <QGroupBox>
#include <QIntValidator>
#include <QLabel>
#include <QLineEdit>
#include <QRadioButton>
#include <QSettings>
#include <QTabWidget>
#include <QVBoxLayout>

#include "lammpswrapper.h"

Preferences::Preferences(LammpsWrapper *_lammps, QWidget *parent) :
    QDialog(parent), tabWidget(new QTabWidget),
    buttonBox(new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel)),
    settings(new QSettings), lammps(_lammps)
{
    tabWidget->addTab(new AcceleratorTab(settings, lammps), "Accelerators");
    tabWidget->addTab(new SnapshotTab(settings), "Snapshot Image");

    connect(buttonBox, &QDialogButtonBox::accepted, this, &Preferences::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);

    auto *layout = new QVBoxLayout;
    layout->addWidget(tabWidget);
    layout->addWidget(buttonBox);
    setLayout(layout);
    setWindowTitle("LAMMPS-GUI - Preferences");
    resize(400, 320);
}

Preferences::~Preferences()
{
    delete buttonBox;
    delete tabWidget;
    delete settings;
}

void Preferences::accept()
{
    // store all data in settings class
    // and then confirm accepting
    QDialog::accept();
}

AcceleratorTab::AcceleratorTab(QSettings *_settings, LammpsWrapper *_lammps, QWidget *parent) :
    QWidget(parent), settings(_settings), lammps(_lammps)
{
    QGroupBox *accelerator = new QGroupBox("Choose Accelerator");
    QRadioButton *none     = new QRadioButton("&None");
    QRadioButton *kokkos   = new QRadioButton("&Kokkos");
    QRadioButton *intel    = new QRadioButton("&Intel");
    QRadioButton *openmp   = new QRadioButton("&OpenMP");
    QRadioButton *opt      = new QRadioButton("O&pt");

    none->setChecked(true);
    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(none);
    layout->addWidget(kokkos);
    layout->addWidget(intel);
    layout->addWidget(openmp);
    layout->addWidget(opt);
    layout->addStretch(1);
    setLayout(layout);
}

SnapshotTab::SnapshotTab(QSettings *_settings, QWidget *parent) :
    QWidget(parent), settings(_settings)
{
    QGridLayout *grid = new QGridLayout;

    auto *xsize = new QLabel("Image width:");
    auto *ysize = new QLabel("Image height:");
    auto *zoom  = new QLabel("Zoom factor:");
    auto *xval  = new QLineEdit(settings->value("xsize", "800").toString());
    auto *yval  = new QLineEdit(settings->value("ysize", "600").toString());
    auto *zval  = new QLineEdit(settings->value("zoom", "1.0").toString());

    auto *intval = new QIntValidator(100, 100000, this);
    xval->setValidator(intval);
    yval->setValidator(intval);
    zval->setValidator(new QDoubleValidator(0.01, 100.0, 100, this));

    grid->addWidget(xsize, 0, 0);
    grid->addWidget(ysize, 1, 0);
    grid->addWidget(zoom, 2, 0);
    grid->addWidget(xval, 0, 1);
    grid->addWidget(yval, 1, 1);
    grid->addWidget(zval, 2, 1);
    setLayout(grid);
}

// Local Variables:
// c-basic-offset: 4
// End:
