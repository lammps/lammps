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

#include "setvariables.h"

#include <QDialogButtonBox>
#include <QGridLayout>
#include <QLabel>
#include <QLineEdit>
#include <QSizePolicy>

SetVariables::SetVariables(QList<QPair<QString, QString>> &vars, QWidget *parent) : QDialog(parent)
{
    auto *layout = new QGridLayout;
    auto *top    = new QLabel("Set Variables:");
    layout->addWidget(top, 0, 0, 1, 2, Qt::AlignHCenter);

    auto *buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(buttonBox, &QDialogButtonBox::accepted, this, &SetVariables::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);

    int i = 1;
    for (const auto &v : vars) {
        auto *name = new QLineEdit(v.first);
        auto *val  = new QLineEdit(v.second);
        name->setObjectName("varname");
        val->setObjectName("varval");
        layout->addWidget(name, i, 0);
        layout->addWidget(val, i, 1);
        ++i;
    }
    layout->addItem(new QSpacerItem(10, 10, QSizePolicy::Expanding, QSizePolicy::Expanding), i, 0);
    layout->addWidget(buttonBox, i + 1, 0, 1, 2);
    setLayout(layout);
    setWindowTitle("LAMMPS-GUI - Set Variables");
    resize(500, 400);
}

void SetVariables::accept()
{
    // store all data in settings class
    // and then confirm accepting

    QDialog::accept();
}
