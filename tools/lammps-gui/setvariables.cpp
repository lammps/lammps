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
#include <QIcon>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSizePolicy>

SetVariables::SetVariables(QList<QPair<QString, QString>> &_vars, QWidget *parent) :
    QDialog(parent), vars(_vars), layout(new QVBoxLayout)
{
    auto *top = new QLabel("Set Variables:");
    layout->addWidget(top, 0, Qt::AlignHCenter);

    int i = 1;
    for (const auto &v : vars) {
        auto *row  = new QHBoxLayout;
        auto *name = new QLineEdit(v.first);
        auto *val  = new QLineEdit(v.second);
        auto *del  = new QPushButton(QIcon(":/icons/edit-delete.png"), "");
        name->setObjectName("varname");
        val->setObjectName("varval");
        del->setObjectName(QString::number(i));
        connect(del, &QPushButton::released, this, &SetVariables::del_row);
        row->addWidget(name);
        row->addWidget(val);
        row->addWidget(del);
        layout->addLayout(row);
        ++i;
    }
    layout->addSpacerItem(new QSpacerItem(10, 10, QSizePolicy::Minimum, QSizePolicy::Expanding));

    auto *buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    auto *add       = new QPushButton("&Add Row");
    add->setObjectName("add_row");
    buttonBox->addButton(add, QDialogButtonBox::ActionRole);
    connect(add, &QPushButton::released, this, &SetVariables::add_row);
    connect(buttonBox, &QDialogButtonBox::accepted, this, &SetVariables::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);

    layout->addWidget(buttonBox);
    setLayout(layout);
    setWindowIcon(QIcon(":/icons/lammps-icon-128x128.png"));
    setWindowTitle("LAMMPS-GUI - Set Variables");
    resize(300, 200);
}

void SetVariables::accept()
{
    // store all data in variables class and then confirm accepting
    vars.clear();
    int nrows = layout->count() - 2;
    for (int i = 1; i < nrows; ++i) {
        auto *row = layout->itemAt(i)->layout();
        auto *var = dynamic_cast<QLineEdit *>(row->itemAt(0)->widget());
        auto *val = dynamic_cast<QLineEdit *>(row->itemAt(1)->widget());
        if (var && val) vars.append(qMakePair(var->text(), val->text()));
    }

    QDialog::accept();
}

void SetVariables::add_row()
{
    int nrows  = layout->count();
    auto *row  = new QHBoxLayout;
    auto *name = new QLineEdit(QString());
    auto *val  = new QLineEdit(QString());
    auto *del  = new QPushButton(QIcon(":/icons/edit-delete.png"), "");
    name->setObjectName("varname");
    val->setObjectName("varval");
    del->setObjectName(QString::number(nrows - 2));
    connect(del, &QPushButton::released, this, &SetVariables::del_row);
    row->addWidget(name);
    row->addWidget(val);
    row->addWidget(del);
    layout->insertLayout(nrows - 2, row);
}

void SetVariables::del_row()
{
    int nrows = layout->count();
    auto *who = sender();
    if (who) {
        // figure out which row was deleted and delete its layout and widgets
        int delrow = who->objectName().toInt();
        auto *row  = layout->takeAt(delrow);
        while (row->layout()->count() > 0) {
            auto *item = row->layout()->takeAt(0);
            if (item) {
                row->layout()->removeItem(item);
                delete item->widget();
                delete item;
            }
        }
        layout->removeItem(row);
        delete row->layout();

        // renumber the delete pushbutton names
        for (int i = delrow; i < nrows - 3; ++i) {
            auto *row    = layout->itemAt(i)->layout();
            auto *widget = row->itemAt(2)->widget();
            widget->setObjectName(QString::number(i));
        }
    }
}

// Local Variables:
// c-basic-offset: 4
// End:
