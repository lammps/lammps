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

#ifndef SET_VARIABLES_H
#define SET_VARIABLES_H

#include <QDialog>
#include <QList>
#include <QPair>
#include <QString>

class SetVariables : public QDialog {
    Q_OBJECT

public:
    explicit SetVariables(QList<QPair<QString, QString>> &vars, QWidget *parent = nullptr);
    ~SetVariables() = default;

private slots:
    void accept() override;
    void add_row();
    void del_row();

private:
    QList<QPair<QString, QString>> &vars;
    class QVBoxLayout *layout;
};

#endif

// Local Variables:
// c-basic-offset: 4
// End:
