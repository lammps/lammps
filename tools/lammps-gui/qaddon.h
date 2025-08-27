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

#ifndef QADDON_H
#define QADDON_H

#include <QCompleter>
#include <QFrame>
#include <QValidator>

// draw horizontal line
class QHline : public QFrame {
public:
    QHline(QWidget *parent = nullptr);
};

// complete color inputs
class QColorCompleter : public QCompleter {
public:
    QColorCompleter(QWidget *parent = nullptr);
};

// validate color inputs
class QColorValidator : public QValidator {
public:
    QColorValidator(QWidget *parent = nullptr);

    void fixup(QString &input) const override;
    QValidator::State validate(QString &input, int &pos) const override;
};

#endif

// Local Variables:
// c-basic-offset: 4
// End:
