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

#ifndef FIND_AND_REPLACE_H
#define FIND_AND_REPLACE_H

#include "codeeditor.h"
#include <QDialog>

class QLineEdit;
class QCheckBox;

class FindAndReplace : public QDialog {
    Q_OBJECT

public:
    explicit FindAndReplace(CodeEditor *_editor, QWidget *parent = nullptr);
    ~FindAndReplace() = default;

private slots:
    void find_next();
    void replace_next();
    void replace_all();
    void quit();

private:
    CodeEditor *editor;
    QLineEdit *search, *replace;
    QCheckBox *withcase, *wrap, *whole;
};

#endif

// Local Variables:
// c-basic-offset: 4
// End:
