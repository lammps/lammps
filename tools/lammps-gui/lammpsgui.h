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

#ifndef LAMMPSGUI_H
#define LAMMPSGUI_H

#include <QMainWindow>
#include <QString>

QT_BEGIN_NAMESPACE
namespace Ui {
class LammpsGui;
}
QT_END_NAMESPACE

class LammpsGui : public QMainWindow {
    Q_OBJECT

public:
    LammpsGui(QWidget *parent = nullptr, const char *filename = nullptr);
    ~LammpsGui();

protected:
    void open_file(const QString &filename);

private slots:
    void new_document();
    void open();
    void save();
    void save_as();
    void quit();
    void copy();
    void cut();
    void paste();
    void undo();
    void redo();
    void clear();
    void run_buffer();
    void run_line();
    void about();
    void about_lammps();

private:
    Ui::LammpsGui *ui;

    QString current_file;
    void *lammps_handle;
};

#endif // LAMMPSGUI_H
