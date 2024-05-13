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

#ifndef LOGWINDOW_H
#define LOGWINDOW_H

#include <QPlainTextEdit>

class LogWindow : public QPlainTextEdit {
    Q_OBJECT

public:
    LogWindow(const QString &filename, QWidget *parent = nullptr);

private slots:
    void quit();
    void save_as();
    void stop_run();

protected:
    void closeEvent(QCloseEvent *event) override;
    void contextMenuEvent(QContextMenuEvent *event) override;
    bool eventFilter(QObject *watched, QEvent *event) override;

private:
    QString filename;
};

#endif
// Local Variables:
// c-basic-offset: 4
// End:
