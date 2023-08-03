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

#include <QList>
#include <QString>
#include <vector>

#include "lammpswrapper.h"

// forward declarations

QT_BEGIN_NAMESPACE
namespace Ui {
class LammpsGui;
}
QT_END_NAMESPACE

class QLabel;
class QPlainTextEdit;
class QProgressBar;
class QTimer;

class Highlighter;
class StdCapture;
class Preferences;
class ImageViewer;

class LammpsGui : public QMainWindow {
    Q_OBJECT

    friend class CodeEditor;

public:
    LammpsGui(QWidget *parent = nullptr, const char *filename = nullptr);
    ~LammpsGui() override;

protected:
    void open_file(const QString &filename);
    void write_file(const QString &filename);
    void start_lammps();
    void run_done();

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
    void stop_run();
    void view_image();
    void about();
    void help();
    void manual();
    void logupdate();
    void modified();
    void preferences();
    void defaults();

private:
    Ui::LammpsGui *ui;
    Highlighter *highlighter;
    StdCapture *capturer;
    QLabel *status;
    QPlainTextEdit *logwindow;
    ImageViewer *imagewindow;
    QTimer *logupdater;
    QLabel *dirstatus;
    QProgressBar *progress;
    Preferences *prefdialog;

    QString current_file;
    QString current_dir;
    LammpsWrapper lammps;
    std::string plugin_path;
    bool is_running;
    QList<QString> recent_files;
    std::vector<char *> lammps_args;
};
#endif // LAMMPSGUI_H

// Local Variables:
// c-basic-offset: 4
// End:
