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

#include <QGridLayout>
#include <QList>
#include <QPair>
#include <QSpacerItem>
#include <QString>
#include <vector>

#include "lammpswrapper.h"

// forward declarations

class GeneralTab;
class LammpsRunner;
class LogWindow;

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
class ChartWindow;
class SlideShow;

class LammpsGui : public QMainWindow {
    Q_OBJECT

    friend class CodeEditor;
    friend class GeneralTab;

public:
    LammpsGui(QWidget *parent = nullptr, const char *filename = nullptr);
    ~LammpsGui() override;

protected:
    void open_file(const QString &filename);
    void write_file(const QString &filename);
    void update_recents(const QString &filename = "");
    void update_variables();
    void do_run(bool use_buffer);
    void start_lammps();
    void run_done();

public slots:
    void quit();
    void stop_run();

private slots:
    void new_document();
    void open();
    void open_recent();
    void start_exe();
    void save();
    void save_as();
    void copy();
    void cut();
    void paste();
    void undo();
    void redo();
    void run_buffer() { do_run(true); }
    void run_file() { do_run(false); }

    void edit_variables();
    void render_image();
    void view_slides();
    void view_image();
    void view_chart();
    void view_log();
    void view_variables();
    void about();
    void help();
    void manual();
    void howto();
    void logupdate();
    void modified();
    void preferences();
    void defaults();

protected:
    Ui::LammpsGui *ui;

private:
    Highlighter *highlighter;
    StdCapture *capturer;
    QLabel *status;
    LogWindow *logwindow;
    ImageViewer *imagewindow;
    ChartWindow *chartwindow;
    SlideShow *slideshow;
    QTimer *logupdater;
    QLabel *dirstatus;
    QProgressBar *progress;
    Preferences *prefdialog;
    QLabel *lammpsstatus;
    QLabel *varwindow;

    QString current_file;
    QString current_dir;
    QList<QString> recent;
    QList<QPair<QString, QString>> variables;

    LammpsWrapper lammps;
    LammpsRunner *runner;
    std::string plugin_path;
    bool is_running;
    int run_counter;
    std::vector<char *> lammps_args;
};
#endif // LAMMPSGUI_H

// Local Variables:
// c-basic-offset: 4
// End:
