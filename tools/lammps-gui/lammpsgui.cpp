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

#include "lammpsgui.h"

#include "highlighter.h"
#include "lammpsrunner.h"
#include "stdcapture.h"
#include "ui_lammpsgui.h"

#include <QFileDialog>
#include <QFileInfo>
#include <QFont>
#include <QLabel>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QShortcut>
#include <QStatusBar>
#include <QTextStream>
#include <QTimer>
#include <cstring>
#include <string>

#if defined(LAMMPS_GUI_USE_PLUGIN)
#include "liblammpsplugin.h"
#else
#include "library.h"
#endif

LammpsGui::LammpsGui(QWidget *parent, const char *filename) :
    QMainWindow(parent), ui(new Ui::LammpsGui), highlighter(nullptr), capturer(nullptr),
    status(nullptr), logwindow(nullptr), logupdater(nullptr), progress(nullptr),
    lammps_handle(nullptr), plugin_handle(nullptr), is_running(false)
{
    ui->setupUi(this);
    this->setCentralWidget(ui->textEdit);
    current_file.clear();
    capturer = new StdCapture;

    QFont text_font;
    text_font.setStyleHint(QFont::TypeWriter);
    ui->textEdit->document()->setDefaultFont(text_font);
    ui->textEdit->setMinimumSize(800, 600);
    highlighter = new Highlighter(ui->textEdit->document());

    connect(ui->actionNew, &QAction::triggered, this, &LammpsGui::new_document);
    connect(ui->actionOpen, &QAction::triggered, this, &LammpsGui::open);
    connect(ui->actionSave, &QAction::triggered, this, &LammpsGui::save);
    connect(ui->actionSave_As, &QAction::triggered, this, &LammpsGui::save_as);
    connect(ui->actionQuit, &QAction::triggered, this, &LammpsGui::quit);
    connect(ui->actionCopy, &QAction::triggered, this, &LammpsGui::copy);
    connect(ui->actionCut, &QAction::triggered, this, &LammpsGui::cut);
    connect(ui->actionPaste, &QAction::triggered, this, &LammpsGui::paste);
    connect(ui->actionUndo, &QAction::triggered, this, &LammpsGui::undo);
    connect(ui->actionRedo, &QAction::triggered, this, &LammpsGui::redo);
    connect(ui->actionRun_Buffer, &QAction::triggered, this, &LammpsGui::run_buffer);
    connect(ui->actionStop_LAMMPS, &QAction::triggered, this, &LammpsGui::stop_run);
    connect(ui->actionAbout_LAMMPS_GUI, &QAction::triggered, this, &LammpsGui::about);

#if !QT_CONFIG(clipboard)
    ui->actionCut->setEnabled(false);
    ui->actionCopy->setEnabled(false);
    ui->actionPaste->setEnabled(false);
#endif

    if (filename)
        open_file(filename);
    else
        setWindowTitle(QString("LAMMPS-GUI - *unknown*"));
    status = new QLabel("Ready.");
    status->setFixedWidth(300);
    ui->statusbar->addWidget(status);
    progress = new QProgressBar();
    progress->setRange(0, 1000);
    progress->setFixedWidth(500);
    ui->statusbar->addWidget(progress);

#if defined(LAMMPS_GUI_USE_PLUGIN)
    liblammpsplugin_t *lammps = liblammpsplugin_load("liblammps.so");
    if (!lammps) lammps = liblammpsplugin_load("liblammps.dylib");
    if (!lammps) lammps = liblammpsplugin_load("liblammps.dll");
    bool do_exit = !lammps || (lammps && lammps->abiversion != LAMMPSPLUGIN_ABI_VERSION);
    if (!lammps) QMessageBox::critical(this, "Warning", "Cannot open LAMMPS shared library file");
    if (lammps && (lammps->abiversion != LAMMPSPLUGIN_ABI_VERSION))
        QMessageBox::critical(this, "Warning",
                              "ERROR: LAMMPS lib plugin ABI version does not match");
    plugin_handle = lammps;
    if (do_exit) exit(1);
#endif
}

LammpsGui::~LammpsGui()
{
    delete ui;
    delete highlighter;
    delete capturer;
    delete status;
    delete logwindow;
}

void LammpsGui::new_document()
{
    current_file.clear();
    ui->textEdit->document()->setPlainText(QString());

#if defined(LAMMPS_GUI_USE_PLUGIN)
    if (lammps_handle) ((liblammpsplugin_t *)plugin_handle)->close(lammps_handle);
#else
    if (lammps_handle) lammps_close(lammps_handle);
#endif
    lammps_handle = nullptr;
    setWindowTitle(QString("LAMMPS-GUI - *unknown*"));
}

void LammpsGui::open()
{
    QString fileName = QFileDialog::getOpenFileName(this, "Open the file");
    open_file(fileName);
}

// open file and switch CWD to path of file
void LammpsGui::open_file(const QString &fileName)
{

    QFileInfo path(fileName);
    current_file = path.fileName();
    current_dir  = path.absolutePath();
    QFile file(path.absoluteFilePath());

    QDir::setCurrent(current_dir);
    if (!file.open(QIODevice::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, "Warning",
                             "Cannot open file " + path.absoluteFilePath() + ": " +
                                 file.errorString());
        return;
    }
    setWindowTitle(QString("LAMMPS-GUI - " + current_file));
    QTextStream in(&file);
    QString text = in.readAll();
    ui->textEdit->document()->setPlainText(text);
    ui->textEdit->moveCursor(QTextCursor::Start, QTextCursor::MoveAnchor);
    file.close();
}

void LammpsGui::write_file(const QString &fileName)
{
    QFile file(fileName);
    QFileInfo path(file);
    current_file = path.fileName();
    if (!file.open(QIODevice::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, "Warning", "Cannot save file: " + file.errorString());
        return;
    }
    setWindowTitle(QString("LAMMPS-GUI - " + current_file));
    QTextStream out(&file);
    QString text = ui->textEdit->toPlainText();
    out << text;
    file.close();
}

void LammpsGui::save()
{
    QString fileName = current_file;
    // If we don't have a filename from before, get one.
    if (fileName.isEmpty()) fileName = QFileDialog::getSaveFileName(this, "Save");

    write_file(fileName);
}

void LammpsGui::save_as()
{
    QString fileName = QFileDialog::getSaveFileName(this, "Save as");
    write_file(fileName);
}

void LammpsGui::quit()
{
#if defined(LAMMPS_GUI_USE_PLUGIN)
    if (lammps_handle) {
        liblammpsplugin_t *lammps = (liblammpsplugin_t *)plugin_handle;
        lammps->close(lammps_handle);
        lammps->mpi_finalize();
        lammps->kokkos_finalize();
        lammps->python_finalize();
    }
#else
    if (lammps_handle) {
        lammps_close(lammps_handle);
        lammps_mpi_finalize();
        lammps_kokkos_finalize();
        lammps_python_finalize();
    }
#endif
    QCoreApplication::quit();
}

void LammpsGui::copy()
{
#if QT_CONFIG(clipboard)
    ui->textEdit->copy();
#endif
}

void LammpsGui::cut()
{
#if QT_CONFIG(clipboard)
    ui->textEdit->cut();
#endif
}

void LammpsGui::paste()
{
#if QT_CONFIG(clipboard)
    ui->textEdit->paste();
#endif
}

void LammpsGui::undo()
{
    ui->textEdit->undo();
}

void LammpsGui::redo()
{
    ui->textEdit->redo();
}

void LammpsGui::stop_run()
{
#if defined(LAMMPS_GUI_USE_PLUGIN)
    ((liblammpsplugin_t *)plugin_handle)->force_timeout(lammps_handle);
#else
    lammps_force_timeout(lammps_handle);
#endif
}

void LammpsGui::logupdate()
{
    double t_elapsed, t_remain, t_total;
    int completed = 1000;

#if defined(LAMMPS_GUI_USE_PLUGIN)
    liblammpsplugin_t *lammps = (liblammpsplugin_t *)plugin_handle;
    if (lammps->is_running(lammps_handle)) {
        t_elapsed = lammps->get_thermo(lammps_handle, "cpu");
        t_remain  = lammps->get_thermo(lammps_handle, "cpuremain");
        t_total   = t_elapsed + t_remain + 1.0e-10;
        completed = t_elapsed / t_total * 1000.0;
    }
#else
    if (lammps_is_running(lammps_handle)) {
        t_elapsed = lammps_get_thermo(lammps_handle, "cpu");
        t_remain  = lammps_get_thermo(lammps_handle, "cpuremain");
        t_total   = t_elapsed + t_remain + 1.0e-10;
        completed = t_elapsed / t_total * 1000.0;
    }
#endif

    progress->setValue(completed);
    if (logwindow) {
        const auto text = capturer->GetChunk();
        if (text.size() > 0) {
            logwindow->insertPlainText(text.c_str());
            logwindow->moveCursor(QTextCursor::End);
            logwindow->textCursor().deleteChar();
        }
    }
}

void LammpsGui::run_done()
{
    logupdater->stop();
    delete logupdater;
    logupdater = nullptr;
    progress->setValue(1000);

    capturer->EndCapture();
    auto log = capturer->GetCapture();
    logwindow->insertPlainText(log.c_str());
    logwindow->moveCursor(QTextCursor::End);

    bool success         = true;
    constexpr int BUFLEN = 1024;
    char errorbuf[BUFLEN];

#if defined(LAMMPS_GUI_USE_PLUGIN)
    liblammpsplugin_t *lammps = (liblammpsplugin_t *)plugin_handle;
    if (lammps->has_error(lammps_handle)) {
        lammps->get_last_error_message(lammps_handle, errorbuf, BUFLEN);
        success = false;
    }
#else
    if (lammps_has_error(lammps_handle)) {
        lammps_get_last_error_message(lammps_handle, errorbuf, BUFLEN);
        success = false;
    }
#endif

    if (success) {
        status->setText("Ready.");
    } else {
        status->setText("Failed.");
        QMessageBox::critical(this, "LAMMPS-GUI Error",
                              QString("Error running LAMMPS:\n\n") + errorbuf);
    }
    is_running = false;
}

void LammpsGui::run_buffer()
{
    status->setText("Running LAMMPS. Please wait...");
    status->repaint();
    start_lammps();
    if (!lammps_handle) return;
    clear();
    capturer->BeginCapture();

    std::string buffer = ui->textEdit->toPlainText().toStdString();
    char *input        = new char[buffer.size() + 1];
    memcpy(input, buffer.c_str(), buffer.size() + 1);

    is_running           = true;
    LammpsRunner *runner = new LammpsRunner(this);
    runner->setup_run(lammps_handle, input, plugin_handle);
    connect(runner, &LammpsRunner::resultReady, this, &LammpsGui::run_done);
    connect(runner, &LammpsRunner::finished, runner, &QObject::deleteLater);
    runner->start();

    logwindow = new QPlainTextEdit();
    logwindow->setReadOnly(true);
    logwindow->setCenterOnScroll(true);
    logwindow->moveCursor(QTextCursor::End);
    logwindow->setWindowTitle("LAMMPS-GUI - Output from running LAMMPS on buffer - " +
                              current_file);
    QFont text_font;
    text_font.setStyleHint(QFont::TypeWriter);
    logwindow->document()->setDefaultFont(text_font);
    logwindow->setLineWrapMode(QPlainTextEdit::NoWrap);
    logwindow->setMinimumSize(800, 600);
    QShortcut *shortcut = new QShortcut(QKeySequence(Qt::CTRL + Qt::Key_W), logwindow);
    QObject::connect(shortcut, &QShortcut::activated, logwindow, &QPlainTextEdit::close);
    shortcut = new QShortcut(QKeySequence(Qt::CTRL + Qt::Key_Slash), logwindow);
    QObject::connect(shortcut, &QShortcut::activated, this, &LammpsGui::stop_run);
    logwindow->show();

    logupdater = new QTimer(this);
    connect(logupdater, &QTimer::timeout, this, &LammpsGui::logupdate);
    logupdater->start(1000);
}

void LammpsGui::clear()
{
    if (lammps_handle) {
#if defined(LAMMPS_GUI_USE_PLUGIN)
        ((liblammpsplugin_t *)plugin_handle)->command(lammps_handle, "clear");
#else
        lammps_command(lammps_handle, "clear");
#endif
    }
    ui->textEdit->moveCursor(QTextCursor::Start, QTextCursor::MoveAnchor);
}

void LammpsGui::about()
{
    std::string version = "This is LAMMPS-GUI version 0.1";
    std::string info    = "LAMMPS is currently running. LAMMPS config info not available.";

    // LAMMPS is not re-entrant, so we can only query LAMMPS when it is not running
    if (!is_running) {
        start_lammps();
        capturer->BeginCapture();
#if defined(LAMMPS_GUI_USE_PLUGIN)
        ((liblammpsplugin_t *)plugin_handle)->commands_string(lammps_handle, "info config");
#else
        lammps_commands_string(lammps_handle, "info config");
#endif
        capturer->EndCapture();
        info       = capturer->GetCapture();
        auto start = info.find("LAMMPS version:");
        auto end   = info.find("Info-Info-Info", start);
        info       = std::string(info, start, end - start);
    }

    QMessageBox msg;
    msg.setWindowTitle("About LAMMPS-GUI");
    msg.setText(version.c_str());
    msg.setInformativeText(info.c_str());
    msg.setIcon(QMessageBox::NoIcon);
    msg.setStandardButtons(QMessageBox::Ok);
    QFont font;
    font.setFamilies(QStringList({"Sans", "Arial", "Helvetica"}));
    font.setFixedPitch(true);
    font.setStyleHint(QFont::TypeWriter);
    font.setPointSize(8);
    font.setWeight(QFont::Medium);
    msg.setFont(font);
    msg.exec();
}

void LammpsGui::start_lammps()
{
    char *args[] = {(char *)"LAMMPS GUI", (char *)"-log", (char *)"none"};
    int nargs    = sizeof(args) / sizeof(char *);

#if defined(LAMMPS_GUI_USE_PLUGIN)
    liblammpsplugin_t *lammps = (liblammpsplugin_t *)plugin_handle;

    // Create LAMMPS instance if not already present
    if (!lammps_handle) lammps_handle = lammps->open_no_mpi(nargs, args, nullptr);
    if (lammps->has_error(lammps_handle)) {
        constexpr int BUFLEN = 1024;
        char errorbuf[BUFLEN];
        lammps->get_last_error_message(lammps_handle, errorbuf, BUFLEN);

        QMessageBox::critical(this, "LAMMPS-GUI Error",
                              QString("Error launching LAMMPS:\n\n") + errorbuf);
    }
#else
    // Create LAMMPS instance if not already present
    if (!lammps_handle) lammps_handle = lammps_open_no_mpi(nargs, args, nullptr);
    if (lammps_has_error(lammps_handle)) {
        constexpr int BUFLEN = 1024;
        char errorbuf[BUFLEN];
        lammps_get_last_error_message(lammps_handle, errorbuf, BUFLEN);

        QMessageBox::critical(this, "LAMMPS-GUI Error",
                              QString("Error launching LAMMPS:\n\n") + errorbuf);
    }
#endif
}

// Local Variables:
// c-basic-offset: 4
// End:
