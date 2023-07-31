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
#include "imageviewer.h"
#include "lammpsrunner.h"
#include "preferences.h"
#include "stdcapture.h"
#include "ui_lammpsgui.h"

#include <QDesktopServices>
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
#include <QUrl>

#include <cstring>
#include <string>

#if defined(LAMMPS_GUI_USE_PLUGIN)
#include "liblammpsplugin.h"
#else
#include "library.h"
#endif

#if defined(_OPENMP)
#include <cstdlib>
#include <omp.h>
#endif

#if defined(_WIN32)
#include <io.h>
#else
#include <unistd.h>
#endif

// duplicate string
static char *mystrdup(const std::string &text)
{
    auto tmp = new char[text.size() + 1];
    memcpy(tmp, text.c_str(), text.size() + 1);
    return tmp;
}

LammpsGui::LammpsGui(QWidget *parent, const char *filename) :
    QMainWindow(parent), ui(new Ui::LammpsGui), highlighter(nullptr), capturer(nullptr),
    status(nullptr), logwindow(nullptr), imagewindow(nullptr), logupdater(nullptr),
    dirstatus(nullptr), progress(nullptr), prefdialog(nullptr), lammps_handle(nullptr),
    plugin_handle(nullptr), plugin_path(nullptr), is_running(false)
{
    ui->setupUi(this);
    this->setCentralWidget(ui->textEdit);
    highlighter = new Highlighter(ui->textEdit->document());
    prefdialog  = new Preferences(this);
    capturer    = new StdCapture;
    current_file.clear();
    current_dir = QDir(".").absolutePath();
    recent_files.clear();

#if defined(_OPENMP)
    // use maximum number of available threads unless OMP_NUM_THREADS was set
    auto nthreads = std::to_string(omp_get_max_threads());
#if _WIN32
    if (!getenv("OMP_NUM_THREADS")) {
        _putenv_s("OMP_NUM_THREADS", nthreads.c_str());
    }
#else
    setenv("OMP_NUM_THREADS", nthreads.c_str(), 0);
#endif
#endif

    const char *tmpdir = getenv("TMPDIR");
    if (!tmpdir) tmpdir = getenv("TMP");
    if (!tmpdir) tmpdir = getenv("TEMPDIR");
    if (!tmpdir) tmpdir = getenv("TEMP");
#if _WIN32
    if (!tmpdir) tmpdir = "C:\\Windows\\Temp";
#else
    if (!tmpdir) tmpdir = "/tmp";
#endif
    temp_dir = tmpdir;

    lammps_args.clear();
    lammps_args.push_back(mystrdup("LAMMPS-GUI"));
    lammps_args.push_back(mystrdup("-log"));
    lammps_args.push_back(mystrdup("none"));
    lammps_args.push_back(mystrdup("-suffix"));
    lammps_args.push_back(mystrdup("omp"));

    setWindowIcon(QIcon(":/lammps-icon-128x128.png"));
#if (__APPLE__)
    QFont text_font("Menlo");
#else
    QFont text_font(":/Monospace.ttf");
#endif
    text_font.setStyleHint(QFont::TypeWriter);
    ui->textEdit->document()->setDefaultFont(text_font);
    ui->textEdit->setMinimumSize(600, 400);

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
    connect(ui->actionImage, &QAction::triggered, this, &LammpsGui::view_image);
    connect(ui->actionAbout_LAMMPS_GUI, &QAction::triggered, this, &LammpsGui::about);
    connect(ui->action_Help, &QAction::triggered, this, &LammpsGui::help);
    connect(ui->actionLAMMPS_Manual, &QAction::triggered, this, &LammpsGui::manual);
    connect(ui->actionEdit_Preferences, &QAction::triggered, this, &LammpsGui::preferences);
    connect(ui->textEdit->document(), &QTextDocument::modificationChanged, this,
            &LammpsGui::modified);

#if !QT_CONFIG(clipboard)
    ui->actionCut->setEnabled(false);
    ui->actionCopy->setEnabled(false);
    ui->actionPaste->setEnabled(false);
#endif

    status = new QLabel("Ready.");
    status->setFixedWidth(250);
    ui->statusbar->addWidget(status);
    dirstatus = new QLabel(QString(" Directory: ") + current_dir);
    dirstatus->setMinimumWidth(500);
    ui->statusbar->addWidget(dirstatus);
    progress = new QProgressBar();
    progress->setRange(0, 1000);
    progress->setMinimumWidth(500);
    progress->hide();
    dirstatus->show();
    ui->statusbar->addWidget(progress);

#if defined(LAMMPS_GUI_USE_PLUGIN)
    liblammpsplugin_t *lammps = nullptr;
    for (const auto libfile : {"liblammps.so", "./liblammps.so", "liblammps.dylib",
                               "./liblammps.dylib", "liblammps.dll"}) {
        if (!lammps) lammps = liblammpsplugin_load(libfile);
        if (lammps) {
            plugin_path = libfile;
            break;
        }
    }
    bool do_exit = !lammps || (lammps && lammps->abiversion != LAMMPSPLUGIN_ABI_VERSION);
    if (!lammps) QMessageBox::critical(this, "Error", "Cannot open LAMMPS shared library file");
    if (lammps && (lammps->abiversion != LAMMPSPLUGIN_ABI_VERSION))
        QMessageBox::critical(this, "Warning",
                              "ERROR: LAMMPS lib plugin ABI version does not match");
    plugin_handle = lammps;
    if (do_exit) exit(1);
#endif

    if (filename) {
        open_file(filename);
    } else {
        setWindowTitle(QString("LAMMPS-GUI - *unknown*"));
    }

}

LammpsGui::~LammpsGui()
{
    delete ui;
    delete highlighter;
    delete capturer;
    delete status;
    delete logwindow;
    delete imagewindow;
    delete dirstatus;
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
    if (ui->textEdit->document()->isModified()) {
        QMessageBox msg;
        msg.setWindowTitle("Unsaved Changes");
        msg.setWindowIcon(windowIcon());
        msg.setText(QString("The buffer ") + current_file + " has changes");
        msg.setInformativeText("Do you want to save the file before opening a new file?");
        msg.setIcon(QMessageBox::Question);
        msg.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
        int rv = msg.exec();
        switch (rv) {
            case QMessageBox::Yes:
                save();
                break;
            case QMessageBox::Cancel:
                return;
                break;
            case QMessageBox::No: // fallthrough
            default:
                // do nothing
                break;
        }
    }

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
    ui->textEdit->document()->setModified(false);
    file.close();
    dirstatus->setText(QString(" Directory: ") + current_dir);
}

void LammpsGui::write_file(const QString &fileName)
{
    QFile file(fileName);
    QFileInfo path(file);
    current_file = path.fileName();
    current_dir  = path.absolutePath();

    if (!file.open(QIODevice::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, "Warning", "Cannot save file: " + file.errorString());
        return;
    }
    setWindowTitle(QString("LAMMPS-GUI - " + current_file));
    QTextStream out(&file);
    QString text = ui->textEdit->toPlainText();
    out << text;
    file.close();
    dirstatus->setText(QString(" Directory: ") + current_dir);
    ui->textEdit->document()->setModified(false);
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
    if (ui->textEdit->document()->isModified()) {
        QMessageBox msg;
        msg.setWindowTitle("Unsaved Changes");
        msg.setWindowIcon(windowIcon());
        msg.setText(QString("The buffer ") + current_file + " has changes");
        msg.setInformativeText("Do you want to save the file before exiting?");
        msg.setIcon(QMessageBox::Question);
        msg.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
        int rv = msg.exec();
        switch (rv) {
            case QMessageBox::Yes:
                save();
                break;
            case QMessageBox::Cancel:
                return;
                break;
            case QMessageBox::No: // fallthrough
            default:
                // do nothing
                break;
        }
    }
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

void LammpsGui::modified()
{
    const QString modflag(" - *modified*");
    auto title = windowTitle().remove(modflag);
    if (ui->textEdit->document()->isModified())
        setWindowTitle(title + modflag);
    else
        setWindowTitle(title);
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
    progress->hide();
    dirstatus->show();
}

void LammpsGui::run_buffer()
{
    progress->setValue(0);
    dirstatus->hide();
    progress->show();
    status->setText("Running LAMMPS. Please wait...");
    status->repaint();
    start_lammps();
    if (!lammps_handle) return;
    clear();
    capturer->BeginCapture();

    char *input = mystrdup(ui->textEdit->toPlainText().toStdString());
    is_running  = true;

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
    logwindow->setWindowIcon(QIcon(":/lammps-icon-128x128.png"));
#if (__APPLE__)
    QFont text_font("Menlo");
#else
    QFont text_font(":/Monospace.ttf");
#endif
    text_font.setStyleHint(QFont::TypeWriter);
    logwindow->document()->setDefaultFont(text_font);
    logwindow->setLineWrapMode(QPlainTextEdit::NoWrap);
    logwindow->setMinimumSize(600, 400);
    QShortcut *shortcut = new QShortcut(QKeySequence(Qt::CTRL + Qt::Key_W), logwindow);
    QObject::connect(shortcut, &QShortcut::activated, logwindow, &QPlainTextEdit::close);
    shortcut = new QShortcut(QKeySequence(Qt::CTRL + Qt::Key_Slash), logwindow);
    QObject::connect(shortcut, &QShortcut::activated, this, &LammpsGui::stop_run);
    logwindow->show();

    logupdater = new QTimer(this);
    connect(logupdater, &QTimer::timeout, this, &LammpsGui::logupdate);
    logupdater->start(1000);
}

void LammpsGui::view_image()
{
    // LAMMPS is not re-entrant, so we can only query LAMMPS when it is not running
    if (!is_running) {
        start_lammps();
        int box = 0;
#if defined(LAMMPS_GUI_USE_PLUGIN)
        box = ((liblammpsplugin_t *)plugin_handle)->extract_setting(lammps_handle, "box_exists");
#else
        box = lammps_extract_setting(lammps_handle, "box_exist");
#endif
        if (!box) {
            QMessageBox::warning(this, "ImageViewer Error",
                                 "Cannot create snapshot image without a system box");
            return;
        }

        std::string dumpcmd = "write_dump all image '";
        QString dumpfile    = temp_dir;
#if defined(_WIN32)
        dumpfile += '\\';
#else
        dumpfile += '/';
#endif
        dumpfile += current_file + ".ppm";

        dumpcmd += dumpfile.toStdString() + "' type type size 800 600";

#if defined(LAMMPS_GUI_USE_PLUGIN)
        ((liblammpsplugin_t *)plugin_handle)->command(lammps_handle, dumpcmd.c_str());
#else
        lammps_command(lammps_handle, dumpcmd.c_str());
#endif
        imagewindow = new ImageViewer(dumpfile);
#if 0
#if defined(_WIN32)
        _unlink(dumpfile.toLocal8Bit());
#else
        unlink(dumpfile.toLocal8Bit());
#endif
#endif
    } else {
        QMessageBox::warning(this, "ImageViewer Error",
                             "Cannot create snapshot image while LAMMPS is running");
        return;
    }
    imagewindow->show();
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
    std::string version = "This is LAMMPS-GUI version " LAMMPS_GUI_VERSION;
#if defined(LAMMPS_GUI_USE_PLUGIN)
    version += " - LAMMPS linked dynamically";
    if (plugin_path) {
        version += " from file ";
        version += plugin_path;
    }
#else
    version += " - LAMMPS linked statically";
#endif
    std::string info = "LAMMPS is currently running. LAMMPS config info not available.";

    // LAMMPS is not re-entrant, so we can only query LAMMPS when it is not running
    if (!is_running) {
        start_lammps();
        capturer->BeginCapture();
#if defined(LAMMPS_GUI_USE_PLUGIN)
        ((liblammpsplugin_t *)plugin_handle)->commands(lammps_handle, "info config");
#else
        lammps_command(lammps_handle, "info config");
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
    msg.setIconPixmap(QPixmap(":/lammps-icon-128x128.png").scaled(64, 64));
    msg.setStandardButtons(QMessageBox::Ok);
    QFont font;
    font.setFixedPitch(true);
    font.setStyleHint(QFont::TypeWriter);
    font.setFamilies(QStringList({"Arial", "Helvetica"}));
    font.setPointSize(8);
    msg.setFont(font);
    msg.exec();
}

void LammpsGui::help()
{
    QString helpmsg = "This is LAMMPS-GUI version " LAMMPS_GUI_VERSION;
    QMessageBox::information(this, "LAMMPS-GUI Help", helpmsg);
}

void LammpsGui::manual()
{
    QDesktopServices::openUrl(QUrl("https://docs.lammps.org/"));
}

void LammpsGui::preferences()
{
    QString helpmsg = "This is LAMMPS-GUI version " LAMMPS_GUI_VERSION;
    QMessageBox::information(this, "LAMMPS-GUI Help", helpmsg);
}

void LammpsGui::start_lammps()
{
    char **args = lammps_args.data();
    int nargs   = lammps_args.size();

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
