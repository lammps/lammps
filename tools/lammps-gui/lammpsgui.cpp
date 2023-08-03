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

#include <QCoreApplication>
#include <QDesktopServices>
#include <QFileDialog>
#include <QFileInfo>
#include <QFont>
#include <QLabel>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QSettings>
#include <QShortcut>
#include <QStatusBar>
#include <QTextStream>
#include <QThread>
#include <QTimer>
#include <QUrl>

#include <cstring>
#include <string>

#if defined(_OPENMP)
#include <cstdlib>
#include <omp.h>
#endif

#if defined(_WIN32)
#include <io.h>
#else
#include <unistd.h>
#endif

static const QString blank(" ");

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
    dirstatus(nullptr), progress(nullptr), prefdialog(nullptr)
{
    ui->setupUi(this);
    this->setCentralWidget(ui->textEdit);
    highlighter = new Highlighter(ui->textEdit->document());
    capturer    = new StdCapture;
    current_file.clear();
    current_dir = QDir(".").absolutePath();
    recent_files.clear();

    QCoreApplication::setOrganizationName("The LAMMPS Developers");
    QCoreApplication::setOrganizationDomain("lammps.org");
    QCoreApplication::setApplicationName("LAMMPS GUI");

    // restorge and initialize settings
    QSettings settings;

    // check and initialize nthreads setting. Default is to use max,
    // but not override OMP_NUM_THREADS and preferences setting.
#if defined(_OPENMP)
    // use maximum number of available threads unless OMP_NUM_THREADS was set
    int nthreads = settings.value("nthreads", omp_get_max_threads()).toInt();
#if _WIN32
    if (!getenv("OMP_NUM_THREADS")) {
        _putenv_s("OMP_NUM_THREADS", std::to_string(nthreads).c_str());
    }
#else
    setenv("OMP_NUM_THREADS", std::to_string(nthreads).c_str(), 0);
#endif
#else
    int nthreads = settings.value("nthreads", 1).toInt();
#endif
    settings.setValue("nthreads", QString::number(nthreads));

    const char *tmpdir = getenv("TMPDIR");
    if (!tmpdir) tmpdir = getenv("TMP");
    if (!tmpdir) tmpdir = getenv("TEMPDIR");
    if (!tmpdir) tmpdir = getenv("TEMP");
#if _WIN32
    if (!tmpdir) tmpdir = "C:\\Windows\\Temp";
#else
    if (!tmpdir) tmpdir = "/tmp";
#endif
    settings.setValue("tempdir", QString(tmpdir));

    lammps_args.clear();
    lammps_args.push_back(mystrdup("LAMMPS-GUI"));
    lammps_args.push_back(mystrdup("-log"));
    lammps_args.push_back(mystrdup("none"));

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
    connect(ui->actionPreferences, &QAction::triggered, this, &LammpsGui::preferences);
    connect(ui->actionDefaults, &QAction::triggered, this, &LammpsGui::defaults);
    connect(ui->textEdit->document(), &QTextDocument::modificationChanged, this,
            &LammpsGui::modified);

#if !QT_CONFIG(clipboard)
    ui->actionCut->setEnabled(false);
    ui->actionCopy->setEnabled(false);
    ui->actionPaste->setEnabled(false);
#endif

    status = new QLabel("Ready.");
    status->setFixedWidth(300);
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
    plugin_path.clear();
    std::string deffile = settings.value("plugin_path", "liblammps.so").toString().toStdString();
    for (const char *libfile : {deffile.c_str(), "./liblammps.so", "liblammps.dylib",
                                "./liblammps.dylib", "liblammps.dll"}) {
        if (lammps.load_lib(libfile)) {
            auto canonical = QFileInfo(libfile).canonicalFilePath();
            plugin_path    = canonical.toStdString();
            settings.setValue("plugin_path", canonical);
            break;
        }
    }

    if (plugin_path.empty()) {
        // none of the plugin paths could load, remove key
        settings.remove("plugin_path");
        QMessageBox::critical(this, "Error", "Cannot open LAMMPS shared library file");
        // QCoreApplication::quit();
        exit(1);
    }
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

    lammps.close();
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
    if (text.back().toLatin1() != '\n') out << "\n"; // add final newline if missing
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
    lammps.close();
    lammps.finalize();

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
    lammps.force_timeout();
}

void LammpsGui::logupdate()
{
    double t_elapsed, t_remain, t_total;
    int completed = 1000;

    if (lammps.is_running()) {
        t_elapsed = lammps.get_thermo("cpu");
        t_remain  = lammps.get_thermo("cpuremain");
        t_total   = t_elapsed + t_remain + 1.0e-10;
        completed = t_elapsed / t_total * 1000.0;
    }

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

    if (lammps.has_error()) {
        lammps.get_last_error_message(errorbuf, BUFLEN);
        success = false;
    }

    if (success) {
        status->setText("Ready.");
    } else {
        status->setText("Failed.");
        QMessageBox::critical(this, "LAMMPS-GUI Error",
                              QString("Error running LAMMPS:\n\n") + errorbuf);
    }
    progress->hide();
    dirstatus->show();
}

void LammpsGui::run_buffer()
{
    QSettings settings;
    progress->setValue(0);
    dirstatus->hide();
    progress->show();
    int nthreads = settings.value("nthreads", 1).toInt();
    status->setText(QString("Running LAMMPS with %1 thread(s)...").arg(nthreads));
    status->repaint();
    start_lammps();
    if (!lammps.is_open()) return;
    clear();
    capturer->BeginCapture();

    // always add final newline since the text edit widget does not
    char *input = mystrdup(ui->textEdit->toPlainText().toStdString() + "\n");
    is_running  = true;

    LammpsRunner *runner = new LammpsRunner(this);
    runner->setup_run(&lammps, input);
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
    if (!lammps.is_running()) {
        start_lammps();
        if (!lammps.extract_setting("box_exists")) {
            QMessageBox::warning(this, "ImageViewer Error",
                                 "Cannot create snapshot image without a system box");
            return;
        }

        QSettings settings;
        QString dumpcmd  = "write_dump all image ";
        QString dumpfile = settings.value("tempdir").toString();
#if defined(_WIN32)
        dumpfile += '\\';
#else
        dumpfile += '/';
#endif
        dumpfile += current_file + ".ppm";
        dumpcmd += dumpfile;

        settings.beginGroup("snapshot");
        dumpcmd += blank + settings.value("color", "type").toString();
        dumpcmd += blank + settings.value("diameter", "type").toString();
        dumpcmd += QString(" size ") + settings.value("xsize", "800").toString();
        dumpcmd += blank + settings.value("ysize", "600").toString();
        dumpcmd += QString(" zoom ") + settings.value("zoom", "1.0").toString();
        settings.endGroup();

        lammps.command(dumpcmd.toLocal8Bit());
        imagewindow = new ImageViewer(dumpfile);
#if defined(_WIN32)
        _unlink(dumpfile.toLocal8Bit());
#else
        unlink(dumpfile.toLocal8Bit());
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
    ui->textEdit->moveCursor(QTextCursor::Start, QTextCursor::MoveAnchor);
}

void LammpsGui::about()
{
    std::string version = "This is LAMMPS-GUI version " LAMMPS_GUI_VERSION;
    if (lammps.has_plugin()) {
        version += " - LAMMPS linked dynamically";
        if (!plugin_path.empty()) {
            version += " from file ";
            version += plugin_path;
        }
    } else {
        version += " - LAMMPS linked statically";
    }
    std::string info = "LAMMPS is currently running. LAMMPS config info not available.";

    // LAMMPS is not re-entrant, so we can only query LAMMPS when it is not running
    if (!lammps.is_running()) {
        start_lammps();
        capturer->BeginCapture();
        lammps.command("info config");
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
    font.setFamily("Arial");
    font.setPointSize(8);
    msg.setFont(font);
    msg.exec();
}

void LammpsGui::help()
{
    QMessageBox msg;
    msg.setWindowTitle("LAMMPS-GUI Quick Help");
    msg.setText("<div>This is LAMMPS-GUI version " LAMMPS_GUI_VERSION "</div>");
    msg.setInformativeText("<b>Overview</b>"
                           "<div align=\"justify\">LAMMPS GUI is a graphical text editor that is "
                           "linked to the LAMMPS library and thus can run LAMMPS directly using "
                           "the contents of the text buffer as input through the LAMMPS C-library "
                           "interface. This makes it convenient to use for beginners and during "
                           "tutorials</div><br><br>"
                           "<b>Features</b>"
                           "<div align=\"justify\">The main window of the LAMMPS GUI is a text "
                           "editor window with syntax highlighting. The output of a LAMMPS run is "
                           "captured and displayed in a log window. The log window is updated "
                           "regularly during the run, as is a progress bar in the main window. "
                           "After the simulation is finished, an image of the simulated system "
                           "can be created and shown (and saved) in image viewer window.  Ongoing "
                           "runs can be stopped at the next run iteration.</div><br>"
                           "<div align=\"justify\">When opening a file, the editor will determine "
                           "the directory where the input file resides and switch its current "
                           "working directory to that same folder. Many LAMMPS inputs contain "
                           "commands that read other files, typically from the folder as the "
                           "input file. The GUI will show its current working directory. "
                           "In addition to using the menu, the editor window also receive files "
                           "as the first command line argument or via drag-n-drop from a "
                           "graphical file manager GUI or a desktop environment.</div><br>"
                           "<div align=\"justify\">Almost all commands are accessible via hotkeys. "
                           "Which those hotkeys are, is shown next to the entries in the menus. "
                           "Log and image viewer windows can be closed with CTRL-W (or Command-W "
                           "on macOS).</div><br>"
                           "<div align=\"justify\">The 'About LAMMPS' dialog will show the "
                           "LAMMPS version and the features included into the LAMMPS library "
                           "linked to the LAMMPS GUI.<br><br>"
                           "Due to its nature as a graphical application, it is <b>not</b> "
                           "possible to use the LAMMPS GUI in parallel with MPI, but OpenMP "
                           "multi-threading is available.</div>");
    msg.setIconPixmap(QPixmap(":/lammps-icon-128x128.png").scaled(64, 64));
    msg.setStandardButtons(QMessageBox::Close);
    msg.exec();
}

void LammpsGui::manual()
{
    QDesktopServices::openUrl(QUrl("https://docs.lammps.org/"));
}

void LammpsGui::defaults()
{
    QSettings settings;
    settings.clear();
    settings.sync();
}

void LammpsGui::preferences()
{
    Preferences prefs(&lammps);
    if (prefs.exec() == QDialog::Accepted) {
        // must delete LAMMPS instance after setting may be changed so we can apply different
        // suffixes
        lammps.close();
    }
}

void LammpsGui::start_lammps()
{
    // temporary extend lammps_args with additional arguments
    int initial_narg = lammps_args.size();
    QSettings settings;
    int nthreads = settings.value("nthreads", 1).toInt();
    int accel    = settings.value("accelerator", AcceleratorTab::None).toInt();
    if (accel == AcceleratorTab::Opt) {
        lammps_args.push_back(mystrdup("-suffix"));
        lammps_args.push_back(mystrdup("opt"));
    } else if (accel == AcceleratorTab::OpenMP) {
        lammps_args.push_back(mystrdup("-suffix"));
        lammps_args.push_back(mystrdup("omp"));
        lammps_args.push_back(mystrdup("-pk"));
        lammps_args.push_back(mystrdup("omp"));
        lammps_args.push_back(mystrdup(std::to_string(nthreads)));
    } else if (accel == AcceleratorTab::Intel) {
        lammps_args.push_back(mystrdup("-suffix"));
        lammps_args.push_back(mystrdup("intel"));
        lammps_args.push_back(mystrdup("-pk"));
        lammps_args.push_back(mystrdup("intel"));
        lammps_args.push_back(mystrdup(std::to_string(nthreads)));
    } else if (accel == AcceleratorTab::Gpu) {
        lammps_args.push_back(mystrdup("-suffix"));
        lammps_args.push_back(mystrdup("gpu"));
        lammps_args.push_back(mystrdup("-pk"));
        lammps_args.push_back(mystrdup("gpu"));
        lammps_args.push_back(mystrdup("0"));
    } else if (accel == AcceleratorTab::Kokkos) {
        lammps_args.push_back(mystrdup("-kokkos"));
        lammps_args.push_back(mystrdup("on"));
        lammps_args.push_back(mystrdup("t"));
        lammps_args.push_back(mystrdup(std::to_string(nthreads)));
        lammps_args.push_back(mystrdup("-suffix"));
        lammps_args.push_back(mystrdup("kk"));
    }
    if (settings.value("echo", "0").toInt()) {
        lammps_args.push_back(mystrdup("-echo"));
        lammps_args.push_back(mystrdup("screen"));
    }
    if (settings.value("cite", "0").toInt()) {
        lammps_args.push_back(mystrdup("-cite"));
        lammps_args.push_back(mystrdup("screen"));
    }

    char **args = lammps_args.data();
    int narg    = lammps_args.size();
    lammps.open(narg, args);

    // delete additional arguments again (3 were there initially
    while (lammps_args.size() > initial_narg) {
        delete lammps_args.back();
        lammps_args.pop_back();
    }

    if (lammps.has_error()) {
        constexpr int BUFLEN = 1024;
        char errorbuf[BUFLEN];
        lammps.get_last_error_message(errorbuf, BUFLEN);

        QMessageBox::critical(this, "LAMMPS-GUI Error",
                              QString("Error launching LAMMPS:\n\n") + errorbuf);
    }
}

// Local Variables:
// c-basic-offset: 4
// End:
