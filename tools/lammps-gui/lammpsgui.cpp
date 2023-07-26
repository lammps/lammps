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
#include "stdcapture.h"
#include "ui_lammpsgui.h"

//#include <QApplication>
#include <QFileDialog>
#include <QFileInfo>
#include <QFont>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QShortcut>
#include <QTextStream>
#include <string>

#include "library.h"

LammpsGui::LammpsGui(QWidget *parent, const char *filename) :
    QMainWindow(parent), ui(new Ui::LammpsGui), lammps_handle(nullptr)
{
    ui->setupUi(this);
    this->setCentralWidget(ui->textEdit);
    current_file.clear();
    capturer = new StdCapture;

    QFont text_font;
    text_font.setFamilies(QStringList({"Consolas", "Monospace", "Sans", "Courier"}));
    text_font.setFixedPitch(true);
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
    connect(ui->actionClear, &QAction::triggered, this, &LammpsGui::clear);
    connect(ui->actionRun_Buffer, &QAction::triggered, this, &LammpsGui::run_buffer);
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
}

LammpsGui::~LammpsGui()
{
    delete ui;
    delete highlighter;
    delete capturer;
    delete status;
}

void LammpsGui::new_document()
{
    current_file.clear();
    ui->textEdit->document()->setPlainText(QString());
    if (lammps_handle) lammps_close(lammps_handle);
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
    if (lammps_handle) {
        lammps_close(lammps_handle);
        lammps_mpi_finalize();
        lammps_kokkos_finalize();
        lammps_python_finalize();
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

void LammpsGui::run_buffer()
{
    status->setText("Running LAMMPS. Please wait...");
    status->repaint();
    start_lammps();
    if (!lammps_handle) return;
    clear();
    capturer->BeginCapture();
    std::string buffer = ui->textEdit->toPlainText().toStdString();
    lammps_commands_string(lammps_handle, buffer.c_str());
    capturer->EndCapture();
    auto log = capturer->GetCapture();
    status->setText("Ready.");

    auto box = new QPlainTextEdit();
    box->document()->setPlainText(log.c_str());
    box->setReadOnly(true);
    box->setWindowTitle("LAMMPS-GUI - Output from running LAMMPS on buffer - " + current_file);
    QFont text_font;
    text_font.setFamilies(QStringList({"Consolas", "Monospace", "Sans", "Courier"}));
    text_font.setFixedPitch(true);
    text_font.setStyleHint(QFont::TypeWriter);
    box->document()->setDefaultFont(text_font);
    box->setLineWrapMode(QPlainTextEdit::NoWrap);
    box->setMinimumSize(800, 600);
    QShortcut *shortcut = new QShortcut(QKeySequence(Qt::CTRL + Qt::Key_W), box);
    QObject::connect(shortcut, &QShortcut::activated, box, &QPlainTextEdit::close);
    box->show();

    if (lammps_has_error(lammps_handle)) {
        constexpr int BUFLEN = 1024;
        char errorbuf[BUFLEN];
        lammps_get_last_error_message(lammps_handle, errorbuf, BUFLEN);

        QMessageBox::critical(this, "LAMMPS-GUI Error",
                              QString("Error running LAMMPS:\n\n") + errorbuf);
    }
}

void LammpsGui::clear()
{
    if (lammps_handle) {
        lammps_command(lammps_handle, "clear");
    }
    ui->textEdit->moveCursor(QTextCursor::Start, QTextCursor::MoveAnchor);
}

void LammpsGui::about()
{
    start_lammps();

    std::string version = "This is LAMMPS-GUI version 0.1";
    capturer->BeginCapture();
    lammps_commands_string(lammps_handle, "info config");
    capturer->EndCapture();
    std::string info = capturer->GetCapture();
    auto start       = info.find("LAMMPS version:");
    auto end         = info.find("Info-Info-Info", start);
    QMessageBox msg;
    msg.setWindowTitle("About LAMMPS-GUI");
    msg.setText(version.c_str());
    msg.setInformativeText(std::string(info, start, end - start).c_str());
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

    // Create LAMMPS instance if not already present
    if (!lammps_handle) lammps_handle = lammps_open_no_mpi(nargs, args, nullptr);
    if (lammps_has_error(lammps_handle)) {
        constexpr int BUFLEN = 1024;
        char errorbuf[BUFLEN];
        lammps_get_last_error_message(lammps_handle, errorbuf, BUFLEN);

        QMessageBox::critical(this, "LAMMPS-GUI Error",
                             QString("Error launching LAMMPS:\n\n") + errorbuf);
    }
}

// Local Variables:
// c-basic-offset: 4
// End:
