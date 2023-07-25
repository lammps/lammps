#include "lammpsgui.h"
#include "ui_lammpsgui.h"

#include <QFileDialog>
#include <QFont>
#include <QMessageBox>
#include <QTextStream>
#include <string>

#include "library.h"

LammpsGui::LammpsGui(QWidget *parent, const char *filename) :
    QMainWindow(parent), ui(new Ui::LammpsGui), lammps_handle(nullptr)
{
    ui->setupUi(this);
    this->setCentralWidget(ui->textEdit);
    current_file.clear();

    QFont text_font;
    text_font.setFamily("monospace");
    text_font.setFixedPitch(true);
    text_font.setStyleHint(QFont::TypeWriter);
    ui->textEdit->document()->setDefaultFont(text_font);

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
    connect(ui->actionRun_Buffer, &QAction::triggered, this,
            &LammpsGui::run_buffer);
    connect(ui->actionExecute_Line, &QAction::triggered, this,
            &LammpsGui::run_line);
    connect(ui->actionAbout_LAMMPS_GUI, &QAction::triggered, this,
            &LammpsGui::about);
    connect(ui->actionLAMMPS_Info, &QAction::triggered, this,
            &LammpsGui::about_lammps);

#if !QT_CONFIG(clipboard)
    ui->actionCut->setEnabled(false);
    ui->actionCopy->setEnabled(false);
    ui->actionPaste->setEnabled(false);
#endif

    if (filename)
        open_file(filename);
    else
        setWindowTitle(QString("LAMMPS-GUI - *unknown*"));
}

LammpsGui::~LammpsGui()
{
    delete ui;
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

void LammpsGui::open_file(const QString &fileName)
{
    QFile file(fileName);
    current_file = fileName;
    if (!file.open(QIODevice::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, "Warning",
                             "Cannot open file: " + file.errorString());
        return;
    }
    setWindowTitle(QString("LAMMPS-GUI - " + fileName));
    QTextStream in(&file);
    QString text = in.readAll();
    ui->textEdit->document()->setPlainText(text);
    ui->textEdit->moveCursor(QTextCursor::Start, QTextCursor::MoveAnchor);
    file.close();
}

void LammpsGui::save()
{
    QString fileName;
    // If we don't have a filename from before, get one.
    if (current_file.isEmpty()) {
        fileName     = QFileDialog::getSaveFileName(this, "Save");
        current_file = fileName;
    } else {
        fileName = current_file;
    }
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, "Warning",
                             "Cannot save file: " + file.errorString());
        return;
    }
    setWindowTitle(QString("LAMMPS-GUI - " + fileName));
    QTextStream out(&file);
    QString text = ui->textEdit->toPlainText();
    out << text;
    file.close();
}

void LammpsGui::save_as()
{
    QString fileName = QFileDialog::getSaveFileName(this, "Save as");
    QFile file(fileName);

    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, "Warning",
                             "Cannot save file: " + file.errorString());
        return;
    }
    current_file = fileName;
    setWindowTitle(QString("LAMMPS-GUI - " + fileName));
    QTextStream out(&file);
    QString text = ui->textEdit->toPlainText();
    out << text;
    file.close();
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
    char *args[] = {(char *)"LAMMPS GUI", (char *)"-log", (char *)"none"};
    int nargs    = sizeof(args) / sizeof(char *);

    clear();
    if (!lammps_handle)
        lammps_handle = lammps_open_no_mpi(nargs, args, nullptr);
    if (!lammps_handle) return;
    std::string buffer = ui->textEdit->toPlainText().toStdString();
    lammps_commands_string(lammps_handle, buffer.c_str());
}

void LammpsGui::run_line()
{
    char *args[] = {(char *)"LAMMPS GUI", (char *)"-log", (char *)"none"};
    int nargs    = sizeof(args) / sizeof(char *);

    if (!lammps_handle)
        lammps_handle = lammps_open_no_mpi(nargs, args, nullptr);
    if (!lammps_handle) return;

    //    std::string buffer = ui->textEdit->toPlainText().toStdString();
    //    lammps_commands_string(lammps_handle, buffer.c_str());
}

void LammpsGui::clear()
{
    if (lammps_handle) {
        lammps_command(lammps_handle, "clear");
    }
    ui->textEdit->moveCursor(QTextCursor::Start, QTextCursor::MoveAnchor);
}

void LammpsGui::about_lammps()
{
    char *args[] = {(char *)"LAMMPS GUI", (char *)"-log", (char *)"none"};
    int nargs    = sizeof(args) / sizeof(char *);

    if (!lammps_handle)
        lammps_handle = lammps_open_no_mpi(nargs, args, nullptr);

    std::string version = "LAMMPS Version " + std::to_string(lammps_version(lammps_handle));
    QString lammps_info(version.c_str());
    QMessageBox::information(this, "About LAMMPS", lammps_info);
}

void LammpsGui::about()
{
    QMessageBox::information(this, "About LAMMPS-GUI",
                             "This is LAMMPS-GUI version 0.1");
}
