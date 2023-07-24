#include "lammpsgui.h"
#include "./ui_lammpsgui.h"
#include "library.h"

#include <string>
#include <QFile>
#include <QFileDialog>
#include <QFont>
#include <QTextStream>
#include <QMessageBox>

LammpsGui::LammpsGui(QWidget *parent)
  : QMainWindow(parent), ui(new Ui::LammpsGui), lammps_handle(nullptr)
{
    ui->setupUi(this);
    this->setCentralWidget(ui->textEdit);
    current_file.clear();
    current_line = 0;
    text_font.setFamily("monospace");
    text_font.setFixedPitch(true);
    text_font.setStyleHint(QFont::TypeWriter);
    ui->textEdit->setCurrentFont(text_font);

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
    connect(ui->actionExecute_Line, &QAction::triggered, this, &LammpsGui::run_line);
//    connect(ui->actionAbout, &QAction::triggered, this, &LammpsGui::about);

#if !QT_CONFIG(clipboard)
    ui->actionCut->setEnabled(false);
    ui->actionCopy->setEnabled(false);
    ui->actionPaste->setEnabled(false);
#endif
}

LammpsGui::~LammpsGui()
{
    delete ui;
}

void LammpsGui::new_document()
{
  current_file.clear();
  current_line = 0;
  ui->textEdit->setText(QString());
  if (lammps_handle) lammps_close(lammps_handle);
  lammps_handle = nullptr;
}

void LammpsGui::open()
{
    QString fileName = QFileDialog::getOpenFileName(this, "Open the file");
    QFile file(fileName);
    current_file = fileName;
    current_line = 0;
    if (!file.open(QIODevice::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, "Warning", "Cannot open file: " + file.errorString());
        return;
    }
    setWindowTitle(fileName);
    QTextStream in(&file);
    QString text = in.readAll();
    ui->textEdit->setText(text);
    file.close();
}

void LammpsGui::save()
{
    QString fileName;
    // If we don't have a filename from before, get one.
    if (current_file.isEmpty()) {
        fileName = QFileDialog::getSaveFileName(this, "Save");
        current_file = fileName;
    } else {
        fileName = current_file;
    }
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, "Warning", "Cannot save file: " + file.errorString());
        return;
    }
    setWindowTitle(fileName);
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
        QMessageBox::warning(this, "Warning", "Cannot save file: " + file.errorString());
        return;
    }
    current_file = fileName;
    setWindowTitle(fileName);
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
  clear();
  if (!lammps_handle) lammps_handle = lammps_open_no_mpi(0, nullptr, nullptr);
  if (!lammps_handle) return;
  std::string buffer = ui->textEdit->toPlainText().toStdString();
  lammps_commands_string(lammps_handle, buffer.c_str());
}

void LammpsGui::run_line()
{
  // dummy
}

void LammpsGui::clear()
{
  if (lammps_handle) {
    lammps_command(lammps_handle, "clear");
    current_line = 0;
  }
}

void LammpsGui::about()
{
  // dummy
}
