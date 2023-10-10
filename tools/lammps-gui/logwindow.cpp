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

#include "logwindow.h"

#include "lammpsgui.h"

#include <QAction>
#include <QApplication>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QIcon>
#include <QKeySequence>
#include <QMenu>
#include <QMessageBox>
#include <QSettings>
#include <QShortcut>
#include <QString>
#include <QTextStream>

LogWindow::LogWindow(const QString &_filename, QWidget *parent) :
    QPlainTextEdit(parent), filename(_filename)
{
    QSettings settings;
    resize(settings.value("logx", 500).toInt(), settings.value("logy", 320).toInt());

    auto action = new QShortcut(QKeySequence::fromString("Ctrl+S"), this);
    connect(action, &QShortcut::activated, this, &LogWindow::save_as);
    action = new QShortcut(QKeySequence::fromString("Ctrl+Q"), this);
    connect(action, &QShortcut::activated, this, &LogWindow::quit);
    action = new QShortcut(QKeySequence(Qt::Key_Slash, Qt::CTRL), this);
    connect(action, &QShortcut::activated, this, &LogWindow::stop_run);
}

void LogWindow::closeEvent(QCloseEvent *event)
{
    QSettings settings;
    if (!isMaximized()) {
        settings.setValue("logx", width());
        settings.setValue("logy", height());
    }
    QPlainTextEdit::closeEvent(event);
}

void LogWindow::quit()
{
    LammpsGui *main;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") main = dynamic_cast<LammpsGui *>(widget);
    main->quit();
}

void LogWindow::stop_run()
{
    LammpsGui *main;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") main = dynamic_cast<LammpsGui *>(widget);
    main->stop_run();
}

void LogWindow::save_as()
{
    QString defaultname = filename + ".log";
    if (filename.isEmpty()) defaultname = "lammps.log";
    QString logFileName = QFileDialog::getSaveFileName(this, "Save Log to File", defaultname,
                                                       "Log files (*.log *.out *.txt)");
    if (logFileName.isEmpty()) return;

    QFileInfo path(logFileName);
    QFile file(path.absoluteFilePath());

    if (!file.open(QIODevice::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, "Warning", "Cannot save file: " + file.errorString());
        return;
    }

    QTextStream out(&file);
    QString text = toPlainText();
    out << text;
    if (text.back().toLatin1() != '\n') out << "\n"; // add final newline if missing
    file.close();
}

void LogWindow::contextMenuEvent(QContextMenuEvent *event)
{
    // show augmented context menu
    auto *menu = createStandardContextMenu();
    menu->addSeparator();
    auto action = menu->addAction(QString("Save Log to File ..."));
    action->setIcon(QIcon(":/icons/document-save-as.png"));
    action->setShortcut(QKeySequence::fromString("Ctrl+S"));
    connect(action, &QAction::triggered, this, &LogWindow::save_as);
    action = menu->addAction("&Close Window", this, &QWidget::close);
    action->setIcon(QIcon(":/icons/window-close.png"));
    action->setShortcut(QKeySequence::fromString("Ctrl+W"));
    menu->exec(event->globalPos());
    delete menu;
}

// Local Variables:
// c-basic-offset: 4
// End:
