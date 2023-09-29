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

#include <QAction>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QIcon>
#include <QMenu>
#include <QMessageBox>
#include <QSettings>
#include <QString>
#include <QTextStream>

LogWindow::LogWindow(QWidget *parent) : QPlainTextEdit(parent)
{
    QSettings settings;
    resize(settings.value("logx", 500).toInt(), settings.value("logy", 320).toInt());
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

void LogWindow::save_as()
{
    QString fileName = QFileDialog::getSaveFileName(this, "Save Log to File");
    if (fileName.isEmpty()) return;

    QFileInfo path(fileName);
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
    connect(action, &QAction::triggered, this, &LogWindow::save_as);
    menu->exec(event->globalPos());
    delete menu;
}

// Local Variables:
// c-basic-offset: 4
// End:
