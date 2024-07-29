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

#include "fileviewer.h"

#include "lammpsgui.h"

#include <QApplication>
#include <QFile>
#include <QIcon>
#include <QKeySequence>
#include <QSettings>
#include <QShortcut>
#include <QString>
#include <QTextCursor>
#include <QTextStream>

FileViewer::FileViewer(const QString &_filename, QWidget *parent) :
    QPlainTextEdit(parent), fileName(_filename)
{
    auto *action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_Q), this);
    connect(action, &QShortcut::activated, this, &FileViewer::quit);
    action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_Slash), this);
    connect(action, &QShortcut::activated, this, &FileViewer::stop_run);

    installEventFilter(this);

    // open and read file. Set editor to read-only.
    QFile file(fileName);
    if (file.open(QIODevice::Text | QIODevice::ReadOnly)) {
        QTextStream in(&file);
        QString content = in.readAll();
        file.close();

        QFont text_font;
        QSettings settings;
        text_font.fromString(settings.value("textfont", text_font.toString()).toString());
        document()->setDefaultFont(text_font);

        document()->setPlainText(content);
        moveCursor(QTextCursor::Start, QTextCursor::MoveAnchor);
        setReadOnly(true);
        setLineWrapMode(NoWrap);
        setMinimumSize(800, 500);
        setWindowIcon(QIcon(":/icons/lammps-icon-128x128.png"));
        setWindowTitle("LAMMPS-GUI - Viewer - " + fileName);
    }
}

void FileViewer::quit()
{
    LammpsGui *main = nullptr;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") main = dynamic_cast<LammpsGui *>(widget);
    if (main) main->quit();
}

void FileViewer::stop_run()
{
    LammpsGui *main = nullptr;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") main = dynamic_cast<LammpsGui *>(widget);
    if (main) main->stop_run();
}

// event filter to handle "Ambiguous shortcut override" issues
bool FileViewer::eventFilter(QObject *watched, QEvent *event)
{
    if (event->type() == QEvent::ShortcutOverride) {
        auto *keyEvent = dynamic_cast<QKeyEvent *>(event);
        if (!keyEvent) return QWidget::eventFilter(watched, event);
        if (keyEvent->modifiers().testFlag(Qt::ControlModifier) && keyEvent->key() == '/') {
            stop_run();
            event->accept();
            return true;
        }
        if (keyEvent->modifiers().testFlag(Qt::ControlModifier) && keyEvent->key() == 'W') {
            close();
            event->accept();
            return true;
        }
    }
    return QWidget::eventFilter(watched, event);
}

// Local Variables:
// c-basic-offset: 4
// End:
