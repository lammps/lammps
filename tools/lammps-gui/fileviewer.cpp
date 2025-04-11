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
#include <QEvent>
#include <QFile>
#include <QFileInfo>
#include <QIcon>
#include <QKeySequence>
#include <QProcess>
#include <QSettings>
#include <QShortcut>
#include <QString>
#include <QStringList>
#include <QTextCursor>
#include <QTextStream>

FileViewer::FileViewer(const QString &_filename, QString title, QWidget *parent) :
    QPlainTextEdit(parent), fileName(_filename)
{
    auto *action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_Q), this);
    connect(action, &QShortcut::activated, this, &FileViewer::quit);
    action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_Slash), this);
    connect(action, &QShortcut::activated, this, &FileViewer::stop_run);

    installEventFilter(this);

    // open and read file. Set editor to read-only.
    QFile file(fileName);
    QFileInfo finfo(file);
    QString command;
    QString content;
    QProcess decomp;
    QStringList args = {"-cdf", fileName};
    bool compressed  = false;

    // match suffix with decompression program
    if (finfo.suffix() == "gz") {
        command    = "gzip";
        compressed = true;
    } else if (finfo.suffix() == "bz2") {
        command    = "bzip2";
        compressed = true;
    } else if (finfo.suffix() == "zst") {
        command    = "zstd";
        compressed = true;
    } else if (finfo.suffix() == "xz") {
        command    = "xz";
        compressed = true;
    } else if (finfo.suffix() == "lzma") {
        command = "xz";
        args.insert(1, "--format=lzma");
        compressed = true;
    } else if (finfo.suffix() == "lz4") {
        command    = "lz4";
        compressed = true;
    }

    // read compressed file from pipe
    if (compressed) {
        decomp.start(command, args, QIODevice::ReadOnly);
        if (decomp.waitForStarted()) {
            while (decomp.waitForReadyRead())
                content += decomp.readAll();
        } else {
            content = "\nCould not open compressed file %1 with decompression program %2\n";
            content = content.arg(fileName).arg(command);
        }
        decomp.close();
    } else if (file.open(QIODevice::Text | QIODevice::ReadOnly)) {
        // read plain text
        QTextStream in(&file);
        content = in.readAll();
        file.close();
    }

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
    if (title.isEmpty())
        setWindowTitle("LAMMPS-GUI - Viewer - " + fileName);
    else
        setWindowTitle(title);
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
