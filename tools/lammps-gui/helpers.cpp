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

#include "helpers.h"

#include <QApplication>
#include <QBrush>
#include <QColor>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QPalette>
#include <QProcess>
#include <QStringList>
#include <QWidget>

// duplicate string, STL version
char *mystrdup(const std::string &text)
{
    auto *tmp = new char[text.size() + 1];
    memcpy(tmp, text.c_str(), text.size() + 1);
    return tmp;
}

// duplicate string, pointer version
char *mystrdup(const char *text)
{
    return mystrdup(std::string(text));
}

// duplicate string, Qt version
char *mystrdup(const QString &text)
{
    return mystrdup(text.toStdString());
}

// Convert string into words on whitespace while handling single and double
// quotes. Adapted from LAMMPS_NS::utils::split_words() to preserve quotes.

std::vector<std::string> split_line(const std::string &text)
{
    std::vector<std::string> list;
    const char *buf = text.c_str();
    std::size_t beg = 0;
    std::size_t len = 0;
    std::size_t add = 0;

    char c = *buf;
    while (c) { // leading whitespace
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') {
            c = *++buf;
            ++beg;
            continue;
        };
        len = 0;

    // handle escaped/quoted text.
    quoted:

        if (c == '\'') { // handle single quote
            add = 0;
            len = 1;
            c   = *++buf;
            while (((c != '\'') && (c != '\0')) || ((c == '\\') && (buf[1] == '\''))) {
                if ((c == '\\') && (buf[1] == '\'')) {
                    ++buf;
                    ++len;
                }
                c = *++buf;
                ++len;
            }
            ++len;
            c = *++buf;

            // handle triple double quotation marks
        } else if ((c == '"') && (buf[1] == '"') && (buf[2] == '"') && (buf[3] != '"')) {
            len = 3;
            add = 1;
            buf += 3;
            c = *buf;

        } else if (c == '"') { // handle double quote
            add = 0;
            len = 1;
            c   = *++buf;
            while (((c != '"') && (c != '\0')) || ((c == '\\') && (buf[1] == '"'))) {
                if ((c == '\\') && (buf[1] == '"')) {
                    ++buf;
                    ++len;
                }
                c = *++buf;
                ++len;
            }
            ++len;
            c = *++buf;
        }

        while (true) { // unquoted
            if ((c == '\'') || (c == '"')) goto quoted;
            // skip escaped quote
            if ((c == '\\') && ((buf[1] == '\'') || (buf[1] == '"'))) {
                ++buf;
                ++len;
                c = *++buf;
                ++len;
            }
            if ((c == ' ') || (c == '\t') || (c == '\r') || (c == '\n') || (c == '\f') ||
                (c == '\0')) {
                if (beg < text.size()) list.push_back(text.substr(beg, len));
                beg += len + add;
                break;
            }
            c = *++buf;
            ++len;
        }
    }
    return list;
}

// get pointer to LAMMPS-GUI main widget

QWidget *get_main_widget()
{
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") return widget;
    return nullptr;
}

// find if executable is in path
// https://stackoverflow.com/a/51041497

bool has_exe(const QString &exe)
{
    QProcess findProcess;
    QStringList arguments;
    arguments << exe;
#if defined(_WIN32)
    findProcess.start("where", arguments);
#else
    findProcess.start("which", arguments);
#endif
    findProcess.setReadChannel(QProcess::ProcessChannel::StandardOutput);

    if (!findProcess.waitForFinished()) return false; // Not found or which does not work

    QString retStr(findProcess.readAll());
    retStr = retStr.trimmed();

    QFile file(retStr);
    QFileInfo check_file(file);
    return (check_file.exists() && check_file.isFile());
}

// recursively remove all contents from a directory

void purge_directory(const QString &dir)
{
    QDir directory(dir);

    directory.setFilter(QDir::AllEntries | QDir::NoDotAndDotDot);
    const auto &entries = directory.entryList();
    for (const auto &entry : entries) {
        if (!directory.remove(entry)) {
            directory.cd(entry);
            directory.removeRecursively();
            directory.cdUp();
        }
    }
}

// compare black level of foreground and background color
bool is_light_theme()
{
    QPalette p;
    int fg = p.brush(QPalette::Active, QPalette::WindowText).color().black();
    int bg = p.brush(QPalette::Active, QPalette::Window).color().black();

    return (fg > bg);
}

// Local Variables:
// c-basic-offset: 4
// End:
