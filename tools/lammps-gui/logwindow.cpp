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
#include <QFile>
#include <QFileDialog>
#include <QIcon>
#include <QKeySequence>
#include <QMenu>
#include <QMessageBox>
#include <QRegularExpression>
#include <QSettings>
#include <QShortcut>
#include <QString>
#include <QTextStream>

const QString LogWindow::yaml_regex =
    QStringLiteral("^(keywords:.*$|data:$|---$|\\.\\.\\.$|  - \\[.*\\]$)");

LogWindow::LogWindow(const QString &_filename, QWidget *parent) :
    QPlainTextEdit(parent), filename(_filename)
{
    QSettings settings;
    resize(settings.value("logx", 500).toInt(), settings.value("logy", 320).toInt());

    auto *action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_S), this);
    connect(action, &QShortcut::activated, this, &LogWindow::save_as);
    action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_Y), this);
    connect(action, &QShortcut::activated, this, &LogWindow::extract_yaml);
    action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_Q), this);
    connect(action, &QShortcut::activated, this, &LogWindow::quit);
    action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_Slash), this);
    connect(action, &QShortcut::activated, this, &LogWindow::stop_run);

    installEventFilter(this);
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
    LammpsGui *main = nullptr;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") main = dynamic_cast<LammpsGui *>(widget);
    if (main) main->quit();
}

void LogWindow::stop_run()
{
    LammpsGui *main = nullptr;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") main = dynamic_cast<LammpsGui *>(widget);
    if (main) main->stop_run();
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

bool LogWindow::check_yaml()
{
    QRegularExpression is_yaml(yaml_regex);
    QStringList lines = toPlainText().split('\n');
    for (const auto &line : lines)
        if (is_yaml.match(line).hasMatch()) return true;
    return false;
}

void LogWindow::extract_yaml()
{
    // ignore if no YAML format lines in buffer
    if (!check_yaml()) return;

    QString defaultname = filename + ".yaml";
    if (filename.isEmpty()) defaultname = "lammps.yaml";
    QString yamlFileName = QFileDialog::getSaveFileName(this, "Save YAML data to File", defaultname,
                                                        "YAML files (*.yaml *.yml)");
    // cannot save without filename
    if (yamlFileName.isEmpty()) return;

    QFileInfo path(yamlFileName);
    QFile file(path.absoluteFilePath());
    if (!file.open(QIODevice::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, "Warning", "Cannot save file: " + file.errorString());
        return;
    }

    QRegularExpression is_yaml(yaml_regex);
    QTextStream out(&file);
    QStringList lines = toPlainText().split('\n');
    for (const auto &line : lines) {
        if (is_yaml.match(line).hasMatch()) out << line << '\n';
    }
    file.close();
}

void LogWindow::contextMenuEvent(QContextMenuEvent *event)
{
    // show augmented context menu
    auto *menu = createStandardContextMenu();
    menu->addSeparator();
    auto *action = menu->addAction(QString("Save Log to File ..."));
    action->setIcon(QIcon(":/icons/document-save-as.png"));
    action->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_S));
    connect(action, &QAction::triggered, this, &LogWindow::save_as);
    // only show export-to-yaml entry if there is YAML format content.
    if (check_yaml()) {
        action = menu->addAction(QString("&Export YAML Data to File ..."));
        action->setIcon(QIcon(":/icons/yaml-file-icon.png"));
        action->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_Y));
        connect(action, &QAction::triggered, this, &LogWindow::extract_yaml);
    }
    action = menu->addAction("&Close Window", this, &QWidget::close);
    action->setIcon(QIcon(":/icons/window-close.png"));
    action->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_W));
    menu->exec(event->globalPos());
    delete menu;
}

// event filter to handle "Ambiguous shortcut override" issues
bool LogWindow::eventFilter(QObject *watched, QEvent *event)
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
