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

#include "flagwarnings.h"
#include "lammpsgui.h"

#include <QAction>
#include <QApplication>
#include <QFile>
#include <QFileDialog>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QIcon>
#include <QKeySequence>
#include <QLabel>
#include <QMenu>
#include <QMessageBox>
#include <QPushButton>
#include <QRegularExpression>
#include <QSettings>
#include <QShortcut>
#include <QSpacerItem>
#include <QString>
#include <QTextStream>

const QString LogWindow::yaml_regex =
    QStringLiteral("^(keywords:.*$|data:$|---$|\\.\\.\\.$|  - \\[.*\\]$)");

LogWindow::LogWindow(const QString &_filename, QWidget *parent) :
    QPlainTextEdit(parent), filename(_filename), warnings(nullptr)
{
    QSettings settings;
    resize(settings.value("logx", 500).toInt(), settings.value("logy", 320).toInt());

    summary = new QLabel("0 Warnings / Errors  -  0 Lines");
    summary->setMargin(1);

    auto *frame = new QFrame;
    frame->setAutoFillBackground(true);
    frame->setFrameStyle(QFrame::Box | QFrame::Plain);
    frame->setLineWidth(2);

    auto *button = new QPushButton(QIcon(":/icons/warning.png"), "");
    button->setToolTip("Jump to next warning");
    connect(button, &QPushButton::released, this, &LogWindow::next_warning);

    auto *spacer = new QSpacerItem(0, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    auto *panel  = new QHBoxLayout(frame);
    auto *grid   = new QGridLayout(this);

    panel->addWidget(summary);
    panel->addWidget(button);
    panel->setStretchFactor(summary, 10);
    panel->setStretchFactor(button, 1);

    grid->addItem(spacer, 0, 0, 1, 3);
    grid->addWidget(frame, 1, 1, 1, 1);
    grid->setColumnStretch(0, 5);
    grid->setColumnStretch(1, 1);
    grid->setColumnStretch(2, 5);

    warnings = new FlagWarnings(summary, document());

    auto *action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_S), this);
    connect(action, &QShortcut::activated, this, &LogWindow::save_as);
    action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_Y), this);
    connect(action, &QShortcut::activated, this, &LogWindow::extract_yaml);
    action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_Q), this);
    connect(action, &QShortcut::activated, this, &LogWindow::quit);
    action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_N), this);
    connect(action, &QShortcut::activated, this, &LogWindow::next_warning);
    action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_Slash), this);
    connect(action, &QShortcut::activated, this, &LogWindow::stop_run);

    installEventFilter(this);
}

LogWindow::~LogWindow()
{
    delete warnings;
    delete summary;
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

void LogWindow::next_warning()
{
#if QT_VERSION < QT_VERSION_CHECK(5, 15, 0)
    auto regex = QRegExp(QStringLiteral("^(ERROR|WARNING).*$"));
#else
    auto regex = QRegularExpression(QStringLiteral("^(ERROR|WARNING).*$"));
#endif

    if (warnings->get_nwarnings() > 0) {
        // wrap around search
        if (!find(regex)) {
            moveCursor(QTextCursor::Start, QTextCursor::MoveAnchor);
            find(regex);
        }
        // move cursor to unselect
        moveCursor(QTextCursor::NextBlock, QTextCursor::MoveAnchor);
    }
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
    action = menu->addAction("&Jump to next warning or error", this, &LogWindow::next_warning);
    action->setIcon(QIcon(":/icons/warning.png"));
    action->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_N));
    menu->addSeparator();
    action = menu->addAction("&Close Window", this, &QWidget::close);
    action->setIcon(QIcon(":/icons/window-close.png"));
    action->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_W));
    action = menu->addAction("&Quit LAMMPS-GUI", this, &LogWindow::quit);
    action->setIcon(QIcon(":/icons/application-exit.png"));
    action->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_Q));
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
