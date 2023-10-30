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

#include "codeeditor.h"
#include "lammpsgui.h"
#include "linenumberarea.h"

#include <QAbstractItemView>
#include <QAction>
#include <QCompleter>
#include <QDesktopServices>
#include <QDir>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QFileInfo>
#include <QIcon>
#include <QKeySequence>
#include <QMenu>
#include <QMimeData>
#include <QPainter>
#include <QRegularExpression>
#include <QScrollBar>
#include <QSettings>
#include <QShortcut>
#include <QStringListModel>
#include <QTextBlock>
#include <QTextDocumentFragment>
#include <QUrl>

#include <string>
#include <vector>

// Convert string into words on whitespace while handling single and double
// quotes. Adapted from LAMMPS_NS::utils::split_words() to preserve quotes.

static std::vector<std::string> split_line(const std::string &text)
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
                list.push_back(text.substr(beg, len));
                beg += len + add;
                break;
            }
            c = *++buf;
            ++len;
        }
    }
    return list;
}

CodeEditor::CodeEditor(QWidget *parent) :
    QPlainTextEdit(parent), current_comp(nullptr), command_comp(new QCompleter(this)),
    fix_comp(new QCompleter(this)), compute_comp(new QCompleter(this)),
    dump_comp(new QCompleter(this)), atom_comp(new QCompleter(this)),
    pair_comp(new QCompleter(this)), bond_comp(new QCompleter(this)),
    angle_comp(new QCompleter(this)), dihedral_comp(new QCompleter(this)),
    improper_comp(new QCompleter(this)), kspace_comp(new QCompleter(this)),
    region_comp(new QCompleter(this)), integrate_comp(new QCompleter(this)),
    minimize_comp(new QCompleter(this)), variable_comp(new QCompleter(this)),
    units_comp(new QCompleter(this)), group_comp(new QCompleter(this)),
    varname_comp(new QCompleter(this)), fixid_comp(new QCompleter(this)),
    compid_comp(new QCompleter(this)), file_comp(new QCompleter(this)), highlight(NO_HIGHLIGHT)
{
    help_action = new QShortcut(QKeySequence::fromString("Ctrl+?"), parent);
    connect(help_action, &QShortcut::activated, this, &CodeEditor::get_help);

    // set up completer class (without a model currently)
#define COMPLETER_SETUP(completer)                                                            \
    completer->setCompletionMode(QCompleter::UnfilteredPopupCompletion);                      \
    completer->setModelSorting(QCompleter::CaseInsensitivelySortedModel);                     \
    completer->setWidget(this);                                                               \
    completer->setMaxVisibleItems(16);                                                        \
    completer->setWrapAround(false);                                                          \
    QObject::connect(completer, QOverload<const QString &>::of(&QCompleter::activated), this, \
                     &CodeEditor::insertCompletedCommand)

    COMPLETER_SETUP(command_comp);
    COMPLETER_SETUP(fix_comp);
    COMPLETER_SETUP(compute_comp);
    COMPLETER_SETUP(dump_comp);
    COMPLETER_SETUP(atom_comp);
    COMPLETER_SETUP(pair_comp);
    COMPLETER_SETUP(bond_comp);
    COMPLETER_SETUP(angle_comp);
    COMPLETER_SETUP(dihedral_comp);
    COMPLETER_SETUP(improper_comp);
    COMPLETER_SETUP(kspace_comp);
    COMPLETER_SETUP(region_comp);
    COMPLETER_SETUP(integrate_comp);
    COMPLETER_SETUP(minimize_comp);
    COMPLETER_SETUP(variable_comp);
    COMPLETER_SETUP(units_comp);
    COMPLETER_SETUP(group_comp);
    COMPLETER_SETUP(varname_comp);
    COMPLETER_SETUP(fixid_comp);
    COMPLETER_SETUP(compid_comp);
    COMPLETER_SETUP(file_comp);
#undef COMPLETER_SETUP

    // initialize help system
    QFile help_index(":/help_index.table");
    if (help_index.open(QIODevice::ReadOnly | QIODevice::Text)) {
        while (!help_index.atEnd()) {
            auto line  = QString(help_index.readLine());
            auto words = line.trimmed().replace('\t', ' ').split(' ');
            if (words.size() > 2) {
                if (words.at(1) == "pair_style") {
                    pair_map[words.at(2)] = words.at(0);
                } else if (words.at(1) == "bond_style") {
                    bond_map[words.at(2)] = words.at(0);
                } else if (words.at(1) == "angle_style") {
                    angle_map[words.at(2)] = words.at(0);
                } else if (words.at(1) == "dihedral_style") {
                    dihedral_map[words.at(2)] = words.at(0);
                } else if (words.at(1) == "improper_style") {
                    improper_map[words.at(2)] = words.at(0);
                } else if (words.at(1) == "fix") {
                    fix_map[words.at(2)] = words.at(0);
                } else if (words.at(1) == "compute") {
                    compute_map[words.at(2)] = words.at(0);
                } else if (words.at(1) == "kspace_style") {
                    cmd_map["kspace_style"] = "kspace_style.html";
                }
                // ignoring: dump, fix_modify ATC
            } else if (words.size() == 2) {
                cmd_map[words.at(1)] = words.at(0);
            } else {
                fprintf(stderr, "unhandled: %s", line.toStdString().c_str());
            }
        }
        help_index.close();
    }

    lineNumberArea = new LineNumberArea(this);
    connect(this, &CodeEditor::blockCountChanged, this, &CodeEditor::updateLineNumberAreaWidth);
    connect(this, &CodeEditor::updateRequest, this, &CodeEditor::updateLineNumberArea);
    updateLineNumberAreaWidth(0);
    setCursorWidth(2);
}

CodeEditor::~CodeEditor()
{
    delete help_action;
    delete lineNumberArea;

    delete command_comp;
    delete fix_comp;
    delete compute_comp;
    delete atom_comp;
    delete pair_comp;
    delete bond_comp;
    delete angle_comp;
    delete dihedral_comp;
    delete improper_comp;
    delete kspace_comp;
    delete region_comp;
    delete integrate_comp;
    delete minimize_comp;
    delete variable_comp;
    delete units_comp;
    delete group_comp;
    delete varname_comp;
    delete fixid_comp;
    delete compid_comp;
    delete file_comp;
}

int CodeEditor::lineNumberAreaWidth()
{
    int digits = 1;
    int max    = qMax(1, blockCount());
    while (max >= 10) {
        max /= 10;
        ++digits;
    }

    int space = 3 + fontMetrics().horizontalAdvance(QLatin1Char('9')) * (digits + 2);
    return space;
}

void CodeEditor::setFont(const QFont &newfont)
{
    lineNumberArea->setFont(newfont);
    document()->setDefaultFont(newfont);
}

void CodeEditor::setCursor(int block)
{
    // move cursor to given position
    auto cursor = textCursor();
    int moves   = block - cursor.blockNumber();
    if (moves < 0)
        cursor.movePosition(QTextCursor::Up, QTextCursor::MoveAnchor, -moves);
    else
        cursor.movePosition(QTextCursor::Down, QTextCursor::MoveAnchor, moves);
    setTextCursor(cursor);
}

void CodeEditor::setHighlight(int block, bool error)
{
    if (error)
        highlight = -block;
    else
        highlight = block;

    // also reset the cursor
    setCursor(block);

    // update graphics
    repaint();
}

// reformat line

QString CodeEditor::reformatLine(const QString &line)
{
    auto words = split_line(line.toStdString());
    QString newtext;
    QSettings settings;
    settings.beginGroup("reformat");
    int cmdsize  = settings.value("command", "16").toInt();
    int typesize = settings.value("type", "4").toInt();
    int idsize   = settings.value("id", "4").toInt();
    int namesize = settings.value("name", "8").toInt();
    settings.endGroup();

    bool rebuildGroupComp     = false;
    bool rebuildVarNameComp   = false;
    bool rebuildComputeIDComp = false;
    bool rebuildFixIDComp     = false;

    if (words.size()) {
        // commented line. do nothing
        if (words[0][0] == '#') return line;

        // start with LAMMPS command plus padding if another word follows
        newtext = words[0].c_str();
        if (words.size() > 1) {
            for (int i = words[0].size() + 1; i < cmdsize; ++i)
                newtext += ' ';
            // new/updated group command -> update completer
            if (words[0] == "group") rebuildGroupComp = true;
            // new/updated variable command -> update completer
            if (words[0] == "variable") rebuildVarNameComp = true;
            // new/updated compute command -> update completer
            if (words[0] == "compute") rebuildComputeIDComp = true;
            // new/updated fix command -> update completer
            if (words[0] == "fix") rebuildFixIDComp = true;
        }

        // append remaining words with just a single blank added.
        for (std::size_t i = 1; i < words.size(); ++i) {
            newtext += ' ';
            newtext += words[i].c_str();

            // special cases

            if (i < 3) {
                // additional space for types or type ranges
                if (words[0] == "pair_coeff")
                    for (int j = words[i].size(); j < typesize; ++j)
                        newtext += ' ';

                // pad 4 for IDs and 8 for groups
                if ((words[0] == "fix") || (words[0] == "compute") || (words[0] == "dump")) {
                    if (i == 1) {
                        for (int j = words[i].size(); j < idsize; ++j)
                            newtext += ' ';
                    } else if (i == 2) {
                        for (int j = words[i].size(); j < namesize; ++j)
                            newtext += ' ';
                    }
                }
            }

            if (i < 2) {
                if ((words[0] == "bond_coeff") || (words[0] == "angle_coeff") ||
                    (words[0] == "dihedral_coeff") || (words[0] == "improper_coeff") ||
                    (words[0] == "mass"))
                    for (int j = words[i].size(); j < typesize; ++j)
                        newtext += ' ';
            }
        }
    }
    if (rebuildGroupComp) setGroupList();
    if (rebuildVarNameComp) setVarNameList();
    if (rebuildComputeIDComp) setComputeIDList();
    if (rebuildFixIDComp) setFixIDList();
    return newtext;
}

#define COMPLETER_INIT_FUNC(keyword, Type)                                     \
    void CodeEditor::set##Type##List(const QStringList &words)                 \
    {                                                                          \
        keyword##_comp->setModel(new QStringListModel(words, keyword##_comp)); \
    }

COMPLETER_INIT_FUNC(command, Command)
COMPLETER_INIT_FUNC(fix, Fix)
COMPLETER_INIT_FUNC(compute, Compute)
COMPLETER_INIT_FUNC(dump, Dump)
COMPLETER_INIT_FUNC(atom, Atom)
COMPLETER_INIT_FUNC(pair, Pair)
COMPLETER_INIT_FUNC(bond, Bond)
COMPLETER_INIT_FUNC(angle, Angle)
COMPLETER_INIT_FUNC(dihedral, Dihedral)
COMPLETER_INIT_FUNC(improper, Improper)
COMPLETER_INIT_FUNC(kspace, Kspace)
COMPLETER_INIT_FUNC(region, Region)
COMPLETER_INIT_FUNC(integrate, Integrate)
COMPLETER_INIT_FUNC(minimize, Minimize)
COMPLETER_INIT_FUNC(variable, Variable)
COMPLETER_INIT_FUNC(units, Units)

#undef COMPLETER_INIT_FUNC

// build completer for groups by parsing through edit buffer

void CodeEditor::setGroupList()
{
    QStringList groups;
#if QT_VERSION < QT_VERSION_CHECK(5, 15, 0)
    QRegExp groupcmd(QStringLiteral("^\\s*group\\s+(\\S+)(\\s+|$)"));
#else
    QRegularExpression groupcmd(QStringLiteral("^\\s*group\\s+(\\S+)(\\s+|$)"));
#endif
    auto saved = textCursor();
    // reposition cursor to beginning of text and search for group commands
    auto cursor = textCursor();
    cursor.movePosition(QTextCursor::Start);
    setTextCursor(cursor);
    while (find(groupcmd)) {
        auto words = split_line(textCursor().block().text().replace('\t', ' ').toStdString());
        if ((words.size() > 1) && !groups.contains(words[1].c_str())) groups << words[1].c_str();
    }
    groups.sort();
    groups.prepend(QStringLiteral("all"));

    setTextCursor(saved);
    group_comp->setModel(new QStringListModel(groups, group_comp));
}

void CodeEditor::setVarNameList()
{
    QStringList vars;

    // variable "gui_run" is always defined by LAMMPS GUI
    vars << QString("${gui_run}");
    vars << QString("v_gui_run");

    LammpsWrapper *lammps = &qobject_cast<LammpsGui *>(parent())->lammps;
    int nvar              = lammps->id_count("variable");
    constexpr int BUFLEN  = 256;
    char buffer[BUFLEN];
    for (int i = 0; i < nvar; ++i) {
        memset(buffer, 0, BUFLEN);
        if (lammps->variable_info(i, buffer, BUFLEN)) {
            if (strlen(buffer) == 1) vars << QString("$%1").arg(buffer);
            vars << QString("${%1}").arg(buffer);
            vars << QString("v_%1").arg(buffer);
        }
    }

#if QT_VERSION < QT_VERSION_CHECK(5, 15, 0)
    QRegExp varcmd(QStringLiteral("^\\s*variable\\s+(\\S+)(\\s+|$)"));
#else
    QRegularExpression varcmd(QStringLiteral("^\\s*variable\\s+(\\S+)(\\s+|$)"));
#endif
    auto saved = textCursor();
    // reposition cursor to beginning of text and search for group commands
    auto cursor = textCursor();
    cursor.movePosition(QTextCursor::Start);
    setTextCursor(cursor);
    while (find(varcmd)) {
        auto words = split_line(textCursor().block().text().replace('\t', ' ').toStdString());
        if ((words.size() > 1)) {
            QString w = QString("$%1").arg(words[1].c_str());
            if ((words[1].size() == 1) && !vars.contains(w)) vars << w;
            w = QString("${%1}").arg(words[1].c_str());
            if (!vars.contains(w)) vars << w;
            w = QString("v_%1").arg(words[1].c_str());
            if (!vars.contains(w)) vars << w;
        }
    }
    vars.sort();

    setTextCursor(saved);
    varname_comp->setModel(new QStringListModel(vars, varname_comp));
}

void CodeEditor::setComputeIDList()
{
    QStringList compid;
#if QT_VERSION < QT_VERSION_CHECK(5, 15, 0)
    QRegExp compcmd(QStringLiteral("^\\s*compute\\s+(\\S+)\\s+"));
#else
    QRegularExpression compcmd(QStringLiteral("^\\s*compute\\s+(\\S+)\\s+"));
#endif
    auto saved = textCursor();
    // reposition cursor to beginning of text and search for group commands
    auto cursor = textCursor();
    cursor.movePosition(QTextCursor::Start);
    setTextCursor(cursor);
    while (find(compcmd)) {
        auto words = split_line(textCursor().block().text().replace('\t', ' ').toStdString());
        if ((words.size() > 1)) {
            QString w = QString("c_%1").arg(words[1].c_str());
            if (!compid.contains(w)) compid << w;
            w = QString("C_%1").arg(words[1].c_str());
            if (!compid.contains(w)) compid << w;
        }
    }
    compid.sort();

    setTextCursor(saved);
    compid_comp->setModel(new QStringListModel(compid, compid_comp));
}

void CodeEditor::setFixIDList()
{
    QStringList fixid;
#if QT_VERSION < QT_VERSION_CHECK(5, 15, 0)
    QRegExp fixcmd(QStringLiteral("^\\s*fix\\s+(\\S+)\\s+"));
#else
    QRegularExpression fixcmd(QStringLiteral("^\\s*fix\\s+(\\S+)\\s+"));
#endif
    auto saved = textCursor();
    // reposition cursor to beginning of text and search for group commands
    auto cursor = textCursor();
    cursor.movePosition(QTextCursor::Start);
    setTextCursor(cursor);
    while (find(fixcmd)) {
        auto words = split_line(textCursor().block().text().replace('\t', ' ').toStdString());
        if ((words.size() > 1)) {
            QString w = QString("f_%1").arg(words[1].c_str());
            if (!fixid.contains(w)) fixid << w;
            w = QString("F_%1").arg(words[1].c_str());
            if (!fixid.contains(w)) fixid << w;
        }
    }
    fixid.sort();

    setTextCursor(saved);
    fixid_comp->setModel(new QStringListModel(fixid, fixid_comp));
}

void CodeEditor::setFileList()
{
    QStringList files;
    QDir dir(".");
    for (const auto &file : dir.entryInfoList(QDir::Files))
        files << file.fileName();
    files.sort();
    file_comp->setModel(new QStringListModel(files, file_comp));
}

void CodeEditor::keyPressEvent(QKeyEvent *event)
{
    const auto key = event->key();

    if (current_comp && current_comp->popup()->isVisible()) {
        // The following keys are forwarded by the completer to the widget
        switch (key) {
            case Qt::Key_Enter:
            case Qt::Key_Return:
            case Qt::Key_Escape:
            case Qt::Key_Tab:
            case Qt::Key_Backtab:
                event->ignore();
                return; // let the completer do default behavior
            default:
                break;
        }
    }

    // reformat current line and consume key event
    if (key == Qt::Key_Tab) {
        reformatCurrentLine();
        return;
    }

    // run command completion and consume key event
    if (key == Qt::Key_Backtab) {
        runCompletion();
        return;
    }

    // automatically reformat when hitting the return or enter key
    if (reformat_on_return && ((key == Qt::Key_Return) || (key == Qt::Key_Enter))) {
        reformatCurrentLine();
    }

    // process key event in parent class
    QPlainTextEdit::keyPressEvent(event);

    // if enabled, try pop up completion automatically after 2 characters
    if (automatic_completion) {
        auto cursor = textCursor();
        auto line   = cursor.block().text();
        if (line.isEmpty()) return;

        // QTextCursor::WordUnderCursor is unusable here since recognizes '/' as word boundary.
        // Work around it by manually searching for the location of the beginning of the word.
        int begin = qMin(cursor.positionInBlock(), line.length() - 1);

        while (begin >= 0) {
            if (line[begin].isSpace()) break;
            --begin;
        }
        if (((cursor.positionInBlock() - begin) > 2) || (line[begin + 1] == '$')) runCompletion();
        if (current_comp && current_comp->popup()->isVisible() &&
            ((cursor.positionInBlock() - begin) < 2)) {
            current_comp->popup()->hide();
        }
    }
}

void CodeEditor::updateLineNumberAreaWidth(int /* newBlockCount */)
{
    setViewportMargins(lineNumberAreaWidth(), 0, 0, 0);
}

void CodeEditor::updateLineNumberArea(const QRect &rect, int dy)
{
    if (dy)
        lineNumberArea->scroll(0, dy);
    else
        lineNumberArea->update(0, rect.y(), lineNumberArea->width(), rect.height());

    if (rect.contains(viewport()->rect())) updateLineNumberAreaWidth(0);
}

void CodeEditor::dragEnterEvent(QDragEnterEvent *event)
{
    event->acceptProposedAction();
}

bool CodeEditor::canInsertFromMimeData(const QMimeData *source) const
{
    return source->hasUrls(); // || source->hasText();
}

void CodeEditor::dropEvent(QDropEvent *event)
{
    if (event->mimeData()->hasUrls()) {
        event->accept();
        auto file = event->mimeData()->urls()[0].toLocalFile();
        auto gui  = dynamic_cast<LammpsGui *>(parent());
        if (gui) {
            moveCursor(QTextCursor::Start, QTextCursor::MoveAnchor);
            gui->open_file(file);
        }
    } else if (event->mimeData()->hasText()) {
        event->accept();
        fprintf(stderr, "Drag - Drop for text block not yet implemented: text=%s\n",
                event->mimeData()->text().toStdString().c_str());
    } else
        event->ignore();
}

void CodeEditor::resizeEvent(QResizeEvent *e)
{
    QPlainTextEdit::resizeEvent(e);

    QRect cr = contentsRect();
    lineNumberArea->setGeometry(QRect(cr.left(), cr.top(), lineNumberAreaWidth(), cr.height()));
}

void CodeEditor::lineNumberAreaPaintEvent(QPaintEvent *event)
{
    QPainter painter(lineNumberArea);
    QTextBlock block = firstVisibleBlock();
    int blockNumber  = block.blockNumber();

    int top    = qRound(blockBoundingGeometry(block).translated(contentOffset()).top());
    int bottom = top + qRound(blockBoundingRect(block).height());
    while (block.isValid() && top <= event->rect().bottom()) {
        if (block.isVisible() && bottom >= event->rect().top()) {
            QString number = QString::number(blockNumber + 1) + " ";
            if ((highlight == NO_HIGHLIGHT) || (blockNumber != std::abs(highlight))) {
                painter.setPen(Qt::black);
            } else {
                number = QString(">") + QString::number(blockNumber + 1) + "<";
                if (highlight < 0)
                    painter.fillRect(0, top, lineNumberArea->width(), fontMetrics().height(),
                                     Qt::darkRed);
                else
                    painter.fillRect(0, top, lineNumberArea->width(), fontMetrics().height(),
                                     Qt::darkGreen);

                painter.setPen(Qt::white);
            }
            painter.drawText(0, top, lineNumberArea->width(), fontMetrics().height(),
                             Qt::AlignRight, number);
        }

        block  = block.next();
        top    = bottom;
        bottom = top + qRound(blockBoundingRect(block).height());
        ++blockNumber;
    }
}

void CodeEditor::contextMenuEvent(QContextMenuEvent *event)
{
    // reposition the cursor here, but only if there is no active selection
    if (!textCursor().hasSelection()) setTextCursor(cursorForPosition(event->pos()));

    QString page, help;
    find_help(page, help);

    auto *menu = createStandardContextMenu();
    menu->addSeparator();
    if (textCursor().hasSelection()) {
        auto action1 = menu->addAction("Comment out selection");
        action1->setIcon(QIcon(":/icons/expand-text.png"));
        connect(action1, &QAction::triggered, this, &CodeEditor::comment_selection);
        auto action2 = menu->addAction("Uncomment selection");
        action2->setIcon(QIcon(":/icons/expand-text.png"));
        connect(action2, &QAction::triggered, this, &CodeEditor::uncomment_selection);
    } else {
        auto action1 = menu->addAction("Comment out line");
        action1->setIcon(QIcon(":/icons/expand-text.png"));
        connect(action1, &QAction::triggered, this, &CodeEditor::comment_line);
        auto action2 = menu->addAction("Uncomment line");
        action2->setIcon(QIcon(":/icons/expand-text.png"));
        connect(action2, &QAction::triggered, this, &CodeEditor::uncomment_line);
    }
    menu->addSeparator();

    // print augmented context menu if an entry was found
    if (!help.isEmpty()) {
        auto action = menu->addAction(QString("Display available completions for '%1'").arg(help));
        action->setIcon(QIcon(":/icons/expand-text.png"));
        connect(action, &QAction::triggered, this, &CodeEditor::runCompletion);
        menu->addSeparator();
    }

    if (!page.isEmpty()) {
        auto action = menu->addAction(QString("Reformat '%1' command").arg(help));
        action->setIcon(QIcon(":/icons/format-indent-less-3.png"));
        connect(action, &QAction::triggered, this, &CodeEditor::reformatCurrentLine);

        menu->addSeparator();
        action = menu->addAction(QString("View Documentation for '%1'").arg(help));
        action->setIcon(QIcon(":/icons/system-help.png"));
        action->setData(page);
        connect(action, &QAction::triggered, this, &CodeEditor::open_help);
        // if we link to help with specific styles (fix, compute, pair, bond, ...)
        // also link to the docs for the primary command
        auto words = help.split(' ');
        if (words.size() > 1) {
            help = words.at(0);
            page = words.at(0);
            page += ".html";
            auto action2 = menu->addAction(QString("View Documentation for '%1'").arg(help));
            action2->setIcon(QIcon(":/icons/system-help.png"));
            action2->setData(page);
            connect(action2, &QAction::triggered, this, &CodeEditor::open_help);
        }
    }
    auto action = menu->addAction(QString("LAMMPS Manual"));
    action->setIcon(QIcon(":/icons/help-browser.png"));
    action->setData(QString());
    connect(action, &QAction::triggered, this, &CodeEditor::open_help);

    menu->exec(event->globalPos());
    delete menu;
}

void CodeEditor::reformatCurrentLine()
{
    auto cursor  = textCursor();
    auto text    = cursor.block().text();
    auto newtext = reformatLine(text);

    // perform edit but only if text has changed
    if (QString::compare(text, newtext)) {
        cursor.beginEditBlock();
        cursor.movePosition(QTextCursor::StartOfLine);
        cursor.movePosition(QTextCursor::EndOfLine, QTextCursor::KeepAnchor, 1);
        cursor.insertText(newtext);
        cursor.endEditBlock();
    }
}

void CodeEditor::comment_line()
{
    auto cursor = textCursor();
    cursor.movePosition(QTextCursor::StartOfLine);
    cursor.insertText("#");
}

void CodeEditor::comment_selection()
{
    auto cursor = textCursor();
    auto text   = cursor.selection().toPlainText();
    auto lines  = text.split('\n');
    QString newtext;
    for (auto line : lines) {
        newtext.append('#');
        newtext.append(line);
        newtext.append('\n');
    }
    if (newtext.isEmpty()) newtext = "#\n";
    cursor.insertText(newtext);
    setTextCursor(cursor);
}

void CodeEditor::uncomment_selection()
{
    auto cursor = textCursor();
    auto text   = cursor.selection().toPlainText();
    auto lines  = text.split('\n');
    QString newtext;
    for (auto line : lines) {
        QString newline;
        bool start = true;
        for (auto letter : line) {
            if (start && (letter == '#')) {
                start = false;
                continue;
            }
            if (start && !letter.isSpace()) start = false;
            newline.append(letter);
        }
        newtext.append(newline);
        newtext.append('\n');
    }
    cursor.insertText(newtext);
    setTextCursor(cursor);
}

void CodeEditor::uncomment_line()
{
    auto cursor = textCursor();
    auto text   = cursor.block().text();
    QString newtext;
    bool start = true;
    for (auto letter : text) {
        if (start && (letter == '#')) {
            start = false;
            continue;
        }
        if (start && !letter.isSpace()) start = false;
        newtext.append(letter);
    }

    // perform edit but only if text has changed
    if (QString::compare(text, newtext)) {
        cursor.beginEditBlock();
        cursor.movePosition(QTextCursor::StartOfLine);
        cursor.movePosition(QTextCursor::EndOfLine, QTextCursor::KeepAnchor, 1);
        cursor.insertText(newtext);
        cursor.endEditBlock();
    }
}

void CodeEditor::runCompletion()
{
    QAbstractItemView *popup = nullptr;
    if (current_comp) popup = current_comp->popup();

    auto cursor = textCursor();
    auto line   = cursor.block().text().trimmed();
    // no completion possible on empty lines
    if (line.isEmpty()) return;
    auto words = split_line(line.toStdString());

    // QTextCursor::WordUnderCursor is unusable here since it recognizes '/' as word boundary.
    // Work around it by manually searching for the beginning and end position of the word
    // under the cursor and then using that substring.
    int begin = qMin(cursor.positionInBlock(), line.length() - 1);
    line      = cursor.block().text();
    while (begin >= 0) {
        if (line[begin].isSpace()) break;
        --begin;
    }
    int end = ++begin;
    while (end < line.length()) {
        if (line[end].isSpace()) break;
        ++end;
    }
    const auto selected = line.mid(begin, end - begin);

    // variable expansion may be anywhere
    if (selected.startsWith("$")) {
        current_comp = varname_comp;
        current_comp->setCompletionPrefix(selected);
        if (popup && (popup != current_comp->popup())) popup->hide();
        popup = current_comp->popup();
        // if the command is already a complete command, remove existing popup
        if (selected == current_comp->currentCompletion()) {
            if (popup->isVisible()) {
                popup->hide();
                current_comp = nullptr;
            }
            return;
        }
        QRect cr = cursorRect();
        cr.setWidth(popup->sizeHintForColumn(0) + popup->verticalScrollBar()->sizeHint().width());
        popup->setAlternatingRowColors(true);
        current_comp->complete(cr);

        // if on first word, try to complete command
    } else if ((words.size() > 0) && (words[0] == selected.toStdString())) {
        // no completion on comment lines
        if (words[0][0] == '#') return;

        current_comp = command_comp;
        current_comp->setCompletionPrefix(words[0].c_str());
        if (popup && (popup != current_comp->popup())) popup->hide();
        popup = current_comp->popup();
        // if the command is already a complete command, remove existing popup
        if (words[0] == current_comp->currentCompletion().toStdString()) {
            if (popup->isVisible()) {
                popup->hide();
                current_comp = nullptr;
            }
            return;
        }
        QRect cr = cursorRect();
        cr.setWidth(popup->sizeHintForColumn(0) + popup->verticalScrollBar()->sizeHint().width());
        popup->setAlternatingRowColors(true);
        current_comp->complete(cr);

        // completions for second word
    } else if ((words.size() > 1) && (words[1] == selected.toStdString())) {
        // no completion on comment lines
        if (words[0][0] == '#') return;

        current_comp = nullptr;
        if (words[0] == "pair_style")
            current_comp = pair_comp;
        else if (words[0] == "bond_style")
            current_comp = bond_comp;
        else if (words[0] == "angle_style")
            current_comp = angle_comp;
        else if (words[0] == "dihedral_style")
            current_comp = dihedral_comp;
        else if (words[0] == "improper_style")
            current_comp = improper_comp;
        else if (words[0] == "kspace_style")
            current_comp = kspace_comp;
        else if (words[0] == "atom_style")
            current_comp = atom_comp;
        else if (words[0] == "run_style")
            current_comp = integrate_comp;
        else if (words[0] == "minimize_style")
            current_comp = minimize_comp;
        else if (words[0] == "units")
            current_comp = units_comp;
        else if ((words[0] == "change_box") || (words[0] == "displace_atoms") ||
                 (words[0] == "velocity") || (words[0] == "write_dump"))
            current_comp = group_comp;
        else if ((words[0] == "fitpod") || (words[0] == "include") || (words[0] == "ndx2group") ||
                 (words[0] == "read_data") || (words[0] == "read_dump") ||
                 (words[0] == "read_restart") || (words[0] == "rerun")) {
            if (selected.contains('/')) {
                if (popup && popup->isVisible()) popup->hide();
            } else
                current_comp = file_comp;
        } else if (selected.startsWith("v_"))
            current_comp = varname_comp;
        else if (selected.startsWith("c_"))
            current_comp = compid_comp;
        else if (selected.startsWith("C_"))
            current_comp = compid_comp;
        else if (selected.startsWith("f_"))
            current_comp = fixid_comp;
        else if (selected.startsWith("F_"))
            current_comp = fixid_comp;

        if (current_comp) {
            current_comp->setCompletionPrefix(words[1].c_str());
            if (popup && (popup != current_comp->popup())) popup->hide();
            popup = current_comp->popup();
            // if the command is already a complete command, remove existing popup
            if (words[1] == current_comp->currentCompletion().toStdString()) {
                if (popup->isVisible()) popup->hide();
                return;
            }
            QRect cr = cursorRect();
            cr.setWidth(popup->sizeHintForColumn(0) +
                        popup->verticalScrollBar()->sizeHint().width());
            popup->setAlternatingRowColors(true);
            current_comp->complete(cr);
        }
        // completions for third word
    } else if ((words.size() > 2) && (words[2] == selected.toStdString())) {
        // no completion on comment lines
        if (words[0][0] == '#') return;

        current_comp = nullptr;
        if (words[0] == "region")
            current_comp = region_comp;
        else if (words[0] == "variable")
            current_comp = variable_comp;
        else if (words[0] == "fix")
            current_comp = group_comp;
        else if (words[0] == "compute")
            current_comp = group_comp;
        else if (words[0] == "dump")
            current_comp = group_comp;
        else if (selected.startsWith("v_"))
            current_comp = varname_comp;
        else if (selected.startsWith("c_"))
            current_comp = compid_comp;
        else if (selected.startsWith("C_"))
            current_comp = compid_comp;
        else if (selected.startsWith("f_"))
            current_comp = fixid_comp;
        else if (selected.startsWith("F_"))
            current_comp = fixid_comp;
        else if ((words[0] == "fitpod") || (words[0] == "molecule")) {
            if (selected.contains('/')) {
                if (popup && popup->isVisible()) popup->hide();
            } else
                current_comp = file_comp;
        }
        if (current_comp) {
            current_comp->setCompletionPrefix(words[2].c_str());
            if (popup && (popup != current_comp->popup())) popup->hide();
            popup = current_comp->popup();
            // if the command is already a complete command, remove existing popup
            if (words[2] == current_comp->currentCompletion().toStdString()) {
                if (popup->isVisible()) popup->hide();
                return;
            }
            QRect cr = cursorRect();
            cr.setWidth(popup->sizeHintForColumn(0) +
                        popup->verticalScrollBar()->sizeHint().width());
            popup->setAlternatingRowColors(true);
            current_comp->complete(cr);
        }
        // completions for fourth word
    } else if ((words.size() > 3) && (words[3] == selected.toStdString())) {
        // no completion on comment lines
        if (words[0][0] == '#') return;

        current_comp = nullptr;
        if (words[0] == "fix")
            current_comp = fix_comp;
        else if (words[0] == "compute")
            current_comp = compute_comp;
        else if (words[0] == "dump")
            current_comp = dump_comp;
        else if ((words[0] == "pair_coeff") && (words[1] == "*") && (words[2] == "*")) {
            if (selected.contains('/')) {
                if (popup && popup->isVisible()) popup->hide();
            } else
                current_comp = file_comp;
        } else if (selected.startsWith("v_"))
            current_comp = varname_comp;
        else if (selected.startsWith("c_"))
            current_comp = compid_comp;
        else if (selected.startsWith("C_"))
            current_comp = compid_comp;
        else if (selected.startsWith("f_"))
            current_comp = fixid_comp;
        else if (selected.startsWith("F_"))
            current_comp = fixid_comp;

        if (current_comp) {
            current_comp->setCompletionPrefix(words[3].c_str());
            if (popup && (popup != current_comp->popup())) popup->hide();
            popup = current_comp->popup();
            // if the command is already a complete command, remove existing popup
            if (words[3] == current_comp->currentCompletion().toStdString()) {
                if (popup->isVisible()) popup->hide();
                return;
            }
            QRect cr = cursorRect();
            cr.setWidth(popup->sizeHintForColumn(0) +
                        popup->verticalScrollBar()->sizeHint().width());
            popup->setAlternatingRowColors(true);
            current_comp->complete(cr);
        }
        // reference located anywhere further right in the line
    } else if (words.size() > 4) {
        current_comp = nullptr;
        if (selected.startsWith("v_"))
            current_comp = varname_comp;
        else if (selected.startsWith("c_"))
            current_comp = compid_comp;
        else if (selected.startsWith("C_"))
            current_comp = compid_comp;
        else if (selected.startsWith("f_"))
            current_comp = fixid_comp;
        else if (selected.startsWith("F_"))
            current_comp = fixid_comp;

        if (current_comp) {
            current_comp->setCompletionPrefix(selected);
            if (popup && (popup != current_comp->popup())) popup->hide();
            popup = current_comp->popup();
            // if the command is already a complete command, remove existing popup
            if (selected == current_comp->currentCompletion()) {
                if (popup->isVisible()) popup->hide();
                return;
            }
            QRect cr = cursorRect();
            cr.setWidth(popup->sizeHintForColumn(0) +
                        popup->verticalScrollBar()->sizeHint().width());
            popup->setAlternatingRowColors(true);
            current_comp->complete(cr);
        }
    }
}

void CodeEditor::insertCompletedCommand(const QString &completion)
{
    auto *completer = qobject_cast<QCompleter *>(sender());
    if (completer->widget() != this) return;

    // select the entire word (non-space text) under the cursor
    // we need to do it in this compicated way, since QTextCursor does not recognize
    // special characters as part of a word.
    auto cursor = textCursor();
    auto line   = cursor.block().text();
    int begin   = cursor.positionInBlock();
    do {
        if (line[begin].isSpace()) break;
        --begin;
    } while (begin >= 0);

    int end = begin + 1;
    while (end < line.length()) {
        if (line[end].isSpace()) break;
        ++end;
    }

    cursor.setPosition(cursor.position() - cursor.positionInBlock() + begin + 1);
    cursor.movePosition(QTextCursor::NextCharacter, QTextCursor::KeepAnchor, end - begin - 1);
    cursor.insertText(completion);
    setTextCursor(cursor);
}

void CodeEditor::get_help()
{
    QString page, help;
    find_help(page, help);
    if (!page.isEmpty())
        QDesktopServices::openUrl(QUrl(QString("https://docs.lammps.org/%1").arg(page)));
}

void CodeEditor::find_help(QString &page, QString &help)
{
    // process line of text where the cursor is
    auto text = textCursor().block().text().replace('\t', ' ').trimmed();
    auto style =
        QRegularExpression("^(pair|bond|angle|dihedral|improper)_style\\s+(\\S+)").match(text);
    help.clear();
    page.clear();
    if (style.hasMatch()) {
        if (style.captured(1) == "pair") {
            page = pair_map.value(style.captured(2), QString());
            help = QString("pair_style %1").arg(style.captured(2));
        } else if (style.captured(1) == "bond") {
            page = bond_map.value(style.captured(2), QString());
            help = QString("bond_style %1").arg(style.captured(2));
        } else if (style.captured(1) == "angle") {
            page = angle_map.value(style.captured(2), QString());
            help = QString("angle_style %1").arg(style.captured(2));
        } else if (style.captured(1) == "dihedral") {
            page = dihedral_map.value(style.captured(2), QString());
            help = QString("dihedral_style %1").arg(style.captured(2));
        } else if (style.captured(1) == "improper") {
            page = improper_map.value(style.captured(2), QString());
            help = QString("improper_style %1").arg(style.captured(2));
        }
    }

    style = QRegularExpression("^(fix|compute)\\s+\\w+\\s+\\w+\\s+(\\S+)").match(text);
    if (style.hasMatch()) {
        help = QString("%1 %2").arg(style.captured(1), style.captured(2));
        if (style.captured(1) == "fix") {
            page = fix_map.value(style.captured(2), QString());
        } else if (style.captured(1) == "compute") {
            page = compute_map.value(style.captured(2), QString());
        }
    }

    // could not find a matching "style", now try the plain command
    if (page.isEmpty() && !text.isEmpty()) {
        auto cmd = text.split(' ').at(0);
        help     = cmd;
        page     = cmd_map.value(cmd, QString());
    }
}

void CodeEditor::open_help()
{
    QAction *act = qobject_cast<QAction *>(sender());
    QDesktopServices::openUrl(
        QUrl(QString("https://docs.lammps.org/%1").arg(act->data().toString())));
}

// Local Variables:
// c-basic-offset: 4
// End:
