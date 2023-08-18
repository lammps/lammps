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

#include <QAction>
#include <QDesktopServices>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QIcon>
#include <QKeySequence>
#include <QMenu>
#include <QMimeData>
#include <QPainter>
#include <QRegularExpression>
#include <QTextBlock>
#include <QUrl>

CodeEditor::CodeEditor(QWidget *parent) : QPlainTextEdit(parent)
{
    help_action = new QShortcut(QKeySequence::fromString("Ctrl+?"), parent);
    connect(help_action, &QShortcut::activated, this, &CodeEditor::get_help);

    // initialize help system
    QFile help_index(":/help_index.table");
    if (help_index.open(QIODevice::ReadOnly | QIODevice::Text)) {
        while (!help_index.atEnd()) {
            auto line  = QString(help_index.readLine());
            auto words = line.trimmed().split(' ');
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
}

int CodeEditor::lineNumberAreaWidth()
{
    int digits = 1;
    int max    = qMax(1, blockCount());
    while (max >= 10) {
        max /= 10;
        ++digits;
    }

    int space = 3 + fontMetrics().horizontalAdvance(QLatin1Char('9')) * digits;

    return space;
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
        auto file = event->mimeData()->urls()[0].url().remove("file://");
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

void CodeEditor::highlightCurrentLine()
{
    QList<QTextEdit::ExtraSelection> extraSelections;

    if (!isReadOnly()) {
        QTextEdit::ExtraSelection selection;

        QColor lineColor = QColor(Qt::yellow).lighter(160);

        selection.format.setBackground(lineColor);
        selection.format.setProperty(QTextFormat::FullWidthSelection, true);
        selection.cursor = textCursor();
        selection.cursor.clearSelection();
        extraSelections.append(selection);
    }

    setExtraSelections(extraSelections);
}

void CodeEditor::lineNumberAreaPaintEvent(QPaintEvent *event)
{
    QPainter painter(lineNumberArea);
    painter.fillRect(event->rect(), Qt::lightGray);
    QTextBlock block = firstVisibleBlock();
    int blockNumber  = block.blockNumber();
    int top          = qRound(blockBoundingGeometry(block).translated(contentOffset()).top());
    int bottom       = top + qRound(blockBoundingRect(block).height());
    while (block.isValid() && top <= event->rect().bottom()) {
        if (block.isVisible() && bottom >= event->rect().top()) {
            QString number = QString::number(blockNumber + 1);
            painter.setPen(Qt::black);
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
    // reposition the cursor here?
    QString page, help;
    find_help(page, help);

    // print augmented context menu if an entry was found
    auto *menu = createStandardContextMenu();
    if (!page.isEmpty()) {
        menu->addSeparator();
        auto action = menu->addAction(QString("View Documentation for '%1'").arg(help));
        action->setIcon(QIcon(":/system-help.png"));
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
            action2->setIcon(QIcon(":/system-help.png"));
            action2->setData(page);
            connect(action2, &QAction::triggered, this, &CodeEditor::open_help);
        }
    }
    auto action3 = menu->addAction(QString("LAMMPS Manual"));
    action3->setIcon(QIcon(":/help-browser.png"));
    action3->setData(QString());
    connect(action3, &QAction::triggered, this, &CodeEditor::open_help);

    menu->exec(event->globalPos());
    delete menu;
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
