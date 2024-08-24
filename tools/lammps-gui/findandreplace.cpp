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

#include "findandreplace.h"

#include "codeeditor.h"
#include "lammpsgui.h"

#include <QApplication>
#include <QCheckBox>
#include <QDialogButtonBox>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QIcon>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QShortcut>
#include <QSizePolicy>
#include <QTextCursor>

/* ---------------------------------------------------------------------- */

FindAndReplace::FindAndReplace(CodeEditor *_editor, QWidget *parent) :
    QDialog(parent), editor(_editor), search(nullptr), replace(nullptr), withcase(nullptr),
    wrap(nullptr), whole(nullptr)
{
    auto *layout  = new QGridLayout;
    search        = new QLineEdit;
    replace       = new QLineEdit;
    withcase      = new QCheckBox("Match case");
    wrap          = new QCheckBox("Wrap around");
    whole         = new QCheckBox("Whole word");
    auto *next    = new QPushButton("&Next");
    auto *replone = new QPushButton("&Replace");
    auto *replall = new QPushButton("Replace &All");
    auto *done    = new QPushButton("&Done");

    layout->addWidget(new QLabel("Find:"), 0, 0, Qt::AlignRight);
    layout->addWidget(search, 0, 1, 1, 2, Qt::AlignLeft);
    layout->addWidget(new QLabel("Replace with:"), 1, 0, Qt::AlignRight);
    layout->addWidget(replace, 1, 1, 1, 2, Qt::AlignLeft);
    layout->addWidget(withcase, 2, 0, Qt::AlignLeft);
    layout->addWidget(wrap, 2, 1, Qt::AlignLeft);
    layout->addWidget(whole, 2, 2, Qt::AlignLeft);
    wrap->setChecked(true);

    auto *buttons = new QHBoxLayout;
    buttons->addWidget(next);
    buttons->addWidget(replone);
    buttons->addWidget(replall);
    buttons->addWidget(done);
    layout->addLayout(buttons, 3, 0, 1, 3, Qt::AlignHCenter);

    connect(next, &QPushButton::released, this, &FindAndReplace::find_next);
    connect(replone, &QPushButton::released, this, &FindAndReplace::replace_next);
    connect(replall, &QPushButton::released, this, &FindAndReplace::replace_all);
    connect(done, &QPushButton::released, this, &QDialog::accept);

    auto action = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_Q), this);
    connect(action, &QShortcut::activated, this, &FindAndReplace::quit);

    setLayout(layout);
    setWindowIcon(QIcon(":/icons/lammps-icon-128x128.png"));
    setWindowTitle("LAMMPS-GUI - Find and Replace");
}

/* ---------------------------------------------------------------------- */

void FindAndReplace::find_next()
{
    auto text = search->text();

    int find_flags = 0;
    if (withcase->isChecked()) find_flags |= QTextDocument::FindCaseSensitively;
    if (whole->isChecked()) find_flags |= QTextDocument::FindWholeWords;

    if (!text.isEmpty()) {
        if (!editor->find(text, (QTextDocument::FindFlag)find_flags) && wrap->isChecked()) {
            // nothing found from the current position to the end, reposition cursor and beginning
            editor->moveCursor(QTextCursor::Start, QTextCursor::MoveAnchor);
            editor->find(text, (QTextDocument::FindFlag)find_flags);
        }
    }
}

/* ---------------------------------------------------------------------- */

void FindAndReplace::replace_next()
{
    auto text = search->text();
    if (text.isEmpty()) return;

    auto cursor = editor->textCursor();
    auto flag   = withcase->isChecked() ? Qt::CaseSensitive : Qt::CaseInsensitive;

    // if selected text at cursor location matches search text, replace
    if (QString::compare(cursor.selectedText(), search->text(), flag) == 0)
        cursor.insertText(replace->text());

    find_next();
}

/* ---------------------------------------------------------------------- */

void FindAndReplace::replace_all()
{
    auto text = search->text();
    if (text.isEmpty()) return;

    // drop selection if we have one
    auto cursor = editor->textCursor();
    if (cursor.hasSelection()) cursor.movePosition(QTextCursor::Left);

    find_next();
    cursor = editor->textCursor();

    // keep replacing until find_next() does not find anything anymore
    while (cursor.hasSelection()) {
        cursor.insertText(replace->text());
        find_next();
        cursor = editor->textCursor();
    }
}

/* ---------------------------------------------------------------------- */

void FindAndReplace::quit()
{
    LammpsGui *main = nullptr;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") main = dynamic_cast<LammpsGui *>(widget);
    if (main) main->quit();
}

// Local Variables:
// c-basic-offset: 4
// End:
