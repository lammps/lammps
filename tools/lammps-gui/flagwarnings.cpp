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

#include "flagwarnings.h"

#include "helpers.h"

#include <QColor>
#include <QFont>
#include <QLabel>
#include <QTextDocument>

// workaround for Qt-5.12
#if QT_VERSION < QT_VERSION_CHECK(5, 14, 0)
namespace QColorConstants {
const QColor Red  = QColor::fromRgb(0xff, 0x00, 0x00);
const QColor Blue = QColor::fromRgb(0x00, 0x00, 0xff);
} // namespace QColorConstants
#endif

FlagWarnings::FlagWarnings(QLabel *label, QTextDocument *parent) :
    QSyntaxHighlighter(parent), isWarning(QStringLiteral("^(ERROR|WARNING).*$")),
    isURL(QStringLiteral("^.*(https://docs.lammps.org/err[0-9]+).*$")), summary(label),
    document(parent)
{
    nwarnings = nlines = 0;

    formatWarning.setForeground(QColorConstants::Red);
    formatWarning.setFontWeight(QFont::Bold);
    formatURL.setForeground(QColorConstants::Blue);
    formatURL.setFontWeight(QFont::Bold);
}

void FlagWarnings::highlightBlock(const QString &text)
{
    // nothing to do for empty lines
    if (text.isEmpty()) return;

    // highlight errors or warnings
    auto match = isWarning.match(text);
    if (match.hasMatch()) {
        ++nwarnings;
        setFormat(match.capturedStart(0), match.capturedLength(0), formatWarning);
    }

    // highlight ErrorURL links
    match = isURL.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatURL);
    }

    // update error summary label
    if (document && summary) {
        summary->setText(
            QString("%1 Warnings / Errors - %2 Lines").arg(nwarnings).arg(document->lineCount()));
        summary->repaint();
    }
}

// Local Variables:
// c-basic-offset: 4
// End:
