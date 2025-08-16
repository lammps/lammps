/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FLAGWARNINGS_H
#define FLAGWARNINGS_H

#include <QRegularExpression>
#include <QSyntaxHighlighter>
#include <QTextCharFormat>

class QLabel;
class QTextDocument;

class FlagWarnings : public QSyntaxHighlighter {
    Q_OBJECT

public:
    explicit FlagWarnings(QLabel *label = nullptr, QTextDocument *parent = nullptr);
    ~FlagWarnings() override = default;

    FlagWarnings()                                = delete;
    FlagWarnings(const FlagWarnings &)            = delete;
    FlagWarnings(FlagWarnings &&)                 = delete;
    FlagWarnings &operator=(const FlagWarnings &) = delete;
    FlagWarnings &operator=(FlagWarnings &&)      = delete;

    int get_nwarnings() const { return nwarnings; }

protected:
    void highlightBlock(const QString &text) override;

private:
    QRegularExpression isWarning;
    QRegularExpression isURL;
    QTextCharFormat formatWarning;
    QTextCharFormat formatURL;
    QLabel *summary;
    QTextDocument *document;
    int nwarnings, oldwarnings, nlines, oldlines;
};
#endif
// Local Variables:
// c-basic-offset: 4
// End:
