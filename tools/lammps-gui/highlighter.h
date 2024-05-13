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

#ifndef HIGHLIGHTER_H
#define HIGHLIGHTER_H

#include <QRegularExpression>
#include <QSyntaxHighlighter>
#include <QTextCharFormat>

class Highlighter : public QSyntaxHighlighter {
    Q_OBJECT

public:
    Highlighter(QTextDocument *parent = 0);

protected:
    void highlightBlock(const QString &text) override;

private:
    QRegularExpression isLattice1, isLattice2, isLattice3;
    QRegularExpression isOutput1, isOutput2, isRead;
    QTextCharFormat formatOutput, formatRead, formatLattice, formatSetup;
    QRegularExpression isStyle, isForce, isDefine, isUndo;
    QRegularExpression isParticle, isRun, isSetup, isSetup1;
    QTextCharFormat formatParticle, formatRun, formatDefine;
    QRegularExpression isVariable, isReference;
    QTextCharFormat formatVariable;
    QRegularExpression isNumber1, isNumber2, isNumber3, isNumber4;
    QTextCharFormat formatNumber;
    QRegularExpression isSpecial, isContinue;
    QTextCharFormat formatSpecial;
    QRegularExpression isComment;
    QRegularExpression isQuotedComment;
    QTextCharFormat formatComment;
    QRegularExpression isTriple;
    QRegularExpression isString;
    QTextCharFormat formatString;

    int in_triple;
    int startIndex;
};
#endif
// Local Variables:
// c-basic-offset: 4
// End:
