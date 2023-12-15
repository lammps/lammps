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

#include "highlighter.h"

Highlighter::Highlighter(QTextDocument *parent) :
    QSyntaxHighlighter(parent),
    isLattice1(QStringLiteral("^\\s*(units|atom_style|change_box|dielectric|dimension)\\s+(\\S+)")),
    isLattice2(QStringLiteral("^\\s*(lattice|region|create_box|create_atoms|delete_atoms|displace_"
                              "atoms)\\s+(\\S+)\\s+(\\S+)")),
    isLattice3(QStringLiteral("^\\s*(boundary|replicate)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)")),
    isOutput1(QStringLiteral("^\\s*(echo|log|write_data|write_coeff|write_restart|restart|info|"
                             "thermo|print|thermo_style|"
                             "timer|pair_write|bond_write|angle_write|dihedral_write)\\s+(\\S+)")),
    isOutput2(QStringLiteral("^\\s*(write_dump|shell|thermo_modify)\\s+(\\S+)\\s+(\\S+)")),
    isRead(QStringLiteral("^\\s*(include|read_restart|read_data|read_dump|molecule)")),
    isStyle(QStringLiteral("^\\s*(fix|compute|dump)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)")),
    isForce(QStringLiteral(
        "^\\s*(pair_style|bond_style|angle_style|dihedral_style|improper_style|kspace_style|pair_"
        "coeff|angle_coeff|bond_coeff|dihedral_coeff|improper_coeff)\\s+(\\S+)")),
    isDefine(QStringLiteral("^\\s*(group|variable|python|set|group2ndx|ndx2group|kim|kim_query|mdi)"
                            "\\s+(\\S+)\\s+(\\S+)")),
    isUndo(QStringLiteral("^\\s*(unfix|uncompute|undump|label|jump|next)\\s+(\\S+)")),
    isParticle(QStringLiteral("^\\s*(pair_modify|mass|velocity|create_bonds|delete_"
                              "bonds|kspace_modify|labelmap|atom_modify)\\s+(\\S+)")),
    isRun(QStringLiteral("^\\s*(minimize|minimize/kk|run|rerun|tad|neb|neb/spin|prd|server|temper/"
                         "npt|temper/grem|temper|message|hyper|dynamical_matrix|dynamical_matrix/"
                         "kk|third_order|third_order/kk|fitpod)")),
    isSetup(QStringLiteral("^\\s*(min_modify|neighbor|neigh_modify|special_bonds|balance|box|clear|"
                           "quit|newton|partition|processors|reset_atoms|reset_ids)")),
    isSetup1(
        QStringLiteral("^\\s*(min_style|run_style|timestep|suffix|plugin|comm_modify|comm_style|"
                       "package|reset_timestep|dump_modify|fix_modify|compute_modify)\\s+(\\S+)")),
    isVariable(QStringLiteral("(\\$[a-z]|\\${[^} ]+}|\\$\\(\\S+\\))")),
    isReference(
        QStringLiteral("\\s+(c_\\S+|C_\\S+|f_\\S+|F_\\S+|i_\\S+|i2_\\S+|d_\\S+|d2_\\S+|v_\\S+)")),
    isNumber1(QStringLiteral("(^|\\s+)[-+]?[0-9:*]+")), // integer and integer ranges
    isNumber2(QStringLiteral("(^|\\s+)[-+]?[0-9]+\\.[0-9]*[edED]?[-+]?[0-9]*")), // floating point 1
    isNumber3(QStringLiteral("(^|\\s+)[-+]?[0-9]*\\.[0-9]+[edED]?[-+]?[0-9]*")), // floating point 2
    isNumber4(QStringLiteral("(^|\\s+)[-+]?[0-9]+([edED][-+]?[0-9]+)?")),        // floating point 3
    isSpecial(QStringLiteral("(\\sINF|\\sEDGE|\\sNULL|\\sSELF|if\\s|then\\s|else\\s|elif\\s)")),
    isContinue(QStringLiteral("&$")), isComment(QStringLiteral("#.*")),
    isQuotedComment(QStringLiteral("(\".*#.*\"|'.*#.*')")),
    isTriple(QStringLiteral("[^\"]*\"\"\"[^\"]*")),
    isString(QStringLiteral("(\".+?\"|'.+?'|\"\"\".*\"\"\")")), in_triple(false)
{
    formatNumber.setForeground(Qt::blue);
    formatString.setForeground(Qt::darkGreen);
    formatString.setFontWeight(QFont::Normal);
    formatComment.setForeground(Qt::red);
    formatSpecial.setForeground(Qt::darkMagenta);
    formatSpecial.setFontWeight(QFont::Bold);
    formatParticle.setForeground(Qt::darkRed);
    formatParticle.setFontWeight(QFont::Bold);
    formatRun.setForeground(Qt::darkBlue);
    formatRun.setFontWeight(QFont::Bold);
    formatVariable.setForeground(Qt::darkGray);
    formatVariable.setFontWeight(QFont::Bold);

    formatOutput.setForeground(Qt::darkYellow);
    formatOutput.setFontWeight(QFont::Bold);
    formatRead.setForeground(Qt::magenta);
    formatRead.setFontWeight(QFont::Bold);
    formatLattice.setForeground(Qt::darkGreen);
    formatLattice.setFontWeight(QFont::Bold);
    formatSetup.setForeground(Qt::darkCyan);
    formatSetup.setFontWeight(QFont::Bold);
}

void Highlighter::highlightBlock(const QString &text)
{
    // nothing to do for empty lines
    if (text.isEmpty()) return;

    auto match = isLattice1.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatLattice);
        setFormat(match.capturedStart(2), match.capturedLength(2), formatRun);
    }

    match = isLattice2.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatLattice);
        setFormat(match.capturedStart(2), match.capturedLength(2), formatString);
        setFormat(match.capturedStart(3), match.capturedLength(3), formatRun);
    }

    match = isLattice3.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatLattice);
        setFormat(match.capturedStart(2), match.capturedLength(2), formatString);
        setFormat(match.capturedStart(3), match.capturedLength(3), formatString);
        setFormat(match.capturedStart(4), match.capturedLength(4), formatString);
    }

    match = isOutput1.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatOutput);
        setFormat(match.capturedStart(2), match.capturedLength(2), formatString);
    }

    match = isOutput2.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatOutput);
        setFormat(match.capturedStart(2), match.capturedLength(2), formatString);
        setFormat(match.capturedStart(3), match.capturedLength(3), formatRun);
    }

    match = isRead.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatRead);
        setFormat(match.capturedStart(2), match.capturedLength(2), formatString);
    }

    match = isStyle.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatParticle);
        setFormat(match.capturedStart(2), match.capturedLength(2), formatNumber);
        setFormat(match.capturedStart(3), match.capturedLength(3), formatString);
        setFormat(match.capturedStart(4), match.capturedLength(4), formatRun);
    }

    match = isForce.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatParticle);
        setFormat(match.capturedStart(2), match.capturedLength(2), formatRun);
    }

    match = isUndo.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatSpecial);
        setFormat(match.capturedStart(2), match.capturedLength(2), formatString);
    }

    match = isDefine.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatParticle);
        setFormat(match.capturedStart(2), match.capturedLength(2), formatString);
        setFormat(match.capturedStart(3), match.capturedLength(3), formatRun);
    }

    match = isParticle.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatParticle);
        setFormat(match.capturedStart(2), match.capturedLength(2), formatString);
    }

    match = isRun.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatRun);
    }

    match = isSetup.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatSetup);
    }

    match = isSetup1.match(text);
    if (match.hasMatch()) {
        setFormat(match.capturedStart(1), match.capturedLength(1), formatSetup);
        setFormat(match.capturedStart(2), match.capturedLength(2), formatString);
    }

    // numbers
    QRegularExpression numbers[] = {isNumber1, isNumber2, isNumber3, isNumber4};
    for (auto &number : numbers) {
        auto num = number.globalMatch(text);
        while (num.hasNext()) {
            auto match = num.next();
            setFormat(match.capturedStart(), match.capturedLength(), formatNumber);
        }
    }

    // variables
    auto vars = isVariable.globalMatch(text);
    while (vars.hasNext()) {
        auto match = vars.next();
        setFormat(match.capturedStart(), match.capturedLength(), formatVariable);
    }

    // references
    auto refs = isReference.globalMatch(text);
    while (refs.hasNext()) {
        auto match = refs.next();
        setFormat(match.capturedStart(), match.capturedLength(), formatVariable);
    }

    // continuation character
    auto multiline = isContinue.match(text);
    if (multiline.hasMatch())
        setFormat(multiline.capturedStart(0), multiline.capturedLength(0), formatSpecial);

    // special keywords
    auto special = isSpecial.globalMatch(text);
    while (special.hasNext()) {
        auto match = special.next();
        setFormat(match.capturedStart(), match.capturedLength(), formatSpecial);
    }

    // comments, must come before strings but after other keywords.
    auto comment = isComment.match(text);
    if (comment.hasMatch() && !isQuotedComment.match(text).hasMatch() && !in_triple) {
        setFormat(comment.capturedStart(0), comment.capturedLength(0), formatComment);
        return;
    }

    // strings, must come last so they can overwrite other formatting
    auto string = isString.globalMatch(text);
    while (string.hasNext()) {
        auto match = string.next();
        setFormat(match.capturedStart(), match.capturedLength(), formatString);
    }

    auto triple = isTriple.match(text);
    if (triple.hasMatch()) {
        if (in_triple) {
            in_triple = false;
            setFormat(0, triple.capturedStart(0) + triple.capturedLength(0), formatString);
        } else {
            in_triple = true;
            setFormat(triple.capturedStart(0), -1, formatString);
        }
    } else {
        if (in_triple) setFormat(0, text.size(), formatString);
    }
}
// Local Variables:
// c-basic-offset: 4
// End:
