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

#ifndef CODEEDITOR_H
#define CODEEDITOR_H

#include <QMap>
#include <QPlainTextEdit>
#include <QShortcut>
#include <QString>

class CodeEditor : public QPlainTextEdit {
    Q_OBJECT

public:
    CodeEditor(QWidget *parent = nullptr);

    void lineNumberAreaPaintEvent(QPaintEvent *event);
    int lineNumberAreaWidth();

protected:
    void resizeEvent(QResizeEvent *event) override;
    void dragEnterEvent(QDragEnterEvent *event) override;
    bool canInsertFromMimeData(const QMimeData *source) const override;
    void dropEvent(QDropEvent *event) override;
    void contextMenuEvent(QContextMenuEvent *event) override;

private slots:
    void updateLineNumberAreaWidth(int newBlockCount);
    void highlightCurrentLine();
    void updateLineNumberArea(const QRect &rect, int dy);
    void get_help();
    void find_help(QString &page, QString &help);
    void open_help();

private:
    QWidget *lineNumberArea;
    QShortcut *help_action;

    QMap<QString, QString> cmd_map;
    QMap<QString, QString> fix_map;
    QMap<QString, QString> compute_map;
    QMap<QString, QString> pair_map;
    QMap<QString, QString> bond_map;
    QMap<QString, QString> angle_map;
    QMap<QString, QString> dihedral_map;
    QMap<QString, QString> improper_map;
    QMap<QString, QString> dump_map;
};

#endif
// Local Variables:
// c-basic-offset: 4
// End:
