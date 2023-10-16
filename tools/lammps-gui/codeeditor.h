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

#include <QFont>
#include <QMap>
#include <QPlainTextEdit>
#include <QString>
#include <QStringList>

class QCompleter;
class QStringListModel;
class QShortcut;

class CodeEditor : public QPlainTextEdit {
    Q_OBJECT

public:
    CodeEditor(QWidget *parent = nullptr);
    ~CodeEditor() override;

    void lineNumberAreaPaintEvent(QPaintEvent *event);
    int lineNumberAreaWidth();
    void setFont(const QFont &newfont);
    void setCursor(int block);
    void setHighlight(int block, bool error);
    void setReformatOnReturn(bool flag) { reformat_on_return = flag; }
    void setAutoComplete(bool flag) { automatic_completion = flag; }
    QString reformatLine(const QString &line);

    void setCommandList(const QStringList &words);
    void setFixList(const QStringList &words);
    void setComputeList(const QStringList &words);
    void setDumpList(const QStringList &words);
    void setAtomList(const QStringList &words);
    void setPairList(const QStringList &words);
    void setBondList(const QStringList &words);
    void setAngleList(const QStringList &words);
    void setDihedralList(const QStringList &words);
    void setImproperList(const QStringList &words);
    void setKspaceList(const QStringList &words);
    void setRegionList(const QStringList &words);
    void setIntegrateList(const QStringList &words);
    void setMinimizeList(const QStringList &words);
    void setVariableList(const QStringList &words);
    void setUnitsList(const QStringList &words);
    void setGroupList();
    void setVarNameList();
    void setComputeIDList();
    void setFixIDList();
    void setFileList();

    static constexpr int NO_HIGHLIGHT = 1 << 30;

protected:
    void resizeEvent(QResizeEvent *event) override;
    void dragEnterEvent(QDragEnterEvent *event) override;
    bool canInsertFromMimeData(const QMimeData *source) const override;
    void dropEvent(QDropEvent *event) override;
    void contextMenuEvent(QContextMenuEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;

private slots:
    void updateLineNumberAreaWidth(int newBlockCount);
    void updateLineNumberArea(const QRect &rect, int dy);
    void get_help();
    void find_help(QString &page, QString &help);
    void open_help();
    void reformatCurrentLine();
    void runCompletion();
    void insertCompletedCommand(const QString &completion);

private:
    QWidget *lineNumberArea;
    QShortcut *help_action;
    QCompleter *current_comp, *command_comp, *fix_comp, *compute_comp, *dump_comp, *atom_comp,
        *pair_comp, *bond_comp, *angle_comp, *dihedral_comp, *improper_comp, *kspace_comp,
        *region_comp, *integrate_comp, *minimize_comp, *variable_comp, *units_comp, *group_comp,
        *varname_comp, *fixid_comp, *compid_comp, *file_comp;

    int highlight;
    bool reformat_on_return;
    bool automatic_completion;

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
