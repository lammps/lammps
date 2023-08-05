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

#ifndef CHARTVIEWER_H
#define CHARTVIEWER_H

#include <QList>
#include <QString>
#include <QWidget>
#include <QtCharts>

class QAction;
class QMenuBar;
class QMenu;
class QComboBox;
class ChartViewer;

class ChartWindow : public QWidget {
    Q_OBJECT

public:
    ChartWindow(const QString &filename, QWidget *parent = nullptr);

    bool has_charts() const { return !charts.isEmpty(); }
    void add_chart(const QString &title, int index);
    void add_data(int step, double data, int index);

private slots:
    void saveAs();
    void change_chart(int index);

    //    void normalSize();

protected:
    void closeEvent(QCloseEvent *event) override;

private:
    QMenuBar *menu;
    QMenu *file;
    QComboBox *columns;
    QAction *saveAsAct;
    QAction *closeAct;

    QString filename;
    int active_chart;
    QList<ChartViewer *> charts;
};

/* -------------------------------------------------------------------- */

class ChartViewer : public QtCharts::QChartView {
    Q_OBJECT

public:
    explicit ChartViewer(const QString &title, int index, QWidget *parent = nullptr);

    void add_data(int step, double data);
    int get_index() const { return index; };

private:
    int last_step, index;
    QtCharts::QChart *chart;
    QtCharts::QLineSeries *series;
    QtCharts::QValueAxis *xaxis;
    QtCharts::QValueAxis *yaxis;
};
#endif

// Local Variables:
// c-basic-offset: 4
// End:
