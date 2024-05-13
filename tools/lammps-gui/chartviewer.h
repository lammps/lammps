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

#include <QComboBox>
#include <QList>
#include <QString>
#include <QWidget>

class QAction;
class QMenuBar;
class QMenu;
namespace QtCharts {
class ChartViewer;
}

class ChartWindow : public QWidget {
    Q_OBJECT

public:
    ChartWindow(const QString &filename, QWidget *parent = nullptr);

    int num_charts() const { return charts.size(); }
    bool has_title(const QString &title, int index) const
    {
        return (columns->itemText(index) == title);
    }
    int get_step() const;
    void reset_charts();
    void add_chart(const QString &title, int index);
    void add_data(int step, double data, int index);

private slots:
    void quit();
    void reset_zoom();
    void stop_run();

    void saveAs();
    void exportDat();
    void exportCsv();

    void change_chart(int index);

protected:
    void closeEvent(QCloseEvent *event) override;
    bool eventFilter(QObject *watched, QEvent *event) override;

private:
    QMenuBar *menu;
    QMenu *file;
    QComboBox *columns;
    QAction *saveAsAct, *exportCsvAct, *exportDatAct;
    QAction *closeAct, *stopAct, *quitAct;

    QString filename;
    QList<QtCharts::ChartViewer *> charts;
};

/* -------------------------------------------------------------------- */

#include <QChart>
#include <QChartView>
#include <QLineSeries>
#include <QValueAxis>

namespace QtCharts {
class ChartViewer : public QChartView {
    Q_OBJECT

public:
    explicit ChartViewer(const QString &title, int index, QWidget *parent = nullptr);

    void add_data(int step, double data);
    void reset_zoom();

    int get_index() const { return index; };
    int get_count() const { return series->count(); }
    const char *get_title() const { return series->name().toLocal8Bit(); }
    double get_step(int index) const { return (index < 0) ? 0.0 : series->at(index).x(); }
    double get_data(int index) const { return (index < 0) ? 0.0 : series->at(index).y(); }

private:
    int last_step, index;
    QChart *chart;
    QLineSeries *series;
    QValueAxis *xaxis;
    QValueAxis *yaxis;
};
} // namespace QtCharts
#endif

// Local Variables:
// c-basic-offset: 4
// End:
