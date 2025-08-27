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
#include <QLineEdit>
#include <QList>
#include <QRectF>
#include <QString>
#include <QTime>
#include <QWidget>

class QAction;
class QCheckBox;
class QCloseEvent;
class QEvent;
class QLabel;
class QMenuBar;
class QMenu;
class QSpinBox;
class RangeSlider;

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
    void set_units(const QString &_units);
    void set_norm(bool norm);

private slots:
    void quit();
    void stop_run();
    void select_smooth(int selection);
    void update_smooth();
    void update_tlabel();
    void update_ylabel();
    void update_xrange(int low, int high);
    void update_yrange(int low, int high);

    void saveAs();
    void exportDat();
    void exportCsv();
    void exportYaml();

    void change_chart(int index);

protected:
    void closeEvent(QCloseEvent *event) override;
    bool eventFilter(QObject *watched, QEvent *event) override;

private:
    bool do_raw, do_smooth;
    QMenuBar *menu;
    QMenu *file;
    QComboBox *columns;
    QAction *saveAsAct, *exportCsvAct, *exportDatAct, *exportYamlAct;
    QAction *closeAct, *stopAct, *quitAct;
    QComboBox *smooth;
    QSpinBox *window, *order;
    QLineEdit *chartTitle, *chartYlabel;
    QLabel *units;
    QCheckBox *norm;
    RangeSlider *xrange, *yrange;

    QString filename;
    QList<QtCharts::ChartViewer *> charts;
};

/* -------------------------------------------------------------------- */

#include <QChartView>
#include <QLineSeries>
#include <QValueAxis>
class QChart;

namespace QtCharts {
class ChartViewer : public QChartView {
    Q_OBJECT

public:
    explicit ChartViewer(const QString &title, int index, QWidget *parent = nullptr);
    ~ChartViewer() override;

    ChartViewer()                               = delete;
    ChartViewer(const ChartViewer &)            = delete;
    ChartViewer(ChartViewer &&)                 = delete;
    ChartViewer &operator=(const ChartViewer &) = delete;
    ChartViewer &operator=(ChartViewer &&)      = delete;

    void add_data(int step, double data);
    QRectF get_minmax() const;
    QList<QAbstractAxis *> get_axes() const { return chart->axes(); }
    void reset_zoom();
    void smooth_param(bool _do_raw, bool _do_smooth, int _window, int _order);
    void update_smooth();

    int get_index() const { return index; };
    int get_count() const { return series->count(); }
    QString get_title() const { return series->name(); }
    double get_step(int index) const { return (index < 0) ? 0.0 : series->at(index).x(); }
    double get_data(int index) const { return (index < 0) ? 0.0 : series->at(index).y(); }
    void set_tlabel(const QString &tlabel);
    void set_ylabel(const QString &ylabel);
    QString get_tlabel() const { return chart->title(); }
    QString get_xlabel() const { return xaxis->titleText(); }
    QString get_ylabel() const { return yaxis->titleText(); }

private:
    int last_step, index;
    int window, order;
    QChart *chart;
    QLineSeries *series, *smooth;
    QValueAxis *xaxis;
    QValueAxis *yaxis;
    QTime last_update;
    bool do_raw, do_smooth;
};
} // namespace QtCharts
#endif

// Local Variables:
// c-basic-offset: 4
// End:
