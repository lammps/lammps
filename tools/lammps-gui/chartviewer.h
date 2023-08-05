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

#include <QtCharts>

class ChartViewer : public QtCharts::QChartView {
    Q_OBJECT

public:
    ChartViewer(QWidget *parent = nullptr);
    bool has_columns() const { return last_step >= 0; }
    void add_column(const QString &title);
    void add_data(int step, int column, double data);
    int get_last_step() const;

protected:
    void closeEvent(QCloseEvent *event) override;

private:
    int last_step;
    QtCharts::QChart *chart;
    QtCharts::QLineSeries *series;
    QtCharts::QValueAxis *xaxis;
    QtCharts::QValueAxis *yaxis;
};
#endif

// Local Variables:
// c-basic-offset: 4
// End:
