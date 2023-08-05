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

#include "chartviewer.h"

#include <QLineSeries>
#include <QSettings>

using namespace QtCharts;

ChartViewer::ChartViewer(QWidget *parent) :
    QChartView(parent), last_step(-1), chart(new QChart), series(new QLineSeries),
    xaxis(new QValueAxis), yaxis(new QValueAxis)
{
    chart->legend()->hide();
    chart->addAxis(xaxis,Qt::AlignBottom);
    chart->addAxis(yaxis,Qt::AlignLeft);
    chart->addSeries(series);
    series->attachAxis(xaxis);
    series->attachAxis(yaxis);
    xaxis->setTitleText("Time step");
    xaxis->setTickCount(5);
    yaxis->setTickCount(5);
    xaxis->setMinorTickCount(5);
    yaxis->setMinorTickCount(5);

    setRenderHint(QPainter::Antialiasing);
    setChart(chart);
    setRubberBand(QChartView::RectangleRubberBand);

    QSettings settings;
    resize(settings.value("chartx", 500).toInt(), settings.value("charty", 320).toInt());
}

void ChartViewer::add_column(const QString &title)
{
    yaxis->setTitleText(title);
    series->setName(title);
}

void ChartViewer::add_data(int step, int column, double data)
{
    if (last_step < step) {
        last_step = step;
        series->append(step, data);
        auto points = series->pointsVector();

        qreal xmin = 1.0e100;
        qreal xmax = -1.0e100;
        qreal ymin = 1.0e100;
        qreal ymax = -1.0e100;
        for (auto &p : points) {
            xmin = qMin(xmin, p.x());
            xmax = qMax(xmax, p.x());
            ymin = qMin(ymin, p.y());
            ymax = qMax(ymax, p.y());
        }
        xaxis->setRange(xmin, xmax);
        yaxis->setRange(ymin, ymax);
    }
}

void ChartViewer::closeEvent(QCloseEvent *event)
{
    QSettings settings;
    if (!isMaximized()) {
        settings.setValue("chartx", width());
        settings.setValue("charty", height());
    }
    QChartView::closeEvent(event);
}

// Local Variables:
// c-basic-offset: 4
// End:
