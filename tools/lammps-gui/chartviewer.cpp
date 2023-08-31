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

#include <QHBoxLayout>
#include <QLineSeries>
#include <QSettings>
#include <QSpacerItem>
#include <QVBoxLayout>

using namespace QtCharts;

ChartWindow::ChartWindow(const QString &_filename, QWidget *parent) :
    QWidget(parent), menu(new QMenuBar), file(new QMenu("&File")), active_chart(-1),
    filename(_filename)
{
    auto *top = new QHBoxLayout;
    menu->addMenu(file);
    menu->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Preferred);

    columns = new QComboBox;
    top->addWidget(menu);
    top->addSpacerItem(new QSpacerItem(1, 1, QSizePolicy::Expanding, QSizePolicy::Minimum));
    top->addWidget(new QLabel("Select data:"));
    top->addWidget(columns);
    saveAsAct = file->addAction("&Save Graph As...", this, &ChartWindow::saveAs);
    saveAsAct->setIcon(QIcon(":/document-save-as.png"));
    exportCsvAct = file->addAction("&Export data to CSV...", this, &ChartWindow::exportCsv);
    exportCsvAct->setIcon(QIcon(":/application-calc.png"));
    exportDatAct = file->addAction("Export data to &Gnuplot...", this, &ChartWindow::exportDat);
    exportDatAct->setIcon(QIcon(":/application-plot.png"));
    file->addSeparator();
    closeAct = file->addAction("&Close", this, &QWidget::close);
    closeAct->setIcon(QIcon(":/window-close.png"));
    auto *layout = new QVBoxLayout;
    layout->addLayout(top);
    setLayout(layout);

    connect(columns, SIGNAL(currentIndexChanged(int)), this, SLOT(change_chart(int)));
    QSettings settings;
    resize(settings.value("chartx", 500).toInt(), settings.value("charty", 320).toInt());
}

int ChartWindow::get_step() const
{
    if (charts.size() > 0) {
        auto *v = charts[0];
        return (int)v->get_step(v->get_count() - 1);
    } else {
        return -1;
    }
}

void ChartWindow::reset_charts()
{
    while (layout()->count() > 1) {
        auto *item = layout()->takeAt(1);
        if (item) {
            layout()->removeItem(item);
            delete item->widget();
            delete item;
        }
    }
    charts.clear();
    columns->clear();
    active_chart = 0;
}

void ChartWindow::add_chart(const QString &title, int index)
{
    auto *chart = new ChartViewer(title, index);
    layout()->addWidget(chart);
    columns->addItem(title, index);
    columns->show();
    // hide all but the first chart added
    if (charts.size() > 0) chart->hide();
    charts.append(chart);
    active_chart = 0;
}

void ChartWindow::add_data(int step, double data, int index)
{
    for (auto &c : charts)
        if (c->get_index() == index) c->add_data(step, data);
}

void ChartWindow::saveAs()
{
    if (charts.empty() || (active_chart < 0)) return;
    QString defaultname = filename + "." + columns->currentText() + ".png";
    if (filename.isEmpty()) defaultname = columns->currentText() + ".png";
    QString fileName = QFileDialog::getSaveFileName(this, "Save Chart as Image", defaultname,
                                                    "Image Files (*.jpg *.png *.bmp *.ppm)");
    if (!fileName.isEmpty()) {
        charts[active_chart]->grab().save(fileName);
    }
}

void ChartWindow::exportDat()
{
    if (charts.empty() || (active_chart < 0)) return;
    QString defaultname = filename + ".dat";
    if (filename.isEmpty()) defaultname = "lammpsdata.dat";
    QString fileName = QFileDialog::getSaveFileName(this, "Save Chart as Gnuplot data", defaultname,
                                                    "Image Files (*.dat)");
    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {

            file.write("# Thermodynamic data from ");
            file.write(filename.toLocal8Bit());
            file.write("\n# Columns:");
            for (auto &c : charts) {
                file.write(" ");
                file.write(c->get_title());
            }
            file.write("\n");

            int lines = charts[0]->get_count();
            for (int i = 0; i < lines; ++i) {
                // timestep
                file.write(QString::number(charts[0]->get_step(i)).toLocal8Bit());
                for (auto &c : charts) {
                    file.write(" ");
                    file.write(QString::number(c->get_data(i)).toLocal8Bit());
                }
                file.write("\n");
            }
            file.close();
        }
    }
}

void ChartWindow::exportCsv()
{
    if (charts.empty() || (active_chart < 0)) return;
    QString defaultname = filename + ".csv";
    if (filename.isEmpty()) defaultname = "lammpsdata.csv";
    QString fileName = QFileDialog::getSaveFileName(this, "Save Chart as CSV data", defaultname,
                                                    "Image Files (*.csv)");
    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {

            file.write("Step");
            for (auto &c : charts) {
                file.write(",");
                file.write(c->get_title());
            }
            file.write("\n");

            int lines = charts[0]->get_count();
            for (int i = 0; i < lines; ++i) {
                // timestep
                file.write(QString::number(charts[0]->get_step(i)).toLocal8Bit());
                for (auto &c : charts) {
                    file.write(",");
                    file.write(QString::number(c->get_data(i)).toLocal8Bit());
                }
                file.write("\n");
            }
            file.close();
        }
    }
}

void ChartWindow::change_chart(int index)
{
    int choice = columns->currentData().toInt();
    for (auto &c : charts) {
        if (choice == c->get_index())
            c->show();
        else
            c->hide();
    }
}

void ChartWindow::closeEvent(QCloseEvent *event)
{
    QSettings settings;
    if (!isMaximized()) {
        settings.setValue("chartx", width());
        settings.setValue("charty", height());
    }
    QWidget::closeEvent(event);
}

/* -------------------------------------------------------------------- */

ChartViewer::ChartViewer(const QString &title, int _index, QWidget *parent) :
    QChartView(parent), last_step(-1), index(_index), chart(new QChart), series(new QLineSeries),
    xaxis(new QValueAxis), yaxis(new QValueAxis)
{
    chart->legend()->hide();
    chart->addAxis(xaxis, Qt::AlignBottom);
    chart->addAxis(yaxis, Qt::AlignLeft);
    chart->addSeries(series);
    series->attachAxis(xaxis);
    series->attachAxis(yaxis);
    xaxis->setTitleText("Time step");
    xaxis->setTickCount(5);
    xaxis->setLabelFormat("%d");
    yaxis->setTickCount(5);
    xaxis->setMinorTickCount(5);
    yaxis->setMinorTickCount(5);
    yaxis->setTitleText(title);
    series->setName(title);

    setRenderHint(QPainter::Antialiasing);
    setChart(chart);
    setRubberBand(QChartView::RectangleRubberBand);
}

/* -------------------------------------------------------------------- */

void ChartViewer::add_data(int step, double data)
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

// Local Variables:
// c-basic-offset: 4
// End:
