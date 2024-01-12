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

#include "lammpsgui.h"

#include <QAction>
#include <QApplication>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QKeySequence>
#include <QLabel>
#include <QLayout>
#include <QLineSeries>
#include <QMenu>
#include <QMenuBar>
#include <QPushButton>
#include <QSettings>
#include <QSpacerItem>
#include <QTextStream>
#include <QVBoxLayout>

#include <cmath>

using namespace QtCharts;

ChartWindow::ChartWindow(const QString &_filename, QWidget *parent) :
    QWidget(parent), menu(new QMenuBar), file(new QMenu("&File")), filename(_filename)
{
    auto *top = new QHBoxLayout;
    menu->addMenu(file);
    menu->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Preferred);

    // workaround for incorrect highlight bug on macOS
    auto *dummy = new QPushButton(QIcon(), "");
    dummy->hide();
    auto *normal = new QPushButton(QIcon(":/icons/gtk-zoom-fit.png"), "");
    normal->setToolTip("Reset zoom to normal");

    columns = new QComboBox;
    top->addWidget(menu);
    top->addSpacerItem(new QSpacerItem(1, 1, QSizePolicy::Expanding, QSizePolicy::Minimum));
    top->addWidget(dummy);
    top->addWidget(normal);
    top->addWidget(new QLabel("Select data:"));
    top->addWidget(columns);
    saveAsAct = file->addAction("&Save Graph As...", this, &ChartWindow::saveAs);
    saveAsAct->setIcon(QIcon(":/icons/document-save-as.png"));
    exportCsvAct = file->addAction("&Export data to CSV...", this, &ChartWindow::exportCsv);
    exportCsvAct->setIcon(QIcon(":/icons/application-calc.png"));
    exportDatAct = file->addAction("Export data to &Gnuplot...", this, &ChartWindow::exportDat);
    exportDatAct->setIcon(QIcon(":/icons/application-plot.png"));
    file->addSeparator();
    stopAct = file->addAction("Stop &Run", this, &ChartWindow::stop_run);
    stopAct->setIcon(QIcon(":/icons/process-stop.png"));
    stopAct->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_Slash));
    closeAct = file->addAction("&Close", this, &QWidget::close);
    closeAct->setIcon(QIcon(":/icons/window-close.png"));
    closeAct->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_W));
    quitAct = file->addAction("&Quit", this, &ChartWindow::quit);
    quitAct->setIcon(QIcon(":/icons/application-exit.png"));
    quitAct->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_Q));
    auto *layout = new QVBoxLayout;
    layout->addLayout(top);
    setLayout(layout);

    connect(normal, &QPushButton::released, this, &ChartWindow::reset_zoom);
    connect(columns, SIGNAL(currentIndexChanged(int)), this, SLOT(change_chart(int)));
    installEventFilter(this);

    QSettings settings;
    resize(settings.value("chartx", 500).toInt(), settings.value("charty", 320).toInt());
}

int ChartWindow::get_step() const
{
    if (charts.size() > 0) {
        auto *v = charts[0];
        if (v)
            return (int)v->get_step(v->get_count() - 1);
        else
            return -1;
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
}

void ChartWindow::add_data(int step, double data, int index)
{
    for (auto &c : charts)
        if (c->get_index() == index) c->add_data(step, data);
}

void ChartWindow::quit()
{
    LammpsGui *main = nullptr;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") main = dynamic_cast<LammpsGui *>(widget);
    if (main) main->quit();
}

void ChartWindow::reset_zoom()
{
    int choice = columns->currentData().toInt();
    charts[choice]->reset_zoom();
}

void ChartWindow::stop_run()
{
    LammpsGui *main = nullptr;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") main = dynamic_cast<LammpsGui *>(widget);
    if (main) main->stop_run();
}

void ChartWindow::saveAs()
{
    if (charts.empty()) return;
    QString defaultname = filename + "." + columns->currentText() + ".png";
    if (filename.isEmpty()) defaultname = columns->currentText() + ".png";
    QString fileName = QFileDialog::getSaveFileName(this, "Save Chart as Image", defaultname,
                                                    "Image Files (*.jpg *.png *.bmp *.ppm)");
    if (!fileName.isEmpty()) {
        int choice = columns->currentData().toInt();
        for (auto &c : charts)
            if (choice == c->get_index()) c->grab().save(fileName);
    }
}

void ChartWindow::exportDat()
{
    if (charts.empty()) return;
    QString defaultname = filename + ".dat";
    if (filename.isEmpty()) defaultname = "lammpsdata.dat";
    QString fileName = QFileDialog::getSaveFileName(this, "Save Chart as Gnuplot data", defaultname,
                                                    "Image Files (*.dat)");
    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QTextStream out(&file);
            constexpr int fw = 16;
            out.setFieldAlignment(QTextStream::AlignRight);
            out.setRealNumberPrecision(8);

            out << "# Thermodynamic data from " << filename << "\n";
            out << "#          Step";
            for (auto &c : charts)
                out << qSetFieldWidth(0) << ' ' << qSetFieldWidth(fw) << c->get_title();
            out << qSetFieldWidth(0) << '\n';

            int lines = charts[0]->get_count();
            for (int i = 0; i < lines; ++i) {
                // timestep
                out << qSetFieldWidth(0) << ' ' << qSetFieldWidth(fw) << charts[0]->get_step(i);
                for (auto &c : charts)
                    out << qSetFieldWidth(0) << ' ' << qSetFieldWidth(fw) << c->get_data(i);
                out << qSetFieldWidth(0) << '\n';
            }
            file.close();
        }
    }
}

void ChartWindow::exportCsv()
{
    if (charts.empty()) return;
    QString defaultname = filename + ".csv";
    if (filename.isEmpty()) defaultname = "lammpsdata.csv";
    QString fileName = QFileDialog::getSaveFileName(this, "Save Chart as CSV data", defaultname,
                                                    "Image Files (*.csv)");
    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QTextStream out(&file);
            out.setRealNumberPrecision(8);

            out << "Step";
            for (auto &c : charts)
                out << ',' << c->get_title();
            out << '\n';

            int lines = charts[0]->get_count();
            for (int i = 0; i < lines; ++i) {
                // timestep
                out << charts[0]->get_step(i);
                for (auto &c : charts)
                    out << ',' << c->get_data(i);
                out << '\n';
            }
            file.close();
        }
    }
}

void ChartWindow::change_chart(int)
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

// event filter to handle "Ambiguous shortcut override" issues
bool ChartWindow::eventFilter(QObject *watched, QEvent *event)
{
    if (event->type() == QEvent::ShortcutOverride) {
        QKeyEvent *keyEvent = dynamic_cast<QKeyEvent *>(event);
        if (!keyEvent) return QWidget::eventFilter(watched, event);
        if (keyEvent->modifiers().testFlag(Qt::ControlModifier) && keyEvent->key() == '/') {
            stop_run();
            event->accept();
            return true;
        }
        if (keyEvent->modifiers().testFlag(Qt::ControlModifier) && keyEvent->key() == 'W') {
            close();
            event->accept();
            return true;
        }
    }
    return QWidget::eventFilter(watched, event);
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
        reset_zoom();
    }
}

/* -------------------------------------------------------------------- */

void ChartViewer::reset_zoom()
{
    auto points = series->points();

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

    // avoid (nearly) empty ranges
    double deltax = xmax - xmin;
    if ((deltax / ((xmax == 0.0) ? 1.0 : xmax)) < 1.0e-10) {
        if ((xmin == 0.0) || (xmax == 0.0)) {
            xmin = -0.025;
            xmax = 0.025;
        } else {
            xmin -= 0.025 * fabs(xmin);
            xmax += 0.025 * fabs(xmax);
        }
    }

    double deltay = ymax - ymin;
    if ((deltay / ((ymax == 0.0) ? 1.0 : ymax)) < 1.0e-10) {
        if ((ymin == 0.0) || (ymax == 0.0)) {
            ymin = -0.025;
            ymax = 0.025;
        } else {
            ymin -= 0.025 * fabs(ymin);
            ymax += 0.025 * fabs(ymax);
        }
    }

    xaxis->setRange(xmin, xmax);
    yaxis->setRange(ymin, ymax);
}

// Local Variables:
// c-basic-offset: 4
// End:
