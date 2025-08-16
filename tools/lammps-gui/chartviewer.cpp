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

#include "helpers.h"
#include "lammpsgui.h"
#include "rangeslider.h"

#include <QAction>
#include <QApplication>
#include <QChart>
#include <QCheckBox>
#include <QCloseEvent>
#include <QComboBox>
#include <QEvent>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QKeySequence>
#include <QLabel>
#include <QLayout>
#include <QLineSeries>
#include <QList>
#include <QMenu>
#include <QMenuBar>
#include <QPushButton>
#include <QSettings>
#include <QSpacerItem>
#include <QSpinBox>
#include <QTextStream>
#include <QTime>
#include <QVBoxLayout>
#include <QValueAxis>
#include <QVariant>

#include <cmath>

using namespace QtCharts;

namespace {

// Set RangeSlider resolution to 1000 steps
constexpr int SLIDER_RANGE       = 1000;
constexpr double SLIDER_FRACTION = 1.0 / (double)SLIDER_RANGE;

// brush color index must be kept in sync with preferences

const QList<QBrush> mybrushes = {
    QBrush(QColor(0, 0, 0)),       // black
    QBrush(QColor(100, 150, 255)), // blue
    QBrush(QColor(255, 125, 125)), // red
    QBrush(QColor(100, 200, 100)), // green
    QBrush(QColor(120, 120, 120)), // grey
};

// convenience class

class QHline : public QFrame {
public:
    QHline(QWidget *parent = nullptr) : QFrame(parent)
    {
        setGeometry(QRect(0, 0, 100, 3));
        setFrameShape(QFrame::HLine);
        setFrameShadow(QFrame::Sunken);
    }
};
} // namespace

ChartWindow::ChartWindow(const QString &_filename, QWidget *parent) :
    QWidget(parent), menu(new QMenuBar), file(new QMenu("&File")), saveAsAct(nullptr),
    exportCsvAct(nullptr), exportDatAct(nullptr), exportYamlAct(nullptr), closeAct(nullptr),
    stopAct(nullptr), quitAct(nullptr), smooth(nullptr), window(nullptr), order(nullptr),
    chartTitle(nullptr), chartYlabel(nullptr), units(nullptr), norm(nullptr), filename(_filename)
{
    QSettings settings;
    auto *top  = new QVBoxLayout;
    auto *row1 = new QHBoxLayout;
    auto *row2 = new QHBoxLayout;
    top->addLayout(row1);
    top->addWidget(new QHline);
    top->addLayout(row2);
    top->addWidget(new QHline);

    menu->addMenu(file);
    menu->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);

    // workaround for incorrect highlight bug on macOS
    auto *dummy = new QPushButton(QIcon(), "");
    dummy->hide();

    // plot title and axis labels
    settings.beginGroup("charts");
    chartTitle =
        new QLineEdit(settings.value("title", "Thermo: %f").toString().replace("%f", filename));
    chartYlabel = new QLineEdit("");

    // plot smoothing
    int smoothchoice = settings.value("smoothchoice", 0).toInt();
    switch (smoothchoice) {
        case 0:
            do_raw    = true;
            do_smooth = false;
            break;
        case 1:
            do_raw    = false;
            do_smooth = true;
            break;
        case 2: // fallthrough
        default:
            do_raw    = true;
            do_smooth = true;
            break;
    }
    // list of choices must be kepy in sync with list in preferences
    smooth = new QComboBox;
    smooth->addItem("Raw");
    smooth->addItem("Smooth");
    smooth->addItem("Both");
    smooth->setCurrentIndex(smoothchoice);
    smooth->show();
    window = new QSpinBox;
    window->setRange(5, 999);
    window->setValue(settings.value("smoothwindow", 10).toInt());
    window->setEnabled(true);
    window->setToolTip("Smoothing Window Size");
    order = new QSpinBox;
    order->setRange(1, 20);
    order->setValue(settings.value("smoothorder", 4).toInt());
    order->setEnabled(true);
    order->setToolTip("Smoothing Order");
    settings.endGroup();

    columns = new QComboBox;
    row1->addWidget(menu);
    row1->addWidget(dummy);
    row2->addWidget(dummy);
    row1->addWidget(new QLabel("Title:"));
    row1->addWidget(chartTitle);
    row1->addWidget(new QLabel("Y-Axis:"));
    row1->addWidget(chartYlabel);

    units = new QLabel("Units:");
    row2->addWidget(units);
    row2->addWidget(new QLabel("Norm:"));
    norm = new QCheckBox("");
    norm->setChecked(Qt::Unchecked);
    norm->setEnabled(false);
    row2->addWidget(norm);
    xrange = new RangeSlider;
    xrange->setMinimum(0);
    xrange->setMaximum(SLIDER_RANGE);
    xrange->setLow(0);
    xrange->setHigh(SLIDER_RANGE);
    xrange->setToolTip("Adjust x-axis data range");
    xrange->setTickPosition(QSlider::TicksBothSides);
    xrange->setTickInterval(100);
    yrange = new RangeSlider;
    yrange->setMinimum(0);
    yrange->setMaximum(SLIDER_RANGE);
    yrange->setLow(0);
    yrange->setHigh(SLIDER_RANGE);
    yrange->setToolTip("Adjust y-axis data range");
    yrange->setTickPosition(QSlider::TicksBothSides);
    yrange->setTickInterval(100);
    row2->addWidget(new QLabel("X:"));
    row2->addWidget(xrange);
    row2->addWidget(new QLabel("Y:"));
    row2->addWidget(yrange);
    row2->addWidget(new QLabel("Plot:"));
    row2->addWidget(smooth);
    row2->addWidget(new QLabel(" Smooth:"));
    row2->addWidget(window);
    row2->addWidget(order);
    row1->addWidget(new QLabel(" Data:"));
    row1->addWidget(columns);
    saveAsAct = file->addAction("&Save Graph As...", this, &ChartWindow::saveAs);
    saveAsAct->setIcon(QIcon(":/icons/document-save-as.png"));
    exportCsvAct = file->addAction("&Export data to CSV...", this, &ChartWindow::exportCsv);
    exportCsvAct->setIcon(QIcon(":/icons/application-calc.png"));
    exportDatAct = file->addAction("Export data to &Gnuplot...", this, &ChartWindow::exportDat);
    exportDatAct->setIcon(QIcon(":/icons/application-plot.png"));
    exportYamlAct = file->addAction("Export data to &YAML...", this, &ChartWindow::exportYaml);
    exportYamlAct->setIcon(QIcon(":/icons/yaml-file-icon.png"));
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

    connect(chartTitle, &QLineEdit::editingFinished, this, &ChartWindow::update_tlabel);
    connect(chartYlabel, &QLineEdit::editingFinished, this, &ChartWindow::update_ylabel);
    connect(smooth, SIGNAL(currentIndexChanged(int)), this, SLOT(select_smooth(int)));
    connect(window, &QAbstractSpinBox::editingFinished, this, &ChartWindow::update_smooth);
    connect(order, &QAbstractSpinBox::editingFinished, this, &ChartWindow::update_smooth);
    connect(window, QOverload<int>::of(&QSpinBox::valueChanged), this, &ChartWindow::update_smooth);
    connect(order, QOverload<int>::of(&QSpinBox::valueChanged), this, &ChartWindow::update_smooth);
    connect(columns, SIGNAL(currentIndexChanged(int)), this, SLOT(change_chart(int)));
    connect(xrange, &RangeSlider::sliderMoved, this, &ChartWindow::update_xrange);
    connect(yrange, &RangeSlider::sliderMoved, this, &ChartWindow::update_yrange);

    installEventFilter(this);
    resize(settings.value("chartx", 640).toInt(), settings.value("charty", 480).toInt());
}

int ChartWindow::get_step() const
{
    if (!charts.empty()) {
        auto *v = charts[0];
        if (v) {
            return (int)v->get_step(v->get_count() - 1);
        }
    }
    return -1;
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
    if (!charts.empty()) {
        chart->hide();
    } else {
        // must initialize QLineEdit with first title
        // will be automatically updated when changing charts.
        chartYlabel->setText(title);
    }
    charts.append(chart);
    update_tlabel();
    select_smooth(0);
}

void ChartWindow::add_data(int step, double data, int index)
{
    for (auto &c : charts)
        if (c->get_index() == index) c->add_data(step, data);
}

void ChartWindow::set_units(const QString &_units)
{
    units->setText(_units);
}

void ChartWindow::set_norm(bool _norm)
{
    norm->setChecked(_norm ? Qt::Checked : Qt::Unchecked);
}

void ChartWindow::quit()
{
    auto *main = dynamic_cast<LammpsGui *>(get_main_widget());
    if (main) main->quit();
}

void ChartWindow::stop_run()
{
    auto *main = dynamic_cast<LammpsGui *>(get_main_widget());
    if (main) main->stop_run();
}

void ChartWindow::select_smooth(int)
{
    switch (smooth->currentIndex()) {
        case 0:
            do_raw    = true;
            do_smooth = false;
            break;
        case 1:
            do_raw    = false;
            do_smooth = true;
            break;
        case 2: // fallthrough
        default:
            do_raw    = true;
            do_smooth = true;
            break;
    }
    window->setEnabled(do_smooth);
    order->setEnabled(do_smooth);
    update_smooth();
}

void ChartWindow::update_smooth()
{
    int wval = window->value();
    int oval = order->value();

    for (auto &c : charts)
        c->smooth_param(do_raw, do_smooth, wval, oval);
}

void ChartWindow::update_tlabel()
{
    for (auto &c : charts)
        c->set_tlabel(chartTitle->text());
}

void ChartWindow::update_ylabel()
{
    for (auto &c : charts) {
        if (c->isVisible()) c->set_ylabel(chartYlabel->text());
    }
}

void ChartWindow::update_xrange(int low, int high)
{
    for (auto &c : charts) {
        if (c->isVisible()) {
            auto axes   = c->get_axes();
            auto ranges = c->get_minmax();
            double xmin = ranges.left() + (double)low * SLIDER_FRACTION * ranges.width();
            double xmax = ranges.left() + (double)high * SLIDER_FRACTION * ranges.width();
            axes[0]->setRange(xmin, xmax);
        }
    }
}

void ChartWindow::update_yrange(int low, int high)
{
    for (auto &c : charts) {
        if (c->isVisible()) {
            constexpr double fraction = 1.0 / (double)SLIDER_RANGE;
            auto axes                 = c->get_axes();
            auto ranges               = c->get_minmax();
            double ymin = ranges.bottom() - (double)low * SLIDER_FRACTION * ranges.height();
            double ymax = ranges.bottom() - (double)high * SLIDER_FRACTION * ranges.height();
            axes[1]->setRange(ymin, ymax);
        }
    }
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
void ChartWindow::exportYaml()
{
    if (charts.empty()) return;
    QString defaultname = filename + ".yaml";
    if (filename.isEmpty()) defaultname = "lammpsdata.yaml";
    QString fileName = QFileDialog::getSaveFileName(this, "Save Chart as YAML data", defaultname,
                                                    "Image Files (*.yaml, *.yml)");
    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QTextStream out(&file);
            out.setRealNumberPrecision(8);
            out << "---\n";

            out << "keywords: ['Step'";
            for (auto &c : charts)
                out << ", " << c->get_title();
            out << "]\n";

            out << "data: \n";
            int lines = charts[0]->get_count();
            for (int i = 0; i < lines; ++i) {
                // timestep
                out << "  - [" << charts[0]->get_step(i);
                // data
                for (auto &c : charts)
                    out << ", " << c->get_data(i);
                out << "]\n";
            }
            out << "...\n";
            file.close();
        }
    }
}

void ChartWindow::change_chart(int)
{
    int choice = columns->currentData().toInt();
    for (auto &c : charts) {
        if (choice == c->get_index()) {
            c->show();
            chartTitle->setText(c->get_tlabel());
            chartYlabel->setText(c->get_ylabel());
        } else {
            c->hide();
        }
    }

    // reset plot range selection
    xrange->setLow(0);
    xrange->setHigh(SLIDER_RANGE);
    yrange->setLow(0);
    yrange->setHigh(SLIDER_RANGE);
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
        auto *keyEvent = dynamic_cast<QKeyEvent *>(event);
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
    QChartView(parent), last_step(-1), index(_index), window(10), order(4), chart(new QChart),
    series(new QLineSeries), smooth(nullptr), xaxis(new QValueAxis), yaxis(new QValueAxis),
    do_raw(true), do_smooth(false)
{
    chart->legend()->hide();
    chart->addAxis(xaxis, Qt::AlignBottom);
    chart->addAxis(yaxis, Qt::AlignLeft);
    chart->setTitle("");
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
    setRubberBand(QChartView::NoRubberBand);
    last_update = QTime::currentTime();
    update_smooth();
}

/* -------------------------------------------------------------------- */

ChartViewer::~ChartViewer()
{
    delete xaxis;
    delete yaxis;
    delete smooth;
    delete series;
    delete chart;
}

/* -------------------------------------------------------------------- */

void ChartViewer::add_data(int step, double data)
{
    if (last_step < step) {
        last_step = step;
        series->append(step, data);

        QSettings settings;
        // update the chart display only after at least updchart milliseconds have passed
        if (last_update.msecsTo(QTime::currentTime()) > settings.value("updchart", "500").toInt()) {
            last_update = QTime::currentTime();
            update_smooth();
            reset_zoom();
        }
    }
}

/* -------------------------------------------------------------------- */

QRectF ChartViewer::get_minmax() const
{
    auto points = series->points();

    // get min/max for plot
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

    // if plotting the smoothed plot, check for its min/max values, too
    if (smooth) {
        auto spoints = smooth->points();
        for (auto &p : spoints) {
            xmin = qMin(xmin, p.x());
            xmax = qMax(xmax, p.x());
            ymin = qMin(ymin, p.y());
            ymax = qMax(ymax, p.y());
        }
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

    return QRectF(xmin, ymax, xmax - xmin, ymin - ymax);
}

/* -------------------------------------------------------------------- */

void ChartViewer::reset_zoom()
{
    auto ranges = get_minmax();
    xaxis->setRange(ranges.left(), ranges.right());
    yaxis->setRange(ranges.bottom(), ranges.top());
}

/* -------------------------------------------------------------------- */

void ChartViewer::smooth_param(bool _do_raw, bool _do_smooth, int _window, int _order)
{
    // turn off raw plot
    if (!_do_raw) {
        if (do_raw) chart->removeSeries(series);
    }
    // turn off smooth plot
    if (!_do_smooth) {
        if (smooth) {
            chart->removeSeries(smooth);
            delete smooth;
            smooth = nullptr;
        }
    }
    do_raw    = _do_raw;
    do_smooth = _do_smooth;
    window    = _window;
    order     = _order;
    update_smooth();
}

/* -------------------------------------------------------------------- */

void ChartViewer::set_tlabel(const QString &tlabel)
{
    chart->setTitle(tlabel);
}

/* -------------------------------------------------------------------- */

void ChartViewer::set_ylabel(const QString &ylabel)
{
    yaxis->setTitleText(ylabel);
}

// local implementation of Savitzky-Golay filter

namespace {

//! array of doubles
using float_vect = std::vector<double>;

//! array of ints;
using int_vect = std::vector<int>;

// forward declaration
float_vect sg_smooth(const float_vect &v, const int w, const int deg);

// savitzky golay smoothing.
QList<QPointF> calc_sgsmooth(const QList<QPointF> &input, int window, int order)
{
    const std::size_t ndat = input.count();
    if (ndat < ((2 * window) + 2)) window = (ndat / 2) - 1;

    if (window > 1) {
        float_vect in(ndat);
        QList<QPointF> rv(input);

        for (int i = 0; i < ndat; ++i)
            in[i] = input[i].y();

        float_vect out = sg_smooth(in, window, order);

        for (int i = 0; i < ndat; ++i)
            rv[i].setY(out[i]);

        return rv;
    }
    return input;
}

/*! matrix class.
 *
 * This is a matrix class derived from a vector of float_vects.  Note that
 * the matrix elements indexed [row][column] with indices starting at 0 (c
 * style). Also note that because of its design looping through rows should
 * be faster than looping through columns.
 *
 * \brief two dimensional floating point array
 */
class float_mat : public std::vector<float_vect> {

public:
    // disable selected default constructors and assignment operators
    float_mat()                             = delete;
    float_mat(float_mat &&)                 = default;
    ~float_mat()                            = default;
    float_mat &operator=(const float_mat &) = delete;
    float_mat &operator=(float_mat &&)      = delete;

    //! constructor with sizes
    float_mat(const std::size_t rows, const std::size_t cols, const double def = 0.0);
    //! copy constructor for matrix
    float_mat(const float_mat &m);
    //! copy constructor for vector
    float_mat(const float_vect &v);

    //! get size
    std::size_t nr_rows() const { return size(); };
    //! get size
    std::size_t nr_cols() const { return front().size(); };
};

// constructor with sizes
float_mat::float_mat(const std::size_t rows, const std::size_t cols, const double defval) :
    std::vector<float_vect>(rows)
{
    for (std::size_t i = 0; i < rows; ++i) {
        (*this)[i].resize(cols, defval);
    }
}

// copy constructor for matrix
float_mat::float_mat(const float_mat &m) : std::vector<float_vect>(m.size())
{

    auto inew = begin();
    auto iold = m.begin();
    for (/* empty */; iold < m.end(); ++inew, ++iold) {
        const auto oldsz = iold->size();
        inew->resize(oldsz);
        const float_vect &oldvec(*iold);
        *inew = oldvec;
    }
}

// copy constructor for vector
float_mat::float_mat(const float_vect &v) : std::vector<float_vect>(1)
{

    const auto oldsz = v.size();
    front().resize(oldsz);
    front() = v;
}

//////////////////////
// Helper functions //
//////////////////////

//! permute() orders the rows of A to match the integers in the index array.
void permute(float_mat &A, int_vect &idx)
{
    int_vect i(idx.size());

    for (int j = 0; j < A.nr_rows(); ++j) {
        i[j] = j;
    }

    // loop over permuted indices
    for (int j = 0; j < A.nr_rows(); ++j) {
        if (i[j] != idx[j]) {

            // search only the remaining indices
            for (int k = j + 1; k < A.nr_rows(); ++k) {
                if (i[k] == idx[j]) {
                    std::swap(A[j], A[k]); // swap the rows and
                    i[k] = i[j];           // the elements of
                    i[j] = idx[j];         // the ordered index.
                    break;                 // next j
                }
            }
        }
    }
}

/*! \brief Implicit partial pivoting.
 *
 * The function looks for pivot element only in rows below the current
 * element, A[idx[row]][column], then swaps that row with the current one in
 * the index map. The algorithm is for implicit pivoting (i.e., the pivot is
 * chosen as if the max coefficient in each row is set to 1) based on the
 * scaling information in the vector scale. The map of swapped indices is
 * recorded in swp. The return value is +1 or -1 depending on whether the
 * number of row swaps was even or odd respectively. */
int partial_pivot(float_mat &A, const std::size_t row, const std::size_t col, float_vect &scale,
                  int_vect &idx)
{
    int swapNum = 1;

    // default pivot is the current position, [row,col]
    std::size_t pivot = row;
    double piv_elem   = fabs(A[idx[row]][col]) * scale[idx[row]];

    // loop over possible pivots below current
    for (std::size_t j = row + 1; j < A.nr_rows(); ++j) {

        const double tmp = fabs(A[idx[j]][col]) * scale[idx[j]];

        // if this elem is larger, then it becomes the pivot
        if (tmp > piv_elem) {
            pivot    = j;
            piv_elem = tmp;
        }
    }

    if (pivot > row) {         // bring the pivot to the diagonal
        int j      = idx[row]; // reorder swap array
        idx[row]   = idx[pivot];
        idx[pivot] = j;
        swapNum    = -swapNum; // keeping track of odd or even swap
    }
    return swapNum;
}

/*! \brief Perform backward substitution.
 *
 * Solves the system of equations A*b=a, ASSUMING that A is upper
 * triangular. If diag==1, then the diagonal elements are additionally
 * assumed to be 1.  Note that the lower triangular elements are never
 * checked, so this function is valid to use after a LU-decomposition in
 * place.  A is not modified, and the solution, b, is returned in a. */
void lu_backsubst(float_mat &A, float_mat &a, bool diag = false)
{
    for (int r = (A.nr_rows() - 1); r >= 0; --r) {
        for (int c = (A.nr_cols() - 1); c > r; --c) {
            for (int k = 0; k < A.nr_cols(); ++k) {
                a[r][k] -= A[r][c] * a[c][k];
            }
        }
        if (!diag) {
            for (int k = 0; k < A.nr_cols(); ++k) {
                a[r][k] /= A[r][r];
            }
        }
    }
}

/*! \brief Perform forward substitution.
 *
 * Solves the system of equations A*b=a, ASSUMING that A is lower
 * triangular. If diag==1, then the diagonal elements are additionally
 * assumed to be 1.  Note that the upper triangular elements are never
 * checked, so this function is valid to use after a LU-decomposition in
 * place.  A is not modified, and the solution, b, is returned in a. */
void lu_forwsubst(float_mat &A, float_mat &a, bool diag = true)
{
    for (int r = 0; r < A.nr_rows(); ++r) {
        for (int c = 0; c < r; ++c) {
            for (int k = 0; k < A.nr_cols(); ++k) {
                a[r][k] -= A[r][c] * a[c][k];
            }
        }
        if (!diag) {
            for (int k = 0; k < A.nr_cols(); ++k) {
                a[r][k] /= A[r][r];
            }
        }
    }
}

/*! \brief Performs LU factorization in place.
 *
 * This is Crout's algorithm (cf., Num. Rec. in C, Section 2.3).  The map of
 * swapped indeces is recorded in idx. The return value is +1 or -1
 * depending on whether the number of row swaps was even or odd
 * respectively.  idx must be preinitialized to a valid set of indices
 * (e.g., {1,2, ... ,A.nr_rows()}). */
int lu_factorize(float_mat &A, int_vect &idx)
{
    float_vect scale(A.nr_rows()); // implicit pivot scaling
    for (int i = 0; i < A.nr_rows(); ++i) {
        double maxval = 0.0;
        for (int j = 0; j < A.nr_cols(); ++j) {
            maxval = std::max(fabs(A[i][j]), maxval);
        }
        if (maxval == 0.0) {
            return 0;
        }
        scale[i] = 1.0 / maxval;
    }

    int swapNum = 1;
    for (int c = 0; c < A.nr_cols(); ++c) {            // loop over columns
        swapNum *= partial_pivot(A, c, c, scale, idx); // bring pivot to diagonal
        for (int r = 0; r < A.nr_rows(); ++r) {        //  loop over rows
            int lim = (r < c) ? r : c;
            for (int j = 0; j < lim; ++j) {
                A[idx[r]][c] -= A[idx[r]][j] * A[idx[j]][c];
            }
            if (r > c) A[idx[r]][c] /= A[idx[c]][c];
        }
    }
    permute(A, idx);
    return swapNum;
}

/*! \brief Solve a system of linear equations.
 * Solves the inhomogeneous matrix problem with lu-decomposition. Note that
 * inversion may be accomplished by setting a to the identity_matrix. */
float_mat lin_solve(const float_mat &A, const float_mat &a)
{
    float_mat B(A);
    float_mat b(a);
    int_vect idx(B.nr_rows());

    for (int j = 0; j < B.nr_rows(); ++j) {
        idx[j] = j; // init row swap label array
    }
    lu_factorize(B, idx); // get the lu-decomp.
    permute(b, idx);      // sort the inhomogeneity to match the lu-decomp
    lu_forwsubst(B, b);   // solve the forward problem
    lu_backsubst(B, b);   // solve the backward problem
    return b;
}

///////////////////////
// related functions //
///////////////////////

//! Returns the inverse of a matrix using LU-decomposition.
float_mat invert(const float_mat &A)
{
    const std::size_t n = A.size();
    float_mat E(n, n, 0.0);
    const float_mat &B(A);

    for (std::size_t i = 0; i < n; ++i) {
        E[i][i] = 1.0;
    }

    return lin_solve(B, E);
}

//! returns the transposed matrix.
float_mat transpose(const float_mat &a)
{
    float_mat res(a.nr_cols(), a.nr_rows());

    for (std::size_t i = 0; i < a.nr_rows(); ++i) {
        for (std::size_t j = 0; j < a.nr_cols(); ++j) {
            res[j][i] = a[i][j];
        }
    }
    return res;
}

//! matrix multiplication.
float_mat operator*(const float_mat &a, const float_mat &b)
{
    float_mat res(a.nr_rows(), b.nr_cols());
    for (std::size_t i = 0; i < a.nr_rows(); ++i) {
        for (std::size_t j = 0; j < b.nr_cols(); ++j) {
            double sum(0.0);
            for (std::size_t k = 0; k < a.nr_cols(); ++k) {
                sum += a[i][k] * b[k][j];
            }
            res[i][j] = sum;
        }
    }
    return res;
}

//! calculate savitzky golay coefficients.
float_vect sg_coeff(const float_vect &b, const std::size_t deg)
{
    const std::size_t rows(b.size());
    const std::size_t cols(deg + 1);
    float_mat A(rows, cols);
    float_vect res(rows);

    // generate input matrix for least squares fit
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            A[i][j] = pow(double(i), double(j));
        }
    }

    float_mat c(invert(transpose(A) * A) * (transpose(A) * transpose(b)));

    for (std::size_t i = 0; i < b.size(); ++i) {
        res[i] = c[0][0];
        for (std::size_t j = 1; j <= deg; ++j) {
            res[i] += c[j][0] * pow(double(i), double(j));
        }
    }
    return res;
}

/*! \brief savitzky golay smoothing.
 *
 * This method means fitting a polynome of degree 'deg' to a sliding window
 * of width 2w+1 throughout the data.  The needed coefficients are
 * generated dynamically by doing a least squares fit on a "symmetric" unit
 * vector of size 2w+1, e.g. for w=2 b=(0,0,1,0,0). evaluating the polynome
 * yields the sg-coefficients.  at the border non symmectric vectors b are
 * used. */
float_vect sg_smooth(const float_vect &v, const int width, const int deg)
{
    float_vect res(v.size(), 0.0);
    const int window = (2 * width) + 1;
    const int endidx = v.size() - 1;

    // do a regular sliding window average
    if (deg == 0) {
        // handle border cases first because we need different coefficients
        for (int i = 0; i < width; ++i) {
            const double scale = 1.0 / double(i + 1);
            const float_vect c1(width, scale);
            for (int j = 0; j <= i; ++j) {
                res[i] += c1[j] * v[j];
                res[endidx - i] += c1[j] * v[endidx - j];
            }
        }

        // now loop over rest of data. reusing the "symmetric" coefficients.
        const double scale = 1.0 / double(window);
        const float_vect c2(window, scale);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t i = 0; i <= (v.size() - window); ++i) {
            for (int j = 0; j < window; ++j) {
                res[i + width] += c2[j] * v[i + j];
            }
        }
        return res;
    }

    // handle border cases first because we need different coefficients
    for (int i = 0; i < width; ++i) {
        float_vect b1(window, 0.0);
        b1[i] = 1.0;

        const float_vect c1(sg_coeff(b1, deg));
        for (int j = 0; j < window; ++j) {
            res[i] += c1[j] * v[j];
            res[endidx - i] += c1[j] * v[endidx - j];
        }
    }

    // now loop over rest of data. reusing the "symmetric" coefficients.
    float_vect b2(window, 0.0);
    b2[width] = 1.0;
    const float_vect c2(sg_coeff(b2, deg));

    for (std::size_t i = 0; i <= (v.size() - window); ++i) {
        for (int j = 0; j < window; ++j) {
            res[i + width] += c2[j] * v[i + j];
        }
    }
    return res;
}
} // namespace

/* -------------------------------------------------------------------- */

// update smooth plot data

void ChartViewer::update_smooth()
{
    QSettings settings;
    settings.beginGroup("charts");
    int rawidx    = settings.value("rawbrush", 1).toInt();
    int smoothidx = settings.value("smoothbrush", 2).toInt();
    if ((rawidx < 0) || (rawidx >= mybrushes.size())) rawidx = 0;
    if ((smoothidx < 0) || (smoothidx >= mybrushes.size())) smoothidx = 0;
    settings.endGroup();

    auto allseries = chart->series();
    if (do_raw) {
        // add raw data if not in chart
        if (!allseries.contains(series)) {
            series->setPen(QPen(mybrushes[rawidx], 3, Qt::SolidLine, Qt::RoundCap));
            chart->addSeries(series);
            series->attachAxis(xaxis);
            series->attachAxis(yaxis);
        }
    }

    if (do_smooth) {
        if (series->count() > (2 * window)) {
            if (!smooth) {
                smooth = new QLineSeries;
                smooth->setPen(QPen(mybrushes[smoothidx], 3, Qt::SolidLine, Qt::RoundCap));
                chart->addSeries(smooth);
                smooth->attachAxis(xaxis);
                smooth->attachAxis(yaxis);
            }
            smooth->replace(calc_sgsmooth(series->points(), window, order));
        }
    }
}

// Local Variables:
// c-basic-offset: 4
// End:
