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
#include <QChart>
#include <QCloseEvent>
#include <QComboBox>
#include <QEvent>
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
#include <QSpinBox>
#include <QTextStream>
#include <QTime>
#include <QVBoxLayout>
#include <QValueAxis>
#include <QVariant>

#include <cmath>

using namespace QtCharts;

ChartWindow::ChartWindow(const QString &_filename, QWidget *parent) :
    QWidget(parent), menu(new QMenuBar), file(new QMenu("&File")), saveAsAct(nullptr),
    exportCsvAct(nullptr), exportDatAct(nullptr), exportYamlAct(nullptr), closeAct(nullptr),
    stopAct(nullptr), quitAct(nullptr), filename(_filename)
{
    auto *top = new QHBoxLayout;
    menu->addMenu(file);
    menu->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Preferred);

    // workaround for incorrect highlight bug on macOS
    auto *dummy = new QPushButton(QIcon(), "");
    dummy->hide();

    do_raw    = true;
    do_smooth = true;
    smooth    = new QComboBox;
    smooth->addItem("Raw");
    smooth->addItem("Smooth");
    smooth->addItem("Both");
    smooth->setCurrentIndex(2);
    smooth->show();
    window = new QSpinBox;
    window->setRange(5, 999);
    window->setValue(10);
    window->setEnabled(true);
    window->setToolTip("Smoothing Window Size");
    order = new QSpinBox;
    order->setRange(1, 20);
    order->setValue(4);
    order->setEnabled(true);
    order->setToolTip("Smoothing Order");

    auto *normal = new QPushButton(QIcon(":/icons/gtk-zoom-fit.png"), "");
    normal->setToolTip("Reset zoom to normal");

    columns = new QComboBox;
    top->addWidget(menu);
    top->addSpacerItem(new QSpacerItem(1, 1, QSizePolicy::Expanding, QSizePolicy::Minimum));
    top->addWidget(dummy);
    top->addWidget(new QLabel("Plot:"));
    top->addWidget(smooth);
    top->addWidget(new QLabel(" Smooth:"));
    top->addWidget(window);
    top->addWidget(order);
    top->addWidget(new QLabel(" "));
    top->addWidget(normal);
    top->addWidget(new QLabel(" Data:"));
    top->addWidget(columns);
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

    connect(smooth, SIGNAL(currentIndexChanged(int)), this, SLOT(select_smooth(int)));
    connect(window, &QAbstractSpinBox::editingFinished, this, &ChartWindow::update_smooth);
    connect(order, &QAbstractSpinBox::editingFinished, this, &ChartWindow::update_smooth);
    connect(window, QOverload<int>::of(&QSpinBox::valueChanged), this, &ChartWindow::update_smooth);
    connect(order, QOverload<int>::of(&QSpinBox::valueChanged), this, &ChartWindow::update_smooth);
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
    if ((choice >= 0) && (choice < charts.size())) {
        charts[choice]->update_smooth();
        charts[choice]->reset_zoom();
    }
}

void ChartWindow::stop_run()
{
    LammpsGui *main = nullptr;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") main = dynamic_cast<LammpsGui *>(widget);
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
    do_raw(true), do_smooth(true)
{
    chart->legend()->hide();
    chart->addAxis(xaxis, Qt::AlignBottom);
    chart->addAxis(yaxis, Qt::AlignLeft);
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

        // do not add data that deviates by more than 4 sigma from the average
        // over the last 5 to 20 data items.  this is a hack to work around
        // getting corrupted data from lammps_get_last_thermo()
        const auto &points = series->points();
        const auto count   = points.count();
        if (count > 4) {
            double ysum   = 0.0;
            double ysumsq = 0.0;
            int first     = count - 20;
            if (first < 0) first = 0;
            for (int i = first; i < count; ++i) {
                double val = points[i].y();
                ysum += val;
                ysumsq += val * val;
            }
            const double num   = count - first;
            const double avg   = ysum / num;
            const double avgsq = ysumsq / num;
            if (fabs(data - avg) > (5.0 * sqrt(avgsq - (avg * avg)))) return;
        }
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

void ChartViewer::reset_zoom()
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

    xaxis->setRange(xmin, xmax);
    yaxis->setRange(ymin, ymax);
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

// update smooth plot data

static QList<QPointF> calc_sgsmooth(const QList<QPointF> &input, const int window, const int order);

void ChartViewer::update_smooth()
{
    auto allseries = chart->series();
    if (do_raw) {
        // add raw data if not in chart
        if (!allseries.contains(series)) {
            series->setPen(QPen(QBrush(QColor(100, 150, 255)), 3, Qt::SolidLine, Qt::RoundCap));
            chart->addSeries(series);
            series->attachAxis(xaxis);
            series->attachAxis(yaxis);
        }
    }

    if (do_smooth) {
        if (series->count() > (2 * window)) {
            if (!smooth) {
                smooth = new QLineSeries;
                smooth->setPen(QPen(QBrush(QColor(255, 125, 125)), 3, Qt::SolidLine, Qt::RoundCap));
                chart->addSeries(smooth);
                smooth->attachAxis(xaxis);
                smooth->attachAxis(yaxis);
            }
            smooth->clear();
            smooth->append(calc_sgsmooth(series->points(), window, order));
        }
    }
}

//! default convergence
static constexpr double TINY_FLOAT = 1.0e-300;

//! comfortable array of doubles
typedef std::vector<double> float_vect;
//! comfortable array of ints;
typedef std::vector<int> int_vect;

// savitzky golay smoothing.
static float_vect sg_smooth(const float_vect &v, const int w, const int deg);

QList<QPointF> calc_sgsmooth(const QList<QPointF> &input, int window, int order)
{
    const int ndat = input.count();
    if (ndat < 2 * window + 2) window = ndat / 2 - 1;

    if (window > 1) {
        float_vect in(ndat);
        QList<QPointF> rv;

        for (int i = 0; i < ndat; ++i) {
            in[i] = input[i].y();
        }
        float_vect out = sg_smooth(in, window, order);

        for (int i = 0; i < ndat; ++i) {
            rv.append(QPointF(input[i].x(), out[i]));
        }
        return rv;
    } else {
        return input;
    }
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
private:
    //! disable the default constructor
    explicit float_mat() {};
    //! disable assignment operator until it is implemented.
    float_mat &operator=(const float_mat &) { return *this; };

public:
    //! constructor with sizes
    float_mat(const std::size_t rows, const std::size_t cols, const double def = 0.0);
    //! copy constructor for matrix
    float_mat(const float_mat &m);
    //! copy constructor for vector
    float_mat(const float_vect &v);

    //! use default destructor
    // ~float_mat() {};

    //! get size
    int nr_rows(void) const { return size(); };
    //! get size
    int nr_cols(void) const { return front().size(); };
};

// constructor with sizes
float_mat::float_mat(const std::size_t rows, const std::size_t cols, const double defval) :
    std::vector<float_vect>(rows)
{
    for (std::size_t i = 0; i < rows; ++i) {
        (*this)[i].resize(cols, defval);
    }
#if 0
    if ((rows < 1) || (cols < 1)) {
        char buffer[1024];

        sprintf(buffer, "cannot build matrix with %d rows and %d columns\n",
                rows, cols);
        sgs_error(buffer);
    }
#endif
}

// copy constructor for matrix
float_mat::float_mat(const float_mat &m) : std::vector<float_vect>(m.size())
{

    float_mat::iterator inew       = begin();
    float_mat::const_iterator iold = m.begin();
    for (/* empty */; iold < m.end(); ++inew, ++iold) {
        const auto oldsz = iold->size();
        inew->resize(oldsz);
        const float_vect oldvec(*iold);
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
static int partial_pivot(float_mat &A, const std::size_t row, const std::size_t col,
                         float_vect &scale, int_vect &idx, double tol)
{
    if (tol <= 0.0) tol = TINY_FLOAT;

    int swapNum = 1;

    // default pivot is the current position, [row,col]
    std::size_t pivot = row;
    double piv_elem   = fabs(A[idx[row]][col]) * scale[idx[row]];

    // loop over possible pivots below current
    for (int j = row + 1; j < A.nr_rows(); ++j) {

        const double tmp = fabs(A[idx[j]][col]) * scale[idx[j]];

        // if this elem is larger, then it becomes the pivot
        if (tmp > piv_elem) {
            pivot    = j;
            piv_elem = tmp;
        }
    }

#if 0
    if(piv_elem < tol) {
        sgs_error("partial_pivot(): Zero pivot encountered.\n")
#endif

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
static void lu_backsubst(float_mat &A, float_mat &a, bool diag = false)
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
static void lu_forwsubst(float_mat &A, float_mat &a, bool diag = true)
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
static int lu_factorize(float_mat &A, int_vect &idx, double tol = TINY_FLOAT)
{
    if (tol <= 0.0) tol = TINY_FLOAT;
#if 0
    if ((A.nr_rows() == 0) || (A.nr_rows() != A.nr_cols())) {
        sgs_error("lu_factorize(): cannot handle empty "
                  "or nonsquare matrices.\n");

        return 0;
    }
#endif
    float_vect scale(A.nr_rows()); // implicit pivot scaling
    for (int i = 0; i < A.nr_rows(); ++i) {
        double maxval = 0.0;
        for (int j = 0; j < A.nr_cols(); ++j) {
            if (fabs(A[i][j]) > maxval) maxval = fabs(A[i][j]);
        }
        if (maxval == 0.0) {
#if 0
            sgs_error("lu_factorize(): zero pivot found.\n");
#endif
            return 0;
        }
        scale[i] = 1.0 / maxval;
    }

    int swapNum = 1;
    for (int c = 0; c < A.nr_cols(); ++c) {                 // loop over columns
        swapNum *= partial_pivot(A, c, c, scale, idx, tol); // bring pivot to diagonal
        for (int r = 0; r < A.nr_rows(); ++r) {             //  loop over rows
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
static float_mat lin_solve(const float_mat &A, const float_mat &a, double tol = TINY_FLOAT)
{
    float_mat B(A);
    float_mat b(a);
    int_vect idx(B.nr_rows());

    for (int j = 0; j < B.nr_rows(); ++j) {
        idx[j] = j; // init row swap label array
    }
    lu_factorize(B, idx, tol); // get the lu-decomp.
    permute(b, idx);           // sort the inhomogeneity to match the lu-decomp
    lu_forwsubst(B, b);        // solve the forward problem
    lu_backsubst(B, b);        // solve the backward problem
    return b;
}

///////////////////////
// related functions //
///////////////////////

//! Returns the inverse of a matrix using LU-decomposition.
static float_mat invert(const float_mat &A)
{
    const int n = A.size();
    float_mat E(n, n, 0.0);
    float_mat B(A);

    for (int i = 0; i < n; ++i) {
        E[i][i] = 1.0;
    }

    return lin_solve(B, E);
}

//! returns the transposed matrix.
static float_mat transpose(const float_mat &a)
{
    float_mat res(a.nr_cols(), a.nr_rows());

    for (int i = 0; i < a.nr_rows(); ++i) {
        for (int j = 0; j < a.nr_cols(); ++j) {
            res[j][i] = a[i][j];
        }
    }
    return res;
}

//! matrix multiplication.
float_mat operator*(const float_mat &a, const float_mat &b)
{
    float_mat res(a.nr_rows(), b.nr_cols());
#if 0
    if (a.nr_cols() != b.nr_rows()) {
        sgs_error("incompatible matrices in multiplication\n");
        return res;
    }
#endif
    for (int i = 0; i < a.nr_rows(); ++i) {
        for (int j = 0; j < b.nr_cols(); ++j) {
            double sum(0.0);
            for (int k = 0; k < a.nr_cols(); ++k) {
                sum += a[i][k] * b[k][j];
            }
            res[i][j] = sum;
        }
    }
    return res;
}

//! calculate savitzky golay coefficients.
static float_vect sg_coeff(const float_vect &b, const std::size_t deg)
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
#if 0
    if ((width < 1) || (deg < 0) || (v.size() < (2 * width + 2))) {
        sgs_error("sgsmooth: parameter error.\n");
        return res;
    }
#endif
    const int window = 2 * width + 1;
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

// Local Variables:
// c-basic-offset: 4
// End:
