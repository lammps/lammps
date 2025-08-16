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

#include "preferences.h"

#include "helpers.h"
#include "lammpsgui.h"
#include "lammpswrapper.h"
#include "ui_lammpsgui.h"

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QCoreApplication>
#include <QDialogButtonBox>
#include <QDoubleValidator>
#include <QFileDialog>
#include <QFontDialog>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHash>
#include <QIcon>
#include <QIntValidator>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QRadioButton>
#include <QSettings>
#include <QSpacerItem>
#include <QSpinBox>
#include <QTabWidget>
#if defined(_OPENMP)
#include <QThread>
#endif
#include <QVBoxLayout>

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <process.h>
#define execl(exe, arg0, arg1) _execl(exe, arg0, arg1)
#else
#include <unistd.h>
#endif

// convenience class
namespace {
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

Preferences::Preferences(LammpsWrapper *_lammps, QWidget *parent) :
    QDialog(parent), tabWidget(new QTabWidget),
    buttonBox(new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel)),
    settings(new QSettings), lammps(_lammps), need_relaunch(false)
{
    tabWidget->addTab(new GeneralTab(settings, lammps), "&General Settings");
    tabWidget->addTab(new AcceleratorTab(settings, lammps), "&Accelerators");
    tabWidget->addTab(new SnapshotTab(settings), "&Snapshot Image");
    tabWidget->addTab(new EditorTab(settings), "&Editor Settings");
    tabWidget->addTab(new ChartsTab(settings), "Cha&rts Settings");

    connect(buttonBox, &QDialogButtonBox::accepted, this, &Preferences::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);

    auto *layout = new QVBoxLayout;
    layout->addWidget(tabWidget);
    layout->addWidget(buttonBox);
    setLayout(layout);
    setWindowIcon(QIcon(":/icons/lammps-icon-128x128.png"));
    setWindowTitle("LAMMPS-GUI - Preferences");
    resize(700, 500);
}

Preferences::~Preferences()
{
    delete buttonBox;
    delete tabWidget;
    delete settings;
}

namespace {
const QHash<QString, int> buttonToChoice = {
    {"none", AcceleratorTab::None},     {"opt", AcceleratorTab::Opt},
    {"openmp", AcceleratorTab::OpenMP}, {"intel", AcceleratorTab::Intel},
    {"kokkos", AcceleratorTab::Kokkos}, {"gpu", AcceleratorTab::Gpu}};

const QHash<QString, int> buttonToPrecision = {{"inteldouble", AcceleratorTab::Double},
                                               {"intelmixed", AcceleratorTab::Mixed},
                                               {"intelsingle", AcceleratorTab::Single}};
} // namespace

void Preferences::accept()
{
    // store all data in settings class
    // and then confirm accepting

    // store selected accelerator and precision settings from radiobuttons
    QList<QRadioButton *> allButtons = tabWidget->findChildren<QRadioButton *>();
    for (const auto &anyButton : allButtons) {
        if (anyButton->isChecked()) {
            const auto &button = anyButton->objectName();
            if (buttonToChoice.contains(button)) {
                settings->setValue("accelerator", buttonToChoice.value(button));
            } else if (buttonToPrecision.contains(button)) {
                settings->setValue("intelprec", buttonToPrecision.value(button));
            }
        }
    }

    QLineEdit *field;

#if defined(_OPENMP)
    // store number of threads, reset to 1 for "None" and "Opt" settings
    auto *mainwidget = dynamic_cast<LammpsGui *>(get_main_widget());
    field            = tabWidget->findChild<QLineEdit *>("nthreads");
    if (field && mainwidget) {
        int accel = settings->value("accelerator", AcceleratorTab::None).toInt();
        if ((accel == AcceleratorTab::None) || (accel == AcceleratorTab::Opt)) {
            mainwidget->nthreads = 1;
        } else if (field->hasAcceptableInput()) {
            settings->setValue("nthreads", field->text());
            mainwidget->nthreads = settings->value("nthreads", 1).toInt();
        }
    }
#endif

    // store setting for GPU package
    auto *box = tabWidget->findChild<QCheckBox *>("gpuneigh");
    if (box) settings->setValue("gpuneigh", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("gpupaironly");
    if (box) settings->setValue("gpupaironly", box->isChecked());

    // store image width, height, zoom, and rendering settings

    settings->beginGroup("snapshot");
    field = tabWidget->findChild<QLineEdit *>("xsize");
    if (field)
        if (field->hasAcceptableInput()) settings->setValue("xsize", field->text());
    field = tabWidget->findChild<QLineEdit *>("ysize");
    if (field)
        if (field->hasAcceptableInput()) settings->setValue("ysize", field->text());
    field = tabWidget->findChild<QLineEdit *>("zoom");
    if (field)
        if (field->hasAcceptableInput()) settings->setValue("zoom", field->text());
    box = tabWidget->findChild<QCheckBox *>("anti");
    if (box) settings->setValue("antialias", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("ssao");
    if (box) settings->setValue("ssao", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("shiny");
    if (box) settings->setValue("shinystyle", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("box");
    if (box) settings->setValue("box", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("axes");
    if (box) settings->setValue("axes", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("vdwstyle");
    if (box) settings->setValue("vdwstyle", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("autobond");
    if (box) settings->setValue("autobond", box->isChecked());
    auto *combo = tabWidget->findChild<QComboBox *>("background");
    if (combo) settings->setValue("background", combo->currentText());
    combo = tabWidget->findChild<QComboBox *>("boxcolor");
    if (combo) settings->setValue("boxcolor", combo->currentText());
    settings->endGroup();

    // general settings
    box = tabWidget->findChild<QCheckBox *>("echo");
    if (box) settings->setValue("echo", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("cite");
    if (box) settings->setValue("cite", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("logreplace");
    if (box) settings->setValue("logreplace", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("chartreplace");
    if (box) settings->setValue("chartreplace", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("imagereplace");
    if (box) settings->setValue("imagereplace", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("viewlog");
    if (box) settings->setValue("viewlog", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("viewchart");
    if (box) settings->setValue("viewchart", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("viewslide");
    if (box) settings->setValue("viewslide", box->isChecked());

    settings->beginGroup("tutorial");
    box = tabWidget->findChild<QCheckBox *>("solution");
    if (box) settings->setValue("solution", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("webpage");
    if (box) settings->setValue("webpage", box->isChecked());
    settings->endGroup();

    auto *spin = tabWidget->findChild<QSpinBox *>("updfreq");
    if (spin) settings->setValue("updfreq", spin->value());
    spin = tabWidget->findChild<QSpinBox *>("updchart");
    if (spin) settings->setValue("updchart", spin->value());

    field = tabWidget->findChild<QLineEdit *>("proxyval");
    if (field) settings->setValue("https_proxy", field->text());

    if (need_relaunch) {
        QMessageBox msg(QMessageBox::Information, QString("Relaunching LAMMPS-GUI"),
                        QString("LAMMPS library plugin path was changed.\n"
                                "LAMMPS-GUI must be relaunched."),
                        QMessageBox::Ok);
        msg.exec();
        const char *path = mystrdup(QCoreApplication::applicationFilePath());
        const char *arg0 = mystrdup(QCoreApplication::arguments().at(0));
        execl(path, arg0, (char *)nullptr);
    }

    // reformatting settings

    settings->beginGroup("reformat");
    spin = tabWidget->findChild<QSpinBox *>("cmdval");
    if (spin) settings->setValue("command", spin->value());
    spin = tabWidget->findChild<QSpinBox *>("typeval");
    if (spin) settings->setValue("type", spin->value());
    spin = tabWidget->findChild<QSpinBox *>("idval");
    if (spin) settings->setValue("id", spin->value());
    spin = tabWidget->findChild<QSpinBox *>("nameval");
    if (spin) settings->setValue("name", spin->value());
    box = tabWidget->findChild<QCheckBox *>("retval");
    if (box) settings->setValue("return", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("autoval");
    if (box) settings->setValue("automatic", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("savval");
    if (box) settings->setValue("autosave", box->isChecked());
    settings->endGroup();

    // chart window settings

    settings->beginGroup("charts");
    field = tabWidget->findChild<QLineEdit *>("title");
    if (field) settings->setValue("title", field->text());
    combo = tabWidget->findChild<QComboBox *>("smoothchoice");
    if (combo) settings->setValue("smoothchoice", combo->currentIndex());
    combo = tabWidget->findChild<QComboBox *>("rawbrush");
    if (combo) settings->setValue("rawbrush", combo->currentIndex());
    combo = tabWidget->findChild<QComboBox *>("smoothbrush");
    if (combo) settings->setValue("smoothbrush", combo->currentIndex());
    spin = tabWidget->findChild<QSpinBox *>("smoothwindow");
    if (spin) settings->setValue("smoothwindow", spin->value());
    spin = tabWidget->findChild<QSpinBox *>("smoothorder");
    if (spin) settings->setValue("smoothorder", spin->value());
    settings->endGroup();
    spin = tabWidget->findChild<QSpinBox *>("chartx");
    if (spin) settings->setValue("chartx", spin->value());
    spin = tabWidget->findChild<QSpinBox *>("charty");
    if (spin) settings->setValue("charty", spin->value());

    QDialog::accept();
}

GeneralTab::GeneralTab(QSettings *_settings, LammpsWrapper *_lammps, QWidget *parent) :
    QWidget(parent), settings(_settings), lammps(_lammps)
{
    auto *layout = new QGridLayout;

    auto *echo = new QCheckBox("Echo input to output buffer");
    echo->setObjectName("echo");
    echo->setCheckState(settings->value("echo", false).toBool() ? Qt::Checked : Qt::Unchecked);
    auto *cite = new QCheckBox("Include citation details");
    cite->setObjectName("cite");
    cite->setCheckState(settings->value("cite", false).toBool() ? Qt::Checked : Qt::Unchecked);
    auto *logv = new QCheckBox("Show Output window by default");
    logv->setObjectName("viewlog");
    logv->setCheckState(settings->value("viewlog", true).toBool() ? Qt::Checked : Qt::Unchecked);
    auto *pltv = new QCheckBox("Show Charts window by default");
    pltv->setObjectName("viewchart");
    pltv->setCheckState(settings->value("viewchart", true).toBool() ? Qt::Checked : Qt::Unchecked);
    auto *sldv = new QCheckBox("Show Slide Show window by default");
    sldv->setObjectName("viewslide");
    sldv->setCheckState(settings->value("viewslide", true).toBool() ? Qt::Checked : Qt::Unchecked);
    auto *logr = new QCheckBox("Replace Output window on new run");
    logr->setObjectName("logreplace");
    logr->setCheckState(settings->value("logreplace", true).toBool() ? Qt::Checked : Qt::Unchecked);
    auto *imgr = new QCheckBox("Replace Image window on new render");
    imgr->setObjectName("imagereplace");
    imgr->setCheckState(settings->value("imagereplace", true).toBool() ? Qt::Checked
                                                                       : Qt::Unchecked);
    auto *pltr = new QCheckBox("Replace Charts window on new run");
    pltr->setObjectName("chartreplace");
    pltr->setCheckState(settings->value("chartreplace", true).toBool() ? Qt::Checked
                                                                       : Qt::Unchecked);

    settings->beginGroup("tutorial");
    auto *solution = new QCheckBox("Download tutorial solutions enabled");
    solution->setObjectName("solution");
    solution->setCheckState(settings->value("solution", false).toBool() ? Qt::Checked
                                                                        : Qt::Unchecked);
    auto *webpage = new QCheckBox("Open tutorial webpage enabled");
    webpage->setObjectName("webpage");
    webpage->setCheckState(settings->value("webpage", true).toBool() ? Qt::Checked : Qt::Unchecked);
    settings->endGroup();

    auto *getallfont =
        new QPushButton(QIcon(":/icons/preferences-desktop-font.png"), "Select Default Font...");
    auto *gettextfont =
        new QPushButton(QIcon(":/icons/preferences-desktop-font.png"), "Select Text Font...");
    connect(getallfont, &QPushButton::released, this, &GeneralTab::newallfont);
    connect(gettextfont, &QPushButton::released, this, &GeneralTab::newtextfont);

    auto *freqlabel = new QLabel("Data update interval (ms):");
    auto *freqval   = new QSpinBox;
    freqval->setRange(1, 1000);
    freqval->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
    freqval->setValue(settings->value("updfreq", "10").toInt());
    freqval->setObjectName("updfreq");

    auto *chartlabel = new QLabel("Charts update interval (ms):");
    auto *chartval   = new QSpinBox;
    chartval->setRange(1, 5000);
    chartval->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
    chartval->setValue(settings->value("updchart", "500").toInt());
    chartval->setObjectName("updchart");

    int nrow = 0;
    layout->addWidget(new QHline, nrow++, 0, 1, 2);
    layout->addWidget(echo, nrow, 0);
    layout->addWidget(cite, nrow++, 1);
    layout->addWidget(new QHline, nrow++, 0, 1, 2);
    layout->addWidget(logv, nrow, 0);
    layout->addWidget(logr, nrow++, 1);
    layout->addWidget(pltv, nrow, 0);
    layout->addWidget(pltr, nrow++, 1);
    layout->addWidget(sldv, nrow, 0);
    layout->addWidget(imgr, nrow++, 1);
    layout->addWidget(new QHline, nrow++, 0, 1, 2);
    layout->addWidget(solution, nrow, 0);
    layout->addWidget(webpage, nrow++, 1);
    layout->addWidget(new QHline, nrow++, 0, 1, 2);
    layout->addWidget(getallfont, nrow, 0);
    layout->addWidget(gettextfont, nrow++, 1);
    layout->addWidget(new QHline, nrow++, 0, 1, 2);
    layout->addWidget(freqlabel, nrow, 0);
    layout->addWidget(freqval, nrow++, 1);
    layout->addWidget(chartlabel, nrow, 0);
    layout->addWidget(chartval, nrow++, 1);
    layout->addWidget(new QHline, nrow++, 0, 1, 2);

    auto *proxylabel = new QLabel("HTTPS proxy setting (empty for no proxy):");
    layout->addWidget(proxylabel, nrow, 0);

    auto https_proxy = QString::fromLocal8Bit(qgetenv("https_proxy"));
    if (https_proxy.isEmpty()) {
        https_proxy     = settings->value("https_proxy", "").toString();
        auto *proxyedit = new QLineEdit(https_proxy);
        proxyedit->setObjectName("proxyval");
        layout->addWidget(proxyedit, nrow++, 1);
    } else {
        layout->addWidget(new QLabel(https_proxy), nrow++, 1);
    }

#if defined(LAMMPS_GUI_USE_PLUGIN)
    layout->addWidget(new QHline, nrow++, 0, 1, 2);
    auto *pluginlabel = new QLabel("Path to LAMMPS Shared Library File:");
    auto *pluginedit =
        new QLineEdit(settings->value("plugin_path", "liblammpsplugin.so").toString());
    auto *pluginbrowse = new QPushButton("Browse...");
    auto *pluginlayout = new QHBoxLayout;
    pluginedit->setObjectName("pluginedit");
    pluginlayout->addWidget(pluginedit);
    pluginlayout->addWidget(pluginbrowse);

    connect(pluginbrowse, &QPushButton::released, this, &GeneralTab::pluginpath);

    layout->addWidget(pluginlabel, nrow++, 0, 1, 2);
    layout->addLayout(pluginlayout, nrow++, 0, 1, 2);
#endif
    layout->addWidget(new QHline, nrow++, 0, 1, 2);

    layout->addItem(new QSpacerItem(10, 10, QSizePolicy::Minimum, QSizePolicy::Expanding), nrow, 0);
    layout->addItem(new QSpacerItem(10, 10, QSizePolicy::Minimum, QSizePolicy::Expanding), nrow++,
                    1);
    setLayout(layout);
}

void GeneralTab::updatefonts(const QFont &all, const QFont &text)
{
    auto *mainwidget = dynamic_cast<LammpsGui *>(get_main_widget());
    if (mainwidget) {
        mainwidget->setFont(all);
        mainwidget->ui->textEdit->document()->setDefaultFont(text);
        if (mainwidget->wizard) mainwidget->wizard->setFont(all);
    }

    Preferences *prefs = nullptr;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "preferences") prefs = dynamic_cast<Preferences *>(widget);
    if (prefs) prefs->setFont(all);
}

void GeneralTab::newallfont()
{
    QSettings settings;
    QFont all, text;
    all.fromString(settings.value("allfont", QFont("Arial", -1).toString()).toString());
    text.fromString(settings.value("textfont", QFont("Monospace", -1).toString()).toString());

    bool ok    = false;
    QFont font = QFontDialog::getFont(&ok, all, this, QString("Select Default Font"));
    if (ok) updatefonts(font, text);

    settings.setValue("allfont", font.toString());
}

void GeneralTab::newtextfont()
{
    QSettings settings;
    QFont all, text;
    all.fromString(settings.value("allfont", QFont("Arial", -1).toString()).toString());
    text.fromString(settings.value("textfont", QFont("Monospace", -1).toString()).toString());

    bool ok    = false;
    QFont font = QFontDialog::getFont(&ok, text, this, QString("Select Text Font"));
    if (ok) updatefonts(all, font);

    settings.setValue("textfont", font.toString());
}

void GeneralTab::pluginpath()
{
    auto *field = findChild<QLineEdit *>("pluginedit");
    QString pluginfile =
        QFileDialog::getOpenFileName(this, "Select Shared LAMMPS Library to Load", field->text(),
                                     "Shared Objects (*.so *.dll *.dylib)");
    if (!pluginfile.isEmpty() && pluginfile.contains("liblammps", Qt::CaseSensitive)) {
        auto canonical = QFileInfo(pluginfile).canonicalFilePath();
        field->setText(pluginfile);
        settings->setValue("plugin_path", canonical);
        // ugly hack
        qobject_cast<Preferences *>(parent()->parent()->parent())->set_relaunch(true);
    }
}

AcceleratorTab::AcceleratorTab(QSettings *_settings, LammpsWrapper *_lammps, QWidget *parent) :
    QWidget(parent), settings(_settings), lammps(_lammps)
{
    auto *mainLayout  = new QHBoxLayout;
    auto *accelerator = new QGroupBox("Choose Accelerator:");
    auto *none        = new QRadioButton("&None");
    auto *opt         = new QRadioButton("O&pt");
    auto *openmp      = new QRadioButton("&OpenMP");
    auto *intel       = new QRadioButton("&Intel");
    auto *kokkos      = new QRadioButton("&Kokkos");
    auto *gpu         = new QRadioButton("GP&U");

    auto *accelframe   = new QFrame;
    auto *accelLayout  = new QVBoxLayout;
    auto *buttonLayout = new QVBoxLayout;
    accelLayout->addWidget(accelerator);
    buttonLayout->addWidget(none);
    buttonLayout->addWidget(opt);
    buttonLayout->addWidget(openmp);
    buttonLayout->addWidget(intel);
    buttonLayout->addWidget(kokkos);
    buttonLayout->addWidget(gpu);
    buttonLayout->addStretch(1);
    accelerator->setLayout(buttonLayout);
    accelframe->setLayout(accelLayout);
    mainLayout->addWidget(accelframe);

    none->setEnabled(true);
    none->setObjectName("none");
    opt->setEnabled(lammps->config_has_package("OPT"));
    opt->setObjectName("opt");
    openmp->setEnabled(lammps->config_has_package("OPENMP"));
    openmp->setObjectName("openmp");
    intel->setEnabled(lammps->config_has_package("INTEL"));
    intel->setObjectName("intel");
    // Kokkos support only works with OpenMP for now.
    kokkos->setEnabled(false);
    if (lammps->config_has_package("KOKKOS")) {
        if (lammps->config_accelerator("KOKKOS", "api", "openmp") &&
            !(lammps->config_accelerator("KOKKOS", "api", "cuda") ||
              lammps->config_accelerator("KOKKOS", "api", "hip") ||
              lammps->config_accelerator("KOKKOS", "api", "sycl")))
            kokkos->setEnabled(true);
    }
    kokkos->setObjectName("kokkos");
    gpu->setEnabled(lammps->config_has_package("GPU") && lammps->has_gpu_device());
    gpu->setObjectName("gpu");

    auto *choices      = new QFrame;
    auto *choiceLayout = new QVBoxLayout;
#if defined(_OPENMP)
    // maximum number of threads is limited half of available threads and no more than 16
    // unless OMP_NUM_THREADS is set to a larger value
    int maxthreads = std::min(QThread::idealThreadCount() / 2, 16);
    maxthreads     = std::max(maxthreads, 1);
    maxthreads     = std::max(maxthreads, qEnvironmentVariable("OMP_NUM_THREADS").toInt());

    auto *ntlabel  = new QLabel(QString("Number of threads (max %1):").arg(maxthreads));
    auto *ntchoice = new QLineEdit(settings->value("nthreads", maxthreads).toString());
    auto *intval   = new QIntValidator(1, maxthreads, this);
    ntchoice->setValidator(intval);
#else
    auto *ntlabel  = new QLabel("Number of threads (OpenMP not available):");
    auto *ntchoice = new QLineEdit("1");
    ntchoice->setEnabled(false);
#endif
    ntchoice->setObjectName("nthreads");

    connect(none, &QRadioButton::released, this, &AcceleratorTab::update_accel);
    connect(opt, &QRadioButton::released, this, &AcceleratorTab::update_accel);
    connect(openmp, &QRadioButton::released, this, &AcceleratorTab::update_accel);
    connect(intel, &QRadioButton::released, this, &AcceleratorTab::update_accel);
    connect(kokkos, &QRadioButton::released, this, &AcceleratorTab::update_accel);
    connect(gpu, &QRadioButton::released, this, &AcceleratorTab::update_accel);

    auto *intelLayout = new QHBoxLayout;
    auto *intelprec   = new QGroupBox("Intel Precision:");
    auto *inteldouble = new QRadioButton("&Double");
    auto *intelmixed  = new QRadioButton("&Mixed");
    auto *intelsingle = new QRadioButton("&Single");
    intelLayout->addWidget(inteldouble);
    inteldouble->setObjectName("inteldouble");
    intelLayout->addWidget(intelmixed);
    intelmixed->setObjectName("intelmixed");
    intelLayout->addWidget(intelsingle);
    intelsingle->setObjectName("intelsingle");
    intelprec->setLayout(intelLayout);
    intelprec->setObjectName("intelprec");
    intelprec->setEnabled(false);

    connect(inteldouble, &QRadioButton::released, this, &AcceleratorTab::update_accel);
    connect(intelmixed, &QRadioButton::released, this, &AcceleratorTab::update_accel);
    connect(intelsingle, &QRadioButton::released, this, &AcceleratorTab::update_accel);

    auto *gpuLayout   = new QHBoxLayout;
    auto *gpuchoice   = new QGroupBox("GPU Settings:");
    auto *gpuneigh    = new QCheckBox("Neighbor&list on GPU");
    auto *gpupaironly = new QCheckBox("Pair st&yles only");
    gpuLayout->addWidget(gpuneigh);
    gpuneigh->setObjectName("gpuneigh");
    gpuneigh->setCheckState(settings->value("gpuneigh", true).toBool() ? Qt::Checked
                                                                       : Qt::Unchecked);
    gpuLayout->addWidget(gpupaironly);
    gpupaironly->setObjectName("gpupaironly");
    gpupaironly->setCheckState(settings->value("gpupaironly", false).toBool() ? Qt::Checked
                                                                              : Qt::Unchecked);
    gpuchoice->setLayout(gpuLayout);
    gpuchoice->setObjectName("gpuchoice");
    gpuchoice->setEnabled(false);

    choiceLayout->addWidget(new QLabel("Settings for accelerator packages:\n"));
    choiceLayout->addWidget(ntlabel);
    choiceLayout->addWidget(ntchoice);
    choiceLayout->addWidget(intelprec);
    choiceLayout->addWidget(gpuchoice);
    choiceLayout->addStretch(1);
    choices->setLayout(choiceLayout);
    mainLayout->addWidget(choices);
    setLayout(mainLayout);

    // trigger update of nthreads line editor field depending on accelerator choice
    // fall back on None, if configured accelerator package is no longer available
    int choice = settings->value("accelerator", AcceleratorTab::None).toInt();
    int iprec  = settings->value("intelprec", AcceleratorTab::Mixed).toInt();
    if (iprec == AcceleratorTab::Double)
        inteldouble->setChecked(true);
    else if (iprec == AcceleratorTab::Mixed)
        intelmixed->setChecked(true);
    else if (iprec == AcceleratorTab::Single)
        intelsingle->setChecked(true);

    switch (choice) {
        case AcceleratorTab::Opt:
            if (opt->isEnabled())
                opt->click();
            else
                none->click();
            break;
        case AcceleratorTab::OpenMP:
            if (openmp->isEnabled())
                openmp->click();
            else
                none->click();
            break;
        case AcceleratorTab::Intel:
            if (intel->isEnabled()) {
                intel->click();
                intelprec->setEnabled(true);
            } else {
                none->click();
            }
            break;
        case AcceleratorTab::Kokkos:
            if (kokkos->isEnabled())
                kokkos->click();
            else
                none->click();
            break;
        case AcceleratorTab::Gpu:
            if (gpu->isEnabled()) {
                gpu->click();
                gpuchoice->setEnabled(true);
            } else
                none->click();
            break;
        case AcceleratorTab::None: // fallthrough
        default:
            none->click();
            break;
    }
}

void AcceleratorTab::update_accel()
{
    // store selected accelerator
    int choice = AcceleratorTab::None;

    QList<QRadioButton *> allButtons = findChildren<QRadioButton *>();
    for (auto &anyButton : allButtons) {
        if (anyButton->isChecked()) {
            const auto &button = anyButton->objectName();
            if (buttonToChoice.contains(button)) {
                choice = buttonToChoice.value(button);
            }
        }
    }

    auto *group = findChild<QGroupBox *>("intelprec");
    if (choice == AcceleratorTab::Intel) {
        group->setEnabled(true);
    } else {
        group->setEnabled(false);
    }

    group = findChild<QGroupBox *>("gpuchoice");
    if (choice == AcceleratorTab::Gpu) {
        group->setEnabled(true);
    } else {
        group->setEnabled(false);
    }

#if defined(_OPENMP)
    // The number of threads field is disabled and the value set to 1 for "None" and "Opt" choice
    auto *field = findChild<QLineEdit *>("nthreads");
    if (field) {
        if ((choice == AcceleratorTab::None) || (choice == AcceleratorTab::Opt)) {
            field->setText("1");
            field->setEnabled(false);
        } else {
            field->setText(settings->value("nthreads", 1).toString());
            field->setEnabled(true);
        }
    }
#endif
}

SnapshotTab::SnapshotTab(QSettings *_settings, QWidget *parent) :
    QWidget(parent), settings(_settings)
{
    auto *grid = new QGridLayout;

    auto *xsize = new QLabel("Image width:");
    auto *ysize = new QLabel("Image height:");
    auto *zoom  = new QLabel("Zoom factor:");
    auto *anti  = new QLabel("Antialias:");
    auto *ssao  = new QLabel("HQ Image mode:");
    auto *shiny = new QLabel("Shiny Image mode:");
    auto *bbox  = new QLabel("Show Box:");
    auto *axes  = new QLabel("Show Axes:");
    auto *vdw   = new QLabel("VDW Style:");
    auto *bond  = new QLabel("Dynamic Bonds:");
    auto *bclbl  = new QLabel("Bond Cutoff:");
    auto *cback = new QLabel("Background Color:");
    auto *cbox  = new QLabel("Box Color:");
    settings->beginGroup("snapshot");
    auto *xval = new QLineEdit(settings->value("xsize", "600").toString());
    auto *yval = new QLineEdit(settings->value("ysize", "600").toString());
    auto *zval = new QLineEdit(settings->value("zoom", "1.0").toString());
    auto *aval = new QCheckBox;
    auto *sval = new QCheckBox;
    auto *hval = new QCheckBox;
    auto *bval = new QCheckBox;
    auto *eval = new QCheckBox;
    auto *vval = new QCheckBox;
    auto *uval = new QCheckBox;
    auto *bcut = new QLineEdit(settings->value("bondcut", "1.6").toString());

    sval->setCheckState(settings->value("ssao", false).toBool() ? Qt::Checked : Qt::Unchecked);
    sval->setObjectName("ssao");
    aval->setCheckState(settings->value("antialias", false).toBool() ? Qt::Checked : Qt::Unchecked);
    aval->setObjectName("anti");
    hval->setCheckState(settings->value("shinystyle", true).toBool() ? Qt::Checked : Qt::Unchecked);
    hval->setObjectName("shiny");
    bval->setCheckState(settings->value("box", true).toBool() ? Qt::Checked : Qt::Unchecked);
    bval->setObjectName("box");
    eval->setCheckState(settings->value("axes", false).toBool() ? Qt::Checked : Qt::Unchecked);
    eval->setObjectName("axes");
    vval->setCheckState(settings->value("vdwstyle", false).toBool() ? Qt::Checked : Qt::Unchecked);
    vval->setObjectName("vdwstyle");
    uval->setCheckState(settings->value("autobond", false).toBool() ? Qt::Checked : Qt::Unchecked);
    uval->setObjectName("autobond");
    bcut->setObjectName("bondcut");

    auto *intval = new QIntValidator(100, 100000, this);
    xval->setValidator(intval);
    xval->setObjectName("xsize");
    yval->setValidator(intval);
    yval->setObjectName("ysize");
    zval->setValidator(new QDoubleValidator(0.01, 100.0, 100, this));
    zval->setObjectName("zoom");

    auto *background = new QComboBox;
    background->setObjectName("background");
    background->addItem("black");
    background->addItem("darkgray");
    background->addItem("gray");
    background->addItem("silver");
    background->addItem("white");
    background->setCurrentText(settings->value("background", "black").toString());

    auto *boxcolor = new QComboBox;
    boxcolor->setObjectName("boxcolor");
    boxcolor->addItem("yellow");
    boxcolor->addItem("silver");
    boxcolor->addItem("gray");
    boxcolor->addItem("darkred");
    boxcolor->addItem("darkgreen");
    boxcolor->addItem("darkblue");
    boxcolor->setCurrentText(settings->value("boxcolor", "yellow").toString());
    settings->endGroup();

    int i = 0;
    grid->addWidget(xsize, i, 0, Qt::AlignTop);
    grid->addWidget(xval, i++, 1, Qt::AlignTop);
    grid->addWidget(ysize, i, 0, Qt::AlignTop);
    grid->addWidget(yval, i++, 1, Qt::AlignTop);
    grid->addWidget(zoom, i, 0, Qt::AlignTop);
    grid->addWidget(zval, i++, 1, Qt::AlignTop);
    grid->addWidget(anti, i, 0, Qt::AlignTop);
    grid->addWidget(aval, i++, 1, Qt::AlignTop);
    grid->addWidget(ssao, i, 0, Qt::AlignTop);
    grid->addWidget(sval, i++, 1, Qt::AlignVCenter);
    grid->addWidget(shiny, i, 0, Qt::AlignTop);
    grid->addWidget(hval, i++, 1, Qt::AlignVCenter);
    grid->addWidget(bbox, i, 0, Qt::AlignTop);
    grid->addWidget(bval, i++, 1, Qt::AlignVCenter);
    grid->addWidget(axes, i, 0, Qt::AlignTop);
    grid->addWidget(eval, i++, 1, Qt::AlignVCenter);
    grid->addWidget(vdw, i, 0, Qt::AlignTop);
    grid->addWidget(vval, i++, 1, Qt::AlignVCenter);
    grid->addWidget(bond, i, 0, Qt::AlignTop);
    grid->addWidget(uval, i++, 1, Qt::AlignVCenter);
    grid->addWidget(bclbl, i, 0, Qt::AlignTop);
    grid->addWidget(bcut, i++, 1, Qt::AlignVCenter);
    grid->addWidget(cback, i, 0, Qt::AlignTop);
    grid->addWidget(background, i++, 1, Qt::AlignVCenter);
    grid->addWidget(cbox, i, 0, Qt::AlignTop);
    grid->addWidget(boxcolor, i++, 1, Qt::AlignVCenter);

    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), i, 0);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), i, 1);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Expanding, QSizePolicy::Expanding), i, 2);
    setLayout(grid);

    connect(vval, &QCheckBox::toggled, this, &SnapshotTab::choose_vdw);
    connect(uval, &QCheckBox::toggled, this, &SnapshotTab::choose_bond);
}

void SnapshotTab::choose_vdw()
{
    auto *vdw = findChild<QCheckBox *>("vdwstyle");
    auto *bnd = findChild<QCheckBox *>("autobond");
    if (vdw && bnd) {
        if (vdw->isChecked()) bnd->setCheckState(Qt::Unchecked);
    }
}

void SnapshotTab::choose_bond()
{
    auto *vdw = findChild<QCheckBox *>("vdwstyle");
    auto *bnd = findChild<QCheckBox *>("autobond");
    if (vdw && bnd) {
        if (bnd->isChecked()) vdw->setCheckState(Qt::Unchecked);
    }
}

EditorTab::EditorTab(QSettings *_settings, QWidget *parent) : QWidget(parent), settings(_settings)
{
    settings->beginGroup("reformat");
    auto *grid     = new QGridLayout;
    auto *reformat = new QLabel("Tab Reformatting settings:");
    auto *cmdlbl   = new QLabel("Command width:");
    auto *typelbl  = new QLabel("Type width:");
    auto *idlbl    = new QLabel("ID width:");
    auto *namelbl  = new QLabel("Name width:");
    auto *retlbl   = new QLabel("Reformat with 'Enter':");
    auto *autolbl  = new QLabel("Automatic completion:");
    auto *savlbl   = new QLabel("Auto-save on 'Run' and 'Quit':");
    auto *cmdval   = new QSpinBox;
    auto *typeval  = new QSpinBox;
    auto *idval    = new QSpinBox;
    auto *nameval  = new QSpinBox;
    auto *retval   = new QCheckBox;
    auto *autoval  = new QCheckBox;
    auto *savval   = new QCheckBox;
    cmdval->setObjectName("cmdval");
    cmdval->setRange(1, 32);
    cmdval->setValue(settings->value("command", "16").toInt());
    typeval->setObjectName("typeval");
    typeval->setRange(1, 32);
    typeval->setValue(settings->value("type", "4").toInt());
    idval->setObjectName("idval");
    idval->setRange(1, 32);
    idval->setValue(settings->value("id", "8").toInt());
    nameval->setObjectName("nameval");
    nameval->setRange(1, 32);
    nameval->setValue(settings->value("name", "8").toInt());
    retval->setObjectName("retval");
    retval->setCheckState(settings->value("return", false).toBool() ? Qt::Checked : Qt::Unchecked);
    autoval->setObjectName("autoval");
    autoval->setCheckState(settings->value("automatic", true).toBool() ? Qt::Checked
                                                                       : Qt::Unchecked);
    savval->setObjectName("savval");
    savval->setCheckState(settings->value("autosave", false).toBool() ? Qt::Checked
                                                                      : Qt::Unchecked);
    settings->endGroup();

    int i = 0;
    grid->addWidget(reformat, i++, 0, 1, 2, Qt::AlignTop | Qt::AlignHCenter);
    grid->addWidget(cmdlbl, i, 0, Qt::AlignTop);
    grid->addWidget(cmdval, i++, 1, Qt::AlignTop);
    grid->addWidget(typelbl, i, 0, Qt::AlignTop);
    grid->addWidget(typeval, i++, 1, Qt::AlignTop);
    grid->addWidget(idlbl, i, 0, Qt::AlignTop);
    grid->addWidget(idval, i++, 1, Qt::AlignTop);
    grid->addWidget(namelbl, i, 0, Qt::AlignTop);
    grid->addWidget(nameval, i++, 1, Qt::AlignTop);
    grid->addWidget(retlbl, i, 0, Qt::AlignTop);
    grid->addWidget(retval, i++, 1, Qt::AlignVCenter);
    grid->addWidget(autolbl, i, 0, Qt::AlignTop);
    grid->addWidget(autoval, i++, 1, Qt::AlignVCenter);
    grid->addWidget(new QLabel(" "), i++, 0);
    grid->addWidget(savlbl, i, 0, Qt::AlignTop);
    grid->addWidget(savval, i++, 1, Qt::AlignVCenter);

    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), i, 0);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), i, 1);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Expanding, QSizePolicy::Expanding), i, 2);
    setLayout(grid);
}

ChartsTab::ChartsTab(QSettings *_settings, QWidget *parent) : QWidget(parent), settings(_settings)
{
    auto *grid     = new QGridLayout;
    auto *chartlbl = new QLabel("Charts default settings:");

    settings->beginGroup("charts");
    auto *titlelbl = new QLabel("Default chart title:");
    auto *titletxt = new QLineEdit(settings->value("title", "Thermo: %f").toString());
    auto *titlehlp = new QLabel("(use %f for current input file)");

    // list of choices must be kepy in sync with list in chartviewer
    auto *smoothlbl = new QLabel("Default plot data choice:");
    auto *smoothval = new QComboBox;
    smoothval->addItem("Raw");
    smoothval->addItem("Smooth");
    smoothval->addItem("Both");
    smoothval->setObjectName("smoothchoice");
    smoothval->setCurrentIndex(settings->value("smoothchoice", 0).toInt());

    auto *rawbrlbl = new QLabel("Raw plot color:");
    auto *rawbrush = new QComboBox;
    rawbrush->addItem("Black");
    rawbrush->addItem("Blue");
    rawbrush->addItem("Red");
    rawbrush->addItem("Green");
    rawbrush->addItem("Gray");
    rawbrush->setObjectName("rawbrush");
    rawbrush->setCurrentIndex(settings->value("rawbrush", 1).toInt());

    auto *smoothbrlbl = new QLabel("Smooth plot color:");
    auto *smoothbrush = new QComboBox;
    smoothbrush->addItem("Black");
    smoothbrush->addItem("Blue");
    smoothbrush->addItem("Red");
    smoothbrush->addItem("Green");
    smoothbrush->addItem("Gray");
    smoothbrush->setObjectName("smoothbrush");
    smoothbrush->setCurrentIndex(settings->value("smoothbrush", 2).toInt());

    auto *smwindlbl = new QLabel("Default smoothing window:");
    auto *smwindval = new QSpinBox;
    smwindval->setRange(5, 999);
    smwindval->setValue(settings->value("smoothwindow", 10).toInt());
    smwindval->setObjectName("smoothwindow");

    auto *smordrlbl = new QLabel("Default smoothing order:");
    auto *smordrval = new QSpinBox;
    smordrval->setRange(1, 20);
    smordrval->setValue(settings->value("smoothorder", 4).toInt());
    smordrval->setObjectName("smoothorder");
    settings->endGroup();

    auto *chartxlbl = new QLabel("Chart default width:");
    auto *chartxval = new QSpinBox;
    chartxval->setRange(400, 40000);
    chartxval->setValue(settings->value("chartx", 500).toInt());
    chartxval->setObjectName("chartx");
    auto *chartylbl = new QLabel("Chart default height:");
    auto *chartyval = new QSpinBox;
    chartyval->setRange(300, 30000);
    chartyval->setValue(settings->value("charty", 320).toInt());
    chartyval->setObjectName("charty");

    int i = 0;
    grid->addWidget(chartlbl, i++, 0, 1, 2, Qt::AlignTop | Qt::AlignHCenter);
    grid->addWidget(titlelbl, i, 0, Qt::AlignTop);
    grid->addWidget(titletxt, i, 1, Qt::AlignTop);
    grid->addWidget(titlehlp, i++, 2, Qt::AlignTop);
    grid->addWidget(smoothlbl, i, 0, Qt::AlignTop);
    grid->addWidget(smoothval, i++, 1, Qt::AlignTop);
    grid->addWidget(rawbrlbl, i, 0, Qt::AlignTop);
    grid->addWidget(rawbrush, i++, 1, Qt::AlignTop);
    grid->addWidget(smoothbrlbl, i, 0, Qt::AlignTop);
    grid->addWidget(smoothbrush, i++, 1, Qt::AlignTop);
    grid->addWidget(smwindlbl, i, 0, Qt::AlignTop);
    grid->addWidget(smwindval, i++, 1, Qt::AlignTop);
    grid->addWidget(smordrlbl, i, 0, Qt::AlignTop);
    grid->addWidget(smordrval, i++, 1, Qt::AlignVCenter);
    grid->addWidget(chartxlbl, i, 0, Qt::AlignTop);
    grid->addWidget(chartxval, i++, 1, Qt::AlignVCenter);
    grid->addWidget(chartylbl, i, 0, Qt::AlignTop);
    grid->addWidget(chartyval, i++, 1, Qt::AlignVCenter);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), i, 0);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), i, 1);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Expanding, QSizePolicy::Expanding), i, 2);
    setLayout(grid);
}

// Local Variables:
// c-basic-offset: 4
// End:
