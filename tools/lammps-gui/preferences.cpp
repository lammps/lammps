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
#include <QDir>
#include <QDoubleValidator>
#include <QFileDialog>
#include <QFontDialog>
#include <QGroupBox>
#include <QHBoxLayout>
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
#include <QThread>
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

Preferences::Preferences(LammpsWrapper *_lammps, QWidget *parent) :
    QDialog(parent), need_relaunch(false), tabWidget(new QTabWidget),
    buttonBox(new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel)),
    settings(new QSettings), lammps(_lammps)
{
    tabWidget->addTab(new GeneralTab(settings, lammps), "&General Settings");
    tabWidget->addTab(new AcceleratorTab(settings, lammps), "&Accelerators");
    tabWidget->addTab(new SnapshotTab(settings), "&Snapshot Image");
    tabWidget->addTab(new EditorTab(settings), "&Editor Settings");

    connect(buttonBox, &QDialogButtonBox::accepted, this, &Preferences::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);

    auto *layout = new QVBoxLayout;
    layout->addWidget(tabWidget);
    layout->addWidget(buttonBox);
    setLayout(layout);
    setWindowIcon(QIcon(":/icons/lammps-icon-128x128.png"));
    setWindowTitle("LAMMPS-GUI - Preferences");
    resize(600, 450);
}

Preferences::~Preferences()
{
    delete buttonBox;
    delete tabWidget;
    delete settings;
}

void Preferences::accept()
{
    // store all data in settings class
    // and then confirm accepting

    // store selected accelerator
    QList<QRadioButton *> allButtons = tabWidget->findChildren<QRadioButton *>();
    for (int i = 0; i < allButtons.size(); ++i) {
        if (allButtons[i]->isChecked()) {
            if (allButtons[i]->objectName() == "none")
                settings->setValue("accelerator", QString::number(AcceleratorTab::None));
            if (allButtons[i]->objectName() == "opt")
                settings->setValue("accelerator", QString::number(AcceleratorTab::Opt));
            if (allButtons[i]->objectName() == "openmp")
                settings->setValue("accelerator", QString::number(AcceleratorTab::OpenMP));
            if (allButtons[i]->objectName() == "intel")
                settings->setValue("accelerator", QString::number(AcceleratorTab::Intel));
            if (allButtons[i]->objectName() == "kokkos")
                settings->setValue("accelerator", QString::number(AcceleratorTab::Kokkos));
            if (allButtons[i]->objectName() == "gpu")
                settings->setValue("accelerator", QString::number(AcceleratorTab::Gpu));
        }
    }

    // store number of threads, reset to 1 for "None" and "Opt" settings
    QLineEdit *field = tabWidget->findChild<QLineEdit *>("nthreads");
    if (field) {
        int accel = settings->value("accelerator", AcceleratorTab::None).toInt();
        if ((accel == AcceleratorTab::None) || (accel == AcceleratorTab::Opt))
            settings->setValue("nthreads", 1);
        else if (field->hasAcceptableInput())
            settings->setValue("nthreads", field->text());
    }

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
    QCheckBox *box = tabWidget->findChild<QCheckBox *>("anti");
    if (box) settings->setValue("antialias", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("ssao");
    if (box) settings->setValue("ssao", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("box");
    if (box) settings->setValue("box", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("axes");
    if (box) settings->setValue("axes", box->isChecked());
    box = tabWidget->findChild<QCheckBox *>("vdwstyle");
    if (box) settings->setValue("vdwstyle", box->isChecked());
    QComboBox *combo = tabWidget->findChild<QComboBox *>("background");
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

    auto spin = tabWidget->findChild<QSpinBox *>("updfreq");
    if (spin) settings->setValue("updfreq", spin->value());

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
    settings->endGroup();

    QDialog::accept();
}

GeneralTab::GeneralTab(QSettings *_settings, LammpsWrapper *_lammps, QWidget *parent) :
    QWidget(parent), settings(_settings), lammps(_lammps)
{
    auto *layout = new QVBoxLayout;

    auto *echo = new QCheckBox("Echo input to output buffer");
    echo->setObjectName("echo");
    echo->setCheckState(settings->value("echo", false).toBool() ? Qt::Checked : Qt::Unchecked);
    auto *cite = new QCheckBox("Include citation details");
    cite->setObjectName("cite");
    cite->setCheckState(settings->value("cite", false).toBool() ? Qt::Checked : Qt::Unchecked);
    auto *logv = new QCheckBox("Show log window by default");
    logv->setObjectName("viewlog");
    logv->setCheckState(settings->value("viewlog", true).toBool() ? Qt::Checked : Qt::Unchecked);
    auto *pltv = new QCheckBox("Show chart window by default");
    pltv->setObjectName("viewchart");
    pltv->setCheckState(settings->value("viewchart", true).toBool() ? Qt::Checked : Qt::Unchecked);
    auto *sldv = new QCheckBox("Show slide show window by default");
    sldv->setObjectName("viewslide");
    sldv->setCheckState(settings->value("viewslide", true).toBool() ? Qt::Checked : Qt::Unchecked);
    auto *logr = new QCheckBox("Replace log window on new run");
    logr->setObjectName("logreplace");
    logr->setCheckState(settings->value("logreplace", true).toBool() ? Qt::Checked : Qt::Unchecked);
    auto *imgr = new QCheckBox("Replace image window on new render");
    imgr->setObjectName("imagereplace");
    imgr->setCheckState(settings->value("imagereplace", true).toBool() ? Qt::Checked
                                                                       : Qt::Unchecked);
    auto *pltr = new QCheckBox("Replace chart window on new run");
    pltr->setObjectName("chartreplace");
    pltr->setCheckState(settings->value("chartreplace", true).toBool() ? Qt::Checked
                                                                       : Qt::Unchecked);

#if defined(LAMMPS_GUI_USE_PLUGIN)
    auto *pluginlabel = new QLabel("Path to LAMMPS Shared Library File:");
    auto *pluginedit =
        new QLineEdit(settings->value("plugin_path", "liblammpsplugin.so").toString());
    auto *pluginbrowse = new QPushButton("Browse...");
    auto *pluginlayout = new QHBoxLayout;
    pluginedit->setObjectName("pluginedit");
    pluginlayout->addWidget(pluginedit);
    pluginlayout->addWidget(pluginbrowse);

    connect(pluginbrowse, &QPushButton::released, this, &GeneralTab::pluginpath);
#endif

    auto *fontlayout = new QHBoxLayout;
    auto *getallfont =
        new QPushButton(QIcon(":/icons/preferences-desktop-font.png"), "Select Default Font...");
    auto *gettextfont =
        new QPushButton(QIcon(":/icons/preferences-desktop-font.png"), "Select Text Font...");
    fontlayout->addWidget(getallfont);
    fontlayout->addWidget(gettextfont);
    connect(getallfont, &QPushButton::released, this, &GeneralTab::newallfont);
    connect(gettextfont, &QPushButton::released, this, &GeneralTab::newtextfont);

    auto *freqlayout = new QHBoxLayout;
    auto *freqlabel  = new QLabel("GUI update interval (ms)");
    auto *freqval    = new QSpinBox;
    freqval->setRange(1, 1000);
    freqval->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
    freqval->setValue(settings->value("updfreq", "100").toInt());
    freqval->setObjectName("updfreq");
    freqlayout->addWidget(freqlabel);
    freqlayout->addWidget(freqval);
    freqlayout->addStretch(1);

    layout->addWidget(echo);
    layout->addWidget(cite);
    layout->addWidget(logv);
    layout->addWidget(pltv);
    layout->addWidget(sldv);
    layout->addWidget(logr);
    layout->addWidget(pltr);
    layout->addWidget(imgr);
#if defined(LAMMPS_GUI_USE_PLUGIN)
    layout->addWidget(pluginlabel);
    layout->addLayout(pluginlayout);
#endif
    layout->addLayout(fontlayout);
    layout->addLayout(freqlayout);
    layout->addStretch(1);
    setLayout(layout);
}

void GeneralTab::updatefonts(const QFont &all, const QFont &text)
{
    LammpsGui *main = nullptr;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") main = dynamic_cast<LammpsGui *>(widget);

    QApplication::setFont(all);
    if (main) main->ui->textEdit->document()->setDefaultFont(text);
}

void GeneralTab::newallfont()
{
    QSettings settings;
    QFont all, text;
    all.fromString(settings.value("allfont", "").toString());
    text.fromString(settings.value("textfont", "").toString());

    bool ok    = false;
    QFont font = QFontDialog::getFont(&ok, all, this, QString("Select Default Font"));
    if (ok) updatefonts(font, text);

    settings.setValue("allfont", font.toString());
}

void GeneralTab::newtextfont()
{
    QSettings settings;
    QFont all, text;
    all.fromString(settings.value("allfont", "").toString());
    text.fromString(settings.value("textfont", "").toString());

    bool ok    = false;
    QFont font = QFontDialog::getFont(&ok, text, this, QString("Select Text Font"));
    if (ok) updatefonts(all, font);

    settings.setValue("textfont", font.toString());
}

void GeneralTab::pluginpath()
{
    QLineEdit *field = findChild<QLineEdit *>("pluginedit");
    QString pluginfile =
        QFileDialog::getOpenFileName(this, "Select Shared LAMMPS Library to Load", field->text(),
                                     "Shared Objects (*.so *.dll *.dylib)");
    if (!pluginfile.isEmpty() && pluginfile.contains("liblammps", Qt::CaseSensitive)) {
        auto canonical = QFileInfo(pluginfile).canonicalFilePath();
        field->setText(pluginfile);
        settings->setValue("plugin_path", canonical);
        // ugly hack
        qobject_cast<Preferences *>(parent()->parent()->parent())->need_relaunch = true;
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
    auto *gpu         = new QRadioButton("&GPU");

    auto *buttonLayout = new QVBoxLayout;
    buttonLayout->addWidget(none);
    buttonLayout->addWidget(opt);
    buttonLayout->addWidget(openmp);
    buttonLayout->addWidget(intel);
    buttonLayout->addWidget(kokkos);
    buttonLayout->addWidget(gpu);
    buttonLayout->addStretch(1);
    accelerator->setLayout(buttonLayout);
    mainLayout->addWidget(accelerator);

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

    int choice = settings->value("accelerator", AcceleratorTab::None).toInt();
    switch (choice) {
        case AcceleratorTab::Opt:
            if (opt->isEnabled()) opt->setChecked(true);
            break;
        case AcceleratorTab::OpenMP:
            if (openmp->isEnabled()) openmp->setChecked(true);
            break;
        case AcceleratorTab::Intel:
            if (intel->isEnabled()) intel->setChecked(true);
            break;
        case AcceleratorTab::Kokkos:
            if (kokkos->isEnabled()) kokkos->setChecked(true);
            break;
        case AcceleratorTab::Gpu:
            if (gpu->isEnabled()) gpu->setChecked(true);
            break;
        case AcceleratorTab::None: // fallthrough
        default:
            none->setChecked(true);
            break;
    }

    int maxthreads = 1;
#if defined(_OPENMP)
    maxthreads = QThread::idealThreadCount();
#endif
    auto *choices      = new QFrame;
    auto *choiceLayout = new QVBoxLayout;
#if defined(_OPENMP)
    auto *ntlabel      = new QLabel(QString("Number of threads (max %1):").arg(maxthreads));
    auto *ntchoice     = new QLineEdit(settings->value("nthreads", maxthreads).toString());
#else
    auto *ntlabel      = new QLabel(QString("Number of threads (OpenMP not available):"));
    auto *ntchoice     = new QLineEdit("1");
#endif
    auto *intval       = new QIntValidator(1, maxthreads, this);
    ntchoice->setValidator(intval);
    ntchoice->setObjectName("nthreads");
#if !defined(_OPENMP)
    ntchoice->setEnabled(false);
#endif

    choiceLayout->addWidget(ntlabel);
    choiceLayout->addWidget(ntchoice);
    choices->setLayout(choiceLayout);
    choiceLayout->addStretch(1);

    mainLayout->addWidget(choices);
    setLayout(mainLayout);
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
    auto *bbox  = new QLabel("Show Box:");
    auto *axes  = new QLabel("Show Axes:");
    auto *vdw   = new QLabel("VDW Style:");
    auto *cback = new QLabel("Background Color:");
    auto *cbox  = new QLabel("Box Color:");
    settings->beginGroup("snapshot");
    auto *xval = new QLineEdit(settings->value("xsize", "800").toString());
    auto *yval = new QLineEdit(settings->value("ysize", "600").toString());
    auto *zval = new QLineEdit(settings->value("zoom", "1.0").toString());
    auto *aval = new QCheckBox;
    auto *sval = new QCheckBox;
    auto *bval = new QCheckBox;
    auto *eval = new QCheckBox;
    auto *vval = new QCheckBox;
    sval->setCheckState(settings->value("ssao", false).toBool() ? Qt::Checked : Qt::Unchecked);
    sval->setObjectName("ssao");
    aval->setCheckState(settings->value("antialias", false).toBool() ? Qt::Checked : Qt::Unchecked);
    aval->setObjectName("anti");
    bval->setCheckState(settings->value("box", true).toBool() ? Qt::Checked : Qt::Unchecked);
    bval->setObjectName("box");
    eval->setCheckState(settings->value("axes", false).toBool() ? Qt::Checked : Qt::Unchecked);
    eval->setObjectName("axes");
    vval->setCheckState(settings->value("vdwstyle", false).toBool() ? Qt::Checked : Qt::Unchecked);
    vval->setObjectName("vdwstyle");

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
    grid->addWidget(bbox, i, 0, Qt::AlignTop);
    grid->addWidget(bval, i++, 1, Qt::AlignVCenter);
    grid->addWidget(axes, i, 0, Qt::AlignTop);
    grid->addWidget(eval, i++, 1, Qt::AlignVCenter);
    grid->addWidget(vdw, i, 0, Qt::AlignTop);
    grid->addWidget(vval, i++, 1, Qt::AlignVCenter);
    grid->addWidget(cback, i, 0, Qt::AlignTop);
    grid->addWidget(background, i++, 1, Qt::AlignVCenter);
    grid->addWidget(cbox, i, 0, Qt::AlignTop);
    grid->addWidget(boxcolor, i++, 1, Qt::AlignVCenter);

    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), i, 0);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), i, 1);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Expanding, QSizePolicy::Expanding), i, 2);
    setLayout(grid);
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
    auto *cmdval   = new QSpinBox;
    auto *typeval  = new QSpinBox;
    auto *idval    = new QSpinBox;
    auto *nameval  = new QSpinBox;
    auto *retval   = new QCheckBox;
    auto *autoval  = new QCheckBox;
    cmdval->setRange(1, 32);
    cmdval->setValue(settings->value("command", "16").toInt());
    cmdval->setObjectName("cmdval");
    typeval->setRange(1, 32);
    typeval->setValue(settings->value("type", "4").toInt());
    typeval->setObjectName("typeval");
    idval->setRange(1, 32);
    idval->setValue(settings->value("id", "8").toInt());
    idval->setObjectName("idval");
    nameval->setRange(1, 32);
    nameval->setValue(settings->value("name", "8").toInt());
    nameval->setObjectName("nameval");
    retval->setCheckState(settings->value("return", true).toBool() ? Qt::Checked : Qt::Unchecked);
    retval->setObjectName("retval");
    autoval->setCheckState(settings->value("automatic", true).toBool() ? Qt::Checked
                                                                       : Qt::Unchecked);
    autoval->setObjectName("autoval");
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

    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), i, 0);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), i, 1);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Expanding, QSizePolicy::Expanding), i, 2);
    setLayout(grid);
}

// Local Variables:
// c-basic-offset: 4
// End:
