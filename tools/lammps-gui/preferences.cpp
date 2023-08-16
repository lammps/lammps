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
#include <QTabWidget>
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

// duplicate string
static char *mystrdup(const std::string &text)
{
    auto tmp = new char[text.size() + 1];
    memcpy(tmp, text.c_str(), text.size() + 1);
    return tmp;
}

Preferences::Preferences(LammpsWrapper *_lammps, QWidget *parent) :
    QDialog(parent), tabWidget(new QTabWidget),
    buttonBox(new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel)),
    settings(new QSettings), lammps(_lammps), need_relaunch(false)
{
    tabWidget->addTab(new GeneralTab(settings, lammps), "&General Settings");
    tabWidget->addTab(new AcceleratorTab(settings, lammps), "&Accelerators");
    tabWidget->addTab(new SnapshotTab(settings), "&Snapshot Image");

    connect(buttonBox, &QDialogButtonBox::accepted, this, &Preferences::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);

    auto *layout = new QVBoxLayout;
    layout->addWidget(tabWidget);
    layout->addWidget(buttonBox);
    setLayout(layout);
    setWindowTitle("LAMMPS-GUI - Preferences");
    resize(500, 400);
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

    // store number of threads
    QLineEdit *field = tabWidget->findChild<QLineEdit *>("nthreads");
    if (field)
        if (field->hasAcceptableInput()) settings->setValue("nthreads", field->text());

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

    if (need_relaunch) {
        QMessageBox msg(QMessageBox::Information, QString("Relaunching LAMMPS-GUI"),
                        QString("LAMMPS library plugin path was changed.\n"
                                "LAMMPS-GUI must be relaunched."),
                        QMessageBox::Ok);
        msg.exec();
        const char *path = mystrdup(QCoreApplication::applicationFilePath().toStdString());
        const char *arg0 = mystrdup(QCoreApplication::arguments().at(0).toStdString());
        execl(path, arg0, (char *)NULL);
    }
    QDialog::accept();
}

GeneralTab::GeneralTab(QSettings *_settings, LammpsWrapper *_lammps, QWidget *parent) :
    QWidget(parent), settings(_settings), lammps(_lammps)
{
    auto *layout = new QVBoxLayout;

    auto *echo = new QCheckBox("Echo input to log");
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
    auto *logr = new QCheckBox("Replace log window on new run");
    logr->setObjectName("logreplace");
    logr->setCheckState(settings->value("logreplace", false).toBool() ? Qt::Checked
                                                                      : Qt::Unchecked);
    auto *imgr = new QCheckBox("Replace image window on new render");
    imgr->setObjectName("imagereplace");
    imgr->setCheckState(settings->value("imagereplace", false).toBool() ? Qt::Checked
                                                                        : Qt::Unchecked);
    auto *pltr = new QCheckBox("Replace chart window on new run");
    pltr->setObjectName("chartreplace");
    pltr->setCheckState(settings->value("chartreplace", false).toBool() ? Qt::Checked
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
        new QPushButton(QIcon(":/preferences-desktop-font.png"), "Select Default Font...");
    auto *gettextfont =
        new QPushButton(QIcon(":/preferences-desktop-font.png"), "Select Text Font...");
    fontlayout->addWidget(getallfont);
    fontlayout->addWidget(gettextfont);
    connect(getallfont, &QPushButton::released, this, &GeneralTab::newallfont);
    connect(gettextfont, &QPushButton::released, this, &GeneralTab::newtextfont);

    layout->addWidget(echo);
    layout->addWidget(cite);
    layout->addWidget(logv);
    layout->addWidget(pltv);
    layout->addWidget(logr);
    layout->addWidget(pltr);
    layout->addWidget(imgr);
#if defined(LAMMPS_GUI_USE_PLUGIN)
    layout->addWidget(pluginlabel);
    layout->addLayout(pluginlayout);
#endif
    layout->addLayout(fontlayout);
    layout->addStretch(1);
    setLayout(layout);
}

void GeneralTab::updatefonts(const QFont &all, const QFont &text)
{
    LammpsGui *main;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") main = dynamic_cast<LammpsGui *>(widget);

    QApplication::setFont(all);
    main->ui->textEdit->document()->setDefaultFont(text);
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
    maxthreads = omp_get_max_threads();
#endif
    auto *choices      = new QFrame;
    auto *choiceLayout = new QVBoxLayout;
    auto *ntlabel      = new QLabel("Number of threads:");
    auto *ntchoice     = new QLineEdit(settings->value("nthreads", maxthreads).toString());
    auto *intval       = new QIntValidator(1, maxthreads, this);
    ntchoice->setValidator(intval);
    ntchoice->setObjectName("nthreads");

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
    sval->setCheckState(settings->value("ssao", false).toBool() ? Qt::Checked : Qt::Unchecked);
    sval->setObjectName("ssao");
    aval->setCheckState(settings->value("antialias", false).toBool() ? Qt::Checked : Qt::Unchecked);
    aval->setObjectName("anti");
    bval->setCheckState(settings->value("box", true).toBool() ? Qt::Checked : Qt::Unchecked);
    bval->setObjectName("box");
    eval->setCheckState(settings->value("axes", false).toBool() ? Qt::Checked : Qt::Unchecked);
    eval->setObjectName("axes");

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

    grid->addWidget(xsize, 0, 0, Qt::AlignTop);
    grid->addWidget(ysize, 1, 0, Qt::AlignTop);
    grid->addWidget(zoom, 2, 0, Qt::AlignTop);
    grid->addWidget(anti, 3, 0, Qt::AlignTop);
    grid->addWidget(ssao, 4, 0, Qt::AlignTop);
    grid->addWidget(bbox, 5, 0, Qt::AlignTop);
    grid->addWidget(axes, 6, 0, Qt::AlignTop);
    grid->addWidget(cback, 7, 0, Qt::AlignTop);
    grid->addWidget(cbox, 8, 0, Qt::AlignTop);
    grid->addWidget(xval, 0, 1, Qt::AlignTop);
    grid->addWidget(yval, 1, 1, Qt::AlignTop);
    grid->addWidget(zval, 2, 1, Qt::AlignTop);
    grid->addWidget(aval, 3, 1, Qt::AlignTop);
    grid->addWidget(sval, 4, 1, Qt::AlignVCenter);
    grid->addWidget(bval, 5, 1, Qt::AlignVCenter);
    grid->addWidget(eval, 6, 1, Qt::AlignVCenter);
    grid->addWidget(background, 7, 1, Qt::AlignVCenter);
    grid->addWidget(boxcolor, 8, 1, Qt::AlignVCenter);

    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), 9, 0);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), 9, 1);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Expanding, QSizePolicy::Expanding), 9, 2);
    setLayout(grid);
}

// Local Variables:
// c-basic-offset: 4
// End:
