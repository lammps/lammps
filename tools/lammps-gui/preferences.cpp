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

#include "lammpswrapper.h"

#include <QCheckBox>
#include <QCoreApplication>
#include <QDialogButtonBox>
#include <QDir>
#include <QDoubleValidator>
#include <QFileDialog>
#include <QGroupBox>
#include <QHBoxLayout>
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

Preferences::Preferences(LammpsWrapper *_lammps, QWidget *parent) :
    QDialog(parent), tabWidget(new QTabWidget),
    buttonBox(new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel)),
    settings(new QSettings), lammps(_lammps)
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

    // store image width, height, and zoom

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
    settings->endGroup();

    // general settings
    QCheckBox *box = tabWidget->findChild<QCheckBox *>("echo");
    if (box) settings->setValue("echo", box->isChecked() ? "1" : "0");
    box = tabWidget->findChild<QCheckBox *>("cite");
    if (box) settings->setValue("cite", box->isChecked() ? "1" : "0");

    QDialog::accept();
}

GeneralTab::GeneralTab(QSettings *_settings, LammpsWrapper *_lammps, QWidget *parent) :
    QWidget(parent), settings(_settings), lammps(_lammps)
{
    auto *layout = new QVBoxLayout;

    auto *echo = new QCheckBox("Echo input to log");
    echo->setCheckState(settings->value("echo", "0").toInt() ? Qt::Checked : Qt::Unchecked);
    echo->setObjectName("echo");
    auto *cite = new QCheckBox("Include Citations");
    cite->setCheckState(settings->value("cite", "0").toInt() ? Qt::Checked : Qt::Unchecked);
    cite->setObjectName("cite");
#if !defined(__APPLE__)
    auto *tmplabel  = new QLabel("Scratch Folder:");
    auto *tmpedit   = new QLineEdit(settings->value("tempdir", ".").toString());
    auto *tmpbrowse = new QPushButton("Browse...");
    auto *tmplayout = new QHBoxLayout;
    tmpedit->setObjectName("tmpedit");
    tmplayout->addWidget(tmplabel);
    tmplayout->addWidget(tmpedit);
    tmplayout->addWidget(tmpbrowse);
#endif

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
    layout->addWidget(echo);
    layout->addWidget(cite);
#if !defined(__APPLE__)
    layout->addLayout(tmplayout);
#endif
#if defined(LAMMPS_GUI_USE_PLUGIN)
    layout->addWidget(pluginlabel);
    layout->addLayout(pluginlayout);
#endif
    layout->addStretch(1);
    setLayout(layout);
}

void GeneralTab::newtmpfolder()
{
    QLineEdit *field = findChild<QLineEdit *>("tmpedit");
    QString tmpdir =
        QFileDialog::getExistingDirectory(this, "Find Folder for Temporary Files", field->text());
    if (!tmpdir.isEmpty()) field->setText(tmpdir);
}

void GeneralTab::pluginpath()
{
    QLineEdit *field = findChild<QLineEdit *>("pluginedit");
    QString pluginfile =
        QFileDialog::getOpenFileName(this, "Select Shared LAMMPS Library to Load", field->text(),
                                     "Shared Objects (*.so *.dll *.dylib)");
    if (!pluginfile.isEmpty() && pluginfile.contains("liblammps", Qt::CaseSensitive)) {
        if (lammps->load_lib(pluginfile.toStdString().c_str())) {
            auto canonical = QFileInfo(pluginfile).canonicalFilePath();
            field->setText(pluginfile);
            settings->setValue("plugin_path", canonical);
        } else {
            // plugin did not load cannot continue
            settings->remove("plugin_path");
            QMessageBox::critical(this, "Error", "Cannot open LAMMPS shared library file");
            QCoreApplication::quit();
        }
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
    settings->beginGroup("snapshot");
    auto *xval = new QLineEdit(settings->value("xsize", "800").toString());
    auto *yval = new QLineEdit(settings->value("ysize", "600").toString());
    auto *zval = new QLineEdit(settings->value("zoom", "1.0").toString());
    settings->endGroup();

    auto *intval = new QIntValidator(100, 100000, this);
    xval->setValidator(intval);
    xval->setObjectName("xsize");
    yval->setValidator(intval);
    yval->setObjectName("ysize");
    zval->setValidator(new QDoubleValidator(0.01, 100.0, 100, this));
    zval->setObjectName("zoom");

    grid->addWidget(xsize, 0, 0, Qt::AlignTop);
    grid->addWidget(ysize, 1, 0, Qt::AlignTop);
    grid->addWidget(zoom, 2, 0, Qt::AlignTop);
    grid->addWidget(xval, 0, 1, Qt::AlignTop);
    grid->addWidget(yval, 1, 1, Qt::AlignTop);
    grid->addWidget(zval, 2, 1, Qt::AlignTop);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), 3, 0);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Minimum, QSizePolicy::Expanding), 3, 1);
    grid->addItem(new QSpacerItem(100, 100, QSizePolicy::Expanding, QSizePolicy::Expanding), 3, 2);
    setLayout(grid);
}

// Local Variables:
// c-basic-offset: 4
// End:
