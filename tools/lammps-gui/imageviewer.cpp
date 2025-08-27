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

#include "imageviewer.h"

#include "helpers.h"
#include "lammpsgui.h"
#include "lammpswrapper.h"
#include "qaddon.h"

#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QDir>
#include <QDoubleValidator>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QFontMetrics>
#include <QGuiApplication>
#include <QHBoxLayout>
#include <QIcon>
#include <QImage>
#include <QImageReader>
#include <QIntValidator>
#include <QKeySequence>
#include <QLabel>
#include <QLineEdit>
#include <QMenu>
#include <QMenuBar>
#include <QPalette>
#include <QPixmap>
#include <QPushButton>
#include <QScrollArea>
#include <QScrollBar>
#include <QSettings>
#include <QSizePolicy>
#include <QSpinBox>
#include <QStringList>
#include <QVBoxLayout>
#include <QVariant>

#include <algorithm>
#include <cmath>
#include <unordered_set>

// clang-format off
/* periodic table of elements for translation of ordinal to atom type */
namespace {
    const char * const pte_label[] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P" , "S",  "Cl", "Ar", "K",  "Ca", "Sc",
    "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc",
    "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
    "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg"
};
constexpr int nr_pte_entries = sizeof(pte_label) / sizeof(char *);

/* corresponding table of masses. */
constexpr double pte_mass[] = {
    /* X  */ 0.00000, 1.00794, 4.00260, 6.941, 9.012182, 10.811,
    /* C  */ 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
    /* Na */ 22.989770, 24.3050, 26.981538, 28.0855, 30.973761,
    /* S  */ 32.065, 35.453, 39.948, 39.0983, 40.078, 44.955910,
    /* Ti */ 47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.9332,
    /* Ni */ 58.6934, 63.546, 65.409, 69.723, 72.64, 74.92160,
    /* Se */ 78.96, 79.904, 83.798, 85.4678, 87.62, 88.90585,
    /* Zr */ 91.224, 92.90638, 95.94, 98.0, 101.07, 102.90550,
    /* Pd */ 106.42, 107.8682, 112.411, 114.818, 118.710, 121.760,
    /* Te */ 127.60, 126.90447, 131.293, 132.90545, 137.327,
    /* La */ 138.9055, 140.116, 140.90765, 144.24, 145.0, 150.36,
    /* Eu */ 151.964, 157.25, 158.92534, 162.500, 164.93032,
    /* Er */ 167.259, 168.93421, 173.04, 174.967, 178.49, 180.9479,
    /* W  */ 183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,
    /* Hg */ 200.59, 204.3833, 207.2, 208.98038, 209.0, 210.0, 222.0,
    /* Fr */ 223.0, 226.0, 227.0, 232.0381, 231.03588, 238.02891,
    /* Np */ 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0,
    /* Md */ 258.0, 259.0, 262.0, 261.0, 262.0, 266.0, 264.0, 269.0,
    /* Mt */ 268.0, 271.0, 272.0
};

/*
 * corresponding table of VDW radii.
 * van der Waals radii are taken from A. Bondi,
 * J. Phys. Chem., 68, 441 - 452, 1964,
 * except the value for H, which is taken from R.S. Rowland & R. Taylor,
 * J.Phys.Chem., 100, 7384 - 7391, 1996. Radii that are not available in
 * either of these publications have RvdW = 2.00 \AA
 * The radii for Ions (Na, K, Cl, Ca, Mg, and Cs are based on the CHARMM27
 * Rmin/2 parameters for (SOD, POT, CLA, CAL, MG, CES) by default.
 */
constexpr double pte_vdw_radius[] = {
    /* X  */ 1.5, 1.2, 1.4, 1.82, 2.0, 2.0,
    /* C  */ 1.7, 1.55, 1.52, 1.47, 1.54,
    /* Na */ 1.36, 1.18, 2.0, 2.1, 1.8,
    /* S  */ 1.8, 2.27, 1.88, 1.76, 1.37, 2.0,
    /* Ti */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Ni */ 1.63, 1.4, 1.39, 1.07, 2.0, 1.85,
    /* Se */ 1.9, 1.85, 2.02, 2.0, 2.0, 2.0,
    /* Zr */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Pd */ 1.63, 1.72, 1.58, 1.93, 2.17, 2.0,
    /* Te */ 2.06, 1.98, 2.16, 2.1, 2.0,
    /* La */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Eu */ 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Er */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* W  */ 2.0, 2.0, 2.0, 2.0, 1.72, 1.66,
    /* Hg */ 1.55, 1.96, 2.02, 2.0, 2.0, 2.0, 2.0,
    /* Fr */ 2.0, 2.0, 2.0, 2.0, 2.0, 1.86,
    /* Np */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Md */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    /* Mt */ 2.0, 2.0, 2.0
};

// clang-format on

int get_pte_from_mass(double mass)
{
    int idx = 0;
    for (int i = 0; i < nr_pte_entries; ++i)
        if (fabs(mass - pte_mass[i]) < 0.65) idx = i;
    if ((mass > 0.0) && (mass < 2.2)) idx = 1;
    // discriminate between Cobalt and Nickel. The loop will detect Nickel
    if ((mass < 61.24) && (mass > 58.8133)) idx = 27;
    return idx;
}

QStringList defaultcolors = {"white", "gray",  "magenta", "cyan",   "yellow",
                             "blue",  "green", "red",     "orange", "brown"};

// constants
const QString blank(" ");
constexpr double VDW_ON           = 1.6;
constexpr double VDW_OFF          = 0.5;
constexpr double VDW_CUT          = 1.0;
constexpr double SHINY_ON         = 0.6;
constexpr double SHINY_OFF        = 0.2;
constexpr double SHINY_CUT        = 0.4;
constexpr double MAX_BOND_CUT     = 99.0;
constexpr int DEFAULT_BUFLEN      = 1024;
constexpr int DEFAULT_NPOINTS     = 100000;
constexpr double DEFAULT_DIAMETER = 0.2;

enum { FRAME, FILLED, POINTS };

} // namespace

class RegionInfo {
public:
    RegionInfo() = delete;
    RegionInfo(bool _enabled, int _style, const std::string &_color, double _diameter,
               int _npoints) :
        enabled(_enabled), style(_style), color(_color), diameter(_diameter), npoints(_npoints)
    {
    }

    bool enabled;
    int style;
    std::string color;
    double diameter;
    int npoints;
};

ImageViewer::ImageViewer(const QString &fileName, LammpsWrapper *_lammps, QWidget *parent) :
    QDialog(parent), menuBar(new QMenuBar), imageLabel(new QLabel), scrollArea(new QScrollArea),
    buttonBox(nullptr), scaleFactor(1.0), atomSize(1.0), saveAsAct(nullptr), copyAct(nullptr),
    cmdAct(nullptr), zoomInAct(nullptr), zoomOutAct(nullptr), normalSizeAct(nullptr),
    lammps(_lammps), group("all"), molecule("none"), filename(fileName), useelements(false),
    usediameter(false), usesigma(false)
{
    imageLabel->setBackgroundRole(QPalette::Base);
    imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    imageLabel->setScaledContents(true);
    imageLabel->minimumSizeHint();

    scrollArea->setBackgroundRole(QPalette::Dark);
    scrollArea->setWidget(imageLabel);
    scrollArea->setVisible(false);

    auto *mainLayout = new QVBoxLayout;

    QSettings settings;
    settings.beginGroup("snapshot");
    xsize       = settings.value("xsize", "600").toInt();
    ysize       = settings.value("ysize", "600").toInt();
    zoom        = settings.value("zoom", 1.0).toDouble();
    hrot        = settings.value("hrot", 60).toInt();
    vrot        = settings.value("vrot", 30).toInt();
    shinyfactor = settings.value("shinystyle", true).toBool() ? SHINY_ON : SHINY_OFF;
    vdwfactor   = settings.value("vdwstyle", false).toBool() ? VDW_ON : VDW_OFF;
    autobond    = settings.value("autobond", false).toBool();
    bondcutoff  = settings.value("bondcutoff", 1.6).toDouble();
    showbox     = settings.value("box", true).toBool();
    showaxes    = settings.value("axes", false).toBool();
    usessao     = settings.value("ssao", false).toBool();
    antialias   = settings.value("antialias", false).toBool();
    xcenter = ycenter = zcenter = 0.5;
    settings.endGroup();

    auto pix   = QPixmap(":/icons/emblem-photos.png");
    auto bsize = QFontMetrics(QApplication::font()).size(Qt::TextSingleLine, "Height:  200");

    auto *renderstatus = new QLabel(QString());
    renderstatus->setPixmap(pix.scaled(22, 22, Qt::KeepAspectRatio));
    renderstatus->setEnabled(false);
    renderstatus->setToolTip("Render status");
    renderstatus->setObjectName("renderstatus");
    auto *asize = new QLineEdit(QString::number(atomSize));
    auto *valid = new QDoubleValidator(1.0e-10, 1.0e10, 10, asize);
    asize->setValidator(valid);
    asize->setObjectName("atomSize");
    asize->setToolTip("Set Atom size");
    asize->setEnabled(false);
    asize->hide();

    auto *xval = new QSpinBox;
    xval->setRange(100, 10000);
    xval->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
    xval->setValue(xsize);
    xval->setObjectName("xsize");
    xval->setToolTip("Set rendered image width");
    xval->setMinimumSize(bsize);
    auto *yval = new QSpinBox;
    yval->setRange(100, 10000);
    yval->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
    yval->setValue(ysize);
    yval->setObjectName("ysize");
    yval->setToolTip("Set rendered image height");
    yval->setMinimumSize(bsize);

    connect(asize, &QLineEdit::editingFinished, this, &ImageViewer::set_atom_size);
    connect(xval, &QAbstractSpinBox::editingFinished, this, &ImageViewer::edit_size);
    connect(yval, &QAbstractSpinBox::editingFinished, this, &ImageViewer::edit_size);

    // workaround for incorrect highlight bug on macOS
    auto *dummy1 = new QPushButton(QIcon(), "");
    dummy1->hide();
    auto *dummy2 = new QPushButton(QIcon(), "");
    dummy2->hide();

    auto *dossao = new QPushButton(QIcon(":/icons/hd-img.png"), "");
    dossao->setCheckable(true);
    dossao->setToolTip("Toggle SSAO rendering");
    dossao->setObjectName("ssao");
    auto *doanti = new QPushButton(QIcon(":/icons/antialias.png"), "");
    doanti->setCheckable(true);
    doanti->setToolTip("Toggle anti-aliasing");
    doanti->setObjectName("antialias");
    auto *doshiny = new QPushButton(QIcon(":/icons/image-shiny.png"), "");
    doshiny->setCheckable(true);
    doshiny->setToolTip("Toggle shininess");
    doshiny->setObjectName("shiny");
    auto *dovdw = new QPushButton(QIcon(":/icons/vdw-style.png"), "");
    dovdw->setCheckable(true);
    dovdw->setToolTip("Toggle VDW style representation");
    dovdw->setObjectName("vdw");
    auto *dobond = new QPushButton(QIcon(":/icons/autobonds.png"), "");
    dobond->setCheckable(true);
    dobond->setToolTip("Toggle dynamic bond representation");
    dobond->setObjectName("autobond");
    auto *bondcut = new QLineEdit(QString::number(bondcutoff));
    bondcut->setMaxLength(5);
    bondcut->setObjectName("bondcut");
    bondcut->setToolTip("Set dynamic bond cutoff");
    QFontMetrics metrics(bondcut->fontMetrics());
    bondcut->setFixedSize(metrics.averageCharWidth() * 6, metrics.height() + 4);
    bondcut->setEnabled(false);
    auto *dobox = new QPushButton(QIcon(":/icons/system-box.png"), "");
    dobox->setCheckable(true);
    dobox->setToolTip("Toggle displaying box");
    dobox->setObjectName("box");
    auto *doaxes = new QPushButton(QIcon(":/icons/axes-img.png"), "");
    doaxes->setCheckable(true);
    doaxes->setToolTip("Toggle displaying axes");
    doaxes->setObjectName("axes");
    auto *zoomin = new QPushButton(QIcon(":/icons/gtk-zoom-in.png"), "");
    zoomin->setToolTip("Zoom in by 10 percent");
    auto *zoomout = new QPushButton(QIcon(":/icons/gtk-zoom-out.png"), "");
    zoomout->setToolTip("Zoom out by 10 percent");
    auto *rotleft = new QPushButton(QIcon(":/icons/object-rotate-left.png"), "");
    rotleft->setToolTip("Rotate left by 15 degrees");
    auto *rotright = new QPushButton(QIcon(":/icons/object-rotate-right.png"), "");
    rotright->setToolTip("Rotate right by 15 degrees");
    auto *rotup = new QPushButton(QIcon(":/icons/gtk-go-up.png"), "");
    rotup->setToolTip("Rotate up by 15 degrees");
    auto *rotdown = new QPushButton(QIcon(":/icons/gtk-go-down.png"), "");
    rotdown->setToolTip("Rotate down by 15 degrees");
    auto *recenter = new QPushButton(QIcon(":/icons/move-recenter.png"), "");
    recenter->setToolTip("Recenter on group");
    auto *reset = new QPushButton(QIcon(":/icons/gtk-zoom-fit.png"), "");
    reset->setToolTip("Reset view to defaults");
    auto *regviz = new QPushButton("Regions");
    regviz->setToolTip("Open dialog for visualizing regions");
    regviz->setObjectName("regions");
    regviz->setEnabled(false);

    constexpr int BUFLEN = 256;
    char gname[BUFLEN];
    auto *combo = new QComboBox;
    combo->setToolTip("Select group to display");
    combo->setObjectName("group");
    int ngroup = lammps->id_count("group");
    for (int i = 0; i < ngroup; ++i) {
        memset(gname, 0, BUFLEN);
        lammps->id_name("group", i, gname, BUFLEN);
        combo->addItem(gname);
    }

    auto *molbox = new QComboBox;
    molbox->setToolTip("Select molecule to display");
    molbox->setObjectName("molecule");
    molbox->addItem("none");
    int nmols = lammps->id_count("molecule");
    for (int i = 0; i < nmols; ++i) {
        memset(gname, 0, BUFLEN);
        lammps->id_name("molecule", i, gname, BUFLEN);
        molbox->addItem(gname);
    }

    auto *menuLayout   = new QHBoxLayout;
    auto *buttonLayout = new QHBoxLayout;
    auto *topLayout    = new QVBoxLayout;
    topLayout->addLayout(menuLayout);
    topLayout->addLayout(buttonLayout);

    menuLayout->addWidget(menuBar);
    menuLayout->addWidget(renderstatus);
    menuLayout->addWidget(new QLabel(" Atom Size: "));
    menuLayout->addWidget(asize);
    // hide item initially
    menuLayout->itemAt(2)->widget()->setObjectName("AtomLabel");
    menuLayout->itemAt(2)->widget()->hide();
    menuLayout->addWidget(new QLabel(" Width: "));
    menuLayout->addWidget(xval);
    menuLayout->addWidget(new QLabel(" Height: "));
    menuLayout->addWidget(yval);
    menuLayout->addWidget(dummy1);
    menuLayout->addWidget(new QLabel(" Group: "));
    menuLayout->addWidget(combo);
    menuLayout->addWidget(new QLabel(" Molecule: "));
    menuLayout->addWidget(molbox);
    buttonLayout->addWidget(dummy2);
    buttonLayout->addWidget(dossao);
    buttonLayout->addWidget(doanti);
    buttonLayout->addWidget(doshiny);
    buttonLayout->addWidget(dovdw);
    buttonLayout->addWidget(dobond);
    buttonLayout->addWidget(bondcut);
    buttonLayout->addWidget(dobox);
    buttonLayout->addWidget(doaxes);
    buttonLayout->addWidget(zoomin);
    buttonLayout->addWidget(zoomout);
    buttonLayout->addWidget(rotleft);
    buttonLayout->addWidget(rotright);
    buttonLayout->addWidget(rotup);
    buttonLayout->addWidget(rotdown);
    buttonLayout->addWidget(recenter);
    buttonLayout->addWidget(reset);
    buttonLayout->addWidget(regviz);
    buttonLayout->addStretch(1);

    connect(dossao, &QPushButton::released, this, &ImageViewer::toggle_ssao);
    connect(doanti, &QPushButton::released, this, &ImageViewer::toggle_anti);
    connect(doshiny, &QPushButton::released, this, &ImageViewer::toggle_shiny);
    connect(dovdw, &QPushButton::released, this, &ImageViewer::toggle_vdw);
    connect(dobond, &QPushButton::released, this, &ImageViewer::toggle_bond);
    connect(bondcut, &QLineEdit::editingFinished, this, &ImageViewer::set_bondcut);
    connect(dobox, &QPushButton::released, this, &ImageViewer::toggle_box);
    connect(doaxes, &QPushButton::released, this, &ImageViewer::toggle_axes);
    connect(zoomin, &QPushButton::released, this, &ImageViewer::do_zoom_in);
    connect(zoomout, &QPushButton::released, this, &ImageViewer::do_zoom_out);
    connect(rotleft, &QPushButton::released, this, &ImageViewer::do_rot_left);
    connect(rotright, &QPushButton::released, this, &ImageViewer::do_rot_right);
    connect(rotup, &QPushButton::released, this, &ImageViewer::do_rot_up);
    connect(rotdown, &QPushButton::released, this, &ImageViewer::do_rot_down);
    connect(recenter, &QPushButton::released, this, &ImageViewer::do_recenter);
    connect(reset, &QPushButton::released, this, &ImageViewer::reset_view);
    connect(regviz, &QPushButton::released, this, &ImageViewer::region_settings);
    connect(combo, SIGNAL(currentIndexChanged(int)), this, SLOT(change_group(int)));
    connect(molbox, SIGNAL(currentIndexChanged(int)), this, SLOT(change_molecule(int)));

    mainLayout->addLayout(topLayout);
    mainLayout->addWidget(scrollArea);
    setWindowIcon(QIcon(":/icons/lammps-icon-128x128.png"));
    setWindowTitle(QString("LAMMPS-GUI - Image Viewer - ") + QFileInfo(fileName).fileName());
    createActions();

    reset_view();
    // layout has not yet be established, so we need to fix up some pushbutton
    // properties directly since lookup in reset_view() will have failed
    dobox->setChecked(showbox);
    doshiny->setChecked(shinyfactor > SHINY_CUT);
    dovdw->setChecked(vdwfactor > VDW_CUT);
    dovdw->setEnabled(useelements || usediameter || usesigma);
    dobond->setChecked(autobond);
    doaxes->setChecked(showaxes);
    dossao->setChecked(usessao);
    doanti->setChecked(antialias);

    scaleFactor = 1.0;
    resize(image.width() + 25, image.height() + 80);

    scrollArea->setVisible(true);
    updateActions();
    setLayout(mainLayout);
    update_regions();
}

void ImageViewer::reset_view()
{
    QSettings settings;
    settings.beginGroup("snapshot");
    xsize       = settings.value("xsize", "600").toInt();
    ysize       = settings.value("ysize", "600").toInt();
    zoom        = settings.value("zoom", 1.0).toDouble();
    hrot        = settings.value("hrot", 60).toInt();
    vrot        = settings.value("vrot", 30).toInt();
    shinyfactor = settings.value("shinystyle", true).toBool() ? SHINY_ON : SHINY_OFF;
    vdwfactor   = settings.value("vdwstyle", false).toBool() ? VDW_ON : VDW_OFF;
    autobond    = settings.value("autobond", false).toBool();
    bondcutoff  = settings.value("bondcutoff", 1.6).toDouble();
    showbox     = settings.value("box", true).toBool();
    showaxes    = settings.value("axes", false).toBool();
    usessao     = settings.value("ssao", false).toBool();
    antialias   = settings.value("antialias", false).toBool();
    xcenter = ycenter = zcenter = 0.5;
    settings.endGroup();

    // reset state of checkable push buttons and combo box (if accessible)

    auto *field = findChild<QSpinBox *>("xsize");
    if (field) field->setValue(xsize);
    field = findChild<QSpinBox *>("ysize");
    if (field) field->setValue(ysize);

    auto *button = findChild<QPushButton *>("ssao");
    if (button) button->setChecked(usessao);
    button = findChild<QPushButton *>("antialias");
    if (button) button->setChecked(antialias);
    button = findChild<QPushButton *>("shiny");
    if (button) button->setChecked(shinyfactor > SHINY_CUT);
    button = findChild<QPushButton *>("vdw");
    if (button) button->setChecked(vdwfactor > VDW_CUT);
    button = findChild<QPushButton *>("autobond");
    if (button) button->setChecked(autobond);
    auto *cutoff = findChild<QLineEdit *>("bondcut");
    if (cutoff) {
        cutoff->setEnabled(autobond);
        cutoff->setText(QString::number(bondcutoff));
    }
    button = findChild<QPushButton *>("box");
    if (button) button->setChecked(showbox);
    button = findChild<QPushButton *>("axes");
    if (button) button->setChecked(showaxes);
    auto *cb = findChild<QComboBox *>("combo");
    if (cb) cb->setCurrentText("all");
    createImage();
}

void ImageViewer::set_atom_size()
{
    auto *field = qobject_cast<QLineEdit *>(sender());
    atomSize    = field->text().toDouble();
    createImage();
}

void ImageViewer::edit_size()
{
    auto *field = qobject_cast<QSpinBox *>(sender());
    if (field->objectName() == "xsize") {
        xsize = field->value();
    } else if (field->objectName() == "ysize") {
        ysize = field->value();
    }
    createImage();
}

void ImageViewer::toggle_ssao()
{
    auto *button = qobject_cast<QPushButton *>(sender());
    usessao      = !usessao;
    button->setChecked(usessao);
    createImage();
}

void ImageViewer::toggle_anti()
{
    auto *button = qobject_cast<QPushButton *>(sender());
    antialias    = !antialias;
    button->setChecked(antialias);
    createImage();
}

void ImageViewer::toggle_shiny()
{
    auto *button = qobject_cast<QPushButton *>(sender());
    if (shinyfactor > SHINY_CUT)
        shinyfactor = SHINY_OFF;
    else
        shinyfactor = SHINY_ON;
    button->setChecked(shinyfactor > SHINY_CUT);
    createImage();
}

void ImageViewer::toggle_vdw()
{
    auto *button = qobject_cast<QPushButton *>(sender());

    if (button->isChecked())
        vdwfactor = VDW_ON;
    else
        vdwfactor = VDW_OFF;

    // when enabling VDW rendering, we must turn off autobond
    bool do_vdw = vdwfactor > VDW_CUT;
    if (do_vdw) {
        autobond   = false;
        auto *bond = findChild<QPushButton *>("autobond");
        if (bond) bond->setChecked(false);
        auto *cutoff = findChild<QLineEdit *>("bondcut");
        if (cutoff) cutoff->setEnabled(false);
    }

    button->setChecked(do_vdw);
    createImage();
}

void ImageViewer::toggle_bond()
{
    auto *button = qobject_cast<QPushButton *>(sender());
    if (button) autobond = button->isChecked();
    auto *cutoff = findChild<QLineEdit *>("bondcut");
    if (cutoff) cutoff->setEnabled(autobond);
    set_bondcut();

    // when enabling autobond, we must turn off VDW
    if (autobond) {
        vdwfactor = VDW_OFF;
        auto *vdw = findChild<QPushButton *>("vdw");
        if (vdw) vdw->setChecked(false);
    }

    button->setChecked(autobond);
    createImage();
}

void ImageViewer::set_bondcut()
{
    auto *cutoff = findChild<QLineEdit *>("bondcut");
    if (cutoff) {
        auto *dptr            = (double *)lammps->extract_global("neigh_cutmax");
        double max_bondcutoff = (dptr) ? *dptr : 0.0;
        double new_bondcutoff = cutoff->text().toDouble();

        if ((max_bondcutoff > 0.1) && (new_bondcutoff > max_bondcutoff))
            new_bondcutoff = max_bondcutoff;
        if (new_bondcutoff > 0.1) bondcutoff = new_bondcutoff;

        cutoff->setText(QString::number(bondcutoff));
    }
    createImage();
}

void ImageViewer::toggle_box()
{
    auto *button = qobject_cast<QPushButton *>(sender());
    showbox      = !showbox;
    button->setChecked(showbox);
    createImage();
}

void ImageViewer::toggle_axes()
{
    auto *button = qobject_cast<QPushButton *>(sender());
    showaxes     = !showaxes;
    button->setChecked(showaxes);
    createImage();
}

void ImageViewer::do_zoom_in()
{
    zoom = zoom * 1.1;
    zoom = std::min(zoom, 10.0);
    createImage();
}

void ImageViewer::do_zoom_out()
{
    zoom = zoom / 1.1;
    zoom = std::max(zoom, 0.25);
    createImage();
}

void ImageViewer::do_rot_left()
{
    vrot -= 10;
    if (vrot < -180) vrot += 360;
    createImage();
}

void ImageViewer::do_rot_right()
{
    vrot += 10;
    if (vrot > 180) vrot -= 360;
    createImage();
}

void ImageViewer::do_rot_down()
{
    hrot -= 10;
    if (hrot < 0) hrot += 360;
    createImage();
}

void ImageViewer::do_rot_up()
{
    hrot += 10;
    if (hrot > 360) hrot -= 360;
    createImage();
}

void ImageViewer::do_recenter()
{
    QString commands = QString("variable LAMMPSGUI_CX delete\n"
                               "variable LAMMPSGUI_CY delete\n"
                               "variable LAMMPSGUI_CZ delete\n"
                               "variable LAMMPSGUI_CX equal (xcm(%1,x)-xlo)/lx\n"
                               "variable LAMMPSGUI_CY equal (xcm(%1,y)-ylo)/ly\n"
                               "variable LAMMPSGUI_CZ equal (xcm(%1,z)-zlo)/lz\n")
                           .arg(group);
    lammps->commands_string(commands);
    xcenter = lammps->extract_variable("LAMMPSGUI_CX");
    ycenter = lammps->extract_variable("LAMMPSGUI_CY");
    zcenter = lammps->extract_variable("LAMMPSGUI_CZ");
    lammps->commands_string("variable LAMMPSGUI_CX delete\n"
                            "variable LAMMPSGUI_CY delete\n"
                            "variable LAMMPSGUI_CZ delete\n");
    createImage();
}

void ImageViewer::cmd_to_clipboard()
{
    auto words = split_line(last_dump_cmd.toStdString());
    int modidx = 0;
    int maxidx = words.size();
    for (int i = 0; i < maxidx; ++i) {
        if (words[i] == "modify") {
            modidx = i;
            break;
        }
    }

    std::string dumpcmd = "dump viz ";
    dumpcmd += words[1];
    dumpcmd += " image 100 myimage-*.ppm";
    for (int i = 4; i < modidx; ++i)
        if (words[i] != "noinit") dumpcmd += " " + words[i];
    dumpcmd += '\n';

    dumpcmd += "dump_modify viz pad 9";
    for (int i = modidx + 1; i < maxidx; ++i)
        dumpcmd += " " + words[i];
    dumpcmd += '\n';
#if QT_CONFIG(clipboard)
    QGuiApplication::clipboard()->setText(dumpcmd.c_str());
#else
    fprintf(stderr, "# customized dump image command:\n%s", dumpcmd.c_str())
#endif
}

void ImageViewer::region_settings()
{
    update_regions();
    if (regions.size() == 0) return;
    QDialog regionview;
    regionview.setWindowTitle(QString("LAMMPS-GUI - Visualize Regions"));
    regionview.setWindowIcon(QIcon(":/icons/lammps-icon-128x128.png"));
    regionview.setMinimumSize(100, 50);
    regionview.setContentsMargins(5, 5, 5, 5);
    regionview.setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);

    auto *title = new QLabel("Visualize Regions:");
    title->setFrameStyle(QFrame::Panel | QFrame::Raised);
    title->setLineWidth(1);

    auto *layout = new QGridLayout;
    layout->addWidget(title, 0, 0, 1, 6, Qt::AlignHCenter);

    layout->addWidget(new QLabel("Region:"), 1, 0);
    layout->addWidget(new QLabel("Show:"), 1, 1, Qt::AlignHCenter);
    layout->addWidget(new QLabel("Style:"), 1, 2, Qt::AlignHCenter);
    layout->addWidget(new QLabel("Color:"), 1, 3, Qt::AlignHCenter);
    layout->addWidget(new QLabel("Size:"), 1, 4, Qt::AlignHCenter);
    layout->addWidget(new QLabel("# Points:"), 1, 5, Qt::AlignHCenter);

    auto *colorcompleter = new QColorCompleter;
    auto *colorvalidator = new QColorValidator;
    auto *framevalidator = new QDoubleValidator(1.0e-10, 1.0e10, 10);
    auto *pointvalidator = new QIntValidator(100, 1000000);
    QFontMetrics metrics(regionview.fontMetrics());

    int idx = 2;
    for (const auto &reg : regions) {
        layout->addWidget(new QLabel(reg.first.c_str()), idx, 0);
        layout->setObjectName(QString(reg.first.c_str()));

        auto *check = new QCheckBox("");
        check->setCheckState(reg.second->enabled ? Qt::Checked : Qt::Unchecked);
        layout->addWidget(check, idx, 1, Qt::AlignHCenter);
        auto *style = new QComboBox;
        style->setEditable(false);
        style->addItem("frame");
        style->addItem("filled");
        style->addItem("points");
        style->setCurrentIndex(reg.second->style);
        layout->addWidget(style, idx, 2);
        auto *color = new QLineEdit(reg.second->color.c_str());
        color->setCompleter(colorcompleter);
        color->setValidator(colorvalidator);
        color->setFixedSize(metrics.averageCharWidth() * 12, metrics.height() + 4);
        color->setText(reg.second->color.c_str());
        layout->addWidget(color, idx, 3);
        auto *frame = new QLineEdit(QString::number(reg.second->diameter));
        frame->setValidator(framevalidator);
        frame->setFixedSize(metrics.averageCharWidth() * 8, metrics.height() + 4);
        frame->setText(QString::number(reg.second->diameter));
        layout->addWidget(frame, idx, 4);
        auto *points = new QLineEdit(QString::number(reg.second->npoints));
        points->setValidator(pointvalidator);
        points->setFixedSize(metrics.averageCharWidth() * 10, metrics.height() + 4);
        points->setText(QString::number(reg.second->npoints));
        layout->addWidget(points, idx, 5);
        ++idx;
    }
    auto *cancel = new QPushButton("&Cancel");
    auto *apply  = new QPushButton("&Apply");
    cancel->setAutoDefault(false);
    apply->setAutoDefault(true);
    layout->addWidget(cancel, idx, 0, 1, 3, Qt::AlignHCenter);
    layout->addWidget(apply, idx, 3, 1, 3, Qt::AlignHCenter);
    connect(cancel, &QPushButton::released, &regionview, &QDialog::reject);
    connect(apply, &QPushButton::released, &regionview, &QDialog::accept);
    regionview.setLayout(layout);

    int rv = regionview.exec();

    // return immediately on cancel
    if (!rv) return;

    // retrieve data from dialog and store in map
    for (int idx = 2; idx < (int)regions.size() + 2; ++idx) {
        auto *item           = layout->itemAtPosition(idx, 0);
        auto *label          = qobject_cast<QLabel *>(item->widget());
        auto id              = label->text().toStdString();
        item                 = layout->itemAtPosition(idx, 1);
        auto *box            = qobject_cast<QCheckBox *>(item->widget());
        regions[id]->enabled = (box->checkState() == Qt::Checked);
        item                 = layout->itemAtPosition(idx, 2);
        auto *combo          = qobject_cast<QComboBox *>(item->widget());
        regions[id]->style   = combo->currentIndex();
        item                 = layout->itemAtPosition(idx, 3);
        auto *line           = qobject_cast<QLineEdit *>(item->widget());
        if (line && line->hasAcceptableInput()) regions[id]->color = line->text().toStdString();
        item = layout->itemAtPosition(idx, 4);
        line = qobject_cast<QLineEdit *>(item->widget());
        if (line && line->hasAcceptableInput()) regions[id]->diameter = line->text().toDouble();
        item = layout->itemAtPosition(idx, 5);
        line = qobject_cast<QLineEdit *>(item->widget());
        if (line && line->hasAcceptableInput()) regions[id]->npoints = line->text().toInt();
    }
    createImage();
}

void ImageViewer::change_group(int)
{
    auto *box = findChild<QComboBox *>("group");
    group     = box ? box->currentText() : "all";

    // reset molecule to "none" when changing group
    box = findChild<QComboBox *>("molecule");
    if (box && (box->currentIndex() > 0)) {
        box->setCurrentIndex(0); // triggers call to createImage()
    } else {
        createImage();
    }
}

void ImageViewer::change_molecule(int)
{
    auto *box = findChild<QComboBox *>("molecule");
    molecule  = box ? box->currentText() : "none";

    box = findChild<QComboBox *>("group");
    if (molecule == "none") {
        box->setEnabled(true);
    } else {
        box->setEnabled(false);
    }

    createImage();
}

// This function creates a visualization of the current system using the
// "dump image" command and reads and displays the renderd image.
// To visualize molecules we create new atoms with create_atoms and
// put them into a new, temporary group and then visualize that group.
// After rendering the image, the atoms and group are deleted.
// to update bond data, we also need to issue a "run 0" command.

void ImageViewer::createImage()
{
    auto *renderstatus = findChild<QLabel *>("renderstatus");
    if (renderstatus) renderstatus->setEnabled(true);
    repaint();

    QString oldgroup = group;

    if (molecule != "none") {

        // get center of box
        double *boxlo, *boxhi, xmid, ymid, zmid;
        boxlo = (double *)lammps->extract_global("boxlo");
        boxhi = (double *)lammps->extract_global("boxhi");
        if (boxlo && boxhi) {
            xmid = 0.5 * (boxhi[0] + boxlo[0]);
            ymid = 0.5 * (boxhi[1] + boxlo[1]);
            zmid = 0.5 * (boxhi[2] + boxlo[2]);
        } else {
            xmid = ymid = zmid = 0.0;
        }

        QString molcreate = "create_atoms 0 single %1 %2 %3 mol %4 312944 group %5 units box";
        group             = "imgviewer_tmp_mol";
        lammps->command(molcreate.arg(xmid).arg(ymid).arg(zmid).arg(molecule).arg(group));
        lammps->command(QString("neigh_modify exclude group all %1").arg(group));
        lammps->command("run 0 post no");
    }

    QSettings settings;
    QString dumpcmd = QString("write_dump ") + group + " image ";
    QDir dumpdir(QDir::tempPath());
    QFile dumpfile(dumpdir.absoluteFilePath(filename + ".ppm"));
    dumpcmd += "'" + dumpfile.fileName() + "'";

    settings.beginGroup("snapshot");
    int hhrot = (hrot > 180) ? 360 - hrot : hrot;

    // determine elements from masses and set their covalent radii
    int ntypes       = lammps->extract_setting("ntypes");
    int nbondtypes   = lammps->extract_setting("nbondtypes");
    auto *masses     = (double *)lammps->extract_atom("mass");
    QString units    = (const char *)lammps->extract_global("units");
    QString elements = "element ";
    QString adiams;
    useelements = false;
    if ((units == "real") || (units == "metal")) {
        useelements = true;
        for (int i = 1; i <= ntypes; ++i) {
            int idx = get_pte_from_mass(masses[i]);
            if (idx == 0) useelements = false;
            elements += QString(pte_label[idx]) + blank;
            adiams += QString("adiam %1 %2 ").arg(i).arg(vdwfactor * pte_vdw_radius[idx]);
        }
    }
    usediameter = lammps->extract_setting("radius_flag") != 0;
    // use Lennard-Jones sigma for radius, if available
    usesigma               = false;
    const char *pair_style = (const char *)lammps->extract_global("pair_style");
    if (!useelements && !usediameter && pair_style && (strncmp(pair_style, "lj/", 3) == 0)) {
        auto **sigma = (double **)lammps->extract_pair("sigma");
        if (sigma) {
            usesigma = true;
            for (int i = 1; i <= ntypes; ++i) {
                if (sigma[i][i] > 0.0)
                    adiams += QString("adiam %1 %2 ").arg(i).arg(vdwfactor * sigma[i][i]);
            }
        }
    }
    // adjust pushbutton state and clear adiams string to disable VDW display, if needed
    if (useelements || usediameter || usesigma) {
        auto *button = findChild<QPushButton *>("vdw");
        if (button) button->setEnabled(true);
        auto *edit = findChild<QLineEdit *>("atomSize");
        if (edit) {
            edit->setEnabled(false);
            edit->hide();
        }
        auto *label = findChild<QLabel *>("AtomLabel");
        if (label) {
            label->setEnabled(false);
            label->hide();
        }

    } else {
        adiams.clear();
        auto *button = findChild<QPushButton *>("vdw");
        if (button) button->setEnabled(false);

        auto *label = findChild<QLabel *>("AtomLabel");
        if (label) {
            label->setEnabled(true);
            label->show();
        }
        auto *edit = findChild<QLineEdit *>("atomSize");
        if (edit) {
            if (!edit->isEnabled()) {
                edit->setEnabled(true);
                edit->show();
                // initialize with lattice spacing
                const auto *xlattice = (const double *)lammps->extract_global("xlattice");
                if (xlattice) atomSize = *xlattice;
                edit->setText(QString::number(atomSize));
            }
            atomSize = edit->text().toDouble();
        }
        if (atomSize != 1.0) {
            for (int i = 1; i <= ntypes; ++i)
                adiams += QString("adiam %1 %2 ").arg(i).arg(atomSize);
        }
    }

    // color
    if (useelements)
        dumpcmd += blank + "element";
    else
        dumpcmd += blank + settings.value("color", "type").toString();

    bool do_vdw = vdwfactor > VDW_CUT;
    // diameter
    if (usediameter && do_vdw)
        dumpcmd += blank + "diameter";
    else
        dumpcmd += blank + settings.value("diameter", "type").toString();
    dumpcmd += QString(" size %1 %2").arg(xsize).arg(ysize);
    dumpcmd += QString(" zoom %1").arg(zoom);
    dumpcmd += QString(" shiny %1 ").arg(shinyfactor);
    dumpcmd += QString(" fsaa %1").arg(antialias ? "yes" : "no");
    if (nbondtypes > 0) {
        if (do_vdw)
            dumpcmd += " bond none none ";
        else
            dumpcmd += " bond atom 0.5 ";
    }
    if (lammps->extract_setting("dimension") == 3) {
        dumpcmd += QString(" view %1 %2").arg(hhrot).arg(vrot);
    }
    if (usessao) dumpcmd += " ssao yes 453983 0.75";
    if (showbox)
        dumpcmd += " box yes 0.025";
    else
        dumpcmd += " box no 0.0";

    if (showaxes)
        dumpcmd += " axes yes 0.5 0.025";
    else
        dumpcmd += " axes no 0.0 0.0";

    if (autobond) dumpcmd += blank + "autobond" + blank + QString::number(bondcutoff) + " 0.5";

    if (regions.size() > 0) {
        for (const auto &reg : regions) {
            if (reg.second->enabled) {
                QString id(reg.first.c_str());
                QString color(reg.second->color.c_str());
                switch (reg.second->style) {
                    case FRAME:
                        dumpcmd += " region " + id + blank + color;
                        dumpcmd += " frame " + QString::number(reg.second->diameter);
                        break;
                    case FILLED:
                        dumpcmd += " region " + id + blank + color + " filled";
                        break;
                    case POINTS:
                    default:
                        dumpcmd += " region " + id + blank + color;
                        dumpcmd += " points " + QString::number(reg.second->npoints);
                        dumpcmd += blank + QString::number(reg.second->diameter);
                        break;
                }
                dumpcmd += blank;
            }
        }
    }

    dumpcmd += QString(" center s %1 %2 %3").arg(xcenter).arg(ycenter).arg(zcenter);
    dumpcmd += " noinit";
    dumpcmd += " modify boxcolor " + settings.value("boxcolor", "yellow").toString();
    dumpcmd += " backcolor " + settings.value("background", "black").toString();
    if (useelements) dumpcmd += blank + elements + blank + adiams + blank;
    if (usesigma) dumpcmd += blank + adiams + blank;
    if (!useelements && !usesigma && (atomSize != 1.0)) dumpcmd += blank + adiams + blank;
    settings.endGroup();

    last_dump_cmd = dumpcmd;
    lammps->command(dumpcmd);

    QImageReader reader(dumpfile.fileName());
    reader.setAutoTransform(true);
    const QImage newImage = reader.read();
    dumpfile.remove();

    // read of new image failed. nothing left to do.
    if (newImage.isNull()) return;

    // show show image
    image = newImage;
    imageLabel->setPixmap(QPixmap::fromImage(image));
    imageLabel->adjustSize();
    if (renderstatus) renderstatus->setEnabled(false);
    repaint();

    if (molecule != "none") {
        lammps->command("neigh_modify exclude none");
        lammps->command(QString("delete_atoms group %1 compress no").arg(group));
        lammps->command(QString("group %1 delete").arg(group));
        group = oldgroup;
    }
}

void ImageViewer::saveAs()
{
    QString fileName = QFileDialog::getSaveFileName(this, "Save Image File As", QString(),
                                                    "Image Files (*.jpg *.png *.bmp *.ppm)");
    saveFile(fileName);
}

void ImageViewer::copy() {}

void ImageViewer::quit()
{
    auto *main = dynamic_cast<LammpsGui *>(get_main_widget());
    if (main) main->quit();
}

void ImageViewer::saveFile(const QString &fileName)
{
    if (!fileName.isEmpty()) image.save(fileName);
}

void ImageViewer::createActions()
{
    QMenu *fileMenu = menuBar->addMenu("&File");

    saveAsAct = fileMenu->addAction("&Save As...", this, &ImageViewer::saveAs);
    saveAsAct->setIcon(QIcon(":/icons/document-save-as.png"));
    saveAsAct->setEnabled(false);
    fileMenu->addSeparator();
    copyAct = fileMenu->addAction("&Copy Image", this, &ImageViewer::copy);
    copyAct->setIcon(QIcon(":/icons/edit-copy.png"));
    copyAct->setShortcut(QKeySequence::Copy);
    copyAct->setEnabled(false);
    cmdAct = fileMenu->addAction("Copy &dump image command", this, &ImageViewer::cmd_to_clipboard);
    cmdAct->setIcon(QIcon(":/icons/file-clipboard.png"));
    cmdAct->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_D));
    fileMenu->addSeparator();
    QAction *exitAct = fileMenu->addAction("&Close", this, &QWidget::close);
    exitAct->setIcon(QIcon(":/icons/window-close.png"));
    exitAct->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_W));
    QAction *quitAct = fileMenu->addAction("&Quit", this, &ImageViewer::quit);
    quitAct->setIcon(QIcon(":/icons/application-exit.png"));
    quitAct->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_Q));
}

void ImageViewer::updateActions()
{
    saveAsAct->setEnabled(!image.isNull());
    copyAct->setEnabled(!image.isNull());
}

void ImageViewer::scaleImage(double factor)
{
    scaleFactor *= factor;
    imageLabel->resize(scaleFactor * imageLabel->pixmap(Qt::ReturnByValue).size());

    adjustScrollBar(scrollArea->horizontalScrollBar(), factor);
    adjustScrollBar(scrollArea->verticalScrollBar(), factor);
}

void ImageViewer::adjustScrollBar(QScrollBar *scrollBar, double factor)
{
    scrollBar->setValue(
        int((factor * scrollBar->value()) + ((factor - 1) * scrollBar->pageStep() / 2)));
}

void ImageViewer::update_regions()
{
    if (!lammps) return;

    // remove any regions that no longer exist. to avoid inconsistencies while looping
    // over the regions, we first collect the list of missing id and then apply it.
    std::unordered_set<std::string> oldkeys;
    for (const auto &reg : regions) {
        if (!lammps->has_id("region", reg.first.c_str())) oldkeys.insert(reg.first);
    }
    for (const auto &id : oldkeys) {
        delete regions[id];
        regions.erase(id);
    }

    // add any new regions
    char buffer[DEFAULT_BUFLEN];
    int nregions = lammps->id_count("region");
    for (int i = 0; i < nregions; ++i) {
        if (lammps->id_name("region", i, buffer, DEFAULT_BUFLEN)) {
            std::string id = buffer;
            if (regions.count(id) == 0) {
                const auto &color = defaultcolors[i % defaultcolors.size()].toStdString();
                auto *reginfo =
                    new RegionInfo(false, FRAME, color, DEFAULT_DIAMETER, DEFAULT_NPOINTS);
                regions[id] = reginfo;
            }
        }
    }

    auto *button = findChild<QPushButton *>("regions");
    if (button) button->setEnabled(regions.size() > 0);
}

// Local Variables:
// c-basic-offset: 4
// End:
