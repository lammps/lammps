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

#include "lammpsgui.h"
#include "lammpswrapper.h"

#include <QAction>
#include <QApplication>
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
#include <QVBoxLayout>
#include <QVariant>

#include <cmath>

// clang-format off
/* periodic table of elements for translation of ordinal to atom type */
static const char *pte_label[] = {
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
static constexpr int nr_pte_entries = sizeof(pte_label) / sizeof(char *);

/* corresponding table of masses. */
static constexpr double pte_mass[] = {
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
static constexpr double pte_vdw_radius[] = {
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

static int get_pte_from_mass(double mass)
{
    int idx = 0;
    for (int i = 0; i < nr_pte_entries; ++i)
        if (fabs(mass - pte_mass[i]) < 0.65) idx = i;
    if ((mass > 0.0) && (mass < 2.2)) idx = 1;
    // discriminate between Cobalt and Nickel. The loop will detect Nickel
    if ((mass < 61.24) && (mass > 58.8133)) idx = 27;
    return idx;
}

static const QString blank(" ");

ImageViewer::ImageViewer(const QString &fileName, LammpsWrapper *_lammps, QWidget *parent) :
    QDialog(parent), menuBar(new QMenuBar), imageLabel(new QLabel), scrollArea(new QScrollArea),
    buttonBox(nullptr), scaleFactor(1.0), atomSize(1.0), saveAsAct(nullptr), copyAct(nullptr),
    cmdAct(nullptr), zoomInAct(nullptr), zoomOutAct(nullptr), normalSizeAct(nullptr),
    lammps(_lammps), group("all"), filename(fileName), useelements(false), usediameter(false),
    usesigma(false)
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

    vdwfactor   = 0.5;
    shinyfactor = 0.6;
    auto pix    = QPixmap(":/icons/emblem-photos.png");
    xcenter = ycenter = zcenter = 0.5;
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
    settings.beginGroup("snapshot");
    auto *xval = new QSpinBox;
    xval->setRange(100, 10000);
    xval->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
    xval->setValue(settings.value("xsize", "600").toInt());
    xval->setObjectName("xsize");
    xval->setToolTip("Set rendered image width");
    xval->setMinimumSize(bsize);
    auto *yval = new QSpinBox;
    yval->setRange(100, 10000);
    yval->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
    yval->setValue(settings.value("ysize", "600").toInt());
    yval->setObjectName("ysize");
    yval->setToolTip("Set rendered image height");
    yval->setMinimumSize(bsize);
    settings.endGroup();
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
    auto *combo = new QComboBox;
    combo->setObjectName("group");
    combo->setToolTip("Select group to display");
    combo->setObjectName("group");
    int ngroup           = lammps->id_count("group");
    constexpr int BUFLEN = 256;
    char gname[BUFLEN];
    for (int i = 0; i < ngroup; ++i) {
        memset(gname, 0, BUFLEN);
        lammps->id_name("group", i, gname, BUFLEN);
        combo->addItem(gname);
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
    buttonLayout->addWidget(dummy2);
    buttonLayout->addWidget(dossao);
    buttonLayout->addWidget(doanti);
    buttonLayout->addWidget(doshiny);
    buttonLayout->addWidget(dovdw);
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
    buttonLayout->addStretch(1);

    connect(dossao, &QPushButton::released, this, &ImageViewer::toggle_ssao);
    connect(doanti, &QPushButton::released, this, &ImageViewer::toggle_anti);
    connect(doshiny, &QPushButton::released, this, &ImageViewer::toggle_shiny);
    connect(dovdw, &QPushButton::released, this, &ImageViewer::toggle_vdw);
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
    connect(combo, SIGNAL(currentIndexChanged(int)), this, SLOT(change_group(int)));

    mainLayout->addLayout(topLayout);
    mainLayout->addWidget(scrollArea);
    setWindowIcon(QIcon(":/icons/lammps-icon-128x128.png"));
    setWindowTitle(QString("LAMMPS-GUI - Image Viewer - ") + QFileInfo(fileName).fileName());
    createActions();

    reset_view();
    // layout has not yet be established, so we need to fix up some pushbutton
    // properties directly since lookup in reset_view() will have failed
    dobox->setChecked(showbox);
    doshiny->setChecked(shinyfactor > 0.4);
    dovdw->setChecked(vdwfactor > 1.0);
    dovdw->setEnabled(useelements || usediameter || usesigma);
    doaxes->setChecked(showaxes);
    dossao->setChecked(usessao);
    doanti->setChecked(antialias);

    scaleFactor = 1.0;
    resize(image.width() + 25, image.height() + 80);

    scrollArea->setVisible(true);
    updateActions();
    setLayout(mainLayout);
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
    shinyfactor = settings.value("shinystyle", true).toBool() ? 0.6 : 0.2;
    vdwfactor   = settings.value("vdwstyle", false).toBool() ? 1.6 : 0.5;
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
    if (button) button->setChecked(shinyfactor > 0.4);
    button = findChild<QPushButton *>("vdw");
    if (button) button->setChecked(vdwfactor > 1.0);
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
    if (shinyfactor > 0.4)
        shinyfactor = 0.2;
    else
        shinyfactor = 0.6;
    button->setChecked(shinyfactor > 0.4);
    createImage();
}

void ImageViewer::toggle_vdw()
{
    auto *button = qobject_cast<QPushButton *>(sender());
    if (vdwfactor > 1.0)
        vdwfactor = 0.5;
    else
        vdwfactor = 1.6;
    button->setChecked(vdwfactor > 1.0);
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
    if (zoom > 5.0) zoom = 5.0;
    createImage();
}

void ImageViewer::do_zoom_out()
{
    zoom = zoom / 1.1;
    if (zoom < 0.5) zoom = 0.5;
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
                               "variable LAMMPSGUI_CZ equal (xcm(%1,z)-zlo)/lz\n").arg(group);
    lammps->commands_string(commands.toLocal8Bit());
    xcenter = lammps->extract_variable("LAMMPSGUI_CX");
    ycenter = lammps->extract_variable("LAMMPSGUI_CZ");
    zcenter = lammps->extract_variable("LAMMPSGUI_CZ");
    lammps->commands_string("variable LAMMPSGUI_CX delete\n"
                            "variable LAMMPSGUI_CY delete\n"
                            "variable LAMMPSGUI_CZ delete\n");
    createImage();
}

void ImageViewer::cmd_to_clipboard()
{
    auto words    = last_dump_cmd.split(" ");
    QString blank = QStringLiteral(" ");
    int modidx    = words.indexOf("modify");
    int maxidx    = words.size();

    QString dumpcmd = "dump viz ";
    dumpcmd += words[1] + " image 100 myimage-*.ppm";
    for (int i = 4; i < modidx; ++i)
        if (words[i] != "noinit") dumpcmd += blank + words[i];
    dumpcmd += '\n';

    dumpcmd += "dump_modify viz pad 9";
    for (int i = modidx + 1; i < maxidx; ++i)
        dumpcmd += blank + words[i];
    dumpcmd += '\n';
#if QT_CONFIG(clipboard)
    QGuiApplication::clipboard()->setText(dumpcmd);
#endif
}

void ImageViewer::change_group(int)
{
    auto *box = findChild<QComboBox *>("group");
    if (box) group = box->currentText();
    createImage();
}

void ImageViewer::createImage()
{
    auto *renderstatus = findChild<QLabel *>("renderstatus");
    if (renderstatus) renderstatus->setEnabled(true);
    repaint();

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
        double **sigma = (double **)lammps->extract_pair("sigma");
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
                auto *xlattice = (const double *)lammps->extract_global("xlattice");
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

    // diameter
    if (usediameter && (vdwfactor > 1.0))
        dumpcmd += blank + "diameter";
    else
        dumpcmd += blank + settings.value("diameter", "type").toString();
    dumpcmd += QString(" size %1 %2").arg(xsize).arg(ysize);
    dumpcmd += QString(" zoom %1").arg(zoom);
    dumpcmd += QString(" shiny %1 ").arg(shinyfactor);
    dumpcmd += QString(" fsaa %1").arg(antialias ? "yes" : "no");
    if (nbondtypes > 0) {
        if (vdwfactor > 1.0)
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

    dumpcmd += QString(" center s %1 %2 %3").arg(xcenter).arg(ycenter).arg(zcenter);
    dumpcmd += " noinit";
    dumpcmd += " modify boxcolor " + settings.value("boxcolor", "yellow").toString();
    dumpcmd += " backcolor " + settings.value("background", "black").toString();
    if (useelements) dumpcmd += blank + elements + blank + adiams + blank;
    if (usesigma) dumpcmd += blank + adiams + blank;
    if (!useelements && !usesigma && (atomSize != 1.0)) dumpcmd += blank + adiams + blank;
    settings.endGroup();

    last_dump_cmd = dumpcmd;
    lammps->command(dumpcmd.toLocal8Bit());

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
    LammpsGui *main = nullptr;
    for (QWidget *widget : QApplication::topLevelWidgets())
        if (widget->objectName() == "LammpsGui") main = dynamic_cast<LammpsGui *>(widget);
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
#if QT_VERSION < QT_VERSION_CHECK(5, 15, 0)
    imageLabel->resize(scaleFactor * imageLabel->pixmap()->size());
#else
    imageLabel->resize(scaleFactor * imageLabel->pixmap(Qt::ReturnByValue).size());
#endif

    adjustScrollBar(scrollArea->horizontalScrollBar(), factor);
    adjustScrollBar(scrollArea->verticalScrollBar(), factor);
}

void ImageViewer::adjustScrollBar(QScrollBar *scrollBar, double factor)
{
    scrollBar->setValue(
        int(factor * scrollBar->value() + ((factor - 1) * scrollBar->pageStep() / 2)));
}

// Local Variables:
// c-basic-offset: 4
// End:
