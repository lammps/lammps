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
#include "lammpswrapper.h"

#include <QAction>
#include <QDialogButtonBox>
#include <QDir>
#include <QFileDialog>
#include <QGuiApplication>
#include <QImage>
#include <QImageReader>
#include <QLabel>
#include <QLineEdit>
#include <QMenuBar>
#include <QMessageBox>
#include <QPalette>
#include <QPoint>
#include <QPushButton>
#include <QScreen>
#include <QScrollArea>
#include <QScrollBar>
#include <QSettings>
#include <QSpinBox>
#include <QStatusBar>
#include <QVBoxLayout>
#include <QWheelEvent>
#include <QWidgetAction>

#include <cmath>

extern "C" {
#include "periodic_table.h"
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
}

static const QString blank(" ");

ImageViewer::ImageViewer(const QString &fileName, LammpsWrapper *_lammps, QWidget *parent) :
    QDialog(parent), imageLabel(new QLabel), scrollArea(new QScrollArea), menuBar(new QMenuBar),
    lammps(_lammps), group("all"), filename(fileName), useelements(false)
{
    imageLabel->setBackgroundRole(QPalette::Base);
    imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    imageLabel->setScaledContents(true);
    imageLabel->minimumSizeHint();

    scrollArea->setBackgroundRole(QPalette::Dark);
    scrollArea->setWidget(imageLabel);
    scrollArea->setVisible(false);

    buttonBox = new QDialogButtonBox(QDialogButtonBox::Close);

    connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);

    QVBoxLayout *mainLayout = new QVBoxLayout;

    QSettings settings;

    vdwfactor = 0.5;
    auto pix  = QPixmap(":/emblem-photos.png");

    auto *renderstatus = new QLabel(QString());
    renderstatus->setPixmap(pix.scaled(22, 22, Qt::KeepAspectRatio));
    renderstatus->setEnabled(false);
    renderstatus->setToolTip("Render status");
    renderstatus->setObjectName("renderstatus");
    settings.beginGroup("snapshot");
    auto *xval = new QSpinBox;
    xval->setRange(100, 10000);
    xval->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
    xval->setValue(settings.value("xsize", "800").toInt());
    xval->setObjectName("xsize");
    xval->setToolTip("Set rendered image width");
    auto *yval = new QSpinBox;
    yval->setRange(100, 10000);
    yval->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
    yval->setValue(settings.value("ysize", "600").toInt());
    yval->setObjectName("ysize");
    yval->setToolTip("Set rendered image height");
    settings.endGroup();
    connect(xval, &QAbstractSpinBox::editingFinished, this, &ImageViewer::edit_size);
    connect(yval, &QAbstractSpinBox::editingFinished, this, &ImageViewer::edit_size);

    // workaround for incorrect highlight bug on macOS
    auto *dummy = new QPushButton(QIcon(), "");
    dummy->hide();

    auto *dossao = new QPushButton(QIcon(":/hd-img.png"), "");
    dossao->setCheckable(true);
    dossao->setToolTip("Toggle SSAO rendering");
    dossao->setObjectName("ssao");
    auto *doanti = new QPushButton(QIcon(":/antialias.png"), "");
    doanti->setCheckable(true);
    doanti->setToolTip("Toggle anti-aliasing");
    doanti->setObjectName("antialias");
    auto *dovdw = new QPushButton(QIcon(":/vdw-style.png"), "");
    dovdw->setCheckable(true);
    dovdw->setToolTip("Toggle VDW style representation");
    dovdw->setObjectName("vdw");
    auto *dobox = new QPushButton(QIcon(":/system-box.png"), "");
    dobox->setCheckable(true);
    dobox->setToolTip("Toggle displaying box");
    dobox->setObjectName("box");
    auto *doaxes = new QPushButton(QIcon(":/axes-img.png"), "");
    doaxes->setCheckable(true);
    doaxes->setToolTip("Toggle displaying axes");
    doaxes->setObjectName("axes");
    auto *zoomin = new QPushButton(QIcon(":/gtk-zoom-in.png"), "");
    zoomin->setToolTip("Zoom in by 10 percent");
    auto *zoomout = new QPushButton(QIcon(":/gtk-zoom-out.png"), "");
    zoomout->setToolTip("Zoom out by 10 percent");
    auto *rotleft = new QPushButton(QIcon(":/object-rotate-left.png"), "");
    rotleft->setToolTip("Rotate left by 15 degrees");
    auto *rotright = new QPushButton(QIcon(":/object-rotate-right.png"), "");
    rotright->setToolTip("Rotate right by 15 degrees");
    auto *rotup = new QPushButton(QIcon(":/gtk-go-up.png"), "");
    rotup->setToolTip("Rotate up by 15 degrees");
    auto *rotdown = new QPushButton(QIcon(":/gtk-go-down.png"), "");
    rotdown->setToolTip("Rotate down by 15 degrees");
    auto *reset = new QPushButton(QIcon(":/gtk-zoom-fit.png"), "");
    reset->setToolTip("Reset view to defaults");
    auto *combo = new QComboBox;
    combo->setObjectName("group");
    combo->setToolTip("Select group to display");
    combo->setObjectName("group");
    int ngroup = lammps->id_count("group");
    char gname[64];
    for (int i = 0; i < ngroup; ++i) {
        lammps->id_name("group", i, gname, 64);
        combo->addItem(gname);
    }

    QHBoxLayout *menuLayout = new QHBoxLayout;
    menuLayout->addWidget(menuBar);
    menuLayout->addWidget(renderstatus);
    menuLayout->addWidget(new QLabel(" Width: "));
    menuLayout->addWidget(xval);
    menuLayout->addWidget(new QLabel(" Height: "));
    menuLayout->addWidget(yval);
    menuLayout->addWidget(dummy);
    menuLayout->addWidget(dossao);
    menuLayout->addWidget(doanti);
    menuLayout->addWidget(dovdw);
    menuLayout->addWidget(dobox);
    menuLayout->addWidget(doaxes);
    menuLayout->addWidget(zoomin);
    menuLayout->addWidget(zoomout);
    menuLayout->addWidget(rotleft);
    menuLayout->addWidget(rotright);
    menuLayout->addWidget(rotup);
    menuLayout->addWidget(rotdown);
    menuLayout->addWidget(reset);
    menuLayout->addWidget(new QLabel(" Group: "));
    menuLayout->addWidget(combo);

    connect(dossao, &QPushButton::released, this, &ImageViewer::toggle_ssao);
    connect(doanti, &QPushButton::released, this, &ImageViewer::toggle_anti);
    connect(dovdw, &QPushButton::released, this, &ImageViewer::toggle_vdw);
    connect(dobox, &QPushButton::released, this, &ImageViewer::toggle_box);
    connect(doaxes, &QPushButton::released, this, &ImageViewer::toggle_axes);
    connect(zoomin, &QPushButton::released, this, &ImageViewer::do_zoom_in);
    connect(zoomout, &QPushButton::released, this, &ImageViewer::do_zoom_out);
    connect(rotleft, &QPushButton::released, this, &ImageViewer::do_rot_left);
    connect(rotright, &QPushButton::released, this, &ImageViewer::do_rot_right);
    connect(rotup, &QPushButton::released, this, &ImageViewer::do_rot_up);
    connect(rotdown, &QPushButton::released, this, &ImageViewer::do_rot_down);
    connect(reset, &QPushButton::released, this, &ImageViewer::reset_view);
    connect(combo, SIGNAL(currentIndexChanged(int)), this, SLOT(change_group(int)));

    mainLayout->addLayout(menuLayout);
    mainLayout->addWidget(scrollArea);
    mainLayout->addWidget(buttonBox);
    setWindowIcon(QIcon(":/lammps-icon-128x128.png"));
    setWindowTitle(QString("Image Viewer: ") + QFileInfo(fileName).fileName());
    createActions();

    reset_view();
    // layout has not yet be established, so we need to fix up some pushbutton
    // properties directly since lookup in reset_view() will have failed
    dobox->setChecked(showbox);
    dovdw->setChecked(vdwfactor > 1.0);
    dovdw->setEnabled(useelements);
    doaxes->setChecked(showaxes);
    dossao->setChecked(usessao);
    doanti->setChecked(antialias);

    scaleFactor = 1.0;
    resize(image.width() + 20, image.height() + 50);

    scrollArea->setVisible(true);
    updateActions();
    setLayout(mainLayout);
}

void ImageViewer::reset_view()
{
    QSettings settings;
    settings.beginGroup("snapshot");
    xsize     = settings.value("xsize", "800").toInt();
    ysize     = settings.value("ysize", "600").toInt();
    zoom      = settings.value("zoom", 1.0).toDouble();
    hrot      = settings.value("hrot", 60).toInt();
    vrot      = settings.value("vrot", 30).toInt();
    vdwfactor = settings.value("vdwstyle", false).toBool() ? 1.6 : 0.5;
    showbox   = settings.value("box", true).toBool();
    showaxes  = settings.value("axes", false).toBool();
    usessao   = settings.value("ssao", false).toBool();
    antialias = settings.value("antialias", false).toBool();
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

void ImageViewer::edit_size()
{
    QSpinBox *field = qobject_cast<QSpinBox *>(sender());
    if (field->objectName() == "xsize") {
        xsize = field->value();
    } else if (field->objectName() == "ysize") {
        ysize = field->value();
    }
    createImage();
}

void ImageViewer::toggle_ssao()
{
    QPushButton *button = qobject_cast<QPushButton *>(sender());
    usessao             = !usessao;
    button->setChecked(usessao);
    createImage();
}

void ImageViewer::toggle_anti()
{
    QPushButton *button = qobject_cast<QPushButton *>(sender());
    antialias           = !antialias;
    button->setChecked(antialias);
    createImage();
}

void ImageViewer::toggle_vdw()
{
    QPushButton *button = qobject_cast<QPushButton *>(sender());
    if (vdwfactor > 1.0)
        vdwfactor = 0.5;
    else
        vdwfactor = 1.6;
    button->setChecked(vdwfactor > 1.0);
    createImage();
}

void ImageViewer::toggle_box()
{
    QPushButton *button = qobject_cast<QPushButton *>(sender());
    showbox             = !showbox;
    button->setChecked(showbox);
    createImage();
}

void ImageViewer::toggle_axes()
{
    QPushButton *button = qobject_cast<QPushButton *>(sender());
    showaxes            = !showaxes;
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

void ImageViewer::change_group(int idx)
{
    QComboBox *box = findChild<QComboBox *>("group");
    if (box) group = box->currentText();
    createImage();
}

void ImageViewer::createImage()
{
    QLabel *renderstatus = findChild<QLabel *>("renderstatus");
    if (renderstatus) renderstatus->setEnabled(true);
    repaint();

    QSettings settings;
    QString dumpcmd = QString("write_dump ") + group + " image ";
    QDir dumpdir(QDir::tempPath());
    QFile dumpfile(dumpdir.absoluteFilePath(filename + ".ppm"));
    dumpcmd += "'" + dumpfile.fileName() + "'";

    settings.beginGroup("snapshot");
    int aa       = antialias ? 2 : 1;
    int tmpxsize = xsize * aa;
    int tmpysize = ysize * aa;
    int hhrot    = (hrot > 180) ? 360 - hrot : hrot;

    // determine elements from masses and set their covalent radii
    int ntypes       = lammps->extract_setting("ntypes");
    int nbondtypes   = lammps->extract_setting("nbondtypes");
    double *masses   = (double *)lammps->extract_atom("mass");
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

    // adjust pushbutton state and clear adiams string to disable VDW display, if needed
    if (useelements) {
        auto *button = findChild<QPushButton *>("vdw");
        if (button) button->setEnabled(true);
    } else {
        adiams.clear();
        auto *button = findChild<QPushButton *>("vdw");
        if (button) button->setEnabled(false);
    }

    if (!adiams.isEmpty())
        dumpcmd += blank + "element";
    else
        dumpcmd += blank + settings.value("color", "type").toString();
    dumpcmd += blank + settings.value("diameter", "type").toString();
    dumpcmd += QString(" size ") + QString::number(tmpxsize) + blank + QString::number(tmpysize);
    dumpcmd += QString(" zoom ") + QString::number(zoom);
    dumpcmd += " shiny 0.5 ";
    if (nbondtypes > 0) {
        if (vdwfactor > 1.0)
            dumpcmd += " bond none none ";
        else
            dumpcmd += " bond atom 0.5 ";
    }
    if (lammps->extract_setting("dimension") == 3) {
        dumpcmd += QString(" view ") + QString::number(hhrot) + blank + QString::number(vrot);
    }
    if (usessao) dumpcmd += QString(" ssao yes 453983 0.75");
    if (showbox)
        dumpcmd += QString(" box yes 0.025");
    else
        dumpcmd += QString(" box no 0.0");

    if (showaxes)
        dumpcmd += QString(" axes yes 0.5 0.025");
    else
        dumpcmd += QString(" axes no 0.0 0.0");

    dumpcmd += " modify boxcolor " + settings.value("boxcolor", "yellow").toString();
    dumpcmd += " backcolor " + settings.value("background", "black").toString();
    if (!adiams.isEmpty()) dumpcmd += blank + elements + blank + adiams + blank;
    settings.endGroup();

    lammps->command(dumpcmd.toLocal8Bit());

    QImageReader reader(dumpfile.fileName());
    reader.setAutoTransform(true);
    const QImage newImage = reader.read();

    if (newImage.isNull()) {
        QMessageBox::warning(
            this, QGuiApplication::applicationDisplayName(),
            QString("Cannot load %1: %2").arg(dumpfile.fileName(), reader.errorString()));
        return;
    }
    dumpfile.remove();

    // scale back to achieve antialiasing
    image = newImage.scaled(xsize, ysize, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
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

void ImageViewer::saveFile(const QString &fileName)
{
    if (!fileName.isEmpty()) image.save(fileName);
}

void ImageViewer::createActions()
{
    QMenu *fileMenu = menuBar->addMenu("&File");

    saveAsAct = fileMenu->addAction("&Save As...", this, &ImageViewer::saveAs);
    saveAsAct->setIcon(QIcon(":/document-save-as.png"));
    saveAsAct->setEnabled(false);
    fileMenu->addSeparator();
    copyAct = fileMenu->addAction("&Copy", this, &ImageViewer::copy);
    copyAct->setIcon(QIcon(":/edit-copy.png"));
    copyAct->setShortcut(QKeySequence::Copy);
    copyAct->setEnabled(false);
    fileMenu->addSeparator();
    QAction *exitAct = fileMenu->addAction("&Close", this, &QWidget::close);
    exitAct->setIcon(QIcon(":/window-close.png"));
    exitAct->setShortcut(QKeySequence::fromString("Ctrl+W"));
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
