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
#include <QMenuBar>
#include <QMessageBox>
#include <QPalette>
#include <QPoint>
#include <QPushButton>
#include <QScreen>
#include <QScrollArea>
#include <QScrollBar>
#include <QSettings>
#include <QStatusBar>
#include <QVBoxLayout>
#include <QWheelEvent>
#include <QWidgetAction>

static const QString blank(" ");

ImageViewer::ImageViewer(const QString &fileName, LammpsWrapper *_lammps, QWidget *parent) :
    QDialog(parent), imageLabel(new QLabel), scrollArea(new QScrollArea), menuBar(new QMenuBar),
    lammps(_lammps), group("all"), filename(fileName)
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

    auto *dossao = new QPushButton(QIcon(":/hd-img.png"), "");
    dossao->setCheckable(true);
    auto *doanti = new QPushButton(QIcon(":/antialias.png"), "");
    doanti->setCheckable(true);
    auto *dobox = new QPushButton(QIcon(":/system-box.png"), "");
    dobox->setCheckable(true);
    auto *doaxes = new QPushButton(QIcon(":/axes-img.png"), "");
    doaxes->setCheckable(true);
    auto *zoomin   = new QPushButton(QIcon(":/gtk-zoom-in.png"), "");
    auto *zoomout  = new QPushButton(QIcon(":/gtk-zoom-out.png"), "");
    auto *rotleft  = new QPushButton(QIcon(":/object-rotate-left.png"), "");
    auto *rotright = new QPushButton(QIcon(":/object-rotate-right.png"), "");
    auto *rotup    = new QPushButton(QIcon(":/gtk-go-up.png"), "");
    auto *rotdown  = new QPushButton(QIcon(":/gtk-go-down.png"), "");
    auto *reset    = new QPushButton(QIcon(":/gtk-zoom-fit.png"), "");
    auto *combo    = new QComboBox;
    combo->setObjectName("group");
    int ngroup = lammps->id_count("group");
    char gname[64];
    for (int i = 0; i < ngroup; ++i) {
        lammps->id_name("group", i, gname, 64);
        combo->addItem(gname);
    }

    QHBoxLayout *menuLayout = new QHBoxLayout;
    menuLayout->addWidget(menuBar);
    menuLayout->addWidget(dossao);
    menuLayout->addWidget(doanti);
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
    setWindowTitle(QString("Image Viewer: ") + QFileInfo(fileName).fileName());
    createActions();

    reset_view();
    dobox->setChecked(showbox);
    doaxes->setChecked(showaxes);
    dossao->setChecked(usessao);
    doanti->setChecked(antialias);

    scaleFactor = 1.0;
    resize(image.width() + 20, image.height() + 50);

    scrollArea->setVisible(true);
    fitToWindowAct->setEnabled(true);
    updateActions();
    if (!fitToWindowAct->isChecked()) imageLabel->adjustSize();
    setLayout(mainLayout);
}

void ImageViewer::reset_view()
{
    QSettings settings;
    settings.beginGroup("snapshot");
    zoom      = settings.value("zoom", 1.0).toDouble();
    hrot      = settings.value("hrot", 60).toInt();
    vrot      = settings.value("vrot", 30).toInt();
    showbox   = settings.value("box", true).toBool();
    showaxes  = settings.value("axes", false).toBool();
    usessao   = settings.value("ssao", false).toBool();
    antialias = settings.value("antialias", false).toBool();
    settings.endGroup();

    // reset state of checkable push buttons and combo box (after main layout is set up)
    auto *lo = layout();
    if (lo) {
        // grab layout manager for the top bar
        lo = lo->itemAt(0)->layout();
        // grab the first 4 buttons after the menu bar
        auto *button = qobject_cast<QPushButton *>(lo->itemAt(1)->widget());
        button->setChecked(usessao);
        button = qobject_cast<QPushButton *>(lo->itemAt(2)->widget());
        button->setChecked(antialias);
        button = qobject_cast<QPushButton *>(lo->itemAt(3)->widget());
        button->setChecked(showbox);
        button = qobject_cast<QPushButton *>(lo->itemAt(4)->widget());
        button->setChecked(showaxes);
        // grab the last entry -> group selector
        auto *cb = qobject_cast<QComboBox *>(lo->itemAt(lo->count() - 1)->widget());
        cb->setCurrentText("all");
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
    vrot -= 15;
    if (vrot < -180) vrot += 360;
    createImage();
}

void ImageViewer::do_rot_right()
{
    vrot += 15;
    if (vrot > 180) vrot -= 360;
    createImage();
}

void ImageViewer::do_rot_down()
{
    hrot -= 15;
    if (hrot < 0) hrot += 360;
    createImage();
}

void ImageViewer::do_rot_up()
{
    hrot += 15;
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
    QSettings settings;
    QString dumpcmd = QString("write_dump ") + group + " image ";
    QDir dumpdir(QDir::tempPath());
    QFile dumpfile(dumpdir.absoluteFilePath(filename + ".ppm"));
    dumpcmd += dumpfile.fileName();

    settings.beginGroup("snapshot");
    int aa    = antialias ? 2 : 1;
    int xsize = settings.value("xsize", 800).toInt() * aa;
    int ysize = settings.value("ysize", 600).toInt() * aa;
    int hhrot = (hrot > 180) ? 360 - hrot : hrot;

    dumpcmd += blank + settings.value("color", "type").toString();
    dumpcmd += blank + settings.value("diameter", "type").toString();
    dumpcmd += QString(" size ") + QString::number(xsize) + blank + QString::number(ysize);
    dumpcmd += QString(" zoom ") + QString::number(zoom);
    lammps->command(dumpcmd.toLocal8Bit());
    if (lammps->extract_setting("dimension") == 3) {
        dumpcmd += QString(" view ") + QString::number(hhrot) + blank + QString::number(vrot);
    }
    if (usessao) dumpcmd += QString(" ssao yes 453983 0.75");
    if (showbox)
        dumpcmd += QString(" box yes 0.025");
    else
        dumpcmd += QString(" box no 0.0");

    if (showaxes)
        dumpcmd += QString(" axes yes 0.2 0.025");
    else
        dumpcmd += QString(" axes no 0.0 0.0");

    dumpcmd += " modify boxcolor " + settings.value("boxcolor", "yellow").toString();
    dumpcmd += " backcolor " + settings.value("background", "black").toString();
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

    settings.beginGroup("snapshot");
    xsize = settings.value("xsize", 800).toInt();
    ysize = settings.value("ysize", 600).toInt();
    settings.endGroup();
    // scale back to achieve antialiasing
    image = newImage.scaled(xsize, ysize, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
    imageLabel->setPixmap(QPixmap::fromImage(image));
}

void ImageViewer::saveAs()
{
    QString fileName = QFileDialog::getSaveFileName(this, "Save Image File As", QString(),
                                                    "Image Files (*.jpg *.png *.bmp *.ppm)");
    saveFile(fileName);
}

void ImageViewer::copy() {}

void ImageViewer::zoomIn()
{
    scaleImage(1.25);
}

void ImageViewer::zoomOut()
{
    scaleImage(0.8);
}

void ImageViewer::normalSize()
{
    imageLabel->adjustSize();
    scaleFactor = 1.0;
}

void ImageViewer::fitToWindow()
{
    bool fitToWindow = fitToWindowAct->isChecked();
    scrollArea->setWidgetResizable(fitToWindow);
    if (!fitToWindow) normalSize();
    updateActions();
}

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

    QMenu *viewMenu = menuBar->addMenu("&View");

    zoomInAct = viewMenu->addAction("Image Zoom &In (25%)", this, &ImageViewer::zoomIn);
    zoomInAct->setShortcut(QKeySequence::ZoomIn);
    zoomInAct->setIcon(QIcon(":/gtk-zoom-in.png"));
    zoomInAct->setEnabled(false);

    zoomOutAct = viewMenu->addAction("Image Zoom &Out (25%)", this, &ImageViewer::zoomOut);
    zoomOutAct->setShortcut(QKeySequence::ZoomOut);
    zoomOutAct->setIcon(QIcon(":/gtk-zoom-out.png"));
    zoomOutAct->setEnabled(false);

    normalSizeAct = viewMenu->addAction("&Reset Image Size", this, &ImageViewer::normalSize);
    normalSizeAct->setShortcut(QKeySequence::fromString("Ctrl+0"));
    normalSizeAct->setIcon(QIcon(":/gtk-zoom-fit.png"));
    normalSizeAct->setEnabled(false);

    viewMenu->addSeparator();

    fitToWindowAct = viewMenu->addAction("&Fit to Window", this, &ImageViewer::fitToWindow);
    fitToWindowAct->setEnabled(false);
    fitToWindowAct->setCheckable(true);
    fitToWindowAct->setShortcut(QKeySequence::fromString("Ctrl+="));
}

void ImageViewer::updateActions()
{
    saveAsAct->setEnabled(!image.isNull());
    copyAct->setEnabled(!image.isNull());
    zoomInAct->setEnabled(!fitToWindowAct->isChecked());
    zoomOutAct->setEnabled(!fitToWindowAct->isChecked());
    normalSizeAct->setEnabled(!fitToWindowAct->isChecked());
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
    zoomInAct->setEnabled(scaleFactor < 3.0);
    zoomOutAct->setEnabled(scaleFactor > 0.333);
}

void ImageViewer::adjustScrollBar(QScrollBar *scrollBar, double factor)
{
    scrollBar->setValue(
        int(factor * scrollBar->value() + ((factor - 1) * scrollBar->pageStep() / 2)));
}

// Local Variables:
// c-basic-offset: 4
// End:
