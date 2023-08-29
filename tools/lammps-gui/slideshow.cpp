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

#include "slideshow.h"
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
#include <QWidgetAction>

#include <cmath>

static const QString blank(" ");

SlideShow::SlideShow(const QString &fileName, LammpsWrapper *_lammps, QWidget *parent) :
    QDialog(parent), imageLabel(new QLabel), scrollArea(new QScrollArea), menuBar(new QMenuBar),
    lammps(_lammps)
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

    mainLayout->addWidget(scrollArea);
    mainLayout->addWidget(buttonBox);
    setWindowTitle(QString("Slide Show: ") + QFileInfo(fileName).fileName());
    createActions();

    scaleFactor = 1.0;
    resize(image.width() + 20, image.height() + 50);

    scrollArea->setVisible(true);
    fitToWindowAct->setEnabled(true);
    updateActions();
    if (!fitToWindowAct->isChecked()) imageLabel->adjustSize();
    setLayout(mainLayout);
}

void SlideShow::saveAs()
{
    QString fileName = QFileDialog::getSaveFileName(this, "Save Image File As", QString(),
                                                    "Image Files (*.jpg *.png *.bmp *.ppm)");
    saveFile(fileName);
}

void SlideShow::copy() {}

void SlideShow::zoomIn()
{
    scaleImage(1.25);
}

void SlideShow::zoomOut()
{
    scaleImage(0.8);
}

void SlideShow::normalSize()
{
    imageLabel->adjustSize();
    scaleFactor = 1.0;
}

void SlideShow::fitToWindow()
{
    bool fitToWindow = fitToWindowAct->isChecked();
    scrollArea->setWidgetResizable(fitToWindow);
    if (!fitToWindow) normalSize();
    updateActions();
}

void SlideShow::saveFile(const QString &fileName)
{
    if (!fileName.isEmpty()) image.save(fileName);
}

void SlideShow::createActions()
{
    QMenu *fileMenu = menuBar->addMenu("&File");

    saveAsAct = fileMenu->addAction("&Save As...", this, &SlideShow::saveAs);
    saveAsAct->setIcon(QIcon(":/document-save-as.png"));
    saveAsAct->setEnabled(false);
    fileMenu->addSeparator();
    copyAct = fileMenu->addAction("&Copy", this, &SlideShow::copy);
    copyAct->setIcon(QIcon(":/edit-copy.png"));
    copyAct->setShortcut(QKeySequence::Copy);
    copyAct->setEnabled(false);
    fileMenu->addSeparator();
    QAction *exitAct = fileMenu->addAction("&Close", this, &QWidget::close);
    exitAct->setIcon(QIcon(":/window-close.png"));
    exitAct->setShortcut(QKeySequence::fromString("Ctrl+W"));

    QMenu *viewMenu = menuBar->addMenu("&View");

    zoomInAct = viewMenu->addAction("Image Zoom &In (25%)", this, &SlideShow::zoomIn);
    zoomInAct->setShortcut(QKeySequence::ZoomIn);
    zoomInAct->setIcon(QIcon(":/gtk-zoom-in.png"));
    zoomInAct->setEnabled(false);

    zoomOutAct = viewMenu->addAction("Image Zoom &Out (25%)", this, &SlideShow::zoomOut);
    zoomOutAct->setShortcut(QKeySequence::ZoomOut);
    zoomOutAct->setIcon(QIcon(":/gtk-zoom-out.png"));
    zoomOutAct->setEnabled(false);

    normalSizeAct = viewMenu->addAction("&Reset Image Size", this, &SlideShow::normalSize);
    normalSizeAct->setShortcut(QKeySequence::fromString("Ctrl+0"));
    normalSizeAct->setIcon(QIcon(":/gtk-zoom-fit.png"));
    normalSizeAct->setEnabled(false);

    viewMenu->addSeparator();

    fitToWindowAct = viewMenu->addAction("&Fit to Window", this, &SlideShow::fitToWindow);
    fitToWindowAct->setEnabled(false);
    fitToWindowAct->setCheckable(true);
    fitToWindowAct->setShortcut(QKeySequence::fromString("Ctrl+="));
}

void SlideShow::updateActions()
{
    saveAsAct->setEnabled(!image.isNull());
    copyAct->setEnabled(!image.isNull());
    zoomInAct->setEnabled(!fitToWindowAct->isChecked());
    zoomOutAct->setEnabled(!fitToWindowAct->isChecked());
    normalSizeAct->setEnabled(!fitToWindowAct->isChecked());
}

void SlideShow::scaleImage(double factor)
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

void SlideShow::adjustScrollBar(QScrollBar *scrollBar, double factor)
{
    scrollBar->setValue(
        int(factor * scrollBar->value() + ((factor - 1) * scrollBar->pageStep() / 2)));
}

// Local Variables:
// c-basic-offset: 4
// End:
