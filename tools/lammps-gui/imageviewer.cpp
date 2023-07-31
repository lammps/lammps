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
#include <QScreen>
#include <QScrollArea>
#include <QScrollBar>
#include <QStatusBar>
#include <QVBoxLayout>
#include <QWheelEvent>

ImageViewer::ImageViewer(const QString &fileName, QWidget *parent) :
    QDialog(parent), imageLabel(new QLabel), scrollArea(new QScrollArea), menuBar(new QMenuBar)
{
    imageLabel->setBackgroundRole(QPalette::Base);
    imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    imageLabel->setScaledContents(true);
    imageLabel->minimumSizeHint();

    scrollArea->setBackgroundRole(QPalette::Dark);
    scrollArea->setWidget(imageLabel);
    scrollArea->setMouseTracking(true);
    scrollArea->installEventFilter(this);
    scrollArea->setVisible(false);

    buttonBox = new QDialogButtonBox(QDialogButtonBox::Close);

    connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(menuBar);
    mainLayout->addWidget(scrollArea);
    mainLayout->addWidget(buttonBox);
    setWindowTitle(QString("Image Viewer: ") + QFileInfo(fileName).completeBaseName());

    createActions();

    QImageReader reader(fileName);
    reader.setAutoTransform(true);
    const QImage newImage = reader.read();
    if (newImage.isNull()) {
        QMessageBox::warning(this, QGuiApplication::applicationDisplayName(),
                             tr("Cannot load %1: %2").arg(fileName, reader.errorString()));
        return;
    }
    image = newImage;
    imageLabel->setPixmap(QPixmap::fromImage(image));
    scaleFactor = 1.0;
    resize(image.width() + 20, image.height() + 50);

    scrollArea->setVisible(true);
    fitToWindowAct->setEnabled(true);
    updateActions();
    if (!fitToWindowAct->isChecked()) imageLabel->adjustSize();
    setLayout(mainLayout);
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
    scaleImage(1.1);
}

void ImageViewer::zoomOut()
{
    scaleImage(0.9);
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
    QMenu *fileMenu = menuBar->addMenu(tr("&File"));

    saveAsAct = fileMenu->addAction(tr("&Save As..."), this, &ImageViewer::saveAs);
    saveAsAct->setEnabled(false);
    fileMenu->addSeparator();
    copyAct = fileMenu->addAction(tr("&Copy"), this, &ImageViewer::copy);
    copyAct->setShortcut(QKeySequence::Copy);
    copyAct->setEnabled(false);
    fileMenu->addSeparator();
    QAction *exitAct = fileMenu->addAction(tr("&Close"), this, &QWidget::close);
    exitAct->setShortcut(tr("Ctrl+W"));

    QMenu *viewMenu = menuBar->addMenu(tr("&View"));

    zoomInAct = viewMenu->addAction(tr("Zoom &In (10%)"), this, &ImageViewer::zoomIn);
    zoomInAct->setShortcut(QKeySequence::ZoomIn);
    zoomInAct->setEnabled(false);

    zoomOutAct = viewMenu->addAction(tr("Zoom &Out (10%)"), this, &ImageViewer::zoomOut);
    zoomOutAct->setShortcut(QKeySequence::ZoomOut);
    zoomOutAct->setEnabled(false);

    normalSizeAct = viewMenu->addAction(tr("&Normal Size"), this, &ImageViewer::normalSize);
    normalSizeAct->setShortcut(tr("Ctrl+S"));
    normalSizeAct->setEnabled(false);

    viewMenu->addSeparator();

    fitToWindowAct = viewMenu->addAction(tr("&Fit to Window"), this, &ImageViewer::fitToWindow);
    fitToWindowAct->setEnabled(false);
    fitToWindowAct->setCheckable(true);
    fitToWindowAct->setShortcut(tr("Ctrl+F"));
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

bool ImageViewer::eventFilter(QObject *, QEvent *event)
{
    if (event->type() == QEvent::Wheel) {
        wheelEvent((QWheelEvent *)event);
        return true;
    }
    return false;
}

void ImageViewer::wheelEvent(QWheelEvent *event)
{
    QPoint num = event->angleDelta();
    if (!num.isNull()) {
        if (num.y() > 0)
            zoomIn();
        else
            zoomOut();
    }
    event->accept();
}

// Local Variables:
// c-basic-offset: 4
// End:
