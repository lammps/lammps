/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include <QDialog>
#include <QImage>

class QAction;
class QMenuBar;
class QDialogButtonBox;
class QLabel;
class QObject;
class QScrollArea;
class QScrollBar;
class QStatusBar;
class QWheelEvent;

class ImageViewer : public QDialog {
    Q_OBJECT

public:
    explicit ImageViewer(const QString &fileName, QWidget *parent = nullptr);

private slots:
    void saveAs();
    void copy();
    void zoomIn();
    void zoomOut();
    void normalSize();
    void fitToWindow();

private:
    void createActions();
    void updateActions();
    void saveFile(const QString &fileName);
    void scaleImage(double factor);
    void adjustScrollBar(QScrollBar *scrollBar, double factor);
    bool eventFilter(QObject *object, QEvent *event);
    void wheelEvent(QWheelEvent *event);

private:
    QImage image;
    QMenuBar *menuBar;
    QLabel *imageLabel;
    QScrollArea *scrollArea;
    QDialogButtonBox *buttonBox;
    double scaleFactor = 1.0;

    QAction *saveAsAct;
    QAction *copyAct;
    QAction *zoomInAct;
    QAction *zoomOutAct;
    QAction *normalSizeAct;
    QAction *fitToWindowAct;
};
#endif

// Local Variables:
// c-basic-offset: 4
// End:
