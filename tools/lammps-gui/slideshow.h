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

#ifndef SLIDESHOW_H
#define SLIDESHOW_H

#include <QDialog>
#include <QImage>
#include <QString>
#include <QStringList>

class QDialogButtonBox;
class QLabel;
class QTimer;

class SlideShow : public QDialog {
    Q_OBJECT

public:
    explicit SlideShow(const QString &fileName, QWidget *parent = nullptr);
    void add_image(const QString &filename);
    void clear();

private slots:
    void first();
    void last();
    void next();
    void prev();
    void play();
    void loop();
    void zoomIn();
    void zoomOut();
    void normalSize();

private:
    void scaleImage(double factor);
    void loadImage(int idx);

private:
    QImage image;
    QTimer *playtimer;
    QLabel *imageLabel, *imageName;
    QDialogButtonBox *buttonBox;
    double scaleFactor = 1.0;

    int current;
    int maxwidth, maxheight;
    bool do_loop;
    QStringList imagefiles;
};
#endif

// Local Variables:
// c-basic-offset: 4
// End:
