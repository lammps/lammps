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

#include <QDialogButtonBox>
#include <QFileInfo>
#include <QGuiApplication>
#include <QHBoxLayout>
#include <QImage>
#include <QImageReader>
#include <QLabel>
#include <QPalette>
#include <QPushButton>
#include <QScreen>
#include <QSettings>
#include <QShortcut>
#include <QSpacerItem>
#include <QTimer>
#include <QVBoxLayout>

SlideShow::SlideShow(const QString &fileName, QWidget *parent) :
    QDialog(parent), playtimer(nullptr), imageLabel(new QLabel), imageName(new QLabel("(none)")),
    do_loop(true)
{
    imageLabel->setBackgroundRole(QPalette::Base);
    imageLabel->setScaledContents(false);
    imageLabel->setMinimumSize(100, 100);

    imageName->setFrameStyle(QFrame::Raised);
    imageName->setFrameShape(QFrame::Panel);
    imageName->setAlignment(Qt::AlignCenter);
    imageName->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    auto *shortcut = new QShortcut(QKeySequence(Qt::CTRL + Qt::Key_W), this);
    QObject::connect(shortcut, &QShortcut::activated, this, &SlideShow::close);

    buttonBox = new QDialogButtonBox(QDialogButtonBox::Close);

    connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);

    auto *mainLayout = new QVBoxLayout;
    auto *navLayout  = new QHBoxLayout;
    auto *botLayout  = new QHBoxLayout;

    // workaround for incorrect highlight bug on macOS
    auto *dummy = new QPushButton(QIcon(), "");
    dummy->hide();

    auto *gofirst = new QPushButton(QIcon(":/go-first.png"), "");
    gofirst->setToolTip("Go to first Image");
    auto *goprev = new QPushButton(QIcon(":/go-previous-2.png"), "");
    goprev->setToolTip("Go to previous Image");
    auto *goplay = new QPushButton(QIcon(":/media-playback-start-2.png"), "");
    goplay->setToolTip("Play animation");
    goplay->setCheckable(true);
    goplay->setChecked(playtimer);
    goplay->setObjectName("play");
    auto *gonext = new QPushButton(QIcon(":/go-next-2.png"), "");
    gonext->setToolTip("Go to next Image");
    auto *golast = new QPushButton(QIcon(":/go-last.png"), "");
    golast->setToolTip("Go to last Image");
    auto *goloop = new QPushButton(QIcon(":/media-playlist-repeat.png"), "");
    goloop->setToolTip("Loop animation");
    goloop->setCheckable(true);
    goloop->setChecked(do_loop);

    auto *zoomin = new QPushButton(QIcon(":/gtk-zoom-in.png"), "");
    zoomin->setToolTip("Zoom in by 10 percent");
    auto *zoomout = new QPushButton(QIcon(":/gtk-zoom-out.png"), "");
    zoomout->setToolTip("Zoom out by 10 percent");
    auto *normal = new QPushButton(QIcon(":/gtk-zoom-fit.png"), "");
    normal->setToolTip("Reset zoom to normal");

    connect(gofirst, &QPushButton::released, this, &SlideShow::first);
    connect(goprev, &QPushButton::released, this, &SlideShow::prev);
    connect(goplay, &QPushButton::released, this, &SlideShow::play);
    connect(gonext, &QPushButton::released, this, &SlideShow::next);
    connect(golast, &QPushButton::released, this, &SlideShow::last);
    connect(goloop, &QPushButton::released, this, &SlideShow::loop);
    connect(zoomin, &QPushButton::released, this, &SlideShow::zoomIn);
    connect(zoomout, &QPushButton::released, this, &SlideShow::zoomOut);
    connect(gofirst, &QPushButton::released, this, &SlideShow::first);
    connect(normal, &QPushButton::released, this, &SlideShow::normalSize);

    navLayout->addSpacerItem(new QSpacerItem(10, 10, QSizePolicy::Expanding, QSizePolicy::Minimum));
    navLayout->addWidget(dummy);
    navLayout->addWidget(gofirst);
    navLayout->addWidget(goprev);
    navLayout->addWidget(goplay);
    navLayout->addWidget(gonext);
    navLayout->addWidget(golast);
    navLayout->addWidget(goloop);

    navLayout->addWidget(zoomin);
    navLayout->addWidget(zoomout);
    navLayout->addWidget(normal);

    mainLayout->addWidget(imageLabel);
    mainLayout->addLayout(navLayout);

    botLayout->addWidget(imageName);
    botLayout->addWidget(buttonBox);
    botLayout->setStretch(0, 3);
    mainLayout->addLayout(botLayout);

    setWindowIcon(QIcon(":/lammps-icon-128x128.png"));
    setWindowTitle(QString("LAMMPS-GUI - Slide Show: ") + QFileInfo(fileName).fileName());

    imagefiles.clear();
    scaleFactor = 1.0;
    current     = 0;

    auto maxsize = QGuiApplication::primaryScreen()->availableSize() * 4 / 5;
    maxheight    = maxsize.height();
    maxwidth     = maxsize.width();

    setLayout(mainLayout);
}

void SlideShow::add_image(const QString &filename)
{
    if (!imagefiles.contains(filename)) {
        int lastidx = imagefiles.size();
        imagefiles.append(filename);
        loadImage(lastidx);
    }
}

void SlideShow::clear()
{
    imagefiles.clear();
    image.fill(Qt::black);
    imageLabel->setPixmap(QPixmap::fromImage(image));
    imageLabel->adjustSize();
    imageName->setText("(none)");
    repaint();
}

void SlideShow::loadImage(int idx)
{
    if ((idx < 0) || (idx >= imagefiles.size())) return;

    do {
        QImageReader reader(imagefiles[idx]);
        reader.setAutoTransform(true);
        const QImage newImage = reader.read();

        // There was an error reading the image file. Try reading the previous image instead.
        if (newImage.isNull()) {
            --idx;
        } else {
            int newheight = (int)newImage.height() * scaleFactor;
            int newwidth  = (int)newImage.width() * scaleFactor;
            image         = newImage.scaled(newwidth, newheight, Qt::IgnoreAspectRatio,
                                            Qt::SmoothTransformation);
            imageLabel->setPixmap(QPixmap::fromImage(image));
            imageLabel->setMinimumSize(newwidth, newheight);
            imageName->setText(QString(" Image %1 / %2 : %3 ")
                                   .arg(idx + 1)
                                   .arg(imagefiles.size())
                                   .arg(imagefiles[idx]));
            adjustSize();
            current = idx;
            break;
        }
    } while (idx >= 0);
}

void SlideShow::first()
{
    current = 0;
    loadImage(current);
}

void SlideShow::last()
{
    current = imagefiles.size() - 1;
    loadImage(current);
}

void SlideShow::play()
{
    // if we do not loop, start animation from beginning
    if (!do_loop) current = 0;

    if (playtimer) {
        playtimer->stop();
        delete playtimer;
        playtimer = nullptr;
    } else {
        playtimer = new QTimer(this);
        connect(playtimer, &QTimer::timeout, this, &SlideShow::next);
        playtimer->start(100);
    }

    // reset push button state. use findChild() if not triggered from button.
    QPushButton *button = qobject_cast<QPushButton *>(sender());
    if (!button) button = findChild<QPushButton *>("play");
    if (button) button->setChecked(playtimer);
}

void SlideShow::next()
{
    ++current;
    if (current >= imagefiles.size()) {
        if (do_loop) {
            current = 0;
        } else {
            // stop animation
            if (playtimer) play();
            --current;
        }
    }
    loadImage(current);
}

void SlideShow::prev()
{
    --current;
    if (current < 0) {
        if (do_loop)
            current = imagefiles.size() - 1;
        else
            current = 0;
    }
    loadImage(current);
}

void SlideShow::loop()
{
    QPushButton *button = qobject_cast<QPushButton *>(sender());
    do_loop             = !do_loop;
    button->setChecked(do_loop);
}

void SlideShow::zoomIn()
{
    scaleImage(1.1);
}

void SlideShow::zoomOut()
{
    scaleImage(0.9);
}

void SlideShow::normalSize()
{
    scaleFactor = 1.0;
    scaleImage(1.0);
}

void SlideShow::scaleImage(double factor)
{
    // compute maxfactor so the image is not scaled beyond 80 of width or height of screen
    double maxfactor = 10.0;
    maxfactor        = qMin((double)maxheight / (double)image.height(), maxfactor);
    maxfactor        = qMin((double)maxwidth / (double)image.width(), maxfactor);

    if (factor > maxfactor) factor = maxfactor;
    scaleFactor *= factor;
    if (scaleFactor < 0.25) scaleFactor = 0.25;

    loadImage(current);
}

// Local Variables:
// c-basic-offset: 4
// End:
