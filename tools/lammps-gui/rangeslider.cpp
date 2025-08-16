/** Copyright Hoyoung Lee (2024-07-14).
 **
 ** hoyoung.yi@gmail.com
 **
 ** This software is a computer program whose purpose is to provide
 ** "Qt widget of range slider with two handles".
 **
 ** This software is governed by the CeCILL-A license under French law
 ** and abiding by the rules of distribution of free software. You can
 ** use, modify and/or redistribute the software under the terms of the
 ** CeCILL-A license as circulated by CEA, CNRS and INRIA at the following
 ** URL: "http://www.cecill.info".
 **
 ** As a counterpart to the access to the source code and rights to copy,
 ** modify and redistribute granted by the license, users are provided
 ** only with a limited warranty and the software's author, the holder of
 ** the economic rights, and the successive licensors have only limited
 ** liability.
 **
 ** In this respect, the user's attention is drawn to the risks associated
 ** with loading, using, modifying and/or developing or reproducing the
 ** software by the user in light of its specific status of free software,
 ** that may mean that it is complicated to manipulate, and that also
 ** therefore means that it is reserved for developers and experienced
 ** professionals having in-depth computer knowledge. Users are therefore
 ** encouraged to load and test the software's suitability as regards
 ** their requirements in conditions enabling the security of their
 ** systems and/or data to be ensured and, more generally, to use and
 ** operate it in the same conditions as regards security.
 **
 ** The fact that you are presently reading this means that you have had
 ** knowledge of the CeCILL-A license and that you accept its terms.
 */

#include "rangeslider.h"

#include <QColor>
#include <QMouseEvent>
#include <QPaintEvent>
#include <QPainter>
#include <QRect>
#include <QStyle>
#include <QStyleOptionSlider>

/*
 * Originated from
 * https://www.mail-archive.com/pyqt@riverbankcomputing.com/msg22889.html
 * Modification refered from
 * https://gist.github.com/Riateche/27e36977f7d5ea72cf4f
 */
RangeSlider::RangeSlider(Qt::Orientation ot, QWidget *parent) : QSlider(ot, parent)
{
    lowLimit        = minimum();
    highLimit       = maximum();
    pressed_control = QStyle::SC_None;
    tick_interval   = 0;
    tick_position   = QSlider::NoTicks;
    hover_control   = QStyle::SC_None;
    click_offset    = 0;
    active_slider   = 0;
}

// based on http://qt.gitorious.org/qt/qt/blobs/master/src/gui/widgets/qslider.cpp
void RangeSlider::paintEvent(QPaintEvent *)
{
    QPainter painter(this);
    QStyleOptionSlider opt;

    // Draw groove
    initStyleOption(&opt);
    opt.sliderValue    = 0;
    opt.sliderPosition = 0;
    opt.subControls    = QStyle::SC_SliderGroove;
    if (tickPosition() != NoTicks) opt.subControls |= QStyle::SC_SliderTickmarks;
    style()->drawComplexControl(QStyle::CC_Slider, &opt, &painter, this);
    QRect groove = style()->subControlRect(QStyle::CC_Slider, &opt, QStyle::SC_SliderGroove, this);

    // Draw span
    initStyleOption(&opt);
    opt.subControls    = QStyle::SC_SliderGroove;
    opt.sliderValue    = 0;
    opt.sliderPosition = lowLimit;
    QRect low_rect =
        style()->subControlRect(QStyle::CC_Slider, &opt, QStyle::SC_SliderHandle, this);
    opt.sliderPosition = highLimit;
    QRect high_rect =
        style()->subControlRect(QStyle::CC_Slider, &opt, QStyle::SC_SliderHandle, this);

    int low_pos  = pick(low_rect.center());
    int high_pos = pick(high_rect.center());
    int min_pos  = std::min(low_pos, high_pos);
    int max_pos  = std::max(low_pos, high_pos);

    QPoint c = QRect(low_rect.center(), high_rect.center()).center();
    QRect span_rect;
    if (opt.orientation == Qt::Horizontal) {
        span_rect = QRect(QPoint(min_pos, c.y() - 2), QPoint(max_pos, c.y() + 1));
        groove.adjust(0, 0, -1, 0);
    } else {
        span_rect = QRect(QPoint(c.x() - 2, min_pos), QPoint(c.x() + 1, max_pos));
        groove.adjust(0, 0, 0, -1);
    }

    QColor highlight = palette().color(QPalette::Highlight);
    painter.setBrush(QBrush(highlight));
    painter.setPen(QPen(highlight, 0));
    painter.drawRect(span_rect.intersected(groove));

    initStyleOption(&opt);
    opt.subControls = QStyle::SC_SliderHandle; // | QStyle::SC_SliderGroove;
    if (tickPosition() != QSlider::NoTicks) opt.subControls |= QStyle::SC_SliderTickmarks;

    if (pressed_control)
        opt.activeSubControls = pressed_control;
    else
        opt.activeSubControls = hover_control;

    opt.sliderPosition = lowLimit;
    opt.sliderValue    = lowLimit;
    style()->drawComplexControl(QStyle::CC_Slider, &opt, &painter, this);

    initStyleOption(&opt);
    opt.subControls = QStyle::SC_SliderHandle;
    if (tickPosition() != QSlider::NoTicks) opt.subControls |= QStyle::SC_SliderTickmarks;

    if (pressed_control)
        opt.activeSubControls = pressed_control;
    else
        opt.activeSubControls = hover_control;

    opt.sliderPosition = highLimit;
    opt.sliderValue    = highLimit;
    style()->drawComplexControl(QStyle::CC_Slider, &opt, &painter, this);
}

void RangeSlider::mousePressEvent(QMouseEvent *ev)
{
    /*  # In a normal slider control, when the user clicks on a point in the
        # slider's total range, but not on the slider part of the control the
        # control would jump the slider value to where the user clicked.
        # For this control, clicks which are not direct hits will slide both
        # slider parts
        */
    if (ev->button() == Qt::LeftButton) {
        ev->accept();
        QStyleOptionSlider opt;
        initStyleOption(&opt);
        active_slider = -1;

        QStyle::SubControl hit;
        opt.sliderPosition = lowLimit;
        hit = style()->hitTestComplexControl(QStyle::CC_Slider, &opt, ev->pos(), this);
        if (hit == QStyle::SC_SliderHandle) {
            active_slider   = 0;
            pressed_control = hit;
            triggerAction(SliderMove);
            setRepeatAction(SliderNoAction);
            setSliderDown(true);
        } else {
            opt.sliderPosition = highLimit;
            hit = style()->hitTestComplexControl(QStyle::CC_Slider, &opt, ev->pos(), this);
            if (hit == QStyle::SC_SliderHandle) {
                active_slider   = 1;
                pressed_control = hit;
                triggerAction(SliderMove);
                setRepeatAction(SliderNoAction);
                setSliderDown(true);
            }
        }

        if (active_slider < 0) {
            pressed_control = QStyle::SC_SliderHandle;
            click_offset    = pixelPosToRangeValue(pick(ev->pos()));
            triggerAction(SliderMove);
            setRepeatAction(SliderNoAction);
        }
    } else {
        ev->ignore();
    }

    QSlider::mousePressEvent(ev);
}

void RangeSlider::mouseMoveEvent(QMouseEvent *ev)
{
    if (pressed_control != QStyle::SC_SliderHandle) {
        ev->ignore();
        return;
    }

    ev->accept();
    int new_pos = pixelPosToRangeValue(pick(ev->pos()));
    QStyleOptionSlider opt;
    initStyleOption(&opt);

    int offset, diff;
    if (active_slider < 0) {
        offset = new_pos - click_offset;
        highLimit += offset;
        lowLimit += offset;
        if (lowLimit < minimum()) {
            diff = minimum() - lowLimit;
            lowLimit += diff;
            highLimit += diff;
        }
        if (highLimit > maximum()) {
            diff = maximum() - highLimit;
            lowLimit += diff;
            highLimit += diff;
        }
    } else if (active_slider == 0) {
        if (new_pos >= highLimit) new_pos = highLimit - 1;
        lowLimit = new_pos;
    } else {
        if (new_pos <= lowLimit) new_pos = lowLimit + 1;
        highLimit = new_pos;
    }

    click_offset = new_pos;
    update();
    emit sliderMoved(lowLimit, highLimit);
}

int RangeSlider::pixelPosToRangeValue(int pos)
{
    QStyleOptionSlider opt;
    initStyleOption(&opt);

    QRect gr = style()->subControlRect(QStyle::CC_Slider, &opt, QStyle::SC_SliderGroove, this);
    QRect sr = style()->subControlRect(QStyle::CC_Slider, &opt, QStyle::SC_SliderHandle, this);

    int slider_length, slider_min, slider_max;
    if (orientation() == Qt::Horizontal) {
        slider_length = sr.width();
        slider_min    = gr.x();
        slider_max    = gr.right() - slider_length + 1;
    } else {
        slider_length = sr.height();
        slider_min    = gr.y();
        slider_max    = gr.bottom() - slider_length + 1;
    }

    return style()->sliderValueFromPosition(minimum(), maximum(), pos - slider_min,
                                            slider_max - slider_min, opt.upsideDown);
}
