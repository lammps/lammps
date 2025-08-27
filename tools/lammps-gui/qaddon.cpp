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

#include "qaddon.h"

#include <QString>
#include <QStringList>
#include <QWidget>

namespace {
// clang-format off
QStringList imagecolors = {
    "aliceblue", "antiquewhite", "aqua", "aquamarine", "azure", "beige", "bisque", "black",
    "blanchedalmond", "blue", "blueviolet", "brown", "burlywood", "cadetblue", "chartreuse",
    "chocolate", "coral", "cornflowerblue", "cornsilk", "crimson", "cyan", "darkblue", "darkcyan",
    "darkgoldenrod", "darkgray", "darkgreen", "darkkhaki", "darkmagenta", "darkolivegreen",
    "darkorange", "darkorchid", "darkred", "darksalmon", "darkseagreen", "darkslateblue",
    "darkslategray", "darkturquoise", "darkviolet", "deeppink", "deepskyblue", "dimgray",
    "dodgerblue", "firebrick", "floralwhite", "forestgreen", "fuchsia", "gainsboro", "ghostwhite",
    "gold", "goldenrod", "gray", "green", "greenyellow", "honeydew", "hotpink", "indianred",
    "indigo", "ivory", "khaki", "lavender", "lavenderblush", "lawngreen", "lemonchiffon",
    "lightblue", "lightcoral", "lightcyan", "lightgoldenrodyellow", "lightgreen", "lightgrey",
    "lightpink", "lightsalmon", "lightseagreen", "lightskyblue", "lightslategray", "lightsteelblue",
    "lightyellow", "lime", "limegreen", "linen", "magenta", "maroon", "mediumaquamarine",
    "mediumblue", "mediumorchid", "mediumpurple", "mediumseagreen", "mediumslateblue",
    "mediumspringgreen", "mediumturquoise", "mediumvioletred", "midnightblue", "mintcream",
    "mistyrose", "moccasin", "navajowhite", "navy", "oldlace", "olive", "olivedrab", "orange",
    "orangered", "orchid", "palegoldenrod", "palegreen", "paleturquoise", "palevioletred",
    "papayawhip", "peachpuff", "peru", "pink", "plum", "powderblue", "purple", "red", "rosybrown",
    "royalblue", "saddlebrown", "salmon", "sandybrown", "seagreen", "seashell", "sienna", "silver",
    "skyblue", "slateblue", "slategray", "snow", "springgreen", "steelblue", "tan", "teal",
    "thistle", "tomato", "turquoise", "violet", "wheat", "white", "whitesmoke", "yellow",
    "yellowgreen"
};
// clang-format on
} // namespace

QHline::QHline(QWidget *parent) : QFrame(parent)
{
    setGeometry(QRect(0, 0, 100, 3));
    setFrameShape(QFrame::HLine);
    setFrameShadow(QFrame::Sunken);
}

QColorValidator::QColorValidator(QWidget *parent) : QValidator(parent) {}

void QColorValidator::fixup(QString &input) const
{
    // remove leading/trailing whitespace and make lowercase
    input = input.trimmed();
    input = input.toLower();
}

QValidator::State QColorValidator::validate(QString &input, int &pos) const
{
    QString match;

    // find if input string is contained in list of colors
    for (auto color : imagecolors) {
        if (color.startsWith(input)) {
            match = color;
            break;
        }
    }

    if (match == input) {
        return QValidator::Acceptable;
    } else if (match.size() > 0) {
        return QValidator::Intermediate;
    }
    return QValidator::Invalid;
}

// complete color inputs
QColorCompleter::QColorCompleter(QWidget *parent) : QCompleter(imagecolors, parent)
{
    setCompletionMode(QCompleter::InlineCompletion);
    setModelSorting(QCompleter::CaseInsensitivelySortedModel);
};

// Local Variables:
// c-basic-offset: 4
// End:
