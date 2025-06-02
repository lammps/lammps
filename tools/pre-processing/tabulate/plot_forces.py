#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   https://www.lammps.org/ Sandia National Laboratories
#   LAMMPS Development team: developers@lammps.org
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------
# Author:  Germain Clavier (Unicaen), germain.clavier at unicaen.fr

"""
Plot LAMMPS tabulated forces.
"""

import argparse
import numpy as np
import os
import logging
import sys
from matplotlib import pyplot as plt

logger = logging.getLogger(__name__)


units = {
    "lj": {
        "Distance": "Reduced units",
        "Energy": "Reduced units",
        "Force": "Reduced units",
        "kb": 1,
    },
    "real": {
        "Distance": "[A]",
        "Energy": "[kcal/mol]",
        "Force": "[kcal/mol/A]",
        "kb": 0.001985875,
    },
    "metal": {
        "Distance": "[A]",
        "Energy": "[eV]",
        "Force": "[eV/A]",
        "kb": 8.6173332e-5,
    },
    "si": {"Distance": "[m]", "Energy": "[J]", "Force": "[N]", "kb": 1.380649e-23},
    "cgs": {
        "Distance": "[cm]",
        "Energy": "[ergs]",
        "Force": "[dynes]",
        "kb": 1.3806504e-16,
    },
    "electron": {
        "Distance": "[Bohr]",
        "Energy": "[Hartrees]",
        "Force": "[Hartree/Bohr]",
        "kb": 3.16681534e-6,
    },
    "micro": {
        "Distance": "[um]",
        "Energy": "[pg·um^2/us^2]",
        "Force": "[pg·um/us^2]",
        "kb": 1.3806504e-8,
    },
    "nano": {
        "Distance": "[nm]",
        "Energy": "[ag·nm^2/ns^2]",
        "Force": "[ag·nm/ns^2]",
        "kb": 0.013806504,
    },
}


def compute_energy(tp):
    r = tp[0]
    fo = tp[2]
    e = np.zeros(r.shape)
    for i, (ri, fi) in enumerate(zip(r, fo)):
        if i == 0:
            continue
        dr = ri - r[i - 1]
        e[i] = e[i - 1] - dr * fo[i - 1]
    e -= e[-1]
    return e


def main():

    parser = argparse.ArgumentParser(
        description="""
        Plots LAMMPS tabulated forces. This script takes a table
        file as an input and plots all the tabulated forces inside with their
        corresponding energy. The forces label is the token used to name the
        force in the file. It can be used to output all the forces in separate
        files and/or recompute the energy from forces through finite difference
        (assuming e(rc)=0). This script requires the matplotlib and numpy
        Python libraries. Bitmap format is not supported.
        """
    )
    parser.add_argument(
        "-u",
        "--units",
        dest="units",
        default="real",
        help="Units of the file (LAMMPS units system)",
    )
    parser.add_argument(
        "-f",
        "--file",
        dest="infile",
        default="",
        help="File to read",
    )
    parser.add_argument(
        "-x",
        dest="xrange",
        default="",
        help="xrange separated by : (for negative values use the '=' sign: -x=-3:10)",
    )
    parser.add_argument(
        "-y",
        dest="yrange",
        default="",
        help="yrange separated by :",
    )
    parser.add_argument(
        "-t",
        dest="temp",
        default=None,
        type=float,
        help="temperature for KbT plot [default none]",
    )
    parser.add_argument(
        "-d",
        "--diff-num",
        dest="recompute",
        action="store_true",
        help="Recompute the energies from forces and distances through finite differences",
    )
    parser.add_argument(
        "-e",
        dest="extract",
        action="store_true",
        help="Extract the forces in separate files",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    ##########
    # Manage arguments

    udistance = units[args.units]["Distance"]
    uenergy = units[args.units]["Energy"]
    uforce = units[args.units]["Force"]
    kb = units[args.units]["kb"]
    rlabel = " ".join(["Rij", udistance])
    elabel = " ".join(["E", uenergy])
    flabel = " ".join(["F", uforce])
    etitle = "Energy"
    ftitle = "Force"
    font = "DejaVu Sans"
    fontsize = 30

    infile = args.infile
    if not os.path.isfile(infile):
        logger.error("Input file not found")
        sys.exit(1)

    toplot = []
    with open(infile, "r") as f:
        lines = iter(f.readlines())
        while True:
            try:
                r = []
                force = []
                ener = []
                tok = []
                while not tok:
                    tok = next(lines).partition("#")[0].rstrip()
                logger.info("Found {} token".format(tok))
                infos = next(lines).split()
                npoints = int(infos[1])
                next(lines)
                if "bitmap" in infos:
                    logger.info("Unsupported bitmap format for token {:s}".format(tok))
                    for _ in range(npoints):
                        continue
                else:
                    for i in range(npoints):
                        line = next(lines).split()
                        r.append(float(line[1]))
                        ener.append(float(line[2]))
                        force.append(float(line[3]))
                    r = np.array(r)
                    ener = np.array(ener)
                    force = np.array(force)
                    toplot.append([r, ener, force, tok])
                tok = []
                next(lines)
            except StopIteration:
                break
    if args.recompute:
        etitle = "Estimated energy"
        for tp in toplot:
            tp[1] = compute_energy(tp)

    fig, axes = plt.subplots(1, 2)

    for tp in toplot:
        axes[0].plot(tp[0], tp[1], label=tp[3], linewidth=3)
        axes[1].plot(tp[0], tp[2], label=tp[3], linewidth=3)
        hmin, hmax = axes[1].get_xlim()
        axes[1].hlines(0, hmin, hmax, color="black", linewidth=3, linestyles="dashdot")

    if args.temp:
        if args.temp > 0:
            hmin, hmax = axes[0].get_xlim()
            axes[0].hlines(
                kb * args.temp,
                hmin,
                hmax,
                color="orange",
                label=r"$k_BT$",
                linewidth=3,
                linestyles="dashdot",
            )
            axes[0].text(hmax / 2.0, kb * args.temp, "KbT", fontsize=0.7 * fontsize)
            logger.info("KbT value= {:e} {:s}".format(kb * args.temp, uenergy))
        else:
            logger.info("Invalid temperature value: {:e}".format(args.temp))

    if args.xrange:
        xmin, xmax = list(map(float, args.xrange.split(":")))
        axes[0].set_xlim(xmin, xmax)
        axes[1].set_xlim(xmin, xmax)
    if args.yrange:
        ymin, ymax = list(map(float, args.yrange.split(":")))
        axes[0].set_ylim(ymin, ymax)
        axes[1].set_ylim(ymin, ymax)

    # Setting axes 0
    axes[0].set_title(etitle, fontsize=fontsize)
    axes[0].set_xlabel(
        rlabel, fontname=font, fontsize=fontsize
    )  # xlabel name, size 30pts
    axes[0].set_ylabel(
        elabel, fontname=font, fontsize=fontsize
    )  # ylabel name, size 30pts
    axes[0].tick_params(
        axis="both", which="major", labelsize=fontsize
    )  # Biggers ticks, bigger tick labels!

    # Setting axes 1
    axes[1].set_title(ftitle, fontsize=fontsize)
    axes[1].legend(frameon=False, fontsize=fontsize)  # Fat font, no frame
    axes[1].set_xlabel(
        rlabel, fontname=font, fontsize=fontsize
    )
    axes[1].set_ylabel(
        flabel, fontname=font, fontsize=fontsize
    )
    axes[1].tick_params(
        axis="both", which="major", labelsize=fontsize
    )

    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.show()

    if args.extract:
        for tp in toplot:
            outfile = "".join([tp[3], ".plot"])
            logger.info("Writing file {}".format(outfile))
            with open(outfile, "w") as f:
                f.write("# {} force extracted from {}\n".format(tp[3], infile))
                f.write("# {:^20} {:^20} {:^20}\n".format("r", "energy", "force"))
                for a, b, c in zip(tp[0], tp[1], tp[2]):
                    f.write("{:>18.16e} {:>18.16e} {:>18.16e}\n".format(a, b, c))
    return


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        raise SystemExit("User interruption.")
