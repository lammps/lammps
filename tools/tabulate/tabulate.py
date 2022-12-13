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
"""
Python classes and support functions to generate tabulated potential files for LAMMPS.
"""

# for python2/3 compatibility
from __future__ import print_function

import sys
import argparse

########################################################################

def numdiff(x, func):
    """ Get the value of the first derivative of a function 'func(x)' at 'x'
        from numerical differentiation."""

    # optimal delta x value for 2-point numerical differentiation of two floating point numbers
    epsilon = 6.05504e-6
    f1 = func(x - epsilon)
    f2 = func(x + epsilon)
    return 0.5 * (f1-f2) / epsilon

########################################################################

def mktable(tstyle, label, num, xmin, xmax, efunc, diff=False, ffunc=None):
    """ Do the tabulation of the provided energy function. Compute force from
    numerical differentiation if no force function is provided.  Also detect
    minimum for use to determine potential shifting in bonded potentials."""

    # must use numerical differentiation if no force function provided
    if not ffunc: diff = True

    print("# Creating %s table %s with %d points from %g to %g" % (tstyle, label, num, xmin, xmax))

    table = []
    delx = (xmax - xmin) / (float(num) - 1.0)
    emin = 999999999.0
    xzero = 0.0

    for i in range(0,num):
        x = xmin + i*delx
        e = efunc(x)
        if e < emin:
            emin = e
            xzero = x
        if diff:
            f = numdiff(x, efunc)
        else:
            f = ffunc(x)
        table.append([i+1, x, e, f])

    return table, xzero


########################################################################
# base class with shared functionality
class Tabulate(object):

    def __init__(self, style, efunc, ffunc=None):
        self.fp = sys.stdout
        self.tstyle = style
        self.energyfunc = efunc
        self.forcefunc = ffunc
        self.eshift = 0.0
        self.args = None
        self.diff = True

        self.parser = argparse.ArgumentParser(description='Tool to generate tabulated '
                                              + self.tstyle + ' potential files for LAMMPS')
        self.parser.add_argument('--num-points', '-n', dest='num', default=1000, type=int,
                                 help="Number of tabulated points")
        self.parser.add_argument('--filename', '-f', dest='filename', default='-',
                                 help="Name of output file")
        self.parser.add_argument('--diff-num', '-d', dest='diff',default=False,
                                 action='store_true',
                                 help="Differentiate energy function numerically")
        self.parser.add_argument('--inner', '-i', dest='xmin', required=True, type=float,
                                 help="Inner cutoff of table")
        self.parser.add_argument('--outer', '-o', dest='xmax', required=True, type=float,
                                 help="Outer cutoff of table")

    def openfile(self, label):
        if self.args and self.args.filename != '-':
            self.fp = open(self.args.filename, 'a', encoding='ascii')
        self.fp.write('\n' + label + '\n')

    def writetable(self, table, offset):
        for i,r,e,f in table:
            self.fp.write("%8d %- 22.15g %- 22.15g %- 22.15g\n" % (i, r, e - offset, f))

################################################################################
# create tabulation for pair styles
class PairTabulate(Tabulate):
    def __init__(self, efunc, ffunc=None):
        super(PairTabulate, self).__init__('pair', efunc, ffunc)
        self.parser.add_argument('--eshift', '-e', dest='eshift', default=False,
                                 action='store_true',
                                 help="Shift potential energy to be zero at outer cutoff")
        try:
            self.args = self.parser.parse_args()
        except argparse.ArgumentError:
            sys.exit()

    def run(self, label):
        # sanity checks
        if self.args.num < 2:
            sys.exit('\nExpect 2 or more points in table for tabulation\n\n' + self.parser.format_help())
        if self.args.xmin <= 0.0:
            sys.exit('\nInner tabulation cutoff must be > 0 for pair style table\n\n' + self.parser.format_help())
        if self.args.xmax <= self.args.xmin:
            sys.exit('\nOuter cutoff must be larger than inner cutoff\n\n' + self.parser.format_help())

        self.diff = self.args.diff
        if not self.forcefunc: self.diff = True
        offset = 0.0
        if self.args.eshift: offset=self.energyfunc(self.args.xmax)

        (table, dummy) = mktable(self.tstyle, label, self.args.num, self.args.xmin, self.args.xmax,
                                 self.energyfunc, self.args.diff, self.forcefunc)

        # open table file
        self.openfile(label)

        # write pair style specific header
        if self.forcefunc:
            diffmin = numdiff(self.args.xmin, self.forcefunc)
            diffmax = numdiff(self.args.xmax, self.forcefunc)
            self.fp.write("N %d R %g %g FPRIME %- 22.15g %- 22.15g\n\n"
                          % (self.args.num, self.args.xmin, self.args.xmax, diffmin, diffmax))
        else:
            self.fp.write("N %d R %g %g\n\n" % (self.args.num, self.args.xmin, self.args.xmax))

        self.writetable(table, offset)
        if self.args.filename != '-': self.fp.close()


################################################################################
# shared functionality to create tabulation for bond or angle styles
class BondAngleTabulate(Tabulate):
    def __init__(self, style, efunc, ffunc=None):
        super(BondAngleTabulate, self).__init__(style, efunc, ffunc)
        self.parser.add_argument('--eshift', '-e', dest='eshift', default=False,
                                 action='store_true',
                                 help="Shift potential energy to be zero at minimum")
        idx = [a.dest for a in self.parser._actions].index('xmin')
        self.parser._actions[idx].required=False
        self.parser._actions[idx].default=0.0
        if style == 'angle':
            idx = [a.dest for a in self.parser._actions].index('xmax')
            self.parser._actions[idx].required=False
            self.parser._actions[idx].default=180.0
        try:
            self.args = self.parser.parse_args()
        except argparse.ArgumentError:
            sys.exit()

    def run(self, label):
        # sanity checks
        if self.args.num < 2:
            sys.exit('\nExpect 2 or more points in table for tabulation\n\n' + self.parser.format_help())
        if self.args.xmin < 0.0:
            sys.exit('\nInner cutoff must not be negative\n\n' + self.parser.format_help())
        if self.tstyle == 'angle' and self.args.xmax > 180.0:
            sys.exit('\nOuter cutoff must not be larger than 180.0 degrees\n\n' + self.parser.format_help())

        self.diff = self.args.diff
        if not self.forcefunc: self.diff = True

        (table, xzero) = mktable(self.tstyle, label, self.args.num, self.args.xmin, self.args.xmax,
                                 self.energyfunc, self.args.diff, self.forcefunc)
        print("# Minimum energy of tabulated potential is at %g" % xzero)
        offset = 0.0
        if self.args.eshift: offset=self.energyfunc(xzero)

        self.openfile(label)

        if self.forcefunc:
            diffmin = numdiff(self.args.xmin, self.forcefunc)
            diffmax = numdiff(self.args.xmax, self.forcefunc)
            self.fp.write("N %d FP %- 22.15g %- 22.15g EQ %g\n\n" %
                          (self.args.num, diffmin, diffmax, xzero))
        else:
            self.fp.write("N %d EQ %g\n\n" % (self.args.num, xzero))

        self.writetable(table, offset)
        if self.args.filename != '-': self.fp.close()

################################################################################
class BondTabulate(BondAngleTabulate):
    def __init__(self, efunc, ffunc=None):
        super(BondTabulate, self).__init__('bond', efunc, ffunc)

################################################################################
class AngleTabulate(BondAngleTabulate):
    def __init__(self, efunc, ffunc=None):
        super(AngleTabulate, self).__init__('angle', efunc, ffunc)

################################################################################
