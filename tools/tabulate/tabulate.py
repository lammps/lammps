#! /usr/bin/env python
"""
Python classes and support functions to generate tabulated potential files for LAMMPS.
"""

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

class Tabulate(object):

    def __init__(self, style, efunc, ffunc=None):
        self.fp = sys.stdout
        self.tstyle = style
        self.energyfunc = efunc
        self.forcefunc = ffunc
        self.eshift = 0.0

        self.parser = argparse.ArgumentParser(description='Tool to generate tabulated '
                                              + self.tstyle + ' potential files for LAMMPS')
        self.parser.add_argument('--num-points', '-n', dest='num', default=1000,
                                 help="Number of tabulated points")
        self.parser.add_argument('--filename', '-f', dest='filename', default='-',
                                 help="Name of output file")
        self.parser.add_argument('--diff-num', '-d', dest='diff',default=False,
                                 action='store_true',
                                 help="Differentiate energy function numerically")
        self.parser.add_argument('--inner', '-i', dest='xmin', required=True,
                                 help="Inner cutoff of table")
        self.parser.add_argument('--outer', '-o', dest='xmax', required=True,
                                 help="Outer cutoff of table")

    def openfile(self, label):
        if self.args.filename != '-':
            self.fp = open(self.args.filename, 'a')

        self.fp.write('\n' + label + '\n')

    def writetable(self, table, offset):
        for i,r,e,f in table:
            self.fp.write("%8d %- 22.15g %- 22.15g %- 22.15g\n" % (i, r, e - offset, f))

################################################################################

class PairTabulate(Tabulate):
    def __init__(self, efunc, ffunc=None):
        super(PairTabulate, self).__init__('pair', efunc, ffunc)
        self.parser.add_argument('--eshift', '-e', dest='eshift', default=False, action='store_true',
                                  help="Shift potential energy to be zero at outer cutoff")
        try:
            self.args = self.parser.parse_args(sys.argv[1:])
        except:
            sys.exit()

    def run(self, label):
        num = int(self.args.num)
        xmin = float(self.args.xmin)
        xmax = float(self.args.xmax)
        self.diff = self.args.diff
        if not self.forcefunc: self.diff = True
        offset = 0.0
        if self.args.eshift: offset=self.energyfunc(xmax)

        (table, dummy) = mktable(self.tstyle, label, num, xmin, xmax,
                                 self.energyfunc, self.args.diff, self.forcefunc)

        # open table file
        self.openfile(label)

        # write pair style specific header
        if self.forcefunc:
            diffmin = numdiff(xmin, self.forcefunc)
            diffmax = numdiff(xmax, self.forcefunc)
            self.fp.write("N %d R %g %g FPRIME %- 22.15g %- 22.15g\n\n"
                          % (num, xmin, xmax, diffmin, diffmax))
        else:
            self.fp.write("N %d R %g %g\n\n" % (num, xmin, xmax))

        self.writetable(table, offset)
        if self.args.filename != '-': self.fp.close()
        

################################################################################

class BondTabulate(Tabulate):
    def __init__(self, efunc, ffunc=None):
        super(BondTabulate, self).__init__('bond', efunc, ffunc)
        self.parser.add_argument('--eshift', '-e', dest='eshift', default=False, action='store_true',
                                  help="Shift potential energy to be zero at minimum")
        try:
            self.args = self.parser.parse_args(sys.argv[1:])
        except:
            sys.exit()

    def run(self, label):
        num = int(self.args.num)
        xmin = float(self.args.xmin)
        xmax = float(self.args.xmax)
        self.diff = self.args.diff
        if not self.forcefunc: self.diff = True

        (table, xzero) = mktable(self.tstyle, label, num, xmin, xmax,
                                 self.energyfunc, self.args.diff, self.forcefunc)
        print("# Minimum energy of tabulated potential is at %g" % xzero)
        offset = 0.0
        if self.args.eshift: offset=self.energyfunc(xzero)

        self.openfile(label)

        if self.forcefunc:
            diffmin = numdiff(xmin, self.forcefunc)
            diffmax = numdiff(xmax, self.forcefunc)
            self.fp.write("N %d FP %- 22.15g %- 22.15g EQ %g\n\n" % (num, diffmin, diffmax, xzero))
        else:
            self.fp.write("N %d EQ %g\n\n" % (num, xzero))

        self.writetable(table, offset)
        if self.args.filename != '-': self.fp.close()

################################################################################
