#!/usr/bin/env python3

from tabulate import PairTabulate
import sys
import argparse
import numpy as np
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

"""
   This script gives an example on how to make tabulated forces from radial
   distribution function using tabulate.py.
   Required: python3, numpy, scipy.
   BI stands for Boltzmann Inversion.
   WARNING: Using BI does not garranty an out of the box working potential for
            your simulation. Check the relevamt literature.
"""
###############################################################################


class BI(PairTabulate):
    def __init__(self, units=None, comment=None):
        super(PairTabulate, self).__init__("pair", self.energy, units, comment)
        self.parser.add_argument(
            "--eshift",
            "-e",
            dest="eshift",
            default=False,
            action="store_true",
            help="Shift potential energy to be zero at outer cutoff",
        )
        self.parser.add_argument("--units", required=True, help="LAMMPS units to use")
        self.parser.add_argument("--rdffile", required=True, help="Rdf file to be read")
        self.parser.add_argument("--temperature", required=True, type=float,
                                 help="Temperature for Boltzman Inversion.")
        try:
            self.args = self.parser.parse_args()
        except argparse.ArgumentError:
            sys.exit()

        if self.args.temperature > 0:
            T = self.args.temperature
        else:
            sys.exit("Invalid temperature provided.")

        kb = 1.0
        # Add more kb units if you need
        self.units = self.args.units
        if self.units == "lj":
            kb = 1.0           # reduced
        elif self.units == "si":
            kb = 1.380649e-23  # J/K
        elif self.units == "metal":
            kb = 8.617333e-5   # eV/K
        elif self.units == "real":
            kb = 1.987204e-3   # kcal/mol/K
        else:
            sys.exit("Unsupported units setting " + self.units)
        self.kbT = kb * T
        self.r, self.e = self.read_rdf(self.args.rdffile)

    # This function assumes LAMMPS format for rdf with a single entry
    def read_rdf(self, rdffile):
        data = np.loadtxt(rdffile, skiprows=4)
        r = data[:, 1]
        g = data[:, 2]

        # savgol_filter is an example of smoothing.
        # Other filters/functions can be used.
        g = savgol_filter(g, 9, 5)
        return self.inversion(r, g)

    def inversion(self, r, g):
        e = -self.kbT * np.log(np.clip(g,1.0e-100,1.0e100))
        e = self.complete_exponential(r, e)
        return r, e,

    def complete_exponential(self, r, e):
        r_temp = r[e != np.inf]
        e_temp = e[e != np.inf]

        # Optimising the parameter for a function for derivation
        # to be continuous.
        # Here a gaussian function, can be anything relevant defined in func.
        popt, pcov = curve_fit(self.func, r_temp[:2], e_temp[:2])
        for i, _ in enumerate(e):
            if e[i] == np.inf:
                e[i] = self.func(r[i], *popt)
        return e

    def func(self, x, K, s):
        return K * np.exp(-0.5 * (x / s) ** 2) / (s * np.sqrt(2 * np.pi))

    def energy(self, x):
        e = self.e
        r = self.r
        # Force estimation at minimum distance.
        # Should not be that useful
        f0 = (e[1] - e[0]) / (r[1] - r[0])

        minr = min(r)
        maxr = max(r)
        # Note that you might want OOB to return an error.
        if x >= maxr:
            return 0
        if x < minr:
            dx = minr - x
            return -f0 * dx
        else:
            # Linear interpolation between points.
            for i, ri in enumerate(r):
                if r[i] < x:
                    r1, e1 = r[i], e[i]
                    r2, e2 = r[i + 1], e[i + 1]
            dr12 = r2 - r1
            dr = x - r1
            de = (e2 - e1) / (r2 - r1)
            return e1 + (de * dr / dr12)


###############################################################################


if __name__ == "__main__":
    ptable = BI()
    ptable.run("BI")
