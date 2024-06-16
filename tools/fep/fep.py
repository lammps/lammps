#!/usr/bin/env python
# fep.py - calculate free energy from compute fep results

import sys
from argparse import ArgumentParser
import math

def compute_fep():

    parser = ArgumentParser(description='A python script to calculate free energy from compute fep results')

    parser.add_argument("units", help="unit system can be lj, real or si")
    parser.add_argument("Temperature", type=float, help="The temperature of the system")
    parser.add_argument("InputFile", help="The name of the file with the data from compute fep.")

    args = parser.parse_args()

    r_value = {'lj': 1.0, 'real': 0.0019872036, 'si': 8.31446}

    if args.units in r_value:
        rt = r_value[args.units] * args.Temperature
    else:
        sys.exit("The provided units keyword is not valid")

    v = 1.0
    mysum = 0.0

    with open(args.InputFile, "r") as input_file:
        for line in input_file:
            if line.startswith("#"):
                continue
            tok = line.split()
            if len(tok) == 4:
                v = float(tok[3])
            mysum += math.log(float(tok[2]) / v)

    print(-rt * mysum)

if __name__ == "__main__":
    compute_fep()

