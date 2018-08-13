#!/usr/bin/env python
# Author: Andrew Jewett (jewett.aij at g mail)
# License: 3-clause BSD License  (See LICENSE.TXT)
# Copyright (c) 2017, California Institute of Technology
# All rights reserved.

"""
   Reorder the integer arguments to the commands in a LAMMPS input
   file if these arguments violate LAMMPS order requirements.
   We have to do this because the moltemplate.sh script will automatically
   assign these integers in a way which may violate these restrictions
   and the user has little control over this.

   This script:
   swaps the I and J integers in    "pair_coeff I J ..." commands when I > J

   Other features may be added later

"""


import sys

def main():
    lines_orig = []
    f = None
    fname = None
    num_lines_ignore = 0


    # Lines from files passed as arguments are read and processed silently.
    # (Why? Sometimes it's necessary to read the contents of previous input scripts
    #  in order to be able to understand a script command which appears later.
    #  I'm assuming these files will be processed by lammps in the same order. So I
    #  must insure that moltemplate.sh passes them to this program in that order.
    #  I'm too lazy to read the "include" commands in input scripts correctly.)
    if len(sys.argv) > 1:
        for fname in sys.argv[1:]:
            f = open(fname, 'r')
            in_stream = f
            lines_orig += in_stream.readlines()
            num_lines_ignore += len(lines_orig)
            f.close()


    # Lines read from the standard input are read, processed, and printed to stdout
    in_stream = sys.stdin
    lines_orig += in_stream.readlines()

    pair_style_list = []
    swap_occured = False
    warn_wildcard = False

    i = 0
    while i < len(lines_orig):
        # Read the next logical line
        # Any lines ending in '&' should be merged with the next line before
        # breaking
        line_orig = ''
        while i < len(lines_orig):
            line_counter = 1 + i - num_lines_ignore
            line_orig += lines_orig[i]
            if ((len(line_orig) < 2) or (line_orig[-2:] != '&\n')):
                break
            i += 1
        line = line_orig.replace('&\n', '\n').rstrip('\n')

        comment = ''
        if '#' in line_orig:
            ic = line.find('#')
            line = line_orig[:ic]
            # keep track of comments (put them back later)
            comment = line_orig[ic:].rstrip()

        tokens = line.strip().split()
        if ((len(tokens) >= 2) and (tokens[0] == 'pair_style')):
            pair_style_list = tokens[1:]

        if ((len(tokens) >= 3) and (tokens[0] == 'pair_coeff')):

            if ((tokens[1].isdigit() and (tokens[2].isdigit())) and
                (int(tokens[1]) > int(tokens[2]))):

                swap_occured = True
                tmp = tokens[2]
                tokens[2] = tokens[1]
                tokens[1] = tmp

                if i >= num_lines_ignore:

                    # polite warning:
                    sys.stderr.write(
                        'swapped pair_coeff order on line ' + str(line_counter))
                    # if (fname != None):
                    #    sys.stderr.write(' of file \"'+fname+'\"')
                    sys.stderr.write('\n')

                    # Deal with the "hbond/" pair coeffs.
                    #
                    # The hbond/dreiding pair style designates one of the two atom types
                    # as a donor, and the other as an acceptor (using the 'i','j' flags)
                    # If swapped atom types eariler, we also need to swap 'i' with 'j'.
                    #
                    # If "hbond/dreiding.." pair style is used with "hybrid" or
                    # "hybrid/overlay" then tokens[3] is the name of the pair style
                    # and tokens[5] is either 'i' or 'j'.
                    if len(pair_style_list) > 0:
                        if ((pair_style_list[0] == 'hybrid') or
                                (pair_style_list[0] == 'hybrid/overlay')):
                            if ((len(tokens) > 5) and (tokens[5] == 'i') and (tokens[3][0:6] == 'hbond/')):
                                tokens[5] = 'j'
                                sys.stderr.write(
                                    '  (and replaced \"i\" with \"j\")\n')
                            elif ((len(tokens) > 5) and (tokens[5] == 'j') and (tokens[3][0:6] == 'hbond/')):
                                tokens[5] = 'i'
                                sys.stderr.write(
                                    '  (and replaced \"j\" with \"i\")\n')
                        elif (pair_style_list[0][0:6] == 'hbond/'):
                            if ((len(tokens) > 4) and (tokens[4] == 'i')):
                                tokens[4] = 'j'
                                sys.stderr.write(
                                    '  (and replaced \"i\" with \"j\")\n')
                            elif ((len(tokens) > 4) and (tokens[4] == 'j')):
                                tokens[4] = 'i'
                                sys.stderr.write(
                                    '  (and replaced \"j\" with \"i\")\n')

                    sys.stdout.write(
                        (' '.join(tokens) + comment).replace('\n', '&\n') + '\n')

            else:
                if ((('*' in tokens[1]) or ('*' in tokens[2]))
                    and
                    (not (('*' == tokens[1]) and ('*' == tokens[2])))):
                    warn_wildcard = True
                if i >= num_lines_ignore:
                    sys.stdout.write(line_orig)
        else:
            if i >= num_lines_ignore:
                sys.stdout.write(line_orig)

        i += 1


    if swap_occured:
        sys.stderr.write('\n'
                         '  WARNING: Atom order in some pair_coeff commands was swapped to pacify LAMMPS.\n'
                         '  For some exotic pair_styles such as hbond/dreiding, this is not enough. If you\n'
                         '  use exotic pair_styles, please verify the \"pair_coeff\" commands are correct.\n')

    if warn_wildcard:
        sys.stderr.write('\n'
                         '  WARNING: The use of wildcard characters (\"*\") in your \"pair_coeff\"\n'
                         '           commands is not recommended.\n'
                         '           (It is safer to specify each interaction pair manually.\n'
                         '            Check every pair_coeff command.  Make sure that every atom type in\n'
                         '            the first group is <= atom types in the second group.\n'
                         '            Moltemplate does NOT do this when wildcards are used.)\n'
                         '        If you are using a many-body pair style then ignore this warning.\n')

    return

if __name__ == '__main__':
    main()
