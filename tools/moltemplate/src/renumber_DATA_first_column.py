#!/usr/bin/env python

"""
   Get rid of lines containing duplicate copies of the same atom in the "Atoms"
   section of a LAMMPS data file.  Duplicate lines which occur later are
   preserved and the earlier lines are erased.
   The file is read from sys.stdin.  This program does not parse the entire 
   data file.  The text from the "Atoms" section of the LAMMPS file must 
   be extracted in advance before it is sent to this program.)

"""

import sys
from operator import itemgetter

in_stream = sys.stdin
f = None
fname = None
if len(sys.argv) == 2:
    fname = sys.argv[1]
    f = open(fname, 'r')
    in_stream = f


lines = in_stream.readlines()
column1_iorig_columnsAfter1 = []

# Start at the end of the file and read backwards.
# If duplicate lines exist, eliminate the ones that occur earlier in the file.
i = 0
while i < len(lines):
    line_orig = lines[i]
    line = line_orig.rstrip('\n')
    if '#' in line_orig:
        ic = line.find('#')
        line = line_orig[:ic]

    tokens = line.strip().split()
    if len(tokens) > 0:
        if str.isdigit(tokens[0]):
            column1 = int(tokens[0])
        else:
            column1 = tokens[0]
        column1_iorig_columnsAfter1.append([column1, i, ' '.join(tokens[1:])])
    i += 1

# Sort the list of lines by the first number on each line

column1_iorig_columnsAfter1.sort(key=itemgetter(0))

# Change all of these numbers so that they are consecutive starting at 1
for i in range(0, len(column1_iorig_columnsAfter1)):
    column1_iorig_columnsAfter1[i][0] = i+1

# Sort the list of lines so they are back in the original order

column1_iorig_columnsAfter1.sort(key=itemgetter(1))

for i in range(0, len(column1_iorig_columnsAfter1)):
    column1       = column1_iorig_columnsAfter1[i][0]
    columnsAfter1 = column1_iorig_columnsAfter1[i][2]
    sys.stdout.write(str(column1) +' ' + columnsAfter1 + '\n')

if f != None:
    f.close()

