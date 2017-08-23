#!/usr/bin/env python

import sys

def main():
    for line_orig in sys.stdin:
        line = line_orig.rstrip('\n')
        comment = ''
        if '#' in line_orig:
            ic = line.find('#')
            line = line_orig[:ic]
            comment = ' '+line_orig[ic:].rstrip('\n')

        tokens = line.strip().split()
        if len(tokens) > 2:
            atomid = -1
            atomtype = -1
            pos_found = False
            for i in range(0,len(tokens)):
                if (tokens[i] == 'part') and (i+1 < len(tokens)):
                    atomid = tokens[i+1]
                elif (tokens[i] == 'type') and (i+1 < len(tokens)):
                    atomtype = tokens[i+1]
                elif (tokens[i] == 'pos') and (i+2 < len(tokens)):
                    pos_found = True
            if (atomid != -1) and (atomtype != -1) and pos_found:
                sys.stdout.write(atomid+' '+atomtype+'\n')

if __name__ == "__main__":
    main()
