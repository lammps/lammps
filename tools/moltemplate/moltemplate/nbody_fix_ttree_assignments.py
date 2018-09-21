#!/usr/bin/env python

"""
nbody_fix_ttree_assignments.py

This is an ugly little script which was not intended to be run by end users.

Typical usage:

nbody_fix_ttree_assignments.py "angles" new_Angles.template \
  < ttree_assignments.txt > ttree_assigmnents_new.txt

What it does:

In this example, this program extracts the first column from
"new_Angles.template", and appends to it the first column from
the lines in ttree_assignments.txt containing "@angle:".
Then it adds a second column which is just a sequence of integers
counting upwards.
Finally it inserts this 2-column text into the appropriate place in a
ttree_assignments.txt file, replacing the original @angle variables
with the new ones (followed by the renumbered original).


   AWK/GREP equivalent
This program is roughly equivalent to the following lines of awk/grep:

awk 'BEGIN{i=-1} {if(substr($0,0,8)=="$/angle") i=NR; if (i==-1){print $0}}'\
    < ttree_assignments.txt > ttree_assignments_new.txt

awk '{print $1}' < new_Angles.template > Angles_column1.txt
grep '$/angle:' ttree_assignments.txt | awk '{print $1}' >> Angles_column1.txt
awk '{print $1 "  " NR}' <  Angles_column1.txt >> ttree_assignments_new.txt

awk 'BEGIN{found=0;passed=0} {if(substr($0,0,8)=="$/angle") found=1;
                              else {if (found) {passed=1}}
                              if (passed) print $0}' \
          < ttree_assignments.txt >> ttree_assignments_new.txt

I wrote this python script (instead of using awk) just to handle quoted stings
(and strings with other fancy characters and escape sequences).

"""

import sys

try:
    from .ttree_lex import SplitQuotedString, InputError
except (ImportError, SystemError, ValueError):
    # not installed as a package
    from ttree_lex import *

g_program_name = __file__.split('/')[-1]


def main():
    try:
        if (len(sys.argv) != 3):
            raise InputError('Error running  \"' + g_program_name + '\"\n'
                             '   Wrong number of arguments.\n'
                             '   (This is likely a programmer error.\n'
                             '    This script was not intended to be run by end users.)\n')

        cat_name = sys.argv[1]
        f = open(sys.argv[2])
        lines_generated = f.readlines()
        f.close()

        # Selections are simply lists of 2-tuples (pairs)
        #f = open('ttree_assignments.txt','r')
        #lines_bindings = f.readlines()
        # f.close()
        lines_bindings = sys.stdin.readlines()

        # Figure out which lines in the 'ttree_assignments.txt' file
        # contain the variables of the type you are looking for.
        # Make note of the relevant line numbers
        i_preexisting_begin = -1
        i_preexisting_end = -1
        in_section = False
        possible_cat_names = set(
            ['$' + cat_name, '$/' + cat_name, '${' + cat_name, '${/' + cat_name])

        preexisting_interaction_list = []
        for i in range(0, len(lines_bindings)):
            line = lines_bindings[i].strip()
            tokens = SplitQuotedString(line)  # strip comments, handle quotes
            if len(tokens) == 2:
                before_colon = tokens[0].split(':')[0]
                if before_colon in possible_cat_names:
                    if i_preexisting_begin == -1:
                        i_preexisting_begin = i
                        in_section = True
                else:
                    if in_section:
                        i_preexisting_end = i
                    in_section = False

        if i_preexisting_end == -1:
            i_preexisting_end = len(lines_bindings)

        if i_preexisting_begin == -1:
            for line in lines_bindings:
                sys.stdout.write(line)
        else:
            # write out all the lines in the original file up until the point where
            # the variables in the category we are looking for were encountered
            for i in range(0, i_preexisting_begin):
                sys.stdout.write(lines_bindings[i])

        sys.stderr.write('  (adding new lines)\n')

        # Now add some new lines (2-column format).
        # As with any ttree_assignment.txt file:
        #   The first column has our generated variable names
        #   The second column has the counter assigned to that variable
        new_counter = 1
        for line_orig in lines_generated:
            line = line_orig.strip()
            if len(line) > 0:
                tokens = SplitQuotedString(line)  # strip comments, handle quotes
                sys.stdout.write(tokens[0] + '  ' + str(new_counter) + '\n')
                new_counter += 1

        sys.stderr.write('  (adding pre-exisiting lines)\n')
        if i_preexisting_begin != -1:
            # Append the original pre-existing interactions of that type, but assign
            # them to higher numbers.  (Hopefully this helps to make sure that these
            # assignments will override any of the automatic/generated assignments.)
            # As with any ttree_assignment.txt file:
            #   The first column has our generated variable names
            #   The second column has the counter assigned to that variable

            # sys.stderr.write('  i_preexisting_begin='+
            #                 str(i_preexisting_begin)+
            #                 ' i_preexisting_end='+str(i_preexisting_end)+'\n')

            for i in range(i_preexisting_begin, i_preexisting_end):
                line = lines_bindings[i].strip()
                tokens = SplitQuotedString(line)  # strip comments, handle quotes
                if len(tokens) == 2:
                    sys.stdout.write(tokens[0] + '  ' + str(new_counter) + '\n')
                    new_counter += 1

            #sys.stderr.write('  (writing pre-exisiting lines)\n')

            # write out all the lines in the original file after this point.
            for i in range(i_preexisting_end, len(lines_bindings)):
                sys.stdout.write(lines_bindings[i])

        sys.exit(0)


    except (ValueError, InputError) as err:
        sys.stderr.write('\n' + str(err) + '\n')
        sys.exit(-1)

if __name__ == '__main__':
    main()
