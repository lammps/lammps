#!/usr/bin/env python

man_page_text = """
Usage (example):

ttree_render.py ttree_assignments.txt < file.template > file.rendered

The argument (ttree_assignments.txt) should be a 2-column file containing
ttree-style variables (1st column), and their values (bindings, 2nd column).

This program reads a text file containing ttree-style variables,
substitutes the corresponding values stored in ttree_assignments.txt,
and prints out the new (rendered) text to the standard-out.

"""


import sys
from ttree_lex import *
from ttree import *
import gc


g_filename    = __file__.split('/')[-1]
g_module_name  = g_filename
if g_filename.rfind('.py') != -1:
    g_module_name = g_filename[:g_filename.rfind('.py')]
g_date_str     = '2012-9-06'
g_version_str  = '0.1'
g_program_name = g_filename
#sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+' ')



try:

    if (len(sys.argv) != 2):
        raise InputError('Error running  \"'+g_program_name+'\"\n'
                         ' Typical usage:\n'
                         ' ttree_render.py ttree_assignments.txt < file.template > file.rendered\n'
                         '\n'
                         '   Missing argument.\n'
                         '   Expected the name of a 2-column file containing\n'
                         '   variable names and their bindings (values).\n'
                         '   (This is likely a programmer error.\n'
                         '    This script was not intended to be run by end users.)\n')

    bindings_filename = sys.argv[1]
    f = open(bindings_filename)
    assignments = {}

    #BasicUIReadBindingsStream(assignments, f, bindings_filename)

    # The line above is robust but it uses far too much memory.
    # This for loop below works for most cases.
    for line in f:
        #tokens = lines.strip().split()
        tokens = SplitQuotedString(line.strip()) # like split but handles quotes
        if len(tokens) < 2:
            continue
        assignments[tokens[0]] = tokens[1]

    f.close()
    gc.collect()

    lex = TemplateLexer(sys.stdin, '__standard_input_for_ttree_render__')
    lex.var_delim = '$@'

    text_block_list = lex.ReadTemplate(simplify_output=True)

    output = []

    for entry in text_block_list:
        assert(isinstance(entry, str))

        if ((len(entry) > 1) and (entry[0] in lex.var_delim)):
            if '.' in entry:
                ic = entry.find('.')
                var_name = entry[:ic]
                var_suffix = entry[ic:]
            else:
                var_name = entry
                var_suffix = ''

            var_name = entry
            if var_name not in assignments:
                raise(InputError('Error('+g_program_name+')'
                                 #' at '+ErrorLeader(var_ref.src_loc.infile,
                                 #                  var_ref.src_loc.lineno)+
                                 ' unknown variable:\n'
                                 '         \"'+var_name+'\"\n'))
            else:
                var_value = assignments[var_name]

            format_fname, args = ExtractFormattingCommands(var_suffix)
            if format_fname == 'ljust':
                if len(args) == 1:
                    var_value = var_value.ljust(int(args[0]))
                else:
                    var_value = var_value.ljust(int(args[0]), args[1])
            elif format_fname == 'rjust':
                if len(args) == 1:
                    var_value = var_value.rjust(int(args[0]))
                else:
                    var_value = var_value.rjust(int(args[0]), args[1])
            output.append(var_value)
        else:
            output += entry

    sys.stdout.write(''.join(output))


except (ValueError, InputError) as err:
    sys.stderr.write('\n'+str(err)+'\n')
    sys.exit(-1)
