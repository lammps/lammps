#!/usr/bin/env python

# Author: Andrew Jewett (jewett.aij at g mail)
#         http://www.chem.ucsb.edu/~sheagroup
# License: 3-clause BSD License  (See LICENSE.TXT)
# Copyright (c) 2011, Regents of the University of California
# All rights reserved.

"""
ettree.py

ettree.py is an extension of the generic ttree.py program.
This version can understand and manipulate ttree-style templates which 
are specialized for storing molecule-specific data for use in ESPresSo/TCL.

The main difference between ettree.py and ttree.py is:
Unlike ttree.py, ettree.py understands rigid-body movement commands like 
"rot()" and "move()" which allows it to reorient and move each copy
of a molecule to a new location.  (ttree.py just ignores these commands.
Consequently ESPresSo/TCL input file (fragments) created with ttree.py have
invalid (overlapping) atomic coordinates and must be modified or aguemted 
later (by loading atomic coordinates from a PDB file or an XYZ file).
ettree.py understands and can manipulate atomic coordinates.

Additional ESPresSo/TCL-specific features may be added in the future.

"""

import sys

# Problem:
# One of the python files I need is in a different git repository
# which is linked to the parent directory using "git subtree".
# The result of this is that he python code I need to access is in
# a directory which is outside the current one (in "../moltemplate/src")
# For now, I'm willing to resort to using a hack to import this file.

import os, inspect
# use this if you want to include modules from a subfolder
# http://stackoverflow.com/questions/279237/import-a-module-from-a-relative-path
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"..","moltemplate","src")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)


try:
    from .ttree import BasicUISettings, BasicUIParseArgs, EraseTemplateFiles, \
        StackableCommand, PopCommand, PopRightCommand, PopLeftCommand, \
        PushCommand, PushLeftCommand, PushRightCommand, ScopeCommand, \
        WriteVarBindingsFile, StaticObj, InstanceObj, \
        BasicUI, ScopeBegin, ScopeEnd, WriteFileCommand, Render
    from .ttree_lex import InputError, TextBlock, DeleteLinesWithBadVars, \
        TemplateLexer
    from .ettree_styles import espt_delim_atom_fields, \
        LinesWSlashes, SplitMultiDelims, SplitAtomLine, \
        iEsptAtomCoords, iEsptAtomVects, iEsptAtomType, iEsptAtomID
    from .ttree_matrix_stack import AffineTransform, MultiAffineStack, \
        LinTransform
except (ImportError, SystemError, ValueError):
    # not installed as a package
    from ttree import *
    from ttree_lex import *
    from ettree_styles import *
    from ttree_matrix_stack import *

try:
    unicode
except NameError:
    # Python 3
    basestring = unicode = str


data_atoms = 'Data Atoms' # <-- The name of the file/section storing Atom data.



class EttreeSettings(BasicUISettings):
    """ Currently EttreeSettings is identical to BasicUISettings.
    Later on, if I find I need to add custom settings specific to ESPresSoTCL,
    I will add them here.
    (See "class LttreeSettings" in "lttree.py" for comparison.)

    """

    def __init__(self,
                 user_bindings_x=None,
                 user_bindings=None,
                 order_method='by_command'):

        BasicUISettings.__init__(self, 
                                 user_bindings_x, 
                                 user_bindings, 
                                 order_method)



def EttreeParseArgs(argv, settings):
    """
    This function currently does nothing except invoke BasicUIParseArgs()
    (and throw an error message if the user forgot to specify an file name).
      Later on, if I think of some command-line arguments specific 
      to ESPresSo(tcl), then I will deal with them here.
    (See the "LttreeParseArgs()" function in "lttree.py" for comparison.)

    """

    BasicUIParseArgs(argv, settings)

    # Loop over the remaining arguments not processed yet.
    # These arguments are specific to the ettree.py program
    # and are not understood by ttree.py:
    i = 1
    while i < len(argv):
        #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')
        if ((argv[i][0] == '-') and (__name__ == "__main__")):
            #elif (__name__ == "__main__"):
            raise InputError('Error('+__file__+'):\n'
                             'Unrecogized command line argument \"'+argv[i]+'\"\n')
        else:
            i += 1


    if __name__ == "__main__":

        # Instantiate the lexer we will be using.
        #  (The lexer's __init__() function requires an openned file.
        #   Assuming __name__ == "__main__", then the name of that file should
        #   be the last remaining (unprocessed) argument in the argument list.
        #   Otherwise, then name of that file will be determined later by the
        #   python script which imports this module, so we let them handle it.)

        if len(argv) == 1:
            raise InputError('Error: This program requires at least one argument\n'
                             '       the name of a file containing ttree template commands\n')
        elif len(argv) == 2:
            try:
                # Parse text from the file named argv[1]
                settings.lex.infile = argv[1]  
                settings.lex.instream = open(argv[1], 'r')
            except IOError: 
                sys.stderr.write('Error: unable to open file\n'
                                 '       \"'+argv[1]+'\"\n'
                                 '       for reading.\n')
                sys.exit(1)
            del(argv[1:2])

        else:
            # if there are more than 2 remaining arguments,
            problem_args = ['\"'+arg+'\"' for arg in argv[1:]]
            raise InputError('Syntax Error('+__file__+'):\n\n'
                             '       Problem with argument list.\n'
                             '       The remaining arguments are:\n\n'
                             '         '+(' '.join(problem_args))+'\n\n'
                             '       (The actual problem may be earlier in the argument list.\n'
                             '       If these arguments are source files, then keep in mind\n'
                             '       that this program can not parse multiple source files.)\n'
                             '       Check the syntax of the entire argument list.\n')





def TransformAtomText(text, matrix):
    """ Apply transformations to the coordinates and other vector degrees 
    of freedom stored in the atom declaration section.
    This is the \"text\" argument.
    The \"matrix\" stores the aggregate sum of combined transformations
    to be applied.

    """

    #sys.stderr.write('matrix_stack.M = \n'+ MatToStr(matrix) + '\n')

    # lines = text.split('\n') <-- this does not work because a backslash at
    #                              the end of a line can merge multiple lines

    lines = [line for line in LinesWSlashes(text)]  # <-- handles backslashes
    for i in range(0, len(lines)):
        line = lines[i]
        tokens = SplitAtomLine(line)
        if len(tokens) > 0:
            x0 = [0.0, 0.0, 0.0]
            x  = [0.0, 0.0, 0.0]
            for icrd in iEsptAtomCoords(tokens):
                coords_str = tokens[icrd].split()
                for d in range(0,3):
                    x0[d] = float(coords_str[d])
                AffineTransform(x, matrix, x0)  # x = matrix * x0 + b
                for d in range(0,3):            # ("b" is part of "matrix")
                    coords_str[d] = str(x[d])
                tokens[icrd] = ' '.join(coords_str)

            for ivect in iEsptAtomVects(tokens):
                coords_str = tokens[ivect].split()
                for d in range(0,3):
                    x0[d] = float(coords_str[d])
                LinearTransform(x, matrix, x0)  # x = matrix * x0
                for d in range(0,3):
                    coords_str[d] = str(x[d])
                tokens[ivect] = ' '.join(coords_str)
        line = ' '.join(tokens)
        lines[i] = line
    return '\n'.join(lines)+'\n'


# NOT IMPLEMENTED YET:
def CalcCM(text_Atoms, 
           text_Masses=None, 
           settings=None):
    # FILL IN THE CONTENTS OF THIS FUNCTION LATER
    xcm = [0.0, 0.0, 0.0]
    return xcm



def _ExecCommands(command_list,
                  index,
                  global_files_content,
                  settings,
                  matrix_stack,
                  current_scope_id=None,
                  substitute_vars=True):
    """ 
    _ExecCommands():
    The argument "commands" is a nested list of lists of 
    "Command" data structures (defined in ttree.py).

    Carry out the write() and write_once() commands (which 
    write out the contents of the templates contain inside them).
    Instead of writing the files, save their contents in a string.

    The argument "global_files_content" should be of type defaultdict(list)
    It is an associative array whose key is a string (a filename)
    and whose value is a lists of strings (of rendered templates).

    """
    files_content = defaultdict(list)
    postprocessing_commands = []

    while index < len(command_list):
        command = command_list[index]
        index += 1

        # For debugging only
        # sys.stderr.write(str(command)+'\n')

        if isinstance(command, PopCommand):
            assert(current_scope_id != None)
            if command.context_node == None:
                command.context_node = current_scope_id
            if isinstance(command, PopRightCommand):
                matrix_stack.PopRight(which_stack = command.context_node)
            elif isinstance(command, PopLeftCommand):
                matrix_stack.PopLeft(which_stack = command.context_node)
            else:
                assert(False)

        elif isinstance(command, PushCommand):
            assert(current_scope_id != None)
            if command.context_node == None:
                command.context_node = current_scope_id
            # Some commands are post-processing commands, and must be
            # carried out AFTER all the text has been rendered.  For example
            # the "movecm(0,0,0)" waits until all of the coordinates have
            # been rendered, calculates the center-of-mass, and then applies
            # a translation moving the center of mass to the origin (0,0,0).
            # We need to figure out which of these commands need to be 
            # postponed, and which commands can be carried out now.
            # ("now"=pushing transformation matrices onto the matrix stack).
            # UNFORTUNATELY POSTPONING SOME COMMANDS MAKES THE CODE UGLY
            transform_list = command.contents.split('.')
            transform_blocks = []
            i_post_process = -1
            #  Example:  Suppose:
            #command.contents = '.rot(30,0,0,1).movecm(0,0,0).rot(45,1,0,0).scalecm(2.0).move(-2,1,0)'
            #  then
            #transform_list = ['rot(30,0,0,1)', 'movecm(0,0,0)', 'rot(45,1,0,0)', 'scalecm(2.0)', 'move(-2,1,0)']
            # Note: the first command 'rot(30,0,0,1)' is carried out now. 
            # The remaining commands are carried out during post-processing,
            # (when processing the "ScopeEnd" command.
            # 
            # We break up the commands into "blocks" separated by center-
            # of-mass transformations ('movecm', 'rotcm', or 'scalecm')
            # 
            # transform_blocks = ['.rot(30,0,0,1)',
            #                     '.movecm(0,0,0).rot(45,1,0,0)',
            #                     '.scalecm(2.0).move(-2,1,0)']

            i = 0
            while i < len(transform_list):
                transform_block = ''
                while i < len(transform_list):
                    transform = transform_list[i]
                    i += 1
                    if transform != '':
                        transform_block += '.' + transform
                        transform = transform.split('(')[0]
                        if ((transform == 'movecm') or
                            (transform == 'rotcm') or
                            (transform == 'scalecm')):

                            #break

                            raise InputError("Error: center-of-mass transformations are not yet implemented\n"
                                             "       Avoid using \""+transform+"()\" transformations.\n")

                transform_blocks.append(transform_block)

            if len(postprocessing_commands) == 0:
                # The first block (before movecm, rotcm, or scalecm)
                # can be executed now by modifying the matrix stack.
                if isinstance(command, PushRightCommand):
                    matrix_stack.PushCommandsRight(transform_blocks[0].strip('.'),
                                                   command.srcloc,
                                                   which_stack=command.context_node)
                elif isinstance(command, PushLeftCommand):
                    matrix_stack.PushCommandsLeft(transform_blocks[0].strip('.'),
                                                  command.srcloc,
                                                  which_stack=command.context_node)
                # Everything else must be saved for later.
                postprocessing_blocks = transform_blocks[1:]
            else:
                # If we already encountered a "movecm" "rotcm" or "scalecm"
                # then all of the command blocks must be handled during
                # postprocessing.
                postprocessing_blocks = transform_blocks

            for transform_block in postprocessing_blocks:
                assert(isinstance(block, basestring))
                if isinstance(command, PushRightCommand):
                    postprocessing_commands.append(PushRightCommand(transform_block,
                                                                    command.srcloc,
                                                                    command.context_node))
                elif isinstance(command, PushLeftCommand):
                    postprocessing_commands.append(PushLeftCommand(transform_block,
                                                                   command.srcloc,
                                                                   command.context_node))


        elif isinstance(command, WriteFileCommand):

            # --- Throw away lines containin references to deleted variables:---

            # First: To edit the content of a template, 
            #        you need to make a deep local copy of it
            tmpl_list = []
            for entry in command.tmpl_list:
                if isinstance(entry, TextBlock):
                    tmpl_list.append(TextBlock(entry.text, 
                                               entry.srcloc)) #, entry.srcloc_end))
                else:
                    tmpl_list.append(entry)


            # --- Now throw away lines with deleted variables ---

            DeleteLinesWithBadVars(tmpl_list)

            # --- Now render the text ---
            text = Render(tmpl_list, 
                          substitute_vars)

            # ---- Coordinates of the atoms, must be rotated 
            # and translated after rendering.
            # In addition, other vectors (dipoles, ellipsoid orientations)
            # must be processed.
            # This requires us to re-parse the contents of this text
            # (after it has been rendered), and apply these transformations
            # before passing them on to the caller.
            if command.filename == data_atoms:
                text = TransformAtomText(text, matrix_stack.M)

            files_content[command.filename].append(text)


        elif isinstance(command, ScopeBegin):

            if isinstance(command.node, InstanceObj):
                if ((command.node.children != None) and 
                    (len(command.node.children) > 0)):
                    matrix_stack.PushStack(command.node)

            # "command_list" is a long list of commands.
            # ScopeBegin and ScopeEnd are (usually) used to demarcate/enclose
            # the commands which are issued for a single class or 
            # class instance.  _ExecCommands() carries out the commands for 
            # a single class/instance.  If we reach a ScopeBegin(), 
            # then recursively process the commands belonging to the child.
            index = _ExecCommands(command_list,
                                  index,
                                  files_content, 
                                  settings,
                                  matrix_stack, 
                                  command.node,
                                  substitute_vars)

        elif isinstance(command, ScopeEnd):
            if 'Data Atoms' in files_content:
                for ppcommand in postprocessing_commands:
                    if 'Data Masses' in files_content:
                        pass
                        #xcm = CalcCM(files_content['Data Atoms'],
                        #             files_content['Data Masses'],
                        #             settings)
                    else:
                        pass
                        #xcm = CalcCM(files_content['Data Atoms'])

                    if isinstance(ppcommand, PushRightCommand):
                        matrix_stack.PushCommandsRight(ppcommand.contents,
                                                       ppcommand.srcloc,
                                                       xcm,
                                                       which_stack=command.context_node)
                    elif isinstance(ppcommand, PushLeftCommand):
                        matrix_stack.PushCommandsLeft(ppcommand.contents,
                                                      ppcommand.srcloc,
                                                      xcm,
                                                      which_stack=command.context_node)
                    files_content['Data Atoms'] = \
                        TransformAtomText(Files_content['Data Atoms'],
                                          matrix_stack.M)

                for ppcommand in postprocessing_commands:
                    matrix_stack.Pop(which_stack = command.context_node)
                    #(same as PopRight())

            if isinstance(command.node, InstanceObj):
                if ((command.node.children != None) and 
                    (len(command.node.children) > 0)):
                    matrix_stack.PopStack()

            # "ScopeEnd"  means we're done with this class/instance.
            break

        else:
            assert(False)
            # no other command types allowed at this point


    # After processing the commands in this list,
    # merge the templates with the callers template list
    for file_name, tmpl_list in files_content.items():
        global_files_content[file_name] += \
            files_content[file_name]
        
    return index



def ExecCommands(commands, 
                 files_content, 
                 settings, 
                 substitute_vars=True):

    matrix_stack = MultiAffineStack()

    index = _ExecCommands(commands,
                          0,
                          files_content, 
                          settings,
                          matrix_stack, 
                          None,
                          substitute_vars)
    assert(index == len(commands))




def WriteFiles(files_content, suffix='', write_to_stdout=True):
    for file_name, str_list in files_content.items():
        if file_name != None:
            out_file = None
            if file_name == '':
                if write_to_stdout:
                    out_file = sys.stdout
            else:
                out_file = open(file_name+suffix, 'a')
            if out_file != None:
                out_file.write(''.join(str_list))
                if file_name != '':
                    out_file.close()



def main():
    """
    This is is a "main module" wrapper for invoking ettree.py
    as a stand alone program.  This program:

    1)reads a ttree file, 
    2)constructs a tree of class definitions (g_objectdefs)
    3)constructs a tree of instantiated class objects (g_objects),
    4)automatically assigns values to the variables,
    5)and carries out the "write" commands to write the templates a file(s).

    """
    g_program_name = 'ettree.py'
    g_date_str     = '2018-6-26'
    g_version_str  = '0.37.0'

    SimpleCounter.default_n0 = 0   # counters in Espresso begin at 0, not 1

    #######  Main Code Below: #######
    sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+' ')
    sys.stderr.write('\n(python version '+str(sys.version)+')\n')
    if sys.version < '2.6':
        raise InputError('Error: Alas, you must upgrade to a newever version of python.')

    try:

        #settings = BasicUISettings()
        #BasicUIParseArgs(sys.argv, settings)
        settings = EttreeSettings()
        EttreeParseArgs(sys.argv, settings)

        # Data structures to store the class definitionss and instances
        g_objectdefs = StaticObj('', None) # The root of the static tree
                                           # has name '' (equivalent to '/')
        g_objects    = InstanceObj('', None) # The root of the instance tree
                                             # has name '' (equivalent to '/')

        # A list of commands to carry out
        g_static_commands = []
        g_instance_commands = []


        BasicUI(settings,
                g_objectdefs, 
                g_objects,
                g_static_commands,
                g_instance_commands)

        # Now, carry out the commands
        # This involves rendering the templates and post-processing them.

        sys.stderr.write(' done\nbuilding templates...')

        files_content = defaultdict(list)

        ExecCommands(g_static_commands, 
                     files_content, 
                     settings, 
                     False)
        ExecCommands(g_instance_commands, 
                     files_content,
                     settings, 
                     False)

        # Erase the files that will be written to:
        sys.stderr.write(' done\nwriting templates...')
        EraseTemplateFiles(g_static_commands)
        EraseTemplateFiles(g_instance_commands)

        # Write the files as templates 
        # (with the original variable names present)
        WriteFiles(files_content, suffix=".template", write_to_stdout=False)

        # Write the files with the variables substituted by values
        sys.stderr.write(' done\nbuilding and rendering templates...')
        files_content = defaultdict(list)
        ExecCommands(g_static_commands, files_content, settings, True)
        ExecCommands(g_instance_commands, files_content, settings, True)
        sys.stderr.write(' done\nwriting rendered templates...\n')
        WriteFiles(files_content)

        # Now write the variable bindings/assignments table.
        sys.stderr.write('writing \"ttree_assignments.txt\" file...')
        open('ttree_assignments.txt', 'w').close() # <-- erase previous version.
        WriteVarBindingsFile(g_objectdefs)
        WriteVarBindingsFile(g_objects)
        sys.stderr.write(' done\n')

    except (ValueError, InputError) as err:
        sys.stderr.write('\n\n'+str(err)+'\n')
        sys.exit(-1)



if __name__ == "__main__":
    main()
