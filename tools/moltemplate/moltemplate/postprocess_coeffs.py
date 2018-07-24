#!/usr/bin/env python

man_page_text = """
Usage (example):

postprocess_coeffs.py ttree_assignments.txt < file.template > file_new.template

Moltemplate users would like to be able to use commands with wildcards like:
  bond_coeff @bond:*/A/* harmonic 533.7 120.0
  pair_coeff @atom:A* @atom:B* lj/cut 0.35 2.7
as shorthand for specifying force field parameters for multiple atom&bond types:
  bond_coeff @bond:a/A/c harmonic 533.7 120.0
  bond_coeff @bond:b/A/d harmonic 533.7 120.0
  pair_coeff @atom:A* @atom:B* lj/cut 0.35 2.7
  pair_coeff @atom:A  @atom:B  lj/cut 0.35 2.7
  pair_coeff @atom:A1 @atom:B1 lj/cut 0.35 2.7
  pair_coeff @atom:A1 @atom:B2 lj/cut 0.35 2.7
  pair_coeff @atom:A2 @atom:B1 lj/cut 0.35 2.7
  pair_coeff @atom:A2 @atom:B2 lj/cut 0.35 2.7
      :         :        :       :     :    :
However LAMMPS does not interpret the * character this way.
Hence, this script will replace the line with the * wildcards above with the 
lines text that follow above.  This script also works for bond_coeff, 
angle_coeff, dihedral_coeff, and improper_coeff commands (expanding wildcard
characters appearing in @bond, @angle, @dihedral, and @improper variables).

This program expects an argument (ttree_assignments.txt) which has a list of all
of the atom types (and bond_types, and angle_types...) which have been defined
and it will substitute those variable names in place of the wildcard expressions
(The ttree_assignments.txt file is a 2-column file variables in the 1st column, 
 and their integer values in the 2nd column.  The 2nd column is ignored.)

The resulting list of commands with explicit variable names will be 
printed to the standard-out.

"""


import sys
import gc

try:
    from .ttree import ExtractFormattingCommands
    from .ttree_lex import *

except (ImportError, SystemError, ValueError):
    # not installed as a package
    from ttree import ExtractFormattingCommands
    from ttree_lex import *


g_filename = __file__.split('/')[-1]
g_module_name = g_filename
if g_filename.rfind('.py') != -1:
    g_module_name = g_filename[:g_filename.rfind('.py')]
g_date_str = '2018-6-09'
g_version_str = '0.2.0'
g_program_name = g_filename
#sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+' ')


def ExtractVarName(text):
    """ Read a string like 'atom:A  '  or  '{/atom:A B/C/../D }ABC '
        and return ('','@atom:A','  ')  or  ('{','atom:A B/C/../D ','}ABC')
        These are 3-tuples containing the portion of the text containing 
        only the variable's name (assumed to be within the text),
        ...in addition to the text on either side of the variable name.
    """
    i_begin = 0
    escape = '\''
    lparen = '{'
    rparen = '}'
    escaped = False
    commenters = '#'
    whitespace = ' \t\r\f\n'
    terminators = whitespace + commenters
    # Ideally, perhaps I should lookup these values from ttree_lex.TtreeLex to 
    # make sure I am being consistent, instead of re-defining them in this file.
    #while ((i_begin < len(text)) and
    #       (text[i_begin] in whitespace)):
    #    i_begin += 1
    in_paren = text[i_begin:i_begin+1] == lparen
    if in_paren:
        terminators = rparen
        i_begin += 1
    i_end = i_begin
    while ((i_end < len(text)) and
           (text[i_end] not in terminators)):
        i_end += 1
    return (text[0: i_begin],
            text[i_begin:i_end],
            text[i_end:])



def main():
    try:
        if (len(sys.argv) != 2):
            raise InputError('Error running  \"' + g_program_name + '\"\n'
                             ' Typical usage:\n'
                             ' postprocess_coeffs.py ttree_assignments.txt < file.template > file.rendered\n'
                             '\n'
                             '   Missing argument.\n'
                             '   Expected the name of a 2-column file containing\n'
                             '   variable names and their bindings (values).\n'
                             '   (This is likely a programmer error.\n'
                             '    This script was not intended to be run by end users.)\n')

        bindings_filename = sys.argv[1]
        f = open(bindings_filename)
        atom_types = set([])
        bond_types = set([])
        angle_types = set([])
        dihedral_types = set([])
        improper_types = set([])

        #BasicUIReadBindingsStream(assignments, f, bindings_filename)

        # The line above is robust but it uses far too much memory.
        # This for loop below works for most cases.
        for line in f:
            #tokens = lines.strip().split()
            # like split but handles quotes
            tokens = SplitQuotedString(line.strip())
            if len(tokens) < 2:
                continue
            if tokens[0].find('@') != 0:
                continue
            if tokens[0][2:].find('atom') == 0:
                atom_types.add(tokens[0][1:])
            elif tokens[0][2:].find('bond') == 0:
                bond_types.add(tokens[0][1:])
            elif tokens[0][2:].find('angle') == 0:
                angle_types.add(tokens[0][1:])
            elif tokens[0][2:].find('dihedral') == 0:
                dihedral_types.add(tokens[0][1:])
            elif tokens[0][2:].find('improper') == 0:
                improper_types.add(tokens[0][1:])

        f.close()
        gc.collect()

        lex = LineLex(sys.stdin, '__standard_input_for_postprocess_coeffs__')
        #lex = LineLex(open('deleteme.template', 'r'), '__standard_input_for_postprocess_coeffs_')
        lex.commenters = ''            #(don't attempt to skip over comments)
        lex.line_extend_chars += '&'   #(because LAMMPS interprets '&' as '\')

        while True:
            line_orig = lex.ReadLine()
            #sys.stderr.write('line_orig = \"'+str(line_orig)+'\"\n')
            if (not line_orig) or (line_orig == ''):
                break
            tokens = line_orig.strip().split('@')

            if ((len(tokens) >= 2) and
                (tokens[0].find('bond_coeff') == 0) and
                #does this token contain '*' or '?'
                HasWildcard(tokens[1]) and
                (tokens[1][0:1] != '{')):  #(don't pattern match * in {})?'
                left_paren, typepattern, text_after = ExtractVarName(tokens[1])
                for btype in bond_types:
                    if MatchesPattern(btype, typepattern):
                        #assert(left_paren == '')
                        tokens[1] = btype + text_after
                        sys.stdout.write('@'.join(tokens) + '\n')

            elif ((len(tokens) >= 2) and
                  (tokens[0].find('angle_coeff') == 0) and
                  #does this token contain '*' or '?'
                  HasWildcard(tokens[1]) and
                  (tokens[1][0:1] != '{')):
                left_paren, typepattern, text_after = ExtractVarName(tokens[1])
                for antype in angle_types:
                    if MatchesPattern(antype, typepattern):
                        #assert(left_paren == '')
                        tokens[1] = antype + text_after
                        sys.stdout.write('@'.join(tokens) + '\n')

            elif ((len(tokens) >= 2) and
                  (tokens[0].find('dihedral_coeff') == 0) and
                  #does this token contain '*' or '?'
                  HasWildcard(tokens[1]) and
                  (tokens[1][0:1] != '{')):  #(don't pattern match * in {})
                left_paren, typepattern, text_after = ExtractVarName(tokens[1])
                for dtype in dihedral_types:
                    if MatchesPattern(dtype, typepattern):
                        #assert(left_paren == '')
                        tokens[1] = dtype + text_after
                        sys.stdout.write('@'.join(tokens) + '\n')

            elif ((len(tokens) >= 2) and
                  (tokens[0].find('improper_coeff') == 0) and
                  #does this token contain '*' or '?'
                  HasWildcard(tokens[1]) and
                  (tokens[1][0:1] != '{')):  #(don't pattern match * in {})
                left_paren, typepattern, text_after = ExtractVarName(tokens[1])
                for itype in improper_types:
                    if MatchesPattern(itype, typepattern):
                        #assert(left_paren == '')
                        tokens[1] = itype + text_after
                        sys.stdout.write('@'.join(tokens) + '\n')

            #elif ((len(tokens) >= 3) and
            #      (tokens[0].find('pair_coeff') == 0) and
            #      (HasWildcard(tokens[1]) or HasWildcard(tokens[2]))):
            elif ((len(tokens) >= 2) and
                  (tokens[0].find('pair_coeff') == 0)):
                # First deal with cases with only one @variable, such as:
                #    pair_coeff @atom:A*   *       ...
                #    pair_coeff    *     @atom:A*  ...
                # (We don't deal with cases like "* *" because LAMMPS interprets
                #  these in a special way:  manybody pair_styles use "* *")
                if len(tokens) == 2:
                    if tokens[0].rstrip()[-1:] == '*':
                        tokens[0] = tokens[0].rstrip()[:-1]
                        tokens.insert(1, '/atom:* ')
                    else:
                        ic = tokens[1].find(' * ')
                        tokens.append('/atom:* '+tokens[1][ic+3:])
                        tokens[1] = tokens[1][:ic]+' '
               
                assert(len(tokens) >= 3)
                left_paren1,typepattern1,text_after1=ExtractVarName(tokens[1])

                # Then deal with cases like this:
                #  pair_coeff @{/atom:r1}*@{/atom:r3} @{/atom:r4}*@{/atom:r6}
                # In this case we should be using ' ' as the delimeter, not '@'
                # to separate the two arguments from eachother, since
                #   @{/atom:r1}*@{/atom:r3} is the first argument, and
                #   @{/atom:r4}*@{/atom:r6} is the second argument

                # Check: Were there any whitespace characters in the text
                #        separating token[1] from token[2]?
                if ((left_paren1 == '{') and
                    (len(SplitQuotedString(text_after1)) == 1)):
                    # If not, then tokens[1] and tokens[2] are both part of
                    # the 1st argument.
                    tokens[1] = tokens[1]+'@'+tokens[2]
                    left_paren1 = ''
                    text_after1 = ''
                    typepattern1 = tokens[1]
                    del tokens[2]

                left_paren2,typepattern2,text_after2=ExtractVarName(tokens[2])
                # Check: Were there any whitespace characters in the text
                #        separating token[2] from what follows?
                if ((len(tokens) > 3) and
                    (len(SplitQuotedString(text_after2)) == 1)):
                    # If not, then tokens[2] and tokens[3] are both part of
                    # the 2nd argument.
                    tokens[2] = tokens[2]+'@'+tokens[3]
                    left_paren2 = ''
                    text_after2 = ''
                    typepattern2 = tokens[2]
                    del tokens[3]
                if (HasWildcard(tokens[1]) and
                    (tokens[1][0:1] != '{')):  #(don't pattern match * in {})
                    atom_types1 = atom_types
                else:
                    atom_types1 = set([typepattern1])
                if (HasWildcard(tokens[2]) and
                    (tokens[2][0:1] != '{')):  #(don't pattern match * in {})
                    atom_types2 = atom_types
                else:
                    atom_types2 = set([typepattern2])
                for atype1 in atom_types1:
                    #sys.stderr.write('atype1 = \"'+str(atype1)+'\"\n')
                    if MatchesPattern(atype1, typepattern1):
                        #assert(left_paren1 == '')
                        tokens[1] = left_paren1 + atype1 + text_after1
                        for atype2 in atom_types2:
                            #sys.stderr.write(' atype2 = \"'+str(atype2)+'\"\n')
                            if MatchesPattern(atype2, typepattern2):
                                #assert(left_paren2 == '')
                                tokens[2] = left_paren2 + atype2 + text_after2
                                sys.stdout.write('@'.join(tokens) + '\n')
            else:
                sys.stdout.write(line_orig)

    except (ValueError, InputError) as err:
        sys.stderr.write('\n' + str(err) + '\n')
        sys.exit(-1)

    return

if __name__ == '__main__':
    main()
