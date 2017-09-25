

espt_delim_atom_fields = set(["pos", "type", "v", "f", 
                              "bond", 
                              "temp", "gamma", 
                              "q", 
                              "quat", "omega", "torque", 
                              "rinertia", "fix", "unfix", "ext_force",
                              "exclude", "delete", 
                              "mass",
                              "dipm", "dip", "virtual",
                              "vs_relative", "distance", "vs_auto_relate_to"])



def LinesWSlashes(text):
    """ 
    Iterate over the lines contained in a string of text.
    Merge lines ending in backslashes.

    """

    current_line = ''
    for line in text.split('\n'):
        current_line += line
        if (len(line) > 0) and (line[-1] != '\\'):
            yield current_line
            current_line = ''
    if len(current_line) > 0:
        yield current_line




def SplitMultiDelims(line, delimiters):
    """
    Split a string into tokens using one or more (multi-character) delimiters.
    (Bug: The current version of this function does not preserve white space,
          but this should not matter.)

    """

    token = ''
    for sub_token in line.strip().split():
        if sub_token in delimiters:
            yield token
            yield sub_token
            token = ''
        elif len(token) > 0:
            token += ' ' + sub_token
        else:
            token += sub_token
    if len(token) > 0:
        yield token



def SplitAtomLine(line):
    l = []
    for token in SplitMultiDelims(line, espt_delim_atom_fields):
        l.append(token)
    return l

    # In this type of TCL command, all of the delimiters 
    #    (like 'pos', 'type', 'q', ...) 
    # are supposed to be followed by an argument.  If the last 
    # token on this line IS a delimiter, then this is a syntax error.

    if token in espt_delim_atom_fields:
        raise InputError("Error: Incomplete line:\n"
                         "\""+line+"\"\n")



def iEsptAtomCoords(tokens):
    #tokens = SplitMultiDelims(line)
    i = 0
    while i < len(tokens):
        if tokens[i] in set(['pos', 'fix', 'unfix']):
            assert(i+1 < len(tokens))
            yield i+1
            i += 1
        i += 1




def iEsptAtomVects(tokens):
    #tokens = SplitMultiDelims(line)
    i = 0
    while i < len(tokens):
        if tokens[i] in set(['dip', 'rinertia', 'v', 'f', 'omega', 'torque']):
            assert(i+1 < len(tokens))
            yield i+1
            i += 1
        i += 1


def iEsptAtomType(tokens):
    #tokens = SplitMultiDelims(line)
    i = 0
    while i < len(tokens):
        if tokens[i] == 'type':
            assert(i+1 < len(tokens))
            yield i+1
            i += 1
        i += 1

def iEsptAtomID(tokens):
    if len(tokens) > 1:
        return 1
    else:
        raise InputError("Error: Incomplete line:\n"
                         "\""+line+"\"\n")

