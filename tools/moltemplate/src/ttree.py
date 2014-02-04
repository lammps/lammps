#!/usr/bin/env python

# Authors: Andrew Jewett (jewett.aij at g mail)
#         http://www.chem.ucsb.edu/~sheagroup
# License: 3-clause BSD License  (See LICENSE.TXT)
# Copyright (c) 2011, Regents of the University of California
# All rights reserved.

"""
ttree Ttree is a simple program for recursively composing and generating  
      large redundant text files from small template files.

      By default, the large number of unique template variables generated 
      in the process are automatically substituted with integers 
      (or other numeric counters, all of which can be overridden),
      rendered, and the rendered templates are written to a file.  

      ttree was designed to be useful for generating input files for
      molecular simulation software like LAMMPS or NAMD.

BasicUI  This section of the code contains the user interface for ttree
         when run as a stand-alone program, as described above.  (This
         section of code contains the "if __name__ == __main__:" code block.)

-- Data Types --

StaticObj   Static nodes are data structures used to store ttree class definitions.
                (Static nodes are useful for defining molecule types or 
                 namespaces in LAMMPS or other molecular simulation programs.)
         The nodes themselves are stored in a tree of nested class definitions.
         Static variables (such as "@atom:C") are also associated with 
         StaticObjs.

InstanceObj   Instance nodes are created when a user creates one (or many) 
         copies of a class, using the "new" command.
         These classes in turn may instantiate other classes.
            (Example: A user may manually instantiate several copies of a
                      molecule, such as a protein, however each of those 
                      molecules may contain molecular subunits, such as
                      amino acids, which are automatically instantiated.)
         Instance variables (such as "$atom:CA") are also associated with 
         InstanceObjs.

"""

import sys
from collections import defaultdict
import operator
import random
#import gc

try:
    unicode
except NameError:
    # Python 3
    basestring = unicode = str


#   -- ttree_lex.py --
# TtreeShlex is a backwards-compatible version of python's standard shlex module.
# It has the additional member: "self.wordterminators", which overrides 
# the "self.wordchars" member.  This enables better handling of unicode 
# characters by allowing a much larger variety of characters to appear 
# in words or tokens parsed by TtreeShlex.  Otherwise it is identical to shlex.
from ttree_lex import *


if sys.version < '2.6':
    raise InputError('Error: Using python '+sys.version+'\n'
                     '       Alas, you must upgrade to a newer version of python (2.7 or later).')
elif sys.version < '2.7':
    sys.stderr.write('--------------------------------------------------------\n'
                     '----------------- WARNING: OLD PYTHON VERSION ----------\n'
                     '  This program is untested on your python version ('+sys.version+').\n'
                     '  PLEASE LET ME KNOW IF THIS PROGRAM CRASHES (and upgrade python).\n'
                     '    -Andrew   2013-10-25\n'
                     '--------------------------------------------------------\n'
                     '--------------------------------------------------------\n')
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict


if sys.version > '3':
    import io
else:
    import cStringIO


# We keep track of the program name and version.
# (This is only used for generating error messages.)
#g_filename = 'ttree.py'
g_filename    = __file__.split('/')[-1]
g_module_name  = g_filename
if g_filename.rfind('.py') != -1:
    g_module_name = g_filename[:g_filename.rfind('.py')]
g_date_str     = '2013-9-12'
g_version_str  = '0.76'







class ClassReference(object):
    """ Every class defined by the user (stored in an StaticObj data structure) 
    may contain references to other classes (ie. other StaticObjs).
    (Note: All of these StaticObjs are stored in the same tree, the
     global static tree.)
    Examples:
       Whenever an instance of a class is created, this may automatically spawn 
    the creation of additional classes (which are instantiated because a 'new' 
    command appeared within the first class's definition).  These are stored in
    the "StaticObj.instance_commands[i].class_ref" attribute.
       Similarly, each class (StaticObj) can optionally inherit some of its 
    traits (consisting of write() and new commands) from one or more 
    "class_parents" (also StaticObjs).  A list of these parents is stored in the 
    "StaticObj.class_parents" attribute.  In both cases (self.instance_commands 
    or self.class_parents) we need to storea pointer to the StaticObj(s) 
    corresponding to the instance-childen or class-parents.
    (This stored in self.statobj).
    However, for the purposes of debugging and interactivity, it is also 
    convenient to permanently keep track of the string that the user used to 
    specify the name/location of that class/StaticObj 
    (stored in self.statobj_str), in addition to the location 
    in the file where that string occurs (stored in self.srcloc)."""

    __slots__=["statobj_str","srcloc","statobj"]

    def __init__(self, 
                 statobj_str=None, 
                 srcloc=None, 
                 statobj=None):
        self.statobj_str  = statobj_str
        if srcloc is None:
            self.srcloc = OSrcLoc('', -1)
        else:
            self.srcloc = srcloc
        self.statobj = statobj

    #def __repr__(self):
    #    return repr((self.statobj_str, self.srcloc))



# "Command"s are tasks to carry out.
# (...either immediately, or later during instantiation)
# Commands are used to write to files, create new instances, delete instances,
# or custom commands to modify an instance of a class.
# (For example "instance = new Class.move(1,0,0).rot(45,0,0,1)"
#  The ".move(1,0,0)" and ".rot(45,0,0,1)" commands are "stackable" and 
#  have similar syntax to member functions in C++, JAVA, Python.)

class Command(object):
    __slots__=["srcloc"]

    def __init__(self, srcloc=None):
        self.srcloc = srcloc

    # COMMENTING OUT: "COUNT" AND "ORDER" ARE NO LONGER NEEDED

    #count = 0
    #def __init__(self, srcloc=None):
    #    self.srcloc = srcloc
    #    # The "order" member is a counter that keeps track of the order 
    #    #  in which the Command data types are created (issued by the user).
    #    Command.count += 1
    #    self.order = Command.count
    #def __lt__(self, x):
    #    return self.order < x.order



class WriteFileCommand(Command):
    """ WriteFileCommand

    filename  This is the name of the file that will be written to
               when the command is executed.
    tmpl_list  This is the contents of what will be written to the file.
               Text strings are often simple strings, however more
               generally, they can be strings which include other variables
               (ie templates).  In general, templates are lists of alternating
               TextBlocks and VarRefs, (with additional tags and data to 
               identify where they occur in in the original user's files).
               
    """
    __slots__=["filename", "tmpl_list"]

    def __init__(self, 
                 filename = None,
                 tmpl_list = None,
                 srcloc = None):
        self.filename = filename
        if tmpl_list is None:
            self.tmpl_list = []
        else:
            Command.__init__(self, srcloc)
            self.tmpl_list = tmpl_list
    def __str__(self):
        if self.filename:
            return 'WriteFileCommand(\"'+self.filename+'\")'
        else:
            return 'WriteFileCommand(NULL)'
    def __copy__(self):
        tmpl_list = []
        CopyTmplList(self.tmpl_list, tmpl_list) #CHECK:IS_MEMORY_WASTED_HERE?
        return WriteFileCommand(self.filename, tmpl_list, self.srcloc)


class InstantiateCommand(Command):
    """ InstantiateCommand is a simple tuple-like datatype used to 
    store pairs of names (strings, stored in self.name), 
    and ClassReferences (see above, stored in self.class_ref).n
    The "suffix" argument is an optional string which may contain 
    additional instructions how to instantiate the object.

    """

    __slots__=["name", "class_ref"]

    def __init__(self,
                 name = None,
                 class_ref = None,
                 srcloc = None):
        Command.__init__(self, srcloc)
        self.name = name
        #if class_ref is None:
        #    self.class_ref = ClassReference()
        #else:
        self.class_ref = class_ref

    def __str__(self):
        return 'InstantiateCommand('+self.name+')'

    def __copy__(self):
        return InstantiateCommand(self.name, self.class_ref, self.srcloc)


class DeleteCommand(Command):
    __slots__=[]
    def __init__(self,
                 srcloc = None):
        Command.__init__(self, srcloc)

    def __str__(self):
        return 'DeleteCommand()'

    def __copy__(self):
        return DeleteCommand(self.srcloc)




class StackableCommand(Command):
    """ StackableCommand is a class for storing commands
    that effect the environment of the object being created.
    The combined effect of these commands can be thought of as a "stack"
    Commands can be pushed on the stack, or popped off

    The actual commands themselves are represented by the "contents" member
    which is usually a text string.
    ttree.py does not attempt to understand the content of these commands.
    That job is left up to the __main___ module.  (IE. whatever script that
    happens to be importing ttree.py.  If there is no script, and
    ttree.py IS the main module, then it simply ignores these commands.)

    """

    __slots__=["context_node"]

    def __init__(self,
                 srcloc,
                 context_node=None):
        Command.__init__(self, srcloc)
        self.context_node = context_node # if multiple stacks are present, then use "context_node"
                             # as a key to identify which stack you want
                             # the command to modify

class PushCommand(StackableCommand):

    __slots__=["contents"]

    def __init__(self,
                 contents,
                 srcloc,
                 context_node=None):
        StackableCommand.__init__(self, srcloc, context_node)
        self.contents = contents

    def __copy__(self):
        return PushCommand(self.contents, self.srcloc, self.context_node)

    def __str__(self):
        return 'PushCommand('+str(self.contents)+')'

class PushRightCommand(PushCommand):
    __slots__=[]
    def __init__(self,
                 contents,
                 srcloc,
                 context_node=None):
        PushCommand.__init__(self, contents, srcloc, context_node)

    def __copy__(self):
        return PushRightCommand(self.contents, self.srcloc, self.context_node)

    def __str__(self):
        return 'PushRightCommand('+str(self.contents)+')'

class PushLeftCommand(PushCommand):
    __slots__=[]
    def __init__(self,
                 contents,
                 srcloc,
                 context_node=None):
        PushCommand.__init__(self, contents, srcloc, context_node)

    def __copy__(self):
        return PushLeftCommand(self.contents, self.srcloc, self.context_node)

    def __str__(self):
        return 'PushLeftCommand('+str(self.contents)+')'

class PopCommand(StackableCommand):
    __slots__=["partner"]
    def __init__(self,
                 partner,
                 srcloc,
                 context_node=None):
        StackableCommand.__init__(self, srcloc, context_node)
        self.partner = partner

    def __copy__(self):
        return PopCommand(self.partner, self.srcloc, self.context_node)

    def __str__(self):
        return 'PopCommand('+str(self.partner.contents)+')'

class PopRightCommand(PopCommand):
    __slots__=[]
    def __init__(self,
                 partner,
                 srcloc,
                 context_node=None):
        PopCommand.__init__(self, partner, srcloc, context_node)
        assert((partner is None) or isinstance(partner, PushRightCommand))

    def __copy__(self):
        return PopRightCommand(self.partner, self.srcloc, self.context_node)

    def __str__(self):
        return 'PopRightCommand('+str(self.partner.contents)+')'

class PopLeftCommand(PopCommand):
    __slots__=[]
    def __init__(self,
                 partner,
                 srcloc,
                 context_node=None):
        PopCommand.__init__(self, partner, srcloc, context_node)
        assert((partner is None) or isinstance(partner, PushLeftCommand))

    def __copy__(self):
        return PopLeftCommand(self.partner, self.srcloc, self.context_node)

    def __str__(self):
        return 'PopLeftCommand('+str(self.partner.contents)+')'





# The ScopeCommand, ScopeBegin, and ScopeEnd commands are useful to designate 
# which commands belong to a particular class definition (or class instance).
# (This is useful later on, when a linear list of commands has been created.)
# They are simply markers an do not do anything.  These classes can be ignored.
class ScopeCommand(Command):
    __slots__=["node"]

    def __init__(self,
                 node,
                 srcloc):
        Command.__init__(self, srcloc)
        self.node = node
        #self.srcloc = srcloc

    def __copy__(self):
        return ScopeCommand(self.node, self.srcloc)

    def __str__(self):
        if self.node:
            return 'ScopeCommand('+self.node.name+')'
        else:
            return 'ScopeCommand(None)'

class ScopeBegin(ScopeCommand):
    __slots__=[]
    def __init__(self, node, srcloc):
        ScopeCommand.__init__(self, node, srcloc)
    def __copy__(self):
        return ScopeBegin(self.node, self.srcloc)
    def __str__(self):
        if self.node:
            return 'ScopeBegin('+NodeToStr(self.node)+')'
        else:
            return 'ScopeBegin(None)'

class ScopeEnd(ScopeCommand):
    __slots__=[]
    def __init__(self, node, srcloc):
        ScopeCommand.__init__(self, node, srcloc)
    def __copy__(self):
        return ScopeEnd(self.node, self.srcloc)
    def __str__(self):
        if self.node:
            return 'ScopeEnd('+NodeToStr(self.node)+')'
        else:
            return 'ScopeEnd(None)'



# COMMENTING OUT: NOT NEEDED AT THE MOMENT
#class VarAssignCommand(Command):
#    """ VarAssignCommand
#
#    This class is used whenever the user makes an explicit request to assign
#    a variable to a value (values are text strings).
#
#    var_ref    The variable name (tecnically speaking, I call this
#                 a variable descriptor string and it includes at least one of 
#                 the following: the name of a leaf node, a category node name,
#                 and category name)
#               the location in the file where variable appears, and (eventually
#               after subsequent lookup), references to the leaf_node, cat_node,
#               "Category", and "VarBinding" data structures associated with it.
#    text_tmpl  Text strings are often simple strings, however more
#               generally, they can be strings which include other variables
#               (ie templates).  In general, templates are lists of alternating
#               TextBlocks and VarRefs, (with additional tags and data to 
#               identify where they occur in in the original user's files).
#               
#    """
#    __slots__=["var_ref","text_tmpl"]
#
#    def __init__(self,
#                 #command_name = '=',  <-- ?!?
#                 var_ref = None,
#                 text_tmpl=None):
#        Command.__init__(self, srcloc)
#        self.var_ref = var_ref
#        self.text_tmpl = text_tmpl



class ModCommand(object):
    __slots__=["command","multi_descr_str"]
    def __init__(self,
                 command,
                 multi_descr_str):
        self.command = command
        self.multi_descr_str = multi_descr_str

    def __str__(self):
        return 'ModCommand('+str(self.command)+')'

    def __copy__(self):
        return ModCommand(self.command.__copy__(), self.multi_descr_str)
        



def CopyTmplList(source_tmpl_list, dest_cpy):
    for entry in source_tmpl_list:
        if isinstance(entry, TextBlock):
            dest_cpy.append(entry) # Then make a shallow copy
                                          # (pointer assignment) to the text
                                          # block (Text blocks do not change
                                          # during instantiation.)
        elif isinstance(entry, VarRef):
            assert(len(entry.prefix)>0)
            if entry.prefix[0] == '@': # '@' vars refer to static data
                dest_cpy.append(entry) # Then make a shallow copy
                                          # pointer assignment) to the static
                                          # variable. (Static variables do
                                          # not change during instantiation.)

            elif entry.prefix[0] == '$': # new '$' vars are created 
                                         # during every instantiation.

                # var_refs do change when you instantiate them.  So 
                # create a new VarRef object, and copy the attributes.
                var_ref = VarRef(entry.prefix,
                                 entry.descr_str,
                                 entry.suffix,
                                 entry.srcloc)
                                 # Note: for instance variables ('$' vars)
                                 #       "entry.nptr" should not contain
                                 #       any data yet, so we just ignore it.
                                 # I assert this below:
                assert((entry.nptr.cat_node is None) and
                       (entry.nptr.leaf_node is None))

                dest_cpy.append(var_ref)
            else:
                assert(False) # prefix[0] should be either '@' or '$'
        else:
            assert(False) # type(entry) should be either TextBlock or VarRef







def RecursiveJoin(tokens_expr, delimiter = ''):
    """ RecursiveJoin() converts a tree-like list/tuple of tokens, for example:
    ['a ', ('tree', '-', ['like', 'container']), [[' '], 'of'], ' strings']
    to an ordinary string, eg: 
    'a tree-like container of strings'
    This behavees similarly to  "reduce(lambda a, b: a+b, tokens)", 
    except that it works with arbitrarily nested lists/tuples."""
    text = ''
    if isinstance(tokens_expr, basestring):
        return tokens_expr
    else:
        text_lstr = []
        for i in range(0, len(tokens_expr)):
            text.append( TokensToStr(tokens_expr[i]) )
        return ''.join(text_lstr, delimiter)






#----------------------------------------------------------
#----------------------------------------------------------
#    The following code is specific to ttree.
#
#    (Up until this point, we have only defined 
#    a few simple general text parsing routines.)
#----------------------------------------------------------
#----------------------------------------------------------



def PtknsToStr(path_tokens):
    """ 
    There are three ways to store paths:
      As a single string: '/Protein/Phe/Ca'  <- the format entered by the user
      As a list of tokens ['Protein', 'Phe', 'Ca'] <- split into tokens
      As a list of nodes in a tree (pointers to nodes in a tree hierarchy)
    This function converts between the first two formats.
    """
    text = ''
    if len(path_tokens) > 0:
        text = path_tokens[0]
        for i in range(1, len(path_tokens)):
            text += '/' + path_tokens[i]
    else:
        text = ''
    return text



def StrToPtkns(path_string):
    """ The inverse of PtknsToStr(), this function splits a string like
        '/usr/local/../bin/awk' into ['usr','local','..','bin','awk'].
        For illustrative purposes only. Use text.split('/') directly instead."""
    return orig_text.split('/')



def FindChild(name, node, dbg_loc):
    """ FindChild looks over the list of node.children to find a child
    which matches the name given in the first argument.
    If it is not found, it returns None.
       Note: I have not yet specified what kind of nodes FindChild() operates
    on.  Both StaticObjs and InstanceObjs have self.children and self.parent.
    However only StaticObjs have "self.class_parents". 
    ("class_parents" are "parents" in the object-oriented sense.)
    If "node" (2nd argument) happens t be an StaticObj, this means it also
    We must search over the children of these class_parents as well.

          Terminology used here differs from Object Oriented Programming

    Children in node.children are not children in the object-oriented 
    programming sense.  However, in OOP, "children" are objects that share all
    of the traits of their ancestors (and may have additionl traits as well).
    I have implemented OOP style children and parents, but this informtion
    is stored in "node.class_parents", instead of "node.parents".
       For comparison, instantiated nodes (InstanceObjs) are different.  Altough 
    instantiated classes (InstanceObjs) have access to the attributes of the 
    class_parents of the StaticObjs that define them, they do not remember the 
    ownership of that data.  (It just gets merged with their own member data,
    including their .children.)
       Hence we must treat StaticObjs carefully because their are two ways we can
    access child data.  We should loop over both of them.  We do that below:
    """
    for child in node.children:
        if child.name == name:
            return child

    # The object-oriented inheritance stuff appears here.
    # If you don't care about OOP or inheritance, 
    # then comment out the loop that follows:

    if isinstance(node, StaticObj):
        #CONTINUEHERE: DEBUG THIS FEATURE AND EXTEND NAMESPACES TO INSTOBJS
        # Search recursively over the "children" (ie attributes or members)
        # belonging to any OOP ancestors of this node.
        for class_parent in node.class_parents:
            for child in class_parent.children:
                child = FindChild(name, class_parent, dbg_loc)
                if child != None:
                    return child
        for namespace_node in node.namespaces:
            for child in namespace_node.children:
                child = FindChild(name, namespace_node, dbg_loc)
                if child != None:
                    return child

    # Otherwise, a child name match was not found
    return None

    


def FollowPath(path_tokens, starting_node, dbg_loc):
    """ FollowPath() returns the "last_node", a node whose position in the
        tree is indicated by a list of path_tokens, describing the names
        of nodes connecting "starting_node" to "last_node".
        If it one of the strings in the list of path_tokens turns out
        not to match then names of classes in the tree, then this function
        returns the last_node that did match before the error occurred,
        as well as an integer which stores the number of tokens in
        the path_tokens list which were successfully processed.
        In other words, the list of node naes is not a full path, but the 
        relative path that takes you from one node (not necessarily the root) 
        to another.  Return Value:
           Ideally, each node in the list should be a parent or a child of the 
        previous node. (See comment for PathTokensToStr(), for more details.)
        This function returns the number of path_tokens successfully 
        parsed.  Under normal termination, this is len(path_tokens).
        If the path can not be followed (because at some point, a child 
        or parent does not exist), then this function returns a number
        smaller than len(path_tokens).
        We let the caller handle undefined paths. """

    #print('    FollowPath() invoked on: ', path_tokens)

    if len(path_tokens) == 0:
        return 0, starting_node

    node = starting_node
    # Is this path a relative path, or a full path?
    # If the path-string began with '/', then it's a full path.  This means 
    # that after processing by split('/'), the first token will be ''
    # Example:       path_tokens='/Prot/Alanine'.split('/')
    #            --> path_tokens[0] == ''
    if path_tokens[0] == '':
        # In that case, then take us to the root node:
        while node.parent != None:
            node = node.parent
            #sys.stdout.write('FollowPath(): Retreating to node \"'+node.name+'\"\n')
        i0 = 1 # <- We've just processed the first token.  Skip over it later.
    else:
        i0 = 0


    i = i0
    while i < len(path_tokens):

        if path_tokens[i] == '..':
            if node.parent is None:
                return i, node # <-return the index into the token list
                               #   Caller will know that something went awry
                               #   if the return value is not equal to the
                               #   length of the token list
            else:
                node = node.parent
            i += 1

        elif path_tokens[i] == '...':

            node_before_ellipsis = node
            if i == len(path_tokens)-1:
                return i, node_before_ellipsis

            search_target = path_tokens[i+1]
            # Now search over the "children" of this node 
            # for one who's name matches path_tokens[i]. 
            # If not found, then move up to the parent node's children. 
            # (This is not an exhaustive tree search. Only the nodes which 
            #  are immediate children of this node's parents are searched.)
            while node != None:
                child = FindChild(search_target, node, dbg_loc)
                if child is None:
                    node = node.parent 
                else:
                    node = child
                    break

            if node is None:
                #   Caller will know that something went awry if the return
                #   value is not equal to the length of the token list.
                return i, node_before_ellipsis

            i += 2

        elif path_tokens[i] in ('','.'):  # <-Note we ignore empty tokens from now on.
            # (Same convention is used in specifying a 
            # directory in a filesystem, eg. using /usr/local
            # or /usr//local or /usr/./local.  These are equivalent.)
            i += 1

        else:

            # Now search over the "children" of this 
            # node for one who's name matches path_tokens[i].
            child = FindChild(path_tokens[i], node, dbg_loc)

            if child is None:
                # In that case, return with the node_list incomplete.
                # Let the caller check to see if something went wrong.
                return i, node # <-return the index into the token list (i)
                          #   Caller will know that something went awry
                          #   if the return value is not equal to the
                          #   length of the token list
            else:
                node = child

            i += 1

            if node.IsDeleted():
                #sys.stderr.write('(debug_msg: encountered deleted node: \"'+node.name+'\")\n')
                break

    return len(path_tokens), node





def PtknsToNode(path_tokens, starting_node, dbg_loc):
    """ PtknsToNode() is identical to def FollowPath() except
    that it raises syntax-error exceptions if the path is undefined."""

    i_last_ptkn, last_node = FollowPath(path_tokens, starting_node, dbg_loc)

    if i_last_ptkn < len(path_tokens):
        #assert(isinstance(last_node,StaticObj)) <--why did I assert this? seems wrong

        if (last_node.parent is None) and (path_tokens[i_last_ptkn] == '..'):
            #In that case, we tried to back out beyond the root of the tree.
            raise InputError('Error('+g_module_name+'.PtknsToNode()):\n'
                             '        Invalid variable/class name:\n'
                             '        \"'+PtknsToStr(path_tokens)+'\" located near '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno)+'\n'
                             '       There are too many \"..\" tokens in the path string.')

        elif path_tokens[i_last_ptkn] == '...':
            if i_last_ptkn+1 == len(path_tokens):
                raise InputError('Error('+g_module_name+'.PtknsToNode()):\n'
                                 '       Error in '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno)+'\n'
                                 '       Expected name following \"...\"\n')
            else:
                search_target = path_tokens[i_last_ptkn+1]
                #In that case, we were unable to find the node referenced by "..." 
                raise InputError('Error('+g_module_name+'.PtknsToNode()):\n'
                                 '       Class or variable \"'+search_target+'\" not found\n'
                                 '       in this context: \"'+PtknsToStr(path_tokens)+'\"\n'
                                 '       located near '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno))

        else:
            #Then the reason is: The string in path_tokens[i_last_ptkn]
            #was supposed to be a child of last_node but a child
            #of that name was not found.
            err_msg = 'Error('+g_module_name+'.PtknsToNode()):\n'+\
                      '    Undefined variable/class name:\n'+\
                      '       \"'+PtknsToStr(path_tokens)+'\",\n'+\
                      '    This occured near or before '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno)+'\n'+\
                      '    (Specifically \"'+path_tokens[i_last_ptkn]+\
                      '\" is not a subordinate of \"'+MaxLenStr(last_node.name,'/')+'\".)\n'+\
                      '    This may be due to a typo located here or earlier.\n'+\
                      '    It may also occur if you deleted the object earlier.  (Referring to a\n'+\
                      '        deleted object is only forgiven when using [0-9] or [0:10] notation.)\n'+\
                      '    If this object refers to an array you must use brackets []\n'+\
                      '    to explicitly specify the element(s) you want from that array.\n'+\
                      '    (To select multiple elements, you can use [*] or [0-9] or [0:10].)\n'
                   
            if (path_tokens[i_last_ptkn] in NodeToPtkns(last_node)):
                err_msg += '\nIn this case:\n'+\
                           '    It seems like you may have omitted a } character somewhere before:\n'+\
                           '    '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno)
            raise InputError(err_msg)
        assert(False) # One of the two conditions above should be true.

    return last_node





def StrToNode(obj_name, starting_node, dbg_loc):
    path_tokens = obj_name.split('/')
    return PtknsToNode(path_tokens, starting_node, dbg_loc)



def NodeListToPtkns(node_list, dbg_loc=None):
    assert(len(node_list) > 0) #The path must contain at least the starting node
    path_tokens = [node_list[0].name]
    for i in range(1, len(node_list)):
        if node_list[i] == node_list[i-1].parent:
            path_tokens.append('..')
        else:
            path_tokens.append(node_list[i].name)
            # Now check to make sure the user supplied consistent information:
            if (node_list[i] not in node_list[i-1].children):
                raise InputError('Error('+g_module_name+'.NodeListToPtkns()):\n'
                                 '        Undefined variable/class name:\n'
                                 '        \"'+PtknsToStr(path_tokens)+'\" located near '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno)+'\n'
                                 '       (\"'+path_tokens[i]+'\" is not subordinate to \"'+MaxLenStr(node_list[i-1].name,'/')+'\")\n'
                                 '       This could be an internal error.')
    return path_tokens



def NodeListToStr(node_list, dbg_loc=None):
    assert(len(node_list) > 0) #The path must contain at least the starting node
    path_str = node_list[0].name
    for i in range(1, len(node_list)):
        if node_list[i] == node_list[i-1].parent:
            path_str += '/..'
        else:
            path_str += '/' + node_list[i].name
            # Now check to make sure the user supplied consistent information:
            if (node_list[i] not in node_list[i-1].children):
                err_msg = 'Error('+g_module_name+'.NodeListToStr()):\n'+\
                          '    Invalid variable/class name:\n' +\
                          '    \"'+PtknsToStr(path_tokens)+'\"'
                if dbg_loc != None:
                    err_msg += ' located near '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno)
                err_msg += '\n' +\
                           '       (\"'+node_list[i].name+'\" is not a subordinate of \"'+MaxLenStr(node_list[i-1].name,'/')+'\")\n'+\
                           '    This could be an internal error.'
                raise InputError(err_msg)
    return path_str



def NodeToPtkns(node):
    ptkns = []
    nd = node
    while nd != None:
        ptkns.append(nd.name)
        nd = nd.parent
    ptkns.reverse()
    return ptkns



def NodeToStr(node):
    ptkns = NodeToPtkns(node)
    assert(len(ptkns) > 0)
    if node.parent is None:
        assert(node.name == '')
        return '/'
    path_str = ptkns[0]
    i = 1
    while i < len(ptkns):
        path_str += '/'+ptkns[i]
        i += 1
    return path_str



def CatLeafNodesToTkns(cat_name, cat_node, leaf_node, dbg_loc):
    assert((cat_node != None) and (leaf_node != None))
    assert((cat_name != None) and  (cat_name != ''))

    # Determine the path of the cat node
    cat_node_ptkns = NodeToPtkns(cat_node)
    cat_node_ptkns.append(cat_name+':')

    # Determine the path of the leaf node (which should inherit from cat)
    deleted = False
    leaf_node_ptkns = []
    if cat_node != leaf_node:
        node = leaf_node
        while node.parent != None:
            if node.IsDeleted():
                deleted = True
                leaf_node_ptkns.append('DELETED_'+node.name)
                break
            leaf_node_ptkns.append(node.name)
            if node.parent == cat_node:
                break
            node = node.parent
        leaf_node_ptkns.reverse()

        if not deleted:
            # Check that leaf inherits from cat.  If not, print error.
            if ((node.parent != cat_node) and (node != cat_node)):
                err_msg = 'Error('+g_module_name+'.CatLeafNodesToPtkns()):\n'+\
                          '      Invalid variable (category:leaf) pair\n'
                if dbg_loc != None:
                    cat_node_str = NodeToStr(cat_node)
                    leaf_node_str = NodeToStr(leaf_node)
                    err_msg += '    located near '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno)+'\n'+\
                               '    (\"'+leaf_node.name+'\" is not in the scope of \"'+cat_node_str+'/'+cat_name+':\")\n'+\
                               '    This will happen if you used the \"category\" command to manually\n'+\
                               '    create a category/counter which is not defined globally.\n'+\
                               '\n'+\
                               '    Note: Using the analogy of a unix style file system, \n'+\
                               '    the problem is that \"'+leaf_node_str+'\"\n'+\
                               '    is not a subdirectory of \"'+cat_node_str+'\".\n'+\
                               '\n'+\
                               '    Note: This often occurs when \".../\" is used. In that case, you may\n'+\
                               '    be able to avoid this error by referring to your variable explicitly\n'+\
                               '    by using chains of \"../\" tokens in the path instead of \".../\".\n'
                               #'       Make sure that your variable you are using is defined in \n'+\
                               #'       an environment (currently \"'+leaf_node_str+'\")\n'+\
                               #'       which lies WITHIN the environment where the category was defined.\n'+\
                               #'       (currently \"'+cat_node_str+'\").\n'
                    raise InputError(err_msg)
    else:
        err_msg = 'Warning: Strange variable path'
        if dbg_loc != None:
            err_msg += ' near '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno)
        err_msg += '\n' +\
            '    The category and leaf nodes for variable \"'+cat_name+':'+leaf_node.name+'\" are the same.\n'+\
            '    Check to see that this variable is behaving the way you intended.\n'+\
            '    (It\'s possible this could be an internal error in the program.)\n'
        sys.stderr.write(err_msg)

    # Merge the list of strings together into a single string:
    return cat_node_ptkns + leaf_node_ptkns



def CanonicalCatName(cat_name, cat_node, dbg_loc=None):
    # Determine the path of the cat node
    tkns = NodeToPtkns(cat_node)
    tkns.append(cat_name)
    #full_cat_name = tkns[0]
    #for i in range(1,len(tkns)):
    #    full_cat_name += '/'+tkns[i]
    #  better way:
    return '/'.join(tkns)


    
def CanonicalDescrStr(cat_name, cat_node, leaf_node, dbg_loc=None):
    tkns = CatLeafNodesToTkns(cat_name, cat_node, leaf_node, dbg_loc)
    descr_str = tkns[0]
    for i in range(1, len(tkns)):
        if (len(descr_str)>0) and (descr_str[-1] == ':'):
            descr_str += tkns[i]
        else:
            descr_str += '/'+tkns[i]
    return descr_str



def CollapsePath(path_tokens):
    """
    CollapsePath() takes a list of Strings argument representing a 
    directory-like path string 
    (for example '/SUB1A/Sub2A/../Sub2B/sub3b/../sub3c/entry'), 
    and replaces it with a version which should contain no '..' patterns. 
    (In the example above, it returns /SUB1A/Sub2B/sub3c/entry')
    """
    new_ptkns = []
    ndelete = 0
    i = len(path_tokens)-1
    while i >= 0:
        if path_tokens[i] == '..':
            ndelete += 1
        else:
            if (ndelete > 0) and (path_tokens[i] != ''):
                # Note: "path_tokens[i] != '')"  means "/a/b//c" <-> "/a/b/c"
                ndelete -= 1
            else:
                if len(path_tokens[i]) > 0: 
                    new_ptkns.append(path_tokens[i])
        i -= 1
    new_ptkns.reverse()

    if ndelete > 0:
        return ndelete  # <-- useful to let caller know an error ocurred

    return new_ptkns





def FindCatNode(category_name, current_node, srcloc):
    """ Search upwards (toward the ancester nodes), looking for a node
    containing a category matching category_name (first argument).
    Useful when the user specifies a category name, but neglects to 
    specify which node it was defined in.
    Note: there is no gaurantee that the category node returned by this function
          contains an entry in it's "categories" list corresponding to this
          category name.  You must check for this condition and handle it."""
    cat_node = None
    node = current_node
    while True:
        if category_name in node.categories:
            cat_node = node
            break
        elif node.parent != None:
            node = node.parent
        else:
            # node.parent is None,    ... we're done
            break

    if cat_node is None:
        assert(node.parent is None)
        #sys.stderr.write('Warning near ' +
        #                 ErrorLeader(srcloc.infile, 
        #                             srcloc.lineno)+'\n'+
        #                 '       no category named \"'+category_name+'\" found.\n'+
        #                 '       Creating a new global category: /'+
        #                 category_name+':\n')
        cat_node = node # the global node

    assert(cat_node != None)

    return cat_node


def RemoveNullTokens(in_ptkns):
    """This function just gets rid of useless empty tokens in the path ('', '.')
       (However if '' appears at the beginning of a path, we leave it alone.)

    """
    out_ptkns = []
    for i in range(0,len(in_ptkns)):
        if ((in_ptkns[i] != '.') and
            ((in_ptkns[i] != '') or (i==0))):
            out_ptkns.append(in_ptkns[i])
    # (I'm sure there are ways to write this in python 
    #  using fewer lines of code.  Sigh.)
    return out_ptkns
            


def DescrToCatLeafPtkns(descr_str, dbg_loc):
    """  
        Review: Variables in this program have three parts:
     1) A variable category name (designating the type of variable).
     2) A variable category path, which consists of a node which is an ancestor 
        of the variable leaf (1) in the tree
     3) A variable name ("leaf"), which refers to a node in the tree
        (either a static type tree or instance tree)

       DescrToCatLeafPtkns() takes a string describing a variable, 
    as it appears in a template     (ie, a write() command, once it has been 
    stripped of it's '$' or '@' prefix, and surrounding {} brackets)
    ...and divides it into strings which specify the location of that leaf in
    a static or instance tree, in addition to the name and location of the
    category node.  Descriptor examples for atoms in water:
              "AtomType:/Water/O",  There are only 2 --types-- of atoms in 
              "AtomType:/Water/H",  a water molecule. We identify them this way.
              "AtomID:O"     However each water molecule has 3 atoms, and we 
              "AtomID:H1"    can give each atom in each water molecule a unique 
              "AtomID:H2"    AtomID number.  "AtomID:H2" is the id number of the
                             second hydrogen atom in the current water molecule.

      ---- Output:  This function returns a 3-tuple: ----

    leaf_ptkns  The name of the variable's leaf node, as well as the list of 
               tokens denoting the path (named list of nodes) which lead to it.
    cat_name   The name of the variable category (no path information)
    cat_ptkns   A --suggestion-- for where to find the node containing the
               category mentioned in "cat_name".  Same format as leaf_ptkns.

    Examples:
    "AtomType:/Water/O"  cat_name='AtomType', cat_path=[], leaf_ptkns=['','Water','O']
    "AtomType:/Water/H"  cat_name='AtomType', cat_path=[], leaf_ptkns=['','Water','H']

    "AtomID:O"       cat_name='AtomID', cat_path=[], leaf_ptkns=['O']
    "AtomID:H1"      cat_name='AtomID', cat_path=[], leaf_ptkns=['H1']
    "AtomID:H2"      cat_name='AtomID', cat_path=[], leaf_ptkns=['H2']
    "mol:/"        cat_name='mol', cat_path=[], leaf_ptkns=['']
    "mol:"         cat_name='mol', cat_path=[], leaf_ptkns=[]
    "mol:../"      cat_name='mol', cat_path=[], leaf_ptkns=['..']
    "../mol"       cat_name='mol', cat_path=[], leaf_ptkns=['..']
    "$/peptide[3]/ResID:res[25]"  cat_name='ResID', cat_path=['', 'peptide[3]'], leaf_ptkns=['res[25]']

    """

    split_colon = descr_str.split(':')
    if len(split_colon) > 2:
        raise InputError('Error('+g_module_name+'.DescrToCatLeafPtkns())\n'
                         '       Error near '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno)+'\n\n'
                         '       Bad variable descriptor: \"'+descr_str+'\"\n'+
                         '       There can be at most one \':\' character in a variable descriptor.\n')
    # ---- Are we using colon syntax (example '$atom:H1')? 
    elif len(split_colon) == 2:
        # The category name = text after the last '/' (if present)and before ':'
        cat_ptkns  = split_colon[0].split('/')
        cat_name   = cat_ptkns[-1]
        # The text before that is the suggested (category) path 
        cat_ptkns  = cat_ptkns[:-1] 
        # if len(cat_ptkns) == 0:
        #    cat_ptkns.append('.')

        # The remaining text is the path leading to the leaf node.
        if split_colon[1] != '':
            leaf_ptkns = split_colon[1].split('/')
        else:
            leaf_ptkns = []

        if (cat_name == ''):
            raise InputError('Error('+g_module_name+'.DescrToCatLeafPtkns()):\n'
                             '       Error near '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno)+'\n\n'
                             '       Bad variable descriptor: \"'+descr_str+'\"\n')
    else:
        # ---- Are we using colon-less syntax (example: "$../mol") ?
        ptkns = split_colon[0].split('/')
        cat_name  = ptkns[-1] # last token (eg. "mol") is the cat_name
        leaf_ptkns = ptkns[:-1] # the rest is the leaf's path ("..")
        if len(leaf_ptkns) == 0:
            leaf_ptkns.append('.')
        #cat_ptkns  = ptkns[:-1] # the same goes for the cat path suggestion
        #if len(cat_ptkns) == 0:
        #    cat_ptkns.append('.')
        cat_ptkns  = []



    #  On 2012-8-22, I commented out this line:
    #return cat_name, RemoveNullTokens(cat_ptkns), RemoveNullTokens(leaf_ptkns)
    #  and replaced it with:

    return cat_name, RemoveNullTokens(cat_ptkns), leaf_ptkns








def DescrToCatLeafNodes(descr_str, 
                        context_node, 
                        dbg_loc, 
                        create_missing_nodes=False):
    """
    Variables in ttree correspond to nodes in a tree 
    (and also categories to which they belong). 
    DescrToCatLeafNodes() reads the name of a variable,
    (its descriptor) and determines where in the tree
    does this variable reside, and what is it's category?
    This function is the heart of ttree because it is
    the function used to interpret ttree variable syntax.
    (It is very messy right now. I will clean up the code later. AJ 2011-9-06)

       Arguments:
    descr_str    The complete name that the user gave 
                 to the variable. (Excluding '$' or '@')

    context_node  The class (node) in which the variable
                  was used.  descr_str is interpeted 
                  relative to this context.  (This argument
                  is similar to the current directory
                  in which a command was issued in unix.)

    dbg_loc      The location in the user's input file(s)
                 where this variable is referred to.

    create_missing_nodes
                 If we lookup a variable whose leaf node
                 does not exist yet, should we create it?
                 Setting this argument to "True" allows
                 us to augment the tree to add nodes 
                 corresponding to variables.



    -- Here is a greatly simplified version of DescrToCatLeafNodes(): --

    def DescrToCatLeafNodes(descr_str, context_node, dbg_loc):
       cat_name, cat_ptkns, leaf_ptkns = DescrToCatLeafPtkns(descr_str, dbg_loc)
       cat_node = PtknsToNode(cat_ptkns, context_node, dbg_loc)
       if len(cat_ptkns) > 0:
           leaf_node = PtknsToNode(leaf_ptkns, cat_node, dbg_loc)
       else:
           leaf_node = PtknsToNode(leaf_ptkns, context_node, dbg_loc)
       return cat_name, cat_node, leaf_node

    (This version works, but it does not handle "..." corectly,
    and it does not create missing nodes when needed.)


    -- Here is a (probably unnecessary) review of terminology: --

    Descriptor String:
    The first argument ("descr_str") is a descriptor string.
    A descriptor string typically contains ":" and "/" 
    characters to to divide the string into pieces in order
    to identify a category name, category node, and leaf node.
    Conceptually, the variable's NAME is the leaf node.
    The variable's TYPE is the category (node and name).

    Node:
    Nodes are used to represent both class objects and variable names
      1) class objects 
    Each type of class objects is represented by an StaticObj.
    Each instantiated object is represented by an InstanceObj.
      2) variable names (leaf nodes)
    However variable names are also represented using either
    StaticObjs (for @ static variables) or
    InstanceObjs (for $ instance variables)
    Again, all variables in ttree are members of a class object.
    In this case, the name of the node corresponds to the variable's
    name, and it's position in the tree refers to the class to which 
    it belongs.
    However "leaf nodes" do not uniquely identify the 
    actual variable itself. A single node can refer to two different 
    variables if they are in different categories.
    All 3 identifiers (leaf node, category node, category name)
    are needed to uniquely identify a ttree variable.  See below.

    Ptkn (Path Token) 
    Strings containing multiple '/' characters are typically used
    to identify the location of the category and leaf nodes in the
    tree (ie the path to the node).  The '/' characters are 
    delimiters which break up the string into small pieces, (which 
    are usually the names of classes).  
    These pieces are called "path tokens" or "ptkns"

    Leaf Node:    
    It exists as a node in a tree (instead of a simple string)
    because, just like member variables in a class in an
    object oriented programming language (or in a C struct)
    language, variables in ttree belong to the class in 
    which they are defined.  The node's location in the 
    tree represents which class it belongs to.
    If a variable's leaf node name 
    refers to a node which does no exist yet, then we create it
    (assuming the "create_missing_nodes" argument is "True").

    Category Node/Name:
    Categories are a peculiar feature of ttree.  Categories
    are groups of variables that share the same counter when
    numeric values are automatically given to each variable.
    So you can think of a category as a counter with a name.
    Variables in different categories have different counters,
    and are assigned numeric values independently. 
    Consequently two variables in different categories 
    may be assigned the same number.  But two variables 
    in the same category are always given unique values.
    Counters are typically global, but can have local scope.
    (ie, only defined within a Class, or an instantiated 
     class, and whatever other classes are nested or
     instantiated beneath it.)  
    Therefore to identify a counter/category you must specify
    both a name AND a node. The node identifies the class where 
    the scope is defined.  It is assumed that the Leaf Node 
    (see above) lies within this scope (ie. somewhere after
    it in the tree).
      Example: local counters are used to keep track of the 
    residues within in a protein chain. If we use a class to
    represent the protein, we can create a local residue-
    counter (category) within that protein.  Then when we 
    instantiate the protein multiple times, this counter 
    is reset for every new instance of of the protein.

    """

    cat_name, cat_ptkns, leaf_ptkns = DescrToCatLeafPtkns(descr_str, dbg_loc)



    # ---- ellipsis hack ----
    #
    # Search for class:
    # Most users expect ttree.py to behave like a 
    # standard programming language: If the class they are 
    # instantiating was not defined in this specific 
    # location, they expect ttree.py to search for
    # it outwards, first in the parent's environment, 
    # and then in the parent's parent's environment, 
    # and so on, until the object is found.
    # For example, most users expect this to work:
    # class Res{
    #   write("Atoms") { 
    #     $atom:CA @atom:CA 0.123 1.234 2.345
    #     $atom:CB @atom:CB 1.234 2.345 3.456
    #   }
    # }
    # class Protein{
    #   write_once("AnglesByType") {
    #     @angle:backbone @atom:Res/CA @atom:Res/CA @atom:Res/CA
    # }
    # Notice that in class Protein, we did not have to specify 
    # where "Res" was defined because it is defined in the parent
    # environment (ie. immediately outside Proteins's environment).
    #    The general way to do this in ttree.py, is to
    # use ellipsis syntax "@atom:.../Res/CA" symbol.  The 
    # ellipsis ".../" tells ttree.py to search upwards
    # for the object to the right of it ("Res")
    #    In order to make ttree.py behave the way 
    # most users are expecting, we artificially insert a 
    # ".../" before the class name here.  (Later on, the
    # code that processes the ".../" symbol will take
    # care of finding A.  We don't have to worry about
    # about doing that now.)
    #
    #   I think we only want to do this for variables with path information 
    # such as "@atom:Res/CA" (which means that leaf_ptkns = ['Res', 'CA']).
    # For simple variables like "@atom:CA", we don't automatically look upwards 
    # unless the user eplicitly requests it.
    #    (That's why we check to make sure that len(leaf_ptkns) > 1 below
    #     before we insert '...' into the leaf_ptkns.)
    # In other words, the two variables "@atom:CA" below are treated differently
    #
    # A {
    #   write("Atoms") {
    #     @atom:CA
    #   }
    #   class B {
    #     write("Atoms") { 
    #       @atom:CA
    #     }
    #   }
    # }
    #
    if ((descr_str.find(':') != -1) and
        #(not ((len(leaf_ptkns) == 1) and
        #      (leaf_ptkns[0] == context_node.name))) and
        #(len(leaf_ptkns) > 0) and 
        (len(leaf_ptkns) > 1) and
        (len(leaf_ptkns[0]) > 0) and 
        (leaf_ptkns[0][0] not in ('.','*','?'))):

        leaf_ptkns.insert(0, '...')
    # ---- Done with "ellipsis hack" -----



    #sys.stderr.write(' DescrToCatLeafNodes(): (cat_ptkns, cat_name, lptkns) = ('+
    #                 str(cat_ptkns)+', \"'+cat_name+'\", '+str(leaf_ptkns)+')\n')

    cat_node   = None
    cat_start_node  = context_node
    leaf_start_node = context_node

    if (len(cat_ptkns) > 0):
        if cat_ptkns[-1] == '...':        
            # The "..." in this position means trace the path from the
            # current node (context_node) up to cat_ptkns[:-1].
            cat_start_node = PtknsToNode(cat_ptkns[:-1], context_node, dbg_loc)
            # Later on, we will search upwards until we find an ancestor
            # node containing a category matching cat_name.  This will
            # be taken care of later.  (See "if cat_node is None:" below.)
        else:
            # In this case, the user supplied an explicit path
            # for the category node.  Find it now.
            cat_node   = PtknsToNode(cat_ptkns, context_node, dbg_loc)
            # Whenever the user supplies an explicit path, then
            # the cat node should be the starting location from
            # which the leaf path is interpreted.  This nearly 
            # insures that the leaf node will be an ancestor 
            # of the category node, which is what we want.
            leaf_start_node = cat_node

    if cat_node is None:
        # Otherwise, the user did not indicate where the category
        # node is defined, but only supplied the category name.
        # (This is the most common scenario.)
        # In this case, climb up the tree to the parent
        # until you find an ancestor with a category whose
        # name matches cat_name.  
        cat_node = FindCatNode(cat_name, cat_start_node, dbg_loc)

    if (cat_name not in cat_node.categories):
        if create_missing_nodes:
            # If this is the first time we encountered a variable in this 
            # category (ie if it's the first time we encountered a variable 
            # with this category's name and node), then we must create a
            # new entry in the cat_node.categories associative container
            # (using cat_name as the dictionary key).
            cat_node.categories[cat_name] = Category(cat_name)
        else:
            raise InputError('Error('+g_module_name+'.DescrToCatLeafNodes()):\n'
                             '       Error near '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno)+'\n'
                             '       Category named \"'+cat_name+'\" not found at\n'
                             '       position '+NodeToStr(cat_node)+'\n')


    # ---------- Now look up the leaf node -----------



    if (len(leaf_ptkns) > 0) and (leaf_ptkns[-1] == 'query()'):
        #   Special case: "query()"
        # Variables named "query()" are not really variables. 
        # (They are a way for users to query a category's counter.)
        # But we treat them as such internally. Consequently we
        # give them unique names to avoid clashes (just in case 
        # "query()" appears multiple times in the same context).
        leaf_ptkns[-1] = '__query__'+str(dbg_loc.order)



    # Lookup the path for the leaf:
    #
    # Often, the leaf that the path refers to does not 
    # exist yet.  For example, it is common for a template to
    # contain a reference to "$atom:CA". If the current context_node
    # is "/protein1/res22", this means that the leaf should be 
    # at "/protein1/res22/CA".  (However in this example, "CA" 
    # is not a class that has been defined yet.  It is the name
    # of a variable which which may not have even been mentioned 
    # before.  Think of "CA" as a variable placeholder. 
    #
    # So we follow the path tokens as far as we can:
    i_last_ptkn, last_node = FollowPath(leaf_ptkns, 
                                        leaf_start_node, 
                                        dbg_loc)

    # Did we find the node?
    if i_last_ptkn == len(leaf_ptkns):

        leaf_node = last_node

    else:

        # If we are here, then we did not find the node.
        # The unrecognized token is stored in
        #   leaf_ptkns[i_last_ptkn]

        if leaf_ptkns[i_last_ptkn] == '...':
            # ----------------------------------------------
            # ----    UGHH  I hate dealing with '...'   ----
            # ----(Messy code to follow in this section)----
            # ----------------------------------------------
            # The "..." means different things depending on
            # whether or not it is the last token in leaf_ptkns.

            if i_last_ptkn+1 < len(leaf_ptkns):
                # If "..." is NOT the last token in leaf_ptkns, we
                # should search for an ancestor of this node who has
                # a child whose name matches a the requested target
                # string (located in leaf_ptkns[i_last_ptkn+1])
                search_target = leaf_ptkns[i_last_ptkn+1]
                # If such an ancestor exists, then FollowPath() 
                # should have already found it for us. 
                #     This means it was not found.
                # So if there is only one more token in the
                # list of tokens, then create the needed node
                if (create_missing_nodes and
                    (i_last_ptkn+1 == len(leaf_ptkns)-1)):
                    # Create a new leaf node and link it:
                    new_leaf_name = leaf_ptkns[-1]
                    parent_node = last_node
                    # Is this parent_node an StaticObj? (..or inherit from StaticObj?)
                    if isinstance(parent_node, StaticObj):
                        parent_node.children.append(StaticObj(new_leaf_name, parent_node))
                    elif isinstance(parent_node, InstanceObj):
                        parent_node.children.append(InstanceObjBasic(new_leaf_name, parent_node))
                    else:
                        assert(False) # (only 2 types of nodes are possible)
                    # Now assign the pointer
                    leaf_node = parent_node.children[-1]
                else:
                    #In that case, we were unable to find the node referenced by "..." 
                    raise InputError('Error('+g_module_name+'.DescrToCatLeafNodes()):\n'
                                     '       Broken path.\n' # containing ellipsis (...)\n'
                                     '       class/variable \"'+search_target+'\" not found in this\n'
                                     '       context: \"'
                                     #+var_ref.prefix + var_ref.descr_str + var_ref.suffix+'\"\n'
                                     +descr_str+'\"\n'
                                     '       located near '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno))


            else:    # if i_last_ptkn+1 < len(leaf_ptkns):


                # If "..." IS the last token, then it means:
                # we want to search for the CATEGORY NAME,
                # This is very different.
                # It means we need to:
                #   search backwards up the ancestor tree until
                #   we find an ancestor variable (of last_node)
                #   which has the right category, (ie until you 
                #   find an ancestor node with a variable (VarRef) 
                #   pointing to it with belonging to the correct
                #   category node and name (determined above).)
                #  If not found, then use the current context_node.

                assert(cat_name in cat_node.categories)
                var_bindings = cat_node.categories[cat_name].bindings

                node = last_node
                while (node != None):
                    # Recall that cat_node.categories[cat_name]
                    # is a dictionary whose keys are leaf nodes
                    # corresponding to the variables in this category.
                    if node in var_bindings:
                        # then we found it, and we're done
                        break
                    else:
                        node = node.parent

                if node != None:
                    leaf_node = node

                else:
                     # If not found, have it point to the 
                     # current (context) node.
                     leaf_node = context_node

            # -----------------------------------------------
            # -- Finished dealing with '...' in leaf_ptkns --
            # -----------------------------------------------


        elif (create_missing_nodes and
              ((i_last_ptkn == len(leaf_ptkns)-1) or
               HasWildCard('/'.join(leaf_ptkns)))):

        #elif (create_missing_nodes and
        #      (i_last_ptkn == len(leaf_ptkns)-1)):

            # Again, another reason the leaf-node was not found is
            # that it refers to a leaf node which has not yet been
            # created.  If the path was valid until up to the last
            # token, then we sould create a new node with this name.
            #  -- This is a common scenario.                 --
            #  -- This is how all new variables are created. --
            # Anyway, we handle that here:

            # Create a new leaf node and link it:
            new_leaf_name = leaf_ptkns[-1]
            new_leaf_name = '/'.join(leaf_ptkns[i_last_ptkn:])
            parent_node = last_node
            # Is this parent_node an StaticObj? (..or does it inherit from StaticObj?)
            if isinstance(parent_node, StaticObj):
                parent_node.children.append(StaticObj(new_leaf_name, parent_node))
            elif isinstance(parent_node, InstanceObj):
                parent_node.children.append(InstanceObjBasic(new_leaf_name, parent_node))
            else:
                assert(False) # (only 2 types of nodes are possible)
            # Now assign the pointer
            leaf_node = parent_node.children[-1]


        else:

            # Otherwise, the user made a mistake in the path.
            # Figure out which kind of mistake and print an error.

            if (last_node.parent is None) and (leaf_ptkns[i_last_ptkn] == '..'):
                #In that case, we tried to back out beyond the root of the tree.
                raise InputError('Error('+g_module_name+'.DescrToCatLeafNodes()):\n'
                                 '       Broken path in variable:\n'
                                 #'       \"'+var_ref.prefix + var_ref.descr_str + var_ref.suffix+'\"\n'
                                 '       \"'+ descr_str + '\"\n'
                                 '    located near '+
                                 ErrorLeader(dbg_loc.infile, 
                                             dbg_loc.lineno)+'\n'
                                 '       There are too many \"..\" tokens in the path string.')


            else:

                #Then the reason is: The string in leaf_ptkns[i_last_ptkn]
                #was supposed to be a child of last_node but a child
                #of that name was not found.
                raise InputError('Error('+g_module_name+'.DescrToCatLeafNodes()):\n'
                                 '       Broken path / Undefined variable:\n'
                                 #'       \"'+var_ref.prefix + var_ref.descr_str + var_ref.suffix+'\"\n'
                                 '       \"'+ descr_str + '\"\n'
                                 '    located near '+
                                 ErrorLeader(dbg_loc.infile, 
                                             dbg_loc.lineno)+'\n'
                                 '       Undefined: \"'+PtknsToStr(leaf_ptkns)+'\"\n'
                                 '       (Specifically \"'+leaf_ptkns[i_last_ptkn]+
                                 '\" is not a subordinate of \"'+MaxLenStr(last_node.name,'/')+'\")')
                                 #'\n    This could be a typo or spelling error.')

    return cat_name, cat_node, leaf_node









def DescrToVarBinding(descr_str, context_node, dbg_loc):
    """ DescrToVarBinding() is identical to LookupVar(), but it has a name 
    that is harder to remember.  See comment for LookupVar() below.

    """
    cat_name, cat_node, leaf_node = DescrToCatLeafNodes(descr_str, 
                                                        context_node, 
                                                        dbg_loc)

    if cat_name in cat_node.categories:
        category = cat_node.categories[cat_name]
        var_bindings = category.bindings
        if leaf_node in var_bindings:
            var_binding = var_bindings[leaf_node]
        else:
            raise InputError('Error('+g_module_name+'.DescrToVarBinding()):\n'
                             '       Error near '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno)+'\n'
                             '       Bad variable reference: \"'+descr_str+'\".  There is\n'
                             '       There no category named \"'+cat_name+'\" defined for "'+NodeToStr(cat_node)+'\".\n')
    else:
        raise InputError('Error('+g_module_name+'.DescrToVarBinding()):\n'
                         '       Error near '+ErrorLeader(dbg_loc.infile, dbg_loc.lineno)+'\n'
                         '       Bad variable reference: \"'+descr_str+'\".  There is\n'
                         '       no category named \"'+cat_name+'\" defined for "'+NodeToStr(cat_node)+'\".\n')

    return var_binding





# Wrappers:

def LookupVar(descr_str, context_node, dbg_loc):
    """ LookupVar() looks up a string (a variable descriptor, which is the
    variable's name, excluding the '$', '@' prefixes and any '{}' brackets.)
    This function returns the variable's "VarBinding" (the variable-name:value 
    pair). This is useful for querying or changing the value of a variable.
        Because nearly all variables are local, you must specify the starting 
    node (ie. the node corresponding to the class in which this class 
    or variable was referred to).  This is typically the global node.

    """
    return DescrToVarBinding(descr_str, context_node, dbg_loc)




def LookupNode(obj_name, starting_node, dbg_loc):
    """ LookupNode() parses through a string like 
          '../ClassA/NestedClassB'
        and returns the corresponding node.
        Nodes are data types used for representing a class or class instance.
        They are also used for storing variables.
          'ClassA/NestedClassB/VariableC'
        Because nearly all variables are local, you must specify the starting 
        node (ie. the node corresponding to the class in which this class 
        or variable was referred to).  This is typically the global node.

        """
    return StrToNode(obj_name, starting_node, dbg_loc)





class SimpleCounter(object):
    __slots__=["n","nincr"]

    def __init__(self, n0 = 1, nincr = 1):
        self.n = n0 - nincr
        self.nincr = nincr

    def query(self):
        return self.n

    def incr(self):
        self.n += self.nincr

    def __copy__(self): #makes a (deep) copy of the counter in its current state
        return SimpleCounter(self.n + self.nincr, self.nincr)
    



class Category(object):
    """
    Category contains a list of all of the variables which belong
    to the same category, as well as some supporting information.

       Attributes:

    name         The name of the category (a string)

    bindings     An OrderedDict() containing leaf_node:VarBinding 
                 (key:value) pairs. Variables are looked up by their leaf node.
                 The actual variable name (which simply refers to the leaf node)
                 and values are both stored in the VarBinding data structure.

    counter      A counter object like "SimpleCounter". Each time counter.incr()
                 is invoked it should return a unique string (typically this is
                 simply a string representing an integer which is incremented).

    """

    __slots__=["name","bindings","counter","manual_assignments","reserved_values"]

    def __init__(self, 
                 name = '', 
                 bindings = None, 
                 counter = None,
                 manual_assignments = None,
                 reserved_values = None):

        self.name = name

        if bindings is None:
            self.bindings = OrderedDict()
        else:
            self.bindings = bindings

        if counter is None:
            self.counter = SimpleCounter(1,1)
        else:
            self.counter = counter

        if manual_assignments is None:
            self.manual_assignments = OrderedDict()
        else:
            self.manual_assignments = manual_assignments

        if reserved_values is None:
            self.reserved_values = OrderedDict()
        else:
            self.reserved_values = reserved_values



class StaticObj(object):
    """  StaticObjs and InstanceObjs:

       The state of the system is stored in two different trees:
    1) The static tree:
       StaticObj trees are similar "class" definitions in an OOP language.
       These trees contains class definitions, and their nested classes,
       and instructions for how to create new instances (copies) of this class.
       Nodes in this tree are stored using StaticObjs:
    2) The instance tree: 
       This tree contains classes that have been instantiated, and any sub-
       classes (members or attributes) that are instantiated as a result.
       This tree is automatically generated by instantiating the root
       StaticObj.  Nodes in this tree are stored using InstanceObjs.

      StaticObjs and InstanceObjs both contain 
         "commands" (commands which usually involve instructions 
                     for writing templates)
         "categories" (local counters used to assign variables. See below.)
         "children" (branches in the tree. -Not- OOP class children. See below.)
      StaticObjs also contain
         "instance_commands"
         "instance_categories"
         These three members contain information to create a new instance/copy
         of this class (how to construct an InstanceObj from an StaticObj).

    StaticObj contains the member function Parse() which builds the global static 
    tree by parsing the contents of a text file supplied by the user.

    The external function BuildInstanceTree(), creates the global instance tree
    from the global static tree (a tree of StaticObjs).

    -----  CLASS MEMBERS OF StaticObj: ----
     0) Name:
        Every class (object type) has a name.  It is stored in self.name.
     To make it easier to distinguish the names of classes from the names of
     individual instances of that class, I recommend using a capital letter
     for the name of a class type (and lower-case letters for instances).

     1) Commands 
        Commands are usually instructions for writing templates.
      Templates are blocks of ordinary text which contain variables.  
      (Variables in this program consist of variable names, categories, 
      and (eventually) bound values (usually generated automatically),
      which will be substituted into the template to generate a text file.)
      A class can contain multiple templates, each one having a unique name
      which also happens to be the name of the file that will be created when
      the template is written.
      Variants:

          self.commands:
      Some templates are written immediate after the class is defined
      (stored in "self.commands").
      Example: The "write_once()" command.

          self.instance_commands:
      Other templates are written when an instance/copy of the class is created
      (stored in "self.instance_commands".
      Example: The "write()" command.

    2) Children
          self.children:
       Class definitions can be defined from within the definition of other
       ("parent") classes.  These nested classes are referred to as "children".
       These sub-classes are not "children" in the OOP sense of the word at
       all (they do not have any of the the traits of their "parents").
       However in the source code I refer to them as "children" because here 
       they are implemented as "child" nodes (branches) in the tree-like
       data-structure used to store class definitions (the static tree).

    3) Categories
        This is a new concept and is difficult to explain.
        Recall that each class contains a list of templates containing raw text,
        interspersed with variables (whose values will determined later).
        In most cases, variables are assigned to integer values which are 
        automatically generated by incrementing a counter. Simply put, 
        "categories" are collections of variables which share the same counter.
        Within a category, the goal is to assign a unique integer (or other 
        symbol) to each distinct variable in this category.
        To avoid name-clashes, variable names have local "scope".
        This scope is the "leaf_token"

        Categories can be specific to a particular class (node), and any of the
        classes (nodes) which are nested within it, but by default are global.
        (This means they "belong" to the global (outermost) node by default.)
        All the various categories which are defined within a particular
        StaticObj are stored in self.categories.
        Static variables (ie. variables with a '@' prefix) are stored this way.

           "self.categories"
        If a class contains a new category, it means that if any nested 
        classes defined within that class definition contain (static, '@') 
        variables whose categories match the category name, their values will 
        be determined by looking up the couter associated with this category
        stored locally (here) in self.categories.  All variables belonging
        to this category are stored in "self.categories[category_name]".
           "self.instance_categories"
        Recall that creating a new copy (instance) of a class automatically 
        creates an InstanceObj in the instance-tree.  InstanceObj's have a 
        ".categories" attribute of their own, the contents of which are
        copied from this StaticObj's "self.instance_categories" attribute.
        Instantiating a new class also spawns the instantiation of any 
        sub-classes.
        If any of these "instance children" contain variables whose category
        names match a category stored in the parent InstanceObj's .categories
        dictionary, then their values will be determined by that InstanceObj's
        counter for that category name.
    4)  Parent: 
        A link to the parent StaticObj is stored in self.parent.

    """

    # sigh

    __slots__=["name",
               "parent",
               "children",
               "categories",
               "commands",
               "srcloc_begin",
               "srcloc_end",
               "deleted",
               "class_parents",
               "namespaces",
               "instname_refs",
               "instance_categories",
               "instance_commands_push",
               "instance_commands",
               "instance_commands_pop"]

    def __init__(self, 
                 name='',
                 parent=None):
        """
        The members/attributes of StaticObj are defined in the comment
        for StaticObj above.  """

        # The following members are shared by both InstanceObj and StaticObj:

        self.name = name

        self.parent = parent    #For traversing the global static template tree

        self.children = []           # Nested class definitions. 

        self.categories=OrderedDict()  #<- new variable categories that are only defined
                                 # in the context of this molecule's type definition
        self.commands=[]         # Commands to carry out (only once)

        ##vb##self.var_bindings=[]     # List of variables assigned to this object.

        self.srcloc_begin = None     # Keep track of location in user files
        self.srcloc_end   = None     # (useful for error message reporting)

        self.deleted      = False    # Users can delete static objects?
                                     # (why not?)

        # The following members are not shared with InstanceObj:

        self.class_parents = []  # classes we inherit traits from (this is
                                 # similar to the parent/child relationship
                                 # in an object-oriented-programming language)

        self.namespaces = []  # A list of classes we also look in when searching
                              # for other static nodes or variables. (similar to
                              # class_parents, but only used for searches.)

        self.instname_refs = {}  # <-- used for debugging to insure that
                                 #     two instances do not have the same name

        self.instance_categories=OrderedDict()#<-new variable categories that are defined
                                   #within the scope of this molecule's instance
        self.instance_commands_push=[]  #1)process these commands first by adding
                                       #  these commands to InstanceObj.commands
                                       #  (before you deal with class_parents)
        self.instance_commands=[]      #2) then add this to InstanceObj.commands
        self.instance_commands_pop=[] #3) finally add these commands



    def DeleteSelf(self):
        self.deleted = True



    def IsDeleted(self):
        return self.deleted


    ##vb##def AddVarBinding(self, var_binding):
    ##vb##    if self.var_bindings is None:
    ##vb##        self.var_bindings = [var_binding]
    ##vb##    else:
    ##vb##        self.var_bindings.append(var_binding)


    def Parse(self, lex):

        """ Parse() builds a static tree of StaticObjs by parsing text file.
        -The "lex" argument is a file or input stream which has been converted 
         to a "TemplateLexer" object (similar to the python's built-in shlex lexer).
        """

        #sys.stdout.write(' -- Parse() invoked --\n')

        # Keep track of the location in the users' input files where this
        # class object is first defined.  (Keep in mind that the user might 
        # augment their original class definition, adding new content to an
        # existing class.  In that case self.srcloc_begin will have already 
        # been assigned.  We don't want to overwrite it in that case.)
        if self.srcloc_begin is None:  # <-- not defined yet?
            self.srcloc_begin = lex.GetSrcLoc()

        while True:

            cmd_token = lex.get_token()

            #print('Parse(): token = \"'+cmd_token+'\", '+lex.error_leader())

            if cmd_token == lex.eof:
                #print('Parse(): EOF encountered\n')
                break


            if ((cmd_token == 'write') or 
                (cmd_token == 'write_once') or
                (cmd_token == 'create_var') or
                (cmd_token == 'create_vars')):
                open_paren  = lex.get_token()

                #print('Parse():     open_paren=\"'+open_paren+'\"')
                if open_paren=='{': 
                    # ..then the user neglected to specify the "dest" file-name
                    # argument.  In that case, supply the default, ''.
                    # (which is shorthand for the standard out in this case)
                    open_curly     = open_paren[0]
                    open_paren     = ''
                    close_paren    = ''
                    tmpl_filename = ''
                    srcloc  = lex.GetSrcLoc()
                else:
                    tmpl_filename = lex.get_token()
                    if tmpl_filename == ')':
                        tmpl_filename = ''
                        close_paren = ')'
                    else:
                        close_paren    = lex.get_token()
                    open_curly     = lex.get_token()
                    srcloc        = lex.GetSrcLoc()

                if ((cmd_token == 'create_var') or
                    (cmd_token == 'create_vars')):
                    tmpl_filename = None
                    # This means: define the template without attaching 
                    # a file name to it. (IE., don't write the contents
                    # of what's enclosed in the curly brackets { } to a file.)
                
                if ((open_curly != '{') or 
                    ((open_paren == '')  and (close_paren != '')) or
                    ((open_paren == '(') and (close_paren != ')'))):
                    raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                     '       Error in '+lex.error_leader()+'\n\n'
                                     'Syntax error at the beginning of the \"'+cmd_token+'\" command.')
                if tmpl_filename != None:
                    tmpl_filename = RemoveOuterQuotes(tmpl_filename, lex.quotes)
                # ( The previous line is similar to:
                #   tmpl_filename = tmpl_filename.strip(lex.quotes) )

                tmpl_contents = lex.ReadTemplate()
                StaticObj.CleanupReadTemplate(tmpl_contents, lex)

                #sys.stdout.write('  Parse() after ReadTemplate, tokens:\n\n')
                #print(tmpl_contents)
                #sys.stdout.write('\n----------------\n')

                        
                if cmd_token == 'write_once':
                    # Check for a particular bug:
                    #    Ordinary instance variables (preceded by a '$')
                    #    should never appear in a write_once() statement.
                    for entry in tmpl_contents:
                        if (isinstance(entry, VarRef) and 
                            (entry.prefix[0]=='$')):
                            raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                             '       Error near '+ErrorLeader(entry.srcloc.infile,
                                                                       entry.srcloc.lineno)+'\n'
                                             '       Illegal variable: \"'+entry.prefix+entry.descr_str+entry.suffix+'\"\n'
                                             '       All variables in a \"write_once()\" statement must be statically\n'
                                             '       defined, and hence they must begin with a \'@\' prefix character.\n'
                                             '       (not a \'$\' character).\n'
                                             '       Suggestion: Use the \"write()\" command instead.\n')



                if cmd_token == 'write':
                    commands = self.instance_commands
                elif cmd_token == 'write_once':
                    commands = self.commands
                elif ((cmd_token == 'create_var') or 
                      (cmd_token == 'create_vars')):
                    commands = self.instance_commands
                else:
                    assert(False) 

                command = WriteFileCommand(tmpl_filename, 
                                           tmpl_contents,
                                           srcloc)
                commands.append(command)

            # end of "if (cmd_token == 'write') or (cmd_token == 'write_once'):"

            elif cmd_token == 'delete':

                instobj_descr_str = lex.get_token()
                instobj_srcloc = lex.GetSrcLoc()
                delete_command = DeleteCommand(instobj_srcloc)
                mod_command = ModCommand(delete_command,
                                         instobj_descr_str)
                self.instance_commands.append(mod_command)

            elif cmd_token == 'using':

                namespacecom_str = lex.get_token()
                if namespacecom_str != 'namespace':
                    raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                     '       Error near '+lex.error_leader()+'\n'
                                     '       The \"'+cmd_token+'\" command must be followed by the \"namespace\" keyword.')
                namespace_str = lex.get_token()

                stnode = StrToNode(namespace_str,
                                   self,
                                   lex.GetSrcLoc())

                self.namespaces.append(stnode)

            elif cmd_token == 'category':

                cat_name = lex.get_token()

                cat_count_start = 1
                cat_count_incr  = 1

                open_paren = lex.get_token()
                if (open_paren == '('):
                    token = lex.get_token()
                    if token == ',':
                       token = lex.get_token()
                    if token != ')':
                        # Interpret token as an integer, float, or string
                        try:
                            cat_count_start = int(token)
                        except ValueError:
                            try:
                                cat_count_start = float(token)
                            except ValueError:
                                cat_count_start = RemoveOuterQuotes(token, '\'\"')
                        token = lex.get_token()
                        if token == ',':
                            token = lex.get_token()
                        if token != ')':
                            # Interpret token as an integer,float,or string
                            try:
                                cat_count_incr = int(token)
                            except ValueError:
                                try:
                                    cat_count_incr = float(token)
                                except ValueError:
                                    cat_count_incr = RemoveOuterQuotes(token, '\'\"')
                            token = lex.get_token()
                            if token != ')':
                                raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                                 '       Error near '+lex.error_leader()+'\n'
                                                 '       \"'+cmd_token+' '+cat_name+'...\" has too many arguments,\n'
                                                 '       or lacks a close-paren \')\'.\n')
                else:
                    lex.push_token(open_paren)

                if (isinstance(cat_count_start, basestring) or 
                    isinstance(cat_count_incr, basestring)):
                    raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                     '       Error near '+lex.error_leader()+'\n'
                                     '       \"'+cmd_token+' '+cat_name+'('+
                                     str(cat_count_start)+','+
                                     str(cat_count_incr)+')\"\n'
                                     '       Only numeric counters are currently supported.\n')

                # check for really stupid and unlikely errors:
                if type(cat_count_start) is not type(cat_count_incr):
                    if ((isinstance(cat_count_start, int) or
                         isinstance(cat_count_start, float))
                        and
                        (isinstance(cat_count_incr, int) or
                         isinstance(cat_count_incr, float))):
                        cat_count_start = float(cat_count_start)
                        cat_count_incr  = float(cat_count_incr)
                    else:
                        raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                         '      Error near '+lex.error_leader()+'\n'
                                         '      Problem with \"'+cmd_token+'\" command.\n')

                prefix = cat_name[0]
                cat_name = cat_name[1:]
                # Add this category to the list.
                if prefix == '@':
                    self.categories[cat_name] = Category(cat_name)
                    self.categories[cat_name].counter=SimpleCounter(cat_count_start,
                                                                    cat_count_incr)
                elif prefix == '$':
                    self.instance_categories[cat_name] = Category(cat_name)
                    self.instance_categories[cat_name].counter=SimpleCounter(cat_count_start,
                                                                             cat_count_incr)
                else:
                    raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                     '       Error near '+lex.error_leader()+'\n'
                                     '       category name = \"'+cat_name+'\" lacks a \'$\' or \'&\' prefix.\n'
                                     '       This one-character prefix indicates whether the variables in this\n'
                                     '       new category will be static or dynamics variables\n')



            elif (cmd_token == '}') or (cmd_token == ''):
                # a '}' character means we have reached the end of our scope.
                # Stop parsing and let the caller deal with the remaining text.
                # (And a '' means we reached the end of the file... I think.)
                break


            #elif (cmd_token == 'include'):
                # "include filename" loads a file (adds it to the file stack)
                # The "TtreeShlex" class (from which "lex" inherits) handles 
                # "include" statements (ie. "source" statements) automatically.

            else:
                # Otherwise, 'cmd_token' is not a command at all.
                # Instead it's the name of an object which needs to be
                # defined or instantiated.
                # First, let's figure out which.

                # (small detail: The "class" keyword is optional
                #                and can be skipped.)
                if cmd_token == 'class':
                    object_name = lex.get_token()
                else:
                    object_name = cmd_token

                next_symbol = lex.get_token()
                #print('Parse(): next_token=\"'+next_symbol+'\"')

                class_parents = []

                if next_symbol == 'inherits':

                    # Then read in the list of classes which are parents of
                    # of this class.  (Multiple inheritance is allowed.)
                    # (We don't yet check to insure that these are valid class 
                    #  names.  We'll do this later in LookupStaticRefs().)

                    syntax_err_inherits = False

                    while True:
                        next_symbol = lex.get_token()
                        if ((next_symbol == '{') or
                            (next_symbol == lex.eof)):
                            break
                        elif (next_symbol == '='):
                            syntax_err_inherits = True
                            break
                        else:
                            class_parents.append(StrToNode(next_symbol,
                                                           self,
                                                           lex.GetSrcLoc()))
                    if len(class_parents) == 0:
                        syntax_err_inherits = True
 
                    if syntax_err_inherits:
                        raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                         '      Error near '+lex.error_leader()+'\n'
                                         '      \"inherits\" should be followed by one or more class names.\n')


                if next_symbol == '{':
                    child_name = object_name

                    # Check to see if this class has already been defined. 
                    # (IE. check if it present in the list of children.) 
                    # If the name (child_name) matches another class (child), 
                    # then the contents of the new class will be appended to 
                    # the old.  This way, class definitions can be augmented 
                    # later.  (This is the way "namespaces" work in C++.)
                    child = None
                    for child in self.children:
                        if child.name == child_name:
                            break
                    # If found, we refer to it as "child". 
                    # If not, then we create a new StaticObj named "child". 
                    if (child is None) or (child.name != child_name):
                        child = StaticObj(child_name, self)
                        self.children.append(child)
                    assert(child in self.children)
                    # Either way we invoke child.Parse(), to 
                    # add contents (class commands) to child. 
                    child.Parse(lex)
                    child.class_parents += class_parents



                elif next_symbol == '=':
                    next_symbol = lex.get_token()
                    if next_symbol == 'new':
                        base_name = object_name
                        base_srcloc = lex.GetSrcLoc()
                        if base_name.find('/') != -1:
                            raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                             '       Error near '+ErrorLeader(base_srcloc.infile,
                                                                        base_srcloc.lineno)+'\n'
                                             '          (You can not instantiate some other class\'s members.)\n'
                                             '       Invalid instance name: \"'+base_name+'\"\n')
                            
                        elif base_name in self.instname_refs:
                            ref_srcloc = self.instname_refs[base_name]
                            raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                             '       Duplicate class/array \"'+base_name+'\"\n'
                                             '       This occurs near:\n'
                                             '         '+ErrorLeader(ref_srcloc.infile,
                                                                     ref_srcloc.lineno)+'\n'
                                             '       and also near:\n'
                                             '         '+ErrorLeader(base_srcloc.infile,
                                                                     base_srcloc.lineno)+'\n')
                        else:
                            self.instname_refs[base_name] = base_srcloc

                        # If the statobj_str token contains a ']' character
                        # then this means the user wants us to make multiple 
                        # copies of this template.  The number of copies 
                        # to instantiate is enclosed in the [] characters
                        # (Example wat = new Water[3000] creates
                        #  3000 instantiations of the Water template
                        #  named wat[1], wat[2], wat[3], ... wat[3000]).

                        # Note: Here '[' and ']' have a special meaning.
                        # So lex.get_token() should not treat them as 
                        # ordinary word characters.  To prevent this:
                        orig_wordterminators = lex.wordterminators
                        lex.wordterminators += '[],'

                        class_name_str = lex.get_token()
                        if ((class_name_str == lex.eof) or
                            (class_name_str == '}')):
                                raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                                 '       Error near '+lex.error_leader()+'\n'
                                                 'Class ends prematurely. (Incomplete \"new\" statement.)')

                        assert(len(class_name_str) > 0)

                        if (class_name_str[0] == '['):
                            raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                             '      Error near '+lex.error_leader()+'\n'
                                             '      new'+class_name_str+'\n'
                                             'Bracketed number should be preceeded by a class name.')
                        class_names = []
                        weights = []
                        if class_name_str == 'random':
                            class_names, weights = self._ParseRandom(lex)
                            tmp_token = lex.get_token()
                            if len(tmp_token)>0:
                                if tmp_token[0]=='.':
                                    raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                                     '      Error near '+lex.error_leader()+'\n'
                                                     '      \"'+tmp_token+'\" should not follow random()\n'
                                                     '\n'
                                                     '      Coordinate transformations and other commands (such as \"'+tmp_token+'\")\n'
                                                     '      should appear after each class name inside the random() statement,\n'
                                                     '      not after it.  For example, do not use:\n'
                                                     '      \"lipids=new random([DPPC,DLPC],[0.5,0.5]).move(0,0,23.6)\"\n'
                                                     '         Use this instead:\n'
                                                     '      \"lipids=new random([DPPC.move(0,0,23.6),DLPC.move(0,0,23.6)],[0.5,0.5])\"\n')
                                lex.push_token(tmp_token)
                        else:
                            class_name, class_suffix, class_suffix_srcloc = \
                                self._ProcessClassName(class_name_str, lex)

                        array_size = []
                        array_suffixes = []
                        array_srclocs = []

                        # A general "new" statement could look like this:
                        # "m = new Mol.scale(3) [2].trans(0,4.5,0).rotate(30,0,0,1)
                        #                       [3].trans(0,0,4.5)"
                        # So far we have processed "m = new Mol.scale(3)".
                        # Now, we need to deal with:
                        #  "[2].trans(0,4.5,0).rotate(30,0,0,1) [3].trans(0,0,4.5)"
                        while True:
                            new_token = lex.get_token()
                            #if ((new_token == '') or (new_token == lex.eof)):
                            #    break
                            if new_token == '[':
                                number_str = lex.get_token()
                                close_bracket = lex.get_token()
                                if ((not str.isdigit(number_str)) or
                                    (close_bracket != ']')):
                                    raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                                     '     Error in new statement near '+lex.error_leader()+'\n'
                                                     '     A \'[\' character should be followed by a number and a \']\' character.')
                                array_size.append(int(number_str))
                                suffix = lex.get_token()

                                if ((suffix == '') or (suffix == lex.eof)):
                                    array_suffixes.append('')
                                    array_srclocs.append(base_srcloc)
                                    break
                                if suffix[0] == '.':
                                    lex.push_token(suffix[1:])
                                    suffix_func = lex.GetParenExpr()
                                    suffix = '.' + suffix_func
                                    array_suffixes.append(suffix)
                                    array_srclocs.append(lex.GetSrcLoc())
                                else:
                                    array_suffixes.append('')
                                    array_srclocs.append(base_srcloc)
                                    lex.push_token(suffix)
                                    if suffix != '[':
                                        break
                            else:
                                lex.push_token(new_token)
                                break
                        srcloc_final = lex.GetSrcLoc()

                        lex.wordterminators = orig_wordterminators

                        assert(len(array_size) == len(array_suffixes))

                        # If the user wants us to instantiate a
                        # multidimensional array of class instances
                        # then we must loop through this multidimensional
                        # array and create a new instance for each entry.
                        # For example fill a 3 dimensional volume 
                        # with 1000 water molecules
                        # Example 1:
                        #    solvent = new Water [10][10][10]
                        #    (The coordinates must be read separately.)
                        #    In this example array_size = [10,10,10]
                        #                      array_suffixes = ['','','']
                        # Example 2:
                        #    solvent = new Water.transcm(0,0,0)
                        #                        [10].trans(0,0,4)
                        #                        [10].trans(0,4,0).rot(45,0,0,1)
                        #                        [10].trans(4,0,0)
                        #    (This command generates a 10x10x10 lattice
                        #     simple cubic lattice of regularly spaced
                        #     water molecules pointing the same direction.)
                        #    In this example array_size = [10,10,10]
                        #      and
                        #    class_suffix = 'transcm(0,0,0)'
                        #      and
                        #    array_suffixes = ['trans(0,0,4)', 
                        #                      'trans(0,4,0).rot(45,0,0,1)', 
                        #                      'trans(4,0,0)']
                        # Note that tree ignores the "trans()" 
                        #     commands, it stores them so that inherited 
                        #     classes can attempt to process them.


                        D = len(array_size)
                        if D > 0:
                            digits = [0 for d in range(0, D)]
                            table_filled = False
                            pushed_commands = []
                            while (not table_filled):
                                instance_name = base_name
                                for d in range(0, D):
                                    i = digits[d]
                                    instance_name += '['+str(i)+']'

                                # Does the user want us to select 
                                # a class at random? 
                                if len(class_names) > 0:
                                    class_name_str = RandomSelect(class_names,
                                                                  weights)
                                    class_name, class_suffix, class_suffix_srcloc= \
                                        self._ProcessClassName(class_name_str, lex)

                                if class_suffix != '':
                                    class_suffix_command = \
                                        PushRightCommand(class_suffix.lstrip('.'),
                                                         class_suffix_srcloc)
                                    self.instance_commands.append(class_suffix_command)
                                command = \
                                 InstantiateCommand(instance_name,
                                                    ClassReference(class_name,
                                                                   base_srcloc),
                                                    base_srcloc)
                                self.instance_commands.append(command)

                                if class_suffix != '':
                                    command = \
                                        PopRightCommand(class_suffix_command, 
                                                        srcloc_final)
                                    self.instance_commands.append(command)


                                # Now go to the next entry in the table.
                                # The indices of this table are similar to
                                # a D-digit integer.  We increment this d-digit number now.
                                d_carry = D-1
                                while True:
                                    digits[d_carry] += 1
                                    if digits[d_carry] >= array_size[d_carry]:
                                        digits[d_carry] = 0
                                        if array_suffixes[d_carry] != '':
                                            for i in range(0, array_size[d_carry]-1):
                                                partner = pushed_commands.pop()
                                                command = PopRightCommand(partner,
                                                                          srcloc_final)
                                                self.instance_commands.append(command)
                                        d_carry -= 1
                                    else:
                                        if array_suffixes[d_carry] != '':
                                            command = PushRightCommand(array_suffixes[d_carry].lstrip('.'),
                                                                       array_srclocs[d_carry])
                                            pushed_commands.append(command)
                                            self.instance_commands.append(command)
                                        break
                                    if d_carry < 0:
                                        table_filled = True
                                        break


                            pass


                        else:
                            if len(class_names) > 0:
                                class_name_str = RandomSelect(class_names,
                                                              weights)
                                class_name, class_suffix, class_suffix_srcloc= \
                                    self._ProcessClassName(class_name_str, lex)
                            if class_suffix != '':
                                class_suffix_command = \
                                    PushRightCommand(class_suffix.lstrip('.'),
                                                     class_suffix_srcloc)
                                self.instance_commands.append(class_suffix_command)
                            command = \
                              InstantiateCommand(base_name,
                                                 ClassReference(class_name,
                                                                base_srcloc),
                                                 base_srcloc)
                            self.instance_commands.append(command)

                            if class_suffix != '':
                                command = \
                                    PopRightCommand(class_suffix_command, 
                                                    srcloc_final)
                                self.instance_commands.append(command)

                    else:

                        # Now check for commands using this syntax:
                        #
                        # "MolNew    =    MolOld.rot(45,1,0,0).scale(100.0)"
                        #    /|\           /|\   `-----------.------------'
                        #     |             |                |
                        # child_name  parent_name    optional suffix
 
                        child_name = object_name
                        parent_name_str = next_symbol
 
                        child = StaticObj(child_name, self)
 
                        parent_name, suffix, suffix_srcloc = \
                            self._ProcessClassName(parent_name_str, lex)

                        child.class_parents.append(StrToNode(parent_name,
                                                             self,
                                                             lex.GetSrcLoc()))

                        if suffix != '':
                            # Assume the command is a StackableCommand. (This
                            # way it will enclose the commands of the parents.)
                            # Stackable commands come in (Push...Pop) pairs.
                            push_command = PushLeftCommand(suffix,
                                                           suffix_srcloc)
                            pop_command  = PopLeftCommand(push_command,
                                                          suffix_srcloc)
                            push_mod_command = ModCommand(push_command, './')
                            pop_mod_command  = ModCommand(pop_command, './')
                            child.instance_commands_push.append(push_mod_command)
                            child.instance_commands_pop.insert(0,pop_mod_command)

                            #sys.stderr.write('child.instance_commands_push = '+str(child.instance_commands_push)+'\n')

                            #sys.stderr.write('child.instance_commands_pop = '+str(child.instance_commands_pop)+'\n')
                         
                        # Check to see if this class has already been defined. 
                        # (IE. check if it present in the list of children.) 
                        for i in range(0, len(self.children)):
                            if self.children[i].name == child_name:
                                if self.children[i].IsDeleted():
                                    del self.children[i]
                                else:
                                    raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                                     '      Error near '+lex.error_leader()+'\n'
                                                     '      The name \"'+child_name+'\" is already in use.')

                        self.children.append(child)



                else:

                    # Otherwise hopefully this is a post-instance command
                    # (a command applied to a class which has been instantiated)
                    #   In that case, the object_name would be followed by
                    #   a dot and a function-call containing a '(' paren (which
                    #   would have ended up stored in the next_symbol variable).

                    open_paren_encountered = False
                    if (next_symbol == '('):
                        open_paren_encountered = True
                        lex.push_token(next_symbol) #put '(' back in the stream

                    i_dot   = object_name.rfind('.')
                    i_slash = object_name.rfind('/')
                    dot_encountered = ((i_dot != -1) and
                                       ((i_slash == -1) or (i_slash < i_dot)))

                    if (open_paren_encountered and dot_encountered and
                        (object_name[:1] != '[')):

                        obj_descr_str, suffix, suffix_srcloc = \
                            self._ExtractSuffix(object_name, lex)

                        path_tokens = obj_descr_str.split('/')
                        i_last_ptkn, staticobj = FollowPath(path_tokens, 
                                                            self, 
                                                            lex.GetSrcLoc())
                        instobj_descr_str = './'+'/'.join(path_tokens[i_last_ptkn:])

                        # I still support the "object_name.delete()" syntax for
                        # backwards compatibility.  (However newer input files
                        # use this equivalent syntax: "delete object_name")
                        if suffix == 'delete()':
                            delete_command = DeleteCommand(suffix_srcloc)
                            mod_command    = ModCommand(delete_command,
                                                        instobj_descr_str)
                            staticobj.instance_commands.append(mod_command)
                        else:
                            push_command = PushLeftCommand(suffix,
                                                           suffix_srcloc,
                                                           '.')
                            pop_command  = PopLeftCommand(push_command,
                                                          suffix_srcloc,
                                                           '.')
                            push_mod_command = ModCommand(push_command, 
                                                          instobj_descr_str)
                            pop_mod_command  = ModCommand(pop_command, 
                                                          instobj_descr_str)
                            if instobj_descr_str != './':
                                #sys.stderr.write('DEBUG: Adding '+str(push_command)+' to '+
                                #                 staticobj.name+'.instance_commands\n')
                                staticobj.instance_commands.append(push_mod_command)
                                staticobj.instance_commands.append(pop_mod_command)
                            else:
                                #sys.stderr.write('DEBUG: Adding '+str(push_command)+' to '+
                                #                 staticobj.name+'.instance_commands_push\n')
                                # CONTINUEHERE: should I make these PushRight commands and
                                #               append them in the opposite order?
                                #               If so I also have to worry about the case above.
                                staticobj.instance_commands_push.append(push_mod_command)
                                staticobj.instance_commands_pop.insert(0,pop_mod_command)

                    else:
                        # Otherwise, the cmd_token is not any of these:
                        # "write", "write_once", "create_vars"
                        # "delete", or "category".
                        # ... and it is ALSO not any of these:
                        # the name of a class (StaticObj), or
                        # the name of an instance (InstanceObj) 
                        #   followed by either a '.' or "= new"
                        #
                        # In that case, it is a syntax error:
                        raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                         '       Syntax error at or before '+lex.error_leader()+'\n'
                                         '       \"'+object_name+' '+next_symbol+'\".')

        # Keep track of the location in the user's input files 
        # where the definition of this object ends.
        self.srcloc_end = lex.GetSrcLoc()



    @staticmethod
    def CleanupReadTemplate(tmpl_contents, lex):
        #1) Remove any newlines at the beginning of the first text block
        # in tmpl_content.(Sometimes they cause ugly extra blank lines)
        assert(len(tmpl_contents) > 0)
        if isinstance(tmpl_contents[0], TextBlock):
            first_token_strip = tmpl_contents[0].text.lstrip(' ')
            if ((len(first_token_strip) > 0) and 
                (first_token_strip[0] in lex.newline)):
                tmpl_contents[0].text = first_token_strip[1:]
                tmpl_contents[0].srcloc.lineno += 1

        #2) Remove any trailing '}' characters, and complain if absent.
        #   The last token
        assert(isinstance(tmpl_contents[-1], TextBlock))
        assert(tmpl_contents[-1].text in ['}',''])
        if tmpl_contents[-1].text == '}':
            del tmpl_contents[-1]
        else:
            tmpl_begin = None
            if isinstance(tmpl_contents[0], TextBlock):
                tmpl_begin = tmpl_contents[0].srcloc
            elif isinstance(tmpl_contents[0], VarRef):
                tmpl_begin = tmpl_contents[0].srcloc
            else:
                assert(False)
            raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                             '      Error near '+lex.error_leader()+'\n\n'
                             '      Premature end to template.\n'
                             '(Missing terminator character, usually a \'}\'.) The\n'
                             'incomplete template begins near '+ErrorLeader(tmpl_begin.infile, tmpl_begin.lineno)+'\n')
        #3) Finally, if there is nothing but whitespace between the
        #   last newline and the end, then strip that off too.
        if isinstance(tmpl_contents[-1], TextBlock):
            i = len(tmpl_contents[-1].text)-1
            if i >= 0:
                while ((i >= 0) and 
                       (tmpl_contents[-1].text[i] in lex.whitespace) and
                       (tmpl_contents[-1].text[i] not in lex.newline)):
                    i -= 1
                if (tmpl_contents[-1].text[i] in lex.newline):
                    tmpl_contents[-1].text = tmpl_contents[-1].text[0:i+1]



    def LookupStaticRefs(self):
        """ Whenever the user requests to instantiate a new copy of a class,
        the name of that class is stored in self.instance_commands.
        This name is stored as a string.  After all of the classes have been
        defined, then we go back through the tree and replace these names
        with pointers to actual StaticObjs which correspond to those classes.
        (This was deferred until all of the classes have been defined so 
        that users can refer to classes that they will define later on.)

        """

        # Now do the same for any children which 
        # are created during instantiation: 
        for command in self.instance_commands:
            # Does this command create/instantiate a new copy of a class?
            if isinstance(command, InstantiateCommand):
                # If so, figure out which statobj is referred to by statobj_str.
                assert(isinstance(command.class_ref.statobj_str, basestring))
                command.class_ref.statobj = StrToNode(command.class_ref.statobj_str,
                                                      self,
                                                      command.class_ref.srcloc)
                
        # Now recursively resolve StaticObj pointers for the "children"
        # (in this case, "children" refers to classes whose definitions 
        #  are nested within this one).
        for child in self.children:
            child.LookupStaticRefs()




    def _ExtractSuffix(self, class_name_str, lex):
        """ 

        This ugly function helps process "new" commands such as:
        mola = new ForceFieldA/../MoleculeA.move(30,0,0).rot(45,0,0,1)
        This function expects a string,
        (such as "ForceFieldA/../MoleculeA.move(30,0,0).rot(45,0,0,1)")
        It extracts the class name "ForceFieldA/../MoleculeA"
        and suffix "move(30,0,0).rot(45,0,0,1)"
        """
        # Dots in class names can appear for 2 reasons:
        #   1) as part of a path like "../" describing the location 
        #      where this class was defined relative to the caller.
        #      In that case it will be preceeded or followed by
        #      either another dot '.', or a slash '/'
        #   2) as part of a "suffix" which appears after the name
        #      containing instructions which modify how to
        #      instantiate that class.
        # Case 1 is handled elsewhere.  Case 2 is handled here.
        i_dot = 0
        while i_dot < len(class_name_str):

            i_dot = class_name_str.find('.', i_dot)

            if i_dot == -1:
                break
            # Is the '.' character followed by another '.', as in ".."?
            # If so, it's part of a path such as "../Parent/Mol', (if 
            # so, it's not what we're looking for, so keep searching)
            if i_dot < len(class_name_str)-1:
                if class_name_str[i_dot+1] == '.':
                    i_dot += 1
                    #otherwise, check to see if it is followed by a '/'?
                elif class_name_str[i_dot+1] != '/':
                    # if not, then it must be part of a function name
                    break;

        class_suffix = ''
        class_name = class_name_str
        class_suffix_srcloc = None

        if ((i_dot != -1) and 
            (i_dot < len(class_name_str))):
            class_suffix = class_name_str[i_dot:]
            class_name   = class_name_str[:i_dot]
            if class_name_str[-1] != ')':
                # If it does not already contains the parenthesis?
                class_suffix += lex.GetParenExpr()
                class_suffix_srcloc = lex.GetSrcLoc()
            #sys.stderr.write('  splitting class name into class_name.suffix\n'
            #                 '  class_name=\"'+class_name+'\"\n'
            #                 '      suffix=\"'+class_suffix+'\"\n')

        return class_name, class_suffix.lstrip('.'), class_suffix_srcloc




    def _ProcessClassName(self, class_name_str, lex):
        """

        This function does some additional 
        processing (occasionaly inserting "..." before class_name).
        """

        class_name, class_suffix, class_suffix_srcloc = \
            self._ExtractSuffix(class_name_str, lex)

        # ---- ellipsis hack ----
        #   (Note-to-self 2012-4-15)
        # Most users expect ttree.py to behave like a 
        # standard programming language: If the class they are 
        # instantiating was not defined in this specific 
        # location, they expect ttree.py to search for
        # it outwards, first in the parent's environment, 
        # and then in the parent's parent's environment, 
        # and so on, until the object is found.
        # For example, most users expect this to work:
        # class A{
        #   <definition_of_a_goes_here...>
        # }
        # class B{
        #   a = new A
        # }
        # Notice in the example above we did not have to specify where "A"
        # was defined, because it is defined in the parent's
        # environment (ie. immediately outside B's environment).
        # 
        #    One can obtain the equivalent behavior in ttree.py 
        # using ellipsis syntax:     "a = new .../A" symbol. 
        # The ellipsis ".../" tells ttree.py to search upwards
        # for the object to the right of it ("A")
        #    In order to make ttree.py behave the way 
        # most users are expecting, we artificially insert a 
        # ".../" before the class name here.  (Later on, the
        # code that processes the ".../" symbol will take
        # care of finding A.)

        if (len(class_name)>0) and (class_name[0] not in ('.','/','*','?')):
            class_name = '.../' + class_name

        return class_name, class_suffix, class_suffix_srcloc





    def _ParseRandom(self, lex):
        bracket1 = lex.get_token()
        bracket2 = lex.get_token()
        if ((bracket1 != '(') and (bracket1 != '[')):
            raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                             '       Error near '+lex.error_leader()+'\n'
                             'Expected a \"([\" following '+class_name+'.')
        class_names = []
        token = ''
        prev_token = '['
        while True:
            token = lex.get_token()
            if (token == '('):
                lex.push_token(token)
                token = lex.GetParenExpr()
                if (prev_token not in (',','[','(')):
                    assert(len(class_names) > 0)
                    class_names[-1] = prev_token + token
                    prev_token = prev_token + token
                else:
                    class_names.append(token)
                    prev_token = token
            else:
                if ((token == ']') or
                    (token == lex.eof) or
                    (token == '}') or
                    ((token in lex.wordterminators) and 
                     (token != ','))):
                    if (prev_token in (',','[','(')):
                        class_names.append('')
                    break
                if token != ',':
                    class_names.append(token)
                elif (prev_token in (',','[','(')):
                    class_names.append('')
                prev_token = token


        token_comma = lex.get_token()
        bracket1 = lex.get_token()
        if ((token != ']') or 
            (token_comma != ',') or
            (bracket1 != '[')):
            raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                             '       Error near '+lex.error_leader()+'\n'
                             'Expected a list of class names enclosed in [] brackets, followed by\n'
                             'a comma, and then a list of probabilities also enclosed in [] brackets.\n'
                             '(A random-seed following another comma is optional.)')

        weights = []
        while True:
            token = lex.get_token()
            if ((token == ']') or
                (token == lex.eof) or
                (token == '}') or
                ((token in lex.wordterminators) and
                 (token != ','))):
                break
            if token != ',':
                try:
                    weight = float(token)
                except ValueError:
                    raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                     '       Error near '+lex.error_leader()+'\n'
                                     '       \"'+token+'\"\n'
                                     'Expected a list of numbers enclosed in [] brackets.')
                if (weight < 0.0):
                    raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                     '       Error near '+lex.error_leader()+'\n'
                                     '       Negative numbers are not allowed in \"random(\" argument list.\n')
                weights.append(weight)

        bracket2 = lex.get_token()
        if ((token != ']') or 
            (bracket2 not in (')',','))):
            raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                             '       Error near '+lex.error_leader()+'\n'
                             'Expected a \")\" or a \",\" following the list of numeric weights.')

        if len(class_names) != len(weights):
            raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                             '       Error near '+lex.error_leader()+'\n'
                             'Unequal number of entries in object list and probability list.\n')

        tot_weight = sum(weights)
        if (tot_weight <= 0.0):
                if (weight < 0.0):
                    raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                 '       Error near '+lex.error_leader()+'\n'
                                     '       The numbers in the \"random(\" argument list can not all be zero.\n')
        for i in range(0,len(weights)):
            weights[i] /= tot_weight


        if bracket2 == ',':
            try:
                token = lex.get_token()
                seed = int(token)
                random.seed(seed)
            except ValueError:
                raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                 '       Error near '+lex.error_leader()+'\n'
                                 '       \"'+token+'\"\n'
                                 'Expected an integer (a seed) following the list of weights.')
            bracket2 = lex.get_token()
            if (bracket2 != ')'):
                raise InputError('Error('+g_module_name+'.StaticObj.Parse()):\n'
                                 '       Error near '+lex.error_leader()+'\n'
                                 '       \"'+token+'\"\n'
                                 'Expected a \")\".')
        else:
            random.seed()


        return (class_names, weights)


    def BuildCommandList(self, command_list):
        """ 
        Search the commands in the tree and make a linear list of commands
        in the order they should be carried out.

        """

        if self.IsDeleted():
            return

        # Add a special note to the list of commands to indicate which object
        # the commands refer to.  (This might be useful one day.)
        # Later we can loop through this command list and still be able to tell 
        # whether or not we are within the scope of a particular class or instance
        # (by seeing if we are between a "ScopeBegin" and "ScopeEnd" pair).
        command_list.append(ScopeBegin(self, self.srcloc_begin))

        # We want to append commands to the command_list in the same order
        # that these commands appear in the user's input files.

        #    Unfortunately the commands may be interspersed with the creation of
        # new StaticObjs which have their own commands which we have to explore
        # recursively. 
        # Fortunately each child (StaticObj) has a srcloc_begin member, so we
        # can infer the correct order of all the commands belonging to the 
        # children and correctly insert them into the correct place in between
        # the commands of the parent.

        srcloc2command_or_child = {}

        for command in self.commands:
            srcloc2command_or_child[command.srcloc] = command

        for child in self.children:
            srcloc = child.srcloc_begin
            # special case: Some children do not have a srcloc because 
            # they were generated automatically.  These children should 
            # not have any commands either so we can ignore them.
            if srcloc != None:
                srcloc2command_or_child[srcloc] = child
            else:
                assert(len(child.commands) == 0)

        for srcloc in sorted(srcloc2command_or_child.keys()):
            entry = srcloc2command_or_child[srcloc]
            if isinstance(entry, StaticObj):
                child = entry
                child_commands = []
                child.BuildCommandList(child_commands)
                command_list += child_commands
            else:
                command_list.append(entry)

        command_list.append(ScopeEnd(self, self.srcloc_end))




    def __str__(self):
        out_str = self.name
        if len(self.children) > 0:
            out_str += '('
            for i in range(0,len(self.children)):
                if i+1 < len(self.children):
                    out_str += str(self.children[i])+', '
                else:
                    out_str += str(self.children[i])+')'
        return out_str






def RandomSelect(entries, weights):
    """ Return an entry from a list at random using
        a (normalized) list of probabilities. """
    assert(len(entries) == len(weights))
    x = random.random()
    i = 0
    tot_probability = 0.0
    while i < len(weights)-1:
        tot_probability += weights[i]
        if x <= tot_probability:
            break
        i += 1
    return entries[i]




class InstanceObjBasic(object):
    """ A simplified version of InstanceObj.
        See the documentation/comments for InstanceObj for more details.
        (Leaf nodes (variables) are typically stored as InstanceObjBasic objects
         More general, non-leaf nodes are stored using InstanceObj objects.)

    """

    __slots__=["name","parent"]

    def __init__(self, 
                 name = '', 
                 parent = None):
        self.parent = parent  # the environment/object which created this object
                              # Example:
                              # Suppose this "molecule" is an amino acid monomer
                              # belonging to a protein.  The "parent" refers to
                              # the InstanceObj for the protein.  ".parent" is
                              # useful for traversing the global instance tree.
                              #   (use InstanceObj.statobj.parent for 
                              #   traversing the global static tree)

        self.name   = name    # A string uniquely identifying this object in
                              # in it's "parent" environment.
                              # (It is always the same for every instance
                              #  of the parent object.  It would save memory to
                              #  get rid of this member.  Andrew 2012/9/13)

        ##vb##self.var_bindings=None  # List of variables assigned to this object
        ##vb##                        # or None (None takes up less space than an
        ##vb##                        # empty list.)


    ##vb##def AddVarBinding(self, var_binding):
    ##vb##    if self.var_bindings is None:
    ##vb##        self.var_bindings = [var_binding]
    ##vb##    else:
    ##vb##        self.var_bindings.append(var_binding)


    def Dealloc(self):
        pass
        ##vb##if self.var_bindings is None:
        ##vb##    return
        ##vb##N = len(self.var_bindings)-1
        ##vb##for i in range(0,len(self.var_bindings)):
        ##vb##    vb = self.var_bindings[N-i]
        ##vb##    cat = vb.category
        ##vb##    assert(self in cat.bindings)
        ##vb##    del cat.bindings[self]
        ##vb##    del self.var_bindings[N-i]
        ##vb##self.var_bindings = None


    def DeleteSelf(self):
        self.Dealloc()
        self.parent = self  # This condition (normally never true)
                            # flags the node as "deleted".  (Nodes are never
                            # actually deleted, just flagged.)
                            # I used to have a separate boolean member variable
                            # which was set True when deleted, but I started 
                            # eliminated unnecessary data members to save space.

    def IsDeleted(self):
        # Return true if self.parent == self
        # for this node (or for any ancestor node).
        node = self
        while node.parent != None:
            if node.parent == node:
                return True
            node = node.parent
        return False



class InstanceObj(InstanceObjBasic):
    """  InstanceObjs are used to store nodes in the global
    "instance tree", the tree of all classes (molecules) which have
    been instantiated.  Recall that whenever a class is instantiated,
    it's members will be instantiated as well.  Since these
    members can also be classes, this relationship is hierarchical, 
    and can be represented as a tree.
    "InstanceObjs" are the data type used to store the nodes in that tree."""

    __slots__=["statobj",
               "children",
               "categories",
               "commands",
               "commands_push",
               "commands_pop",
               "srcloc_begin",
               "srcloc_end"]


    def __init__(self, 
                 name = '', 
                 parent = None):
        
        InstanceObjBasic.__init__(self, name, parent)

        self.statobj = None       # The statobj node refered to by this instance
        self.children = []       # A list of statobjs corresponding to
                                 # constituent parts (members) of the
                                 # current class instance.
                                 # The typical example is to consider the 
                                 # multiple amino acids (child-molecules)
                                 # which must be created in order to create a 
                                 # new protein (instance) to which they belong 
                                 # (which would be "self" in this example)
        
        self.categories = OrderedDict()         # This member stores the same data as the 
                                     # Instance variables (ie. variables 
                                     # with a '$' prefix) are stored in a
                                     # category belonging to node.categories
                                     # where "node" is of type InstanceObj.
                                     # (There is a long explanation of 
                                     # "categories" in the comments 
                                     # of class StaticObj.)
        
                                     # always equal self.statobj.categories.
                                     # It may seem redundant to make a 
                                     # new member to store the same data, 
                                     # but I do it because then the two 
                                     # types of trees (StaticObjs, InstanceObjs) have
                                     # the same member names (although different
                                     # connectivity).  This way I can write
                                     # simple functions (eg FindCatNode()) 
                                     # which can traverse through either 
                                     # type of tree.

        self.commands = []           # An ordered list of commands to carry out
                                     # during instantiation

        self.commands_push = []      # Stackable commands to carry out (first, before children)
        self.commands_pop = []       # Stackable commands to carry out (last, after children)


        self.srcloc_begin = None     # Keep track of location in user files
        self.srcloc_end   = None     # (useful for error message reporting)




    def LookupMultiDescrStr(self, 
                            multi_descr_str, 
                            srcloc,
                            null_list_warning=False,
                            null_list_error=False):
        """ 
        Post-Instance (PI) modifiers/commands are commands which modify
        an instance of a class after it has already been instantiated.

        Simple Example:

        class A {
           ...
        }
        class B {
           a = new A.command_1()
           a.command_2()
        }

        In the example above "command_2()" is a ModCommand, and
        "a" is the multi_descr_str (string describing the correspond InstanceObj).
        The "command_2()" command will be retroactively pushed onto the
        list of commands to execute once "a" is instantiated. 
        (This is somewhat counter-intuitive.)

        When array brackets [] and wildcards are used, a single ModCommand
        can modify many different instances, for example suppose:

        a = new A [2][5][3]

           then "a[1][2][*].command_3()" is equivalent to

        a[0][2][0].command_3()
        a[0][2][1].command_3()
        a[0][2][2].command_3()

        In this example "a[1][2][*]" is the multi_descr_str

           "a[*][3][*].command_4()" is equivalent to

        a[0][3][0].command_4()
        a[0][3][1].command_4()
        a[1][3][0].command_4()
        a[1][3][1].command_4()

        In this function, we interpret strings like "a" and "a[*][3][*]"
        in the examples above, and figure out which InstanceObjs they refer to,
        and push the corresponding command into that InstanceObjs instance 
        command stack retroactively.

           In addition to [*], you can use [a-b] and [a:b] syntax. For example:
           "a[0][1-2][0-1].command_3()" and
           "a[0][1:3][0:2].command_3()" are both equivalent to:

        a[0][1][0].command_3()
        a[0][1][1].command_3()
        a[0][2][0].command_3()
        a[0][2][1].command_3()

        """

        pattern_str = multi_descr_str

        # Suppose pattern_str = 'a[1][*][3]/b[**][2]'
        # We want to split this string into a list of string fragments
        # which omits the '*' characters:  [ 'a[',  '][3]/b',  '][2]' ]
        # However, we only want to do this when * is enclosed in [].
        pattern_fragments = []
        ranges_ab         = []
        i_close_prev = 0
        i_close = 0
        i_open = 0
        while True:
            i_open = pattern_str.find('[', i_open+1)
            if i_open == -1:
                pattern_fragments.append(pattern_str[i_close_prev:])
                break
            else:
                i_close = pattern_str.find(']', i_open+1)
                if i_close == -1:
                    pattern_fragments.append(pattern_str[i_close_prev:])
                    break

                # If there is a '*' or a ':' character between
                # the [] brackets, then split the string at '['
                # (at i_open) and resume reading again at ']' 
                # (at i_close) (and create a new entry in the
                #  pattern_fragments[] and ranges_ab[] lists)
                wildcard_here = True
                range_ab = [0,-1]
                for j in range(i_open+1, i_close):
                    if ((pattern_str[j] == ':') or 
                        ((pattern_str[j] == '-') and (j > i_open+1)) or
                        (pattern_str[j] == '*')):
                        i_wildcard = len(pattern_fragments)
                        range_a_str = pattern_str[i_open+1 : j]
                        range_b_str = pattern_str[j+1 : i_close]
                        if (range_a_str != ''):
                            if str.isdigit(range_a_str):
                                range_ab[0] = int(range_a_str)
                            else:
                                raise InputError('Error near '+
                                                 ErrorLeader(srcloc.infile,
                                                             srcloc.lineno)+'\n'
                                                 '   Expected colon-separated integers.\n')
                        if (range_b_str != ''):
                            if str.isdigit(range_b_str):
                                range_ab[1] = int(range_b_str)
                                # special case: When [a-b] type syntax is 
                                #   used, it selects from a to b inclusive.
                                #   (IE. b is not a strict upper bound.)
                                if pattern_str[j] == '-':
                                    range_ab[1] += 1
                            else:
                                raise InputError('Error near '+
                                                 ErrorLeader(srcloc.infile,
                                                             srcloc.lineno)+'\n'
                                                 '   Expected colon-separated integers.\n')
                        break
                    elif j == i_close-1:
                        wildcard_here = False

                if wildcard_here:
                    pattern_fragments.append(pattern_str[i_close_prev:i_open+1])
                    ranges_ab.append(range_ab)
                    i_close_prev = i_close

        assert(len(pattern_fragments)-1==len(ranges_ab))
        # Now figure out which InstanceObj or InstanceObjs correspond to 
        # the name or set of names suggested by the multi_descr_str, 
        # (after wildcard characters have been substituted with integers).

        instobj_list = []
        if len(pattern_fragments) == 1:
            # commenting out:
            # instobj_list.append(StrToNode(pattern_str, self, srcloc))
            #
            # Line above will print an error message if the node is not found.
            # However sometimes we don't want this.  Use this code instead:
            path_tokens = pattern_str.split('/')
            i_last_ptkn, instobj = FollowPath(path_tokens, 
                                              self, 
                                              srcloc)
            # If found add to instobj_list
            if ((i_last_ptkn == len(path_tokens)) 
                and (not instobj.IsDeleted())): # (make sure not "deleted")
                instobj_list.append(instobj)
        else:
            # num_counters equals the number of bracket-enclosed wildcards
            num_counters= len(pattern_fragments)-1 
            multi_counters = [ranges_ab[i][0] for i in range(0, num_counters)]
            all_matches_found = False
            d_carry = 0            
            while d_carry < num_counters:
                    
                # Find the next InstanceObj in the set of InstanceObjs which
                # satisfy the wild-card pattern in pattern_fragments.

                while d_carry < num_counters:

                    candidate_descr_str = ''.join([pattern_fragments[i] +
                                                   str(multi_counters[i])
                                                   for i in range(0,num_counters)] \
                                                      + \
                                                      [pattern_fragments[num_counters]])

                    #sys.stderr.write('DEBUG: /'+self.name+
                    #                 '.LookupMultiDescrStr()\n'
                    #                 '       looking up \"'+
                    #                 candidate_descr_str+'\"\n')

                    path_tokens = candidate_descr_str.split('/')
                    i_last_ptkn, instobj = FollowPath(path_tokens, 
                                                      self, 
                                                      srcloc)

                    # If there is an InstanceObj with that name,
                    # then add it to the list of InstanceObjs to
                    # which we will apply this modifier function,
                    # and increment the counters

                    # If found (and if the counter is within the range)...
                    if ((i_last_ptkn == len(path_tokens)) and
                        ((ranges_ab[d_carry][1] == -1) or
                         (multi_counters[d_carry]<ranges_ab[d_carry][1]))):
                        # (make sure it has not yet been "deleted")
                        if (not instobj.IsDeleted()):
                            instobj_list.append(instobj)
                        d_carry = 0
                        multi_counters[0] += 1
                        #sys.stderr.write('DEBUG: InstanceObj found.\n')
                        break


                    # If there is no InstanceObj with that name,
                    # then perhaps it is because we have incremented
                    # the counter too high.  If there are multiple 
                    # counters, increment the next most significant
                    # counter, and reset this counter to 0.
                    # Keep looking
                    # (We only do this if the user neglected to explicitly
                    #  specify an upper bound --> ranges_ab[d_carry[1]==-1)
                    elif ((ranges_ab[d_carry][1] == -1) or
                          (multi_counters[d_carry]>=ranges_ab[d_carry][1])):
                        #sys.stderr.write('DEBUG: InstanceObj not found.\n')
                        multi_counters[d_carry] = ranges_ab[d_carry][0]
                        d_carry += 1
                        if d_carry >= num_counters:
                            break
                        multi_counters[d_carry] += 1
                    else:
                        # Object was not found but we keep going.  Skip 
                        # to the next entry in the multi-dimensional list.
                        d_carry = 0
                        multi_counters[0] += 1
                        break

        if (null_list_warning and (len(instobj_list) == 0)):
            sys.stderr.write('WARNING('+g_module_name+'.LookupMultiDescrStr()):\n'
                             '       Potential problem near '+
                             ErrorLeader(srcloc.infile,
                                         srcloc.lineno)+'\n'
                             '       No objects (yet) matching name \"'+pattern_str+'\".\n')
        if (null_list_error and 
            (len(instobj_list) == 0)):
            if len(pattern_fragments) == 1:
                raise InputError('Error('+g_module_name+'.LookupMultiDescrStr()):\n'
                                 '       Syntax error near '+
                                 ErrorLeader(srcloc.infile,
                                             srcloc.lineno)+'\n'
                                 '       No objects matching name \"'+pattern_str+'\".')
            else:
                sys.stderr.write('WARNING('+g_module_name+'.LookupMultiDescrStr()):\n'
                                 '       Potential problem near '+
                                 ErrorLeader(srcloc.infile,
                                             srcloc.lineno)+'\n'
                                 '       No objects (yet) matching name \"'+pattern_str+'\".\n')
                

        return instobj_list





    def __str__(self):
        if len(self.children) > 0:
            out_str = self.name + '('
            for i in range(0,len(self.children)):
                if i+1 < len(self.children):
                    out_str += str(self.children[i])+', '
                else:
                    out_str += str(self.children[i])+')'
        else:
            out_str = self.name

        return out_str



    def Dealloc(self):
        InstanceObjBasic.Dealloc(self)
        # Trying to remove pointers and variables.
        # Hope I'm doing it correctly
        #sys.stderr.write(self.name+'.Dealloc() invoked.\n')
        N = len(self.commands)-1
        for i in range(0, len(self.commands)):
            del self.commands[N-i]
        N = len(self.children)-1
        for i in range(0, len(self.children)):
            child = self.children[N-i]
            child.Dealloc()
            del self.children[N-i]



    #def __del__(self):
    #    sys.stderr.write(self.name+'.__del__() invoked.\n')


    def DeleteSelf(self):
        #  Modification:  Don't get rid of pointers to yourself.  Knowing which
        #                 objects you instantiated and destroyed might be useful
        #                 in case you want to apply multiple delete [*] commands
        #  Commenting out:
        #assert(parent != None)
        #i = 0
        #while i < len(parent.children):
        #    if command.instobj == parent.children[i]:
        #        del parent.children[i]
        #    else:
        #        i += 1

        self.Dealloc()
        InstanceObjBasic.DeleteSelf(self)




    def BuildInstanceTree(self,
                          statobj,
                          class_parents_in_use):
        """
        This takes care of the details of copying relevant data from an StaticObj
        into a newly-created InstanceObj.  It allocates space for and performs 
        a deep-copy of any instance variables (and new instance categories), but
        it performs a shallow copy of everything else (template text, etc..).
        This is done recursively for every child that this class instantiates.

        """

        if self.IsDeleted():
            return

        #sys.stderr.write('  DEBUG: '+self.name+
        #                 '.BuildInstanceTree('+statobj.name+')\n')

        #instance_refs = {}
        # Keep track of which line in the file (and which file) we were 
        # in when we began parsing the class which defines this instance, 
        # as well as when we stopped parsing.
        # (Don't do this if you are recusively searching class_parents because
        #  in that case you would be overwritting .statobj with with the parent.)
        if len(class_parents_in_use) == 0:
            self.statobj      = statobj
            self.srcloc_begin = statobj.srcloc_begin
            self.srcloc_end   = statobj.srcloc_end

        # Make copies of the class_parents' StaticObj data.

        # First deal with the "self.instance_commands_push"
        # These commands should be carried out before any of the commands 
        # in "self.instance_commands".
        for command in statobj.instance_commands_push:
            self.ProcessCommand(command)

        # Then deal with class parents
        for class_parent in statobj.class_parents:
            # Avoid the "Diamond of Death" multiple inheritance problem
            if class_parent not in class_parents_in_use:
                #sys.stderr.write('  DEBUG: '+self.name+'.class_parent = '+
                #                 class_parent.name+'\n')
                self.BuildInstanceTree(class_parent,
                                       class_parents_in_use)
            class_parents_in_use.add(class_parent)

        # Now, deal with the data in THIS object and its children
        assert((self.commands != None) and (self.categories != None))

        # "instance_categories" contains a list of new "categories" (ie new 
        # types of variables) to create whenever this class is instantiated.
        # (This is used whenever we create a local counter variable: Suppose we
        #  want to count the residues within a particular protein, when there
        #  are multiple copies of the same protein in the simulation.)
        for cat_name, cat in statobj.instance_categories.items():
            assert(len(cat.bindings) == 0)
            self.categories[cat_name] = Category(cat_name)
            self.categories[cat_name].counter = cat.counter.__copy__()
            # Note: Later on we will generate leaf nodes corresponding to
            #       variables, and put references to them in this category.

        # Deal with the "instance_commands",
        for command in statobj.instance_commands:
            self.ProcessCommand(command)


        # Finally deal with the "self.instance_commands_pop"
        # These commands should be carried out after all of the commands 
        # in "self.instance_commands".
        for command in statobj.instance_commands_pop:
            self.ProcessCommand(command)




    def ProcessCommand(self, command):

        if isinstance(command, ModCommand):
            sys.stderr.write('  processing command \"'+str(command)+'\"\n')
            mod_command = command
            instobj_list = self.LookupMultiDescrStr(mod_command.multi_descr_str,
                                                    mod_command.command.srcloc)
            if len(instobj_list) == 0:
                if isinstance(mod_command.command, DeleteCommand):
                    # It's possible that some nodes we want to delete have 
                    # not yet been created yet.  Deal with these later.
                    # (See InvokeAllDeletes().)
                    self.commands.append(mod_command.__copy__())
                else:
                    raise InputError('Error('+g_module_name+'.ProcessCommand()):\n'
                                     '       Syntax error at or before '+
                                     ErrorLeader(mod_command.command.srcloc.infile,
                                                 mod_command.command.srcloc.lineno)+'\n'
                                     '       No objects matching name \"'+
                                     mod_command.multi_descr_str+'\"\n'
                                     '          (If the object is an array, include brackets. Eg. \"molecules[*][*][*]\")')
            else:
                for instobj in instobj_list:
                    if isinstance(mod_command.command, DeleteCommand):
                        # We need the instance tree to be fully defined after
                        # this function is called, so carry out the deletion now
                        instobj.DeleteSelf()
                    else:
                        command = mod_command.command.__copy__()
                        self.ProcessContextNodes(command)
                        if isinstance(command, PushCommand):
                            instobj.commands_push.append(command)
                        elif isinstance(mod_command.command, PopCommand):
                            instobj.commands_pop.insert(0, command)
                        else:
                            # I don't know if any other types commands will ever
                            # occur but I handle them below, just in case...
                            assert(not isinstance(command, InstantiateCommand))
                            instobj.commands.append(command.__copy__())

            return

        # Not needed:
        #if isinstance(command, DeleteCommand):
        #    self.DeleteSelf()
        #    return  # (Note: we do not append this command to self.commands)

        # Otherwise:
        command = command.__copy__()
        self.ProcessContextNodes(command)
        if isinstance(command, InstantiateCommand):
            sys.stderr.write('  processing command \"'+str(command)+'\"\n')
            self.commands.append(command) # <- useful later to keep track of the
                                          #    order that children were created

            child = InstanceObj(command.name, self)
            if command.class_ref.statobj_str == '':
                child.DeleteSelf()
                # Why?  This if-then check handles the case when the user
                # wants to create an array of molecules with random vacancies.
                # When this happens, some of the instance commands will
                # contain instructions to create a copy of a molecule with
                # an empty molecule-type-string (class_ref.statobj_str).
                # Counter-intuitively, ...
                #  ...we DO want to create something here so that the user can
                # safely loop over the array without generating an error.
                # (Such as to delete elements, or move the remaining
                #  members in the array.)  We just want to mark it as
                # 'deleted'.  (That's what "DeleteSelf()" does.)
            else:
                # This is the heart of "BuildInstanceTree()"
                # (which implements object composition)
                new_class_parents_in_use = set([])
                child.BuildInstanceTree(command.class_ref.statobj,
                                        new_class_parents_in_use)
            self.children.append(child)
        else:
            # Otherwise, we don't know what this command is yet.
            # Append it to the list of commands and process it/ignore it later.
            self.commands.append(command)


    def ProcessContextNodes(self, command):
        if hasattr(command, 'context_node'):
            # Lookup any nodes pointers to instobjs
            if command.context_node != None:
                if type(command.context_node) is str:
                    command.context_node = StrToNode(command.context_node,
                                                     self,
                                                     command.srcloc)
            # (Otherwise, just leave it as None)



    def BuildCommandList(self, command_list):
        """ 
        Search the commands in the tree and make a linear list of commands
        in the order they should be carried out.

        """

        if self.IsDeleted():
            return

        if (len(self.commands) == 0):
            # To save memory don't generate any commands
            # for trivial (leaf) nodes
            assert(len(self.children) == 0)
            return

        # Add a special note to the list of commands to indicate which object
        # the commands refer to.  (This might be useful one day.)
        # Later we can loop through this command list and still be able to tell 
        # whether or not we are within the scope of a particular class or instance
        # (by seeing if we are between a "ScopeBegin" and "ScopeEnd" pair).

        command_list.append(ScopeBegin(self, self.srcloc_begin))
        # Note:
        # The previous version looped over all commands in this node, and then
        # recursively invoke BuildCommandList() on all the children of this node
        # We don't do that anymore because it does not take into account the
        # order that various child objects were created/instantiated
        # which potentially could occur in-between other commands.  Instead,
        # now we loop through the command_list and recursively visit child
        # nodes only when we encounter them in the command list.

        for command in self.commands_push:
            assert(isinstance(command, InstantiateCommand) == False)
            command_list.append(command)

        i_child = 0
        for command in self.commands:
            if isinstance(command, InstantiateCommand):
                child = self.children[i_child]
                child.BuildCommandList(command_list)
                i_child += 1
            else:
                command_list.append(command)

        for command in self.commands_pop:
            assert(isinstance(command, InstantiateCommand) == False)
            command_list.append(command)

        command_list.append(ScopeEnd(self, self.srcloc_begin))



def AssignVarPtrs(context_node, search_instance_commands = False):
    """ 
       Now scan through all the variables within the templates defined 
    for this context_node (either static or dynamic depending on var_filter). 
    Each reference to a variable in the template has a descriptor which 
    indicates the variable's type, and in which molecule it is defined (ie 
    where it is located in the molecule instance tree or type definition tree).
    (See comments for "class VarNPtr(object):" above for details.)
       Also check to see the nodes in the tree that are refered to by the 
    descriptor exist, and, if so, assign pointers to these nodes. 
    If not, then (in some cases), create new node whose names match the 
    variable descriptor names (and assign pointers to them). 
    Also figure out the category node and name which specify the type of 
    variable.  These three pieces of information 
    This node (it's name and location in the tree) essentially identifies the 
    variable, (like a variable name) but it does not store anything.
         Also: Eventually we want to assign a value to each variable. 
    This same variable (node) may appear multiple times in diffent templates. 
    So we also create a place to store this variable's value, and also assign 
    (two-way) pointers from the VarRef in the template, to this storage area
    so that later on when we write out the contents of the template to a file,
    we can substitute this variable with it's value.
    """

    #sys.stdout.write('AssignVarPtrs() invoked on node: \"'+NodeToStr(context_node)+'\"\n')

    if search_instance_commands:
        assert(isinstance(context_node, StaticObj))
        commands = context_node.instance_commands
    else:
        # Note: Leaf nodes contain no commands, so skip them
        if (not hasattr(context_node, 'commands')):
            return
        # Otherwise process their commands
        commands = context_node.commands

    for command in commands:
        if isinstance(command, WriteFileCommand):

            for var_ref in command.tmpl_list:
                # Process the VarRef entries in the tmpl_list, 
                #   (and check they have the correct prefix: either '$' or '@')
                # Ignore other entries (for example, ignore TextBlocks).

                if (isinstance(var_ref, VarRef) and
                    ((isinstance(context_node, StaticObj) and (var_ref.prefix[0] == '@')) or
                     (isinstance(context_node, InstanceObjBasic) and (var_ref.prefix[0] == '$')))):

                    var_ref.nptr.cat_name, var_ref.nptr.cat_node, var_ref.nptr.leaf_node = \
                        DescrToCatLeafNodes(var_ref.descr_str,
                                            context_node,
                                            var_ref.srcloc,
                                            True)

                    categories = var_ref.nptr.cat_node.categories

                    # "categories" is a dictionary storing "Category" objects
                    # indexed by category names.

                    # Note to self:  Always use the ".categories" member, 
                    #  (never the ".instance_categories" member.
                    #  ".instance_categories" are only used temporarilly before 
                    #  we instantiate, ie. before we build the tree of InstanceObjs.)


                    category = categories[var_ref.nptr.cat_name]
                    # "category" is a Category object containing a
                    # dictionary of VarBinding objects, and an internal counter.
                    
                    var_bindings = category.bindings
                    # "var_bindings" is a dictionary storing "VarBinding" 
                    # objects, indexed by leaf nodes.  Each leaf node 
                    # corresponds to a unique variable in this category.

                    # --- Now update "var_bindings" ---

                    # Search for the "VarBinding" object that
                    # corresponds to this leaf node.
                    # If not found, then create one.

                    if var_ref.nptr.leaf_node in var_bindings:
                        var_binding = var_bindings[var_ref.nptr.leaf_node]
                        # "var_binding" stores the information for a variable,
                        # including pointers to all of the places the variable 
                        # is rerefenced, the variable's (full) name, and value.
                        #
                        # Keep track of all the places that varible is
                        # referenced by updating the ".refs" member
                        var_binding.refs.append(var_ref)
                    else:
                        # Not found, so we create a new binding.
                        var_binding       = VarBinding()

                        # var_binding.refs contains a list of all the places
                        # this variable is referenced. Start with this var_ref:
                        var_binding.refs  = [var_ref]

                        # keep track of the cat_node, cat_name, leaf_node:
                        var_binding.nptr = var_ref.nptr

                        # "var_binding.full_name" stores a unique string like 
                        #   '@/atom:Water/H' or '$/atom:water[1423]/H2', 
                        # which contains the full path for the category and leaf
                        # nodes, and uniquely identifies this variable globally.
                        # Thus these strings correspond uniquely (ie. in a 
                        # one-to-one fashion) with the nodes they represent.


                        var_binding.full_name = var_ref.prefix[0] + \
                            CanonicalDescrStr(var_ref.nptr.cat_name,
                                              var_ref.nptr.cat_node,
                                              var_ref.nptr.leaf_node,
                                              var_ref.srcloc)
                        # (These names can always be generated later when needed
                        #  but it doesn't hurt to keep track of it here too.)


                        # Now add this binding to the other
                        # bindings in this category:
                        var_bindings[var_ref.nptr.leaf_node] = var_binding

                        ##vb## var_ref.nptr.leaf_node.AddVarBinding(var_binding)

                        var_binding.category = category

                    # It's convenient to add a pointer in the opposite direction
                    # so that later if we find the var_ref, we can find its
                    # binding and visa-versa. (Ie. two-way pointers)
                    var_ref.binding = var_binding

                    assert(var_ref.nptr.leaf_node in var_bindings)



    # Recursively invoke AssignVarPtrs() on all (non-leaf) child nodes:
    for child in context_node.children:
        AssignVarPtrs(child, search_instance_commands)





def AssignVarOrderByFile(context_node, search_instance_commands=False):
    """
    For each category in context_node, and each variable in that category,
    set the order of each variable equal to the position of that variable 
    in the user's input file.

    """

    if search_instance_commands:
        assert(isinstance(context_node, StaticObj))
        commands = context_node.instance_commands_push + \
                   context_node.instance_commands + \
                   context_node.instance_commands_pop
    else:
        commands = context_node.commands
    for command in commands:
        if isinstance(command, WriteFileCommand):
            tmpl_list = command.tmpl_list
            for var_ref in tmpl_list:
                if isinstance(var_ref, VarRef):
                    if (((var_ref.prefix == '@') and
                         isinstance(context_node, StaticObj)) or
                        ((var_ref.prefix == '$') and
                         isinstance(context_node, InstanceObjBasic))):
                    #if ((var_ref.prefix == '@') or
                    #    (not search_instance_commands)):
                        if ((var_ref.binding.order is None) or 
                            (var_ref.binding.order > var_ref.srcloc.order)):
                            var_ref.binding.order = var_ref.srcloc.order

    for child in context_node.children:
        AssignVarOrderByFile(child, search_instance_commands)




def AssignVarOrderByCommand(command_list, prefix_filter):
    """
    For each category in context_node, and each variable in that category,
    set the order of each variable according to the position of the
    write(), write_once(), or other command that created it.
    Only variables with the correct prefix ('$' or '@') are affected.

    """
    count = 0
    for command in command_list:
        if isinstance(command, WriteFileCommand):
            tmpl_list = command.tmpl_list
            for var_ref in tmpl_list:
                if isinstance(var_ref, VarRef):
                    if var_ref.prefix in prefix_filter:
                        count += 1
                        if ((var_ref.binding.order is None) or 
                            (var_ref.binding.order > count)):
                            var_ref.binding.order = count






def AutoAssignVals(cat_node, 
                   sort_variables,
                   reserved_values = None,
                   ignore_prior_values = False):
    """
    This function automatically assigns all the variables 
    belonging to all the categories in cat_node.categories.
    Each category has its own internal counter.  For every variable in that
    category, query the counter (which usually returns an integer), 
    and assign the variable to it.  Exceptions can be made if the integer 
    is reserved by some other variable, or if it has been already assigned.
    Afterwards, we recursively search the child nodes recursively 
    (in a depth-first-search order).

    sort_variables: Sorting the variables according to their "binding.order" 
                    counters is optional.

    """    

    if (not hasattr(cat_node, 'categories')):
        # (sometimes leaf nodes lack a 'categories' member, to save memory)
        return

    # Search the tree in a depth-first-search manner.
    # For each node, examine the "categories" associated with that node
    # (ie the list of variables whose counters lie within that node's scope).
    for cat_name, cat in cat_node.categories.items():

        # Loop through all the variables in this category.

        if sort_variables:

            # Sort the list of variables according to var_binding.order

            # First, print a progress indicator (this could be slow)
            prefix = '$'
            # Is this parent_node an StaticObj? (..or inherit from StaticObj?)
            if isinstance(cat_node, StaticObj):
                prefix = '@'
            sys.stderr.write('  sorting variables in category: '+prefix+
                             CanonicalCatName(cat_name, cat_node)+':\n')
            
            var_bind_iter = iter(sorted(cat.bindings.items(), 
                                        key=operator.itemgetter(1)))
        else:
            # Just iterate through them in the order that they were added
            # to the category list.  (This happens to be the same order as 
            # we found it earlier when searching the tree.)
            var_bind_iter = iter(cat.bindings.items())


        for leaf_node,var_binding in var_bind_iter:

            if ((var_binding.value is None) or ignore_prior_values):

                if var_binding.nptr.leaf_node.name[:9] == '__query__':
                    #   -- THE "COUNT" HACK --
                    # '__query__...' variables are not really variables.
                    # They are a mechanism to allow the user to query the
                    # category counter without incrementing it.
                    var_binding.value = str(cat.counter.query())

                elif HasWildCard(var_binding.full_name):
                    #   -- The wildcard hack ---
                    # Variables containing * or ? characters in their names
                    # are not allowed.  These are not variables, but patterns
                    # to match with other variables.  Represent them by the
                    # (full-path-expanded) string containing the * or ?.
                    var_binding.value = var_binding.full_name

                else:

                    if (not var_binding.nptr.leaf_node.IsDeleted()):

                        # For each (regular) variable, query this category's counter
                        # (convert it to a string), and see if it is already in use
                        # (in this category). If not, then set this variable's value
                        # to the counter's value. Either way, increment the counter.
                        while True:
                            cat.counter.incr()
                            value = str(cat.counter.query())
                            if ((reserved_values is None) or 
                                ((cat, value) not in reserved_values)):
                                break

                        var_binding.value = value

    # Recursively invoke AssignVarValues() on all child nodes
    for child in cat_node.children:
        AutoAssignVals(child,
                       sort_variables,
                       reserved_values,
                       ignore_prior_values)








def InvokeAllDeletes(node, 
                     null_list_warning=False,
                     null_list_error=True):
    """ This is a short simple function to carry out any remaining delete
    statements which have not yet been processed. """

    i_child = 0
    i = 0
    while i < len(node.commands):
        command = node.commands[i]
        if isinstance(command, ModCommand):
            mod_command = command
            if isinstance(mod_command.command, DeleteCommand):
                instobj_list = node.LookupMultiDescrStr(mod_command.multi_descr_str,
                                                        mod_command.command.srcloc,
                                                        null_list_warning,
                                                        null_list_error)
                for instobj in instobj_list:
                    instobj.DeleteSelf()
                del node.commands[i]
            else:
                assert(False)
        elif isinstance(command, DeleteCommand):
            node.DeleteSelf()
            del node.commands[i]
        elif isinstance(command, InstantiateCommand):
            child = node.children[i_child]
            InvokeAllDeletes(child, 
                             null_list_warning,
                             null_list_error)
            i_child += 1
            i += 1
        else:
            i += 1



# Did the user ask us to reformat the output string?
# This information is encoded in the variable's suffix.
def ExtractFormattingCommands(suffix):
    if (len(suffix) <= 1):
        return None, None
    if suffix[-1] == '}': # Get rid of any trailing '}' characters
        suffix = suffix[:-1]
    if suffix[-1] != ')': # Format functions are always followed by parens
        return None, None
    else:
        idot = suffix.find('.') # Format functions usually preceeded by '.'
        ioparen = suffix.find('(')
        icparen = suffix.find(')')
        format_fname = suffix[idot+1:ioparen]
        args = suffix[ioparen+1:icparen]
        args = args.split(',')
        for i in range(0, len(args)):
            args[i] = RemoveOuterQuotes(args[i].strip(), '\"\'')
        return format_fname, args




def Render(tmpl_list, substitute_vars=True):
    """ 
    This function converts a TextBlock,VarRef list into a string.
    It is invoked by WriteTemplatesValue() in order to print
    out the templates stored at each node of the tree.

    """

    out_str_list = []
    i = 0
    while i < len(tmpl_list):
        entry = tmpl_list[i]
        if isinstance(entry, VarRef):
            var_ref = entry
            var_bindings = var_ref.nptr.cat_node.categories[var_ref.nptr.cat_name].bindings
            #if var_ref.nptr.leaf_node not in var_bindings:
            #assert(var_ref.nptr.leaf_node in var_bindings)
            if var_ref.nptr.leaf_node.IsDeleted():
                raise InputError('Error near '+
                                 ErrorLeader(var_ref.srcloc.infile,
                                             var_ref.srcloc.lineno)+'\n'
                                 '   The variable you referred to does not exist:\n'
                                 '     '+var_ref.prefix+var_ref.descr_str+var_ref.suffix+'\n'
                                 '   (You probably deleted it or something it belonged to earlier.)\n')
            else:
                if substitute_vars:
                    value = var_bindings[var_ref.nptr.leaf_node].value
                    format_fname, args = ExtractFormattingCommands(var_ref.suffix)
                    if format_fname == 'ljust':
                        if len(args) == 1:
                            value = value.ljust(int(args[0]))
                        else:
                            value = value.ljust(int(args[0]), args[1])
                    elif format_fname == 'rjust':
                        if len(args) == 1:
                            value = value.rjust(int(args[0]))
                        else:
                            value = value.rjust(int(args[0]), args[1])
                    out_str_list.append(value)

                else:
                    out_str_list.append(var_ref.prefix +
                        SafelyEncodeString(var_bindings[var_ref.nptr.leaf_node].full_name[1:]) +
                        var_ref.suffix)

        else:
            assert(isinstance(entry, TextBlock))
            out_str_list.append(entry.text)
        i += 1

    return ''.join(out_str_list)






def MergeWriteCommands(command_list):
    """ Write commands are typically to the same file.
    We can improve performance by appending all of
    commands that write to the same file together before
    carrying out the write commands.

    """
    file_templates = defaultdict(list)
    for command in command_list:
        if isinstance(command, WriteFileCommand):
            if command.filename != None:
                file_templates[command.filename] += \
                    command.tmpl_list

    return file_templates



def WriteTemplatesValue(file_templates):
    """ Carry out the write() and write_once() commands (which 
    write out the contents of the templates contain inside them).

    """
    for filename, tmpl_list in file_templates.items():
        if filename == '':
            out_file = sys.stdout
        else:
            out_file = open(filename, 'a')

        out_file.write(Render(tmpl_list, substitute_vars=True))
        if filename != '':
            out_file.close()

    # Alternate (old method):
    #for command in command_list:
    #    if isinstance(command, WriteFileCommand):
    #        if command.filename != None:
    #            if command.filename == '':
    #                out_file = sys.stdout
    #            else:
    #                out_file = open(command.filename, 'a')
    #
    #            out_file.write(Render(command.tmpl_list))
    #
    #           if command.filename != '':
    #               out_file.close()



def WriteTemplatesVarName(file_templates):
    """ Carry out the write() and write_once() commands (which 
    write out the contents of the templates contain inside them).
    However variables within the templates are represented by their 
    full name instead of their assigned value.

    """

    for filename, tmpl_list in file_templates.items():
        if filename != '':
            out_file = open(filename + '.template', 'a')
            out_file.write(Render(tmpl_list, substitute_vars=False))
            out_file.close()



def EraseTemplateFiles(command_list):
    filenames = set([])
    for command in command_list:
        if isinstance(command, WriteFileCommand):
            if (command.filename != None) and (command.filename != ''):
                if command.filename not in filenames:
                    filenames.add(command.filename)
                    # Openning the files (in mode 'w') and closing them again
                    # erases their contents.
                    out_file = open(command.filename, 'w')
                    out_file.close()
                    out_file = open(command.filename + '.template', 'w')
                    out_file.close()

#def ClearTemplates(file_templates):
#    for filename in file_templates:
#        if filename != '':
#            out_file = open(filename, 'w')
#            out_file.close()
#            out_file = open(filename + '.template', 'w')
#            out_file.close()



def WriteVarBindingsFile(node):
    """ Write out a single file which contains a list of all 
    of the variables defined (regardless of which class they 
    were defined in).  Next to each variable name is the corresponding
    information stored in that variable (a number) that variable. 

    """
    if (not hasattr(node, 'categories')):
        # (sometimes leaf nodes lack a 'categories' member, to save memory)
        return

    out = open('ttree_assignments.txt', 'a')
    for cat_name in node.categories:
        var_bindings = node.categories[cat_name].bindings
        for nd, var_binding in var_bindings.items():
            if nd.IsDeleted():
                continue   # In that case, skip this variable

            #if type(node) is type(nd):
            if ((isinstance(node, InstanceObjBasic) and isinstance(nd, InstanceObjBasic))
                or
                (isinstance(node, StaticObj) and isinstance(nd, StaticObj))):

                # Now omit variables whos names contain "*" or "?"
                # (these are actually not variables, but wildcard patterns)
                if not HasWildCard(var_binding.full_name):
                    if len(var_binding.refs) > 0:
                        usage_example = '       #'+\
                            ErrorLeader(var_binding.refs[0].srcloc.infile, \
                                        var_binding.refs[0].srcloc.lineno)
                    else:
                        usage_example = ''
                    out.write(SafelyEncodeString(var_binding.full_name) +'   '+
                              SafelyEncodeString(var_binding.value)
                              +usage_example+'\n')
    out.close()
    for child in node.children:
        WriteVarBindingsFile(child)



def CustomizeBindings(bindings,
                      g_objectdefs,
                      g_objects):

    var_assignments = set()

    for name,vlpair in bindings.items():

        prefix = name[0]
        var_descr_str = name[1:]

        value   = vlpair.val
        dbg_loc = vlpair.loc

        if prefix == '@':
            var_binding = LookupVar(var_descr_str,
                                    g_objectdefs, 
                                    dbg_loc)
        elif prefix == '$':
            var_binding = LookupVar(var_descr_str, 
                                    g_objects,
                                    dbg_loc)
        else:
            # If the user neglected a prefix, this should have generated
            # an error earlier on.
                assert(False) 

        # Change the assignment:
        var_binding.value = value
        var_assignments.add((var_binding.category, value))

        #sys.stderr.write('  CustomizeBindings: descr=' + var_descr_str +
        #                 ', value=' + value + '\n')

    return var_assignments




##############################################################
#####################  BasicUI functions #####################
# These functions are examples of how to use the StaticObj 
# and InstanceObj data structures above, and to read a ttree file.
# These are examples only.  New programs based on ttree_lib.py 
# will probably require their own settings and functions.
##############################################################


def BasicUIReadBindingsFile(bindings_so_far, filename):
    try:
        f = open(filename, 'r')
    except IOError: 
        sys.stderr.write('Error('+g_filename+'):\n''       : unable to open file\n'
                         '\n'
                         '       \"'+bindings_filename+'\"\n'
                         '       for reading.\n'
                         '\n'
                         '       (If you were not trying to open a file with this name, then this could\n'
                         '        occur if you forgot to enclose your command-line-argument in quotes,\n'
                         '        For example, use: \'$atom:wat[2]/H1 20\' or "\$atom:wat[2]/H1 to 20"\n'
                         '        to set the variable $atom:wat[2]/H1 to 20.)\n')
        sys.exit(1)

    BasicUIReadBindingsStream(bindings_so_far, f, filename)
    f.close()




def BasicUIReadBindingsText(bindings_so_far, text, source_name=''):
    if sys.version > '3':
        in_stream = io.StringIO(text)
    else:
        in_stream = cStringIO.StringIO(text)
    return BasicUIReadBindingsStream(bindings_so_far, in_stream, source_name)



class ValLocPair(object):
    __slots__=["val","loc"]
    def __init__(self,
                 val = None,
                 loc = None):
        self.val = val
        self.loc = loc


def BasicUIReadBindingsStream(bindings_so_far, in_stream, source_name=''):

    # EXAMPLE (simple version)
    # The simple version of this function commented out below
    # does not handle variable whose names or values
    # contain strange or escaped characters, quotes or whitespace.
    # But I kept it in for illustrative purposes:
    #
    #for line in f:
    #    line = line.strip()
    #    tokens = line.split()
    #    if len(tokens) == 2:
    #        var_name  = tokens[0]
    #        var_value = tokens[1]
    #        var_assignments[var_name] = var_value
    #f.close()


    lex = TemplateLexer(in_stream, source_name)
    tmpllist = lex.ReadTemplate()
    i = 0
    if isinstance(tmpllist[0], TextBlock):
        i += 1
    while i+1 < len(tmpllist):
        # process one line at a time (2 entries per line)
        var_ref = tmpllist[i]
        text_block = tmpllist[i+1]
        assert(isinstance(var_ref, VarRef))
        if (not isinstance(text_block, TextBlock)):
            raise InputError('Error('+g_filename+'):\n'
                             '       This is not a valid name-value pair:\n'
                             '          \"'+var_ref.prefix+var_ref.descr_str+' '+text_block.text.rstrip()+'\"\n'
                             '       Each variable asignment should contain a variable name (beginning with\n'
                             '       @ or $) followed by a space, and then a string you want to assign to it.\n'
                             '       (Surrounding quotes are optional and will be removed.)\n')

        # Variables in the ttree_assignments.txt file use "full-path" style.
        # In other words, the full name of the variable, (including all
        # path information) is stored var_ref.descr_str,
        # and the first character of the prefix stores either a @ or $
        var_name  = var_ref.prefix[:1] + var_ref.descr_str
        text = SplitQuotedString(text_block.text.strip())
        var_value = EscCharStrToChar(RemoveOuterQuotes(text, '\'\"'))
        bindings_so_far[var_name] = ValLocPair(var_value, lex.GetSrcLoc())
        i += 2




class BasicUISettings(object):
    """
    BasicUISettings() contains several run-time user customisations 
    for ttree. (These effect the order and values assigned to variables 
    in a ttreee file). 
    This object, along with the other "UI" functions below are examples only.
    (New programs based on ttree_lib.py will probably have their own settings
     and functions.)

    Members:
        user_bindings 
        user_bindings_x
    These are lists containing pairs of variable names, 
    and the string values they are bound to (which are typically numeric).
    Values specified in the "user_bindings_x" list are "exclusive".
    This means their values are reserved, so that later on, when other 
    variables (in the same category) are automatically assigned to values, care
    care will be taken to avoid duplicating the values in user_bindings_x.
    However variables in the "user_bindings" list are assigned without regard
    to the values assigned to other variables, and may or may not be unique.

        order_method
    The order_method specifies the order in which values will be automatically 
    assigned to variables.  (In the context of building molecular simulation
    input files, this helps the user to insure that the order of the atoms
    created by the ttree file matches the order they appear in other files
    created by other programs.)

    """
    def __init__(self,
                 user_bindings_x=None,
                 user_bindings=None,
                 order_method='by_command',
                 lex=None):
        if user_bindings_x:
            self.user_bindings_x = user_bindings_x
        else:
            self.user_bindings_x = OrderedDict()
        if user_bindings:
            self.user_bindings = user_bindings
        else:
            self.user_bindings = OrderedDict()
        self.order_method = order_method
        self.lex = lex




def BasicUIParseArgs(argv, settings):
    """
    BasicUIParseArgs()
    The following function contains part of the user interface for a 
    typical ttree-based program.  This function processes an argument list
    and extracts the common ttree user settings.
    This function, along with the other "UI" functions below are examples only.
    (New programs based on ttree_lib.py will probably have their own UI.)
    
    """

    #argv = [arg for arg in orig_argv] # (make a deep copy of "orig_argv")

    # This error message is used in multiple places:

    bind_err_msg = 'should either be followed by a 2-column\n'+\
        '      file (containing variable-value pairs on each line).\n'+\
        '   --OR-- a quoted string (such as \"@atom:x 2\")\n'+\
        '      with the full variable name and its desired value.'
    bind_err_msg_var = 'Missing value, or space needed separating variable\n'+\
        '      and value. (Remember to use quotes to surround the argument\n'+\
        '      containing the variable name, and it\'s assigned value.)'

    i = 1


    while i < len(argv):
        #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')
        if argv[i] == '-a':
            if ((i+1 >= len(argv)) or (argv[i+1][:1] == '-')):
                raise InputError('Error('+g_filename+'):\n'
                                 '      Error in -a \"'+argv[i+1]+' argument.\"\n'
                                 '      The -a flag '+bind_err_msg)
            if (argv[i+1][0] in '@$'):
                #tokens = argv[i+1].strip().split(' ')
                tokens = SplitQuotedString(argv[i+1].strip())
                if len(tokens) < 2:
                    raise InputError('Error('+g_filename+'):\n'
                                     '      Error in -a \"'+argv[i+1]+'\" argument.\n'
                                     '      '+bind_err_msg_var)
                BasicUIReadBindingsText(settings.user_bindings_x,
                                        argv[i+1], 
                                        '__command_line_argument__')
            else:
                BasicUIReadBindingsFile(settings.user_bindings_x,
                                        argv[i+1])
            #i += 2
            del(argv[i:i+2])
        elif argv[i] == '-b':
            if ((i+1 >= len(argv)) or (argv[i+1][:1] == '-')):
                raise InputError('Error('+g_filename+'):\n'
                                 '      Error in -b \"'+argv[i+1]+' argument.\"\n'
                                 '      The -b flag '+bind_err_msg)
            if (argv[i+1][0] in '@$'):
                #tokens = argv[i+1].strip().split(' ')
                tokens = SplitQuotedString(argv[i+1].strip())
                if len(tokens) < 2:
                    raise InputError('Error('+g_filename+'):\n'
                                     '      Error in -b \"'+argv[i+1]+'\" argument.\n'
                                     '      '+bind_err_msg_var)
                BasicUIReadBindingsText(settings.user_bindings,
                                        argv[i+1], 
                                        '__command_line_argument__')
            else:
                BasicUIReadBindingsFile(settings.user_bindings,
                                        argv[i+1])
            #i += 2
            del(argv[i:i+2])
        elif argv[i] == '-order-command':
            settings.order_method = 'by_command'
            #i += 1
            del(argv[i:i+1])
        elif argv[i] == '-order-file':
            settings.order_method = 'by_file'
            #i += 1
            del(argv[i:i+1])
        elif ((argv[i] == '-order-tree') or (argv[i] == '-order-dfs')):
            settings.order_method = 'by_tree'
            del(argv[i:i+1])
        elif ((argv[i] == '-importpath') or
              (argv[i] == '-import-path') or
              (argv[i] == '-import_path')):
            if ((i+1 >= len(argv)) or (argv[i+1][:1] == '-')):
                raise InputError('Error('+g_filename+'):\n'
                                 '     Error in \"'+argv[i]+'\" argument.\"\n'
                                 '     The \"'+argv[i]+'\" argument should be followed by the name of\n'
                                 '     an environment variable storing a path for including/importing files.\n')
            TtreeShlex.custom_path = RemoveOuterQuotes(argv[i+1])
            del(argv[i:i+2])

        elif ((argv[i][0] == '-') and (__name__ == "__main__")):
            #elif (__name__ == "__main__"):
            raise InputError('Error('+g_filename+'):\n'
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
            raise InputError('Error('+g_filename+'):\n'
                             '       This program requires at least one argument\n'
                             '       the name of a file containing ttree template commands\n')
        elif len(argv) == 2:
            try:
                settings.lex = TemplateLexer(open(argv[1], 'r'), argv[1]) # Parse text from file
            except IOError: 
                sys.stderr.write('Error('+g_filename+'):\n'
                                 '       unable to open file\n'
                                 '       \"'+argv[1]+'\"\n'
                                 '       for reading.\n')
                sys.exit(1)
            del(argv[1:2])

        else:
            # if there are more than 2 remaining arguments,
            problem_args = ['\"'+arg+'\"' for arg in argv[1:]]
            raise InputError('Syntax Error ('+g_filename+'):\n'
                             '       Problem with argument list.\n'
                             '       The remaining arguments are:\n\n'
                             '         '+(' '.join(problem_args))+'\n\n'
                             '       (The actual problem may be earlier in the argument list.\n'
                             '       If these arguments are source files, then keep in mind\n'
                             '       that this program can not parse multiple source files.)\n'
                             '       Check the syntax of the entire argument list.\n')






def BasicUI(settings,
            static_tree_root, 
            instance_tree_root,
            static_commands,
            instance_commands):
    """
    BasicUI()
    This function loads a ttree file and optional custom bindings for it,
    creates a "static" tree (of defined ttree classes),
    creates an "instance" tree (of instantiated ttree objects),
    automatically assigns values to unbound variables,
    substitutes them into text templates (renders the template).
    The actual writing of the templates to a file is not handled here.

    """

    # Parsing, and compiling is a multi-pass process.

    # Step 1: Read in the StaticObj (class) defintions, without checking
    # whether or not the instance_children refer to valid StaticObj types.

    sys.stderr.write('parsing the class definitions...')
    static_tree_root.Parse(settings.lex)
    #gc.collect()

    #sys.stderr.write('static = ' + str(static_tree_root) + '\n')

    # Step 2: Now that the static tree has been constructed, lookup 
    #         any references to classes (StaticObjs), contained within
    #         the instance_children or class_parents of each node in
    #         static_tree_root.  Replace them with (pointers to)
    #         the StaticObjs they refer to (and check validity).
    #         (Note: Variables stored within the templates defined by write() 
    #                and write_once() statements may also refer to StaticObjs in
    #                the tree, but we leave these references alone.  We handle
    #                these assignments later using "AssignVarPtrs()" below.)
    sys.stderr.write(' done\nlooking up classes...')
    static_tree_root.LookupStaticRefs()
    #gc.collect()

    # Step 3: Now scan through all the (static) variables within the templates
    #         and replace the (static) variable references to pointers
    #         to nodes in the StaticObj tree:
    sys.stderr.write(' done\nlooking up @variables...')
    # Here we assign pointers for variables in "write_once(){text}" templates:
    AssignVarPtrs(static_tree_root, search_instance_commands=False)
    # Here we assign pointers for variables in "write(){text}" templates:
    AssignVarPtrs(static_tree_root, search_instance_commands=True)
    sys.stderr.write(' done\nconstructing the tree of class definitions...')
    sys.stderr.write(' done\n\nclass_def_tree = ' + str(static_tree_root) + '\n\n')
    #gc.collect()

    # Step 4: Construct the instance tree (the tree of instantiated
    #         classes) from the static tree of type definitions.
    sys.stderr.write('constructing the instance tree...\n')
    class_parents_in_use = set([])
    instance_tree_root.BuildInstanceTree(static_tree_root, class_parents_in_use)
    #sys.stderr.write('done\n  garbage collection...')
    #gc.collect()
    sys.stderr.write(' done\n')
    #sys.stderr.write('instance_tree = ' + str(instance_tree_root) + '\n')

    # Step 5: Now scan through all the (instance) variables within the templates
    #         and replace the (instance) variable references to pointers
    #         to nodes in the InstanceObj tree:
    sys.stderr.write(' done\nlooking up $variables...')
    AssignVarPtrs(instance_tree_root, search_instance_commands=False)
    #sys.stderr.write('done\n  garbage collection...')
    #gc.collect()

    # Step 6: Now carry out all of the "delete" commands (deferred earlier).
    # (These were deferred because the instance tree must be complete before any
    #  references to target nodes (with non-trivial paths) can be understood.)
    InvokeAllDeletes(instance_tree_root, 
                     null_list_warning=False,
                     null_list_error=True)
    sys.stderr.write(' done\n')
    #sys.stderr.write('instance_v = ' + str(instance_tree_root) + '\n')
    #gc.collect()

    # Step 7: The commands must be carried out in a specific order.
    #         (for example, the "write()" and "new" commands).
    #         Search through the tree, and append commands to a command list.
    #         Then re-order the list in the order the commands should have
    #         been executed in.  (We don't carry out the commands yet, 
    #         we just store them and sort them.)
    class_parents_in_use = set([])
    static_tree_root.BuildCommandList(static_commands)
    instance_tree_root.BuildCommandList(instance_commands)
    #sys.stderr.write('static_commands = '+str(static_commands)+'\n')
    #sys.stderr.write('instance_commands = '+str(instance_commands)+'\n')

    # Step 8: We are about to assign numbers to the variables.
    #         We need to decide the order in which to assign them.
    #         By default static variables (@) are assigned in the order 
    #         they appear in the file.
    #         And, by default instance variables ($) 
    #         are assigned in the order they are created during instantiation.
    #sys.stderr.write(' done\ndetermining variable count order...')
    AssignVarOrderByFile(static_tree_root, search_instance_commands=False)
    AssignVarOrderByFile(static_tree_root, search_instance_commands=True)
    AssignVarOrderByCommand(instance_commands, '$')

    # Step 9: Assign the variables.
    #         (If the user requested any customized variable bindings,
    #          load those now.)
    if len(settings.user_bindings_x) > 0:
        reserved_values = CustomizeBindings(settings.user_bindings_x,
                                            static_tree_root,
                                            instance_tree_root)
    else:
        reserved_values = None

    sys.stderr.write('sorting variables...\n')
    AutoAssignVals(static_tree_root,
                   (settings.order_method != 'by_tree'),
                   reserved_values)

    AutoAssignVals(instance_tree_root,
                   (settings.order_method != 'by_tree'),
                   reserved_values)
                        
    if len(settings.user_bindings) > 0:
        CustomizeBindings(settings.user_bindings,
                          static_tree_root,
                          instance_tree_root)

    sys.stderr.write(' done\n')




if __name__ == "__main__":

    """
    This is is a "main module" wrapper for invoking ttree.py
    as a stand alone program.  This program:

    1)reads a ttree file, 
    2)constructs a tree of class definitions (g_objectdefs)
    3)constructs a tree of instantiated class objects (g_objects),
    4)automatically assigns values to the variables,
    5)and carries out the "write" commands to write the templates a file(s).

    """

    #######  Main Code Below: #######
    g_program_name = g_filename
    sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+' ')
    sys.stderr.write('\n(python version '+str(sys.version)+')\n')

    try:

        settings = BasicUISettings()
        BasicUIParseArgs(sys.argv, settings)


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

        # Now write the files
        # (Finally carry out the "write()" and "write_once()" commands.)

        # Optional: Multiple commands to write to the same file can be merged to
        #           reduce the number of times the file is openned and closed.
        sys.stderr.write('writing templates...\n')
        # Erase the files that will be written to:
        EraseTemplateFiles(g_static_commands)
        EraseTemplateFiles(g_instance_commands)

        g_static_commands   = MergeWriteCommands(g_static_commands)
        g_instance_commands = MergeWriteCommands(g_instance_commands)

        # Write the files with the original variable names present
        WriteTemplatesVarName(g_static_commands)
        WriteTemplatesVarName(g_instance_commands)
        # Write the files with the variable names substituted by values
        WriteTemplatesValue(g_static_commands)
        WriteTemplatesValue(g_instance_commands)
    
        sys.stderr.write(' done\n')

        # Step 11: Now write the variable bindings/assignments table.
        sys.stderr.write('writing \"ttree_assignments.txt\" file...')
        open('ttree_assignments.txt', 'w').close() # <-- erase previous version.
        WriteVarBindingsFile(g_objectdefs)
        WriteVarBindingsFile(g_objects)

    except (ValueError, InputError) as err:
        sys.stderr.write('\n\n'+str(err)+'\n')
        sys.exit(-1)

