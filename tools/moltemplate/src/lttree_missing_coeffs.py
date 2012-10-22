if __name__ == "__main__":

    try:

        # Data structures to store the class definitionss and instances
        static_tree_root = StaticObj('', None) # The root of the static tree
                                               # has name ''

        # Parsing, and compiling is a multi-pass process.

        # Step 1: Read in the StaticObj (class) defintions, without checking
        # whether or not the instance_children refer to valid StaticObj types.

        #sys.stderr.write('parsing the class definitions...')
        static_tree_root.Parse(settings.lex)

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
        #sys.stderr.write(' done\nlooking up classes...')
        static_tree_root.LookupStaticRefs()

        # Step 3: Now scan through all the (static) variables within the templates
        #         and replace the (static) variable references to pointers
        #         to nodes in the StaticObj tree:
        #sys.stderr.write(' done\nlooking up @variables...')
        AssignVarPtrs(static_tree_root, search_instance_commands=False)
        AssignVarPtrs(static_tree_root, search_instance_commands=True)
        #sys.stderr.write(' done\nconstructing the tree of class definitions...')
        sys.stderr.write(' done\n\nclass_def_tree = '+str(static_tree_root)+'\n\n')
