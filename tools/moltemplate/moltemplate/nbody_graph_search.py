# Author: Andrew Jewett (jewett.aij at g mail)
#         http://www.chem.ucsb.edu/~sheagroup
# License: 3-clause BSD License  (See LICENSE.TXT)
# Copyright (c) 2012, Regents of the University of California
# All rights reserved.


#__all__ = ['Ugraph', 'GraphMatcher', 'DFS', 'GenError', 'GraphError', 'Disconnected', 'NotUndirected']


import sys
import copy
from operator import itemgetter


class GenError(Exception):
    """
    An exception class containing string for error reporting.

    """

    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return self.err_msg

    def __repr__(self):
        return str(self)


class GraphError(GenError):
    """
    An exception class containing a graph and a string for error reporting.

    """

    def __init__(self, g, err_msg):
        GenError.__init__(self, err_msg)
        self.g = g

    def __str__(self):
        g_str = str(g)
        # If the string representation of the graph is too
        # large to fit in one screen, truncate it
        g_str_lines = g_str.split('\n')
        if (len(g_str_lines) > 12):
            g_str_lines = g_str_lines[0:12] + \
                [' ...(additional lines not shown)]']
            g_str = '\n'.join(g_str_lines)
        return 'Problem with graph:\n' + g_str + '\n' + self.err_msg

    def __repr__(self):
        return str(self)


class Disconnected(GraphError):

    def __init__(self, g, err_msg):
        GraphError.__init__(self, g, err_msg)


class NotUndirected(GraphError):

    def __init__(self, g, err_msg):
        GraphError.__init__(self, g, err_msg)


class Edge(object):
    __slots__ = ["start", "stop", "attr"]

    def __init__(self,
                 iv_start,  # edge starts here (index into vertex list)
                 iv_stop,   # edge ends here (index into vertex list)
                 attr=None):  # edges have an optional type attribute
        self.start = iv_start
        self.stop = iv_stop
        self.attr = attr

    def __str__(self):
        return '(' + str(self.start) + ',' + str(self.stop) + ')'

    def __repr__(self):
        return str(self)


class Vertex(object):
    __slots__ = ["attr"]

    def __init__(self, attr=None):
        self.attr = attr


class Dgraph(object):
    """
    This class is a minimal implementation of a directed graph.
    Vertices and edges are accessed by integer index only (beginning at 0).
    Multiple edges connecting the same pair of vertices are allowed.
    (One would use the AddEdge() member function to accomplish this.)
    Both vertices and edges have an optional "attr" attribute.

    """

    NULL = -1  # forbidden vertex id number (used several places)

    def __init__(self, edgelist=None):
        """
        The constructor accepts an optional neighborlist argument.
        This is a simple list of neighbors for every vertex in the graph
        and it completely defines the topology of the graph.
        (Vertex and edge attributes can be specified later.)

        Alternatley, you can leave the neighborlist argument blank,
        and build the graph one vertex at a time later
        using the "AddVertex()" and "AddEdge()" commands.
        (AddEdge() commands must be issued strictly after
         all vertices have been defined.)

        """

        if edgelist == None:
            self.verts = []
            self.edges = []
            # integer keeps track of # of vertices = len(self.verts)
            self.nv = 0
            # integer keeps track of # of edges    = len(self.edges)
            self.ne = 0
            self.neighbors = []  # The adjacency list.

        else:
            # Parse the edge-list format:
            iv_max = 0  # <-- what's the vertex with the maximum id number?
            for i in range(0, len(edgelist)):
                iv = edgelist[i][0]
                jv = edgelist[i][1]
                if ((iv < 0) or (jv < 0)):
                    raise(GenError(
                        'Error in Dgraph.__init__: Negative vertex number pair encountered: (' + str(iv) + ',' + str(jv) + ')'))
                if iv > iv_max:
                    iv_max = iv
                if jv > iv_max:
                    iv_max = jv

            self.nv = iv_max + 1
            self.verts = [Vertex() for iv in range(0, self.nv)]
            self.edges = []
            self.ne = 0
            self.neighbors = [[] for iv in range(0, self.nv)]

            for i in range(0, len(edgelist)):
                iv = edgelist[i][0]
                jv = edgelist[i][1]
                self.neighbors[iv].append(self.ne)
                self.edges.append(Edge(iv, jv))
                self.ne += 1
            assert(self.ne == len(self.edges))

            self.SortNeighborLists()

    def AddVertex(self, iv=-1, attr=None):
        """
        Add a vertex to the graph.
        (Edges connected to this vertex must be added later using "AddEdge()"
         All vertices should be added before "AddEdge()" is ever invoked.)

        Optional "attr" argument allows you to set the attribute of this vertex.
        (for example, in a molecule this might correspond to the type of atom
         in the molecule).

        Optional "iv" argument allows you to specify the index of that vertex.
        Vertices can be added in any order, but thei vertex id numbers
        should eventually fill the range from 0 to self.nv-1.

        """
        if iv == -1:  # if iv unspecified, put the vertex at the end of the list
            iv = self.nv

        if iv < self.nv:
            self.verts[iv].attr = attr
        else:
            # In case there is a gap between iv and nv, fill it with blanks
            self.verts += ([Vertex()] * ((1 + iv) - self.nv))
            self.neighbors += ([[]] * ((1 + iv) - self.nv))
            self.verts[iv].attr = attr
            self.nv = iv + 1
            assert(self.nv == len(self.verts))
            assert(self.nv == len(self.neighbors))

    def AddEdge(self, iv, jv, attr=None, remove_duplicates=False):
        """
        Add an edge to graph connecting vertex iv to jv.
        (both are integers from 0 to self.nv-1)
        This function must not be called until all vertices have been added.
        If the edge is already present (and remove_duplicates==True),
        no new edge will be added.

        """
        if remove_duplicates:
            for je in self.neighbors[iv]:
                if jv == self.edges[je].stop:
                    return  # In that case, do nothing, the edge is already present
        self.edges.append(Edge(iv, jv, attr))
        self.neighbors[iv].append(self.ne)
        self.ne += 1
        assert(self.ne == len(self.edges))

    def ReorderVerts(self, vpermutation, invert=False):
        """
        This function allows the user to re-order (relabel) the vertices
        in a graph, making the necessary changes to the
        self.verts, self.edges, and self.neighbors lists.
        By default (invert=False).  The vpermutation is a list
        from 1 to self.nv which is interpreted this way:
            iv = vpermutation[iv_orig]
        where "iv" and "iv_orig" are the vertex id numbers before
        and after the mapping (which also corresponds to its
        position in the self.verts and self.neighbors arrays).

        """
        assert(len(self.verts) == self.nv)
        assert(len(self.edges) == self.ne)
        assert(len(vpermutation) == self.nv)

        if (invert):
            vperm = [-1 for iv in vpermutation]
            for iv in range(0, self.nv):
                vperm[vpermutation[iv]] = iv
        else:
            vperm = vpermutation

        orig_verts = [vert for vert in self.verts]
        for iv_old in range(0, self.nv):
            iv = vperm[iv_old]
            self.verts[iv] = orig_verts[iv_old]

        for ie in range(0, self.ne):
            self.edges[ie].start = vperm[self.edges[ie].start]
            self.edges[ie].stop = vperm[self.edges[ie].stop]

        orig_neighbors = [nlist for nlist in self.neighbors]
        # self.neighbors is a 2-d array.
        # We need to re-sort "self.neighbors" because the first index is
        # a vertex id number, and these id numbers have been permuted.
        # However, there's no need to sort the contents of each sub-array
        # (self.neighbors[iv]), because these are edge id numbers (indices into
        # the self.edges[] array). These edge index numbers are never altered.
        # (However the entries stored in self.edges were modified earlier.)
        for iv_old in range(0, self.nv):
            iv = vperm[iv_old]
            self.neighbors[iv] = orig_neighbors[iv_old]

        # Optional:
        self.SortNeighborLists()

    def ReorderEdges(self, epermutation, invert=False):
        """
        This function allows the user to re-order (relabel) the
        edges in a graph, making the necessary changes to the
        self.edges and self.neighbors lists.
        By default (invert=False).  The epermutation is a list
        from 1 to self.ne which is interpreted this way:
            ie = epermutation[ie_orig]
        where "ie" and "ie_orig" are the edge id numbers before
        and after the mapping (which also corresponds to that edge's
        position in the self.edges array).
            (Minor detail: Recall that in this code, Ugraphs
        are implemented by placing two (directed) edges between each pair of
        connected, adjacent vertices, which point back-and-forth between them.
        Consequently the list of edges in self.edges is often typically
        twice as large you might expect.)

        """
        assert(len(self.verts) == self.nv)
        assert(len(self.edges) == self.ne)
        assert(len(epermutation) == self.ne)

        if (invert):
            eperm = [-1 for ie in epermutation]
            for ie in range(0, self.ne):
                eperm[epermutation[ie]] = ie
        else:
            eperm = epermutation

        orig_edges = [edge for edge in self.edges]
        for ie_old in range(0, self.ne):
            ie = eperm[ie_old]
            self.edges[ie] = orig_edges[ie_old]

        for iv in range(0, self.nv):
            for j in range(0, len(self.neighbors[iv])):
                je_old = self.neighbors[iv][j]
                self.neighbors[iv][j] = eperm[je_old]

    def SortNeighborLists(self):
        assert(self.nv == len(self.neighbors))
        for iv in range(0, self.nv):
            # Back when self.neighbors was just a 2-dimensional list of
            # vertex id numbers, then the following line would have worked:
            #  self.neighbors[iv].sort()

            # ugly python code alert:
            # Unfortunately, we had to change the format of self.neighbors. Now
            # it is a list of indices into the self.edges array ("ie" numbers).
            # We want to sort the "ie" numbers by the vertices they point to.
            # self.edge[ie].start should point to the current vertex (hopefully).
            # self.edge[ie].stop should point to the vertex it's attached to.
            # So we want to sort the ie's in self.neighbors by self.edge[ie].stop
            # Create a temporary array of 2-tuples (ie, jv)
            nlist = [(ie, self.edges[ie].stop)
                     for ie in self.neighbors[iv]]
            self.neighbors[iv] = [ie for ie, jv in sorted(nlist,
                                                          key=itemgetter(1))]

    def FindEdge(self, istart, istop):
        """
        A simple function looks up the edge id number
        corresponding to an edge connecting vertex istart to istop.
        If not present returns Dgraph.NULL.

        """
        iv = istart
        for je in self.neighbors[iv]:
            jv = self.edges[je].stop
            if jv == istop:
                return je
        return Dgraph.NULL

    def GetVert(self, iv):
        return self.verts[iv]

    def GetEdge(self, ie):
        return self.edges[ie]

    def GetNumVerts(self):
        return self.nv

    def GetNumEdges(self):
        return self.ne

    # Commenting out.  I think it's clearer to use python's deepcopy instead
    # def makecopy(self):
    #    new_copy = Ugraph()
    #    new_copy.verts = [vertex for vertex in self.verts]
    #    new_copy.edges = [  edge for edge   in self.edges]
    #    new_copy.neighbors = [nlist for nlist in self.neighbors]
    #    new_copy.nv = self.nv
    #    new_copy.ne = self.ne
    #    return new_copy

    def __str__(self):
        # Print the graph as a list of neighbor-lists.
        # (Note: This is the same format as the first argument to __init__().
        #        The Vertex.attr and Edge.attr attributes are not printed.)
        l = ['([']
        for iv in range(0, self.nv):
            l.append('[')
            for j in range(0, len(self.neighbors[iv])):
                je = self.neighbors[iv][j]
                jv = self.edges[je].stop
                l.append(str(jv))
                if j < len(self.neighbors[iv]) - 1:
                    l.append(', ')
                else:
                    l.append(']')
            if iv < self.nv - 1:
                l.append(',\n  ')
            else:
                l.append(']')
        l.append(',\n [')
        for ie in range(0, self.ne):
            l.append(str(self.edges[ie]))
            if ie < self.ne - 1:
                l.append(', ')
            else:
                l.append('])\n')
        return ''.join(l)

    def __repr__(self):
        return str(self)


class Ugraph(Dgraph):
    """
    This class is a minimal implementation of an undirected graph.
    Vertices and edges are accessed by integer index only (beginning at 0).
    Multiple edges connecting the same pair of vertices are allowed.
    (One would use the AddEdge() member function to accomplish this.)
    Both vertices and edges have an optional "attr" attribute.

        Undirected graphs (Ugraphs) are represented internally as
        directed graphs.  This means that for every edge in the Ugraph,
        connecting vertex 2 to 3, for example, two edges are stored
        internally, (2 -> 3,   and   3 -> 2),
        Edges which begin and end at the same vertex are stored only once.)

    """

    def __init__(self, edgelist=None):
        Dgraph.__init__(self, edgelist)

        # Now add the extra edges which point in the reverse direction.
        neu = self.ne
        ned = self.ne
        for ieu in range(0, self.ne):
            iv = self.edges[ieu].start
            jv = self.edges[ieu].stop
            if iv != jv:
                ned += 1

        self.ieu_to_ied = [Dgraph.NULL for ieu in range(0, neu)]
        self.ied_to_ieu = [Dgraph.NULL for ied in range(0, ned)]

        ied_redundant = neu
        for ie in range(0, neu):
            iv = self.edges[ie].start
            jv = self.edges[ie].stop
            attr = self.edges[ie].attr
            self.ieu_to_ied[ie] = ie
            self.ied_to_ieu[ie] = ie

            if iv != jv:
               # Then create another edge which points in the reverse direction
                # <--this increments self.ne
                Dgraph.AddEdge(self, jv, iv, attr)
                self.ied_to_ieu[ied_redundant] = ie
                ied_redundant += 1

        self.neu = neu
        assert(self.ne == ned)

    def AddEdge(self, iv, jv, attr=None, remove_duplicates=False):
        """
        Add an edge to an undirected graph connecting vertices iv and jv.
        If the edge is already present (and remove_duplicates==True),
        no new edge will be added.

        Note: Undirected Ugraphs are implemented by creating two separate
              digraph edges that conect iv->jv  and jv->iv.

        """

        self.ieu_to_ied.append(len(self.edges))
        Dgraph.AddEdge(self, iv, jv, attr, remove_duplicates)
        self.ied_to_ieu.append(self.neu)
        if jv != iv:
            Dgraph.AddEdge(self, jv, iv, attr, remove_duplicates)
            self.ied_to_ieu.append(self.neu)
        self.neu += 1

        assert(len(self.ieu_to_ied) == self.neu)
        assert(len(self.ied_to_ieu) == len(self.edges))

    def ReorderEdges(self, epermutation, invert=False):
        Dgraph.ReorderEdges(self, epermutation, invert)

        # Now update the
        # self.ieu_to_ied and
        # self.ied_to_ieu lookup tables:

        if (invert):  # (first invert the permutation if necessary)
            eperm = [-1 for ie in epermutation]
            for ie in range(0, self.ne):
                eperm[epermutation[ie]] = ie
        else:
            eperm = epermutation
        # epermutation.reverse()

        ieu_to_ied_orig = [ied for ied in self.ieu_to_ied]
        ied_to_ieu_orig = [ieu for ieu in self.ied_to_ieu]

        for ieu in range(0, self.neu):
            ied_old = ieu_to_ied_orig[ieu]
            ied = eperm[ied_old]
            self.ieu_to_ied[ieu] = ied
        for ied_old in range(0, self.ne):
            ieu = ied_to_ieu_orig[ied_old]
            ied = eperm[ied_old]
            self.ied_to_ieu[ied] = ieu

        eperm = epermutation

    def LookupDirectedEdgeIdx(self, ieu):
        return self.ieu_to_ied[ieu]

    def LookupUndirectedEdgeIdx(self, ied):
        return self.ied_to_ieu[ied]

    # def GetVert(self, iv):    <-- (inherited from parent)
    #    return self.verts[iv]

    def GetEdge(self, ieu):
        ied = self.ieu_to_ied[ieu]
        return self.edges[ied]

    # def GetNumVerts(self):    <-- (inherited from parent)
    #    return self.nv

    def GetNumEdges(self):
        return self.neu

    def FindEdge(self, istart, istop):
        """
        A simple function looks up the (undirected) edge id number
        corresponding to an edge connecting vertices istart and istop.
        If not present returns Dgraph.NULL.

        To find the corresponding entry in the self.edges[] list,
        you can either:
            use the LookupDirectedEdge() lookup function
             or
            you can use the parent-class' version of this function
            Dgraph.FindEdge(self, istart, istop) which returns
            this number by default.

        """
        ied = Dgraph.FindEdge(self, istart, istop)
        ieu = self.LookupUndirectedEdgeIdx(ied)
        return ieu

    def CalcEdgeLookupTable(self):
        """
        COMMENT: THIS NEXT FUNCTION IS PROBABLY NOT NECESSARY AND MIGHT BE
                 REMOVED AT A LATER TIME WHEN I FIGURE OUT A BETTER WAY.

        Because undirected graphs (Ugraphs) are implemented as directed graphs
        (Dgraphs) with redundant edges, they may have some extra edges which
        the user never explicitly asked for.
        There is some confusion about whether the i'th edge refers to
        the i'th undirected edge that the user explicitly added, or
        the i'th directed edge which is stored internally.

           (The number of directed edges is usually twice the number of
           edges that the user asked for.  But not always, because edges
           wich start and end at the same vertex are only represented once.)

        This function calculates lookup tables to translate between
        the two edge numbering systems:

        self.ieu_to_ied[ieu] returns a directed edge id number,
                             (which is an index into the self.edges list)
                             corresponding to the ieu'th undirected edge
                             which was explicitly added by the caller.

        self.ied_to_ieu[ied] takes a directed edge id number (ied,
                             an index into the self.edges list)
                             and returns the undirected edge number,
                             which is allways <= ied

        """

        self.ieu_to_ied = []
        self.ied_to_ieu = [Ugraph.NULL for ied in range(0, self.ne)]
        for ied in range(0, self.ne):
            iv = self.edges[ied].start
            jv = self.edges[ied].stop
            ieu = len(self.ieu_to_ied)
            self.ied_to_ieu[ied] = ieu
            if iv <= jv:
                self.ieu_to_ied.append(ied)


def SortVertsByDegree(g):
    vert_numneighbors = [(iv, len(g.neighbors[iv])) for iv in range(0, g.nv)]
    vert_numneighbors.sort(key=itemgetter(1))
    order = [vert_numneighbors[iv][0] for iv in range(0, g.nv)]
    g.ReorderVerts(order, invert=True)


class DFS(object):
    """
    This class contains a member function (Order()) calculates the order
    of vertices visited in a depth-first-search over a connected graph.

    """

    def __init__(self, g):
        self.g = g
        self.sv = 0  # integer sv keeps track of how many vertices visited so far
        self.se = 0  # integer se keeps track of how many edges visited so far
        self.vvisited = [False for iv in range(0, self.g.nv)]  # verts visited
        self.vorder = [Dgraph.NULL for iv in range(
            0, self.g.nv)]  # search order
        self.evisited = [False for ie in range(0, self.g.ne)]  # edges visited
        self.eorder = [Dgraph.NULL for ie in range(
            0, self.g.ne)]  # search order

    def Reset(self):
        self.sv = 0
        self.se = 0
        for iv in range(0, self.g.nv):
            self.vvisited[iv] = False
            self.vorder[iv] = Dgraph.NULL
        for ie in range(0, self.g.ne):
            self.evisited[ie] = False
            self.eorder[ie] = Dgraph.NULL

    def Order(self, starting_node=0):
        """
        VisitOrder(starting_node)
        generates a list of integers from 0 to self.g.nv-1 (=#vertices minus 1)
        which represents the order in which the vertices would be visited
        during a Depth-First-Search.

        The first vertex visited is specified by the "starting_node" argument
        (an integer (from 0 to g.nv-1)).

        """
        self.Reset()
        # The first vertex to be visited should be the starting_node
        self.vorder[0] = starting_node
        self.vvisited[starting_node] = True
        self.sv = 1
        self._Order(starting_node)
        if self.sv != self.g.nv:
            raise(Disconnected(self.g, "Error(Order): " +
                               "The input graph is not connected."))
        assert(self.se == self.g.ne)
        return ([iv for iv in self.vorder], [ie for ie in self.eorder])
        # return self.order

    def _Order(self, iv):
        """
        _Order() is a recursive function which carries out a
        Depth-First-Search over the graph "self.g", starting with vertex iv.

        """
        for je in self.g.neighbors[iv]:
            jv = self.g.edges[je].stop
            if not self.evisited[je]:
                self.eorder[self.se] = je
                self.se += 1
                self.evisited[je] = True
                if not self.vvisited[jv]:
                    self.vorder[self.sv] = jv
                    self.sv += 1
                    self.vvisited[jv] = True
                    self._Order(jv)

    def IsConnected(self):
        self.Reset()
        self._Order(0)
        return (self.sv == self.g.nv)

    def IsCyclic(self):
        """
        IsCyclic() returns True if the graph is cyclic (and connected).
        (An exception is raised on disconnected graphs.)
        This function quits early as soon as a cycle is found.

        """
        self.Reset()
        if (type(self.g) is Ugraph):
            is_cyclic = self._IsCyclicUgraph(0, Dgraph.NULL)
        else:
            is_cyclic = self._IsCyclic(0)
        if ((self.sv != self.g.nv) and (not is_cyclic)):
            raise(Disconnected(self.g, "Error(IsCyclic): " +
                               "The input graph is not connected."))
        return is_cyclic

    def _IsCyclicUgraph(self, iv, ivprev):
        """
        _IsCyclicUgraph() is a recursive function which carries out a
        Depth-First-Search over the graph "self.g" to determine whether the
        graph is cyclic.  This function works on undirected graphs (Ugraphs).

        Indirected graphs (Ugraphs) are a special case.
        Ugraphs are implemented by using two (redundant) forward/backward edges
        connecting each pair of adjacent vertices.  This creates trivial loops.
        This version of _IsCyclicUgraph() only counts loops between more
        distantly connected vertices.

        """
        self.sv += 1
        self.vvisited[iv] = True
        for je in self.g.neighbors[iv]:
            jv = self.g.edges[je].stop
            if self.vvisited[jv]:
                if jv != ivprev:
                    return True
            elif self._IsCyclicUgraph(jv, iv):
                return True

        return False

    def _IsCyclic(self, iv):
        """
        _IsCyclic() is a recursive function which carries out a
        Depth-First-Search over the graph "self.g" to determine whether
        the graph is cyclic.
        This function works on directed graphs.

        """
        self.sv += 1
        self.vvisited[iv] = True
        for je in self.g.neighbors[iv]:
            jv = self.g.edges[je].stop
            if self.vvisited[jv]:
                return True
            elif self._IsCyclic(jv):
                return True

        return False


class GraphMatcher(object):
    """
    This class is a variant of the VF2 algorithm for searching
    for small connected subgraphs (g) within a larger graph (G).
    GraphMatcher works on directed or underected graphs (Dgraph or Ugraph).
    This particular version is better optimized for detecting subgraph
    isomorphisms between two graphs of highly unequal size.  It should be
    faster in these situations because, the computation required for
    each step is independent of the number of vertices in the larger graph
    In the original VF2 algorithm, the computation time for each step
    is proportional to the number of vertices in the larger graph.
    (The distinction matters when one graph is much smaller than the other.)

    Limitations: At the moment, the matching process uses a simple
    depth-first-search to search the vertices of the small graph "g".
    Hence this approach fails when the smaller graph g is disconnected.
    (but it can probably be fixed by picking a different algorithm to search
     the small graph).

    """

    def __init__(self,
                 G,  # The "big" graph
                 g):  # The little graph (number of vertices in g must be <= G)

        self.G = G
        self.g = copy.deepcopy(g)

        if (type(self.G) is Ugraph):
            assert(type(self.g) is Ugraph)
        #    self.G.CalcEdgeLookupTable() <-- not needed anymore

        self.sv = 0
        self.se = 0
        self.voccupiedG = [False for iv in range(0, G.nv)]
        self.eoccupiedG = [False for ie in range(0, G.ne)]
        self.G_is_too_small = False
        if ((g.nv > G.nv) or
                (g.ne > G.ne)):
            self.G_is_too_small = True
            # raise GenErr('Error: The first argument of GraphMatcher(G,g),\n'+
            #             '       must be at least as large as the second.')

        # The list self.iv_to_Iv is the mapping between the graph vertices.
        # Iv is an index into the large graph's list of vertices.
        # iv is an index into the small graph's list of vertices.
        # The mapping is stored in the iv_to_Iv list.
        self.iv_to_Iv = [Dgraph.NULL for Iv in range(0, self.g.nv)]
        self.ie_to_Ie = [Dgraph.NULL for Ie in range(0, self.g.ne)]
        #  (This used to be called "core_2" in the VF2 algorithm)

        # Due to the large number of recursion limit
        self.old_recursion_limit = sys.getrecursionlimit()
        expected_max_recursion = self.g.nv

        if self.old_recursion_limit < 1.5 * expected_max_recursion:
            # Give some breathing room.
            sys.setrecursionlimit(int(1.5 * expected_max_recursion))

        subgraph_searcher = DFS(self.g)
        # Perform a Depth-First-Search on the small graph.
        self.vorder_g, self.eorder_g = subgraph_searcher.Order()
        # Then re-order the vertices and edgers to
        # match the order they were visited.

        # Note on permutation order:
        # (The DFS.Order() function returns the permutation in this format
        #       old_index[ new_index ]
        #  where new_index is the DFS iteration when the vertex/edge was visited
        #    and old_index is the original vertex/edge order.
        #  However the ReorderVerts() and ReorderEdges() functions expect
        #  the permutation to have the opposite order: new_index[ old_index ]
        #  Hence we set "invert=True", when we invoke these functions.)
        self.g.ReorderVerts(self.vorder_g, invert=True)
        self.g.ReorderEdges(self.eorder_g, invert=True)

        # Initialize state
        self.Reset()

    def Reset(self):
        """Reinitializes the state of the match-search algorithm.

        """
        for iv in range(0, self.g.nv):
            self.iv_to_Iv[iv] = Dgraph.NULL
        for ie in range(0, self.g.ne):
            self.ie_to_Ie[ie] = Dgraph.NULL
        for Iv in range(0, self.G.nv):
            self.voccupiedG[Iv] = False
        for Ie in range(0, self.G.ne):
            self.eoccupiedG[Ie] = False

        self.se = 0
        self.sv = 0

        # OPTIONAL: First, do a partial sort for the vertices in the graphs
        # based on number of edges emanating from each vertex.
        # (This is probably unnecessary for small subgraphs.)
        # SortVertsByDegree(self.g)

    def Matches(self):
        """
        Iterator over all matches between G and g.
        Each "match" corresponds to a subgraph of G which is isomorphic to g.
        Matches is formatted as a 2-tuple of lists:
           (list of vertex ids from G, list of edge ids from G)
        The vertex ids in the list are a subset of the integers from 0 to G.nv.
        The edge   ids in the list are a subset of the integers from 0 to G.ne.

        (The corresponding vertices and edges from g are indicated by the order)

        """

        self.Reset()
        if self.G_is_too_small:
            # Then there are fewer verts and edges in G than in g.
            # Thus it is impossible for a subgraph of G to be isomorphic to g.
            return  # return no matches

        for Iv in range(0, self.G.nv):

            # match vertex Iv from G with vertex 0 from graph g
            self.iv_to_Iv[0] = Iv
            self.voccupiedG[Iv] = True

            # Implementation:
            # In this loop we begin the search process
            # starting with a different vertex (Iv) from big graph G,
            # and matching it with the first vertex (iv=0) from small graph g.
            # In this way the match "begins" from vertex Iv in G.
            #
            # Any matches found which begin from vertex Iv are distinct
            # from matches beginning from any other vertex in G.
            # Looping over all Iv in G is necessary and sufficient
            # to insure that all possible subgraphs of G
            # (which are isomorphic to g) are considered.

            self.sv = 1  # we have matched one vertex already
            self.se = 0  # we haven't matched any edges yet
            for match in self.Match():
                yield match
            self.voccupiedG[Iv] = False

    def Match(self):

        # self.se represents how many vertices have been matched so far.
        # We are done searching if all of the edges from 0 to self.se-1
        # from graph g have been selected (matched with edges from graph G).
        if self.se == self.g.ne:
            # Note: This also gaurantees that all vertices have been visited.
            assert(self.sv == self.g.nv)
            yield self.ReformatMatch()

        else:
            # VF2-style recursive loop:

            # We know the next edge to be matched is connected to at least
            # one previously visited vertex from g which has already been
            # been added to the the current match-in-progress.
            iv = self.g.edges[self.se].start
            Iv = self.iv_to_Iv[iv]
            assert(iv < self.sv)  # <-- check to verify this is so

            # The other vertex may or may not have been visited (matched) yet.
            iv_neighbor = self.g.edges[self.se].stop

            # Two cases:
            # Case 1: edge self.se points to a previously visited vertex from g
            #         This means we have a loop.
            if iv_neighbor < self.sv:
                # In that case, then the corresponding edge in G must
                # connect the corresponding pair of vertices from G.
                # (Which we know have already been assigned to vertices in g
                #  because both iv and iv_neighbor are < self.sv)
                Iv_neighbor = self.iv_to_Iv[iv_neighbor]
                # Loop over all of the edges in G which connect this pair
                # of vertices (Iv --> Iv_neighbor)
                for Je in self.G.neighbors[Iv]:
                    Jv = self.G.edges[Je].stop
                    if ((Jv == Iv_neighbor) and
                            (not self.eoccupiedG[Je])):

                        # Match edge Je from big   graph G with
                        #  edge self.se from small graph g
                        self.ie_to_Ie[self.se] = Je
                        self.se += 1
                        self.eoccupiedG[Je] = True
                        for match in self.Match():
                            yield match
                        self.eoccupiedG[Je] = False
                        self.se -= 1
                        self.ie_to_Ie[self.se] = Dgraph.NULL

            # Case 2:
            else:  # this would mean that iv_neighbor >= self.sv

                # If iv_neighbor>=self.sv, then this edge points to to a vertex
                # in g which has not yet been paired with a vertex from G.

                # Loop over all of the edges in G which connect vertex
                # Iv from G to new (unvisited) vertices in G
                for Je in self.G.neighbors[Iv]:
                    Jv = self.G.edges[Je].stop
                    if (not self.voccupiedG[Jv]):

                        assert(not self.eoccupiedG[Je])
                        # Match both edge Je with je
                        #      AND vertex Jv with jv
                        self.ie_to_Ie[self.se] = Je
                        self.se += 1
                        self.eoccupiedG[Je] = True
                        self.iv_to_Iv[self.sv] = Jv
                        self.sv += 1
                        self.voccupiedG[Jv] = True
                        # Then continue the recursion
                        for match in self.Match():
                            yield match
                        self.voccupiedG[Jv] = False
                        self.sv -= 1
                        self.iv_to_Iv[self.sv] = Dgraph.NULL
                        self.eoccupiedG[Je] = False
                        self.se -= 1
                        self.ie_to_Ie[self.se] = Dgraph.NULL

    def ReformatMatch(self):
        #   (This is because we are assuming g is connected.
        #    IT should not have any orphanned vertices.)
        # Now return the match:
        #
        # There are different ways of doing this
        # version 1:
        # match = (self.iv_to_Iv, self.ie_to_Ie) <-return a pointer to array
        # version 2:
        # match = ([Iv for Iv in self.iv_to_Iv], <-return a copy of the array
        #         [Ie for Ie in self.ie_to_Ie])
        # version 3:
        #   Recall that the vertices and edges and g have been re-ordered,
        #   so sort the list of Iv indices in the order they would be
        #   matched with the original vertices from the original graph g:
        # match = ([self.iv_to_Iv[self.vorder_g[iv]]
        #          for iv in range(0,self.g.nv)],
        #         [self.ie_to_Ie[self.eorder_g[ie]]
        #          for ie in range(0,self.g.ne)])
        # version 4: Similar to version 3 above, but we also translate
        #            the directed edge id list into a shorter undirected
        #            edge id list.
        match_verts = [self.iv_to_Iv[self.vorder_g[iv]]
                       for iv in range(0, self.g.nv)]

        if type(self.g) is Dgraph:
            match_edges = [self.ie_to_Ie[self.eorder_g[ie]]
                           for ie in range(0, self.g.ne)]
        else:
            #assert(atype(self.g) is Ugraph)
            match_edges = [Dgraph.NULL for ieu in range(0, self.g.neu)]

            for ie in range(0, self.g.ne):
                iv = self.g.edges[ie].start
                jv = self.g.edges[ie].stop
                if iv <= jv:  # <-- avoid duplicating edges (iv,jv) and (jv,iv)
                    ieu = self.g.LookupUndirectedEdgeIdx(ie)
                    Ie = self.ie_to_Ie[ie]
                    Ieu = self.G.LookupUndirectedEdgeIdx(Ie)
                    match_edges[ieu] = Ieu

        return (tuple(match_verts), tuple(match_edges))
