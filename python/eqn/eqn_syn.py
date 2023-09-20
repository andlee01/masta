from enum import Enum
import sys, getopt, re, shutil, os, pprint, networkx as nx
import numpy as np

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement

class Circuit:

    def __init__(self):
        self.num_edges = 0
        self.G = nx.MultiGraph()
        self.internal_node_count = -1

    # Add an edge to the Graph.
    # Note that explicit pos_node and neg_node attributes are stored as part of the edge dict
    # as edge iterators return a (u,v) node pair which are not always in the order of (pos,neg).
    # Note that each edge has a unique key. For some types of graph, an edge iterator does not return
    # the edge key. Therefore, store a copy in the dict (ref).
    def add_edge(self, pos_node, neg_node, edge_info):

        edge_info.set_ref(self.num_edges)
        edge_weight = edge_info.get_weight()

        self.G.add_edge(pos_node,
                        neg_node,
                        key=self.num_edges,
                        weight=edge_weight,
                        pos_node=pos_node,
                        neg_node=neg_node,
                        ref=self.num_edges,
                        sys_var_ref=-1,
                        info=edge_info)
        self.num_edges += 1

        return self.num_edges - 1

    def get_edge_info(self, edge_ref):

       for (u,v,d) in self.G.edges(data=True):
           if d['info'].ref == edge_ref:
               return d['info']

    def get_ref_from_instance(self, instance):

        for (u,v,d) in self.G.edges(data=True):
           if d['info'].instance == instance:
               return d['info'].ref

        return -1

    def get_im(self, x, t):

        I = np.zeros(self.num_edges)

        for (u,v,d) in self.G.edges(data=True):
            I[d['info'].get_ref()] = d['info'].get_current(x=x, t=t)

        return I

    def get_vm(self, x, t):

        V = np.zeros(self.num_edges)

        for (u,v,d) in self.G.edges(data=True):
            V[d['info'].get_ref()] = d['info'].get_voltage(x=x, t=t)

        return V

    def get_linalg_mtrx(self, x, sys, t):

        linalg_qf = self.qf.copy()
        linalg_bf = self.bf.copy()
        linalg_b  = np.zeros(self.num_edges)

        #print(linalg_qf)

        for (u,v,d) in self.G.edges(data=True):
            [linalg_qf, linalg_bf, linalg_b] = d['info'].upd_linalg_mtrx(x=x, sys=sys, t=t, linalg_qf=linalg_qf, linalg_bf=linalg_bf, linalg_b=linalg_b)

        return [linalg_qf, linalg_bf, linalg_b]

    def set_value(self, ref, value):

        for (u,v,d) in self.G.edges(data=True):
            if d['info'].get_ref() == ref:
                d['info'].set_value(value)
                return

    def consistency_check(self):

        consistency_failed = False

        if not nx.is_connected(self.G):
            print ("Warning: The graph is not connected")
            consistency_failed = True

        for node in list(self.G.nodes()):
            connected_nodes = 0
            for c_node in nx.all_neighbors(self.G,node):
                connected_nodes = connected_nodes + 1

            if connected_nodes < 2:
                print ("Warning: Node " + str(node) + " is only connected to node " + str(c_node))
                consistency_failed = True

        return consistency_failed

    # Find nodes that are a connected only to current sources and inductors. Such nodes exhibit
    # linear dependence between inductor currents, i.e.:
    #
    # Is + IL1 + IL2 = 0
    #
    # Arbitary values (or initial conditions) for IL1 and IL2 cannot be made. One of these inductors
    # cannot form part of the system vector. An inductor that is a member of the spanning-tree is not
    # a system variable.
    #
    # Current sources and inductors have the two lowest priorities (highest edge weight) for inclusion into the
    # sapnning-tree. By default, inductors have the highest edge weight so that they are forced into the co-tree
    # and therefore form part of the system vector. The exception is that described above. To force one of the
    # inductors into the spanning-tree ahead of a connected current source, its weight is reduced to below that
    # of the current source. This is done for only 1 of the connected inductors.
    #
    # Before:
    # Is  - weight = 3
    # IL1 - weight = 4
    # IL2 - weight = 4
    #
    # After:
    # Is  - weight = 3
    # IL1 - weight = 4
    # IL2 - weight = 2.5
    def chk_linear_ind_dep(self):

        # Iterate through all nodes
        for node in list(self.G.nodes()):

            # Get all edges connected to this node
            edges = self.G.edges(node, data=True)

            # Indicate if the node is connected only to current sources and inductors
            isrc_ind_node = True

            for (u,v,d) in edges:
                edge_info = d['info']

                inductor   = edge_info.get_type() == ElementType.inductor

                indep_isrc = (edge_info.get_type()      == ElementType.current_src) and \
                             (edge_info.get_dependant() == False)

                if not inductor and not indep_isrc:
                    isrc_ind_node = False

            # Current src/inductor node found, reduced weight of first connected inductor
            if isrc_ind_node == True:
                update_made = False
                for (u,v,d) in edges:
                    edge_info = d['info']

                    if edge_info.get_type() == ElementType.inductor and \
                       update_made == False:
                        d['weight'] = 2.5
                        update_made = True

    #def chk_dep_vsrc_cap_loop(self, src_edge_ref):
    #
    #    chk_tree = self.G.copy()
    #
    #    # Iterate through all edges
    #    for (u,v,d) in self.G.edges(data=True):
    #        edge_ref  = d['ref']
    #        edge_info = d['info']
    #
    #        capacitor  = edge_info.get_type() == ElementType.capacitor
    #        src_edge   = edge_ref == src_edge_ref
    #
    #        if not capacitor and not src_edge:
    #           chk_tree.remove_edge(u,v,edge_ref)
    #
    #    # Convert to Graph type as cycle_basic not supported for MultiGraph
    #    loop_tree = nx.Graph()
    #    for (x,y) in chk_tree.edges():
    #        loop_tree.add_edge(x,y)
    #
    #    # Check for loops
    #    loops = nx.cycle_basis(loop_tree)
    #    print (loops)
    #
    #    # Iterate through all edges
    #    for (u,v,d) in chk_tree.edges(data=True):
    #        print (d)

    def get_spanning_tree(self):
        self.spanning_tree = nx.minimum_spanning_tree(self.G)

    def get_co_tree(self):

        self.co_tree = self.G.copy()

        # Remove all edges that are in the spanning-tree from the co-tree
        for (u,v,d) in self.spanning_tree.edges(data=True):
            key = d['ref']
            self.co_tree.remove_edge(u,v,key)

    # Calculates the Qf matrix relating to the network. The algorithm is as follows:
    # 1. For each edge of the spanning-tree, set Qf[ref][ref] = 1
    # 2. For the current spanning-tree edge, add a single co-tree edge.
    # 3. If the addition of the co-tree edge causes a loop to form, then the co-tree edge
    #    forms a basic cut with the current spanning-tree edge (see html description).
    # 4. Remove the added co-tree edge and current spanning-tree edge. The graph will now be
    #    disconnected and will be split into two sets of connected nodes.
    # 5. Determine which set of nodes contains the postive node of the spanning-tree edge.
    #    Positive current will be directed towards the postive node of the spanning-tree edge
    #    (this corresponds to the blue nodes in the html description). If the postive node of the
    #    co-tree edge is also contained in this set of nodes, the Qf entry is positive.
    def get_qf_matrix(self):

        self.qf = np.zeros([self.num_edges, self.num_edges])

        for (u,v,d) in self.spanning_tree.edges(data=True):

            spanning_tree_edge_ref = d['ref']
            spanning_tree_pos_node = d['pos_node']

            # For non-vsrc edges, positive current flows from the positive edge to the negative edge.
            # For vsrc edges, this direction is reversed. Therefore, the defined blue and green node
            # sets need swapped.
            if (d['info'].get_type() == ElementType.voltage_src):
                spanning_tree_pos_node = d['neg_node']

            self.qf[spanning_tree_edge_ref][spanning_tree_edge_ref] = 1

            for (i,j,e) in self.co_tree.edges(data=True):

                # create copy of self.spanning_tree as Graph type as cycle_basis only supports Graph
                tree = nx.Graph()
                for (x,y) in self.spanning_tree.edges():
                    tree.add_edge(x,y)

                parallel_edge = (u == i and v == j) or (u == j and v == i)

                if not parallel_edge:
                    tree.add_edge(i,j)

                loops = nx.cycle_basis(tree)
                loop_found = 0
                if len(loops) != 0:
                    if (u in loops[0] and v in loops[0]):
                        loop_found = 1

                if loop_found or parallel_edge:
                    if not parallel_edge:
                        tree.remove_edge(i,j)
                    tree.remove_edge(u,v)

                    co_tree_edge_ref = e['ref']
                    co_tree_pos_node = e['pos_node']

                    green_surface = nx.node_connected_component(tree, spanning_tree_pos_node)

                    if co_tree_pos_node in green_surface:
                        self.qf[spanning_tree_edge_ref][co_tree_edge_ref] = 1
                    else:
                        self.qf[spanning_tree_edge_ref][co_tree_edge_ref] = -1

    def get_bf_matrix(self):

        self.bf = np.zeros([self.num_edges, self.num_edges])

        for (u,v,d) in self.spanning_tree.edges(data=True):

            spanning_tree_edge_ref = d['ref']
            spanning_tree_pos_node = d['pos_node']
            spanning_tree_neg_node = d['neg_node']

            for (i,j,e) in self.co_tree.edges(data=True):

                # create copy of self.spanning_tree as Graph type as cycle_basis only supports Graph
                tree = nx.Graph()
                for (x,y) in self.spanning_tree.edges():
                    tree.add_edge(x,y)

                parallel_edge = (u == i and v == j) or (u == j and v == i)

                co_tree_edge_ref = e['ref']
                co_tree_pos_node = e['pos_node']
                co_tree_neg_node = e['neg_node']

                self.bf[co_tree_edge_ref][co_tree_edge_ref] = 1

                if parallel_edge:
                    if co_tree_pos_node == spanning_tree_pos_node:
                        self.bf[co_tree_edge_ref][spanning_tree_edge_ref] = -1
                    else:
                        self.bf[co_tree_edge_ref][spanning_tree_edge_ref] = 1
                else:
                    tree.add_edge(i,j)
                    loops = nx.cycle_basis(tree, co_tree_pos_node)
                    loop_found = len(loops) != 0

                    if loop_found:
                        if spanning_tree_pos_node in loops[0] and spanning_tree_neg_node in loops[0]:
                            spanning_tree_pos_node_index = loops[0].index(spanning_tree_pos_node)
                            spanning_tree_neg_node_index = loops[0].index(spanning_tree_neg_node)

                            co_tree_pos_node_index = loops[0].index(co_tree_pos_node)
                            co_tree_neg_node_index = loops[0].index(co_tree_neg_node)

                            spanning_tree_pos_traversal = ((loops[0][len(loops[0])-1] == spanning_tree_pos_node and
                                                            loops[0][0]               == spanning_tree_neg_node) or
                                                           (spanning_tree_pos_node_index == (spanning_tree_neg_node_index - 1)))

                            co_tree_pos_traversal = ((loops[0][len(loops[0])-1] == co_tree_pos_node and
                                                      loops[0][0]               == co_tree_neg_node) or
                                                     (co_tree_pos_node_index == (co_tree_neg_node_index - 1)))

                            if spanning_tree_pos_traversal == co_tree_pos_traversal:
                                self.bf[co_tree_edge_ref][spanning_tree_edge_ref] = 1
                            else:
                                self.bf[co_tree_edge_ref][spanning_tree_edge_ref] = -1

    def init_circuit(self):
        self.consistency_check()
        self.chk_linear_ind_dep()
        self.get_spanning_tree()
        self.get_co_tree()
        self.get_qf_matrix()
        self.get_bf_matrix()
