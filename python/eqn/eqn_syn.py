from enum import Enum
import sys, getopt, shutil, os, pprint, networkx as nx
import numpy as np
from sympy import *

from numpy.linalg import inv

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement

import re as regex

class Scoreboard:

    def __init__(self, num_edges, num_sys_vars, degen_mtrx):
        self.v          = np.zeros(num_edges)
        self.i          = np.zeros(num_edges)
        self.x          = np.zeros(num_edges)
        self.sys        = np.zeros(num_sys_vars)
        self.t          = 0
        self.dv         = np.zeros(num_edges)
        self.di         = np.zeros(num_edges)
        self.degen_mtrx = degen_mtrx


class Circuit:

    def __init__(self, lti=False):
        self.num_edges = 0
        self.G = nx.MultiGraph()
        self.internal_node_count = -1

        # Indication if the Circuit is LTI
        #  - When true disables insertion of depedency variables for dependent sources
        self.lti = lti
        self.num_outputs_lti = 0

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

    def get_edge_info_from_sys_var_ref(self, sys_var_ref):

       for (u,v,d) in self.G.edges(data=True):
           if d['info'].sys_var_ref == sys_var_ref:
               return d['info']

    def get_ref_from_instance(self, instance):

        for (u,v,d) in self.G.edges(data=True):
           if d['info'].instance == instance:
               return d['info'].ref

        return -1

    def get_im(self, x, sys, t):

        self.scb.i   = np.zeros(self.num_edges)
        self.scb.x   = x
        self.scb.sys = sys
        self.scb.t   = t

        for (u,v,d) in self.G.edges(data=True):
            d['info'].get_current(self.scb)

        return self.scb.i

    def get_im_dep(self):

        for (u,v,d) in self.G.edges(data=True):
            d['info'].get_dependent_current(self.scb)

        return self.scb.i

    def get_degen(self):

        for (u,v,d) in self.G.edges(data=True):

            if d['info'].get_degen_ref() != -1:
                d['info'].get_degen_current(self.scb)

        return self.scb.i

    def get_vm(self, x, sys, t):

        self.scb.v = np.zeros(self.num_edges)

        for (u,v,d) in self.G.edges(data=True):
            d['info'].get_voltage(self.scb)

        self.get_ddt_vm()

        return self.scb.v

    def get_ddt_vm(self):

        self.scb.dv = np.zeros(self.num_edges)

        for (u,v,d) in self.spanning_tree.edges(data=True):
            edge = d['info']

            if edge.get_type() == ElementType.capacitor or \
               edge.get_type() == ElementType.voltage_src:
               edge.get_ddt(self.scb)
            else:
                # Other types of element can be members of the spanning tree, but only for
                # the case of voltage-capacitor loops is the edge voltage derivative required.
                self.scb.dv[d['info'].ref] = 0.0

    def get_dy(self, x):

        dy = np.zeros(self.num_sys_vars)

        for (u,v,d) in self.G.edges(data=True):
            dy = d['info'].get_dy(x=x, sys=dy)

        return dy

    def get_linalg_mtrx(self, x, sys, t):

        linalg_qf = self.qf.copy()
        linalg_bf = self.bf.copy()
        linalg_b  = np.zeros(self.num_edges)

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

                #indep_isrc = (edge_info.get_type()      == ElementType.current_src) and \
                #             (edge_info.get_dependant() == False)

                indep_isrc = (edge_info.get_type()      == ElementType.current_src)

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

    # Each edge has a unique ref. Capacitors and inductors that are members of the spanning-tree
    # and co-tree respectively are system variables. Assign these branches a unique sys_var_ref.
    # Since capacitors and inductors which are system variables are evaluated as part of distinct
    # matrices (Qf and Bf), the references both start at zero.
    def set_sys_var_ref(self):

        sys_var_cap_ref = 0
        sys_var_ind_ref = 0

        for (u,v,d) in self.spanning_tree.edges(data=True):
            edge = d['info']

            if edge.get_type() == ElementType.capacitor:
                edge.set_sys_var(ref=sys_var_cap_ref)
                sys_var_cap_ref  += 1

        for (u,v,d) in self.co_tree.edges(data=True):
            edge = d['info']

            if edge.get_type() == ElementType.inductor:
                edge.set_sys_var(ref=(sys_var_cap_ref + sys_var_ind_ref) )
                sys_var_ind_ref  += 1

        self.num_sys_var_cap = sys_var_cap_ref
        self.num_sys_var_ind = sys_var_ind_ref

        self.num_sys_vars = self.num_sys_var_cap + self.num_sys_var_ind


    def set_const_sys_var_ref(self):

        for (u,v,d) in self.G.edges(data=True):
            elem = d['info']

            if elem.get_is_const():
                elem.set_sys_var(self.num_sys_vars)
                self.num_sys_vars += 1

    def set_lti_outputs(self):

        self.num_outputs_lti = 0

        for (u,v,d) in self.G.edges(data=True):
            elem = d['info']
            if elem.get_is_output():
                elem.set_output_ref(self.num_outputs_lti)
                self.num_outputs_lti += 1

    def set_denerate_ref(self):

        degen_ref = 0

        self.degen_mtrx = np.zeros((self.num_edges, self.num_edges))

        for (u,v,d) in self.co_tree.edges(data=True):
            edge = d['info']

            if edge.get_type() == ElementType.capacitor:
                edge.set_degen_ref(ref=self.num_edges)
                degen_ref += 1
                edge_ref   = edge.get_ref()

                # Add additional row/column to Qf, Bf and degen
                self.qf = np.vstack([self.qf, np.zeros( self.qf.shape[0])])
                self.qf = np.hstack((self.qf, np.zeros((self.qf.shape[0], 1))))

                self.bf = np.vstack([self.bf, np.zeros( self.bf.shape[0])])
                self.bf = np.hstack((self.bf, np.zeros((self.bf.shape[0], 1))))

                self.degen_mtrx = np.vstack([self.degen_mtrx, np.zeros( self.degen_mtrx.shape[0])])
                self.degen_mtrx = np.hstack((self.degen_mtrx, np.zeros((self.degen_mtrx.shape[0], 1))))

                # Update Qf for degen
                self.qf[self.num_edges][edge_ref]       = 1
                self.qf[self.num_edges][self.num_edges] = -1

                # Update degen matrix
                self.degen_mtrx[edge_ref,:] = self.bf[edge_ref,:]

                # Mask Bf equation
                self.bf[edge_ref,:] = np.zeros(self.bf.shape[0])

                self.num_edges += 1

    def get_num_sys_vars(self):
        return self.num_sys_vars

    # Determine if a particular edge ref is a resistor in the spanning-tree.
    def is_spanning_tree_res(self, ref):

        for (u,v,d) in self.spanning_tree.edges(data=True):
            edge     = d['info']
            edge_ref = d['ref']

            if edge.get_type() == ElementType.resistor and ref == edge_ref:
                return True

        return False

    # Determine if a particular edge ref is a resistor in the co-tree.
    def is_co_tree_res(self, ref):

        for (u,v,d) in self.co_tree.edges(data=True):
            edge     = d['info']
            edge_ref = d['ref']

            if edge.get_type() == ElementType.resistor and ref == edge_ref:
                return True

        return False

    def is_sys_var_cap(self, ref):

        for (u,v,d) in self.spanning_tree.edges(data=True):
            edge     = d['info']
            edge_ref = d['ref']

            if edge.get_type() == ElementType.capacitor and ref == edge_ref:
                return True

        return False

    def is_non_sys_var_cap(self, ref):

        for (u,v,d) in self.co_tree.edges(data=True):
            edge     = d['info']
            edge_ref = d['ref']

            if edge.get_type() == ElementType.capacitor and ref == edge_ref:
                return True

        return False

    # Determine if a particular edge ref is a system variable inductor.
    # That is, the inductor is a member of the co-tree.
    def is_sys_var_ind(self, ref):

        for (u,v,d) in self.co_tree.edges(data=True):
            edge     = d['info']
            edge_ref = d['ref']

            if edge.get_type() == ElementType.inductor and ref == edge_ref:
                return True

        return False

    def is_non_sys_var_ind(self, ref):

        for (u,v,d) in self.spanning_tree.edges(data=True):
            edge     = d['info']
            edge_ref = d['ref']

            if edge.get_type() == ElementType.inductor and ref == edge_ref:
                return True

        return False

    def is_current_source(self, ref):

        for (u,v,d) in self.co_tree.edges(data=True):
            edge     = d['info']
            edge_ref = d['ref']

            if edge.get_type() == ElementType.current_src and ref == edge_ref:
                return True

        for (u,v,d) in self.spanning_tree.edges(data=True):
            edge     = d['info']
            edge_ref = d['ref']

            if edge.get_type() == ElementType.current_src and ref == edge_ref:
                return True

        return False

    def is_voltage_source(self, ref):

        for (u,v,d) in self.co_tree.edges(data=True):
            edge     = d['info']
            edge_ref = d['ref']

            if edge.get_type() == ElementType.voltage_src and ref == edge_ref:
                return True

        for (u,v,d) in self.spanning_tree.edges(data=True):
            edge     = d['info']
            edge_ref = d['ref']

            if edge.get_type() == ElementType.voltage_src and ref == edge_ref:
                return True

        return False

    def is_co_tree_elem(self, ref):

        for (u,v,d) in self.co_tree.edges(data=True):
            edge     = d['info']
            edge_ref = d['ref']

            if ref == edge_ref:
                return True

    def is_output_lti(self, ref):

        for (u,v,d) in self.co_tree.edges(data=True):
            edge     = d['info']
            edge_ref = d['ref']

            if edge_ref == ref and edge.get_is_output():
                return True

        for (u,v,d) in self.spanning_tree.edges(data=True):
            edge     = d['info']
            edge_ref = d['ref']

            if edge_ref == ref and edge.get_is_output():
                return True

        return False

    def extract_numeric_values(self, expression):
        numeric_values = regex.findall(r'[-+]?\d+\.\d+', expression)
        numeric_values = [float(value) for value in numeric_values]
        return numeric_values

    def separate_numeric_non_numeric(self, input_string):
        non_numeric_list = regex.findall(r'\b(?<!\d)(?<![a-zA-Z_])[a-zA-Z_]+[a-zA-Z0-9_]*\b', input_string)
        return non_numeric_list

    def get_ss(self):

        A_ss, B_ss, C_ss = self.get_ss_A()
        var_string, var_list, sys_var_list, inp_var_list = self.get_ss_sym_matrices()

        x    = Matrix(sys_var_list)
        u    = Matrix(inp_var_list)

        A    = Matrix(A_ss)
        B    = Matrix(B_ss)
        C    = Matrix(C_ss)

        # Vaux = vector of circuit variables
        #        For capacitors, the variable is ddt(V_C<n>)
        #        For inductors, the variable is ddt(i_L<n>)

        # A.Vaux + B.x + C.u = 0
        # Vaux = -inv(A).B.x - inv(A).C.u

        Vaux = -A.inv().multiply(B).multiply(x) - A.inv().multiply(C).multiply(u)

        # Default state space system matrices
        self.A = np.zeros((self.num_sys_vars, self.num_sys_vars))
        self.B = np.zeros((self.num_sys_vars, 1))

        if self.num_outputs_lti != 0:
            self.C = np.zeros((self.num_outputs_lti, self.num_sys_vars))
            self.D = np.zeros((self.num_outputs_lti, 1))
        else:
            self.C = np.eye(self.num_sys_vars)
            self.D = np.zeros((self.num_sys_vars, 1))

        # Loop through resulting equations (Vaux = ...)
        #  - Select the equations that correspond to a system variable (cap or ind)
        #  - Format the equations into the matrices, e.g.
        #
        #    1000.0*u0 + 1000.001*x0 + 1000.0*x1 - 1000.001*x2 + 1.0*x3
        #
        #    is formatted into A[row] = [1000.001 1000.0 -1000.001 1.0]
        #                      B[row] = [1000.0]
        #
        #    Assuming that the system has 4 state variables and 1 input

        # Offset for non sys var capacitors or inductors
        #  - Non sys var capacitors and inductors are removed from the system of equations
        #  - Therefore, the element index of any later edges in the list needs corrected by the number
        #    of non sys var edges before it
        #
        #  Before: Vs - 0                After: Vs - 0   elem_offset = 0
        #          C1 - 1                       C1 - 1   elem_offset = 0
        #          C2 - 2 non-sys               R2 - 2   elem_offset = 1
        #          R2 - 3                       L1 - 3   elem_offset = 2
        #          C3 - 4 non-sys
        #          L1 - 5
        elem_offset = 0

        for elem in range(self.num_edges):

            if self.is_non_sys_var_cap(elem) or self.is_non_sys_var_ind(elem):
                elem_offset += 1

            elif self.is_sys_var_cap(elem) or self.is_sys_var_ind(elem):

                elem_sel = elem - elem_offset

                # Determine the sys var index
                sys_var_ref = self.get_edge_info(elem).get_sys_var_ref()

                str_expr = str(Vaux[elem_sel])
                str_expr = regex.sub(r"\s+", "", str_expr)

                # Parse the symbolic equation string into 2 lists
                # coeffs = [1000.0 1000.001 1000.0 -1000.001 1.0]
                # vect   = [u0 x0 x1 x2 x3]
                coeffs = self.extract_numeric_values(str_expr)
                vect   = self.separate_numeric_non_numeric(str_expr)

                # Parse the above lists and update the state space A and B matrices
                self.parse_sym_list(coeffs, vect, sys_var_ref, state=True)

            if self.is_output_lti(elem):

                if self.is_sys_var_cap(elem):
                    output_ref  = self.get_edge_info(elem).get_output_ref()
                    sys_var_ref = self.get_edge_info(elem).get_sys_var_ref()
                    self.C[output_ref][sys_var_ref] = 1

    def parse_sym_list(self, coeffs, vect, sys_var_ref, state):

        for i in range(len(coeffs)):

            # Seperate each vect element into matrix and element:
            #  - x0  --> x, 0
            #  - x12 --> x, 12
            #  - u0  --> u, 0

            # Get Matrix name (x or u)
            mtrx = regex.findall(r'([a-zA-Z]+)', str(vect[i]))

            # Get Index
            index = regex.findall(r'(\d+)', str(vect[i]))

            if state:
                if mtrx[0] == 'x':
                    self.A[sys_var_ref][int(index[0])] = coeffs[i]
                else:
                    self.B[sys_var_ref][int(index[0])] = coeffs[i]
            else:
                if mtrx[0] == 'x':
                    self.C[sys_var_ref][int(index[0])] = coeffs[i]
                else:
                    self.D[sys_var_ref][int(index[0])] = coeffs[i]

    def get_ss_A(self):

        qf_A_ss = self.qf.copy()
        bf_A_ss = self.bf.copy()

        # For numerical precision reasons, use the long double type.
        # Element values can have large ranges. For example, resistors could be 10M and capacitor values
        # could be 20f. Values with such large differences in magnitude represent a problem for double precision
        # values. That is, they may lack the required precision to adequately represent the numbers, especially when added.
        qf_A_ss = qf_A_ss.astype(np.longdouble)
        bf_A_ss = bf_A_ss.astype(np.longdouble)
        B_ss    = np.zeros([self.num_edges, self.num_sys_vars], dtype=np.longdouble)

        non_sys = np.zeros([self.num_edges, self.num_edges], dtype=np.longdouble)

        # REVISIT (needs to be num_inputs)
        C_ss    = np.zeros([self.num_edges, 1], dtype=np.longdouble)

        for elem in range(self.num_edges):

            # Remove sys var ind currents from A into B
            #  - i_L<n> is a system variable
            #
            # Qf
            # ----
            #               +---+
            #        Vs  C2 |L1 |  C1  C3  R2  RS
            # Vs     1.  0. |0. |  0.  0.  0. -1
            # C2     0.  1. |1. |  0.  0.  0. -1
            # L1     0.  0. |0. |  0.  0.  0.  0
            # C1     0.  0. |0. |  1.  0.  1. -1
            # C3     0.  0. |-1.|  0.  1. -1.  1
            # R2     0.  0. |0. |  0.  0.  0.  0
            # RS     0.  0. |0. |  0.  0.  0.  0
            #               +---+
            #
            # B_ss
            # -----
            #                         +----+
            #        v_C2  v_C1  v_C3 |i_L1|
            #          0     0     0  |  0 |
            #          0     0     0  |  1 |
            #          0     0     0  |  0 |
            #          0     0     0  |  0 |
            #          0     0     0  |  -1|
            #          0     0     0  |  0 |
            #          0     0     0  |  0 |
            #                         +----+
            #
            # Co-tree Resistors require V_R in their Bf equation but the system
            # variable associated with a resistor is i_R.
            #  - Therefore, Bf equation must specify i_R * R (to equal V_R)
            #  - Modify references in associated columns to multiply by R
            # Bf
            # ----
            #                           +------+
            #        Vs  C2  L1  C1  C3 |R2  RS|
            # Vs     0.  0.  0.  0.  0. |0.  0.|
            # C2     0.  0.  0.  0.  0. |0.  0.|
            # L1     0. -1.  1.  0.  1. |0.  0.|
            # C1     0.  0.  0.  0.  0. |0.  0.|
            # C3     0.  0.  0.  0.  0. |0.  0.|
            # R2     0.  0.  0. -1.  1. |1.  0.|
            # RS    -1.  1.  0.  1. -1. |0.  1.|
            #                           +------+
            #
            #                           +------+
            #        Vs  C2  L1  C1  C3 |R2  RS|
            # Vs     0.  0.  0.  0.  0. |0.  0.|
            # C2     0.  0.  0.  0.  0. |0.  0.|
            # L1     0. -1.  1.  0.  1. |0.  0.|
            # C1     0.  0.  0.  0.  0. |0.  0.|
            # C3     0.  0.  0.  0.  0. |0.  0.|
            # R2     0.  0.  0. -1.  1. |R2  0.|
            # RS    -1.  1.  0.  1. -1. |0.  RS|
            #                           +------+

            if self.is_sys_var_ind(elem):

                # Get the Qf matrix column
                qf_column = qf_A_ss[:,elem]

                # Determine the column of B to update (i.e. sys var index)
                sys_var_ref = self.get_edge_info(elem).get_sys_var_ref()

                # Update B column
                B_ss[:,sys_var_ref] = qf_column

                # Mask A column
                qf_A_ss[:,elem] = np.zeros(self.num_edges)

                # Change references to V_L<n> in bf to L<n> * ddt(i_L<n>)
                #  - This implies that the variable inside Vaux is ddt(i_L<n>)
                #  - Therefore, the bf column related to this inductor should be multiplied
                #    by the inductor value
                bf_A_ss[:,elem] = bf_A_ss[:,elem] * self.get_edge_info(elem).get_value()

            elif self.is_non_sys_var_ind(elem):
                assert False, "Non sys var inductors not supported"

            elif self.is_sys_var_cap(elem):

                # Get the Bf matrix column
                bf_column = bf_A_ss[:,elem]

                # Determine the column of B to update (i.e. sys var index)
                sys_var_ref = self.get_edge_info(elem).get_sys_var_ref()

                # Update B column
                B_ss[:,sys_var_ref] = bf_column

                # Mask A column
                bf_A_ss[:,elem] = np.zeros(self.num_edges)

                # Change references to i_C<n> in qf to C<n> * ddt(V_C<n>)
                #  - This implies that the variable inside Vaux is ddt(V_C<n>)
                #  - Therefore, the qf column related to this capacitor should be multiplied
                #    by the capacitor value
                qf_A_ss[:,elem] = qf_A_ss[:,elem] * self.get_edge_info(elem).get_value()

            elif self.is_non_sys_var_cap(elem):

                # Get the Bf matrix row
                bf_row = self.bf[elem ,:].copy()

                for elem_chk in range(self.num_edges):
                    if bf_row[elem_chk] and elem_chk != elem:
                        assert self.is_sys_var_cap(elem_chk), "Only capacitors allowed in  loops with co-tree capacitors " + str(elem_chk)

                # Mask references to self
                #  i.e. VC3 - VC2 + VC1 = 0
                #       If the element being processed is VC3, remove from equation
                #       VC2 + VC1 = 0

                # Assume element being process is VC3 and VC3 - VC2 + VC1 = 0.
                #  - To determine current through C3:
                #       ddt(VC3) - ddt(VC2) + ddt(VC1) = 0
                #       ddt(VC3) = ddt(VC2) - ddt(VC1)
                #       (iC3/C3) = ddt(VC2) - ddt(VC1)
                #        iC3     = C3 * (ddt(VC2) - ddt(VC1))
                #
                #       VC2 VC3    VC1
                #           +--+
                #  1 -1  0  |0 | 0  0
                #  0  0  0  |0 | 0  0
                #  0 -1  1  |1 | 0  0
                #  0  0  0  |0 | 0  0
                #  0  0  0  |0 | 0  0
                #  0  0  0  |-1| 1  1
                #           +--+
                #
                # In the above example, Qf equations 2 and 6 require iC3 as a variable.
                # The following matrix should therefore be added to qf_A_ss:
                #
                # 0 0   0 0 0 0
                # 0 0   0 0 0 0
                # 0 0  C3 0 0 -C3
                # 0 0   0 0 0 0
                # 0 0   0 0 0 0
                # 0 0 -C3 0 0 C3
                #

                # Form the base bf row to add to qf_A_ss
                bf_row[elem] = 0
                bf_row = -1 * (bf_row * self.get_edge_info(elem).get_value())

                # Multiply the base_bf_row by the qf column corresponding to the element to form the addition
                # matrix
                #
                # [ 0 0  C3 0 0 -C3 ] * [ 0] =  [0 0   0 0 0 0  ]
                #                       [ 0]    [0 0   0 0 0 0  ]
                #                       [ 1]    [0 0  C3 0 0 -C3]
                #                       [ 0]    [0 0   0 0 0 0  ]
                #                       [ 0]    [0 0   0 0 0 0  ]
                #                       [-1]    [0 0 -C3 0 0 C3 ]
                non_sys += bf_row * self.qf[:,[elem]]

                # Remove any references to this element:
                #  - bf_A_ss remove the equation i.e. delete VC3 - VC2 + VC1 = 0
                #  - qf_A_ss remove references to column iC3 in qf
                bf_A_ss[elem ,:] = 0
                qf_A_ss[:, elem] = 0

            elif self.is_co_tree_res(elem) or self.is_spanning_tree_res(elem):

                # Get the Bf matrix column
                bf_column = bf_A_ss[:,elem]

                # Get the resistance
                res_val = self.get_edge_info(elem).get_value()

                # Multiply by resistance
                bf_column = bf_column * res_val

                # write back to column
                bf_A_ss[:,elem] = bf_column

            elif self.is_voltage_source(elem):

                if self.get_edge_info(elem).get_is_input():

                    # Get the Bf matrix column
                    bf_column = bf_A_ss[:,elem]

                    # Get the input index
                    input_ref = 0

                    # Update C column
                    C_ss[:,input_ref] = bf_column

                    # Mask Bf column
                    bf_A_ss[:,elem] = np.zeros(self.num_edges)

                elif self.get_edge_info(elem).get_is_const():

                    # Get the Bf matrix column
                    bf_column = bf_A_ss[:,elem]

                    # Determine the column of B to update (i.e. sys var index)
                    sys_var_ref = self.get_edge_info(elem).get_sys_var_ref()

                    # Update B column
                    B_ss[:,sys_var_ref] = bf_column

                    # Mask A column
                    bf_A_ss[:,elem] = np.zeros(self.num_edges)

                else:
                    assert False, "Dependent voltage sources not yet supported"

            elif self.is_current_source(elem):

                if self.get_edge_info(elem).get_is_input():

                    # Get the Qf matrix column
                    qf_column = qf_A_ss[:,elem]

                    # Get the input index
                    input_ref = 0

                    # Update C column
                    C_ss[:,input_ref] = qf_column

                    # Mask Qf column
                    qf_A_ss[:,elem] = np.zeros(self.num_edges)

                elif self.get_edge_info(elem).get_is_const():

                    # Get the Qf matrix column
                    qf_column = self.qf[:,elem].copy()

                    # Determine the column of B to update (i.e. sys var index)
                    sys_var_ref = self.get_edge_info(elem).get_sys_var_ref()

                    # Update B column
                    B_ss[:,sys_var_ref] = qf_column

                    # Mask A column
                    qf_A_ss[:,elem] = qf_A_ss[:,elem] - qf_column

                else:

                    # Dependent current source (i1)
                    #
                    # i1 is a dependent current source with value gm * control_voltage
                    #
                    #                           +--+
                    #        Vs  C2  L1  C1  C3 |i1| RS  i2
                    #        -1  1   0   0   0  |1 | 0   0
                    #        0   0   1   1   0  |-1| 1   0
                    #        0   0   1   0   1  |1 | 1   1
                    #                           +--+
                    # To formulate the qf equations correctly, the gm must appear in the column of the control voltage.
                    # Depending on the control element, the control voltage column differs.
                    #
                    # If the control voltage is a current source, then the Vaux variable is the element voltage.
                    #   - Assume that the control voltage is the element voltage of current source i2. The qf equations
                    #     therefore become:
                    #
                    #                           +--+    +****+
                    #        Vs  C2  L1  C1  C3 |i1| RS |i2  |
                    #        -1  1   0   0   0  |1 | 0  |gm  |
                    #        0   0   1   1   0  |-1| 1  |-gm |
                    #        0   0   1   0   1  |1 | 1  |1+gm|
                    #                           +--+    +****+
                    #
                    # The next stage of the process involves removing i1 references from the qf equations.
                    # This cannot simply be forced to zero in case another dependent element uses the element voltage of
                    # i1 as its control voltage. Therefore, the original qf column should be subtracted. This preserves
                    # any additional gm terms added there.
                    #
                    #                           +--+    +****+
                    #        Vs  C2  L1  C1  C3 |i1| RS |i2  |
                    #        -1  1   0   0   0  |0 | 0  |gm  |
                    #        0   0   1   1   0  |0 | 1  |-gm |
                    #        0   0   1   0   1  |0 | 1  |1+gm|
                    #                           +--+    +****+
                    # Control voltage edge
                    dep_ref = self.get_edge_info(elem).get_src_dep_ref()

                    assert self.is_current_source(dep_ref), "Only current sources as dep ref are currently supported"

                    # Control voltage is from an element whose circuit variable is the element voltage

                    # Get the Qf matrix column
                    #  - Must use the original qf column here in case another dependent source has added its gm
                    qf_column = self.qf[:,elem].copy()

                    # Get the Qf matrix column corresponding to the control voltage element
                    #  - This column may not be clear if other qf equations require the current through this element
                    #  - Therefore, the dependent current source contribution must be added to the column values
                    #    already there
                    qf_A_ss[:,dep_ref] = qf_A_ss[:,dep_ref] + (qf_column * self.get_edge_info(elem).get_value())

                    # The original references to the current source in the qf matrix must be masked out
                    qf_A_ss[:,elem] = qf_A_ss[:,elem] - qf_column
            else:
                assert False, "Unmapped element"

        # Combine the final qf and bf equations
        A_ss = qf_A_ss + bf_A_ss + non_sys

        elem_offset = 0
        for elem in range(self.num_edges):
            if self.is_non_sys_var_ind(elem) or self.is_non_sys_var_cap(elem):
                A_ss = np.delete(A_ss, (elem - elem_offset), axis=0)
                A_ss = np.delete(A_ss, (elem - elem_offset), axis=1)

                B_ss = np.delete(B_ss, (elem - elem_offset), axis=0)
                C_ss = np.delete(C_ss, (elem - elem_offset), axis=0)

                elem_offset += 1

        return A_ss, B_ss, C_ss

    def get_ss_sym_matrices(self):

        # Create sympy variable string
        var_string = ""

        # Lists to form sympy Matrix
        var_list = []
        sys_var_list = []
        inp_var_list = []

        # 1 variable for each edge
        for elem in range(self.num_edges):
            var_string += "v" + str(elem) + " "
            var_list.append("v" + str(elem))

        # 1 variable for each sys var
        for elem in range(self.num_sys_vars):
            var_string += "x" + str(elem) + " "
            sys_var_list.append("x" + str(elem))

        # 1 variable for each input
        var_string += "u" + str(0) + " "
        inp_var_list.append("u" + str(0))

        # Add remaining symbolic matrices
        var_string += "A B C Vaux x u"

        return var_string, var_list, sys_var_list, inp_var_list

    def get_lti_const_x0(self, x0):

        for elem in range(self.num_edges):

            if (self.is_current_source(elem) and self.get_edge_info(elem).get_is_const()) or \
               (self.is_voltage_source(elem) and self.get_edge_info(elem).get_is_const()):

                sys_var_ref = self.get_edge_info(elem).get_sys_var_ref()

                x0[sys_var_ref] = self.get_edge_info(elem).get_value()

        return x0

    def init_circuit(self):
        self.consistency_check()
        self.chk_linear_ind_dep()
        self.get_spanning_tree()
        self.get_co_tree()
        self.get_qf_matrix()
        self.get_bf_matrix()
        self.set_sys_var_ref()

        self.degen_mtrx = 0

        if not self.lti:
            self.set_denerate_ref()
            for (u,v,d) in self.G.edges(data=True):
                d['info'].set_dependencies(self)

        if self.lti:
            self.set_const_sys_var_ref()
            self.set_lti_outputs()

        self.scb = Scoreboard(num_edges=self.num_edges, num_sys_vars=self.num_sys_vars, degen_mtrx=self.degen_mtrx)
