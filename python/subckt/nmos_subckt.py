import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

class nmos_subckt(subckt):

    def __init__(self, **nodes, **params):

        self.g = nodes["g"]
        self.d = nodes["d"]
        self.s = nodes["s"]

        self.KP       = params["KP"]
        self.vth      = params["vth"]
        self.l_lambda = params["l_lambda"]
        self.L        = params["L"]
        self.W        = params["W"]

    def add(self, ckt):
        self.ids = vccs_l1_mosfet()
        self.ids.set_type(ElementType.current_src)
        self.ids.set_instance("IDS")
        self.ids.set_params(KP=self.KP, \
                            vth=self.vth, \
                            l_lambda=self.l_lambda, \
                            L=self.L, \
                            W=self.W)

        self.vgs = current_src()
        self.vgs.set_type(ElementType.current_src)
        self.vgs.set_instance("VGS")

        self.ids_ref = ckt.add_edge(self.d, self.s, self.ids)
        self.vgs_ref = ckt.add_edge(self.g, self.s, self.vgs)

        ids.set_vgs_ref(vgs_ref=self.vgs_ref)

    def add_small(self, op, ckt_sml, **nodes):

        g = nodes["g"]
        d = nodes["d"]
        s = nodes["s"]

        ids = vccs_l1_mosfet_small()
        ids.set_instance("IDS")
        ids.set_params(KP=self.KP, \
                       vth=self.vth, \
                       l_lambda=self.l_lambda, \
                       L=self.L, \
                       W=self.W)

        i_ro = resistor()
        i_ro.set_type(ElementType.resistor)
        i_ro.set_instance("Ro")

        vgs = current_src()
        vgs.set_type(ElementType.current_src)
        vgs.set_is_const()
        vgs.set_value(0.0)
        vgs.set_instance("VGS")

        ids_ref = ckt_sml.add_edge(d, s, ids)
        vgs_ref = ckt_sml.add_edge(g, s, vgs)
        ckt_sml.add_edge(d, s, i_ro)

        ids.set_vgs_ref(vgs_ref=vgs_ref)
        ids.set_src_dep_ref(ref=vgs_ref)

        [gm, ro] = ids.get_op_t(op=op, i_x_ref=self.ids.i_x_ref, vgs_ref=self.vgs_ref)
        i_ro.set_value(ro)
