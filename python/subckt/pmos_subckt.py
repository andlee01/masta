import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from subckt import *

class pmos_subckt(subckt):

    def __init__(self, **nodes):

        self.g = nodes["g"]
        self.d = nodes["d"]
        self.s = nodes["s"]

    def set_params(self, **params):
        self.KP       = params["KP"]
        self.vth      = params["vth"]
        self.l_lambda = params["lambda"]
        self.L        = params["L"]
        self.W        = params["W"]

    def add(self, ckt):
        self.isd = vccs_l1_mosfet()
        self.isd.set_type(ElementType.current_src)
        self.isd.set_instance("IDS")
        self.isd.set_params(KP=self.KP, \
                            vth=self.vth, \
                            l_lambda=self.l_lambda, \
                            L=self.L, \
                            W=self.W)

        self.vsg = current_src()
        self.vsg.set_type(ElementType.current_src)
        self.vsg.set_instance("VGS")

        self.isd_ref = ckt.add_edge(self.s, self.d, self.isd)
        self.vsg_ref = ckt.add_edge(self.s, self.g, self.vsg)

        self.isd.set_vgs_ref(vgs_ref=self.vsg_ref)

    def add_small(self, op, ckt_sml, **nodes):

        g = nodes["g"]
        d = nodes["d"]
        s = nodes["s"]

        isd = vccs_l1_mosfet_small()
        isd.set_instance("IDS")
        isd.set_params(KP=self.KP, \
                       vth=self.vth, \
                       l_lambda=self.l_lambda, \
                       L=self.L, \
                       W=self.W)

        i_ro = resistor()
        i_ro.set_type(ElementType.resistor)
        i_ro.set_instance("Ro")

        vsg = current_src()
        vsg.set_type(ElementType.current_src)
        vsg.set_is_const()
        vsg.set_value(0.0)
        vsg.set_instance("VGS")

        isd_ref = ckt_sml.add_edge(s, d, isd)
        vsg_ref = ckt_sml.add_edge(s, g, vsg)
        ckt_sml.add_edge(s, d, i_ro)

        isd.set_vgs_ref(vgs_ref=vsg_ref)
        isd.set_src_dep_ref(ref=vsg_ref)

        [gm, ro] = isd.get_op_t(op=op, i_x_ref=self.isd.i_x_ref, vgs_ref=self.vsg_ref)
        i_ro.set_value(ro)
