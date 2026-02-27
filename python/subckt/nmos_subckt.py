import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from subckt import *

class nmos_subckt(subckt):

    def __init__(self, instance: str = "M_{1}", **nodes):

        self.g = nodes["g"]
        self.d = nodes["d"]
        self.s = nodes["s"]
        self.instance = instance

    def set_params(self, **params):
        self.KP       = params["KP"]
        self.vth      = params["vth"]
        self.l_lambda = params["lambda"]
        self.L        = params["L"]
        self.W        = params["W"]

        self.CGDO     = 200e-12
        self.scale    = 1e-6
        self.cox_dash = 1.75e-15
        #self.cox      = 35e-15

    def add(self, ckt):
        self.ids = vccs_l1_mosfet()
        self.ids.set_type(ElementType.current_src)
        self.ids.set_instance("i_{DS_{" +self.instance + "}}")
        self.ids.set_params(KP=self.KP, \
                            vth=self.vth, \
                            l_lambda=self.l_lambda, \
                            L=self.L, \
                            W=self.W)

        self.vgs = current_src()
        self.vgs.set_type(ElementType.current_src)
        self.vgs.set_instance("i_{GS_{" +self.instance + "}}")

        self.ids_ref = ckt.add_edge(self.d, self.s, self.ids)
        self.vgs_ref = ckt.add_edge(self.g, self.s, self.vgs)

        self.ids.set_vgs_ref(vgs_ref=self.vgs_ref)

        return self.ids_ref, self.vgs_ref

    def add_small(self, op, ckt_sml, dyn=False, output_resistance=True, **nodes):

        g = nodes["g"]
        d = nodes["d"]
        s = nodes["s"]

        ids = vccs_l1_mosfet_small()
        ids.set_instance("i_{DS_{" +self.instance + "}}")
        ids.set_params(KP=self.KP, \
                       vth=self.vth, \
                       l_lambda=self.l_lambda, \
                       L=self.L, \
                       W=self.W)

        if output_resistance:
            i_ro = resistor()
            i_ro.set_instance("R_{o_{" + self.instance + "}}")
        vgs = current_src()
        vgs.set_is_const()
        vgs.set_value(0.0)
        vgs.set_instance("i_{GS_{" +self.instance + "}}")

        ids_ref = ckt_sml.add_edge(d, s, ids)
        vgs_ref = ckt_sml.add_edge(g, s, vgs)

        if output_resistance:
            ckt_sml.add_edge(d, s, i_ro)

        ids.set_vgs_ref(vgs_ref=vgs_ref)
        ids.set_src_dep_ref(ref=vgs_ref)

        [gm, ro] = ids.get_op_t(op=op, i_x_ref=self.ids.i_x_ref, vgs_ref=self.vgs_ref)

        if output_resistance:
            i_ro.set_value(ro)

        self.ids_sml = ids

        if dyn:
            cox = self.cox_dash * self.W * self.L * self.scale**2

            cox = 35e-15
            Cgs_val = (2/3) * cox
            Cgd_val = self.CGDO * self.W * self.scale

            Cgs = capacitor()
            Cgs.set_value(Cgs_val)
            ckt_sml.add_edge(g, s, Cgs)

            Cgd = capacitor()
            Cgd.set_value(Cgd_val)
            ckt_sml.add_edge(g, d, Cgd)

        return ids_ref, vgs_ref
    
    def get_op_sml(self, op):
        [gm, ro] = self.ids_sml.get_op_t(op=op, i_x_ref=self.ids.i_x_ref, vgs_ref=self.vgs_ref)

        return gm, ro
