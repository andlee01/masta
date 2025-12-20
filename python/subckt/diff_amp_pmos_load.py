import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from subckt import *
from nmos_subckt import *
from pmos_subckt import *

class diff_amp_pmos_load(subckt):

    def __init__(self, neg_feedback=False, pos_feedback=False, **nodes):

        self.VCC  = nodes["vcc"]
        self.GND  = nodes["gnd"]
        self.vlop = nodes["vlop"]
        self.vlon = nodes["vlon"]
        self.vfp  = nodes["vfp"]
        self.vfn  = nodes["vfn"]

        self.neg_feedback = neg_feedback
        self.pos_feedback = pos_feedback

    def set_params(self, **params):

        self.iref_n.set_value(params["Iref"])

        if self.neg_feedback:
            self.Rf.set_value(params["Rf"])
            self.Rin.set_value(params["Rin"])

    def set_params_sml(self, **params):

        if self.neg_feedback:
            self.Rf_sml.set_value(params["Rf"])
            self.Rin_sml.set_value(params["Rin"])

    def add(self, ckt):

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 10}

        pmos_params = {"KP"     : 40e-6,
                       "vth"    : 0.9,
                       "lambda" : 0.0125,
                       "L"      : 2,
                       "W"      : 30}

        n1     = ckt.get_internal_node()

        if self.neg_feedback or self.pos_feedback:
            n2 = ckt.get_internal_node()

        # Add nmos
        nodes = {"g": n2 if self.neg_feedback else self.vlon, "d": self.vfn, "s": n1}
        self.nmos_m2 = nmos_subckt(**nodes)
        self.nmos_m2.set_params(**nmos_params)
        self.ids_ref_m2, self.vgs_ref_m2 = self.nmos_m2.add(ckt)

        # Add nmos
        nodes = {"g": self.vlop, "d": self.vfp, "s": n1}
        self.nmos_m1 = nmos_subckt(**nodes)
        self.nmos_m1.set_params(**nmos_params)
        self.ids_ref_m1, self.vgs_ref_m1 = self.nmos_m1.add(ckt)

        # Add Iref
        self.iref_n = current_src()
        self.iref_n.set_instance("iref_n")
        ckt.add_edge(n1, self.GND, self.iref_n)

        # Add pmos
        nodes = {"g": self.vfp, "d": self.vfp, "s": self.VCC}
        self.pmos_m3 = pmos_subckt(**nodes)
        self.pmos_m3.set_params(**pmos_params)
        self.isd_ref_m3, self.vsg_ref_m3 = self.pmos_m3.add(ckt)

        # Add pmos
        nodes = {"g": self.vfp, "d": self.vfn, "s": self.VCC}
        self.pmos_m4 = pmos_subckt(**nodes)
        self.pmos_m4.set_params(**pmos_params)
        self.isd_ref_m4, self.vsg_ref_m4 = self.pmos_m4.add(ckt)

        if self.neg_feedback:
            self.Rf = resistor()
            ckt.add_edge(self.vfn, n2, self.Rf)

            self.Rin = resistor()
            ckt.add_edge(n2, self.vlon, self.Rin)

    def add_small(self, op, ckt_sml, dyn=False, **nodes):

        VCC  = nodes["vcc"]
        GND  = nodes["gnd"]
        vlop = nodes["vlop"]
        vlon = nodes["vlon"]
        vfp  = nodes["vfp"]
        vfn  = nodes["vfn"]

        n1 = ckt_sml.get_internal_node()

        if self.neg_feedback or self.pos_feedback:
            n2 = ckt_sml.get_internal_node()

        # Add nmos
        nodes = {"g": n2 if self.neg_feedback else vlon, "d": vfn, "s": n1}
        self.ids_ref_sml_m2, self.vgs_ref_sml_m2 = self.nmos_m2.add_small(op=op, ckt_sml=ckt_sml, **nodes)

        # Add nmos
        nodes = {"g": vlop, "d": vfp, "s": n1}
        self.ids_ref_sml_m1, self.vgs_ref_sml_m1 = self.nmos_m1.add_small(op=op, ckt_sml=ckt_sml, **nodes)
        
        # Add pmos
        nodes = {"g": vfp, "d": vfp, "s": VCC}
        self.isd_ref_sml_m3, self.vsg_ref_sml_m3 = self.pmos_m3.add_small(op=op, ckt_sml=ckt_sml, **nodes)

        # Add pmos
        nodes = {"g": vfp, "d": vfn, "s": VCC}
        self.isd_ref_sml_m4, self.vsg_ref_sml_m4 = self.pmos_m4.add_small(op=op, ckt_sml=ckt_sml, **nodes)

        if self.neg_feedback:
            self.Rf_sml = resistor()
            ckt_sml.add_edge(vfn, n2, self.Rf_sml)

            self.Rin_sml = resistor()
            ckt_sml.add_edge(n2, vlon, self.Rin_sml)
