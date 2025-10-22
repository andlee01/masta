import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from subckt import *
from nmos_subckt import *
from pmos_subckt import *

class beta_multiplier(subckt):

    def __init__(self, **nodes):

        self.VCC    = nodes["vcc"]
        self.GND    = nodes["gnd"]
        self.vbiasp = nodes["vbiasp"]
        self.vbiasn = nodes["vbiasn"]


    def set_params(self, **params):

        self.Rref_val   = params["Rref"]

        self.Rref.set_value(self.Rref_val)

    def set_sml_params(self, **params):

        self.Rref_sml_val   = params["Rref"]

        self.Rref_sml.set_value(self.Rref_sml_val)

    def add(self, ckt, output_resistance=True):

        n1     = ckt.get_internal_node()

        if output_resistance:
            lambda_nmos = 0.0125
            lambda_pmos = 0.01
        else:
            lambda_nmos = 0
            lambda_pmos = 0

        # pmos current mirror
        # -------------------

        m3_params = {"KP"     : 40e-6,
                     "vth"    : 0.9,
                     "lambda" : lambda_pmos,
                     "L"      : 2,
                     "W"      : 30}

        nodes = {"g": self.vbiasp, "d": self.vbiasn, "s": self.VCC}
        self.pmos_m3 = pmos_subckt(**nodes)
        self.pmos_m3.set_params(**m3_params)
        self.pmos_m3.add(ckt)

        m4_params = {"KP"     : 40e-6,
                     "vth"    : 0.9,
                     "lambda" : lambda_pmos,
                     "L"      : 2,
                     "W"      : 30}

        nodes = {"g": self.vbiasp, "d": self.vbiasp, "s": self.VCC}
        self.pmos_m4 = pmos_subckt(**nodes)
        self.pmos_m4.set_params(**m4_params)
        self.pmos_m4.add(ckt)

        # nmos current mirror
        # -------------------

        m1_params = {"KP"     : 120e-6,
                     "vth"    : 0.8,
                     "lambda" : lambda_nmos,
                     "L"      : 2,
                     "W"      : 10}

        nodes = {"g": self.vbiasn, "d": self.vbiasn, "s": self.GND}
        self.nmos_m1 = nmos_subckt(**nodes)
        self.nmos_m1.set_params(**m1_params)
        self.nmos_m1.add(ckt)

        m2_params = {"KP"     : 120e-6,
                     "vth"    : 0.8,
                     "lambda" : lambda_nmos,
                     "L"      : 2,
                     "W"      : 40}

        nodes = {"g": self.vbiasn, "d": self.vbiasp, "s": n1}
        self.nmos_m2 = nmos_subckt(**nodes)
        self.nmos_m2.set_params(**m2_params)
        self.nmos_m2.add(ckt)

        # Add Rn
        self.Rref_idx = ckt.num_edges
        self.Rref = resistor()
        ckt.add_edge(n1, self.GND, self.Rref)

    def add_small(self, op, ckt_sml, output_resistance=True, **nodes):

        n1     = ckt_sml.get_internal_node()

        # pmos current mirror
        # -------------------

        nodes = {"g": self.vbiasp, "d": self.vbiasn, "s": self.VCC}
        self.isd_ref_sml_m3, self.vsg_ref_sml_m3 = self.pmos_m3.add_small(op=op, ckt_sml=ckt_sml, output_resistance=output_resistance, **nodes)

        nodes = {"g": self.vbiasp, "d": self.vbiasp, "s": self.VCC}
        self.isd_ref_sml_m4, self.vsg_ref_sml_m4 = self.pmos_m4.add_small(op=op, ckt_sml=ckt_sml, output_resistance=output_resistance, **nodes)

        # nmos current mirror
        # -------------------

        nodes = {"g": self.vbiasn, "d": self.vbiasn, "s": self.GND}
        self.ids_ref_sml_m1, self.vgs_ref_sml_m1 = self.nmos_m1.add_small(op=op, ckt_sml=ckt_sml, output_resistance=output_resistance, **nodes)

        nodes = {"g": self.vbiasn, "d": self.vbiasp, "s": n1}
        self.ids_ref_sml_m2, self.vgs_ref_sml_m2 = self.nmos_m2.add_small(op=op, ckt_sml=ckt_sml, output_resistance=output_resistance, **nodes)

        # Add Rn
        self.Rref_sml_idx = ckt_sml.num_edges
        self.Rref_sml = resistor()
        ckt_sml.add_edge(n1, self.GND, self.Rref_sml)
