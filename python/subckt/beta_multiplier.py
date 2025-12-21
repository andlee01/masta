import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from subckt import *
from nmos_subckt import *
from pmos_subckt import *

class beta_multiplier(subckt):

    def __init__(self, instance: str = "B_{1}", **nodes):

        self.VCC    = nodes["vcc"]
        self.GND    = nodes["gnd"]
        self.vbiasp = nodes["vbiasp"]
        self.vbiasn = nodes["vbiasn"]

        self.instance = instance

    def set_params(self, **params):

        self.Rref_val   = params["Rref"]

        self.Rref.set_value(self.Rref_val)

    def set_sml_params(self, **params):

        self.Rref_sml_val   = params["Rref"]

        self.Rref_sml.set_value(self.Rref_sml_val)

    def add(self, ckt, output_resistance=True, res_m2=True):

        n1     = ckt.get_internal_node()
        n2     = ckt.get_internal_node()

        if output_resistance:
            lambda_nmos = 0.0125
            lambda_pmos = 0.01
        else:
            lambda_nmos = 0
            lambda_pmos = 0

        self.res_m2 = res_m2
        self.output_resistance = output_resistance

        # pmos current mirror
        # -------------------

        m3_params = {"KP"     : 40e-6,
                     "vth"    : 0.9,
                     "lambda" : lambda_pmos,
                     "L"      : 2,
                     "W"      : 30}

        nodes = {"g": self.vbiasp, "d": self.vbiasn, "s": self.VCC}
        self.pmos_m3 = pmos_subckt(instance = self.instance + "M_3", **nodes)
        self.pmos_m3.set_params(**m3_params)
        self.isd_ref_m3, self.vsg_ref_m3 = self.pmos_m3.add(ckt)

        m4_params = {"KP"     : 40e-6,
                     "vth"    : 0.9,
                     "lambda" : lambda_pmos,
                     "L"      : 2,
                     "W"      : 30}

        nodes = {"g": self.vbiasp, "d": self.vbiasp, "s": self.VCC}
        self.pmos_m4 = pmos_subckt(instance = self.instance + "M_4", **nodes)
        self.pmos_m4.set_params(**m4_params)
        self.isd_ref_m4, self.vsg_ref_m4 = self.pmos_m4.add(ckt)

        # nmos current mirror
        # -------------------

        m1_params = {"KP"     : 120e-6,
                     "vth"    : 0.8,
                     "lambda" : lambda_nmos,
                     "L"      : 2,
                     "W"      : 10 if res_m2 else 40}

        nodes = {"g": self.vbiasn, "d": self.vbiasn, "s": self.GND if res_m2 else n1}
        self.nmos_m1 = nmos_subckt(instance = self.instance + "M_1", **nodes)
        self.nmos_m1.set_params(**m1_params)
        self.ids_ref_m1, self.vgs_ref_m1 = self.nmos_m1.add(ckt)

        m2_params = {"KP"     : 120e-6,
                     "vth"    : 0.8,
                     "lambda" : lambda_nmos,
                     "L"      : 2,
                     "W"      : 40 if res_m2 else 10}

        nodes = {"g": self.vbiasn, "d": self.vbiasp, "s": n1 if res_m2 else self.GND}
        self.nmos_m2 = nmos_subckt(instance = self.instance + "M_2", **nodes)
        self.nmos_m2.set_params(**m2_params)
        self.ids_ref_m2, self.vgs_ref_m2 = self.nmos_m2.add(ckt)

        # Add Rn
        self.Rref_idx = ckt.num_edges
        self.Rref = resistor()
        self.Rref.set_instance(self.instance + "R_{ref}")
        ckt.add_edge(n1, self.GND, self.Rref)

        # Start-up
        return
        msu3_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.0125,
                       "L"      : 1,
                       "W"      : 10}

        nodes = {"g": n2, "d": self.vbiasp, "s": self.vbiasn}
        self.nmos_msu3 = nmos_subckt(**nodes)
        self.nmos_msu3.set_params(**msu3_params)
        self.ids_ref_msu3, self.vgs_ref_msu3 = self.nmos_msu3.add(ckt)

        msu1_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.0125,
                       "L"      : 2,
                       "W"      : 10}

        nodes = {"g": self.vbiasn, "d": n2, "s": self.GND}
        self.nmos_msu1 = nmos_subckt(**nodes)
        self.nmos_msu1.set_params(**msu1_params)
        self.ids_ref_msu1, self.vgs_ref_msu1 = self.nmos_msu1.add(ckt)

        msu2_params = {"KP"     : 40e-6,
                       "vth"    : 0.9,
                       "lambda" : 0.01,
                       "L"      : 100,
                       "W"      : 10}

        nodes = {"g": n2, "d": n2, "s": self.VCC}
        self.pmos_msu2 = pmos_subckt(**nodes)
        self.pmos_msu2.set_params(**msu2_params)
        self.isd_ref_msu2, self.vsg_ref_msu2 = self.pmos_msu2.add(ckt)

    def add_small(self, op, ckt_sml, **nodes):

        n1     = ckt_sml.get_internal_node()

        # pmos current mirror
        # -------------------

        nodes = {"g": self.vbiasp, "d": self.vbiasn, "s": self.VCC}
        self.isd_ref_sml_m3, self.vsg_ref_sml_m3 = self.pmos_m3.add_small(op=op, ckt_sml=ckt_sml, output_resistance=self.output_resistance, **nodes)

        nodes = {"g": self.vbiasp, "d": self.vbiasp, "s": self.VCC}
        self.isd_ref_sml_m4, self.vsg_ref_sml_m4 = self.pmos_m4.add_small(op=op, ckt_sml=ckt_sml, output_resistance=self.output_resistance, **nodes)

        # nmos current mirror
        # -------------------

        nodes = {"g": self.vbiasn, "d": self.vbiasn, "s": self.GND if self.res_m2 else n1}
        self.ids_ref_sml_m1, self.vgs_ref_sml_m1 = self.nmos_m1.add_small(op=op, ckt_sml=ckt_sml, output_resistance=self.output_resistance, **nodes)

        nodes = {"g": self.vbiasn, "d": self.vbiasp, "s": n1 if self.res_m2 else self.GND}
        self.ids_ref_sml_m2, self.vgs_ref_sml_m2 = self.nmos_m2.add_small(op=op, ckt_sml=ckt_sml, output_resistance=self.output_resistance, **nodes)

        # Add Rn
        self.Rref_sml_idx = ckt_sml.num_edges
        self.Rref_sml = resistor()
        self.Rref_sml.set_instance("R_{{ref}_{" + self.instance + "}}")
        ckt_sml.add_edge(n1, self.GND, self.Rref_sml)
