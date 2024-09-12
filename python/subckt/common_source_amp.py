import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from subckt import *
from nmos_subckt import *
from pmos_subckt import *

class common_source_amp(subckt):

    def __init__(self, **nodes):

        self.VCC  = nodes["vcc"]
        self.GND  = nodes["gnd"]
        self.vin  = nodes["vin"]
        self.vout = nodes["vout"]

    def set_params(self, **params):

        self.rs_val = params["Rs"]
        self.rd_val = params["Rd"]

        if self.src_degen:
            self.Rs.set_value(self.rs_val)
        self.Rd.set_value(self.rd_val)

    def add(self, ckt, src_degen=False):

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 10,
                       "W"      : 20}

        self.src_degen = src_degen

        if src_degen:
            vdegen     = ckt.get_internal_node()

        self.rd_ref = ckt.num_edges

        # Add Rd
        self.Rd = resistor()
        self.Rd.set_instance("Rd")
        ckt.add_edge(self.VCC, self.vout, self.Rd)

        # Add nmos
        if src_degen:
            nodes = {"g": self.vin, "d": self.vout, "s": vdegen}
        else:
            nodes = {"g": self.vin, "d": self.vout, "s": self.GND}
        self.nmos_m1 = nmos_subckt(**nodes)
        self.nmos_m1.set_params(**nmos_params)
        self.nmos_m1.add(ckt)

        if src_degen:
            # Add Rs
            self.Rs = resistor()
            self.Rs.set_instance("Rs")
            ckt.add_edge(vdegen, self.GND, self.Rs)

    def add_small(self, op, ckt_sml, dyn=False, **nodes):

        pass
