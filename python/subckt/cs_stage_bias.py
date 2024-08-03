import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from subckt import *
from nmos_subckt import *
from pmos_subckt import *

class cs_stage_bias(subckt):

    def __init__(self, **nodes):

        self.VCC  = nodes["vcc"]
        self.GND  = nodes["gnd"]
        self.vin  = nodes["vin"]
        self.vx   = nodes["vx"]

    def set_params(self, **params):

        self.Cb_val = params["Cb"]
        self.Rb_val = params["Rb"]
        self.Ib_val = params["Ib"]

        self.ib.set_value(self.Ib_val)
        self.Cb.set_value(self.Cb_val)
        self.Rb.set_value(self.Rb_val)

    def add(self, ckt):

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 10}

        vb     = ckt.get_internal_node()

        self.nmos_vds_ref = ckt.num_edges

        # Add nmos
        nodes = {"g": vb, "d": vb, "s": self.GND}
        self.nmos_m2 = nmos_subckt(**nodes)
        self.nmos_m2.set_params(**nmos_params)
        self.nmos_m2.add(ckt)

        # Add Ib
        self.ib = current_src()
        self.ib.set_instance("iref_n")
        ckt.add_edge(self.VCC, vb, self.ib)

        self.Rb_ref = ckt.num_edges

        # Add Rb
        self.Rb = resistor()
        self.Rb.set_instance("Rb")
        ckt.add_edge(self.vx, vb, self.Rb)

        # Add Cb
        self.Cb = capacitor()
        self.Cb.set_instance("Cb")
        ckt.add_edge(self.vin, self.vx, self.Cb)

    def add_small(self, op, ckt_sml, dyn=False, **nodes):

        pass
