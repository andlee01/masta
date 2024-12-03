import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from subckt import *
from nmos_subckt import *
from pmos_subckt import *

class common_source_lc_tank(subckt):

    def __init__(self, **nodes):

        self.VCC  = nodes["vcc"]
        self.GND  = nodes["gnd"]
        self.vin  = nodes["vin"]
        self.vout = nodes["vout"]

    def set_params(self, **params):

        self.L_val = params["L"]
        self.C_val = params["C"]
        self.R_val = params["R"]

        self.L.set_value(self.L_val)
        self.C.set_value(self.C_val)
        self.R.set_value(self.R_val)

    def add(self, ckt):

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 10}

        n1     = ckt.get_internal_node()
        n2     = ckt.get_internal_node()

        # Add L
        self.L = inductor()
        self.L.set_instance("L")
        ckt.add_edge(self.VCC, n1, self.L)

        # Add R
        self.R = resistor()
        self.R.set_instance("R")
        ckt.add_edge(n1, self.vout, self.R)

        # Add C
        self.C = capacitor()
        self.C.set_instance("C")
        ckt.add_edge(self.VCC, self.vout, self.C)

        # Add nmos
        nodes = {"g": self.vin, "d": self.vout, "s": n2}
        self.nmos_m1 = nmos_subckt(**nodes)
        self.nmos_m1.set_params(**nmos_params)
        self.nmos_m1.add(ckt)

        # Add Ibias
        self.Ibias = current_src()
        self.Ibias.set_value(100e-6)
        ckt.add_edge(n2, self.GND, self.Ibias)

    def add_op(self, ckt_op):

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 10}

        n2     = ckt_op.get_internal_node()

        # Add R
        self.R_op = resistor()
        self.R_op.set_instance("R")
        self.R_op.set_value(self.R.get_value())
        ckt_op.add_edge(self.VCC, self.vout, self.R_op)

        # Add nmos
        nodes = {"g": self.vin, "d": self.vout, "s": n2}
        self.nmos_m1_op = nmos_subckt(**nodes)
        self.nmos_m1_op.set_params(**nmos_params)
        self.nmos_m1_op.add(ckt_op)

        # Add Ibias
        self.Ibias_op = current_src()
        self.Ibias_op.set_value(100e-6)
        ckt_op.add_edge(n2, self.GND, self.Ibias_op)


    def add_small(self, op, ckt_sml, **nodes):

        VCC  = nodes["vcc"]
        GND  = nodes["gnd"]
        vin  = nodes["vin"]
        vout = nodes["vout"]

        n1     = ckt_sml.get_internal_node()

        # Add L
        self.L_sml = inductor()
        self.L_sml.set_value(self.L.get_value())
        ckt_sml.add_edge(VCC, n1, self.L_sml)

        # Add R
        self.R_sml = resistor()
        self.R_sml.set_value(self.R.get_value())
        ckt_sml.add_edge(n1, vout, self.R_sml)

        # Add C
        self.C_sml = capacitor()
        self.C_sml.set_value(self.C.get_value())
        ckt_sml.add_edge(VCC, vout, self.C_sml)

        # Add nmos
        nodes = {"g": self.vin, "d": self.vout, "s": self.GND}
        self.nmos_m1_op.add_small(op=op, ckt_sml=ckt_sml, **nodes)
