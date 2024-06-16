import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from subckt import *
from nmos_subckt import *

class gilbert_cell(subckt):

    def __init__(self, **nodes):

        self.VCC  = nodes["vcc"]
        self.GND  = nodes["gnd"]
        self.vlop = nodes["vlop"]
        self.vlon = nodes["vlon"]
        self.vrfp = nodes["vrfp"]
        self.vrfn = nodes["vrfn"]
        self.vfp  = nodes["vfp"]
        self.vfn  = nodes["vfn"]

    def set_params(self, **params):
        pass

    def add(self, ckt):

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 10}

        n1 = ckt.get_internal_node()
        n2 = ckt.get_internal_node()
        n3 = ckt.get_internal_node()

        # Add nmos
        nodes = {"g": self.vlop, "d": self.vfp, "s": n1}
        self.nmos_m4 = nmos_subckt(**nodes)
        self.nmos_m4.set_params(**nmos_params)
        self.nmos_m4.add(ckt)

        # Add nmos
        nodes = {"g": self.vlon, "d": self.vfn, "s": n1}
        self.nmos_m5 = nmos_subckt(**nodes)
        self.nmos_m5.set_params(**nmos_params)
        self.nmos_m5.add(ckt)

        # Add nmos
        nodes = {"g": self.vlon, "d": self.vfp, "s": n2}
        self.nmos_m6 = nmos_subckt(**nodes)
        self.nmos_m6.set_params(**nmos_params)
        self.nmos_m6.add(ckt)

        # Add nmos
        nodes = {"g": self.vlop, "d": self.vfn, "s": n2}
        self.nmos_m7 = nmos_subckt(**nodes)
        self.nmos_m7.set_params(**nmos_params)
        self.nmos_m7.add(ckt)

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 20}

        # Add nmos
        nodes = {"g": self.vrfp, "d": n1, "s": n3}
        self.nmos_m2 = nmos_subckt(**nodes)
        self.nmos_m2.set_params(**nmos_params)
        self.nmos_m2.add(ckt)

        # Add nmos
        nodes = {"g": self.vrfn, "d": n2, "s": n3}
        self.nmos_m3 = nmos_subckt(**nodes)
        self.nmos_m3.set_params(**nmos_params)
        self.nmos_m3.add(ckt)

        # Add Iref
        self.iref = current_src()
        self.iref.set_instance("iref")
        ckt.add_edge(n3, self.GND, self.iref)

        # Add Rd1
        Rd1 = resistor()
        Rd1.set_instance("Rd1")
        Rd1.set_value(value=10e3)
        self.rd1_ref = ckt.add_edge(self.VCC, self.vfp, Rd1)

        # Add Rd2
        Rd2 = resistor()
        Rd2.set_instance("Rd2")
        Rd2.set_value(value=10e3)
        self.rd2_ref = ckt.add_edge(self.VCC, self.vfn, Rd2)

    def add_small(self, op, ckt_sml, dyn=False, **nodes):

        pass
