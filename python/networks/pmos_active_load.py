import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from pmos_subckt import *

class pmos_active_load():

    def __init__(self, ckt):

        # Nodes
        GND = 0
        VCC = 1
        n1  = 2
        v1  = 3
        n2  = 4

        pmos_params = {"KP"     : 40e-6,
                       "vth"    : 0.9,
                       "lambda" : 0.0125,
                       "L"      : 2,
                       "W"      : 30}

        # Add pmos
        nodes = {"g": n2, "d": n2, "s": VCC}
        self.pmos_m1 = pmos_subckt(**nodes)
        self.pmos_m1.set_params(**pmos_params)
        self.pmos_m1.add(ckt)

        # Add Iref
        iref_p = current_src()
        iref_p.set_type(ElementType.current_src)
        iref_p.set_instance("iref_p")
        iref_p.set_value(value=20e-6)
        ckt.add_edge(n2, GND, iref_p)

        # Add pmos
        nodes = {"g": n2, "d": n1, "s": VCC}
        self.pmos_m2 = pmos_subckt(**nodes, **pmos_params)
        self.pmos_m2.set_params(**pmos_params)
        self.pmos_m2.add(ckt)

        # Add Iref
        iref_n = current_src()
        iref_n.set_type(ElementType.current_src)
        iref_n.set_instance("iref_n")
        iref_n.set_value(value=20e-6)
        ckt.add_edge(n1, GND, iref_n)

        # Add Vs
        vs = voltage_src()
        vs.set_type(ElementType.voltage_src)
        vs.set_instance("VS")
        vs.set_value(value=5.0)
        ckt.add_edge(VCC, GND, vs)

    def add_sml_ckt(self, ckt_sml, op):

        # Nodes
        GND = 0
        VCC = 1
        n1  = 2
        v1  = 3
        n2  = 4

        pmos_params = {"KP"     : 40e-6,
                       "vth"    : 0.9,
                       "lambda" : 0.0125,
                       "L"      : 2,
                       "W"      : 30}

        # Add pmos
        nodes = {"g": n2, "d": n2, "s": GND}
        self.pmos_m1.add_small(op=op, ckt_sml=ckt_sml, **nodes)

        # Add Iref
        iref_p = current_src()
        iref_p.set_type(ElementType.current_src)
        iref_p.set_instance("iref_p")
        iref_p.set_value(value=0)
        ckt_sml.add_edge(n2, GND, iref_p)

        # Add pmos
        nodes = {"g": n2, "d": n1, "s": GND}
        self.pmos_m2.add_small(op=op, ckt_sml=ckt_sml, **nodes)

        # Add Iref
        iref_n = current_src()
        iref_n.set_type(ElementType.current_src)
        iref_n.set_instance("iref_n")
        iref_n.set_value(value=0)
        ckt_sml.add_edge(n1, GND, iref_n)
