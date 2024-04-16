import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from nmos_subckt import *

class sense_amp_cross_coupled_nmos():

    def __init__(self, ckt, trans=False):

        # Nodes
        GND   = 0
        Vcol0 = 1
        Vcol1 = 2

        self.nmos_params = {"KP"     : 120e-6,
                            "vth"    : 0.8,
                            "lambda" : 0.01,
                            "L"      : 2,
                            "W"      : 10}

        # Add nmos
        nodes = {"g": Vcol1, "d": Vcol0, "s": GND}
        self.nmos_m0 = nmos_subckt(**nodes)
        self.nmos_m0.set_params(**self.nmos_params)
        self.nmos_m0.add(ckt)

        # Add nmos
        nodes = {"g": Vcol0, "d": Vcol1, "s": GND}
        self.nmos_m1 = nmos_subckt(**nodes)
        self.nmos_m1.set_params(**self.nmos_params)
        self.nmos_m1.add(ckt)

        if not trans:
            # Add Vcol0
            Vc0 = voltage_src()
            Vc0.set_instance("Vcol0")
            Vc0.set_value(value=2.5)
            ckt.add_edge(Vcol0, GND, Vc0)

            # Add Vcol1
            Vc1 = voltage_src()
            Vc1.set_instance("Vcol1")
            Vc1.set_value(value=2.5)
            ckt.add_edge(Vcol1, GND, Vc1)

        if trans:
            # Add Vcol0
            Vc0 = capacitor()
            Vc0.set_instance("Vcol0")
            Vc0.set_value(value=100e-15)
            ckt.add_edge(Vcol0, GND, Vc0)

            # Add Vcol1
            Vc1 = capacitor()
            Vc1.set_instance("Vcol1")
            Vc1.set_value(value=100e-15)
            ckt.add_edge(Vcol1, GND, Vc1)

    def print_op_regions(self, x):

        print (self.nmos_m0.ids.get_region(x))
        print (self.nmos_m1.ids.get_region(x))

    def add_sml_ckt(self, ckt_sml, op, dyn=False):

        # Nodes
        GND     = 0
        Vcol0   = 1
        Vcol1   = 2
        Vcol0_R = 3
        Vcol1_R = 4

        # Add nmos
        nodes = {"g": Vcol1, "d": Vcol0, "s": GND}
        self.nmos_m0.add_small(op=op, ckt_sml=ckt_sml, dyn=dyn, **nodes)

        # Add nmos
        nodes = {"g": Vcol0, "d": Vcol1, "s": GND}
        self.nmos_m1.add_small(op=op, ckt_sml=ckt_sml, dyn=dyn, **nodes)

        # Add Vcol0
        Vc0 = capacitor()
        Vc0.set_instance("Vcol0")
        Vc0.set_value(value=100e-15)
        ckt_sml.add_edge(Vcol0, GND, Vc0)

        # Add Vcol1
        Vc1 = capacitor()
        Vc1.set_instance("Vcol1")
        Vc1.set_value(value=100e-15)
        ckt_sml.add_edge(Vcol1, GND, Vc1)
