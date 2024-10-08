import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from subckt import *
from nmos_subckt import *
from pmos_subckt import *

class diff_amp_pmos_load(subckt):

    def __init__(self, **nodes):

        self.VCC  = nodes["vcc"]
        self.GND  = nodes["gnd"]
        self.vlop = nodes["vlop"]
        self.vlon = nodes["vlon"]
        self.vfp  = nodes["vfp"]
        self.vfn  = nodes["vfn"]
        self.vout = nodes["vout"]

    def set_params(self, **params):
        pass

    def add(self, ckt):

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 10}

        nmos_params_out = {"KP"     : 120e-6,
                           "vth"    : 0.8,
                           "lambda" : 0.01,
                           "L"      : 2,
                           "W"      : 20}

        pmos_params = {"KP"     : 40e-6,
                       "vth"    : 0.9,
                       "lambda" : 0.0125,
                       "L"      : 2,
                       "W"      : 30}

        n1     = ckt.get_internal_node()
        vbias3 = ckt.get_internal_node()

        # Add nmos
        nodes = {"g": self.vlop, "d": self.vfp, "s": n1}
        self.nmos_m2 = nmos_subckt(**nodes)
        self.nmos_m2.set_params(**nmos_params)
        self.nmos_m2.add(ckt)

        # Add nmos
        nodes = {"g": self.vlon, "d": self.vfn, "s": n1}
        self.nmos_m1 = nmos_subckt(**nodes)
        self.nmos_m1.set_params(**nmos_params)
        self.nmos_m1.add(ckt)

        # Add Iref
        #self.iref_n = current_src()
        #self.iref_n.set_instance("iref_n")
        #ckt.add_edge(n1, self.GND, self.iref_n)

        nodes = {"g": vbias3, "d": n1, "s": self.GND}
        self.nmos_m6t = nmos_subckt(**nodes)
        self.nmos_m6t.set_params(**nmos_params_out)
        self.nmos_m6t.add(ckt)

        # Add pmos
        nodes = {"g": self.vfn, "d": self.vfn, "s": self.VCC}
        self.pmos_m3 = pmos_subckt(**nodes)
        self.pmos_m3.set_params(**pmos_params)
        self.pmos_m3.add(ckt)

        # Add pmos
        nodes = {"g": self.vfn, "d": self.vfp, "s": self.VCC}
        self.pmos_m4 = pmos_subckt(**nodes)
        self.pmos_m4.set_params(**pmos_params)
        self.pmos_m4.add(ckt)

        self.pmos_m4_vsd_ref = ckt.num_edges

        # Add pmos
        nodes = {"g": self.vfp, "d": self.vout, "s": self.VCC}
        self.pmos_m7 = pmos_subckt(**nodes)
        self.pmos_m7.set_params(**pmos_params)
        self.pmos_m7.add(ckt)

        # Add nmos
        nodes = {"g": vbias3, "d": self.vout, "s": self.GND}
        self.nmos_m8t = nmos_subckt(**nodes)
        self.nmos_m8t.set_params(**nmos_params)
        self.nmos_m8t.add(ckt)

        self.v_bias3 = voltage_src()
        self.v_bias3.set_value(value=1.5)
        ckt.add_edge(vbias3, self.GND, self.v_bias3)

        ## Add nmos
        #nodes = {"g": self.vlon, "d": self.vfn, "s": n1}
        #self.nmos_m1 = nmos_subckt(**nodes)
        #self.nmos_m1.set_params(**nmos_params)
        #self.nmos_m1.add(ckt)


        ## Add Iref
        #self.iref_out = current_src()
        #self.iref_out.set_instance("iref_n")
        #ckt.add_edge(vout, self.GND, self.iref_out)

    def add_small(self, op, ckt_sml, dyn=False, **nodes):

        VCC  = nodes["vcc"] # revisit - should use vcc in here
        GND  = nodes["gnd"]
        vlop = nodes["vlop"]
        vlon = nodes["vlon"]
        vfp  = nodes["vfp"]
        vfn  = nodes["vfn"]

        n1 = ckt_sml.get_internal_node()

        # Add nmos
        nodes = {"g": vlop, "d": vfp, "s": n1}
        self.nmos_m1.add_small(op=op, ckt_sml=ckt_sml, **nodes)

        # Add nmos
        nodes = {"g": vlon, "d": vfn, "s": n1}
        self.nmos_m2.add_small(op=op, ckt_sml=ckt_sml, **nodes)

        # Add Iref
        self.iref_n_sml = current_src()
        self.iref_n_sml.set_instance("iref_n")
        self.iref_n_sml.set_value(value=0)
        ckt_sml.add_edge(n1, GND, self.iref_n_sml)

        # Add Rd1
        Rd1 = resistor()
        Rd1.set_instance("Rd1")
        Rd1.set_value(value=1e4)
        self.rd1_ref_sml = ckt_sml.add_edge(GND, vfp, Rd1)

        # Add Rd2
        Rd2 = resistor()
        Rd2.set_instance("Rd2")
        Rd2.set_value(value=1e4)
        self.rd2_ref_sml = ckt_sml.add_edge(GND, vfn, Rd2)
