import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from beta_multiplier import *

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

class beta_multiplier_ntwrk():

    def __init__(self, output_resistance=True, res_m2=True, stray_capacitance=False):

        # Nodes
        GND    = 0
        VCC    = 1
        vbiasn = 2
        vbiasp = 3

        self.ckt = Circuit()

        # Add beta multiplier
        nodes = {"vcc": VCC, "gnd": GND, "vbiasp": vbiasp, "vbiasn": vbiasn}
        self.beta_mult = beta_multiplier(instance = "B_{1}", **nodes)
        self.beta_mult.add(ckt=self.ckt, output_resistance=output_resistance, res_m2=res_m2, stray_capacitance=stray_capacitance)

        # Add Vs
        self.v_vs = voltage_src()
        self.v_vs.set_instance("V_{DD}")
        self.v_vs.set_value(5.0)
        self.ckt.add_edge(VCC, GND, self.v_vs)

        self.ckt.init_circuit()

        self.scb = Scoreboard(num_edges=self.ckt.num_edges, \
                              num_sys_vars=self.ckt.num_sys_vars, \
                              degen_mtrx=self.ckt.degen_mtrx)

    def add_sml_ckt(self, op, lti=False, gate_topology=None, degen_topology=None, output_topology=None):

        # Nodes
        GND    = 0
        VCC    = 0 #if lti else 1
        vbiasn = 2
        vbiasp = 3

        self.ckt_sml = Circuit(lti=lti)

        # Add beta multiplier
        nodes = {"vcc": VCC, "gnd": GND, "vbiasp": vbiasp, "vbiasn": vbiasn}
        self.beta_mult.add_small(op=op, \
                                 ckt_sml=self.ckt_sml, \
                                    gate_topology=gate_topology, \
                                        degen_topology=degen_topology, \
                                            output_topology=output_topology, \
                                                **nodes)

        # Add Vs
        #self.v_vs_sml = voltage_src()
        #self.v_vs_sml.set_instance("V_S")
        #self.v_vs_sml.set_value(0.0)
        #self.ckt_sml.add_edge(VCC, GND, self.v_vs_sml)

        self.ckt_sml.init_circuit()

    def gate_break_m1_m2(self):
        return self.beta_mult._broken_gate_m1_m2
    
    def gate_break_m1_m2_closed(self):
        return self.beta_mult._broken_gate_m1_m2_closed

    def _capacitive_degen_topology(self):
        return self.beta_mult._capacitive_degen_topology
    
    def _output_vds_m2_topology(self):
        return self.beta_mult._output_vds_m2_topology

    def get_Rref_current(self, x):
        return x[self.beta_mult.Rref_idx]
    
    def get_rref_idx(self):
        return self.beta_mult.Rref_idx

    def get_Rref_sml_current(self, x):
        return x[self.beta_mult.Rref_sml_idx]

    def get_mos_vltg_sml(self, x):
        return x[self.beta_mult.isd_ref_sml_m3], x[self.beta_mult.vsg_ref_sml_m3], \
               x[self.beta_mult.isd_ref_sml_m4], x[self.beta_mult.vsg_ref_sml_m4], \
               x[self.beta_mult.ids_ref_sml_m1], x[self.beta_mult.vgs_ref_sml_m1], \
               x[self.beta_mult.ids_ref_sml_m2], x[self.beta_mult.vgs_ref_sml_m2]

    def get_mos_vltg(self, x):
        return x[self.beta_mult.isd_ref_m3], x[self.beta_mult.vsg_ref_m3], \
               x[self.beta_mult.isd_ref_m4], x[self.beta_mult.vsg_ref_m4], \
               x[self.beta_mult.ids_ref_m1], x[self.beta_mult.vgs_ref_m1], \
               x[self.beta_mult.ids_ref_m2], x[self.beta_mult.vgs_ref_m2]

    def set_beta_multiplier_params(self, **params):
        self.beta_mult.set_params(**params)
        self.v_vs.set_value(params["Vs"])

    def set_beta_multiplier_sml_params(self, **params):
        self.beta_mult.set_sml_params(**params)
        #self.v_vs_sml.set_value(params["Vs"])

    def set_source(self, source, val):
        self.beta_mult.set_source(name=source, value=val)

    def get_source_idx(self, source):
        return self.beta_mult.get_source_idx(name=source)