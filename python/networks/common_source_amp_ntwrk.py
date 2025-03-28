import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from common_source_amp import *

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

class common_source_amp_ntwrk():

    def __init__(self, src_degen=False, **params):

        # Nodes
        GND  = 0
        VCC  = 1
        vin  = 2
        vout = 3

        self.ckt = Circuit()

        # Add common source amp
        nodes = {"vcc": VCC, "gnd": GND, "vin": vin, "vout": vout}
        self.common_source_amp = common_source_amp(**nodes)
        self.common_source_amp.add(self.ckt, src_degen=src_degen)

        self.common_source_amp.set_params(**params)

        # Add Vin
        self.v_vin = voltage_src()
        self.v_vin.set_instance("VS")
        self.ckt.add_edge(vin, GND, self.v_vin)

        # Add Vs
        self.v_vs = voltage_src()
        self.v_vs.set_instance("VS")
        self.v_vs.set_value(5.0)
        self.ckt.add_edge(VCC, GND, self.v_vs)

        self.ckt.init_circuit()

        self.scb = Scoreboard(num_edges=self.ckt.num_edges, \
                              num_sys_vars=self.ckt.num_sys_vars, \
                              degen_mtrx=self.ckt.degen_mtrx)

    def add_sml_ckt(self, op):
        # Nodes
        GND  = 0
        vin  = 1
        vout = 2

        self.ckt_sml = Circuit()

        nodes = {"gnd": GND, "vin": vin, "vout": vout}
        self.common_source_amp.add_small(ckt_sml=self.ckt_sml, op=op, **nodes)

        # Add Vs
        self.v_vin_sml = voltage_src()
        self.v_vin_sml.set_instance("VS")
        self.ckt_sml.add_edge(vin, GND, self.v_vin_sml)

        self.ckt_sml.init_circuit()

    def set_vin(self, vin):
        self.v_vin.set_value(value=vin)

    def get_vout(self, x):

        self.scb.x = x

        rd_edge = self.ckt.get_edge_info(self.common_source_amp.rd_ref)
        rd_edge.get_voltage(scb=self.scb)
        v_rd = self.scb.v[self.common_source_amp.rd_ref]

        return self.v_vs.get_value() - v_rd

    def get_sml_vout(self, x):

        self.scb.x = x

        rd_edge = self.ckt_sml.get_edge_info(self.common_source_amp.rd_sml_ref)
        rd_edge.get_voltage(scb=self.scb)
        v_rd = self.scb.v[self.common_source_amp.rd_sml_ref]

        return -v_rd

    def check_saturation_region(self, x):

        if self.common_source_amp.nmos_m1.ids.get_region(x) == 0:
            return True
        return False

    def check_triode_region(self, x):

        if self.common_source_amp.nmos_m1.ids.get_region(x) == 1:
            return True
        return False

    def get_op(self, x):

        return self.common_source_amp.nmos_m1.ids.get_op(x)

    def get_op_sml(self):

        return self.common_source_amp.nmos_m1.ids_sml.get_op()
