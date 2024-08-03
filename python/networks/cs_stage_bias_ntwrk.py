import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from cs_stage_bias import *

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

class cs_stage_bias_ntwrk():

    def __init__(self, **params):

        # Nodes
        GND  = 0
        VCC  = 1
        vin  = 2
        vx   = 3

        self.ckt = Circuit()

        self.omega = params["omega"]
        self.mag   = params["mag"]

        # Add cs_stage_bias
        nodes = {"vcc": VCC, "gnd": GND, "vin": vin, "vx": vx}
        self.cs_stage_bias = cs_stage_bias(**nodes)
        self.cs_stage_bias.add(self.ckt)

        params = {"Cb": 1e-12, "Rb":1e3, "Ib":20e-6}
        self.cs_stage_bias.set_params(**params)

        # Add Vin
        self.v_vin = sine_voltage_src(omega=self.omega, mag=self.mag, phi=0, bias=0)
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

    #def add_sml_ckt(self, op, sine_src=False, **params):
    #
    #    pass
    #
    #def set_iref(self, iref):
    #    return
    #    self.diff_amp.iref_n.set_value(value=iref)
    #    #self.diff_amp.iref_out.set_value(value=iref)
    #
    #def set_iref_sml(self, iref):
    #    self.diff_amp.iref_n_sml.set_value(value=iref)
    #
    #def set_vlo(self, vlo):
    #    #self.v_vlon.set_value(value=self.v_lo_bias - vlo)
    #    self.v_vlop.set_value(value=self.v_lo_bias + vlo)
    #
    #def set_vlo_sml(self, vlo):
    #    self.v_vlon_sml.set_value(value=-vlo)
    #    self.v_vlop_sml.set_value(value= vlo)
    #
    #def set_vlo_bias(self, bias):
    #    self.v_vlon.set_value(value=bias)
    #    self.v_vlop.set_value(value=bias)
    #
    #    self.v_lo_bias = bias
    #
    #def check_saturation_region(self, x):
    #
    #    if self.diff_amp.nmos_m1.ids.get_region(x) == 0 and self.diff_amp.nmos_m2.ids.get_region(x) == 0:
    #        return True
    #    return False
    #
    #def get_vf_output(self, x):
    #
    #    self.scb.x = x
    #
    #    rd1_edge = self.ckt.get_edge_info(self.diff_amp.pmos_m4_vsd_ref)
    #    rd1_edge.get_voltage(scb=self.scb)
    #    v_rd1 = self.scb.v[self.diff_amp.pmos_m4_vsd_ref]
    #
    #    return (5.0 - v_rd1)
    #
    #def get_vf_output_base(self, x):
    #
    #    self.scb.x = x
    #
    #    rd1_edge = self.ckt.get_edge_info(self.diff_amp.rd1_ref)
    #    rd1_edge.get_voltage(scb=self.scb)
    #    v_rd1 = self.scb.v[self.diff_amp.rd1_ref]
    #
    #    rd2_edge = self.ckt.get_edge_info(self.diff_amp.rd2_ref)
    #    rd2_edge.get_voltage(scb=self.scb)
    #    v_rd2 = self.scb.v[self.diff_amp.rd2_ref]
    #
    #    return v_rd1, v_rd2
    #
    #def get_vf_output_sml(self, ckt, scb, v_rd1_base, v_rd2_base):
    #
    #    rd1_edge = ckt.get_edge_info(self.rd1_ref_sml)
    #    rd1_edge.get_voltage(scb=scb)
    #    v_rd1 = scb.v[self.rd1_ref_sml]
    #
    #    rd2_edge = ckt.get_edge_info(self.rd2_ref_sml)
    #    rd2_edge.get_voltage(scb=scb)
    #    v_rd2 = scb.v[self.rd2_ref_sml]
    #
    #    return (5.0 - (v_rd1 + v_rd1_base)) - (5.0 - (v_rd2 + v_rd2_base))

    def get_vx(self, sys, t):
        return self.v_vin.get_value() - sys[0]
