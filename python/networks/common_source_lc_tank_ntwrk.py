import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from common_source_lc_tank import *

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

class common_source_lc_tank_ntwrk():

    def __init__(self, src_degen=False, **params):

        # Nodes
        GND  = 0
        VCC  = 1
        vin  = 2
        vout = 3

        self.ckt = Circuit()

        # Add common source amp
        nodes = {"vcc": VCC, "gnd": GND, "vin": vin, "vout": vout}
        self.common_source_lc_tank = common_source_lc_tank(**nodes)
        self.common_source_lc_tank.add(self.ckt)

        self.common_source_lc_tank.set_params(**params)

        # Add Vin
        self.v_vin = sine_voltage_src(omega=(50.34e3 * math.pi * 2), mag=0.2, phi=0, bias=4.0)
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

    def add_op_ckt(self):

        # Nodes
        GND  = 0
        VCC  = 1
        vin  = 2
        vout = 3

        self.ckt_op = Circuit()

        self.common_source_lc_tank.add_op(self.ckt_op)

        # Add Vin
        self.v_vin_op = voltage_src()
        self.v_vin_op.set_instance("VS")
        self.v_vin_op.set_value(5.0)
        self.ckt_op.add_edge(vin, GND, self.v_vin_op)

        # Add Vs
        self.ckt_op.add_edge(VCC, GND, self.v_vs)

        self.ckt_op.init_circuit()

        self.scb_op = Scoreboard(num_edges=self.ckt_op.num_edges, \
                                 num_sys_vars=self.ckt_op.num_sys_vars, \
                                 degen_mtrx=self.ckt_op.degen_mtrx)


    def add_sml_ckt(self, op):

        # Nodes
        GND  = 0
        VCC  = 0
        vin  = 2
        vout = 3

        self.ckt_sml = Circuit(lti=True)

        # Add common source amp
        nodes = {"vcc": VCC, "gnd": GND, "vin": vin, "vout": vout}

        self.common_source_lc_tank.add_small(op=op, ckt_sml=self.ckt_sml, **nodes)

        # Add Vin
        self.v_vin_sml = sine_voltage_src(omega=(50.34e3 * math.pi * 2), mag=0.2, phi=0, bias=0.0)
        self.v_vin_sml.set_is_input()
        self.ckt_sml.add_edge(vin, GND, self.v_vin_sml)

        self.vmeas = current_src()
        self.vmeas.set_is_const()
        self.vmeas.set_value(0.0)
        self.vmeas.set_is_output()
        self.vmeas.set_instance("VMEAS")
        self.ckt_sml.add_edge(vout, GND, self.vmeas)

        self.ckt_sml.init_circuit()


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

    #def set_vin(self, vin):
     #   self.v_vin.set_value(value=vin)

    def set_vin(self, omega, mag, phi, bias):
        self.v_vin.set_params(omega=omega, mag=mag, phi=phi, bias=bias)

    def set_vin_op(self, vin):
        self.v_vin_op.set_value(value=vin)

    def get_vout(self, x):
        return 0
        self.scb.x = x

        rd_edge = self.ckt.get_edge_info(self.common_source_amp.rd_ref)
        rd_edge.get_voltage(scb=self.scb)
        v_rd = self.scb.v[self.common_source_amp.rd_ref]

        return self.v_vs.get_value() - v_rd

    def check_saturation_region(self, x):

        if self.common_source_amp.nmos_m1.ids.get_region(x) == 0:
            return True
        return False

    def check_triode_region(self, x):

        if self.common_source_amp.nmos_m1.ids.get_region(x) == 1:
            return True
        return False

    def get_ckt_op_op(self, x):

        return self.common_source_lc_tank.nmos_m1_op.ids.get_op_t(x)

    def get_vout_op(self, x):

        self.scb_op = x
        self.scb_op.v = self.ckt_op.get_vm(x=self.scb_op.x, sys=0, t=0)

        rd_edge = self.ckt_op.get_edge_info(self.common_source_lc_tank.R_op_ref)
        rd_edge.get_voltage(scb=self.scb_op)
        v_rd = self.scb_op.v[self.common_source_lc_tank.R_op_ref]

        print (self.common_source_lc_tank.R_op_ref)

        v_rd = x.x[0] * 3e3

        return self.v_vs.get_value() - v_rd
