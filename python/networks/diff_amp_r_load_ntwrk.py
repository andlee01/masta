import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from diff_amp_r_load import *

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

class diff_amp_r_load_ntwrk():

    def __init__(self, sine_src=False, two_sine_src=False, **params):

        # Nodes
        GND  = 0
        VCC  = 1
        n1   = 2
        vlop = 3
        vlon = 4
        vfp  = 5
        vfn  = 6

        self.ckt = Circuit()

        if sine_src:
            self.omega = params["omega"]
            self.mag   = params["mag"]
            self.bias  = params["bias"]

        if two_sine_src:
            self.omega_0 = params["omega_0"]
            self.omega_1 = params["omega_1"]
            self.mag     = params["mag"]
            self.bias    = params["bias"]

        # Add diff_amp
        nodes = {"vcc": VCC, "gnd": GND, "vlop": vlop, "vlon": vlon, "vfp": vfp, "vfn": vfn}
        self.diff_amp = diff_amp_r_load(**nodes)
        self.diff_amp.add(self.ckt)

        if sine_src:
            # Add VLop
            self.v_vlop = sine_voltage_src(omega=self.omega, mag=self.mag, phi=0, bias=self.bias)
            self.v_vlop.set_instance("VS")
            self.ckt.add_edge(vlop, GND, self.v_vlop)

            # Add VLon
            self.v_vlon = sine_voltage_src(omega=self.omega, mag=self.mag, phi=math.pi, bias=self.bias)
            self.v_vlon.set_instance("VS")
            self.ckt.add_edge(vlon, GND, self.v_vlon)

        elif two_sine_src:

            v_p1 = self.ckt.get_internal_node()
            v_n1 = self.ckt.get_internal_node()

            # Add VLop
            self.v_vlop_0 = sine_voltage_src(omega=self.omega_0, mag=self.mag, phi=0, bias=self.bias)
            self.v_vlop_0.set_instance("VS")
            self.ckt.add_edge(v_p1, GND, self.v_vlop_0)

            self.v_vlop_1 = sine_voltage_src(omega=self.omega_1, mag=self.mag, phi=0, bias=self.bias)
            self.v_vlop_1.set_instance("VS")
            self.ckt.add_edge(vlop, v_p1, self.v_vlop_1)

            # Add VLon
            self.v_vlon_0 = sine_voltage_src(omega=self.omega_0, mag=self.mag, phi=math.pi, bias=self.bias)
            self.v_vlon_0.set_instance("VS")
            self.ckt.add_edge(v_n1, GND, self.v_vlon_0)

            self.v_vlon_1 = sine_voltage_src(omega=self.omega_1, mag=self.mag, phi=math.pi, bias=self.bias)
            self.v_vlon_1.set_instance("VS")
            self.ckt.add_edge(vlon, v_n1, self.v_vlon_1)

        else:
            # Add VLop
            self.v_vlop = voltage_src()
            self.v_vlop.set_instance("VS")
            self.v_vlop.set_value(value=2.3)#2.25
            self.ckt.add_edge(vlop, GND, self.v_vlop)

            # Add VLon
            self.v_vlon = voltage_src()
            self.v_vlon.set_instance("VS")
            self.v_vlon.set_value(value=1.7)#1.75
            self.ckt.add_edge(vlon, GND, self.v_vlon)


        if sine_src or two_sine_src:
            # Add C1
            C1 = capacitor()
            C1.set_instance("C1")
            C1.set_value(value=1e-12)
            self.ckt.add_edge(vfp, GND, C1)

            # Add C2
            C2 = capacitor()
            C2.set_instance("C2")
            C2.set_value(value=1e-12)
            self.ckt.add_edge(vfn, GND, C2)

        # Add Vcc
        vs = voltage_src()
        vs.set_instance("VS")
        vs.set_value(value=5.0)
        self.ckt.add_edge(VCC, GND, vs)

        self.ckt.init_circuit()

        self.scb = Scoreboard(num_edges=self.ckt.num_edges, \
                              num_sys_vars=self.ckt.num_sys_vars, \
                              degen_mtrx=self.ckt.degen_mtrx)

    def add_sml_ckt(self, op, sine_src=False, **params):

        # Nodes
        GND  = 0
        VCC  = 1
        n1   = 2
        vlop = 3
        vlon = 4
        vfp  = 5
        vfn  = 6

        self.ckt_sml = Circuit()

        nodes = {"vcc": GND, "gnd": GND, "vlop": vlop, "vlon": vlon, "vfp": vfp, "vfn": vfn}
        self.diff_amp.add_small(op=op, ckt_sml=self.ckt_sml, **nodes)

        if sine_src:
            self.omega = params["omega"]
            self.mag   = params["mag"]
            self.bias  = params["bias"]

            self.vlop_bias = -(self.mag)
            self.vlon_bias =  (self.mag)

        if sine_src:
            # Add VLop
            self.v_vlop_sml = sine_voltage_src(omega=self.omega, mag=self.mag, phi=0, bias=self.vlop_bias)
            self.v_vlop_sml.set_instance("VS")
            self.ckt_sml.add_edge(vlop, GND, self.v_vlop_sml)

            # Add VLon
            self.v_vlon_sml = sine_voltage_src(omega=self.omega, mag=self.mag, phi=math.pi, bias=self.vlon_bias)
            self.v_vlon_sml.set_instance("VS")
            self.ckt_sml.add_edge(vlon, GND, self.v_vlon_sml)

        else:
            # Add VLop
            self.v_vlop_sml = voltage_src()
            self.v_vlop_sml.set_instance("VS")
            self.v_vlop_sml.set_value(value=0)
            self.ckt_sml.add_edge(vlop, GND, self.v_vlop_sml)

            # Add VLon
            self.v_vlon_sml = voltage_src()
            self.v_vlon_sml.set_instance("VS")
            self.v_vlon_sml.set_value(value=0)
            self.ckt_sml.add_edge(vlon, GND, self.v_vlon_sml)

        if sine_src:
            # Add C1
            C1 = capacitor()
            C1.set_instance("C1")
            C1.set_value(value=1e-12)
            self.ckt_sml.add_edge(vfp, GND, C1)

            # Add C2
            C2 = capacitor()
            C2.set_instance("C2")
            C2.set_value(value=1e-12)
            self.ckt_sml.add_edge(vfn, GND, C2)

        self.ckt_sml.init_circuit()

    def set_iref(self, iref):
        self.diff_amp.iref_n.set_value(value=iref)

    def set_iref_sml(self, iref):
        self.diff_amp.iref_n_sml.set_value(value=iref)

    def set_vlo(self, vlo):
        self.v_vlon.set_value(value=self.v_lo_bias - vlo)
        self.v_vlop.set_value(value=self.v_lo_bias + vlo)

    def set_vlo_sml(self, vlo):
        self.v_vlon_sml.set_value(value=-vlo)
        self.v_vlop_sml.set_value(value= vlo)

    def set_vlo_bias(self, bias):
        self.v_vlon.set_value(value=bias)
        self.v_vlop.set_value(value=bias)

        self.v_lo_bias = bias

    def check_saturation_region(self, x):

        if self.diff_amp.nmos_m1.ids.get_region(x) == 0 and self.diff_amp.nmos_m2.ids.get_region(x) == 0:
            return True
        return False

    def get_vf_output(self, x):

        self.scb.x = x

        rd1_edge = self.ckt.get_edge_info(self.diff_amp.rd1_ref)
        rd1_edge.get_voltage(scb=self.scb)
        v_rd1 = self.scb.v[self.diff_amp.rd1_ref]

        rd2_edge = self.ckt.get_edge_info(self.diff_amp.rd2_ref)
        rd2_edge.get_voltage(scb=self.scb)
        v_rd2 = self.scb.v[self.diff_amp.rd2_ref]

        return (5.0 - v_rd1) - (5.0 - v_rd2)

    def get_vf_output_base(self, x):

        self.scb.x = x

        rd1_edge = self.ckt.get_edge_info(self.diff_amp.rd1_ref)
        rd1_edge.get_voltage(scb=self.scb)
        v_rd1 = self.scb.v[self.diff_amp.rd1_ref]

        rd2_edge = self.ckt.get_edge_info(self.diff_amp.rd2_ref)
        rd2_edge.get_voltage(scb=self.scb)
        v_rd2 = self.scb.v[self.diff_amp.rd2_ref]

        return v_rd1, v_rd2

    def get_vf_output_sml(self, ckt, scb, v_rd1_base, v_rd2_base):

        rd1_edge = ckt.get_edge_info(self.rd1_ref_sml)
        rd1_edge.get_voltage(scb=scb)
        v_rd1 = scb.v[self.rd1_ref_sml]

        rd2_edge = ckt.get_edge_info(self.rd2_ref_sml)
        rd2_edge.get_voltage(scb=scb)
        v_rd2 = scb.v[self.rd2_ref_sml]

        return (5.0 - (v_rd1 + v_rd1_base)) - (5.0 - (v_rd2 + v_rd2_base))
