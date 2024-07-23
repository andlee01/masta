import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from gilbert_cell import *

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

class gilbert_cell_ntwrk():

    def __init__(self, sine_src=False, pulse_src=False, **params):

        # Nodes
        GND   = 0
        VCC   = 1
        n1    = 2
        vlop  = 3
        vlon  = 4
        vrfp  = 5
        vrfn  = 6
        vfp   = 7
        vfn   = 8

        self.ckt = Circuit()

        if sine_src or pulse_src:
            self.vlo_omega = params["vlo_omega"]
            self.vlo_mag   = params["vlo_mag"]
            self.vlo_bias  = params["vlo_bias"]

        self.vrf_omega = params["vrf_omega"]
        self.vrf_mag   = params["vrf_mag"]
        self.vrf_bias  = params["vrf_bias"]

        # Add gilbert cell
        nodes = {"vcc": VCC, "gnd": GND, "vlop": vlop, "vlon": vlon, "vrfp": vrfp, "vrfn": vrfn, "vfp": vfp, "vfn": vfn}
        self.glbrt_cell = gilbert_cell(**nodes)
        self.glbrt_cell.add(self.ckt)

        if sine_src:
            # Add VLop
            self.v_vlop = sine_voltage_src(omega=self.vlo_omega, mag=self.vlo_mag, phi=0, bias=self.vlo_bias)
            self.v_vlop.set_instance("VS")
            self.ckt.add_edge(vlop, GND, self.v_vlop)

            # Add VLon
            self.v_vlon = sine_voltage_src(omega=self.vlo_omega, mag=self.vlo_mag, phi=math.pi, bias=self.vlo_bias)
            self.v_vlon.set_instance("VS")
            self.ckt.add_edge(vlon, GND, self.v_vlon)

        elif pulse_src:

            # Add VLop
            self.v_vlop = pulse_voltage_src(T=1e-6, t_on=0.6e-6, t_r=100e-12, t_f=100e-12, t_del=0.25e-6, mag=4.0, bias=0.0)
            self.v_vlop.set_instance("VS")
            self.ckt.add_edge(vlop, GND, self.v_vlop)

            # Add VLon
            self.v_vlon = pulse_voltage_src(T=1e-6, t_on=0.6e-6, t_r=100e-12, t_f=100e-12, t_del=0.75e-6, mag=4.0, bias=0.0)
            self.v_vlon.set_instance("VS")
            self.ckt.add_edge(vlon, GND, self.v_vlon)

        # Add Vrfp
        self.v_vrfp = sine_voltage_src(omega=self.vrf_omega, mag=self.vrf_mag, phi=0, bias=self.vrf_bias)
        self.v_vrfp.set_instance("VS")
        self.ckt.add_edge(vrfp, GND, self.v_vrfp)

        # Add Vrfn
        self.v_vrfn = sine_voltage_src(omega=self.vrf_omega, mag=self.vrf_mag, phi=math.pi, bias=self.vrf_bias)
        self.v_vrfn.set_instance("VS")
        self.ckt.add_edge(vrfn, GND, self.v_vrfn)

        if sine_src or pulse_src:
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

        pass

    def set_iref(self, iref):
        self.glbrt_cell.iref.set_value(value=iref)

    def set_iref_sml(self, iref):
        self.glbrt_cell.iref_n_sml.set_value(value=iref)

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

        if self.glbrt_cell.nmos_m1.ids.get_region(x) == 0 and self.glbrt_cell.nmos_m2.ids.get_region(x) == 0:
            return True
        return False

    def check_switching_cutoff(self, x):
        if self.glbrt_cell.nmos_m4.ids.get_region(x) == -1 and \
           self.glbrt_cell.nmos_m5.ids.get_region(x) == -1 and \
           self.glbrt_cell.nmos_m6.ids.get_region(x) == -1 and \
           self.glbrt_cell.nmos_m7.ids.get_region(x) == -1:
            return True

        return False

    def get_vf_output(self, x):

        self.scb.x = x

        rd1_edge = self.ckt.get_edge_info(self.glbrt_cell.rd1_ref)
        rd1_edge.get_voltage(scb=self.scb)
        v_rd1 = self.scb.v[self.glbrt_cell.rd1_ref]

        rd2_edge = self.ckt.get_edge_info(self.glbrt_cell.rd2_ref)
        rd2_edge.get_voltage(scb=self.scb)
        v_rd2 = self.scb.v[self.glbrt_cell.rd2_ref]

        return (5.0 - v_rd1) - (5.0 - v_rd2)

    def get_vf_output_base(self, x):

        self.scb.x = x

        rd1_edge = self.ckt.get_edge_info(self.glbrt_cell.rd1_ref)
        rd1_edge.get_voltage(scb=self.scb)
        v_rd1 = self.scb.v[self.glbrt_cell.rd1_ref]

        rd2_edge = self.ckt.get_edge_info(self.glbrt_cell.rd2_ref)
        rd2_edge.get_voltage(scb=self.scb)
        v_rd2 = self.scb.v[self.glbrt_cell.rd2_ref]

        return v_rd1, v_rd2

    def get_vf_output_sml(self, ckt, scb, v_rd1_base, v_rd2_base):

        rd1_edge = ckt.get_edge_info(self.rd1_ref_sml)
        rd1_edge.get_voltage(scb=scb)
        v_rd1 = scb.v[self.rd1_ref_sml]

        rd2_edge = ckt.get_edge_info(self.rd2_ref_sml)
        rd2_edge.get_voltage(scb=scb)
        v_rd2 = scb.v[self.rd2_ref_sml]

        return (5.0 - (v_rd1 + v_rd1_base)) - (5.0 - (v_rd2 + v_rd2_base))
