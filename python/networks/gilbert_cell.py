import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from nmos_subckt import *

class gilbert_cell_switching_stage():

    def __init__(self, ckt):

        # Nodes
        GND  = 0
        VCC  = 1
        n1   = 2
        vlop = 3
        vlon = 4
        vfp  = 5
        vfn  = 6

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 10}

        # Add nmos
        nodes = {"g": vlop, "d": vfp, "s": n1}
        self.nmos_m1 = nmos_subckt(**nodes)
        self.nmos_m1.set_params(**nmos_params)
        self.nmos_m1.add(ckt)

        # Add nmos
        nodes = {"g": vlon, "d": vfn, "s": n1}
        self.nmos_m2 = nmos_subckt(**nodes)
        self.nmos_m2.set_params(**nmos_params)
        self.nmos_m2.add(ckt)

        # Add Iref
        self.iref_n = current_src()
        self.iref_n.set_instance("iref_n")
        self.iref_n.set_value(value=20e-6)
        ckt.add_edge(n1, GND, self.iref_n)

        # Add VLop
        self.v_vlop = voltage_src()
        self.v_vlop.set_instance("VS")
        self.v_vlop.set_value(value=2.6)
        ckt.add_edge(vlop, GND, self.v_vlop)

        # Add VLon
        self.v_vlon = voltage_src()
        self.v_vlon.set_instance("VS")
        self.v_vlon.set_value(value=2.5)
        ckt.add_edge(vlon, GND, self.v_vlon)

        # Add Rd1
        Rd1 = resistor()
        Rd1.set_instance("Rd1")
        Rd1.set_value(value=1e4)
        self.rd1_ref = ckt.add_edge(VCC, vfp, Rd1)

        # Add Rd2
        Rd2 = resistor()
        Rd2.set_instance("Rd2")
        Rd2.set_value(value=1e4)
        self.rd2_ref = ckt.add_edge(VCC, vfn, Rd2)

        # Add Vcc
        vs = voltage_src()
        vs.set_instance("VS")
        vs.set_value(value=5.0)
        ckt.add_edge(VCC, GND, vs)

    def add_sml_ckt(self, ckt_sml, op):

        # Nodes
        GND  = 0
        VCC  = 1
        n1   = 2
        vlop = 3
        vlon = 4
        vfp  = 5
        vfn  = 6

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 10}

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

        # Add VLop
        self.v_vlop_sml = voltage_src()
        self.v_vlop_sml.set_instance("VS")
        self.v_vlop_sml.set_value(value=0)
        ckt_sml.add_edge(vlop, GND, self.v_vlop_sml)

        # Add VLon
        self.v_vlon_sml = voltage_src()
        self.v_vlon_sml.set_instance("VS")
        self.v_vlon_sml.set_value(value=0)
        ckt_sml.add_edge(vlon, GND, self.v_vlon_sml)

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

    def set_iref(self, iref):
        self.iref_n.set_value(value=iref)

    def set_iref_sml(self, iref):
        self.iref_n_sml.set_value(value=iref)

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

    def check_saturation_region(self, scb):

        if self.nmos_m1.ids.get_region(scb.x) == 0 and self.nmos_m2.ids.get_region(scb.x) == 0:
            return True
        return False

    def get_vf_output(self, ckt, scb):

        rd1_edge = ckt.get_edge_info(self.rd1_ref)
        rd1_edge.get_voltage(scb=scb)
        v_rd1 = scb.v[self.rd1_ref]

        rd2_edge = ckt.get_edge_info(self.rd2_ref)
        rd2_edge.get_voltage(scb=scb)
        v_rd2 = scb.v[self.rd2_ref]

        return (5.0 - v_rd1) - (5.0 - v_rd2)

    def get_vf_output_base(self, ckt, scb):

        rd1_edge = ckt.get_edge_info(self.rd1_ref)
        rd1_edge.get_voltage(scb=scb)
        v_rd1 = scb.v[self.rd1_ref]

        rd2_edge = ckt.get_edge_info(self.rd2_ref)
        rd2_edge.get_voltage(scb=scb)
        v_rd2 = scb.v[self.rd2_ref]

        return v_rd1, v_rd2

    def get_vf_output_sml(self, ckt, scb, v_rd1_base, v_rd2_base):

        rd1_edge = ckt.get_edge_info(self.rd1_ref_sml)
        rd1_edge.get_voltage(scb=scb)
        v_rd1 = scb.v[self.rd1_ref_sml]

        rd2_edge = ckt.get_edge_info(self.rd2_ref_sml)
        rd2_edge.get_voltage(scb=scb)
        v_rd2 = scb.v[self.rd2_ref_sml]

        return (5.0 - (v_rd1 + v_rd1_base)) - (5.0 - (v_rd2 + v_rd2_base))

class gilbert_cell_switching_stage_dual():

    def __init__(self, ckt):

        # Nodes
        GND  = 0
        VCC  = 1
        n1   = 2
        n2   = 3
        vlop = 4
        vlon = 5
        vfp  = 6
        vfn  = 7

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 10}

        # Add nmos
        nodes = {"g": vlop, "d": vfp, "s": n1}
        self.nmos_m4 = nmos_subckt(**nodes)
        self.nmos_m4.set_params(**nmos_params)
        self.nmos_m4.add(ckt)

        # Add nmos
        nodes = {"g": vlon, "d": vfn, "s": n1}
        self.nmos_m5 = nmos_subckt(**nodes)
        self.nmos_m5.set_params(**nmos_params)
        self.nmos_m5.add(ckt)

        # Add nmos
        nodes = {"g": vlon, "d": vfp, "s": n2}
        self.nmos_m6 = nmos_subckt(**nodes)
        self.nmos_m6.set_params(**nmos_params)
        self.nmos_m6.add(ckt)

        # Add nmos
        nodes = {"g": vlop, "d": vfn, "s": n2}
        self.nmos_m7 = nmos_subckt(**nodes)
        self.nmos_m7.set_params(**nmos_params)
        self.nmos_m7.add(ckt)

        # Add Iref
        self.iref_n = current_src()
        self.iref_n.set_instance("iref_n")
        self.iref_n.set_value(value=20e-6)
        ckt.add_edge(n1, GND, self.iref_n)

        # Add Iref
        self.iref_p = current_src()
        self.iref_p.set_instance("iref_p")
        self.iref_p.set_value(value=20e-6)
        ckt.add_edge(n2, GND, self.iref_p)

        # Add VLop
        self.v_vlop = voltage_src()
        self.v_vlop.set_instance("VS")
        self.v_vlop.set_value(value=2.6)
        ckt.add_edge(vlop, GND, self.v_vlop)

        # Add VLon
        self.v_vlon = voltage_src()
        self.v_vlon.set_instance("VS")
        self.v_vlon.set_value(value=2.5)
        ckt.add_edge(vlon, GND, self.v_vlon)

        # Add Rd1
        Rd1 = resistor()
        Rd1.set_instance("Rd1")
        Rd1.set_value(value=1e4)
        self.rd1_ref = ckt.add_edge(VCC, vfp, Rd1)

        # Add Rd2
        Rd2 = resistor()
        Rd2.set_instance("Rd2")
        Rd2.set_value(value=1e4)
        self.rd2_ref = ckt.add_edge(VCC, vfn, Rd2)

        # Add Vcc
        vs = voltage_src()
        vs.set_instance("VS")
        vs.set_value(value=5.0)
        ckt.add_edge(VCC, GND, vs)

    def set_iref(self, iref_n, iref_p):
        self.iref_n.set_value(value=iref_n)
        self.iref_p.set_value(value=iref_p)

    def set_vlo(self, vlo):
        self.v_vlon.set_value(value=self.v_lo_bias - vlo)
        self.v_vlop.set_value(value=self.v_lo_bias + vlo)


    def set_vlo_bias(self, bias):
        self.v_vlon.set_value(value=bias)
        self.v_vlop.set_value(value=bias)

        self.v_lo_bias = bias

    def check_saturation_region(self, scb):

        if self.nmos_m4.ids.get_region(scb.x) == 0 and self.nmos_m5.ids.get_region(scb.x) == 0 and \
           self.nmos_m6.ids.get_region(scb.x) == 0 and self.nmos_m7.ids.get_region(scb.x) == 0:
            return True
        return False

    def get_vf_output(self, ckt, scb):

        rd1_edge = ckt.get_edge_info(self.rd1_ref)
        rd1_edge.get_voltage(scb=scb)
        v_rd1 = scb.v[self.rd1_ref]

        rd2_edge = ckt.get_edge_info(self.rd2_ref)
        rd2_edge.get_voltage(scb=scb)
        v_rd2 = scb.v[self.rd2_ref]

        return (5.0 - v_rd1) - (5.0 - v_rd2)

    def add_sml_ckt(self, ckt_sml, op):
        pass


class gilbert_cell_full():

    def __init__(self, ckt):

        # Nodes
        GND  = 0
        VCC  = 1
        n1   = 2
        n2   = 3
        n3   = 4
        vlop = 5
        vlon = 6
        vrfp = 7
        vrfn = 8
        vfp  = 9
        vfn  = 10
        vref = 11

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 10}

        # Add nmos
        nodes = {"g": vlop, "d": vfp, "s": n1}
        self.nmos_m4 = nmos_subckt(**nodes)
        self.nmos_m4.set_params(**nmos_params)
        self.nmos_m4.add(ckt)

        # Add nmos
        nodes = {"g": vlon, "d": vfn, "s": n1}
        self.nmos_m5 = nmos_subckt(**nodes)
        self.nmos_m5.set_params(**nmos_params)
        self.nmos_m5.add(ckt)

        # Add nmos
        nodes = {"g": vlon, "d": vfp, "s": n2}
        self.nmos_m6 = nmos_subckt(**nodes)
        self.nmos_m6.set_params(**nmos_params)
        self.nmos_m6.add(ckt)

        # Add nmos
        nodes = {"g": vlop, "d": vfn, "s": n2}
        self.nmos_m7 = nmos_subckt(**nodes)
        self.nmos_m7.set_params(**nmos_params)
        self.nmos_m7.add(ckt)

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 20}

        # Add nmos
        nodes = {"g": vrfp, "d": n1, "s": n3}
        self.nmos_m2 = nmos_subckt(**nodes)
        self.nmos_m2.set_params(**nmos_params)
        self.nmos_m2.add(ckt)

        # Add nmos
        nodes = {"g": vrfn, "d": n2, "s": n3}
        self.nmos_m3 = nmos_subckt(**nodes)
        self.nmos_m3.set_params(**nmos_params)
        self.nmos_m3.add(ckt)

        # Add Iref
        self.iref = current_src()
        self.iref.set_instance("iref")
        self.iref.set_value(value=20e-6)
        ckt.add_edge(n3, GND, self.iref)

        # Add nmos
        #nodes = {"g": vref, "d": n3, "s": GND}
        #self.nmos_m1 = nmos_subckt(**nodes)
        #self.nmos_m1.set_params(**nmos_params)
        #self.nmos_m1.add(ckt)

        # Add VLop
        self.v_vlop = voltage_src()
        self.v_vlop.set_instance("VS")
        self.v_vlop.set_value(value=2.6)
        ckt.add_edge(vlop, GND, self.v_vlop)

        # Add VLon
        self.v_vlon = voltage_src()
        self.v_vlon.set_instance("VS")
        self.v_vlon.set_value(value=2.5)
        ckt.add_edge(vlon, GND, self.v_vlon)

        # Add Vrfp
        self.v_vrfp = voltage_src()
        self.v_vrfp.set_instance("VS")
        self.v_vrfp.set_value(value=2.6)
        ckt.add_edge(vrfp, GND, self.v_vrfp)

        # Add Vrfn
        self.v_vrfn = voltage_src()
        self.v_vrfn.set_instance("VS")
        self.v_vrfn.set_value(value=2.5)
        ckt.add_edge(vrfn, GND, self.v_vrfn)

        # Add Rd1
        Rd1 = resistor()
        Rd1.set_instance("Rd1")
        Rd1.set_value(value=10e3)
        self.rd1_ref = ckt.add_edge(VCC, vfp, Rd1)

        # Add Rd2
        Rd2 = resistor()
        Rd2.set_instance("Rd2")
        Rd2.set_value(value=10e3)
        self.rd2_ref = ckt.add_edge(VCC, vfn, Rd2)

        # Add Vcc
        vs = voltage_src()
        vs.set_instance("VS")
        vs.set_value(value=5.0)
        ckt.add_edge(VCC, GND, vs)

        # Add Vref
        #v_vref = voltage_src()
        #v_vref.set_instance("VS")
        #v_vref.set_value(value=2.8)
        #ckt.add_edge(vref, GND, v_vref)

    def set_iref(self, iref):
        self.iref.set_value(value=iref)

    def set_vlo(self, vlo):
        self.v_vlop.set_value(value=vlo)

    def set_vrf(self, vrf):
        self.v_vrfp.set_value(value=vrf)

    def set_vlo_bias(self, bias):
        self.v_vlon.set_value(value=bias)
        self.v_vlop.set_value(value=bias)

    def set_vrf_bias(self, bias):
        self.v_vrfn.set_value(value=bias)
        self.v_vrfp.set_value(value=bias)

    def check_saturation_region(self, scb):

        if self.nmos_m4.ids.get_region(scb.x) == 0 and self.nmos_m5.ids.get_region(scb.x) == 0 and \
           self.nmos_m6.ids.get_region(scb.x) == 0 and self.nmos_m7.ids.get_region(scb.x) == 0 and \
           self.nmos_m2.ids.get_region(scb.x) == 0 and self.nmos_m3.ids.get_region(scb.x) == 0:
            return True

        self.nmos_m3.ids.print_stats()
        self.nmos_m2.ids.print_stats()

        return False

    def get_vf_output(self, ckt, scb):

        rd1_edge = ckt.get_edge_info(self.rd1_ref)
        rd1_edge.get_voltage(scb=scb)
        v_rd1 = scb.v[self.rd1_ref]

        rd2_edge = ckt.get_edge_info(self.rd2_ref)
        rd2_edge.get_voltage(scb=scb)
        v_rd2 = scb.v[self.rd2_ref]

        return (5.0 - v_rd1) - (5.0 - v_rd2)

    def add_sml_ckt(self, ckt_sml, op):
        pass


class gilbert_cell_full_dyn():

    def __init__(self, ckt):

        # Nodes
        GND  = 0
        VCC  = 1
        n1   = 2
        n2   = 3
        n3   = 4
        vlop = 5
        vlon = 6
        vrfp = 7
        vrfn = 8
        vfp  = 9
        vfn  = 10
        vref = 11

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 10}

        # Add nmos
        nodes = {"g": vlop, "d": vfp, "s": n1}
        self.nmos_m4 = nmos_subckt(**nodes)
        self.nmos_m4.set_params(**nmos_params)
        self.nmos_m4.add(ckt)

        # Add nmos
        nodes = {"g": vlon, "d": vfn, "s": n1}
        self.nmos_m5 = nmos_subckt(**nodes)
        self.nmos_m5.set_params(**nmos_params)
        self.nmos_m5.add(ckt)

        # Add nmos
        nodes = {"g": vlon, "d": vfp, "s": n2}
        self.nmos_m6 = nmos_subckt(**nodes)
        self.nmos_m6.set_params(**nmos_params)
        self.nmos_m6.add(ckt)

        # Add nmos
        nodes = {"g": vlop, "d": vfn, "s": n2}
        self.nmos_m7 = nmos_subckt(**nodes)
        self.nmos_m7.set_params(**nmos_params)
        self.nmos_m7.add(ckt)

        nmos_params = {"KP"     : 120e-6,
                       "vth"    : 0.8,
                       "lambda" : 0.01,
                       "L"      : 2,
                       "W"      : 20}

        # Add nmos
        nodes = {"g": vrfp, "d": n1, "s": n3}
        self.nmos_m2 = nmos_subckt(**nodes)
        self.nmos_m2.set_params(**nmos_params)
        self.nmos_m2.add(ckt)

        # Add nmos
        nodes = {"g": vrfn, "d": n2, "s": n3}
        self.nmos_m3 = nmos_subckt(**nodes)
        self.nmos_m3.set_params(**nmos_params)
        self.nmos_m3.add(ckt)

        # Add Iref
        self.iref = current_src()
        self.iref.set_instance("iref")
        self.iref.set_value(value=20e-6)
        ckt.add_edge(n3, GND, self.iref)

        # Add nmos
        #nodes = {"g": vref, "d": n3, "s": GND}
        #self.nmos_m1 = nmos_subckt(**nodes)
        #self.nmos_m1.set_params(**nmos_params)
        #self.nmos_m1.add(ckt)

        ## Add VLop
        #self.v_vlop = voltage_src()
        #self.v_vlop.set_instance("VS")
        #self.v_vlop.set_value(value=4.2)
        #ckt.add_edge(vlop, GND, self.v_vlop)
        #
        ## Add VLon
        #self.v_vlon = voltage_src()
        #self.v_vlon.set_instance("VS")
        #self.v_vlon.set_value(value=3.8)
        #ckt.add_edge(vlon, GND, self.v_vlon)

        # Add VLop
        self.v_vlop = sine_voltage_src(omega=1e6, mag=0.15, phi=0, bias=4.0)
        self.v_vlop.set_instance("VS")
        self.v_vlop.set_value(value=4.2)
        ckt.add_edge(vlop, GND, self.v_vlop)

        # Add VLon
        self.v_vlon = sine_voltage_src(omega=1e6, mag=0.15, phi=math.pi, bias=4.0)
        self.v_vlon.set_instance("VS")
        self.v_vlon.set_value(value=3.8)
        ckt.add_edge(vlon, GND, self.v_vlon)

        # Add Vrfp
        self.v_vrfp = sine_voltage_src(omega=5e6, mag=0.15, phi=0, bias=1.5)
        self.v_vrfp.set_instance("VS")
        self.v_vrfp.set_value(value=2.6)
        ckt.add_edge(vrfp, GND, self.v_vrfp)

        # Add Vrfn
        self.v_vrfn = sine_voltage_src(omega=5e6, mag=0.15, phi=math.pi, bias=1.5)
        self.v_vrfn.set_instance("VS")
        self.v_vrfn.set_value(value=2.5)
        ckt.add_edge(vrfn, GND, self.v_vrfn)

        # Add Rd1
        Rd1 = resistor()
        Rd1.set_instance("Rd1")
        Rd1.set_value(value=10e3)
        self.rd1_ref = ckt.add_edge(VCC, vfp, Rd1)

        # Add Rd2
        Rd2 = resistor()
        Rd2.set_instance("Rd2")
        Rd2.set_value(value=10e3)
        self.rd2_ref = ckt.add_edge(VCC, vfn, Rd2)

        # Add Vcc
        vs = voltage_src()
        vs.set_instance("VS")
        vs.set_value(value=5.0)
        ckt.add_edge(VCC, GND, vs)

        # Add C1
        C1 = capacitor()
        C1.set_instance("C1")
        C1.set_value(value=1e-12)
        ckt.add_edge(vfp, GND, C1)

        # Add C2
        C2 = capacitor()
        C2.set_instance("C2")
        C2.set_value(value=1e-12)
        ckt.add_edge(vfn, GND, C2)

        # Add Vref
        #v_vref = voltage_src()
        #v_vref.set_instance("VS")
        #v_vref.set_value(value=2.8)
        #ckt.add_edge(vref, GND, v_vref)

    def set_iref(self, iref):
        self.iref.set_value(value=iref)

    def set_vlo(self, vlo):
        self.v_vlop.set_value(value=vlo)

    def set_vrf(self, vrf):
        self.v_vrfp.set_value(value=vrf)

    def set_vlo_bias(self, bias):
        self.v_vlon.set_value(value=bias)
        self.v_vlop.set_value(value=bias)

    def set_vrf_bias(self, bias):
        self.v_vrfn.set_value(value=bias)
        self.v_vrfp.set_value(value=bias)

    def set_vlo(self, vlop, vlon):
        self.v_vlon.set_value(value=vlon)
        self.v_vlop.set_value(value=vlop)

    def set_vrf(self, vrfp, vrfn):
        self.v_vrfn.set_value(value=vrfn)
        self.v_vrfp.set_value(value=vrfp)

    def check_saturation_region(self, scb):

        if self.nmos_m4.ids.get_region(scb.x) == 0 and self.nmos_m5.ids.get_region(scb.x) == 0 and \
           self.nmos_m6.ids.get_region(scb.x) == 0 and self.nmos_m7.ids.get_region(scb.x) == 0 and \
           self.nmos_m2.ids.get_region(scb.x) == 0 and self.nmos_m3.ids.get_region(scb.x) == 0:
            return True

        self.nmos_m3.ids.print_stats()
        self.nmos_m2.ids.print_stats()

        return False

    def get_vf_output(self, ckt, scb):

        rd1_edge = ckt.get_edge_info(self.rd1_ref)
        rd1_edge.get_voltage(scb=scb)
        v_rd1 = scb.v[self.rd1_ref]

        rd2_edge = ckt.get_edge_info(self.rd2_ref)
        rd2_edge.get_voltage(scb=scb)
        v_rd2 = scb.v[self.rd2_ref]

        return (5.0 - v_rd1) - (5.0 - v_rd2)

    def add_sml_ckt(self, ckt_sml, op):
        pass
