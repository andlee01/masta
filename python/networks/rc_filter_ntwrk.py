import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from rc_filter import *

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

class rc_filter_ntwrk():

    def __init__(self, lti=False, output_topology=None, input_type=None):

        # Nodes
        GND    = 0
        VIN    = 1
        VOUT   = 2

        self.ckt = Circuit(lti=lti)

        # Add beta multiplier
        nodes = {"vin": VIN, "gnd": GND, "vout": VOUT}
        self.rc = rc_filter(instance = "B_{1}", **nodes)
        self.rc.add(ckt=self.ckt, output_topology=output_topology)

        self.sources = {}

        if input_type is None:
            input_type = self._default_input_source
        else:
            input_type = input_type.__get__(self, type(self))

        input_type(VIN, GND)

        # # Add Vs
        # self.v_vs = voltage_src()
        # self.v_vs.set_instance("V_S")
        # self.v_vs.set_is_input()
        # self.ckt.add_edge(VIN, GND, self.v_vs)

        self.ckt.init_circuit()

        self.scb = Scoreboard(num_edges=self.ckt.num_edges, \
                              num_sys_vars=self.ckt.num_sys_vars, \
                              degen_mtrx=self.ckt.degen_mtrx)

    def _sine_source(self, vin=0, gnd=0):
        self.v_vs = sine_voltage_src()
        self.v_vs.set_instance("V_S")
        self.v_vs.set_is_input()
        self.ckt.add_edge(vin, gnd, self.v_vs)

        self.sources["sine_input_src"] = self.v_vs

    def _default_input_source(self, vin, gnd):
        self.v_vs = voltage_src()
        self.v_vs.set_instance("V_S")
        self.v_vs.set_is_input()
        self.ckt.add_edge(vin, gnd, self.v_vs)

        self.sources["dc_input_src"] = self.v_vs

    def get_sine_time_series_input(self, tr=0):
        return self.sources["sine_input_src"].get_time_series_input(tr)

    def set_sine_input_source(self, omega=0, mag=0, phi=0, bias=0):
        self.sources["sine_input_src"].set_params(omega=omega, \
                                                  mag=mag, \
                                                    phi=phi, \
                                                        bias=bias)

    def add_sml_ckt(self, op, lti=False, gate_topology=None, degen_topology=None, output_topology=None):
        pass

    def set_rc_filter_params(self, **params):
        self.rc.set_params(**params)
        self.v_vs.set_value(params["VIN"])

    def set_source(self, source, val):
        self.rc.set_source(name=source, value=val)

    def get_source_idx(self, source):
        return self.rc.get_source_idx(name=source)