import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from parallel_lrc import *

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

class parallel_lrc_ntwrk():

    def __init__(self, lti=False, output_topology=None, input_type=None):

        # Nodes
        GND  = 0
        VDD   = 1
        N2   = 2

        self.ckt = Circuit(lti=lti)

        # Add beta multiplier
        nodes = {"n1": VDD, "n2": N2}
        self.rc = parallel_lrc(instance = "B_{1}", **nodes)
        self.rc.add(ckt=self.ckt, output_topology=output_topology)

        self.sources = {}

        if input_type is None:
            input_type = self._default_input_source
        else:
            input_type = input_type.__get__(self, type(self))

        input_type(N2, GND)

        # Add Vdd
        self.v_vdd = voltage_src()
        self.v_vdd.set_instance("V_DD")
        self.v_vdd.set_is_const()
        self.ckt.add_edge(VDD, GND, self.v_vdd)

        self.sources["vdd"] = self.v_vdd

        self.ckt.init_circuit()

        self.scb = Scoreboard(num_edges=self.ckt.num_edges, \
                              num_sys_vars=self.ckt.num_sys_vars, \
                              degen_mtrx=self.ckt.degen_mtrx)

    def _sine_source(self, vin=0, gnd=0):
        self.i_is = sine_current_src()
        self.i_is.set_instance("I_S")
        self.i_is.set_is_input()
        self.ckt.add_edge(vin, gnd, self.i_is)

        self.sources["sine_input_src"] = self.i_is

    def _default_input_source(self, vin, gnd):
        self.i_is = current_src()
        self.i_is.set_instance("I_S")
        self.i_is.set_is_input()
        self.ckt.add_edge(vin, gnd, self.i_is)

        self.sources["dc_input_src"] = self.i_is

    def get_sine_time_series_input(self, tr=0):
        return self.sources["sine_input_src"].get_time_series_input(tr)

    def set_sine_input_source(self, omega=0, mag=0, phi=0, bias=0):
        self.sources["sine_input_src"].set_params(omega=omega, \
                                                  mag=mag, \
                                                    phi=phi, \
                                                        bias=bias)

    def add_sml_ckt(self, op, lti=False, gate_topology=None, degen_topology=None, output_topology=None):
        pass

    def set_parallel_lrc_params(self, **params):
        self.rc.set_params(**params)
        self.v_vdd.set_value(params["VDD"])

    def set_source(self, source, val):
        self.rc.set_source(name=source, value=val)

    def get_source_idx(self, source):
        return self.rc.get_source_idx(name=source)