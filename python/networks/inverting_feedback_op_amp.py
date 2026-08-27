import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from op_amp_model import *

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

class inverting_feedback_op_amp():

    def __init__(self, lti=False, output_topology=None, input_type=None, gain_type="linear", feedback_type="resistor"):

        # Nodes
        GND   = 0
        VIN_P = 1
        VIN_N = 2
        VOUT  = 3
        INPUT = 4

        self.ckt = Circuit(lti=lti)

        self.instance = "{B_{1}}"

        if feedback_type == "None":
            VIN_N = INPUT

        # Add beta multiplier
        nodes = {"in_p": GND, "in_n": VIN_N, "out": VOUT, "gnd": GND}
        self.amp = op_amp_model(instance = self.instance, **nodes)
        self.amp.add(ckt=self.ckt, output_topology=output_topology, gain_type=gain_type)

        self.feedback_type = "None"

        if feedback_type == "resistor":
            # Add RF
            self.RF_idx = self.ckt.num_edges
            self.RF = resistor()
            self.RF.set_instance("R_F")
            self.ckt.add_edge(VOUT, VIN_N, self.RF)

            # Add RIN
            self.RIN_idx = self.ckt.num_edges
            self.RIN = resistor()
            self.RIN.set_instance("R_{IN}")
            self.ckt.add_edge(VIN_N, INPUT, self.RIN)

            self.feedback_type = "resistor"

        self.sources = {}

        if feedback_type == "None":
            # Add ROUT
            self.ROUT_idx = self.ckt.num_edges
            self.ROUT = resistor()
            self.ROUT.set_value(1e6)
            self.ROUT.set_instance("R_{OUT}")
            self.ckt.add_edge(VOUT, GND, self.ROUT)          

        if input_type is None:
            input_type = self._default_input_source
        else:
            input_type = input_type.__get__(self, type(self))

        input_type(INPUT, GND)

        # # Add Vdd
        # self.v_vdd = voltage_src()
        # self.v_vdd.set_instance("V_DD")
        # self.v_vdd.set_is_const()
        # self.ckt.add_edge(VDD, GND, self.v_vdd)

        # self.sources["vdd"] = self.v_vdd

        self.ckt.init_circuit()

        self.scb = Scoreboard(num_edges=self.ckt.num_edges, \
                              num_sys_vars=self.ckt.num_sys_vars, \
                              degen_mtrx=self.ckt.degen_mtrx)

    def _sine_source(self, vin=0, gnd=0):
        self.v_in = sine_voltage_src()
        self.v_in.set_instance("I_S")
        self.v_in.set_is_input()
        self.ckt.add_edge(vin, gnd, self.v_in)

        self.sources["sine_input_src"] = self.v_in

    def _default_input_source(self, vin, gnd):
        self.v_in = voltage_src()
        self.v_in.set_instance("I_S")
        self.v_in.set_is_input()
        self.ckt.add_edge(vin, gnd, self.v_in)

        self.sources["dc_input_src"] = self.v_in

    def get_sine_time_series_input(self, tr=0):
        return self.sources["sine_input_src"].get_time_series_input(tr)

    def set_sine_input_source(self, omega=0, mag=0, phi=0, bias=0):
        self.sources["sine_input_src"].set_params(omega=omega, \
                                                  mag=mag, \
                                                    phi=phi, \
                                                        bias=bias)

    def add_sml_ckt(self, op, lti=False, gate_topology=None, degen_topology=None, output_topology=None):
        pass

    def set_params(self, **params):
        op_amp_params = {"R_in": 100e6, 
                         "R_out": 1e-3, 
                         "linear_gain": params["linear_gain"], 
                         "a2": params["a2"], 
                         "alpha3": params["alpha3"],
                         "vclip": params["vclip"]}
        self.amp.set_params(**op_amp_params)

        if self.feedback_type == "resistor":
            self.RF.set_value(params["RF"])
            self.RIN.set_value(params["RIN"])

    def set_source(self, source, val):
        if source not in self.sources:
            raise KeyError(f"Unknown source: {source}")
        self.sources[source].set_value(val)

    def get_source_idx(self, source):
        return self.rc.get_source_idx(name=source)
    
    def get_vout(self):
        return self.amp.get_vout(self.ckt.scb)