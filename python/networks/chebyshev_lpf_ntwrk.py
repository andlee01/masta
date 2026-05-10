import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from chebyshev_lpf import *

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

class chebyshev_lpf_ntwrk():

    def __init__(self, 
                 filter_reqs: dict, 
                 lti=False, 
                 input_type=None,
                 termination_type=None):
        """
        filter_reqs: dict containing Rs, Rl, fc, f_stop, atten_db, and topology
        """
        # Nodes
        GND = 0
        VIN_NODE = 1
        FILTER_IN = 2
        FILTER_OUT = 3

        self.ckt = Circuit(lti=lti)

        # 1. Add the Subcircuit Filter
        # We pass the full dict to .add() so it can calculate the order 
        # and build the correct number of L/C edges.
        filter_nodes = {"n_in": FILTER_IN, "n_out": FILTER_OUT, "gnd": GND}
        self.filter = chebyshev_lpf(instance="F_{1}", **filter_nodes)
        
        # This now triggers the chebyshev_order() calculation internally
        self.filter.add(ckt=self.ckt, params=filter_reqs)

        # 2. Add Source Termination Resistor (Rs)
        self.rs = resistor()
        self.rs.set_instance("R_S")
        self.ckt.add_edge(VIN_NODE, FILTER_IN, self.rs)

        # 3. Add Load Termination Resistor (Rl / LNA Input)
        self.terminations = {}

        if termination_type is None:
            termination_type = self._default_termination
        else:
            termination_type = termination_type.__get__(self, type(self))

        TERMINATION_OUT_MODE = termination_type(FILTER_OUT, GND)

        # 4. Add Input Source
        self.sources = {}
        if input_type is None:
            input_type = self._default_input_source
        else:
            input_type = input_type.__get__(self, type(self))
            
        input_type(VIN_NODE, GND)

        self.vmeas = current_src()
        self.vmeas.set_is_const()
        self.vmeas.set_value(0.0)
        self.vmeas.set_is_output()
        self.vmeas.set_instance("VMEAS")
        self.ckt.add_edge(TERMINATION_OUT_MODE, GND, self.vmeas)

        # Initialize Circuit and Scoreboard
        # The num_edges and sys_vars are now correct because .add() 
        # already created the ladder based on the required order.
        self.ckt.init_circuit()
        self.scb = Scoreboard(num_edges=self.ckt.num_edges, 
                              num_sys_vars=self.ckt.num_sys_vars, 
                              degen_mtrx=self.ckt.degen_mtrx)

    def _default_input_source(self, vin, gnd):
        self.v_in = voltage_src()
        self.v_in.set_instance("V_ANT")
        self.v_in.set_is_input()
        self.ckt.add_edge(vin, gnd, self.v_in)
        self.sources["input_src"] = self.v_in

    def _sine_source(self, vin, gnd):
        self.v_in = sine_voltage_src()
        self.v_in.set_instance("V_ANT_SINE")
        self.v_in.set_is_input()
        self.ckt.add_edge(vin, gnd, self.v_in)
        self.sources["input_src"] = self.v_in

    def _default_termination(self, node, gnd):

        self.rl = resistor()
        self.rl.set_instance("R_L")

        self.ckt.add_edge(node, gnd, self.rl)

        self.termination_setter = (
            self._set_default_termination_params
        )

        return node

    def _set_default_termination_params(self, **params):

        self.rl.set_value(params["Rl"])

    def _ac_rebias_termination(self, node, gnd):

        REBIAS_NODE = self.ckt.get_internal_node()
        VBIAS_NODE  = self.ckt.get_internal_node()

        # Termination resistance
        self.rl = resistor()
        self.rl.set_instance("R_L")
        self.ckt.add_edge(node, gnd, self.rl)

        # Rebias DC supplu
        self.v_bias = voltage_src()
        self.v_bias.set_instance("V_BIAS")
        self.v_bias.set_is_input()
        self.ckt.add_edge(VBIAS_NODE, gnd, self.v_bias)

        # AC coupling capacitor
        self.c_bias = capacitor()
        self.c_bias.set_instance("C_BIAS")
        self.ckt.add_edge(node, REBIAS_NODE, self.c_bias)

        # upper divider resistor
        self.r_bias_top = resistor()
        self.r_bias_top.set_instance("R_BIAS_TOP")
        self.ckt.add_edge(VBIAS_NODE, REBIAS_NODE, self.r_bias_top)

        # lower divider resistor
        self.r_bias_bot = resistor()
        self.r_bias_bot.set_instance("R_BIAS_BOT")
        self.ckt.add_edge(REBIAS_NODE, gnd, self.r_bias_bot)

        self.termination_setter = (
            self._set_ac_rebias_termination_params
        )

        return REBIAS_NODE

    def _set_ac_rebias_termination_params(self, **params):

        self.rl.set_value(params["Rl"])

        self.c_bias.set_value(
            params["C_BIAS"]
        )

        self.r_bias_top.set_value(
            params["R_BIAS_TOP"]
        )

        self.r_bias_bot.set_value(
            params["R_BIAS_BOT"]
        )

    def set_filter_params(self, **params):
        """
        Updates the numerical values (L, C, R) of the components.
        """
        # Set values for the termination resistors
        self.rs.set_value(params["Rs"])

        self.termination_setter(**params)

        # Pass dict to subcircuit to perform the scaling math and update component values
        self.filter.set_params(**params)

    # ... (rest of methods remain the same)