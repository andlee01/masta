import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from subckt import *

class op_amp_model(subckt):

    def __init__(self, instance: str = "B_{1}", **nodes):

        self.IN_P   = nodes["in_p"]
        self.IN_N   = nodes["in_n"]
        self.OUT    = nodes["out"]
        self.GND    = nodes["gnd"]        

        self.instance = instance

        self.gain = 1e6  # Your a1 linear gain
            
        # Non-linear coefficients
        self.a2 = 0.0    # Set > 0 for 2nd harmonic distortion (asymmetry)
            
        # Alpha controls how fast the 3rd-order compression kicks in.
        # If alpha = 1e10, compression is noticeable when Vin approaches 10 microvolts.
        self.alpha3 = 1e10

        # Clipping limits
        self.vsat = 5.0

        self.test_sources = {}

    def set_params(self, **params):

        self.RIN_val  = params["R_in"]
        self.ROUT_val = params["R_out"]

        self.gain   = params["linear_gain"]
        self.a2     = params["a2"]
        self.alpha3 = params["alpha3"]
        self.vsat   = params["vclip"]

        self.RIN.set_value(self.RIN_val)
        self.ROUT.set_value(self.ROUT_val)

    def linear_gain(self, scb):
        return self.gain * scb.i[self.RIN_idx] * self.RIN.get_value()

    def nonlinear_gain(self, scb):
        v_in = scb.i[self.RIN_idx] * self.RIN.get_value()

        # 1. Calculate the EXACT alpha3 required to make the curve peak at vsat
        # This overrides self.alpha3 to guarantee mathematical continuity
        alpha3_ideal = (4.0 * self.gain**2) / (27.0 * self.vsat**2)

        # 2. Find the input voltage where this ideal curve peaks
        v_peak_in = np.sqrt(1.0 / (3.0 * alpha3_ideal))

        # 3. Clip the input to prevent the mathematical fold-back
        v_in_eff = np.clip(v_in, -v_peak_in, v_peak_in)

        # 4. Calculate the polynomial (it will now naturally peak at exactly vsat)
        v_lin = self.gain * v_in_eff
        v_cub = -1.0 * (self.gain * alpha3_ideal * (v_in_eff ** 3))

        return v_lin + v_cub

    def add(self, ckt, output_topology=None, gain_type="linear"):

        n1 = ckt.get_internal_node()

        # 1. Dictionary Mapping (Replaces if/else)
        # This is cleaner and easier to extend later
        gain_map = {
            "linear": self.linear_gain,
            "nonlinear": self.nonlinear_gain
        }
        
        # Use .get() to provide a fallback or raise a helpful error
        func = gain_map.get(gain_type)
        if not func:
            raise ValueError(f"Unknown gain_type: {gain_type}. Choose from {list(gain_map.keys())}")

        # Add Rin
        self.RIN_idx = ckt.num_edges
        self.RIN = resistor()
        self.RIN.set_instance(f"R_{{IN_{{{self.instance}}}}}")
        self.RIN.set_value(100e6)
        ckt.add_edge(self.IN_P, self.IN_N, self.RIN)

        # Add Vout
        self.VOUT_idx = ckt.num_edges
        self.VOUT = vcvs_src(func=func)
        self.VOUT.set_instance(f"V_{{OUT_{{{self.instance}}}}}")
        ckt.add_edge(n1, self.GND, self.VOUT)        

        # Add Rout
        self.ROUT_idx = ckt.num_edges
        self.ROUT = resistor()
        self.ROUT.set_instance(f"R_{{OUT_{{{self.instance}}}}}")
        self.ROUT.set_value(1e-3)
        ckt.add_edge(n1, self.OUT, self.ROUT)

        if output_topology is None:
            output_topology = self._default_output_topology

        output_topology(ckt)

    def set_source(self, name, value):
        if name not in self.test_sources:
            raise KeyError(f"Unknown source: {name}")
        self.test_sources[name].set_value(value)

    def get_source_idx(self, name):
        if name not in self.test_sources:
            raise KeyError(f"Unknown source: {name}")

        return self.test_sources[name].get_input_ref()

    def _default_output_topology(self, ckt):
        # Voltage across capacitor
        pass

    def _output_r_topology(self, ckt):
        
        self.vmeas = current_src()
        self.vmeas.set_is_const()
        self.vmeas.set_value(0.0)
        self.vmeas.set_is_output()
        self.vmeas.set_instance("VMEAS")
        ckt_sml.add_edge(self.VIN, self.VOUT, self.vmeas)

    def get_vout(self, scb):
        return scb.v[self.VOUT.i_x_ref] - scb.v[self.ROUT_idx]

    def add_small(self, op, ckt_sml, lti=False, gate_topology=None, degen_topology=None, output_topology=None, **nodes):
        pass