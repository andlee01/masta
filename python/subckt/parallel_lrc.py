import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from subckt import *

class parallel_lrc(subckt):

    def __init__(self, instance: str = "B_{1}", **nodes):

        self.N1   = nodes["n1"]
        self.N2   = nodes["n2"]

        self.instance = instance

        self.test_sources = {}

    def set_params(self, **params):

        self.R_val = params["R"]
        self.C_val = params["C"]
        self.L_val = params["L"]

        self.R.set_value(self.R_val)
        self.C.set_value(self.C_val)
        self.L.set_value(self.L_val)        

    def add(self, ckt, output_topology=None):

        # Add Cstray
        self.C_idx = ckt.num_edges
        self.C = capacitor()
        self.C.set_instance(f"C_{self.instance}")
        ckt.add_edge(self.N1, self.N2, self.C)

        # Add Rn
        self.R_idx = ckt.num_edges
        self.R = resistor()
        self.R.set_instance(f"R_{self.instance}")
        ckt.add_edge(self.N1, self.N2, self.R)

        # Add L
        self.L_idx = ckt.num_edges
        self.L = inductor()
        self.L.set_instance(f"L_{self.instance}")
        ckt.add_edge(self.N1, self.N2, self.L)

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

    def add_small(self, op, ckt_sml, lti=False, gate_topology=None, degen_topology=None, output_topology=None, **nodes):
        pass