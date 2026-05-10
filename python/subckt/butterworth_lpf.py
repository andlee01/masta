import sys

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from subckt import *

sys.path.append("../../../filter")
from butterworth import *

class butterworth_lpf(subckt):
    def __init__(self, instance: str = "F1", **nodes):
        self.N_IN = nodes["n_in"]
        self.N_OUT = nodes["n_out"]
        self.GND = nodes.get("gnd", 0)
        self.instance = instance
        self.components = [] 

    def add(self, ckt, params: dict):
        """
        Refactored LPF add method using LadderElement objects.
        """
        # ========================================================
        # 1. DESIGN LOWPASS FILTER
        # ========================================================
        result = design_lowpass_filter(
            Rl=params["Rl"],
            Rs=params["Rs"],
            wc=2 * math.pi * params["fc"],
            w=2 * math.pi * params["f_stop"],
            A=params["atten_db"],
        )

        topo = params.get("topology", "inductor_first")
        key = f"lowpass_{'inductor' if topo == 'inductor_first' else 'cap'}_first"
        elements = result[key]

        # ========================================================
        # 2. STORE INSTANCE VARIABLES
        # ========================================================
        self.order = result["order"]
        self.result = result
        self.components = []

        # --------------------------------------------------------
        # 3. COUNT SERIES ELEMENTS
        # --------------------------------------------------------
        series_count = sum(1 for e in elements if e.topology == "series")
        num_internal_nodes = max(0, series_count - 1)
        
        internal_nodes = [
            f"{self.instance}_n{i}"
            for i in range(num_internal_nodes)
        ]
        series_path = [self.N_IN] + internal_nodes + [self.N_OUT]
        series_idx = 0

        # ========================================================
        # 4. PLACE ELEMENTS (LADDER CONSTRUCTION)
        # ========================================================
        for elem in elements:

            # ----------------------------------------------------
            # SHUNT ELEMENT
            # ----------------------------------------------------
            if elem.topology == "shunt":
                node = series_path[series_idx]
                
                dev = capacitor() if elem.component_type == "C" else inductor()
                dev.set_instance(f"{elem.component_type}_{self.instance}_{elem.name}")
                dev.set_value(elem.value)

                ckt.add_edge(node, self.GND, dev)
                print(f"{node} -> {elem.component_type} -> {self.GND} (Shunt)")
                self.components.append(dev)

            # ----------------------------------------------------
            # SERIES ELEMENT
            # ----------------------------------------------------
            else: # topology == "series"
                node_a = series_path[series_idx]
                node_b = series_path[series_idx + 1]

                dev = inductor() if elem.component_type == "L" else capacitor()
                dev.set_instance(f"{elem.component_type}_{self.instance}_{elem.name}")
                dev.set_value(elem.value)

                ckt.add_edge(node_a, node_b, dev)
                print(f"{node_a} -> {elem.component_type} -> {node_b} (Series)")
                
                series_idx += 1
                self.components.append(dev)

    def set_params(self, **params):
        # Implementation remains the same as previous (lookup + scaling)
        wc = 2 * math.pi * params["fc"]
        result = design_lowpass_filter(
            Rl=params["Rl"], Rs=params["Rs"], wc=wc,
            w=2*math.pi*params["f_stop"], A=params["atten_db"]
        )
        
        topo_key = params.get("topology", "inductor_first")
        design_data = result[f"lowpass_{topo_key}"]

        for dev_obj, spec in zip(self.components, design_data):
            dev_obj.set_value(spec.value)

    def print_lowpass_filter(self):
        print_lowpass_filter(self.result)

    def add_small():
        pass