import sys
import math

sys.path.append("../devices")
from base_devices import *
from base_devices import ElementType, TwoPortElement

sys.path.append("../subckt")
from subckt import *

sys.path.append("../../../filter")
from butterworth import design_bandpass_filter, print_bandpass_filter


class butterworth_bpf(subckt):

    def __init__(self, instance: str = "F1", **nodes):
        self.N_IN = nodes["n_in"]
        self.N_OUT = nodes["n_out"]
        self.GND = nodes.get("gnd", 0)
        self.instance = instance

        self.components = []
        self.resonators = []

    # ============================================================
    # Build circuit
    # ============================================================

    def add(self, ckt, params: dict):

        result = design_bandpass_filter(
            Rs=params["Rs"],
            Rl=params["Rl"],
            f0=params["f0"],
            BW_pass=params["BW_pass"],
            BW_stop=params["BW_stop"],
            attenuation_db=params["atten_db"],
        )

        topo = params.get("topology", "inductor_first")
        resonators = result[f"bandpass_{'inductor' if topo=='inductor_first' else 'cap'}_first"]

        self.order = result["order"]
        self.result = result

        self.components = []

        # --------------------------------------------------------
        # 1. COUNT ONLY SERIES RESONATORS
        # --------------------------------------------------------
        series_count = sum(1 for r in resonators if r.topology == "series")

        num_internal_nodes = max(0, series_count - 1)

        internal_nodes = [
            f"{self.instance}_n{i}"
            for i in range(num_internal_nodes)
        ]

        series_path = [self.N_IN] + internal_nodes + [self.N_OUT]

        series_idx = 0

        # --------------------------------------------------------
        # 2. PLACE ELEMENTS
        # --------------------------------------------------------
        for i, r in enumerate(resonators):

            # ====================================================
            # SERIES RESONATOR
            # ====================================================
            if r.topology == "series":

                node_a = series_path[series_idx]
                node_b = series_path[series_idx + 1]

                mid = f"{self.instance}_r{i}_mid"

                L = inductor()
                L.set_instance(f"L_{self.instance}_{r.name}")
                L.set_value(r.L)

                C = capacitor()
                C.set_instance(f"C_{self.instance}_{r.name}")
                C.set_value(r.C)

                ckt.add_edge(node_a, mid, L)
                ckt.add_edge(mid, node_b, C)

                print(f"{node_a} -> L -> {mid} -> C -> {node_b} (Series)")

                series_idx += 1
                self.components.extend([L, C])

            # ====================================================
            # PARALLEL RESONATOR
            # ====================================================
            else:

                node = series_path[series_idx]

                L = inductor()
                L.set_instance(f"L_{self.instance}_{r.name}")
                L.set_value(r.L)

                C = capacitor()
                C.set_instance(f"C_{self.instance}_{r.name}")
                C.set_value(r.C)

                ckt.add_edge(node, self.GND, L)
                ckt.add_edge(node, self.GND, C)

                print(f"{node} (Parallel resonator)")

                self.components.extend([L, C])

    # ============================================================
    # Parameter update (post-synthesis scaling hook)
    # ============================================================
    def set_params(self, **params):
            
        result = design_bandpass_filter(
            Rs=params["Rs"],
            Rl=params["Rl"],
            f0=params["f0"],
            BW_pass=params["BW_pass"],
            BW_stop=params["BW_stop"],
            attenuation_db=params["atten_db"],
        )
            
        topo = params.get("topology", "inductor_first")
        design_data = result[f"bandpass_{'inductor' if topo=='inductor_first' else 'cap'}_first"]

        # Flatten the L and C values so they map 1:1 with self.components
        flat_values = []
        for r in design_data:
            flat_values.extend([r.L, r.C])

        # Now they perfectly align
        for dev, val in zip(self.components, flat_values):
            if hasattr(dev, "set_value"):
                dev.set_value(val)

    def print_bandpass_filter(self):
        print_bandpass_filter(self.result)

    def add_small(self):
        pass