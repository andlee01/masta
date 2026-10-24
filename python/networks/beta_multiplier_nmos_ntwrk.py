from pathlib import Path
import sys
import numpy as np
from scipy.optimize import brentq

# Local imports
ROOT = Path(__file__).resolve().parents[1]
for directory in ("devices", "eqn", "subckt"):
    sys.path.insert(0, str(ROOT / directory))

from base_devices import current_src, resistor
from eqn_syn import Circuit
from nmos_subckt import nmos_subckt


class beta_multiplier_nmos_ntwrk:
    """Lower beta-multiplier pair driven by equal ideal drain currents."""

    def __init__(self, Rref=6.5e3, current=1e-6, l_lambda=0.0125):
        GND, gate, drain_m2, source_m2 = 0, 1, 2, 3
        self.ckt = Circuit()

        params = {"KP": 120e-6, "vth": 0.8, "lambda": l_lambda, "L": 2}
        self.m1 = nmos_subckt(instance="M_{1}", g=gate, d=gate, s=GND)
        self.m1.set_params(**params, W=10)
        self.m1.add(self.ckt)
        self.m2 = nmos_subckt(instance="M_{2}", g=gate, d=drain_m2, s=source_m2)
        self.m2.set_params(**params, W=40)
        self.m2.add(self.ckt)

        self.Rref = resistor()
        self.Rref.set_instance("R_{ref}")
        self.Rref.set_value(Rref)
        self.rref_idx = self.ckt.add_edge(source_m2, GND, self.Rref)

        self.sources = []
        for index, node in enumerate((gate, drain_m2), start=1):
            source = current_src()
            source.set_instance(f"I_{{test,{index}}}")
            source.set_value(current)
            self.ckt.add_edge(GND, node, source)
            self.sources.append(source)

        self.ckt.init_circuit()
        self.scb = self.ckt.scb

    def set_current(self, current):
        for source in self.sources:
            source.set_value(current)

    def current_driven_guess(self, current):
        """Seed from the existing device law, including its numerical gmin.

        At lambda=0 a forced current can require enormous drain voltages;
        a bounded Newton voltage step is unsuitable for reaching that seed.
        The full circuit residual is still solved and validated by the caller.
        """
        def positive_root(function):
            upper = 1.0
            while function(upper) < 0:
                upper *= 2
                if upper > 1e12:
                    raise ValueError("could not bracket current-driven device voltage")
            return brentq(function, 0, upper, xtol=1e-13)

        vg = positive_root(lambda v: self.m1.ids.drain_current(v, v) - current)
        vs = current * self.Rref.get_value()
        vds = positive_root(lambda v: self.m2.ids.drain_current(v, vg - vs) - current)
        guess = np.zeros(self.ckt.num_edges)
        for device, drain_voltage, gate_voltage in (
                (self.m1, vg, vg), (self.m2, vds, vg - vs)):
            guess[device.ids_ref] = drain_voltage
            guess[device.vgs_ref] = gate_voltage
            guess[device.ids.i_x_ref] = current
        guess[self.rref_idx] = current
        guess[self.sources[0].ref] = -vg
        guess[self.sources[1].ref] = -(vds + vs)
        return guess

    def get_operating_point(self):
        v = self.ckt.scb.v
        i = self.ckt.scb.i
        vs = v[self.rref_idx]
        return {
            "vg": float(v[self.m1.vgs_ref]),
            "vs2": float(vs),
            "vd2": float(v[self.m2.ids_ref] + vs),
            "vgs2": float(v[self.m2.vgs_ref]),
            "vds2": float(v[self.m2.ids_ref]),
            "overdrive2": float(v[self.m2.vgs_ref] - self.m2.vth),
            "i1": float(i[self.m1.ids_ref]),
            "i2": float(i[self.m2.ids_ref]),
        }
