from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for directory in ("circuits", "devices", "eqn", "networks", "subckt"):
    sys.path.insert(0, str(ROOT / directory))

from base_devices import vccs_l1_mosfet
from beta_multiplier_circuit_tb import (
    _damped_newton,
    _seeded_guess,
    circuit_eqn,
    circuit_jacobian,
    op_solve,
)
from beta_multiplier_ntwrk import beta_multiplier_ntwrk


def make_mosfet():
    device = vccs_l1_mosfet()
    device.set_params(KP=120e-6, vth=0.8, l_lambda=0.0125, L=2, W=10)
    return device


class NumericalRobustnessTests(unittest.TestCase):
    def test_mos_dependencies_add_one_unknown_and_no_empty_equations(self):
        network = beta_multiplier_ntwrk()
        circuit = network.ckt
        self.assertEqual(circuit.num_elems, 10)
        self.assertEqual(circuit.num_edges, 14)
        point = _seeded_guess(circuit, np.zeros(circuit.num_edges), 0.5)
        jacobian = circuit_jacobian(point, 0, circuit, 0)
        self.assertEqual(np.linalg.matrix_rank(jacobian), circuit.num_edges)

    def test_mosfet_current_is_finite_and_continuous_at_region_boundaries(self):
        device = make_mosfet()
        points = [
            (0.0, 0.8),       # source/drain reversal and threshold
            (0.2, 1.0),       # classical triode/saturation boundary
            (-0.2, 1.0),
        ]
        for vds, vgs in points:
            values = np.array([
                device.drain_current(vds - 1e-8, vgs),
                device.drain_current(vds, vgs),
                device.drain_current(vds + 1e-8, vgs),
            ])
            self.assertTrue(np.all(np.isfinite(values)))
            self.assertLess(np.max(np.abs(np.diff(values))), 1e-9)

    def test_small_signal_parameters_match_dc_finite_differences(self):
        device = make_mosfet()
        vds, vgs = 1.5, 1.2
        gm, gds = device.small_signal_params(vds, vgs)
        h = 2e-6
        gm_check = (device.drain_current(vds, vgs + h) -
                    device.drain_current(vds, vgs - h)) / (2 * h)
        gds_check = (device.drain_current(vds + h, vgs) -
                     device.drain_current(vds - h, vgs)) / (2 * h)
        self.assertTrue(np.isclose(gm, gm_check, rtol=1e-7, atol=1e-12))
        self.assertTrue(np.isclose(gds, gds_check, rtol=1e-6, atol=1e-12))

    def test_beta_multiplier_selects_valid_nonzero_branch_from_zero_guess(self):
        network = beta_multiplier_ntwrk()
        network.set_beta_multiplier_params(Rref=6.5e3, Vs=5.0)
        solution = op_solve(
            network.ckt,
            np.zeros(network.ckt.num_edges),
            branch_ref=network.get_rref_idx(),
        )
        residual = circuit_eqn(solution, 0, network.ckt, 0)
        self.assertGreater(abs(solution[network.get_rref_idx()]), 1e-6)
        self.assertLess(np.max(np.abs(residual)), 1e-9)

    def test_companion_jacobian_matches_numerical_derivative(self):
        network = beta_multiplier_ntwrk()
        network.set_beta_multiplier_params(Rref=6.5e3, Vs=5.0)
        circuit = network.ckt
        point = _seeded_guess(circuit, np.zeros(circuit.num_edges), 0.5)
        analytic = circuit_jacobian(point, 0, circuit, 0)
        numerical = np.zeros_like(analytic)
        for column in range(circuit.num_edges):
            step = 1e-7
            plus = point.copy()
            minus = point.copy()
            plus[column] += step
            minus[column] -= step
            numerical[:, column] = (
                circuit_eqn(plus, 0, circuit, 0) -
                circuit_eqn(minus, 0, circuit, 0)
            ) / (2 * step)
        self.assertTrue(np.allclose(analytic, numerical, rtol=1e-7, atol=1e-7))

    def test_companion_newton_converges_without_nonlinear_fallback(self):
        network = beta_multiplier_ntwrk()
        network.set_beta_multiplier_params(Rref=6.5e3, Vs=5.0)
        initial = _seeded_guess(
            network.ckt, np.zeros(network.ckt.num_edges), voltage_seed=0.5
        )
        solution, converged, iterations, _ = _damped_newton(network.ckt, initial)
        self.assertTrue(converged)
        self.assertLess(iterations, 15)
        self.assertGreater(abs(solution[network.get_rref_idx()]), 1e-6)

    def test_ideal_output_resistance_uses_conducting_homotopy_branch(self):
        network = beta_multiplier_ntwrk(output_resistance=False, res_m2=True)
        network.set_beta_multiplier_params(Rref=6.5e3, Vs=4.0)
        solution = op_solve(
            network.ckt,
            np.ones(network.ckt.num_edges),
            branch_ref=network.get_rref_idx(),
        )
        self.assertGreater(abs(solution[network.get_rref_idx()]), 1e-6)
        self.assertLess(
            np.max(np.abs(circuit_eqn(solution, 0, network.ckt, 0))), 1e-9
        )

    def test_m1_degenerated_topology_continues_from_m2_solution(self):
        m2_network = beta_multiplier_ntwrk(output_resistance=True, res_m2=True)
        m2_network.set_beta_multiplier_params(Rref=6.5e3, Vs=4.0)
        m2_solution = op_solve(
            m2_network.ckt,
            np.zeros(m2_network.ckt.num_edges),
            branch_ref=m2_network.get_rref_idx(),
        )
        m1_network = beta_multiplier_ntwrk(output_resistance=True, res_m2=False)
        m1_network.set_beta_multiplier_params(Rref=6.5e3, Vs=4.0)
        m1_solution = op_solve(
            m1_network.ckt,
            m2_solution,
            branch_ref=m1_network.get_rref_idx(),
        )
        self.assertGreater(abs(m1_solution[m1_network.get_rref_idx()]), 1e-6)


if __name__ == "__main__":
    unittest.main()
