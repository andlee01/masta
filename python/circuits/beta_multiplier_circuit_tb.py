import sys
import os
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

from scipy import optimize
from scipy.integrate import ode

# Local imports
sys.path.extend(["../eqn", "../devices", "../networks"])
from eqn_syn import Circuit, Scoreboard
from base_devices import *
from beta_multiplier_ntwrk import *

import control as ct
from scipy.signal import freqresp, StateSpace

from pathlib import Path

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

DEFAULT_RREF = 6.5e3
DEFAULT_VS = 5.0
SWEEP_RANGE = (-0.2, 0.2, 100)
LRG_SWEEP_RANGE = (4, 6, int(1000))
RREF_SWEEP_RANGE = (2.5e3, 6.5e3, int(10))
TEND = 50e-12
TSTEP = 10000
TOL = 1e-6
M2_SWEEP_RANGE = (0, 1, 100)

CURRENT_SCALE = 100e-6
VOLTAGE_SCALE = 1.0
CURRENT_ABSTOL = 1e-9
VOLTAGE_ABSTOL = 1e-8


class OperatingPointError(RuntimeError):
    pass

# ------------------------------------------------------------------------------
# Core circuit equations
# ------------------------------------------------------------------------------

def circuit_eqn(x, sys, ckt, t):
    """Return Kirchhoff equations for given circuit state x."""
    ckt.get_im(x=x, sys=sys, t=t)
    V = ckt.get_vm(x=x, sys=sys, t=t)
    I = ckt.get_im_dep()
    return ckt.qf @ I + ckt.bf @ V


def _equation_scales(ckt):
    """Scale KCL rows as currents and KVL rows as voltages."""
    voltage_rows = np.any(np.abs(ckt.bf) > 0.0, axis=1)
    return np.where(voltage_rows, VOLTAGE_SCALE, CURRENT_SCALE)


def _equation_tolerances(ckt):
    voltage_rows = np.any(np.abs(ckt.bf) > 0.0, axis=1)
    return np.where(voltage_rows, VOLTAGE_ABSTOL, CURRENT_ABSTOL)


def _variable_scales(ckt):
    scales = np.full(ckt.num_edges, CURRENT_SCALE)
    # Original current-source edge variables are voltages.  Dependency
    # variables appended during init_circuit are currents.
    for ref in range(ckt.num_elems):
        edge = ckt.get_edge_info(ref)
        if edge is not None and edge.get_type() == ElementType.current_src:
            scales[ref] = VOLTAGE_SCALE
    return scales


def _seeded_guess(ckt, guess, voltage_seed=1.0):
    """Create a conducting seed without treating current unknowns as 1 A."""
    seeded = np.asarray(guess, dtype=float).copy()
    x_scale = _variable_scales(ckt)
    current_vars = x_scale == CURRENT_SCALE
    voltage_vars = ~current_vars
    seeded[current_vars & (np.abs(seeded) < 1e-15)] = CURRENT_SCALE
    seeded[current_vars & (np.abs(seeded) > 0.1)] = CURRENT_SCALE
    seeded[voltage_vars & (np.abs(seeded) < 1e-15)] = voltage_seed
    return seeded


def circuit_jacobian(x, sys, ckt, t):
    """Stamp the analytic Jacobian of the circuit residual.

    This is the matrix form of the usual nonlinear-device companion model:
    each MOSFET contributes its present gm and gds to the linearized circuit.
    """
    size = ckt.num_edges
    dcurrent = np.zeros((size, size))
    dvoltage = np.zeros((size, size))

    # Populate the scoreboard voltages used by controlled devices.
    ckt.get_im(x=x, sys=sys, t=t)
    ckt.get_vm(x=x, sys=sys, t=t)

    for _, _, data in ckt.G.edges(data=True):
        edge = data["info"]
        ref = edge.ref

        if isinstance(edge, vccs_l1_mosfet):
            vds = ckt.scb.v[ref]
            vgs = ckt.scb.v[edge.vgs_ref]
            _, gds, gm = edge.current_and_derivatives(vds, vgs)
            dvoltage[ref, ref] = 1.0
            dcurrent[edge.i_x_ref, edge.i_x_ref] = 1.0
            dcurrent[ref, ref] = gds
            dcurrent[ref, edge.vgs_ref] = gm
        elif isinstance(edge, vccs_l1_mosfet_small):
            dvoltage[ref, ref] = 1.0
            dcurrent[edge.i_x_ref, edge.i_x_ref] = 1.0
            dcurrent[ref, edge.vgs_ref] = edge.gm
        elif isinstance(edge, resistor):
            dcurrent[ref, ref] = 1.0
            dvoltage[ref, ref] = edge.get_value()
        elif edge.get_type() == ElementType.voltage_src:
            dcurrent[ref, ref] = 1.0
        elif edge.get_type() == ElementType.current_src:
            dvoltage[ref, ref] = 1.0
        elif edge.get_type() == ElementType.capacitor:
            dcurrent[ref, ref] = 1.0
        elif edge.get_type() == ElementType.inductor:
            dvoltage[ref, ref] = 1.0
        else:
            raise NotImplementedError(
                f"no analytic Jacobian stamp for {type(edge).__name__}"
            )

    return ckt.qf @ dcurrent + ckt.bf @ dvoltage


def _damped_newton(ckt, guess, sys=0, t=0, max_iterations=40):
    """Solve one companion-model/Newton sequence in scaled coordinates."""
    x = np.asarray(guess, dtype=float).copy()
    eq_scale = _equation_scales(ckt)
    eq_tol = _equation_tolerances(ckt)
    x_scale = _variable_scales(ckt)

    for iteration in range(max_iterations):
        residual = circuit_eqn(x, sys, ckt, t)
        if (np.all(np.isfinite(residual)) and
                np.max(np.abs(residual) / eq_tol) <= 1.0):
            return x, True, iteration, "converged"

        try:
            jacobian = circuit_jacobian(x, sys, ckt, t)
        except (NotImplementedError, AttributeError):
            return x, False, iteration, "analytic stamp unavailable"

        # Solve in normalized variables/equations to avoid mixing volts and
        # microamps in the Newton linear system.
        scaled_jacobian = ((jacobian / eq_scale[:, None]) *
                           x_scale[None, :])
        scaled_residual = residual / eq_scale
        try:
            delta_scaled = np.linalg.solve(scaled_jacobian, -scaled_residual)
        except np.linalg.LinAlgError:
            # The fundamental cut-set/loop formulation can retain redundant
            # zero or dependent rows.  A rank-revealing least-norm linear
            # solve is still a Newton step; unlike nonlinear least_squares it
            # does not reevaluate the whole circuit to estimate a Jacobian.
            delta_scaled, _, rank, _ = np.linalg.lstsq(
                scaled_jacobian, -scaled_residual, rcond=None
            )
            if rank == 0:
                return x, False, iteration, "singular companion matrix"

        delta = delta_scaled * x_scale
        if not np.all(np.isfinite(delta)):
            return x, False, iteration, "non-finite Newton step"

        # Device voltage limiting and current limiting prevent a poor seed
        # from crossing several operating regions in one iteration.
        voltage_vars = x_scale == VOLTAGE_SCALE
        delta[voltage_vars] = np.clip(delta[voltage_vars], -0.5, 0.5)
        delta[~voltage_vars] = np.clip(
            delta[~voltage_vars], -5 * CURRENT_SCALE, 5 * CURRENT_SCALE
        )

        old_norm = np.linalg.norm(scaled_residual, ord=np.inf)
        alpha = 1.0
        accepted = False
        while alpha >= 2.0**-12:
            trial = x + alpha * delta
            trial_residual = circuit_eqn(trial, sys, ckt, t)
            trial_norm = np.linalg.norm(trial_residual / eq_scale, ord=np.inf)
            if np.isfinite(trial_norm) and trial_norm < old_norm:
                x = trial
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            return x, False, iteration, "line search stalled"

    return x, False, max_iterations, "iteration limit"


def op_solve(ckt, root_start, branch_ref=None, min_branch_current=1e-9,
             max_nfev=3000):
    """Find and validate a scaled DC operating point.

    If ``branch_ref`` is supplied, both the caller's guess and a conducting
    seed are tried and the finite-current branch is selected explicitly.
    """
    root_start = np.asarray(root_start, dtype=float)
    if root_start.shape != (ckt.num_edges,):
        raise ValueError(f"expected {ckt.num_edges} initial values, got {root_start.shape}")

    eq_scale = _equation_scales(ckt)
    eq_tol = _equation_tolerances(ckt)
    x_scale = _variable_scales(ckt)

    def scaled_residual(x):
        residual = circuit_eqn(x, 0, ckt, 0)
        if not np.all(np.isfinite(residual)):
            return np.full_like(residual, 1e30)
        return residual / eq_scale

    guesses = [root_start]
    if branch_ref is not None:
        # The all-off solution has a wide basin of attraction.  Try modest
        # and strongly conducting voltage seeds, then select a validated
        # finite-current solution rather than trusting one arbitrary seed.
        for voltage_seed in (0.5, 1.0, 5.0):
            seeded = _seeded_guess(ckt, np.zeros_like(root_start), voltage_seed)
            if not any(np.array_equal(seeded, old) for old in guesses):
                guesses.append(seeded)

        # With ideal output resistance (lambda == 0), the conducting and
        # all-off branches can have nearly disjoint attraction basins.  Use a
        # temporary finite-gds companion model to obtain a nearby conducting
        # seed, then restore the requested ideal devices before the real solve.
        ideal_mosfets = []
        for _, _, data in ckt.G.edges(data=True):
            edge = data["info"]
            if isinstance(edge, vccs_l1_mosfet) and edge.l_lambda == 0.0:
                ideal_mosfets.append(edge)
        branch_guess = abs(root_start[branch_ref])
        cold_or_unphysical = (
            branch_guess < min_branch_current or
            branch_guess > 10.0 * CURRENT_SCALE
        )
        if ideal_mosfets and cold_or_unphysical:
            for edge in ideal_mosfets:
                edge.l_lambda = 0.01
            try:
                homotopy_start = _seeded_guess(
                    ckt, np.zeros_like(root_start), voltage_seed=0.1
                )
                homotopy_x, homotopy_ok, _, _ = _damped_newton(
                    ckt, homotopy_start
                )
            finally:
                for edge in ideal_mosfets:
                    edge.l_lambda = 0.0
            if (homotopy_ok and
                    abs(homotopy_x[branch_ref]) >= min_branch_current):
                guesses.insert(0, homotopy_x)

    diagnostics = []
    for guess in guesses:
        newton_x, newton_ok, iterations, newton_message = _damped_newton(
            ckt, guess
        )
        newton_residual = circuit_eqn(newton_x, 0, ckt, 0)
        newton_error = np.max(np.abs(newton_residual) / eq_tol)
        newton_branch_ok = (branch_ref is None or
                            abs(newton_x[branch_ref]) >= min_branch_current)
        diagnostics.append(("newton: " + newton_message, newton_error,
                            newton_x[branch_ref] if branch_ref is not None else None,
                            iterations))
        if newton_ok and newton_branch_ok:
            circuit_eqn(newton_x, 0, ckt, 0)
            return newton_x

        # Retain nonlinear least squares as a robust fallback for unsupported
        # elements, singular companion matrices, or unusually poor seeds.
        result = optimize.least_squares(
            scaled_residual,
            newton_x,
            x_scale=x_scale,
            jac="3-point",
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
            max_nfev=max_nfev,
        )
        residual = circuit_eqn(result.x, 0, ckt, 0)
        normalized_error = np.max(np.abs(residual) / eq_tol)
        valid = (result.success and np.all(np.isfinite(result.x)) and
                 np.all(np.isfinite(residual)) and normalized_error <= 1.0)
        branch_ok = (branch_ref is None or
                     abs(result.x[branch_ref]) >= min_branch_current)
        diagnostics.append(("least_squares: " + result.message, normalized_error,
                            result.x[branch_ref] if branch_ref is not None else None,
                            result.nfev))
        if valid and branch_ok:
            circuit_eqn(result.x, 0, ckt, 0)
            return result.x

    raise OperatingPointError(
        "operating-point solve failed residual/branch validation: " +
        repr(diagnostics)
    )

def dypc_litmus0(t, state, ckt, algebraic_state):
    """ODE callback with a persistent warm start for the algebraic solve."""
    eq_scale = _equation_scales(ckt)
    eq_tol = _equation_tolerances(ckt)
    x_scale = _variable_scales(ckt)

    solution, converged, _, message = _damped_newton(
        ckt, algebraic_state["x"], sys=state, t=t, max_iterations=20
    )
    if not converged:
        result = optimize.least_squares(
            lambda x: circuit_eqn(x, state, ckt, t) / eq_scale,
            solution,
            x_scale=x_scale,
            jac="3-point",
            ftol=1e-10,
            xtol=1e-10,
            gtol=1e-10,
            max_nfev=1000,
        )
        solution = result.x
        message = result.message
        converged = result.success

    residual = circuit_eqn(solution, state, ckt, t)
    if (not converged or not np.all(np.isfinite(solution)) or
            np.max(np.abs(residual) / eq_tol) > 1.0):
        raise OperatingPointError(
            f"transient algebraic solve failed at t={t:g}: {message}"
        )
    algebraic_state["x"] = solution.copy()
    return ckt.get_dy(x=solution)

def ode_solve(ckt, tend=50e-12, tstep=10000, x0=0):

    num_sys_vars    = ckt.get_num_sys_vars()

    t     = np.linspace(0, tend, tstep)

    algebraic_state = {
        "x": _seeded_guess(ckt, np.zeros(ckt.num_edges), voltage_seed=0.5)
    }
    r = ode(dypc_litmus0).set_integrator('lsoda', method='bdf', atol=1e-9, rtol=1e-9)
    r.set_initial_value(x0, 0.0)
    r.set_f_params(ckt, algebraic_state)

    y    = np.empty((tstep, num_sys_vars))
    y[0] = x0

    k = 1
    while r.successful() and k < tstep:
        r.integrate(t[k])

        y[k] = r.y

        k += 1

    return t, y

# ------------------------------------------------------------------------------
# Large-signal analysis
# ------------------------------------------------------------------------------

def sweep(nw, vs_sweep, root_start=None):
    """Perform small-signal sweep of beta-multiplier network."""
    N = len(vs_sweep)

    # Preallocate arrays
    isd_m3_sweep = np.zeros(N)
    vsg_m3_sweep = np.zeros(N)
    isd_m4_sweep = np.zeros(N)
    vsg_m4_sweep = np.zeros(N)
    ids_m1_sweep = np.zeros(N)
    vgs_m1_sweep = np.zeros(N)
    ids_m2_sweep = np.zeros(N)
    vgs_m2_sweep = np.zeros(N)
    irref_sweep  = np.zeros(N)

    if root_start is None:
        root_start = np.ones(nw.ckt.num_edges)
    else:
        root_start = np.asarray(root_start, dtype=float).copy()

    for vs_idx, vs in enumerate(vs_sweep):
        params = {"Rref": 6.5e3, "Vs": vs}
        nw.set_beta_multiplier_params(**params)

        x = op_solve(nw.ckt, root_start, branch_ref=nw.get_rref_idx())
        if vs_idx == 0:
            nw._sweep_first_solution = x.copy()
        root_start = x.copy()

        irref_sweep[vs_idx] = nw.get_Rref_current(x)

        # flatten each returned value to scalar
        vals = [np.atleast_1d(v).item() for v in nw.get_mos_vltg(x)]
        (isd_m3_sweep[vs_idx],
         vsg_m3_sweep[vs_idx],
         isd_m4_sweep[vs_idx],
         vsg_m4_sweep[vs_idx],
         ids_m1_sweep[vs_idx],
         vgs_m1_sweep[vs_idx],
         ids_m2_sweep[vs_idx],
         vgs_m2_sweep[vs_idx]) = vals

    return {
        "$v_{i_{SD_{B_{1}M_{3}}}}$": isd_m3_sweep,
        "$v_{i_{SG_{B_{1}M_{3}}}}$": vsg_m3_sweep,
        "$v_{i_{SD_{B_{1}M_{4}}}}$": isd_m4_sweep,
        "$v_{i_{SG_{B_{1}M_{4}}}}$": vsg_m4_sweep,
        "$v_{i_{DS_{B_{1}M_{1}}}}$": ids_m1_sweep,
        "$v_{i_{GS_{B_{1}M_{1}}}}$": vgs_m1_sweep,
        "$v_{i_{DS_{B_{1}M_{2}}}}$": ids_m2_sweep,
        "$v_{i_{GS_{B_{1}M_{2}}}}$": vgs_m2_sweep,
        "$i_{B_{1}R_{ref}}$":        irref_sweep
    }

# ------------------------------------------------------------------------------
# Small-signal analysis
# ------------------------------------------------------------------------------

def sml_sweep(nw, vs_sml_sweep):
    """Perform small-signal sweep of beta-multiplier network."""
    N = len(vs_sml_sweep)

    # Preallocate arrays
    isd_m3_sml_sweep = np.zeros(N)
    vsg_m3_sml_sweep = np.zeros(N)
    isd_m4_sml_sweep = np.zeros(N)
    vsg_m4_sml_sweep = np.zeros(N)
    ids_m1_sml_sweep = np.zeros(N)
    vgs_m1_sml_sweep = np.zeros(N)
    ids_m2_sml_sweep = np.zeros(N)
    vgs_m2_sml_sweep = np.zeros(N)
    irref_sml_sweep  = np.zeros(N)

    root_start = np.zeros(nw.ckt_sml.num_edges)
    for vs_idx, vs in enumerate(vs_sml_sweep):
        params = {"Rref": 6.5e3, "Vs": vs}
        nw.set_beta_multiplier_sml_params(**params)

        x = op_solve(nw.ckt_sml, root_start)
        root_start = x.copy()

        irref_sml_sweep[vs_idx] = nw.get_Rref_sml_current(x)

        # flatten each returned value to scalar
        vals = [np.atleast_1d(v).item() for v in nw.get_mos_vltg_sml(x)]
        (isd_m3_sml_sweep[vs_idx],
         vsg_m3_sml_sweep[vs_idx],
         isd_m4_sml_sweep[vs_idx],
         vsg_m4_sml_sweep[vs_idx],
         ids_m1_sml_sweep[vs_idx],
         vgs_m1_sml_sweep[vs_idx],
         ids_m2_sml_sweep[vs_idx],
         vgs_m2_sml_sweep[vs_idx]) = vals

    return {
        "$v_{i_{SD_{B_{1}M_{3}}}}$": isd_m3_sml_sweep,
        "$v_{i_{SG_{B_{1}M_{3}}}}$": vsg_m3_sml_sweep,
        "$v_{i_{SD_{B_{1}M_{4}}}}$": isd_m4_sml_sweep,
        "$v_{i_{SG_{B_{1}M_{4}}}}$": vsg_m4_sml_sweep,
        "$v_{i_{DS_{B_{1}M_{1}}}}$": ids_m1_sml_sweep,
        "$v_{i_{GS_{B_{1}M_{1}}}}$": vgs_m1_sml_sweep,
        "$v_{i_{DS_{B_{1}M_{2}}}}$": ids_m2_sml_sweep,
        "$v_{i_{GS_{B_{1}M_{2}}}}$": vgs_m2_sml_sweep,
        "$i_{B_{1}R_{ref}}$":        irref_sml_sweep
    }

def sml_sweep_m2(nw, vs_sml_sweep):
    """Perform small-signal sweep of beta-multiplier network."""
    N = len(vs_sml_sweep)

    # Preallocate arrays
    isd_m3_sml_sweep = np.zeros(N)
    vsg_m3_sml_sweep = np.zeros(N)
    isd_m4_sml_sweep = np.zeros(N)
    vsg_m4_sml_sweep = np.zeros(N)
    ids_m1_sml_sweep = np.zeros(N)
    vgs_m1_sml_sweep = np.zeros(N)
    ids_m2_sml_sweep = np.zeros(N)
    vgs_m2_sml_sweep = np.zeros(N)
    irref_sml_sweep  = np.zeros(N)

    root_start = np.zeros(nw.ckt_sml.num_edges)
    for vs_idx, vs in enumerate(vs_sml_sweep):
        nw.set_source(source="gate_m2", val=vs)
        params = {"Rref": DEFAULT_RREF, "Vs": 0}
        nw.set_beta_multiplier_sml_params(**params)
        nw.set_source("gate_m2", vs)

        x = op_solve(nw.ckt_sml, root_start)
        root_start = x.copy()

        irref_sml_sweep[vs_idx] = nw.get_Rref_sml_current(x)

        # flatten each returned value to scalar
        vals = [np.atleast_1d(v).item() for v in nw.get_mos_vltg_sml(x)]
        (isd_m3_sml_sweep[vs_idx],
         vsg_m3_sml_sweep[vs_idx],
         isd_m4_sml_sweep[vs_idx],
         vsg_m4_sml_sweep[vs_idx],
         ids_m1_sml_sweep[vs_idx],
         vgs_m1_sml_sweep[vs_idx],
         ids_m2_sml_sweep[vs_idx],
         vgs_m2_sml_sweep[vs_idx]) = vals

    return {
        "$v_{i_{SD_{B_{1}M_{3}}}}$": isd_m3_sml_sweep,
        "$v_{i_{SG_{B_{1}M_{3}}}}$": vsg_m3_sml_sweep,
        "$v_{i_{SD_{B_{1}M_{4}}}}$": isd_m4_sml_sweep,
        "$v_{i_{SG_{B_{1}M_{4}}}}$": vsg_m4_sml_sweep,
        "$v_{i_{DS_{B_{1}M_{1}}}}$": ids_m1_sml_sweep,
        "$v_{i_{GS_{B_{1}M_{1}}}}$": vgs_m1_sml_sweep,
        "$v_{i_{DS_{B_{1}M_{2}}}}$": ids_m2_sml_sweep,
        "$v_{i_{GS_{B_{1}M_{2}}}}$": vgs_m2_sml_sweep,
        "$i_{B_{1}R_{ref}}$":        irref_sml_sweep
    }

# ------------------------------------------------------------------------------
# Plotting Functions
# ------------------------------------------------------------------------------
def plot_tr_new(
    x,
    y,
    filename="plot.html",
    title="System Response",
    x_label="Time (s)",
    labels=None,
    quantity_types=None,
    extensions=None,
):
    import numpy as np
    import plotly.graph_objects as go

    # --- Ensure proper array shapes ---
    x = np.asarray(x)
    y = np.asarray(y)

    if y.ndim == 1:
        y = y.reshape(-1, 1)

    num_signals = y.shape[1]

    labels = labels or [f"Signal {i}" for i in range(num_signals)]
    quantity_types = quantity_types or ["voltage"] * num_signals

    context = {}

    # --- Run preprocess extensions ---
    if extensions:
        for ext in extensions:
            if getattr(ext, "preprocess", False):
                ext(None, x, y, labels, quantity_types, context=context)

    x_plot = context.get("x_crop", x)
    y_plot = context.get("y_crop", y)

    # --- Helper: compute padded limits ---
    def padded_limits(data, pad=0.05, include_zero=False):

        if data.size == 0:
            return [-1, 1]

        dmin = np.min(data)
        dmax = np.max(data)

        if include_zero:
            dmin = min(dmin, 0)
            dmax = max(dmax, 0)

        span = dmax - dmin

        if span == 0:
            span = abs(dmax) if dmax != 0 else 1

        return [
            dmin - pad * span,
            dmax + pad * span
        ]

    # --- Separate voltage and current signals safely ---
    v_idxs = [
        i for i, qt in enumerate(quantity_types[:y_plot.shape[1]])
        if qt.lower() == "voltage"
    ]

    i_idxs = [
        i for i, qt in enumerate(quantity_types[:y_plot.shape[1]])
        if qt.lower() == "current"
    ]

    v_range = padded_limits(y_plot[:, v_idxs]) if v_idxs else [-1, 1]
    i_range = padded_limits(y_plot[:, i_idxs]) if i_idxs else [-1, 1]

    # --- Create figure ---
    fig = go.Figure()

    for i in range(num_signals):

        is_voltage = quantity_types[i].lower() == "voltage"

        fig.add_trace(
            go.Scatter(
                x=x_plot,
                y=y_plot[:, i],
                name=labels[i],
                yaxis="y" if is_voltage else "y2",
                mode="lines",
            )
        )

    # --- Layout ---
    fig.update_layout(

        title=dict(
            text=title,
            x=0.5,
            xanchor="center"
        ),

        xaxis=dict(
            title=x_label,
            domain=[0.0, 0.85],
            showgrid=True,
        ),

        yaxis=dict(
            title="Voltage (V)",
            side="left",
            range=v_range,
            zeroline=True,
            zerolinecolor="gray",
            zerolinewidth=1,
        ),

        yaxis2=dict(
            title="Current (A)",
            side="right",
            overlaying="y",
            anchor="x",
            range=i_range,
            showgrid=False,
            zeroline=True,
            zerolinecolor="gray",
            zerolinewidth=1,
        ),

        legend=dict(
            x=1.02,
            y=1,
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="lightgray",
            borderwidth=1
        ),

        template="plotly_white",

        margin=dict(
            l=70,
            r=140,
            t=60,
            b=60
        ),
    )

    # --- Run post extensions ---
    if extensions:
        for ext in extensions:
            if not getattr(ext, "preprocess", False):
                ext(fig, x_plot, y_plot, labels, quantity_types, context=context)

    # --- Save ---
    fig.write_html(filename, include_mathjax="cdn")

    print(f"✅ Plot saved to {filename}")


def plot_sweep(vs_sweep, results, output_dir, suffix="ro"):
    """Generate interactive Plotly plots with LaTeX-rendered titles and legends."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    def save_plot(x, ys, labels, filename, ylabel, title):
        fig = go.Figure()

        for y, label in zip(ys, labels):
            hover_label = label.replace("$$", "").replace("{", "").replace("}", "")
            hover_label = hover_label.replace("_", "<sub>") 
            fig.add_trace(go.Scatter(
                x=x,
                y=y,
                mode="lines",
                name=label, 
                hovertemplate=f"{hover_label}: %{{y:.6g}}<extra></extra>"
            ))

        fig.update_layout(
            title=title,
            xaxis_title="$$V_{sweep}$$ (V)",
            yaxis_title=ylabel,
            template="plotly_white",
            hovermode="x unified",
            legend=dict(title="Legend", orientation="h", y=-0.25),
            font=dict(size=14)
        )

        file_path = output_path / filename
        html_str = fig.to_html(
            include_plotlyjs="cdn",
            full_html=True,
            include_mathjax="cdn",  
            config={"responsive": True}
        )
        with open(file_path, "w") as f:
            f.write(html_str)

        print(f"✅ Saved interactive Plotly figure: {file_path}")

    # Drain/source currents
    save_plot(
        vs_sweep,
        [results["isd_m3"], results["isd_m4"], results["$v_{i_{DS_{B_{1}M_{1}}}}$"], results["ids_m2"]],
        ["$$I_{SD,M3}$$", "$$I_{SD,M4}$$", "$$I_{DS,M1}$$", "$$I_{DS,M2}$$"],
        f"beta_multiplier_ds_{suffix}.html",
        "$$I$$ (A)",
        "Drain/Source Currents vs Sweep Voltage"
    )

    # Gate/source voltages
    save_plot(
        vs_sweep,
        [results["vsg_m3"], results["vsg_m4"], results["vgs_m1"], results["vgs_m2"]],
        ["$$V_{SG,M3}$$", "$$V_{SG,M4}$$", "$$V_{GS,M1}$$", "$$V_{GS,M2}$$"],
        f"beta_multiplier_gs_{suffix}.html",
        "$$V$$ (V)",
        "Gate/Source Voltages vs Sweep Voltage"
    )

    # Reference current
    save_plot(
        vs_sweep,
        [results["irref"]],
        ["$$I_{R_{ref}}$$"],
        f"beta_multiplier_iref_{suffix}.html",
        "$$I$$ (A)",
        "Reference Current vs Sweep Voltage"
    )

def plot_tr(tr, yr, output_dir, filename="plot.html", title="Plot of yr vs tr", x_label="tr", y_label="yr"):
    """Plots yr vs tr using Plotly and saves the plot as an HTML file."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    tr = np.array(tr).flatten()
    yr = np.array(yr).flatten()
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=tr, y=yr, mode='lines+markers', name='yr vs tr'))
    
    fig.update_layout(
        title=title,
        xaxis_title=x_label,
        yaxis_title=y_label,
        template="plotly_white"
    )
    
    file_path = output_path / filename
    fig.write_html(str(file_path))
    print(f"Plot saved as {file_path}")

def plot_bode_plotly(mag, phase, omega, output_dir, *,
                     single_out=True,
                     output_idx=(0, 0),
                     filename="bode_output.html"):
    """Plot magnitude and phase using Plotly and save to a standalone HTML file."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if single_out:
        mag_vec = mag
        phase_vec = phase
    else:
        o, i = output_idx
        mag_vec = mag[o][i]
        phase_vec = phase[o][i]

    mag_db = 20 * np.log10(mag_vec)
    phase_deg = np.degrees(phase_vec)

    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.08,
        subplot_titles=("Magnitude (dB)", "Phase (deg)")
    )

    fig.add_trace(go.Scatter(x=omega, y=mag_db, mode="lines", name="Magnitude"), row=1, col=1)
    fig.add_trace(go.Scatter(x=omega, y=phase_deg, mode="lines", name="Phase"), row=2, col=1)

    fig.update_xaxes(type="log", title_text="Frequency (rad/s)", row=1, col=1)
    fig.update_xaxes(type="log", title_text="Frequency (rad/s)", row=2, col=1)
    fig.update_yaxes(title_text="Magnitude (dB)", row=1, col=1)
    fig.update_yaxes(title_text="Phase (deg)", row=2, col=1)

    fig.update_layout(
        height=700,
        width=900,
        showlegend=False,
        title_text=f"Bode Plot ({filename})"
    )

    file_path = output_path / filename
    fig.write_html(str(file_path), include_plotlyjs="cdn")
    print(f"Plot saved as {file_path}")

def plot_pz_plotly(sys, poles, zeros, output_dir, filename="pz_plot.html"):
    """Plot Pole-Zero map and save to standalone HTML file."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=np.real(poles), y=np.imag(poles),
        mode="markers", marker=dict(symbol="x", size=10), name="Poles"
    ))

    fig.add_trace(go.Scatter(
        x=np.real(zeros), y=np.imag(zeros),
        mode="markers", marker=dict(symbol="circle", size=10), name="Zeros"
    ))

    fig.add_shape(
        type="line",
        x0=min(np.real(poles.min()), -1) if len(poles) else -1,
        x1=max(np.real(poles.max()), 1) if len(poles) else 1,
        y0=0, y1=0, line=dict(dash="dash", color="gray")
    )

    fig.add_shape(
        type="line",
        x0=0, x1=0,
        y0=min(np.imag(poles.min()), -1) if len(poles) else -1,
        y1=max(np.imag(poles.max()), 1) if len(poles) else 1,
        line=dict(dash="dash", color="gray")
    )

    fig.update_layout(
        title="Pole–Zero Map",
        xaxis_title="Real Axis (σ)",
        yaxis_title="Imaginary Axis (jω)",
        width=700,
        height=700,
        showlegend=True
    )

    file_path = output_path / filename
    fig.write_html(str(file_path), include_plotlyjs="cdn")
    print(f"Plot saved as {file_path}")


def quarto_math_block(equations):
    r"""
    Convert an array of equation strings into a Quarto-ready MathJax aligned block.

    Example output:

    $$
    \begin{aligned}
    a &= b + c \\
    d &= e + f
    \end{aligned}
    $$
    """

    if not equations:
        return "$$\n\\begin{aligned}\n\\end{aligned}\n$$"

    lines = ["$$", r"\begin{aligned}"]

    # Add equations with proper line breaks except the last
    for i, eq in enumerate(equations):
        if i < len(equations) - 1:
            lines.append(eq + r" \\")
        else:
            lines.append(eq)

    lines.append(r"\end{aligned}")
    lines.append("$$")

    return "\n".join(lines)



# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():
    # Define output directory once and ensure it exists
    output_dir = Path(os.environ.get(
        "BETA_MULTIPLIER_OUTPUT_DIR", "../../doc/beta_multiplier"
    ))
    output_dir.mkdir(parents=True, exist_ok=True)

    vs_sweep = np.linspace(*SWEEP_RANGE)
    m2_sweep = np.linspace(*M2_SWEEP_RANGE)
    lrg_vs_sweep = np.linspace(*LRG_SWEEP_RANGE)
    Rref_sweep = np.linspace(*RREF_SWEEP_RANGE)
    params = {"Rref": DEFAULT_RREF, "Vs": DEFAULT_VS, "Cstray": 100e-12}

    fstep    = 10000
    fstart   = 2 * math.pi * 1e3
    fstop    = 2 * math.pi * 100e6
    omega    = np.linspace(fstart, fstop, fstep)

    nw_m2             = beta_multiplier_ntwrk(output_resistance=True,  res_m2=True)
    nw_m2_ss          = beta_multiplier_ntwrk(output_resistance=True,  res_m2=True)

    nw_no_ro_m2       = beta_multiplier_ntwrk(output_resistance=False, res_m2=True)
    nw_no_ro_m2_stray = beta_multiplier_ntwrk(output_resistance=True,  res_m2=True, stray_capacitance=True)
    nw_m1             = beta_multiplier_ntwrk(output_resistance=True,  res_m2=False)
    nw_no_ro_m1       = beta_multiplier_ntwrk(output_resistance=False, res_m2=False)

    eqs = nw_no_ro_m2.ckt.build_mathjax_equations()
    math_block = quarto_math_block(eqs)
    print(math_block)

    # Solve OPs and build small-signal circuits
    sweep_cases = [
        (nw_m2, "ro_m2", None),
        (nw_no_ro_m2, "no_ro_m2", None),
        # Moving the degeneration resistor creates another separated solution
        # basin.  Continue from the corresponding M2 topology at the same
        # first supply point instead of restarting from an artificial vector.
        (nw_m1, "ro_m1", nw_m2),
        (nw_no_ro_m1, "no_ro_m1", nw_no_ro_m2),
    ]
    for net, label, topology_seed in sweep_cases:
        initial = (None if topology_seed is None else
                   topology_seed._sweep_first_solution)
        results = sweep(net, lrg_vs_sweep, root_start=initial)
        # Pass output_dir directly
        #plot_sweep(lrg_vs_sweep, results, output_dir, suffix=label)

        cycle_filename = output_dir / f"beta_multiplier_{label}.html"

        # Labels (dict keys)
        labels = list(results.keys())

        # Stack signals into 2D array
        yr = np.column_stack([
            np.asarray(results[k]).reshape(-1)
            for k in labels
        ])

        # generate quantity types FROM labels (guaranteed same length)
        quantity_types = [
            "current" if label.startswith("$i_") else "voltage"
            for label in labels
        ]

        plot_tr_new(
            lrg_vs_sweep,
            yr,
            labels=labels,
            quantity_types=quantity_types,
            x_label="$V_{DD}$",
            filename=str(cycle_filename),
            extensions=[]
        )

    # Solve OPs and build small-signal circuits
    for net, base_label in [(nw_m2, "ro_m2_sml")]:
            
        # Update parameters (including Rref)
        net.set_beta_multiplier_params(**params)

        # Solve operating point
        root_start = np.zeros(net.ckt.num_edges)
        x = op_solve(net.ckt, root_start, branch_ref=net.get_rref_idx())
        op = net.ckt.scb

        # Build small-signal circuit with loop break
        net.add_sml_ckt(op=op, gate_topology=net.gate_break_m1_m2())

        # Extract operating-point small-signal parameters
        gm_m1, ro_m1 = net.beta_mult.nmos_m1.get_op_sml(op)
        gm_m2, ro_m2 = net.beta_mult.nmos_m2.get_op_sml(op)
        gm_m3, ro_m3 = net.beta_mult.pmos_m3.get_op_sml(op)
        gm_m4, ro_m4 = net.beta_mult.pmos_m4.get_op_sml(op)

        print(f"m1 gm = {gm_m1:.6f}, ro = {ro_m1:.6f}")
        print(f"m2 gm = {gm_m2:.6f}, ro = {ro_m2:.6f}")
        print(f"m3 gm = {gm_m3:.6f}, ro = {ro_m3:.6f}")
        print(f"m4 gm = {gm_m4:.6f}, ro = {ro_m4:.6f}")

        # Run small-signal sweep
        results = sml_sweep_m2(net, m2_sweep)

        label = f"{base_label}"
        
        # Pass output_dir directly
        #plot_sweep(m2_sweep, results, output_dir, suffix=label)

        cycle_filename = output_dir / f"beta_multiplier_open_loop_{label}.html"

        # Labels (dict keys)
        labels = list(results.keys())

        # Stack signals into 2D array
        yr = np.column_stack([
            np.asarray(results[k]).reshape(-1)
            for k in labels
        ])

        # generate quantity types FROM labels (guaranteed same length)
        quantity_types = [
            "current" if label.startswith("$i_") else "voltage"
            for label in labels
        ]

        plot_tr_new(
            m2_sweep,
            yr,
            labels=labels,
            quantity_types=quantity_types,
            x_label="$V_{test}$",
            filename=str(cycle_filename),
            extensions=[]
        )

    # Solve OPs and build SS small-signal circuits
    for net, base_label in [(nw_m2_ss, "ro_m2_sml")]:

        # -----------------------------
        # Parameter setup (OP only once)
        # -----------------------------
        Rref = 6.5e3

        params_with_rref = dict(params)
        params_with_rref["Rref"] = Rref
        net.set_beta_multiplier_params(**params_with_rref)

        # Solve operating point
        root_start = np.zeros(net.ckt.num_edges)
        x = op_solve(net.ckt, root_start, branch_ref=net.get_rref_idx())
        op = net.ckt.scb

        # -----------------------------
        # Gate-topology configurations
        # -----------------------------
        gate_configs = [
            {
                "label": "open_loop",
                "gate_topology": net.gate_break_m1_m2()
            },
            {
                "label": "closed_loop",
                "gate_topology": net.gate_break_m1_m2_closed()
            }
        ]

        # -----------------------------
        # Small-signal analyses
        # -----------------------------
        for cfg in gate_configs:

            label = f"{base_label}_{cfg['label']}"

            # Build small-signal circuit
            net.add_sml_ckt(
                op=op,
                lti=True,
                gate_topology=cfg["gate_topology"],
                degen_topology=net._capacitive_degen_topology(),
                output_topology=net._output_vds_m2_topology()
            )

            # -----------------------------
            # State-space construction
            # -----------------------------
            net.ckt_sml.get_ss()
            sys = ct.ss(
                net.ckt_sml.A,
                net.ckt_sml.B,
                net.ckt_sml.C,
                net.ckt_sml.D
            )

            sys_scypi = StateSpace(net.ckt_sml.A, net.ckt_sml.B, net.ckt_sml.C, net.ckt_sml.D)

            w, H = freqresp(sys_scypi, omega)

            mag_scypi = np.abs(H)
            phase_scypi = np.angle(H, deg=True)

            input_idx = net.get_source_idx("gate_m2")

            # Force SISO: output 0, input 0
            sys_siso = sys[0, input_idx]

            # Poles and zeros
            poles = ct.pole(sys_siso)
            zeros = ct.zero(sys_siso)

            # Frequency response
            mag, phase, omega_out = ct.freqresp(sys_siso, omega=omega)

            # -----------------------------
            # Plotting
            # -----------------------------
            
            # Pass output_dir directly
            plot_bode_plotly(
               mag,
               phase,
               omega_out,
               output_dir,
               single_out=True,
               output_idx=(0, 0),
               filename=f"bode_{label}.html"
            )

            plot_pz_plotly(
                sys_siso,
                poles=poles,
                zeros=zeros,
                output_dir=output_dir,
                filename=f"pz_{label}.html"
            )

    exit(1)

if __name__ == "__main__":
    main()
