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
from parallel_lrc_ntwrk import *

import control as ct
from scipy.signal import freqresp, StateSpace

from pathlib import Path

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
DEFAULT_R = 31.6e3
DEFAULT_C = 1e-9
DEFAULT_L = 10e-3
DEFAULT_VDD = 5.0

DEFAULT_RREF = 6.5e3
DEFAULT_VS = 5.0
SWEEP_RANGE = (-0.2, 0.2, 100)
LRG_SWEEP_RANGE = (4, 6, int(1000))
RREF_SWEEP_RANGE = (2.5e3, 6.5e3, int(10))
TEND = 50e-12
TSTEP = 10000
TOL = 1e-6
M2_SWEEP_RANGE = (0, 1, 100)

# ------------------------------------------------------------------------------
# Core circuit equations
# ------------------------------------------------------------------------------

def circuit_eqn(x, sys, ckt, t):
    """Return Kirchhoff equations for given circuit state x."""
    I = ckt.get_im(x=x, sys=sys, t=t)
    V = ckt.get_vm(x=x, sys=sys, t=t)
    I = ckt.get_im_dep()
    qf_num, bf_num = ckt.qf.copy(), ckt.bf.copy()
    return np.dot(qf_num, I) + np.dot(bf_num, V)


def op_solve(ckt, root_start):
    """Find DC operating point."""
    root = optimize.root(circuit_eqn, root_start, args=(0, ckt, 0), tol=1e-9, method='lm')
    if not root.success:
        print("⚠️ OP solver did not converge, retrying with updated initial guess")
        return op_solve(ckt, root.x)
    return root.x

def dypc_litmus0(t, sys, ckt):

    root_start = np.ones(ckt.num_edges)
    vs = 1

    root = optimize.root(circuit_eqn, root_start, args=(sys, ckt, t), tol=1e-9, method='lm')
    root_start = root.x
    if not root.success:
        print (root.x)
        sys.exit(0)

    return ckt.get_dy(x=root.x)

def ode_solve(ckt, tend=50e-12, tstep=10000, x0=0):

    num_sys_vars    = ckt.get_num_sys_vars()

    t     = np.linspace(0, tend, tstep)

    r = ode(dypc_litmus0).set_integrator('lsoda', method='bdf', atol=1e-9, rtol=1e-9)
    r.set_initial_value(x0, 0.0)
    r.set_f_params(ckt)

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

def sweep(nw, vs_sweep):
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

    root_start = np.ones(nw.ckt.num_edges)

    for vs_idx, vs in enumerate(vs_sweep):
        params = {"Rref": 6.5e3, "Vs": vs}
        nw.set_beta_multiplier_params(**params)

        x = op_solve(nw.ckt, root_start)
        root_start = x.copy()

        #rref_idx = nw.get_rref_idx()
        #root_start[rref_idx] = 20e-3

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
        "isd_m3": isd_m3_sweep,
        "vsg_m3": vsg_m3_sweep,
        "isd_m4": isd_m4_sweep,
        "vsg_m4": vsg_m4_sweep,
        "ids_m1": ids_m1_sweep,
        "vgs_m1": vgs_m1_sweep,
        "ids_m2": ids_m2_sweep,
        "vgs_m2": vgs_m2_sweep,
        "irref":  irref_sweep
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

    for vs_idx, vs in enumerate(vs_sml_sweep):
        params = {"Rref": 6.5e3, "Vs": vs}
        nw.set_beta_multiplier_sml_params(**params)

        root_start = np.ones(nw.ckt_sml.num_edges)
        x = op_solve(nw.ckt_sml, root_start)

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
        "isd_m3": isd_m3_sml_sweep,
        "vsg_m3": vsg_m3_sml_sweep,
        "isd_m4": isd_m4_sml_sweep,
        "vsg_m4": vsg_m4_sml_sweep,
        "ids_m1": ids_m1_sml_sweep,
        "vgs_m1": vgs_m1_sml_sweep,
        "ids_m2": ids_m2_sml_sweep,
        "vgs_m2": vgs_m2_sml_sweep,
        "irref":  irref_sml_sweep
    }

def sml_sweep_m2(nw, vs_sml_sweep, Rref):
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

    for vs_idx, vs in enumerate(vs_sml_sweep):
        nw.set_source(source="gate_m2", val=vs)
        params = {"Rref": Rref, "Vs": 0}
        nw.set_beta_multiplier_sml_params(**params)
        nw.set_source("gate_m2", vs)

        root_start = np.ones(nw.ckt_sml.num_edges)
        x = op_solve(nw.ckt_sml, root_start)

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
        "isd_m3": isd_m3_sml_sweep,
        "vsg_m3": vsg_m3_sml_sweep,
        "isd_m4": isd_m4_sml_sweep,
        "vsg_m4": vsg_m4_sml_sweep,
        "ids_m1": ids_m1_sml_sweep,
        "vgs_m1": vgs_m1_sml_sweep,
        "ids_m2": ids_m2_sml_sweep,
        "vgs_m2": vgs_m2_sml_sweep,
        "irref":  irref_sml_sweep
    }



import os
import plotly.graph_objects as go

def plot_sweep(vs_sweep, results, suffix="ro", output_dir="../../doc"):
    """Generate interactive Plotly plots with LaTeX-rendered titles and legends,
    and hover labels formatted with subscripts using HTML."""
    os.makedirs(output_dir, exist_ok=True)

    def save_plot(x, ys, labels, filename, ylabel, title):
        fig = go.Figure()

        for y, label in zip(ys, labels):
            # Convert LaTeX labels to simple HTML for hover (remove $$, replace _ with <sub>)
            hover_label = label.replace("$$", "").replace("{", "").replace("}", "")
            hover_label = hover_label.replace("_", "<sub>")  # simple subscript
            fig.add_trace(go.Scatter(
                x=x,
                y=y,
                mode="lines",
                name=label,  # legend still renders LaTeX
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

        # Save HTML with MathJax for titles/legend
        file_path = os.path.join(output_dir, filename)
        html_str = fig.to_html(
            include_plotlyjs="cdn",
            full_html=True,
            include_mathjax="cdn",  # ensures LaTeX in titles and legend
            config={"responsive": True}
        )
        with open(file_path, "w") as f:
            f.write(html_str)

        print(f"✅ Saved interactive Plotly figure: {file_path}")

    # Drain/source currents
    save_plot(
        vs_sweep,
        [results["isd_m3"], results["isd_m4"], results["ids_m1"], results["ids_m2"]],
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

import numpy as np
import plotly.graph_objects as go
from scipy.signal import find_peaks

def ext_crop_cycles(fig, tr, yr, labels, quantity_types, cycles=5.5, skip_cycles=2, context=None):
    """PREPROCESS: Crops exactly 2.5 cycles, starting cleanly on a zero-crossing."""
    y = yr[:, 0]
    signs = np.sign(y)
    rising = np.where(np.diff(signs) > 0)[0]
    
    if len(rising) < skip_cycles + 1: return context
    
    period_samples = np.mean(np.diff(rising))
    start_idx = int(rising[skip_cycles])
    end_idx = int(start_idx + (cycles * period_samples))
    
    mask = slice(max(0, start_idx), min(len(tr)-1, end_idx))
    
    if context is not None:
        context["tr_crop"] = tr[mask]
        context["yr_crop"] = yr[mask, :]
        context["fundamental_period"] = (tr[1] - tr[0]) * period_samples
    return context

def format_engineering(value, unit="C"):
    """Formats a value into engineering notation with MathJax SI prefixes."""
    prefixes = {
        0:  "",
        -3: "m",
        -6: "\\mu ", # LaTeX for micro
        -9: "n",
        -12:"p"
    }
    
    # Calculate the exponent in steps of 3
    if value == 0:
        return f"0\\,{unit}"
    
    exp = int(np.floor(np.log10(abs(value)) / 3.0) * 3)
    # Clamp to supported prefixes
    exp = max(-12, min(0, exp))
    
    scaled_val = value / (10**exp)
    prefix = prefixes.get(exp, "")
    
    return f"{scaled_val:.2f}\\,{prefix}{unit}"

def ext_integrate_current(fig, tr, yr, labels, quantity_types, context=None):
    """POSTPROCESS: Targets the positive lobe 1 cycle before the center."""
    current_indices = [i for i, qt in enumerate(quantity_types) if qt == "current"]
    if not current_indices: return context

    y = yr[:, current_indices[0]]
    period = context.get("fundamental_period", (tr[-1]-tr[0])/2.5)
    mid_time = tr[len(tr)//2]
    target_t = mid_time - period 
    
    signs = np.sign(y)
    rising = np.where(np.diff(signs) > 0)[0]
    falling = np.where(np.diff(signs) < 0)[0]
    
    lobes = []
    for r_idx in rising:
        f_idxs = falling[falling > r_idx]
        if len(f_idxs) > 0:
            lobes.append((r_idx, f_idxs[0]))
            
    if not lobes: return context
    
    lobe_centers = [tr[(s+e)//2] for s, e in lobes]
    best_lobe = np.argmin(np.abs(np.array(lobe_centers) - target_t))
    
    idx_s, idx_e = lobes[best_lobe]
    tr_int, y_int = tr[idx_s:idx_e], y[idx_s:idx_e]
    
    # Calculate integral
    val = np.trapz(y_int, tr_int)
    # Format with our new engineering helper
    formatted_val = format_engineering(val, "C")
    
    fig.add_trace(go.Scatter(
        x=np.concatenate([tr_int, tr_int[::-1]]),
        y=np.concatenate([y_int, np.zeros_like(y_int)]),
        fill="toself", fillcolor="rgba(0, 200, 100, 0.25)",
        line=dict(width=0), yaxis="y2", showlegend=False, hoverinfo="skip"
    ))

    fig.add_annotation(
        x=np.mean(tr_int), y=np.max(y_int), yref="y2",
        # Use the formatted string in MathJax
        text=f"$\\int i(t)dt = {formatted_val}$", 
        showarrow=True, arrowhead=2, ay=-35, bgcolor="white"
    )
    return context

def ext_phase_voltage(fig, tr, yr, labels, quantity_types, context=None):
    """POSTPROCESS: Anchors to the absolute CENTER and draws a horizontal arrow."""
    v_idxs = [i for i, qt in enumerate(quantity_types) if qt == "voltage"]
    if len(v_idxs) < 2: return context

    y1, y2 = yr[:, v_idxs[0]], yr[:, v_idxs[1]]
    period = context.get("fundamental_period", (tr[-1]-tr[0])/2.5)
    mid_time = tr[len(tr)//2]
    
    pks1, _ = find_peaks(y1, distance=len(y1)//5, prominence=np.max(y1)*0.2)
    pks2, _ = find_peaks(y2, distance=len(y2)//5, prominence=np.max(y2)*0.2)

    if len(pks1) == 0 or len(pks2) == 0: return context

    t1 = tr[pks1[np.argmin(np.abs(tr[pks1] - mid_time))]]
    t2 = tr[pks2[np.argmin(np.abs(tr[pks2] - t1))]]
    
    dt = t2 - t1
    phase = -(dt / period) * 360
    phase = ((phase + 180) % 360) - 180 

    ymax = max(np.max(y1), np.max(y2))
    
    # Vertical dashed lines
    fig.add_vline(x=t1, line_dash="dash", line_color="black", opacity=0.5)
    fig.add_vline(x=t2, line_dash="dash", line_color="red", opacity=0.5)

    # 1. NEW: The Horizontal Double-Headed Arrow connecting t1 and t2
    fig.add_annotation(
        x=t2, y=ymax,           # Arrow Head
        ax=t1, ay=ymax,         # Arrow Tail
        xref="x", yref="y", axref="x", ayref="y",
        showarrow=True, arrowhead=2, arrowside="start+end",
        arrowwidth=1.5, arrowcolor="black"
    )

    # 2. NEW: The Text floating just above the arrow
    fig.add_annotation(
        x=(t1+t2)/2, y=ymax, 
        text=f"$\\Delta\\phi = {phase:.1f}^\\circ$",
        showarrow=False, # We don't need a pointer line for the text anymore
        yshift=15,       # Pushes the text 15 pixels UP from the arrow
        bgcolor="rgba(255, 255, 255, 0.8)", bordercolor="black"
    )
    return context

def make_extrema_ext(signal_index=0, prominence=0.0):
    def ext_mark_extrema(fig, tr, yr, labels, quantity_types,
                        signal_index=signal_index,
                        prominence=prominence,
                        context=None):
        """
        Draw vertical dashed lines at peaks and troughs of a selected signal.
        
        Parameters:
        - signal_index: which column of yr to analyze
        - prominence: minimum peak prominence (filters noise)
        """
        import numpy as np
        from scipy.signal import find_peaks

        y = yr[:, signal_index]

        # Find peaks
        peaks, _ = find_peaks(y, prominence=prominence)

        # Find troughs (invert signal)
        troughs, _ = find_peaks(-y, prominence=prominence)

        extrema = np.sort(np.concatenate([peaks, troughs]))

        # Add vertical lines
        for idx in extrema:
            fig.add_vline(
                x=tr[idx],
                line=dict(color="gray", width=1, dash="dash"),
                opacity=0.5
            )
    return ext_mark_extrema

def ext_color_abs_slope(
    fig, tr, yr, labels, quantity_types,
    signal_index=0,
    color_increasing="rgba(0, 200, 0, 0.25)",
    color_decreasing="rgba(200, 0, 0, 0.25)",
    context=None
):
    import numpy as np
    import plotly.graph_objects as go

    y = yr[:, signal_index]
    abs_y = np.abs(y)

    # Compute |y| slope
    d_abs = np.diff(abs_y)
    state = np.sign(d_abs)
    for i in range(1, len(state)):
        if state[i] == 0:
            state[i] = state[i-1]
    if len(state) > 0 and state[0] == 0:
        state[0] = 1
    state = np.concatenate([[state[0]], state])

    # Segment the waveform
    segments = []
    start = 0
    for i in range(1, len(state)):
        if state[i] != state[i-1]:
            segments.append((start, i))
            start = i
    segments.append((start, len(state)))

    # Add filled polygons
    for s, e in segments:
        if e - s < 2:
            continue
        tr_seg = tr[s:e]
        y_seg = y[s:e]
        color = color_increasing if state[s] > 0 else color_decreasing

        # Build closed polygon
        x_poly = np.concatenate([tr_seg, tr_seg[::-1]])
        y_poly = np.concatenate([y_seg, np.zeros_like(y_seg)])

        fig.add_trace(go.Scatter(
            x=x_poly,
            y=y_poly,
            mode="lines",
            fill="toself",
            fillcolor=color,
            line=dict(width=0),
            yaxis="y2",
            showlegend=False,
            hoverinfo="skip"
        ))

    return context

def make_abs_slope_ext(
    signal_index=0,
    color_increasing="rgba(0, 200, 0, 0.25)",
    color_decreasing="rgba(200, 0, 0, 0.25)"
):
    def ext(fig, tr, yr, labels, quantity_types, context=None):
        return ext_color_abs_slope(
            fig, tr, yr, labels, quantity_types,
            signal_index=signal_index,
            color_increasing=color_increasing,
            color_decreasing=color_decreasing,
            context=context
        )
    return ext

ext_crop_cycles.preprocess = True
ext_integrate_current.preprocess = False
ext_phase_voltage.preprocess = False

def plot_tr(
    tr, yr,
    filename="plot.html",
    title="System Response",
    labels=None,
    quantity_types=None,
    extensions=None,
    v_bound=None,   # NEW: allow fixed bounds
    i_bound=None    # NEW: allow fixed bounds
):
    import numpy as np
    import plotly.graph_objects as go

    tr, yr = np.asarray(tr), np.asarray(yr)
    if yr.ndim == 1:
        yr = yr.reshape(-1, 1)

    num_signals = yr.shape[1]
    labels = labels or [f"Sig {i}" for i in range(num_signals)]
    quantity_types = quantity_types or ["voltage"] * num_signals

    context = {}

    # --- Run preprocess extensions ---
    if extensions:
        for ext in extensions:
            if getattr(ext, "preprocess", False):
                ext(None, tr, yr, labels, quantity_types, context=context)

    # --- Use cropped data if provided ---
    tr_p = context.get("tr_crop", tr)
    yr_p = context.get("yr_crop", yr)

    # --- Identify signal types ---
    v_idxs = [i for i, qt in enumerate(quantity_types) if qt.lower() == "voltage"]
    i_idxs = [i for i, qt in enumerate(quantity_types) if qt.lower() == "current"]

    # --- Compute bounds ONLY if not provided ---
    if v_bound is None:
        v_bound = np.max(np.abs(yr[:, v_idxs])) * 1.1 if v_idxs else 1

    if i_bound is None:
        i_bound = np.max(np.abs(yr[:, i_idxs])) * 1.1 if i_idxs else 1

    # --- Create figure ---
    fig = go.Figure()

    for i in range(num_signals):
        is_v = quantity_types[i].lower() == "voltage"
        fig.add_trace(go.Scatter(
            x=tr_p,
            y=yr_p[:, i],
            name=labels[i],
            yaxis="y" if is_v else "y2"
        ))

    # --- Layout with FIXED scaling ---
    fig.update_layout(
        title=title,
        xaxis_title="Time (s)",
        yaxis=dict(
            title="Voltage (V)",
            side="left",
            range=[-v_bound, v_bound],
            zeroline=True,
            zerolinecolor="lightgray",
            zerolinewidth=1
        ),
        yaxis2=dict(
            title="Current (A)",
            side="right",
            overlaying="y",
            showgrid=False,
            range=[-i_bound, i_bound],
            zeroline=True,
            zerolinecolor="lightgray",
            zerolinewidth=1
        ),
        template="plotly_white",
        legend=dict(orientation="h", y=1.1, x=1)
    )

    # --- Run postprocess extensions ---
    if extensions:
        for ext in extensions:
            if not getattr(ext, "preprocess", False):
                ext(fig, tr_p, yr_p, labels, quantity_types, context=context)

    fig.write_html(filename, include_mathjax="cdn")
    print(f"Plot saved to {filename}")

def ext_lc_visualization(
    fig, tr, yr, labels, quantity_types,
    signal_index=0,               # Which signal to analyze
    color_increasing="rgba(0, 255, 0, 0.3)",
    color_decreasing="rgba(255, 0, 0, 0.3)",
    draw_extrema=True,            # Whether to draw vertical lines
    prominence=0.0,               # Peak/trough detection threshold
    context=None
):
    """
    Full LC visualization extension:
    - Colors under waveform based on |y| increasing/decreasing
    - Draws vertical dashed lines at peaks and troughs
    """

    import numpy as np
    from scipy.signal import find_peaks
    import plotly.graph_objects as go

    y = yr[:, signal_index]
    abs_y = np.abs(y)

    # 1️⃣ Color under curve based on |y| slope
    d_abs = np.diff(abs_y)
    state = np.sign(d_abs)
    for i in range(1, len(state)):
        if state[i] == 0:
            state[i] = state[i-1]
    if len(state) > 0 and state[0] == 0:
        state[0] = 1
    state = np.concatenate([[state[0]], state])

    # Segment continuous slope regions
    segments = []
    start = 0
    for i in range(1, len(state)):
        if state[i] != state[i-1]:
            segments.append((start, i))
            start = i
    segments.append((start, len(state)))

    # Draw filled polygons under curve
    for s, e in segments:
        if e - s < 2:
            continue
        tr_seg = tr[s:e]
        y_seg = y[s:e]
        color = color_increasing if state[s] > 0 else color_decreasing

        # Closed polygon
        x_poly = np.concatenate([tr_seg, tr_seg[::-1]])
        y_poly = np.concatenate([y_seg, np.zeros_like(y_seg)])

        fig.add_trace(go.Scatter(
            x=x_poly,
            y=y_poly,
            mode="lines",
            fill="toself",
            fillcolor=color,
            line=dict(width=0),
            yaxis="y2",
            showlegend=False,
            hoverinfo="skip"
        ))

    # 2️⃣ Draw vertical dashed lines at peaks & troughs
    if draw_extrema:
        peaks, _ = find_peaks(y, prominence=prominence)
        troughs, _ = find_peaks(-y, prominence=prominence)
        extrema = np.sort(np.concatenate([peaks, troughs]))
        for idx in extrema:
            fig.add_vline(
                x=tr[idx],
                line=dict(color="gray", width=1, dash="dash"),
                opacity=0.5
            )

    return context

def make_lc_visualization_ext(
    signal_index=0,
    color_increasing="rgba(0, 255, 0, 0.3)",
    color_decreasing="rgba(255, 0, 0, 0.3)",
    draw_extrema=True,
    prominence=0.0
):
    """
    Factory wrapper for ext_lc_visualization.
    Returns a configured extension function compatible with your plot_tr pipeline.
    """
    def ext(fig, tr, yr, labels, quantity_types, context=None):
        return ext_lc_visualization(
            fig, tr, yr, labels, quantity_types,
            signal_index=signal_index,
            color_increasing=color_increasing,
            color_decreasing=color_decreasing,
            draw_extrema=draw_extrema,
            prominence=prominence,
            context=context
        )
    return ext


def quarto_math_block(equations):
    """
    Convert an array of equation strings into a Quarto-ready MathJax block.
    """

    lines = ["$$"]

    for eq in equations:
        lines.append(eq + r" \\")  # new line in MathJax

    # Remove trailing \\ from last equation
    if len(lines) > 1:
        lines[-1] = lines[-1].rstrip(r" \\")

    lines.append("$$")

    return "\n".join(lines)

def plot_bode_plotly(mag, phase, omega, *,
                     single_out=True,
                     output_idx=(0, 0),
                     filename="bode_output.html"):
    """
    Plot magnitude and phase using Plotly and save to a standalone HTML file.

    Parameters
    ----------
    mag : ndarray
        Magnitude array from ct.freqresp
    phase : ndarray
        Phase array from ct.freqresp (radians)
    omega : ndarray
        Frequency vector (rad/s)
    single_out : bool
        True for SISO systems
    output_idx : tuple
        (output, input) index for MIMO systems
    filename : str
        Output HTML filename
    """

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

    # Magnitude plot
    fig.add_trace(
        go.Scatter(
            x=omega,
            y=mag_db,
            mode="lines",
            name="Magnitude"
        ),
        row=1, col=1
    )

    # Phase plot
    fig.add_trace(
        go.Scatter(
            x=omega,
            y=phase_deg,
            mode="lines",
            name="Phase"
        ),
        row=2, col=1
    )

    fig.update_xaxes(
        type="log",
        title_text="Frequency (rad/s)",
        row=1, col=1
    )

    fig.update_xaxes(
        type="log",
        title_text="Frequency (rad/s)",
        row=2, col=1
    )

    fig.update_yaxes(
        title_text="Magnitude (dB)",
        row=1, col=1
    )

    fig.update_yaxes(
        title_text="Phase (deg)",
        row=2, col=1
    )

    fig.update_layout(
        height=700,
        width=900,
        showlegend=False,
        title_text=f"Bode Plot ({filename})"
    )

    fig.write_html(filename, include_plotlyjs="cdn")

def plot_pz_plotly(sys, poles, zeros, filename="pz_plot.html"):

    fig = go.Figure()

    # Poles
    fig.add_trace(go.Scatter(
        x=np.real(poles),
        y=np.imag(poles),
        mode="markers",
        marker=dict(symbol="x", size=10),
        name="Poles"
    ))

    # Zeros
    fig.add_trace(go.Scatter(
        x=np.real(zeros),
        y=np.imag(zeros),
        mode="markers",
        marker=dict(symbol="circle", size=10),
        name="Zeros"
    ))

    # Axes
    fig.add_shape(
        type="line",
        x0=min(np.real(poles.min()), -1),
        x1=max(np.real(poles.max()), 1),
        y0=0, y1=0,
        line=dict(dash="dash", color="gray")
    )

    fig.add_shape(
        type="line",
        x0=0, x1=0,
        y0=min(np.imag(poles.min()), -1),
        y1=max(np.imag(poles.max()), 1),
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

    fig.write_html(filename, include_plotlyjs="cdn")

NUM_CYCLES = 25
DEFAULT_MAG = 100e-6

def run_rc_transient(net, omega, num_cycles=NUM_CYCLES, mag=DEFAULT_MAG, tstep_per_cycle=100):
    """
    Run transient simulation for a single omega.

    Automatically sets simulation duration to cover ~num_cycles cycles.
    """

    # Derived values
    f = omega / (2 * math.pi)
    period = 1 / f
    tend = num_cycles * period

    # Choose timestep so each cycle has sufficient resolution
    tstep = int(num_cycles * tstep_per_cycle)

    # Configure source
    net.set_sine_input_source(omega=omega, mag=mag)

    # Initial conditions
    x0 = np.zeros(net.ckt.get_num_sys_vars())

    # Run solver
    tr, yr = ode_solve(net.ckt, tend=tend, tstep=tstep, x0=x0)

    # Ensure yr is 2D with correct orientation
    yr = np.atleast_2d(yr)
    if yr.shape[0] != len(tr):
        yr = yr.T

    # Input waveform
    i_input = net.get_sine_time_series_input(tr=tr)

    # Resistor voltage and current
    v_res = yr[:, 0]
    i_res = v_res / DEFAULT_R

    # Inductor current
    i_ind = yr[:, 1]

    # Capacitor current
    i_cap = i_input - i_res - i_ind

    # DIfference
    i_diff = i_cap + i_ind

    # Append signals
    yr = np.column_stack((yr, i_input, i_res, i_cap, i_diff))

    return tr, yr


def plot_rc_results(tr, yr, omega, prefix="parallel_lrc_m1", base_path="."):
    """
    Generate standard plots and save to a specified directory.
    """
    # 1. Ensure the directory exists
    # parents=True creates intermediate folders if missing
    # exist_ok=True prevents errors if the folder already exists
    output_dir = Path(base_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    label = f"{prefix}_{omega}"

    # 2. Construct full file paths
    transient_filename = output_dir / f"transient_{label}.html"
    cycle_filename = output_dir / f"transient_cycle_{label}.html"

    # Compute global bounds ONCE
    v_idxs = [0]  # or however you define voltage indices
    i_idxs = [1,2,3,4,5]

    v_bound = np.max(np.abs(yr[:, v_idxs])) * 1.1
    i_bound = np.max(np.abs(yr[:, i_idxs])) * 1.1

    # 3. Standard transient plot
    plot_tr(
        tr=tr,
        yr=yr,
        v_bound=v_bound,
        i_bound=i_bound,
        labels=[r"$V_C$", r"$i_{L}$", r"$i_{IN}$", r"$i_{R}$", r"$i_{C}$", r"$i_{diff}$"],
        quantity_types=["voltage", "current", "current", "current", "current", "current"],
        filename=str(transient_filename) # Passing as string for compatibility
    )

    # 4. Cycle-processed plot
    plot_tr(
        tr,
        yr,
        v_bound=v_bound,
        i_bound=i_bound,
        labels=[r"$V_C$", r"$i_{L}$", r"$i_{IN}$", r"$i_{R}$", r"$i_{C}$", r"$i_{diff}$"],
        quantity_types=["voltage", "current", "current", "current", "current", "current"],
        filename=str(cycle_filename),
        extensions=[
            ext_crop_cycles,

            make_lc_visualization_ext(
                signal_index=5,                    # e.g., i_L
                color_increasing="rgba(0,255,0,0.3)",
                color_decreasing="rgba(255,0,0,0.3)",
                draw_extrema=True,
                prominence=0.01 * np.max(np.abs(yr[:,5]))
            ),
        ]
    )


def run_rc_sweep(net, omega_values, base_path="."):
    """
    Run simulation for multiple omega values.
    """

    results = {}

    for omega in omega_values:
        print(f"Running simulation for omega = {omega:.3e}")

        tr, yr = run_rc_transient(net, omega)

        plot_rc_results(tr, yr, omega, base_path=base_path)

        results[omega] = (tr, yr)

    return results

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():

    fstep    = 100000
    fstart   = 2 * math.pi * 1e3
    fstop    = 2 * math.pi * 100e6
    omega    = np.linspace(fstart, fstop, fstep)

    params = {"R": DEFAULT_R, "C": DEFAULT_C, "L": DEFAULT_L, "VDD": DEFAULT_VDD}

    nw_m0             = parallel_lrc_ntwrk(lti=True)
    nw_m0.set_parallel_lrc_params(**params)
    nw_m1             = parallel_lrc_ntwrk(lti=False, input_type=parallel_lrc_ntwrk._sine_source)
    nw_m1.set_parallel_lrc_params(**params)

    net = nw_m0

    eqs = net.ckt.build_mathjax_equations()
    math_block = quarto_math_block(eqs)
    print(math_block)

    net.ckt.get_ss()
    sys = ct.ss(
                net.ckt.A,
                net.ckt.B,
                net.ckt.C,
                net.ckt.D
            )

    # Force SISO: output 0, input 1
    #  - Output 0 = Vc
    #  - Input 1  = i_in
    sys_siso = sys[0, 1]


    # Taking the transfer function of sys_siso gives:
    # Z = Vc / i_in
    #sys_impedance = sys[0, 0]
    sys_Z = ct.ss2tf(sys_siso)

    # Poles and zeros
    poles = ct.pole(sys_siso)
    zeros = ct.zero(sys_siso)

    # Frequency response
    mag, phase, omega_out = ct.freqresp(sys_Z, omega=omega)

    label = "parallel_lrc_m0"

    # 1. Define the target directory
    output_dir = Path("../../doc/parallel_lrc")

    # 2. Create the directory (and any necessary parent directories) if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # 3. Create the full file path
    file_path = output_dir / f"bode_{label}.html"

    plot_bode_plotly(
        mag,
        phase,
        omega_out,
        single_out=True,
        output_idx=(0, 0),
        filename=str(file_path)
    )

    net = nw_m1

    # Example 1: Explicit list
    omega_values = [
        50.3e3 * (2*math.pi),
        60.3e3 * (2*math.pi)
    ]

    results = run_rc_sweep(net, omega_values, base_path=output_dir)

    exit(1)


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

    # Solve OPs and build small-signal circuits
    # for net, label in [(nw_m2, "ro_m2"), (nw_no_ro_m2, "no_ro_m2"), (nw_m1, "ro_m1"), (nw_no_ro_m1, "no_ro_m1")]:

    #     results = sweep(net, lrg_vs_sweep)
    #     plot_sweep(lrg_vs_sweep, results, suffix=label)

    # Solve OPs and build small-signal circuits
    # for net, base_label in [(nw_m2, "ro_m2_sml")]:
    #     for Rref in Rref_sweep:

    #         # Update parameters (including Rref)
    #         params_with_rref = dict(params)
    #         params_with_rref["Rref"] = Rref
    #         net.set_beta_multiplier_params(**params_with_rref)

    #         # Solve operating point
    #         root_start = np.zeros(net.ckt.num_edges)
    #         x = op_solve(net.ckt, root_start)
    #         op = net.ckt.scb

    #         # Build small-signal circuit with loop break
    #         net.add_sml_ckt(op=op, gate_topology=net.gate_break_m1_m2())

    #         # Extract operating-point small-signal parameters
    #         gm_m1, ro_m1 = net.beta_mult.nmos_m1.get_op_sml(op)
    #         gm_m2, ro_m2 = net.beta_mult.nmos_m2.get_op_sml(op)
    #         gm_m3, ro_m3 = net.beta_mult.pmos_m3.get_op_sml(op)
    #         gm_m4, ro_m4 = net.beta_mult.pmos_m4.get_op_sml(op)

    #         #print(f"\n[{cfg['label']}]")
    #         print(f"m1 gm = {gm_m1:.6f}, ro = {ro_m1:.6f}")
    #         print(f"m2 gm = {gm_m2:.6f}, ro = {ro_m2:.6f}")
    #         print(f"m3 gm = {gm_m3:.6f}, ro = {ro_m3:.6f}")
    #         print(f"m4 gm = {gm_m4:.6f}, ro = {ro_m4:.6f}")

    #         # Run small-signal sweep
    #         results = sml_sweep_m2(net, m2_sweep, Rref)

    #         # Generate a unique label per Rref
    #         label = f"{base_label}_Rref={Rref:g}"

    #         plot_sweep(m2_sweep, results, suffix=label)

    #Rref_sweep = np.linspace(6.5e3)

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
        x = op_solve(net.ckt, root_start)
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

            A = net.ckt_sml.A
            B = net.ckt_sml.B[:, [0]]   # keep column shape (n,1)
            C = net.ckt_sml.C[[0], :]   # keep row shape (1,n)
            D = net.ckt_sml.D[[0], [0]] # shape (1,1)

            # # Assume A, B, C, D are your matrices (n x n, n x 1, 1 x n, 1 x 1)
            # omega_num = np.logspace(0, 8, 500)  # example: 1 Hz to 100 MHz
            # H = np.zeros_like(omega, dtype=complex)

            # for k, w in enumerate(omega_num):
            #     # Solve (A + j*w*I) V = -B
            #     # Use pseudo-inverse to avoid singularity issues
            #     A_freq = A + 1j*w*np.eye(A.shape[0])
            #     Vaux = - np.linalg.pinv(A_freq) @ B
            #     # Output response
            #     H[k] = C @ Vaux + D

            # # Magnitude and phase
            # mag_num = np.abs(H)
            # phase_num = np.angle(H, deg=True)

            # # Optional: convert magnitude to dB
            # mag_dB_num = 20*np.log10(mag_num)

            # # Example: print DC and high-frequency gain
            # print("DC gain (linear):", mag_num[0])
            # print("DC gain (dB):", mag_dB_num[0])
            # print("HF gain (linear):", mag_num[-1])
            # print("HF gain (dB):", mag_dB_num[-1])

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
            plot_bode_plotly(
               mag,
               phase,
               omega_out,
               single_out=True,
               output_idx=(0, 0),
              filename=f"bode_{label}.html"
            )

            # plot_bode_plotly(
            #    mag_num,
            #    phase_num,
            #    omega_num,
            #    single_out=True,
            #    output_idx=(0, 0),
            #   filename=f"bode_{label}_num.html"
            # )

            plot_pz_plotly(
                sys_siso,
                poles=poles,
                zeros=zeros,
                filename=f"pz_{label}.html"
            )

            plot_bode_plotly(
                mag_scypi,
                phase_scypi,
                w,
                single_out=True,
                output_idx=(0, 0),
                filename=f"bode_{label}_scipy.html"
            )

    exit(1)

    # # Solve OPs and build small-signal circuits
    # for net, label in [(nw_m2, "ro_m2_sml"), (nw_no_ro_m2, "no_ro_m2_sml"), (nw_m1, "ro_m1_sml"), (nw_no_ro_m1, "no_ro_m1_sml")]:
    #     net.set_beta_multiplier_params(**params)
    #     root_start = np.zeros(net.ckt.num_edges)
    #     x = op_solve(net.ckt, root_start)
    #     op = net.ckt.scb
    #     net.add_sml_ckt(op=op)

    #     results = sml_sweep(net, vs_sweep)
    #     plot_sweep(vs_sweep, results, suffix=label)



    # Stability analysis
    for net, label in [(nw_m2, "ro_m2"), (nw_no_ro_m2, "no_ro_m2"), (nw_no_ro_m2_stray, "no_ro_m2_stray")]:
        net.set_beta_multiplier_params(**params)
        root_start = np.zeros(net.ckt.num_edges)
        x = op_solve(net.ckt, root_start)
        op = net.ckt.scb
        net.add_sml_ckt(op=op)

        iref = net.get_Rref_current(x)
        #print (iref)

        #eqs = net.ckt_sml.build_mathjax_equations()
        #math_block = quarto_math_block(eqs)
        #print(math_block)

        print(f"iref = {iref:.6f}")

        gm_m1, ro_m1 = net.beta_mult.nmos_m1.get_op_sml(op)
        gm_m2, ro_m2 = net.beta_mult.nmos_m2.get_op_sml(op)
        gm_m3, ro_m3 = net.beta_mult.pmos_m3.get_op_sml(op)
        gm_m4, ro_m4 = net.beta_mult.pmos_m4.get_op_sml(op)

        print(f"m1 gm = {gm_m1:.6f}, ro = {ro_m1:.6f}")
        print(f"m2 gm = {gm_m2:.6f}, ro = {ro_m2:.6f}")
        print(f"m3 gm = {gm_m3:.6f}, ro = {ro_m3:.6f}")
        print(f"m4 gm = {gm_m4:.6f}, ro = {ro_m4:.6f}")

    params = {"Rref": DEFAULT_RREF, "Vs": 0.05, "Cstray": 100e-12}
    nw_no_ro_m2_stray.set_beta_multiplier_sml_params(**params)
    x0    = np.zeros (nw_no_ro_m2_stray.ckt_sml.get_num_sys_vars())
    tr, yr = ode_solve(nw_no_ro_m2_stray.ckt_sml, tend=500e-9, tstep=10000, x0=x0)
    print (yr)
    plot_tr(tr, yr, filename="plot_beta_dyn.html")

    exit (1)
    for net, label in [(nw_m2, "ro_m2"), (nw_no_ro_m2, "no_ro_m2"), (nw_no_ro_m2_stray, "no_ro_m2_stray")]:
        net.set_beta_multiplier_params(**params)
        root_start = np.zeros(net.ckt.num_edges)
        #x = op_solve(net.ckt, root_start)
        #op = net.ckt.scb
        net.add_sml_ckt(op=root_start)

        #iref = net.get_Rref_current(x)
        #print (iref)

        eqs = net.ckt_sml.build_mathjax_equations()
        math_block = quarto_math_block(eqs)
        print(math_block)




if __name__ == "__main__":
    main()