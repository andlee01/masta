"""Sweep equal NMOS drain currents using the existing smooth MOS model.

Saturation is classified by the conventional Vds >= Vgs - Vth boundary.
The existing model smoothly transitions between regions; this boundary is
a diagnostic, rather than a switch in its current law.
"""

from pathlib import Path
import argparse
import sys

import numpy as np
import plotly.graph_objects as go
from scipy import optimize

# Local imports
ROOT = Path(__file__).resolve().parents[1]
for directory in ("devices", "eqn", "networks", "subckt"):
    sys.path.insert(0, str(ROOT / directory))

from beta_multiplier_nmos_ntwrk import beta_multiplier_nmos_ntwrk
from beta_multiplier_circuit_tb import (
    OperatingPointError, circuit_eqn, circuit_jacobian, op_solve,
    _equation_scales, _variable_scales,
)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

DEFAULT_RREF = 6.5e3
CURRENT_SWEEP = (0.1e-6, 30e-6, 600)
DEFAULT_DRAIN_LIMIT = 5.0
DEFAULT_OUTPUT = ROOT.parent / "doc" / "beta_multiplier" / "beta_multiplier_nmos.html"

# ------------------------------------------------------------------------------
# Large-signal analysis
# ------------------------------------------------------------------------------

def sweep(network, currents, drain_limit):
    """Continue DC solutions without clipping drain voltages to compliance."""
    results = []
    guess = np.zeros(network.ckt.num_edges)
    for current in currents:
        network.set_current(current)
        try:
            if network.m2.l_lambda == 0:
                guess = network.current_driven_guess(current)
            solution = op_solve(network.ckt, guess)
            # The shared OP solver allows 1 nA KCL error. Refine in scaled
            # coordinates so the equal-current sweep is accurate below 1 µA.
            ckt = network.ckt
            eq_scale = _equation_scales(ckt)
            x_scale = _variable_scales(ckt)
            refined = optimize.least_squares(
                lambda z: circuit_eqn(z * x_scale, 0, ckt, 0) / eq_scale,
                solution / x_scale,
                jac=lambda z: (circuit_jacobian(z * x_scale, 0, ckt, 0) *
                               x_scale[None, :] / eq_scale[:, None]),
                ftol=1e-13, xtol=1e-13, gtol=1e-13, max_nfev=200,
            )
            solution = refined.x * x_scale
            residual = circuit_eqn(solution, 0, network.ckt, 0)
            point = network.get_operating_point()
            if (np.max(np.abs(residual)) > 1e-9 or
                    not np.allclose([point["i1"], point["i2"]], current,
                                    rtol=1e-5, atol=1e-10)):
                raise OperatingPointError("current/residual validation failed")
            guess = solution.copy()
            boundary = point["vg"] - network.m2.vth
            saturated = (point["vgs2"] > network.m2.vth and
                         point["vd2"] >= boundary)
            point.update(boundary=boundary, saturated=saturated,
                         acceptable=saturated and point["vd2"] <= drain_limit,
                         solved=True)
        except OperatingPointError as error:
            print(f"No validated solution at {current * 1e6:.6g} uA: {error}")
            point = dict.fromkeys(("vg", "vs2", "vd2", "vgs2", "vds2",
                                   "i1", "i2", "boundary", "overdrive2"), np.nan)
            point.update(saturated=False, acceptable=False, solved=False)
        results.append(point)
    return {key: np.array([point[key] for point in results])
            for key in results[0]}

# ------------------------------------------------------------------------------
# Plotting functions
# ------------------------------------------------------------------------------

def main_current_interval(currents, results, drain_limit):
    """Interpolate the boundaries of the widest sampled acceptable interval.

    The smooth model can also produce a separate near-threshold interval;
    this annotation describes the main conducting interval only.
    """
    indices = np.flatnonzero(results["acceptable"])
    if not indices.size:
        return None
    groups = np.split(indices, np.where(np.diff(indices) != 1)[0] + 1)
    group = max(groups, key=lambda g: currents[g[-1]] - currents[g[0]])
    margin = np.minimum.reduce([
        results["vd2"] - results["boundary"],
        drain_limit - results["vd2"],
        results["overdrive2"],
    ])

    def crossing(inside, outside):
        if (outside < 0 or outside >= len(currents) or
                not np.isfinite(margin[outside]) or margin[outside] >= 0):
            return currents[inside]
        fraction = margin[inside] / (margin[inside] - margin[outside])
        return currents[inside] + fraction * (currents[outside] - currents[inside])

    return crossing(group[0], group[0] - 1), crossing(group[-1], group[-1] + 1)


def plot_sweep(currents, results, drain_limit, filename):
    x = currents * 1e6
    boundary = results["boundary"]
    # Shade only where the geometrical voltage window exists and Vov > 0.
    window = (boundary <= drain_limit) & (results["overdrive2"] > 0)
    lower = np.where(window, boundary, np.nan)
    upper = np.where(window, drain_limit, np.nan)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=lower, mode="lines", line=dict(width=0),
                             showlegend=False, hoverinfo="skip"))
    fig.add_trace(go.Scatter(
        x=x, y=upper, mode="lines", line=dict(width=0), fill="tonexty",
        fillcolor="rgba(40,170,90,0.16)", name="Acceptable drain-voltage window",
        hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=x, y=boundary, name="Saturation boundary: V<sub>G</sub> − V<sub>TH</sub>",
                             line=dict(dash="dash")))
    fig.add_trace(go.Scatter(x=x, y=np.full_like(x, drain_limit),
                             name="Available drain-voltage limit", line=dict(dash="dot")))
    fig.add_trace(go.Scatter(x=x, y=results["vd2"], name="M<sub>2</sub> drain voltage, λ = 0.0125 V⁻¹",
                             mode="lines+markers", marker=dict(color="royalblue", size=6),
                             line=dict(color="royalblue"), connectgaps=False))
    fig.update_layout(
        title=dict(text=r"$\text{NMOS beta multiplier: } I_{D1}=I_{D2}$",
                   x=0.5, xanchor="center", font=dict(size=28)),
        xaxis_title=r"$\text{Common drain current, } I_D\;(\mu\mathrm{A})$",
        yaxis_title=r"$M_2\text{ drain voltage, } V_{D2}\;(\mathrm{V})$",
        template="plotly_white", hovermode="x unified",
        font=dict(size=22), hoverlabel=dict(font_size=22),
        height=900, autosize=True, showlegend=False,
        margin=dict(l=115, r=45, t=95, b=100),
    )
    fig.update_xaxes(title_font=dict(size=25), tickfont=dict(size=22),
                     title_standoff=24, automargin=True)
    fig.update_yaxes(title_font=dict(size=25), tickfont=dict(size=22),
                     title_standoff=24, automargin=True)
    interval = main_current_interval(currents, results, drain_limit)
    if interval is not None:
        imin, imax = np.asarray(interval) * 1e6
        for value in (imin, imax):
            fig.add_vline(x=value, line=dict(color="darkgreen", width=2, dash="dash"))
        fig.add_annotation(
            x=(imin + imax) / 2, y=0.65 * drain_limit,
            text=(rf"$\text{{Range for }}\lambda=0.0125:\;"
                  rf"{imin:.2f}\lesssim I_D\lesssim {imax:.2f}\,\mu\mathrm{{A}}$"),
            showarrow=True, ax=-210, ay=-65, arrowhead=2,
            font=dict(size=22), bgcolor="rgba(255,255,255,0.95)",
            bordercolor="darkgreen", borderpad=8,
        )
        print(f"Interpolated main acceptable current range: {imin:.4f}–{imax:.4f} µA")
    fig.add_annotation(
        x=0.98, xref="paper", y=drain_limit, yref="y",
        text=rf"$V_{{D2,\max}}={drain_limit:g}\,\mathrm{{V}}\;\text{{(drain limit)}}$", showarrow=False,
        xanchor="right", yanchor="bottom", yshift=8,
        font=dict(size=22), bgcolor="rgba(255,255,255,0.9)",
    )
    finite = np.flatnonzero(np.isfinite(boundary))
    if finite.size:
        index = finite[-1]
        fig.add_annotation(
            x=x[index], y=boundary[index],
            text=r"$\text{Saturation boundary: } V_{D2}=V_G-V_{TH}$", showarrow=True,
            ax=-180, ay=-80, arrowhead=2,
            font=dict(size=22), bgcolor="rgba(255,255,255,0.9)",
        )
    # Retain out-of-compliance points in the HTML; zoom out to inspect them.
    fig.update_yaxes(range=[0, drain_limit * 1.2])
    filename.parent.mkdir(parents=True, exist_ok=True)
    # Give the standalone graph a viewport-sized container. Plotly's default
    # percentage height otherwise collapses inside the document's iframe.
    html = fig.to_html(include_mathjax="cdn", full_html=True,
                       default_width="100%", default_height="100vh",
                       config={"responsive": True})
    html = html.replace("<head>", """<head>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
html, body { margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; }
.plotly-graph-div { width: 100% !important; height: 100vh !important; }
</style>""", 1)
    filename.write_text(html, encoding="utf-8")
    print(f"Plot saved to {filename}")

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--drain-limit", type=float, default=DEFAULT_DRAIN_LIMIT)
    args = parser.parse_args()
    if args.drain_limit <= 0:
        parser.error("--drain-limit must be positive")
    network = beta_multiplier_nmos_ntwrk(Rref=DEFAULT_RREF)
    currents = np.linspace(*CURRENT_SWEEP)
    results = sweep(network, currents, args.drain_limit)
    print(f"Validated {results['solved'].sum()}/{len(currents)} operating points")
    indices = np.flatnonzero(results["acceptable"])
    if indices.size:
        for group in np.split(indices, np.where(np.diff(indices) != 1)[0] + 1):
            print(f"Acceptable sampled interval: {currents[group[0]]*1e6:.6g} to "
                  f"{currents[group[-1]]*1e6:.6g} µA")
        print(f"Grid spacing: {np.diff(currents)[0]*1e6:.6g} µA. Near threshold, "
              "the conventional region label reflects the model's smoothing.")
    else:
        print("No acceptable operating points in this sweep")
    plot_sweep(currents, results, args.drain_limit, args.output)


if __name__ == "__main__":
    main()
