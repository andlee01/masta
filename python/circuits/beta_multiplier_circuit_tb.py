import sys
import os
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from scipy import optimize
from scipy.integrate import ode

# Local imports
sys.path.extend(["../eqn", "../devices", "../networks"])
from eqn_syn import Circuit, Scoreboard
from base_devices import *
from beta_multiplier_ntwrk import *

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
M2_SWEEP_RANGE = (0, 1e-3, 100)

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

def plot_tr(tr, yr, filename="plot.html", title="Plot of yr vs tr", x_label="tr", y_label="yr"):
    """
    Plots yr vs tr using Plotly and saves the plot as an HTML file.

    Parameters:
    - tr: list or array-like, values for x-axis
    - yr: list or array-like (can be 2D), values for y-axis
    - filename: str, name of output HTML file
    - title: str, plot title
    - x_label: str, label for x-axis
    - y_label: str, label for y-axis
    """
    # Convert to numpy array and flatten yr
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
    
    fig.write_html(filename)
    print(f"Plot saved as {filename}")



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

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():
    vs_sweep = np.linspace(*SWEEP_RANGE)
    m2_sweep = np.linspace(*M2_SWEEP_RANGE)
    lrg_vs_sweep = np.linspace(*LRG_SWEEP_RANGE)
    Rref_sweep = np.linspace(*RREF_SWEEP_RANGE)
    params = {"Rref": DEFAULT_RREF, "Vs": DEFAULT_VS, "Cstray": 100e-12}

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

    #         # Run small-signal sweep
    #         results = sml_sweep_m2(net, m2_sweep, Rref)

    #         # Generate a unique label per Rref
    #         label = f"{base_label}_Rref={Rref:g}"

    #         plot_sweep(m2_sweep, results, suffix=label)

    #Rref_sweep = np.linspace(6.5e3)

    # Solve OPs and build SS small-signal circuits
    for net, base_label in [(nw_m2_ss, "ro_m2_sml")]:
        #for Rref in Rref_sweep:
            
            Rref = 6.5e3

            # Update parameters (including Rref)
            params_with_rref = dict(params)
            params_with_rref["Rref"] = Rref
            net.set_beta_multiplier_params(**params_with_rref)

            # Solve operating point
            root_start = np.zeros(net.ckt.num_edges)
            x = op_solve(net.ckt, root_start)
            op = net.ckt.scb

            # Build small-signal circuit with loop break
            net.add_sml_ckt(op=op, \
                            gate_topology=net.gate_break_m1_m2(), \
                            degen_topology=net._capacitive_degen_topology(),
                            output_topology=net._output_vds_m2_topology())

            # State Space
            net.ckt_sml.get_ss()

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