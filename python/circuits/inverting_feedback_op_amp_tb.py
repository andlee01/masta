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
from inverting_feedback_op_amp import *

import control as ct
from scipy.signal import freqresp, StateSpace

from pathlib import Path

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
DEFAULT_RF = 1000
DEFAULT_RIN = 100
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
    V = ckt.get_vm_dep()
    qf_num, bf_num = ckt.qf.copy(), ckt.bf.copy()
    return np.dot(qf_num, I) + np.dot(bf_num, V)


def op_solve(ckt, root_start):
    """Find DC operating point."""
    root = optimize.root(circuit_eqn, root_start, args=(0, ckt, 0), tol=1e-9, method='lm')
    if not root.success:
        print("⚠️ OP solver did not converge, retrying with updated initial guess")
        return op_solve(ckt, root.x)
    return root.x

def op_solve_trial(ckt, root_start):
    # Define physical bounds: no node can exceed Vsat or a reasonable margin
    # Assuming vsat is 15, we might bound between -20 and 20
    v_limit = 100e-6
    lower_bounds = np.full_like(root_start, -v_limit)
    upper_bounds = np.full_like(root_start, v_limit)

    # Note: 'lm' and 'dogleg' don't support bounds, so we use 'hybr' or 'trf'
    # 'least_squares' is often more robust for bounded circuit problems
    res = optimize.least_squares(
        circuit_eqn, 
        root_start, 
        args=(0, ckt, 0), 
        bounds=(lower_bounds, upper_bounds),
        xtol=1e-9
    )
    
    if not res.success:
        print(f"⚠️ Solver failed: {res.message}")
        return root_start # Or handle retry logic
        
    return res.x

# def dypc_litmus0(t, sys, ckt):

#     root_start = np.ones(ckt.num_edges)
#     vs = 1

#     root = optimize.root(circuit_eqn, root_start, args=(sys, ckt, t), tol=1e-9, method='lm')
#     root_start = root.x
#     if not root.success:
#         print (root.x)
#         sys.exit(0)

#     return ckt.get_dy(x=root.x)

# def ode_solve(ckt, tend=50e-12, tstep=10000, x0=0):

#     num_sys_vars    = ckt.get_num_sys_vars()

#     t     = np.linspace(0, tend, tstep)

#     r = ode(dypc_litmus0).set_integrator('lsoda', method='bdf', atol=1e-9, rtol=1e-9)
#     r.set_initial_value(x0, 0.0)
#     r.set_f_params(ckt)

#     y    = np.empty((tstep, num_sys_vars))
#     y[0] = x0

#     k = 1
#     while r.successful() and k < tstep:
#         r.integrate(t[k])

#         y[k] = r.y

#         k += 1

#     return t, y

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


def plot_tr(
    tr, yr,
    filename="plot.html",
    title="System Response",
    labels=None,
    quantity_types=None,
    v_bound=None,   # Fixed bounds for Voltage axis
    i_bound=None    # Fixed bounds for Current axis
):
    import numpy as np
    import plotly.graph_objects as go

    # Ensure inputs are numpy arrays
    tr, yr = np.asarray(tr), np.asarray(yr)
    if yr.ndim == 1:
        yr = yr.reshape(-1, 1)

    num_signals = yr.shape[1]
    
    print (labels)

    # Fallback for labels and types if not provided
    labels = labels or [f"Sig {i}" for i in range(num_signals)]
    quantity_types = quantity_types or ["voltage"] * num_signals

    # Identify indices for scaling and axis assignment
    v_idxs = [i for i, qt in enumerate(quantity_types) if qt.lower() == "voltage"]
    i_idxs = [i for i, qt in enumerate(quantity_types) if qt.lower() == "current"]

    # Compute bounds for axes if not explicitly provided
    if v_bound is None:
        v_bound = np.max(np.abs(yr[:, v_idxs])) * 1.1 if v_idxs else 1

    if i_bound is None:
        i_bound = np.max(np.abs(yr[:, i_idxs])) * 1.1 if i_idxs else 1

    # Create figure
    fig = go.Figure()

    for i in range(num_signals):
        is_v = quantity_types[i].lower() == "voltage"
        fig.add_trace(go.Scatter(
            x=tr,
            y=yr[:, i],
            name=labels[i],
            yaxis="y" if is_v else "y2",
            mode='lines'
        ))

    # Layout configuration
    fig.update_layout(
        title=title,
        xaxis_title="Time (s)",
        yaxis=dict(
            title="Voltage (V)",
            side="left",
            range=[-v_bound, v_bound],
            zeroline=True,
            zerolinecolor="lightgray"
        ),
        yaxis2=dict(
            title="Current (A)",
            side="right",
            overlaying="y",
            showgrid=False,
            range=[-i_bound, i_bound],
            zeroline=True,
            zerolinecolor="lightgray"
        ),
        template="plotly_white",
        legend=dict(
            orientation="v",    # Stacked vertically
            yanchor="top",
            y=1,                # Top of the plot
            xanchor="left",
            x=1.05,             # Just to the right of the y2-axis
            font=dict(size=12), # Fixed font size helps prevent clipping
            itemsizing='constant' # Keeps the legend icons a consistent size
        ),
        # Increase the right margin significantly for the legend labels
        margin=dict(l=80, r=180, t=80, b=80)
    )
    # Save output
    fig.write_html(filename, include_mathjax="cdn")
    print(f"Plot saved to {filename}")

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

import numpy as np
from scipy.fft import fft, fftfreq

def compute_fft(tr, yr, crop_cycles=0, omega_input=None):
    """
    Compute the FFT of a time-domain signal.
    
    Parameters:
    -----------
    tr : ndarray (Time vector)
    yr : ndarray (Signal values)
    crop_cycles : float (Number of cycles to ignore from the start)
    omega_input : float (Optional: used to calculate how many points to crop)
    """
    # 1. Handle cropping for steady-state analysis
    if crop_cycles > 0 and omega_input is not None:
        period = 2 * np.pi / omega_input
        start_time = crop_cycles * period
        mask = tr >= start_time
        tr, yr = tr[mask], yr[mask]

    n = len(tr)
    # Calculate average timestep to find sampling frequency
    dt = (tr[-1] - tr[0]) / (n - 1)
    fs = 1.0 / dt

    # 2. Perform FFT
    # We use rfft (Real FFT) because circuit signals are real-valued
    yf = np.fft.rfft(yr)
    xf = np.fft.rfftfreq(n, d=dt) # Frequency in Hz

    # Convert Hz to rad/s to match your Bode function
    omega = 2 * np.pi * xf

    # 3. Get Magnitude and Phase
    # Multiply magnitude by 2/n to normalize (except DC component)
    mag = 2.0/n * np.abs(yf)
    phase = np.angle(yf)

    return mag, phase, omega

import plotly.graph_objects as go
from plotly.subplots import make_subplots

def plot_fft_plotly(mag, phase, omega, filename="fft_output.html", title="FFT Analysis", ref_val=1.0, log_x=True):
    """
    Plot FFT Magnitude and Phase with optional log scale on the x-axis.
    """
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import os

    # Calculate relative magnitude in dB
    mag_db = 20 * np.log10((mag / ref_val) + 1e-12)
    phase_deg = np.degrees(phase)

    # Filter out 0 for log scale to avoid Plotly rendering issues
    # (Log of 0 is undefined)
    if log_x:
        mask = omega > 0
        omega = omega[mask]
        mag_db = mag_db[mask]
        phase_deg = phase_deg[mask]

    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.1,
        subplot_titles=(f"Magnitude (dB relative to {ref_val})", "Phase (deg)")
    )

    # Magnitude Plot
    fig.add_trace(
        go.Scatter(x=omega, y=mag_db, mode="lines", name="Magnitude"),
        row=1, col=1
    )

    # Phase Plot
    fig.add_trace(
        go.Scatter(x=omega, y=phase_deg, mode="lines", name="Phase"),
        row=2, col=1
    )

    # --- Update Axes to Log Scale ---
    x_axis_config = dict(
        title_text="Frequency (rad/s)",
        type="log" if log_x else "linear",
        gridcolor="lightgray"
    )

    fig.update_xaxes(x_axis_config, row=1, col=1)
    fig.update_xaxes(x_axis_config, row=2, col=1)
    
    fig.update_yaxes(title_text="Gain (dB)", row=1, col=1)
    fig.update_yaxes(title_text="Phase (deg)", row=2, col=1)

    fig.update_layout(
        height=700,
        width=900,
        template="plotly_white",
        title_text=title,
        showlegend=False
    )

    # Save and Print
    fig.write_html(filename, include_plotlyjs="cdn")
    print(f"✅ FFT Plot successfully saved to: {os.path.abspath(filename)}")

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
DEFAULT_MAG = 350e-3

def run_rc_transient_manual(net, omega, num_cycles=NUM_CYCLES, mag=DEFAULT_MAG, tstep_per_cycle=100):
    """
    Run simulation for a single omega by stepping through time manually
    and solving the operating point (OP) at each step.
    """
    # Derived time parameters
    f = omega / (2 * math.pi)
    period = 1 / f
    tend = num_cycles * period
    num_steps = int(num_cycles * tstep_per_cycle)
    
    # Generate time vector
    tr = np.linspace(0, tend, num_steps)
    
    # Initialize container for results and the initial guess for the solver
    yr_list = np.zeros((num_steps, net.ckt.num_edges))
    i_input = np.zeros(num_steps)
    v_out   = np.zeros(num_steps)
    x_guess = 1e-3 * np.ones(net.ckt.num_edges)
    
    # Manual time-stepping loop
    for i, t in enumerate(tr):
        # 1. Calculate the instantaneous value of the source
        inst_val = mag * math.cos(omega * t)
        
        # 2. Update the circuit's source to this specific DC value
        # NOTE: You will need to replace this method call with the actual method 
        # your 'net' or 'ckt' object uses to set an instantaneous/DC source value.
        net.set_source(source="dc_input_src", val=inst_val) 
        
        # 3. Solve the operating point
        x_sol = op_solve(net.ckt, x_guess)
        
        I = net.ckt.get_im(x=x_sol, sys=sys, t=t)
        V = net.ckt.get_vm(x=x_sol, sys=sys, t=t)
        I = net.ckt.get_im_dep()
        V = net.ckt.get_vm_dep()

        net.ckt.scb.x = x_sol
        net.ckt.scb.v = V
        net.ckt.scb.i = I

        #print (x_sol)

        # print (inst_val)
        # print (x_sol)

        # 4. Store the result and use it as the starting guess for the next step 
        # (This "warm start" helps the LM solver converge much faster)
        yr_list[i] = x_sol
        i_input[i] = inst_val
        v_out[i]   = net.get_vout(x_sol)
        x_guess = x_sol
        
    # Convert list of results to a 2D numpy array
    yr = np.array(yr_list)
    
    # Ensure yr is 2D with correct orientation
    if yr.shape[0] != len(tr):
        yr = yr.T

    # Append signals
    yr = np.column_stack((yr, i_input))
    yr = np.column_stack((yr, v_out))

    return tr, yr

def plot_rc_results_1(net, tr, yr, omega, prefix="inverting_feedback_op_amp", base_path="."):
    """
    Generate standard plots using dynamic labels from the circuit instances.
    """
    output_dir = Path(base_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    label = f"{prefix}_{omega}"
    transient_filename = output_dir / f"transient_{label}.html"

    # 1. Integrate over net.get_instance(idx) for the first N labels
    # These correspond to the system variables (v_c, i_l, etc.)
    num_edges = net.ckt.num_edges
    dynamic_labels = []
    quantity_types = []
    
    for idx in range(num_edges):
        # Retrieve the instance name (e.g., 'C1', 'L1', 'R1')
        try:
            elem = net.ckt.get_edge_info(idx)   
            instance_name = elem.get_instance()

            if elem.get_type() == ElementType.resistor:
                dynamic_labels.append(f"$i_{{{instance_name}}}$")
                quantity_types.append("current")
            if elem.get_type() == ElementType.voltage_src:
                dynamic_labels.append(f"$i_{{{instance_name}}}$")
                quantity_types.append("current")

        except:
            instance_name = "dep"
            dynamic_labels.append(instance_name)
            quantity_types.append("current")
        
        # We need to decide if the state variable is V or I.
        # Usually, if it's a capacitor in the tree, it's V; if it's an inductor link, it's I.
        # You can check the element type to format the LaTeX correctly:
        #
        # if net.ckt.is_co_tree_elem(idx):
        #     dynamic_labels.append(f"$V_{{{instance_name}}}$")
        #     quantity_types.append("voltage")
        # else:
        #     dynamic_labels.append(f"$i_{{{instance_name}}}$")
        #     quantity_types.append("current")

        # if elem.get_type() == ElementType.resistor:
        #     dynamic_labels.append(f"$i_{{{instance_name}}}$")
        #     quantity_types.append("current")
        # if elem.get_type() == ElementType.voltage_src:
        #     dynamic_labels.append(f"$i_{{{instance_name}}}$")
        #     quantity_types.append("current")

    # 2. Handle the "appended" signals
    # Your run_rc_transient function appends: i_input, i_res, i_cap, i_diff
    # You may want to keep these hardcoded or map them as well:
    extra_labels = [r"$v_{IN}$", r"$v_{OUT}$"]
    all_labels = dynamic_labels + extra_labels
    all_quantities = quantity_types + ["voltage", "voltage"]

    # 3. Define quantity types for the plotter
    # This assumes yr[:, 0] is voltage and the rest are currents per your previous snippet
    # Adjust the slice if your num_edges > 2
    #quantity_types = ["voltage"] * 1 + ["current"] * (yr.shape[1] - 1)

    # 4. Bounds calculation
    v_bound = np.max(np.abs(yr[:, 0])) * 1.1 # Assuming index 0 is the voltage
    i_bound = np.max(np.abs(yr[:, 1:])) * 1.1 # All others are currents

    # 5. Standard transient plot
    plot_tr(
        tr=tr,
        yr=yr,
        v_bound=v_bound,
        i_bound=i_bound,
        labels=all_labels,
        quantity_types=all_quantities,
        filename=str(transient_filename)
    )


def run_rc_sweep(net, omega_values, base_path="."):
    """
    Run simulation for multiple omega values.
    """

    results = {}

    for omega in omega_values:
        print(f"Running simulation for omega = {omega:.3e}")

        tr, yr = run_rc_transient_manual(net, omega, mag=DEFAULT_MAG, tstep_per_cycle=100)

        plot_rc_results_1(net, tr, yr, omega, base_path=base_path)

        # 1. Use -1 to grab the last column of yr
        # 2. Use os.path.join to combine base_path with the filename
        fft_filename = os.path.join(base_path, f"fft_results_{omega}.html")

        mag, phase, omega_vec = compute_fft(tr, yr[:, -1], crop_cycles=2, omega_input=omega)

        plot_fft_plotly(mag, phase, omega_vec, filename=fft_filename, ref_val=DEFAULT_MAG)

        results[omega] = (tr, yr)

    return results

def plot_transfer_curve(net_element, v_limit, base_path, num_points=500):
    """
    Evaluates and plots the DC transfer characteristic of the non-linear amplifier.
    
    Parameters
    ----------
    net_element : object
        The instance of your amplifier class (containing nonlinear_gain).
    v_limit : float
        The range for the input voltage sweep (-v_limit to +v_limit).
    base_path : Path or str
        Directory to save the HTML output.
    """
    # 1. Create the input voltage sweep
    v_in_vec = np.linspace(-v_limit, v_limit, num_points)
    v_out_vec = []

    # 2. Mock a simple SCB (Simulation Control Block)
    # We only need it to have an .i array that nonlinear_gain can read
    class MockSCB:
        def __init__(self, current_val, idx):
            # We assume current * RIN_value = v_in
            # So we set i[idx] such that i * R = v_in
            self.i = {idx: current_val}

    # 3. Evaluate the non-linear function
    r_val = net_element.amp.RIN.get_value()
    for v_in in v_in_vec:
        # Calculate the current that would produce this v_in across RIN
        mock_i = v_in / r_val
        scb = MockSCB(mock_i, net_element.amp.RIN_idx)
        
        v_out = net_element.amp.nonlinear_gain(scb)
        v_out_vec.append(v_out)

    # 4. Create the Plotly Figure
    fig = go.Figure()

    # Add the non-linear curve
    fig.add_trace(go.Scatter(
        x=v_in_vec, y=v_out_vec,
        mode='lines',
        name='Non-linear Transfer',
        line=dict(color='firebrick', width=3)
    ))

    # # Add an ideal linear reference for comparison
    fig.add_trace(go.Scatter(
        x=v_in_vec, y=net_element.amp.gain * v_in_vec,
        mode='lines',
        name='Ideal Linear',
        line=dict(color='royalblue', dash='dash')
    ))

    fig.update_layout(
        title=f"Amplifier Transfer Characteristic (Gain={net_element.amp.gain:.1e})",
        xaxis_title="Input Voltage (V)",
        yaxis_title="Output Voltage (V)",
        template="plotly_white",
        hovermode="x unified"
    )

    # 5. Save the output
    os.makedirs(base_path, exist_ok=True)
    filename = os.path.join(base_path, "amplifier_transfer_curve.html")
    fig.write_html(filename, include_plotlyjs="cdn")
    
    print(f"✅ Transfer curve saved to: {os.path.abspath(filename)}")

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():

    fstep    = 100000
    fstart   = 2 * math.pi * 1e3
    fstop    = 2 * math.pi * 100e6
    omega    = np.linspace(fstart, fstop, fstep)

    params_fb    = {"RF": DEFAULT_RF, "RIN": DEFAULT_RIN, "linear_gain": 100e6, "a2": 0, "alpha3": 1e10, "vclip": 5.0}
    params_no_fb = {"linear_gain": 10, "a2": 0, "alpha3": 1e1, "vclip": 5.0}


    nw_m1             = inverting_feedback_op_amp(lti=False, input_type=inverting_feedback_op_amp._default_input_source, gain_type="linear")
    nw_m1.set_params(**params_fb)
    nw_m2             = inverting_feedback_op_amp(lti=False, input_type=inverting_feedback_op_amp._default_input_source, gain_type="nonlinear")
    nw_m2.set_params(**params_fb)
    nw_m3             = inverting_feedback_op_amp(lti=False, 
                                                  input_type=inverting_feedback_op_amp._default_input_source, 
                                                  gain_type="nonlinear",
                                                  feedback_type="None")
    nw_m3.set_params(**params_no_fb)

    eqs = nw_m1.ckt.build_mathjax_equations()
    math_block = quarto_math_block(eqs)
    print(math_block)

    # 1. Define your networks and their corresponding directory names
    networks = [nw_m1, nw_m2, nw_m3]
    dir_names = ["inverting_feedback_op_amp_linear", 
                 "inverting_feedback_op_amp_nonlinear",
                 "inverting_open_loop_op_amp_nonlinear"]

    v_limits = [120e-9, 120e-9, 100e-2]

    # 2. Define frequencies (defined once if they are the same for both runs)
    omega_values = [
        50.3e3 * (2 * math.pi),
        60.3e3 * (2 * math.pi)
    ]

    # 3. Loop through them simultaneously using zip()
    for net, dir_name, v_limit in zip(networks, dir_names, v_limits):
        print(f"Processing: {dir_name}...")
        
        # Define the target directory dynamically using f-strings
        output_dir = Path(f"../../doc/{dir_name}")

        # Create the directory (and parent directories) if it doesn't exist
        output_dir.mkdir(parents=True, exist_ok=True)

        # Run the sweep and save to the specific output directory
        results = run_rc_sweep(net, omega_values, base_path=output_dir)
        
        plot_transfer_curve(net, v_limit=v_limit, base_path=output_dir)

    print("All sweeps completed!")

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