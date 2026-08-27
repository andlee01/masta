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

        # 4. Store the result and use it as the starting guess for the next step 
        # (This "warm start" helps the LM solver converge much faster)
        yr_list[i] = x_sol
        i_input[i] = inst_val
        v_out[i]   = net.get_vout()
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

def plot_rc_results(net, tr, yr, omega, prefix="inverting_feedback_op_amp", base_path="."):
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
        
    # 2. Handle the "appended" signals
    # Your run_rc_transient function appends: i_input, i_res, i_cap, i_diff
    # You may want to keep these hardcoded or map them as well:
    extra_labels = [r"$v_{IN}$", r"$v_{OUT}$"]
    all_labels = dynamic_labels + extra_labels
    all_quantities = quantity_types + ["voltage", "voltage"]

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

        plot_rc_results(net, tr, yr, omega, base_path=base_path)

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

if __name__ == "__main__":
    main()