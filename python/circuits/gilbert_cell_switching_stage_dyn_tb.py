import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../networks")
from gilbert_cell import *
from diff_amp_r_load_ntwrk import *

from scipy import optimize
from scipy.integrate import ode

import numpy as np

import itertools

import matplotlib.pyplot as plt

import math

import control as ct

def dypc_litmus0(t, sys, ckt):

    root_start = np.ones(ckt.num_edges)
    vs = 1

    root = optimize.root(circuit_eqn_dyn, root_start, args=(sys, ckt, t), tol=1e-9)
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

def circuit_eqn_dyn(x, sys, ckt, t):

    # current Im
    I = ckt.get_im(x=x, sys=sys, t=t)

    # voltage Vm
    V = ckt.get_vm(x=x, sys=sys, t=t)

    # get dependent currents
    I = ckt.get_im_dep()

    # Copy matrices
    qf_num = ckt.qf.copy()
    bf_num = ckt.bf.copy()

    eqn_qf = np.dot(qf_num, I)
    eqn_bf = np.dot(bf_num, V)

    return eqn_qf + eqn_bf

def trim(tstart, s, t):

    i = 0

    while t[i] < tstart:
        i += 1

    t_trim_start = t[i]

    t_trimmed = t[i:] - t_trim_start
    s_trimmed = s[i:]

    t_end_trimmed = t_trimmed[-1]

    return t_end_trimmed, s_trimmed, t_trimmed

def sweep_vlo(nw, vlo_sweep, vlo_bias, iref_val):

    vf_out = np.zeros(len(vlo_sweep))

    root_start = np.ones(nw.ckt.num_edges)

    nw.set_iref(iref_val)
    nw.set_vlo_bias(vlo_bias)

    for idx_vlo, vlo in enumerate(vlo_sweep):
        nw.set_vlo(vlo)

        # Operating Point
        root       = optimize.root(circuit_eqn_dyn, root_start, args=(0, nw.ckt,0 ), tol=1e-9)

        # Start point for next iteration
        root_start = root.x

        saturation_region = nw.check_saturation_region(root.x)

        if not root.success:
            print ("Failed to calculate small signal circuit")

        if not saturation_region:
            print ("Not in saturation region " + str(vlo))

        vf_out[idx_vlo] = nw.get_vf_output(root.x)

    return vf_out

def trimmed_fft(vf_out, tr, N, tstart):

    # Trim first 1us of transient to remove waveform until circuit in steady state
    t_end_trimmed, vf_out_trimmed, tr_trimmed = trim(tstart, vf_out, tr)
    timestep_trimmed = t_end_trimmed / N

    # New equally spaced time points for trimmed transient time
    tr = np.linspace(0, t_end_trimmed, N)

    # Resample the trimmed Vf_out at the equally spaced time points
    vf_out = np.interp(tr, tr_trimmed, vf_out_trimmed)

    # Calculate FFT
    sp   = np.fft.fft(vf_out)
    freq = np.fft.fftfreq(N, d=timestep_trimmed)

    # FFT magnitude in dB
    X =  20 * np.log10(np.abs(sp))

    return freq, X

def main():

    vlo_bias = 2.0
    vrf_bias = 1.5
    iref_val = 200e-6

    freq  = 1e6
    omega = math.pi * 2 * freq





    tend = 50e-6
    N    = 5500
    timestep = tend / N

    # Intermodulation
    # ---------------
    vlo_diff  = 0.3
    two_lrg_ac_params = {"omega_0": omega, "omega_1": 0.8 * omega, "mag": vlo_diff/2, "bias": vlo_bias}
    nw = diff_amp_r_load_ntwrk(two_sine_src=True, **two_lrg_ac_params)
    #nw.set_vlo_bias(bias=vlo_bias)
    nw.set_iref(iref_val)

    # Transient for large signal
    x0    = 4 * np.ones (nw.ckt.get_num_sys_vars())
    tr, yr = ode_solve(nw.ckt, tend=tend, tstep=10000, x0=x0)

    # Calculate Vf_out
    vf_out = yr[:,0] - yr[:,1]

    # Perform FFT
    freq, X = trimmed_fft(vf_out, tr, N, 1e-6)

    # Plot Vf_out
    fig, ax1 = plt.subplots()
    ax1.plot(tr, vf_out)
    plt.xlabel("t (s)")
    plt.ylabel("$V_f$")
    plt.title("$V_f$ Transient For Large Signal Circuit Intermodulation Test")
    plt.grid()
    plt.show(block=False)

    # Plot FFT
    fig, ax1 = plt.subplots()
    ax1.plot(freq[:N//2], X[:N//2])
    plt.xscale('log')
    plt.xlabel("f (Hz)")
    plt.ylabel("$V_f$ (db)")
    plt.title("FFT of $V_f$ For Large Signal Circuit Intermodulation Test")
    plt.grid()
    plt.show(block=False)

    #sys.exit(0)

    # Non-linear Gain
    # ----------------

    # Sweep Vlo and plot Vf
    #  - Demonstrates non-linear relationship of Vlo vs Vf

    # Higher degree of non-linearity with increasing Vlo_diff
    vlo_diff  = 0.6
    vlo_sweep = np.linspace(-vlo_diff/2, vlo_diff/2, 20)

    nw_nonlin     = diff_amp_r_load_ntwrk()
    vf_out_nonlin = sweep_vlo(nw_nonlin, vlo_sweep, vlo_bias, iref_val)

    # Plot Vlo vs Vf
    fig, ax1 = plt.subplots()
    ax1.plot(vlo_sweep, vf_out_nonlin)
    plt.xlabel("$V_{lo}$")
    plt.ylabel("$V_f$")
    plt.title("$V_{lo}$ vs $V_f$")
    plt.grid()
    plt.show(block=False)

    # Operating Point
    # ---------------

    # Apply the maximum differential input and calculate operating point.
    #  - Use the operating point to calculate the small signal circuit

    vlo_diff  = 0.6
    sml_ac_params = {"omega": omega, "mag": vlo_diff/2, "bias": 0}
    nw_nonlin.set_iref(iref_val)
    nw_nonlin.set_vlo_bias(bias=vlo_bias)
    nw_nonlin.set_vlo(vlo=vlo_diff/2)

    # Calculate small signal
    root_start = np.ones(nw_nonlin.ckt.num_edges)
    root       = optimize.root(circuit_eqn_dyn, root_start, args=(0, nw_nonlin.ckt, 0), tol=1e-9)
    if not root.success:
        print ("Failed to calculate small signal circuit")

    op = nw_nonlin.ckt.scb
    nw_nonlin.add_sml_ckt(op=op, sine_src=True, **sml_ac_params)

    # Small Signal Transient
    # -----------------------

    # Calculate small signal transient and perform FFT on result

    # Transient for small signal
    x0    = np.zeros (nw_nonlin.ckt_sml.get_num_sys_vars())
    tr_sml, yr_sml = ode_solve(nw_nonlin.ckt_sml, tend=tend, tstep=10000, x0=x0)

    # Calculate Vf_out
    vf_out_sml = yr_sml[:,0] - yr_sml[:,1]

    # Perform FFT
    freq_sml, X_sml = trimmed_fft(vf_out_sml, tr_sml, N, 1e-6)

    # Plot FFT
    fig, ax1 = plt.subplots()
    ax1.plot(freq_sml[:N//2], X_sml[:N//2])
    plt.xscale('log')
    plt.xlabel("f (Hz)")
    plt.ylabel("$V_f$ (db)")
    plt.title("FFT of $V_f$ For Small Signal Circuit")
    plt.grid()
    plt.show(block=False)

    # Operating Point
    # ---------------

    # Calculate large signal (non-linear) transient and perform FFT on result

    lrg_ac_params = {"omega": omega, "mag": vlo_diff/2, "bias": vlo_bias}
    nw = diff_amp_r_load_ntwrk(sine_src=True, **lrg_ac_params)
    nw.set_vlo_bias(bias=vlo_bias)
    nw.set_iref(iref_val)

    # Transient for large signal
    x0    = 4 * np.ones (nw.ckt.get_num_sys_vars())
    tr, yr = ode_solve(nw.ckt, tend=tend, tstep=10000, x0=x0)

    # Calculate Vf_out
    vf_out = yr[:,0] - yr[:,1]

    # Perform FFT
    freq, X = trimmed_fft(vf_out, tr, N, 1e-6)

    # Plot Vf_out
    fig, ax1 = plt.subplots()
    ax1.plot(tr, vf_out)
    plt.xlabel("t (s)")
    plt.ylabel("$V_f$")
    plt.title("$V_f$ Transient For Large Signal Circuit")
    plt.grid()
    plt.show(block=False)

    # Plot FFT
    fig, ax1 = plt.subplots()
    ax1.plot(freq[:N//2], X[:N//2])
    plt.xscale('log')
    plt.xlabel("f (Hz)")
    plt.ylabel("$V_f$ (db)")
    plt.title("FFT of $V_f$ For Large Signal Circuit")
    plt.grid()
    plt.show()

if __name__ == "__main__":

    main()
