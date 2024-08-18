import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../networks")
from diff_amp_pmos_load_cs_bias_input_ntwrk import *

from scipy import optimize
from scipy.integrate import ode

import numpy as np

import itertools

import matplotlib.pyplot as plt

import math

import control as ct

def dypc_litmus0(t, sys, nw):

    #if t == 0:
    #    root_start = np.ones(nw.ckt.num_edges)
    #else:
    root_start = nw.ckt.scb.x

    root = optimize.root(circuit_eqn, root_start, args=(sys, nw.ckt, t), tol=1e-9)
    #root_start = root.x
    #if not root.success:
     #   print (t)
        #nw.ckt.scb.x = np.ones(nw.ckt.num_edges)
        #sys.exit(0)

    return nw.ckt.get_dy(x=root.x)

def ode_solve(nw, tend=50e-12, tstep=10000, x0=0):

    num_sys_vars    = nw.ckt.get_num_sys_vars()

    t     = np.linspace(0, tend, tstep, dtype=np.float64)
    #t     = np.linspace(0, tend, tstep)

    r = ode(dypc_litmus0).set_integrator('lsoda', method='bdf', atol=1e-6, rtol=1e-6)
    r.set_initial_value(x0, 0.0)
    r.set_f_params(nw)

    y    = np.empty((tstep, num_sys_vars))
    y[0] = x0

    vx_out = np.zeros(len(t))

    k = 1
    while r.successful() and k < tstep:
        r.integrate(t[k])

        y[k] = r.y

        vx_out[k] = nw.get_vout(nw.ckt.scb.x, y[k], t[k])

        k += 1

    return t, y, vx_out

def circuit_eqn(x, sys, ckt, t):

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

def trimmed_fft(vf_out, tr, tstart, A=1):

    # Trim first 1us of transient to remove waveform until circuit in steady state
    t_end_trimmed, vf_out_trimmed, tr_trimmed = trim(tstart, vf_out, tr)
    N = len(tr_trimmed)
    timestep_trimmed = t_end_trimmed / N

    vf_out = vf_out_trimmed

    # Window
    window = np.hanning(len(vf_out))
    #vf_out_windowed = window * vf_out
    vf_out_windowed = vf_out

    # Calculate FFT
    sp   = np.fft.fft(vf_out_windowed)
    freq = np.fft.fftfreq(N, d=timestep_trimmed)

    mag = np.abs(sp) / (N/2)

    # FFT magnitude in dB
    X =  20 * np.log10(mag / A)
    return freq, X, N

def main():

    np.random.seed(42)

    vin_omega = 1e6
    vin_mag   = 0.1

    ac_params = {"omega_vp": vin_omega, "mag_vp": 0.0, "bias_vp": 0.0,
                 "omega_vn": vin_omega, "mag_vn": vin_mag,     "bias_vn": 0.0}

    nw = diff_amp_pmos_load_cs_bias_input_ntwrk(**ac_params)

    nw.ckt.scb.x = np.ones(nw.ckt.num_edges)
    #x0    = 4.0 * np.zeros (nw.ckt.get_num_sys_vars())
    x0    = np.zeros (nw.ckt.get_num_sys_vars(), dtype=np.float64)

    tr, yr, v_vx_out = ode_solve(nw, tend=300e-6, tstep=5000, x0=x0)

    vout = yr[:,2]

    # Perform FFT
    freq, X, N = trimmed_fft(vout, tr, 100e-6, vin_mag/2)

    fig, ax1 = plt.subplots()
    plt.grid()

#    ax1.plot(tr, yr[:,2],       color='blue',  label="$V_{OUT}$")

    # Plot FFT
    ax1.plot(freq[:N//2], X[:N//2], label="$Vlo$ = "+ str(vin_mag))

    ax1.set_ylabel('v')
    ax1.legend()
    plt.show(block=False)

    fig, ax1 = plt.subplots()
    plt.grid()

    ax1.plot(tr, yr[:,2],       color='blue',  label="$V_{OUT}$")

    ax1.set_ylabel('v')
    ax1.legend()
    plt.show()

if __name__ == "__main__":

    main()
