import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../networks")
from gilbert_cell_ntwrk import *

from scipy import optimize
from scipy.integrate import ode

import numpy as np

import itertools

import matplotlib.pyplot as plt

import math

import control as ct

def dypc_litmus0(t, sys, nw):

    root_start = np.ones(nw.ckt.num_edges)
    vs = 1

    root = optimize.root(circuit_eqn_dyn, root_start, args=(sys, nw.ckt, t), tol=1e-9)
    root_start = root.x
    if not root.success:
        print (t)

        if nw.check_switching_cutoff(nw.scb.x):
            print ("cut-off")
        sys.exit(0)

    return nw.ckt.get_dy(x=root.x)

def ode_solve(nw, tend=50e-12, tstep=10000, x0=0):

    num_sys_vars    = nw.ckt.get_num_sys_vars()

    t     = np.linspace(0, tend, tstep)

    r = ode(dypc_litmus0).set_integrator('lsoda', method='bdf', atol=1e-9, rtol=1e-9)
    r.set_initial_value(x0, 0.0)
    r.set_f_params(nw)

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

    # Operating Point
    # ---------------

    vlo_omega = 1e6 * math.pi * 2
    vlo_mag   = 0.15
    vlo_bias  = 4.0

    vrf_omega = 5e6 * math.pi * 2
    vrf_mag   = 0.15
    vrf_bias  = 1.5

    iref_val = 200e-6

    #tend = 10e-6
    #tstep = 1000
    #t     = np.linspace(0, tend, tstep)
    #
    #v = np.zeros(len(t))
    #x = np.zeros(len(t))
    ##v = vlo_bias + (vlo_mag * (0.5 * signal.square(vlo_omega * t + math.pi) + 0.5) )
    ##x = vlo_bias + (vlo_mag * (0.5 * signal.square(vlo_omega * t          ) + 0.5) )
    #
    #for t_idx, t_step in enumerate(t):
    #    v[t_idx] = cunt(t_step, 1e-6, 100e-12, 100e-12, 0.6e-6, 0.25e-6)
    #    x[t_idx] = cunt(t_step, 1e-6, 100e-12, 100e-12, 0.6e-6, 0.75e-6)
    #
    #
    #
    #
    #
    #fig, ax1 = plt.subplots()
    #ax1.plot(t, v)
    #ax1.plot(t, x)
    #plt.show()

    ac_params = {"vlo_omega": vlo_omega, "vlo_mag": vlo_mag, "vlo_bias": vlo_bias, \
                 "vrf_omega": vrf_omega, "vrf_mag": vrf_mag, "vrf_bias": vrf_bias}
    nw = gilbert_cell_ntwrk(sine_src=True, **ac_params)
    nw.set_iref(iref_val)


    x0    = 4 * np.ones (nw.ckt.get_num_sys_vars())


    tr, yr = ode_solve(nw, tend=50e-6, tstep=10000, x0=x0)

    # Calculate Vf_out
    vf_out = yr[:,0] - yr[:,1]

    # Perform FFT
    freq, X, N = trimmed_fft(vf_out, tr, 1e-6, vrf_mag/2)

    fig, ax1 = plt.subplots()
    #
    #for idx_vlo, vlo in enumerate(vlo_sweep):
    #    ax1.plot(vrf_sweep, vf_out[idx_vlo])
    #ax1.plot(T, yout[1][0])
    ax1.plot(tr, vf_out)
    plt.show(block=False)

    fig, ax1 = plt.subplots()

    # Plot FFT
    ax1.plot(freq[:N//2], X[:N//2])
    plt.show()

if __name__ == "__main__":

    main()
