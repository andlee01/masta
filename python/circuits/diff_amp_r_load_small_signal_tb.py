import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../networks")
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

def intermodulation_test(omega, omega_1, mag, bias, iref_val, tend):

    two_lrg_ac_params = {"omega_0": omega, "omega_1": omega_1, "mag": mag/2, "bias": bias}
    nw = diff_amp_r_load_ntwrk(two_sine_src=True, **two_lrg_ac_params)
    nw.set_iref(iref_val)

    # Transient for large signal
    x0    = 4 * np.ones (nw.ckt.get_num_sys_vars())
    tr, yr = ode_solve(nw.ckt, tend=tend, tstep=10000, x0=x0)

    # Calculate Vf_out
    vf_out = yr[:,0] - yr[:,1]

    # Perform FFT
    freq, X, N = trimmed_fft(vf_out, tr, 1e-6, mag/2)

    return tr, yr, freq, X, vf_out, N

def op_solve(nw, vlop_mag, vlon_mag):

    nw.set_vlo_dc(vlop=vlop_mag, vlon=vlon_mag)

    root_start = np.ones(nw.ckt.num_edges)
    root = optimize.root(circuit_eqn_dyn, root_start, args=(0, nw.ckt, 0), tol=1e-9)
    if not root.success:
        print (root.x)
        sys.exit(0)

    vds_vref_m1 = nw.get_nmos_m1_vltg(root.x)[0]
    vgs_vref_m1 = nw.get_nmos_m1_vltg(root.x)[1]

    vds_vref_m2 = nw.get_nmos_m2_vltg(root.x)[0]
    vgs_vref_m2 = nw.get_nmos_m2_vltg(root.x)[1]

    return vds_vref_m1, vgs_vref_m1, vds_vref_m2, vgs_vref_m2

def plot_mos(vgs_ref, vds_vref, vds_sweep, ids_copy, block, ref):

    fig, ax1 = plt.subplots()
    plt.grid()

    scb = Scoreboard(num_edges=3, num_sys_vars=0, degen_mtrx=0)

    i_ds    = np.zeros([len(vgs_ref), len(vds_sweep)])

    for vgs_idx, vgs in enumerate(vgs_ref):
        for vds_idx, vds in enumerate(vds_sweep):

            x = [vds, vgs]
            scb.v = x
            ids_copy.get_dependent_current(scb=scb)
            i_ds[vgs_idx][vds_idx] = scb.i[2]

        if vgs_idx == 0:
            color = "red"
        else:
            color = "green"

        if ref == 1:
            base_str_vgs = "$V_{GS}M_1 = $"
            base_str_vds = "$V_{DS}M_1$"
        else:
            base_str_vgs = "$V_{GS}M_2 = $"
            base_str_vds = "$V_{DS}M_2$"


        l_vgs = "{:.5f}".format(vgs)
        ax1.plot(vds_sweep, i_ds[vgs_idx], color=color, label=base_str_vgs + l_vgs)
        ax1.axvline(x=vds_vref[vgs_idx], color=color, linestyle="dotted", label=base_str_vds)
        ax1.legend()

    plt.show(block=block)

def main():

    iref_val = 40e-6

    vlop_mag = 2.5
    vlon_mag = 2.4

    dc_params = {"vlop_mag": 2.5, "vlon_mag": 2.5}
    nw = diff_amp_r_load_ntwrk(dc_src=True, **dc_params)
    nw.set_iref(iref_val)

    vgs_ref_m1 = np.zeros(2)
    vds_vref_m1 = np.zeros(2)
    vgs_ref_m2 = np.zeros(2)
    vds_vref_m2 = np.zeros(2)
    vds_sweep  = np.linspace(0.0, 4.8, 1000)
    i_ds_m1    = np.zeros([len(vgs_ref_m1), len(vds_sweep)])

    vds_vref_m1[0], vgs_ref_m1[0], vds_vref_m2[0], vgs_ref_m2[0] = op_solve(nw, vlop_mag=2.45, vlon_mag=2.5)
    vds_vref_m1[1], vgs_ref_m1[1], vds_vref_m2[1], vgs_ref_m2[1] = op_solve(nw, vlop_mag=2.55, vlon_mag=2.5)

    ids_copy = nw.get_nmos_m1_copy()

    scb = Scoreboard(num_edges=3, num_sys_vars=0, degen_mtrx=0)

    plot_mos(vgs_ref_m1, vds_vref_m1 ,vds_sweep, ids_copy, False, 1)
    plot_mos(vgs_ref_m2, vds_vref_m2 ,vds_sweep, ids_copy, True, 2)

if __name__ == "__main__":

    main()
