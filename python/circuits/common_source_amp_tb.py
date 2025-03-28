import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../networks")
from common_source_amp_ntwrk import *

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

    root = optimize.root(circuit_eqn, root_start, args=(sys, nw.ckt, t), tol=1e-6)
    root_start = root.x
    if not root.success:
        print (t)
        #sys.exit(0)

    return nw.ckt.get_dy(x=root.x)

def ode_solve(nw, tend=50e-12, tstep=10000, x0=0):

    num_sys_vars    = nw.ckt.get_num_sys_vars()

    t     = np.linspace(0, tend, tstep)

    r = ode(dypc_litmus0).set_integrator('lsoda', method='bdf', atol=1e-9, rtol=1e-9)
    r.set_initial_value(x0, 0.0)
    r.set_f_params(nw)

    y    = np.empty((tstep, num_sys_vars))
    y[0] = x0

    vx_out = np.zeros(len(t))

    k = 1
    while r.successful() and k < tstep:
        r.integrate(t[k])

        y[k] = r.y

        vx_out[k] = nw.get_vx(y[k], t[k])

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

def op_solve(nw, vin, root_start):

    nw.set_vin(vin=vin)

    root = optimize.root(circuit_eqn, root_start, args=(0, nw.ckt, 0), tol=1e-9)
    if not root.success:
        root_start = root.x
        op_solve(nw, vin, root_start)

    vout = nw.get_vout(x=root.x)

    return vout, root.x

def plot_mos(vgs, vds_sweep, ids_copy):

    scb = Scoreboard(num_edges=3, num_sys_vars=0, degen_mtrx=0)

    i_ds    = np.zeros(len(vds_sweep))

    for vds_idx, vds in enumerate(vds_sweep):

        x = [vds, vgs]
        scb.v = x
        ids_copy.get_dependent_current(scb=scb)
        i_ds[vds_idx] = scb.i[2]

    return i_ds

def plot_rd(vds_sweep, Rd, Rs, Vs):

    i_rd = np.zeros(len(vds_sweep))

    fact = Rd / (Rd+Rs)

    for vds_idx, vds in enumerate(vds_sweep):
        i_rd[vds_idx] = ((Vs - vds)*fact) / Rd

    return i_rd

def plot_vin_reponse(nw_degen, vin_sweep, vin_op_sweep, degen):

    vout_degen = np.zeros(len(vin_sweep))
    ids_degen  = np.zeros(len(vin_sweep))
    degen_op   = np.zeros([len(vin_op_sweep), 3])

    # Solve for swept Vin
    off_region = True
    sat_region = False
    tri_region = False

    root_start = np.ones(nw_degen.ckt.num_edges)
    for vin_idx, vin in enumerate(vin_sweep):
        vout_degen[vin_idx], x = op_solve(nw_degen, vin, root_start)
        ids_degen[vin_idx] = nw_degen.get_op(x)[2]

        if nw_degen.check_saturation_region(x):
            if not sat_region:
                sat_op = nw_degen.get_op(x)
                sat_region_vin = vin
            sat_region = True

        if nw_degen.check_triode_region(x):
            if not tri_region:
                tri_op = nw_degen.get_op(x)
                tri_region_vin = vin
            tri_region = True

        root_start = x

    # Solve for Vin operating points (selected)
    vin_op_sweep = np.linspace(0.0, 5.0, 11)
    degen0_op    = np.zeros([len(vin_op_sweep), 3])

    root_start = np.ones(nw_degen.ckt.num_edges)
    for vin_idx, vin in enumerate(vin_op_sweep):
        _, x = op_solve(nw_degen, vin, root_start)
        degen_op[vin_idx] = nw_degen.get_op(x)

    fig, ax1 = plt.subplots()
    plt.grid()

    ax1.plot(vin_sweep, vout_degen,       color='blue',  label="$V_{OUT}$")
    ax1.set_ylabel("$V_{out}$")
    ax1.set_xlabel("$V_{in}$")

    plt.fill_between(vin_sweep, 6, where=(vin_sweep < sat_region_vin), color='lightgreen', alpha=0.5, label="Off")
    plt.fill_between(vin_sweep, 6, where=(vin_sweep > tri_region_vin), color='lightblue', alpha=0.5, label="Triode")

    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.2, 0.8))

    ax2 = ax1.twinx()
    ax2.plot(vin_sweep, ids_degen,       color='red',  label="$i_{DS}$")
    ax2.set_ylabel("$i_{DS}$")
    ax2.legend(loc="upper left", bbox_to_anchor=(1.2, 1))

    if degen:
        plt.savefig("../../doc/common_source_amp_degen_transfer.svg", bbox_inches = 'tight')
    else:
        plt.savefig("../../doc/common_source_amp_transfer.svg", bbox_inches = 'tight')
    #plt.show(block=False)

    return sat_op, tri_op, degen_op

def plot_load_line(ids, vds_sweep, vin_op_sweep, degen_op, tri_op, Rd, Rs, block):

    i_ds_m1     = np.zeros([len(vin_op_sweep), len(vds_sweep)])
    i_rd        = np.zeros(len(vds_sweep))

    fig, ax1 = plt.subplots()
    plt.grid()

    for vin_idx, vin in enumerate(vin_op_sweep):

        vgs = degen_op[vin_idx][1]

        if vgs < tri_op[1]:
            i_ds_m1[vin_idx] = plot_mos(vgs, vds_sweep, ids)
            l_vgs = "{:.5f}".format(vgs)
            l_vin = "{:.5f}".format(vin)
            ax1.plot(vds_sweep, i_ds_m1[vin_idx],  label="$V_{GS} = $" + l_vgs + ", $V_{in} = $" + l_vin)

    i_ds_m1_tri_op = plot_mos(tri_op[1], vds_sweep, ids)
    l_vgs = "{:.5f}".format(tri_op[1])
    ax1.plot(vds_sweep, i_ds_m1_tri_op,  label="$V_{GS_{tri}} = $" + l_vgs)
    ax1.legend(loc='upper left', bbox_to_anchor=(1, 1))


    i_rd = plot_rd(vds_sweep=vds_sweep, Rd=Rd, Rs=Rs, Vs=5.0)
    ax1.plot(vds_sweep, i_rd,  label="$i_{RD}$")
    ax1.set_ylabel("$i_{DS}$")
    ax1.set_xlabel("$V_{DS}$")

    plt.fill_between(vds_sweep, 2.5e-3, where=(vds_sweep < tri_op[0]), color='lightblue', alpha=0.5)

    for i in range(len(degen_op)):

        if degen_op[i][0] > tri_op[0]:
            if i < (len(degen_op) - 1):
                plt.vlines(x=degen_op[i][0], ymin=0, ymax=degen_op[i+1][2] * 1.2, color='red', linestyle='--')
            else:
                plt.vlines(x=degen_op[i][0], ymin=0, ymax=degen_op[i][2] * 1.2, color='red', linestyle='--')

        if i >= 3 and degen_op[i][0] > tri_op[0]:

            plt.ylim(0, degen_op[i][2] * 1.8)

            # Annotate with a horizontal arrow between the two vertical lines
            # Set `xy` at one line and `xytext` at the other line
            y_pos = degen_op[i][2]

            plt.annotate('', xy=(degen_op[i][0], y_pos), xytext=(degen_op[i-1][0], y_pos),
                         arrowprops=dict(arrowstyle='<->', color='black', lw=2))  # <-> creates a bidirectional arrow

            # Add the distance label between the lines
            midpoint = (degen_op[i][0] + degen_op[i-1][0]) / 2  # Midpoint between the two vertical lines
            plt.text(midpoint, y_pos + 0.05e-4, "$\\Delta V_{{DS}} = {:.5f}$".format(degen_op[i-1][0] - degen_op[i][0]), ha='center')

    if Rs != 0:
        plt.savefig("../../doc/common_source_amp_degen_loadline.svg", bbox_inches = 'tight')
    else:
        plt.savefig("../../doc/common_source_amp_loadline.svg", bbox_inches = 'tight')
    #plt.show(block=block)

def calc_small_signal_gain(nw, vin, vin_delta):

    vin_sweep    = np.linspace(vin - vin_delta/2, vin + vin_delta/2, 50)
    vout         = np.zeros(len(vin_sweep))

    root_start = np.ones(nw.ckt.num_edges)
    _, x = op_solve(nw, vin, root_start)
    op = nw.ckt.scb
    nw.add_sml_ckt(op=op)

    for i in range(len(vin_sweep)):
        root_start = np.ones(nw.ckt_sml.num_edges)
        _, x = op_solve(nw, vin_sweep[i], root_start)

        vout[i] = nw.get_sml_vout(x)

    gain = (vout[-1] - vout[0]) / vin_delta

    return vout, vin_sweep, gain

def calc_small_signal_gm(nw, vin_sweep, num):

    gm = np.zeros(len(vin_sweep))

    for i in range(len(vin_sweep)):
        root_start = np.ones(nw.ckt.num_edges)
        _, x = op_solve(nw, vin_sweep[i], root_start)
        op = nw.ckt.scb
        op.x = x
        op.i = nw.ckt.get_im(x=x, sys=0, t=0)
        op.v = nw.ckt.get_vm(x=x, sys=0, t=0)
        nw.add_sml_ckt(op=op)

        gm[i], _ = nw.get_op_sml()

        KP = 120e-6
        L = 10
        W = 20

        beta = KP * (W / L)

        #if num == 0:
        #gm[i] = math.sqrt(2 * gm[i] * beta)

    return gm

def plot_small_signal_gain(vin, vout_degen0, vout_degen1):

    fig, ax1 = plt.subplots()
    plt.grid()

    ax1.plot(vin, vout_degen0,  label="$V_{OUT_0}$")
    ax1.plot(vin, vout_degen1,  label="$V_{OUT_1}$")

    ax1.legend(loc='upper left', bbox_to_anchor=(1, 1))

    plt.show(block=False)

def plot_small_signal_gm(vin, gm0, gm1):

    fig, ax1 = plt.subplots()
    plt.grid()

    ax1.plot(vin, gm0,  label="$g_{m_0}$")
    ax1.plot(vin, gm1,  label="$g_{m_1}$")

    ax1.legend(loc='upper left', bbox_to_anchor=(1, 1))

    ax1.set_ylabel("$gm$")
    ax1.set_xlabel("$V_{in}$")

    plt.savefig("../../doc/common_source_amp_gm.svg", bbox_inches = 'tight')

def main():


    # HERE
    params = {"Rs": 5e3, "Rd":5e3}

    vin_sweep    = np.linspace(0.0, 5.0, 1000)
    vin_op_sweep = np.linspace(0.0, 5.0, 11)
    vds_sweep    = np.linspace(0.0, 5.0, 1000)

    nw_degen0 = common_source_amp_ntwrk(src_degen=False, **params)
    nw_degen1 = common_source_amp_ntwrk(src_degen=True,  **params)

    sat0_op, tri0_op, degen0_op = plot_vin_reponse(nw_degen=nw_degen0, vin_sweep=vin_sweep, vin_op_sweep=vin_op_sweep, degen=False)
    sat0_op, tri1_op, degen1_op = plot_vin_reponse(nw_degen=nw_degen1, vin_sweep=vin_sweep, vin_op_sweep=vin_op_sweep, degen=True)

    vds_sweep   = np.linspace(0.0, 5.0, 1000)


    nmos_params = {"KP"     : 120e-6,
                   "vth"    : 0.8,
                   "lambda" : 0.01,
                   "L"      : 10,
                   "W"      : 20}

    ids = vccs_l1_mosfet()
    ids.set_params(KP=nmos_params["KP"], \
                   vth=nmos_params["vth"], \
                   l_lambda=nmos_params["lambda"], \
                   L=nmos_params["L"], \
                   W=nmos_params["W"])

    ids.set_ref(0)
    ids.set_vgs_ref(1)
    ids.i_d_ref = 2


    plot_load_line(ids=ids, vds_sweep=vds_sweep, vin_op_sweep=vin_op_sweep, degen_op=degen0_op, tri_op=tri0_op, Rd=5e3, Rs=0, block=False)

    plot_load_line(ids=ids, vds_sweep=vds_sweep, vin_op_sweep=vin_op_sweep, degen_op=degen1_op, tri_op=tri1_op, Rd=5e3, Rs=5e3, block=True)

    # Small signal
    # ------------
    #vout_degen0, vin, gain0 = calc_small_signal_gain(nw=nw_degen0, vin=2.5, vin_delta=0.2)
    #vout_degen1, _  , gain1 = calc_small_signal_gain(nw=nw_degen1, vin=2.5, vin_delta=0.2)
    #
    #plot_small_signal_gain(vin=vin, vout_degen0=vout_degen0, vout_degen1=vout_degen1)

    gm0 = calc_small_signal_gm(nw=nw_degen0, vin_sweep=vin_sweep, num=0)
    gm1 = calc_small_signal_gm(nw=nw_degen1, vin_sweep=vin_sweep, num=1)

    plot_small_signal_gm(vin=vin_sweep, gm0=gm0, gm1=gm1)

if __name__ == "__main__":

    main()
