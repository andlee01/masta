import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../networks")
from cross_coupled_lc_osc_ntwrk import *

from scipy import optimize
from scipy.integrate import ode
from scipy.signal import freqresp, StateSpace

import numpy as np

import itertools

import matplotlib.pyplot as plt

import math

import control as ct

def dypc_litmus0(t, sys, ckt):

    root_start = np.ones(ckt.num_edges)
    vs = 1

    root = optimize.root(circuit_eqn, root_start, args=(sys, ckt, t), tol=1e-6)
    root_start = root.x
    if not root.success:
        print (t)
        #sys.exit(0)

    return ckt.get_dy(x=root.x)

def ode_solve(ckt, tend=50e-12, tstep=10000, x0=0):

    num_sys_vars    = ckt.get_num_sys_vars()

    t     = np.linspace(0, tend, tstep)

    r = ode(dypc_litmus0).set_integrator('lsoda', method='bdf', atol=1e-9, rtol=1e-9)
    r.set_initial_value(x0, 0.0)
    r.set_f_params(ckt)

    y    = np.empty((tstep, num_sys_vars))
    y[0] = x0

    vx_out = np.zeros(len(t))

    k = 1
    while r.successful() and k < tstep:
        r.integrate(t[k])

        y[k] = r.y

        k += 1

    return t, y

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

def op_solve(nw, root_start):

    root = optimize.root(circuit_eqn, root_start, args=(0, nw.ckt_op, 0), tol=1e-9)
    if not root.success:
        root_start = root.x
        op_solve(nw, root_start)

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


def main():


    R = [2e-3, 2, 20, 200, 1e3, 1.2e3, 2e3, 3e3]

    fstep    = 10000
    fstart   = 2 * math.pi * 1e3
    fstop    = 2 * math.pi * 10e6
    omega    = np.linspace(fstart, fstop, fstep)

    mag_lst   = []
    phase_lst = []
    omega_lst = []

    for R_idx, R_osc in enumerate(R):

        # HERE
        params = {"L": 10e-3, "R":R_osc, "C":1e-9}

        nw_lc_tank = cross_coupled_lc_osc_ntwrk(**params)

        # Operating point
        nw_lc_tank.add_op_ckt()
        root_start = np.ones(nw_lc_tank.ckt_op.num_edges)
        _, g = op_solve(nw=nw_lc_tank, root_start=root_start)

        op = nw_lc_tank.ckt_op.scb
        op.x = g

        # Small signal
        nw_lc_tank.add_sml_ckt(op=op)

        # State Space
        nw_lc_tank.ckt_sml.get_ss()

        x0    = np.zeros (nw_lc_tank.ckt_sml.get_num_sys_vars())
        x0 = nw_lc_tank.ckt_sml.get_lti_const_x0(x0)
        sys = ct.ss(nw_lc_tank.ckt_sml.A, nw_lc_tank.ckt_sml.B, nw_lc_tank.ckt_sml.C, nw_lc_tank.ckt_sml.D)

        single_out = True

        mag, phase, omega_out = ct.freqresp(sys, omega=omega)

        if single_out:
            mag_out = 20 * np.log10(np.abs(mag))
            phase_out = phase
        else:
            mag_out   = 20 * np.log10(np.abs(mag[1, 0, :]))
            phase_out = phase[1, 0, :]

        mag_lst.append(mag_out)
        phase_lst.append(phase_out)
        omega_lst.append(omega_out)

    fig, ax1 = plt.subplots()
    plt.grid()
    for R_idx, R_osc in enumerate(R):
        ax1.semilogx(omega_lst[R_idx], mag_lst[R_idx], label="$R = $" + str(R_osc))
    ax1.legend()
    plt.show(block=False)

    fig, ax1 = plt.subplots()
    plt.grid()
    for R_idx, R_osc in enumerate(R):
        ax1.semilogx(omega_lst[R_idx], phase_lst[R_idx], label="$R = $" + str(R_osc))
    ax1.legend()
    plt.show(block=True)


    #plt.xscale("log")
    ##plt.yscale("log")
    ##ax1.plot(omega_out, mag_out)
    #if single_out:
    #    ax1.plot(omega_out, phase)
    #else:
    #    ax1.plot(omega_out, phase[0][0])
    #plt.show(block=True)



if __name__ == "__main__":

    main()
