import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../networks")
from gilbert_cell import *

from scipy import optimize
from scipy.integrate import ode

import numpy as np

import itertools

import matplotlib.pyplot as plt

import math

import control as ct

def dypc_litmus0(t, sys, ckt):

    root_start = np.zeros(ckt.num_edges)
    vs = 1

    root = optimize.root(circuit_eqn_dyn, root_start, args=(sys, ckt), tol=1e-9)
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

def circuit_eqn_dyn(x, sys, ckt):

    # current Im
    I = ckt.get_im(x=x, sys=sys, t=0)

    # voltage Vm
    V = ckt.get_vm(x=x, sys=sys, t=0)

    # get dependent currents
    I = ckt.get_im_dep()

    # Copy matrices
    qf_num = ckt.qf.copy()
    bf_num = ckt.bf.copy()

    eqn_qf = np.dot(qf_num, I)
    eqn_bf = np.dot(bf_num, V)

    return eqn_qf + eqn_bf

def circuit_eqn(x, ckt):

    # current Im
    I = ckt.get_im(x=x, sys=0, t=0)

    # voltage Vm
    V = ckt.get_vm(x=x, sys=0, t=0)

    # get dependent currents
    I = ckt.get_im_dep()

    # Copy matrices
    qf_num = ckt.qf.copy()
    bf_num = ckt.bf.copy()

    eqn_qf = np.dot(qf_num, I)
    eqn_bf = np.dot(bf_num, V)

    return eqn_qf + eqn_bf

def main():

    # Operating Point
    # ---------------
    ckt_op = Circuit()
    nw = gilbert_cell_full(ckt_op)

    ckt_op.init_circuit()

    scb = Scoreboard(num_edges=ckt_op.num_edges, num_sys_vars=ckt_op.num_sys_vars, degen_mtrx=ckt_op.degen_mtrx)

    vrf_sweep = np.linspace(1.2, 1.8, 20)
    vlo_sweep = np.linspace(3.7, 4.3, 20)

    vf_out     = np.zeros([len(vlo_sweep), len(vrf_sweep)])

    vlo_bias = 4.0
    vrf_bias = 1.5
    iref_val = 200e-6

    nw.set_vlo_bias(bias=vlo_bias)
    nw.set_vrf_bias(bias=vrf_bias)
    nw.set_iref(iref_val)

    for idx_vlo, vlo in enumerate(vlo_sweep):
        nw.set_vlo(vlo)

        root_start = np.ones(ckt_op.num_edges)
        for idx, vrf in enumerate(vrf_sweep):

            nw.set_vrf(vrf=vrf)

            # Operating Point
            root       = optimize.root(circuit_eqn, root_start, args=(ckt_op), tol=1e-9)

            scb.x = root.x
            scb.t = 0

            root_start = root.x

            saturation_region = nw.check_saturation_region(scb)

            ref_vout = (vlo - vlo_bias)

            if ref_vout != 0:
                vf_out[idx_vlo][idx] = nw.get_vf_output(ckt_op, scb) / ref_vout
            else:
                vf_out[idx_vlo][idx] = nw.get_vf_output(ckt_op, scb)

            if not root.success:
                print ("Failed to calculate small signal circuit")

            if not saturation_region:
                print ("Not in saturation region " + str(vlo) + " " + str(vrf))

    fig, ax1 = plt.subplots()
    #
    for idx_vlo, vlo in enumerate(vlo_sweep):
        ax1.plot(vrf_sweep, vf_out[idx_vlo])
    #ax1.plot(T, yout[1][0])
    #ax2.plot(tr, yr)
    plt.show()


if __name__ == "__main__":

    main()
