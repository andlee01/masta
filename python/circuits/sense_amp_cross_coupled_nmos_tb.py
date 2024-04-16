import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../networks")
from sense_amp_cross_coupled_nmos import *

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
    nw = sense_amp_cross_coupled_nmos(ckt_op)

    ckt_op.init_circuit()

    # Operating Point
    root_start = np.ones(ckt_op.num_edges)
    root       = optimize.root(circuit_eqn, root_start, args=(ckt_op), tol=1e-9)
    op         = ckt_op.scb

    if not root.success:
        print ("Failed to calculate small signal circuit")

    # Large Signal
    # ---------------
    ckt_lrg_dyn = Circuit()
    nw_lrg_dyn = sense_amp_cross_coupled_nmos(ckt_lrg_dyn, trans=True)

    ckt_lrg_dyn.init_circuit()

    x0    = np.zeros (ckt_lrg_dyn.get_num_sys_vars())
    for x in range(len(x0)):
        x0[x] = 2.5
    x0[0] -= 0.05

    tr, yr = ode_solve(ckt_lrg_dyn, tend=50e-12, tstep=10000, x0=x0)

    # Small signal
    # -------------
    ckt_small = Circuit(lti=True)
    nw.add_sml_ckt(ckt_sml=ckt_small, op=op, dyn=True)
    ckt_small.init_circuit()
    ckt_small.get_ss()

    sys = ct.ss(ckt_small.A, ckt_small.B, ckt_small.C, ckt_small.D)

    tstep = 10000
    t     = np.linspace(0, 50e-12, tstep)
    x0    = np.zeros (ckt_small.get_num_sys_vars())
    x0[1] = -0.05

    T, yout = ct.step_response(sys, T=t, X0=x0)

    fig, (ax1, ax2) = plt.subplots(2, 1)

    ax1.plot(T, yout[0][0])
    ax1.plot(T, yout[1][0])
    ax2.plot(tr, yr)
    plt.show()


if __name__ == "__main__":

    main()
