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

def main():

    # Operating Point
    # ---------------

    vlo_omega = 1e6
    vlo_mag   = 0.15
    vlo_bias  = 4.0

    vrf_omega = 5e6
    vrf_mag   = 0.15
    vrf_bias  = 1.5

    iref_val = 200e-6

    ac_params = {"vlo_omega": vlo_omega, "vlo_mag": vlo_mag, "vlo_bias": vlo_bias, \
                 "vrf_omega": vrf_omega, "vrf_mag": vrf_mag, "vrf_bias": vrf_bias}
    nw = gilbert_cell_ntwrk(sine_src=True, **ac_params)
    nw.set_iref(iref_val)


    x0    = 4 * np.ones (nw.ckt.get_num_sys_vars())


    tr, yr = ode_solve(nw.ckt, tend=10e-6, tstep=10000, x0=x0)

    fig, ax1 = plt.subplots()
    #
    #for idx_vlo, vlo in enumerate(vlo_sweep):
    #    ax1.plot(vrf_sweep, vf_out[idx_vlo])
    #ax1.plot(T, yout[1][0])
    ax1.plot(tr, yr)
    plt.show()


if __name__ == "__main__":

    main()
