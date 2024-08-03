import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../networks")
from cs_stage_bias_ntwrk import *

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

    root = optimize.root(circuit_eqn, root_start, args=(sys, nw.ckt, t), tol=1e-8)
    root_start = root.x
    if not root.success:
        print (t)
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

def main():

    vin_omega = 1e6
    vin_mag   = 0.1

    ac_params = {"omega": vin_omega, "mag": vin_mag}

    nw = cs_stage_bias_ntwrk(sine_src=True, **ac_params)

    print(nw.ckt.num_sys_vars)

    x0    = np.zeros (nw.ckt.get_num_sys_vars())
    tr, yr, v_vx_out = ode_solve(nw, tend=500e-6, tstep=10000, x0=x0)

    fig, ax1 = plt.subplots()
    plt.grid()

    ax1.plot(tr, v_vx_out,       color='blue',  label="$V_{OUT}$")
    ax1.set_ylabel('v')
    ax1.legend()
    plt.show()

if __name__ == "__main__":

    main()
