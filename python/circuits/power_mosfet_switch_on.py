import sys, getopt, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

from diode import diode_id
from power_mosfet import power_mosfet

from scipy import optimize
from scipy.integrate import ode

import numpy as np

import itertools

import matplotlib.pyplot as plt

import math

import control as ct

from sympy import *

import re as regex

def add_power_nmos(g, d, s, KP, VTO, L_LAMBDA, L, W, Cgs_val, Cgd_val, ckt):

    d_int = -1
    s_int = -2

    ids = power_mosfet()
    ids.set_instance("IDS")
    ids.set_params(KP=KP, VTO=VTO, L_LAMBDA=L_LAMBDA, L=L, W=W)

    vgs = current_src()
    vgs.set_type(ElementType.current_src)
    vgs.set_instance("VGS")

    ids_ref = ckt.add_edge(d_int, s_int, ids)
    vgs_ref = ckt.add_edge(g,     s_int, vgs)

    ids.set_vgs_ref(vgs_ref=vgs_ref)

    # Add Cgs
    Cgs = capacitor()
    Cgs.set_value(Cgs_val)
    Cgs.set_instance("Cgs")
    ckt.add_edge(g, s_int, Cgs)

    # Add Cgd
    Cgd = capacitor()
    Cgd.set_value(Cgd_val)
    Cgd.set_instance("Cgd")
    ckt.add_edge(d_int, g, Cgd)

    # Add Rd
    Rd = resistor()
    Rd.set_value(0.0017997)
    Rd.set_instance("Rd")
    ckt.add_edge(d, d_int, Rd)

    # Add Rs
    Rs = resistor()
    Rs.set_value(0.014066)
    Rs.set_instance("Rs")
    ckt.add_edge(s, s_int, Rs)

def add_circuit(ckt, rser):

    Vs_val     = 12.0
    Vdrive_val = 12.0
    R2_val     = 100

    # Add Vs
    Vs = voltage_src(ramp=False, ramp_ddt=1e3, delay=5e-6)
    Vs.set_value(Vs_val)
    Vs.set_instance("Vs")
    Vs.set_is_input()
    ckt.add_edge(1, 0, Vs)

    # Add Iload
    Iload = current_src()
    Iload.set_is_const()
    Iload.set_value(100e-3)
    Iload.set_instance("Iload")
    ckt.add_edge(1, 2, Iload)

    # Add freewheeling diode
    D1 = diode_id()
    D1.set_instance("D1")
    D1.set_params(IS=75e-12, N=1.0, EG=0.7, CJO=26e-12, M=0.5, IBV=5e-6, BV=400, TT=4.2e-6)
    ckt.add_edge(2, 1, D1)

    # Add power nmos M1
    add_power_nmos(g=3, d=2, s=0, KP=67.9211, VTO=2.08819, L_LAMBDA=0.0038193, L=100e-6, W=100e-6, Cgs_val=5e-9, Cgd_val=20e-9, ckt=ckt)

    # Add R2
    R2 = resistor()
    R2.set_value(R2_val)
    R2.set_instance("R2")
    ckt.add_edge(4, 3, R2)

    # Add Vdrive
    Vdrive = voltage_src(ramp=True, ramp_ddt=1e5, delay=50e-6)
    Vdrive.set_value(Vdrive_val)
    Vdrive.set_instance("Vdrive")
    Vdrive.set_is_input()
    ckt.add_edge(4, 0, Vdrive)

def circuit_eqn(x, sys, ckt, vs, t):

    # current Im
    I = ckt.get_im(x=x, sys=sys, t=t)

    # voltage Vm
    V = ckt.get_vm(x=x, sys=sys, t=t)

    # get dependent currents
    I = ckt.get_im_dep()

    #I = ckt.get_degen()

    # Copy matrices
    qf_num = ckt.qf.copy()
    bf_num = ckt.bf.copy()

    eqn_qf = np.dot(qf_num, I)
    eqn_bf = np.dot(bf_num, V)

    return eqn_qf + eqn_bf

def dypc_litmus0(t, sys, ckt):

    root_start = np.ones(ckt.num_edges)
    vs = 1

    root = optimize.root(circuit_eqn, root_start, args=(sys, ckt, vs, t), tol=1e-3)
    if not root.success:
        print (root.x)
        sys.exit(0)

    return ckt.get_dy(x=root.x)

def ode_solve(ckt):

    num_sys_vars    = ckt.get_num_sys_vars()

    tstep = 10000

    t     = np.linspace(0, 150e-6, tstep)
    x0    = np.zeros (num_sys_vars)

    r = ode(dypc_litmus0).set_integrator('lsoda', method='bdf', atol=1e-3, rtol=1e-3)
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


def main():

    ckt = Circuit()

    rser = False

    add_circuit(ckt, rser)

    # Initilise the Circuit
    ckt.init_circuit()

    # Plot numerical solver response
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    # Calculate the step response using a numerical ODE solver
    t, y = ode_solve(ckt)

    for sys_var in range(len(y[0,:])):

        inst = str(ckt.get_edge_info_from_sys_var_ref(sys_var).instance)

        if inst[0] == "C":
            l = "$V_{" + inst + "}$"
            ax1.plot(t, y[:,sys_var], label=l)
        else:
            l = "$i_{" + inst + "}$"
            ax2.plot(t, y[:,sys_var], label=l, linestyle="dotted")

    ax1.legend()
    ax2.legend(loc="lower right")

    ax1.set_ylabel("$V (V)$")
    ax1.set_xlabel("$t (S)$")

    ax2.set_ylabel("$i (A)$")
    ax2.set_xlabel("$t (S)$")

    if rser == False:
        plt.savefig("litmus0_rser0.png", bbox_inches = 'tight')
    else:
        plt.savefig("litmus0_rser1.svg", bbox_inches = 'tight')

    plt.show()

if __name__ == "__main__":

    main()
