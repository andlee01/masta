import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

from scipy import optimize
from scipy.integrate import ode

import numpy as np

import itertools

import matplotlib.pyplot as plt

import math

nmos_params = {"KP"     : 120e-6,
               "vth"    : 0.8,
               "lambda" : 0.01,
               "L"      : 2,
               "W"      : 10}

pmos_params = {"KP"     : 40e-6,
               "vth"    : 0.9,
               "lambda" : 0.0125,
               "L"      : 2,
               "W"      : 30}

vgs_1_ref = 0

def add_nmos(g, d, s, KP, vth, l_lambda, L, W, ckt):
    ids = vccs_l1_mosfet()
    ids.set_type(ElementType.current_src)
    ids.set_instance("IDS")
    ids.set_params(KP=KP, \
                   vth=vth, \
                   l_lambda=l_lambda, \
                   L=L, \
                   W=W)

    vgs = current_src()
    vgs.set_type(ElementType.current_src)
    vgs.set_instance("VGS")

    ids_ref = ckt.add_edge(d, s, ids)
    vgs_ref = ckt.add_edge(g, s, vgs)

    ids.set_vgs_ref(vgs_ref=vgs_ref)

def add_pmos(g, d, s, KP, vth, l_lambda, L, W, ckt):
    isd = vccs_l1_mosfet()
    isd.set_type(ElementType.current_src)
    isd.set_instance("IDS")
    isd.set_params(KP=KP, \
                   vth=vth, \
                   l_lambda=l_lambda, \
                   L=L, \
                   W=W)

    vsg = current_src()
    vsg.set_type(ElementType.current_src)
    vsg.set_instance("VGS")

    ids_ref = ckt.add_edge(s, d, isd)
    vgs_ref = ckt.add_edge(s, g, vsg)

    isd.set_vgs_ref(vgs_ref=vgs_ref)

def add_circuit(ckt, rser):

    Vs_val = 1.0
    C2_val = 1e-9
    L1_val = 1e-3
    C1_val = 2e-9
    C3_val = 3e-9
    R2_val = 1e3
    R3_val = 1e-3

    if (rser):
        v_node = 4
    else:
        v_node = 1

    #Add Vs
    Vs = voltage_src(ramp=True, ramp_ddt=1e6)
    Vs.set_type(ElementType.voltage_src)
    Vs.set_value(Vs_val)
    Vs.set_instance("Vs")
    ckt.add_edge(v_node, 0, Vs)

    # Add C2
    C2 = capacitor()
    C2.set_type(ElementType.capacitor)
    C2.set_value(C2_val)
    C2.set_instance("C2")
    ckt.add_edge(1, 3, C2)

    # Add L1
    L1 = inductor()
    L1.set_type(ElementType.inductor)
    L1.set_value(L1_val)
    L1.set_instance("L1")
    ckt.add_edge(1, 2, L1)

    # Add C1
    C1 = capacitor()
    C1.set_type(ElementType.capacitor)
    C1.set_value(C1_val)
    C1.set_instance("C1")
    ckt.add_edge(2, 0, C1)

    # Add C3
    C3 = capacitor()
    C3.set_type(ElementType.capacitor)
    C3.set_value(C3_val)
    C3.set_instance("C3")
    ckt.add_edge(2, 3, C3)

    # Add R2
    R2 = resistor()
    R2.set_type(ElementType.resistor)
    R2.set_value(R2_val)
    R2.set_instance("R2")
    ckt.add_edge(3, 0, R2)

    if (rser):
        R3 = resistor()
        R3.set_type(ElementType.resistor)
        R3.set_value(R3_val)
        R3.set_instance("R3")
        ckt.add_edge(v_node, 1, R3)

def circuit_eqn(x, sys, ckt, vs, t):

    # set Vs
    #vs_ref = ckt.get_ref_from_instance("Vs")
    #ckt.set_value(ref=vs_ref, value=vs)

    # current Im
    I = ckt.get_im(x=x, sys=sys, t=t)

    # voltage Vm
    V = ckt.get_vm(x=x, sys=sys, t=t)

    I = ckt.get_degen()

    # Copy matrices
    qf_num = ckt.qf.copy()
    bf_num = ckt.bf.copy()

    eqn_qf = np.dot(qf_num, I)
    eqn_bf = np.dot(bf_num, V)

    return eqn_qf + eqn_bf

def dypc_litmus0(t, sys, ckt):

    root_start = np.ones(ckt.num_edges)
    vs = 1

    root = optimize.root(circuit_eqn, root_start, args=(sys, ckt, vs, t), tol=1e-7)
    root_start = root.x
    if not root.success:
        print (root.x)
        sys.exit(0)

    return ckt.get_dy(x=root.x)

def ode_solve(ckt):

    num_sys_vars    = ckt.get_num_sys_vars()

    tstep = 10000

    t     = np.linspace(0, 50e-6, tstep)
    x0    = np.zeros (num_sys_vars)

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

    return y

def main():

    ckt = Circuit()

    rser = True

    add_circuit(ckt, rser)

    ckt.init_circuit()

    y = ode_solve(ckt)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    tstep = 10000
    t     = np.linspace(0, 50e-6, tstep)

    for sys_var in range(len(y[0,:])):

        inst = str(ckt.get_edge_info_from_sys_var_ref(sys_var).instance)

        if inst[0] == "C":
            l = "$V_{" + inst + "}$"
            ax1.plot(t, y[:,sys_var], label=l)
        else:
            l = "$i_{" + inst + "}$"
            ax2.plot(t, y[:,sys_var], label=l, linestyle="dotted")

    if rser == False:

        VC3   = np.zeros(len(t))

        for idx in range(len(t)):
            for sys_var in range(len(y[0,:])):
                ckt.scb.v[ckt.get_edge_info_from_sys_var_ref(sys_var).ref] = y[idx, sys_var]
            ckt.scb.v[0] = 1.0

            C3 = ckt.get_edge_info(ckt.get_ref_from_instance("C3"))
            C3.get_voltage(ckt.scb)
            VC3[idx] = ckt.scb.v[ckt.get_edge_info(ckt.get_ref_from_instance("C3")).ref]

        l = "$V_{C3}$"
        ax1.plot(t, VC3, label=l)

    ax1.legend()
    ax2.legend(loc="lower right")

    ax1.set_ylabel("$V (V)$")
    ax1.set_xlabel("$t (S)$")

    ax2.set_ylabel("$i (A)$")
    ax2.set_xlabel("$t (S)$")

    if rser == False:
        plt.savefig("litmus0_rser0.png", bbox_inches = 'tight')
    else:
        plt.savefig("litmus0_rser1.png", bbox_inches = 'tight')

    plt.show()

if __name__ == "__main__":

    main()
