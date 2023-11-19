import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

from scipy import optimize

import numpy as np

import itertools

import matplotlib.pyplot as plt

import math

import control as ct

from scipy import optimize
from scipy.integrate import ode

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

def add_nmos_small(g, d, s, KP, vth, l_lambda, L, W, ckt, op):
    ids = vccs_l1_mosfet_small()
    ids.set_instance("IDS")
    ids.set_params(KP=KP, \
                   vth=vth, \
                   l_lambda=l_lambda, \
                   L=L, \
                   W=W)

    i_ro = resistor()
    i_ro.set_type(ElementType.resistor)
    i_ro.set_instance("Ro")

    vgs = current_src()
    vgs.set_type(ElementType.current_src)
    vgs.set_is_const()
    vgs.set_value(0.0)
    vgs.set_instance("VGS")

    ids_ref = ckt.add_edge(d, s, ids)
    vgs_ref = ckt.add_edge(g, s, vgs)
    ckt.add_edge(d, s, i_ro)

    ids.set_vgs_ref(vgs_ref=vgs_ref)
    ids.set_src_dep_ref(ref=vgs_ref)

    [gm, ro] = ids.get_op(op=op)
    i_ro.set_value(ro)

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

def add_pmos_small(g, d, s, KP, vth, l_lambda, L, W, ckt, op):
    isd = vccs_l1_mosfet_small()
    isd.set_type(ElementType.current_src)
    isd.set_instance("IDS")
    isd.set_params(KP=KP, \
                   vth=vth, \
                   l_lambda=l_lambda, \
                   L=L, \
                   W=W)

    i_ro = resistor()
    i_ro.set_type(ElementType.resistor)
    i_ro.set_instance("Ro")

    vsg = current_src()
    vsg.set_type(ElementType.current_src)
    vsg.set_is_const()
    vsg.set_value(0.0)
    vsg.set_instance("VGS")

    ids_ref = ckt.add_edge(s, d, isd)
    vgs_ref = ckt.add_edge(s, g, vsg)
    ckt.add_edge(s, d, i_ro)

    isd.set_vgs_ref(vgs_ref=vgs_ref)
    isd.set_src_dep_ref(ref=vgs_ref)

    [gm, ro] = isd.get_op(op=op)
    i_ro.set_value(ro)

def add_circuit(ckt):

    # Nodes
    GND = 0
    VCC = 1
    n1  = 2
    n2  = 3
    n3  = 4
    n4  = 5
    n5  = 6
    v1  = 7
    v2  = 8
    v3  = 9
    v4  = 10
    nvout = 11

    # Add Vs
    vs = voltage_src()
    vs.set_type(ElementType.voltage_src)
    vs.set_instance("VS")
    vs.set_value(value=5.0)
    ckt.add_edge(VCC, GND, vs)

    # Add pmos M3
    add_pmos(g=n1, d=n1, s=VCC, KP=40e-6, vth=0.9, l_lambda=0.0125, L=2, W=30, ckt=ckt)

    # Add pmos M4
    add_pmos(g=n1, d=nvout, s=VCC, KP=40e-6, vth=0.9, l_lambda=0.0125, L=2, W=30, ckt=ckt)

    # Add nmos M1
    add_nmos(g=v1, d=n1, s=n4, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=10, ckt=ckt)

    # Add nmos M2
    add_nmos(g=v2, d=nvout, s=n5, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=10, ckt=ckt)

    # Add M1 Tail Current
    m1_bias = current_src()
    m1_bias.set_type(ElementType.current_src)
    m1_bias.set_instance("m1_bias")
    m1_bias.set_value(value=10e-6)
    ckt.add_edge(n4, GND, m1_bias)

    # Add M2 Tail Current
    m2_bias = current_src()
    m2_bias.set_type(ElementType.current_src)
    m2_bias.set_instance("m1_bias")
    m2_bias.set_value(value=10e-6)
    ckt.add_edge(n5, GND, m2_bias)

    # Add source degen resistor
    r_src_dgen = resistor()
    r_src_dgen.set_instance("R_SRC_DEGEN")
    r_src_dgen.set_value(value=1e3)
    ckt.add_edge(n5, n4, r_src_dgen)

    # Add Vgs V1
    vgs_1 = voltage_src()
    vgs_1.set_type(ElementType.current_src)
    vgs_1.set_instance("VGS_1")
    ckt.add_edge(v1, GND, vgs_1)

    # Add Vgs V2
    vgs_2 = voltage_src()
    vgs_2.set_type(ElementType.voltage_src)
    vgs_2.set_instance("VGS_2")
    vgs_2.set_value(value=2.0)
    ckt.add_edge(v2, GND, vgs_2)

def add_circuit_large_dyn(ckt):

    # Nodes
    GND = 0
    VCC = 1
    n1  = 2
    n2  = 3
    n3  = 4
    n4  = 5
    n5  = 6
    v1  = 7
    v2  = 8
    v3  = 9
    v4  = 10
    nvout = 11

    # Add Vs
    vs = voltage_src()
    vs.set_type(ElementType.voltage_src)
    vs.set_instance("VS")
    vs.set_value(value=5.0)
    ckt.add_edge(VCC, GND, vs)

    # Add pmos M3
    add_pmos(g=n1, d=n1, s=VCC, KP=40e-6, vth=0.9, l_lambda=0.0125, L=2, W=30, ckt=ckt)

    # Add pmos M4
    add_pmos(g=n1, d=nvout, s=VCC, KP=40e-6, vth=0.9, l_lambda=0.0125, L=2, W=30, ckt=ckt)

    # Add nmos M1
    add_nmos(g=v1, d=n1, s=n4, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=10, ckt=ckt)

    # Add nmos M2
    add_nmos(g=v2, d=nvout, s=n5, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=10, ckt=ckt)

    # Add M1 Tail Current
    m1_bias = current_src()
    m1_bias.set_type(ElementType.current_src)
    m1_bias.set_instance("m1_bias")
    m1_bias.set_value(value=10e-6)
    ckt.add_edge(n4, GND, m1_bias)

    # Add M2 Tail Current
    m2_bias = current_src()
    m2_bias.set_type(ElementType.current_src)
    m2_bias.set_instance("m1_bias")
    m2_bias.set_value(value=10e-6)
    ckt.add_edge(n5, GND, m2_bias)

    # Add source degen resistor
    r_src_dgen = resistor()
    r_src_dgen.set_instance("R_SRC_DEGEN")
    r_src_dgen.set_value(value=1e3)
    ckt.add_edge(n5, n4, r_src_dgen)

    c_src_dgen = capacitor()
    c_src_dgen.set_instance("C_SRC_DEGEN")
    c_src_dgen.set_value(value=1e-9)
    ckt.add_edge(n5, n4, c_src_dgen)

    # Add Vgs V1
    vgs_1 = voltage_src()
    vgs_1.set_type(ElementType.current_src)
    vgs_1.set_instance("VGS_1")
    ckt.add_edge(v1, GND, vgs_1)

    # Add Vgs V2
    vgs_2 = voltage_src()
    vgs_2.set_type(ElementType.voltage_src)
    vgs_2.set_instance("VGS_2")
    vgs_2.set_value(value=2.0)
    ckt.add_edge(v2, GND, vgs_2)

def add_circuit_small(ckt, op):

    # Nodes
    GND = 0
    VCC = 1
    n1  = 2
    n2  = 3
    n3  = 4
    n4  = 5
    n5  = 6
    v1  = 7
    v2  = 8
    v3  = 9
    v4  = 10
    nvout = 11

    # Add pmos M3
    add_pmos_small(g=n1, d=n1, s=GND, KP=40e-6, vth=0.9, l_lambda=0.0125, L=2, W=30, ckt=ckt, op=op)

    # Add pmos M4
    add_pmos_small(g=n1, d=nvout, s=GND, KP=40e-6, vth=0.9, l_lambda=0.0125, L=2, W=30, ckt=ckt, op=op)

    # Add nmos M1
    add_nmos_small(g=v1, d=n1, s=n4, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=10, ckt=ckt, op=op)

    # Add nmos M2
    add_nmos_small(g=v2, d=nvout, s=n5, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=10, ckt=ckt, op=op)

    # Add source degen resistor (12)
    r_src_dgen = resistor()
    r_src_dgen.set_instance("R_SRC_DEGEN")
    r_src_dgen.set_value(value=1e3)
    ckt.add_edge(n5, n4, r_src_dgen)

    # (13)
    c_src_dgen = capacitor()
    c_src_dgen.set_instance("C_SRC_DEGEN")
    c_src_dgen.set_value(value=1e-9)
    ckt.add_edge(n5, n4, c_src_dgen)

    # Add Vgs V1 (14)
    vgs_1 = voltage_src()
    vgs_1.set_type(ElementType.voltage_src)
    vgs_1.set_instance("VGS_1")
    vgs_1.set_value(value=0.0)
    vgs_1.set_is_input()
    ckt.add_edge(v1, GND, vgs_1)

    # Add Vgs V2 (15)
    vgs_2 = voltage_src()
    vgs_2.set_type(ElementType.voltage_src)
    vgs_2.set_instance("VGS_2")
    vgs_2.set_is_const()
    vgs_2.set_value(value=0.0)
    ckt.add_edge(v2, GND, vgs_2)

    vmeas = current_src()
    vmeas.set_is_const()
    vmeas.set_value(0.0)
    vmeas.set_instance("VMEAS")
    ckt.add_edge(nvout, GND, vmeas)

def circuit_eqn(x, ckt, vgs_1):

    # set Vgs_1
    vgs_1_ref = ckt.get_ref_from_instance("VGS_1")
    ckt.set_value(ref=vgs_1_ref, value=vgs_1)

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

def circuit_eqn_dyn(x, sys, ckt, vgs_1):

    # set Vgs_1
    vgs_1_ref = ckt.get_ref_from_instance("VGS_1")
    ckt.set_value(ref=vgs_1_ref, value=vgs_1)

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

def dypc_litmus0(t, sys, ckt):

    root_start = np.ones(ckt.num_edges)
    vs = 1

    root = optimize.root(circuit_eqn_dyn, root_start, args=(sys, ckt, vs), tol=1e-6)
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

    # Large signal
    ckt = Circuit()
    add_circuit(ckt)
    ckt.init_circuit()

    r_src_dgen_sweep  = np.arange(1e3, 10e3, 1e3)

    v_cm   = 2.0
    v_diff = 0.1

    # set Vgs_2
    vgs_2_ref = ckt.get_ref_from_instance("VGS_2")
    ckt.set_value(ref=vgs_2_ref, value=v_cm)

    # set Vgs_1
    vgs_1_ref = ckt.get_ref_from_instance("VGS_1")
    ckt.set_value(ref=vgs_1_ref, value=v_cm)

    # Operating Point
    root_start = np.ones(ckt.num_edges)
    root       = optimize.root(circuit_eqn, root_start, args=(ckt, v_cm), tol=1e-9)
    op         = ckt.scb

    if not root.success:
        print ("Failed to calculate small signal circuit")

    # Small signal
    ckt_small = Circuit(lti=True)
    add_circuit_small(ckt=ckt_small, op=op)
    ckt_small.init_circuit()

    print (ckt_small.qf)

    ckt_small.get_ss()

    x0 = np.zeros(ckt_small.num_sys_vars)
    x0 = ckt_small.get_lti_const_x0(x0)

    sys = ct.ss(ckt_small.A, ckt_small.B, ckt_small.C, ckt_small.D)

    print (ckt_small.B)

    tstep = 10000
    t     = np.linspace(0, 50e-6, tstep)

    T, yout = ct.step_response(sys, output=0, T=t, X0=x0)

    outputs = np.zeros(len(yout[0][0]))
    for i in range(len(outputs)):
        outputs[i] = yout[0][0][i]

    plt.plot(T, outputs)
    plt.show(block=False)

    # Large signal
    ckt_dyn = Circuit()
    add_circuit_small(ckt_dyn, op)
    ckt_dyn.init_circuit()

    y = ode_solve(ckt_dyn)

    plt.plot(t,y)
    plt.show()

    print (ckt_small.num_sys_vars)


if __name__ == "__main__":

    main()
