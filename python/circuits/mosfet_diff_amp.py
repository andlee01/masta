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
    add_nmos(g=v2, d=nvout, s=n4, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=10, ckt=ckt)

    # Add nmos M6T
    add_nmos(g=v3, d=n4, s=n5, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=20, ckt=ckt)

    # Add nmos M6B
    add_nmos(g=v4, d=n5, s=GND, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=20, ckt=ckt)

    # Add Vgs V1
    vgs_1_ref = ckt.num_edges
    vgs_1 = voltage_src()
    vgs_1.set_type(ElementType.voltage_src)
    vgs_1.set_instance("VGS_1")
    ckt.add_edge(v1, GND, vgs_1)

    # Add Vgs V2
    vgs_2 = voltage_src()
    vgs_2.set_type(ElementType.voltage_src)
    vgs_2.set_instance("VGS_2")
    vgs_2.set_value(value=4.0)
    ckt.add_edge(v2, GND, vgs_2)

    # Add Vbias3
    vb3 = voltage_src()
    vb3.set_type(ElementType.voltage_src)
    vb3.set_instance("Vbias3")
    vb3.set_value(value=1.49)
    ckt.add_edge(v3, GND, vb3)

    # Add Vbias4
    vb4 = voltage_src()
    vb4.set_type(ElementType.voltage_src)
    vb4.set_instance("Vbias4")
    vb4.set_value(value=1.086)
    ckt.add_edge(v4, GND, vb4)

def circuit_eqn(x, ckt, vgs_1):

    # set Vgs_1
    vgs_1_ref = ckt.get_ref_from_instance("VGS_1")
    ckt.set_value(ref=vgs_1_ref, value=vgs_1)

    # current Im
    I = ckt.get_im(x=x, t=0)

    # voltage Vm
    V = ckt.get_vm(x=x, t=0)

    # Copy matrices
    qf_num = ckt.qf.copy()
    bf_num = ckt.bf.copy()

    eqn_qf = np.dot(qf_num, I)
    eqn_bf = np.dot(bf_num, V)

    return eqn_qf + eqn_bf

def main():

    # Diode Connected
    # ---------------
    ckt = Circuit()
    add_circuit(ckt)

    ckt.init_circuit()

    root_start = np.ones(ckt.num_edges)
    vs         = 5.0
    vgs_sweep  = np.arange(3.9, 4.1, 0.001)
    i_vs       = np.zeros([len(vgs_sweep)])
    v_vds_n1   = np.zeros([len(vgs_sweep)])

    for idx, vgs in enumerate(vgs_sweep):
        root       = optimize.root(circuit_eqn, root_start, args=(ckt, vgs), tol=1e-9)
        root_start = root.x
        if not root.success:
            sys.exit()

        i_vs[idx]     = root.x[0]
        v_vds_n1[idx] = 5.0 - root.x[3]

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    plt.grid()

    ax1.plot(vgs_sweep, v_vds_n1,       color='blue',  label="V_DS_N1")
    ax1.set_ylabel('v')
    ax1.legend()

    ax2.plot(vgs_sweep, i_vs,    color='red', label="i_vs")
    ax2.set_ylabel('i')

    ax2.legend()

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    plt.show()

if __name__ == "__main__":

    main()
