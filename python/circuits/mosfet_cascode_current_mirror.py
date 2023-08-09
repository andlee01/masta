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
    n4  = 5
    vout  = 7

    # Add Vs
    vs = voltage_src()
    vs.set_type(ElementType.voltage_src)
    vs.set_instance("VS")
    vs.set_value(value=5.0)
    ckt.add_edge(VCC, GND, vs)

    # Add nmos M4
    add_nmos(g=n1, d=vout, s=n2, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=10, ckt=ckt)

    # Add nmos M2
    add_nmos(g=n4, d=n2, s=GND, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=10, ckt=ckt)

    # Add nmos M3
    add_nmos(g=n1, d=n1, s=n4, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=10, ckt=ckt)

    # Add nmos M1
    add_nmos(g=n4, d=n4, s=GND, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=10, ckt=ckt)

    # Add Vo
    vo = voltage_src()
    vo.set_type(ElementType.voltage_src)
    vo.set_instance("Vo")
    ckt.add_edge(vout, GND, vo)

    # Add Iref
    iref = current_src()
    iref.set_type(ElementType.current_src)
    iref.set_instance("Iref")
    iref.set_value(20e-6)
    ckt.add_edge(VCC, n1, iref)

def circuit_eqn(x, ckt, vgs_1):

    # set Vo
    Vo_ref = ckt.get_ref_from_instance("Vo")
    ckt.set_value(ref=Vo_ref, value=vgs_1)

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

    ckt = Circuit()
    add_circuit(ckt)

    ckt.init_circuit()

    root_start = np.ones(ckt.num_edges)
    vs         = 5.0
    vo_sweep   = np.arange(0.15, 5.0, 0.01)
    i_vs       = np.zeros([len(vo_sweep)])
    v_vds_n2   = np.zeros([len(vo_sweep)])

    v_vgs_n1   = np.zeros([len(vo_sweep)])
    v_vgs_n2   = np.zeros([len(vo_sweep)])
    v_vgs_n3   = np.zeros([len(vo_sweep)])
    v_vgs_n4   = np.zeros([len(vo_sweep)])

    for idx, vo in enumerate(vo_sweep):
        root       = optimize.root(circuit_eqn, root_start, args=(ckt, vo), tol=1e-9)
        root_start = root.x
        if not root.success:
            sys.exit()

        ids_n1 = ckt.get_edge_info(7)
        igs_n1 = ckt.get_edge_info(8)
        ids_n2 = ckt.get_edge_info(3)
        igs_n2 = ckt.get_edge_info(4)
        ids_n3 = ckt.get_edge_info(5)
        igs_n3 = ckt.get_edge_info(6)
        ids_n4 = ckt.get_edge_info(1)
        igs_n4 = ckt.get_edge_info(2)

        i_vs[idx]     = root.x[0]
        v_vds_n2[idx] = ids_n2.get_voltage(x=root.x, t=0)

        v_vgs_n1[idx] = igs_n1.get_voltage(x=root.x, t=0)
        v_vgs_n2[idx] = igs_n2.get_voltage(x=root.x, t=0)
        v_vgs_n3[idx] = igs_n3.get_voltage(x=root.x, t=0)
        v_vgs_n4[idx] = igs_n4.get_voltage(x=root.x, t=0)


    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax1.grid()

    ax1.plot(vo_sweep, v_vds_n2,       color='blue',  label="V_DS_N2")
    ax1.set_ylabel('v')
    ax1.legend()

    ax2.plot(vo_sweep, i_vs,    color='red', label="i_vs")
    ax2.set_ylabel('i')

    ax2.legend()

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    plt.show(block=False)
    plt.grid()

    fig, ax1 = plt.subplots()
    ax1.plot(vo_sweep, v_vgs_n1,       color='red',   label="V_GS_N1")
    ax1.plot(vo_sweep, v_vgs_n2,       color='green', label="V_GS_N2")
    ax1.plot(vo_sweep, v_vgs_n3,       color='orange',label="V_GS_N3")
    ax1.plot(vo_sweep, v_vgs_n4,       color='purple',label="V_GS_N4")

    ax1.grid()
    ax1.legend()

    plt.show()

if __name__ == "__main__":

    main()
