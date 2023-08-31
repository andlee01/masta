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

    # Add Vs
    vs = voltage_src()
    vs.set_type(ElementType.voltage_src)
    vs.set_instance("VS")
    vs.set_value(value=5.0)
    ckt.add_edge(VCC, GND, vs)

    # Add nmos M1
    add_nmos(g=n1, d=n1, s=GND, KP=120e-6, vth=0.8, l_lambda=0.01, L=2, W=10, ckt=ckt)

    # Add Iref
    iref = current_src()
    iref.set_type(ElementType.current_src)
    iref.set_instance("Iref")
    iref.set_value(20e-6)
    ckt.add_edge(VCC, n1, iref)

def circuit_eqn(x, ckt, Iref):

    # set Vo
    iref_ref = ckt.get_ref_from_instance("Iref")
    ckt.set_value(ref=iref_ref, value=Iref)

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

    iref_sweep   = np.arange(1e-6, 150e-6, 0.1e-6)
    v_vds_m1     = np.zeros([len(iref_sweep)])
    i_ids_m1     = np.zeros([len(iref_sweep)])

    root_start = np.ones(ckt.num_edges)

    for idx, iref in enumerate(iref_sweep):
        root       = optimize.root(circuit_eqn, root_start, args=(ckt, iref), tol=1e-9)
        root_start = root.x
        if not root.success:
            sys.exit()

        ids_m1 = ckt.get_edge_info(1)

        v_vds_m1[idx] = ids_m1.get_voltage(x=root.x, t=0)
        i_ids_m1[idx] = ids_m1.get_current(x=root.x, t=0)


    vgs_sweep  = np.arange(0.8, 1.5, 0.05)
    vds_sweep  = np.arange(0.0, 2.5, 0.01)

    i_ds_n1       = np.zeros([len(vgs_sweep), len(vds_sweep)])

    ids = vccs_l1_mosfet()
    ids.set_type(ElementType.current_src)
    ids.set_instance("IDS")
    ids.set_params(KP=120e-6, \
                   vth=0.8, \
                   l_lambda=0.01, \
                   L=2, \
                   W=10)

    ids.set_ref(0)
    ids.set_vgs_ref(1)

    fig, ax1 = plt.subplots()
    ax1.grid()

    for vgs_idx, vgs in enumerate(vgs_sweep):
        for vds_idx, vds in enumerate(vds_sweep):

            x = [vds, vgs]
            i_ds_n1[vgs_idx][vds_idx] = ids.get_current(x=x, t=0)

        l = "Vgs = {:.2f}".format(vgs)
        plt.plot(vds_sweep, i_ds_n1[vgs_idx], label=l)

    ax1.plot(v_vds_m1, i_ids_m1,       color='blue',  label="i_DS_M1")
    ax1.set_ylabel("Ids (A)")
    ax1.set_xlabel("Vds (V)")
    ax1.legend(loc="right")

    plt.savefig("../../doc/plt_mosfet_bias_ex1.svg", bbox_inches = 'tight')

    #plt.show()

if __name__ == "__main__":

    main()
