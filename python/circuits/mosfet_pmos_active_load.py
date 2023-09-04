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

def add_nmos(g, d, s, ckt):
    ids = vccs_l1_mosfet()
    ids.set_type(ElementType.current_src)
    ids.set_instance("IDS")
    ids.set_params(KP=nmos_params["KP"], \
                   vth=nmos_params["vth"], \
                   l_lambda=nmos_params["lambda"], \
                   L=nmos_params["L"], \
                   W=nmos_params["W"])

    vgs = current_src()
    vgs.set_type(ElementType.current_src)
    vgs.set_instance("VGS")

    ids_ref = ckt.add_edge(d, s, ids)
    vgs_ref = ckt.add_edge(g, s, vgs)

    ids.set_vgs_ref(vgs_ref=vgs_ref)

def add_pmos(g, d, s, ckt):
    isd = vccs_l1_mosfet()
    isd.set_type(ElementType.current_src)
    isd.set_instance("IDS")
    isd.set_params(KP=pmos_params["KP"], \
                   vth=pmos_params["vth"], \
                   l_lambda=pmos_params["lambda"], \
                   L=pmos_params["L"], \
                   W=pmos_params["W"])

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
    v1  = 3
    n2  = 4

    # Add Vs
    vs = voltage_src()
    vs.set_type(ElementType.voltage_src)
    vs.set_instance("VS")
    vs.set_value(value=5.0)
    ckt.add_edge(VCC, GND, vs)

    # Add pmos
    add_pmos(g=n2, d=n2, s=VCC, ckt=ckt)

    # Add Iref
    iref_p = current_src()
    iref_p.set_type(ElementType.current_src)
    iref_p.set_instance("iref_p")
    iref_p.set_value(value=20e-6)
    ckt.add_edge(n2, GND, iref_p)

    # Add pmos
    add_pmos(g=n2, d=n1, s=VCC, ckt=ckt)

    # Add Iref
    iref_n = current_src()
    iref_n.set_type(ElementType.current_src)
    iref_n.set_instance("iref_n")
    iref_n.set_value(value=20e-6)
    ckt.add_edge(n1, GND, iref_n)

def circuit_eqn(x, ckt, iref):

    # set Iref
    iref_p_ref = ckt.get_ref_from_instance("iref_p")
    ckt.set_value(ref=iref_p_ref, value=(20e-6 + iref))

    # set Iref
    iref_n_ref = ckt.get_ref_from_instance("iref_n")
    ckt.set_value(ref=iref_n_ref, value=(20e-6 - iref))

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

    root_start = np.ones(7)
    vs          = 5.0
    iref_delta  = 0.1e-6
    iref_sweep  = np.linspace(-iref_delta, iref_delta, 1000)
    v_iref_p    = np.zeros([len(iref_sweep)])
    v_iref_n    = np.zeros([len(iref_sweep)])
    iref_limits = [20e-6 + iref_delta, 20e-6 - iref_delta]

    isg_pmos_vref = np.zeros(2)
    isd_pmos_vref = np.zeros(2)

    for idx, iref in enumerate(iref_sweep):
        root       = optimize.root(circuit_eqn, root_start, args=(ckt, iref), tol=1e-9)
        root_start = root.x

        if not root.success:
            sys.exit()

        iref_p_ref = ckt.get_ref_from_instance("iref_p")
        iref_edge = ckt.get_edge_info(iref_p_ref)
        v_iref_p[idx] = iref_edge.get_voltage(x=root.x, t=0)

        iref_n_ref = ckt.get_ref_from_instance("iref_n")
        iref_edge = ckt.get_edge_info(iref_n_ref)
        v_iref_n[idx] = iref_edge.get_voltage(x=root.x, t=0)

        # Limits
        if idx == 0:
            isg_pmos_vref[0] = ckt.get_edge_info(5).get_voltage(x=root.x, t=0)
            isd_pmos_vref[0] = ckt.get_edge_info(4).get_voltage(x=root.x, t=0)
        elif idx == len(iref_sweep)-1:
            isg_pmos_vref[1] = ckt.get_edge_info(5).get_voltage(x=root.x, t=0)
            isd_pmos_vref[1] = ckt.get_edge_info(4).get_voltage(x=root.x, t=0)

    fig, ax1 = plt.subplots()

    plt.grid()

    ax1.plot(iref_sweep, v_iref_p,       color='blue',   label="$V_{i_{REF+}}$")
    ax1.plot(iref_sweep, v_iref_n,       color='green',  label="$V_{i_{REF-}}$")
    ax1.set_ylabel("$V (V)$")
    ax1.set_xlabel("$i_{REF} (A)$")
    ax1.legend()

    plt.savefig("../../doc/mosfet_pmos_active_load_viref.svg", bbox_inches = 'tight')

    # Limits
    # ------
    isd = vccs_l1_mosfet()
    isd.set_type(ElementType.current_src)
    isd.set_instance("IDS")
    isd.set_params(KP=pmos_params["KP"], \
                   vth=pmos_params["vth"], \
                   l_lambda=pmos_params["lambda"], \
                   L=pmos_params["L"], \
                   W=pmos_params["W"])

    isd.set_ref(0)
    isd.set_vgs_ref(1)

    vds_sweep  = np.linspace(0.2, 4.8, 1000)
    i_sd_p1    = np.zeros([len(isg_pmos_vref), len(vds_sweep)])

    fig, ax1 = plt.subplots()
    plt.grid()

    for vsg_idx, vsg in enumerate(isg_pmos_vref):
        for vds_idx, vds in enumerate(vds_sweep):

            x = [vds, vsg]
            i_sd_p1[vsg_idx][vds_idx] = isd.get_current(x=x, t=0)

        if vsg_idx == 0:
            color = "red"
        else:
            color = "green"

        l_vgs = "{:.5f}".format(vsg)
        ax1.plot(vds_sweep, i_sd_p1[vsg_idx], color=color, label="$V_{GS}M_2 = $" + l_vgs)
        ax1.axhline(y=iref_limits[vsg_idx],   color=color, linestyle="--")
        ax1.axvline(x=isd_pmos_vref[vsg_idx], color=color, linestyle="dotted", label="$V_{SD}M_2$")
        ax1.set_ylabel('$i_{REF-} (A)$')
        ax1.set_xlabel('$V_{SD}M_2 (V)$')
        ax1.legend()

    plt.savefig("../../doc/mosfet_pmos_active_load_iref.svg", bbox_inches = 'tight')

if __name__ == "__main__":

    main()
