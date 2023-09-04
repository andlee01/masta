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

def ro(params, vds, vgs):

    vdssat = vgs - params["vth"]

    idsat = (params["KP"] / 2) * (params["W"] / params["L"]) * vdssat**2

    ro = 1 / (params["lambda"] * idsat)

    return ro

# equation 9.22
def gm(params, Ids):
    beta = params["KP"] * (params["W"] / params["L"])

    gm = math.sqrt(2 * beta * Ids)

    return gm

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

def add_nmos_small_signal(g, d, s, gm, ro, ckt):
    i_gm = TwoPortElement()
    i_gm.set_type(ElementType.current_src)
    i_gm.set_instance("IGM")
    i_gm.set_value(gm)

    i_ro = TwoPortElement()
    i_ro.set_type(ElementType.resistor)
    i_ro.set_instance("Ro")
    i_ro.set_value(ro)

    vgs = TwoPortElement()
    vgs.set_type(ElementType.voltage_src)
    vgs.set_instance("VGS")

    ckt.add_edge(d, s, i_gm)
    ckt.add_edge(d, s, i_ro)
    ckt.add_edge(g, s, vgs)

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

    # Add nmos
    add_nmos(g=v1, d=n2, s=GND, ckt=ckt)

    # Add Vgs
    vgs = voltage_src()
    vgs.set_type(ElementType.voltage_src)
    vgs.set_instance("VGS")
    ckt.add_edge(v1, GND, vgs)

def add_circuit_vref(ckt):

    # Nodes
    GND = 0
    VCC = 1
    v1  = 2
    v2  = 3
    n2  = 4

    # Add Vs
    vs = voltage_src()
    vs.set_type(ElementType.voltage_src)
    vs.set_instance("VS")
    vs.set_value(value=5.0)
    ckt.add_edge(VCC, GND, vs)

    # Add pmos
    add_pmos(g=v2, d=n2, s=VCC, ckt=ckt)

    # Add nmos
    add_nmos(g=v1, d=n2, s=GND, ckt=ckt)

    # Add Vgs
    vgs = voltage_src()
    vgs.set_type(ElementType.voltage_src)
    vgs.set_instance("VGS")
    ckt.add_edge(v1, GND, vgs)

    # Add Vref
    vref = voltage_src()
    vref.set_type(ElementType.voltage_src)
    vref.set_instance("VREF")
    ckt.add_edge(v2, GND, vref)

def circuit_eqn(x, ckt, vs, vgs):

    # set Vgs
    ckt.set_value(ref=5, value=vgs)

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

def circuit_eqn_vref(x, ckt, vs, vgs, vref):

    # set Vgs
    ckt.set_value(ref=5, value=vgs)

    # set Vref
    ckt.set_value(ref=6, value=vref)

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

    # Plots
    # ------
    vgs_sweep  = np.arange(1.0, 1.2, 0.01)
    vds_sweep  = np.arange(0.0, 5.0, 0.01)

    i_ds_n1       = np.zeros([len(vgs_sweep), len(vds_sweep)])
    i_sd_p1       = np.zeros([                len(vds_sweep)])

    ids = vccs_l1_mosfet()
    ids.set_type(ElementType.current_src)
    ids.set_instance("IDS")
    ids.set_params(KP=nmos_params["KP"], \
                   vth=nmos_params["vth"], \
                   l_lambda=nmos_params["lambda"], \
                   L=nmos_params["L"], \
                   W=nmos_params["W"])

    ids.set_ref(0)
    ids.set_vgs_ref(1)

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

    ax = plt.subplot()

    for vgs_idx, vgs in enumerate(vgs_sweep):
        for vds_idx, vds in enumerate(vds_sweep):

            x = [vds, vgs]
            i_ds_n1[vgs_idx][vds_idx] = ids.get_current(x=x, t=0)

        l = "Vgs = {:.2f}".format(vgs)
        plt.plot(vds_sweep, i_ds_n1[vgs_idx], label=l)

    vsg = 5.0 - 3.84
    for vds_idx, vds in enumerate(vds_sweep):

        vsd = 5.0 - vds
        x = [vsd, vsg]
        i_sd_p1[vds_idx] = isd.get_current(x=x, t=0)

    l = "isd"
    plt.plot(vds_sweep, i_sd_p1, label=l)

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show(block=False)

    # Diode Connected
    # ---------------
    ckt = Circuit()
    add_circuit(ckt)

    ckt.init_circuit()

    root_start = np.ones(6)
    vs         = 5.0
    vgs_sweep  = np.arange(1.0, 1.2, 0.0001)
    i_vs       = np.zeros([len(vgs_sweep)])
    v_vds_n1   = np.zeros([len(vgs_sweep)])

    for idx, vgs in enumerate(vgs_sweep):
        root       = optimize.root(circuit_eqn, root_start, args=(ckt, vs, vgs), tol=1e-9)
        root_start = root.x

        if not root.success:
            sys.exit()

        i_vs[idx]     = root.x[0]
        v_vds_n1[idx] = root.x[3]

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

    plt.show(block=False)

    # Vref Load
    # ---------
    ckt_vref = Circuit()
    add_circuit_vref(ckt_vref)

    ckt_vref.init_circuit()

    root_start = np.ones(7)
    vref       = 3.84
    i_vs       = np.zeros([len(vgs_sweep)])
    i_ids_n1   = np.zeros([len(vgs_sweep)])
    i_isd_p1   = np.zeros([len(vgs_sweep)])
    v_vds_n1   = np.zeros([len(vgs_sweep)])

    r_n1       = np.zeros([len(vgs_sweep)])
    r_p1       = np.zeros([len(vgs_sweep)])

    for idx, vgs in enumerate(vgs_sweep):
        root       = optimize.root(circuit_eqn_vref, root_start, args=(ckt_vref, vs, vgs, vref), tol=1e-8)
        root_start = root.x

        if not root.success:
            sys.exit()

        v_vds_n1[idx] = root.x[3]
        i_vs[idx]     = root.x[0]

        isd = ckt_vref.get_edge_info(1)
        ids = ckt_vref.get_edge_info(3)

        i_ids_n1[idx] = ids.get_current(x=root.x, t=0)
        i_isd_p1[idx] = isd.get_current(x=root.x, t=0)

        r_n1[idx] = ids.get_region(x=root.x)
        r_p1[idx] = isd.get_region(x=root.x)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    plt.grid()

    ax1.plot(vgs_sweep, v_vds_n1,       color='blue',    label="V_DS_N1")

    # Add an invisible plot and use to fill region
    ax1.plot(vgs_sweep, r_n1, alpha=0.0)
    ax1.fill_between(vgs_sweep, 0, 1, where=r_n1 > 0.5,
                     color='green', alpha=0.2, transform=ax1.get_xaxis_transform())
    ax1.text(1.1, 2.5, r"$n_1$ in triode", fontsize=10)

    ax1.plot(vgs_sweep, r_p1, alpha=0.0)
    ax1.fill_between(vgs_sweep, 0, 1, where=r_p1 > 0.5,
                     color='blue', alpha=0.2, transform=ax1.get_xaxis_transform())
    ax1.text(1.0, 2.5, r"$p_1$ in triode", fontsize=10)


    ax1.set_ylabel("$V_{DSM1} (V)$")
    ax1.set_xlabel("$V_{GSM1} (V)$")
    ax1.legend(loc="center right")

    ax2.plot(vgs_sweep, i_vs,     color='red',   label="$i_{vs}$")
    ax2.plot(vgs_sweep, i_ids_n1, color='blue',  label="$i_{n1}$")
    ax2.plot(vgs_sweep, i_isd_p1, color='green', label="$i_{p1}$")
    ax2.set_ylabel("$i (A)$")

    ax2.legend()

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    plt.show()

if __name__ == "__main__":

    main()
