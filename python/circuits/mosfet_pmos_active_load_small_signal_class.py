import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../networks")
from pmos_active_load import *

from scipy import optimize

import numpy as np

import itertools

import matplotlib.pyplot as plt

import math

def circuit_eqn(x, ckt, iref_p, iref_n):

    # set Iref
    iref_p_ref = ckt.get_ref_from_instance("iref_p")
    ckt.set_value(ref=iref_p_ref, value=(iref_p))

    # set Iref
    iref_n_ref = ckt.get_ref_from_instance("iref_n")
    ckt.set_value(ref=iref_n_ref, value=(iref_n))

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

def sweep_sources(ckt, iref_sweep, iref_base, tol):

    root_start = np.ones(ckt.num_edges)

    scb = Scoreboard(num_edges=ckt.num_edges, num_sys_vars=ckt.num_sys_vars, degen_mtrx=ckt.degen_mtrx)

    v_iref_p    = np.zeros([len(iref_sweep)])
    v_iref_n    = np.zeros([len(iref_sweep)])

    for idx, iref in enumerate(iref_sweep):

        iref_p     = iref_base + iref
        iref_n     = iref_base - iref
        #print (iref_p)
        root       = optimize.root(circuit_eqn, root_start, args=(ckt, iref_p, iref_n), tol=tol)
        root_start = root.x

        if not root.success:
            print (root.x)
            sys.exit()

        scb.x = root.x
        scb.t = 0

        iref_p_ref = ckt.get_ref_from_instance("iref_p")
        iref_edge = ckt.get_edge_info(iref_p_ref)
        iref_edge.get_voltage(scb=scb)
        v_iref_p[idx] = scb.v[iref_p_ref]

        iref_n_ref = ckt.get_ref_from_instance("iref_n")
        iref_edge = ckt.get_edge_info(iref_n_ref)
        iref_edge.get_voltage(scb=scb)
        v_iref_n[idx] = scb.v[iref_n_ref]

    return v_iref_p, v_iref_n


def main():

    # Large Signal
    # ------------
    ckt_lrg = Circuit()
    nw = pmos_active_load(ckt_lrg)
    ckt_lrg.init_circuit()

    # Sweep Parameters
    # ----------------
    iref_delta  = 0.1e-6
    iref_sweep  = np.linspace(-iref_delta, iref_delta, 1000)

    # Operating Point
    # ---------------
    iref_p     = 20e-6
    iref_n     = 20e-6
    root_start = np.ones(ckt_lrg.num_edges)
    root       = optimize.root(circuit_eqn, root_start, args=(ckt_lrg, iref_p, iref_n), tol=1e-9)
    op         = ckt_lrg.scb

    if not root.success:
        print ("Failed to calculate small signal circuit")

    # Small signal
    # ------------
    ckt_small = Circuit()
    nw.add_sml_ckt(ckt_sml=ckt_small, op=op)
    ckt_small.init_circuit()

    # Sweep Large Signal
    # ------------------
    v_iref_p_lrg, v_iref_n_lrg = sweep_sources(ckt_lrg, iref_sweep, 20e-6, 1e-9)

    # Sweep Small Signal
    # ------------------
    v_iref_p_sml, v_iref_n_sml = sweep_sources(ckt_small, iref_sweep, 0, 1e-3)

    fig, (ax1, ax2) = plt.subplots(2, 1)

    ax1.grid()
    ax1.plot(iref_sweep, v_iref_p_sml,       color='blue',   label="$V_{i_{REF+}}$")
    ax1.plot(iref_sweep, v_iref_n_sml,       color='green',  label="$V_{i_{REF-}}$")
    ax1.set_ylabel("$V (V)$")
    ax1.set_xlabel("$i_{REF} (A)$")
    ax1.legend()
    ax1.title.set_text("Small Signal Sweep")

    ax2.grid()
    ax2.plot(iref_sweep, v_iref_p_lrg,       color='blue',   label="$V_{i_{REF+}}$")
    ax2.plot(iref_sweep, v_iref_n_lrg,       color='green',  label="$V_{i_{REF-}}$")
    ax2.set_ylabel("$V (V)$")
    ax2.set_xlabel("$i_{REF} (A)$")
    ax2.legend()
    ax2.title.set_text("Large Signal Sweep")

    plt.show()

if __name__ == "__main__":

    main()
