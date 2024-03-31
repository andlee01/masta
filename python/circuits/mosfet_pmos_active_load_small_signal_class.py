import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../subckt")
from mos import *

sys.path.append("../networks")
from pmos_active_load import *

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

def add_circuit(ckt):

    # Nodes
    GND = 0
    VCC = 1
    n1  = 2
    v1  = 3
    n2  = 4

    # Add pmos
    add_pmos(g=n2, d=n2, s=VCC, \
             KP=pmos_params["KP"], \
             vth=pmos_params["vth"], \
             l_lambda=pmos_params["lambda"], \
             L=pmos_params["L"], \
             W=pmos_params["W"], \
             ckt=ckt)

    # Add Iref
    iref_p = current_src()
    iref_p.set_type(ElementType.current_src)
    iref_p.set_instance("iref_p")
    iref_p.set_value(value=20e-6)
    ckt.add_edge(n2, GND, iref_p)

    # Add pmos
    add_pmos(g=n2, d=n1, s=VCC, \
             KP=pmos_params["KP"], \
             vth=pmos_params["vth"], \
             l_lambda=pmos_params["lambda"], \
             L=pmos_params["L"], \
             W=pmos_params["W"], \
             ckt=ckt)

    # Add Iref
    iref_n = current_src()
    iref_n.set_type(ElementType.current_src)
    iref_n.set_instance("iref_n")
    iref_n.set_value(value=20e-6)
    ckt.add_edge(n1, GND, iref_n)

    # Add Vs
    vs = voltage_src()
    vs.set_type(ElementType.voltage_src)
    vs.set_instance("VS")
    vs.set_value(value=5.0)
    ckt.add_edge(VCC, GND, vs)

def add_circuit_small(ckt, op):

    # Nodes
    GND = 0
    VCC = 1
    n1  = 2
    v1  = 3
    n2  = 4

    # Add pmos
    add_pmos_small(g=n2, d=n2, s=GND, \
                   KP=pmos_params["KP"], \
                   vth=pmos_params["vth"], \
                   l_lambda=pmos_params["lambda"], \
                   L=pmos_params["L"], \
                   W=pmos_params["W"], \
                   ckt=ckt, \
                   op=op, \
                   strt=0)

    # Add Iref
    iref_p = current_src()
    iref_p.set_type(ElementType.current_src)
    iref_p.set_instance("iref_p")
    iref_p.set_value(value=0)
    ckt.add_edge(n2, GND, iref_p)

    # Add pmos
    add_pmos_small(g=n2, d=n1, s=GND, \
                   KP=pmos_params["KP"], \
                   vth=pmos_params["vth"], \
                   l_lambda=pmos_params["lambda"], \
                   L=pmos_params["L"], \
                   W=pmos_params["W"], \
                   ckt=ckt, \
                   op=op, \
                   strt=3)

    # Add Iref
    iref_n = current_src()
    iref_n.set_type(ElementType.current_src)
    iref_n.set_instance("iref_n")
    iref_n.set_value(value=0)
    ckt.add_edge(n1, GND, iref_n)

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

def main():

    # Diode Connected
    # ---------------
    ckt_lrg = Circuit()
    nw = pmos_active_load(ckt_lrg)

    #add_circuit(ckt_lrg)

    ckt_lrg.init_circuit()

    # Operating Point
    iref_p     = 20e-6
    iref_n     = 20e-6
    root_start = np.ones(ckt_lrg.num_edges)
    root       = optimize.root(circuit_eqn, root_start, args=(ckt_lrg, iref_p, iref_n), tol=1e-9)
    op         = ckt_lrg.scb

    if not root.success:
        print ("Failed to calculate small signal circuit")

    #sys.exit(1)

    # Small signal
    ckt_small = Circuit()
    #add_circuit_small(ckt=ckt_small, op=op)
    nw.add_sml_ckt(ckt_sml=ckt_small, op=op)
    ckt_small.init_circuit()

    # Sweep
    # -----
    iref_delta  = 0.1e-6
    iref_sweep  = np.linspace(-iref_delta, iref_delta, 1000)
    v_iref_p    = np.zeros([len(iref_sweep)])
    v_iref_n    = np.zeros([len(iref_sweep)])

    root_start = np.ones(ckt_small.num_edges)

    scb = Scoreboard(num_edges=ckt_small.num_edges, num_sys_vars=ckt_small.num_sys_vars, degen_mtrx=ckt_small.degen_mtrx)

    for idx, iref in enumerate(iref_sweep):

        iref_p     =  iref
        iref_n     = -iref

        root       = optimize.root(circuit_eqn, root_start, args=(ckt_small, iref_p, iref_n), tol=1e-9)
        #root_start = root.x

        if not root.success:
            print (root.x)
            sys.exit()

        scb.x = root.x
        scb.t = 0

        iref_p_ref = ckt_small.get_ref_from_instance("iref_p")
        iref_edge = ckt_small.get_edge_info(iref_p_ref)
        iref_edge.get_voltage(scb=scb)
        v_iref_p[idx] = scb.v[iref_p_ref]

        iref_n_ref = ckt_small.get_ref_from_instance("iref_n")
        iref_edge = ckt_small.get_edge_info(iref_n_ref)
        iref_edge.get_voltage(scb=scb)
        v_iref_n[idx] = scb.v[iref_n_ref]

    fig, ax1 = plt.subplots()

    plt.grid()

    ax1.plot(iref_sweep, v_iref_p,       color='blue',   label="$V_{i_{REF+}}$")
    ax1.plot(iref_sweep, v_iref_n,       color='green',  label="$V_{i_{REF-}}$")
    ax1.set_ylabel("$V (V)$")
    ax1.set_xlabel("$i_{REF} (A)$")
    ax1.legend()

    plt.show()

if __name__ == "__main__":

    main()
