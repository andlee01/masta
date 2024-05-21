import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../networks")
from gilbert_cell import *

from scipy import optimize
from scipy.integrate import ode

import numpy as np

import itertools

import matplotlib.pyplot as plt

import math

import control as ct

def circuit_eqn(x, ckt):

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

def sweep_lrg(iref_sweep, vlo_sweep, vlo_bias, nw, ckt_op):

    # Vector to hold results
    vf_out = np.zeros([len(vlo_sweep), len(iref_sweep)])

    # Define Vlo bias voltage
    nw.set_vlo_bias(bias=vlo_bias)

    scb = Scoreboard(num_edges=ckt_op.num_edges, num_sys_vars=ckt_op.num_sys_vars, degen_mtrx=ckt_op.degen_mtrx)

    for idx_vlo, vlo in enumerate(vlo_sweep):
        nw.set_vlo(vlo)

        root_start = np.ones(ckt_op.num_edges)
        for idx, iref in enumerate(iref_sweep):
            nw.set_iref(iref)

            # Operating Point
            root       = optimize.root(circuit_eqn, root_start, args=(ckt_op), tol=1e-9)

            scb.x = root.x
            scb.t = 0

            root_start = root.x

            saturation_region = nw.check_saturation_region(scb)

            if vlo != 0:
                vf_out[idx_vlo][idx] = nw.get_vf_output(ckt_op, scb) / vlo
            else:
                vf_out[idx_vlo][idx] = nw.get_vf_output(ckt_op, scb)

            if not root.success:
                print ("Failed to calculate small signal circuit")

            if not saturation_region:
                print ("Not in saturation region " + str(vlo) + " " + str(iref))

    return vf_out

def smer(iref_sweep, iref_op, vlo_sweep, vlo_bias, nw):

    # Instantiate large signal circuit
    ckt_op = Circuit()
    nw = gilbert_cell_switching_stage(ckt_op)
    ckt_op.init_circuit()

    # Define Vlo bias voltage
    nw.set_vlo_bias(bias=vlo_bias)

    scb = Scoreboard(num_edges=ckt_op.num_edges, num_sys_vars=ckt_op.num_sys_vars, degen_mtrx=ckt_op.degen_mtrx)

    # Vector to hold results
    vf_out_sml     = np.zeros([len(vlo_sweep), len(iref_sweep)])

    for idx_vlo, vlo in enumerate(vlo_sweep):

        # Define large signal operating point conditions
        nw.set_iref(iref_op)
        nw.set_vlo(vlo=vlo)

        # Calculate small signal
        root_start = np.ones(ckt_op.num_edges)
        root       = optimize.root(circuit_eqn, root_start, args=(ckt_op), tol=1e-9)
        if not root.success:
            print ("Failed to calculate small signal circuit")

        op = ckt_op.scb
        scb.x = root.x

        # Get Vrd1 and Vrd2 from large signal
        #  - Used to combine with small signal output values to create representative outputs from small signal
        v_rd1_base, v_rd2_base = nw.get_vf_output_base(ckt=ckt_op, scb=scb)

        # Create the small signal equivalent
        ckt_small = Circuit()
        nw.add_sml_ckt(ckt_sml=ckt_small, op=op)
        ckt_small.init_circuit()

        # Vlo for small signal is always zero
        #  - Value of Vlo used during creation of small signal circuit and Vlo is not varied
        nw.set_vlo_sml(0)

        # Solve for all values of iref
        root_start = np.ones(ckt_small.num_edges)
        for idx, iref in enumerate(iref_sweep):
            nw.set_iref_sml(iref)

            # Operating Point
            root       = optimize.root(circuit_eqn, root_start, args=(ckt_small), tol=1e-9)

            scb.x = root.x
            scb.t = 0

            root_start = root.x

            if vlo != 0:
                vf_out_sml[idx_vlo][idx] = nw.get_vf_output_sml(ckt_small, scb, v_rd1_base, v_rd2_base) / vlo
            else:
                vf_out_sml[idx_vlo][idx] = nw.get_vf_output_sml(ckt_small, scb, v_rd1_base, v_rd2_base)

            if not root.success:
                print ("Failed to calculate small signal circuit")

    return vf_out_sml

def main():

    # Define Circuit
    ckt_op = Circuit()
    nw = gilbert_cell_switching_stage(ckt_op)

    ckt_op.init_circuit()

    # Define Sweep parameters
    vlo_bias   = 1.5
    vlo_diff   = 0.3
    vlo_sweep  = np.linspace(vlo_diff/2, -vlo_diff/2, 20)

    iref_mid   = 150e-6
    iref_diff  = 50e-6

    # Small Signal
    # ------------

    iref_sweep_sml = np.linspace(-iref_diff, iref_diff, 100)

    vf_out_sml = smer(iref_sweep_sml, iref_mid, vlo_sweep, vlo_bias, nw)

    fig, ax1 = plt.subplots()
    #
    for idx_vlo, vlo in enumerate(vlo_sweep):
        ax1.plot(iref_sweep_sml, vf_out_sml[idx_vlo])
    plt.show(block=False)

    # Large Signal
    # ------------

    iref_sweep_lrg = np.linspace(iref_mid - iref_diff, iref_mid + iref_diff, 100)

    vf_out = sweep_lrg(iref_sweep_lrg,  vlo_sweep, vlo_bias, nw, ckt_op)

    fig, ax1 = plt.subplots()
    #
    for idx_vlo, vlo in enumerate(vlo_sweep):
        ax1.plot(iref_sweep_lrg, vf_out[idx_vlo])
    plt.show()

if __name__ == "__main__":

    main()
