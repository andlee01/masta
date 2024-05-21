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

def main():

    # Operating Point
    # ---------------
    ckt_op = Circuit()
    nw = gilbert_cell_switching_stage_dual(ckt_op)

    ckt_op.init_circuit()

    scb = Scoreboard(num_edges=ckt_op.num_edges, num_sys_vars=ckt_op.num_sys_vars, degen_mtrx=ckt_op.degen_mtrx)

    vlo_bias = 1.5
    vlo_diff = 0.3
    iref_min   = 100e-6
    iref_delta = 100e-6
    iref_max   = iref_min + iref_delta

    iref_sweep = np.linspace(0, iref_delta, 100)

    vlo_sweep  = np.linspace(-vlo_diff/2, vlo_diff/2, 20)

    vf_out     = np.zeros([len(vlo_sweep), len(iref_sweep)])

    nw.set_vlo_bias(bias=vlo_bias)

    for idx_vlo, vlo in enumerate(vlo_sweep):
        nw.set_vlo(vlo)

        root_start = np.ones(ckt_op.num_edges)
        for idx, iref in enumerate(iref_sweep):

            iref_n_val = 200e-6 - iref
            iref_p_val = 100e-6 + iref

            nw.set_iref(iref_n_val, iref_p_val)

            # Operating Point
            root       = optimize.root(circuit_eqn, root_start, args=(ckt_op), tol=1e-9)

            scb.x = root.x
            scb.t = 0

            root_start = root.x

            saturation_region = nw.check_saturation_region(scb)

            ref_vout = (vlo - vlo_bias)

            if vlo != 0:
                vf_out[idx_vlo][idx] = nw.get_vf_output(ckt_op, scb) / vlo
            else:
                vf_out[idx_vlo][idx] = nw.get_vf_output(ckt_op, scb)

            if not root.success:
                print ("Failed to calculate small signal circuit")

            if not saturation_region:
                print ("Not in saturation region " + str(vlo) + " " + str(iref))

    fig, ax1 = plt.subplots()
    #
    for idx_vlo, vlo in enumerate(vlo_sweep):
        ax1.plot(iref_sweep, vf_out[idx_vlo])
    #ax1.plot(T, yout[1][0])
    #ax2.plot(tr, yr)
    plt.show()


if __name__ == "__main__":

    main()
