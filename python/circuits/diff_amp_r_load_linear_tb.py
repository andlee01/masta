import sys, getopt, re, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

sys.path.append("../networks")
from diff_amp_pmos_load_ntwrk import *

from scipy import optimize
from scipy.integrate import ode

import numpy as np

import itertools

import matplotlib.pyplot as plt

import math

import control as ct

def circuit_eqn(x, sys, ckt, t):

    # current Im
    I = ckt.get_im(x=x, sys=sys, t=t)

    # voltage Vm
    V = ckt.get_vm(x=x, sys=sys, t=t)

    # get dependent currents
    I = ckt.get_im_dep()

    # Copy matrices
    qf_num = ckt.qf.copy()
    bf_num = ckt.bf.copy()

    eqn_qf = np.dot(qf_num, I)
    eqn_bf = np.dot(bf_num, V)

    return eqn_qf + eqn_bf

def saturation_test():

    mag  = 0.0
    bias = 3.0
    iref_val = 10e-6

    dc_params = {"mag": mag/2, "bias": bias}
    nw = diff_amp_pmos_load_ntwrk(dc_src=True, **dc_params)
    nw.set_iref(iref_val)
    nw.set_vlo_bias(bias=bias)

    vdiff_sweep  = np.arange(-0.1, 0.1, 0.001)
    v_vout       = np.zeros([len(vdiff_sweep)])

    root_start = np.ones(nw.ckt.num_edges)

    for idx, vdiff in enumerate(vdiff_sweep):

        nw.set_vlo(vdiff)
        root       = optimize.root(circuit_eqn, root_start, args=(0, nw.ckt, 0), tol=1e-6)
        root_start = root.x
        if not root.success:
            print (root.x)
            print (vdiff)
            sys.exit()

        if not nw.check_saturation_region(root.x):
            print ("not in saturation region at vdiff: " + str(vdiff))

        v_vout[idx] = nw.get_vf_output(root.x)

    return vdiff_sweep, v_vout


def main():
    vdiff_sweep, v_vout = saturation_test()

    fig, ax1 = plt.subplots()
    plt.grid()

    ax1.plot(vdiff_sweep, v_vout,       color='blue',  label="$V_{OUT}$")
    ax1.set_ylabel('v')
    ax1.legend()
    plt.show()

if __name__ == "__main__":

    main()
